#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <fluidsim.h>

#include <time.h>

#include <omp.h>

void key_value_desc(char* key, char* valuetype, char* desc) {
    printf("-%s\t\t<%s>\t\t%s\n", key, valuetype, desc);
}
void print_help() {
    key_value_desc("prefix", "string", "Attaches this prefix to output file names.");
    key_value_desc("grid-extents", "float", "Describes the extents of the square grid.");
    key_value_desc("grid-subdiv", "int", "Describes the finite element subdivisions of the square grid.");
    key_value_desc("parallel-count", "int", "Describes the parallelization of the grid. Must be a power of 2.");
    key_value_desc("steps", "int", "describes how far to run the simulation.");
    key_value_desc("num_particles", "int", "describes the number of particles in fluid system.");
    exit(0);
    
}
void print_convergence(FILE* f, double* conv_errs, long long* conv_iters, int nt) {
    for (int n = 0; n < nt; n++) {
        fprintf(f, "%lf,%ld\n", conv_errs[n], conv_iters[n]);
    }
}

void print_particles(FILE* f, particle_system* p, int nt) {
    for (int n = 0; n < p->num_particles; n++) {
        fprintf(f, "%d,%d,%lf,%lf\n", nt, p->ids[n], p->x[n], p->y[n]);
    }
}
void parse_args(int argc, char** argv, char* prefix, float* grid_extents, int* grid_subdiv, int* parallel_count, int* timesteps, int* num_particles) {
    for (int i = 0; i < argc; i++) {
       
        if (strcmp("-help", argv[i]) == 0) {
            
            print_help();
        } 
        char kvpair[64];
        if (strcmp(argv[i], "-prefix") == 0) {
            memcpy(prefix, argv[i+1], 64);
        }
        if (strcmp(argv[i], "-grid-extents") == 0) {
            *grid_extents = atof(argv[i+1]);   
        }
        if (strcmp(argv[i], "-grid-subdiv") == 0) {
            *grid_subdiv = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-parallel-count") == 0) {
            *parallel_count = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-steps") == 0) {
            *timesteps = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-num_particles") == 0) {
            *num_particles = atoi(argv[i+1]);
        }

    }
}

int main(int argc, char** argv) {

    char prefix[64]; float extent = 32.0; int parallel_count = 1, subdiv = 64, steps = 100, num_particles = 10;
    parse_args(argc, argv, prefix, &extent, &subdiv, &parallel_count, &steps, &num_particles);

    if (parallel_count == 1) {
        char fname[64], fname_conv[64], fname_part[64];

        sprintf(fname, "%s.csv", prefix);
        sprintf(fname_conv, "%s_convergence.csv", prefix);
        sprintf(fname_part, "%s_particles.csv", prefix);
#ifndef TIMING
        FILE* fp = fopen(fname, "w");
        FILE* fp_p = fopen(fname_part, "w");
#endif
        FILE* fpconv = fopen(fname_conv, "w");

        int T = steps;
        double* conv_errs = malloc(sizeof(double) * T);
        long long* conv_iters = malloc(sizeof(long long) * T);


        float X[2] = {0.0, extent}, Y[2] = {0.0, extent};
        
        fluid_grid* fl1 = fluid_grid_new(X, Y, subdiv, subdiv);
        double* bc1 = (double*)malloc(sizeof(double) * fl1->_ny * fl1->_nx);
        for (int y = 0; y < fl1->_ny; y++) {
            for (int x = 0; x < fl1->_nx; x++) {
                int I = y * fl1->_nx + x;
                bc1[I] = 1.0;
                if (y == 0 || y == fl1->_ny - 1)
                    bc1[I] = 0.0;
                if (x == 0 || x == fl1->_nx - 1)
                    bc1[I] = 0.0; 
            }
        }
        set_properties(fl1, 1.0, 0.3, 0.0001, bc1);

        particle_system* part = particle_system_new(extent, extent, extent - 0.1, 0.1, num_particles, 100, 1.0, 0.2, 2, 0.0001);

        int nt = 0;

        clock_t c = clock();

        while (nt < T) {
            enforce_bc(fl1);

            ld_cavity(fl1, 20.0, 0.0, 0);

            exchange_momentum_f2p(fl1, part);

            contain_particles_single(part);

            exchange_force_p2f(fl1, part);

            int nit = 0;
            double err = 1.0;
            
            poisson_rhs(fl1);

            while (err > 0.1 && nit < 200) {
                
                err = pressure_poisson_single(fl1, 1, 1, 1, 1);
                poisson_swap(fl1);
                nit++;
            }
            conv_errs[nt] = err;
            conv_iters[nt] = nit;

            advance_fluid(fl1);

            accelerate_fluid(fl1);

    #ifndef TIMING
            print_fluid_grid(fp, fl1, nt);
            if (num_particles > 0)
                print_particles(fp_p, part, nt);
    #endif

            printf("%d\n", nt);
            nt++;
        }

        c = (double)(clock() - c) / (double)(CLOCKS_PER_SEC / 1000.0);

        printf("Total time taken for %d steps:\t%lf ms", nt - 1, (double)c);

        print_convergence(fpconv, conv_errs, conv_iters, nt);
    }
    else {
        float fine_extent = extent / parallel_count;
        int fine_subdiv = subdiv / parallel_count;

        int nparallel_count = parallel_count + 2;

        FILE** fx = malloc(sizeof(FILE*) * nparallel_count * nparallel_count);
        memset(fx, 0, sizeof(FILE*) * nparallel_count * nparallel_count);

#ifndef TIMING
        for (int y = 1; y < nparallel_count - 1; y++) {
            for (int x = 1; x < nparallel_count - 1; x++) {
                char fname[64];
                sprintf(fname, "%s_%d_%d.csv", prefix, x, y);
                fx[y * (nparallel_count) + x] = fopen(fname, "w");
                memset(fname, 0, 64);
            }

        }
#endif
        char fname_coarse[64], fname_conv[64], fname_part[64];
        sprintf(fname_coarse, "%s_coarse.csv", prefix);
        sprintf(fname_conv, "%s_convergence.csv", prefix);
        sprintf(fname_part, "%s_particles.csv", prefix);
        
        FILE* convergence_graph = fopen(fname_conv, "w");
        FILE* fx_coarse = fopen(fname_coarse, "w");
        FILE* fx_p = fopen(fname_part, "w");

        int T = steps;
        double* conv_errs = malloc(sizeof(double) * T);
        long long* conv_iters = malloc(sizeof(long long) * T);

        float Xe[2] = {0.0, fine_extent}, Ye[2] = {0.0, fine_extent};

        fluid_grid** stacked_grid = malloc(sizeof(fluid_grid*) * nparallel_count * nparallel_count);
        memset(stacked_grid, 0, sizeof(fluid_grid*) * nparallel_count * nparallel_count);

        int* stacked_grid_lrtb = malloc(sizeof(int) * nparallel_count * nparallel_count * 4);

        for (int Y = 1; Y < nparallel_count - 1; Y++) {
            for (int X = 1; X < nparallel_count - 1; X++) {

                stacked_grid[Y * nparallel_count + X] = fluid_grid_new(Xe, Ye, fine_subdiv, fine_subdiv);
                fluid_grid* fl = stacked_grid[Y * nparallel_count + X];
                int* lrtb = &stacked_grid_lrtb[Y * nparallel_count + X * 4];
                double* bc = (double*)malloc(sizeof(double) * fl->_ny * fl->_nx);
                int top = 0, bottom = 0, left = 0, right = 0;
                if (X - 1 == 0) 
                    left = 1;
                if (X + 1 == nparallel_count - 1)
                    right = 1;
                if (Y - 1 == 0)
                    top = 1;
                if (Y + 1 == nparallel_count - 1)
                    bottom = 1;
                
                lrtb[0] = left;
                lrtb[1] = right;
                lrtb[2] = top;
                lrtb[3] = bottom;

                printf("%d %d: %d %d %d %d\n", X, Y, left, right, top, bottom);

                for (int y = 0; y < fl->_ny; y++) {
                    for (int x = 0; x < fl->_nx; x++) {
                        int I = y * fl->_nx + x;
                        bc[I] = 1.0;
                        if (top && y == 0)
                            bc[I] = 0.0;
                        if (left && x == 0)
                            bc[I] = 0.0; 
                        if (bottom && y == fl->_ny - 1)
                            bc[I] = 0.0; 
                        if (right && x == fl->_nx - 1)
                            bc[I] = 0.0; 
                    }
                }
                set_properties(fl, 1.0, 0.3, 0.0001, bc);
            }

        }
        

        Xe[0] = 0.0;
        Ye[0] = 0.0;
        Xe[1] = extent;
        Ye[1] = extent;

        fluid_grid* flc = fluid_grid_new(Xe, Ye, nparallel_count - 1, nparallel_count - 1);
        double* bcc = (double*)malloc(sizeof(double) * flc->_ny * flc->_nx);
        for (int y = 0; y < flc->_ny; y++) {
            for (int x = 0; x < flc->_nx; x++) {
                int I = y * flc->_nx + x;
                bcc[I] = 1.0;
                if (y == 0 || y == flc->_ny - 1)
                    bcc[I] = 0.0;
                if (x == 0 || x == flc->_nx - 1)
                    bcc[I] = 0.0; 
            }
        }
        set_properties(flc, 1.0, 0.3, 0.0001, bcc);

        particle_system** stacked_p_grid = malloc(sizeof(particle_system*) * nparallel_count * nparallel_count);
        memset(stacked_p_grid, 0, sizeof(particle_system*) * nparallel_count * nparallel_count);

        if (num_particles > 0) {
            for (int y = 1; y < nparallel_count - 1; y++) {
                for (int x = 1; x < nparallel_count - 1; x++) {
                    int nump = num_particles;
                    if (y != 1 || x != 1)
                        nump = 0;
                    stacked_p_grid[y * nparallel_count + x] = particle_system_new(fine_extent, fine_extent, fine_extent - 0.1, 0.1, nump, 100, 0.2, 0.2, 2, 0.0001);
                }
            }
        }
        
        omp_set_dynamic(0);  
        omp_set_num_threads(nparallel_count * nparallel_count);

        int nt = 0;
        float* errp = malloc(sizeof(float) * nparallel_count * nparallel_count);

        clock_t c = clock();

        while (nt < T) {

            #pragma omp parallel
            {
                int rownum = omp_get_thread_num() / nparallel_count;
                int rowofs = omp_get_thread_num() % nparallel_count;
                fluid_grid* G = stacked_grid[rownum * nparallel_count + rowofs];
                if (G)
                    enforce_bc(G);
            }

            // #pragma omp parallel for
            // for (int y = 1; y < nparallel_count - 1; y++) {
            //     for (int x = 1; x < nparallel_count - 1; x++) {
                    
            //         enforce_bc(stacked_grid[y * nparallel_count + x]);
            //     }
            // }

            for (int x = 1; x < nparallel_count - 1; x++) {
                ld_cavity(stacked_grid[nparallel_count + x], 20.0, 0.0, 0);
            }

            if (num_particles > 0) {

                #pragma omp parallel for
                for (int y = 1; y < nparallel_count - 1; y++) {
                    for (int x = 1; x < nparallel_count - 1; x++) {
                        exchange_momentum_f2p(stacked_grid[y * nparallel_count + x], stacked_p_grid[y * nparallel_count + x]);
                    }
                }

                exchange_inter_particles(stacked_p_grid, nparallel_count, nparallel_count);

                #pragma omp parallel for
                for (int y = 1; y < nparallel_count - 1; y++) {
                    for (int x = 1; x < nparallel_count - 1; x++) {
                        exchange_force_p2f(stacked_grid[y * nparallel_count + x], stacked_p_grid[y * nparallel_count + x]);
                    }
                }
            }   

            
            // ld_cavity(fl2, 0.0, 20.0, 0);
            
            double err = 1.0;
            int nit = 0;

            #pragma omp parallel
            {
                int rownum = omp_get_thread_num() / nparallel_count;
                int rowofs = omp_get_thread_num() % nparallel_count;
                fluid_grid* G = stacked_grid[rownum * nparallel_count + rowofs];
                if (G)
                poisson_rhs(G);
            }

            // #pragma omp parallel for
            // for (int y = 1; y < nparallel_count - 1; y++) {
            //     for (int x = 1; x < nparallel_count - 1; x++) {
            //     poisson_rhs(stacked_grid[y * nparallel_count + x]);
            //     }
            // }


            while (err > 0.1 && nit < 200) {
        
                
                for (int i = 0; i < nparallel_count * nparallel_count; i++)
                    errp[i] = 0.0;

                #pragma omp parallel
                {
                    int rownum = omp_get_thread_num() / nparallel_count;
                    int rowofs = omp_get_thread_num() % nparallel_count;
                    int I = rownum * nparallel_count + rowofs;
                    fluid_grid* G = stacked_grid[rownum * nparallel_count + rowofs];
                    int* lrtb = &stacked_grid_lrtb[rownum * nparallel_count + 4 * rowofs];
                    if (G) {
                        errp[I] = 
                            pressure_poisson_single(
                                G, lrtb[0], lrtb[1], lrtb[2], lrtb[3]);
                    }
                        
                }
                
                // #pragma omp parallel for
                // for (int y = 1; y < nparallel_count - 1; y++) {
                //     for (int x = 1; x < nparallel_count - 1; x++) {
                //         errp[y * nparallel_count + x] = 
                //             pressure_poisson_single(
                //                 stacked_grid[y * nparallel_count + x],
                //                 stacked_grid_lrtb[y * nparallel_count + 4 * x + 0],
                //                 stacked_grid_lrtb[y * nparallel_count + 4 * x + 1],
                //                 stacked_grid_lrtb[y * nparallel_count + 4 * x + 2],
                //                 stacked_grid_lrtb[y * nparallel_count + 4 * x + 3]);
                //     }
                // }
                
                err = 0.0;
                for (int i = 0; i < nparallel_count * nparallel_count; i++) {
                    err += errp[i];
                }
                err /= (double)(nparallel_count * nparallel_count);
                
                match_interface_p(stacked_grid, nparallel_count, nparallel_count);

                fl_restrict_p(flc->_p, stacked_grid, nparallel_count, nparallel_count);
                
                poisson_rhs(flc);
                double err3 = pressure_poisson_single(flc, 1, 1, 1, 1);

                diff(flc->_p_diff, flc->_p_next, flc->_p, flc->_nx, flc->_ny);

                #pragma omp parallel
                {
                    int rownum = omp_get_thread_num() / nparallel_count;
                    int rowofs = omp_get_thread_num() % nparallel_count;
                    int I = rownum * nparallel_count + rowofs;
                    fluid_grid* G = stacked_grid[I];
                    
                    if (G) {
                        fl_prolongate_p(flc->_p_diff[I], G);
                    }
                        
                }

                // #pragma omp parallel for
                // for (int y = 1; y < nparallel_count - 1; y++) {
                //     for (int x = 1; x < nparallel_count - 1; x++) {
                //         fl_prolongate_p(flc->_p_diff[y * nparallel_count + x], stacked_grid[y * nparallel_count + x]);
                //     }
                // }

                #pragma omp parallel
                {
                    int rownum = omp_get_thread_num() / nparallel_count;
                    int rowofs = omp_get_thread_num() % nparallel_count;
                    int I = rownum * nparallel_count + rowofs;
                    fluid_grid* G = stacked_grid[I];
                    
                    if (G) {
                        poisson_swap(G);
                    }
                        
                }

                // #pragma omp parallel for
                // for (int y = 1; y < nparallel_count - 1; y++) {
                //     for (int x = 1; x < nparallel_count - 1; x++) {
                //         poisson_swap(stacked_grid[y * nparallel_count + x]);
                //     }
                // }            
                nit++;
            }

            conv_iters[nt] = nit;
            conv_errs[nt] = err;

            nit = 0;
            while (nit < 20) {

                #pragma omp parallel
                {
                    int rownum = omp_get_thread_num() / nparallel_count;
                    int rowofs = omp_get_thread_num() % nparallel_count;
                    int I = rownum * nparallel_count + rowofs;
                    fluid_grid* G = stacked_grid[I];
                    
                    if (G) {
                        advance_fluid(G);
                    }
                        
                }

                // #pragma omp parallel for
                // for (int y = 1; y < nparallel_count - 1; y++) {
                //     for (int x = 1; x < nparallel_count - 1; x++) {
                //         advance_fluid(stacked_grid[y * nparallel_count + x]);
                //     }
                // }
                
                match_interface_uv(stacked_grid, nparallel_count, nparallel_count);

                fl_restrict_u(flc->_u, stacked_grid, nparallel_count, nparallel_count);
                fl_restrict_v(flc->_v, stacked_grid, nparallel_count, nparallel_count);
                fl_restrict_F(flc->Fx, flc->Fy, stacked_grid, nparallel_count, nparallel_count);

                enforce_bc(flc);
                advance_fluid(flc);

                diff(flc->_u_diff, flc->_u_next, flc->_u, nparallel_count, nparallel_count);
                diff(flc->_v_diff, flc->_v_next, flc->_v, nparallel_count, nparallel_count);

                #pragma omp parallel
                {
                    int rownum = omp_get_thread_num() / nparallel_count;
                    int rowofs = omp_get_thread_num() % nparallel_count;
                    int I = rownum * nparallel_count + rowofs;
                    fluid_grid* G = stacked_grid[I];
                    
                    if (G) {
                        fl_prolongate_u(flc->_u_diff[I], G);
                        fl_prolongate_v(flc->_v_diff[I], G);
                    }
                        
                }
                // #pragma omp parallel for
                // for (int y = 1; y < nparallel_count - 1; y++) {
                //     for (int x = 1; x < nparallel_count - 1; x++) {
                //         fl_prolongate_u(flc->_u_diff[y * nparallel_count + x], stacked_grid[y * nparallel_count + x]);
                //         fl_prolongate_v(flc->_v_diff[y * nparallel_count + x], stacked_grid[y * nparallel_count + x]);

                //     }
                // }

                nit++;
            }

            #pragma omp parallel
            {
                int rownum = omp_get_thread_num() / nparallel_count;
                int rowofs = omp_get_thread_num() % nparallel_count;
                int I = rownum * nparallel_count + rowofs;
                fluid_grid* G = stacked_grid[I];
                
                if (G) {
                    accelerate_fluid(G);
                }
                    
            }

            // #pragma omp parallel for
            // for (int y = 1; y < nparallel_count - 1; y++) {
            //     for (int x = 1; x < nparallel_count - 1; x++) {
            //         accelerate_fluid(stacked_grid[y * nparallel_count + x]);
            //     }
            // }

    #ifndef TIMING
            for (int y = 1; y < nparallel_count - 1; y++) {
                for (int x = 1; x < nparallel_count - 1; x++) {
                    print_fluid_grid(fx[y * nparallel_count + x], stacked_grid[y * nparallel_count + x], nt);
                    if (num_particles > 0)
                        print_particles(fx_p, stacked_p_grid[y * nparallel_count + x], nt);
                }
            }
            print_fluid_grid(fx_coarse, flc, nt);
    #endif
            printf("%d\n", nt);
            nt++;

        }

        c = (double)(clock() - c) / (double)(CLOCKS_PER_SEC / 1000.0);

        printf("Total time taken for %d steps:\t%lf ms", nt - 1, (double)c);

        print_convergence(convergence_graph, conv_errs, conv_iters, nt);
    }
}
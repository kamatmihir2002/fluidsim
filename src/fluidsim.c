#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fluidsim.h>

double diffx(double* a, int i, double dx) {
    return (a[i + 1] - a[i - 1]) / (2 * dx);
}

double diffy(double* a, int i, double dy, int nx) {
    return (a[i + nx] - a[i - nx]) / (2 * dy);
}

double diff2x(double* a, int i, double dx2) {
    return (a[i + 1] - 2 * a[i] + a[i - 1]) / (dx2);
}

double diff2y(double* a, int i, double dy2, int nx) {
    return (a[i + nx] - 2 * a[i] + a[i - nx]) / (dy2);
}

// particle_system particle_system_new(int num_particles, double* initial_x, double* initial_y) {

//     particle_system sys;
//     sys.num_particles = num_particles;
//     sys.ax = malloc(sizeof(double) * num_particles);
//     sys.ay = malloc(sizeof(double) * num_particles);
//     sys.x = malloc(sizeof(double) * num_particles);
//     sys.y = malloc(sizeof(double) * num_particles);
//     sys.u = malloc(sizeof(double) * num_particles);
//     sys.v = malloc(sizeof(double) * num_particles);
//     for (int i = 0; i < num_particles; i++) {
//         sys.x[i] = initial_x[i];
//         sys.y[i] = initial_y[i];
//         sys.ax[i] = 0.0;
//         sys.ay[i] = 0.0;
//         sys.u[i] = 0.0;
//         sys.v[i] = 0.0;
//     }
//     return sys;
// }


// void integrate_particle_system(particle_system** pgrid, fluid_grid** fgrid, int Nx, int Ny)  {
//     for (int y = 1; y < Ny - 1; y++) {
//         for (int x = 1; x < Nx - 1; x++) {
//             particle_system* p = pgrid[y * Nx + x];
//             particle_system* p_above = pgrid[(y - 1) * Nx + x];
//             particle_system* p_below = pgrid[(y + 1) * Nx + x];
//             fluid_grid* f = fgrid[y * Nx + x];
//             for (int i = 0; i < p->num_particles; i++) {
//                 size_t cellx = p->x[i] / f->_dx;
//                 size_t celly = p->y[i] / f->_dy;
//                 size_t fI = celly * f->_nx + cellx;
//                 if (cellx == 0 || cellx == f->_nx - 1) {
//                     //reflect 
//                 }
//                 else if (celly == 0) {
//                     // put_particle_in_system(p, p_above, i);
//                 }
//                 else if (celly == f->_ny - 1) {
//                     // put_particle_in_system(p, p_below, i);
//                 }
//                 else {
//                     p->ax[i] = p->drag[i] * (f->_u[fI] - p->u[i]) - (diffx(f->_p, fI, f->_dx) + diffy(f->_p, fI, f->_dy, f->_nx)) * 0.5 * p->inv_mass[i];
//                     p->ay[i] = p->drag[i] * (f->_v[fI] - p->v[i]) - (diffx(f->_p, fI, f->_dx) + diffy(f->_p, fI, f->_dy, f->_nx)) * 0.5 * p->inv_mass[i];
//                     p->u[i] += p->ax[i] * p->dt;
//                     p->v[i] += p->ay[i] * p->dt;
//                 }
//                 p->x[i] += p->u[i] * p->dt;
//                 p->y[i] += p->v[i] * p->dt;
//             }
//         }
//     }
    
// }


fluid_grid* fluid_grid_new(float x_extents[2], float y_extents[2], size_t x_subdiv, size_t y_subdiv) {
    fluid_grid* f = malloc(sizeof(fluid_grid));
    f->x_domain_extents[0] = x_extents[0];
    f->x_domain_extents[1] = x_extents[1];
    f->y_domain_extents[0] = y_extents[0];
    f->y_domain_extents[1] = y_extents[1];
    f->_nx = x_subdiv + 1;
    f->_ny = y_subdiv + 1;

    int nx = f->_nx, ny = f->_ny;

    f->_dx = (x_extents[1] - x_extents[0]) / x_subdiv;
    f->_dy = (y_extents[1] - y_extents[0]) / y_subdiv;
    f->Fx = (double*)malloc(sizeof(double) * nx * ny);
    
    f->Fy = (double*)malloc(sizeof(double) * nx * ny);

    f->b = (double*)malloc(sizeof(double) * nx * ny);

    f->_u = (double*)malloc(sizeof(double) * nx * ny);
    f->accel_x = (double*)malloc(sizeof(double) * nx * ny);
    f->accel_y = (double*)malloc(sizeof(double) * nx * ny);
    f->_u_next = (double*)malloc(sizeof(double) * nx * ny);
    f->_v = (double*)malloc(sizeof(double) * nx * ny);
    f->_v_next = (double*)malloc(sizeof(double) * nx * ny);
    f->_p = (double*)malloc(sizeof(double) * nx * ny);
    f->_p_next = (double*)malloc(sizeof(double) * nx * ny);
    f->_p_diff = (double*)malloc(sizeof(double) * ny * nx);
    f->_u_diff = (double*)malloc(sizeof(double) * ny * nx);
    f->_v_diff = (double*)malloc(sizeof(double) * ny * nx);
    
    f->_p_prol = (double*)malloc(sizeof(double) * ny * nx);
    f->_u_prol = (double*)malloc(sizeof(double) * ny * nx);
    f->_v_prol = (double*)malloc(sizeof(double) * ny * nx);
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            f->Fy[I] = 0.0;
            f->_u[I] = 0.0;
            f->b[I] = 0.0;
            f->_u_next[I] = 0.0;
            f->_v[I] = 0.0;
            f->_v_next[I] = 0.0;
            f->_p[I] = 0.0;
            f->_p_next[I] = 0.0;
            f->accel_x[I] = 0.0;
            f->accel_y[I] = 0.0;
            f->_p_diff[I] = 0.0;
            f->_u_diff[I] = 0.0;
            f->_v_diff[I] = 0.0;
            f->_p_prol[I] = 0.0;
            f->_u_prol[I] = 0.0;
            f->_v_prol[I] = 0.0;
        }
    }
    return f;
}

void enforce_bc(fluid_grid* f) {

    double *u = f->_u, *v = f->_v;
    
    size_t nx = f->_nx, ny = f->_ny;

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            u[I] *= f->bc[I];
            v[I] *= f->bc[I];
        }
    }
}

void enforce_neumann_bc(fluid_grid* f, float delx, float dely, int include_left, int include_right, int include_top, int include_bottom) {
    size_t nx = f->_nx, ny = f->_ny;
    double* p = f->_p_next;
    double* po = f->_p;

    for (int y = 0; y < ny; y++) {
        int I = y * ny;
        int Ie = y * ny + (nx - 1);
        if (include_left) {
            po[I] = po[I + 1] + delx * f->_dx;
        }
        if (include_right) {
            po[Ie] = po[Ie - 1] + delx * f->_dx;
        }
        if (include_left){
            p[I] = p[I + 1] + delx * f->_dx;
        }
        if (include_right) {
            p[Ie] = p[Ie - 1] + delx * f->_dx;
        }
        
    }

    for (int x = 0; x < nx; x++) {
        int I = x;
        int Ie = (ny - 1) * nx;
        if (include_top)
            po[I] = po[I + nx] + dely * f->_dy;
        if (include_bottom)
            po[Ie] = po[Ie - nx] + dely * f->_dy;
        if (include_top)
            p[I] = p[I + nx] + dely * f->_dy;
        if (include_bottom)
            p[Ie] = p[Ie - nx] + dely * f->_dy;
    }

    
}

void ld_cavity(fluid_grid* f, float velx, float vely, int col) {

    double *u = f->_u, *v = f->_v;
    size_t nx = f->_nx, ny = f->_ny;
    for (int x = 0; x < nx; x++) {
        u[nx * col + x] = velx;
        v[nx * col + x] = vely;
    
    }

}


void set_properties(fluid_grid* f, float rho, float nu, float dt, double* bc) {
    f->rho = rho;
    f->nu = nu;
    f->_dt = dt;
    f->bc = bc;
}


void build_up_b(fluid_grid* f) {
    
    double *u = f->_u,
            *v = f->_v,
            dx = f->_dx, dy = f->_dy,
            rho = f->rho, dt = f->_dt,
            *b = f->b;
    size_t nx = f->_nx, ny = f->_ny;

    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            int I = y * nx + x;
            b[I] = rho * (1.0 / dt * (diffx(u, I, dx) + diffy(v, I, dy, nx)) 
            - pow(diffx(u, I, dx), 2.0) 
            - pow(diffy(v, I, dy, nx), 2.0) 
            - 2 * diffy(u, I, dy, nx) * diffx(v, I, dx)
            );

        }
    }
}

double pressure_poisson_single(fluid_grid* f, int include_left, int include_right, int include_top, int include_bottom) {
    double dx = f->_dx,
            dy = f->_dy,
            *p_old = f->_p,
            *p = f->_p_next,
            *b = f->b;

    size_t ny = f->_ny,
            nx = f->_nx;

    double dy2 = pow(dy, 2.0);
    double dx2 = pow(dx, 2.0);

    f->p_avg = 0.0;
    double E = 0.0;
    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            int I = y * nx + x;
            p[I] = 
                (((p_old[I + 1] + p_old[I - 1]) * dy2 + 
                (p_old[I + nx] + p_old[I - nx]) * dx2) / 
                (2.0 * (dx2 + dy2)) - 
                dx2 * dy2 / (2.0 * (dx2 + dy2)) * b[I]);
            f->p_avg += p[I];
            E += fabs(p[I] - p_old[I]);
        }
    }
    
    
    f->p_avg /= (ny * nx);
    
    
    enforce_neumann_bc(f, 0.0, 0.0, include_left, include_right, include_top, include_bottom);
    return E / (ny * nx);
        
        
}


void pressure_poisson(fluid_grid* f, double* b, int include_top, int include_bottom) {
    double dx = f->_dx,
            dy = f->_dy,
            *p_old = f->_p,
            *p = f->_p_next;

    size_t ny = f->_ny,
            nx = f->_nx;

    double err = 1.0;
    double dy2 = pow(dy, 2.0);
    double dx2 = pow(dx, 2.0);
    int iter = 0;
    while (err >= 0.02 && iter <= 50) {
        
        memcpy(p_old, p, sizeof(double) * ny * nx);

        err = 0.0;
        f->p_avg = 0.0;
        for (int y = 1; y < ny - 1; y++) {
            for (int x = 1; x < nx - 1; x++) {
                int I = y * nx + x;
                p[I] = 
                    (((p_old[I + 1] + p_old[I - 1]) * dy2 + 
                    (p_old[I + nx] + p_old[I - nx]) * dx2) / 
                    (2.0 * (dx2 + dy2)) - 
                    dx2 * dy2 / (2.0 * (dx2 + dy2)) * b[I]);
                f->p_avg += p[I];
                err += fabs(p_old[I] - p[I]);
            }
        }
        f->p_avg /= (ny * nx);
        
        
        enforce_neumann_bc(f, 0.0, 0.0, 1, 1, include_top, include_bottom);

        err /= (double)(ny * nx);
        
        iter++;

    }
}

void poisson_swap(fluid_grid* f) {

    double *p_old = f->_p,
    *p = f->_p_next;

    size_t ny = f->_ny,
            nx = f->_nx;

    memcpy(p_old, p, sizeof(double) * ny * nx);

}

void accelerate_fluid(fluid_grid* f) {
    double *un = f->_u_next,
    *u = f->_u,
    *vn = f->_v_next,
    *v = f->_v;
    

    size_t ny = f->_ny, nx = f->_nx;

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            u[I] = un[I];
            v[I] = vn[I];
            
            
        }
    }

    

}

void advance_fluid(fluid_grid* f) {
    double *un = f->_u_next,
        *accel_x = f->accel_x,
        *accel_y = f->accel_y,
        *u = f->_u,
        *vn = f->_v_next,
        *v = f->_v,
        *p = f->_p,
        dx = f->_dx,
        dy = f->_dy,
        nu = f->nu,
        *Fx = f->Fx,
        *Fy = f->Fy,
        dt = f->_dt,
        rho = f->rho;

    size_t ny = f->_ny, nx = f->_nx;

    double dx2 = pow(dx, 2.0);
    double dy2 = pow(dy, 2.0);

    f->u_avg = 0.0;
    f->v_avg = 0.0;

    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            int I = y * nx + x;
            double udx = diffx(u, I, dx);
            double vdx = diffx(v, I, dx);
            double pdx = diffx(p, I, dx);

            double udy = diffy(u, I, dy, nx);
            double vdy = diffy(v, I, dy, nx);
            double pdy = diffy(p, I, dy, nx);

            double ustress = diff2x(u, I, dx2) + diff2y(u, I, dy2, nx);
            double vstress = diff2x(v, I, dx2) + diff2y(v, I, dy2, nx);

            accel_x[I] = -u[I] * udx - v[I] * udy - 1.0 / rho * pdx + nu / rho * ustress + Fx[I] / rho; 
            accel_y[I] = -u[I] * vdx - v[I] * vdy - 1.0 / rho * pdy + nu / rho * vstress + Fy[I] / rho; 

            un[I] = u[I] + accel_x[I] * dt;
            vn[I] = v[I] + accel_y[I] * dt;
            f->u_avg += un[I];
            f->v_avg += vn[I];
            f->Fx_avg += Fx[I];
            f->Fy_avg += Fy[I];
        }
    }

    for (int y = 0; y < ny; y++) {
        int I = y * nx;
        int Ie = y * nx + (nx - 1);
        un[I] = un[I + 1];
        un[Ie] = un[Ie - 1];

        vn[I] = vn[I + 1];
        vn[Ie] = vn[Ie - 1];
        
    }

    for (int x = 0; x < nx; x++) {
        int I = x;
        int Ie = (ny - 1) * nx + x;
        un[I] = un[I + nx];
        un[Ie] = un[Ie - nx];

        vn[I] = vn[I + nx];
        vn[Ie] = vn[Ie - nx];
    }


    f->u_avg /= (ny * nx);
    f->v_avg /= (ny * nx);
    f->Fx_avg /= (ny * nx);
    f->Fy_avg /= (ny * nx);
    // printf("%f ", f->u_avg);
    // printf("%f\n", f->v_avg);
}

void print_fluid_grid(FILE* f, fluid_grid* fl, size_t timestep) {
    size_t ny = fl->_ny, nx = fl->_nx;
    double *u = fl->_u, *v = fl->_v, *p = fl->_p;

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            if (x == 0)
                fprintf(f, "%d,x,%lf", timestep, u[I]);
            else
                fprintf(f, ",%lf", u[I]);
        }
        fprintf(f, "\n");
    }
    
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            if (x == 0)
                fprintf(f, "%d,y,%lf", timestep, v[I]);
            else
                fprintf(f, ",%lf", v[I]);
        }
        fprintf(f, "\n");
    }
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int I = y * nx + x;
            if (x == 0)
                fprintf(f, "%d,p,%lf", timestep, p[I]);
            else
                fprintf(f, ",%lf", p[I]);
        }
        fprintf(f, "\n");
    }
}

double fl_restrict_p(double* p, fluid_grid** grids, int Nx, int Ny) {
    for (int Y = 1; Y < Ny - 1; Y++) {
        for (int X = 1; X < Nx - 1; X++) {
            int I = Y * Nx + X;
            if (grids[I])
                p[I] = grids[I]->p_avg;  
        } 
    }
    
}

double fl_restrict_u(double* u, fluid_grid** grids, int Nx, int Ny) {
    for (int Y = 1; Y < Ny - 1; Y++) {
        for (int X = 1; X < Nx - 1; X++) {
            int I = Y * Nx + X;
            if (grids[I])
                u[I] = grids[I]->u_avg;
        } 
    }
    
}

double fl_restrict_v(double* v, fluid_grid** grids, int Nx, int Ny) {
    for (int Y = 1; Y < Ny - 1; Y++) {
        for (int X = 1; X < Nx - 1; X++) {
            int I = Y * Nx + X;
            if (grids[I])
                v[I] = grids[I]->v_avg;
        } 
    }
    
}

double fl_restrict_F(double* Fx, double* Fy, fluid_grid** grids, int Nx, int Ny) {
    for (int Y = 1; Y < Ny - 1; Y++) {
        for (int X = 1; X < Nx - 1; X++) {
            int I = Y * Nx + X;
            if (grids[I]) {
                Fx[I] = grids[I]->Fx_avg;
                Fy[I] = grids[I]->Fy_avg;
            }
        } 
    }
    
}


double fl_prolongate_p(double p, fluid_grid* grid) {
    fluid_grid* G = grid;
    
    for (int y = 0; y < G->_ny; y++) {
        for (int x = 0; x < G->_nx; x++) {
            int IF = y * G->_nx + x;
            G->_p_next[IF] += p;
            
        }
    }
}

double fl_prolongate_u(double u, fluid_grid* grid) {
    fluid_grid* G = grid;
    
    for (int y = 0; y < G->_ny; y++) {
        for (int x = 0; x < G->_nx; x++) {
            int IF = y * G->_nx + x;
            G->_u_next[IF] += u;
            
        }
    }

}

double fl_prolongate_v(double v, fluid_grid* grid) {
    fluid_grid* G = grid;
    
    for (int y = 0; y < G->_ny; y++) {
        for (int x = 0; x < G->_nx; x++) {
            int IF = y * G->_nx + x;
            G->_v_next[IF] += v;
            
        }
    }
}

void diff(double* result, double* A, double* B, int nx, int ny) {

    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            int I = y * nx + x;
            result[I] = A[I] - B[I];
        }
    }
}

void match_interface_p(fluid_grid** grid, int nx, int ny) {
    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            fluid_grid* F = grid[y * nx + x];
            fluid_grid* Fleft = grid[y * nx + x - 1];
            fluid_grid* Fright = grid[y * nx + x + 1];
            fluid_grid* Ftop = grid[(y - 1) * nx + x];
            fluid_grid* Fbottom = grid[(y + 1) * nx + x];
            if (F && Fleft) {
                double* Fleft_end = &Fleft->_p_next[Fleft->_nx - 1];
                double* F_start = &F->_p_next[0];
                double avg = 0.0;
                for (int Y = 0; Y < Fleft->_ny; Y++) {
                    avg = (F_start[Y * (F->_nx)] + Fleft_end[Y * Fleft->_nx]) * 0.5;
                    Fleft_end[Y * Fleft->_nx] = avg;
                    F_start[Y * F->_nx] = avg;
                }
            }
            if (F && Fright) {
                double* F_end = &F->_p_next[F->_nx - 1];
                double* Fright_start = &Fright->_p_next[0];
                double avg = 0.0;
                for (int Y = 0; Y < Fright->_ny; Y++) {
                    avg = (Fright_start[Y * Fright->_nx] + F_end[Y * F->_nx]) * 0.5;
                    F_end[Y * Fright->_nx] = avg;
                    Fright_start[Y * F->_nx] = avg;
                }

            }
            if (F && Ftop) {
                double* Ftop_end = &Ftop->_p_next[(Ftop->_ny - 1) * Ftop->_nx];
                double* F_start = &F->_p_next[0];
                double avg = 0.0;
                for (int X = 0; X < Ftop->_nx; X++) {
                    avg = (F_start[X] + Ftop_end[X]) * 0.5;
                    Ftop_end[X] = avg;
                    F_start[X] = avg;
                }
            }
            if (F && Fbottom) {
                double* F_end = &F->_p_next[(F->_ny - 1) * F->_nx];
                double* Fbottom_start = &Fbottom->_p_next[0];
                double avg = 0.0;
                for (int X = 0; X < Fbottom->_nx; X++) {
                    avg = (Fbottom_start[X] + F_end[X]) * 0.5;
                    F_end[X] = avg;
                    Fbottom_start[X] = avg;
                }

            }

        }
    }
}

void match_interface_uv(fluid_grid** grid, int nx, int ny) {
    for (int y = 1; y < ny - 1; y++) {
        for (int x = 1; x < nx - 1; x++) {
            fluid_grid* F = grid[y * nx + x];
            fluid_grid* Ftop = grid[(y - 1) * nx + x];
            fluid_grid* Fbottom = grid[(y + 1) * nx + x];
            fluid_grid* Fleft = grid[y * nx + x - 1];
            fluid_grid* Fright = grid[y * nx + x + 1];
            if (F && Fleft) {
                double* Fleft_end = &Fleft->_u_next[Fleft->_nx - 1];
                double* F_start = &F->_u_next[0];
                double avg = 0.0;
                for (int Y = 0; Y < Fleft->_ny; Y++) {
                    avg = (F_start[Y * F->_nx] + Fleft_end[Y * Fleft->_nx]) * 0.5;
                    Fleft_end[Y * Fleft->_nx] = avg;
                    F_start[Y * F->_nx] = avg;
                }
                Fleft_end = &Fleft->_v_next[Fleft->_nx - 1];
                F_start = &F->_v_next[0];
                avg = 0.0;
                for (int Y = 0; Y < Fleft->_ny; Y++) {
                    avg = (F_start[Y * F->_nx] + Fleft_end[Y * Fleft->_nx]) * 0.5;
                    Fleft_end[Y * Fleft->_nx] = avg;
                    F_start[Y * F->_nx] = avg;
                }
            }
            if (F && Fright) {
                double* F_end = &F->_u_next[F->_nx - 1];
                double* Fright_start = &Fright->_u_next[0];
                double avg = 0.0;
                for (int Y = 0; Y < Fright->_ny; Y++) {
                    avg = (Fright_start[Y * Fright->_nx] + F_end[Y * F->_nx]) * 0.5;
                    F_end[Y * Fright->_nx] = avg;
                    Fright_start[Y * F->_nx] = avg;
                }
                F_end = &F->_v_next[F->_nx - 1];
                Fright_start = &Fright->_v_next[0];
                avg = 0.0;
                for (int Y = 0; Y < Fright->_ny; Y++) {
                    avg = (Fright_start[Y * Fright->_nx] + F_end[Y * F->_nx]) * 0.5;
                    F_end[Y * Fright->_nx] = avg;
                    Fright_start[Y * F->_nx] = avg;
                }
            }
            if (F && Ftop) {
                double* Ftop_end = &Ftop->_u_next[(Ftop->_ny - 1) * Ftop->_nx];
                double* F_start = &F->_u_next[0];
                double avg = 0.0;
                for (int X = 0; X < Ftop->_nx; X++) {
                    avg = (F_start[X] + Ftop_end[X]) * 0.5;
                    Ftop_end[X] = avg;
                    F_start[X] = avg;
                }
                Ftop_end = &Ftop->_v_next[(Ftop->_ny - 1) * Ftop->_nx];
                F_start = &F->_v_next[0];
                avg = 0.0;
                for (int X = 0; X < Ftop->_nx; X++) {
                    avg = (F_start[X] + Ftop_end[X]) * 0.5;
                    Ftop_end[X] = avg;
                    F_start[X] = avg;
                }
            }
            if (F && Fbottom) {
                double* F_end = &F->_u_next[(F->_ny - 1) * F->_nx];
                double* Fbottom_start = &Fbottom->_u_next[0];
                double avg = 0.0;
                for (int X = 0; X < Fbottom->_nx; X++) {
                    avg = (Fbottom_start[X] + F_end[X]) * 0.5;
                    F_end[X] = avg;
                    Fbottom_start[X] = avg;
                }

                F_end = &F->_v_next[(F->_ny - 1) * F->_nx];
                Fbottom_start = &Fbottom->_v_next[0];
                avg = 0.0;
                for (int X = 0; X < Fbottom->_nx; X++) {
                    avg = (Fbottom_start[X] + F_end[X]) * 0.5;
                    F_end[X] = avg;
                    Fbottom_start[X] = avg;
                }

            }

        }
    }
}
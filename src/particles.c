#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fluidsim.h>

void push_particle(particle_system* p, double mass, double x, double y, double u, double v, double ax, double ay, double drag) {
    int top = p->num_particles;

    p->m_inv[top] = 1.0 / mass;
    p->x[top] = x;
    p->y[top] = y;
    p->u[top] = u;
    p->v[top] = v;
    p->ax[top] = ax;
    p->ay[top] = ay;
    p->ids[top] = top;
    p->num_particles++;
}

particle_system* particle_system_new(double x_extents, double y_extents, double init_x, double init_y, int num_particles, int capacity, double mass, double drag, int aoe, double dt) {
    srand(time(0));
    particle_system* p = malloc(sizeof(particle_system));
    p->particle_capacity = capacity;
    p->aoe = aoe;
    p->dt = dt;
    p->num_particles = num_particles;
    p->m_inv = malloc(sizeof(double) * capacity);
    p->drag = malloc(sizeof(double) * capacity);
    p->x = malloc(sizeof(double) * capacity);
    p->y = malloc(sizeof(double) * capacity);
    p->u = malloc(sizeof(double) * capacity);
    p->v = malloc(sizeof(double) * capacity);
    p->ax = malloc(sizeof(double) * capacity);
    p->ay = malloc(sizeof(double) * capacity);
    p->ids = malloc(sizeof(double) * capacity);
    p->nx = x_extents;
    p->ny = y_extents;
    for (int i = 0; i < capacity; i++) {
        p->ids[i] = i;
        float mul = 1.0;
        if (i >= num_particles)
            mul = 0.0;
        p->m_inv[i] =  mul * 1.0 / clamp((mass + 0.02 * (((float)rand() / (float)RAND_MAX) * 2.0 - 1.0)), 0.0, 2.0);
        p->drag[i] =  mul * clamp(drag + 0.05 * (((float)rand() / (float)RAND_MAX) * 2.0 - 1.0), 0.0, 1.0);
        p->x[i]  = mul *  clamp(init_x + 0.05 * (((float)rand() / (float)RAND_MAX) * 2.0 - 1.0), 0.0, p->nx);
        p->y[i]  = mul *  clamp(init_y + 0.05 * (((float)rand() / (float)RAND_MAX) * 2.0 - 1.0), 0.0, p->ny);
        p->u[i] = 0.0;
        p->v[i] = 0.0;
        p->ax[i] = 0.0;
        p->ay[i] = 0.0;
    }
    return p;
}


void exchange_momentum_f2p(fluid_grid* f, particle_system* p) {
    double dt = p->dt;
    p->avg_u = 0.0;
    p->avg_v = 0.0;
    p->avg_x = 0.0;
    p->avg_y = 0.0;
    p->avg_ax = 0.0;
    p->avg_ay = 0.0;
    p->avg_mass = 0.0;
    p->avg_drag = 0.0;

    for (int i = 0; i < p->num_particles; i++) {
        int cellpos_x = p->x[i] / f->_dx;
        int cellpos_y = p->y[i] / f->_dy;
        int I = cellpos_y * f->_nx + cellpos_x;

        float dxmul = 1.0, dymul = 1.0;

        if (cellpos_x <= 0 || cellpos_x >= f->_nx - 1)
            dxmul = 0.0;
        else
            dxmul = diffx(f->_p, I, f->_dx);
        
        if (cellpos_y <= 0 || cellpos_y >= f->_ny - 1)
            dymul = 0.0;
        else
            dymul = diffy(f->_p, I, f->_dy, f->_nx);
        

        
        
        p->ax[i] = p->drag[i] * (f->_u[I] - p->u[i]) - 0.5 * p->m_inv[i] * (dxmul + dymul);
        p->ay[i] = p->drag[i] * (f->_v[I] - p->v[i]) - 0.5 * p->m_inv[i] * (dxmul + dymul);
        p->u[i] += p->ax[i] * dt;
        p->v[i] += p->ay[i] * dt;
        p->x[i] += p->u[i] * dt;
        p->y[i] += p->v[i] * dt;
        p->avg_ax += p->ax[i];
        p->avg_ay += p->ay[i];
        p->avg_x += p->x[i];
        p->avg_y += p->y[i];
        p->avg_u += p->u[i];
        p->avg_v += p->v[i];
        p->avg_mass += 1.0 / p->m_inv[i];
        p->avg_drag += p->drag[i];
    }
    p->avg_ax /= (double)p->num_particles;
    p->avg_ay /= (double)p->num_particles;
    p->avg_x /= (double)p->num_particles;
    p->avg_y /= (double)p->num_particles;
    p->avg_u /= (double)p->num_particles;
    p->avg_v /= (double)p->num_particles;
    p->avg_mass /= (double)p->num_particles;
    p->avg_drag /= (double)p->num_particles;
}

double clamp(double x, double l, double u) {
    return (x > l)?((x < u)?x:u):l;
}

int iclamp(int x, int l, int u) {
    return (x > l)?((x < u)?x:u):l;
}

void contain_particles_single(particle_system* p) {
    for (int i = 0; i < p->num_particles; i++) {
        
        if (p->x[i] > p->nx || p->x[i] < 0.0) {
            p->u[i] *= -1.0;

        }
        if (p->y[i] > p->ny || p->y[i] < 0.0) {
            p->v[i] *= -1.0;
        }
        p->x[i] = clamp(p->x[i], 0.0, p->nx);
        p->y[i] = clamp(p->y[i], 0.0, p->ny);
    }

}

void exchange_inter_particles(particle_system** stacked_p_grid, int nx, int ny) {
    for (int Y = 1; Y < ny - 1; Y++) {
        for (int X = 1; X < nx - 1; X++) {
            particle_system* p = stacked_p_grid[Y * nx + X];
            for (int pn = 0; pn < p->num_particles; pn++) {
                int toleft = 0, toright = 0, totop = 0, tobottom = 0;

                if (p->x[pn] > (double)p->nx) {
                    toright = 1;
                }

                if (p->y[pn] > (double)p->ny) {
                    tobottom = 1;
                }

                if (p->x[pn] < 0.0) {
                    toleft = 1;
                }

                if (p->y[pn] < 0.0) {
                    totop = 1;
                }

                // swap top with current
                int id = p->ids[pn];
                double u = p->u[pn];
                double v = p->v[pn];
                double x = p->x[pn];
                double y = p->y[pn];
                double ax = p->ax[pn];
                double ay = p->ay[pn];
                double m_inv = p->m_inv[pn];
                double drag = p->drag[pn];
                
                p->ids[pn] = p->ids[p->num_particles - 1];
                p->u[pn] = p->u[p->num_particles - 1];
                p->v[pn] = p->v[p->num_particles - 1];
                p->x[pn] = p->x[p->num_particles - 1];
                p->y[pn] = p->y[p->num_particles - 1];
                p->ax[pn] = p->ax[p->num_particles - 1];
                p->ay[pn] = p->ay[p->num_particles - 1];
                p->m_inv[pn] = p->m_inv[p->num_particles - 1];
                p->drag[pn] = p->drag[p->num_particles - 1];

                p->ids[p->num_particles - 1] = id;
                p->u[p->num_particles - 1] = u;
                p->v[p->num_particles - 1] = v;
                p->x[p->num_particles - 1] = x;
                p->y[p->num_particles - 1] = y;
                p->ax[p->num_particles - 1] = ax;
                p->ay[p->num_particles - 1] = ay;
                p->m_inv[p->num_particles - 1] = m_inv;
                p->drag[p->num_particles - 1] = drag;

                // now curr particle is on top
                
                particle_system* pdest = NULL;
                // find correct destination system

                if (toleft) {
                    pdest = stacked_p_grid[Y * nx + X - 1];
                    if (totop)
                        pdest = stacked_p_grid[(Y - 1) * nx + X - 1];    
                    if (tobottom)
                        pdest = stacked_p_grid[(Y + 1) * nx + X - 1];    
                }
                if (toright) {
                    pdest = stacked_p_grid[Y * nx + X + 1];
                    if (totop)
                        pdest = stacked_p_grid[(Y - 1) * nx + X + 1];    
                    if (tobottom)
                        pdest = stacked_p_grid[(Y + 1) * nx + X + 1];
                }

                // if dest is not null

                if (pdest != NULL) {

                    // pop from current stack
                    p->num_particles--;
                    // push into dest stack
                    pdest->u[pdest->num_particles] = p->u[p->num_particles];
                    pdest->v[pdest->num_particles] = p->v[p->num_particles];
                    pdest->ids[pdest->num_particles] = p->ids[p->num_particles];
                    if (toleft) {
                        pdest->x[pdest->num_particles] = pdest->nx + p->x[p->num_particles];
                        if (totop)
                            pdest->y[pdest->num_particles] = pdest->ny + p->y[p->num_particles];
                        if (tobottom)
                            pdest->y[pdest->num_particles] = p->y[p->num_particles] - p->ny;
                        
                    }
                    if (toright) {
                        pdest->x[pdest->num_particles] = p->x[p->num_particles] - p->nx;
                        if (totop)
                            pdest->y[pdest->num_particles] = pdest->ny + p->y[p->num_particles];
                        if (tobottom)
                            pdest->y[pdest->num_particles] = p->y[p->num_particles] - p->ny;
                    }
                    pdest->ax[pdest->num_particles] = p->ax[p->num_particles];
                    pdest->ay[pdest->num_particles] = p->ay[p->num_particles];
                    pdest->m_inv[pdest->num_particles] = p->m_inv[p->num_particles];
                    pdest->drag[pdest->num_particles] = p->drag[p->num_particles];
                    pdest->num_particles++;
                }
                else {
                    // reflect velocity
                    if (toleft) {
                        p->x[p->num_particles - 1] = 0.001;
                    }
                    if (toright) {
                        p->x[p->num_particles - 1] = p->nx - 0.001;
                    }
                    if (totop) {
                        p->y[p->num_particles - 1] = 0.001;
                    }
                    if (tobottom) {
                        p->y[p->num_particles - 1] = p->ny - 0.001;
                    }
                    p->u[p->num_particles - 1] *= -(float)(toleft || toright);
                    p->v[p->num_particles - 1] *= -(float)(totop || tobottom);

                }
            }
        
        }
    }
}

int max(int a, int b) {
    return (a > b)?a:b;
}
int min(int a, int b) {
    return (a <= b)?a:b;
}
void exchange_force_p2f(fluid_grid* f, particle_system* p) {
    double dt = p->dt;
    for (int i = 0; i < p->num_particles; i++) {
        double cx = ((double)p->x[i] / (double)f->_dx);
        double cy = ((double)p->y[i] / (double)f->_dy);
        long long cellpos_x = (long long)cx;
        long long cellpos_y = (long long)cy;
        long long I = cellpos_y * f->_nx + cellpos_x;
        f->Fx[I] = p->m_inv[i] * p->ax[i];
        f->Fy[I] = p->m_inv[i] * p->ay[i];
        for (int y = - p->aoe / 2 ; y < p->aoe / 2; y++) {
            for (int x = - p->aoe / 2; x < p->aoe / 2; x++) {
                int YY = iclamp(cellpos_y + y, 0, f->_ny - 1);
                int XX = iclamp(cellpos_x + x, 0, f->_nx - 1);
                int N = YY * f->_nx + XX;
                f->_u[N] = 0.0;
                f->_v[N] = 0.0;
                
            }
        }
        
    }
}

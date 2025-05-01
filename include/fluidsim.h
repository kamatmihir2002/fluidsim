#ifndef FLUID_SIM_H
#define FLUID_SIM_H

double diffx(double* a, int i, double dx);

double diffy(double* a, int i, double dy, int nx);

double diff2x(double* a, int i, double dx2);

double diff2y(double* a, int i, double dy2, int nx);

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

typedef struct fluid_grid {
    double* _u;
    double* _v;
    double* _u_next;
    double* _v_next;
    double* accel_x;
    double* accel_y;
    double* Fx;
    double* Fy;
    double* bc;

    double* b;

    double u_avg;
    double v_avg;
    double p_avg;
    
    double* _p_diff;
    double* _u_diff;
    double* _v_diff;

    double* _p_prol;
    double* _u_prol;
    double* _v_prol;

    double* _p_next;
    double* _p;

    size_t _nx;
    size_t _ny;

    float _dx;
    float _dy;
    float _dt;

    float rho;
    float nu;

    float x_domain_extents[2];
    float y_domain_extents[2];

    double Fx_avg;
    double Fy_avg;

} fluid_grid;


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

// void particle_system_fluid_momentum_exchange(particle_system* p, fluid_grid* f) {
    
// }


fluid_grid* fluid_grid_new(float x_extents[2], float y_extents[2], size_t x_subdiv, size_t y_subdiv);

void enforce_bc(fluid_grid* f);

void enforce_neumann_bc(fluid_grid* f, float delx, float dely, int include_left, int include_right, int include_top, int include_bottom);

void ld_cavity(fluid_grid* f, float velx, float vely, int col);


void set_properties(fluid_grid* f, float rho, float nu, float dt, double* bc);


void poisson_rhs(fluid_grid* f);

double pressure_poisson_single(fluid_grid* f, int include_left, int include_right, int include_top, int include_bottom);

void pressure_poisson(fluid_grid* f, double* b, int include_top, int include_bottom);

void poisson_swap(fluid_grid* f);

void accelerate_fluid(fluid_grid* f);

void advance_fluid(fluid_grid* f);

void print_fluid_grid(FILE* f, fluid_grid* fl, size_t timestep);

double fl_restrict_p(double* p, fluid_grid** grids, int Nx, int Ny);

double fl_restrict_u(double* u, fluid_grid** grids, int Nx, int Ny);

double fl_restrict_v(double* v, fluid_grid** grids, int Nx, int Ny);

double fl_restrict_F(double* Fx, double* Fy, fluid_grid** grids, int Nx, int Ny);

double fl_prolongate_p(double p, fluid_grid* grid);

double fl_prolongate_u(double u, fluid_grid* grid);

double fl_prolongate_v(double v, fluid_grid* grid);


void diff(double* result, double* A, double* B, int nx, int ny);

void match_interface_p(fluid_grid** grid, int nx, int ny);

void match_interface_uv(fluid_grid** grid, int nx, int ny);

typedef struct particle_system {
    double* x;
    double* y;

    double* u;
    double* v;

    double* ax;
    double* ay;

    double* m_inv;
    double* drag;
    size_t* ids;
    int aoe;


    double dt;

    double nstart_x;
    double nstart_y;
    double nx;
    double ny;

    double avg_mass;
    double avg_drag;
    double avg_x;
    double avg_y;
    double avg_u;
    double avg_v;
    double avg_ax;
    double avg_ay;

    size_t num_particles;

    size_t particle_capacity;
} particle_system;

particle_system* particle_system_new(double x_extents, double y_extents, double init_x, double init_y, int num_particles, int capacity, double mass, double drag, int aoe, double dt);


void exchange_momentum_f2p(fluid_grid* f, particle_system* p);

double clamp(double x, double l, double u);

int iclamp(int x, int l, int u);

void contain_particles_single(particle_system* p);

void exchange_inter_particles(particle_system** stacked_p_grid, int nx, int ny);

int max(int a, int b);

int min(int a, int b);

void exchange_force_p2f(fluid_grid* f, particle_system* p);


#endif
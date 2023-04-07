#include <mdlib.h>
#include <omp.h>


/* first part: propagate velocities by half and positions by full step */
// NB: not static since we need to call it from the test
static void velverlet_first_half(mdsys_t *sys)
{
    int i;
    // execute the loop in parallel
    #pragma omp parallel for
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}


/* second part: propagate velocities by another half step */
static void velverlet_second_half(mdsys_t *sys)
{
    int i;
    // execute the loop in parallel
    #pragma omp parallel for
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}


/* velocity verlet */
void velverlet(mdsys_t *sys, cells_t *cells)
{
    /* first part: propagate velocities by half and positions by full step */
    if (sys->mpirank == 0) {
        velverlet_first_half(sys);
    }
    /* compute forces and potential energy */
    force(sys, cells);
    /* second part: propagate velocities by another half step */
    if (sys->mpirank == 0) {
        velverlet_second_half(sys);
    }
}

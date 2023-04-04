#include <mdlib.h>


/* first part: propagate velocities by half and positions by full step */
// NB: not static since we need to call it from the test
void velverlet_first_half(mdsys_t *sys)
{
    int i;
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
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}


/* velocity verlet */
void velverlet(mdsys_t *sys)
{
    /* first part: propagate velocities by half and positions by full step */
    velverlet_first_half(sys);
    /* compute forces and potential energy */
    force(sys);
    /* second part: propagate velocities by another half step */
    velverlet_second_half(sys);
}


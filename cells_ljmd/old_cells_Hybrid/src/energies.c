#include <mdlib.h>

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    int i;
    double ekin_tmp = 0.0;
    sys->ekin=0.0;

    // execute the loop in parallel
    #pragma omp parallel for reduction(+:ekin_tmp)
    for (i=0; i<sys->natoms; ++i) {
        ekin_tmp += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }

    // update the global variable
    sys->ekin = ekin_tmp;
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

#ifndef MDLIB_MPI_H
#define MDLIB_MPI_H

// for MPI
#include <mpi.h>

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    // for MPI parallelization we define extra sys variables
    MPI_Comm mpicomm;
    int nsize, mpirank;
    // as well as force vectors for parallelization
    double *fx_mpi, *fy_mpi, *fz_mpi;
};

typedef struct _mdsys mdsys_t;

// compute forces
void force(mdsys_t *sys);

#endif

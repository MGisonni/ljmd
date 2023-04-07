#ifndef MDLIB_SERIAL_H
#define MDLIB_SERIAL_H

/* structure to hold the complete information
 * about the MD system */
struct _mdsys{
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};

typedef struct _mdsys mdsys_t;

// compute forces
void force(mdsys_t *sys);

#endif

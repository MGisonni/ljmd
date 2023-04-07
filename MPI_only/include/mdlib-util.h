#ifndef MDLIB_UTIL_H
#define MDLIB_UTIL_H


// for input and output declarations
#include <stdio.h>
// for struct mdsys_t
#ifdef MPIYES
    #include <mdlib-mpi.h>
#else
    #include <mdlib-serial.h>
#endif


/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */


// compute kinetic energy
void ekin(mdsys_t *sys);
// velocity verlet
void velverlet(mdsys_t *sys);
void velverlet_first_half(mdsys_t *sys);
// set vector elements to zero
void azzero(double *d, const int n);
// walltime
double wallclock();
// input function
int get_a_line(FILE *fp, char *buf);
// output function
void output(mdsys_t *sys, FILE *erg, FILE *traj);
// reading files
int readinput (mdsys_t *sys, int * nprint, char restfile[BLEN], char trajfile[BLEN], char ergfile[BLEN]);
int readrest (mdsys_t *sys, char restfile[BLEN]);

#endif


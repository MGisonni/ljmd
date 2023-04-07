#ifndef MDLIB_H
#define MDLIB_H

// for input and output declarations
#include <stdio.h>
#include <mpi.h>



/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */



/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps,nthreads;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    // for MPI parallelization we define extra sys variables
    MPI_Comm mpicomm;
    int mpisize, mpirank;
    // as well as force vectors for parallelization
    double *fx_mpi, *fy_mpi, *fz_mpi;
    // number of cells in each direction and total number of cells
    int ncell_x, ncell; 
    double cell_size;
    // number of pairs of cells and list of pairs
    int npairs;
    int *pairs;
};
typedef struct _mdsys mdsys_t;



// cells_t structure
struct _cells {
    // identifying number of the cell
    int id;
    // boundary of x, y, z directions
    double x_min, x_max, y_min, y_max, z_min, z_max;
    // identifying numbers of 27 neighbor cells (including itself, with periodic boundary conditions)
    int *n_list;
    // list of atoms in the cell
    int catoms;
    int *a_list;
};
typedef struct _cells cells_t;



// cells functions
void initialize_cells(mdsys_t *sys, cells_t *cells);
void initialize_pairs(mdsys_t *sys, cells_t *cells);
void distribute_atoms(mdsys_t *sys, cells_t *cells);

// compute forces
void force(mdsys_t *sys, cells_t *cells);

// compute kinetic energy
void ekin(mdsys_t *sys);

// velocity verlet
void velverlet(mdsys_t *sys, cells_t *cells);

// set vector elements to zero
void azzero(double *d, const int n);
void iazzero(int *d, const int n);

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

/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <ljmd.h>
#include <mdlib-util.h>
#include <mdlib.h>


/* main */
int main(int argc, char **argv)
{
    int nprint, return_value;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *traj=NULL,*erg=NULL;
    mdsys_t sys;
    double t_start=0.0;

    
    /* initialize MPI communicator */
    MPI_Init(&argc, &argv);
    sys.mpicomm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &sys.nsize);
    // DEBUG: let each process print its rank
    printf("\n");
    printf("Rank %d of %d\n", sys.mpirank, sys.nsize);

    // NB: kinda useless since this code runs only if MPICH is defined
    #ifdef MPIYES
        printf("MPICH is defined\n");
    #else
        printf("MPICH is not defined\n");
    #endif

    // Initialization operations for the master (input reading, print info, etc.)
    if (sys.mpirank == 0) {

        // print version
        printf("LJMD version %3.1f\n", LJMD_VERSION);

        // only master process keeps track of time
        t_start = wallclock();

        /* read input file */
        return_value = readinput(&sys, &nprint, restfile, trajfile, ergfile);
        if (return_value != 0) {
            printf("Error reading input file\n");
            return return_value;
        }

        // print info
        printf("Number of particles: %d   Number of steps: %d  Time step: %f  Print frequency: %d Rest file: %s  Trajectory file: %s  Energy file: %s Box size: %f    Temperature: %f   Cutoff: %f    Mass: %f    Epsilon: %f    Sigma: %f    Number of processes: %d  \n", sys.natoms, sys.nsteps, sys.dt, nprint, restfile, trajfile, ergfile, sys.box, sys.temp, sys.rcut, sys.mass, sys.epsilon, sys.sigma, sys.nsize);

        /* allocate memory */
        // NB: other processes are used for forces computations only, 
        // they don't need velocity vectors and have f_mpi vectors for forces
        sys.vx=(double *)malloc(sys.natoms*sizeof(double));
        sys.vy=(double *)malloc(sys.natoms*sizeof(double));
        sys.vz=(double *)malloc(sys.natoms*sizeof(double));
        sys.fx=(double *)malloc(sys.natoms*sizeof(double));
        sys.fy=(double *)malloc(sys.natoms*sizeof(double));
        sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    }

    // Broadcast necessary system parameters to all processes (to be reviewed)
    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, sys.mpicomm);
    MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, sys.mpicomm);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.nfi, 1, MPI_INT, 0, sys.mpicomm);


    // Allocate memory for positions
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    // Allocate memory for forces
    sys.fx_mpi=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy_mpi=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz_mpi=(double *)malloc(sys.natoms*sizeof(double));


    /* read restfile */
    if (sys.mpirank == 0) {
        return_value = readrest(&sys, restfile);
        if (return_value != 0) {
            printf("Error reading restart file\n");
            return return_value;
        }
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);


    if (sys.mpirank == 0) {
        ekin(&sys);
    }

    
    /* master manages output files */
    if (sys.mpirank == 0) {
    
        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        printf("Startup time: %10.3fs\n", wallclock()-t_start);
        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);

        // print informations about the system
        printf("Number of atoms: %d \t Number of steps: %d \t Number of processes: %d \t Number of atoms per    process: %d \t Number of steps per process: %d\n", sys.natoms, sys.nsteps, sys.nsize, sys.natoms/sys.nsize, sys.nsteps/sys.nsize);   

        /* reset timer */
        t_start = wallclock();
    }



    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* master writes output, if requested */
        if (sys.mpirank == 0) {
            if ((sys.nfi % nprint) == 0)
                output(&sys, erg, traj);
        }

        /* propagate system and recompute energies */
        velverlet(&sys);

        /* compute kinectic energy */
        if (sys.mpirank == 0) {
            ekin(&sys);
        }
    }
    /**************************************************/



    /* clean up: close files, free memory */
    if (sys.mpirank == 0) {
        
        // print final results
        printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
        
        // close files
        fclose(erg);
        fclose(traj);

        // free velocity and force vectors
        free(sys.vx);
        free(sys.vy);
        free(sys.vz);
        free(sys.fx);
        free(sys.fy);
        free(sys.fz);
    }

    /* clean up positions and force-mpi vectors as well */
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.fx_mpi);
    free(sys.fy_mpi);
    free(sys.fz_mpi);

    /* finalize MPI */
    MPI_Finalize();

    return 0;
}

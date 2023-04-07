/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <ljmd.h>
#include <mdlib.h>
#include <omp.h>


/* main */
int main(int argc, char **argv)
{
    int i, nprint, return_value=0;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *traj=NULL,*erg=NULL;
    mdsys_t sys;
    double t_start=0.0;

    /* initialize MPI communicator */
    MPI_Init(&argc, &argv);
    sys.mpicomm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &sys.mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &sys.mpisize);

    // set number of threads
    if (argc>1) {
        omp_set_num_threads(atoi(argv[1]));
        sys.nthreads = atoi(argv[1]);
    } else {
        sys.nthreads = omp_get_max_threads();
    }

    // Initialization operations for the master (input reading, print info, etc.)
    if (sys.mpirank == 0) {

        // print version
        printf("\nLJMD version %3.1f\n", LJMD_VERSION);

        // print number of processes and threads
        printf("Number of processes: %d, number of threads: %d\n\n", sys.mpisize, sys.nthreads);

        // only master process keeps track of time
        t_start = wallclock();

        /* read input file */
        return_value = readinput(&sys, &nprint, restfile, trajfile, ergfile);
        if (return_value != 0) {
            printf("Error reading input file\n");
            return return_value;
        }

        /* allocate memory */
        // NB: other processes are used for forces computations only, 
        // they don't need velocity vectors and have f_mpi vectors for forces
        sys.vx=(double *)malloc(sys.natoms*sizeof(double));
        sys.vy=(double *)malloc(sys.natoms*sizeof(double));
        sys.vz=(double *)malloc(sys.natoms*sizeof(double));
        sys.fx=(double *)malloc(sys.natoms*sizeof(double));
        sys.fy=(double *)malloc(sys.natoms*sizeof(double));
        sys.fz=(double *)malloc(sys.natoms*sizeof(double));


        /* divide the 3D box in cells of size rcut*rcut*rcut */
        sys.cell_size = sys.box/((int) (sys.box/sys.rcut));  // surely greater or equal than cutoff
        sys.ncell_x   = (int) (sys.box/sys.rcut);            // number of cells in each direction
        sys.ncell     = sys.ncell_x*sys.ncell_x*sys.ncell_x;

        // exit if number of cells is too small (less than 27)
        if (sys.ncell < 27) {
            printf("Error: number of cells is too small for cell implementation (less than 27)\n");
            return 1;
        }
    }

    
    // Broadcast system parameters to all processes
    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, sys.mpicomm);
    MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, sys.mpicomm);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.nfi, 1, MPI_INT, 0, sys.mpicomm);
    MPI_Bcast(&sys.cell_size, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.ncell_x, 1, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(&sys.ncell, 1, MPI_DOUBLE, 0, sys.mpicomm);
  

    // Allocate memory for positions
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    // Allocate memory for forces
    sys.fx_mpi=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy_mpi=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz_mpi=(double *)malloc(sys.natoms*sizeof(double));




    /* read restfile: get initial positions */
    if (sys.mpirank == 0) {
        return_value = readrest(&sys, restfile);
        if (return_value != 0) {
            printf("Error reading restart file\n");
            return return_value;
        }
    }

    // Broadcast positions to all processes
    MPI_Bcast(sys.rx, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(sys.ry, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);
    MPI_Bcast(sys.rz, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);

    
    /* Cells initialization: all processes have all informations on cells */
    /* (could be improved: tracking carefully the cells, using ghosts cells communication btw processes...) */

    /* allocate ncell cells_t structures */
    cells_t *cells = (cells_t *) malloc(sys.ncell*sizeof(cells_t));
    /* initialize cells */
    initialize_cells(&sys, cells);
    /* compute cells pairs of neighbours */
    initialize_pairs(&sys, cells);
    /* distribute atoms amongst cells */
    distribute_atoms(&sys, cells);


    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys, cells);

    // Let master process compute kinetic energy
    if (sys.mpirank == 0) {
        ekin(&sys);
    }



    /* master manages output files */
    if (sys.mpirank == 0) {
    
        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");

        printf("Startup time: %10.3fs\n", wallclock()-t_start);
        printf("Starting simulation with %d atoms for %d steps.\n\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);

        /* reset timer */
        t_start = wallclock();
    }



    // /**************************************************/
    // /* main MD loop */
    // for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

    //     // DEBUG print
    //     printf("About to start step %d\n", sys.nfi);

    //     /* write output, if requested */
    //     if ((sys.nfi % nprint) == 0)
    //         output(&sys, erg, traj);

    //     /* propagate system and recompute energies */
    //     velverlet(&sys, cells);
    //     ekin(&sys);

    //     // DEBUG print
    //     printf("Evolved step %d\n", sys.nfi);

    //     /* redistribute atoms amongst cells */
    //     for(int i=0; i<ncell; i++) {
    //         cells[i].catoms = 0;
    //         iazzero(cells[i].a_list,sys.natoms);
    //     }

    //     // DEBUG print
    //     printf("Azzeroed step %d\n", sys.nfi);

    //     for (int i=0; i<sys.natoms; i++) {
    //         // compute cell index, with periodic boundary conditions
    //         int a_x = ((int) (sys.rx[i]/cell_size))%ncell_x;
    //         int a_y = ((int) (sys.ry[i]/cell_size))%ncell_x;
    //         int a_z = ((int) (sys.rz[i]/cell_size))%ncell_x;
    //         int cell_index = a_x + a_y*ncell_x + a_z*ncell_x*ncell_x;
    //         // add atom to the list of atoms in the cell
    //         cells[cell_index].a_list[cells[cell_index].catoms] = i;
    //         cells[cell_index].catoms++;
    //     }

    //     // DEBUG print
    //     printf("Ridistributed step %d\n\n", sys.nfi);

    // }
    // /**************************************************/

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* master writes output */
        if (sys.mpirank == 0) {
            if ((sys.nfi % nprint) == 0)
                output(&sys, erg, traj);
        }

        /* propagate system and recompute energies */
        velverlet(&sys, cells);
        
        /* compute kinectic energy */
        if (sys.mpirank == 0) {
            ekin(&sys);
        }

        /* broadcast new position of the atoms */
        MPI_Bcast(sys.rx, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);
        MPI_Bcast(sys.ry, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);
        MPI_Bcast(sys.rz, sys.natoms, MPI_DOUBLE, 0, sys.mpicomm);

        /* redistribute atoms amongst cells */
        distribute_atoms(&sys, cells);

    }


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

    // free cells allocations
    for(i=0; i<sys.ncell; i++) {
        free(cells[i].a_list);
    }
    for(i=0; i<sys.ncell; i++) {
        free(cells[i].n_list);
    }
    free(cells);
    free(sys.pairs);

    /* finalize MPI */
    MPI_Finalize();

    return 0;
}

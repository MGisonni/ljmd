/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <ljmd.h>
#include <mdlib.h>


/* main */
int main(int argc, char **argv)
{
    int i, nprint, return_value=0;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *traj,*erg;
    mdsys_t sys;
    double t_start;
     
    printf("LJMD version %3.1f\n", LJMD_VERSION);

    // global timer
    t_start = wallclock();

    /* read input file */
    return_value = readinput(&sys, &nprint, restfile, trajfile, ergfile);
    if (return_value != 0) {
        printf("Error reading input file\n");
        return return_value;
    }

    // print parameters
    printf("nsteps = %d, dt = %f, mass = %f, epsilon = %f, sigma = %f, box = %f, rcut = %f, nprint = %d\n", sys.nsteps, sys.dt, sys.mass, sys.epsilon, sys.sigma, sys.box, sys.rcut, nprint);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* read restfile */
    return_value = readrest(&sys, restfile);
    if (return_value != 0) {
        printf("Error reading restart file\n");
        return return_value;
    }


    /* divide the 3D box in cells of size rcut*rcut*rcut */
    sys.cell_size = sys.box/((int) (sys.box/sys.rcut));  // surely greater or equal than cutoff
    sys.ncell_x   = (int) (sys.box/sys.rcut);            // number of cells in each direction
    sys.ncell     = sys.ncell_x*sys.ncell_x*sys.ncell_x;

    /* print parameters vs box and rcut */ 
    printf("\nbox = %f, rcut = %f\n", sys.box, sys.rcut);
    printf("cell_size = %f, ncell_x = %d, ncell = %d\n\n", sys.cell_size, sys.ncell_x, sys.ncell);

    // exit if number of cells is too small (less than 27)
    if (sys.ncell < 27) {
        printf("Error: number of cells is too small for cell implementation (less than 27)\n");
        return 1;
    }
    

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
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /* reset timer */
    t_start = wallclock();


    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {


        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys, cells);
        ekin(&sys);

        /* redistribute atoms amongst cells */
        distribute_atoms(&sys, cells);


    }
    /**************************************************/

    // print times
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);

    /* clean up: close files, free memory */
    fclose(erg);
    fclose(traj);
    /* clean up: close files, free memory */
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    // cells allocations
    for(i=0; i<sys.ncell; i++) {
        free(cells[i].a_list);
    }
    for(i=0; i<sys.ncell; i++) {
        free(cells[i].n_list);
    }
    free(cells);
    free(sys.pairs);

    return 0;
}

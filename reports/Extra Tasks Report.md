# Cells implementation

As extra individual task I implemented the cells subdivision for the system. I adopted the following strategy:

## Vanilla implementation (`cells_ljmd/cells_Vanilla/` folder)

Most of the following can be found in the main `ljmd.c` code as well as in the `cells_utilities.c` one.

+ Create a structure `cells_t` holding informations for each cell, such as `id`, its list of neighbours, the number and the list of atoms belonging to her

+ Subdivide the `box` in the maximum number of (cubic) cells of radius greater or equal than `rcut` as

	    sys.cell_size = sys.box/((int) (sys.box/sys.rcut));  // surely greater or equal than cutoff
	    sys.ncell_x   = (int) (sys.box/sys.rcut);            // number of cells in each direction
	    sys.ncell     = sys.ncell_x*sys.ncell_x*sys.ncell_x; // total number of cells

+ Allocate an array of structures `cells_t` of size the total number of cells in the box. Hence initialize every cell in a 3D indices fashion as

	    int cell_index = i + j*ncell_x + k*ncell_x*ncell_x;

+ Distribute atoms in each cell according to their position in the box, accounting for boundary conditions (atoms might get out of the box during time evoluton) and for negative coordinates (which I rescale to positive in order to apply `%` remainders operations)

+ Compute (once) the list of neighbours of each cell and store them into the relative structure variable

+ Using their unique ids and the just created neighbouring list, store the (unique) pairs of neighbouring cells in a struct variable of `sys_t`. A cell is not considered a neighbour of itself.

+ At every time step, we first compute the forces between different cells through three nested loops as

	    /* loop over all pairs of cells */
	    for (i=0; i<sys->npairs; i++) {
		/* loop over all particles in first cell */
		for (j=0; j<cells[sys->pairs[2*i]].catoms; j++) {
		    /* loop over all particles in second cell */
		    for (k=0; k<cells[sys->pairs[2*i+1]].catoms; k++) {

    notice that in this way we can directly apply Newton 3rd law.

+ Hence we compute forces inside every cell as in the basic `ljmd` implementation

+ Once the computation of forces is done, we evolve positions of atoms via `second_half_velverlet` and hence redistribute atoms amongst cells based on the new configuration



### Side note on allocation

There are two obvious ways to allocate space for each cell to gather positions of atoms under their domain: either allocate for every cell a big `sys->natoms` size array at the beginning of the simulation, or allocate and deallocate at every time step an array with the exact number of atoms falling in the cells.

The former procedure favours computation time above memory, the latter the other way around. I tried to implement both and got very similar results on time comparison, with a sligh advantage for the first one. Nonetheless, I preferred the first one, as for the `argon_*` codes proposed I computed to have enough space to run the simulations.

Precisely, for the `argon_78732` case, `1728` cells would be needed, each allocating a `78732`array of integers. This should amount to roughly `544MB` which I found still reasonable.

The second procedure can still be found at the `cells_ljmd/cells_vanilla_v1/` folder.



## MPI implementation (`cells_ljmd/cells_MPI/` folder)

Introducing MPI I adopted the following:

+ *Let every MPI process know the positions of every atom at every step.* In this way the redistribution of atoms at each time step becomes much easier. This is kinda of a simplification and not very efficient, it would be nice to implement a reasonable geometric subdivision between cells as well as *ghost cells* communications between different processes.

+ The above implies an `MPI_Bcast` to be performed from the master at each time step, commmunicating the updated atoms configuration, as the master is the delegate process to execute the velocity Verlet calls.

+ The split of forces computation amongst different processes is performed assigning to each process the (roughly) same amount of cells pairs to compute forces on, as in

	    /* loop over all pairs of cells, distribute among processes */
	    for (i_mpi=0; i_mpi<sys->npairs; i_mpi=i_mpi+sys->mpisize) {
		i = i_mpi + sys->mpirank;
		if (i >= sys->npairs) break;
		int cell1 = sys->pairs[2*i];
		int cell2 = sys->pairs[2*i+1];
		/* loop over all particles in first cell */
		for (j=0; j<cells[sys->pairs[2*i]].catoms; j++) {
		    /* loop over all particles in second cell */
		    for (k=0; k<cells[sys->pairs[2*i+1]].catoms; k++) {
                

+ For the *inter* cell forces, the redistribution is performed on the total number of cells as in

	    // loop over all cells, distribute among processes */ 
	    for (i_mpi=0; i_mpi<sys->ncell; i_mpi=i_mpi+sys->mpisize) {
		i = i_mpi + sys->mpirank;
		if (i >= sys->ncell) break;

+ Once all processes are done with forces computation, we perform a reduction of forces and potential energy to the master.


## OpenMP implementation (`cells_ljmd/cells_Hybrid/` folder)

+ `#pragma omp` clauses were implemented in the master process for the computation of kinectic energy as well as in the velocity Verlet functions.

+ In the forces computation, I applied a `#pragma omp parallel` call to the second `for` loop, precisely the `j-loop`. In order to avoid cache conflicts, as well as allocation/deallocation routines amongst threads, I decided to drop the Newton 3rd law implementation and just perform the inner `j` and `k` loops twice, switching the roles of the first and second cell in the pair. I actually tested the dropped implementation as well and it was performing slightly worse.

+ At the end of the parallel region a `reduction` is performed on the potential energy

+ At the end of the main loop, a `MPI_reduce` is performed for forces and potential energy






## Benchmarking


### Specifics

The simulation were run on 1,2,4,6,8 nodes of the `regular1` partition in `Ulysses` cluster. The following modules were loaded

    	module load cmake/3.15.4
	module load intel/2021.1
	module load openmpi3

+ The Benchmark is executed on the provided `argon_2916.*` and `argon_78732.*` files with, respectively, 1000 and 20 timesteps
+ All simulations are executed on nodes of `regular1` partition in Ulysses ***ADD SPECIFICS***
+ compiler is `GNU 8.3.0` 
+ flags are `-O3 -ffast-math -fomit-frame-pointer`



### Vanilla code

|         | argon_2916 | argon_78732 |
|---------|------------|-------------|
| serial  | 86.358s    | 46.594s   |


### MPI implementation

| MPI     | argon_2916 | argon_78732 |
|---------|------------|-------------|
| 2 nodes | 8.360s     | 4.937s      |
| 4 nodes | 4.680s     | 2.902s      |
| 6 nodes | 2.824s     | 1.812s      |



### Hybrid OpenMP+MPI implementation

| Hybrid   | argon_2916 | argon_78732 |
|----------|------------|-------------|
| 3 nodes  | 4.895s     | 3.107s      |
| 6 nodes  | 3.036s     | 1.854s      |
| 9 nodes  | 2.368s     | 1.441s      |
| 12 nodes | 2.085s     | 1.297s      |

### Report for the *Best Practices in Scientific Software Development*
The following report contains my contribution to the group project https://github.com/Project-MD-GPS/ljmd by GPS organization.
A report about my cell implementations is found in this same folder. 

# MPI implementation

For the overall MPI implementation I adopted the following:
+ parallelize forces computation
+ delegate to master the rest (reading and output, computation of kinectic energy, updating of positions and velocities)


### Version1: *Alternating split*

In the very first implementation, the main `i-loop` was splitted amongst threads in an alternating way, namely

    // main loop
    for(i=0; i < (sys->natoms); i += sys->nsize) {
        i_mpi = i + sys->mpirank;
        if (i_mpi >= sys->natoms) break;
            for(j=i_mpi+1; j < (sys->natoms); ++j) { 

Inside the loop, `forces` and `potential_energy` are computed; to this end the positions of the involved atoms must be known. Since the `j-loop` can run over all atoms, all processes need to be aware of the positions of all atoms. This requires a broadcast to be performed before entering the main loop, as

        MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

At the end of the loop, we need all MPI processes to communicate back to master the forces and potential energy they just computed, as

    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);


### Version2: *Balanced chunks split*

As an alternative to the above, and in hindsight with my group to use a `#pragma omp for` clause to split the `i-loop` on each process, I implemented the second basic version as well, which is: divide equally the `i-loop` in consecutive chunks amongst all processes, letting each process computing its `istart` and `iend`. 

    //main loop
    for(i_mpi = istart; i_mpi < iend; ++i_mpi) {
        for(j=i_mpi+1; j < (sys->natoms); ++j) {

For example:

	natoms=80, mpisize=4, then processes 0,1,2,3 are respectively assigned the:
		 0-19, 20-39, 40-59, 60-79 chunks of the i-loop

However, I realized this would actually cause an unbalanced workload amongst processes, as the `j-loop`s for the first processes would be much longer than those of the last ones (in the above example, process 0 will be loaded with more or less 3 times the work of process 3). Since we are calling a broadcast at the end of the `forces` code, a barrier will be in place and all the processes would have to wait for the slowest one.

To put a remedy on this, one can split the `i-loop` in half and then, on each half, assign the chunks of atoms to processes in reverse order. In the above example:

	natoms=80, mpisize=4, then processes 0,1,2,3 are respectively assigned the:
		 0-9 & 70-79, 10-19 & 60-69, 20-29 & 50-59, 30-39 & 40-49 chunks of the i-loop

This is almost evenly balanced and still retains the asset of each process working on contiguous memory addresses. Start and end for each process are stored in arrays `istart[2]` and `iend[2]` and computed as in the `mpi-bounds.c` code.

**NB** *This is not necessarily better than the "Version1: Alternating split", arguably **it is slightly worse** (see timings below). It was implemented mostly with the idea of better adapting to the OpenMP implementation, offering contiguous memory in each MPI process, and thus being in principle more cache friendly*





## MPI benchmarking


The code already comes equipped **with basic optimization** (removal of `pow` and `sqrt` functions, avoid recomputing constant values in loops) and implementation of **Newton's 3rd law**. I report the timings for both the implementations described above.


### Specifics

The simulation were run on 1,2,4,6,8 nodes of the `regular1` partition in `Ulysses` cluster. The following modules were loaded

    module load cmake/3.15.4
	module load intel/2021.1
	module load openmpi3

+ The Benchmark is executed on the provided `argon_108.*`, `argon_2916.*` and `argon_78732.*` files
+ All simulations are executed on nodes of `regular1` partition in Ulysses ***ADD SPECIFICS***
+ compiler is `GNU 8.3.0` 
+ flags are `-O3 -ffast-math -fomit-frame-pointer`




### *Alternating* MPI implementation

In the following we can see an almost perfect scaling from 2 to 8 nodes. Of course, we do pay a toll passing from the serial to the parallel version, due to the introduced communications of positions and forces at each time-step.

|         | argon_108 | argon_2916 | argon_78732 |
|---------|-----------|------------|-------------|
| serial  | 1.386s    | 47.269s    | 365.354s    |
| 2 nodes | 0.790s    | 26.885s    | 224.017s    |
| 3 nodes | 0.537s    | 18.430s    | 148.097s    |
| 4 nodes | 0.470s    | 14.236s    | 112.337s    |
| 6 nodes | 0.408s    | 9.663s     | 75.033s     |
| 8 nodes | 0.381s    | 7.746s     | 56.442s     |



### *Balanced Chunks* MPI implementation

The same considerations apply. As anticipated, this is slightly worse than the above.

|         | argon_108 | argon_2916 | argon_78732 |
|---------|-----------|------------|-------------|
| serial  | 1.386s    | 47.269s    | 365.354s    |
| 2 nodes | 0.762s    | 28.595s    | 235.094s    |
| 3 nodes | 0.527s    | 20.067s    | 153.769s    |
| 4 nodes | 0.428s    | 15.430s    | 118.349s    |
| 6 nodes | 0.372s    | 11.034s    | 79.521s     |
| 8 nodes | 0.378s    | 8.743s     | 59.928s     |

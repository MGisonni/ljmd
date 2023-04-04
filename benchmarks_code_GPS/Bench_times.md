# Benchmarking times for ljmd code

+ All simulations are executed on (full) nodes of `regular1` partition in Ulysses ***ADD SPECIFICS***

+ The Benchmark is executed on the provided `argon_108.*` and `argon_2916.*` files

+ compiler is `gnu 6.3.1` *(the reason being: is the one that comes with `intel/18.0.3.222`, the necessary mpodule in Ulysses to run mpi)*

+ **?** should we also change `set(CMAKE_BUILD_TYPE RelWithDebInfo Release)` **?** 

+ flags are `-O3 -ffast-math -fomit-frame-pointer`





## Vanilla code

+ Time for `argon_108`: 11.932s (with `-Wall -g` flags)
+ Time for `argon_108`: 3.251s

+ Time for `argon_2916`: 182.250s (with `-Wall -g` flags)
+ Time for `argon_2916`: 103.459s





## Code with forces optimization (with Newton 3rd law implemented)

+ Time for `argon_108`: 1.284s
+ Time for `argon_2916`: 40.677s




## MPI only, with optimized forces (no Newton 3rd law implementation)

We ran the simulation on 4 nodes of the `regular1` partition in `Ulysses` cluster. The following modules were loaded

    module load cmake/3.15.4
	module load intel/18.0.3.222
	module load openmpi/2.1.3

+ Time for `argon_108`: 1.235s
+ Time for `argon_2916`: 25.976s

In the last and heavier job, one can appreciate a perfect scaling with respect to the vanilla code.


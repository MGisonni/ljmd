rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cd examples
mpirun -np 4 ../ljmd.x < argon_108.inp 

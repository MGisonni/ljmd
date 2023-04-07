rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cd examples
# ./../ljmd-serial.x < argon_108.inp 
mpirun -np 4 ../ljmd.x 3 < argon_2916.inp 
mpirun -np 4 ../ljmd.x 3 < argon_78732.inp 
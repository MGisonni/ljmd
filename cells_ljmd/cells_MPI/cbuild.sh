rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cd examples
# ./../ljmd-serial.x < argon_108.inp 
mpirun -np 4 ../ljmd.x < argon_2916.inp 
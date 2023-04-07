rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cd examples
# ./../ljmd-serial.x < argon_108.inp 
./../ljmd-serial.x < argon_2916.inp 
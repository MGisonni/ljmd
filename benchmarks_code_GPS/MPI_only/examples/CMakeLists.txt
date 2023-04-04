# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)

# Copy files from the current directory to the build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/argon_108.inp
    ${CMAKE_CURRENT_BINARY_DIR}/argon_108.inp COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/argon_2916.inp
    ${CMAKE_CURRENT_BINARY_DIR}/argon_2916.inp COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/argon_108.rest
    ${CMAKE_CURRENT_BINARY_DIR}/argon_108.rest COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/argon_2916.rest
    ${CMAKE_CURRENT_BINARY_DIR}/argon_2916.rest COPYONLY)

# Also copy the reference files
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../reference/argon_108.dat
    ${CMAKE_CURRENT_BINARY_DIR}/../reference/argon_108.dat COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../reference/argon_2916.dat
    ${CMAKE_CURRENT_BINARY_DIR}/../reference/argon_2916.dat COPYONLY)

# add check target
add_custom_target(check_example
COMMAND ../ljmd.x < argon_108.inp
COMMAND head -10 argon_108.dat > a.dat
COMMAND head -10 ../reference/argon_108.dat > b.dat
COMMAND cmp a.dat b.dat || exit 1
COMMAND rm -f a.dat b.dat
)

# add clean target
add_custom_target(clean_example COMMAND rm -f *.xyz *.dat)
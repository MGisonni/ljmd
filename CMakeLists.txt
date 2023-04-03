# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C CXX) 
set(CMAKE_BUILD_TYPE RelWithDebInfo Release) 
# set(CMAKE_CXX_STANDARD 14)


# include directory for header files
include_directories(include)

find_package(MPI REQUIRED)
# if found we built the parallel version
if (MPI_FOUND)
    message(MPI_YES " This application does compile with MPI")
    # set(MPI_CXX_SKIP_MPICXX TRUE)
    # Add executable 
    add_executable(ljmd-mpi.x src/ljmd_mpi.c)
    # compile options
    target_compile_options(ljmd-mpi.x PRIVATE -Wall -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
    # add library for force computation
    add_library(mdlib src/utilities.c src/verlet_ti.c src/forces_mpi.c src/energies.c src/input.c src/output.c src/reading.c)
    target_link_libraries(mdlib m MPI::MPI_C)
    # target_link_libraries(mdlib m mpi)
    # link libraries to executable
    target_link_libraries(ljmd-mpi.x mdlib)
# else we build the serial version
else (MPI_FOUND)
    message(MPI_NO " This application cannot compile without MPI")
    # Add executable 
    add_executable(ljmd-serial.x src/ljmd.c)
    # compile options
    target_compile_options(ljmd-serial.x PRIVATE -Wall -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
    # add library for force computation
    add_library(mdlib src/utilities.c src/verlet_ti.c src/forces.c src/energies.c src/input.c src/output.c src/reading.c)
    target_link_libraries(mdlib m)
    # link libraries
    target_link_libraries(ljmd-serial.x m)
    target_link_libraries(ljmd-serial.x mdlib)
endif (MPI_FOUND)


# add subdirectories
add_subdirectory(examples)


# set target for serial version
add_custom_target(serial
    COMMAND make
)

# set target for cleaning
add_custom_target(clean_data_files
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples clean_example
)

# set target for checking
add_custom_target(check
    DEPENDS serial
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples check_example
)

enable_testing()

add_subdirectory(test)


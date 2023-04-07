# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C) 
set(CMAKE_BUILD_TYPE RelWithDebInfo Release) 


# include directory for header files
include_directories(include)
# add subdirectories
add_subdirectory(examples)


# MP is mandatory
find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
	message(OpenMP_FOUND " OpenMP found")
	set(BUILD_OMP_DEFAULT ON)
else (OpenMP_FOUND)
    message(OpenMP_NOT_FOUND " ERROR: This application needs OpenMP")
endif (OpenMP_FOUND)



# MPI is mandatory as well
find_package(MPI REQUIRED)

# if found we built the parallel version
if (MPI_FOUND)

    message(MPI_FOUND " This application will compile with MPI")
    # Add executable 

    add_executable(ljmd.x src/ljmd.c)
    # compile options
    target_compile_options(ljmd.x PRIVATE -Wall -g -O3 -fopenmp  -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
    
    # add library for auxiliary functions
    add_library(mdlib src/forces.c src/utilities.c src/verlet_ti.c src/energies.c src/input.c src/output.c src/reading.c src/cells_utilities.c)
    target_link_libraries(mdlib m MPI::MPI_C)
    # compile options for library
    target_compile_options(mdlib PRIVATE -Wall -g -O3 -fopenmp -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
    
    # link libraries to executable
    target_link_options(ljmd.x PRIVATE -fopenmp)
    target_link_libraries(ljmd.x mdlib)

# else we build the serial version
else (MPI_FOUND)
    message(MPI_NOT_FOUND " ERROR: This application needs MPI")
endif (MPI_FOUND)

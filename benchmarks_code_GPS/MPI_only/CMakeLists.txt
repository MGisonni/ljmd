# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C) 
set(CMAKE_BUILD_TYPE RelWithDebInfo Release) 
# set(CMAKE_CXX_STANDARD 14)


# include directory for header files
include_directories(include)



find_package(MPI REQUIRED)
# if found we built the parallel version
if (MPI_FOUND)
    message(MPI_FOUND " This application will compile with MPI")
    # Add executable 
    add_executable(ljmd.x src/ljmd-mpi.c)
    # compile options
    target_compile_options(ljmd.x PRIVATE -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0 -DMPIYES)
    # add library for force computation
    add_library(mdlib-mpi src/forces_mpi.c)
    target_link_libraries(mdlib-mpi m MPI::MPI_C)
    # compile options for library
    target_compile_options(mdlib-mpi PRIVATE -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0 -DMPIYES)
# add common library to serial and parallel version
add_library(mdlib-util src/utilities.c src/verlet_ti.c src/energies.c src/input.c src/output.c src/reading.c)
# compile options for library
target_link_libraries(mdlib-util m MPI::MPI_C)
target_compile_options(mdlib-util PRIVATE -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0 -DMPIYES)
    # link libraries to executable
    target_link_libraries(ljmd.x mdlib-util mdlib-mpi)
# else we build the serial version
else (MPI_FOUND)
    message(MPI_NOT_FOUND " This application cannot find MPI, will compile serial version")
    # Add executable 
    add_executable(ljmd.x src/ljmd-serial.c)
    # compile options
    target_compile_options(ljmd.x PRIVATE -Wall -pedantic -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
    # add library for force computation
    add_library(mdlib-serial src/forces-serial.c)
    target_link_libraries(mdlib-serial m)
# add common library to serial and parallel version
add_library(mdlib-util src/utilities.c src/verlet_ti.c src/energies.c src/input.c src/output.c src/reading.c)
# compile options for library
target_compile_options(mdlib-util PRIVATE -Wall -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)    # link libraries
    target_link_libraries(ljmd.x m)
    target_link_libraries(ljmd.x mdlib-util mdlib-serial)
endif (MPI_FOUND)


# add subdirectories
add_subdirectory(examples)


# # set target for serial version
# add_custom_target(serial
#     COMMAND make
# )

# # set target for cleaning
# add_custom_target(clean_data_files
#     COMMAND ${CMAKE_MAKE_PROGRAM} -C examples clean_example
# )

# # set target for checking
# add_custom_target(check
#     DEPENDS serial
#     COMMAND ${CMAKE_MAKE_PROGRAM} -C examples check_example
# )

# enable_testing()

# add_subdirectory(test)


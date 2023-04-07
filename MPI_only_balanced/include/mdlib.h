#ifndef MDLIB_H
#define MDLIB_H

#ifdef MPIYES
    #include <mdlib-mpi.h>
#else
    #include <mdlib-serial.h>
#endif

#endif


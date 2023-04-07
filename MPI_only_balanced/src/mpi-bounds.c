#include <mdlib.h>

void mpi_bounds(mdsys_t *sys) {

    /* compute sys->istart and sys->iend for each mpi process, dividing main i-loop */
    int isize[2], chunk[2];

    /* NB: splitting the i-loop in chunks is not optimal, as the first process
    will have to compute more forces than the others (its j-loop will be longer)
    To make even it out we split sys_natoms in two chunks and assign to the first 
    process the first chunk of the first half and the last of the second half,
    to the second process the second chunk of the first and the second to last of the second half etc..*/

    chunk[0] = sys->natoms/2 + sys->natoms%2;
    chunk[1] = sys->natoms - chunk[0];

    isize[0] = chunk[0]/sys->nsize;
    if ( sys->mpirank < chunk[0]%sys->nsize) {
        isize[0]++;
        sys->istart[0] = sys->mpirank*isize[0];
        sys->iend[0]   = sys->istart[0] + isize[0];
    } else {
        sys->istart[0] = sys->mpirank*isize[0] + chunk[0]%sys->nsize;
        sys->iend[0]   = sys->istart[0] + isize[0];
    }
    isize[1] = chunk[1]/sys->nsize;
    if ( (sys->nsize - sys->mpirank - 1) < chunk[1]%sys->nsize) {
        isize[1]++;
        sys->istart[1] = (sys->nsize - sys->mpirank - 1)*isize[1] + chunk[0];
        sys->iend[1]   = sys->istart[1] + isize[1];
    } else {
        sys->istart[1] = (sys->nsize - sys->mpirank - 1)*isize[1] + chunk[1]%sys->nsize + chunk[0];
        sys->iend[1]   = sys->istart[1] + isize[1];
    }
}
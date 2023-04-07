#include <mdlib.h>
#include <stdlib.h>

/* cell related functions */

/* initialize cells: get cell index and list of neighbours*/
void initialize_cells(mdsys_t *sys, cells_t *cells)
{   
    int i,j,k,ii,jj,kk;
    int ncell_x = sys->ncell_x;
    int ncell   = sys->ncell;

    // get cells index and list of neighbors
    for (i=0; i<ncell_x; i++) {
        for (j=0; j<ncell_x; j++) {
            for (k=0; k<ncell_x; k++) {
                // compute cell index
                int cell_index = i + j*ncell_x + k*ncell_x*ncell_x;
                cells[cell_index].id = cell_index;
                // initialize the list of neighbors
                int ncount = 0;
                cells[cell_index].n_list = (int *) malloc(27*sizeof(int));
                // loop over all neighbors (including itself, with periodic boundary conditions)
                for (ii=-1; ii<=1; ii++) {
                    for (jj=-1; jj<=1; jj++) {
                        for (kk=-1; kk<=1; kk++) {
                            // compute neighbor cell index
                            int n_x = (i+ii+ncell_x)%ncell_x;
                            int n_y = (j+jj+ncell_x)%ncell_x;
                            int n_z = (k+kk+ncell_x)%ncell_x;
                            int cell_index_neighbor = n_x + n_y*ncell_x + n_z*ncell_x*ncell_x;
                            // add neighbor to the list of neighbors
                            cells[cell_index].n_list[ncount] = cell_index_neighbor;
                            ncount++;
                        }
                    }
                }// end loop over neighbors
            }
        }
    }// end loop over cells

    // allocate memory for the list of atoms in each cell
    for (i=0; i<ncell; i++) {
        cells[i].a_list = (int *) malloc(sys->natoms*sizeof(int));
    }
}




/* initialize pairs: get the list of neighbouring cells */
void initialize_pairs(mdsys_t *sys, cells_t *cells) {
    
    int i,j;
    int ncell   = sys->ncell;
    
    sys->npairs = 0;
    sys->pairs  = (int *) malloc(2*(ncell*14)*sizeof(int));  // exactly 14 unique neighbors per cell (including self)
    
    // loop over all cells and add neighbor pairs
    for (i=0; i<ncell; i++) {
        for (j=0; j<27; j++) {
            int cell_index_neighbor = cells[i].n_list[j];
            // condition ensures that each pair is added only once (i.e. (i,j) and not (j,i)), including self
            if (cell_index_neighbor > i) {
                sys->pairs[2*sys->npairs] = i;
                sys->pairs[2*sys->npairs+1] = cell_index_neighbor;
                sys->npairs++;
            }
        }
    }
}


/* distribute atoms amongst cells */
void distribute_atoms(mdsys_t *sys, cells_t *cells) {

    int i;
    int ncell_x      = sys->ncell_x;
    double cell_size = sys->cell_size;
    double box_size  = sys->box;
    double rx,ry,rz;

    /* reinitialize cells */
    for(i=0; i<sys->ncell; i++) {
        cells[i].catoms = 0;
    }


    /* distribute atoms in each cell */     // this cannot really be parallelized
    for (i=0; i<sys->natoms; i++) {

        // get atom coordinates
        rx = sys->rx[i];    ry = sys->ry[i];    rz = sys->rz[i];
        // take care of negative coordinates
        while (rx < 0) rx += box_size;
        while (ry < 0) ry += box_size;
        while (rz < 0) rz += box_size;

        // compute cell index, with periodic boundary conditions (atoms might go out of the box...)
        int a_x = ((int) (rx/cell_size))%ncell_x;
        int a_y = ((int) (ry/cell_size))%ncell_x;
        int a_z = ((int) (rz/cell_size))%ncell_x;
        // compute cell index
        int cell_index = a_x + a_y*ncell_x + a_z*ncell_x*ncell_x;
        
        // add atom to the list of atoms in the cell, update counter
        cells[cell_index].a_list[cells[cell_index].catoms] = i;
        cells[cell_index].catoms++;
    }

}
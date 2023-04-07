#include <mdlib.h>
#include <math.h>


/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}


/* compute forces */
void force(mdsys_t *sys, cells_t *cells)
{   
    int i,j,k,i_mpi;
    double r,ffac;
    double rx,ry,rz;

    //new variables for the optimization
    double c12,c6,rcsq, rsq, rm6, rm2;
    // variable for potential energy amongst mpi processes
    double epot_tmp=0.0;

    /* set master energy and forces to zero */
    if (sys->mpirank == 0) {
        sys->epot=0.0;
        azzero(sys->fx_mpi,sys->natoms);
        azzero(sys->fy_mpi,sys->natoms);
        azzero(sys->fz_mpi,sys->natoms);
    }

    /* set mpi forces to zero */
    azzero(sys->fx_mpi,sys->natoms);
    azzero(sys->fy_mpi,sys->natoms);
    azzero(sys->fz_mpi,sys->natoms);

    /* broadcast positions from master */
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

    //remove the expensive sqrt(), pow() and division functions from the inner loop
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0); //c12 is the 12th power of sigma
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0); //c6 is the 6th power of sigma
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius


    /* INTER CELLS FORCES AND ENERGY */

    /* loop over all pairs of cells, distribute among processes */
    for (i_mpi=0; i_mpi<sys->npairs; i_mpi=i_mpi+sys->mpisize) {

        i = i_mpi + sys->mpirank;
        if (i >= sys->natoms) break;

        /* loop over all particles in first cell */
        for (j=0; j<cells[sys->pairs[2*i]].catoms; j++) {
            /* loop over all particles in second cell */
            for (k=0; k<cells[sys->pairs[2*i+1]].catoms; k++) {
                
                /* particles have no interactions with themselves */
                if (cells[sys->pairs[2*i]].a_list[j] == cells[sys->pairs[2*i+1]].a_list[k]) continue;
            
                /* compute distance between particles j and k */
                rx=pbc(sys->rx[cells[sys->pairs[2*i]].a_list[j]]-sys->rx[cells[sys->pairs[2*i+1]].a_list[k]], 0.5*sys->box);
                ry=pbc(sys->ry[cells[sys->pairs[2*i]].a_list[j]]-sys->ry[cells[sys->pairs[2*i+1]].a_list[k]], 0.5*sys->box);
                rz=pbc(sys->rz[cells[sys->pairs[2*i]].a_list[j]]-sys->rz[cells[sys->pairs[2*i+1]].a_list[k]], 0.5*sys->box);
                rsq=rx*rx+ry*ry+rz*rz;
                
                /* compute force and energy if within cutoff */
                if (rsq<rcsq) {
                    rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                    ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;

                    // compute potential energy
                    epot_tmp += rm6*(c12*rm6 - c6);

                    // compute forces (applying Newton's third law)
                    sys->fx_mpi[cells[sys->pairs[2*i]].a_list[j]] += ffac*rx;
                    sys->fy_mpi[cells[sys->pairs[2*i]].a_list[j]] += ffac*ry;
                    sys->fz_mpi[cells[sys->pairs[2*i]].a_list[j]] += ffac*rz;
                    sys->fx_mpi[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*rx;
                    sys->fy_mpi[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*ry;
                    sys->fz_mpi[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*rz;
                }
            }
        }       
    }
    


    /* INTRA CELLS FORCES AND ENERGY */
    
    // loop over all cells, distribute among processes */ 
    for (i_mpi=0; i_mpi<sys->ncell; i_mpi=i_mpi+sys->mpisize) {

        i = i_mpi + sys->mpirank;
        if (i >= sys->natoms) break;

        // loop over all particles in cell i
        for (j=0; j<cells[i].catoms; j++) {    
            for (k=j+1; k<cells[i].catoms; k++) {    
            
                /* compute distance between particles j and k */
                rx=pbc(sys->rx[cells[i].a_list[j]]-sys->rx[cells[i].a_list[k]], 0.5*sys->box);
                ry=pbc(sys->ry[cells[i].a_list[j]]-sys->ry[cells[i].a_list[k]], 0.5*sys->box);
                rz=pbc(sys->rz[cells[i].a_list[j]]-sys->rz[cells[i].a_list[k]], 0.5*sys->box);
                rsq=rx*rx+ry*ry+rz*rz;
                
                /* compute force and energy if within cutoff */
                if (rsq<rcsq) {
                    rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                    ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;

                    // compute potential energy
                    epot_tmp += rm6*(c12*rm6 - c6);

                    // compute forces (applying Newton's third law)
                    sys->fx_mpi[cells[i].a_list[j]] += ffac*rx;
                    sys->fy_mpi[cells[i].a_list[j]] += ffac*ry;
                    sys->fz_mpi[cells[i].a_list[j]] += ffac*rz;
                    sys->fx_mpi[cells[i].a_list[k]] -= ffac*rx;
                    sys->fy_mpi[cells[i].a_list[k]] -= ffac*ry;
                    sys->fz_mpi[cells[i].a_list[k]] -= ffac*rz;
                }
            }
        }       
    }

    /* reduce the forces and epot across processes */
    /* master needs them for velocity verlet */
    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fy_mpi, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fz_mpi, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot_tmp, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}
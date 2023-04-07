#include <mdlib.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


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
    int i,j,k,l,i_mpi,j_tds;
    double rm6,rm2,rsq,ffac;
    double rx,ry,rz;

    //new variables for the optimization
    double c12,c6,rcsq;
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
        if (i >= sys->npairs) break;
        int cell1 = sys->pairs[2*i];
        int cell2 = sys->pairs[2*i+1];

        /* open parallel region for first cell forces */
#pragma omp parallel default(shared) firstprivate(i,cell1,cell2) private(j,j_tds,k,l,rx,ry,rz,rm6,rm2,rsq,ffac) reduction(+:epot_tmp)
        {
            
            /* Thread-distribute loop over all particles in first cell */
            for (j_tds=0; j_tds<cells[cell1].catoms; j_tds=j_tds+sys->nthreads) {

                j = j_tds + omp_get_thread_num();
                if(j >= cells[cell1].catoms) continue;

                /* loop over all particles in second cell */
                for (k=0; k<cells[cell2].catoms; k++) {
                    
                    /* particles have no interactions with themselves */
                    if (cells[cell1].a_list[j] == cells[cell2].a_list[k]) continue;
                
                    /* compute distance between particles j and k */
                    rx=pbc(sys->rx[cells[cell1].a_list[j]]-sys->rx[cells[cell2].a_list[k]], 0.5*sys->box);
                    ry=pbc(sys->ry[cells[cell1].a_list[j]]-sys->ry[cells[cell2].a_list[k]], 0.5*sys->box);
                    rz=pbc(sys->rz[cells[cell1].a_list[j]]-sys->rz[cells[cell2].a_list[k]], 0.5*sys->box);
                    rsq=rx*rx+ry*ry+rz*rz;
                    
                    /* compute force and energy if within cutoff */
                    if (rsq<rcsq) {
                        rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                        ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;

                        // compute potential energy
                        epot_tmp += 0.5*rm6*(c12*rm6 - c6);

                        // // distribute forces to first cell
                        sys->fx_mpi[cells[cell1].a_list[j]] += rx*ffac;
                        sys->fy_mpi[cells[cell1].a_list[j]] += ry*ffac;
                        sys->fz_mpi[cells[cell1].a_list[j]] += rz*ffac;
                    }
                }
            }//end of thread-distributed loop

        }//end of parallel region


        /* open parallel region for second cell forces */
#pragma omp parallel default(shared) firstprivate(i,cell1,cell2) private(j,j_tds,k,l,rx,ry,rz,rm6,rm2,rsq,ffac) reduction(+:epot_tmp)
        {
            
            /* Thread-distribute loop over all particles in first cell */
            for (j_tds=0; j_tds<cells[cell2].catoms; j_tds=j_tds+sys->nthreads) {

                j = j_tds + omp_get_thread_num();
                if(j >= cells[cell2].catoms) continue;

                /* loop over all particles in second cell */
                for (k=0; k<cells[cell1].catoms; k++) {
                    
                    /* particles have no interactions with themselves */
                    if (cells[cell2].a_list[j] == cells[cell1].a_list[k]) continue;
                
                    /* compute distance between particles j and k */
                    rx=pbc(sys->rx[cells[cell2].a_list[j]]-sys->rx[cells[cell1].a_list[k]], 0.5*sys->box);
                    ry=pbc(sys->ry[cells[cell2].a_list[j]]-sys->ry[cells[cell1].a_list[k]], 0.5*sys->box);
                    rz=pbc(sys->rz[cells[cell2].a_list[j]]-sys->rz[cells[cell1].a_list[k]], 0.5*sys->box);
                    rsq=rx*rx+ry*ry+rz*rz;
                    
                    /* compute force and energy if within cutoff */
                    if (rsq<rcsq) {
                        rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                        ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;

                        // compute potential energy
                        epot_tmp += 0.5*rm6*(c12*rm6 - c6);

                        // // distribute forces to first cell
                        sys->fx_mpi[cells[cell2].a_list[j]] += rx*ffac;
                        sys->fy_mpi[cells[cell2].a_list[j]] += ry*ffac;
                        sys->fz_mpi[cells[cell2].a_list[j]] += rz*ffac;
                    }
                }
            }//end of thread-distributed loop

        }//end of parallel region

    } // end of loop over all pairs of cells


// NB: we had to close the first parallel region and open a new one:
// the first one was for the loop over all pairs of cells, the second
// one is for the loop over all particles in the first cell
// as we have no smart way to control pairs and single cells dependencies
// it could very well be a "slow" process from above and a "fast" one from below
// would try to access the same memory location at the same time



    /* INTRA CELLS FORCES AND ENERGY */
    
    // loop over all cells, distribute among processes */ 
    for (i_mpi=0; i_mpi<sys->ncell; i_mpi=i_mpi+sys->mpisize) {

        i = i_mpi + sys->mpirank;
        if (i >= sys->ncell) break;

        /* open parallel region */
#pragma omp parallel default(shared) firstprivate(i) private(j,j_tds,k,l,rx,ry,rz,rm6,rm2,rsq,ffac) reduction(+:epot_tmp)
        {
            // loop over all particles in cell i, distribute among threads
            for (j_tds=0; j_tds<cells[i].catoms; j_tds=j_tds+sys->nthreads) {

                j = j_tds + omp_get_thread_num();
                if(j >= cells[i].catoms) continue;

                for (k=0; k<cells[i].catoms; k++) {    

                    /* particles have no interactions with themselves */
                    if (j == k) continue;

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
                        epot_tmp += 0.5*rm6*(c12*rm6 - c6);

                        sys->fx_mpi[3*j]   += rx*ffac;
                        sys->fy_mpi[3*j+1] += ry*ffac;
                        sys->fz_mpi[3*j+2] += rz*ffac;
                    }
                }
            }       

        }//end of parallel region

    }// end of loop over cells


    /* reduce the forces and epot across processes */
    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fy_mpi, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fz_mpi, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot_tmp, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}
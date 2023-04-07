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

        /* open parallel region */
#pragma omp parallel default(shared) firstprivate(i,cell1,cell2) private(j,j_tds,k,l,rx,ry,rz,rm6,rm2,rsq,ffac) reduction(+:epot_tmp)
        {
            // each thread owns a copy of forces vector
            double *omp_forces1, *omp_forces2;
            omp_forces1 = (double*)malloc(3*(cells[cell1].catoms)*sizeof(double));
            omp_forces2 = (double*)malloc(3*(cells[cell2].catoms)*sizeof(double));
            azzero(omp_forces1, 3*(cells[cell1].catoms));
            azzero(omp_forces2, 3*(cells[cell2].catoms));

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
                        epot_tmp += rm6*(c12*rm6 - c6);

                        // distribute forces to first cell
                        omp_forces1[3*j]     += rx*ffac;
                        omp_forces1[3*j+1]   += ry*ffac;
                        omp_forces1[3*j+2]   += rz*ffac;
                        // distribute forces to second cell
                        omp_forces2[3*k]   -= rx*ffac;
                        omp_forces2[3*k+1] -= ry*ffac;
                        omp_forces2[3*k+2] -= rz*ffac;
                    }
                }
            }//end of thread-distributed loop

            // reduce forces vector to the master thread for first cell
            #pragma omp critical
            for (l=0; l<cells[cell1].catoms; l++) {
                sys->fx_mpi[cells[cell1].a_list[l]] += omp_forces1[3*l];
                sys->fy_mpi[cells[cell1].a_list[l]] += omp_forces1[3*l+1];
                sys->fz_mpi[cells[cell1].a_list[l]] += omp_forces1[3*l+2];
            }

            // reduce forces vector to the master thread for second cell
            #pragma omp critical
            for (l=0; l<cells[cell2].catoms; l++) {
                sys->fx_mpi[cells[cell2].a_list[l]] += omp_forces2[3*l];
                sys->fy_mpi[cells[cell2].a_list[l]] += omp_forces2[3*l+1];
                sys->fz_mpi[cells[cell2].a_list[l]] += omp_forces2[3*l+2];
            }

            // free thread's forces vector
            free(omp_forces1);
            free(omp_forces2);

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
            // each thread owns a copy of forces vector for the cell i
            double *omp_forces;
            omp_forces = (double*)malloc(3*(cells[i].catoms)*sizeof(double));
            azzero(omp_forces, 3*(cells[i].catoms));

            // loop over all particles in cell i, distribute among threads
            for (j_tds=0; j_tds<cells[i].catoms; j_tds=j_tds+sys->nthreads) {

                j = j_tds + omp_get_thread_num();
                if(j >= cells[i].catoms) continue;

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

                        omp_forces[3*j]   += rx*ffac;
                        omp_forces[3*j+1] += ry*ffac;
                        omp_forces[3*j+2] += rz*ffac;
                        omp_forces[3*k]   -= rx*ffac;
                        omp_forces[3*k+1] -= ry*ffac;
                        omp_forces[3*k+2] -= rz*ffac;
                    }
                }
            }       

            // reduce forces vector to the master thread
            #pragma omp critical
            for (l=0; l<cells[i].catoms; l++) {
                sys->fx_mpi[cells[i].a_list[l]] += omp_forces[3*l];
                sys->fy_mpi[cells[i].a_list[l]] += omp_forces[3*l+1];
                sys->fz_mpi[cells[i].a_list[l]] += omp_forces[3*l+2];
            }
            // free thread's forces vector
            free(omp_forces); 

        }//end of parallel region

    }// end of loop over cells


    /* reduce the forces and epot across processes */
    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fy_mpi, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fz_mpi, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot_tmp, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}
#include <mdlib-util.h>
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
void force(mdsys_t *sys)
{
    double ffac;
    double rx,ry,rz;
    int i,j,i_mpi;
    //new variables for the optimization
    double c12,c6,rcsq, rsq, rm6, rm2;
    // variable for potential energy amongst mpi processes
    double epot=0.0;

    /* set master energy and forces to zero */
    if (sys->mpirank == 0) {
        sys->epot=0.0;
        azzero(sys->fx,sys->natoms);
        azzero(sys->fy,sys->natoms);
        azzero(sys->fz,sys->natoms);
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

    for(i=0; i < (sys->natoms); i += sys->nsize) {

        i_mpi = i + sys->mpirank;
        if (i_mpi >= sys->natoms) break;

        for(j=0; j < (sys->natoms); ++j) { 
            
            /* particles have no interactions with themselves */
            if (i_mpi==j) continue; 

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i_mpi] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i_mpi] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i_mpi] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz; 

            /* compute force and energy if within cutoff */
           if (rsq < rcsq) {
                //remove the expensive pow() and division functions from the inner loop
                double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                epot += 0.5*rm6*(c12*rm6 - c6); 

                sys->fx_mpi[i_mpi] += rx*ffac;
                sys->fy_mpi[i_mpi] += ry*ffac;
                sys->fz_mpi[i_mpi] += rz*ffac;
            }
        }
    }

    MPI_Reduce(sys->fx_mpi, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fy_mpi, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->fz_mpi, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);

}
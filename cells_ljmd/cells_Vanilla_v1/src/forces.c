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
    int i,j,k;
    double r,ffac;
    double rx,ry,rz;

    //new variables for the optimization
    double c12,c6,rcsq, rsq, rm6, rm2;
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    //remove the expensive sqrt(), pow() and division functions from the inner loop
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0); //c12 is the 12th power of sigma
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0); //c6 is the 6th power of sigma
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius


    /* INTER CELLS FORCES AND ENERGY */

    /* loop over all pairs of cells */
    for (i=0; i<sys->npairs; i++) {
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
                    sys->epot += rm6*(c12*rm6 - c6);

                    // compute forces (applying Newton's third law)
                    sys->fx[cells[sys->pairs[2*i]].a_list[j]] += ffac*rx;
                    sys->fy[cells[sys->pairs[2*i]].a_list[j]] += ffac*ry;
                    sys->fz[cells[sys->pairs[2*i]].a_list[j]] += ffac*rz;
                    sys->fx[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*rx;
                    sys->fy[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*ry;
                    sys->fz[cells[sys->pairs[2*i+1]].a_list[k]] -= ffac*rz;
                }
            }
        }       
    }
    


    /* INTRA CELLS FORCES AND ENERGY */
    
    // loop over all cells 
    for (i=0; i<sys->ncell; i++) {

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
                    sys->epot += rm6*(c12*rm6 - c6);

                    // compute forces (applying Newton's third law)
                    sys->fx[cells[i].a_list[j]] += ffac*rx;
                    sys->fy[cells[i].a_list[j]] += ffac*ry;
                    sys->fz[cells[i].a_list[j]] += ffac*rz;
                    sys->fx[cells[i].a_list[k]] -= ffac*rx;
                    sys->fy[cells[i].a_list[k]] -= ffac*ry;
                    sys->fz[cells[i].a_list[k]] -= ffac*rz;
                }
            }
        }       
    }
}
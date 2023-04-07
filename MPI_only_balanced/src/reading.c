#include <stdlib.h>
#include <mdlib-util.h>
#include <mdlib.h>

int readinput (mdsys_t *sys, int * nprint, char restfile[BLEN], char trajfile[BLEN], char ergfile[BLEN]) {

    char line[BLEN];
    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys->natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys->mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys->epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys->sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys->rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys->box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys->nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys->dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    *nprint=atoi(line);

    return 0;
}


/* read restfile */
int readrest (mdsys_t *sys, char restfile[BLEN]) {

    int i;
    FILE *fp;

    fp=fopen(restfile,"r");

    if(fp) {
        for (i=0; i<sys->natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys->rx+i, sys->ry+i, sys->rz+i);
        }
        for (i=0; i<sys->natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys->vx+i, sys->vy+i, sys->vz+i);
        }
        fclose(fp);

        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    return 0;
}
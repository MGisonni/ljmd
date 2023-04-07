#include <sys/time.h>

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: get current time in seconds since epoch */
double wallclock()
{
    struct timeval t;
    gettimeofday(&t,0);
    return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
}

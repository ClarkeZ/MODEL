#include "algo.h"

double wtime(){
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* generate a random floating point number from min to max */
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

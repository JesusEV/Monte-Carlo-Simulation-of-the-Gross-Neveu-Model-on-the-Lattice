#ifndef GAUSS_RAN_GEN_H
#define GAUSS_RAN_GEN_H

#include <cmath>

extern "C" 
{
    #include "ranlxd.h"
}

void gaussian_random(double* gauss_rands, const int M_dim, const int seed);

#endif
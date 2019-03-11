#ifndef DET_WRAPP_H
#define DET_WRAPP_H

#include <complex>
#include <cmath>

extern "C" 
{
    #include "suiteSparseDet.h"
}

#include <sparseDiracMatrixClass.h>

typedef std::complex<double> dcmplx;

dcmplx determinant(SparseDiracMatrix& diracMatrix);

#endif
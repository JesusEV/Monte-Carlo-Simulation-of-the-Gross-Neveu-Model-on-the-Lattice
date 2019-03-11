#ifndef IMP_SAMP_INT_H
#define IMP_SAMP_INT_H

#include <string>
#include <iostream>
#include <chrono>
#include <sparseDiracMatrixClass.h>
#include <gaussianRandomGenerator.h>
#include <determinantWrapper.h>
#include <observables.h>

using namespace std::chrono;
using time_point = high_resolution_clock::time_point;

typedef std::complex<double> dcmplx;

dcmplx importSamplIntegrator(SparseDiracMatrix& diracMatrix, 
                             const std::string& observable,
                             const int MC_steps,
                             const int npes, const int rank);

#endif
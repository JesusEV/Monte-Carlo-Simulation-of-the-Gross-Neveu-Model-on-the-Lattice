#ifndef SUI_SPA_DET_H
#define SUI_SPA_DET_H

#include <math.h>
#include <suitesparse/umfpack.h>
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

void* suiteSparseDet(double* determinant,
                     const int M_dim, 
                     const int non_zero_entries,
                     const double* const M_valsRe, 
                     const double* const M_valsIm,
                     const int* const M_xs, 
                     const int* const M_ys, 
                     double* Ax, double* Az,
                     int* Ap, int* Ai);

#endif
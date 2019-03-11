#include <determinantWrapper.h>

dcmplx determinant(SparseDiracMatrix& diracMatrix)
{
    dcmplx det;
    double determinant_pair[2] = {0,0};

    suiteSparseDet(determinant_pair,
                   diracMatrix.get_M_dim(),
                   diracMatrix.getSize(),
                   diracMatrix.get_vals_Re(),
                   diracMatrix.get_vals_Im(),
                   diracMatrix.get_xs(),
                   diracMatrix.get_ys(),
                   diracMatrix.get_Ax(),
                   diracMatrix.get_Az(),
                   diracMatrix.get_Ap(),
                   diracMatrix.get_Ai());
        
    det = {determinant_pair[0], determinant_pair[1]};

    return det;
}


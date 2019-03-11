#include <suiteSparseDet.h>

void* suiteSparseDet(double* determinant,
                     const int M_dim, 
                     const int non_zero_entries,
                     const double* const M_valsRe, 
                     const double* const M_valsIm,
                     const int* const M_xs, 
                     const int* const M_ys, 
                     double* Ax, double* Az,
                     int* Ap, int* Ai)
{
    static double                       x[5], z[5], r[5];
    double                              determinantRe, determinantIm;
    double                              Control [UMFPACK_CONTROL];
    double                              Info [UMFPACK_INFO];
    void                                *Symbolic, *Numeric;
    int                                 status;

    // unpack triplets
    status = umfpack_zi_triplet_to_col (M_dim, M_dim, 
                                        non_zero_entries, 
                                        M_xs, 
                                        M_ys, 
                                        M_valsRe, 
                                        M_valsIm, 
                                        Ap, 
                                        Ai, 
                                        Ax, 
                                        Az, 
                                        (int *) NULL);

    // symbolic factorization 
    status = umfpack_zi_symbolic(M_dim, M_dim, Ap, Ai, Ax, Az, &Symbolic, Control, Info);  
    // // numeric factorization 
    status = umfpack_zi_numeric(Ap, Ai, Ax, Az, Symbolic, &Numeric, Control, Info);
    // // compute the determinant
    status = umfpack_zi_get_determinant(x, z, r, Numeric, Info);
   
    determinantRe = x[0]*pow(10.0, r[0]);
    determinantIm = z[0]*pow(10.0, r[0]);
    
    determinant[0] = determinantRe;
    determinant[1] = determinantIm;

    umfpack_di_free_symbolic (&Symbolic) ;
    umfpack_di_free_numeric (&Numeric) ;
}

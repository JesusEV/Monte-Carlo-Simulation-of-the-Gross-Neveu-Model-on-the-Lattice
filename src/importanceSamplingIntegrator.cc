#include <importanceSamplingIntegrator.h>

dcmplx importSamplIntegrator(SparseDiracMatrix& diracMatrix, 
                             const std::string& observable,
                             const int MC_steps,
                             const int npes, const int rank)
{
    int                         MCstart, MCend; 
    int                         nloc, rem; 
    int                         M_dim = diracMatrix.getSize();
    int                         TMP = M_dim/2;
    double                      (*obsrv)(const int, const double* const);
    dcmplx                      tot_sum = {0., 0.};
    dcmplx                      det = {0., 0.};
    std::unique_ptr<double[]>   gauss_rands{new double[TMP]{}};

    // selection of observable function pointer
    if (observable == "mean_scalar") { obsrv = &mean_scalar; }
    else if (observable == "mean_sq_scalar") { obsrv = &mean_sq_scalar; }
    else {obsrv = &non;}

    nloc = MC_steps / npes;
    rem  = MC_steps % npes;

    // fair distribution of the remainder among ranks:
    MCstart = rank * nloc  + ((rem-rank+npes-1)/npes)*(rank-rem) + rem;
    MCend   = MCstart + nloc + (rem-rank+npes-1)/npes;

    for (int step{MCstart}; step < MCend; ++step)
    {
        gaussian_random(gauss_rands.get(), TMP, step+1);
        diracMatrix.update_diagonal(gauss_rands.get());
        det = determinant(diracMatrix);     
        tot_sum += det*obsrv(TMP, gauss_rands.get());
    }

    return tot_sum;
}   

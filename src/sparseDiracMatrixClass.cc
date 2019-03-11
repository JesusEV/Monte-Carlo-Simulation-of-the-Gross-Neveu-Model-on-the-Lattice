#include <sparseDiracMatrixClass.h>

void SparseDiracMatrix::build_dirac_triplets()
{

    int                         i, j, k, cont = 0;
    int                         a_p_b_c;
    int                         N, M, n_1, n_2, m_1, m_2;
    int                         r_0_1, r_0_2, r_1_1, r_1_2;
    int                         r_2_1, r_2_2, r_3_1, r_3_2;
    int                         r_4_1, r_4_2;
    int                         P0, Q0, P1, Q1, P2, Q2;
    int                         P3, Q3, P4, Q4; 
    dcmplx                      I(0.0,1.0);
    dcmplx                      S[20];
    dcmplx                      A00, A01, A10, A11;

    for (n_1 = 0; n_1 < num_cols; n_1++)
    {
        for (n_2 = 0; n_2 < num_rows; n_2++)
        {
            for (m_1 = 0; m_1 < num_cols; m_1++)
            {
                for (m_2 = 0; m_2 < num_rows; m_2++)
                {
                    for(i=0; i < 20; i++) {S[i] = 0;} 

                    N = n_2*num_cols+n_1;
                    M = m_2*num_cols+m_1;

                    //  % operators to take into account boundary conds
                    r_0_1 = (n_1)%num_cols;
                    if (r_0_1 < 0) r_0_1 += num_cols;                   
                    r_0_2 = (n_2)%num_rows;
                    if (r_0_2 < 0) r_0_2 += num_rows;

                    r_1_1 = (n_1+1)%num_cols;
                    if (r_1_1 < 0) r_1_1 += num_cols;
                    r_1_2 = (n_2)%num_rows;
                    if (r_1_2 < 0) r_1_2 += num_rows;

                    r_2_1 = (n_1-1)%num_cols;
                    if (r_2_1 < 0) r_2_1 += num_cols;                   
                    r_2_2 = (n_2)%num_rows;
                    if (r_2_2 < 0) r_2_2 += num_rows;

                    r_3_1 = (n_1)%num_cols;
                    if (r_3_1 < 0) r_3_1 += num_cols;                   
                    r_3_2 = (n_2+1)%num_rows;
                    if (r_3_2 < 0) r_3_2 += num_rows;

                    r_4_1 = (n_1)%num_cols;
                    if (r_4_1 < 0) r_4_1 += num_cols;                   
                    r_4_2 = (n_2-1)%num_rows;
                    if (r_4_2 < 0) r_4_2 += num_rows;

                    P0 = r_0_2*num_cols+r_0_1;
                    Q0 = m_2*num_cols+m_1;
                    
                    P1 = r_1_2*num_cols+r_1_1;
                    Q1 = m_2*num_cols+m_1;

                    P2 = r_2_2*num_cols+r_2_1;
                    Q2 = m_2*num_cols+m_1;

                    P3 = r_3_2*num_cols+r_3_1;
                    Q3 = m_2*num_cols+m_1;

                    P4 = r_4_2*num_cols+r_4_1;
                    Q4 = m_2*num_cols+m_1;


                    if(P0 == Q0){
                        S[0] = 0.; /* 00 entry */
                        S[1] = 0.; /* 01 entry */
                        S[2] = 0.; /* 10 entry */
                        S[3] = 0.; /* 11 entry */ 
                    } /* P0 == Q0 */

                    if(P1 == Q1){
                        a_p_b_c = 1;
                        if (n_1 == num_cols-1) a_p_b_c = -1; 
                        S[4] = 0.5*a_p_b_c; /* 00 entry */
                        S[5] = 0.5*a_p_b_c; /* 01 entry */
                        S[6] = 0.5*a_p_b_c; /* 10 entry */
                        S[7] = 0.5*a_p_b_c; /* 11 entry */ 
                    } /* P0 == Q0 */

                    if(P2 == Q2){
                        a_p_b_c = 1;
                        if (n_1 == 0) a_p_b_c = -1;                         
                        S[8] = 0.5*a_p_b_c; /* 00 entry */
                        S[9] = -0.5*a_p_b_c; /* 01 entry */
                        S[10] = -0.5*a_p_b_c; /* 10 entry */
                        S[11] = 0.5*a_p_b_c; /* 11 entry */ 
                    } /* P0 == Q0 */

                    if(P3 == Q3){
                        S[12] = 0.5; /* 00 entry */
                        S[13] = -0.5*I; /* 01 entry */
                        S[14] = 0.5*I; /* 10 entry */
                        S[15] = 0.5; /* 11 entry */ 
                    } /* P0 == Q0 */

                    if(P4 == Q4){
                        S[16] = 0.5; /* 00 entry */
                        S[17] = 0.5*I; /* 01 entry */
                        S[18] = -0.5*I; /* 10 entry */
                        S[19] = 0.5; /* 11 entry */ 
                    } /* P0 == Q0 */


                    A00 = S[0] - S[4] - S[8] - S[12] - S[16];
                    A01 = S[1] - S[5] - S[9] - S[13] - S[17];
                    A10 = S[2] - S[6] - S[10] - S[14] - S[18];
                    A11 = S[3] - S[7] - S[11] - S[15] - S[19];


                    /* Create Triplet (COO) representation of the sparse
                     Matrix */
                    if (A00 != ZERO)
                    {
                        M_valsRe[cont] = std::real(A00);
                        M_valsIm[cont] = std::imag(A00);
                        M_xs[cont] = 2*N;
                        M_ys[cont] = 2*M;
                        cont++;
                    }

                    if (A01 != ZERO)
                    {
                        M_valsRe[cont] = std::real(A01);
                        M_valsIm[cont] = std::imag(A01);
                        M_xs[cont] = 2*N;
                        M_ys[cont] = 2*M+1;
                        cont++;
                    }

                    if (A10 != ZERO)
                    {
                        M_valsRe[cont] = std::real(A10);
                        M_valsIm[cont] = std::imag(A10);
                        M_xs[cont] = 2*N+1;
                        M_ys[cont] = 2*M;   
                        cont++;
                    }

                    if (A11 != ZERO)
                    {
                        M_valsRe[cont] = std::real(A11);
                        M_valsIm[cont] = std::imag(A11);
                        M_xs[cont] = 2*N+1 ;
                        M_ys[cont] = 2*M+1; 
                        cont++;
                    }

                    A00 = 0.;
                    A01 = 0.;
                    A10 = 0.;
                    A11 = 0.;

                    
                } /* m_2 on num_rows */ 
            } /* m_1 on num_cols */
        } /* n_2 on num_rows */
    } /* n_1 on num_cols */
} /* Triplets construction */

void SparseDiracMatrix::print_dirac_matrix(bool print)
{
    if (print)
    {
        double xv, yv;
        for (std::size_t i{0}; i <  M_dim; ++i)
        {
            for (std::size_t j{0}; j <  M_dim; ++j)
            {
                xv = 0;
                yv = 0;
                for (std::size_t k{0}; k <  non_zero_entries; ++k)
                {
                    if (M_xs[k] == i && M_ys[k] == j) 
                    {
                        xv = M_valsRe[k];               
                        yv = M_valsIm[k];               
                        break;
                    } // if
                } // for j
                std::cout << xv << "+" << yv << "i\t";      
            } // for i
            std::cout << std::endl;
        }   //for k 
        std::cout << std::endl;
    } // if
}

void SparseDiracMatrix::update_diagonal(double* gauss_rands)
{
    for (std::size_t cont{0}; cont < M_dim/2; cont++)
    {
        M_valsRe[(9-1)*M_dim + 2*cont] = 2 + mass +  g_cc * gauss_rands[cont];
        M_valsIm[(9-1)*M_dim + 2*cont] = 0.0;
        M_xs[(9-1)*M_dim + 2*cont] = 2*cont;
        M_ys[(9-1)*M_dim + 2*cont] = 2*cont;

        M_valsRe[(9-1)*M_dim + 2*cont + 1] = 2 + mass +  g_cc * gauss_rands[cont];
        M_valsIm[(9-1)*M_dim + 2*cont + 1] = 0.0;
        M_xs[(9-1)*M_dim + 2*cont + 1] = 2*cont + 1;
        M_ys[(9-1)*M_dim + 2*cont + 1] = 2*cont + 1;    

        // M_valsRe[(9-1)*M_dim + 2*cont] = 1;
        // M_valsIm[(9-1)*M_dim + 2*cont] = 0.0;
        // M_xs[(9-1)*M_dim + 2*cont] = 2*cont;
        // M_ys[(9-1)*M_dim + 2*cont] = 2*cont;

        // M_valsRe[(9-1)*M_dim + 2*cont + 1] = 1;
        // M_valsIm[(9-1)*M_dim + 2*cont + 1] = 0.0;
        // M_xs[(9-1)*M_dim + 2*cont + 1] = 2*cont + 1;
        // M_ys[(9-1)*M_dim + 2*cont + 1] = 2*cont + 1; 
    }
}

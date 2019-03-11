#ifndef DIRAC_MATRIX_H
#define DIRAC_MATRIX_H

#include <algorithm> 
#include <iostream>
#include <memory>
#include <complex>
#include <cmath>


typedef std::complex<double> dcmplx;
static dcmplx ZERO = {0.,0.};


class SparseDiracMatrix 
{
    private:

        int num_rows;
        int num_cols;
        double mass;
        double g_cc;
        int M_dim;
        int non_zero_entries;

        std::unique_ptr<double[]> M_valsRe;
        std::unique_ptr<double[]> M_valsIm;
        std::unique_ptr<int[]> M_xs;
        std::unique_ptr<int[]> M_ys;

        std::unique_ptr<double[]> Ax;
        std::unique_ptr<double[]> Az;
        std::unique_ptr<int[]> Ap;
        std::unique_ptr<int[]> Ai;

    public:
    // custom ctor
        explicit SparseDiracMatrix(const int rows, const int cols, \
                                  const double m, const double g) noexcept
                :   num_rows{rows}, num_cols{cols}, mass{m}, g_cc{g},
                    M_dim{(2*rows*cols)}, 
                    non_zero_entries{(2*rows*cols)*9}, 
                    M_valsRe{new double[(2*rows*cols)*9]{}}, 
                    M_valsIm{new double[(2*rows*cols)*9]{}}, 
                    M_xs{new int[(2*rows*cols)*9]{}}, 
                    M_ys{new int[(2*rows*cols)*9]{}},
                    Ax{new double[(2*rows*cols)*9]{}},               
                    Az{new double[(2*rows*cols)*9]{}},               
                    Ap{new int[(2*rows*cols) + 1]{}},                
                    Ai{new int[(2*rows*cols)*9]{}}               
        { build_dirac_triplets(); }

        ~SparseDiracMatrix() {} 

        const int getSize() const noexcept {return non_zero_entries;}
        int get_M_dim() const noexcept {return M_dim;}
        int* get_xs()  {return M_xs.get();}
        int* get_ys() {return M_ys.get();}
        double* get_vals_Re() {return M_valsRe.get();}
        double* get_vals_Im() {return M_valsIm.get();}

        int* get_Ap() {return Ap.get();}
        int* get_Ai() {return Ai.get();}
        double* get_Ax() {return Ax.get();}
        double* get_Az() {return Az.get();}
        
        void build_dirac_triplets();
        void update_diagonal(double* gauss_rands);
        void print_dirac_matrix(bool print);
};

#endif
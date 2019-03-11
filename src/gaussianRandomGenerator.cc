#include <gaussianRandomGenerator.h>

void gaussian_random(double* gauss_rands, 
                     const int M_dim, 
                     const int seed)
{
    int                     i{0};
    double                  uniform_rands[2];
    double                  u1, u2;
    double                  gamma0, gamma1;
    
    rlxd_init(1,seed);
    
    while(i < M_dim/2)
    {
        ranlxd(uniform_rands, 2);
        u1 = uniform_rands[0];
        u2 = uniform_rands[1];
        u1 = 2*u1-1;
        u2 = 2*u2-1;
        gamma0 = pow(u1, 2) + pow(u2, 2);
        
        if (gamma0 < 1 && gamma0 != 0)
        { 
            gamma1 = sqrt(-2*log(gamma0)/gamma0);    
            u1 = gamma1*u1;
            u2 = gamma1*u2;
            gauss_rands[2*i] = u1;
            gauss_rands[2*i+1] = u2;
            i += 1;
        }
    }
}   

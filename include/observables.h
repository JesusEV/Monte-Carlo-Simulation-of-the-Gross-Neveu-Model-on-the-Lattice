#ifndef OBSERVABLES_H
#define OBSERVABLES_H

inline double non(const int size, const double* const samples) { return 1.0; }


inline double mean_scalar(const int size, const double* const samples)
{
    double mean = 0.;
    for (std::size_t i{0}; i < size; ++i) 
    {
        mean += samples[i];
    }   
    return mean/size;
}


inline double mean_sq_scalar(const int size, const double* const samples)
{
    double mean = 0.;
    for (std::size_t i{0}; i < size; ++i) 
    {
        mean += samples[i]*samples[i];
    }   
    return mean/size;
}


#endif
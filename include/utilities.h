#ifndef UTILITIES_H
#define UTILITIES_H

#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>

// read input parameters
void read_input_parameters(const char* const file_name, 
                           int& rows, int& cols,
                           double& mass, double& g_cc, 
                           std::string& observable_name,  
                           int& MC_steps, bool& print_matrix, 
                           int& matrix_dim);

// check on input parameters
void check_parameters(const int argc);

#endif 
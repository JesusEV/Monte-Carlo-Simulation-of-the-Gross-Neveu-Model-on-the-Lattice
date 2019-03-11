#include <utilities.h>

void read_input_parameters(const char* const file_name, 
                           int& rows, int& cols,
                           double& mass, double& g_cc, 
                           std::string& observable_name,  
                           int& MC_steps, bool& print_matrix, 
                           int& matrix_dim)
{
    char obsrv[30];
    std::stringstream ss;

    std::ifstream input_file{file_name};
    input_file >> rows;
    input_file >> cols;
    input_file >> mass;
    input_file >> g_cc;
    input_file >> obsrv;
    input_file >> MC_steps;
    input_file >> print_matrix;
    matrix_dim = 2*rows*cols;
    ss << obsrv;
    ss >> observable_name; 
}

void check_parameters(const int argc)
{
    if(argc != 2) 
    {
        std::cout << "wrong number of arguments. ";
        std::cout << "Usage: ./gn input.inp" << std::endl;
        exit(1);
    }
}

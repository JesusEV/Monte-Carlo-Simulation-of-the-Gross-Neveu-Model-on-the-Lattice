#include <utilities.h>
#include <sparseDiracMatrixClass.h>
#include <importanceSamplingIntegrator.h>
#include <conditionalOStreamClass.h>
#include <mpiCommManager.h>

int main(int argc, char* argv[])
{
    int                             rows;
    int                             columns;
    double                          g_cc{0};
    double                          mass{0};
    int                             matrix_dim ;
    int                             MC_steps;
    bool                            print_matrix = false;
    std::string                     observable_name;
    std::string                     no_observable{"non"};
    dcmplx                          locNumerator, locDenominator;
    dcmplx                          numerator{0}, denominator{1};
    dcmplx                          result{0};
    extern bool                     thisIsWorld;
    extern bool                     thisIsNotWorld;
    int                             npes, rank;
    double                          t0, t1, t2, t3, t4;

    MPICommManager world(&argc, &argv, thisIsWorld);
    npes = world.getNprocs();
    rank = world.getRank();


    ConditionalOStream pcout{std::cout, rank == 0};


    check_parameters(argc);
    read_input_parameters(argv[1],  rows,                       \
                                    columns,                    \
                                    mass,                       \
                                    g_cc,                       \
                                    observable_name,            \
                                    MC_steps,                   \
                                    print_matrix,               \
                                    matrix_dim);        
    
    t0 = world.time();
    
    SparseDiracMatrix diracMatrix{rows, columns, mass, g_cc};

    t1 = world.time();

    locNumerator   = importSamplIntegrator(diracMatrix, observable_name, 
                                        MC_steps, npes, rank); 
    locDenominator = importSamplIntegrator(diracMatrix, no_observable, 
                                        MC_steps, npes, rank); 

    t2 = world.time();

    world.dcSumReduce(&locNumerator,   &numerator,   1, 0);
    world.dcSumReduce(&locDenominator, &denominator, 1, 0);

    t3 = world.time();

    result = numerator/denominator;
    
    pcout << "Result: " << result.real();
    pcout << " + " << result.imag() << "i ";
    pcout << "Processes: " << npes << " ";
    pcout << "M-constr-time: "   << t1 - t0 << " ";     
    pcout << "Comp-time: "   << t2 - t1 << " ";     
    pcout << "Comm-time: "   << t3 - t2 << " ";     
    pcout << "Total-time: "   << t3 - t0 << std::endl;     

    return 0;    
}


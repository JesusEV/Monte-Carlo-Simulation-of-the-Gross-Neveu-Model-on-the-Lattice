#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#include <mpi.h>

static bool thisIsWorld = true;
static bool thisIsNotWorld = false;

typedef std::complex<double> dcmplx;

class MPICommManager
{
    private:
        MPI_Comm comm;
        bool CommIsWorld;
        int npes = 0;
        int rank = 0;

    public:

        MPICommManager(int* pargc, char** pargv[], const bool isItWorld)
        : CommIsWorld{isItWorld}                // isItWorld ? T/F 
        {
            // MPI Initialization
            if (CommIsWorld)
            {
                MPI_Init(pargc, pargv);             
                comm = MPI_COMM_WORLD;
            }
            MPI_Comm_size(comm, &npes);
            MPI_Comm_rank( comm, &rank);    
        }

        ~MPICommManager() 
        {   
            if (CommIsWorld) {MPI_Finalize();}
        }


        const int getNprocs () const noexcept {return npes;} 
        const int getRank () const noexcept {return rank;}
        

        inline const double time() const {return MPI_Wtime();}


        void dSendRecv(const double* sendThis, double* recvHere,\
                       const int buffSize, const int sendTo, const int recvFrom) const
        {
            MPI_Status status;
            int error;
            error = MPI_Sendrecv(sendThis, buffSize, MPI_DOUBLE, sendTo  , 0,  \
                                 recvHere, buffSize, MPI_DOUBLE, recvFrom, 0,  \
                                 comm, &status);
        }   


        void dcSumReduce(const dcmplx* sendThis, dcmplx* recvHere, 
                         const int buffSize, const int sendTo) const
        {
            int error;
            error = MPI_Reduce(sendThis, recvHere, buffSize, MPI::DOUBLE_COMPLEX, 
                               MPI_SUM, sendTo, comm);
        }
};


#endif
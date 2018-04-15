//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MpiFunc.hpp"

#include <iostream>

namespace mpi { // mpi namespace

bool init(int argc, char** argv)
{
    if (ENABLE_MPI)
    {
        int32_t mlevel = 0;

        auto merror = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                                      &mlevel);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "init()");
        
            return false;
        }
    }

    return true;
}

bool finalize()
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Finalize();

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "finalize()"); 

            return false;
        }
    }

    return true;
}
    
int32_t rank(MPI_Comm comm)
{
    int32_t mrank = 0;

    if (ENABLE_MPI) MPI_Comm_rank(comm, &mrank);

    return mrank;
}

int32_t nodes(MPI_Comm comm)
{
    int32_t mnodes = 1;

    if (ENABLE_MPI) MPI_Comm_size(comm, &mnodes);

    return mnodes;
}
    
bool compare(MPI_Comm comm1, MPI_Comm comm2)
{
    if (ENABLE_MPI)
    {
        int32_t mcvalue = 0;

        auto merror = MPI_Comm_compare(comm1, comm2, &mcvalue);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "compare()");

            return false;
        }

        if (mcvalue == MPI_IDENT) return true;

        return false;
    }

    return true;
}
    
void bcast(int32_t& value, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "bcast(int32_t)");
    }
}

void bcast(double& value, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "bcast(double)");
    }
}

void bcast(bool& value, int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        int32_t mvalue = 0;
            
        if (rank == mpi::master()) mvalue = (value) ? 1 : 0;
            
        mpi::bcast(mvalue, comm);
            
        value = (mvalue == 1) ? true : false;
    }
}

void bcast(std::vector<int32_t>& vector, int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        int32_t veclen = 0;

        if (rank == mpi::master()) veclen = static_cast<int32_t>(vector.size());

        mpi::bcast(veclen, comm);

        for (int32_t i = 0; i < veclen; i++)
        {
            int32_t mvalue = 0;

            if (rank == mpi::master()) mvalue = vector[i];

            mpi::bcast(mvalue, comm);

            if (rank != mpi::master()) vector.push_back(mvalue);
        }
    }
}
    
// TODO: Add other broadcast methods...
    
void abort(const int errorcode, const char* label)
{
    if (ENABLE_MPI)
    {
        std::cerr << "MPI ERROR - " << label << " : " << errorcode << std::endl;
        
        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}
    
} // mpi namespace

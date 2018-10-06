//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MpiFunc.hpp"

#include <iostream>
#include <sstream>

namespace mpi { // mpi namespace

bool
init(int    argc,
     char** argv)
{
    if (ENABLE_MPI)
    {
        if (initialized())
        {
            return true;
        }

        int32_t mlevel = 0;

        auto merror = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                                      &mlevel);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::init()");
        
            return false;
        }
    }

    return true;
}

bool
initialized()
{
    if (ENABLE_MPI)
    {
        int32_t minit = 0;

        MPI_Initialized(&minit);
        
        if (minit == 1) return true;
    }

    return false;
}

bool
finalize()
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Finalize();

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::finalize()"); 

            return false;
        }
    }

    return true;
}
    
int32_t
rank(MPI_Comm comm)
{
    int32_t mrank = 0;

    if (ENABLE_MPI) MPI_Comm_rank(comm, &mrank);

    return mrank;
}

int32_t
nodes(MPI_Comm comm)
{
    int32_t mnodes = 1;

    if (ENABLE_MPI) MPI_Comm_size(comm, &mnodes);

    return mnodes;
}
    
bool
compare(MPI_Comm comm1,
        MPI_Comm comm2)
{
    if (ENABLE_MPI)
    {
        int32_t mcvalue = 0;

        auto merror = MPI_Comm_compare(comm1, comm2, &mcvalue);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::compare()");

            return false;
        }

        if (mcvalue == MPI_IDENT) return true;

        return false;
    }

    return true;
}
    
void
bcast(int32_t& value,
      MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcast(int32_t)");
    }
}

void
bcast(double&  value,
      MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcast(double)");
    }
}

void
bcast(bool&    value,
      int32_t  rank,
      MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        int32_t mvalue = 0;
            
        if (rank == mpi::master()) mvalue = (value) ? 1 : 0;
            
        mpi::bcast(mvalue, comm);
            
        value = (mvalue == 1) ? true : false;
    }
}

void
bcast(char&    value,
      MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, MPI_CHAR, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::bcast(char)");
    }
}
    
void
bcast(std::vector<int32_t>& vector,
      int32_t               rank,
      MPI_Comm              comm)
{
    if (ENABLE_MPI)
    {
        int32_t veclen = 0;

        if (rank == mpi::master()) veclen = static_cast<int32_t>(vector.size());

        mpi::bcast(veclen, comm);
        
        if (rank != mpi::master()) vector.clear();

        for (int32_t i = 0; i < veclen; i++)
        {
            int32_t mvalue = 0;

            if (rank == mpi::master()) mvalue = vector[i];

            mpi::bcast(mvalue, comm);

            if (rank != mpi::master()) vector.push_back(mvalue);
        
            MPI_Barrier(comm);
        }
    }
}
    
void
bcast(std::vector<double>& vector,
      int32_t              rank,
      MPI_Comm             comm)
{
    if (ENABLE_MPI)
    {
        int32_t veclen = 0;

        if (rank == mpi::master()) veclen = static_cast<int32_t>(vector.size());

        mpi::bcast(veclen, comm);
        
        if (rank != mpi::master()) vector.clear();

        for (int32_t i = 0; i < veclen; i++)
        {
            double mvalue = 0;

            if (rank == mpi::master()) mvalue = vector[i];

            mpi::bcast(mvalue, comm);

            if (rank != mpi::master()) vector.push_back(mvalue);
        
            MPI_Barrier(comm);
        }
    }
}
 
void
bcast(std::string& str,
      int32_t      rank,
      MPI_Comm     comm)
{
    if (ENABLE_MPI)
    {
        int32_t strwidth = 0;

        if (rank == mpi::master()) strwidth = static_cast<int32_t>(str.size());

        mpi::bcast(strwidth, comm);

        if (rank != mpi::master()) str.clear();
        
        for (int32_t i = 0; i < strwidth; i++)
        {
            char symbol;

            if (rank == mpi::master()) symbol = str[i];

            mpi::bcast(symbol, comm);

            if (rank != mpi::master()) str.append(1, symbol);
            
            MPI_Barrier(comm);
        }
    }
}
    
void
bcast(std::vector<std::string>& vector,
      int32_t                   rank,
      MPI_Comm                  comm)
{
    if (ENABLE_MPI)
    {
        int32_t veclen = 0;

        if (rank == mpi::master()) veclen = static_cast<int32_t>(vector.size());

        mpi::bcast(veclen, comm);
        
        if (rank != mpi::master()) vector.clear(); 

        for (int32_t i = 0; i < veclen; i++)
        {
            std::string mstr;

            if (rank == mpi::master()) mstr = vector[i];

            mpi::bcast(mstr, rank, comm);

            if (rank != mpi::master()) vector.push_back(mstr);
        
            MPI_Barrier(comm);
        }
    }
}

void
send(      double&  value,
     const int32_t  rank,
           MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Send(&value, 1, MPI_DOUBLE, rank, 0, comm);
            
        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::send(double)");
    }
}
  
void
receive(      double&  value,
        const int32_t  rank,
              MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        MPI_Status mstat;
        
        auto merror = MPI_Recv(&value, 1, MPI_DOUBLE, rank, 0, comm, &mstat);
        
        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::receive(double)");
    }
}
    
int32_t
batch_size(const int32_t nElements,
           const int32_t rank,
           const int32_t nodes)
{
    int32_t numelem = nElements / nodes;

    int32_t nremind = nElements % nodes;

    if ((nremind != 0) && (rank < nremind)) numelem++;

    return numelem;
}

int32_t
batch_offset(const int32_t nElements,
             const int32_t rank,
             const int32_t nodes)
{
    int32_t index = 0;

    for (int32_t i = 0; i < rank; i++)
    {
        index += mpi::batch_size(nElements, i, nodes);
    }

    return index;
}

void
batches_pattern(      int32_t* pattern,
                const int32_t  nElements,
                const int32_t  nodes)
{
    for (int32_t i = 0; i < nodes; i++)
    {
        pattern[i] = mpi::batch_size(nElements, i, nodes);
    }
}

void
gather(int32_t* vector,
       int32_t  value,
       int32_t  rank,
       MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Gather(&value, 1, MPI_INT32_T, vector, 1, MPI_INT32_T,
                                 mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::gather(integer)");
    }
    else
    {
        vector[0] = value; 
    }
}
    
void
gather(double*  vector,
       double   value,
       int32_t  rank,
       MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        auto merror = MPI_Gather(&value, 1, MPI_DOUBLE, vector, 1, MPI_DOUBLE,
                                 mpi::master(), comm);
        
        if (merror != MPI_SUCCESS) mpi::abort(merror, "mpi::gather(double)");
    }
    else
    {
        vector[0] = value;
    }
}
    
// TODO: Add other MPI functions for generic types
    
void
abort(const int   errorcode,
      const char* label)
{
    if (ENABLE_MPI)
    {
        int32_t errclass = 0;

        MPI_Error_class(errorcode, &errclass);
        
        int32_t errlen =0;
        
        char errstr[MPI_MAX_ERROR_STRING];
        
        MPI_Error_string(errorcode, errstr, &errlen);

        std::stringstream sst;

        sst << "**** Critical Error in " << label << " ****" << std::endl;

        sst << "MPI ERROR " << errclass << ": " << errstr << std::endl;

        std::cerr << sst.str();

        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}
    
} // mpi namespace

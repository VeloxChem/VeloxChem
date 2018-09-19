//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MpiFunc_hpp
#define MpiFunc_hpp

#include <cstdint>
#include <vector>
#include <string>

#include "mpi.h"

namespace mpi { // mpi namespace
    
/**
 Defines a default rank of master MPI process.

 @return the rank of master MPI process.
*/
inline int32_t master() { return 0;}
  
/**
 Initializes parallel execution mode driven by MPI.

 @param argc the number of command line arguments.
 @param argv the array of command line arguments.
 @return true if success, false otherwise.
 */
bool init(int    argc,
          char** argv);
    
/**
 Exits parallel execution mode driven by MPI.

 @return true if success, false otherwise.
 */
bool finalize();
    
/**
 Determines a rank of MPI process within MPI communicator.

 @param comm the MPI communicator.
 @return the rank of MPI process.
*/
int32_t rank(MPI_Comm comm);
    
/**
 Determines a number of MPI processes within MPI communicator.

 @param comm the MPI communicator.
 @return the number of MPI processes.
 */
int32_t nodes(MPI_Comm comm);

/**
 Compares two MPI communicators.

 @param comm1 the first MPI communicator.
 @param comm2 the second MPI communicator.
 @return true if MPI communicators is equal, false otherwise.
 */
bool compare(MPI_Comm comm1,
             MPI_Comm comm2);
    
/**
 Broadcasts an integer number within MPI communicator.

 @param value the integer number.
 @param comm the MPI communicator.
 */
void bcast(int32_t& value,
           MPI_Comm comm);

/**
 Broadcasts a real number within MPI communicator.

 @param value the real number.
 @param comm the MPI communicator.
 */
void bcast(double&  value,
           MPI_Comm comm);
    
/**
 Broadcasts a boolean within MPI communicator.

 @param value the boolean value.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(bool&    value,
           int32_t  rank,
           MPI_Comm comm);
    
/**
 Broadcasts a symbol within domain of MPI communicator.

 @param value the symbol.
 @param comm the MPI communicator.
 */
void bcast(char&    value,
           MPI_Comm comm);
    
/**
 Broadcasts vector of integer numbers within domain of MPI communicator.

 @param vector the vector of integer numbers.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::vector<int32_t>& vector,
           int32_t               rank,
           MPI_Comm              comm);
    
/**
 Broadcasts vector of real numbers within domain of MPI communicator.

 @param vector the vector of real numbers.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::vector<double>& vector,
           int32_t              rank,
           MPI_Comm             comm);
  
/**
 Broadcasts a string within domain of MPI communicator.

 @param str the string.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::string& str,
           int32_t      rank,
           MPI_Comm     comm);
    
    
/**
 Broadcasts vector of stringd within domain of MPI communicator.

 @param vector the vector of strings.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::vector<std::string>& vector,
           int32_t                   rank,
           MPI_Comm                  comm);
    
/**
 Sends a real number to destination MPI process.
 
 @param value the real number.
 @param rank the rank of destination MPI process.
 @param comm the MPI communicator.
*/
void
send(      double&  value,
     const int32_t  rank,
           MPI_Comm comm);
    
/**
 Receives a real number to source MPI process.
     
 @param value the real number.
 @param rank the rank of source MPI process.
 @param comm the MPI communicator.
*/
void
receive(      double&  value,
        const int32_t  rank,
              MPI_Comm comm);
   
/**
 Determines batch size associated with MPI process for data vector within domain
 of MPI communicator.

 @param nElements the size of data vector.
 @param rank the rank of MPI process.
 @param nodes the number of nodes in MPI communicator domain.
 @return the size of data batch.
 */
int32_t batch_size(const int32_t nElements,
                   const int32_t rank,
                   const int32_t nodes);
    
/**
 Determines offset of batch associated with MPI process within indexing
 space of MPI communicator domain.

 @param nElements the number of elements in data vector.
 @param rank the rank of MPI process.
 @param nodes the number of nodes in MPI domain.
 @return the offset of batch.
 */
int32_t batch_offset(const int32_t nElements,
                     const int32_t rank,
                     const int32_t nodes);

/**
 Creates batches distribution pattern for data vector with given number of
 elements.

 @param pattern the batches distribution pattern.
 @param nElements the number of elements in data vector.
 @param nodes the number of nodes in MPI communicator domain.
 */
void batches_pattern(      int32_t* pattern,
                     const int32_t  nElements,
                     const int32_t  nodes);
/**
 Gathers vector of integers on master MPI process by taking single integer from
 all MPI processes within domain of MPI communicator.

 @param vector the vector of integers.
 @param value the integer value.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void gather(int32_t* vector,
            int32_t  value,
            int32_t  rank,
            MPI_Comm comm);
    
/**
 Gathers vector of real numbers on master MPI process by taking single real
 number from all MPI processes within domain of MPI communicator.

 @param vector the vector of real numbers.
 @param value the real number.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void gather(double*  vector,
            double   value,
            int32_t  rank,
            MPI_Comm comm);

/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
void abort(const int   errorcode,
           const char* label);
    
} // mpi namespace


#endif /* MpiFunc_hpp */

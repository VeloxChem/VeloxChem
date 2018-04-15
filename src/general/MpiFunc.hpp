//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef MpiFunc_hpp
#define MpiFunc_hpp

#include <cstdint>
#include <vector>

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
bool init(int argc, char** argv);
    
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
bool compare(MPI_Comm comm1, MPI_Comm comm2);
    
/**
 Broadcasts an integer number within MPI communicator.

 @param value the integer number.
 @param comm the MPI communicator.
 */
void bcast(int32_t& value, MPI_Comm comm);

/**
 Broadcasts a real number within MPI communicator.

 @param value the real number.
 @param comm the MPI communicator.
 */
void bcast(double& value, MPI_Comm comm);
    
/**
 Broadcasts a boolean within MPI communicator.

 @param value the boolean value.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(bool& value, int32_t rank, MPI_Comm comm);
    
/**
 Broadcasts vector of integer numbers within domain of MPI communicator.

 @param vector the vector of integer numbers.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::vector<int32_t>& vector, int32_t rank, MPI_Comm comm);
    
/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
void abort(const int errorcode, const char* label);
    
} // mpi namespace


#endif /* MpiFunc_hpp */

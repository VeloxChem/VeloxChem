//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef MpiFunc_hpp
#define MpiFunc_hpp

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "DenseMatrix.hpp"

namespace mpi {

/**
 Gets default rank of master MPI process.

 @return the rank of master MPI process.
*/
inline constexpr auto
master() -> int64_t
{
    return 0;
}

/**
 Initializes parallel execution mode driven by MPI.

 @param argc the number of command line arguments.
 @param argv the array of command line arguments.
 @return true if success, false otherwise.
 */
auto init(int argc, char** argv) -> bool;

/**
 Check if MPI has been initialized.

 @return true if MPI has been initialized, false otherwise.
 */
auto initialized() -> bool;

/**
 Exits parallel execution mode driven by MPI.

 @return true if success, false otherwise.
 */
auto finalize() -> bool;

/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
auto abort(const int errorcode, const char* label) -> void;

/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
auto abort(const int errorcode, const std::string& label) -> void;

/**
 Determines a rank of MPI process within MPI communicator.

 @param comm the MPI communicator.
 @return the rank of MPI process.
*/
auto rank(MPI_Comm comm) -> int64_t;

/**
 Determines a number of MPI processes within MPI communicator.

 @param comm the MPI communicator.
 @return the number of MPI processes.
 */
auto nodes(MPI_Comm comm) -> int64_t;

/**
 Broadcasts scalar.

 @param val the scalar.
 @param comm the MPI communicator.
 @return the broadcast scalar.
 */
auto bcastScalar(const int64_t val, MPI_Comm comm) -> int64_t;

/**
 Broadcasts scalar.

 @param val the scalar.
 @param comm the MPI communicator.
 @return the broadcast scalar.
 */
auto bcastScalar(const double val, MPI_Comm comm) -> double;

auto bcastStdVectorInt(const std::vector<int64_t>& vector, MPI_Comm comm) -> std::vector<int64_t>;

/**
 Broadcasts dense matrix.

 @param matrix the dense matrix.
 @param comm the MPI communicator.
 @return the broadcast dense matrix.
 */
auto bcastDenseMatrix(const CDenseMatrix& matrix, MPI_Comm comm) -> CDenseMatrix;

/**
 Scatters std vector.

 @param vec the std vector.
 @param comm the MPI communicator.
 @return the scattered vector.
 */
auto scatterStdVector(const std::vector<int64_t>& vec, MPI_Comm comm) -> std::vector<int64_t>;

/**
 Gathers a dense matrix by columns.

 @param matrix the dense matrix.
 @param comm the MPI communicator.
 @return the gathered dense matrix.
 */
auto gatherDenseMatricesByColumns(const CDenseMatrix& matrix, MPI_Comm comm) -> CDenseMatrix;

}  // namespace mpi

#endif /* MpiFunc_hpp */

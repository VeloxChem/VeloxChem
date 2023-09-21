//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef MpiFunc_hpp
#define MpiFunc_hpp

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

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

}  // namespace mpi

#endif /* MpiFunc_hpp */

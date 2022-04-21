//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include "MetaUtils.hpp"

namespace mpi {
template <typename T>
inline MPI_Datatype type_v = metautils::type_to_mpi_datatype<T>::value;

template <typename T>
inline std::string
error_message(std::string op)
{
    return op + "(" + metautils::type_to_string<T>::name + ")";
}

/**
 Gets default rank of master MPI process.

 @return the rank of master MPI process.
*/
inline constexpr int32_t
master()
{
    return 0;
}

/**
 Initializes parallel execution mode driven by MPI.

 @param argc the number of command line arguments.
 @param argv the array of command line arguments.
 @return true if success, false otherwise.
 */
bool init(int argc, char** argv);

/**
 Check if MPI has been initialized.

 @return true if MPI has been initialized, false otherwise.
 */
bool initialized();

/**
 Exits parallel execution mode driven by MPI.

 @return true if success, false otherwise.
 */
bool finalize();

/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
void abort(const int errorcode, const char* label);

/**
 Terminates all MPI processes and prints error message to standard error stream.

 @param errorcode the MPI error code.
 @param label the label of function in which MPI error occured.
 */
void abort(const int errorcode, const std::string& label);

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
 Duplicates an existing MPI communicator.

 @param comm1 the MPI communicator.
 @param comm2 the copy of MPI communicator.
 */
void duplicate(MPI_Comm comm1, MPI_Comm* comm2);

/**
 Destroys data in MPI communicator.

 @param comm the pointer to MPI communicator.
 */
void destroy(MPI_Comm* comm);

/**
 Compares two MPI communicators.

 @param comm1 the first MPI communicator.
 @param comm2 the second MPI communicator.
 @return true if MPI communicators is equal, false otherwise.
 */
bool compare(MPI_Comm comm1, MPI_Comm comm2);

/**
 * Broadcasts a scalar value within MPI communicator.
 *
 * @tparam T type of value.
 * @param value the value.
 * @param comm the MPI communicator.
 */
template <typename T>
auto
bcast(T& value, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Bcast(&value, 1, mpi::type_v<T>, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, error_message<T>("mpi::bcast"));
    }
}

/**
 Broadcasts a string within domain of MPI communicator.

 @param str the string.
 @param rank the rank of MPI process.
 @param comm the MPI communicator.
 */
void bcast(std::string& str, int32_t rank, MPI_Comm comm);

/**
 * Broadcasts vector within domain of MPI communicator.
 *
 * @tparam T type of underlying scalar.
 * @param vector the vector.
 * @param comm the MPI communicator.
 */
template <typename T>
void
bcast(std::vector<T>& vector, int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        auto veclen = (rank == mpi::master()) ? static_cast<int32_t>(vector.size()) : int32_t{0};

        mpi::bcast(veclen, comm);

        if (rank != mpi::master()) vector.clear();

        // a range-based for loop makes this broadcast hang!
        for (int32_t i = 0; i < veclen; ++i)
        {
            auto mvalue = (rank == mpi::master()) ? vector[i] : T{0};

            if constexpr (std::is_same_v<T, std::string>)
            {
                mpi::bcast(mvalue, rank, comm);
            }
            else
            {
                mpi::bcast(mvalue, comm);
            }

            if (rank != mpi::master()) vector.push_back(mvalue);

            MPI_Barrier(comm);
        }
    }
}

/**
 * Sends a scalar value to destination MPI process.
 *
 * @tparam T type of scalar value.
 * @param value the scalar.
 * @param rank the rank of destination MPI process.
 * @param comm the MPI communicator.
 */
template <typename T>
auto
send(T& value, const int32_t rank, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Send(&value, 1, mpi::type_v<T>, rank, 0, comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, error_message<T>("mpi::send"));
    }
}

/**
 * Receives a scalar value from source MPI process.
 *
 * @tparam T type of scalar value.
 * @param value the real number.
 * @param rank the rank of source MPI process.
 * @param comm the MPI communicator.
 */
template <typename T>
auto
receive(T& value, const int32_t rank, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    if constexpr (ENABLE_MPI)
    {
        MPI_Status mstat;

        auto merror = MPI_Recv(&value, 1, mpi::type_v<T>, rank, 0, comm, &mstat);

        if (merror != MPI_SUCCESS) mpi::abort(merror, error_message<T>("mpi::receive"));
    }
}

/**
 Determines batch size associated with MPI process for data vector within domain
 of MPI communicator.

 @param nElements the size of data vector.
 @param rank the rank of MPI process.
 @param nodes the number of nodes in MPI communicator domain.
 @return the size of data batch.
 */
int32_t batch_size(const int32_t nElements, const int32_t rank, const int32_t nodes);

/**
 Determines offset of batch associated with MPI process within indexing
 space of MPI communicator domain.

 @param nElements the number of elements in data vector.
 @param rank the rank of MPI process.
 @param nodes the number of nodes in MPI domain.
 @return the offset of batch.
 */
int32_t batch_offset(const int32_t nElements, const int32_t rank, const int32_t nodes);

/**
 Creates batches distribution pattern for data vector with given number of
 elements.

 @param pattern the batches distribution pattern.
 @param nElements the number of elements in data vector.
 @param nodes the number of nodes in MPI communicator domain.
 */
void batches_pattern(int32_t* pattern, const int32_t nElements, const int32_t nodes);

/**
 * Gathers vector of arithmetic scalars on master MPI process by taking single scalar from
 * all MPI processes within domain of MPI communicator.
 *
 * @tparam T type of scalar value.
 * @param vector the vector of integers.
 * @param value the integer value.
 * @param comm the MPI communicator.
 */
template <typename T>
auto
gather(T* vector, T value, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Gather(&value, 1, mpi::type_v<T>, vector, 1, mpi::type_v<T>, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, error_message<T>("mpi::gather"));
    }
    else
    {
        vector[0] = value;
    }
}

/**
 Sums vector of real numbers on master MPI process from vectors of real numbers
 accross MPI processes.

 @param source the vector of real numbers associated with MPI process.
 @param destination the summed vector of real numbers.
 @param nElements the number of elements in vector.
 @param comm the MPI communicator.
 */
template <typename T>
auto
reduce_sum(const T* source, T* destination, const int32_t nElements, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), void())
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Reduce(source, destination, nElements, mpi::type_v<T>, MPI_SUM, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, error_message<T>("mpi::reduce_sum"));
    }
    else
    {
        std::copy_n(source, nElements, destination);
    }
}

/**
 Performs reduction operation over given real numbers accross MPI processes.

 @param value the partial value of real number.
 @param comm the MPI communicator
 @return the sum of partial real numbes on master node, 0.0 on other nodes.
 */
template <typename T>
auto
reduce_sum(const T value, MPI_Comm comm) -> decltype((void)(std::is_arithmetic_v<T>), T())
{
    if constexpr (ENABLE_MPI)
    {
        auto dval = T{0.0};

        mpi::reduce_sum(&value, &dval, 1, comm);

        MPI_Barrier(comm);

        return dval;
    }

    return value;
}
}  // namespace mpi

#endif /* MpiFunc_hpp */

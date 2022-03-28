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

#include "MpiFunc.hpp"

#include <mpi.h>

#include <iostream>
#include <sstream>

namespace mpi {  // mpi namespace
bool
init(int argc, char** argv)
{
    if constexpr (ENABLE_MPI)
    {
        if (initialized())
        {
            return true;
        }

        int32_t mlevel = 0;

        auto merror = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mlevel);

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
    if constexpr (ENABLE_MPI)
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
    if constexpr (ENABLE_MPI)
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

void
abort(const int errorcode, const char* label)
{
    if constexpr (ENABLE_MPI)
    {
        int32_t errclass = 0;

        MPI_Error_class(errorcode, &errclass);

        int32_t errlen = 0;

        char errstr[MPI_MAX_ERROR_STRING];

        MPI_Error_string(errorcode, errstr, &errlen);

        std::stringstream sst;

        sst << "**** Critical Error in " << label << " ****" << std::endl;

        sst << "MPI ERROR " << errclass << ": " << errstr << std::endl;

        std::cerr << sst.str();

        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}

void
abort(const int errorcode, const std::string& label)
{
    mpi::abort(errorcode, label.c_str());
}

int32_t
rank(MPI_Comm comm)
{
    int32_t mrank = 0;

    if constexpr (ENABLE_MPI) MPI_Comm_rank(comm, &mrank);

    return mrank;
}

int32_t
nodes(MPI_Comm comm)
{
    int32_t mnodes = 1;

    if constexpr (ENABLE_MPI) MPI_Comm_size(comm, &mnodes);

    return mnodes;
}

void
duplicate(MPI_Comm comm1, MPI_Comm* comm2)
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Comm_dup(comm1, comm2);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::duplicate");
        }
    }
}

void
destroy(MPI_Comm* comm)
{
    if constexpr (ENABLE_MPI)
    {
        auto merror = MPI_Comm_free(comm);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::destroy");
        }
    }
}

bool
compare(MPI_Comm comm1, MPI_Comm comm2)
{
    if constexpr (ENABLE_MPI)
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
bcast(std::string& str, int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        auto strwidth = (rank == mpi::master()) ? static_cast<int32_t>(str.size()) : int32_t{0};

        mpi::bcast(strwidth, comm);

        if (rank != mpi::master()) str.clear();

        // a range-based for loop makes this broadcast hang!
        for (int32_t i = 0; i < strwidth; ++i)
        {
            auto symbol = (rank == mpi::master()) ? str[i] : char{};

            mpi::bcast(symbol, comm);

            if (rank != mpi::master()) str.append(1, symbol);

            MPI_Barrier(comm);
        }
    }
}

int32_t
batch_size(const int32_t nElements, const int32_t rank, const int32_t nodes)
{
    int32_t numelem = nElements / nodes;

    int32_t nremind = nElements % nodes;

    if ((nremind != 0) && (rank < nremind)) numelem++;

    return numelem;
}

int32_t
batch_offset(const int32_t nElements, const int32_t rank, const int32_t nodes)
{
    int32_t index = 0;

    for (int32_t i = 0; i < rank; i++)
    {
        index += mpi::batch_size(nElements, i, nodes);
    }

    return index;
}

void
batches_pattern(int32_t* pattern, const int32_t nElements, const int32_t nodes)
{
    for (int32_t i = 0; i < nodes; i++)
    {
        pattern[i] = mpi::batch_size(nElements, i, nodes);
    }
}

// TODO: Add other MPI functions for generic types
}  // namespace mpi

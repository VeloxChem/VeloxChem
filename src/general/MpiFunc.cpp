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

#include "MpiFunc.hpp"

#include <mpi.h>

#include <iostream>
#include <sstream>

namespace mpi {  // mpi namespace

auto
init(int argc, char** argv) -> bool
{
    if constexpr (ENABLE_MPI)
    {
        if (initialized())
        {
            return true;
        }

        int32_t mlevel_int32 = 0;

        auto merror = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mlevel_int32);

        if (merror != MPI_SUCCESS)
        {
            mpi::abort(merror, "mpi::init()");

            return false;
        }
    }

    return true;
}

auto
initialized() -> bool
{
    if constexpr (ENABLE_MPI)
    {
        int32_t minit_int32 = 0;

        MPI_Initialized(&minit_int32);

        if (minit_int32 == 1) return true;
    }

    return false;
}

auto
finalize() -> bool
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

auto
abort(const int errorcode, const char* label) -> void
{
    if constexpr (ENABLE_MPI)
    {
        int32_t errclass_int32 = 0;

        MPI_Error_class(errorcode, &errclass_int32);

        int32_t errlen_int32 = 0;

        char errstr[MPI_MAX_ERROR_STRING];

        MPI_Error_string(errorcode, errstr, &errlen_int32);

        std::stringstream sst;

        sst << "**** Critical Error in " << label << " ****" << std::endl;

        sst << "MPI ERROR " << errclass_int32 << ": " << errstr << std::endl;

        std::cerr << sst.str();

        MPI_Abort(MPI_COMM_WORLD, errorcode);
    }
}

auto
abort(const int errorcode, const std::string& label) -> void
{
    mpi::abort(errorcode, label.c_str());
}

auto
rank(MPI_Comm comm) -> int64_t
{
    int32_t mrank_int32 = 0;

    if constexpr (ENABLE_MPI) MPI_Comm_rank(comm, &mrank_int32);

    return static_cast<int64_t>(mrank_int32);
}

auto
nodes(MPI_Comm comm) -> int64_t
{
    int32_t mnodes_int32 = 1;

    if constexpr (ENABLE_MPI) MPI_Comm_size(comm, &mnodes_int32);

    return static_cast<int64_t>(mnodes_int32);
}

}  // namespace mpi

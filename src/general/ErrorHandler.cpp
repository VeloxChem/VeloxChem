//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "ErrorHandler.hpp"

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "MpiFunc.hpp"

namespace errors {  // errors namespace

auto
msg(const std::string& message, const std::string& header) -> void
{
    std::stringstream sst;

    sst << "**** " << header << " ****" << std::endl;

    sst << "     " << message << std::endl << std::endl;

    std::cerr << sst.str();
}

auto
assertMsgCritical(const bool condition, const std::string& message) -> void
{
    if (!condition)
    {
        msg(message, "Critical Error");

        if (mpi::initialized() && (mpi::nodes(MPI_COMM_WORLD) > 1))
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }

        std::abort();
    }
}

}  // namespace errors

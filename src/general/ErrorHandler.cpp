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

#include "ErrorHandler.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>

#include "MpiFunc.hpp"

namespace errors {  // errors namespace

void
assertMsgCritical(const bool condition, const std::string& label)
{
    if (!condition)
    {
	    msgCritical(label);
    }
}

void
msgCritical(const std::string& label)
{
    std::stringstream sst;

    sst << "**** Critical Error";

    if (mpi::initialized() && mpi::nodes(MPI_COMM_WORLD) > 1)
    {
        sst << " (process " << mpi::rank(MPI_COMM_WORLD) << ")";
    }

    sst << " ****" << std::endl;

    sst << "     " <<  label << std::endl << std::endl;

    std::cerr << sst.str();

    if (mpi::initialized() && mpi::nodes(MPI_COMM_WORLD) > 1)
    {
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
    }

    std::abort();
}

}  // namespace errors

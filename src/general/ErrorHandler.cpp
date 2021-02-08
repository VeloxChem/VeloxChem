//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
}

}  // namespace errors

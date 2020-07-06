//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace errors {  // errors namespace

void
assertMsgCritical(const bool condition, const std::string& label)
{
    if (!condition)
    {
        bool masternode = false;

        if (!mpi::initialized() || mpi::rank(MPI_COMM_WORLD) == mpi::master())
        {
            masternode = true;
        }

        if (masternode)
        {
            std::stringstream sst;

            sst << "**** Critical Error in " << label << " ****" << std::endl;

            std::cerr << sst.str();
        }

        if (mpi::initialized() && mpi::nodes(MPI_COMM_WORLD) > 1)
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        }

        std::abort();
    }
}

}  // namespace errors

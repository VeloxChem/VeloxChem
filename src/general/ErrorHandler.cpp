//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi::initialized() && mpi::nodes(MPI_COMM_WORLD) > 1)
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
        }

        std::abort();
    }
}

}  // namespace errors

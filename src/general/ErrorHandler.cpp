//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"

#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>

namespace errors { // errors namespace

void
assertMsgCritical(const bool         condition,
                  const std::string& label)
{
    if (! condition)
    {
        if (mpi::initialized())
        {
            mpi::abort(MPI_ERR_OTHER, label.c_str());
        }

        std::stringstream sst;

        sst << "**** Critical Error in " << label << " ****" << std::endl;

        std::cerr << sst.str();

        std::abort();
    }
}

} // errors namespace

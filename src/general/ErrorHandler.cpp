//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ErrorHandler.hpp"

#include <string>
#include <cstdlib>
#include <iostream>

namespace errors { // errors namespace

void
assertMsgCritical(const bool         condition,
                  const std::string& message)
{
    if (! condition)
    {
        std::cerr << std::endl;
        
        std::cerr << "Critical Error: " << message << std::endl << std::endl;

        std::abort();
    }
}

} // errors namespace

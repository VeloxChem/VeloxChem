//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DeviceProp_hpp
#define DeviceProp_hpp

#include "OutputStream.hpp"

namespace gpu {

    #ifdef ENABLE_GPU

    void get_device_prop(COutputStream& oStream);

    #endif

}

#endif

//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DeviceFunc.hpp"

#include "ErrorHandler.hpp"

namespace gpu {  // gpu namespace

void setDevices(const int32_t iDevice);
{
#ifdef ENABLE_GPU
    auto cerr = cudaSetDevice(iDevice);

    errors::assertMsgCritical(cerr == cudaSuccess, {"setCudaDevice"});
#endif
}

}  // namespace gpu

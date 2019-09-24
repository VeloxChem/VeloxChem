//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CudaGenFunc.hpp"

namespace gpu {  // gpu namespace

int32_t
getNumberOfGridBlocks(const int32_t nElements,
                      const int32_t gridBlockSize)
{
    auto nblocks = nElements / gridBlockSize;
    
    if ((nElements % gridBlockSize) != 0) nblocks++;
    
    return nblocks;
}
    
}  // namespace gpu

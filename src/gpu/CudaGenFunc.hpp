//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CudaGenFunc_hpp
#define CudaGenFunc_hpp

#include <cstdint>

namespace gpu {

/**
 Gets number of grid blocks used to launch CUDA kernel for given
 size of linear data array.

 @param nElements the number of elements in linear data array.
 @param gridBlockSize the grid block size.
 @return the number of grid blocks. 
 */
int32_t getNumberOfGridBlocks(const int32_t nElements,
                               const int32_t gridBlockSize); 
    
}
#endif /* CudaGenFunc_hpp */

//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MemAlloc.hpp"

#include "mkl.h"

namespace mem { // mem namespace

void*
malloc(const size_t size)
{
    return MKL_malloc(size, VLX_ALIGN);
}

void
free(void* pointer)
{
    if (pointer != nullptr) MKL_free(pointer);
}

} // mem namespace

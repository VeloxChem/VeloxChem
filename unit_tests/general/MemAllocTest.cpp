//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MemAllocTest.hpp"

#include "MemAlloc.hpp"

TEST_F(CMemAllocTest, Malloc)
{
    double* ptr = (double*) mem::malloc(53 * sizeof(double));

    ASSERT_EQ(0, ((size_t) ptr) % VLX_ALIGN);

    mem::free(ptr);
}

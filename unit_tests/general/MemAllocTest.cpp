//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MemAllocTest.hpp"

#include "MemAlloc.hpp"

TEST_F(CMemAllocTest, Malloc)
{
    double* ptr = (double*) mem::malloc(53 * sizeof(double));

    ASSERT_EQ(0u, ((size_t) ptr) % VLX_ALIGN);

    mem::free(ptr);
}

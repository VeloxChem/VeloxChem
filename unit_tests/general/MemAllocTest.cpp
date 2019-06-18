//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MemAllocTest.hpp"

#include "MemAlloc.hpp"

TEST_F(CMemAllocTest, Malloc)
{
    double* ptr = (double*)mem::malloc(53 * sizeof(double));

    ASSERT_EQ(0u, ((size_t)ptr) % VLX_ALIGN);

    mem::free(ptr);
}

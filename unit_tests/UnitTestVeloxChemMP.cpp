//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "gtest/gtest.h"
#include "MpiFunc.hpp"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    mpi::init(argc, argv);

    return RUN_ALL_TESTS();
}

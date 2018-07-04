//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GenFuncTest.hpp"

#include "MemBlock2D.hpp"
#include "GenFunc.hpp"

TEST_F(CGenFuncTest, Contract)
{
    int32_t sposvec[3] = {0, 2, 3};
    
    int32_t eposvec[3] = {2, 3, 6};
    
    CMemBlock2D<double> pdat(CMemBlock2D<double>({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                                  2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                                  2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                                 6, 3));
    
    CMemBlock2D<double> cdat(3, 2);
    
    genfunc::contract(cdat, 0, pdat, 1, sposvec, eposvec, 2, 3);
    
    CMemBlock2D<double> tdat({5.0, 6.0, 16.0, 8.1, -1.0, 17.0}, 3, 2);
    
    ASSERT_EQ(cdat, tdat);
}

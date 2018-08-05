//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef BoysFunctionTest_hpp
#define BoysFunctionTest_hpp

#include "gtest/gtest.h"

#include "DenseMatrix.hpp"

class CBoysFunctionTest : public ::testing::Test
{
protected:

    CBoysFunctionTest() {};

    virtual ~CBoysFunctionTest() {};

    CMemBlock<double> getSmallArguments() const;

    CMemBlock<double> getMediumArguments() const;

    CMemBlock<double> getLargeArguments() const;
};


#endif /* BoysFunctionTest_hpp */

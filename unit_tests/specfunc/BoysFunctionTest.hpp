//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef BoysFunctionTest_hpp
#define BoysFunctionTest_hpp

#include "gtest/gtest.h"

#include "MemBlock.hpp"

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

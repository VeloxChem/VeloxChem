//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MpiFuncTest_hpp
#define MpiFuncTest_hpp

#include "gtest/gtest.h"

class CMpiFuncTest : public ::testing::Test
{
   protected:
    CMpiFuncTest(){};

    virtual ~CMpiFuncTest(){};
};

#endif /* MpiFuncTest_hpp */

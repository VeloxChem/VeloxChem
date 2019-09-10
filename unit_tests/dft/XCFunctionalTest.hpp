//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef XCFunctionalTest_hpp
#define XCFunctionalTest_hpp

#include "gtest/gtest.h"

class CXCFunctionalTest : public ::testing::Test
{
protected:
    CXCFunctionalTest() {};
    
    virtual ~CXCFunctionalTest() {};
};

#endif /* XCFunctionalTest_hpp */

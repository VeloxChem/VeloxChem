//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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

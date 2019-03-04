//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ErrorHandlerTest.hpp"

#include "ErrorHandler.hpp"

TEST_F(CErrorHandlerTest, Critical)
{
    const bool condition = true;

    errors::assertMsgCritical(condition, "");

    ASSERT_TRUE(condition);
}

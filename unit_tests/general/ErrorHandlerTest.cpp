//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ErrorHandlerTest.hpp"

#include "ErrorHandler.hpp"

TEST_F(CErrorHandlerTest, Critical)
{
    const bool condition = true;

    errors::assertMsgCritical(condition, "");

    ASSERT_TRUE(condition);
}

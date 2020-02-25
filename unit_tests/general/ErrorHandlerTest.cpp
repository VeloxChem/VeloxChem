//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ErrorHandlerTest.hpp"

#include "ErrorHandler.hpp"

TEST_F(CErrorHandlerTest, Critical)
{
    const bool condition = true;

    errors::assertMsgCritical(condition, "");

    ASSERT_TRUE(condition);
}

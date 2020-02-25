//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "MathConstTest.hpp"

#include "MathConst.hpp"

TEST_F(CMathConstTest, GetPiValue)
{
    ASSERT_NEAR(3.14159265358979323846, mathconst::getPiValue(), 1.0e-13);
}

//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CodataTest.hpp"

#include "Codata.hpp"

TEST_F(CCodataTest, GetBohrValueInAngstroms)
{
    ASSERT_NEAR(0.52917721092, units::getBohrValueInAngstroms(), 1.0e-13);
}

TEST_F(CCodataTest, GetHatreeValueInElectronVolts)
{
    ASSERT_NEAR(27.21138505, units::getHatreeValueInElectronVolts(), 1.0e-13);
}

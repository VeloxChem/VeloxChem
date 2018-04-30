//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AngularMomentumTest.hpp"

#include "AngularMomentum.hpp"

TEST_F(CAngularMomentumTest, To_SphericalComponents)
{
    ASSERT_EQ(1, angmom::to_SphericalComponents(0));

    ASSERT_EQ(3, angmom::to_SphericalComponents(1));

    ASSERT_EQ(5, angmom::to_SphericalComponents(2));

    ASSERT_EQ(7, angmom::to_SphericalComponents(3));

    ASSERT_EQ(9, angmom::to_SphericalComponents(4));

    ASSERT_EQ(11, angmom::to_SphericalComponents(5));

    ASSERT_EQ(13, angmom::to_SphericalComponents(6));
}

TEST_F(CAngularMomentumTest, To_CartesianComponents)
{
    ASSERT_EQ(1, angmom::to_CartesianComponents(0));

    ASSERT_EQ(3, angmom::to_CartesianComponents(1));

    ASSERT_EQ(6, angmom::to_CartesianComponents(2));

    ASSERT_EQ(10, angmom::to_CartesianComponents(3));

    ASSERT_EQ(15, angmom::to_CartesianComponents(4));

    ASSERT_EQ(21, angmom::to_CartesianComponents(5));

    ASSERT_EQ(28, angmom::to_CartesianComponents(6));
}

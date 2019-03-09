//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CubicGridTest.hpp"

#include <vector>

#include "CubicGrid.hpp"

TEST_F(CCubicGridTest, Getters)
{
    CCubicGrid grid({0.1, 0.2, 0.3}, {1.0, 2.0, 3.0}, {1, 2, 3});

    ASSERT_EQ(0.1, grid.originX());

    ASSERT_EQ(0.2, grid.originY());

    ASSERT_EQ(0.3, grid.originZ());

    ASSERT_EQ(1.0, grid.stepSizeX());

    ASSERT_EQ(2.0, grid.stepSizeY());

    ASSERT_EQ(3.0, grid.stepSizeZ());

    ASSERT_EQ(1, grid.numPointsX());

    ASSERT_EQ(2, grid.numPointsY());

    ASSERT_EQ(3, grid.numPointsZ());
}

TEST_F(CCubicGridTest, CopyMoveConstructor)
{
    CCubicGrid grid({0.1, 0.2, 0.3}, {1.0, 2.0, 3.0}, {1, 2, 3});

    CCubicGrid grid2(grid);

    ASSERT_EQ(grid, grid2);

    CCubicGrid grid3(std::move(grid));

    ASSERT_EQ(grid2, grid3);
}

TEST_F(CCubicGridTest, CopyMoveAssignment)
{
    CCubicGrid grid({0.1, 0.2, 0.3}, {1.0, 2.0, 3.0}, {1, 2, 3});

    CCubicGrid grid2;

    grid2 = grid;

    ASSERT_EQ(grid, grid2);

    CCubicGrid grid3;

    grid3 = std::move(grid);

    ASSERT_EQ(grid2, grid3);
}

//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

#include "LebedevLaikovQuadratureTest.hpp"

#include <vector>

#include "LebedevLaikovQuadrature.hpp"
#include "CheckFunctions.hpp"

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith6Points)
{
    CLebedevLaikovQuadrature lquad(6);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(6, grid.size(0));
    
    ASSERT_EQ(6, grid.size(1));
    
    ASSERT_EQ(6, grid.size(2));
    
    ASSERT_EQ(6, grid.size(3));

    std::vector<double> coordx({1.0, -1.0, 0.0, 0.0, 0.0, 0.0});

    vlxtest::compare(coordx, grid.data(0));

    std::vector<double> coordy({0.0, 0.0, 1.0, -1.0, 0.0, 0.0});

    vlxtest::compare(coordy, grid.data(1));

    std::vector<double> coordz({0.0, 0.0, 0.0, 0.0, 1.0, -1.0});

    vlxtest::compare(coordz, grid.data(2));

    std::vector<double> weights({0.1666666666666667, 0.1666666666666667,
                                 0.1666666666666667, 0.1666666666666667,
                                 0.1666666666666667, 0.1666666666666667});

    vlxtest::compare(weights, grid.data(3));

    vlxtest::checkNorm(grid.data(3), 6);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith50Points)
{
    CLebedevLaikovQuadrature lquad(50);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(50, grid.size(0));
    
    ASSERT_EQ(50, grid.size(1));
    
    ASSERT_EQ(50, grid.size(2));
    
    ASSERT_EQ(50, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 50);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith110Points)
{
    CLebedevLaikovQuadrature lquad(110);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(110, grid.size(0));
    
    ASSERT_EQ(110, grid.size(1));
    
    ASSERT_EQ(110, grid.size(2));
    
    ASSERT_EQ(110, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 110);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith194Points)
{
    CLebedevLaikovQuadrature lquad(194);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(194, grid.size(0));
    
    ASSERT_EQ(194, grid.size(1));
    
    ASSERT_EQ(194, grid.size(2));
    
    ASSERT_EQ(194, grid.size(3));
    
    vlxtest::checkNorm(grid.data(3), 194);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith302Points)
{
    CLebedevLaikovQuadrature lquad(302);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(302, grid.size(0));
    
    ASSERT_EQ(302, grid.size(1));
    
    ASSERT_EQ(302, grid.size(2));
    
    ASSERT_EQ(302, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 302);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith434Points)
{
    CLebedevLaikovQuadrature lquad(434);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(434, grid.size(0));
    
    ASSERT_EQ(434, grid.size(1));
    
    ASSERT_EQ(434, grid.size(2));
    
    ASSERT_EQ(434, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 434);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith590Points)
{
    CLebedevLaikovQuadrature lquad(590);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(590, grid.size(0));
    
    ASSERT_EQ(590, grid.size(1));
    
    ASSERT_EQ(590, grid.size(2));
    
    ASSERT_EQ(590, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 590);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith770Points)
{
    CLebedevLaikovQuadrature lquad(770);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(770, grid.size(0));
    
    ASSERT_EQ(770, grid.size(1));
    
    ASSERT_EQ(770, grid.size(2));
    
    ASSERT_EQ(770, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 770);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith974Points)
{
    CLebedevLaikovQuadrature lquad(974);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(974, grid.size(0));
    
    ASSERT_EQ(974, grid.size(1));
    
    ASSERT_EQ(974, grid.size(2));
    
    ASSERT_EQ(974, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 974);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith2030Points)
{
    CLebedevLaikovQuadrature lquad(2030);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(4, grid.blocks());

    ASSERT_EQ(2030, grid.size(0));
    
    ASSERT_EQ(2030, grid.size(1));
    
    ASSERT_EQ(2030, grid.size(2));
    
    ASSERT_EQ(2030, grid.size(3));

    vlxtest::checkNorm(grid.data(3), 2030);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith100Points)
{
    CLebedevLaikovQuadrature lquad(100);

    auto grid = lquad.getGridPoints();

    ASSERT_EQ(0, grid.blocks());

    ASSERT_EQ(0, grid.size(0));
    
    ASSERT_EQ(0, grid.size(1));
    
    ASSERT_EQ(0, grid.size(2));
    
    ASSERT_EQ(0, grid.size(3));
}





//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DenseLinearAlgebraTest.hpp"

#include "DenseLinearAlgebra.hpp"

TEST_F(CDenseLinearAlgebraTest, MultAB)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                       4, 3);
    
    CDenseMatrix matb({ 1.0, 2.0, 5.0,  2.0, 4.0,
                       -3.0, 3.0, 1.0, -1.0, 2.0,
                        3.0, 0.5, 7.8,  1.0, 2.1},
                      3, 5);
    
    auto matab = denblas::multAB(mata, matb);
    
    CDenseMatrix refab({ 5.0, 15.0, 44.2,  5.0, 22.4,
                        -9.0,  3.5, -4.2, -8.0, -3.9,
                        20.1, 22.4, 86.9, 16.7, 43.3,
                         7.0, 10.0, 38.2,  4.0, 16.4},
                       4, 5);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultAtB)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                       4, 3);
    
    CDenseMatrix matb({ 1.0, 2.0,
                        5.0,  2.0,
                       -3.0, 3.0,
                        1.0, -1.0},
                      4, 2);
    
    auto matab = denblas::multAtB(mata, matb);
    
    CDenseMatrix refab({ -30.0, 15.0,
                          13.1, 16.9,
                          -8.0, 27.0},
                       3, 2);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultABt)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                      4, 3);
    
    CDenseMatrix matb({ 1.0,  2.0, 5.0,
                        2.0, -3.0, 3.0},
                      2, 3);
    
    auto matab = denblas::multABt(mata, matb);
    
    CDenseMatrix refab({ 28.0,   7.0,
                          8.0, -12.0,
                         45.6,  26.1,
                         25.0,   8.0},
                       4, 2);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultDiagByA)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                       6.0, 2.3, 7.0,
                       1.0, 2.0, 4.0},
                      4, 3);
    
    CMemBlock<double> diagb({1.0, 2.0, -2.0, 5.0});
    
    auto matab = denblas::multDiagByA(diagb, mata);
    
    CDenseMatrix refab({  2.0,  3.0,   4.0,
                         -6.0,  6.0,   2.0,
                        -12.0, -4.6, -14.0,
                          5.0, 10.0,  20.0},
                       4, 3);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultDiagByAt)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                       4, 3);
    
    CMemBlock<double> diagb({1.0, 2.0, -2.0});
    
    auto matab = denblas::multDiagByAt(diagb, mata);
    
    CDenseMatrix refab({ 2.0, -3.0,  6.0,   1.0,
                         6.0,  6.0,  4.6,   4.0,
                        -8.0, -2.0, -14.0, -8.0},
                       3, 4);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, SubAB)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                      4, 3);
    
    CDenseMatrix matb({ 0.8, 1.2, 3.2,
                        2.1, 2.5, 1.7,
                        4.9, 2.3, 8.1,
                        0.2, 1.3, 5.1},
                      4, 3);
    
    auto matab = denblas::subAB(mata, matb);
    
    CDenseMatrix refab({ 1.2, 1.8,  0.8,
                        -5.1, 0.5, -0.7,
                         1.1, 0.0, -1.1,
                         0.8, 0.7, -1.1},
                       4, 3);
    
    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, AddAB)
{
    CDenseMatrix mata({ 2.0, 3.0, 4.0,
                       -3.0, 3.0, 1.0,
                        6.0, 2.3, 7.0,
                        1.0, 2.0, 4.0},
                      4, 3);
    
    CDenseMatrix matb({ 0.8, 1.2, 3.2,
                        2.1, 2.5, 1.7,
                        4.9, 2.3, 8.1,
                        0.2, 1.3, 5.1},
                      4, 3);
    
    auto matab = denblas::addAB(mata, matb, 2.0);
    
    CDenseMatrix refab({ 3.6, 5.4, 10.4,
                         1.2, 8.0,  4.4,
                        15.8, 6.9, 23.2,
                         1.4, 4.6, 14.2},
                       4, 3);
    
    ASSERT_EQ(matab, refab);
}

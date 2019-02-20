//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DenseMatrixTest.hpp"

#include "DenseMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CDenseMatrixTest, DefaultConstructor)
{
    CDenseMatrix ma;
    
    CDenseMatrix mb(std::vector<double>{}, 0, 0);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, ConstructorWithDimensions)
{
    CDenseMatrix ma(0, 0);
    
    CDenseMatrix mb(0);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    CDenseMatrix mb(ma);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    CDenseMatrix mb(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3));
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    CDenseMatrix mb = ma;
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    CDenseMatrix mb = CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma.getNumberOfRows(), 2);
}

TEST_F(CDenseMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma.getNumberOfColumns(), 3);
}

TEST_F(CDenseMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma.getNumberOfElements(), 6);
}

TEST_F(CDenseMatrixTest, ValuesConstant)
{
    const CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, ma.values());
}

TEST_F(CDenseMatrixTest, Values)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, ma.values());
    
    auto mdat = ma.values();
    
    mdat[0] = 0.0; mdat[3] = 9.0;
    
    CDenseMatrix mb({0.0, -1.0, -3.0, 9.0, 5.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, RowConstant)
{
    const CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    vlxtest::compare({1.0, -1.0, -3.0}, ma.row(0));
    
    vlxtest::compare({-2.0, 5.0, 4.0}, ma.row(1));
    
    ASSERT_TRUE(ma.row(2) == nullptr);
    
    ASSERT_TRUE(ma.row(3) == nullptr);
}

TEST_F(CDenseMatrixTest, Row)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    vlxtest::compare({1.0, -1.0, -3.0}, ma.row(0));
    
    vlxtest::compare({-2.0, 5.0, 4.0}, ma.row(1));
    
    ASSERT_TRUE(ma.row(2) == nullptr);
    
    ASSERT_TRUE(ma.row(3) == nullptr);
    
    auto r0dat = ma.row(0);
    
    r0dat[0] = 8.0; r0dat[2] = 9.0;
    
    auto r1dat = ma.row(1);
    
    r1dat[1] = -6.0;
    
    CDenseMatrix mb({8.0, -1.0, 9.0, -2.0, -6.0, 4.0}, 2, 3);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, Zero)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 2, 3);
    
    ma.zero();
    
    CDenseMatrix mb({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 3);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, Symmetrize)
{
    CDenseMatrix ma({1.0, -1.0, -3.0,
                    -2.0,  5.0,  4.0,
                     1.0,  2.0,  3.0},
                    3, 3);
    
    ma.symmetrize();
    
    CDenseMatrix mb({2.0, -3.0, -2.0,
                    -3.0, 10.0,  6.0,
                    -2.0,  6.0,  6.0},
                    3, 3);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CDenseMatrixTest, Transpose)
{
    CDenseMatrix ma({ 1.0, -1.0,
                     -3.0, -2.0,
                      5.0,  4.0},
                      3, 2);
    
    CDenseMatrix mb({ 1.0, -3.0, 5.0,
                     -1.0, -2.0, 4.0},
                      2, 3);
    
    ASSERT_EQ(mb, ma.transpose());
}

TEST_F(CDenseMatrixTest, Slice)
{
    CDenseMatrix ma({1.0, -1.0, -3.0,
                    -2.0,  5.0,  4.0,
                     1.0,  2.0,  3.0,
                     2.0,  1.0,  4.0},
                    4, 3);
    
    CDenseMatrix refm14({-1.0, 5.0, 2.0, 1.0}, 4, 1);
    
    ASSERT_EQ(refm14, ma.slice(0, 1, 4, 1));
    
    CDenseMatrix refm32({ 5.0,  4.0,
                          2.0,  3.0,
                          1.0,  4.0},
                        3, 2);
    
    ASSERT_EQ(refm32, ma.slice(1, 1, 3, 2));
    
    CDenseMatrix refm21({1.0, -2.0}, 2, 1);
    
    ASSERT_EQ(refm21, ma.slice(0, 0, 2, 1));
}

TEST_F(CDenseMatrixTest, SelectByColumn)
{
    CDenseMatrix ma({ 1.0, -1.0, 3.0, 5.0,
                     -3.0, -2.0, 1.0, 2.0,
                      5.0,  4.0, 7.0, 5.0},
                    3, 4);
    
    CDenseMatrix refx({ 1.0, -1.0, 5.0,
                       -3.0, -2.0, 2.0,
                        5.0,  4.0, 5.0},
                      3, 3);
    
    ASSERT_EQ(refx, ma.selectByColumn({0, 1, 3}));
    
    CDenseMatrix refy({ -1.0, 5.0,
                        -2.0, 2.0,
                         4.0, 5.0},
                    3, 2);
    
    ASSERT_EQ(refy, ma.selectByColumn({1, 3}));
    
    ASSERT_EQ(CDenseMatrix(), ma.selectByColumn({})); 
}

TEST_F(CDenseMatrixTest, SelectByRow)
{
    CDenseMatrix ma({1.0, -1.0, -3.0,
                    -2.0,  5.0,  4.0,
                     1.0,  2.0,  3.0,
                     2.0,  1.0,  4.0},
                    4, 3);
    
    CDenseMatrix refx({1.0, -1.0, -3.0,
                      -2.0,  5.0,  4.0,
                       2.0,  1.0,  4.0},
                    3, 3);
    
    ASSERT_EQ(refx, ma.selectByRow({0, 1, 3}));
    
    CDenseMatrix refy({-2.0,  5.0,  4.0,
                        2.0,  1.0,  4.0},
                    2, 3);
    
    ASSERT_EQ(refy, ma.selectByRow({1, 3}));
    
    ASSERT_EQ(CDenseMatrix(), ma.selectByRow({}));
}

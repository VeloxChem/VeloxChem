//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SparseMatrixTest.hpp"

#include "SparseMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CSparseMatrixTest, DefaultConstructor)
{
    CSparseMatrix ma;
    
    CSparseMatrix mb(std::vector<double>{}, std::vector<int32_t>{},
                     std::vector<int32_t>{}, 0, 0, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, ConstructorWithDimensions)
{
    CSparseMatrix ma(2, 2, 4.0e-10);
   
    ASSERT_EQ(ma.getNumberOfRows(), 2);
    
    ASSERT_EQ(ma.getNumberOfColumns(), 2);
    
    ASSERT_EQ(ma.getNumberOfElements(), 0);
    
    ASSERT_NEAR(ma.getThreshold(), 4.0e-10, 1.0e-13);
}

TEST_F(CSparseMatrixTest, CopyConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                      {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb(ma);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, MoveConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb(CSparseMatrix({ 1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0,
                                    -4.0, 2.0, 7.0, 8.0, -5.0}, {0, 0, 0, 1, 1,
                                     2, 2, 2, 3, 3, 3, 4, 4}, {0, 1, 2, 0, 1, 2,
                                     3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13));
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, CopyAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb = ma;
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, MoveAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb = CSparseMatrix({ 1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0,
                                      -4.0, 2.0, 7.0, 8.0, -5.0}, {0, 0, 0, 1, 1,
                                       2, 2, 2, 3, 3, 3, 4, 4}, {0, 1, 2, 0, 1, 2,
                                       3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, AppendWithWidth)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                     8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5, 5},
                     {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6, 1.0e-13);
    
  
    CSparseMatrix mb(6, 6, 1.0e-13);
    
    CMemBlock<double> vecdat(6);
    
    CMemBlock<int32_t> vecidx(6);

    // 0-th row
    
    vecdat.zero();
    
    vecdat.at(0) = 1.0; vecdat.at(1) = -1.0; vecdat.at(2) = -3.0;
    
    vecidx.zero();
    
    vecidx.at(0) = 0; vecidx.at(1) = 1; vecidx.at(2) = 2;
    
    mb.append(vecdat, vecidx, 3, 0);
    
    // 1-th row
    
    vecdat.zero();
    
    vecdat.at(0) = -2.0; vecdat.at(1) = 5.0;
    
    vecidx.zero();
    
    vecidx.at(0) = 0; vecidx.at(1) = 1;
    
    mb.append(vecdat, vecidx, 2, 1);
    
    // 2-th row
    
    vecdat.zero();
    
    vecdat.at(0) = 4.0; vecdat.at(1) = 6.0; vecdat.at(2) = 4.0;
    
    vecidx.zero();
    
    vecidx.at(0) = 2; vecidx.at(1) = 3; vecidx.at(2) = 4;
    
    mb.append(vecdat, vecidx, 3, 2);
    
    // 3-th row
    
    vecdat.zero();
    
    vecdat.at(0) = -4.0; vecdat.at(1) = 2.0; vecdat.at(2) = 7.0;
    
    vecidx.zero();
    
    vecidx.at(0) = 0; vecidx.at(1) = 2; vecidx.at(2) = 3;
    
    mb.append(vecdat, vecidx, 3, 3);
    
    // 5-th row
    
    vecdat.zero();
    
    vecdat.at(0) = 8.0; vecdat.at(1) = -5.0;
    
    vecidx.zero();
    
    vecidx.at(0) = 1; vecidx.at(1) = 4;
    
    mb.append(vecdat, vecidx, 2, 5);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, Append)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                     8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5, 5},
                     {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6, 1.0e-13);
    
    
    CSparseMatrix mb(6, 6, 1.0e-13);
    
    // 0-th row
    
    CMemBlock<double> r0dat({1.0, -1.0, -3.0});
    
    CMemBlock<int32_t> r0idx({0, 1, 2});
    
    mb.append(r0dat, r0idx, 0);
    
    // 1-th row
    
    CMemBlock<double> r1dat({-2.0, 5.0});
    
    CMemBlock<int32_t> r1idx({0, 1});
    
    mb.append(r1dat, r1idx, 1);
    
    // 2-th row
    
    CMemBlock<double> r2dat({4.0, 6.0, 4.0});
    
    CMemBlock<int32_t> r2idx({2, 3, 4});
    
    mb.append(r2dat, r2idx, 2);
    
    // 3-th row
    
    CMemBlock<double> r3dat({-4.0, 2.0, 7.0});
    
    CMemBlock<int32_t> r3idx({0, 2, 3});
    
    mb.append(r3dat, r3idx, 3);
    
    // 5-th row
    
    CMemBlock<double> r5dat({8.0, -5.0});
    
    CMemBlock<int32_t> r5idx({1, 4});
    
    mb.append(r5dat, r5idx, 5);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, AppendWithReallocation)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     3, 3, 1.0e-13);
    
    
    CSparseMatrix mb(3, 3, 1.0e-13);
    
    // 0-th row
    
    CMemBlock<double> r0dat({1.0, -1.0, -3.0});
    
    CMemBlock<int32_t> r0idx({0, 1, 2});
    
    mb.append(r0dat, r0idx, 0);
    
    // 1-th row
    
    CMemBlock<double> r1dat({-2.0, 5.0, 4.0});
    
    CMemBlock<int32_t> r1idx({0, 1, 2});
    
    mb.append(r1dat, r1idx, 1);
    
    // 2-th row
    
    CMemBlock<double> r2dat({6.0, 4.0, -4.0});
    
    CMemBlock<int32_t> r2idx({0, 1, 2});
    
    mb.append(r2dat, r2idx, 2);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, Optimize_storage)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CSparseMatrix mb = ma;
    
    ma.optimize_storage();
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, IsOptimizedStorage)
{
    CSparseMatrix ma(6000, 6000, 1.0e-13);
    
    // 0-th row
    
    CMemBlock<double> r0dat({1.0, -1.0, -3.0});
    
    CMemBlock<int32_t> r0idx({0, 1, 2});
    
    ma.append(r0dat, r0idx, 3);
    
    ASSERT_FALSE(ma.isOptimizedStorage());
    
    ma.optimize_storage();
    
    CSparseMatrix mb({1.0, -1.0, -3.0},
                     {3, 3, 3}, {0, 1, 2},
                     6000, 6000, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
    
    ASSERT_TRUE(ma.isOptimizedStorage());
}

TEST_F(CSparseMatrixTest, GetNumberOfRows)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     3, 4, 1.0e-13);
    
    ASSERT_EQ(ma.getNumberOfRows(), 3);
}

TEST_F(CSparseMatrixTest, GetNumberOfColumns)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     3, 4, 1.0e-13);
    
    ASSERT_EQ(ma.getNumberOfColumns(), 4);
}

TEST_F(CSparseMatrixTest, GetNumberOfElements)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     3, 4, 1.0e-13);
    
    ASSERT_EQ(ma.getNumberOfElements(), 9);
}

TEST_F(CSparseMatrixTest, GetNumberOfElementsForRow)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5, 5},
                     {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6, 1.0e-13);

    ASSERT_EQ(ma.getNumberOfElements(0), 3);
    
    ASSERT_EQ(ma.getNumberOfElements(1), 2);
    
    ASSERT_EQ(ma.getNumberOfElements(2), 3);
    
    ASSERT_EQ(ma.getNumberOfElements(3), 3);
    
    ASSERT_EQ(ma.getNumberOfElements(4), 0);
    
    ASSERT_EQ(ma.getNumberOfElements(5), 2);
    
    ASSERT_EQ(ma.getNumberOfElements(6), 0);
    
    ASSERT_EQ(ma.getNumberOfElements(7), 0);
}

TEST_F(CSparseMatrixTest, GetThreshold)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     3, 4, 1.0e-8);
    
    ASSERT_NEAR(ma.getThreshold(), 1.0e-8, 1.0e-13);
}

TEST_F(CSparseMatrixTest, GetSparsity)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     {0, 0, 0, 1, 1, 1, 2, 2, 2}, {0, 1, 2, 0, 1, 2, 0, 1, 2},
                     6, 6, 1.0e-8);
    
    ASSERT_NEAR(ma.getSparsity(), 0.25, 1.0e-13);
}

TEST_F(CSparseMatrixTest,RowConstant)
{
    const CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                            7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                            5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                            1.0e-13);
    
    vlxtest::compare({1.0, -1.0, -3.0}, ma.row(0));
    
    vlxtest::compare({-2.0, 5.0}, ma.row(1));
    
    vlxtest::compare({4.0, 6.0, 4.0}, ma.row(2));
    
    vlxtest::compare({-4.0, 2.0, 7.0}, ma.row(3));
    
    ASSERT_TRUE(ma.row(4) == nullptr);
    
    vlxtest::compare({8.0, -5.0}, ma.row(5));
    
    ASSERT_TRUE(ma.row(6) == nullptr);
    
    ASSERT_TRUE(ma.row(7) == nullptr);
}

TEST_F(CSparseMatrixTest, Row)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                      1.0e-13);
    
    vlxtest::compare({1.0, -1.0, -3.0}, ma.row(0));
    
    vlxtest::compare({-2.0, 5.0}, ma.row(1));
    
    vlxtest::compare({4.0, 6.0, 4.0}, ma.row(2));
    
    vlxtest::compare({-4.0, 2.0, 7.0}, ma.row(3));
    
    ASSERT_TRUE(ma.row(4) == nullptr);
    
    vlxtest::compare({8.0, -5.0}, ma.row(5));
    
    ASSERT_TRUE(ma.row(6) == nullptr);
    
    ASSERT_TRUE(ma.row(7) == nullptr);
    
    auto r2dat = ma.row(2);
    
    r2dat[1] = 7.0;
    
    CSparseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 7.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                     1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, IndexesConstant)
{
    const CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                            7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                            5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                            1.0e-13);
    
    vlxtest::compare({0, 1, 2}, ma.indexes(0));
    
    vlxtest::compare({0, 1}, ma.indexes(1));
    
    vlxtest::compare({2, 3, 4}, ma.indexes(2));
    
    vlxtest::compare({0, 2, 3}, ma.indexes(3));
    
    ASSERT_TRUE(ma.indexes(4) == nullptr);
    
    vlxtest::compare({1, 4}, ma.indexes(5));
    
    ASSERT_TRUE(ma.indexes(6) == nullptr);
    
    ASSERT_TRUE(ma.indexes(7) == nullptr);
}

TEST_F(CSparseMatrixTest, GetRowIndexes)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6, 1.0e-13);
    
    vlxtest::compare({0, 1, 2}, ma.indexes(0));
    
    vlxtest::compare({0, 1}, ma.indexes(1));
    
    vlxtest::compare({2, 3, 4}, ma.indexes(2));
    
    vlxtest::compare({0, 2, 3}, ma.indexes(3));
    
    ASSERT_TRUE(ma.indexes(4) == nullptr);
    
    vlxtest::compare({1, 4}, ma.indexes(5));
    
    ASSERT_TRUE(ma.indexes(6) == nullptr);
    
    ASSERT_TRUE(ma.indexes(7) == nullptr);
    
    auto r2idx = ma.indexes(2);
    
    r2idx[2] = 5;
    
    CSparseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 5, 0, 2, 3, 1, 4}, 6, 6, 1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, ValuesConstant)
{
    const CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                            7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                            5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                            1.0e-13);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, ma.values());
}

TEST_F(CSparseMatrixTest, Values)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                      1.0e-13);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, ma.values());
    
    auto mdat = ma.values();
    
    mdat[1] = 2.0; mdat[6] = -1.0;
    
    CSparseMatrix mb({1.0, 2.0, -3.0, -2.0, 5.0, 4.0, -1.0, 4.0, -4.0, 2.0,
                      7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                      5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                      1.0e-13);
    
    ASSERT_EQ(ma, mb);
}

TEST_F(CSparseMatrixTest, Rows)
{
    const CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                            7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                            5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                            1.0e-13);
    
    vlxtest::compare({0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5, 5}, ma.rows());
}

TEST_F(CSparseMatrixTest, Columns)
{
    const CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0,
                            7.0, 8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 5,
                            5}, {0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, 6, 6,
                            1.0e-13);
    
    vlxtest::compare({0, 1, 2, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4}, ma.columns());
}



//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GenFuncTest.hpp"

#include "MemBlock2D.hpp"
#include "GenFunc.hpp"

TEST_F(CGenFuncTest, Contract)
{
    int32_t sposvec[3] = {0, 2, 3};
    
    int32_t eposvec[3] = {2, 3, 6};
    
    CMemBlock2D<double> pdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                              2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                              2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                             6, 3);
    
    CMemBlock2D<double> cdat(3, 2);
    
    genfunc::contract(cdat, 0, pdat, 1, sposvec, eposvec, 3, 2);
    
    CMemBlock2D<double> tdat({5.0, 6.0, 16.0, 8.1, -1.0, 17.0}, 3, 2);
    
    ASSERT_EQ(cdat, tdat);
}

TEST_F(CGenFuncTest, TransformForS)
{
    CMemBlock2D<double> cartdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 6, 3);
    
    CMemBlock2D<double> spherdat(6, 3);
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(0), 6, 3);
    
    ASSERT_EQ(spherdat, cartdat);
}

TEST_F(CGenFuncTest, TransformForP)
{
    CMemBlock2D<double> cartdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 6, 3);
    
    CMemBlock2D<double> spherdat(6, 3);
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(1), 6, 1);
    
    CMemBlock2D<double> tdat({2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                              2.4, 5.7, -1.0, 8.0,  9.0, 0.0,
                              1.0, 2.0,  3.0, 6.0, -3.0, 4.0},
                             6, 3);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, TransformForD)
{
    CMemBlock2D<double> cartdat({1.0, 2.0, 3.0, 4.0,
                                 1.0, 3.0, 3.0, 5.0,
                                 0.0, 1.0, 4.0, 2.0,
                                 2.0, 1.0, 2.0, 1.0,
                                 4.0, 3.0, 4.0, 5.0,
                                 1.0, 1.0, 2.0, 3.0},
                                2, 12);
    
    CMemBlock2D<double> spherdat(2, 10);
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(2), 2, 2);
    
    CMemBlock2D<double> tdat({ 3.46410161513775, 10.392304845413260,
                              10.39230484541326, 17.320508075688770,
                              13.85640646055100, 10.392304845413300,
                              13.85640646055100, 17.320508075688800,
                              -1.00000000000000, -1.000000000000000,
                              -1.00000000000000,  1.000000000000000,
                               0.00000000000000,  3.464101615137754,
                             13.856406460551018,  6.928203230275509,
                              -1.73205080756888,   1.73205080756888,
                               1.73205080756888,   5.19615242270663},
                             2, 10);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, IsInVectorForTwoIndexes)
{
    CVecTwoIndexes vec{{0, 1}, {2, 4}, {7, -1}, {8, 2}};
 
    ASSERT_TRUE(genfunc::isInVector(vec, {0, 1}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {2, 4}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {7, -1}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {8,  2}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {0, 0}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {1, 0}));

    ASSERT_FALSE(genfunc::isInVector(vec, {1, 1}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {7, 1}));
}

TEST_F(CGenFuncTest, AddValidAndUniquePair)
{
    CVecTwoIndexes vec{{0, 1}, {2, 4}};
    
    genfunc::addValidAndUniquePair(vec, {7, -1});
    
    genfunc::addValidAndUniquePair(vec, {2, 1});
    
    genfunc::addValidAndUniquePair(vec, {2, 4});
    
    genfunc::addValidAndUniquePair(vec, {2, 1});
    
    ASSERT_EQ(vec.size(), 3);
    
    ASSERT_EQ(vec[0], CTwoIndexes(0, 1));
    
    ASSERT_EQ(vec[1], CTwoIndexes(2, 4));
    
    ASSERT_EQ(vec[2], CTwoIndexes(2, 1));
}

TEST_F(CGenFuncTest, FindPairIndex)
{
    CVecTwoIndexes vec{{0, 1}, {2, 4}, {7, -1}, {8, 2}};
    
    std::vector<int32_t> idx{1, 3, 7, 8};
    
    ASSERT_EQ(8, genfunc::findPairIndex(idx, vec, {8, 2}));
    
    ASSERT_EQ(7, genfunc::findPairIndex(idx, vec, {7, -1}));
    
    ASSERT_EQ(3, genfunc::findPairIndex(idx, vec, {2, 4}));
    
    ASSERT_EQ(1, genfunc::findPairIndex(idx, vec, {0, 1}));
    
    ASSERT_EQ(-1, genfunc::findPairIndex(idx, vec, {0, 0}));
    
    ASSERT_EQ(-1, genfunc::findPairIndex(idx, vec, {2, 3}));
}

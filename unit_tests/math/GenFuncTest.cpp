//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GenFuncTest.hpp"

#include "MemBlock2D.hpp"
#include "GenFunc.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MemBlock.hpp"
#include "GtoContainer.hpp"

TEST_F(CGenFuncTest, Contract)
{
    int32_t sposvec[3] = {0, 2, 3};
    
    int32_t eposvec[3] = {2, 3, 6};
    
    CMemBlock2D<double> pdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                              2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                              2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                             6, 3);
    
    CMemBlock2D<double> cdat(3, 3);
    
    cdat.zero();
    
    genfunc::contract(cdat, pdat, 1,  1, sposvec, eposvec, 3, 2);
    
    CMemBlock2D<double> tdat({0.0, 0.0, 0.0, 5.0, 6.0, 16.0, 8.1, -1.0, 17.0}, 3, 3);
    
    ASSERT_EQ(cdat, tdat);
}

TEST_F(CGenFuncTest, ContractWithGtoBlocks)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 1);
    
    CGtoBlock kgtos(lih, bas, 1);
    
    CMemBlock2D<double> pdat({ 1.0,  2.0,  3.0, 6.0,
                              -3.0,  4.0,  2.0, 3.0,
                               6.0,  7.0,  8.0, 1.0,
                               2.4,  5.7, -1.0, 8.0,
                               9.0,  0.0,  2.0, 3.0,
                               3.0,  4.0,  5.0, 6.0,
                               1.0,  2.0,  3.0, 7.0,
                               0.0, -1.0,  2.0, 2.0,
                               0.1, -1.2,  0.0, 2.6,
                               // second data block
                               0.2,  0.8,  1.3, 2.0,
                               0.9,  0.3,  2.5, 1.0,
                               3.0,  1.0,  4.0, 3.0,
                              -1.4,  2.7,  3.0, 0.5,
                               2.0,  3.0,  4.0, 5.0,
                               3.0,  2.0,  1.0, 2.0,
                               0.0,  4.0,  2.0, 0.0,
                               4.0,  0.0,  5.0, 1.0,
                               0.9,  1.5,  2.0, 2.6},
                             4, 18);
    
    CMemBlock2D<double> cdat(3, 9);
    
    genfunc::contract(cdat, pdat, 0, bgtos, kgtos, 0);
    
    CMemBlock2D<double> tdat({ 4.0,  4.3, 8.0,
                               2.2,  4.5, 4.0,
                              17.0, 12.0, 4.0,
                               9.4,  2.0, 8.5,
                              14.0,  6.0, 8.0,
                              12.0,  6.0, 8.0,
                               7.0,  5.0, 7.0,
                               3.0,  7.0, 3.0,
                               1.3,  2.0, 5.2},
                             3, 9);
    
    ASSERT_EQ(cdat, tdat);
}

TEST_F(CGenFuncTest, ContractWithGtoBlockAndGtoPairsBlock)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> pdat({1.0, 0.2, 3.0, 2.3, 0.0, 1.0, 0.7, 0.2, 3.1, 1.2, 0.8,
                              0.3, 0.5, 2.1, 0.5, 0.9, 2.1, 1.7, 3.4, 3.9, 2.7, 1.2,
                              3.1, 0.8, 0.4, 2.1, 2.7, 0.6, 2.0, 2.7, 1.7, 0.6, 0.3,
                              1.2, 2.6, 0.8, 1.3, 2.9, 2.5, 1.7, 2.2, 2.1, 2.0, 1.5,
                              0.9, 3.2, 2.1, 2.6, 2.7, 0.3, 2.2, 0.9, 1.1, 2.7, 0.2,
                              0.8, 1.5, 0.8, 1.5, 1.4, 1.9, 2.7, 4.2, 0.5, 2.3, 2.8,
                              0.1, 1.8, 2.2, 0.5, 0.8, 2.2, 3.7, 1.2, 2.9, 1.7, 1.3,
                              4.1, 3.6, 0.7, 1.5, 3.2, 4.0, 1.7, 3.2, 0.7, 1.6, 2.1},
                              11, 8);
    
    CVecThreeIndexes cvec({{0, 0, 0}, {0, 1, 0}});
    
    std::vector<int32_t> cidx({3, 0});
    
    CVecThreeIndexes pvec({{0, 2, 3}, {0, 0, 0}, {0, 0, 1},  {0, 1, 0}});
    
    std::vector<int32_t> pidx({7, 0, 8, 2});
    
    CMemBlock2D<double> cdat(6, 4);
    
    genfunc::contract(cdat, pdat, cvec, cidx, pvec, pidx, pbra, ppairs, 0);
    
    CMemBlock2D<double> tdat({11.0,  6.6, 11.6, 2.2, 2.9, 3.1,
                              10.5,  8.4,  8.8, 5.0, 3.7, 2.8,
                              18.7, 10.2,  8.0, 1.8, 4.3, 2.3,
                               9.9,  4.0,  6.0, 7.0, 3.9, 2.0},
                              6, 4);
    
    ASSERT_EQ(cdat, tdat);
}

TEST_F(CGenFuncTest, ContractWithGtoPairsBlocks)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> pdat({1.0, 0.2, 3.0, 2.3, 0.0, 1.0, 0.7, 0.2, 3.1, 1.2, 0.8,
                              0.3, 0.5, 2.1, 0.5, 0.9, 2.1, 1.7, 3.4, 3.9, 2.7, 1.2,
                              3.1, 0.8, 0.4, 2.1, 2.7, 0.6, 2.0, 2.7, 1.7, 0.6, 0.3,
                              1.2, 2.6, 0.8, 1.3, 2.9, 2.5, 1.7, 2.2, 2.1, 2.0, 1.5,
                              0.9, 3.2, 2.1, 2.6, 2.7, 0.3, 2.2, 0.9, 1.1, 2.7, 0.2,
                              0.8, 1.5, 0.8, 1.5, 1.4, 1.9, 2.7, 4.2, 0.5, 2.3, 2.8,
                              0.1, 1.8, 2.2, 0.5, 0.8, 2.2, 3.7, 1.2, 2.9, 1.7, 1.3,
                              4.1, 3.6, 0.7, 1.5, 3.2, 4.0, 1.7, 3.2, 0.7, 1.6, 2.1},
                             11, 8);
    
    CVecThreeIndexes cvec({{0, 0, 0}, {0, 3, 1}, {0, 0, 1}});
    
    std::vector<int32_t> cidx({3, 2, 0});
    
    CVecThreeIndexes pvec({{0, 2, 3}, {0, 0, 0}, {0, 0, 1},  {0, 1, 0}});
    
    std::vector<int32_t> pidx({7, 0, 8, 2});
    
    auto p2dat = pdat; 
    
    CMemBlock2D<double> c2dat(2, 4);
    
    genfunc::contract(c2dat, p2dat, cvec, cidx, pvec, pidx, ppairs, ppairs,
                      6, 2, 1);
    
    CMemBlock2D<double> t2dat({11.0,  6.6,
                              10.5,  8.4,
                              18.7, 10.2,
                               9.9,  4.0},
                              2, 4);
    
    ASSERT_EQ(c2dat, t2dat);
    
    CMemBlock2D<double> cdat(6, 4);
    
    genfunc::contract(cdat, pdat, cvec, cidx, pvec, pidx, ppairs, ppairs,
                      11, 6, 1);
    
    CMemBlock2D<double> tdat({11.0,  6.6, 11.6, 2.2, 2.9, 3.1,
                              10.5,  8.4,  8.8, 5.0, 3.7, 2.8,
                              18.7, 10.2,  8.0, 1.8, 4.3, 2.3,
                               9.9,  4.0,  6.0, 7.0, 3.9, 2.0},
                             6, 4);
    
    ASSERT_EQ(cdat, tdat);
}

TEST_F(CGenFuncTest, TransformForS)
{
    CMemBlock2D<double> cartdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 6, 3);
    
    CMemBlock2D<double> spherdat(6, 3);
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(0), 0, 0, 6, 3);
    
    ASSERT_EQ(spherdat, cartdat);
}

TEST_F(CGenFuncTest, TransformForP)
{
    CMemBlock2D<double> cartdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 6, 3);
    
    CMemBlock2D<double> spherdat(6, 3);
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(1), 0, 0, 6, 1);
    
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
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(2), 0, 0, 2, 2);
    
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

TEST_F(CGenFuncTest, TransformForSPAndPS)
{
    CMemBlock2D<double> cartdat({1.0, 2.0,  3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 6, 3);
    
    CMemBlock2D<double> spherdat(6, 3);
    
    CSphericalMomentum bmom(0);
    
    CSphericalMomentum kmom(1);
    
    genfunc::transform(spherdat, cartdat, bmom, kmom, 0, 0, 6);
    
    CMemBlock2D<double> tdat({2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                              2.4, 5.7, -1.0, 8.0,  9.0, 0.0,
                              1.0, 2.0,  3.0, 6.0, -3.0, 4.0},
                              6, 3);
    
    ASSERT_EQ(spherdat, tdat);
    
    genfunc::transform(spherdat, cartdat, kmom, bmom, 0, 0, 6);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, TransformForSDAndDS)
{
    CMemBlock2D<double> cartdat({1.0, 2.0, 3.0, 4.0,
                                 1.0, 3.0, 3.0, 5.0,
                                 0.0, 1.0, 4.0, 2.0,
                                 2.0, 1.0, 2.0, 1.0,
                                 4.0, 3.0, 4.0, 5.0,
                                 1.0, 1.0, 2.0, 3.0},
                                4, 6);
    
    CMemBlock2D<double> spherdat(4, 5);
    
    CSphericalMomentum bmom(0);
    
    CSphericalMomentum kmom(2);
    
    genfunc::transform(spherdat, cartdat, bmom, kmom, 0, 0, 4);
    
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
                             4, 5);
    
    ASSERT_EQ(spherdat, tdat);
    
    genfunc::transform(spherdat, cartdat, kmom, bmom, 0, 0, 4);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, TransformForDD)
{
    CMemBlock2D<double> cartdat(1, 36);
    
    for (int32_t i = 0; i < 36; i++)
    {
        auto vec = cartdat.data(i);
        
        vec[0] = 1.0;
    }
    
    CMemBlock2D<double> spherdat(1, 25);
    
    CSphericalMomentum bmom(2);
    
    CSphericalMomentum kmom(2);
    
    genfunc::transform(spherdat, cartdat, bmom, kmom, 0, 0, 1);
    
    CMemBlock2D<double> tdat({12.0, 12.0, 0.0, 12.0, 0.0,
                              12.0, 12.0, 0.0, 12.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0,
                              12.0, 12.0, 0.0, 12.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0},
                             1, 25);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, TransformHalfIntegrals)
{
    CMemBlock2D<double> cartdat({1.0, 2.0, 3.0, 6.0, -3.0, 4.0,
                                 2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                                 2.4, 5.7, -1.0, 8.0,  9.0, 0.0},
                                 2, 9);
    
    CMemBlock2D<double> spherdat(2, 12);
    
    spherdat.zero();
    
    CVecThreeIndexes sphervec({{1, 0, 1}, {0, 2, 2}});
    
    std::vector<int32_t> spheridx({3, 0});
    
    CVecThreeIndexes cartvec({{1, 1, 0}});
    
    std::vector<int32_t> cartidx({0});
    
    genfunc::transform(spherdat, cartdat, CSphericalMomentum(1), sphervec,
                       spheridx, cartvec, cartidx, 2);
    
    CMemBlock2D<double> tdat({0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
                              2.0, 3.0,  6.0, 7.0,  8.0, 1.0,
                              2.4, 5.7, -1.0, 8.0,  9.0, 0.0,
                              1.0, 2.0,  3.0, 6.0, -3.0, 4.0},
                             2, 12);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, TransformFullIntegralsForDD)
{
    CMemBlock2D<double> cartdat(1, 72);
    
    for (int32_t i = 0; i < 36; i++)
    {
        auto vec = cartdat.data(i);
        
        vec[0] = 1.0;
    }
    
    for (int32_t i = 36; i < 72; i++)
    {
        auto vec = cartdat.data(i);
        
        vec[0] = 2.0;
    }
    
    CMemBlock2D<double> spherdat(1, 50);
    
    CSphericalMomentum bmom(2);
    
    CSphericalMomentum kmom(2);
    
    genfunc::transform(spherdat, cartdat, bmom, kmom, 0, 0, 1, 2);
    
    CMemBlock2D<double> tdat({12.0, 12.0, 0.0, 12.0, 0.0,
                              12.0, 12.0, 0.0, 12.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0,
                              12.0, 12.0, 0.0, 12.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0,
                              24.0, 24.0, 0.0, 24.0, 0.0,
                              24.0, 24.0, 0.0, 24.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0,
                              24.0, 24.0, 0.0, 24.0, 0.0,
                               0.0,  0.0, 0.0,  0.0, 0.0},
                             1, 50);
    
    ASSERT_EQ(spherdat, tdat);
}

TEST_F(CGenFuncTest, Compress)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 1);
    
    CGtoBlock kgtos(lih, bas, 1);
    
    CMemBlock2D<double> adat({ 1.0,  2.0,  0.0,
                              -3.0,  4.0,  2.0,
                               6.0,  7.0,  8.0,
                               1.0,  5.7, -1.0,
                               0.0,  0.0,  2.0,
                               0.0,  0.0,  0.0,
                               1.0,  2.0,  3.0,
                               0.0, -1.0,  2.0,
                               0.1, -1.2,  0.0},
                             3, 9);
    
    CSparseMatrix spmat(6, 9, 1.0e-13);
    
    CMemBlock<double> rvals(9);
    
    CMemBlock<int32_t> ridx(9);
    
    genfunc::compress(spmat, rvals, ridx, adat, bgtos, kgtos, 0);
    
    CSparseMatrix amat({1.0,  2.0, -3.0,  4.0, 2.0, 6.0,  7.0,  8.0,
                        1.0,  5.7, -1.0, 2.0,
                        1.0,  2.0,  3.0,  -1.0,  2.0, 0.1, -1.2},
                       {0, 0, 0, 0, 0, 0, 0, 0,
                        1, 1, 1, 1,
                        2, 2, 2, 2, 2, 2, 2},
                       {5, 6, 8, 9, 10, 11, 12, 13,
                        5, 6, 7, 10,
                        5, 6, 7, 9, 10, 11, 12},
                        6, 9, 1.0e-13);
    
    ASSERT_EQ(spmat, amat);
    
    CMemBlock2D<double> bdat({ 1.0,  2.0,  0.0,
                               0.0,  4.0,  0.0,
                               0.0,  7.0,  8.0,
                               0.0,  0.0, -1.0,
                               0.0,  0.0,  2.0,
                               0.0,  0.0,  0.0,
                               1.0,  2.0,  3.0,
                               0.0, -1.0,  2.0,
                               0.1, -1.2,  0.0},
                             3, 9);
    
    genfunc::compress(spmat, rvals, ridx, bdat, bgtos, kgtos, 1);
    
    CSparseMatrix bmat({ 1.0,  2.0, -3.0,  4.0, 2.0, 6.0,  7.0,  8.0,
                         1.0,  5.7, -1.0,  2.0,
                         1.0,  2.0,  3.0, -1.0,  2.0, 0.1, -1.2,
                         1.0,  2.0,  4.0,  7.0,  8.0,
                        -1.0,  2.0,
                         1.0,  2.0,  3.0, -1.0,  2.0, 0.1, -1.2},
                       {0, 0, 0, 0, 0, 0, 0, 0,
                        1, 1, 1, 1,
                        2, 2, 2, 2, 2, 2, 2,
                        3, 3, 3, 3, 3,
                        4, 4,
                        5, 5, 5, 5, 5, 5, 5},
                       {5, 6, 8, 9, 10, 11, 12, 13,
                        5, 6, 7, 10,
                        5, 6, 7, 9, 10, 11, 12,
                        5, 6, 9, 12, 13,
                        7, 10,
                        5, 6, 7, 9, 10, 11, 12},
                       6, 9, 1.0e-13);
    
    ASSERT_EQ(spmat, bmat);
}

TEST_F(CGenFuncTest, Distribute)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer bcont(lih, bas);
    
    CGtoContainer kcont(lih, bas);
    
    // set up (s|s) integrals matrix
    
    CSparseMatrix matss({1.0, 2.0,  5.0,
                         1.0, 7.2, -1.0, 3.0,
                         2.0, 3.0, 4.0,
                         1.0, 2.0},
                        { 0, 0, 0,
                          1, 1, 1, 1,
                          3, 3, 3,
                          4, 4},
                        { 0, 2, 4,
                          0, 1, 3, 4,
                          1, 2, 3,
                          3, 4},
                        5, 5, 1.0e-13);
    
    // set up (s|p) integrals matrix
    
    CSparseMatrix matsp({2.0, 3.0, 4.0, 5.0, 6.0,
                         2.0, 1.2, 1.0, 2.0,
                         9.0, 7.0, 8.0,
                         3.0, 6.0, 7.0, 1.0},
                        { 0, 0, 0, 0, 0,
                          2, 2, 2, 2,
                          3, 3, 3,
                          4, 4, 4, 4},
                        { 4, 7, 10, 11, 12,
                          4, 5, 6, 7,
                          8, 9, 10,
                          8, 9, 10, 11},
                         5, 9, 1.0e-13);
    
    // set up (p|s) integrals matrix
    
    CSparseMatrix matps({2.0, 3.0, 4.0,
                         2.0, 1.2, 2.0,
                         9.0, 7.0, 8.0,
                         3.0, 6.0, 7.0, 1.0,
                         1.0, 4.0, 6.0, 1.0,
                         2.0, 2.0, 4.0, 5.0,
                         1.2, 3.0,
                         4.0, 5.0,
                         6.7, 1.9},
                        { 0, 0, 0,
                          1, 1, 1,
                          2, 2, 2,
                          3, 3, 3, 3,
                          4, 4, 4, 4,
                          5, 5, 5, 5,
                          6, 6,
                          7, 7,
                          8, 8},
                        { 0, 1, 2,
                          2, 3, 4,
                          0, 3, 4,
                          0, 1, 2, 3,
                          0, 2, 3, 4,
                          1, 2, 3, 4,
                          0, 4,
                          1, 2,
                          3, 4},
                        9, 5, 1.0e-13);
    
    // set up (p|p) integrals matrix
    
    CSparseMatrix matpp({1.0, 2.0, 3.0,
                         4.0, 5.0, 6.0,
                         7.0, 8.0, 9.0},
                        {0, 1, 2,
                         3, 4, 5,
                         6, 7, 8},
                        { 5,  6,  7,
                          8,  9, 10,
                         11, 12, 13},
                        9, 9, 1.0e-13);
    
    // set up list of matrices
    
    CSparseMatrix* matlst = new CSparseMatrix[4];
    
    matlst[0] = matss; matlst[1] = matsp;
    
    matlst[2] = matps; matlst[3] = matpp;
    
    auto spmat = genfunc::distribute(matlst, &bcont, &kcont);
    
    delete [] matlst;
    
    CSparseMatrix tmat({1.0, 2.0,  5.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                        1.0, 7.2, -1.0, 3.0,
                        2.0, 1.2,  1.0, 2.0,
                        2.0, 3.0,  4.0, 9.0, 7.0, 8.0,
                        1.0, 2.0,  3.0, 6.0, 7.0, 1.0,
                        2.0, 3.0, 4.0, 1.0,
                        3.0, 6.0, 7.0, 1.0, 4.0,
                        1.2, 3.0, 7.0,
                        2.0, 1.2, 2.0, 2.0,
                        1.0, 4.0, 6.0, 1.0, 5.0,
                        4.0, 5.0, 8.0,
                        9.0, 7.0, 8.0, 3.0,
                        2.0, 2.0, 4.0, 5.0, 6.0,
                        6.7, 1.9, 9.0},
                       {0, 0, 0, 0, 0, 0, 0, 0,
                        1, 1, 1, 1,
                        2, 2, 2, 2,
                        3, 3, 3, 3, 3, 3,
                        4, 4, 4, 4, 4, 4,
                        5, 5, 5, 5,
                        6, 6, 6, 6, 6,
                        7, 7, 7,
                        8, 8, 8, 8,
                        9, 9, 9, 9, 9,
                        10, 10, 10,
                        11, 11, 11, 11,
                        12, 12, 12, 12, 12,
                        13, 13, 13},
                       {0, 2, 4, 4, 7, 10, 11, 12,
                        0, 1, 3, 4,
                        4, 5, 6, 7,
                        1, 2, 3, 8, 9, 10,
                        3, 4, 8, 9, 10, 11,
                        0, 1, 2, 5,
                        0, 1, 2, 3, 8,
                        0, 4, 11,
                        2, 3, 4, 6,
                        0, 2, 3, 4, 9,
                        1, 2, 12,
                        0, 3, 4, 7,
                        1, 2, 3, 4, 10,
                        3, 4, 13},
                        14, 14, 1.0e-13);
    
    ASSERT_EQ(spmat, tmat); 
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

TEST_F(CGenFuncTest, IsInVectorForThreeIndexes)
{
    CVecThreeIndexes vec{{0, 1, 3}, {2, 4, 2}, {7, -1, 3}, {8, 2, 4}};
    
    ASSERT_TRUE(genfunc::isInVector(vec, {0, 1, 3}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {2, 4, 2}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {7, -1, 3}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {8,  2, 4}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {0, 0, 0}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {1, 0, 1}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {1, 1, 2}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {7, 1, 2}));
}

TEST_F(CGenFuncTest, IsInVectorForFourIndexes)
{
    CVecFourIndexes vec{{0, 1, 3, 4}, {2, 4, 2, 1}, {7, -1, 3, 4}, {8, 2, 4, 5}};
    
    ASSERT_TRUE(genfunc::isInVector(vec, {0, 1, 3, 4}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {2, 4, 2, 1}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {7, -1, 3, 4}));
    
    ASSERT_TRUE(genfunc::isInVector(vec, {8,  2, 4, 5}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {0, 0, 0, 0}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {1, 0, 1, 1}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {1, 1, 2, 1}));
    
    ASSERT_FALSE(genfunc::isInVector(vec, {7, 1, 2, 3}));
}

TEST_F(CGenFuncTest, AddValidAndUniquePair)
{
    CVecTwoIndexes vec{{0, 1}, {2, 4}};
    
    genfunc::addValidAndUniquePair(vec, {7, -1});
    
    genfunc::addValidAndUniquePair(vec, {2, 1});
    
    genfunc::addValidAndUniquePair(vec, {2, 4});
    
    genfunc::addValidAndUniquePair(vec, {2, 1});
    
    ASSERT_EQ(vec.size(), 3u);
    
    ASSERT_EQ(vec[0], CTwoIndexes(0, 1));
    
    ASSERT_EQ(vec[1], CTwoIndexes(2, 4));
    
    ASSERT_EQ(vec[2], CTwoIndexes(2, 1));
}

TEST_F(CGenFuncTest, AddValidAndUniqueTriple)
{
    CVecThreeIndexes vec{{0, 1, 1}, {2, 4, 0}};
    
    genfunc::addValidAndUniqueTriple(vec, {7, -1, 1});
    
    genfunc::addValidAndUniqueTriple(vec, {2, 1, 1});
    
    genfunc::addValidAndUniqueTriple(vec, {2, 4, 0});
    
    genfunc::addValidAndUniqueTriple(vec, {2, 1, -1});
    
    genfunc::addValidAndUniqueTriple(vec, {2, 1, 1});
    
    ASSERT_EQ(vec.size(), 3u);
    
    ASSERT_EQ(vec[0], CThreeIndexes(0, 1, 1));
    
    ASSERT_EQ(vec[1], CThreeIndexes(2, 4, 0));
    
    ASSERT_EQ(vec[2], CThreeIndexes(2, 1, 1));
}

TEST_F(CGenFuncTest, AddValidAndUniqueQuadruple)
{
    CVecFourIndexes vec{{0, 1, 1, 2}, {2, 4, 0, 7}};
    
    genfunc::addValidAndUniqueQuadruple(vec, {7, -1, 1, 2});
    
    genfunc::addValidAndUniqueQuadruple(vec, {2, 1, 1, 3});
    
    genfunc::addValidAndUniqueQuadruple(vec, {2, 4, 0, 7});
    
    genfunc::addValidAndUniqueQuadruple(vec, {2, 1, -1, 1});
    
    genfunc::addValidAndUniqueQuadruple(vec, {2, 1, 1, 2});
    
    genfunc::addValidAndUniqueQuadruple(vec, {-2, 1, 1, 2});
    
    genfunc::addValidAndUniqueQuadruple(vec, {2, 1, 1, -1});
    
    ASSERT_EQ(vec.size(), 4u);
    
    ASSERT_EQ(vec[0], CFourIndexes(0, 1, 1, 2));
    
    ASSERT_EQ(vec[1], CFourIndexes(2, 4, 0, 7));
    
    ASSERT_EQ(vec[2], CFourIndexes(2, 1, 1, 3));
    
    ASSERT_EQ(vec[3], CFourIndexes(2, 1, 1, 2));
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

TEST_F(CGenFuncTest, FindTripleIndex)
{
    CVecThreeIndexes vec{{0, 1, 1}, {2, 4, 3}, {7, -1, 2}, {8, 2, 4}};
    
    std::vector<int32_t> idx{1, 3, 7, 8};
    
    ASSERT_EQ(8, genfunc::findTripleIndex(idx, vec, {8, 2, 4}));
    
    ASSERT_EQ(7, genfunc::findTripleIndex(idx, vec, {7, -1, 2}));
    
    ASSERT_EQ(3, genfunc::findTripleIndex(idx, vec, {2, 4, 3}));
    
    ASSERT_EQ(1, genfunc::findTripleIndex(idx, vec, {0, 1, 1}));
    
    ASSERT_EQ(-1, genfunc::findTripleIndex(idx, vec, {0, 0, 0}));
    
    ASSERT_EQ(-1, genfunc::findTripleIndex(idx, vec, {2, 4, 2}));
}

TEST_F(CGenFuncTest, FindQuadrupleIndex)
{
    CVecFourIndexes vec{{0, 1, 1, 2}, {2, 4, 3, 1}, {7, -1, 2, 0}, {8, 2, 4, 7}};
    
    std::vector<int32_t> idx{1, 3, 7, 8};
    
    ASSERT_EQ(8, genfunc::findQuadrupleIndex(idx, vec, {8, 2, 4, 7}));
    
    ASSERT_EQ(7, genfunc::findQuadrupleIndex(idx, vec, {7, -1, 2, 0}));
    
    ASSERT_EQ(3, genfunc::findQuadrupleIndex(idx, vec, {2, 4, 3, 1}));
    
    ASSERT_EQ(1, genfunc::findQuadrupleIndex(idx, vec, {0, 1, 1, 2}));
    
    ASSERT_EQ(-1, genfunc::findQuadrupleIndex(idx, vec, {0, 0, 0, 0}));
    
    ASSERT_EQ(-1, genfunc::findQuadrupleIndex(idx, vec, {2, 4, 2, 2}));
}

TEST_F(CGenFuncTest, MaxOrderOfPair)
{
    CVecThreeIndexes vec{{0, 1, 1}, {2, 4, 3}, {0, 1, 0}, {0, 1, 4}};
    
    ASSERT_EQ(-1, genfunc::maxOrderOfPair(vec, 0, 0));
    
    ASSERT_EQ( 4, genfunc::maxOrderOfPair(vec, 0, 1));
    
    ASSERT_EQ(-1, genfunc::maxOrderOfPair(vec, 4, 2));
    
    ASSERT_EQ( 3, genfunc::maxOrderOfPair(vec, 2, 4));
}

TEST_F(CGenFuncTest, MaxOrderOfPairForFourIndexes)
{
    CVecFourIndexes vec{{0, 1, 1, 2}, {0, 1, 3, 1}, {7, -1, 2, 0}, {0, 1, 4, 7}};
    
    ASSERT_EQ(-1, genfunc::maxOrderOfPair(vec, 0, 0));
    
    ASSERT_EQ( 4, genfunc::maxOrderOfPair(vec, 0, 1));
    
    ASSERT_EQ(-1, genfunc::maxOrderOfPair(vec, 4, 2));
    
    ASSERT_EQ( 2, genfunc::maxOrderOfPair(vec, 7, -1));
}

TEST_F(CGenFuncTest, GetPairsFromTripleIndexes)
{
    CVecThreeIndexes vec{{0, 0, 1}, {2, 4, 3}, {7, 0, 2}, {8, -2, 4}};
    
    auto xyvec = genfunc::getPairsFromTripleIndexes(vec);
    
    ASSERT_EQ(2u, xyvec.size());
    
    ASSERT_EQ(xyvec[0], CThreeIndexes(0, 1, 0));
    
    ASSERT_EQ(xyvec[1], CThreeIndexes(7, 2, 0));
}

TEST_F(CGenFuncTest, GetTriplesFromQuadrupleIndexes)
{
    CVecFourIndexes vec{{0, 0, 1, 2}, {2, 4, 3, 2}, {1, 7, 0, 2}, {0, 8, -2, 4}};
    
    auto xyzvec = genfunc::getTriplesFromQuadrupleIndexes(vec);
    
    ASSERT_EQ(2u, xyzvec.size());
    
    ASSERT_EQ(xyzvec[0], CThreeIndexes(0, 1, 2));
    
    ASSERT_EQ(xyzvec[1], CThreeIndexes(8, -2, 4));
}

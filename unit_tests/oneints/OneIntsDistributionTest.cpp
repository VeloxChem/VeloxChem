//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OneIntsDistributionTest.hpp"

#include "OneIntsDistributor.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "DenseMatrix.hpp"

TEST_F(COneIntsDistributionTest, DefaultConstructor)
{
    COneIntsDistribution dista;
    
    COneIntsDistribution distb(nullptr, 0, 0, dist1e::batch);
    
    ASSERT_EQ(dista, distb);
}

TEST_F(COneIntsDistributionTest, CopyConstructor)
{
    double val = 1.0;
    
    COneIntsDistribution dista(&val, 10, 20, dist1e::rect);
    
    COneIntsDistribution distb(dista);
    
    ASSERT_EQ(dista, distb);
}

TEST_F(COneIntsDistributionTest, CopyAssignment)
{
    double val = 1.0;
    
    COneIntsDistribution dista(&val, 10, 20, dist1e::rect);
    
    COneIntsDistribution distb = dista;
    
    ASSERT_EQ(dista, distb);
}


TEST_F(COneIntsDistributionTest, DistributeWithBatchPattern)
{
    CDenseMatrix spmat(9, 9);
    
    COneIntsDistribution dist(spmat.values(), 9, 9, dist1e::batch);
    
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 1);
    
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
    
    dist.distribute(adat, bgtos, bgtos, true, 0);
    
    dist.distribute(adat, bgtos, bgtos, true, 2);
    
    CMemBlock2D<double> bdat({ 1.0,  2.0,  0.0,
                               0.0,  4.0,  0.0,
                               0.0,  7.0,  8.0,
                               0.0,  0.0, -1.0,
                               0.0,  0.0,  2.0,
                               0.0,  0.0,  0.0,
                               2.0,  2.0,  3.0,
                               0.0, -1.0,  2.0,
                               0.1, -1.2,  0.0},
                             3, 9);
    
    dist.distribute(bdat, bgtos, bgtos, true, 1);
    
    CDenseMatrix bmat({1.0,  2.0,  0.0, -3.0,  4.0,  2.0, 6.0,  7.0,  8.0,
                       1.0,  2.0,  0.0,  0.0,  4.0,  0.0, 0.0,  7.0,  8.0,
                       1.0,  2.0,  0.0, -3.0,  4.0,  2.0, 6.0,  7.0,  8.0,
                       1.0,  5.7, -1.0,  0.0,  0.0,  2.0, 0.0,  0.0,  0.0,
                       0.0,  0.0, -1.0,  0.0,  0.0,  2.0, 0.0,  0.0,  0.0,
                       1.0,  5.7, -1.0,  0.0,  0.0,  2.0, 0.0,  0.0,  0.0,
                       1.0,  2.0,  3.0,  0.0, -1.0,  2.0, 0.1, -1.2,  0.0,
                       2.0,  2.0,  3.0,  0.0, -1.0,  2.0, 0.1, -1.2,  0.0,
                       1.0,  2.0,  3.0,  0.0, -1.0,  2.0, 0.1, -1.2,  0.0},
                      9, 9);
    
    ASSERT_EQ(spmat, bmat);
}

TEST_F(COneIntsDistributionTest, DistributeWithSymMatrixPattern)
{
    CDenseMatrix spmat(14, 14);
    
    COneIntsDistribution dist(spmat.values(), 14, 14, dist1e::symsq);
    
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 0);
    
    CMemBlock2D<double> s0dat({1.1, 1.2, 1.3, 1.4, 1.5}, 5, 1);
    
    dist.distribute(s0dat, bgtos, bgtos, true, 0);
    
    CMemBlock2D<double> s1dat({2.1, 2.2, 2.3, 2.4, 2.5}, 5, 1);
    
    dist.distribute(s1dat, bgtos, bgtos, true, 1);
    
    CMemBlock2D<double> s2dat({3.1, 3.2, 3.3, 3.4, 3.5}, 5, 1);
    
    dist.distribute(s2dat, bgtos, bgtos, true, 2);
    
    CMemBlock2D<double> s3dat({4.1, 4.2, 4.3, 4.4, 4.5}, 5, 1);
    
    dist.distribute(s3dat, bgtos, bgtos, true, 3);
    
    CMemBlock2D<double> s4dat({5.1, 5.2, 5.3, 5.4, 5.5}, 5, 1);
    
    dist.distribute(s4dat, bgtos, bgtos, true, 4);
    
    auto kgtos = CGtoBlock(lih, bas, 1);
    
    CMemBlock2D<double> sp0dat({2.2, 2.4, 2.6,
                                2.3, 2.6, 2.9,
                                2.0, 2.5, 3.0},
                                3, 3);
    
    dist.distribute(sp0dat, bgtos, kgtos, false, 0);
    
    CMemBlock2D<double> sp1dat({3.2, 3.4, 3.6,
                                3.3, 3.6, 3.9,
                                3.0, 3.5, 4.0},
                                3, 3);
    
    dist.distribute(sp1dat, bgtos, kgtos, false, 1);
    
    CMemBlock2D<double> sp2dat({4.2, 4.4, 4.6,
                                4.3, 4.6, 4.9,
                                4.0, 4.5, 5.0},
                                3, 3);
    
    dist.distribute(sp2dat, bgtos, kgtos, false, 2);
    
    CMemBlock2D<double> sp3dat({5.2, 5.4, 5.6,
                                5.3, 5.6, 5.9,
                                5.0, 5.5, 6.0},
                                3, 3);
    
    dist.distribute(sp3dat, bgtos, kgtos, false, 3);
    
    CMemBlock2D<double> sp4dat({6.2, 6.4, 6.6,
                                6.3, 6.6, 6.9,
                                6.0, 6.5, 7.0},
                                3, 3);
    
    dist.distribute(sp4dat, bgtos, kgtos, false, 4);
    
    CMemBlock2D<double> pp0dat({ 1.0,  2.0,  0.0,
                                -3.0,  4.0,  2.0,
                                 6.0,  7.0,  8.0,
                                 1.0,  5.7, -1.0,
                                 0.0,  0.0,  2.0,
                                 0.0,  0.0,  0.0,
                                 1.0,  2.0,  3.0,
                                 0.0, -1.0,  2.0,
                                 0.1, -1.2,  0.0},
                                3, 9);
    
    dist.distribute(pp0dat, kgtos, kgtos, true, 0);
    
    dist.distribute(pp0dat, kgtos, kgtos, true, 2);
    
    CMemBlock2D<double> pp1dat({ 1.0,  2.0,  0.0,
                                 0.0,  4.0,  0.0,
                                 0.0,  7.0,  8.0,
                                 0.0,  0.0, -1.0,
                                 0.0,  0.0,  2.0,
                                 0.0,  0.0,  0.0,
                                 2.0,  2.0,  3.0,
                                 0.0, -1.0,  2.0,
                                 0.1, -1.2,  0.0},
                                3, 9);
    
    dist.distribute(pp1dat, kgtos, kgtos, true, 1);
   
    CDenseMatrix bmat({1.1, 1.2, 1.3, 1.4, 1.5, 2.2, 2.4,  2.6,  2.3,  2.6, 2.9, 2.0,  2.5, 3.0,
                       2.1, 2.2, 2.3, 2.4, 2.5, 3.2, 3.4,  3.6,  3.3,  3.6, 3.9, 3.0,  3.5, 4.0,
                       3.1, 3.2, 3.3, 3.4, 3.5, 4.2, 4.4,  4.6,  4.3,  4.6, 4.9, 4.0,  4.5, 5.0,
                       4.1, 4.2, 4.3, 4.4, 4.5, 5.2, 5.4,  5.6,  5.3,  5.6, 5.9, 5.0,  5.5, 6.0,
                       5.1, 5.2, 5.3, 5.4, 5.5, 6.2, 6.4,  6.6,  6.3,  6.6, 6.9, 6.0,  6.5, 7.0,
                       2.2, 3.2, 4.2, 5.2, 6.2, 1.0, 2.0,  0.0, -3.0,  4.0, 2.0, 6.0,  7.0, 8.0,
                       2.4, 3.4, 4.4, 5.4, 6.4, 1.0, 2.0,  0.0,  0.0,  4.0, 0.0, 0.0,  7.0, 8.0,
                       2.6, 3.6, 4.6, 5.6, 6.6, 1.0, 2.0,  0.0, -3.0,  4.0, 2.0, 6.0,  7.0, 8.0,
                       2.3, 3.3, 4.3, 5.3, 6.3, 1.0, 5.7, -1.0,  0.0,  0.0, 2.0, 0.0,  0.0, 0.0,
                       2.6, 3.6, 4.6, 5.6, 6.6, 0.0, 0.0, -1.0,  0.0,  0.0, 2.0, 0.0,  0.0, 0.0,
                       2.9, 3.9, 4.9, 5.9, 6.9, 1.0, 5.7, -1.0,  0.0,  0.0, 2.0, 0.0,  0.0, 0.0,
                       2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0,  3.0,  0.0, -1.0, 2.0, 0.1, -1.2, 0.0,
                       2.5, 3.5, 4.5, 5.5, 6.5, 2.0, 2.0,  3.0,  0.0, -1.0, 2.0, 0.1, -1.2, 0.0,
                       3.0, 4.0, 5.0, 6.0, 7.0, 1.0, 2.0,  3.0,  0.0, -1.0, 2.0, 0.1, -1.2, 0.0},
                      14, 14);
    
    ASSERT_EQ(spmat, bmat);
}

TEST_F(COneIntsDistributionTest, DistributeWithAntiSymMatrixPattern)
{
    CDenseMatrix spmat(14, 14);
    
    COneIntsDistribution dist(spmat.values(), 14, 14, dist1e::antisq);
    
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 0);
    
    CMemBlock2D<double> s0dat({1.1, 1.2, 1.3, 1.4, 1.5}, 5, 1);
    
    dist.distribute(s0dat, bgtos, bgtos, true, 0);
    
    CMemBlock2D<double> s1dat({2.1, 2.2, 2.3, 2.4, 2.5}, 5, 1);
    
    dist.distribute(s1dat, bgtos, bgtos, true, 1);
    
    CMemBlock2D<double> s2dat({3.1, 3.2, 3.3, 3.4, 3.5}, 5, 1);
    
    dist.distribute(s2dat, bgtos, bgtos, true, 2);
    
    CMemBlock2D<double> s3dat({4.1, 4.2, 4.3, 4.4, 4.5}, 5, 1);
    
    dist.distribute(s3dat, bgtos, bgtos, true, 3);
    
    CMemBlock2D<double> s4dat({5.1, 5.2, 5.3, 5.4, 5.5}, 5, 1);
    
    dist.distribute(s4dat, bgtos, bgtos, true, 4);
    
    auto kgtos = CGtoBlock(lih, bas, 1);
    
    CMemBlock2D<double> sp0dat({2.2, 2.4, 2.6,
                                2.3, 2.6, 2.9,
                                2.0, 2.5, 3.0},
                                3, 3);
    
    dist.distribute(sp0dat, bgtos, kgtos, false, 0);
    
    CMemBlock2D<double> sp1dat({3.2, 3.4, 3.6,
                                3.3, 3.6, 3.9,
                                3.0, 3.5, 4.0},
                                3, 3);
    
    dist.distribute(sp1dat, bgtos, kgtos, false, 1);
    
    CMemBlock2D<double> sp2dat({4.2, 4.4, 4.6,
                                4.3, 4.6, 4.9,
                                4.0, 4.5, 5.0},
                                3, 3);
    
    dist.distribute(sp2dat, bgtos, kgtos, false, 2);
    
    CMemBlock2D<double> sp3dat({5.2, 5.4, 5.6,
                                5.3, 5.6, 5.9,
                                5.0, 5.5, 6.0},
                                3, 3);
    
    dist.distribute(sp3dat, bgtos, kgtos, false, 3);
    
    CMemBlock2D<double> sp4dat({6.2, 6.4, 6.6,
                                6.3, 6.6, 6.9,
                                6.0, 6.5, 7.0},
                                3, 3);
    
    dist.distribute(sp4dat, bgtos, kgtos, false, 4);
    
    CMemBlock2D<double> pp0dat({ 1.0,  2.0,  0.0,
                                -3.0,  4.0,  2.0,
                                 6.0,  7.0,  8.0,
                                 1.0,  5.7, -1.0,
                                 0.0,  0.0,  2.0,
                                 0.0,  0.0,  0.0,
                                 1.0,  2.0,  3.0,
                                 0.0, -1.0,  2.0,
                                 0.1, -1.2,  0.0},
                                3, 9);
    
    dist.distribute(pp0dat, kgtos, kgtos, true, 0);
    
    dist.distribute(pp0dat, kgtos, kgtos, true, 2);
    
    CMemBlock2D<double> pp1dat({ 1.0,  2.0,  0.0,
                                 0.0,  4.0,  0.0,
                                 0.0,  7.0,  8.0,
                                 0.0,  0.0, -1.0,
                                 0.0,  0.0,  2.0,
                                 0.0,  0.0,  0.0,
                                 2.0,  2.0,  3.0,
                                 0.0, -1.0,  2.0,
                                 0.1, -1.2,  0.0},
                                3, 9);
    
    dist.distribute(pp1dat, kgtos, kgtos, true, 1);
   
    CDenseMatrix bmat({1.1, 1.2, 1.3, 1.4, 1.5, -2.2, -2.4,  -2.6,  -2.3,  -2.6, -2.9, -2.0,  -2.5, -3.0,
                       2.1, 2.2, 2.3, 2.4, 2.5, -3.2, -3.4,  -3.6,  -3.3,  -3.6, -3.9, -3.0,  -3.5, -4.0,
                       3.1, 3.2, 3.3, 3.4, 3.5, -4.2, -4.4,  -4.6,  -4.3,  -4.6, -4.9, -4.0,  -4.5, -5.0,
                       4.1, 4.2, 4.3, 4.4, 4.5, -5.2, -5.4,  -5.6,  -5.3,  -5.6, -5.9, -5.0,  -5.5, -6.0,
                       5.1, 5.2, 5.3, 5.4, 5.5, -6.2, -6.4,  -6.6,  -6.3,  -6.6, -6.9, -6.0,  -6.5, -7.0,
                       2.2, 3.2, 4.2, 5.2, 6.2,  1.0,  2.0,   0.0,  -3.0,   4.0,  2.0,  6.0,   7.0,  8.0,
                       2.4, 3.4, 4.4, 5.4, 6.4,  1.0,  2.0,   0.0,   0.0,   4.0,  0.0,  0.0,   7.0,  8.0,
                       2.6, 3.6, 4.6, 5.6, 6.6,  1.0,  2.0,   0.0,  -3.0,   4.0,  2.0,  6.0,   7.0,  8.0,
                       2.3, 3.3, 4.3, 5.3, 6.3,  1.0,  5.7,  -1.0,   0.0,   0.0,  2.0,  0.0,   0.0,  0.0,
                       2.6, 3.6, 4.6, 5.6, 6.6,  0.0,  0.0,  -1.0,   0.0,   0.0,  2.0,  0.0,   0.0,  0.0,
                       2.9, 3.9, 4.9, 5.9, 6.9,  1.0,  5.7,  -1.0,   0.0,   0.0,  2.0,  0.0,   0.0,  0.0,
                       2.0, 3.0, 4.0, 5.0, 6.0,  1.0,  2.0,   3.0,   0.0,  -1.0,  2.0,  0.1,  -1.2,  0.0,
                       2.5, 3.5, 4.5, 5.5, 6.5,  2.0,  2.0,   3.0,   0.0,  -1.0,  2.0,  0.1,  -1.2,  0.0,
                       3.0, 4.0, 5.0, 6.0, 7.0,  1.0,  2.0,   3.0,   0.0,  -1.0,  2.0,  0.1,  -1.2,  0.0},
                      14, 14);
    
    ASSERT_EQ(spmat, bmat);
}

TEST_F(COneIntsDistributionTest, DistributeRectSymMatrixPattern)
{
    CDenseMatrix spmat(5, 14);
    
    spmat.zero();
    
    COneIntsDistribution dist(spmat.values(), 5, 14, dist1e::rect);
    
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(lih, bas, 0);
    
    CGtoBlock kgtos(lih, bas, 1);
    
    CMemBlock2D<double> sp0dat({2.2, 2.4, 2.6,
                                2.3, 2.6, 2.9,
                                2.0, 2.5, 3.0},
                                3, 3);
    
    dist.distribute(sp0dat, bgtos, kgtos, false, 0);
    
    CMemBlock2D<double> sp1dat({3.2, 3.4, 3.6,
                                3.3, 3.6, 3.9,
                                3.0, 3.5, 4.0},
                                3, 3);
    
    dist.distribute(sp1dat, bgtos, kgtos, false, 1);
    
    CMemBlock2D<double> sp2dat({4.2, 4.4, 4.6,
                                4.3, 4.6, 4.9,
                                4.0, 4.5, 5.0},
                                3, 3);
    
    dist.distribute(sp2dat, bgtos, kgtos, false, 2);
    
    CMemBlock2D<double> sp3dat({5.2, 5.4, 5.6,
                                5.3, 5.6, 5.9,
                                5.0, 5.5, 6.0},
                                3, 3);
    
    dist.distribute(sp3dat, bgtos, kgtos, false, 3);
    
    CMemBlock2D<double> sp4dat({6.2, 6.4, 6.6,
                                6.3, 6.6, 6.9,
                                6.0, 6.5, 7.0},
                                3, 3);
    
    dist.distribute(sp4dat, bgtos, kgtos, false, 4);
    
   
    CDenseMatrix bmat({0.0, 0.0, 0.0, 0.0, 0.0, 2.2, 2.4, 2.6, 2.3, 2.6, 2.9, 2.0, 2.5, 3.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 3.2, 3.4, 3.6, 3.3, 3.6, 3.9, 3.0, 3.5, 4.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 4.2, 4.4, 4.6, 4.3, 4.6, 4.9, 4.0, 4.5, 5.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 5.2, 5.4, 5.6, 5.3, 5.6, 5.9, 5.0, 5.5, 6.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 6.2, 6.4, 6.6, 6.3, 6.6, 6.9, 6.0, 6.5, 7.0},
                      5, 14);
    
    ASSERT_EQ(spmat, bmat);
}


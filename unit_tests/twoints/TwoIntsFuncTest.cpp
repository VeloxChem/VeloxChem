//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "TwoIntsFuncTest.hpp"

#include "TwoIntsFunc.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CTwoIntsFuncTest, CompDistancesAQ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> raq(11, 3);

    twointsfunc::compDistancesAQ(raq, pbra, ppairs, 2);
    
    CMemBlock2D<double> refaq({0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               1.2000, 1.2000, 1.2000, 1.2000, 1.2000, 1.2000,
                               1.7400 / 2.250, 0.360 / 1.100,
                               1.2000, 0.0984 / 0.882, 0.0000},  11, 3);
    
    ASSERT_EQ(raq, refaq);
}

TEST_F(CTwoIntsFuncTest, CompFactorsForThreeCenterElectronRepulsion)
{
    CMolecularBasis bas = vlxbas::getReducedTestBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 2);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(7, 10);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 0);
    
    CMemBlock2D<double> reffacts({1.0000 / 9.0000, 1.0000 / 8.0000, 1.0000 / 8.0000,
                                  1.0000 / 7.0000, 1.0000 / 6.8000, 1.0000 / 5.8000,
                                  1.0000 / 4.6000,
                                  18.000 / 9.0000, 15.000 / 8.0000, 15.000 / 8.0000,
                                  12.000 / 7.0000, 11.400 / 6.8000, 8.4000 / 5.8000,
                                  4.8000 / 4.6000,
                                  1.0000 / 3.0000, 1.0000 / 3.0000, 1.0000 / 3.0000,
                                  1.0000 / 3.0000, 1.0000 / 3.0000, 1.0000 / 3.0000,
                                  1.0000 / 3.0000,
                                  6.0000 / 9.0000, 5.0000 / 8.0000, 5.0000 / 8.0000,
                                  4.0000 / 7.0000, 3.8000 / 6.8000, 2.8000 / 5.8000,
                                  1.6000 / 4.6000,
                                  3.0000 / 9.0000, 3.0000 / 8.0000, 3.0000 / 8.0000,
                                  3.0000 / 7.0000, 3.0000 / 6.8000, 3.0000 / 5.8000,
                                  3.0000 / 4.6000,
                                  1.0000 / 8.0000, 1.0000 / 7.0000, 1.0000 / 7.0000,
                                  1.0000 / 6.0000, 1.0000 / 5.8000, 1.0000 / 4.8000,
                                  1.0000 / 3.6000,
                                  12.000 / 8.0000, 10.000 / 7.0000, 10.000 / 7.0000,
                                  8.0000 / 6.0000, 7.6000 / 5.8000, 5.6000 / 4.8000,
                                  3.2000 / 3.6000,
                                  1.0000 / 2.0000, 1.0000 / 2.0000, 1.0000 / 2.0000,
                                  1.0000 / 2.0000, 1.0000 / 2.0000, 1.0000 / 2.0000,
                                  1.0000 / 2.0000,
                                  6.0000 / 8.0000, 5.0000 / 7.0000, 5.0000 / 7.0000,
                                  4.0000 / 6.0000, 3.8000 / 5.8000, 2.8000 / 4.8000,
                                  1.6000 / 3.6000,
                                  2.0000 / 8.0000, 2.0000 / 7.0000, 2.0000 / 7.0000,
                                  2.0000 / 6.0000, 2.0000 / 5.8000, 2.0000 / 4.8000,
                                  2.0000 / 3.6000}, 7, 10);
    
    ASSERT_EQ(facts, reffacts);
}

TEST_F(CTwoIntsFuncTest, CompCoordinatesForW)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(11, 10);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 0);
    
    CMemBlock2D<double> rw(11, 6);
    
    twointsfunc::compCoordinatesForW(rw, facts, 5, pbra, ppairs, 0);
    
    CMemBlock2D<double> refrw({0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.000 / 4.350, 0.000 / 3.200, 0.000 / 3.200, 0.000 / 2.050,
                               0.000 / 2.982, 0.000 / 1.832, 0.960 / 3.700, 0.960 / 2.550,
                               0.000 / 1.614, 0.960 / 2.332, 1.920 / 3.050,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.000 / 3.200, 0.000 / 2.050, 0.000 / 2.050, 0.000 / 0.900,
                               0.000 / 1.832, 0.000 / 0.682, 0.960 / 2.550, 0.960 / 1.400,
                               0.000 / 0.464, 0.960 / 1.182, 1.920 / 1.900}, 11, 6);
    
    ASSERT_EQ(rw, refrw);
}

TEST_F(CTwoIntsFuncTest, CompDistancesWA)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(11, 5);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 2);
    
    CMemBlock2D<double> rw(11, 3);
    
    twointsfunc::compCoordinatesForW(rw, facts, 5, pbra, ppairs, 2);
    
    CMemBlock2D<double> rwa(11, 3);
    
    twointsfunc::compDistancesWA(rwa, rw, pbra, ppairs, 2);
    
    CMemBlock2D<double> refrwa({0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,
                               -3.4800 / 3.700, -2.1000 / 2.550, -2.1000 / 2.550, -0.7200 / 1.400,
                               -1.8384 / 2.332, -0.4584 / 1.182, -1.7400 / 3.050, -0.3600 / 1.900,
                               -0.1968 / 0.964, -0.0984 / 1.682,  0.0000 / 2.400}, 11, 3);
    
    ASSERT_EQ(rwa, refrwa);
}

TEST_F(CTwoIntsFuncTest, CompDistancesForWD)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(11, 10);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 0);
    
    CMemBlock2D<double> rw(11, 6);
    
    twointsfunc::compCoordinatesForW(rw, facts, 5, pbra, ppairs, 0);
    
    CMemBlock2D<double> rwd(11, 6);
    
    twointsfunc::compDistancesWD(rwd, rw, pbra, ppairs, 0);
    
    CMemBlock2D<double> refrwd({0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000,
                                0.000 / 4.350,  0.0000 / 3.200,  0.000 / 3.200,  0.000 / 2.050,
                                0.000 / 2.982,  0.0000 / 1.832, -3.480 / 3.700, -2.100 / 2.550,
                                0.000 / 1.614, -1.8384 / 2.332, -1.740 / 3.050,
                                0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                0.0000, 0.0000, 0.0000, 0.0000,
                                0.000 / 3.200,  0.000  / 2.050,  0.000 / 2.050,  0.000 / 0.900,
                                0.000 / 1.832,  0.000  / 0.682, -2.100 / 2.550, -0.720 / 1.400,
                                0.000 / 0.464, -0.4584 / 1.182, -0.360 / 1.900}, 11, 6);
    
    ASSERT_EQ(rwd, refrwd);
}

TEST_F(CTwoIntsFuncTest, CompDistancesForPQ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> t0rpq(4, 12);
    
    twointsfunc::compDistancesPQ(t0rpq, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> ref0rpq(4, 12);
    
    ref0rpq.zero();
    
    ASSERT_EQ(t0rpq, ref0rpq);
    
    CMemBlock2D<double> t2rpq(8, 6);
    
    twointsfunc::compDistancesPQ(t2rpq, ppairs, ppairs, 8, 2);
    
    CMemBlock2D<double> ref2rpq({0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000,
                                 0.960 / 2.250, 0.960 / 2.250, 0.960 / 2.250,
                                 0.960 / 2.250, 0.960 / 2.250, 0.960 / 2.250,
                                 0.000, -1.104 / 2.475,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000,
                                 0.960 / 1.100, 0.960 / 1.100, 0.960 / 1.100,
                                 0.960 / 1.100, 0.960 / 1.100, 0.960 / 1.100,
                                 1.104 / 2.475, 0.000},
                                8, 6);
    
    ASSERT_EQ(t2rpq, ref2rpq); 
    
    CMemBlock2D<double> t4rpq(11, 3);
    
    twointsfunc::compDistancesPQ(t4rpq, ppairs, ppairs, 11, 4);
    
    CMemBlock2D<double> ref4rpq({0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000, 0.000,
                                 0.960 / 0.882, 0.960 / 0.882, 0.960 / 0.882,
                                 0.960 / 0.882, 0.960 / 0.882, 0.960 / 0.882,
                                 1.31328 / 1.9845, 0.20928 / 0.9702,
                                 0.960 / 0.882, 0.000, -0.15744 / 1.4112},
                                11, 3);
    
    ASSERT_EQ(t4rpq, ref4rpq);
}

TEST_F(CTwoIntsFuncTest, CompFactorsForElectronRepulsion)
{
    CMolecularBasis bas = vlxbas::getReducedTestBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 2);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> r0facts(4, 16);
    
    twointsfunc::compFactorsForElectronRepulsion(r0facts, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> ref0facts({ 1.0 / 12.0,  1.0 / 11.0,  1.0 / 11.0,  1.0 / 10.0,
                                   36.0 / 12.0, 30.0 / 11.0, 30.0 / 11.0, 24.0 / 10.0,
                                    6.0 / 12.0,  5.0 / 11.0,  5.0 / 11.0,  4.0 / 10.0,
                                    6.0 / 12.0,  6.0 / 11.0,  6.0 / 11.0,  6.0 / 10.0,
                                    1.0 / 11.0,  1.0 / 10.0,  1.0 / 10.0,  1.0 /  9.0,
                                   30.0 / 11.0, 25.0 / 10.0, 25.0 / 10.0, 20.0 /  9.0,
                                    6.0 / 11.0,  5.0 / 10.0,  5.0 / 10.0,  4.0 /  9.0,
                                    5.0 / 11.0,  5.0 / 10.0,  5.0 / 10.0,  5.0 /  9.0,
                                    1.0 / 11.0,  1.0 / 10.0,  1.0 / 10.0,  1.0 /  9.0,
                                   30.0 / 11.0, 25.0 / 10.0, 25.0 / 10.0, 20.0 /  9.0,
                                    6.0 / 11.0,  5.0 / 10.0,  5.0 / 10.0,  4.0 /  9.0,
                                    5.0 / 11.0,  5.0 / 10.0,  5.0 / 10.0,  5.0 /  9.0,
                                    1.0 / 10.0,  1.0 /  9.0,  1.0 /  9.0,  1.0 /  8.0,
                                   24.0 / 10.0, 20.0 /  9.0, 20.0 /  9.0, 16.0 /  8.0,
                                    6.0 / 10.0,  5.0 /  9.0,  5.0 /  9.0,  4.0 /  8.0,
                                    4.0 / 10.0,  4.0 /  9.0,  4.0 /  9.0,  4.0 /  8.0},
                                  4, 16);
    
    ASSERT_EQ(r0facts, ref0facts);
    
    CMemBlock2D<double> r2facts(7, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r2facts, ppairs, ppairs, 7, 2);
    
    CMemBlock2D<double> ref2facts({ 1.00 / 7.60, 1.00 / 6.60, 1.00 / 6.60, 1.00 / 5.60,
                                    1.00 / 5.40, 1.00 / 4.40, 1.00 / 3.20,
                                    9.60 / 7.60, 8.00 / 6.60, 8.00 / 6.60, 6.40 / 5.60,
                                    6.08 / 5.40, 4.48 / 4.40, 2.56 / 3.20,
                                    6.00 / 7.60, 5.00 / 6.60, 5.00 / 6.60, 4.00 / 5.60,
                                    3.80 / 5.40, 2.80 / 4.40, 1.60 / 3.20,
                                    1.60 / 7.60, 1.60 / 6.60, 1.60 / 6.60, 1.60 / 5.60,
                                    1.60 / 5.40, 1.60 / 4.40, 1.60 / 3.20},
                                  7, 4);
    
    ASSERT_EQ(r2facts, ref2facts);
    
    CMemBlock2D<double> r3facts(7, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r3facts, ppairs, ppairs, 7, 2);
    
    ASSERT_EQ(r3facts, ref2facts);
}

TEST_F(CTwoIntsFuncTest, CompCoordinatesForWUsingPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> r0facts(4, 16);
    
    twointsfunc::compFactorsForElectronRepulsion(r0facts, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> r0w(4, 12);
    
    twointsfunc::compCoordinatesForW(r0w, r0facts, 4, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> refr0w(4, 12);
    
    refr0w.zero();
    
    ASSERT_EQ(r0w, refr0w);
    
    CMemBlock2D<double> r5facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r5facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r5w(11, 3);
    
    twointsfunc::compCoordinatesForW(r5w, r5facts, 4, ppairs, ppairs, 11, 5);
    
    
    CMemBlock2D<double> refr5w({0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                1.920 / 4.500, 1.920 / 3.350, 1.920 / 3.350, 1.920 / 2.200,
                                1.920 / 3.132, 1.920 / 1.982, 2.880 / 3.850, 2.880 / 2.700,
                                1.920 / 1.764, 2.880 / 2.482, 3.840 / 3.200},
                               11, 3);
    
    ASSERT_EQ(r5w, refr5w);
    
    CMemBlock2D<double> r6facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r6facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r6w(11, 3);
    
    twointsfunc::compCoordinatesForW(r6w, r6facts, 4, ppairs, ppairs, 11, 5);
    
    ASSERT_EQ(r6w, refr5w);
}

TEST_F(CTwoIntsFuncTest, CompCoordinatesForWP)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> r0facts(4, 16);
    
    twointsfunc::compFactorsForElectronRepulsion(r0facts, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> r0w(4, 12);
    
    twointsfunc::compCoordinatesForW(r0w, r0facts, 4, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> r0wp(4, 12);
    
    twointsfunc::compDistancesWP(r0wp, r0w, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> refr0wp(4, 12);
    
    refr0wp.zero();
    
    ASSERT_EQ(r0wp, refr0wp);
    
    CMemBlock2D<double> r5facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r5facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r5w(11, 3);
    
    twointsfunc::compCoordinatesForW(r5w, r5facts, 4, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r5wp(11, 3);
    
    twointsfunc::compDistancesWP(r5wp, r5w, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> refr5wp({0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                -3.4800 / 4.500, -2.1000 / 3.350, -2.100 / 3.350, -0.720 / 2.200,
                                -1.8384 / 3.132, -0.4584 / 1.982, -1.740 / 3.850, -0.360 / 2.700,
                                -0.1968 / 1.764, -0.0984 / 2.482, 0.000},
                               11, 3);
    
    ASSERT_EQ(r5wp, refr5wp);
    
    CMemBlock2D<double> r6facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r6facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r6w(11, 3);
    
    twointsfunc::compCoordinatesForW(r6w, r6facts, 4, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r6wp(11, 3);
    
    twointsfunc::compDistancesWP(r6wp, r6w, ppairs, ppairs, 11, 5);
    
    ASSERT_EQ(r6wp, refr5wp);
}

TEST_F(CTwoIntsFuncTest, CompCoordinatesForWQUsingPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> r0facts(4, 16);
    
    twointsfunc::compFactorsForElectronRepulsion(r0facts, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> r0w(4, 12);
    
    twointsfunc::compCoordinatesForW(r0w, r0facts, 4, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> r0wq(4, 12);
    
    twointsfunc::compDistancesWQ(r0wq, r0w, ppairs, ppairs, 4, 0);
    
    CMemBlock2D<double> refr0wq(4, 12);
    
    refr0wq.zero();
    
    ASSERT_EQ(r0wq, refr0wq);
    
    CMemBlock2D<double> r5facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r5facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r5w(11, 3);
    
    twointsfunc::compCoordinatesForW(r5w, r5facts, 4, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r5wq(11, 3);
    
    twointsfunc::compDistancesWQ(r5wq, r5w, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> refr5wq({0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                 0.000, 0.000, 0.000,
                                 1.920 / 4.500, 1.92000 / 3.350000, 1.9200 / 3.3500, 1.920 / 2.200,
                                 1.920 / 3.132, 1.92000 / 1.982000, 2.7840 / 8.6625, 0.576 / 2.970,
                                 1.920 / 1.764, 0.15744 / 2.189124, 0.000},
                                11, 3);
        
    ASSERT_EQ(r5wq, refr5wq);
    
    CMemBlock2D<double> r6facts(11, 4);
    
    twointsfunc::compFactorsForElectronRepulsion(r6facts, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r6w(11, 3);
    
    twointsfunc::compCoordinatesForW(r6w, r6facts, 4, ppairs, ppairs, 11, 5);
    
    CMemBlock2D<double> r6wq(11, 3);
    
    twointsfunc::compDistancesWQ(r6wq, r6w, ppairs, ppairs, 11, 5);
    
    ASSERT_EQ(r6wq, refr5wq);
}

TEST_F(CTwoIntsFuncTest, CompEffectiveDistancesForPQ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock<double> t0rpq(6);
    
    t0rpq.zero();
   
    twointsfunc::compEffectiveDistancesPQ(t0rpq, ppairs, ppairs, true, 0);
    
    CMemBlock<double> ref00pq({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    
    ASSERT_EQ(t0rpq, ref00pq);
    
    twointsfunc::compEffectiveDistancesPQ(t0rpq, ppairs, ppairs, false, 0);
    
    CMemBlock<double> ref01pq({0.0000, 0.0000, 0.781076809730870, 0.0000,
                               1.08843537414966, 1.2000});
    
    ASSERT_EQ(t0rpq, ref01pq);
    
    CMemBlock<double> t2rpq(6);
    
    t2rpq.zero();
    
    twointsfunc::compEffectiveDistancesPQ(t2rpq, ppairs, ppairs, true, 2);
    
    CMemBlock<double> ref20pq({0.781076809730870, 0.781076809730870, 0.0,
                               0.0, 0.0, 0.0});
    
    ASSERT_EQ(t2rpq, ref20pq);
    
    twointsfunc::compEffectiveDistancesPQ(t2rpq, ppairs, ppairs, false, 2);
    
    CMemBlock<double> ref21pq({0.781076809730870, 0.781076809730870, 0.0,
                               0.781076809730870, .307358564418790,
                               0.418923190269130});
    
    ASSERT_EQ(t2rpq, ref21pq);
}

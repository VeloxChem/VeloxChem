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

#include "DiagEriRecFactsTest.hpp"

#include "DiagEriRecFacts.hpp"
#include "BinnedGtoBlock.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CDiagEriRecFactsTest, CompHostDistancesPQ)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> s21pairs(s2gtos, s1gtos);
    
    EXPECT_EQ(s21pairs.getNumberOfPrimPairs(), 2);
    
    EXPECT_EQ(s21pairs.getNumberOfContrPairs(), 4);
    
    BufferHostMY<double, 3> rpq(4);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 0, 0);
    
    const auto rpq_00 = BufferHostMY<double, 3>::Zero(4);
    
    EXPECT_EQ(rpq_00, rpq);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 1, 1);
    
    EXPECT_EQ(rpq_00, rpq);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 0, 1);
    
    const BufferHostMY<double, 3> rpq_01(4, { 0.00000000000000000000,
                                             -0.05024154589371980677,
                                              0.03571428571428571428,
                                              0.00000000000000000000,
                                              0.00000000000000000000,
                                             -0.07536231884057971014,
                                              0.05357142857142857142,
                                              0.00000000000000000000,
                                              0.00000000000000000000,
                                             -0.13816425120772946860,
                                              0.09821428571428571428,
                                              0.00000000000000000000,});
    
    EXPECT_EQ(rpq_01, rpq);
}

TEST_F(CDiagEriRecFactsTest, CompHostFactorRho)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> s21pairs(s2gtos, s1gtos);
    
    EXPECT_EQ(s21pairs.getNumberOfPrimPairs(), 2);
    
    EXPECT_EQ(s21pairs.getNumberOfContrPairs(), 4);
    
    BufferHostX<double> frho(4);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 0, 0);
    
    const auto frho_00 = BufferHostX<double>(4, {2.00000000000000000000,
                                                 1.80000000000000000000,
                                                 2.10000000000000000000,
                                                 1.90000000000000000000});

    EXPECT_EQ(frho_00, frho);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 0, 1);
    
    const auto frho_01 = BufferHostX<double>(4, {1.61194029850746268656,
                                                 1.40338983050847457627,
                                                 1.81621621621621621621,
                                                 1.61212121212121212121});

    EXPECT_EQ(frho_01, frho);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 1, 1);
    
    const auto frho_11 = BufferHostX<double>(4, {1.35000000000000000000,
                                                 1.15000000000000000000,
                                                 1.60000000000000000000,
                                                 1.40000000000000000000});
    
    EXPECT_EQ(frho_11, frho);
}

TEST_F(CDiagEriRecFactsTest, CompHostFactorNorm)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> s21pairs(s2gtos, s1gtos);
    
    EXPECT_EQ(s21pairs.getNumberOfPrimPairs(), 2);
    
    EXPECT_EQ(s21pairs.getNumberOfContrPairs(), 4);
    
    BufferHostX<double> fnorm(4);
    
    derirec::compHostFactorNorm(fnorm.data(), &s21pairs, 0, 4, 0, 0);
    
    const auto fnorm_00 = BufferHostX<double>(4, {0.19938331312932041708,
                                                  0.01729165758259379585,
                                                  0.00964583723283758917,
                                                  0.13759209114032758783});

    EXPECT_EQ(fnorm_00, fnorm);
    
    derirec::compHostFactorNorm(fnorm.data(), &s21pairs, 0, 4, 0, 1);
    
    const auto fnorm_01 = BufferHostX<double>(4, {0.45025968744901134950,
                                                  0.05045735105393916422,
                                                  0.02576024760043138316,
                                                  0.32099138041350482142});

    EXPECT_EQ(fnorm_01, fnorm);
    
    derirec::compHostFactorNorm(fnorm.data(), &s21pairs, 0, 4, 1, 1);
    
    const auto fnorm_11 = BufferHostX<double>(4, {0.25420104491165196955,
                                                  0.03680885223437551286,
                                                  0.01719887917495777140,
                                                  0.18721182563226586746});
    
    EXPECT_EQ(fnorm_11, fnorm);
}

TEST_F(CDiagEriRecFactsTest, CompHostBoysArguments)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> s21pairs(s2gtos, s1gtos);
    
    EXPECT_EQ(s21pairs.getNumberOfPrimPairs(), 2);
    
    EXPECT_EQ(s21pairs.getNumberOfContrPairs(), 4);
    
    BufferHostMY<double, 3> rpq(4);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 0, 0);
    
    BufferHostX<double> frho(4);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 0, 0);
    
    BufferHostX<double> fargs(4);
    
    derirec::compHostBoysArguments(fargs, rpq, frho.data(), 4);
    
    const auto fargs_00 = BufferHostX<double>::Zero(4);
    
    EXPECT_EQ(fargs_00, fargs);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 1, 1);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 1, 1);
    
    derirec::compHostBoysArguments(fargs, rpq, frho.data(), 4);
    
    EXPECT_EQ(fargs_00, fargs);
    
    derirec::compHostDistancesPQ(rpq, &s21pairs, 0, 4, 0, 1);
    
    derirec::compHostFactorRho(frho.data(), &s21pairs, 0, 4, 0, 1);
    
    derirec::compHostBoysArguments(fargs, rpq, frho.data(), 4);
    
    const auto fargs_01 = BufferHostX<double>(4, {0.00000000000000000000,
                                                  0.03830279210677147299,
                                                  0.02504826254826254822,
                                                  0.00000000000000000000});
    
    EXPECT_EQ(fargs_01, fargs);
}

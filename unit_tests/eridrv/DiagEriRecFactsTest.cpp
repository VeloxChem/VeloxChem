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

TEST_F(CDiagEriRecFactsTest, CompDistancesPQ)
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

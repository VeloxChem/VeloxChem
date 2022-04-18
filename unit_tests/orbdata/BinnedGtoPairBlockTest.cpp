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

#include "BinnedGtoPairBlockTest.hpp"

#include "BinnedGtoPairBlock.hpp"
#include "BinnedGtoBlock.hpp"
#include "CheckFunctions.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"


TEST_F(CBinnedGtoPairBlockTest, ConstructorWithTwoBinnedGtoBlock)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    const CBinnedGtoPairBlock<double> apairs(gto, gto);
    
    const BufferHostMY<int32_t, 2> atmidx(4, {0, 0, 1, 1, 0, 1, 0, 1});
    
    const BufferHostXY<int32_t> angidx(6, 4,  {6, 6, 7, 7, 9,  9, 10, 10, 12, 12, 13, 13,
                                               6, 7, 6, 7, 9, 10,  9, 10, 12, 13, 12, 13});
    
    const BufferHostMY<double, 24> ppdata(4, {0.16400000000000000000,  0.88200000000000000000,  0.88200000000000000000, 1.6000000000000000000,
                                              6.09756097560975609756,  1.13378684807256235827,  1.13378684807256235827, 0.6250000000000000000,
                                              0.04100000000000000000,  0.07437641723356009070,  0.07437641723356009070, 0.4000000000000000000,
                                              83.8414994866341000000,  6.03959902057394000000,  6.03959902057394000000, 2.7513436295111000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.08843537414965986394,  1.08843537414965986394, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.08843537414965986394, -0.11156462585034013605, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000, -0.11156462585034013605,  1.08843537414965986394, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  1.20000000000000000000, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.20000000000000000000,  0.00000000000000000000, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000, -1.20000000000000000000,  1.20000000000000000000, 0.0000000000000000000,
                                              0.08200000000000000000,  0.08200000000000000000,  0.80000000000000000000, 0.8000000000000000000,
                                              0.08200000000000000000,  0.80000000000000000000,  0.08200000000000000000, 0.8000000000000000000});
    
    CBinnedGtoPairBlock<double> bpairs(1, 1, 1, angidx, atmidx, ppdata);

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, ConstructorWithSingleBinnedGtoBlock)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock apairs(gto);
    
    const BufferHostMY<int32_t, 2> atmidx(3, {0, 0, 1, 0, 1, 1});
    
    const BufferHostXY<int32_t> angidx(6, 3, {6, 6, 7, 9,  9, 10, 12, 12, 13,
                                              6, 7, 7, 9, 10, 10, 12, 13, 13});
    
    const BufferHostMY<double, 24> ppdata(3, {0.16400000000000000000,  0.88200000000000000000, 1.6000000000000000000,
                                              6.09756097560975609756,  1.13378684807256235827, 0.6250000000000000000,
                                              0.04100000000000000000,  0.07437641723356009070, 0.4000000000000000000,
                                              83.8414994866341000000,  6.03959902057394000000, 2.7513436295111000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.08843537414965986394, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.08843537414965986394, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000, -0.11156462585034013605, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  1.20000000000000000000, 1.2000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000,  0.00000000000000000000, 0.0000000000000000000,
                                              0.00000000000000000000, -1.20000000000000000000, 0.0000000000000000000,
                                              0.08200000000000000000,  0.08200000000000000000, 0.8000000000000000000,
                                              0.08200000000000000000,  0.80000000000000000000, 0.8000000000000000000});
    
    CBinnedGtoPairBlock<double> bpairs(1, 1, 1, angidx, atmidx, ppdata);

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    CBinnedGtoBlock<double> bgto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);

    CBinnedGtoPairBlock<double> bpairs(apairs);

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    CBinnedGtoBlock<double> bgto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);

    CBinnedGtoPairBlock<double> bpairs(CBinnedGtoPairBlock(agto, bgto));

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    CBinnedGtoBlock<double> bgto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);

    CBinnedGtoPairBlock<double> bpairs = apairs;

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    CBinnedGtoBlock<double> bgto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);

    CBinnedGtoPairBlock<double> bpairs = CBinnedGtoPairBlock<double>(agto, bgto);

    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CBinnedGtoPairBlockTest, Compress)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto, gto);
    
    CBinnedGtoPairBlock<double> bpairs(gto);
    
    const BufferHostX<int32_t> idx(4, {2, 2, 3, 2});
    
    ASSERT_EQ(bpairs, apairs.compress(idx, 2));
}

TEST_F(CBinnedGtoPairBlockTest, GetBraAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 1, 2);
    
    CBinnedGtoBlock<double> bgto(lih, bas, 0, 5);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);
    
    ASSERT_EQ(apairs.getBraAngularMomentum(), 1);
}

TEST_F(CBinnedGtoPairBlockTest, GetKetAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 1, 2);
    
    CBinnedGtoBlock<double> bgto(lih, bas, 0, 5);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);
    
    ASSERT_EQ(apairs.getKetAngularMomentum(), 0);
}

TEST_F(CBinnedGtoPairBlockTest, GetNumberOfPrimPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 1, 2);
    
    CBinnedGtoBlock<double> bgto(lih, bas, 0, 5);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);
    
    ASSERT_EQ(apairs.getNumberOfPrimPairs(), 10);
}

TEST_F(CBinnedGtoPairBlockTest, GetNumberOfContrPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 1, 1);
    
    CBinnedGtoBlock<double> bgto(lih, bas, 0, 1);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);
    
    ASSERT_EQ(apairs.getNumberOfContrPairs(), 6);
}

TEST_F(CBinnedGtoPairBlockTest, GetNumberOfComponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> agto(lih, bas, 1, 2);
    
    CBinnedGtoBlock<double> bgto(lih, bas, 0, 5);

    CBinnedGtoPairBlock<double> apairs(agto, bgto);
    
    ASSERT_EQ(apairs.getNumberOfComponents(), 3);
}

TEST_F(CBinnedGtoPairBlockTest, GetFactorsXi)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.164, 0.882, 1.600}, apairs.getFactorsXi());
}

TEST_F(CBinnedGtoPairBlockTest, GetFactorsOneOverXi)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({6.09756097560975609756,  1.13378684807256235827, 0.6250000000000000000},
                     apairs.getFactorsOneOverXi());
}

TEST_F(CBinnedGtoPairBlockTest, GetFactorsZeta)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.04100000000000000000,  0.07437641723356009070, 0.4000000000000000000},
                     apairs.getFactorsZeta());
}

TEST_F(CBinnedGtoPairBlockTest, GetOverlaps)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({83.8414994866341000000,  6.03959902057394000000, 2.7513436295111000000},
                     apairs.getOverlaps());
}

TEST_F(CBinnedGtoPairBlockTest, GetCoordinatesPX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesPX());
}

TEST_F(CBinnedGtoPairBlockTest, GetCoordinatesPY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesPY());
}

TEST_F(CBinnedGtoPairBlockTest, GetCoordinatesPZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.00000000000000000000,  1.08843537414965986394, 1.2000000000000000000},
                     apairs.getCoordinatesPZ());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPAX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesPAX());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPAY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesPAY());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPAZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.00000000000000000000,  1.08843537414965986394, 0.0000000000000000000},
                     apairs.getDistancesPAZ());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPBX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesPBX());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPBY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesPBY());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesPBZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.00000000000000000000, -0.11156462585034013605, 0.0000000000000000000},
                     apairs.getDistancesPBZ());
}
                             
TEST_F(CBinnedGtoPairBlockTest, getCoordinatesAX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesAX());
}

TEST_F(CBinnedGtoPairBlockTest, getCoordinatesAY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesAY());
}

TEST_F(CBinnedGtoPairBlockTest, getCoordinatesAZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 1.2}, apairs.getCoordinatesAZ());
}

TEST_F(CBinnedGtoPairBlockTest, getCoordinatesBX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesBX());
}

TEST_F(CBinnedGtoPairBlockTest, getCoordinatesBY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getCoordinatesBY());
}

TEST_F(CBinnedGtoPairBlockTest, getCoordinatesBZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 1.2, 1.2}, apairs.getCoordinatesBZ());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesABX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesABX());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesABY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, 0.0, 0.0}, apairs.getDistancesABY());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesABZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.0, -1.2, 0.0}, apairs.getDistancesABZ());
}

TEST_F(CBinnedGtoPairBlockTest, GetBraExponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.082, 0.082, 0.8}, apairs.getBraExponents());
}

TEST_F(CBinnedGtoPairBlockTest, GetKetExponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0.082, 0.8, 0.8}, apairs.getKetExponents());
}

TEST_F(CBinnedGtoPairBlockTest, GetBraIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({6, 6, 7}, apairs.getBraIdentifiers(0));
    
    vlxtest::compare({9, 9, 10}, apairs.getBraIdentifiers(1));
    
    vlxtest::compare({12, 12, 13}, apairs.getBraIdentifiers(2));
}

TEST_F(CBinnedGtoPairBlockTest, GetKetIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({6, 7, 7}, apairs.getKetIdentifiers(0));
    
    vlxtest::compare({9, 10, 10}, apairs.getKetIdentifiers(1));
    
    vlxtest::compare({12, 13, 13}, apairs.getKetIdentifiers(2));
}

TEST_F(CBinnedGtoPairBlockTest, GetBraAtomicIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0, 0, 1}, apairs.getBraAtomicIdentifiers());
}

TEST_F(CBinnedGtoPairBlockTest, GetKetAtomicIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    vlxtest::compare({0, 1, 1}, apairs.getKetAtomicIdentifiers());
}

TEST_F(CBinnedGtoPairBlockTest, GetDistancesAB)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> gto(lih, bas, 1, 1);

    CBinnedGtoPairBlock<double> apairs(gto);
    
    auto r01rab = BufferHostMY<double, 3>(2, {0.0,  0.0,
                                              0.0,  0.0,
                                              0.0, -1.2});
    
    auto r12rab = BufferHostMY<double, 3>(2, {0.0,  0.0,
                                              0.0,  0.0,
                                             -1.2, 0.0});
    
    auto r02rab = BufferHostMY<double, 3>(3, {0.0,  0.0, 0.0,
                                              0.0,  0.0, 0.0,
                                              0.0, -1.2, 0.0});
    
    ASSERT_EQ(r01rab, apairs.getDistancesAB(0, 2));
    
    ASSERT_EQ(r12rab, apairs.getDistancesAB(1, 3));
    
    ASSERT_EQ(r02rab, apairs.getDistancesAB(0, 3));
}

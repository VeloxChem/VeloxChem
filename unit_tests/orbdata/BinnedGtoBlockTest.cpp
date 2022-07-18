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

#include "BinnedGtoBlockTest.hpp"

#include <vector>

#include "BinnedGtoBlock.hpp"
#include "CheckFunctions.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CBinnedGtoBlockTest, ConstructorWithMoleculeForSP5)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 5);

    const BufferHostXY<int32_t> sidx(std::vector<int32_t>({0}), 1, 1);

    const BufferHostX<int32_t> aidx(std::vector<int32_t>({0}));

    const BufferHostMY<double, 5> sprim(
        {2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00, 2.450300905100e+00, 7.220957185500e-01, 6.492015032500e-03, 4.774786321500e-02,
         2.026879611100e-01, 4.860657481700e-01, 4.362697795500e-01, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
         0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
        5);

    const CBinnedGtoBlock<double> sdat(0, 5, sidx, aidx, sprim);

    ASSERT_EQ(sorb, sdat);
}

TEST_F(CBinnedGtoBlockTest, ConstructorWithMoleculeForSP3)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 3);

    const BufferHostXY<int32_t> sidx(std::vector<int32_t>({3}), 1, 1);

    const BufferHostX<int32_t> aidx(std::vector<int32_t>({1}));

    const BufferHostMY<double, 5> sprim({1.301070100000e+01,
                                         1.962257200000e+00,
                                         4.445379600000e-01,
                                         1.968215800000e-02,
                                         1.379652400000e-01,
                                         4.783193500000e-01,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         1.200000000000e+00,
                                         1.200000000000e+00,
                                         1.200000000000e+00},
                                        3);

    const CBinnedGtoBlock<double> sdat(0, 3, sidx, aidx, sprim);

    ASSERT_EQ(sorb, sdat);
}

TEST_F(CBinnedGtoBlockTest, ConstructorWithMoleculeForSP1)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    const BufferHostXY<int32_t> sidx(
        {
            1,
            2,
            4,
        },
        1,
        3);

    const BufferHostX<int32_t> aidx(
        {
            0,
            0,
            1,
        },
        3);

    const BufferHostMY<double, 5> sprim({5.281088472100e-02,
                                         2.096094879800e-02,
                                         1.219496200000e-01,
                                         1.000000000000e+00,
                                         1.000000000000e+00,
                                         1.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         1.200000000000e+00},
                                        3);

    const CBinnedGtoBlock<double> sdat(0, 1, sidx, aidx, sprim);

    ASSERT_EQ(sorb, sdat);
}

TEST_F(CBinnedGtoBlockTest, ConstructorWithMoleculeForPP2)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> porb(lih, bas, 1, 2);

    const BufferHostXY<int32_t> pidx(
        {
            5,
            8,
            11,
        },
        3,
        1);

    const BufferHostX<int32_t> aidx(std::vector<int32_t>({0}));

    const BufferHostMY<double, 5> pprim({1.450000000000e+00,
                                         3.000000000000e-01,
                                         2.586000000000e-01,
                                         1.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00},
                                        2);

    const CBinnedGtoBlock<double> pdat(1, 2, pidx, aidx, pprim);

    ASSERT_EQ(porb, pdat);
}

TEST_F(CBinnedGtoBlockTest, ConstructorWithAtomListForSP3)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb11(lih, bas, 1, 1, 0, 3);

    const BufferHostXY<int32_t> sidx(

        std::vector<int32_t>({
            3,
        }),
        1,
        1);

    const BufferHostX<int32_t> aidx(std::vector<int32_t>({1}));

    const BufferHostMY<double, 5> sprim({1.301070100000e+01,
                                         1.962257200000e+00,
                                         4.445379600000e-01,
                                         1.968215800000e-02,
                                         1.379652400000e-01,
                                         4.783193500000e-01,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         0.000000000000e+00,
                                         1.200000000000e+00,
                                         1.200000000000e+00,
                                         1.200000000000e+00},
                                        3);

    const CBinnedGtoBlock<double> sdat(0, 3, sidx, aidx, sprim);

    ASSERT_EQ(sorb11, sdat);

    const CBinnedGtoBlock<double> sorb02(lih, bas, 0, 2, 0, 3);

    ASSERT_EQ(sorb02, sdat);

    const CBinnedGtoBlock<double> sorb01(lih, bas, 0, 1, 0, 3);

    ASSERT_EQ(sorb01, CBinnedGtoBlock<double>());
}

TEST_F(CBinnedGtoBlockTest, ConstructorWithAtomListForPP1)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> porb02(lih, bas, 0, 2, 1, 1);

    const BufferHostXY<int32_t> pidx02({6, 7, 9, 10, 12, 13}, 3, 2);

    const BufferHostX<int32_t> aidx02({0, 1}, 2);

    const BufferHostMY<double, 5> pprim02({8.200000000000e-02,
                                           8.000000000000e-01,
                                           1.000000000000e+00,
                                           1.000000000000e+00,
                                           0.000000000000e+00,
                                           0.000000000000e+00,
                                           0.000000000000e+00,
                                           0.000000000000e+00,
                                           0.000000000000e+00,
                                           1.200000000000e+00},
                                          2);

    CBinnedGtoBlock<double> pdat02(1, 1, pidx02, aidx02, pprim02);

    ASSERT_EQ(porb02, pdat02);

    CBinnedGtoBlock<double> porb01(lih, bas, 0, 1, 1, 1);

    const BufferHostXY<int32_t> pidx01({6, 9, 12}, 3, 1);

    const BufferHostX<int32_t> aidx01(std::vector<int32_t>({0}));

    const BufferHostMY<double, 5> pprim01({8.200000000000e-02, 1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00}, 1);

    const CBinnedGtoBlock<double> pdat01(1, 1, pidx01, aidx01, pprim01);

    ASSERT_EQ(porb01, pdat01);

    const CBinnedGtoBlock<double> porb11(lih, bas, 1, 1, 1, 1);

    const BufferHostXY<int32_t> pidx11({7, 10, 13}, 3, 1);

    const BufferHostX<int32_t> aidx11(std::vector<int32_t>({1}), 1);

    const BufferHostMY<double, 5> pprim11({8.000000000000e-01, 1.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 1.200000000000e+00}, 1);

    CBinnedGtoBlock<double> pdat11(1, 1, pidx11, aidx11, pprim11);

    ASSERT_EQ(porb11, pdat11);
}

TEST_F(CBinnedGtoBlockTest, CopyConstructor)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    CBinnedGtoBlock<double> bgto(agto);

    ASSERT_EQ(agto, bgto);
}

TEST_F(CBinnedGtoBlockTest, MoveConstructor)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    const CBinnedGtoBlock<double> bgto(CBinnedGtoBlock<double>(lih, bas, 0, 1));

    ASSERT_EQ(agto, bgto);
}

TEST_F(CBinnedGtoBlockTest, CopyAssignment)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    const CBinnedGtoBlock<double> bgto = agto;

    ASSERT_EQ(agto, bgto);
}

TEST_F(CBinnedGtoBlockTest, MoveAssignment)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> agto(lih, bas, 0, 1);

    const CBinnedGtoBlock<double> bgto = CBinnedGtoBlock<double>(lih, bas, 0, 1);

    ASSERT_EQ(agto, bgto);
}

TEST_F(CBinnedGtoBlockTest, GetAngularMomentum)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    ASSERT_EQ(sorb.getAngularMomentum(), 0);
}

TEST_F(CBinnedGtoBlockTest, Empty)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> porb(lih, bas, 1, 1);

    const CBinnedGtoBlock<double> dorb(lih, bas, 2, 1);

    ASSERT_TRUE(dorb.empty());

    ASSERT_FALSE(porb.empty());
}

TEST_F(CBinnedGtoBlockTest, GetNumberOfPrimGtos)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    ASSERT_EQ(sorb.getNumberOfPrimGtos(), 1);
}

TEST_F(CBinnedGtoBlockTest, GetNumberOfContrGtos)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    ASSERT_EQ(sorb.getNumberOfContrGtos(), 3);
}

TEST_F(CBinnedGtoBlockTest, GetAtomicIdentifiers)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({0, 0, 1}, sorb.getAtomicIdentifiers());
}

TEST_F(CBinnedGtoBlockTest, GetIdentifiers)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({1, 2, 4}, sorb.getIdentifiers(0));
}

TEST_F(CBinnedGtoBlockTest, GetExponents)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({5.281088472100e-02, 2.096094879800e-02, 1.219496200000e-01}, sorb.getExponents());
}

TEST_F(CBinnedGtoBlockTest, GetNormFactors)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({1.000000000000e+00, 1.000000000000e+00, 1.000000000000e+00}, sorb.getNormFactors());
}

TEST_F(CBinnedGtoBlockTest, GetCoordinatesX)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00}, sorb.getCoordinatesX());
}

TEST_F(CBinnedGtoBlockTest, GetCoordinatesY)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00}, sorb.getCoordinatesY());
}

TEST_F(CBinnedGtoBlockTest, GetCoordinatesZ)
{
    const CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorb(lih, bas, 0, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 1.200000000000e+00}, sorb.getCoordinatesZ());
}

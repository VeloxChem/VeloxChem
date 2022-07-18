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

#include "DiagEriRecForSSSSTest.hpp"

#include "BinnedGtoBlock.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "DiagEriRecForSSSS.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CDiagEriRecForSSSSTest, CompBatchPrimIntsForSSSS)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);

    const CBinnedGtoPairBlock<double> s11pairs(s1gtos);

    // [ss|ss] integrals

    auto s11ints = BufferHostX<double>::Zero(3);

    derirec::compHostSSSS(s11ints.data(), &s11pairs, 0, 3);

    const BufferHostX<double> refs11({0.123607744647421e+01, 0.201644152944101e+00, 0.100925300880806e+01}, 3);

    ASSERT_EQ(s11ints, refs11);
}

TEST_F(CDiagEriRecForSSSSTest, CompBatchContrIntsForSSSS)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoPairBlock<double> s22pairs(s2gtos);

    // (ss|ss) integrals

    auto s22ints = BufferHostX<double>::Zero(3);

    derirec::compHostSSSS(s22ints.data(), &s22pairs, 0, 3);

    const BufferHostX<double> refs22({0.162923244046733e+01, 0.379837636124064e-01, 0.177264071094880e+01}, 3);

    ASSERT_EQ(s22ints, refs22);
}

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

#include "BinnedGtoContainerTest.hpp"

#include <vector>

#include "BinnedGtoBlock.hpp"
#include "BinnedGtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CBinnedGtoContainerTest, DefaultConstructor)
{
    CBinnedGtoContainer<double> acont;

    ASSERT_EQ(acont, CBinnedGtoContainer(std::vector<CBinnedGtoBlock<double>>()));
}

TEST_F(CBinnedGtoContainerTest, ConstructorWithMolecule)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorbp1(lih, bas, 0, 1);
    
    const CBinnedGtoBlock<double> sorbp3(lih, bas, 0, 3);
    
    const CBinnedGtoBlock<double> sorbp5(lih, bas, 0, 5);

    const CBinnedGtoBlock<double> porbp1(lih, bas, 1, 1);
    
    const CBinnedGtoBlock<double> porbp2(lih, bas, 1, 2);

    CBinnedGtoContainer<double> acont({sorbp1, sorbp3, sorbp5, porbp1, porbp2});

    CBinnedGtoContainer<double> bcont(lih, bas);

    ASSERT_EQ(acont, bcont);
}

TEST_F(CBinnedGtoContainerTest, ConstructorWithAtomsList)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    const CBinnedGtoBlock<double> sorbp1a02(lih, bas, 0, 2, 0, 1);
    
    const CBinnedGtoBlock<double> sorbp3a02(lih, bas, 0, 2, 0, 3);
    
    const CBinnedGtoBlock<double> sorbp5a02(lih, bas, 0, 2, 0, 5);

    const CBinnedGtoBlock<double> porbp1a02(lih, bas, 0, 2, 1, 1);
    
    const CBinnedGtoBlock<double> porbp2a02(lih, bas, 0, 2, 1, 2);

    CBinnedGtoContainer<double> acont02({sorbp1a02, sorbp3a02, sorbp5a02, porbp1a02, porbp2a02});

    CBinnedGtoContainer<double> bcont02(lih, bas, 0, 2);

    ASSERT_EQ(acont02, bcont02);
    
    const CBinnedGtoBlock<double> sorbp1a01(lih, bas, 0, 1, 0, 1);
    
    const CBinnedGtoBlock<double> sorbp5a01(lih, bas, 0, 1, 0, 5);

    const CBinnedGtoBlock<double> porbp1a01(lih, bas, 0, 1, 1, 1);
    
    const CBinnedGtoBlock<double> porbp2a01(lih, bas, 0, 1, 1, 2);

    CBinnedGtoContainer<double> acont01({sorbp1a01, sorbp5a01, porbp1a01, porbp2a01});

    CBinnedGtoContainer<double> bcont01(lih, bas, 0, 1);

    ASSERT_EQ(acont01, bcont01);
    
    const CBinnedGtoBlock<double> sorbp1a11(lih, bas, 1, 1, 0, 1);
    
    const CBinnedGtoBlock<double> sorbp3a11(lih, bas, 1, 1, 0, 3);
    
    const CBinnedGtoBlock<double> porbp1a11(lih, bas, 1, 1, 1, 1);

    CBinnedGtoContainer<double> acont11({sorbp1a11, sorbp3a11, porbp1a11});

    CBinnedGtoContainer<double> bcont11(lih, bas, 1, 1);

    ASSERT_EQ(acont11, bcont11);
}

TEST_F(CBinnedGtoContainerTest, CopyConstructor)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    CBinnedGtoContainer<double> bcont(acont);

    ASSERT_EQ(acont, bcont);
}

TEST_F(CBinnedGtoContainerTest, MoveConstructor)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    CBinnedGtoContainer<double> bcont(CBinnedGtoContainer<double>(lih, bas));

    ASSERT_EQ(acont, bcont);
}

TEST_F(CBinnedGtoContainerTest, CopyAssignment)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    CBinnedGtoContainer<double> bcont = acont;

    ASSERT_EQ(acont, bcont);
}

TEST_F(CBinnedGtoContainerTest, MoveAssignment)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    CBinnedGtoContainer<double> bcont = CBinnedGtoContainer<double>(lih, bas);

    ASSERT_EQ(acont, bcont);
}

TEST_F(CBinnedGtoContainerTest, GetPointer)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);
    
    ASSERT_TRUE(acont.getPointer() == &acont);
}

TEST_F(CBinnedGtoContainerTest, GetMaxAngularMomentum)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getMaxAngularMomentum(), 1);
}

TEST_F(CBinnedGtoContainerTest, GetAngularMomentum)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getAngularMomentum(0), 0);
    
    ASSERT_EQ(acont.getAngularMomentum(1), 0);
    
    ASSERT_EQ(acont.getAngularMomentum(2), 0);
    
    ASSERT_EQ(acont.getAngularMomentum(3), 1);
    
    ASSERT_EQ(acont.getAngularMomentum(4), 1);
}

TEST_F(CBinnedGtoContainerTest, GetMaxNumberOfContrGtos)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getMaxNumberOfContrGtos(), 3);
}

TEST_F(CBinnedGtoContainerTest, GetNumberOfAtomicOrbitals)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getNumberOfAtomicOrbitals(), 14);
}

TEST_F(CBinnedGtoContainerTest, GetMaxNumberOfContrGtosWithAngularMomentum)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getMaxNumberOfContrGtos(0), 3);
    
    ASSERT_EQ(acont.getMaxNumberOfContrGtos(1), 2);
}

TEST_F(CBinnedGtoContainerTest, GetNumberOfBinnedGtoBlocks)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getNumberOfBinnedGtoBlocks(), 5);
}

TEST_F(CBinnedGtoContainerTest, GetBinnedGtoBlock)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    ASSERT_EQ(acont.getBinnedGtoBlock(0), CBinnedGtoBlock<double>(lih, bas, 0, 1));
    
    ASSERT_EQ(acont.getBinnedGtoBlock(1), CBinnedGtoBlock<double>(lih, bas, 0, 3));
    
    ASSERT_EQ(acont.getBinnedGtoBlock(2), CBinnedGtoBlock<double>(lih, bas, 0, 5));
    
    ASSERT_EQ(acont.getBinnedGtoBlock(3), CBinnedGtoBlock<double>(lih, bas, 1, 1));
    
    ASSERT_EQ(acont.getBinnedGtoBlock(4), CBinnedGtoBlock<double>(lih, bas, 1, 2));
}

TEST_F(CBinnedGtoContainerTest, PrintSummary)
{
    const auto bas = vlxbas::getMolecularBasisForLiH();

    const auto lih = vlxmol::getMoleculeLiH();

    CBinnedGtoContainer<double> acont(lih, bas);

    std::string str("Binned GTOs Container\n");
         
    str.append("=====================\n\n");
    
    str.append("Size: 5\n\n");

    str.append("Binned GTOs Block: L = S Number of GTOs:        3 Contraction Depth:    1\n");
    
    str.append("Binned GTOs Block: L = S Number of GTOs:        1 Contraction Depth:    3\n");
    
    str.append("Binned GTOs Block: L = S Number of GTOs:        1 Contraction Depth:    5\n");
    
    str.append("Binned GTOs Block: L = P Number of GTOs:        2 Contraction Depth:    1\n");
    
    str.append("Binned GTOs Block: L = P Number of GTOs:        1 Contraction Depth:    2\n");
    
    ASSERT_EQ(str, acont.printSummary());
}

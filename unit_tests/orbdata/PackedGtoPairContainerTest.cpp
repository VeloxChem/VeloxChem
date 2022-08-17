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

#include "PackedGtoPairContainerTest.hpp"

#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "PackedGtoPairContainer.hpp"

TEST_F(CPackedGtoPairContainerTest, DefaultConstructor)
{
    CPackedGtoPairContainer<double> acont;

    const auto bcont = CPackedGtoPairContainer<double>(VPackedGtoPairBlocks<double>());

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, CopyConstructor)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    const CPackedGtoPairContainer<double> bcont(acont);

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, MoveConstructor)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    const CPackedGtoPairContainer<double> bcont(CPackedGtoPairContainer<double>({ssblock, spblock, ppblock}));

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, CopyAssignment)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    const CPackedGtoPairContainer<double> bcont = acont;

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, MoveAssignment)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    const CPackedGtoPairContainer<double> bcont = CPackedGtoPairContainer<double>({ssblock, spblock, ppblock});

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, Add)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    CPackedGtoPairContainer<double> acont;

    ASSERT_EQ(acont.getNumberOfBlocks(), 0);

    acont.add(spblock);

    ASSERT_EQ(acont.getNumberOfBlocks(), 1);

    acont.add(ssblock);

    ASSERT_EQ(acont.getNumberOfBlocks(), 2);

    acont.add(ppblock);

    ASSERT_EQ(acont.getNumberOfBlocks(), 3);

    const CPackedGtoPairContainer<double> bcont({spblock, ssblock, ppblock});

    ASSERT_EQ(acont, bcont);
}

TEST_F(CPackedGtoPairContainerTest, To_Pointer)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    ASSERT_EQ(&acont, acont.to_pointer());
}

TEST_F(CPackedGtoPairContainerTest, GetNumberOfPackedGtoPairBlocks)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock});

    const CPackedGtoPairContainer<double> bcont({ssblock, spblock, ppblock});

    ASSERT_EQ(acont.getNumberOfBlocks(), 2);

    ASSERT_EQ(bcont.getNumberOfBlocks(), 3);
}

TEST_F(CPackedGtoPairContainerTest, GetPackedGtoPairBlock)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    ASSERT_EQ(acont.getBlock(0), ssblock);

    ASSERT_EQ(acont.getBlock(1), spblock);

    ASSERT_EQ(acont.getBlock(2), ppblock);
}

TEST_F(CPackedGtoPairContainerTest, PrintSummary)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);

    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);

    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);

    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);

    const BufferHostX<double> refssints({0.1, 0.003, 1.2}, 3);

    const BufferHostX<double> refspints({0.9, 0.003, 1.4, 0.7}, 4);

    const BufferHostX<double> refppints({0.9, 0.003, 1.9}, 3);

    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);

    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);

    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);

    const CPackedGtoPairContainer<double> acont({ssblock, spblock, ppblock});

    std::string str("Packed GTOs Pair Container\n");

    str.append("==========================\n\n");

    str.append("Size: 3\n\n");

    str.append("Packed GTOs Pair Block:\n");

    str.append("=======================\n\n");

    str.append("Angular Momentum: |SS>\n");

    str.append("Number of Contracted Pairs: 4\n\n");

    str.append("Index = 0 No. Pairs = 2\n");

    str.append("Index = 1 No. Pairs = 1\n\n");

    str.append("Packed GTOs Pair Block:\n");

    str.append("=======================\n\n");

    str.append("Angular Momentum: |SP>\n");

    str.append("Number of Contracted Pairs: 2\n\n");

    str.append("Index = 0 No. Pairs = 3\n");

    str.append("Index = 1 No. Pairs = 1\n\n");

    str.append("Packed GTOs Pair Block:\n");

    str.append("=======================\n\n");

    str.append("Angular Momentum: |PP>\n");

    str.append("Number of Contracted Pairs: 1\n\n");

    str.append("Index = 0 No. Pairs = 2\n");

    str.append("Index = 1 No. Pairs = 1\n\n");

    ASSERT_EQ(str, acont.printSummary());
}

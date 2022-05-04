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

#include "PackedGtoPairBlockTest.hpp"

#include "PackedGtoPairBlock.hpp"
#include "BinnedGtoBlock.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CPackedGtoPairBlockTest, ConstructorWithIntegrals)
{
    const auto mlih = vlxmol::getTestLiH();

    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);

    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);

    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);

    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});

    const CPackedGtoPairBlock<double> ablock(sspairs, refints);

    const BufferHostX<int32_t> indexes(4, {0, 1, 0, 0});

    const auto s0pairs = sspairs.compress(indexes, 0);

    const auto s1pairs = sspairs.compress(indexes, 1);

    const CPackedGtoPairBlock<double> bblock({s0pairs, s1pairs}, {0, 1});

    ASSERT_EQ(ablock, bblock);
}

TEST_F(CPackedGtoPairBlockTest, CopyConstructor)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);

    CPackedGtoPairBlock<double> bblock(ablock);

    ASSERT_EQ(ablock, bblock);
}

TEST_F(CPackedGtoPairBlockTest, MoveConstructor)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);

    CPackedGtoPairBlock<double> bblock(CPackedGtoPairBlock<double>(sspairs, refints));

    ASSERT_EQ(ablock, bblock);
}

TEST_F(CPackedGtoPairBlockTest, CopyAssignment)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);

    CPackedGtoPairBlock<double> bblock = ablock;

    ASSERT_EQ(ablock, bblock);
}

TEST_F(CPackedGtoPairBlockTest, MoveAssignment)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();

    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);

    CPackedGtoPairBlock<double> bblock = CPackedGtoPairBlock<double>(sspairs, refints);

    ASSERT_EQ(ablock, bblock);
}

TEST_F(CPackedGtoPairBlockTest, Size)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);
    
    ASSERT_EQ(ablock.size(), 16);
}

TEST_F(CPackedGtoPairBlockTest, SizeWithIndex)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);
    
    ASSERT_EQ(ablock.size(0), 3);
    
    ASSERT_EQ(ablock.size(1), 1);
    
    ASSERT_EQ(ablock.size(2), 0);
}

TEST_F(CPackedGtoPairBlockTest, GetGtoPairBlock)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);
    
    const BufferHostX<int32_t> indexes(4, {0, 1, 0, 0});
    
    const auto s0pairs = sspairs.compress(indexes, 0);
    
    const auto s1pairs = sspairs.compress(indexes, 1);
    
    ASSERT_EQ(ablock.getGtoPairBlock(0), s0pairs);
    
    ASSERT_EQ(ablock.getGtoPairBlock(1), s1pairs);
    
    ASSERT_EQ(ablock.getGtoPairBlock(2), CBinnedGtoPairBlock<double>());
}

TEST_F(CPackedGtoPairBlockTest, GetBraAngularMomentum)
{
    const auto mlih = vlxmol::getTestLiH();
       
    const auto mbas = vlxbas::getTestBasisForLiH();
       
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
       
    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);
       
    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);
       
    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);
       
    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);
    
    const BufferHostX<double> refssints(3, {0.1, 0.003, 1.2});
    
    const BufferHostX<double> refspints(4, {0.9, 0.003, 1.4, 0.7});
    
    const BufferHostX<double> refppints(3, {0.9, 0.003, 1.9});
    
    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);
    
    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);
    
    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);
    
    ASSERT_EQ(ssblock.getBraAngularMomentum(), 0);
    
    ASSERT_EQ(spblock.getBraAngularMomentum(), 0);
    
    ASSERT_EQ(ppblock.getBraAngularMomentum(), 1);
}

TEST_F(CPackedGtoPairBlockTest, GetKetAngularMomentum)
{
    const auto mlih = vlxmol::getTestLiH();
       
    const auto mbas = vlxbas::getTestBasisForLiH();
       
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
       
    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);
       
    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);
       
    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);
       
    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);
    
    const BufferHostX<double> refssints(3, {0.1, 0.003, 1.2});
    
    const BufferHostX<double> refspints(4, {0.9, 0.003, 1.4, 0.7});
    
    const BufferHostX<double> refppints(3, {0.9, 0.003, 1.9});
    
    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);
    
    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);
    
    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);
    
    ASSERT_EQ(ssblock.getKetAngularMomentum(), 0);
    
    ASSERT_EQ(spblock.getKetAngularMomentum(), 1);
    
    ASSERT_EQ(ppblock.getKetAngularMomentum(), 1);
}

TEST_F(CPackedGtoPairBlockTest, GetNumberOfPrimGtoPairs)
{
    const auto mlih = vlxmol::getTestLiH();
       
    const auto mbas = vlxbas::getTestBasisForLiH();
       
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
       
    const CBinnedGtoBlock<double> p1gtos(mlih, mbas, 1, 1);
       
    const CBinnedGtoPairBlock<double> s2s2pairs(s2gtos);
       
    const CBinnedGtoPairBlock<double> s2p1pairs(s2gtos, p1gtos);
       
    const CBinnedGtoPairBlock<double> p1p1pairs(p1gtos);
    
    const BufferHostX<double> refssints(3, {0.1, 0.003, 1.2});
    
    const BufferHostX<double> refspints(4, {0.9, 0.003, 1.4, 0.7});
    
    const BufferHostX<double> refppints(3, {0.9, 0.003, 1.9});
    
    const CPackedGtoPairBlock<double> ssblock(s2s2pairs, refssints);
    
    const CPackedGtoPairBlock<double> spblock(s2p1pairs, refspints);
    
    const CPackedGtoPairBlock<double> ppblock(p1p1pairs, refppints);
    
    ASSERT_EQ(ssblock.getNumberOfPrimGtoPairs(), 4);
    
    ASSERT_EQ(spblock.getNumberOfPrimGtoPairs(), 2);
    
    ASSERT_EQ(ppblock.getNumberOfPrimGtoPairs(), 1);
}

TEST_F(CPackedGtoPairBlockTest, GetNumberOfSubBlocks)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(0, 3), 1);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(0, 2), 2);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(0, 1), 3);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(1, 3), 1);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(1, 2), 1);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(1, 1), 1);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(2, 3), 0);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(2, 2), 0);
    
    ASSERT_EQ(ablock.getNumberOfSubBlocks(2, 1), 0);
}

TEST_F(CPackedGtoPairBlockTest, PrintSummary)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getTestBasisForLiH();
    
    const CBinnedGtoBlock<double> s2gtos(mlih, mbas, 0, 2);
    
    const CBinnedGtoBlock<double> s1gtos(mlih, mbas, 0, 1);
    
    const CBinnedGtoPairBlock<double> sspairs(s2gtos, s1gtos);
    
    const BufferHostX<double> refints(4, {0.9, 0.003, 1.4, 0.7});
    
    const CPackedGtoPairBlock<double> ablock(sspairs, refints);
    
    std::string str("Packed GTOs Pair Block:\n");
         
    str.append("=======================\n\n");
    
    str.append("Angular Momentum: |SS>\n");
    
    str.append("Number of Contracted Pairs: 2\n\n");

    str.append("Index = 0 No. Pairs = 3\n");
    
    str.append("Index = 1 No. Pairs = 1\n");

    ASSERT_EQ(str, ablock.printSummary());
}

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

#include "CommonNeighborsTest.hpp"

#include "CommonNeighbors.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CCommonNeighborsTest, DefaultConstructor)
{
    CCommonNeighbors cna;
    
    CCommonNeighbors cnb(0.0,
                         CMemBlock<int32_t>(),
                         CDenseMatrix(),
                         std::vector<CFourIndexes>(),
                         std::vector<int32_t>(),
                         std::vector<int32_t>(),
                         std::set<int32_t>());

    ASSERT_EQ(cna, cnb);
}

TEST_F(CCommonNeighborsTest, CopyConstructor)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 2.0);

    CCommonNeighbors cnb(cna);

    ASSERT_EQ(cna, cnb);
}

TEST_F(CCommonNeighborsTest, MoveConstructor)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 2.0);

    CCommonNeighbors cnb(CCommonNeighbors(vlxmol::getMoleculeNH3CH4(), 2.0));

    ASSERT_EQ(cna, cnb);
}

TEST_F(CCommonNeighborsTest, CopyAssignment)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 2.0);

    CCommonNeighbors cnb = cna;

    ASSERT_EQ(cna, cnb);
}

TEST_F(CCommonNeighborsTest, MoveAssignment)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 2.0);

    CCommonNeighbors cnb = CCommonNeighbors(vlxmol::getMoleculeNH3CH4(), 2.0);

    ASSERT_EQ(cna, cnb);
}

TEST_F(CCommonNeighborsTest, GetBonds)
{
    CCommonNeighbors cna(vlxmol::getMoleculeH2O(), 2.0);
    
    const CDenseMatrix refbonds({0.00000000000000000000, 1.78044938147648555945, 1.78044938147648555945,
                                 1.78044938147648555945, 0.00000000000000000000, 2.80000000000000000000,
                                 1.78044938147648555945, 2.80000000000000000000, 0.00000000000000000000},
                                3, 3);

    ASSERT_EQ(refbonds, cna.getBonds());
}

TEST_F(CCommonNeighborsTest, GetAdjacencies)
{
    CCommonNeighbors cna(vlxmol::getMoleculeH2O(), 2.0);
    
    const CMemBlock<int32_t> refadjs({0, 1, 1, 1, 0, 0, 1, 0, 0});

    ASSERT_EQ(refadjs, cna.getAdjacencies());
}

TEST_F(CCommonNeighborsTest, GetBondPairs)
{
    CCommonNeighbors cna(vlxmol::getMoleculeH2O(), 2.0);
    
    const std::vector<CTwoIndexes> refpairs({{0, 1}, {0, 2}});

    ASSERT_EQ(refpairs, cna.getBondPairs());
}

TEST_F(CCommonNeighborsTest, Generate)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);

    ASSERT_EQ(cna.getRepetitions(), std::vector<int32_t>({3, 4}));
    
    ASSERT_EQ(cna.getSignatures(), std::vector<CFourIndexes>({CFourIndexes(2, 7, 4, 2),
                                                              CFourIndexes(1, 7, 3, 2)}));
}

TEST_F(CCommonNeighborsTest, GetRepetitions)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);

    ASSERT_EQ(cna.getRepetitions(), std::vector<int32_t>({3, 4}));
}

TEST_F(CCommonNeighborsTest, GetSignatures)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);
    
    ASSERT_EQ(cna.getSignatures(), std::vector<CFourIndexes>({CFourIndexes(2, 7, 4, 2),
                                                              CFourIndexes(1, 7, 3, 2)}));
}

TEST_F(CCommonNeighborsTest, CompJaccardIndex)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    CCommonNeighbors cnb(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);
    
    cnb.generate(3.0);

    ASSERT_NEAR(cna.compJaccardIndex(cna), 1.0, 1.0e-13);
    
    ASSERT_NEAR(cnb.compJaccardIndex(cnb), 1.0, 1.0e-13);
    
    ASSERT_NEAR(cna.compJaccardIndex(cnb), 0.0, 1.0e-13);
    
    ASSERT_NEAR(cnb.compJaccardIndex(cna), 0.0, 1.0e-13);
}

TEST_F(CCommonNeighborsTest, Find)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);
    
    ASSERT_EQ(cna.find(CFourIndexes(2, 7, 4, 2)), 0);
    
    ASSERT_EQ(cna.find(CFourIndexes(1, 7, 3, 2)), 1);
    
    ASSERT_EQ(cna.find(CFourIndexes(1, 7, 3, 1)), -1);
}

TEST_F(CCommonNeighborsTest, GetSignaturesRepr)
{
    CCommonNeighbors cna(vlxmol::getMoleculeNH3CH4(), 3.0);
    
    cna.generate(15.0);
    
    const auto str = cna.getSignaturesRepr();
    
    const std::string refstr("2 : (7,4,2) : 3\n1 : (7,3,2) : 4\n");
    
    ASSERT_EQ(str, refstr);
}

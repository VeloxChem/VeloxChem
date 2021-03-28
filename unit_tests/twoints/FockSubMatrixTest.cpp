//
//                           VELOXCHEM 1.0-RC
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

#include "FockSubMatrixTest.hpp"

#include "FockSubMatrix.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CFockSubMatrixTest, DefaultConstructor)
{
    CFockSubMatrix subma;
    
    CMemBlock<int32_t> midx;
    
    CFockSubMatrix submb({}, midx, midx, midx, midx, 0, 0, 0, 0);
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, ConstructorWithGtoPairs)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjk);
    
    CMemBlock<int32_t> sidx(1); sidx.zero();
    
    CMemBlock<int32_t> pidx({5, 8, 11});
    
    CMemBlock2D<double> mss(25, 1); mss.zero();
    
    CMemBlock2D<double> msp(15, 3); msp.zero();
    
    CFockSubMatrix submb({mss, msp, mss, msp, mss, msp},
                         sidx, sidx, sidx, pidx, 5, 5, 5, 3);
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, CopyConstructor)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restk);
    
    CFockSubMatrix submb(subma);
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, MoveConstructor)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restk);
    
    CFockSubMatrix submb(CFockSubMatrix(bpairs, kpairs, fockmat::restk));
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, CopyAssignment)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restj);
    
    CFockSubMatrix submb = subma;
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, MoveAssignment)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restj);
    
    CFockSubMatrix submb = CFockSubMatrix(bpairs, kpairs, fockmat::restj);
    
    ASSERT_EQ(subma, submb);
}

TEST_F(CFockSubMatrixTest, GetSubMatrixData)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, bpairs, fockmat::restjkx);
    
    std::vector<double> vecsp(15, 0.0);
    
    for (int32_t i = 0; i < 3; i++)
    {
        vlxtest::compare(vecsp, subma.getSubMatrixData(0, i));
        
        vlxtest::compare(vecsp, subma.getSubMatrixData(1, i));
        
        vlxtest::compare(vecsp, subma.getSubMatrixData(3, i));
        
        vlxtest::compare(vecsp, subma.getSubMatrixData(4, i));
    }
    
    std::vector<double> vecss(25, 0.0);
    
    vlxtest::compare(vecsp, subma.getSubMatrixData(2, 0));
    
    std::vector<double> vecpp(9, 0.0);
    
    for (int32_t i = 0; i < 9; i++)
    {
        vlxtest::compare(vecpp, subma.getSubMatrixData(5, i));
    }
}

TEST_F(CFockSubMatrixTest, GetStartPositionsA)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    vlxtest::compare({0}, subma.getStartPositionsA());
    
    vlxtest::compare({5, 8, 11}, submb.getStartPositionsA());
}

TEST_F(CFockSubMatrixTest, GetStartPositionsB)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    vlxtest::compare({5, 8, 11}, subma.getStartPositionsB());
    
    vlxtest::compare({5, 8, 11}, submb.getStartPositionsB());
}

TEST_F(CFockSubMatrixTest, GetStartPositionsC)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    vlxtest::compare({5, 8, 11}, subma.getStartPositionsC());
    
    vlxtest::compare({0}, submb.getStartPositionsC());
}

TEST_F(CFockSubMatrixTest, GetStartPositionsD)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    vlxtest::compare({5, 8, 11}, subma.getStartPositionsD());
    
    vlxtest::compare({5, 8, 11}, submb.getStartPositionsD());
}

TEST_F(CFockSubMatrixTest, GetDimensionsA)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    ASSERT_EQ(5, subma.getDimensionsA());
    
    ASSERT_EQ(3, submb.getDimensionsA());
}

TEST_F(CFockSubMatrixTest, GetDimensionsB)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    ASSERT_EQ(3, subma.getDimensionsB());
    
    ASSERT_EQ(3, submb.getDimensionsB());
}

TEST_F(CFockSubMatrixTest, GetDimensionsC)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    ASSERT_EQ(3, subma.getDimensionsC());
    
    ASSERT_EQ(5, submb.getDimensionsC());
}

TEST_F(CFockSubMatrixTest, GetDimensionsD)
{
    CMolecularBasis mbas = vlxbas::getMolecularBasisForLiH();
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, kgtos, 1.0e-13);
    
    CFockSubMatrix subma(bpairs, kpairs, fockmat::restjkx);
    
    CFockSubMatrix submb(kpairs, bpairs, fockmat::restj);
    
    ASSERT_EQ(3, subma.getDimensionsD());
    
    ASSERT_EQ(3, submb.getDimensionsD());
}

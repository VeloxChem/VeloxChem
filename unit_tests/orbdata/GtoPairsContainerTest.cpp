//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoPairsContainerTest.hpp"

#include "GtoPairsBlock.hpp"
#include "GtoPairsContainer.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(CGtoPairsContainerTest, DefaultConstructor)
{
    std::vector<CGtoPairsBlock> ppbloks;
    
    CGtoPairsContainer acont;
    
    CGtoPairsContainer bcont(ppbloks);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, ConstructorWithMolecule)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    CGtoBlock porb(lih, bas, 1);
    
    CGtoPairsBlock bss(sorb, 1.0e-13);
    
    CGtoPairsBlock bsp(sorb, porb, 1.0e-13);
    
    CGtoPairsBlock bpp(porb, 1.0e-13);
    
    CGtoPairsContainer acont({bss, bsp, bpp});
    
    CGtoPairsContainer bcont(lih, bas, 1.0e-13);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont(acont);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont(CGtoPairsContainer(lih, bas, 1.0e-13));
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont = acont;
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont = CGtoPairsContainer(lih, bas, 1.0e-13);
    
    ASSERT_EQ(acont, bcont);
}

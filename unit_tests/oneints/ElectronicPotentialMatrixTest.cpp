//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronicPotentialMatrixTest.hpp"

#include "ElectronicPotentialMatrix.hpp"

TEST_F(CElectronicPotentialMatrixTest, DefaultConstructor)
{
    CElectronicPotentialMatrix smata;
    
    CElectronicPotentialMatrix smatb(CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

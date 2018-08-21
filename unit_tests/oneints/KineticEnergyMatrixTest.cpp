//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KineticEnergyMatrixTest.hpp"

#include "KineticEnergyMatrix.hpp"

TEST_F(CKineticEnergyMatrixTest, DefaultConstructor)
{
    CKineticEnergyMatrix smata;
    
    CKineticEnergyMatrix smatb(CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CKineticEnergyMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CKineticEnergyMatrix smata(ma);
    
    CKineticEnergyMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CKineticEnergyMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CKineticEnergyMatrix smata(ma);
    
    CKineticEnergyMatrix smatb(CKineticEnergyMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CKineticEnergyMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CKineticEnergyMatrix smata(ma);
    
    CKineticEnergyMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CKineticEnergyMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
        
    CKineticEnergyMatrix smata(ma);
    
    CKineticEnergyMatrix smatb(CKineticEnergyMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

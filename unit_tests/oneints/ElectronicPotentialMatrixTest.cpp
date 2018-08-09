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
    
    CElectronicPotentialMatrix smatb(CSparseMatrix(0, 0, 1.0e-13));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveConstructor)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveAssignment)
{
    CSparseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0,
                      8.0, -5.0}, {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4},
                     {0, 1, 2, 0, 1, 2, 3, 4, 9, 2, 3, 1, 4}, 5, 5, 1.0e-13);
    
    CElectronicPotentialMatrix smata(ma);
    
    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

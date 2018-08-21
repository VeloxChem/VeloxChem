//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "NuclearPotentialMatrixTest.hpp"

#include "NuclearPotentialMatrix.hpp"

TEST_F(CNuclearPotentialMatrixTest, DefaultConstructor)
{
    CNuclearPotentialMatrix smata;
    
    CNuclearPotentialMatrix smatb(CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb(CNuclearPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CNuclearPotentialMatrix smata(ma);
    
    CNuclearPotentialMatrix smatb(CNuclearPotentialMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

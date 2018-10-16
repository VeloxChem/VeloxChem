//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularOrbitalsTest.hpp"

#include "MolecularOrbitals.hpp"

TEST_F(CMolecularOrbitalsTest, DefaultConstructor)
{
    CMolecularOrbitals moa;
    
    CMolecularOrbitals mob({ }, molorb::rest);
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMolecularOrbitals moa({ma, mb}, molorb::unrest);
    
    CMolecularOrbitals mob(moa);
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CMolecularOrbitals moa({ma, mb}, molorb::unrest);
    
    CMolecularOrbitals mob(CMolecularOrbitals({ma, mb}, molorb::unrest));
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CMolecularOrbitals moa({ma}, molorb::rest);
    
    CMolecularOrbitals mob = moa;
    
    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CMolecularOrbitals moa({ma}, molorb::rest);
    
    CMolecularOrbitals mob = CMolecularOrbitals({ma}, molorb::rest);
    
    ASSERT_EQ(moa, mob);
}

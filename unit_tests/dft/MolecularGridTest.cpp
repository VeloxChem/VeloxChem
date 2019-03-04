//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MolecularGridTest.hpp"

#include "MolecularGrid.hpp"

TEST_F(CMolecularGridTest, DefaultConstructor)
{
    CMolecularGrid mgrida;
    
    CMemBlock2D<double> mblock;
    
    CMolecularGrid mgridb(mblock);
    
    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, CopyConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    CMolecularGrid mgridb(mgrida);
    
    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    CMolecularGrid mgridb((CMolecularGrid(mblock)));
    
    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, CopyAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    CMolecularGrid mgridb = mgrida;
    
    ASSERT_EQ(mgrida, mgridb);
}

TEST_F(CMolecularGridTest, MoveAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CMolecularGrid mgrida(mblock);
    
    CMolecularGrid mgridb = CMolecularGrid(mblock);
    
    ASSERT_EQ(mgrida, mgridb);
}

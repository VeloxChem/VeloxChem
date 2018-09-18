//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OverlapMatrixTest.hpp"

#include "OverlapMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(COverlapMatrixTest, DefaultConstructor)
{
    COverlapMatrix smata;
    
    COverlapMatrix smatb(CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(COverlapMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    COverlapMatrix smata(ma);
    
    COverlapMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(COverlapMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    COverlapMatrix smata(ma);
    
    COverlapMatrix smatb(COverlapMatrix({ma}));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(COverlapMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    COverlapMatrix smata(ma);
    
    COverlapMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(COverlapMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    COverlapMatrix smata(ma);
    
    COverlapMatrix smatb = COverlapMatrix({ma});
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(COverlapMatrixTest, GetString)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    COverlapMatrix smata(ma);
    
    ASSERT_EQ(ma.getString(), smata.getString());
}

TEST_F(COverlapMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    COverlapMatrix smata(ma);
    
    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(COverlapMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    COverlapMatrix smata(ma);
    
    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(COverlapMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    COverlapMatrix smata(ma);
    
    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(COverlapMatrixTest, Values)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const COverlapMatrix smata(ma);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, smata.values());
}

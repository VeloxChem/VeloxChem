//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "LinearMomentumMatrixTest.hpp"

#include "LinearMomentumMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CLinearMomentumMatrixTest, DefaultConstructor)
{
    CLinearMomentumMatrix smata;
    
    CLinearMomentumMatrix smatb(CDenseMatrix(0, 0), CDenseMatrix(0, 0),
                                CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CLinearMomentumMatrixTest, CopyConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    CLinearMomentumMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CLinearMomentumMatrixTest, MoveConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    CLinearMomentumMatrix smatb(CLinearMomentumMatrix(mx, my, mz));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CLinearMomentumMatrixTest, CopyAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    CLinearMomentumMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CLinearMomentumMatrixTest, MoveAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    CLinearMomentumMatrix smatb = CLinearMomentumMatrix(mx, my, mz);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CLinearMomentumMatrixTest, GetStringForComponentX)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    ASSERT_EQ(mx.getString(), smata.getStringForComponentX());
}

TEST_F(CLinearMomentumMatrixTest, GetStringForComponentY)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    ASSERT_EQ(my.getString(), smata.getStringForComponentY());
}

TEST_F(CLinearMomentumMatrixTest, GetStringForComponentZ)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CLinearMomentumMatrix smata(mx, my, mz);
    
    ASSERT_EQ(mz.getString(), smata.getStringForComponentZ());
}

TEST_F(CLinearMomentumMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CLinearMomentumMatrix smata(ma, ma, ma);
    
    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CLinearMomentumMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CLinearMomentumMatrix smata(ma, ma, ma);
    
    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CLinearMomentumMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CLinearMomentumMatrix smata(ma, ma, ma);
    
    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CLinearMomentumMatrixTest, XValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    const CLinearMomentumMatrix smata(mx, my, mz);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     smata.xvalues());
}

TEST_F(CLinearMomentumMatrixTest, YValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    const CLinearMomentumMatrix smata(mx, my, mz);
    
    vlxtest::compare({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0},
                     smata.yvalues());
}

TEST_F(CLinearMomentumMatrixTest, ZValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    const CLinearMomentumMatrix smata(mx, my, mz);
    
    vlxtest::compare({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0},
                     smata.zvalues());
}

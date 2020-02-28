//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AngularMomentumMatrixTest.hpp"

#include "AngularMomentumMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAngularMomentumMatrixTest, DefaultConstructor)
{
    CAngularMomentumMatrix smata;

    CAngularMomentumMatrix smatb(CDenseMatrix(0, 0), CDenseMatrix(0, 0), CDenseMatrix(0, 0), 0.0, 0.0, 0.0);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, CopyConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.0, 2.0, 4.0);

    CAngularMomentumMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, MoveConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 2.0, 3.0, 5.0);

    CAngularMomentumMatrix smatb(CAngularMomentumMatrix(mx, my, mz, 2.0, 3.0, 5.0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, CopyAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 2.9);

    CAngularMomentumMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, MoveAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 0.3, -2.1, 1.8);

    CAngularMomentumMatrix smatb = CAngularMomentumMatrix(mx, my, mz, 0.3, -2.1, 1.8);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, SetOriginCoordinates)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    smata.setOriginCoordinates(2.3, 1.6, 8.1);

    CAngularMomentumMatrix smatb(mx, my, mz, 2.3, 1.6, 8.1);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentX)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(mx.getString(), smata.getStringForComponentX());
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentY)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(my.getString(), smata.getStringForComponentY());
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentZ)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(mz.getString(), smata.getStringForComponentZ());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 1.2, 3.2, 8.9);

    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 0.2, 3.1, 4.5);

    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 1.2, 6.7, 8.1);

    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CAngularMomentumMatrixTest, XValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, smata.xvalues());
}

TEST_F(CAngularMomentumMatrixTest, YValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, smata.yvalues());
}

TEST_F(CAngularMomentumMatrixTest, ZValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, smata.zvalues());
}

TEST_F(CAngularMomentumMatrixTest, GetOriginCoordinateX)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    ASSERT_NEAR(1.2, smata.getOriginCoordinateX(), 1.0e-13);
}

TEST_F(CAngularMomentumMatrixTest, GetOriginCoordinateY)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    ASSERT_NEAR(3.4, smata.getOriginCoordinateY(), 1.0e-13);
}

TEST_F(CAngularMomentumMatrixTest, GetOriginCoordinateZ)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    ASSERT_NEAR(6.1, smata.getOriginCoordinateZ(), 1.0e-13);
}

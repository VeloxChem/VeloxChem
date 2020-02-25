//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricFieldMatrixTest.hpp"

#include "CheckFunctions.hpp"
#include "ElectricFieldMatrix.hpp"

TEST_F(CElectricFieldMatrixTest, DefaultConstructor)
{
    CElectricFieldMatrix smata;

    CElectricFieldMatrix smatb(CDenseMatrix(0, 0), CDenseMatrix(0, 0), CDenseMatrix(0, 0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldMatrixTest, CopyConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    CElectricFieldMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldMatrixTest, MoveConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    CElectricFieldMatrix smatb(CElectricFieldMatrix(mx, my, mz));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldMatrixTest, CopyAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    CElectricFieldMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldMatrixTest, MoveAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    CElectricFieldMatrix smatb = CElectricFieldMatrix(mx, my, mz);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldMatrixTest, GetStringForComponentX)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    ASSERT_EQ(mx.getString(), smata.getStringForComponentX());
}

TEST_F(CElectricFieldMatrixTest, GetStringForComponentY)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    ASSERT_EQ(my.getString(), smata.getStringForComponentY());
}

TEST_F(CElectricFieldMatrixTest, GetStringForComponentZ)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CElectricFieldMatrix smata(mx, my, mz);

    ASSERT_EQ(mz.getString(), smata.getStringForComponentZ());
}

TEST_F(CElectricFieldMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CElectricFieldMatrix smata(ma, ma, ma);

    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CElectricFieldMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CElectricFieldMatrix smata(ma, ma, ma);

    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CElectricFieldMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CElectricFieldMatrix smata(ma, ma, ma);

    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CElectricFieldMatrixTest, XValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CElectricFieldMatrix smata(mx, my, mz);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, smata.xvalues());
}

TEST_F(CElectricFieldMatrixTest, YValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CElectricFieldMatrix smata(mx, my, mz);

    vlxtest::compare({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, smata.yvalues());
}

TEST_F(CElectricFieldMatrixTest, ZValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CElectricFieldMatrix smata(mx, my, mz);

    vlxtest::compare({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, smata.zvalues());
}

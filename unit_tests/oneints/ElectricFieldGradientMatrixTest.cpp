//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldGradientMatrixTest.hpp"

#include "ElectricFieldGradientMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CElectricFieldGradientMatrixTest, DefaultConstructor)
{
    CElectricFieldGradientMatrix smata;
    
    CElectricFieldGradientMatrix smatb(CDenseMatrix(0, 0), CDenseMatrix(0, 0), CDenseMatrix(0, 0),
                                       CDenseMatrix(0, 0), CDenseMatrix(0, 0), CDenseMatrix(0, 0));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldGradientMatrixTest, CopyConstructor)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    CElectricFieldGradientMatrix smatb(smata);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldGradientMatrixTest, MoveConstructor)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    CElectricFieldGradientMatrix smatb(CElectricFieldGradientMatrix(mxx, mxy, mxz, myy, myz, mzz));
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldGradientMatrixTest, CopyAssignment)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    CElectricFieldGradientMatrix smatb = smata;
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldGradientMatrixTest, MoveAssignment)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    CElectricFieldGradientMatrix smatb = CElectricFieldGradientMatrix(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentXX)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(mxx.getString(), smata.getStringForComponentXX());
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentXY)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(mxy.getString(), smata.getStringForComponentXY());
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentXZ)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(mxz.getString(), smata.getStringForComponentXZ());
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentYY)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(myy.getString(), smata.getStringForComponentYY());
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentYZ)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(myz.getString(), smata.getStringForComponentYZ());
}

TEST_F(CElectricFieldGradientMatrixTest, GetStringForComponentZZ)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    ASSERT_EQ(mzz.getString(), smata.getStringForComponentZZ());
}


TEST_F(CElectricFieldGradientMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CElectricFieldGradientMatrix smata(ma, ma, ma, ma, ma, ma);
    
    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CElectricFieldGradientMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CElectricFieldGradientMatrix smata(ma, ma, ma, ma, ma, ma);
    
    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CElectricFieldGradientMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CElectricFieldGradientMatrix smata(ma, ma, ma, ma, ma, ma);
    
    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CElectricFieldGradientMatrixTest, XXValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     smata.xxvalues());
}

TEST_F(CElectricFieldGradientMatrixTest, XYValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0},
                     smata.xyvalues());
}

TEST_F(CElectricFieldGradientMatrixTest, XZValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0},
                     smata.xzvalues());
}

TEST_F(CElectricFieldGradientMatrixTest, YYValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     smata.yyvalues());
}

TEST_F(CElectricFieldGradientMatrixTest, YZValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0},
                     smata.yzvalues());
}

TEST_F(CElectricFieldGradientMatrixTest, ZZValues)
{
    CDenseMatrix mxx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mxy({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);
    
    CDenseMatrix mxz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);
    
    CDenseMatrix myy({0.0, -1.0, -1.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix myz({1.0, -2.0, -8.0, -6.0, 6.0, 1.0, 3.0, 4.0, -7.0}, 3, 3);
    
    CDenseMatrix mzz({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0}, 3, 3);
    
    CElectricFieldGradientMatrix smata(mxx, mxy, mxz, myy, myz, mzz);
    
    vlxtest::compare({7.0, -1.0, -2.0, -3.0, 5.0, 8.0, 1.0, 3.0, -4.0},
                     smata.zzvalues());
}





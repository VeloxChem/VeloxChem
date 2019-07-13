//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AOFockMatrixTest.hpp"

#include "AOFockMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAOFockMatrixTest, DefaultConstructor)
{
    CAOFockMatrix fmata;
    
    CAOFockMatrix fmatb({}, {}, {}, {});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, ConstructorWithDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAODensityMatrix dmat({ma, mb}, denmat::rest);
    
    CAOFockMatrix fmata(dmat);
    
    ASSERT_EQ(2, fmata.getNumberOfFockMatrices());

    ASSERT_EQ(3, fmata.getNumberOfRows(0));
    
    ASSERT_EQ(3, fmata.getNumberOfRows(1));
    
    ASSERT_EQ(3, fmata.getNumberOfColumns(0));
    
    ASSERT_EQ(2, fmata.getNumberOfColumns(1));
    
    ASSERT_EQ(0, fmata.getDensityIdentifier(0));
    
    ASSERT_EQ(1, fmata.getDensityIdentifier(1));
    
    ASSERT_TRUE(fockmat::restjk == fmata.getFockType(0));
    
    ASSERT_TRUE(fockmat::restjk == fmata.getFockType(1));
}

TEST_F(CAOFockMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    CAOFockMatrix fmatb(fmata);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    CAOFockMatrix fmatb(CAOFockMatrix({ma, ma, mb, mb}, {fockmat::restj,
                                      fockmat::restk, fockmat::restj,
                                      fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                                      {0, 2, 3, 6}));
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    CAOFockMatrix fmatb(fmata);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                         fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                         {0, 2, 3, 6});
    
    CAOFockMatrix fmatb = CAOFockMatrix({ma, ma, mb, mb}, {fockmat::restj,
                                        fockmat::restk, fockmat::restj,
                                        fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, SetFockType)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    fmata.setFockType(fockmat::restjkx, 0);
    
    fmata.setFockType(fockmat::restkx, 1);
    
    fmata.setFockType(fockmat::restjk, 3);
    
    CAOFockMatrix fmatb({ma, ma, mb, mb}, {fockmat::restjkx, fockmat::restkx,
                        fockmat::restj, fockmat::restjk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, SetFockScaleFactor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    fmata.setFockScaleFactor(3.2, 0);
    
    fmata.setFockScaleFactor(2.1, 1);
    
    fmata.setFockScaleFactor(0.5, 3);
    
    CAOFockMatrix fmatb({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {3.2, 2.1, 0.0, 0.5},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, Zero)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    fmata.zero();
    
    ma.zero();
    
    mb.zero();
    
    CAOFockMatrix fmatb({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, Symmetrize)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    fmata.symmetrize();

    CDenseMatrix mc({2.0, -3.0, 3.0, -3.0, 10.0, 8.0, 3.0, 8.0, -8.0}, 3, 3);
    
    CAOFockMatrix fmatb({mc, mc, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, Add)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    CAOFockMatrix fmatb({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    fmatb.add(fmata);
    
    CDenseMatrix mc({2.0, -2.0, -6.0, -4.0, 10.0, 8.0, 12.0, 8.0, -8.0}, 3, 3);
    
    CDenseMatrix md({2.0, -2.0, -6.0, -4.0, 10.0, 8.0}, 3, 2);
    
    CAOFockMatrix fmatc({mc, mc, md, md}, {fockmat::restj, fockmat::restk,
                         fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmatb, fmatc);
}

TEST_F(CAOFockMatrixTest, GetNumberOfFockMatrices)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(4, fmata.getNumberOfFockMatrices());
}

TEST_F(CAOFockMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(3, fmata.getNumberOfRows(0));
    
    ASSERT_EQ(3, fmata.getNumberOfRows(1));
    
    ASSERT_EQ(3, fmata.getNumberOfRows(2));
    
    ASSERT_EQ(3, fmata.getNumberOfRows(3));
}

TEST_F(CAOFockMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                         fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                         {0, 2, 3, 6});
    
    ASSERT_EQ(3, fmata.getNumberOfColumns(0));
    
    ASSERT_EQ(3, fmata.getNumberOfColumns(1));
    
    ASSERT_EQ(2, fmata.getNumberOfColumns(2));
    
    ASSERT_EQ(2, fmata.getNumberOfColumns(3));
}

TEST_F(CAOFockMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(9, fmata.getNumberOfElements(0));
    
    ASSERT_EQ(9, fmata.getNumberOfElements(1));
    
    ASSERT_EQ(6, fmata.getNumberOfElements(2));
    
    ASSERT_EQ(6, fmata.getNumberOfElements(3));
}

TEST_F(CAOFockMatrixTest, GetFockConstant)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                              fockmat::restj, fockmat::restk},
                              {1.0, 2.0, 0.0, 1.0} , {0, 2, 3, 6});
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getFock(0));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getFock(1));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getFock(2));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getFock(3));
}

TEST_F(CAOFockMatrixTest, GetFock)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getFock(0));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getFock(1));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getFock(2));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getFock(3));
    
    auto pdat = fmata.getFock(1);
    
    pdat[2] = 2.0;  pdat[5] = 0.0;
    
    CDenseMatrix mc({1.0, -1.0, 2.0, -2.0, 5.0, 0.0, 6.0, 4.0, -4.0}, 3, 3);
    
    ASSERT_EQ(fmata, CAOFockMatrix({ma, mc, mb, mb}, {fockmat::restj,
                                   fockmat::restk, fockmat::restj,
                                   fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                                   {0, 2, 3, 6}));
}

TEST_F(CAOFockMatrixTest, GetReferenceToFock)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                              fockmat::restj, fockmat::restk},
                              {1.0, 2.0, 0.0, 1.0} , {0, 2, 3, 6});
    
    ASSERT_EQ(ma, fmata.getReferenceToFock(0));
    
    ASSERT_EQ(ma, fmata.getReferenceToFock(1));
    
    ASSERT_EQ(mb, fmata.getReferenceToFock(2));
    
    ASSERT_EQ(mb, fmata.getReferenceToFock(3));
}

TEST_F(CAOFockMatrixTest, GetFockType)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_TRUE(fmata.getFockType(0) == fockmat::restj);
    
    ASSERT_TRUE(fmata.getFockType(1) == fockmat::restk);
    
    ASSERT_TRUE(fmata.getFockType(2) == fockmat::restj);
    
    ASSERT_TRUE(fmata.getFockType(3) == fockmat::restk);
}

TEST_F(CAOFockMatrixTest, GetScaleFactor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_NEAR(1.0, fmata.getScaleFactor(0), 1.0e-13);
    
    ASSERT_NEAR(2.0, fmata.getScaleFactor(1), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmata.getScaleFactor(2), 1.0e-13);
    
    ASSERT_NEAR(1.0, fmata.getScaleFactor(3), 1.0e-13);
}

TEST_F(CAOFockMatrixTest, GetDensityIdentifier)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(0, fmata.getDensityIdentifier(0));
    
    ASSERT_EQ(2, fmata.getDensityIdentifier(1));
    
    ASSERT_EQ(3, fmata.getDensityIdentifier(2));
    
    ASSERT_EQ(6, fmata.getDensityIdentifier(3));
}

TEST_F(CAOFockMatrixTest, AddOneElectronProperty)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                         fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    CDenseMatrix mc({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
   
    fmata.addOneElectronMatrix(mc, 0);
    
    CDenseMatrix md({2.0, -2.0, -6.0, -4.0, 10.0, 8.0, 12.0, 8.0, -8.0}, 3, 3);
    
    CAOFockMatrix fmatb({md, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, IsSymmetric)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOFockMatrix fmata({ma, ma, mb, mb}, {fockmat::restj, fockmat::restk,
                        fockmat::restj, fockmat::restk}, {1.0, 2.0, 0.0, 1.0},
                        {0, 2, 3, 6});
    
    ASSERT_TRUE(fmata.isSymmetric(0));
    
    ASSERT_TRUE(fmata.isSymmetric(1));
    
    ASSERT_FALSE(fmata.isSymmetric(2));
    
    ASSERT_FALSE(fmata.isSymmetric(3));
}

TEST_F(CAOFockMatrixTest, GetElectronicEnergy)
{
    CDenseMatrix ma({ 1.0, -1.0, -3.0,
                     -2.0,  5.0,  4.0,
                      6.0,  4.0, -4.0},
                    3, 3);
    
    CDenseMatrix mb({ 1.0, -1.0, -3.0,
                     -2.0,  5.0,  4.0,
                      2.0,  1.0,  3.0},
                    3, 3);
    
    CAOFockMatrix fmat({ma, mb}, {fockmat::restj, fockmat::restk}, {1.0, 2.0}, {0, 2});
    
    CDenseMatrix da({ 1.0, -1.0, -3.0,
                     -2.0,  2.0,  3.0,
                      6.0,  3.0, -4.0},
                     3, 3);
    
    CDenseMatrix db({ 1.0, -1.0, -3.0,
                     -2.0,  1.0,  4.0,
                      3.0,  7.0,  3.0},
                    3, 3);
    
    CAODensityMatrix dmat({da, db}, denmat::unrest);
    
    ASSERT_NEAR(19.0, fmat.getElectronicEnergy(0 , dmat, 0), 1.0e-13);
    
    ASSERT_NEAR(15.0, fmat.getElectronicEnergy(0, dmat, 1), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmat.getElectronicEnergy(0, dmat, 2), 1.0e-13);
    
    ASSERT_NEAR(-6.0, fmat.getElectronicEnergy(1 , dmat, 0), 1.0e-13);
    
    ASSERT_NEAR(36.0, fmat.getElectronicEnergy(1, dmat, 1), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmat.getElectronicEnergy(1, dmat, 2), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmat.getElectronicEnergy(2 , dmat, 0), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmat.getElectronicEnergy(2, dmat, 1), 1.0e-13);
    
    ASSERT_NEAR(0.0, fmat.getElectronicEnergy(2, dmat, 2), 1.0e-13);
}


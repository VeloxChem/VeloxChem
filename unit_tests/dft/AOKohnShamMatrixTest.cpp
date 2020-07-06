//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AOKohnShamMatrixTest.hpp"

#include "AOKohnShamMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAOKohnShamMatrixTest, DefaultConstructor)
{
    CAOKohnShamMatrix fmata;
    
    CAOKohnShamMatrix fmatb(std::vector<CDenseMatrix>(), true, 0.0, 0.0);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 2.0);
    
    CAOKohnShamMatrix fmatb(fmata);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    CAOKohnShamMatrix fmatb(CAOKohnShamMatrix({ma, mb}, false, 3.0, 5.0));
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 2.0);
    
    CAOKohnShamMatrix fmatb = fmata;
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    CAOKohnShamMatrix fmatb = CAOKohnShamMatrix({ma, mb}, false, 3.0, 5.0);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, Zero)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 2.0);
    
    fmata.zero();
    
    ma.zero();
    
    mb.zero();
    
    CAOKohnShamMatrix fmatb({ma, mb}, true, 1.0, 2.0);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, Symmetrize)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 2.0, 3.0);
    
    fmata.symmetrize();
    
    CDenseMatrix mc({2.0, -3.0, 3.0, -3.0, 10.0, 8.0, 3.0, 8.0, -8.0}, 3, 3);
    
    CAOKohnShamMatrix fmatb({mc, mb}, false, 2.0, 3.0);
    
    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOKohnShamMatrixTest, SetNumberOfElectrons)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    fmata.setNumberOfElectrons(2.0);
    
    ASSERT_EQ(fmata, CAOKohnShamMatrix({ma, mb}, false, 2.0, 5.0));
}

TEST_F(CAOKohnShamMatrixTest, SetExchangeCorrelationEnergy)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    fmata.setExchangeCorrelationEnergy(1.0);
    
    ASSERT_EQ(fmata, CAOKohnShamMatrix({ma, mb}, false, 3.0, 1.0));
}

TEST_F(CAOKohnShamMatrixTest, IsRestricted)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 1.0, 2.0);
    
    ASSERT_FALSE(fmata.isRestricted());
}

TEST_F(CAOKohnShamMatrixTest, GetNumberOfElectrons)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    ASSERT_NEAR(fmata.getNumberOfElectrons(), 3.0, 1.0e-13);
}

TEST_F(CAOKohnShamMatrixTest, GetExchangeCorrelationEnergy)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, false, 3.0, 5.0);
    
    ASSERT_NEAR(fmata.getExchangeCorrelationEnergy(), 5.0, 1.0e-13);
}

TEST_F(CAOKohnShamMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 2.0);
    
    ASSERT_EQ(3, fmata.getNumberOfRows());
}

TEST_F(CAOKohnShamMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({mb, ma}, true, 2.0, 4.0);
    
    ASSERT_EQ(2, fmata.getNumberOfColumns());
}

TEST_F(CAOKohnShamMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 3.0);
    
    ASSERT_EQ(9, fmata.getNumberOfElements());
}

TEST_F(CAOKohnShamMatrixTest, GetNumberOfMatrices)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb}, true, 1.0, 2.0);
    
    ASSERT_EQ(2, fmata.getNumberOfMatrices());
}

TEST_F(CAOKohnShamMatrixTest, GetKohnShamConstant)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOKohnShamMatrix fmata({ma, mb}, false, 1.0, 3.0);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getKohnSham());
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getKohnSham(false));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getKohnSham(true));
}

TEST_F(CAOKohnShamMatrixTest, GetKohnSham)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb},false, 1.0, 2.0);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getKohnSham());
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getKohnSham(false));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getKohnSham(true));
    
    auto pdat = fmata.getKohnSham();
    
    pdat[2] = 2.0;  pdat[5] = 0.0;
    
    CDenseMatrix mc({1.0, -1.0, 2.0, -2.0, 5.0, 0.0, 6.0, 4.0, -4.0}, 3, 3);
    
    ASSERT_EQ(fmata, CAOKohnShamMatrix({mc, mb}, false, 1.0, 2.0));
}

TEST_F(CAOKohnShamMatrixTest, GetReferenceToKohnSham)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOKohnShamMatrix fmata({ma, mb}, false, 1.0, 2.0);
    
    ASSERT_EQ(ma, fmata.getReferenceToKohnSham());
    
    ASSERT_EQ(ma, fmata.getReferenceToKohnSham(false));
    
    ASSERT_EQ(mb, fmata.getReferenceToKohnSham(true));
}

TEST_F(CAOKohnShamMatrixTest, GetMatrixConstant)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOKohnShamMatrix fmata({ma, mb}, false, 1.0, 3.0);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getMatrix(0));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getMatrix(1));
}

TEST_F(CAOKohnShamMatrixTest, GetMatrix)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    CAOKohnShamMatrix fmata({ma, mb},false, 1.0, 2.0);
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0},
                     fmata.getMatrix(0));
    
    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, fmata.getMatrix(1));
    
    auto pdat = fmata.getMatrix(0);
    
    pdat[2] = 2.0;  pdat[5] = 0.0;
    
    CDenseMatrix mc({1.0, -1.0, 2.0, -2.0, 5.0, 0.0, 6.0, 4.0, -4.0}, 3, 3);
    
    ASSERT_EQ(fmata, CAOKohnShamMatrix({mc, mb}, false, 1.0, 2.0));
}

TEST_F(CAOKohnShamMatrixTest, GetReferenceToMatrix)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);
    
    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);
    
    const CAOKohnShamMatrix fmata({ma, mb}, false, 1.0, 2.0);
    
    ASSERT_EQ(ma, fmata.getReferenceToMatrix(0));
    
    ASSERT_EQ(mb, fmata.getReferenceToMatrix(1));
}

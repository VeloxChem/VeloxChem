//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

TEST_F(COverlapMatrixTest, GetOrthogonalizationMatrix)
{
    CDenseMatrix ma({10.0, 3.0, 4.0,
                     3.0, 1.2, 1.0,
                     4.0, 1.0, 4.0},
                     3, 3);
    
    COverlapMatrix smata(ma);
    
    CDenseMatrix refmat({ 0.5227351824567613, -0.4976329506165474, -0.1947637155486871,
                         -0.4976329506165478,  1.8097622833751084,  0.0808312367804720,
                         -0.1947637155486870,  0.0808312367804720,  0.6298490905403010},
                        3, 3);
    
    ASSERT_EQ(refmat, smata.getOrthogonalizationMatrix(0.0));
    
    CDenseMatrix submat({{ 0.2545280076864466, 0.2423193899561403,
                           0.1578026238399430, 0.0725559748124375,
                          -0.6206955232198551, 0.1178139554973027}},
                          3, 2);
    
    auto omat = smata.getOrthogonalizationMatrix(0.5);
    
    auto prefmat = submat.values();
    
    auto psolmat = omat.values();
    
    for (int32_t i = 0; i < omat.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(prefmat[i]), std::fabs(psolmat[i]), 1.0e-13);
    }
}



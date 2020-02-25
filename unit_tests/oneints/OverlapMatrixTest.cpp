//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "OverlapMatrixTest.hpp"

#include "CheckFunctions.hpp"
#include "OverlapMatrix.hpp"

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
    CDenseMatrix ma({10.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, 4.0}, 3, 3);

    COverlapMatrix smata(ma);

    CDenseMatrix refmat({0.659788099291142,
                         0.2545280076864466,
                         0.2423193899561403,
                         -1.870627094622050,
                         0.1578026238399430,
                         0.0725559748124375,
                         -0.205020511005960,
                         -0.6206955232198551,
                         0.1178139554973027},
                        3,
                        3);

    auto o3mat = smata.getOrthogonalizationMatrix(0.0);

    auto pref3mat = refmat.values();

    auto psol3mat = o3mat.values();

    for (int32_t i = 0; i < o3mat.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(pref3mat[i]), std::fabs(psol3mat[i]), 1.0e-13);
    }

    CDenseMatrix submat(
        {{0.2545280076864466, 0.2423193899561403, 0.1578026238399430, 0.0725559748124375, -0.6206955232198551, 0.1178139554973027}}, 3, 2);

    auto o2mat = smata.getOrthogonalizationMatrix(0.5);

    auto pref2mat = submat.values();

    auto psol2mat = o2mat.values();

    for (int32_t i = 0; i < o2mat.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(pref2mat[i]), std::fabs(psol2mat[i]), 1.0e-13);
    }
}

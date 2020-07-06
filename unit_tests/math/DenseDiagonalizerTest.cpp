//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "DenseDiagonalizerTest.hpp"

#include "DenseDiagonalizer.hpp"

TEST_F(CDenseDiagonalizerTest, DefaultConstructor)
{
    CDenseDiagonalizer diagdrv;

    ASSERT_TRUE(diagdrv.getState());

    ASSERT_EQ(CDenseMatrix(), diagdrv.getEigenVectors());

    ASSERT_EQ(CMemBlock<double>(), diagdrv.getEigenValues());
}

TEST_F(CDenseDiagonalizerTest, Diagonalize)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CMemBlock<double> evals({-4.583556800624598, -0.551057815834388, 6.334614616458989});

    ASSERT_EQ(evals, diagdrv.getEigenValues());

    CDenseMatrix refvecs({0.562733472929688,
                          0.356917927155759,
                          0.745614264696786,
                          -0.151384697403866,
                          -0.842233286839683,
                          0.517422229838614,
                          -0.812658422608437,
                          0.404045398209235,
                          0.419920950119791},
                         3,
                         3);

    auto evecs = diagdrv.getEigenVectors();

    auto prefvec = refvecs.values();

    auto psolvec = evecs.values();

    for (int32_t i = 0; i < refvecs.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(prefvec[i]), std::fabs(psolvec[i]), 1.0e-13);
    }
}

TEST_F(CDenseDiagonalizerTest, GetEigenVectors)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CDenseMatrix refvecs({0.562733472929688,
                          0.356917927155759,
                          0.745614264696786,
                          -0.151384697403866,
                          -0.842233286839683,
                          0.517422229838614,
                          -0.812658422608437,
                          0.404045398209235,
                          0.419920950119791},
                         3,
                         3);

    auto evecs = diagdrv.getEigenVectors();

    auto prefvec = refvecs.values();

    auto psolvec = evecs.values();

    for (int32_t i = 0; i < refvecs.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(prefvec[i]), std::fabs(psolvec[i]), 1.0e-13);
    }
}

TEST_F(CDenseDiagonalizerTest, GetEigenVectorsWithThreshold)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CDenseMatrix refveca({0.562733472929688,
                          0.356917927155759,
                          0.745614264696786,
                          -0.151384697403866,
                          -0.842233286839683,
                          0.517422229838614,
                          -0.812658422608437,
                          0.404045398209235,
                          0.419920950119791},
                         3,
                         3);

    auto evecs = diagdrv.getEigenVectors(-5.0);

    auto prefvec = refveca.values();

    auto psolvec = evecs.values();

    for (int32_t i = 0; i < refveca.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(prefvec[i]), std::fabs(psolvec[i]), 1.0e-13);
    }

    CDenseMatrix refvecb({0.745614264696786, 0.517422229838614, 0.419920950119791}, 3, 1);

    evecs = diagdrv.getEigenVectors(0.0);

    prefvec = refvecb.values();

    psolvec = evecs.values();

    for (int32_t i = 0; i < refvecb.getNumberOfElements(); i++)
    {
        ASSERT_NEAR(std::fabs(prefvec[i]), std::fabs(psolvec[i]), 1.0e-13);
    }
}

TEST_F(CDenseDiagonalizerTest, GetEigenValues)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CMemBlock<double> evals({-4.583556800624598, -0.551057815834388, 6.334614616458989});

    ASSERT_EQ(evals, diagdrv.getEigenValues());
}

TEST_F(CDenseDiagonalizerTest, GetEigenValuesWithThreshold)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CMemBlock<double> evala({-4.583556800624598, -0.551057815834388, 6.334614616458989});

    ASSERT_EQ(evala, diagdrv.getEigenValues(-5.0));

    CMemBlock<double> evalb({-0.551057815834388, 6.334614616458989});

    ASSERT_EQ(evalb, diagdrv.getEigenValues(-1.0));
}

TEST_F(CDenseDiagonalizerTest, GetState)
{
    CDenseMatrix mata({10.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, 4.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());
}

TEST_F(CDenseDiagonalizerTest, isLinearlyDependentBasis)
{
    CDenseMatrix mata({10.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, 4.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_FALSE(diagdrv.isLinearlyDependentBasis(0.0));

    ASSERT_TRUE(diagdrv.isLinearlyDependentBasis(0.5));
}

TEST_F(CDenseDiagonalizerTest, GetInvertedSqrtMatrix)
{
    CDenseMatrix mata({10.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, 4.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CDenseMatrix refmat({0.5227351824567613,
                         -0.4976329506165474,
                         -0.1947637155486871,
                         -0.4976329506165478,
                         1.8097622833751084,
                         0.0808312367804720,
                         -0.1947637155486870,
                         0.0808312367804720,
                         0.6298490905403010},
                        3,
                        3);

    ASSERT_EQ(refmat, diagdrv.getInvertedSqrtMatrix());
}

TEST_F(CDenseDiagonalizerTest, GetInvertedMatrix)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, -2.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    CDenseMatrix refmat({-0.212500000000000,
                         0.625000000000000,
                         -0.112500000000000,
                         0.625000000000000,
                         -1.250000000000000,
                         0.625000000000000,
                         -0.112500000000000,
                         0.625000000000000,
                         -0.412500000000000},
                        3,
                        3);

    ASSERT_EQ(refmat, diagdrv.getInvertedMatrix());
}

TEST_F(CDenseDiagonalizerTest, GetNumberOfEigenValues)
{
    CDenseMatrix mata({10.0, 3.0, 4.0, 3.0, 1.2, 1.0, 4.0, 1.0, 4.0}, 3, 3);

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(mata);

    ASSERT_TRUE(diagdrv.getState());

    ASSERT_EQ(3, diagdrv.getNumberOfEigenValues(0.0));

    ASSERT_EQ(2, diagdrv.getNumberOfEigenValues(0.5));

    ASSERT_EQ(1, diagdrv.getNumberOfEigenValues(3.0));

    ASSERT_EQ(0, diagdrv.getNumberOfEigenValues(20.0));
}

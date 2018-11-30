//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "BasisFunctionTest.hpp"

#include "BasisFunction.hpp"
#include "CheckFunctions.hpp"

TEST_F(CBasisFunctionTest, DefaultConstructor)
{
    CBasisFunction bfa;

    CBasisFunction bfb(std::vector<double>(), std::vector<double>(), 0, -1);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, CopyConstructor)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    CBasisFunction bfb(bfa);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, MoveConstructor)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    CBasisFunction bfb(CBasisFunction({0.3, 0.7}, {1.5, 5.6}, 1));

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, CopyAssignment)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    CBasisFunction bfb = bfa;

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, MoveAssignment)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    CBasisFunction bfb = CBasisFunction({0.3, 0.7}, {1.5, 5.6}, 1);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, SetExponents)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    bfa.setExponents({1.2, 3.2});

    CBasisFunction bfb({1.2, 3.2}, {1.5, 5.6}, 1);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, SetNormalizationFactors)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    bfa.setNormalizationFactors({2.0, 3.4});

    CBasisFunction bfb({0.3, 0.7}, {2.0, 3.4}, 1);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, SetAngularMomentum)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    bfa.setAngularMomentum(3);

    CBasisFunction bfb({0.3, 0.7}, {1.5, 5.6}, 3);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, Normalize)
{
    CBasisFunction bfa({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 7);

    bfa.normalize();

    CBasisFunction bfb({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 7);

    ASSERT_EQ(bfa, bfb);

    bfa.setAngularMomentum(0);

    bfa.normalize();

    bfb = CBasisFunction({0.800000000000000, 1.500000000000000,
                          2.700000000000000},
                         {0.542250667462678, 0.253418226705534,
                         -0.281296410762365},
                          0);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 1);

    bfa.normalize();

    bfb = CBasisFunction({0.800000000000000, 1.500000000000000,
                          2.700000000000000},
                         {0.958298883615608, 0.613252561079874,
                         -0.913275834565931},
                          1);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 2);

    bfa.normalize();

    bfb = CBasisFunction({0.800000000000000, 1.500000000000000,
                          2.700000000000000},
                         {0.490395036432180, 0.429719528168190,
                         -0.858586268491259},
                         2);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 3);

    bfa.normalize();

    bfb = CBasisFunction({0.800000000000000, 1.500000000000000,
                          2.70000000000000},
                         {0.389718603429863, 0.467617545540924,
                         -1.25350450414293},
                          3);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 4);

    bfa.normalize();

    bfb = CBasisFunction({0.8000000000000000, 1.500000000000000,
                          2.700000000000000},
                         {0.0655648873703003, 0.107723787859937,
                         -0.387420831892870},
                         4);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 5);

    bfa.normalize();

    bfb = CBasisFunction({0.8000000000000000, 1.5000000000000000,
                          2.700000000000000},
                         {0.0389717954232193, 0.0876781437334342,
                         -0.423057065413809},
                          5);

    ASSERT_EQ(bfa, bfb);

    bfa = CBasisFunction({0.8, 1.5, 2.7}, {2.4, 0.7, -0.5}, 6);

    bfa.normalize();

    bfb = CBasisFunction({0.8000000000000000, 1.5000000000000000,
                          2.700000000000000},
                         {0.0104896922997734, 0.0323150117007786,
                         -0.209193495104709},
                          6);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, NormalizeWithGeneralContraction)
{
    CBasisFunction bfa({1.5, 0.8}, {1.0, 0.8, 1.4, 1.2}, 2, 0);
    
    bfa.normalize();
    
    CBasisFunction bfb = CBasisFunction({1.5, 0.8},
                                        {0.546234791546083, 0.272720881091111,
                                         0.529489551834861, 0.283243303598867},
                                        2, 0);
    
    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, GetExponents)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    vlxtest::compare({0.3, 0.7}, bfa.getExponents());
}

TEST_F(CBasisFunctionTest, GetNormalizationFactors)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    vlxtest::compare({1.5, 5.6}, bfa.getNormalizationFactors());
}

TEST_F(CBasisFunctionTest, GetNormalizationFactorsWithGeneralContraction)
{
    CBasisFunction bfa({1.5, 0.8}, {1.0, 0.8, 1.4, 1.2}, 2, 0);
    
    vlxtest::compare({1.0, 0.8}, bfa.getNormalizationFactors(0));
    
    vlxtest::compare({1.4, 1.2}, bfa.getNormalizationFactors(1));
    
    auto vec = bfa.getNormalizationFactors(2);
    
    ASSERT_TRUE(vec.empty());
}

TEST_F(CBasisFunctionTest, GetAngularMomentum)
{
    CBasisFunction bfa({0.3, 0.7}, {1.5, 5.6}, 1);

    ASSERT_EQ(1, bfa.getAngularMomentum());
}

TEST_F(CBasisFunctionTest, GetNumberOfContractedFunctions)
{
    CBasisFunction bfa({1.5, 0.8}, {1.0, 0.8, 1.4, 1.2}, 2, 0);
    
    ASSERT_EQ(2, bfa.getNumberOfContractedFunctions());
    
    CBasisFunction bfb({0.3, 0.7}, {1.5, 5.6}, 1);
    
    ASSERT_EQ(1, bfb.getNumberOfContractedFunctions());
}

TEST_F(CBasisFunctionTest, Add)
{
    CBasisFunction bfa;

    bfa.setAngularMomentum(3);

    bfa.add(0.3, 1.5);

    bfa.add(0.7, 5.6);

    CBasisFunction bfb({0.3, 0.7}, {1.5, 5.6}, 3);

    ASSERT_EQ(bfa, bfb);
}

TEST_F(CBasisFunctionTest, GetNumberOfPrimitiveFunctions)
{
    CBasisFunction bfa({0.3, 0.7, 0.9}, {1.5, 5.6, 1.3}, 0);

    ASSERT_EQ(3, bfa.getNumberOfPrimitiveFunctions());
}

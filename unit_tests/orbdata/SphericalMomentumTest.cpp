//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "SphericalMomentumTest.hpp"

#include <cmath>

#include "CheckFunctions.hpp"
#include "SphericalMomentum.hpp"

TEST_F(CSphericalMomentumTest, DefaultConstructor)
{
    CSphericalMomentum ma;

    ASSERT_EQ(ma.getAngularMomentum(), -1);

    ASSERT_EQ(ma.getNumberOfComponents(), 0);
}

TEST_F(CSphericalMomentumTest, ConstructorForS)
{
    CSphericalMomentum ma(0);

    ASSERT_EQ(ma.getAngularMomentum(), 0);

    ASSERT_EQ(ma.getNumberOfComponents(), 1);

    vlxtest::compare(std::vector<double>({1.0}), ma.getFactors(0));

    ASSERT_EQ(1, ma.getNumberOfFactors(0));

    vlxtest::compare(std::vector<int32_t>({0}), ma.getIndexes(0));
}

TEST_F(CSphericalMomentumTest, ConstructorForP)
{
    CSphericalMomentum ma(1);

    ASSERT_EQ(ma.getAngularMomentum(), 1);

    ASSERT_EQ(ma.getNumberOfComponents(), 3);

    vlxtest::compare(std::vector<double>({1.0}), ma.getFactors(0));

    ASSERT_EQ(1, ma.getNumberOfFactors(0));

    vlxtest::compare(std::vector<int32_t>({1}), ma.getIndexes(0));

    vlxtest::compare(std::vector<double>({1.0}), ma.getFactors(1));

    ASSERT_EQ(1, ma.getNumberOfFactors(1));

    vlxtest::compare(std::vector<int32_t>({2}), ma.getIndexes(1));

    vlxtest::compare(std::vector<double>({1.0}), ma.getFactors(2));

    ASSERT_EQ(1, ma.getNumberOfFactors(2));

    vlxtest::compare(std::vector<int32_t>({0}), ma.getIndexes(2));
}

TEST_F(CSphericalMomentumTest, ConstructorForD)
{
    CSphericalMomentum ma(2);

    ASSERT_EQ(ma.getAngularMomentum(), 2);

    ASSERT_EQ(ma.getNumberOfComponents(), 5);

    double f3 = 2.0 * std::sqrt(3.0);

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(0));

    ASSERT_EQ(1, ma.getNumberOfFactors(0));

    vlxtest::compare(std::vector<int32_t>({1}), ma.getIndexes(0));

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(1));

    ASSERT_EQ(1, ma.getNumberOfFactors(1));

    vlxtest::compare(std::vector<int32_t>({4}), ma.getIndexes(1));

    vlxtest::compare({-1.0, -1.0, 2.0}, ma.getFactors(2));

    ASSERT_EQ(3, ma.getNumberOfFactors(2));

    vlxtest::compare({0, 3, 5}, ma.getIndexes(2));

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(3));

    ASSERT_EQ(1, ma.getNumberOfFactors(3));

    vlxtest::compare(std::vector<int32_t>({2}), ma.getIndexes(3));

    vlxtest::compare({0.5 * f3, -0.5 * f3}, ma.getFactors(4));

    ASSERT_EQ(2, ma.getNumberOfFactors(4));

    vlxtest::compare({0, 3}, ma.getIndexes(4));
}

TEST_F(CSphericalMomentumTest, ConstructorForF)
{
    CSphericalMomentum ma(3);

    ASSERT_EQ(ma.getAngularMomentum(), 3);

    ASSERT_EQ(ma.getNumberOfComponents(), 7);

    double f5 = std::sqrt(2.5);

    double f15 = std::sqrt(15.0);

    double f3 = std::sqrt(1.50);

    vlxtest::compare({3.0 * f5, -f5}, ma.getFactors(0));

    ASSERT_EQ(2, ma.getNumberOfFactors(0));

    vlxtest::compare({1, 6}, ma.getIndexes(0));

    vlxtest::compare(std::vector<double>({2.0 * f15}), ma.getFactors(1));

    ASSERT_EQ(1, ma.getNumberOfFactors(1));

    vlxtest::compare(std::vector<int32_t>({4}), ma.getIndexes(1));

    vlxtest::compare({4.0 * f3, -f3, -f3}, ma.getFactors(2));

    ASSERT_EQ(3, ma.getNumberOfFactors(2));

    vlxtest::compare({8, 1, 6}, ma.getIndexes(2));

    vlxtest::compare({2.0, -3.0, -3.0}, ma.getFactors(3));

    ASSERT_EQ(3, ma.getNumberOfFactors(3));

    vlxtest::compare({9, 2, 7}, ma.getIndexes(3));

    vlxtest::compare({4.0 * f3, -f3, -f3}, ma.getFactors(4));

    ASSERT_EQ(3, ma.getNumberOfFactors(4));

    vlxtest::compare({5, 0, 3}, ma.getIndexes(4));

    vlxtest::compare({f15, -f15}, ma.getFactors(5));

    ASSERT_EQ(2, ma.getNumberOfFactors(5));

    vlxtest::compare({2, 7}, ma.getIndexes(5));

    vlxtest::compare({f5, -3.0 * f5}, ma.getFactors(6));

    ASSERT_EQ(2, ma.getNumberOfFactors(6));

    vlxtest::compare({0, 3}, ma.getIndexes(6));
}

TEST_F(CSphericalMomentumTest, ConstructorForG)
{
    CSphericalMomentum ma(4);

    ASSERT_EQ(ma.getAngularMomentum(), 4);

    ASSERT_EQ(ma.getNumberOfComponents(), 9);

    double f35 = 4.0 * std::sqrt(35);

    double f17 = 4.0 * std::sqrt(17.5);

    double f5 = 4.0 * std::sqrt(5.0);

    double f2 = 4.0 * std::sqrt(2.5);

    vlxtest::compare({f35, -f35}, ma.getFactors(0));

    ASSERT_EQ(2, ma.getNumberOfFactors(0));

    vlxtest::compare({1, 6}, ma.getIndexes(0));

    vlxtest::compare({3.0 * f17, -f17}, ma.getFactors(1));

    ASSERT_EQ(2, ma.getNumberOfFactors(1));

    vlxtest::compare({4, 11}, ma.getIndexes(1));

    vlxtest::compare({6.0 * f5, -f5, -f5}, ma.getFactors(2));

    ASSERT_EQ(3, ma.getNumberOfFactors(2));

    vlxtest::compare({8, 1, 6}, ma.getIndexes(2));

    vlxtest::compare({4.0 * f2, -3.0 * f2, -3.0 * f2}, ma.getFactors(3));

    ASSERT_EQ(3, ma.getNumberOfFactors(3));

    vlxtest::compare({13, 4, 11}, ma.getIndexes(3));

    vlxtest::compare({8.0, 3.0, 3.0, 6.0, -24.0, -24.0}, ma.getFactors(4));

    ASSERT_EQ(6, ma.getNumberOfFactors(4));

    vlxtest::compare({14, 0, 10, 3, 5, 12}, ma.getIndexes(4));

    vlxtest::compare({4.0 * f2, -3.0 * f2, -3.0 * f2}, ma.getFactors(5));

    ASSERT_EQ(3, ma.getNumberOfFactors(5));

    vlxtest::compare({9, 2, 7}, ma.getIndexes(5));

    vlxtest::compare({3.0 * f5, -3.0 * f5, -0.5 * f5, 0.5 * f5}, ma.getFactors(6));

    ASSERT_EQ(4, ma.getNumberOfFactors(6));

    vlxtest::compare({5, 12, 0, 10}, ma.getIndexes(6));

    vlxtest::compare({f17, -3.0 * f17}, ma.getFactors(7));

    ASSERT_EQ(2, ma.getNumberOfFactors(7));

    vlxtest::compare({2, 7}, ma.getIndexes(7));

    vlxtest::compare({0.25 * f35, 0.25 * f35, -1.50 * f35}, ma.getFactors(8));

    ASSERT_EQ(3, ma.getNumberOfFactors(8));

    vlxtest::compare({0, 10, 3}, ma.getIndexes(8));
}

TEST_F(CSphericalMomentumTest, CopyConstructor)
{
    CSphericalMomentum ma(3);

    CSphericalMomentum mb(ma);

    ASSERT_EQ(ma, mb);
}

TEST_F(CSphericalMomentumTest, MoveConstructor)
{
    CSphericalMomentum ma(3);

    CSphericalMomentum mb(CSphericalMomentum(3));

    ASSERT_EQ(ma, mb);
}

TEST_F(CSphericalMomentumTest, CopyAssignment)
{
    CSphericalMomentum ma(3);

    CSphericalMomentum mb(ma);

    ASSERT_EQ(ma, mb);
}

TEST_F(CSphericalMomentumTest, MoveAssignment)
{
    CSphericalMomentum ma(3);

    CSphericalMomentum mb = CSphericalMomentum(3);

    ASSERT_EQ(ma, mb);
}

TEST_F(CSphericalMomentumTest, GetAngularMomentum)
{
    CSphericalMomentum ma(2);

    ASSERT_EQ(ma.getAngularMomentum(), 2);
}

TEST_F(CSphericalMomentumTest, GetNumberOfComponents)
{
    CSphericalMomentum ma(2);

    ASSERT_EQ(ma.getNumberOfComponents(), 5);
}

TEST_F(CSphericalMomentumTest, GetFactors)
{
    CSphericalMomentum ma(2);

    double f3 = 2.0 * std::sqrt(3.0);

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(0));

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(1));

    vlxtest::compare({-1.0, -1.0, 2.0}, ma.getFactors(2));

    vlxtest::compare(std::vector<double>({f3}), ma.getFactors(3));

    vlxtest::compare({0.5 * f3, -0.5 * f3}, ma.getFactors(4));
}

TEST_F(CSphericalMomentumTest, getNumberOfFactors)
{
    CSphericalMomentum ma(2);

    ASSERT_EQ(1, ma.getNumberOfFactors(0));

    ASSERT_EQ(1, ma.getNumberOfFactors(1));

    ASSERT_EQ(3, ma.getNumberOfFactors(2));

    ASSERT_EQ(1, ma.getNumberOfFactors(3));

    ASSERT_EQ(2, ma.getNumberOfFactors(4));
}

TEST_F(CSphericalMomentumTest, GetIndexes)
{
    CSphericalMomentum ma(2);

    vlxtest::compare(std::vector<int32_t>({1}), ma.getIndexes(0));

    vlxtest::compare(std::vector<int32_t>({4}), ma.getIndexes(1));

    vlxtest::compare({0, 3, 5}, ma.getIndexes(2));

    vlxtest::compare(std::vector<int32_t>({2}), ma.getIndexes(3));

    vlxtest::compare({0, 3}, ma.getIndexes(4));
}

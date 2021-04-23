//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "MathFuncTest.hpp"

#include "CheckFunctions.hpp"
#include "MathFunc.hpp"

TEST_F(CMathFuncTest, ZeroReal)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {0.0, 0.0, 0.0, 0.0, 0.0};

    mathfunc::zero(veca, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, ZeroInteger)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, -1, 3, 4};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN))) = {0, 0, 0, 0, 0};

    mathfunc::zero(veca, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, Set_ToReal)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {2.0, 2.0, 2.0, 2.0, 2.0};

    mathfunc::set_to(veca, 2.0, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, Set_ToInteger)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, -1, 3, 4};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN))) = {2, 2, 2, 2, 2};

    mathfunc::set_to(veca, 2, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, SumReal)
{
    double vec[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    ASSERT_NEAR(mathfunc::sum(vec, 5), 19.0, 1.0e-13);
}

TEST_F(CMathFuncTest, SumInteger)
{
    int32_t vec[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, -1, 3, 4};

    ASSERT_EQ(mathfunc::sum(vec, 5), 19);
}

TEST_F(CMathFuncTest, Scale)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {10.0, 16.0, -2.0, 6.0, 8.0};

    mathfunc::scale(veca, 2.0, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, Add_Scaled)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {1.0, 3.0, 2.0, 4.0, 5.0};

    double vecc[5] __attribute__((aligned(VLX_ALIGN))) = {7.0, 14.0, 3.0, 11.0, 14.0};

    mathfunc::add_scaled(veca, vecb, 2.0, 5);

    vlxtest::compare(veca, vecc, 5);
}

TEST_F(CMathFuncTest, MaxReal)
{
    double vec[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    ASSERT_NEAR(mathfunc::max(vec, 5), 8.0, 1.0e-13);
}

TEST_F(CMathFuncTest, MaxInteger)
{
    int32_t vec[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, -1, 3, 41};

    ASSERT_EQ(mathfunc::max(vec, 5), 41);
}

TEST_F(CMathFuncTest, Normalize)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 9.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {0.25, 0.45, -0.05, 0.15, 0.20};

    mathfunc::normalize(veca, 5);

    vlxtest::compare(veca, vecb, 5);
}

TEST_F(CMathFuncTest, Distance)
{
    ASSERT_NEAR(mathfunc::distance(4.0, 5.0, 2.0, 3.0, 3.0, 0.0), 3.0, 1.0e-13);
}

TEST_F(CMathFuncTest, Distances)
{
    double rbx[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double rby[5] __attribute__((aligned(VLX_ALIGN))) = {0.0, 1.0, -2.0, 7.0, 3.0};

    double rbz[5] __attribute__((aligned(VLX_ALIGN))) = {1.0, -2.0, 0.0, 5.0, 2.0};

    double abx[4] __attribute__((aligned(VLX_ALIGN)));

    double aby[4] __attribute__((aligned(VLX_ALIGN)));

    double abz[4] __attribute__((aligned(VLX_ALIGN)));

    mathfunc::distances(abx, aby, abz, 6.0, 4.0, -1.0, rbx, rby, rbz, 4);

    double rdx[4] __attribute__((aligned(VLX_ALIGN))) = {1.0, -2.0, 7.0, 3.0};

    vlxtest::compare(abx, rdx, 4);

    double rdy[4] __attribute__((aligned(VLX_ALIGN))) = {4.0, 3.0, 6.0, -3.0};

    vlxtest::compare(aby, rdy, 4);

    double rdz[4] __attribute__((aligned(VLX_ALIGN))) = {-2.0, 1.0, -1.0, -6.0};

    vlxtest::compare(abz, rdz, 4);
}

TEST_F(CMathFuncTest, Indexes)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, 4, 3, 4};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN)));

    int32_t vecc[5] __attribute__((aligned(VLX_ALIGN))) = {0, 5, 13, 17, 20};

    mathfunc::indexes(vecb, veca, 5);

    vlxtest::compare(vecb, vecc, 5);
}

TEST_F(CMathFuncTest, IndexesWithOffset)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, 4, 3, 4};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN)));

    int32_t vecc[5] __attribute__((aligned(VLX_ALIGN))) = {0, 5, 13, 17, 20};

    int32_t vecd[5] __attribute__((aligned(VLX_ALIGN))) = {2, 7, 15, 19, 22};

    mathfunc::indexes(vecb, veca, 0, 5);

    vlxtest::compare(vecb, vecc, 5);

    mathfunc::indexes(vecb, veca, 2, 5);

    vlxtest::compare(vecb, vecd, 5);
}

TEST_F(CMathFuncTest, Ordering)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {0, 1, 0, 1, 0};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN))) = {0, 0, 0, 0, 0};
    ;

    int32_t vecc[5] __attribute__((aligned(VLX_ALIGN))) = {1, 3, 0, 0, 0};

    mathfunc::ordering(vecb, veca, 5);

    vlxtest::compare(vecb, vecc, 5);
}

TEST_F(CMathFuncTest, QuadChebyshevOfKindTwo)
{
    double coords[5] __attribute__((aligned(VLX_ALIGN)));

    double weights[5] __attribute__((aligned(VLX_ALIGN)));

    mathfunc::quadChebyshevOfKindTwo(coords, weights, 5);

    double refcrd[5] __attribute__((aligned(VLX_ALIGN))) = {
        8.66025403784439e-01, 5.00000000000000e-01, 6.12323399573677e-17, -5.00000000000000e-01, -8.66025403784439e-01};

    vlxtest::compare(coords, refcrd, 5);

    double refwgt[5] __attribute__((aligned(VLX_ALIGN))) = {
        1.30899693899575e-01, 3.92699081698724e-01, 5.23598775598299e-01, 3.92699081698724e-01, 1.30899693899575e-01};

    vlxtest::compare(weights, refwgt, 5);
}

TEST_F(CMathFuncTest, CopyReal)
{
    double veca[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 8.0, -1.0, 3.0, 4.0};

    double vecb[5] __attribute__((aligned(VLX_ALIGN))) = {3.0, 1.0, 4.0, 5.0, 6.0};

    mathfunc::copy(veca, 1, vecb, 2, 2);

    double vecc[5] __attribute__((aligned(VLX_ALIGN))) = {5.0, 4.0, 5.0, 3.0, 4.0};

    vlxtest::compare(veca, vecc, 5);
}

TEST_F(CMathFuncTest, CopyInteger)
{
    int32_t veca[5] __attribute__((aligned(VLX_ALIGN))) = {5, 8, -1, 3, 4};

    int32_t vecb[5] __attribute__((aligned(VLX_ALIGN))) = {1, 2, 4, 6, 3};

    mathfunc::copy(veca, 0, vecb, 3, 2);

    int32_t vecc[5] __attribute__((aligned(VLX_ALIGN))) = {6, 3, -1, 3, 4};

    vlxtest::compare(veca, vecc, 5);
}

TEST_F(CMathFuncTest, MaxTensorComponents)
{
    ASSERT_EQ(0, mathfunc::maxTensorComponents(-1));

    ASSERT_EQ(1, mathfunc::maxTensorComponents(0));

    ASSERT_EQ(3, mathfunc::maxTensorComponents(1));

    ASSERT_EQ(9, mathfunc::maxTensorComponents(2));

    ASSERT_EQ(27, mathfunc::maxTensorComponents(3));

    ASSERT_EQ(81, mathfunc::maxTensorComponents(4));
}

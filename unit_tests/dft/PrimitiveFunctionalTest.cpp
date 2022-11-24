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

#include "PrimitiveFunctionalTest.hpp"

#include "DummyFunctionals.hpp"
#include "PrimitiveFunctional.hpp"

TEST_F(CPrimitiveFunctionalTest, DefaultConstructor)
{
    CPrimitiveFunctional rfa({}, xcfun::undefined, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

    CPrimitiveFunctional rfb;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, CopyConstructor)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    CPrimitiveFunctional rfb(rfa);

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, MoveConstructor)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    CPrimitiveFunctional rfb(CPrimitiveFunctional({"Slater"},
                                                  xcfun::lda,
                                                  &vlxtest::dummy_fvxc_ab,
                                                  &vlxtest::dummy_fvxc_a,
                                                  &vlxtest::dummy_fvxc_b,
                                                  &vlxtest::dummy_fvxc2_ab,
                                                  &vlxtest::dummy_fvxc2_a,
                                                  &vlxtest::dummy_fvxc2_b,
                                                  &vlxtest::dummy_fvxc3_ab,
                                                  &vlxtest::dummy_fvxc3_a,
                                                  &vlxtest::dummy_fvxc3_b));

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, CopyAssignment)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    CPrimitiveFunctional rfb = rfa;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, MoveAssignment)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    CPrimitiveFunctional rfb = CPrimitiveFunctional({"Slater"},
                                                    xcfun::lda,
                                                    &vlxtest::dummy_fvxc_ab,
                                                    &vlxtest::dummy_fvxc_a,
                                                    &vlxtest::dummy_fvxc_b,
                                                    &vlxtest::dummy_fvxc2_ab,
                                                    &vlxtest::dummy_fvxc2_a,
                                                    &vlxtest::dummy_fvxc2_b,
                                                    &vlxtest::dummy_fvxc3_ab,
                                                    &vlxtest::dummy_fvxc3_a,
                                                    &vlxtest::dummy_fvxc3_b);

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CPrimitiveFunctionalTest, GetLabel)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    ASSERT_EQ(std::string("Slater"), rfa.getLabel());
}

TEST_F(CPrimitiveFunctionalTest, GetFunctionalType)
{
    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    ASSERT_EQ(xcfun::lda, rfa.getFunctionalType());
}

TEST_F(CPrimitiveFunctionalTest, ComputeWithGradient)
{
    CDensityGrid dgrid(2, 1, xcfun::lda, dengrid::ab);

    CXCGradientGrid xcgrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 2, 3), dengrid::ab, xcfun::lda);

    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    rfa.compute(xcgrid, 1.0, dgrid);

    ASSERT_EQ(xcgrid, CXCGradientGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 3), dengrid::ab, xcfun::lda));
}

TEST_F(CPrimitiveFunctionalTest, ComputeWithHessian)
{
    CDensityGrid dgrid(2, 1, xcfun::lda, dengrid::ab);

    CXCHessianGrid xcgrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 2, 3), dengrid::ab, xcfun::lda);

    CPrimitiveFunctional rfa({"Slater"},
                             xcfun::lda,
                             &vlxtest::dummy_fvxc_ab,
                             &vlxtest::dummy_fvxc_a,
                             &vlxtest::dummy_fvxc_b,
                             &vlxtest::dummy_fvxc2_ab,
                             &vlxtest::dummy_fvxc2_a,
                             &vlxtest::dummy_fvxc2_b,
                             &vlxtest::dummy_fvxc3_ab,
                             &vlxtest::dummy_fvxc3_a,
                             &vlxtest::dummy_fvxc3_b);

    rfa.compute(xcgrid, 1.0, dgrid);

    ASSERT_EQ(xcgrid, CXCHessianGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 3), dengrid::ab, xcfun::lda));
}

TEST(LibXCFunctionalTest, Constructor)
{
    std::cout << getLibXCDescription() << std::endl;

    auto fun = Functional({"LDA_X", "GGA_C_LYP", "MGGA_C_SCAN", "GGA_X_B88"}, {0.5, 0.0323849, 1.3, 0.12});

    std::cout << fun << std::endl;
}

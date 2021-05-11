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

#include "DenseLinearAlgebraTest.hpp"

#include "DenseLinearAlgebra.hpp"

TEST_F(CDenseLinearAlgebraTest, MultAB)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, 4.0, -3.0, 3.0, 1.0, -1.0, 2.0, 3.0, 0.5, 7.8, 1.0, 2.1}, 3, 5);

    auto matab = denblas::multAB(mata, matb);

    CDenseMatrix refab({5.0, 15.0, 44.2, 5.0, 22.4, -9.0, 3.5, -4.2, -8.0, -3.9, 20.1, 22.4, 86.9, 16.7, 43.3, 7.0, 10.0, 38.2, 4.0, 16.4}, 4, 5);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultAtB)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0, 1.0, -1.0}, 4, 2);

    auto matab = denblas::multAtB(mata, matb);

    CDenseMatrix refab({-30.0, 15.0, 13.1, 16.9, -8.0, 27.0}, 3, 2);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultABt)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0}, 2, 3);

    auto matab = denblas::multABt(mata, matb);

    CDenseMatrix refab({28.0, 7.0, 8.0, -12.0, 45.6, 26.1, 25.0, 8.0}, 4, 2);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultDiagByA)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CMemBlock<double> diagb({1.0, 2.0, -2.0, 5.0});

    auto matab = denblas::multDiagByA(diagb, mata);

    CDenseMatrix refab({2.0, 3.0, 4.0, -6.0, 6.0, 2.0, -12.0, -4.6, -14.0, 5.0, 10.0, 20.0}, 4, 3);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultDiagByAt)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CMemBlock<double> diagb({1.0, 2.0, -2.0});

    auto matab = denblas::multDiagByAt(diagb, mata);

    CDenseMatrix refab({2.0, -3.0, 6.0, 1.0, 6.0, 6.0, 4.6, 4.0, -8.0, -2.0, -14.0, -8.0}, 3, 4);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, SubAB)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({0.8, 1.2, 3.2, 2.1, 2.5, 1.7, 4.9, 2.3, 8.1, 0.2, 1.3, 5.1}, 4, 3);

    auto matab = denblas::subAB(mata, matb);

    CDenseMatrix refab({1.2, 1.8, 0.8, -5.1, 0.5, -0.7, 1.1, 0.0, -1.1, 0.8, 0.7, -1.1}, 4, 3);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, AddAB)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({0.8, 1.2, 3.2, 2.1, 2.5, 1.7, 4.9, 2.3, 8.1, 0.2, 1.3, 5.1}, 4, 3);

    auto matab = denblas::addAB(mata, matb, 2.0);

    CDenseMatrix refab({3.6, 5.4, 10.4, 1.2, 8.0, 4.4, 15.8, 6.9, 23.2, 1.4, 4.6, 14.2}, 4, 3);

    ASSERT_EQ(matab, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultABtWithAddition)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0}, 2, 3);

    CDenseMatrix matc({1.0, 4.0, 2.0, 6.0, 8.0, 1.0, 2.0, 5.0}, 4, 2);

    denblas::multABt(matc, 2.0, 1.0, mata, matb);

    CDenseMatrix refc({57.0, 18.0, 18.0, -18.0, 99.2, 53.2, 52.0, 21.0}, 4, 2);

    ASSERT_EQ(matc, refc);
}

TEST_F(CDenseLinearAlgebraTest, MultABtWithAdditionForPointerForm)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);
    
    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0}, 2, 3);
    
    CDenseMatrix matc({1.0, 4.0, 2.0, 6.0, 8.0, 1.0, 2.0, 5.0}, 4, 2);
    
    denblas::multABt(matc.values(), 2.0, 1.0, mata, matb);
    
    CDenseMatrix refc({57.0, 18.0, 18.0, -18.0, 99.2, 53.2, 52.0, 21.0}, 4, 2);
    
    ASSERT_EQ(matc, refc);
}

TEST_F(CDenseLinearAlgebraTest, MultAtBWithAdditionForPointerForm)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);
    
    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0, 1.0, -1.0}, 4, 2);
    
    CDenseMatrix matc({2.0, 1.0, 3.0, 4.0, 5.0, 6.0}, 3, 2);
    
    denblas::multAtB(matc.values(), 1.0, 1.0, mata, matb);
    
    CDenseMatrix refab({-28.0, 16.0, 16.1, 20.9, -3.0, 33.0}, 3, 2);
    
    ASSERT_EQ(matc, refab);
}

TEST_F(CDenseLinearAlgebraTest, MultABWithAdditionForPointerForm)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);
    
    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, 4.0, -3.0, 3.0, 1.0, -1.0, 2.0, 3.0, 0.5, 7.8, 1.0, 2.1}, 3, 5);
    
    CDenseMatrix matc({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, 4, 5);
    
    denblas::multAB(matc.values(), 1.0, 1.0, mata, matb);
    
    CDenseMatrix refab({6.0, 16.0, 45.2, 6.0, 23.4, -8.0, 4.5, -3.2, -7.0, -2.9, 21.1, 23.4, 87.9, 17.7, 44.3, 8.0, 11.0, 39.2, 5.0, 17.4}, 4, 5);
    
    ASSERT_EQ(matc, refab);
}

TEST_F(CDenseLinearAlgebraTest, Dot)
{
    CMemBlock<double> veca({2.0, 3.0, 4.0});

    CMemBlock<double> vecb({1.0, 2.0, -1.5});

    ASSERT_NEAR(2.0, denblas::dot(veca, vecb), 1.0e-13);

    ASSERT_NEAR(2.0, denblas::dot(vecb, veca), 1.0e-13);

    ASSERT_NEAR(29.0, denblas::dot(veca, veca), 1.0e-13);

    ASSERT_NEAR(7.25, denblas::dot(vecb, vecb), 1.0e-13);
}

TEST_F(CDenseLinearAlgebraTest, DotVectorMatrix)
{
    CMemBlock<double> veca({2.0, 3.0, 4.0});

    CDenseMatrix matb({1.0, 2.0, -1.5}, 3, 1);

    ASSERT_NEAR(2.0, denblas::dot(veca, matb), 1.0e-13);
}

TEST_F(CDenseLinearAlgebraTest, Trace)
{
    CDenseMatrix mat({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0}, 3, 3);

    ASSERT_NEAR(12.0, denblas::trace(mat), 1.0e-13);
}

TEST_F(CDenseLinearAlgebraTest, TraceAB)
{
    CDenseMatrix mata({2.0, 3.0, 4.0, -3.0, 3.0, 1.0, 6.0, 2.3, 7.0, 1.0, 2.0, 4.0}, 4, 3);

    CDenseMatrix matb({1.0, 2.0, 5.0, 2.0, -3.0, 3.0, 1.0, -1.0, 3.0, 0.5, 7.8, 1.0}, 3, 4);

    ASSERT_NEAR(99.40, denblas::trace(mata, matb), 1.0e-13);
}

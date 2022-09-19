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

#include <gtest/gtest.h>

#include <cstddef>

#include "IntegerSequence.hpp"

using namespace buffer;

TEST(IntegerSequenceTest, seq)
{
    auto start = 2;
    auto end   = 5;
    auto N     = end - start + 1;

    auto x = seq(start, end);

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], start + i);
    }

    ASSERT_EQ(x[3], end);
}

TEST(IntegerSequenceTest, seqWithIncrement)
{
    auto start = 2;
    auto end   = 8;
    auto step  = 2;
    auto N     = 1 + (end - start) / step;

    auto x = seq(start, end, step);

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], start + i * step);
    }

    ASSERT_EQ(x[3], end);
}

TEST(IntegerSequenceTest, seqN)
{
    auto start = 2;
    auto N     = 5;
    auto end   = start + N - 1;

    auto x = seqN(start, N);

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], start + i);
    }

    ASSERT_EQ(x[4], end);
}

TEST(IntegerSequenceTest, seqNWithIncrement)
{
    auto start = 2;
    auto N     = 3;
    auto step  = 3;
    auto end   = start + (N - 1) * step;

    auto x = seqN(start, N, step);

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], start + i * step);
    }

    ASSERT_EQ(x[2], end);
}

TEST(IntegerSequenceTest, seqReverse)
{
    auto start = 2;
    auto end   = 5;
    auto N     = end - start + 1;

    auto x = seq(start, end).reverse();

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], end - i);
    }

    ASSERT_EQ(x[3], start);
}

TEST(IntegerSequenceTest, seqWithIncrementReverse)
{
    auto start = 2;
    auto end   = 8;
    auto step  = 2;
    auto N     = 1 + (end - start) / step;

    auto x = seq(start, end, step).reverse();

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], end - i * step);
    }

    ASSERT_EQ(x[3], start);
}

TEST(IntegerSequenceTest, seqNReverse)
{
    auto start = 2;
    auto N     = 5;
    auto end   = start + N - 1;

    auto x = seqN(start, N).reverse();

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], end - i);
    }

    ASSERT_EQ(x[4], start);
}

TEST(IntegerSequenceTest, seqNWithIncrementReverse)
{
    auto start = 2;
    auto N     = 3;
    auto step  = 3;
    auto end   = (start + (N - 1) * step);

    auto x = seqN(start, N, step).reverse();

    ASSERT_EQ(x.size(), N);

    for (auto i = 0; i < x.size(); ++i)
    {
        ASSERT_EQ(x[i], end - i * step);
    }

    ASSERT_EQ(x[2], start);
}

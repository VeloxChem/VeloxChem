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

#include "CBufferTest.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "Buffer.hpp"
#include "MemAlloc.hpp"
#include "MpiFunc.hpp"

using namespace buffer;
using namespace detail;

TYPED_TEST_SUITE(CBufferTest, detail::implementations);

TYPED_TEST(CBufferTest, DefaultConstructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    ASSERT_EQ(buf.nRows(), NRows);
    ASSERT_EQ(buf.nColumns(), NCols);
}

TYPED_TEST(CBufferTest, ConstructorWithDimensions)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
    }
}

TYPED_TEST(CBufferTest, ConstructorWithDimensionsAndRawPointer)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto data = new Scalar[100];

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), NRows * 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.size(), 10 * NCols);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_EQ(buf.size(), 50);
    }
    else if constexpr (kind == Kind::N)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_EQ(buf.size(), NRows * NCols);
    }

    delete[] data;
}

TYPED_TEST(CBufferTest, ConstructorWithDimensionsAndStdVector)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto data = std::vector<Scalar>(100);

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), NRows * 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.size(), 10 * NCols);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_EQ(buf.size(), 50);
    }
    else if constexpr (kind == Kind::N)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_EQ(buf.size(), NRows * NCols);
    }
}

TYPED_TEST(CBufferTest, ConstructorWithStdVector)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto data = std::vector<Scalar>(10);

        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else
    {
        // do nothing: this CTOR is only valid for Kind::X
    }
}

TYPED_TEST(CBufferTest, CopyConstructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> src{};
    CBuffer<Scalar, Backend, NRows, NCols> dst{src};

    // the two objects are exactly equal, but independent from each other
    ASSERT_EQ(dst, src);
}

TYPED_TEST(CBufferTest, CopyAssignment)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> src{};
    CBuffer<Scalar, Backend, NRows, NCols> dst{};

    dst = src;

    // the two objects are exactly equal, but independent from each other
    ASSERT_EQ(dst, src);
}

TYPED_TEST(CBufferTest, MoveConstructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> src{};

    auto nrows = src.nRows();
    auto ncols = src.nColumns();

    CBuffer<Scalar, Backend, NRows, NCols> dst{std::move(src)};

    // dst "steals" representation of src
    ASSERT_EQ(dst.nRows(), nrows);
    ASSERT_EQ(dst.nColumns(), ncols);

    // src is still in a valid state
    ASSERT_EQ(src.nRows(), NRows);
    ASSERT_EQ(src.nColumns(), NCols);
    ASSERT_EQ(src.data(), nullptr);
}

TYPED_TEST(CBufferTest, MoveAssignment)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> src{};

    auto nrows = src.nRows();
    auto ncols = src.nColumns();

    CBuffer<Scalar, Backend, NRows, NCols> dst{};

    dst = std::move(src);

    // dst "steals" representation of src
    ASSERT_EQ(dst.nRows(), nrows);
    ASSERT_EQ(dst.nColumns(), ncols);

    // src is still in a valid state
    ASSERT_EQ(src.nRows(), NRows);
    ASSERT_EQ(src.nColumns(), NCols);
    ASSERT_EQ(src.data(), nullptr);
}

TYPED_TEST(CBufferTest, DefaultConstructorAndResize)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.resize(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else
    {
        buf.resize();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorAndSetZero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.setZero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.setZero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        buf.setZero();

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{0}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, ConstructorWithDimensionsAndSetZero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setZero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{0}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorResizeAndSetZero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.resize(10, 5);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.resize();
        buf.setZero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        buf.resize();
        buf.setZero();

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{0}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorAndSetConstant)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.setConstant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{42}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, ConstructorWithDimensionsAndSetConstant)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{42}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorResizeAndSetConstant)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.resize(10, 5);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.resize();
        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        buf.resize();
        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{42}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorAndSetRandom)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == Kind::X)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.setRandom(LB, UB, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4));
    }
}

TYPED_TEST(CBufferTest, ConstructorWithDimensionsAndSetRandom)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4));
    }
}

TYPED_TEST(CBufferTest, DefaultConstructorResizeAndSetRandom)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == Kind::X)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == Kind::MY)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == Kind::XN)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == Kind::XY)
    {
        buf.resize(10, 5);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == Kind::MN)
    {
        buf.resize();

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        buf.resize();

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4));
    }
}

TYPED_TEST(CBufferTest, Zero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Zero();

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{0}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, Constant)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4), Scalar{42}, 1.0e-14);
    }
}

TYPED_TEST(CBufferTest, Random)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB);

        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4));
    }
}

TYPED_TEST(CBufferTest, Broadcast)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf, ref);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf, ref);
    }
}

TYPED_TEST(CBufferTest, ReduceSum)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto sz   = mpi::nodes(MPI_COMM_WORLD);
    auto rank = mpi::rank(MPI_COMM_WORLD);

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42}, 10);

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42}, 10);

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42}, 10);

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42}, 10, 5);

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42});

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        buf.reduce_sum(0, sz, MPI_COMM_WORLD);

        auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(sz * Scalar{42});

        if (rank == mpi::master())
        {
            ASSERT_EQ(buf, ref);
        }
    }
}

TYPED_TEST(CBufferTest, Scatter)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto sz   = mpi::nodes(MPI_COMM_WORLD);
    auto rank = mpi::rank(MPI_COMM_WORLD);

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a chunk_sz-sized piece of the data at the head
        auto                chunk_sz = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10);
        std::fill_n(foo.begin(), chunk_sz, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        for (auto i = 0; i < chunk_sz; ++i)
        {
            ASSERT_NEAR(buf(i), ref(i), 1.0e-14);
        }
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a piece of the data as first row
        auto                rows_in_chunk = mpi::batch_size(NRows, rank, sz);
        std::vector<Scalar> foo(NRows * 10);
        std::fill_n(foo.begin(), rows_in_chunk * 10, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10);

        for (auto i = 0; i < rows_in_chunk; ++i)
        {
            for (decltype(buf.nColumns()) j = 0; j < buf.nColumns(); ++j)
            {
                ASSERT_NEAR(buf(i, j), ref(i, j), 1.0e-14);
            }
        }
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a piece of the data as first row
        auto                rows_in_chunk = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10 * NCols);
        std::fill_n(foo.begin(), rows_in_chunk * NCols, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10);

        for (auto i = 0; i < rows_in_chunk; ++i)
        {
            for (decltype(buf.nColumns()) j = 0; j < buf.nColumns(); ++j)
            {
                ASSERT_NEAR(buf(i, j), ref(i, j), 1.0e-14);
            }
        }
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a piece of the data as first row
        auto                rows_in_chunk = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10 * 5);
        std::fill_n(foo.begin(), rows_in_chunk * 5, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10, 5);

        for (auto i = 0; i < rows_in_chunk; ++i)
        {
            for (decltype(buf.nColumns()) j = 0; j < buf.nColumns(); ++j)
            {
                ASSERT_NEAR(buf(i, j), ref(i, j), 1.0e-14);
            }
        }
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a piece of the data as first row
        auto                rows_in_chunk = mpi::batch_size(NRows, rank, sz);
        std::vector<Scalar> foo(NRows * NCols);
        std::fill_n(foo.begin(), rows_in_chunk * NCols, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        for (auto i = 0; i < rows_in_chunk; ++i)
        {
            for (decltype(buf.nColumns()) j = 0; j < buf.nColumns(); ++j)
            {
                ASSERT_NEAR(buf(i, j), ref(i, j), 1.0e-14);
            }
        }
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        buf.scatter(rank, sz, MPI_COMM_WORLD);

        // after scattering, each process has a chunk_sz-sized piece of the data at the head
        auto                chunk_sz = mpi::batch_size(NCols, rank, sz);
        std::vector<Scalar> foo(NCols);
        std::fill_n(foo.begin(), chunk_sz, Scalar{42});
        auto ref = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        for (auto i = 0; i < chunk_sz; ++i)
        {
            ASSERT_NEAR(buf(i), ref(i), 1.0e-14);
        }
    }
}

TYPED_TEST(CBufferTest, Gather)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    auto sz   = mpi::nodes(MPI_COMM_WORLD);
    auto rank = mpi::rank(MPI_COMM_WORLD);

    if constexpr (kind == Kind::X)
    {
        // on each rank, generate a buffer which contains a chunk_sz part of the data
        auto                chunk_sz = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10);
        std::fill_n(foo.begin(), chunk_sz, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        // after gathering, we have a buffer with 10 elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
            ASSERT_EQ(buf_1, ref);
        }
    }
    else if constexpr (kind == Kind::MY)
    {
        // on each rank, generate a buffer which contains rows_in_chunk rows of the data
        auto                rows_in_chunk = mpi::batch_size(NRows, rank, sz);
        std::vector<Scalar> foo(NRows * 10);
        std::fill_n(foo.begin(), rows_in_chunk * 10, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10);

        // after gathering, we have a buffer with NRows*10 elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
            ASSERT_EQ(buf_1, ref);
        }
    }
    else if constexpr (kind == Kind::XN)
    {
        // on each rank, generate a buffer which contains rows_in_chunk rows of the data
        auto                rows_in_chunk = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10 * NCols);
        std::fill_n(foo.begin(), rows_in_chunk * NCols, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10);

        // after gathering, we have a buffer with 10*NCols elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
            ASSERT_EQ(buf_1, ref);
        }
    }
    else if constexpr (kind == Kind::XY)
    {
        // on each rank, generate a buffer which contains rows_in_chunk rows of the data
        auto                rows_in_chunk = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10 * 5);
        std::fill_n(foo.begin(), rows_in_chunk * 5, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo, 10, 5);

        // after gathering, we have a buffer with 10*5 elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);
            ASSERT_EQ(buf_1, ref);
        }
    }
    else if constexpr (kind == Kind::MN)
    {
        // on each rank, generate a buffer which contains rows_in_chunk rows of the data
        auto                rows_in_chunk = mpi::batch_size(NRows, rank, sz);
        std::vector<Scalar> foo(NRows * NCols);
        std::fill_n(foo.begin(), rows_in_chunk * NCols, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        // after gathering, we have a buffer with NRows*NCols elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
            ASSERT_EQ(buf_1, ref);
        }
    }
    else
    {
        // on each rank, generate a buffer which contains a chunk_sz part of the data
        auto                chunk_sz = mpi::batch_size(10, rank, sz);
        std::vector<Scalar> foo(10);
        std::fill_n(foo.begin(), chunk_sz, Scalar{42});
        auto buf_0 = CBuffer<Scalar, Backend, NRows, NCols>(foo);

        // after gathering, we have a buffer with 10 elements equal to 42
        auto buf_1 = buf_0.gather(rank, sz, MPI_COMM_WORLD);

        if (rank == mpi::master())
        {
            auto ref = CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
            ASSERT_EQ(buf_1, ref);
        }
    }
}

TYPED_TEST(CBufferTest, Slice)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    // take a slice containing rows {0, 1} and columns {2, 3, 4}
    auto rows = seqN(0, 2);
    auto cols = seqN(2, 3);

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 10, 10);

        auto ref = CBuffer<Scalar, Backend, 1, Dynamic>(std::vector<Scalar>{2, 3, 4}, cols.size());

        ASSERT_EQ(buf.slice(cols), ref);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NRows * 5, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{2, 3, 4, 7, 8, 9}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 5 * NCols, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{2, 3, 4, 12, 13, 14}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 15, 3, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{2, 3, 4, 7, 8, 9}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NRows * NCols);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{2, 3, 4, 7, 8, 9}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NCols);

        auto ref = CBuffer<Scalar, Backend, 1, Dynamic>(std::vector<Scalar>{2, 3, 4}, cols.size());

        ASSERT_EQ(buf.slice(cols), ref);
    }
}

TYPED_TEST(CBufferTest, SliceReverse)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = getKind<NRows, NCols>();

    // take a slice containing rows {1, 0} and columns {4, 3, 2}
    auto rows = seqN(0, 2).reverse();
    auto cols = seqN(2, 3).reverse();

    if constexpr (kind == Kind::X)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 10, 10);

        auto ref = CBuffer<Scalar, Backend, 1, Dynamic>(std::vector<Scalar>{4, 3, 2}, cols.size());

        ASSERT_EQ(buf.slice(cols), ref);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NRows * 5, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{9, 8, 7, 4, 3, 2}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 5 * NCols, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{14, 13, 12, 4, 3, 2}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, 15, 3, 5);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{9, 8, 7, 4, 3, 2}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else if constexpr (kind == Kind::MN)
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NRows * NCols);

        auto ref = CBuffer<Scalar, Backend, Dynamic, Dynamic>(std::vector<Scalar>{9, 8, 7, 4, 3, 2}, rows.size(), cols.size());

        ASSERT_EQ(buf.slice(rows, cols), ref);
    }
    else
    {
        auto buf = CBuffer<Scalar, Backend, NRows, NCols>::Linspace(0, NCols);

        auto ref = CBuffer<Scalar, Backend, 1, Dynamic>(std::vector<Scalar>{4, 3, 2}, cols.size());

        ASSERT_EQ(buf.slice(cols), ref);
    }
}

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

using namespace buffer;

namespace detail {
using implementations = ::testing::Types<
    /* @{ BufferX tests */
    CBufferParameters<int32_t, mem::Host, 1, Dynamic>,
    CBufferParameters<float, mem::Host, 1, Dynamic>,
    CBufferParameters<double, mem::Host, 1, Dynamic>,
    /* @} */
    /* @{ BufferN tests */
    CBufferParameters<int32_t, mem::Host, 1, 10>,
    CBufferParameters<float, mem::Host, 1, 10>,
    CBufferParameters<double, mem::Host, 1, 10>,
    /* @} */
    /* @{ BufferXY tests */
    CBufferParameters<int32_t, mem::Host, Dynamic, Dynamic>,
    CBufferParameters<float, mem::Host, Dynamic, Dynamic>,
    CBufferParameters<double, mem::Host, Dynamic, Dynamic>,
    /* @} */
    /* @{ BufferXN tests */
    CBufferParameters<int32_t, mem::Host, Dynamic, 10>,
    CBufferParameters<float, mem::Host, Dynamic, 10>,
    CBufferParameters<double, mem::Host, Dynamic, 10>,
    /* @} */
    /* @{ BufferMY tests */
    CBufferParameters<int32_t, mem::Host, 10, Dynamic>,
    CBufferParameters<float, mem::Host, 10, Dynamic>,
    CBufferParameters<double, mem::Host, 10, Dynamic>,
    /* @} */
    /* @{ BufferMN tests */
    CBufferParameters<int32_t, mem::Host, 10, 5>,
    CBufferParameters<float, mem::Host, 10, 5>,
    CBufferParameters<double, mem::Host, 10, 5>
    /* @} */
    >;
}  // namespace detail

TYPED_TEST_SUITE(CBufferTest, detail::implementations);

TYPED_TEST(CBufferTest, DefaultConstructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

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
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto data = new Scalar[100];

    if constexpr (kind == Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), NRows * 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.size(), 10 * NCols);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_EQ(buf.size(), 50);
    }
    else if constexpr (kind == Kind::N)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data);

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto data = std::vector<Scalar>(100);

    if constexpr (kind == Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), NRows * 10);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.size(), 10 * NCols);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_EQ(buf.size(), 50);
    }
    else if constexpr (kind == Kind::N)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_EQ(buf.size(), 10);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data);

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == Kind::X)
    {
        auto data = std::vector<Scalar>(10);

        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(data);

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> src{};
    buffer::CBuffer<Scalar, Backend, NRows, NCols> dst{src};

    // the two objects are exactly equal, but independent from each other
    ASSERT_EQ(dst, src);
}

TYPED_TEST(CBufferTest, CopyAssignment)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    buffer::CBuffer<Scalar, Backend, NRows, NCols> src{};
    buffer::CBuffer<Scalar, Backend, NRows, NCols> dst{};

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> src{};

    auto nrows = src.nRows();
    auto ncols = src.nColumns();

    buffer::CBuffer<Scalar, Backend, NRows, NCols> dst{std::move(src)};

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> src{};

    auto nrows = src.nRows();
    auto ncols = src.nColumns();

    buffer::CBuffer<Scalar, Backend, NRows, NCols> dst{};

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.resize(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == buffer::Kind::XY)
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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.setZero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
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
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setZero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.resize(10);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.resize(10, 5);

        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.setConstant(Scalar{42}, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.setConstant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
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
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.resize(10);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.resize(10, 5);

        buf.setConstant(Scalar{42});

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == buffer::Kind::X)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.setRandom(LB, UB, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.setRandom(LB, UB, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == buffer::Kind::MN)
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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>(10, 5);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>();

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

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == buffer::Kind::X)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.resize(10);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.resize(10, 5);

        buf.setRandom(LB, UB);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == buffer::Kind::MN)
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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero();

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{0}, 1.0e-14);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero();

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_NEAR(buf(2, 7), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_NEAR(buf(7, 2), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_NEAR(buf(4, 3), Scalar{42}, 1.0e-14);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    auto LB = Scalar{1};
    auto UB = Scalar{5};

    // check that the value in the buffer is between the bounds of the RNG
    auto pred = [LB, UB](Scalar x) { return ((x - UB) * (x - LB) <= 0); };

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(7));
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nColumns(), 10);
        ASSERT_PRED1(pred, buf(2, 7));
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_PRED1(pred, buf(7, 2));
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB, 10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB);

        ASSERT_EQ(buf.nRows(), NRows);
        ASSERT_EQ(buf.nColumns(), NCols);
        ASSERT_PRED1(pred, buf(4, 3));
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Random(LB, UB);

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

    constexpr auto kind = buffer::getKind<NRows, NCols>();

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42}, 10, 5);

        ASSERT_EQ(buf, ref);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf, ref);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});
        buf.broadcast(0, MPI_COMM_WORLD);

        auto ref = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(Scalar{42});

        ASSERT_EQ(buf, ref);
    }
}

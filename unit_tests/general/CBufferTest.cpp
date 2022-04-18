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
#include <tuple>
#include <utility>
#include <vector>

#include "Buffer.hpp"
#include "MemAlloc.hpp"

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

TYPED_TEST(CBufferTest, DefaultContructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    ASSERT_EQ(buf.nRows(), NRows);
    ASSERT_EQ(buf.nColumns(), NCols);
}

TYPED_TEST(CBufferTest, DefaultContructorAndResize)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = decltype(buf)::kind;

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
        // do nothing for Kind::N and Kind::MN
    }
}

/* Test default construction and `setZero`
 */
TYPED_TEST(CBufferTest, DefaultContructorSetZero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    buffer::CBuffer<Scalar, Backend, NRows, NCols> buf;

    constexpr auto kind = decltype(buf)::kind;

    if constexpr (kind == buffer::Kind::X)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        buf.setZero(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        buf.setZero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        buf.setZero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else
    {
        buf.setZero();

        ASSERT_EQ(buf.nColumns(), 10);
    }
}

TYPED_TEST(CBufferTest, CopyConstructor)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = buffer::CBuffer<Scalar, Backend, NRows, NCols>::kind;

    auto val = static_cast<Scalar>(3);

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val, 10);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val, 10);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val, 10);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val, 10, 5);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
    else
    {
        auto buf_src = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Constant(val);
        auto buf_dst{buf_src};

        ASSERT_EQ(buf_dst, buf_src);
    }
}

/* Test static generators: `Zero`
 */
TYPED_TEST(CBufferTest, Zero)
{
    using Scalar         = typename TypeParam::value_type;
    using Backend        = typename TypeParam::backend_type;
    constexpr auto NRows = TypeParam::NRows;
    constexpr auto NCols = TypeParam::NCols;

    constexpr auto kind = buffer::CBuffer<Scalar, Backend, NRows, NCols>::kind;

    if constexpr (kind == buffer::Kind::X)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::MY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nColumns(), 10);
    }
    else if constexpr (kind == buffer::Kind::XN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10);

        ASSERT_EQ(buf.nRows(), 10);
    }
    else if constexpr (kind == buffer::Kind::XY)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero(10, 5);

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else if constexpr (kind == buffer::Kind::MN)
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero();

        ASSERT_EQ(buf.nRows(), 10);
        ASSERT_EQ(buf.nColumns(), 5);
    }
    else
    {
        auto buf = buffer::CBuffer<Scalar, Backend, NRows, NCols>::Zero();

        ASSERT_EQ(buf.nColumns(), 10);
    }
}

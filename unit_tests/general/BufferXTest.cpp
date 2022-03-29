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

#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

#include "Buffer.hpp"
#include "BufferTest.hpp"
#include "MemAlloc.hpp"

namespace detail {
using implementations = ::testing::Types<std::tuple<int32_t, mem::Host>, std::tuple<float, mem::Host>, std::tuple<double, mem::Host>>;
}

TYPED_TEST_CASE(BufferXTest, detail::implementations);

TYPED_TEST(BufferXTest, DefaultContructor)
{
    auto type_param = TypeParam();
    using Scalar    = typename std::tuple_element<0, decltype(type_param)>::type;
    using Backend   = typename std::tuple_element<1, decltype(type_param)>::type;

    BufferX<Scalar, Backend> buf_a;

    ASSERT_EQ(buf_a.size(), Dynamic);
}

TYPED_TEST(BufferXTest, DefaultContructorAndResize)
{
    auto type_param = TypeParam();
    using Scalar    = typename std::tuple_element<0, decltype(type_param)>::type;
    using Backend   = typename std::tuple_element<1, decltype(type_param)>::type;

    BufferX<Scalar, Backend> buf_a;
    buf_a.resize(10);

    ASSERT_EQ(buf_a.size(), 10);
}

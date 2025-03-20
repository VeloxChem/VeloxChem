//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef BoysFuncGPU_hpp
#define BoysFuncGPU_hpp

#include "GpuConstants.hpp"

namespace gpu {  // gpu namespace

template<typename TPtr>
using CharPtr = std::conditional_t<std::is_const_v<std::remove_pointer_t<TPtr>>, const char*, char*>;

template<typename ValueType, typename IndexType, std::enable_if_t<std::is_integral<IndexType>::value, bool> = true>
static inline __device__ __attribute__((always_inline)) IndexType calculateOffset(IndexType index)
{
    __builtin_assume(index >= 0);
    return index * static_cast<IndexType>(sizeof(ValueType));
}

template<typename PointerType, typename IndexType, std::enable_if_t<std::is_integral<IndexType>::value, bool> = true>
static inline __device__ __attribute__((always_inline)) PointerType indexedAddress(PointerType address,
                                                                                          IndexType idx)
{
    return reinterpret_cast<PointerType>(reinterpret_cast<CharPtr<decltype(address)>>(address)
                                         + calculateOffset<std::remove_pointer_t<PointerType>>(idx));
}

template<typename ValueType>
class AmdFastBuffer
{
private:
    ValueType* buffer;

public:
    __device__ AmdFastBuffer(ValueType* buffer) : buffer(buffer) {}
    template<typename IndexType, std::enable_if_t<std::is_integral<IndexType>::value && std::is_const_v<std::remove_pointer_t<ValueType>>, bool> = true>
    inline __device__ __attribute__((always_inline)) const ValueType& operator[](IndexType idx) const
    {
        return *indexedAddress(buffer, idx);
    }
    template<typename IndexType, std::enable_if_t<std::is_integral<IndexType>::value && !std::is_const_v<std::remove_pointer_t<ValueType>>, bool> = true>
    inline __device__ __attribute__((always_inline)) ValueType& operator[](IndexType idx)
    {
        return *indexedAddress(buffer, idx);
    }
};

__device__ void
computeBoysFunction(double* values_in, const double fa, const uint32_t N, const double* bf_table, const double* ft_in)
{
    // Note: 847 = 121 * 7
    AmdFastBuffer<const double> bf_data{bf_table + N * 847};
    AmdFastBuffer<const double> ft{ft_in};
    AmdFastBuffer<double> values{values_in};

    uint32_t pnt = (fa > 1.0e5) ? 1000000 : static_cast<uint32_t>(10.0 * fa + 0.5);

    if (pnt < 121)
    {
        const double w = fa - 0.1 * pnt;

        const double w2 = w * w;

        const double w4 = w2 * w2;

        values[N] = bf_data[pnt * 7 + 0] + bf_data[pnt * 7 + 1] * w + bf_data[pnt * 7 + 2] * w2 + bf_data[pnt * 7 + 3] * w2 * w

                    + bf_data[pnt * 7 + 4] * w4 + bf_data[pnt * 7 + 5] * w4 * w + bf_data[pnt * 7 + 6] * w4 * w2;

        const double f2a = fa + fa;

        const double fx = exp(-fa);

        for (uint32_t j = 0; j < N; j++)
        {
            values[N - j - 1] = ft[N - j - 1] * (f2a * values[N - j]+ fx);
        }
    }
    else
    {
        const double fia = 1.0 / fa;

        double pf = 0.5 * fia;

        values[0] = MATH_CONST_HALF_SQRT_PI * sqrt(fia);

        if (pnt < 921)
        {
            const double fia2 = fia * fia;

            const double f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

            const double fx = exp(-fa);

            values[0] -= f * fx;

            const double rterm = pf * fx;

            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1] - rterm;

                pf += fia;
            }
        }
        else
        {
            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1];

                pf += fia;
            }
        }
    }
}

}  // namespace gpu

#endif /* BoysFuncGPU_hpp */

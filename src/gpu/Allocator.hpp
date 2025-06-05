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

#ifndef Allocator_hpp
#define Allocator_hpp

#include <hip/hip_runtime.h>

#include <memory>
#include <new>
#include <cstdint>
#include <vector>
#include <string>

#include "GpuWrapper.hpp"
#include "GpuSafeChecks.hpp"

template<class T, typename VeloxPolicy>
class VeloxAllocator : public VeloxPolicy
{
public:
    typedef T value_type;
    VeloxAllocator() = default;
    VeloxAllocator(const VeloxPolicy& p) : VeloxPolicy(p) {}

    value_type* allocate(std::size_t n)
    {
        void* p = VeloxPolicy::malloc(n * sizeof(T));
        if (p == nullptr)
        {
            throw std::bad_alloc();
        }
        else
        {
            return static_cast<value_type*>(p);
        }
    }

    void deallocate(value_type* p, std::size_t n)
    {
        VeloxPolicy::free(p);
    }
    using is_always_equal = std::true_type;
};

class VeloxHostPolicy;

template<class T>
using VeloxHostAllocator = VeloxAllocator<T, VeloxHostPolicy>;

template<class T>
using VeloxHostVector = std::vector<T, VeloxHostAllocator<T>>;

class VeloxHostPolicy
{
public:
    VeloxHostPolicy() = default;
    std::size_t alignment() const noexcept
    {
        // we just guess, no need to overcomplicate things
        long pageSize = 4096;
        return static_cast<std::size_t>(pageSize);
    }
    void* malloc(std::size_t bytes) const noexcept
    {
        int flag = gpuHostMallocDefault;
        if (bytes == 0)
        {
            return nullptr;
        }
        void* p;
        gpuSafe(gpuHostMalloc(&p, bytes, flag));
        return p;
    }
    void free(void* buffer) const noexcept
    {
        if (buffer == nullptr)
        {
            return;
        }
        gpuSafe(gpuHostFree(buffer));
    }
};

#endif /* Allocator_hpp */

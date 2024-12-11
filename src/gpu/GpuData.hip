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

#include <hip/hip_runtime.h>

#include "GpuData.hpp"

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuDevices.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

#define MATH_CONST_PI 3.14159265358979323846

GpuData::GpuData(const size_t num_gpus_per_node)
	: _num_gpus_per_node(num_gpus_per_node)
{
     // we can allocate the boys function tables right at construction time
     // as they never change
     allocateBoysFunctionData();
}

GpuData::~GpuData()
{
	freeStorage();
}


auto GpuData::allocateBoysFunctionData() -> void
{
     const auto boys_func_table = boysfunc::getFullBoysFuncTable();
     const auto boys_func_ft = boysfunc::getBoysFuncFactors();
     _storage.resize(_num_gpus_per_node);

     // each device gets its own copy of the data
     for (size_t gpuId = 0; gpuId < _num_gpus_per_node; ++gpuId)
     {
         hipSafe(hipSetDevice(gpuId));

         hipSafe(hipMalloc(&_storage[gpuId].d_boys_func_table, boys_func_table.size() * sizeof(double)));
         hipSafe(hipMemcpy(_storage[gpuId].d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

         hipSafe(hipMalloc(&_storage[gpuId].d_boys_func_ft, boys_func_ft.size() * sizeof(double)));
         hipSafe(hipMemcpy(_storage[gpuId].d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));
     }
}

auto GpuData::allocatePrimInfo(size_t gpuId, size_t numS, size_t numP, size_t numD) -> void
{
    assert(gpuId < _storage.size());
    if (_storage[gpuId].num_s_prim < numS)
    {
        hipSafe(hipFree(_storage[gpuId].d_s_prim_info));
	hipSafe(hipMalloc(&_storage[gpuId].d_s_prim_info, numS * sizeof(double)));
	_storage[gpuId].num_s_prim = numS;
    }
    if (_storage[gpuId].num_p_prim < numP)
    {
        hipSafe(hipFree(_storage[gpuId].d_p_prim_info));
	hipSafe(hipMalloc(&_storage[gpuId].d_p_prim_info, numP * sizeof(double)));
	_storage[gpuId].num_p_prim = numP;
    }
    if (_storage[gpuId].num_d_prim < numD)
    {
        hipSafe(hipFree(_storage[gpuId].d_d_prim_info));
	hipSafe(hipMalloc(&_storage[gpuId].d_d_prim_info, numD * sizeof(double)));
	_storage[gpuId].num_d_prim = numD;
    }
}

auto GpuData::accessStorage(size_t gpuId) const -> const Storage&
{
	return _storage[gpuId];
}

auto GpuData::accessStorage(size_t gpuId) -> Storage&
{
	return _storage[gpuId];
}

auto GpuData::freeStorage() -> void
{
	for (size_t gpuId = 0; gpuId < _num_gpus_per_node; ++gpuId)
	{
	     hipSafe(hipFree(_storage[gpuId].d_boys_func_table));
	     hipSafe(hipFree(_storage[gpuId].d_boys_func_ft));
	     hipSafe(hipFree(_storage[gpuId].d_s_prim_info));
	     hipSafe(hipFree(_storage[gpuId].d_p_prim_info));
	     hipSafe(hipFree(_storage[gpuId].d_d_prim_info));
	}
}

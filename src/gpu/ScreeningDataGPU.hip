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

#include "ScreeningData.hpp"

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


auto CScreeningData::allocateGpuBoysFunctionData() -> void
{
     const auto boys_func_table = boysfunc::getFullBoysFuncTable();
     const auto boys_func_ft = boysfunc::getBoysFuncFactors();
     _gpuData.resize(_num_gpus_per_node);

     // each device gets its own copy of the data
     for (size_t gpuId = 0; gpuId < _num_gpus_per_node; ++gpuId)
     {
         hipSafe(hipSetDevice(gpuId));

         hipSafe(hipMalloc(&_gpuData[gpuId].d_boys_func_table, boys_func_table.size() * sizeof(double)));
         hipSafe(hipMemcpy(_gpuData[gpuId].d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), hipMemcpyHostToDevice));

         hipSafe(hipMalloc(&_gpuData[gpuId].d_boys_func_ft, boys_func_ft.size() * sizeof(double)));
         hipSafe(hipMemcpy(_gpuData[gpuId].d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), hipMemcpyHostToDevice));
     }
}

auto CScreeningData::accessGpuData(size_t gpuId) const -> const GpuData&
{
	return _gpuData[gpuId];
}

auto CScreeningData::freeGpuData() -> void
{
	for (size_t gpuId = 0; gpuId < _num_gpus_per_node; ++gpuId)
	{
	     hipSafe(hipFree(_gpuData[gpuId].d_boys_func_table));
	     hipSafe(hipFree(_gpuData[gpuId].d_boys_func_ft));
	}
}

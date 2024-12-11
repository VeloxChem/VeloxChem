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

#ifndef GpuData_hpp
#define GpuData_hpp

#include <cstdint>
#include <vector>
#include <string>

/**
 Class GpuData stores pointers to GPU resident data structures.

 @author Paul Bauer <paul.bauer@amd.com>
 */
class GpuData
{
    struct Storage
    {
	double* d_boys_func_table = nullptr;
	double* d_boys_func_ft = nullptr;
	double* d_s_prim_info = nullptr;
	double* d_p_prim_info = nullptr;
	double* d_d_prim_info = nullptr;
	double* d_mat_Q = nullptr;
	double* d_mat_S = nullptr;
	double* d_mat_T = nullptr;
	                                                          
	uint32_t *d_ss_first_inds_local = nullptr;
	uint32_t *d_sp_first_inds_local = nullptr;
	uint32_t *d_sd_first_inds_local = nullptr;
	uint32_t *d_pp_first_inds_local = nullptr;
	uint32_t *d_pd_first_inds_local = nullptr;
	uint32_t *d_dd_first_inds_local = nullptr;

	uint32_t  *d_ss_second_inds_local = nullptr;
	uint32_t  *d_sp_second_inds_local = nullptr;
	uint32_t  *d_sd_second_inds_local = nullptr;
	uint32_t  *d_pp_second_inds_local = nullptr;
	uint32_t  *d_pd_second_inds_local = nullptr;
	uint32_t  *d_dd_second_inds_local = nullptr;

	size_t num_s_prim = 0;
	size_t num_d_prim = 0;
	size_t num_p_prim = 0;
    };
    std::vector<Storage> _storage;

    size_t _num_gpus_per_node;

   public:
    GpuData(const size_t num_gpu_per_node);

    ~GpuData();

    auto allocateBoysFunctionData() -> void;

    auto allocatePrimInfo(size_t gpuId, size_t numS, size_t numP, size_t numD) -> void;

    auto accessStorage(size_t gpuId) const -> const Storage&;

    auto accessStorage(size_t gpuId) -> Storage&;

    auto freeStorage() -> void;
};

#endif /* GpuData_hpp */

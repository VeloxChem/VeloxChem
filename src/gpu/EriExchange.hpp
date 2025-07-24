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

#ifndef EriExchange_hpp
#define EriExchange_hpp

#include <cstdint>

#include "GpuConstants.hpp"
#include "BoysFuncGPU.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSSS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double    ss_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_counts_K_ss,
                        const double*   pair_data_K_ss,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSSP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_sp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_sp,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_sp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPSS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    ss_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_ss,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_ss,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_ss,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPSP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_counts_K_sp,
                        const double*   pair_data_K_sp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_ps,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_ps,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_pp,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_pp,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSSD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_sd,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_sd,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_sd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDSS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ss_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_ss,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_ss,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_ss,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPSD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_sd,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_sd,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_sd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDSP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_sp,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_sp,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_sp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDSD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_ss,
                        const uint32_t* pair_inds_k_for_K_ss,
                        const uint32_t  pair_inds_count_for_K_ss,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_counts_K_sd,
                        const double*   pair_data_K_sd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sp,
                        const uint32_t* pair_inds_k_for_K_sp,
                        const uint32_t  pair_inds_count_for_K_sp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSSDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_ps,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_ps,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDPS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_ps,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_ps,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_ps,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_ps,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_ps,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_ps,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_ps,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_pp,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_pp,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDSDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD0(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD1(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD2(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD3(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD4(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD5(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDD6(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP0(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP1(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP2(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP3(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP4(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP5(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDP6(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

template<uint32_t part>
__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold)
{
    // each thread block scans over [i?|k?] and sum up to a primitive K matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM_Y_K][TILE_DIM_X_K + 1];
    __shared__ uint32_t skip_thread_block, i, k, count_i, count_k, displ_i, displ_k;
    __shared__ double   a_i, r_i[3], a_k, r_k[3];
    __shared__ uint32_t a0, a1, c0, c1;
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    const uint32_t ik = blockIdx.x;

    double K_ik = 0.0, d2 = 1.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {

        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        if (ik < pair_inds_count_for_K_dd)
        {
            i = pair_inds_i_for_K_dd[ik];
            k = pair_inds_k_for_K_dd[ik];

            count_i = pair_counts_K_dd[i];
            count_k = pair_counts_K_dd[k];

            displ_i = pair_displs_K_dd[i];
            displ_k = pair_displs_K_dd[k];

            a_i = d_prim_info[i / 6 + d_prim_count * 0];

            r_i[0] = d_prim_info[i / 6 + d_prim_count * 2];
            r_i[1] = d_prim_info[i / 6 + d_prim_count * 3];
            r_i[2] = d_prim_info[i / 6 + d_prim_count * 4];

            a_k = d_prim_info[k / 6 + d_prim_count * 0];

            r_k[0] = d_prim_info[k / 6 + d_prim_count * 2];
            r_k[1] = d_prim_info[k / 6 + d_prim_count * 3];
            r_k[2] = d_prim_info[k / 6 + d_prim_count * 4];

            a0 = d_cart_inds[i % 6][0];
            a1 = d_cart_inds[i % 6][1];
            c0 = d_cart_inds[k % 6][0];
            c1 = d_cart_inds[k % 6][1];

        }
        else
        {
            count_i = 0;
            count_k = 0;
        }

    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // to avoid memory hazard
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_0, PA_1, PB_0, PB_1;
        uint32_t j_prim, j_cgto, b0, b1;

        if ((ik < pair_inds_count_for_K_dd) && (j < count_i))
        {
            Q_ij   = Q_K_dd[displ_i + j];

            j_prim = D_inds_K_dd[displ_i + j];

            j_cgto = d_prim_aoinds[(j_prim / 6) + d_prim_count * (j_prim % 6)];

            a_j = d_prim_info[j_prim / 6 + d_prim_count * 0];

            r_j[0] = d_prim_info[j_prim / 6 + d_prim_count * 2];
            r_j[1] = d_prim_info[j_prim / 6 + d_prim_count * 3];
            r_j[2] = d_prim_info[j_prim / 6 + d_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_dd[displ_i + j];

            b0 = d_cart_inds[j_prim % 6][0];
            b1 = d_cart_inds[j_prim % 6][1];

            PA_0 = (a_j  * inv_S1) * (r_j[a0] - r_i[a0]);
            PA_1 = (a_j  * inv_S1) * (r_j[a1] - r_i[a1]);
            PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);
            PB_1 = (-a_i * inv_S1) * (r_j[b1] - r_i[b1]);

        }

        if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = 0;

        __syncthreads();

        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            if ((ik < pair_inds_count_for_K_dd) && (j < count_i) && (l < count_k))
            {
                const auto Q_kl = Q_K_dd[displ_k + l];

                if (fabs(Q_ij * Q_kl * dd_max_D) > eri_threshold)
                {
                    const auto l_prim = D_inds_K_dd[displ_k + l];

                    const auto l_cgto = d_prim_aoinds[(l_prim / 6) + d_prim_count * (l_prim % 6)];

                    const auto a_l = d_prim_info[l_prim / 6 + d_prim_count * 0];

                    const double r_l[3] = {d_prim_info[l_prim / 6 + d_prim_count * 2],
                                           d_prim_info[l_prim / 6 + d_prim_count * 3],
                                           d_prim_info[l_prim / 6 + d_prim_count * 4]};

                    const auto S_kl_00 = pair_data_K_dd[displ_k + l];

                    const auto d0 = d_cart_inds[l_prim % 6][0];
                    const auto d1 = d_cart_inds[l_prim % 6][1];

                    // J. Chem. Phys. 84, 3963-3974 (1986)

                    const auto S2 = a_k + a_l;

                    const auto inv_S2 = 1.0 / S2;
                    const auto inv_S4 = 1.0 / (S1 + S2);

                    const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                          (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                          (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

                    const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

                    const auto rho = S1 * S2 * inv_S4;

                    if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

                    const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

                    double F8_t[9];

                    gpu::computeBoysFunction(F8_t, rho * d2 * r2_PQ, 8, boys_func_table, boys_func_ft);

                    if (omega != 0.0)
                    {
                        F8_t[1] *= d2;
                        F8_t[2] *= d2 * d2;
                        F8_t[3] *= d2 * d2 * d2;
                        F8_t[4] *= d2 * d2 * d2 * d2;
                        F8_t[5] *= d2 * d2 * d2 * d2 * d2;
                        F8_t[6] *= d2 * d2 * d2 * d2 * d2 * d2;
                        F8_t[7] *= d2 * d2 * d2 * d2 * d2 * d2 * d2;
                        F8_t[8] *= d2 * d2 * d2 * d2 * d2 * d2 * d2 * d2;
                    }

                    const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
                    const auto QC_1 = (a_l * inv_S2) * (r_l[c1] - r_k[c1]);
                    const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);
                    const auto QD_1 = (-a_k * inv_S2) * (r_l[d1] - r_k[d1]);

                    const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * [&]()
                    {
                        if constexpr(part == 0)
                        {
                            return
                            F8_t[0] * (

                                0.125 * inv_S1 * inv_S1 * inv_S2 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.125 * inv_S1 * inv_S2 * inv_S2 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.25 * inv_S1 * inv_S1 * (
                                    (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.25 * inv_S1 * inv_S2 * (
                                    delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * QD_0 * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * QD_0 * QD_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * QD_0 * QD_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.25 * inv_S2 * inv_S2 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.5 * inv_S1 * (
                                    delta[b0][b1] * (PA_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PB_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.5 * inv_S2 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[0] * (

                                + (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[0] * (

                                + 0.0625 * inv_S1 * inv_S1 * inv_S2 * inv_S2 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            )
                            +
                            F8_t[1] * (

                                (-0.125) * inv_S1 * inv_S1 * inv_S2 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.125) * inv_S1 * inv_S2 * inv_S2 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.25) * inv_S1 * inv_S1 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.25) * inv_S2 * inv_S2 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.5) * S1 * inv_S2 * inv_S2 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.5) * S2 * inv_S1 * inv_S1 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + (-0.5) * S1 * inv_S2 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + 0.5 * S2 * inv_S1 * inv_S4 * (
                                    delta[b0][b1] * (PA_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 + PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 + PA_0 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 + PA_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 + PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 + PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PB_1 * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 + PB_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            );
                        }
                        else if constexpr(part == 1)
                        {
                            return
                            F8_t[1] * (

                                + 0.125 * inv_S1 * inv_S2 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1 * (-1.0) + PA_0 * PQ[a1] + PA_1 * PQ[a0])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0 * (-1.0) + PB_0 * PQ[a0] + PA_0 * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * (-1.0) + PB_1 * PQ[a0] + PA_0 * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_0 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * QD_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * QD_1)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1 * (-1.0) + PB_0 * PQ[a1] + PA_1 * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1 * (-1.0) + PB_1 * PQ[a1] + PA_1 * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_1 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_1 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * QD_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * QD_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * (-1.0) + PB_0 * PQ[b1] + PB_1 * PQ[b0])
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PB_0 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PB_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PB_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PB_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * QC_1 * (-1.0) + PQ[c1] * QC_0 * (-1.0) + QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * QD_0 * (-1.0) + PQ[d0] * QC_0 * (-1.0) + QD_0 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * QD_1 * (-1.0) + PQ[d1] * QC_0 * (-1.0) + QD_1 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * QD_0 * (-1.0) + PQ[d0] * QC_1 * (-1.0) + QD_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * QD_1 * (-1.0) + PQ[d1] * QC_1 * (-1.0) + QD_1 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * QD_1 * (-1.0) + PQ[d1] * QD_0 * (-1.0) + QD_0 * QD_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * QD_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * QD_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * QD_1)
                                )

                            )
                            +
                            F8_t[1] * (

                                + 0.25 * inv_S1 * inv_S4 * (
                                    delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[a0] * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[a0] * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[a0] * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[a0] * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[a0] * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[a0] * QD_0 * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b0] * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b0] * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b1] * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b1] * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b1] * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b1] * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b1] * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * QD_0 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * QD_1 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * QD_0 * QD_1 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b0] * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * QD_0 * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * QC_0 * QC_1 + PB_1 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QC_0 + PB_1 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[b1] * QD_1 * QC_0 + PB_1 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QC_1 + PB_1 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[b1] * QD_1 * QC_1 + PB_1 * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QD_1 + PB_1 * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                )

                            )
                            +
                            F8_t[1] * (

                                + 0.25 * inv_S2 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PQ[a1] + PB_0 * PB_1 * PA_1 * PQ[a0] + PB_0 * PA_0 * PA_1 * PQ[b1] + PB_1 * PA_0 * PA_1 * PQ[b0])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PA_1 * QC_0)
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PA_1 * QC_1)
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * QD_0)
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * QD_1)
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PA_1 * QC_0)
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PA_1 * QC_1)
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * QD_0)
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * QD_1)
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PA_0 * PA_1 * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PA_0 * PA_1 * QD_0 * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PA_0 * PA_1 * QD_1 * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PA_0 * PA_1 * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PA_0 * PA_1 * QD_1 * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PA_0 * PA_1 * QD_0 * QD_1 * (-1.0))
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * QC_0)
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_0 * QC_1)
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * QD_0)
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PA_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PA_0 * QD_0 * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_0 * QD_1 * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PA_0 * QD_0 * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_0 * QD_1 * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_0 * QD_0 * QD_1 * (-1.0))
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QC_0 * (-1.0) + PB_1 * PA_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_0 * (-1.0) + PB_1 * PA_0 * QD_0 * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_0 * QD_1 * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_1 * (-1.0) + PB_1 * PA_0 * QD_0 * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_0 * QD_1 * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_0 * QD_0 * QD_1 * (-1.0))
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_1 * QC_0)
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_1 * QC_1)
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * QD_0)
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * QD_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PA_1 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PA_1 * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_1 * QD_1 * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PA_1 * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_1 * QD_1 * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_1 * QD_0 * QD_1 * (-1.0))
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PB_1 * PA_1 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PB_1 * PA_1 * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_1 * QD_1 * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PB_1 * PA_1 * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_1 * QD_1 * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_1 * QD_0 * QD_1 * (-1.0))
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PB_1 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PB_1 * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PB_1 * QD_1 * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PB_1 * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PB_1 * QD_1 * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PB_1 * QD_0 * QD_1 * (-1.0))
                                )

                            );
                        }
                        else if constexpr(part == 2)
                        {
                            return
                            F8_t[1] * (

                                + 0.5 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QC_0 * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * QC_0 * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * QC_0 * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QC_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QC_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QC_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QD_1 * QC_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_1 * QC_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_1 * QC_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QD_1 * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_1 * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_1 * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1)
                                    + delta[b1][d1] * (PB_0 * PA_0 * PA_1 * QD_0 * QC_0 * QC_1)
                                    + delta[b1][d0] * (PB_0 * PA_0 * PA_1 * QD_1 * QC_0 * QC_1)
                                    + delta[b1][c1] * (PB_0 * PA_0 * PA_1 * QD_0 * QD_1 * QC_0)
                                    + delta[b1][c0] * (PB_0 * PA_0 * PA_1 * QD_0 * QD_1 * QC_1)
                                    + delta[b0][d1] * (PB_1 * PA_0 * PA_1 * QD_0 * QC_0 * QC_1)
                                    + delta[b0][d0] * (PB_1 * PA_0 * PA_1 * QD_1 * QC_0 * QC_1)
                                    + delta[b0][c1] * (PB_1 * PA_0 * PA_1 * QD_0 * QD_1 * QC_0)
                                    + delta[b0][c0] * (PB_1 * PA_0 * PA_1 * QD_0 * QD_1 * QC_1)
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][d1] * (PB_0 * PB_1 * PA_0 * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d0] * (PB_0 * PB_1 * PA_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][c1] * (PB_0 * PB_1 * PA_0 * QD_0 * QD_1 * QC_0)
                                    + delta[a1][c0] * (PB_0 * PB_1 * PA_0 * QD_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][d1] * (PB_0 * PB_1 * PA_1 * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d0] * (PB_0 * PB_1 * PA_1 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][c1] * (PB_0 * PB_1 * PA_1 * QD_0 * QD_1 * QC_0)
                                    + delta[a0][c0] * (PB_0 * PB_1 * PA_1 * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0))
                                )

                            )
                            +
                            F8_t[1] * (

                                + S1 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )

                            )
                            +
                            F8_t[1] * (

                                + S2 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[2] * (

                                0.125 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.125 * S2 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.125 * inv_S1 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PQ[a1] * (-1.0) + PA_1 * PQ[a0] * (-1.0) + PQ[a0] * PQ[a1])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * (-1.0) + PA_0 * PQ[b0] * (-1.0) + PQ[a0] * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * (-1.0) + PA_0 * PQ[b1] * (-1.0) + PQ[a0] * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_0 * QC_0 * (-1.0) + PQ[a0] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_0 * QC_1 * (-1.0) + PQ[a0] * QC_1)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * QD_0 * (-1.0) + PQ[a0] * QD_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * QD_1 * (-1.0) + PQ[a0] * QD_1)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * (-1.0) + PA_1 * PQ[b0] * (-1.0) + PQ[a1] * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a1] * (-1.0) + PA_1 * PQ[b1] * (-1.0) + PQ[a1] * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_1 * QC_0 * (-1.0) + PQ[a1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_1 * QC_1 * (-1.0) + PQ[a1] * QC_1)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * QD_0 * (-1.0) + PQ[a1] * QD_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * QD_1 * (-1.0) + PQ[a1] * QD_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[b1] * (-1.0) + PB_1 * PQ[b0] * (-1.0) + PQ[b0] * PQ[b1])
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PB_0 * QC_0 * (-1.0) + PQ[b0] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PB_0 * QC_1 * (-1.0) + PQ[b0] * QC_1)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * QD_0 * (-1.0) + PQ[b0] * QD_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * QD_1 * (-1.0) + PQ[b0] * QD_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PB_1 * QC_0 * (-1.0) + PQ[b1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PB_1 * QC_1 * (-1.0) + PQ[b1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * QD_0 * (-1.0) + PQ[b1] * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * QD_1 * (-1.0) + PQ[b1] * QD_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * QC_1 * 2.0 + PQ[c1] * QC_0 * 2.0 + QC_0 * QC_1 * 2.0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * QD_0 * 2.0 + PQ[d0] * QC_0 * 2.0 + QD_0 * QC_0 * 2.0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * QD_0 * 2.0 + PQ[d0] * QC_1 * 2.0 + QD_0 * QC_1 * 2.0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * QD_1 * 2.0 + PQ[d1] * QC_0 * 2.0 + QD_1 * QC_0 * 2.0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * QD_1 * 2.0 + PQ[d1] * QC_1 * 2.0 + QD_1 * QC_1 * 2.0)
                                    + (delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * QD_1 * 2.0 + PQ[d1] * QD_0 * 2.0 + QD_0 * QD_1 * 2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (QD_0 * QD_1)
                                )

                            );
                        }
                        else if constexpr(part == 13)
                        {
                            return
                            F8_t[2] * (

                                + 0.125 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PQ[a1] * (-2.0) + PA_1 * PQ[a0] * (-2.0) + PA_0 * PA_1 * 2.0)
                                    + (delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PA_0 * PA_1)
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * (-2.0) + PA_0 * PQ[b0] * (-2.0) + PB_0 * PA_0 * 2.0)
                                    + (delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0)
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * (-2.0) + PA_0 * PQ[b1] * (-2.0) + PB_1 * PA_0 * 2.0)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0)
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_0 * PQ[c0] * (-1.0) + PA_0 * QC_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_0 * PQ[c1] * (-1.0) + PA_0 * QC_1 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d0] * (-1.0) + PA_0 * QD_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d1] * (-1.0) + PA_0 * QD_1 * (-1.0))
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * (-2.0) + PA_1 * PQ[b0] * (-2.0) + PB_0 * PA_1 * 2.0)
                                    + (delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_1)
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a1] * (-2.0) + PA_1 * PQ[b1] * (-2.0) + PB_1 * PA_1 * 2.0)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_1)
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_1 * PQ[c0] * (-1.0) + PA_1 * QC_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_1 * PQ[c1] * (-1.0) + PA_1 * QC_1 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d0] * (-1.0) + PA_1 * QD_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d1] * (-1.0) + PA_1 * QD_1 * (-1.0))
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[b1] * (-2.0) + PB_1 * PQ[b0] * (-2.0) + PB_0 * PB_1 * 2.0)
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PB_0 * PQ[c0] * (-1.0) + PB_0 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PB_0 * PQ[c1] * (-1.0) + PB_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d0] * (-1.0) + PB_0 * QD_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d1] * (-1.0) + PB_0 * QD_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PB_1 * PQ[c0] * (-1.0) + PB_1 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PB_1 * PQ[c1] * (-1.0) + PB_1 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d0] * (-1.0) + PB_1 * QD_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d1] * (-1.0) + PB_1 * QD_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * PQ[c1] + PQ[c0] * QC_1 + PQ[c1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * PQ[d0] + PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * PQ[d0] + PQ[c1] * QD_0 + PQ[d0] * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * PQ[d1] + PQ[c0] * QD_1 + PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * PQ[d1] + PQ[c1] * QD_1 + PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * PQ[d1] + PQ[d0] * QD_1 + PQ[d1] * QD_0)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.25 * S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.25 * S1 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PQ[a1] * (-2.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * (-2.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * (-2.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * (-2.0))
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PA_1 * PQ[c0] * (-1.0) + PB_0 * PA_0 * PA_1 * QC_0 * (-1.0))
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PA_1 * PQ[c1] * (-1.0) + PB_0 * PA_0 * PA_1 * QC_1 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * PQ[d0] * (-1.0) + PB_0 * PA_0 * PA_1 * QD_0 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * PQ[d1] * (-1.0) + PB_0 * PA_0 * PA_1 * QD_1 * (-1.0))
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PA_1 * PQ[c0] * (-1.0) + PB_1 * PA_0 * PA_1 * QC_0 * (-1.0))
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PA_1 * PQ[c1] * (-1.0) + PB_1 * PA_0 * PA_1 * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * PQ[d0] * (-1.0) + PB_1 * PA_0 * PA_1 * QD_0 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * PQ[d1] * (-1.0) + PB_1 * PA_0 * PA_1 * QD_1 * (-1.0))
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] + PA_0 * PA_1 * PQ[c0] * QC_1 + PA_0 * PA_1 * PQ[c1] * QC_0)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[d0] + PA_0 * PA_1 * PQ[c0] * QD_0 + PA_0 * PA_1 * PQ[d0] * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * PQ[c0] * PQ[d1] + PA_0 * PA_1 * PQ[c0] * QD_1 + PA_0 * PA_1 * PQ[d1] * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * PQ[c1] * PQ[d0] + PA_0 * PA_1 * PQ[c1] * QD_0 + PA_0 * PA_1 * PQ[d0] * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * PQ[c1] * PQ[d1] + PA_0 * PA_1 * PQ[c1] * QD_1 + PA_0 * PA_1 * PQ[d1] * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * PQ[d0] * PQ[d1] + PA_0 * PA_1 * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[d1] * QD_0)
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PQ[c0] * (-1.0) + PB_0 * PB_1 * PA_0 * QC_0 * (-1.0))
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_0 * PQ[c1] * (-1.0) + PB_0 * PB_1 * PA_0 * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * PQ[d0] * (-1.0) + PB_0 * PB_1 * PA_0 * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * PQ[d1] * (-1.0) + PB_0 * PB_1 * PA_0 * QD_1 * (-1.0))
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] + PB_0 * PA_0 * PQ[c0] * QC_1 + PB_0 * PA_0 * PQ[c1] * QC_0)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[d0] + PB_0 * PA_0 * PQ[c0] * QD_0 + PB_0 * PA_0 * PQ[d0] * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * PQ[c0] * PQ[d1] + PB_0 * PA_0 * PQ[c0] * QD_1 + PB_0 * PA_0 * PQ[d1] * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * PQ[c1] * PQ[d0] + PB_0 * PA_0 * PQ[c1] * QD_0 + PB_0 * PA_0 * PQ[d0] * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * PQ[c1] * PQ[d1] + PB_0 * PA_0 * PQ[c1] * QD_1 + PB_0 * PA_0 * PQ[d1] * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * PQ[d0] * PQ[d1] + PB_0 * PA_0 * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[d1] * QD_0)
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] + PB_1 * PA_0 * PQ[c0] * QC_1 + PB_1 * PA_0 * PQ[c1] * QC_0)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[d0] + PB_1 * PA_0 * PQ[c0] * QD_0 + PB_1 * PA_0 * PQ[d0] * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * PQ[c0] * PQ[d1] + PB_1 * PA_0 * PQ[c0] * QD_1 + PB_1 * PA_0 * PQ[d1] * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * PQ[c1] * PQ[d0] + PB_1 * PA_0 * PQ[c1] * QD_0 + PB_1 * PA_0 * PQ[d0] * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * PQ[c1] * PQ[d1] + PB_1 * PA_0 * PQ[c1] * QD_1 + PB_1 * PA_0 * PQ[d1] * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * PQ[d0] * PQ[d1] + PB_1 * PA_0 * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[d1] * QD_0)
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_1 * PQ[c0] * (-1.0) + PB_0 * PB_1 * PA_1 * QC_0 * (-1.0))
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_1 * PQ[c1] * (-1.0) + PB_0 * PB_1 * PA_1 * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * PQ[d0] * (-1.0) + PB_0 * PB_1 * PA_1 * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * PQ[d1] * (-1.0) + PB_0 * PB_1 * PA_1 * QD_1 * (-1.0))
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] + PB_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PA_1 * PQ[c1] * QC_0)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[d0] + PB_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PA_1 * PQ[d0] * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * PQ[c0] * PQ[d1] + PB_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PA_1 * PQ[d1] * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * PQ[c1] * PQ[d0] + PB_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PA_1 * PQ[d0] * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * PQ[c1] * PQ[d1] + PB_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PA_1 * PQ[d1] * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * PQ[d0] * PQ[d1] + PB_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[d1] * QD_0)
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] + PB_1 * PA_1 * PQ[c0] * QC_1 + PB_1 * PA_1 * PQ[c1] * QC_0)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[d0] + PB_1 * PA_1 * PQ[c0] * QD_0 + PB_1 * PA_1 * PQ[d0] * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * PQ[c0] * PQ[d1] + PB_1 * PA_1 * PQ[c0] * QD_1 + PB_1 * PA_1 * PQ[d1] * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * PQ[c1] * PQ[d0] + PB_1 * PA_1 * PQ[c1] * QD_0 + PB_1 * PA_1 * PQ[d0] * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * PQ[c1] * PQ[d1] + PB_1 * PA_1 * PQ[c1] * QD_1 + PB_1 * PA_1 * PQ[d1] * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * PQ[d0] * PQ[d1] + PB_1 * PA_1 * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[d1] * QD_0)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PQ[c1] * QC_0)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PQ[d0] * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PQ[d0] * QC_1)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PQ[d1] * QC_0)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PQ[d1] * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[d1] * QD_0)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.25 * S2 * S2 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            );
                        }
                        else if constexpr(part == 3)
                        {
                            return
                            F8_t[2] * (

                                + 0.25 * S2 * inv_S1 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * delta[d0][d1] * (PA_0 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * QD_0 * QC_0 * QC_1 * (-1.0) + PQ[a0] * QD_0 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * QD_1 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a0] * QD_0 * QD_1 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a0] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * QD_0 * QC_0 * QC_1 * (-1.0) + PQ[a1] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a1] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a1] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * QD_0 * QC_0 * QC_1 * (-1.0) + PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * QD_0 * QC_0 * QC_1 * (-1.0) + PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * QD_0 * QD_1 * QC_1 * 2.0 + PQ[c1] * QD_0 * QD_1 * QC_0 * 2.0 + PQ[d0] * QD_1 * QC_0 * QC_1 * 2.0 + PQ[d1] * QD_0 * QC_0 * QC_1 * 2.0)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.25 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[a1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] + PB_1 * PA_1 * PQ[a0] * PQ[b0] + PA_0 * PA_1 * PQ[b0] * PQ[b1])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PQ[a1] * QC_0 + PB_0 * PA_1 * PQ[a0] * QC_0 + PA_0 * PA_1 * PQ[b0] * QC_0)
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PQ[a1] * QC_1 + PB_0 * PA_1 * PQ[a0] * QC_1 + PA_0 * PA_1 * PQ[b0] * QC_1)
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * QD_0 + PB_0 * PA_1 * PQ[a0] * QD_0 + PA_0 * PA_1 * PQ[b0] * QD_0)
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * QD_1 + PB_0 * PA_1 * PQ[a0] * QD_1 + PA_0 * PA_1 * PQ[b0] * QD_1)
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PQ[a1] * QC_0 + PB_1 * PA_1 * PQ[a0] * QC_0 + PA_0 * PA_1 * PQ[b1] * QC_0)
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PQ[a1] * QC_1 + PB_1 * PA_1 * PQ[a0] * QC_1 + PA_0 * PA_1 * PQ[b1] * QC_1)
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * QD_0 + PB_1 * PA_1 * PQ[a0] * QD_0 + PA_0 * PA_1 * PQ[b1] * QD_0)
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * QD_1 + PB_1 * PA_1 * PQ[a0] * QD_1 + PA_0 * PA_1 * PQ[b1] * QD_1)
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c0] * QC_1 + PA_0 * PA_1 * PQ[c1] * QC_0 + PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[c0] * QD_0 + PA_0 * PA_1 * PQ[d0] * QC_0 + PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[c0] * QD_1 + PA_0 * PA_1 * PQ[d1] * QC_0 + PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QD_0 + PA_0 * PA_1 * PQ[d0] * QC_1 + PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QD_1 + PA_0 * PA_1 * PQ[d1] * QC_1 + PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[d1] * QD_0 + PA_0 * PA_1 * QD_0 * QD_1)
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PA_1 * QC_0 * QC_1)
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PA_1 * QD_0 * QC_0)
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PA_1 * QD_1 * QC_0)
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PA_1 * QD_0 * QC_1)
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PA_1 * QD_1 * QC_1)
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PA_1 * QD_0 * QD_1)
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * QC_0 + PB_0 * PA_0 * PQ[b1] * QC_0 + PB_1 * PA_0 * PQ[b0] * QC_0)
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a0] * QC_1 + PB_0 * PA_0 * PQ[b1] * QC_1 + PB_1 * PA_0 * PQ[b0] * QC_1)
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * QD_0 + PB_0 * PA_0 * PQ[b1] * QD_0 + PB_1 * PA_0 * PQ[b0] * QD_0)
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * QD_1 + PB_0 * PA_0 * PQ[b1] * QD_1 + PB_1 * PA_0 * PQ[b0] * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c0] * QC_1 + PB_0 * PA_0 * PQ[c1] * QC_0 + PB_0 * PA_0 * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[c0] * QD_0 + PB_0 * PA_0 * PQ[d0] * QC_0 + PB_0 * PA_0 * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[c0] * QD_1 + PB_0 * PA_0 * PQ[d1] * QC_0 + PB_0 * PA_0 * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QD_0 + PB_0 * PA_0 * PQ[d0] * QC_1 + PB_0 * PA_0 * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QD_1 + PB_0 * PA_0 * PQ[d1] * QC_1 + PB_0 * PA_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[d1] * QD_0 + PB_0 * PA_0 * QD_0 * QD_1)
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PA_0 * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PA_0 * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PA_0 * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PA_0 * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PA_0 * QD_1 * QC_1)
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PA_0 * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c0] * QC_1 + PB_1 * PA_0 * PQ[c1] * QC_0 + PB_1 * PA_0 * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[c0] * QD_0 + PB_1 * PA_0 * PQ[d0] * QC_0 + PB_1 * PA_0 * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[c0] * QD_1 + PB_1 * PA_0 * PQ[d1] * QC_0 + PB_1 * PA_0 * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QD_0 + PB_1 * PA_0 * PQ[d0] * QC_1 + PB_1 * PA_0 * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QD_1 + PB_1 * PA_0 * PQ[d1] * QC_1 + PB_1 * PA_0 * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[d1] * QD_0 + PB_1 * PA_0 * QD_0 * QD_1)
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PA_0 * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PA_0 * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PA_0 * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PA_0 * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PA_0 * QD_1 * QC_1)
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PA_0 * QD_0 * QD_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a1] * QC_0 + PB_0 * PA_1 * PQ[b1] * QC_0 + PB_1 * PA_1 * PQ[b0] * QC_0)
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a1] * QC_1 + PB_0 * PA_1 * PQ[b1] * QC_1 + PB_1 * PA_1 * PQ[b0] * QC_1)
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * QD_0 + PB_0 * PA_1 * PQ[b1] * QD_0 + PB_1 * PA_1 * PQ[b0] * QD_0)
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * QD_1 + PB_0 * PA_1 * PQ[b1] * QD_1 + PB_1 * PA_1 * PQ[b0] * QD_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PA_1 * PQ[c1] * QC_0 + PB_0 * PA_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PA_1 * PQ[d0] * QC_0 + PB_0 * PA_1 * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PA_1 * PQ[d1] * QC_0 + PB_0 * PA_1 * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PA_1 * PQ[d0] * QC_1 + PB_0 * PA_1 * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PA_1 * PQ[d1] * QC_1 + PB_0 * PA_1 * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[d1] * QD_0 + PB_0 * PA_1 * QD_0 * QD_1)
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PA_1 * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PA_1 * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PA_1 * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PA_1 * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PA_1 * QD_1 * QC_1)
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PA_1 * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c0] * QC_1 + PB_1 * PA_1 * PQ[c1] * QC_0 + PB_1 * PA_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[c0] * QD_0 + PB_1 * PA_1 * PQ[d0] * QC_0 + PB_1 * PA_1 * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[c0] * QD_1 + PB_1 * PA_1 * PQ[d1] * QC_0 + PB_1 * PA_1 * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QD_0 + PB_1 * PA_1 * PQ[d0] * QC_1 + PB_1 * PA_1 * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QD_1 + PB_1 * PA_1 * PQ[d1] * QC_1 + PB_1 * PA_1 * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[d1] * QD_0 + PB_1 * PA_1 * QD_0 * QD_1)
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PA_1 * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PA_1 * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PA_1 * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PA_1 * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PA_1 * QD_1 * QC_1)
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PA_1 * QD_0 * QD_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PQ[c1] * QC_0 + PB_0 * PB_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PQ[d0] * QC_0 + PB_0 * PB_1 * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PQ[d0] * QC_1 + PB_0 * PB_1 * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PQ[d1] * QC_0 + PB_0 * PB_1 * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PQ[d1] * QC_1 + PB_0 * PB_1 * QD_1 * QC_1)
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PB_1 * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PB_1 * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PB_1 * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PB_1 * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PB_1 * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[d1] * QD_0 + PB_0 * PB_1 * QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PB_1 * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.5 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.5 * S1 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QC_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QC_0 * QC_1 * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QC_0 * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_1 * QC_0 * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QC_1 * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_1 * QC_1 * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1 * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PA_0 * PA_1 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PA_0 * PA_1 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][d1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PB_1 * PA_0 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PB_1 * PA_0 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_0 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PA_0 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PA_0 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PA_0 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PA_0 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PA_0 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_1 * PA_0 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_1 * PA_0 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_1 * PA_0 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_1 * PA_0 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_1 * PA_0 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][d1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PB_1 * PA_1 * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PB_1 * PA_1 * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_1 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_1 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_1 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_1 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_1 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PB_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PB_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PB_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PB_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PB_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.5 * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * (PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.5 * S2 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QC_0 * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QC_0 * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QC_0 * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QC_0 * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QC_0 * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QC_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QC_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_1 * QC_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_1 * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_1 * QC_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_1 * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_1 * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_1 * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_1 * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_1 * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_1 * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_1 * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + delta[b1][d1] * (PB_0 * PA_0 * PQ[a1] * QD_0 * QC_0 * QC_1 + PB_0 * PA_1 * PQ[a0] * QD_0 * QC_0 * QC_1 + PA_0 * PA_1 * PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + delta[b1][d0] * (PB_0 * PA_0 * PQ[a1] * QD_1 * QC_0 * QC_1 + PB_0 * PA_1 * PQ[a0] * QD_1 * QC_0 * QC_1 + PA_0 * PA_1 * PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + delta[b1][c1] * (PB_0 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 + PB_0 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 + PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + delta[b1][c0] * (PB_0 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_1 + PB_0 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_1 + PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][d1] * (PB_1 * PA_0 * PQ[a1] * QD_0 * QC_0 * QC_1 + PB_1 * PA_1 * PQ[a0] * QD_0 * QC_0 * QC_1 + PA_0 * PA_1 * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[b0][d0] * (PB_1 * PA_0 * PQ[a1] * QD_1 * QC_0 * QC_1 + PB_1 * PA_1 * PQ[a0] * QD_1 * QC_0 * QC_1 + PA_0 * PA_1 * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[b0][c1] * (PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 + PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 + PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[b0][c0] * (PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_1 + PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_1 + PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][b1] * (PA_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d1] * (PB_0 * PB_1 * PQ[a0] * QD_0 * QC_0 * QC_1 + PB_0 * PA_0 * PQ[b1] * QD_0 * QC_0 * QC_1 + PB_1 * PA_0 * PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d0] * (PB_0 * PB_1 * PQ[a0] * QD_1 * QC_0 * QC_1 + PB_0 * PA_0 * PQ[b1] * QD_1 * QC_0 * QC_1 + PB_1 * PA_0 * PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + delta[a1][c1] * (PB_0 * PB_1 * PQ[a0] * QD_0 * QD_1 * QC_0 + PB_0 * PA_0 * PQ[b1] * QD_0 * QD_1 * QC_0 + PB_1 * PA_0 * PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + delta[a1][c0] * (PB_0 * PB_1 * PQ[a0] * QD_0 * QD_1 * QC_1 + PB_0 * PA_0 * PQ[b1] * QD_0 * QD_1 * QC_1 + PB_1 * PA_0 * PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PA_0 * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PA_0 * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PA_0 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_1 * PA_0 * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_1 * PA_0 * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_1 * PA_0 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d1] * (PB_0 * PB_1 * PQ[a1] * QD_0 * QC_0 * QC_1 + PB_0 * PA_1 * PQ[b1] * QD_0 * QC_0 * QC_1 + PB_1 * PA_1 * PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d0] * (PB_0 * PB_1 * PQ[a1] * QD_1 * QC_0 * QC_1 + PB_0 * PA_1 * PQ[b1] * QD_1 * QC_0 * QC_1 + PB_1 * PA_1 * PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + delta[a0][c1] * (PB_0 * PB_1 * PQ[a1] * QD_0 * QD_1 * QC_0 + PB_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 + PB_1 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + delta[a0][c0] * (PB_0 * PB_1 * PQ[a1] * QD_0 * QD_1 * QC_1 + PB_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_1 + PB_1 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_1 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_1 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_1 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PB_1 * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PB_1 * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PB_1 * PQ[d1] * QD_0 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[2] * (

                                + S1 * S1 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[2] * (

                                + S1 * S2 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )

                            );
                        }
                        else if constexpr(part == 4)
                        {
                            return
                            F8_t[2] * (

                                + S2 * S2 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.0625 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.0625 * inv_S1 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * 4.0
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1])
                                )

                            )
                            +
                            F8_t[2] * (

                                + 0.0625 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            )
                            +
                            F8_t[3] * (

                                0.0625 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (-2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (-1.0)
                                )

                            )
                            +
                            F8_t[3] * (

                                + 0.0625 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (-2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (-1.0)
                                )

                            );
                        }
                        else if constexpr(part == 16)
                        {
                            return
                            F8_t[3] * (

                                + 0.125 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1 * (-1.0) + PA_0 * PQ[a1] + PA_1 * PQ[a0])
                                    + (delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PA_0 * PA_1 * (-1.0))
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0 * (-1.0) + PB_0 * PQ[a0] + PA_0 * PQ[b0])
                                    + (delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * (-1.0))
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * (-1.0) + PB_1 * PQ[a0] + PA_0 * PQ[b1])
                                    + (delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_0 * PQ[c0])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_0 * PQ[c1])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d0])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d1])
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1 * (-1.0) + PB_0 * PQ[a1] + PA_1 * PQ[b0])
                                    + (delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_1 * (-1.0))
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1 * (-1.0) + PB_1 * PQ[a1] + PA_1 * PQ[b1])
                                    + (delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_1 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_1 * PQ[c0])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_1 * PQ[c1])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d0])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d1])
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * (-1.0) + PB_0 * PQ[b1] + PB_1 * PQ[b0])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * PQ[c1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PB_0 * PQ[c0])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PB_0 * PQ[c1])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d0])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PB_1 * PQ[c0])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PB_1 * PQ[c1])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d1])
                                )

                            )
                            +
                            F8_t[3] * (

                                + (-0.125) * S2 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a0] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a0] * QC_1)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * QD_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * QD_1)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a1] * QC_1)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * QD_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * QD_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PQ[b0] * PQ[b1])
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PQ[b0] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PQ[b0] * QC_1)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * QD_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * QD_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PQ[b1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PQ[b1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * QD_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * QC_1 + PQ[c1] * QC_0 + QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * QD_0 + PQ[d0] * QC_0 + QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * QD_1 + PQ[d1] * QC_0 + QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * QD_0 + PQ[d0] * QC_1 + QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * QD_1 + PQ[d1] * QC_1 + QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * QD_1 + PQ[d1] * QD_0 + QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[3] * (

                                + 0.125 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * (-2.0) + PA_0 * PQ[a1] * 2.0 + PA_1 * PQ[a0] * 2.0)
                                    + (delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PA_0 * PQ[a1] + PA_1 * PQ[a0])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0] * (-2.0) + PB_0 * PQ[a0] * 2.0 + PA_0 * PQ[b0] * 2.0)
                                    + (delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] + PA_0 * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b1] * (-2.0) + PB_1 * PQ[a0] * 2.0 + PA_0 * PQ[b1] * 2.0)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a0] + PA_0 * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a0] * PQ[c0] * (-1.0) + PQ[a0] * QC_0 * (-1.0) + PA_0 * PQ[c0] + PA_0 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a0] * PQ[c1] * (-1.0) + PQ[a0] * QC_1 * (-1.0) + PA_0 * PQ[c1] + PA_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * PQ[d0] * (-1.0) + PQ[a0] * QD_0 * (-1.0) + PA_0 * PQ[d0] + PA_0 * QD_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * PQ[d1] * (-1.0) + PQ[a0] * QD_1 * (-1.0) + PA_0 * PQ[d1] + PA_0 * QD_1)
                                    + (delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c1]) * (PA_0 * QC_0)
                                    + (delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c0]) * (PA_0 * QC_1)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][c0]) * (PA_0 * QD_0)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a1][d0] * delta[b0][c1] * delta[b1][c0]) * (PA_0 * QD_1)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0] * (-2.0) + PB_0 * PQ[a1] * 2.0 + PA_1 * PQ[b0] * 2.0)
                                    + (delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a1] + PA_1 * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b1] * (-2.0) + PB_1 * PQ[a1] * 2.0 + PA_1 * PQ[b1] * 2.0)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a1] + PA_1 * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a1] * PQ[c0] * (-1.0) + PQ[a1] * QC_0 * (-1.0) + PA_1 * PQ[c0] + PA_1 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a1] * PQ[c1] * (-1.0) + PQ[a1] * QC_1 * (-1.0) + PA_1 * PQ[c1] + PA_1 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * PQ[d0] * (-1.0) + PQ[a1] * QD_0 * (-1.0) + PA_1 * PQ[d0] + PA_1 * QD_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * PQ[d1] * (-1.0) + PQ[a1] * QD_1 * (-1.0) + PA_1 * PQ[d1] + PA_1 * QD_1)
                                    + (delta[a0][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c1]) * (PA_1 * QC_0)
                                    + (delta[a0][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c0]) * (PA_1 * QC_1)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][c0]) * (PA_1 * QD_0)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[b0][c1] * delta[b1][c0]) * (PA_1 * QD_1)
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PQ[b0] * PQ[c0] * (-1.0) + PQ[b0] * QC_0 * (-1.0) + PB_0 * PQ[c0] + PB_0 * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c1]) * (PB_0 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PQ[b0] * PQ[c1] * (-1.0) + PQ[b0] * QC_1 * (-1.0) + PB_0 * PQ[c1] + PB_0 * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c0]) * (PB_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PQ[b1] * PQ[c0] * (-1.0) + PQ[b1] * QC_0 * (-1.0) + PB_1 * PQ[c0] + PB_1 * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1]) * (PB_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PQ[b1] * PQ[c1] * (-1.0) + PQ[b1] * QC_1 * (-1.0) + PB_1 * PQ[c1] + PB_1 * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0]) * (PB_1 * QC_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PQ[b0] * PQ[b1] * (-2.0) + PB_0 * PQ[b1] * 2.0 + PB_1 * PQ[b0] * 2.0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * PQ[d0] * (-1.0) + PQ[b0] * QD_0 * (-1.0) + PB_0 * PQ[d0] + PB_0 * QD_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * PQ[d1] * (-1.0) + PQ[b0] * QD_1 * (-1.0) + PB_0 * PQ[d1] + PB_0 * QD_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * PQ[d0] * (-1.0) + PQ[b1] * QD_0 * (-1.0) + PB_1 * PQ[d0] + PB_1 * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * PQ[d1] * (-1.0) + PQ[b1] * QD_1 * (-1.0) + PB_1 * PQ[d1] + PB_1 * QD_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * PQ[c1] * (-2.0) + PQ[c0] * QC_1 * (-2.0) + PQ[c1] * QC_0 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * PQ[d0] * (-2.0) + PQ[c0] * QD_0 * (-2.0) + PQ[d0] * QC_0 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * PQ[d1] * (-2.0) + PQ[c0] * QD_1 * (-2.0) + PQ[d1] * QC_0 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (PQ[c0] * QC_1 * (-1.0) + PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * QD_0 * (-1.0) + PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * QD_1 * (-1.0) + PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * PQ[d0] * (-2.0) + PQ[c1] * QD_0 * (-2.0) + PQ[d0] * QC_1 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * PQ[d1] * (-2.0) + PQ[c1] * QD_1 * (-2.0) + PQ[d1] * QC_1 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * QD_0 * (-1.0) + PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * QD_1 * (-1.0) + PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * PQ[d1] * (-2.0) + PQ[d0] * QD_1 * (-2.0) + PQ[d1] * QD_0 * (-2.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (PQ[d0] * QD_1 * (-1.0) + PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PB_0 * PQ[b1] + PB_1 * PQ[b0])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][c0]) * (PB_0 * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b1][c0]) * (PB_0 * QD_1)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0]) * (PB_1 * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0]) * (PB_1 * QD_1)
                                )

                            )
                            +
                            F8_t[3] * (

                                + 0.25 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PQ[a1] + PB_0 * PB_1 * PA_1 * PQ[a0] + PB_0 * PA_0 * PA_1 * PQ[b1] + PB_1 * PA_0 * PA_1 * PQ[b0])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PA_1 * PQ[c0])
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PA_1 * PQ[c1])
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * PQ[d0])
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PA_1 * PQ[d1])
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PA_1 * PQ[c0])
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PA_1 * PQ[c1])
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * PQ[d0])
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PA_1 * PQ[d1])
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PQ[c0])
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_0 * PQ[c1])
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * PQ[d0])
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_0 * PQ[d1])
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_1 * PQ[c0])
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PA_1 * PQ[c1])
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * PQ[d0])
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PA_1 * PQ[d1])
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * PQ[d0] * PQ[d1] * (-1.0))
                                )

                            );
                        }
                        else if constexpr(part == 17)
                        {
                            return
                            F8_t[3] * (

                                + 0.25 * S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * (-2.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * (-2.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * (-2.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * (-2.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * (-2.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * (-2.0))
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * (-1.0) + PB_0 * PA_0 * PQ[a1] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * (-1.0) + PB_0 * PA_1 * PQ[a0] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * (-1.0) + PA_0 * PA_1 * PQ[b0] * QC_0 * (-1.0))
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PQ[a1] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * QC_1 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[a1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[a0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b0] * QD_0 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * QD_1 * (-1.0))
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * (-1.0) + PB_1 * PA_0 * PQ[a1] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * (-1.0) + PB_1 * PA_1 * PQ[a0] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * (-1.0) + PA_0 * PA_1 * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PQ[a1] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c1] * (-1.0) + PA_0 * PA_1 * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[a1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[a0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b1] * QD_1 * (-1.0))
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] + PA_0 * PQ[a1] * PQ[c0] * QC_1 + PA_0 * PQ[a1] * PQ[c1] * QC_0 + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] + PA_1 * PQ[a0] * PQ[c0] * QC_1 + PA_1 * PQ[a0] * PQ[c1] * QC_0)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PA_1 * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] + PA_0 * PQ[a1] * PQ[c0] * QD_0 + PA_0 * PQ[a1] * PQ[d0] * QC_0 + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] + PA_1 * PQ[a0] * PQ[c0] * QD_0 + PA_1 * PQ[a0] * PQ[d0] * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PA_1 * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d1] + PA_0 * PQ[a1] * PQ[c0] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[c0] * PQ[d1] + PA_1 * PQ[a0] * PQ[c0] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QC_0)
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PA_1 * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] + PA_0 * PQ[a1] * PQ[c1] * QD_0 + PA_0 * PQ[a1] * PQ[d0] * QC_1 + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] + PA_1 * PQ[a0] * PQ[c1] * QD_0 + PA_1 * PQ[a0] * PQ[d0] * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PA_1 * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d1] + PA_0 * PQ[a1] * PQ[c1] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[c1] * PQ[d1] + PA_1 * PQ[a0] * PQ[c1] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QC_1)
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PA_1 * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * PQ[d1] + PA_0 * PQ[a1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[d0] * PQ[d1] + PA_1 * PQ[a0] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QD_0)
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * (-1.0) + PB_0 * PB_1 * PQ[a0] * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * (-1.0) + PB_0 * PA_0 * PQ[b1] * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * (-1.0) + PB_1 * PA_0 * PQ[b0] * QC_0 * (-1.0))
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[c1] * (-1.0) + PB_0 * PB_1 * PQ[a0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PQ[b1] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PQ[b0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[a0] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[b1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[b0] * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[a0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[b1] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[b0] * QD_1 * (-1.0))
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] + PB_0 * PQ[a0] * PQ[c0] * QC_1 + PB_0 * PQ[a0] * PQ[c1] * QC_0 + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] + PA_0 * PQ[b0] * PQ[c0] * QC_1 + PA_0 * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PA_0 * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] + PB_0 * PQ[a0] * PQ[c0] * QD_0 + PB_0 * PQ[a0] * PQ[d0] * QC_0 + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] + PA_0 * PQ[b0] * PQ[c0] * QD_0 + PA_0 * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PA_0 * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d1] + PB_0 * PQ[a0] * PQ[c0] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QC_0 + PA_0 * PQ[b0] * PQ[c0] * PQ[d1] + PA_0 * PQ[b0] * PQ[c0] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QC_0)
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PA_0 * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] + PB_0 * PQ[a0] * PQ[c1] * QD_0 + PB_0 * PQ[a0] * PQ[d0] * QC_1 + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] + PA_0 * PQ[b0] * PQ[c1] * QD_0 + PA_0 * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PA_0 * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d1] + PB_0 * PQ[a0] * PQ[c1] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QC_1 + PA_0 * PQ[b0] * PQ[c1] * PQ[d1] + PA_0 * PQ[b0] * PQ[c1] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QC_1)
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PA_0 * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * PQ[d1] + PB_0 * PQ[a0] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QD_0 + PA_0 * PQ[b0] * PQ[d0] * PQ[d1] + PA_0 * PQ[b0] * PQ[d0] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QD_0)
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] + PB_1 * PQ[a0] * PQ[c0] * QC_1 + PB_1 * PQ[a0] * PQ[c1] * QC_0 + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] + PA_0 * PQ[b1] * PQ[c0] * QC_1 + PA_0 * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PA_0 * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] + PB_1 * PQ[a0] * PQ[c0] * QD_0 + PB_1 * PQ[a0] * PQ[d0] * QC_0 + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] + PA_0 * PQ[b1] * PQ[c0] * QD_0 + PA_0 * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PA_0 * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d1] + PB_1 * PQ[a0] * PQ[c0] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QC_0 + PA_0 * PQ[b1] * PQ[c0] * PQ[d1] + PA_0 * PQ[b1] * PQ[c0] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QC_0)
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PA_0 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PA_0 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PA_0 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PA_0 * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] + PB_1 * PQ[a0] * PQ[c1] * QD_0 + PB_1 * PQ[a0] * PQ[d0] * QC_1 + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] + PA_0 * PQ[b1] * PQ[c1] * QD_0 + PA_0 * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PA_0 * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d1] + PB_1 * PQ[a0] * PQ[c1] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QC_1 + PA_0 * PQ[b1] * PQ[c1] * PQ[d1] + PA_0 * PQ[b1] * PQ[c1] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QC_1)
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PA_0 * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * PQ[d1] + PB_1 * PQ[a0] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QD_0 + PA_0 * PQ[b1] * PQ[d0] * PQ[d1] + PA_0 * PQ[b1] * PQ[d0] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QD_0)
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * (-1.0) + PB_0 * PB_1 * PQ[a1] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * (-1.0) + PB_0 * PA_1 * PQ[b1] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * (-1.0) + PB_1 * PA_1 * PQ[b0] * QC_0 * (-1.0))
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a1] * PQ[c1] * (-1.0) + PB_0 * PB_1 * PQ[a1] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c1] * (-1.0) + PB_0 * PA_1 * PQ[b1] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c1] * (-1.0) + PB_1 * PA_1 * PQ[b0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[a1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[b1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[b0] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[a1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[b1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[b0] * QD_1 * (-1.0))
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] + PB_0 * PQ[a1] * PQ[c0] * QC_1 + PB_0 * PQ[a1] * PQ[c1] * QC_0 + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] + PA_1 * PQ[b0] * PQ[c0] * QC_1 + PA_1 * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PA_1 * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] + PB_0 * PQ[a1] * PQ[c0] * QD_0 + PB_0 * PQ[a1] * PQ[d0] * QC_0 + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] + PA_1 * PQ[b0] * PQ[c0] * QD_0 + PA_1 * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PA_1 * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d1] + PB_0 * PQ[a1] * PQ[c0] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QC_0 + PA_1 * PQ[b0] * PQ[c0] * PQ[d1] + PA_1 * PQ[b0] * PQ[c0] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PA_1 * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] + PB_0 * PQ[a1] * PQ[c1] * QD_0 + PB_0 * PQ[a1] * PQ[d0] * QC_1 + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] + PA_1 * PQ[b0] * PQ[c1] * QD_0 + PA_1 * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PA_1 * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d1] + PB_0 * PQ[a1] * PQ[c1] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QC_1 + PA_1 * PQ[b0] * PQ[c1] * PQ[d1] + PA_1 * PQ[b0] * PQ[c1] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PA_1 * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * PQ[d1] + PB_0 * PQ[a1] * PQ[d0] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QD_0 + PA_1 * PQ[b0] * PQ[d0] * PQ[d1] + PA_1 * PQ[b0] * PQ[d0] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] + PB_1 * PQ[a1] * PQ[c0] * QC_1 + PB_1 * PQ[a1] * PQ[c1] * QC_0 + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] + PA_1 * PQ[b1] * PQ[c0] * QC_1 + PA_1 * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PA_1 * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] + PB_1 * PQ[a1] * PQ[c0] * QD_0 + PB_1 * PQ[a1] * PQ[d0] * QC_0 + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] + PA_1 * PQ[b1] * PQ[c0] * QD_0 + PA_1 * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PA_1 * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d1] + PB_1 * PQ[a1] * PQ[c0] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QC_0 + PA_1 * PQ[b1] * PQ[c0] * PQ[d1] + PA_1 * PQ[b1] * PQ[c0] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PA_1 * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PA_1 * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PA_1 * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PA_1 * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] + PB_1 * PQ[a1] * PQ[c1] * QD_0 + PB_1 * PQ[a1] * PQ[d0] * QC_1 + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] + PA_1 * PQ[b1] * PQ[c1] * QD_0 + PA_1 * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PA_1 * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d1] + PB_1 * PQ[a1] * PQ[c1] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QC_1 + PA_1 * PQ[b1] * PQ[c1] * PQ[d1] + PA_1 * PQ[b1] * PQ[c1] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PA_1 * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * PQ[d1] + PB_1 * PQ[a1] * PQ[d0] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QD_0 + PA_1 * PQ[b1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b1] * PQ[d0] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * QD_0 + PA_1 * PQ[c0] * PQ[d0] * QC_1 + PA_1 * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * QD_1 + PA_1 * PQ[c0] * PQ[d1] * QC_1 + PA_1 * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[d0] * QD_1 + PA_1 * PQ[c0] * PQ[d1] * QD_0 + PA_1 * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PB_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] + PB_0 * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PQ[b1] * PQ[c1] * QC_0 + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] + PB_1 * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PB_1 * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] + PB_0 * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PQ[b1] * PQ[d0] * QC_0 + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] + PB_1 * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PB_1 * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d1] + PB_0 * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QC_0 + PB_1 * PQ[b0] * PQ[c0] * PQ[d1] + PB_1 * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PB_1 * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PB_1 * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PB_1 * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PB_1 * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] + PB_0 * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PQ[b1] * PQ[d0] * QC_1 + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] + PB_1 * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PB_1 * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d1] + PB_0 * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QC_1 + PB_1 * PQ[b0] * PQ[c1] * PQ[d1] + PB_1 * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PB_1 * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * PQ[d1] + PB_0 * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QD_0 + PB_1 * PQ[b0] * PQ[d0] * PQ[d1] + PB_1 * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                )

                            )
                            +
                            F8_t[3] * (

                                + (-0.25) * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * delta[d0][d1] * (PQ[a0] * PQ[a1] * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PQ[a0] * PQ[a1] * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PQ[a0] * PQ[a1] * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PQ[a0] * PQ[a1] * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PQ[a0] * PQ[a1] * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PQ[a0] * PQ[a1] * QD_0 * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PQ[a0] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PQ[a0] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PQ[a0] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PQ[a0] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PQ[a0] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PQ[a0] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PQ[a0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PQ[a0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PQ[a0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PQ[a0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PQ[a0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PQ[a0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * QD_0 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * QD_1 * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * QD_0 * QD_1 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PQ[a1] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PQ[a1] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PQ[a1] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PQ[a1] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PQ[a1] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PQ[a1] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PQ[a1] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PQ[a1] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PQ[a1] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PQ[a1] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PQ[a1] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PQ[a1] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[d1] * QD_0 * QC_0 * QC_1)
                                )

                            );
                        }
                        else if constexpr(part == 5)
                        {
                            return
                            F8_t[3] * (

                                + 0.25 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * QC_0)
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[a1] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * QC_1)
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * QD_0)
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * QD_1)
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * PQ[a1] * QC_0 + PA_0 * PQ[a1] * PQ[b1] * QC_0 + PA_1 * PQ[a0] * PQ[b1] * QC_0)
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PQ[a0] * PQ[a1] * QC_1 + PA_0 * PQ[a1] * PQ[b1] * QC_1 + PA_1 * PQ[a0] * PQ[b1] * QC_1)
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * QD_0 + PA_0 * PQ[a1] * PQ[b1] * QD_0 + PA_1 * PQ[a0] * PQ[b1] * QD_0)
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * QD_1 + PA_0 * PQ[a1] * PQ[b1] * QD_1 + PA_1 * PQ[a0] * PQ[b1] * QD_1)
                                    + delta[b0][b1] * delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QC_1 + PA_0 * PQ[a1] * PQ[c1] * QC_0 + PA_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[c0] * QC_1 + PA_1 * PQ[a0] * PQ[c1] * QC_0 + PA_1 * PQ[a0] * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QD_0 + PA_0 * PQ[a1] * PQ[d0] * QC_0 + PA_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[a0] * PQ[c0] * QD_0 + PA_1 * PQ[a0] * PQ[d0] * QC_0 + PA_1 * PQ[a0] * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[c0] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QD_0 + PA_0 * PQ[a1] * PQ[d0] * QC_1 + PA_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[a0] * PQ[c1] * QD_0 + PA_1 * PQ[a0] * PQ[d0] * QC_1 + PA_1 * PQ[a0] * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[c1] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * QD_1 * QC_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[a0] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * QD_0 * QD_1)
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[a0] * QC_0 * QC_1)
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[a0] * QD_0 * QC_0)
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[a0] * QD_1 * QC_0)
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[a0] * QD_0 * QC_1)
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[a0] * QD_1 * QC_1)
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[a0] * QD_0 * QD_1)
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[b1] * QC_0 + PB_1 * PQ[a0] * PQ[b0] * QC_0 + PA_0 * PQ[b0] * PQ[b1] * QC_0)
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[b1] * QC_1 + PB_1 * PQ[a0] * PQ[b0] * QC_1 + PA_0 * PQ[b0] * PQ[b1] * QC_1)
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * QD_0 + PB_1 * PQ[a0] * PQ[b0] * QD_0 + PA_0 * PQ[b0] * PQ[b1] * QD_0)
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * QD_1 + PB_1 * PQ[a0] * PQ[b0] * QD_1 + PA_0 * PQ[b0] * PQ[b1] * QD_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QC_1 + PB_0 * PQ[a0] * PQ[c1] * QC_0 + PB_0 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[c0] * QC_1 + PA_0 * PQ[b0] * PQ[c1] * QC_0 + PA_0 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QD_0 + PB_0 * PQ[a0] * PQ[d0] * QC_0 + PB_0 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b0] * PQ[c0] * QD_0 + PA_0 * PQ[b0] * PQ[d0] * QC_0 + PA_0 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QC_0 + PB_0 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b0] * PQ[c0] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QC_0 + PA_0 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QD_0 + PB_0 * PQ[a0] * PQ[d0] * QC_1 + PB_0 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b0] * PQ[c1] * QD_0 + PA_0 * PQ[b0] * PQ[d0] * QC_1 + PA_0 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QC_1 + PB_0 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b0] * PQ[c1] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QC_1 + PA_0 * PQ[b0] * QD_1 * QC_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b0] * PQ[d0] * QD_1 + PA_0 * PQ[b0] * PQ[d1] * QD_0 + PA_0 * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b0] * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b0] * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b0] * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b0] * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b0] * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QC_1 + PB_1 * PQ[a0] * PQ[c1] * QC_0 + PB_1 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b1] * PQ[c0] * QC_1 + PA_0 * PQ[b1] * PQ[c1] * QC_0 + PA_0 * PQ[b1] * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QD_0 + PB_1 * PQ[a0] * PQ[d0] * QC_0 + PB_1 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b1] * PQ[c0] * QD_0 + PA_0 * PQ[b1] * PQ[d0] * QC_0 + PA_0 * PQ[b1] * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QC_0 + PB_1 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b1] * PQ[c0] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QC_0 + PA_0 * PQ[b1] * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QD_0 + PB_1 * PQ[a0] * PQ[d0] * QC_1 + PB_1 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b1] * PQ[c1] * QD_0 + PA_0 * PQ[b1] * PQ[d0] * QC_1 + PA_0 * PQ[b1] * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b1] * PQ[c1] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QC_1 + PA_0 * PQ[b1] * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b1] * PQ[d0] * QD_1 + PA_0 * PQ[b1] * PQ[d1] * QD_0 + PA_0 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PQ[a0] * QC_0 * QC_1 + PA_0 * PQ[b1] * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PQ[a0] * QD_0 * QC_0 + PA_0 * PQ[b1] * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PQ[a0] * QD_1 * QC_0 + PA_0 * PQ[b1] * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PQ[a0] * QD_0 * QC_1 + PA_0 * PQ[b1] * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PQ[a0] * QD_1 * QC_1 + PA_0 * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PQ[a0] * QD_0 * QD_1 + PA_0 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PQ[a0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PQ[a0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[c0] * QD_0 * QC_1 + PA_0 * PQ[c1] * QD_0 * QC_0 + PA_0 * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[c0] * QD_1 * QC_1 + PA_0 * PQ[c1] * QD_1 * QC_0 + PA_0 * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PQ[a0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[c0] * QD_0 * QD_1 + PA_0 * PQ[d0] * QD_1 * QC_0 + PA_0 * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PQ[a0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[c1] * QD_0 * QD_1 + PA_0 * PQ[d0] * QD_1 * QC_1 + PA_0 * PQ[d1] * QD_0 * QC_1)
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * PQ[b1] * QC_0 + PB_1 * PQ[a1] * PQ[b0] * QC_0 + PA_1 * PQ[b0] * PQ[b1] * QC_0)
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PQ[a1] * PQ[b1] * QC_1 + PB_1 * PQ[a1] * PQ[b0] * QC_1 + PA_1 * PQ[b0] * PQ[b1] * QC_1)
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * QD_0 + PB_1 * PQ[a1] * PQ[b0] * QD_0 + PA_1 * PQ[b0] * PQ[b1] * QD_0)
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * QD_1 + PB_1 * PQ[a1] * PQ[b0] * QD_1 + PA_1 * PQ[b0] * PQ[b1] * QD_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QC_1 + PB_0 * PQ[a1] * PQ[c1] * QC_0 + PB_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[c0] * QC_1 + PA_1 * PQ[b0] * PQ[c1] * QC_0 + PA_1 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QD_0 + PB_0 * PQ[a1] * PQ[d0] * QC_0 + PB_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b0] * PQ[c0] * QD_0 + PA_1 * PQ[b0] * PQ[d0] * QC_0 + PA_1 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QC_0 + PB_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b0] * PQ[c0] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QC_0 + PA_1 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QD_0 + PB_0 * PQ[a1] * PQ[d0] * QC_1 + PB_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b0] * PQ[c1] * QD_0 + PA_1 * PQ[b0] * PQ[d0] * QC_1 + PA_1 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QC_1 + PB_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b0] * PQ[c1] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QC_1 + PA_1 * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QD_1 + PB_0 * PQ[a1] * PQ[d1] * QD_0 + PB_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b0] * PQ[d0] * QD_1 + PA_1 * PQ[b0] * PQ[d1] * QD_0 + PA_1 * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b0] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b0] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b0] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b0] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b0] * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QC_1 + PB_1 * PQ[a1] * PQ[c1] * QC_0 + PB_1 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b1] * PQ[c0] * QC_1 + PA_1 * PQ[b1] * PQ[c1] * QC_0 + PA_1 * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QD_0 + PB_1 * PQ[a1] * PQ[d0] * QC_0 + PB_1 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b1] * PQ[c0] * QD_0 + PA_1 * PQ[b1] * PQ[d0] * QC_0 + PA_1 * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QC_0 + PB_1 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b1] * PQ[c0] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QC_0 + PA_1 * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QD_0 + PB_1 * PQ[a1] * PQ[d0] * QC_1 + PB_1 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b1] * PQ[c1] * QD_0 + PA_1 * PQ[b1] * PQ[d0] * QC_1 + PA_1 * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QC_1 + PB_1 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b1] * PQ[c1] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QC_1 + PA_1 * PQ[b1] * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QD_1 + PB_1 * PQ[a1] * PQ[d1] * QD_0 + PB_1 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b1] * PQ[d0] * QD_1 + PA_1 * PQ[b1] * PQ[d1] * QD_0 + PA_1 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PQ[a1] * QC_0 * QC_1 + PA_1 * PQ[b1] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PQ[a1] * QD_0 * QC_0 + PA_1 * PQ[b1] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PQ[a1] * QD_1 * QC_0 + PA_1 * PQ[b1] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PQ[a1] * QD_0 * QC_1 + PA_1 * PQ[b1] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PQ[a1] * QD_1 * QC_1 + PA_1 * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PQ[a1] * QD_0 * QD_1 + PA_1 * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[c0] * QD_0 * QC_1 + PA_1 * PQ[c1] * QD_0 * QC_0 + PA_1 * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[c0] * QD_1 * QC_1 + PA_1 * PQ[c1] * QD_1 * QC_0 + PA_1 * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[c0] * QD_0 * QD_1 + PA_1 * PQ[d0] * QD_1 * QC_0 + PA_1 * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[c1] * QD_0 * QD_1 + PA_1 * PQ[d0] * QD_1 * QC_1 + PA_1 * PQ[d1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PQ[b1] * PQ[c1] * QC_0 + PB_0 * PQ[b1] * QC_0 * QC_1 + PB_1 * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PQ[b0] * PQ[c1] * QC_0 + PB_1 * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PQ[b1] * PQ[d0] * QC_0 + PB_0 * PQ[b1] * QD_0 * QC_0 + PB_1 * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PQ[b0] * PQ[d0] * QC_0 + PB_1 * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PQ[b1] * PQ[d0] * QC_1 + PB_0 * PQ[b1] * QD_0 * QC_1 + PB_1 * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PQ[b0] * PQ[d0] * QC_1 + PB_1 * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d0] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QC_0 + PB_0 * PQ[b1] * QD_1 * QC_0 + PB_1 * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QC_0 + PB_1 * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d0] * (PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QC_1 + PB_0 * PQ[b1] * QD_1 * QC_1 + PB_1 * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QC_1 + PB_1 * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PQ[b1] * QC_0 * QC_1 + PB_1 * PQ[b0] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PQ[b1] * QD_0 * QC_0 + PB_1 * PQ[b0] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PQ[b1] * QD_1 * QC_0 + PB_1 * PQ[b0] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PQ[b1] * QD_0 * QC_1 + PB_1 * PQ[b0] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PQ[b1] * QD_1 * QC_1 + PB_1 * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[c0] * QD_0 * QC_1 + PB_0 * PQ[c1] * QD_0 * QC_0 + PB_0 * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[c0] * QD_1 * QC_1 + PB_0 * PQ[c1] * QD_1 * QC_0 + PB_0 * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[c0] * QD_0 * QD_1 + PB_0 * PQ[d0] * QD_1 * QC_0 + PB_0 * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[c1] * QD_0 * QD_1 + PB_0 * PQ[d0] * QD_1 * QC_1 + PB_0 * PQ[d1] * QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[c0] * QD_0 * QC_1 + PB_1 * PQ[c1] * QD_0 * QC_0 + PB_1 * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[c0] * QD_1 * QC_1 + PB_1 * PQ[c1] * QD_1 * QC_0 + PB_1 * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[c0] * QD_0 * QD_1 + PB_1 * PQ[d0] * QD_1 * QC_0 + PB_1 * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[c1] * QD_0 * QD_1 + PB_1 * PQ[d0] * QD_1 * QC_1 + PB_1 * PQ[d1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PQ[b1] * PQ[d1] * QD_0 + PB_0 * PQ[b1] * QD_0 * QD_1 + PB_1 * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PQ[b0] * PQ[d1] * QD_0 + PB_1 * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-2.0) + PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-2.0) + PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-2.0) + PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-2.0) + PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-2.0) + PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-2.0))
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PQ[b1] * QD_0 * QD_1 + PB_1 * PQ[b0] * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[3] * (

                                + (-0.5) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1])
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0])
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1])
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0])
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1])
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1])
                                )

                            );
                        }
                        else if constexpr(part == 6)
                        {
                            return
                            F8_t[3] * (

                                + 0.5 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QC_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QC_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QC_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QC_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QC_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QC_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QC_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QC_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QC_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QC_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QC_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d1] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d1] * QD_0)
                                    + delta[b1][d1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b1][d0] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b1][c1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b1][c0] * (PB_0 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][d1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b0][d0] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b0][c1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b0][c0] * (PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][d1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PB_1 * PA_0 * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a1][d0] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PB_1 * PA_0 * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a1][c1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_0 * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][c0] * (PB_0 * PB_1 * PA_0 * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_0 * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][d1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PB_1 * PA_1 * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a0][d0] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PB_1 * PA_1 * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a0][c1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][c0] * (PB_0 * PB_1 * PA_1 * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_1 * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_1 * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                )

                            )
                            +
                            F8_t[3] * (

                                + 0.5 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QC_0 * QC_1 * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QC_0 * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_1 * QC_0 * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QC_1 * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_1 * QC_1 * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PA_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PB_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_0 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_0 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][d1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PB_1 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PB_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_0 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_0 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PB_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PB_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PB_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[3] * (

                                + (-0.5) * S2 * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * (PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[3] * (

                                + 0.5 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QC_0 * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QC_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * QC_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + delta[b1][d1] * (PB_0 * PQ[a0] * PQ[a1] * QD_0 * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + delta[b1][d0] * (PB_0 * PQ[a0] * PQ[a1] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + delta[b1][c1] * (PB_0 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + delta[b1][c0] * (PB_0 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][d1] * (PB_1 * PQ[a0] * PQ[a1] * QD_0 * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[b1] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[b0][d0] * (PB_1 * PQ[a0] * PQ[a1] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[b1] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[b0][c1] * (PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[b0][c0] * (PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d1] * (PB_0 * PQ[a0] * PQ[b1] * QD_0 * QC_0 * QC_1 + PB_1 * PQ[a0] * PQ[b0] * QD_0 * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d0] * (PB_0 * PQ[a0] * PQ[b1] * QD_1 * QC_0 * QC_1 + PB_1 * PQ[a0] * PQ[b0] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[a1][c1] * (PB_0 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 + PB_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[a1][c0] * (PB_0 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_1 + PB_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_1 * PQ[a0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_1 * PQ[a0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_1 * PQ[a0] * PQ[d1] * QD_0 * QC_0 * QC_1 + PA_0 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_0 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_0 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_0 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d1] * (PB_0 * PQ[a1] * PQ[b1] * QD_0 * QC_0 * QC_1 + PB_1 * PQ[a1] * PQ[b0] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d0] * (PB_0 * PQ[a1] * PQ[b1] * QD_1 * QC_0 * QC_1 + PB_1 * PQ[a1] * PQ[b0] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[a0][c1] * (PB_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 + PB_1 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[a0][c0] * (PB_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_1 + PB_1 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_1 * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_1 * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_1 * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 + PA_1 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PA_1 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PA_1 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PA_1 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_0 * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_0 * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_0 * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 + PB_1 * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PB_1 * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PB_1 * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PB_1 * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                )

                            );
                        }
                        else if constexpr(part == 7)
                        {
                            return
                            F8_t[3] * (

                                + S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                )

                            )
                            +
                            F8_t[3] * (

                                + S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[3] * (

                                + S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )

                            )
                            +
                            F8_t[3] * (

                                + S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[4] * (

                                + 0.5 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][d1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PB_1 * PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PB_1 * PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][d1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PB_1 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PB_1 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[5] * (

                                + (-0.5) * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + delta[b1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * QC_1)
                                    + delta[b1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[b1][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0)
                                    + delta[b1][c0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_1)
                                    + delta[b0][d1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QC_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QC_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * QC_1)
                                    + delta[b0][d0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * QC_1)
                                    + delta[b0][c1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0)
                                    + delta[b0][c0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_1)
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][d1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * QC_1)
                                    + delta[a1][d0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][c1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0)
                                    + delta[a1][c0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][d1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * QC_1)
                                    + delta[a0][d0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][c1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0)
                                    + delta[a0][c0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[5] * (

                                + S1 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )

                            )
                            +
                            F8_t[5] * (

                                + S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                )

                            )
                            +
                            F8_t[5] * (

                                + S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )

                            );
                        }
                        else if constexpr(part == 8)
                        {
                            return
                            F8_t[4] * (

                                0.125 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PQ[a1] * (-1.0) + PA_1 * PQ[a0] * (-1.0) + PQ[a0] * PQ[a1])
                                    + (delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PA_0 * PQ[a1] * (-1.0) + PA_1 * PQ[a0] * (-1.0))
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * (-1.0) + PA_0 * PQ[b0] * (-1.0) + PQ[a0] * PQ[b0])
                                    + (delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * (-1.0) + PA_0 * PQ[b0] * (-1.0))
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * (-1.0) + PA_0 * PQ[b1] * (-1.0) + PQ[a0] * PQ[b1])
                                    + (delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a0] * (-1.0) + PA_0 * PQ[b1] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_0 * PQ[c0] * (-1.0) + PQ[a0] * PQ[c0])
                                    + (delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c1]) * (PA_0 * PQ[c0] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_0 * PQ[c1] * (-1.0) + PQ[a0] * PQ[c1])
                                    + (delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c0]) * (PA_0 * PQ[c1] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d0] * (-1.0) + PQ[a0] * PQ[d0])
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][c0]) * (PA_0 * PQ[d0] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_0 * PQ[d1] * (-1.0) + PQ[a0] * PQ[d1])
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a1][d0] * delta[b0][c1] * delta[b1][c0]) * (PA_0 * PQ[d1] * (-1.0))
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * (-1.0) + PA_1 * PQ[b0] * (-1.0) + PQ[a1] * PQ[b0])
                                    + (delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a1] * (-1.0) + PA_1 * PQ[b0] * (-1.0))
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PQ[a1] * (-1.0) + PA_1 * PQ[b1] * (-1.0) + PQ[a1] * PQ[b1])
                                    + (delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a1] * (-1.0) + PA_1 * PQ[b1] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PA_1 * PQ[c0] * (-1.0) + PQ[a1] * PQ[c0])
                                    + (delta[a0][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c1]) * (PA_1 * PQ[c0] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PA_1 * PQ[c1] * (-1.0) + PQ[a1] * PQ[c1])
                                    + (delta[a0][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c0]) * (PA_1 * PQ[c1] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d0] * (-1.0) + PQ[a1] * PQ[d0])
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][c0]) * (PA_1 * PQ[d0] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PA_1 * PQ[d1] * (-1.0) + PQ[a1] * PQ[d1])
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[b0][c1] * delta[b1][c0]) * (PA_1 * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[b1] * (-1.0) + PB_1 * PQ[b0] * (-1.0) + PQ[b0] * PQ[b1])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PB_0 * PQ[b1] * (-1.0) + PB_1 * PQ[b0] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PB_0 * PQ[c0] * (-1.0) + PQ[b0] * PQ[c0])
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c1]) * (PB_0 * PQ[c0] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PB_0 * PQ[c1] * (-1.0) + PQ[b0] * PQ[c1])
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c0]) * (PB_0 * PQ[c1] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d0] * (-1.0) + PQ[b0] * PQ[d0])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][c0]) * (PB_0 * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PB_0 * PQ[d1] * (-1.0) + PQ[b0] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b1][c0]) * (PB_0 * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PB_1 * PQ[c0] * (-1.0) + PQ[b1] * PQ[c0])
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1]) * (PB_1 * PQ[c0] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PB_1 * PQ[c1] * (-1.0) + PQ[b1] * PQ[c1])
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0]) * (PB_1 * PQ[c1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d0] * (-1.0) + PQ[b1] * PQ[d0])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0]) * (PB_1 * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PB_1 * PQ[d1] * (-1.0) + PQ[b1] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0]) * (PB_1 * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * PQ[c1] * 2.0)
                                    + (delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (PQ[c0] * PQ[c1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * PQ[d0] * 2.0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * PQ[d1] * 2.0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * PQ[d0] * 2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * PQ[d1] * 2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * PQ[d1] * 2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[4] * (

                                + 0.125 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a0] * PQ[c0] + PQ[a0] * QC_0)
                                    + (delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c1]) * (PQ[a0] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a0] * PQ[c1] + PQ[a0] * QC_1)
                                    + (delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c0]) * (PQ[a0] * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0]) * (PQ[a1] * PQ[c0] + PQ[a1] * QC_0)
                                    + (delta[a0][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c1]) * (PQ[a1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0]) * (PQ[a1] * PQ[c1] + PQ[a1] * QC_1)
                                    + (delta[a0][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c0]) * (PQ[a1] * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0]) * (PQ[b0] * PQ[c0] + PQ[b0] * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c1]) * (PQ[b0] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0]) * (PQ[b0] * PQ[c1] + PQ[b0] * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c0]) * (PQ[b0] * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0]) * (PQ[b1] * PQ[c0] + PQ[b1] * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1]) * (PQ[b1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0]) * (PQ[b1] * PQ[c1] + PQ[b1] * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0]) * (PQ[b1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (PQ[c0] * PQ[c1] + PQ[c0] * QC_1 + PQ[c1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (PQ[c0] * QC_1 + PQ[c1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (PQ[c0] * PQ[d0] + PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (PQ[c1] * PQ[d0] + PQ[c1] * QD_0 + PQ[d0] * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * QD_0 + PQ[d0] * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (PQ[c0] * PQ[d1] + PQ[c0] * QD_1 + PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * QD_1 + PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (PQ[c1] * PQ[d1] + PQ[c1] * QD_1 + PQ[d1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * QD_1 + PQ[d1] * QC_1)
                                    + (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * 2.0)
                                    + (delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0] * 2.0)
                                    + (delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b1] * 2.0)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * PQ[d0] + PQ[a0] * QD_0)
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a0] * PQ[d1] + PQ[a0] * QD_1)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * QD_0)
                                    + (delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a1][d0] * delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * QD_1)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0] * 2.0)
                                    + (delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b1] * 2.0)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * PQ[d0] + PQ[a1] * QD_0)
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1]) * (PQ[a1] * PQ[d1] + PQ[a1] * QD_1)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][c0]) * (PQ[a1] * QD_0)
                                    + (delta[a0][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[b0][c1] * delta[b1][c0]) * (PQ[a1] * QD_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PQ[b0] * PQ[b1] * 2.0)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PQ[b0] * PQ[b1])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * PQ[d0] + PQ[b0] * QD_0)
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1]) * (PQ[b0] * PQ[d1] + PQ[b0] * QD_1)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][c0]) * (PQ[b0] * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b1][c0]) * (PQ[b0] * QD_1)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * PQ[d0] + PQ[b1] * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1]) * (PQ[b1] * PQ[d1] + PQ[b1] * QD_1)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0]) * (PQ[b1] * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0]) * (PQ[b1] * QD_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (PQ[d0] * PQ[d1] + PQ[d0] * QD_1 + PQ[d1] * QD_0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (PQ[d0] * QD_1 + PQ[d1] * QD_0)
                                )

                            );
                        }
                        else if constexpr(part == 15)
                        {
                            return
                            F8_t[4] * (

                                + 0.25 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[a1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] + PB_1 * PA_1 * PQ[a0] * PQ[b0] + PA_0 * PA_1 * PQ[b0] * PQ[b1])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PA_0 * PQ[a1] * PQ[c0] + PB_0 * PA_1 * PQ[a0] * PQ[c0] + PA_0 * PA_1 * PQ[b0] * PQ[c0])
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PA_0 * PQ[a1] * PQ[c1] + PB_0 * PA_1 * PQ[a0] * PQ[c1] + PA_0 * PA_1 * PQ[b0] * PQ[c1])
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * PQ[d0] + PB_0 * PA_1 * PQ[a0] * PQ[d0] + PA_0 * PA_1 * PQ[b0] * PQ[d0])
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PA_0 * PQ[a1] * PQ[d1] + PB_0 * PA_1 * PQ[a0] * PQ[d1] + PA_0 * PA_1 * PQ[b0] * PQ[d1])
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PA_0 * PQ[a1] * PQ[c0] + PB_1 * PA_1 * PQ[a0] * PQ[c0] + PA_0 * PA_1 * PQ[b1] * PQ[c0])
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PA_0 * PQ[a1] * PQ[c1] + PB_1 * PA_1 * PQ[a0] * PQ[c1] + PA_0 * PA_1 * PQ[b1] * PQ[c1])
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * PQ[d0] + PB_1 * PA_1 * PQ[a0] * PQ[d0] + PA_0 * PA_1 * PQ[b1] * PQ[d0])
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PA_0 * PQ[a1] * PQ[d1] + PB_1 * PA_1 * PQ[a0] * PQ[d1] + PA_0 * PA_1 * PQ[b1] * PQ[d1])
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[c1])
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PA_1 * PQ[c0] * PQ[c1])
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[d0])
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PA_1 * PQ[c0] * PQ[d0])
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[d1])
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PA_1 * PQ[c0] * PQ[d1])
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[c1] * PQ[d0])
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[c1] * PQ[d0])
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[c1] * PQ[d1])
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[c1] * PQ[d1])
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[d0] * PQ[d1])
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PA_1 * PQ[d0] * PQ[d1])
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[c0] + PB_0 * PA_0 * PQ[b1] * PQ[c0] + PB_1 * PA_0 * PQ[b0] * PQ[c0])
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a0] * PQ[c1] + PB_0 * PA_0 * PQ[b1] * PQ[c1] + PB_1 * PA_0 * PQ[b0] * PQ[c1])
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * PQ[d0] + PB_0 * PA_0 * PQ[b1] * PQ[d0] + PB_1 * PA_0 * PQ[b0] * PQ[d0])
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a0] * PQ[d1] + PB_0 * PA_0 * PQ[b1] * PQ[d1] + PB_1 * PA_0 * PQ[b0] * PQ[d1])
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[c1])
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PA_0 * PQ[c0] * PQ[c1])
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[d0])
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PA_0 * PQ[c0] * PQ[d0])
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[d1])
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PA_0 * PQ[c0] * PQ[d1])
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[c1] * PQ[d0])
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[c1] * PQ[d0])
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[c1] * PQ[d1])
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[c1] * PQ[d1])
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[d0] * PQ[d1])
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PA_0 * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[c1])
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PA_0 * PQ[c0] * PQ[c1])
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[d0])
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PA_0 * PQ[c0] * PQ[d0])
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[d1])
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PA_0 * PQ[c0] * PQ[d1])
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[c1] * PQ[d0])
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[c1] * PQ[d0])
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[c1] * PQ[d1])
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[c1] * PQ[d1])
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[d0] * PQ[d1])
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PA_0 * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PQ[a1] * PQ[c0] + PB_0 * PA_1 * PQ[b1] * PQ[c0] + PB_1 * PA_1 * PQ[b0] * PQ[c0])
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PB_1 * PQ[a1] * PQ[c1] + PB_0 * PA_1 * PQ[b1] * PQ[c1] + PB_1 * PA_1 * PQ[b0] * PQ[c1])
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * PQ[d0] + PB_0 * PA_1 * PQ[b1] * PQ[d0] + PB_1 * PA_1 * PQ[b0] * PQ[d0])
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PB_1 * PQ[a1] * PQ[d1] + PB_0 * PA_1 * PQ[b1] * PQ[d1] + PB_1 * PA_1 * PQ[b0] * PQ[d1])
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[c1])
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PA_1 * PQ[c0] * PQ[c1])
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PA_1 * PQ[c0] * PQ[d0])
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[d1])
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PA_1 * PQ[c0] * PQ[d1])
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[c1] * PQ[d0])
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[c1] * PQ[d1])
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[d0] * PQ[d1])
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PA_1 * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[c1])
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PA_1 * PQ[c0] * PQ[c1])
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PA_1 * PQ[c0] * PQ[d0])
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[d1])
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PA_1 * PQ[c0] * PQ[d1])
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[c1] * PQ[d0])
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[c1] * PQ[d1])
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[d0] * PQ[d1])
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PA_1 * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[c1])
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[d0])
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[d1])
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PB_1 * PQ[c1] * PQ[d0])
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[c1] * PQ[d1])
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PB_1 * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PB_1 * PQ[c0] * PQ[c1])
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PB_1 * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PB_1 * PQ[c0] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PB_1 * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[4] * (

                                + 0.25 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * (-2.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * (-2.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * (-2.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * (-2.0))
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * QC_0 * (-1.0))
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * QC_1 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * QD_0 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * QD_1 * (-1.0))
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[b0][b1] * delta[d0][d1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] + PQ[a0] * PQ[a1] * PQ[c0] * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * QC_0)
                                    + delta[b0][b1] * delta[c1][d1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] + PQ[a0] * PQ[a1] * PQ[c0] * QD_0 + PQ[a0] * PQ[a1] * PQ[d0] * QC_0)
                                    + delta[b0][b1] * delta[c1][d0] * (PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[c0] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QC_0)
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d1] * (PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] + PQ[a0] * PQ[a1] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[d0] * QC_1)
                                    + delta[b0][b1] * delta[c0][d0] * (PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QC_1)
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][c1] * (PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QD_0)
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[a1][b1] * delta[d0][d1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] + PQ[a0] * PQ[b0] * PQ[c0] * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[a1][b1] * delta[c1][d1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] + PQ[a0] * PQ[b0] * PQ[c0] * QD_0 + PQ[a0] * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[a1][b1] * delta[c1][d0] * (PB_0 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] + PQ[a0] * PQ[b0] * PQ[c0] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QC_0)
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d1] * (PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] + PQ[a0] * PQ[b0] * PQ[c1] * QD_0 + PQ[a0] * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[a1][b1] * delta[c0][d0] * (PB_0 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] + PQ[a0] * PQ[b0] * PQ[c1] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QC_1)
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][c1] * (PB_0 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] + PQ[a0] * PQ[b0] * PQ[d0] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QD_0)
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a1][b0] * delta[d0][d1] * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] + PQ[a0] * PQ[b1] * PQ[c0] * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[a1][b0] * delta[c1][d1] * (PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] + PQ[a0] * PQ[b1] * PQ[c0] * QD_0 + PQ[a0] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[a1][b0] * delta[c1][d0] * (PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] + PQ[a0] * PQ[b1] * PQ[c0] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QC_0)
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PQ[a0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PQ[a0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PQ[a0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d1] * (PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] + PQ[a0] * PQ[b1] * PQ[c1] * QD_0 + PQ[a0] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[a1][b0] * delta[c0][d0] * (PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] + PQ[a0] * PQ[b1] * PQ[c1] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QC_1)
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][c1] * (PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] + PQ[a0] * PQ[b1] * PQ[d0] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QD_0)
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[c1] * PQ[d0] * QC_0 * (-1.0) + PQ[a0] * PQ[c0] * PQ[c1] * QD_0 + PQ[a0] * PQ[c0] * PQ[d0] * QC_1 + PQ[a0] * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[c1] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[c0] * PQ[c1] * QD_1 + PQ[a0] * PQ[c0] * PQ[d1] * QC_1 + PQ[a0] * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PA_0 * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PQ[a0] * PQ[c0] * PQ[d0] * QD_1 + PQ[a0] * PQ[c0] * PQ[d1] * QD_0 + PQ[a0] * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PA_0 * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[a0][b1] * delta[d0][d1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] + PQ[a1] * PQ[b0] * PQ[c0] * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * QC_0)
                                    + delta[a0][b1] * delta[c1][d1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] + PQ[a1] * PQ[b0] * PQ[c0] * QD_0 + PQ[a1] * PQ[b0] * PQ[d0] * QC_0)
                                    + delta[a0][b1] * delta[c1][d0] * (PB_0 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] + PQ[a1] * PQ[b0] * PQ[c0] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d1] * (PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] + PQ[a1] * PQ[b0] * PQ[c1] * QD_0 + PQ[a1] * PQ[b0] * PQ[d0] * QC_1)
                                    + delta[a0][b1] * delta[c0][d0] * (PB_0 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] + PQ[a1] * PQ[b0] * PQ[c1] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][c1] * (PB_0 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] + PQ[a1] * PQ[b0] * PQ[d0] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a0][b0] * delta[d0][d1] * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] + PQ[a1] * PQ[b1] * PQ[c0] * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[a0][b0] * delta[c1][d1] * (PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] + PQ[a1] * PQ[b1] * PQ[c0] * QD_0 + PQ[a1] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[a0][b0] * delta[c1][d0] * (PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] + PQ[a1] * PQ[b1] * PQ[c0] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d1] * (PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] + PQ[a1] * PQ[b1] * PQ[c1] * QD_0 + PQ[a1] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[a0][b0] * delta[c0][d0] * (PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] + PQ[a1] * PQ[b1] * PQ[c1] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][c1] * (PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] + PQ[a1] * PQ[b1] * PQ[d0] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[c1] * PQ[d0] * QC_0 * (-1.0) + PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PQ[a1] * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[c1] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PQ[a1] * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PA_1 * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PQ[a1] * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PA_1 * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PQ[a1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] + PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[a0][a1] * delta[c1][d1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] + PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PQ[b0] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[a0][a1] * delta[c1][d0] * (PB_0 * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] + PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QC_0)
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d1] * (PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] + PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PQ[b0] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[a0][a1] * delta[c0][d0] * (PB_0 * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] + PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QC_1)
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][c1] * (PB_0 * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] + PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QD_0)
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[c1] * PQ[d0] * QC_0 * (-1.0) + PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[c1] * PQ[d1] * QC_0 * (-1.0) + PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PB_0 * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PB_0 * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[c1] * PQ[d0] * QC_0 * (-1.0) + PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[c1] * PQ[d1] * QC_0 * (-1.0) + PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PB_1 * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PB_1 * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * 2.0 + PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * 2.0 + PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * 2.0 + PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * 2.0)
                                )

                            )
                            +
                            F8_t[4] * (

                                + 0.25 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * QC_0)
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * QC_1)
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * QC_0)
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * QC_1)
                                    + delta[b0][b1] * delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * QC_0 + PQ[a0] * PQ[a1] * QC_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_0 + PQ[a0] * PQ[a1] * PQ[d0] * QC_0 + PQ[a0] * PQ[a1] * QD_0 * QC_0)
                                    + delta[b0][b1] * delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[d0] * QC_1 + PQ[a0] * PQ[a1] * QD_0 * QC_1)
                                    + delta[b0][b1] * delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QC_0 + PQ[a0] * PQ[a1] * QD_1 * QC_0)
                                    + delta[b0][b1] * delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QC_1 + PQ[a0] * PQ[a1] * QD_1 * QC_1)
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PQ[a0] * PQ[a1] * QC_0 * QC_1)
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * QD_0 * QC_0)
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * QD_1 * QC_0)
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * QD_0 * QC_1)
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * QD_1 * QC_1)
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * QC_0)
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * QC_1)
                                    + delta[a1][b1] * delta[d0][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * QC_0 + PQ[a0] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_0 + PQ[a0] * PQ[b0] * PQ[d0] * QC_0 + PQ[a0] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a1][b1] * delta[c0][d1] * (PQ[a0] * PQ[b0] * PQ[c1] * QD_0 + PQ[a0] * PQ[b0] * PQ[d0] * QC_1 + PQ[a0] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a1][b1] * delta[c1][d0] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QC_0 + PQ[a0] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a1][b1] * delta[c0][d0] * (PQ[a0] * PQ[b0] * PQ[c1] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QC_1 + PQ[a0] * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PQ[a0] * PQ[b0] * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a1][b0] * delta[d0][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * QC_0 + PQ[a0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_0 + PQ[a0] * PQ[b1] * PQ[d0] * QC_0 + PQ[a0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a1][b0] * delta[c0][d1] * (PQ[a0] * PQ[b1] * PQ[c1] * QD_0 + PQ[a0] * PQ[b1] * PQ[d0] * QC_1 + PQ[a0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a1][b0] * delta[c1][d0] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QC_0 + PQ[a0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a1][b0] * delta[c0][d0] * (PQ[a0] * PQ[b1] * PQ[c1] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QC_1 + PQ[a0] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PQ[a0] * PQ[b1] * QC_0 * QC_1)
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * QD_0 * QC_0)
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * QD_1 * QC_0)
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * QD_0 * QC_1)
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_0 * QC_1 + PQ[a0] * PQ[c1] * QD_0 * QC_0 + PQ[a0] * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_1 * QC_1 + PQ[a0] * PQ[c1] * QD_1 * QC_0 + PQ[a0] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * QD_0 * QD_1 + PQ[a0] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[d1] * QD_0 * QC_1)
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * QC_0)
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * QC_1)
                                    + delta[a0][b1] * delta[d0][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * QC_0 + PQ[a1] * PQ[b0] * QC_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_0 + PQ[a1] * PQ[b0] * PQ[d0] * QC_0 + PQ[a1] * PQ[b0] * QD_0 * QC_0)
                                    + delta[a0][b1] * delta[c0][d1] * (PQ[a1] * PQ[b0] * PQ[c1] * QD_0 + PQ[a1] * PQ[b0] * PQ[d0] * QC_1 + PQ[a1] * PQ[b0] * QD_0 * QC_1)
                                    + delta[a0][b1] * delta[c1][d0] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QC_0 + PQ[a1] * PQ[b0] * QD_1 * QC_0)
                                    + delta[a0][b1] * delta[c0][d0] * (PQ[a1] * PQ[b0] * PQ[c1] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QC_1 + PQ[a1] * PQ[b0] * QD_1 * QC_1)
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PQ[a1] * PQ[b0] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * QD_1 * QC_1)
                                    + delta[a0][b0] * delta[d0][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * QC_0 + PQ[a1] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_0 + PQ[a1] * PQ[b1] * PQ[d0] * QC_0 + PQ[a1] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][b0] * delta[c0][d1] * (PQ[a1] * PQ[b1] * PQ[c1] * QD_0 + PQ[a1] * PQ[b1] * PQ[d0] * QC_1 + PQ[a1] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][b0] * delta[c1][d0] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QC_0 + PQ[a1] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][b0] * delta[c0][d0] * (PQ[a1] * PQ[b1] * PQ[c1] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QC_1 + PQ[a1] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PQ[a1] * PQ[b1] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_0 * QC_1 + PQ[a1] * PQ[c1] * QD_0 * QC_0 + PQ[a1] * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_1 * QC_1 + PQ[a1] * PQ[c1] * QD_1 * QC_0 + PQ[a1] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * QD_0 * QD_1 + PQ[a1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[d1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[d0][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * QC_0 + PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PQ[b0] * PQ[b1] * PQ[d0] * QC_0 + PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[a0][a1] * delta[c0][d1] * (PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PQ[b0] * PQ[b1] * PQ[d0] * QC_1 + PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[a0][a1] * delta[c1][d0] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QC_0 + PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[a0][a1] * delta[c0][d0] * (PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QC_1 + PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_0 * QC_1 + PQ[b0] * PQ[c1] * QD_0 * QC_0 + PQ[b0] * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_1 * QC_1 + PQ[b0] * PQ[c1] * QD_1 * QC_0 + PQ[b0] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * QD_0 * QD_1 + PQ[b0] * PQ[d0] * QD_1 * QC_0 + PQ[b0] * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * PQ[c1] * QD_0 * QD_1 + PQ[b0] * PQ[d0] * QD_1 * QC_1 + PQ[b0] * PQ[d1] * QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_0 * QC_1 + PQ[b1] * PQ[c1] * QD_0 * QC_0 + PQ[b1] * PQ[d0] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_1 * QC_1 + PQ[b1] * PQ[c1] * QD_1 * QC_0 + PQ[b1] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * QD_0 * QD_1 + PQ[b1] * PQ[d0] * QD_1 * QC_0 + PQ[b1] * PQ[d1] * QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * PQ[c1] * QD_0 * QD_1 + PQ[b1] * PQ[d0] * QD_1 * QC_1 + PQ[b1] * PQ[d1] * QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1])
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * QD_0)
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * QD_1)
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * QD_0)
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * QD_1)
                                    + delta[b0][b1] * delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * QD_0 * QD_1)
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * QD_0 * QD_1)
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * QD_0)
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * QD_1)
                                    + delta[a1][b1] * delta[c0][c1] * (PQ[a0] * PQ[b0] * PQ[d0] * QD_1 + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 + PQ[a0] * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a1][b0] * delta[c0][c1] * (PQ[a0] * PQ[b1] * PQ[d0] * QD_1 + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 + PQ[a0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * QD_0)
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * QD_1)
                                    + delta[a0][b1] * delta[c0][c1] * (PQ[a1] * PQ[b0] * PQ[d0] * QD_1 + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 + PQ[a1] * PQ[b0] * QD_0 * QD_1)
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * QD_0 * QD_1)
                                    + delta[a0][b0] * delta[c0][c1] * (PQ[a1] * PQ[b1] * PQ[d0] * QD_1 + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 + PQ[a1] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * QD_0 * QD_1)
                                    + delta[a0][a1] * delta[c0][c1] * (PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 + PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                )

                            );
                        }
                        else if constexpr(part == 9)
                        {
                            return
                            F8_t[4] * (

                                + 0.5 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QC_1 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QC_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QC_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QC_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QC_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QC_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QC_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0)
                                    + delta[b1][d1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b1][d0] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b1][c1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b1][c0] * (PB_0 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][d1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QC_0 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b0][d0] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QC_0 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b0][c1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b0][c0] * (PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][b1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * QC_0 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * QC_0 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * QC_0 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * QC_0 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a1][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][c0] * (PB_0 * PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_0 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PA_0 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][d1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * QC_0 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QC_0 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a0][d0] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * QC_0 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QC_0 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a0][c1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][c0] * (PB_0 * PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PB_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PB_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                )

                            );
                        }
                        else if constexpr(part == 10)
                        {
                            return
                            F8_t[4] * (

                                + 0.5 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QC_0 * QC_1 * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QC_0 * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * QC_0 * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QC_1 * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * QC_1 * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][d1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][d1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * QC_1 * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_1 * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PA_1 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0 * (-1.0) + PB_1 * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0 + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0 + PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1)
                                )

                            )
                            +
                            F8_t[4] * (

                                + 0.5 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QC_0 * QC_1)
                                    + delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_0)
                                    + delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_0)
                                    + delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_1)
                                    + delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_1)
                                    + delta[b1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QC_0 * QC_1)
                                    + delta[b1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * QC_0 * QC_1)
                                    + delta[b1][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0)
                                    + delta[b1][c0] * (PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][d1] * (PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[b0][d0] * (PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[b0][c1] * (PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[b0][c0] * (PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d1] * (PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][d0] * (PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[a1][c1] * (PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[a1][c0] * (PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d1] * (PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][d0] * (PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * QC_0 * QC_1)
                                    + delta[a0][c1] * (PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0)
                                    + delta[a0][c0] * (PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 + PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1)
                                    + delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1)
                                )

                            )
                            +
                            F8_t[4] * (

                                + S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )

                            );
                        }
                        else if constexpr(part == 11)
                        {
                            return
                            F8_t[4] * (

                                + S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                )

                            )
                            +
                            F8_t[4] * (

                                + S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[4] * (

                                + S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )

                            )
                            +
                            F8_t[4] * (

                                + S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            );
                        }
                        else if constexpr(part == 12)
                        {
                            return
                            F8_t[4] * (

                                + 0.0625 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] * delta[b1][c0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1] * delta[b1][c0])
                                )

                            )
                            +
                            F8_t[5] * (

                                (-0.125) * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[b0][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1])
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a1][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b0])
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0] + delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a1][d1] * delta[b0][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b1])
                                    + (delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c1]) * (PQ[a0] * PQ[c0])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a1][d1] * delta[b0][d0] * delta[b1][c0]) * (PQ[a0] * PQ[c1])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a1][d1] * delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * PQ[d0])
                                    + (delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a1][d0] * delta[b0][b1] * delta[c0][c1] + delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a1][d0] * delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * PQ[d1])
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[b1][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b0])
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0] + delta[a0][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[b0][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b1])
                                    + (delta[a0][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c1]) * (PQ[a1] * PQ[c0])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[b0][d0] * delta[b1][c0]) * (PQ[a1] * PQ[c1])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[b0][c1] * delta[b1][c0]) * (PQ[a1] * PQ[d0])
                                    + (delta[a0][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[b0][b1] * delta[c0][c1] + delta[a0][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[b0][c1] * delta[b1][c0]) * (PQ[a1] * PQ[d1])
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[c0][d0] + delta[a0][d0] * delta[a1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[c0][c1]) * (PQ[b0] * PQ[b1])
                                    + (delta[a0][a1] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b1][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c1]) * (PQ[b0] * PQ[c0])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][d0] + delta[a0][d0] * delta[a1][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][d0] * delta[b1][c0]) * (PQ[b0] * PQ[c1])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b1][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b1][c1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d1] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b1][c0]) * (PQ[b0] * PQ[d0])
                                    + (delta[a0][a1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][b1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b1][c0]) * (PQ[b0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] + delta[a0][d1] * delta[a1][b0] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1]) * (PQ[b1] * PQ[c0])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][d0] + delta[a0][c0] * delta[a1][b0] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[c0][d1] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][b0] * delta[c0][d0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0]) * (PQ[b1] * PQ[c1])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d1] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d1] + delta[a0][c0] * delta[a1][c1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] + delta[a0][d1] * delta[a1][b0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0]) * (PQ[b1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] + delta[a0][d0] * delta[a1][b0] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0]) * (PQ[b1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1]) * (PQ[c0] * PQ[c1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1]) * (PQ[c0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1]) * (PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][a1] * delta[b0][c0] * delta[b1][c1] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1]) * (PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[5] * (

                                + 0.25 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] + PA_1 * PQ[a0] * PQ[b0] * PQ[c0])
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[c1] + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] + PA_1 * PQ[a0] * PQ[b0] * PQ[c1])
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[d0] + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] + PA_1 * PQ[a0] * PQ[b0] * PQ[d0])
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[a1] * PQ[d1] + PA_0 * PQ[a1] * PQ[b0] * PQ[d1] + PA_1 * PQ[a0] * PQ[b0] * PQ[d1])
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] + PA_1 * PQ[a0] * PQ[b1] * PQ[c0])
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[c1] + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] + PA_1 * PQ[a0] * PQ[b1] * PQ[c1])
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[d0] + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] + PA_1 * PQ[a0] * PQ[b1] * PQ[d0])
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PB_1 * PQ[a0] * PQ[a1] * PQ[d1] + PA_0 * PQ[a1] * PQ[b1] * PQ[d1] + PA_1 * PQ[a0] * PQ[b1] * PQ[d1])
                                    + delta[b0][b1] * delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[a0] * PQ[c0] * PQ[c1])
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[a0] * PQ[c0] * PQ[c1])
                                    + delta[b0][b1] * delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[a0] * PQ[c0] * PQ[d0])
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PA_0 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[a0] * PQ[c0] * PQ[d0])
                                    + delta[b0][b1] * delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[a0] * PQ[c0] * PQ[d1])
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PA_0 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[a0] * PQ[c0] * PQ[d1])
                                    + delta[b0][b1] * delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[a0] * PQ[c1] * PQ[d0])
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[a0] * PQ[c1] * PQ[d0])
                                    + delta[b0][b1] * delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[a0] * PQ[c1] * PQ[d1])
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[a0] * PQ[c1] * PQ[d1])
                                    + delta[b0][b1] * delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[a0] * PQ[d0] * PQ[d1])
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PA_0 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[a0] * PQ[d0] * PQ[d1])
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] + PA_0 * PQ[b0] * PQ[b1] * PQ[c0])
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[c1] + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] + PA_0 * PQ[b0] * PQ[b1] * PQ[c1])
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[d0] + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] + PA_0 * PQ[b0] * PQ[b1] * PQ[d0])
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PB_0 * PQ[a0] * PQ[b1] * PQ[d1] + PB_1 * PQ[a0] * PQ[b0] * PQ[d1] + PA_0 * PQ[b0] * PQ[b1] * PQ[d1])
                                    + delta[a1][b1] * delta[d0][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] + PA_0 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] + PA_0 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + delta[a1][b1] * delta[c1][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] + PA_0 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PB_0 * PQ[a0] * PQ[c0] * PQ[d0] + PA_0 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + delta[a1][b1] * delta[c1][d0] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[d1] + PA_0 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PB_0 * PQ[a0] * PQ[c0] * PQ[d1] + PA_0 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + delta[a1][b1] * delta[c0][d1] * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] + PA_0 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[c1] * PQ[d0] + PA_0 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + delta[a1][b1] * delta[c0][d0] * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[c1] * PQ[d1] + PA_0 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[c1] * PQ[d1] + PA_0 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + delta[a1][b1] * delta[c0][c1] * (PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[d0] * PQ[d1] + PA_0 * PQ[b0] * PQ[d0] * PQ[d1])
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PB_0 * PQ[a0] * PQ[d0] * PQ[d1] + PA_0 * PQ[b0] * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * delta[d0][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] + PA_0 * PQ[b1] * PQ[c0] * PQ[c1])
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] + PA_0 * PQ[b1] * PQ[c0] * PQ[c1])
                                    + delta[a1][b0] * delta[c1][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] + PA_0 * PQ[b1] * PQ[c0] * PQ[d0])
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PB_1 * PQ[a0] * PQ[c0] * PQ[d0] + PA_0 * PQ[b1] * PQ[c0] * PQ[d0])
                                    + delta[a1][b0] * delta[c1][d0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[d1] + PA_0 * PQ[b1] * PQ[c0] * PQ[d1])
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PB_1 * PQ[a0] * PQ[c0] * PQ[d1] + PA_0 * PQ[b1] * PQ[c0] * PQ[d1])
                                    + delta[a1][b0] * delta[c0][d1] * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] + PA_0 * PQ[b1] * PQ[c1] * PQ[d0])
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[c1] * PQ[d0] + PA_0 * PQ[b1] * PQ[c1] * PQ[d0])
                                    + delta[a1][b0] * delta[c0][d0] * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[c1] * PQ[d1] + PA_0 * PQ[b1] * PQ[c1] * PQ[d1])
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[c1] * PQ[d1] + PA_0 * PQ[b1] * PQ[c1] * PQ[d1])
                                    + delta[a1][b0] * delta[c0][c1] * (PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[d0] * PQ[d1] + PA_0 * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PB_1 * PQ[a0] * PQ[d0] * PQ[d1] + PA_0 * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] + PA_1 * PQ[b0] * PQ[b1] * PQ[c0])
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[c1] + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] + PA_1 * PQ[b0] * PQ[b1] * PQ[c1])
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[d0] + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] + PA_1 * PQ[b0] * PQ[b1] * PQ[d0])
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PB_0 * PQ[a1] * PQ[b1] * PQ[d1] + PB_1 * PQ[a1] * PQ[b0] * PQ[d1] + PA_1 * PQ[b0] * PQ[b1] * PQ[d1])
                                    + delta[a0][b1] * delta[d0][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + delta[a0][b1] * delta[c1][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PB_0 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + delta[a0][b1] * delta[c1][d0] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PB_0 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + delta[a0][b1] * delta[c0][d1] * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + delta[a0][b1] * delta[c0][d0] * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + delta[a0][b1] * delta[c0][c1] * (PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PB_0 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b0] * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * delta[d0][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[b1] * PQ[c0] * PQ[c1])
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] + PA_1 * PQ[b1] * PQ[c0] * PQ[c1])
                                    + delta[a0][b0] * delta[c1][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[b1] * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PB_1 * PQ[a1] * PQ[c0] * PQ[d0] + PA_1 * PQ[b1] * PQ[c0] * PQ[d0])
                                    + delta[a0][b0] * delta[c1][d0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[b1] * PQ[c0] * PQ[d1])
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PB_1 * PQ[a1] * PQ[c0] * PQ[d1] + PA_1 * PQ[b1] * PQ[c0] * PQ[d1])
                                    + delta[a0][b0] * delta[c0][d1] * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[b1] * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[c1] * PQ[d0] + PA_1 * PQ[b1] * PQ[c1] * PQ[d0])
                                    + delta[a0][b0] * delta[c0][d0] * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[b1] * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[c1] * PQ[d1] + PA_1 * PQ[b1] * PQ[c1] * PQ[d1])
                                    + delta[a0][b0] * delta[c0][c1] * (PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PB_1 * PQ[a1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][a1] * delta[d0][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] + PB_1 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + delta[a0][a1] * delta[c1][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] + PB_1 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + delta[a0][a1] * delta[c1][d0] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[d1] + PB_1 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + delta[a0][a1] * delta[c0][d1] * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] + PB_1 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + delta[a0][a1] * delta[c0][d0] * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[c1] * PQ[d1] + PB_1 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + delta[a0][a1] * delta[c0][c1] * (PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[d0] * PQ[d1] + PB_1 * PQ[b0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-2.0))
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] + PB_1 * PQ[b0] * PQ[c0] * PQ[c1])
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PB_0 * PQ[b1] * PQ[c0] * PQ[d0] + PB_1 * PQ[b0] * PQ[c0] * PQ[d0])
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PB_0 * PQ[b1] * PQ[c0] * PQ[d1] + PB_1 * PQ[b0] * PQ[c0] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[c1] * PQ[d0] + PB_1 * PQ[b0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[c1] * PQ[d1] + PB_1 * PQ[b0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PB_0 * PQ[b1] * PQ[d0] * PQ[d1] + PB_1 * PQ[b0] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[5] * (

                                + 0.25 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * (-2.0))
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b0] * QC_0 * (-1.0))
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b0] * QC_1 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * (-1.0))
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b0] * QD_1 * (-1.0))
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[b0][b1] * delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PQ[a0] * PQ[a1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[b0][b1] * delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * (-1.0) + PQ[a0] * PQ[b0] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * (-1.0) + PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[a1][b1] * delta[d0][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c1][d1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c1][d0] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PQ[a0] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d1] * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][d0] * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b1] * delta[c0][c1] * (PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a1][b0] * delta[d0][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c1][d1] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c1][d0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PQ[a0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d1] * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][d0] * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a1][b0] * delta[c0][c1] * (PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PQ[a0] * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PQ[a0] * PQ[c1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PQ[a0] * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[c1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * (-1.0) + PQ[a1] * PQ[b0] * PQ[b1] * QC_0 * (-1.0))
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[b1] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * (-1.0) + PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[b1] * QD_1 * (-1.0))
                                    + delta[a0][b1] * delta[d0][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c1][d1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c1][d0] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PQ[a1] * PQ[b0] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d1] * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][d0] * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b1] * delta[c0][c1] * (PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[d1] * QD_0 * (-1.0))
                                    + delta[a0][b0] * delta[d0][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c1][d1] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c1][d0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PQ[a1] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d1] * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][d0] * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][b0] * delta[c0][c1] * (PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PQ[a1] * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PQ[a1] * PQ[c1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PQ[a1] * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[c1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[d0] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[d0][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c1][d1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c1][d0] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PQ[b0] * PQ[b1] * PQ[c0] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * PQ[c0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_0 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d1] * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][d0] * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QC_1 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[c1] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QC_1 * (-1.0))
                                    + delta[a0][a1] * delta[c0][c1] * (PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PQ[b0] * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PQ[b0] * PQ[c1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PQ[b0] * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PQ[b0] * PQ[c1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PQ[b0] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[b0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * (-1.0) + PQ[b1] * PQ[c0] * PQ[d0] * QC_1 * (-1.0) + PQ[b1] * PQ[c1] * PQ[d0] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * QD_1 * (-1.0) + PQ[b1] * PQ[c0] * PQ[d1] * QC_1 * (-1.0) + PQ[b1] * PQ[c1] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * (-1.0) + PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * (-1.0) + PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[b1] * PQ[d0] * PQ[d1] * QC_1 * (-1.0))
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0))
                                )

                            )
                            +
                            F8_t[5] * (

                                + 0.5 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][d1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * (PB_1 * PA_0 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][d1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * (PB_1 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][a1] * (PB_0 * PB_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[5] * (

                                + 0.5 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QC_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QC_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_0 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QC_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[c1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QC_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0)
                                    + delta[c0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_0 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[c0][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1)
                                    + delta[c0][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0)
                                    + delta[b1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b1][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b1][c0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][d1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QC_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b0][d0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b0][c1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b0][c0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][d1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QC_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a1][d0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QC_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a1][c1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][c0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_0 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_0 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][d1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QC_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a0][d0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QC_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a0][c1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][c0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PA_1 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PA_1 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0) + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_0 * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_0 * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PB_1 * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PB_1 * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                )

                            );
                        }
                        else if constexpr(part == 14)
                        {
                            return
                            F8_t[5] * (

                                + S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )

                            )
                            +
                            F8_t[6] * (

                                0.5 * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * (-1.0))
                                    + delta[c1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * (-1.0))
                                    + delta[c1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * (-1.0))
                                    + delta[c0][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[c0][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[c0][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][d1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b1][d0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b1][c1] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b1][c0] * (PB_0 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][d1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[b0][d0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[b0][c1] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][c0] * (PB_1 * PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[b0][b1] * (PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][d1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a1][d0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a1][c1] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][c0] * (PB_0 * PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][d1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * (-1.0))
                                    + delta[a0][d0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * (-1.0))
                                    + delta[a0][c1] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][c0] * (PB_0 * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0))
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PB_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1] * (-1.0) + PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[6] * (

                                + 0.5 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_0)
                                    + delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QC_1)
                                    + delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_0)
                                    + delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QC_1)
                                    + delta[b1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QC_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b1][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b1][c0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][d1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[b0][d0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[b0][c1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[b0][c0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][d1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a1][d0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a1][c1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][c0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][d1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QC_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QC_0)
                                    + delta[a0][d0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QC_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QC_0)
                                    + delta[a0][c1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][c0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 + PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 + PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 + PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0)
                                    + delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0)
                                )

                            )
                            +
                            F8_t[6] * (

                                + S1 * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )

                            )
                            +
                            F8_t[6] * (

                                + S1 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                )

                            )
                            +
                            F8_t[6] * (

                                + S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )

                            )
                            +
                            F8_t[6] * (

                                + 0.25 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1])
                                    + (delta[b1][c1] * delta[d0][d1] + delta[b1][d0] * delta[c1][d1] + delta[b1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0])
                                    + (delta[b1][c0] * delta[d0][d1] + delta[b1][d0] * delta[c0][d1] + delta[b1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1])
                                    + (delta[b1][c0] * delta[c1][d1] + delta[b1][c1] * delta[c0][d1] + delta[b1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[d0])
                                    + (delta[b1][c0] * delta[c1][d0] + delta[b1][c1] * delta[c0][d0] + delta[b1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[d1])
                                    + (delta[b0][c1] * delta[d0][d1] + delta[b0][d0] * delta[c1][d1] + delta[b0][d1] * delta[c1][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0])
                                    + (delta[b0][c0] * delta[d0][d1] + delta[b0][d0] * delta[c0][d1] + delta[b0][d1] * delta[c0][d0]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1])
                                    + (delta[b0][c0] * delta[c1][d1] + delta[b0][c1] * delta[c0][d1] + delta[b0][d1] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[d0])
                                    + (delta[b0][c0] * delta[c1][d0] + delta[b0][c1] * delta[c0][d0] + delta[b0][d0] * delta[c0][c1]) * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[d1])
                                    + (delta[b0][b1] * delta[d0][d1] + delta[b0][d0] * delta[b1][d1] + delta[b0][d1] * delta[b1][d0]) * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1])
                                    + (delta[b0][b1] * delta[c1][d1] + delta[b0][c1] * delta[b1][d1] + delta[b0][d1] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d0])
                                    + (delta[b0][b1] * delta[c1][d0] + delta[b0][c1] * delta[b1][d0] + delta[b0][d0] * delta[b1][c1]) * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[d1])
                                    + (delta[b0][b1] * delta[c0][d1] + delta[b0][c0] * delta[b1][d1] + delta[b0][d1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d0])
                                    + (delta[b0][b1] * delta[c0][d0] + delta[b0][c0] * delta[b1][d0] + delta[b0][d0] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[c1] * PQ[d1])
                                    + (delta[b0][b1] * delta[c0][c1] + delta[b0][c0] * delta[b1][c1] + delta[b0][c1] * delta[b1][c0]) * (PQ[a0] * PQ[a1] * PQ[d0] * PQ[d1])
                                    + (delta[a1][c1] * delta[d0][d1] + delta[a1][d0] * delta[c1][d1] + delta[a1][d1] * delta[c1][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0])
                                    + (delta[a1][c0] * delta[d0][d1] + delta[a1][d0] * delta[c0][d1] + delta[a1][d1] * delta[c0][d0]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1])
                                    + (delta[a1][c0] * delta[c1][d1] + delta[a1][c1] * delta[c0][d1] + delta[a1][d1] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[d0])
                                    + (delta[a1][c0] * delta[c1][d0] + delta[a1][c1] * delta[c0][d0] + delta[a1][d0] * delta[c0][c1]) * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[d1])
                                    + (delta[a1][b1] * delta[d0][d1] + delta[a1][d0] * delta[b1][d1] + delta[a1][d1] * delta[b1][d0]) * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1])
                                    + (delta[a1][b1] * delta[c1][d1] + delta[a1][c1] * delta[b1][d1] + delta[a1][d1] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d0])
                                    + (delta[a1][b1] * delta[c1][d0] + delta[a1][c1] * delta[b1][d0] + delta[a1][d0] * delta[b1][c1]) * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[d1])
                                    + (delta[a1][b1] * delta[c0][d1] + delta[a1][c0] * delta[b1][d1] + delta[a1][d1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d0])
                                    + (delta[a1][b1] * delta[c0][d0] + delta[a1][c0] * delta[b1][d0] + delta[a1][d0] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[c1] * PQ[d1])
                                    + (delta[a1][b1] * delta[c0][c1] + delta[a1][c0] * delta[b1][c1] + delta[a1][c1] * delta[b1][c0]) * (PQ[a0] * PQ[b0] * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[d0][d1] + delta[a1][d0] * delta[b0][d1] + delta[a1][d1] * delta[b0][d0]) * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1])
                                    + (delta[a1][b0] * delta[c1][d1] + delta[a1][c1] * delta[b0][d1] + delta[a1][d1] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d0])
                                    + (delta[a1][b0] * delta[c1][d0] + delta[a1][c1] * delta[b0][d0] + delta[a1][d0] * delta[b0][c1]) * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[d1])
                                    + (delta[a1][b0] * delta[c0][d1] + delta[a1][c0] * delta[b0][d1] + delta[a1][d1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d0])
                                    + (delta[a1][b0] * delta[c0][d0] + delta[a1][c0] * delta[b0][d0] + delta[a1][d0] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[c1] * PQ[d1])
                                    + (delta[a1][b0] * delta[c0][c1] + delta[a1][c0] * delta[b0][c1] + delta[a1][c1] * delta[b0][c0]) * (PQ[a0] * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][d1] + delta[a1][b1] * delta[b0][d1] + delta[a1][d1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a1][b0] * delta[b1][d0] + delta[a1][b1] * delta[b0][d0] + delta[a1][d0] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][c1] + delta[a1][b1] * delta[b0][c1] + delta[a1][c1] * delta[b0][b1]) * (PQ[a0] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a1][b0] * delta[b1][c0] + delta[a1][b1] * delta[b0][c0] + delta[a1][c0] * delta[b0][b1]) * (PQ[a0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][c1] * delta[d0][d1] + delta[a0][d0] * delta[c1][d1] + delta[a0][d1] * delta[c1][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0])
                                    + (delta[a0][c0] * delta[d0][d1] + delta[a0][d0] * delta[c0][d1] + delta[a0][d1] * delta[c0][d0]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1])
                                    + (delta[a0][c0] * delta[c1][d1] + delta[a0][c1] * delta[c0][d1] + delta[a0][d1] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0])
                                    + (delta[a0][c0] * delta[c1][d0] + delta[a0][c1] * delta[c0][d0] + delta[a0][d0] * delta[c0][c1]) * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1])
                                    + (delta[a0][b1] * delta[d0][d1] + delta[a0][d0] * delta[b1][d1] + delta[a0][d1] * delta[b1][d0]) * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1])
                                    + (delta[a0][b1] * delta[c1][d1] + delta[a0][c1] * delta[b1][d1] + delta[a0][d1] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0])
                                    + (delta[a0][b1] * delta[c1][d0] + delta[a0][c1] * delta[b1][d0] + delta[a0][d0] * delta[b1][c1]) * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[d1])
                                    + (delta[a0][b1] * delta[c0][d1] + delta[a0][c0] * delta[b1][d1] + delta[a0][d1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][b1] * delta[c0][d0] + delta[a0][c0] * delta[b1][d0] + delta[a0][d0] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][b1] * delta[c0][c1] + delta[a0][c0] * delta[b1][c1] + delta[a0][c1] * delta[b1][c0]) * (PQ[a1] * PQ[b0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[d0][d1] + delta[a0][d0] * delta[b0][d1] + delta[a0][d1] * delta[b0][d0]) * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1])
                                    + (delta[a0][b0] * delta[c1][d1] + delta[a0][c1] * delta[b0][d1] + delta[a0][d1] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0])
                                    + (delta[a0][b0] * delta[c1][d0] + delta[a0][c1] * delta[b0][d0] + delta[a0][d0] * delta[b0][c1]) * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[d1])
                                    + (delta[a0][b0] * delta[c0][d1] + delta[a0][c0] * delta[b0][d1] + delta[a0][d1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0])
                                    + (delta[a0][b0] * delta[c0][d0] + delta[a0][c0] * delta[b0][d0] + delta[a0][d0] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[c1] * PQ[d1])
                                    + (delta[a0][b0] * delta[c0][c1] + delta[a0][c0] * delta[b0][c1] + delta[a0][c1] * delta[b0][c0]) * (PQ[a1] * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][d1] + delta[a0][b1] * delta[b0][d1] + delta[a0][d1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][b0] * delta[b1][d0] + delta[a0][b1] * delta[b0][d0] + delta[a0][d0] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][c1] + delta[a0][b1] * delta[b0][c1] + delta[a0][c1] * delta[b0][b1]) * (PQ[a1] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][b0] * delta[b1][c0] + delta[a0][b1] * delta[b0][c0] + delta[a0][c0] * delta[b0][b1]) * (PQ[a1] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[d0][d1] + delta[a0][d0] * delta[a1][d1] + delta[a0][d1] * delta[a1][d0]) * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1])
                                    + (delta[a0][a1] * delta[c1][d1] + delta[a0][c1] * delta[a1][d1] + delta[a0][d1] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0])
                                    + (delta[a0][a1] * delta[c1][d0] + delta[a0][c1] * delta[a1][d0] + delta[a0][d0] * delta[a1][c1]) * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1])
                                    + (delta[a0][a1] * delta[c0][d1] + delta[a0][c0] * delta[a1][d1] + delta[a0][d1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[c0][d0] + delta[a0][c0] * delta[a1][d0] + delta[a0][d0] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[c0][c1] + delta[a0][c0] * delta[a1][c1] + delta[a0][c1] * delta[a1][c0]) * (PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][d1] + delta[a0][b1] * delta[a1][d1] + delta[a0][d1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b1][d0] + delta[a0][b1] * delta[a1][d0] + delta[a0][d0] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][c1] + delta[a0][b1] * delta[a1][c1] + delta[a0][c1] * delta[a1][b1]) * (PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b1][c0] + delta[a0][b1] * delta[a1][c0] + delta[a0][c0] * delta[a1][b1]) * (PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][d1] + delta[a0][b0] * delta[a1][d1] + delta[a0][d1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + (delta[a0][a1] * delta[b0][d0] + delta[a0][b0] * delta[a1][d0] + delta[a0][d0] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c1] + delta[a0][b0] * delta[a1][c1] + delta[a0][c1] * delta[a1][b0]) * (PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][c0] + delta[a0][b0] * delta[a1][c0] + delta[a0][c0] * delta[a1][b0]) * (PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[7] * (

                                (-0.5) * S1 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1])
                                    + delta[c1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0])
                                    + delta[c1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1])
                                    + delta[c0][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0])
                                    + delta[c0][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1])
                                    + delta[c0][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1])
                                    + delta[b1][d1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + delta[b1][d0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + delta[b1][c1] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + delta[b1][c0] * (PQ[a0] * PQ[a1] * PQ[b0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[b0][d1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + delta[b0][d0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + delta[b0][c1] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + delta[b0][c0] * (PQ[a0] * PQ[a1] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[b0][b1] * (PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][d1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + delta[a1][d0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + delta[a1][c1] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + delta[a1][c0] * (PQ[a0] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][d1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0])
                                    + delta[a0][d0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1])
                                    + delta[a0][c1] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1])
                                    + delta[a0][c0] * (PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1])
                                )

                            )
                            +
                            F8_t[7] * (

                                + S1 * S1 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )

                            )
                            +
                            F8_t[7] * (

                                + (-1.0) * S1 * S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0
                                )

                            )
                            +
                            F8_t[8] * (

                                S1 * S1 * S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )

                            );
                        }
                        else
                        {
                            static_assert(false);
                        }

                    }();

                    ERIs[threadIdx.y][threadIdx.x] = eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];
                }
                else
                {
                    ERIs[threadIdx.y][threadIdx.x] = 0.0;

                    // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
                    if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = m + 1;
                }
            }
            else
            {
                ERIs[threadIdx.y][threadIdx.x] = 0.0;

                // indicator for early exit for thread block (ERIs[0][0] has the largest upper bound)
                if ((threadIdx.y == 0) && (threadIdx.x == 0)) skip_thread_block = m + 1;
            }

            __syncthreads();

            // early exit for thread block
            if (skip_thread_block == m + 1) break;

            if ((threadIdx.y == 0) && (threadIdx.x == 0))
            {
                for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
                {
                    for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
                    {
                        K_ik += ERIs[y][x];
                    }
                }
            }

            __syncthreads();
        }
    }

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_dd))
    {
        mat_K[ik] += K_ik;
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDSDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ds,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_ds,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_ds,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dd,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_dd,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_dd,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDSDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ds,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_ds,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_ds,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDPDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ps_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_dp,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_dp,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_dp,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD0(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD1(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD2(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD3(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD4(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD5(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDD6(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_pp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_pp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_pp,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_pp,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPSDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ps,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_ps,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_ps,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_ps,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_ps,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDDS(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pd,
                        const uint32_t* pair_inds_k_for_K_pd,
                        const uint32_t  pair_inds_count_for_K_pd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    ds_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_ds,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_ds,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_ds,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_ds,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_ds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSPDD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sp,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_sp,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_sp,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_sp,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_sp,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockSDDP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_sd,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_sd,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_sd,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_sd,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_sd,
                        const double*   pair_data_K_dp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPPPD(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    pd_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pp,
                        const double*   Q_K_pd,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_counts_K_pp,
                        const uint32_t* pair_counts_K_pd,
                        const double*   pair_data_K_pp,
                        const double*   pair_data_K_pd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockPDPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_pd,
                        const double*   Q_K_pp,
                        const uint32_t* D_inds_K_pd,
                        const uint32_t* D_inds_K_pp,
                        const uint32_t* pair_displs_K_pd,
                        const uint32_t* pair_displs_K_pp,
                        const uint32_t* pair_counts_K_pd,
                        const uint32_t* pair_counts_K_pp,
                        const double*   pair_data_K_pd,
                        const double*   pair_data_K_pp,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);


}  // namespace gpu

#endif

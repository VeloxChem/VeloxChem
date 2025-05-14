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

#include <optional>

#include "GpuConstants.hpp"

namespace gpu {  // gpu namespace

template<uint32_t size>
struct GpuIntermediateBlockData
{
    std::array<double, size> boysFuncData_;
    std::array<double, 3> PQ_;
    double Lambda_;
    double S_ij_00_;
    double S_kl_00_;
    double S1_;
    double S2_;
    double inv_S1_;
    double inv_S2_;
    double inv_S4_;
    double PA_0_;
    double PA_1_;
    double PB_0_;
    double PB_1_;
    double QC_0_;
    double QC_1_;
    double QD_0_;
    double QD_1_;
    uint32_t a0_;
    uint32_t a1_;
    uint32_t b0_;
    uint32_t b1_;
    uint32_t c0_;
    uint32_t c1_;
    uint32_t d0_;
    uint32_t d1_;
    uint32_t j_cgto_;
    uint32_t l_cgto_;
    bool valid_;

    __device__ __forceinline__
    GpuIntermediateBlockData(std::array<double, size> boysFuncData, std::array<double,3> PQ, double Lambda,
                    double S_ij_00, double S_kl_00, double S1, double S2, double inv_S1, double inv_S2,
                    double inv_S4, double PA_0, double PA_1, double PB_0, double PB_1, double QC_0, double QC_1,
                    double QD_0, double QD_1, uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1, uint32_t c0,
                    uint32_t c1, uint32_t d0, uint32_t d1, uint32_t j_cgto, uint32_t l_cgto)
            : boysFuncData_(boysFuncData), PQ_(PQ), Lambda_(Lambda), S_ij_00_(S_ij_00), S_kl_00_(S_kl_00),
              S1_(S1), S2_(S2), inv_S1_(inv_S1), inv_S2_(inv_S2), inv_S4_(inv_S4), PA_0_(PA_0), PA_1_(PA_1), PB_0_(PB_0), PB_1_(PB_1),
              QC_0_(QC_0), QC_1_(QC_1), QD_0_(QD_0), QD_1_(QD_1), a0_(a0), a1_(a1), b0_(b0), b1_(b1), c0_(c0),
              c1_(c1), d0_(d0), d1_(d1), j_cgto_(j_cgto), l_cgto_(l_cgto), valid_(true)
        {}

    __device__ __forceinline__
    GpuIntermediateBlockData()
        {}

    __device__ __forceinline__
    GpuIntermediateBlockData(bool valid) : valid_(valid)
        {}

};

template<uint32_t size>
using DeviceStore = GpuIntermediateBlockData<size>;

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

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDDPrepare(
                        DeviceStore<9>* data,
                        const uint32_t* pair_inds_i_for_K_dd,
                        const uint32_t* pair_inds_k_for_K_dd,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    dd_max_D,
                        const double*   Q_K_dd,
                        const uint32_t* D_inds_K_dd,
                        const uint32_t* pair_displs_K_dd,
                        const uint32_t* pair_counts_K_dd,
                        const double*   pair_data_K_dd,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold);

static inline constexpr uint32_t d_cart_inds[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
static inline constexpr double   delta[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

template<uint32_t part>
__device__ double __attribute__((always_inline)) 
computeExchangeFockDDDDInternal(const double* F8_t,
                         const double* PQ,
                         const double Lambda,
                         const double S_ij_00,
                         const double S_kl_00,
                         const double S1,
                         const double S2,
                         const double inv_S1,
                         const double inv_S2,
                         const double inv_S4,
                         const double PA_0,
                         const double PA_1,
                         const double PB_0,
                         const double PB_1,
                         const double QC_0,
                         const double QC_1,
                         const double QD_0,
                         const double QD_1,
                         const uint32_t a0,
                         const uint32_t a1,
                         const uint32_t b0,
                         const uint32_t b1,
                         const uint32_t c0,
                         const uint32_t c1,
                         const uint32_t d0,
                         const uint32_t d1)
{
						if constexpr (part == 0)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.25 * inv_S1 * inv_S1 * (
                                    (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
                                )

                            );
						}
						else if constexpr(part == 1)
						{
                            return Lambda * S_ij_00 * S_kl_00 *F8_t[0] * (

                                + 0.25 * inv_S2 * inv_S2 * (
                                    (delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
                                )

                            );
						}
						else if constexpr(part == 2)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1
                                )

                            );
						}
						else if constexpr(part == 3)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.0625 * inv_S1 * inv_S1 * inv_S2 * inv_S2 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            );
						}
						else if constexpr(part == 4)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

                                (-0.125) * inv_S1 * inv_S1 * inv_S2 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            );
						}
						else if constexpr(part == 5)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

                                + (-0.125) * inv_S1 * inv_S2 * inv_S2 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )

                            );
						}
						else if constexpr(part == 6)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                0.125 * inv_S1 * inv_S1 * inv_S2 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )

                            );
						}
						else if constexpr(part == 7)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.125 * inv_S1 * inv_S2 * inv_S2 * (
                                    (delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
                                    + (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
                                    + (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
                                    + (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
                                    + (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
                                    + (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
                                )

                            );
						}
						else if constexpr(part == 8)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

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
                            );
						}
						else if constexpr(part == 9)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.5 * inv_S1 * (
                                    delta[b0][b1] * (PA_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PA_0 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PA_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PB_1 * QD_0 * QD_1 * QC_0 * QC_1)
                                )
                            );
						}
						else if constexpr(part == 10)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.5 * inv_S2 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1)
                                )
                            );
						}
						else if constexpr(part == 11)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

                                + (-0.25) * inv_S1 * inv_S1 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )
                            );
						}
						else if constexpr(part == 12)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[0] * (

                                + 0.5 * inv_S2 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QC_0 * QC_1)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1)
                                )
                            );
						}
						else if constexpr(part == 13)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

                                + (-0.25) * inv_S1 * inv_S1 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
                                    + (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
                                )
                            );
						}
						else if constexpr(part == 14)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

								+ (-0.5) * S1 * inv_S2 * inv_S2 * inv_S4 * (
									(delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
								)
							);
						}
						else if constexpr(part == 15)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

								+ (-0.5) * S2 * inv_S1 * inv_S1 * inv_S4 * (
									(delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
								)
							);
						}
						else if constexpr(part == 16)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

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
							);
						}
						else if constexpr(part == 17)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

									+ (-0.25) * inv_S2 * inv_S2 * inv_S4 * (
									(delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
									+ (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
									+ (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
									+ (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
									+ (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
									+ (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
								)
							);
						}
						else if constexpr(part == 18)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

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
							);
						}
						else if constexpr(part == 19)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (

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
						else if constexpr(part == 20)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (
									+ (-0.5) * S1 * inv_S2 * inv_S4 * (
									delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QC_0 * QC_1)
									+ delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_0)
									+ delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_0)
									+ delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QC_1)
									+ delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * QD_1 * QC_1)
									+ delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * QD_0 * QD_1)
								)
							);
						}
						else if constexpr(part == 21)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (
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
						else if constexpr(part == 22)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
								+ 0.25 * S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * (
									(delta[c0][c1] * delta[d0][d1] + delta[c0][d0] * delta[c1][d1] + delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1 * PA_0 * PA_1)
								)
							);
						}
						else if constexpr(part == 23)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
								+ 0.25 * S2 * S2 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
									(delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (QD_0 * QD_1 * QC_0 * QC_1)
								)
							);
						}
						else if constexpr(part == 24)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (
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
							);
						}
						else if constexpr(part == 25)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (
								+ S1 * inv_S4 * (
									+ PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
									+ PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
									+ PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
									+ PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
								)
							);
						}
						else if constexpr(part == 26)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[1] * (
								+ S2 * inv_S4 * (                    
									+ PB_0 * PB_1 * PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1
									+ PB_0 * PB_1 * PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1
									+ PB_0 * PA_0 * PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
									+ PB_1 * PA_0 * PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
								)
							);
						}
						else if constexpr(part == 27)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
								0.125 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * (
									(delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[b0][b1] * delta[c0][d1] * delta[c1][d0]) * (PA_0 * PA_1)
									+ (delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a1][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_0)
									+ (delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_0)
									+ (delta[a0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PA_1)
									+ (delta[a0][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[c0][d1] * delta[c1][d0]) * (PB_1 * PA_1)
									+ (delta[a0][a1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[c0][d1] * delta[c1][d0]) * (PB_0 * PB_1)
								)
							);
						}
						else if constexpr(part == 28)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
								+ 0.125 * S2 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
									(delta[a0][a1] * delta[b0][b1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[d0][d1]) * (QC_0 * QC_1)
									+ (delta[a0][a1] * delta[b0][b1] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c1][d1]) * (QD_0 * QC_0)
									+ (delta[a0][a1] * delta[b0][b1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c1][d0]) * (QD_1 * QC_0)
									+ (delta[a0][a1] * delta[b0][b1] * delta[c0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1]) * (QD_0 * QC_1)
									+ (delta[a0][a1] * delta[b0][b1] * delta[c0][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0]) * (QD_1 * QC_1)
									+ (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1]) * (QD_0 * QD_1)
								)
							);
						}
						else if constexpr(part == 29)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
						else if constexpr(part == 30)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
							);
						}
						else if constexpr(part == 31)
						{
							return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
							);
						}
						else if constexpr(part == 32)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
                            );
						}
						else if constexpr(part == 33)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
                            );
						}
						else if constexpr(part == 34)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
                                + 0.5 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * (
                                    delta[d0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QC_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QC_0)
                                    + delta[c1][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_0)
                                    + delta[c1][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_0)
                                    + delta[c0][d1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_0 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QC_1)
                                    + delta[c0][d0] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QC_1)
                                    + delta[c0][c1] * (PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * QD_1 + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d1] * QD_0)
                                )
                            );
						}
						else if constexpr(part == 35)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
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
                            );
						}
						else if constexpr(part == 36)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (
                                + 0.5 * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * (PA_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PB_0 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PB_1 * PQ[a0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PB_0 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PB_1 * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PA_1 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PB_0 * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PB_1 * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1 * (-1.0) + PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                )
                            );
						}
						else if constexpr(part == 37)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

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
                            );
						}
						else if constexpr(part == 38)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

                                + S1 * S1 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )
                            );
						}
						else if constexpr(part == 39)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

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
						else if constexpr(part == 40)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

                                + S2 * S2 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )
                            );
						}
						else if constexpr(part == 41)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

                                + 0.0625 * inv_S1 * inv_S1 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0])
                                )
                            );
						}
						else if constexpr(part == 42)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[2] * (

                                + 0.0625 * inv_S1 * inv_S2 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * 4.0
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1])
                                )
                            );
						}
						else if constexpr(part == 43)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

                                0.0625 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (-2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (-1.0)
                                )
                            );
						}
						else if constexpr(part == 44)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

                                + 0.0625 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0]) * (-2.0)
                                    + (delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1]) * (-1.0)
                                )
                            );
						}
						else if constexpr(part == 45)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 46)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 47)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 48)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
						else if constexpr(part == 49)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 50)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
						else if constexpr(part == 51)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 52)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
						else if constexpr(part == 53)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 54)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 55)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

                                + (-0.5) * S2 * S2 * S2 * inv_S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    delta[b0][b1] * (PQ[a0] * PQ[a1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b1] * (PQ[a0] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a1][b0] * (PQ[a0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b1] * (PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][b0] * (PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                    + delta[a0][a1] * (PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1)
                                )
                            );
						}
						else if constexpr(part == 56)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
						else if constexpr(part == 57)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

                                + S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * QD_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d1] * QD_0 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[d0] * PQ[d1] * QC_1 * (-1.0)
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c1] * PQ[d0] * PQ[d1] * QC_0 * (-1.0)
                                )
                            );
						}
						else if constexpr(part == 58)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 59)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

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
                            );
						}
						else if constexpr(part == 60)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[3] * (

                                + S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )
                            );
						}
						else if constexpr(part == 61)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 62)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
						else if constexpr(part == 63)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 64)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 65)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
						else if constexpr(part == 66)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 67)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
						else if constexpr(part == 68)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 69)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 70)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

                                + S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PA_1 * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )
                            );
						}
						else if constexpr(part == 71)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

                                + S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * QD_0 * QD_1 * QC_0 * QC_1
                                )
                            );
						}
						else if constexpr(part == 72)
						{	
							return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

                                + 0.0625 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    (delta[a0][a1] * delta[b0][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][a1] * delta[b0][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][a1] * delta[b0][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][a1] * delta[b0][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][a1] * delta[b0][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][a1] * delta[b0][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][a1] * delta[b0][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b0] * delta[a1][b1] * delta[c0][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][b1] * delta[c0][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c0] * delta[b1][c1] * delta[d0][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d0] * delta[c1][d1] + delta[a0][b0] * delta[a1][c0] * delta[b1][d1] * delta[c1][d0] + delta[a0][b0] * delta[a1][c1] * delta[b1][c0] * delta[d0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d0] * delta[c0][d1] + delta[a0][b0] * delta[a1][c1] * delta[b1][d1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d0] * delta[b1][c0] * delta[c1][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][c1] * delta[c0][d1] + delta[a0][b0] * delta[a1][d0] * delta[b1][d1] * delta[c0][c1] + delta[a0][b0] * delta[a1][d1] * delta[b1][c0] * delta[c1][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][c1] * delta[c0][d0] + delta[a0][b0] * delta[a1][d1] * delta[b1][d0] * delta[c0][c1] + delta[a0][b1] * delta[a1][b0] * delta[c0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][b0] * delta[c0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c0] * delta[b0][c1] * delta[d0][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d0] * delta[c1][d1] + delta[a0][b1] * delta[a1][c0] * delta[b0][d1] * delta[c1][d0] + delta[a0][b1] * delta[a1][c1] * delta[b0][c0] * delta[d0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d0] * delta[c0][d1] + delta[a0][b1] * delta[a1][c1] * delta[b0][d1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d0] * delta[b0][c0] * delta[c1][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][c1] * delta[c0][d1] + delta[a0][b1] * delta[a1][d0] * delta[b0][d1] * delta[c0][c1] + delta[a0][b1] * delta[a1][d1] * delta[b0][c0] * delta[c1][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][c1] * delta[c0][d0] + delta[a0][b1] * delta[a1][d1] * delta[b0][d0] * delta[c0][c1] + delta[a0][c0] * delta[a1][b0] * delta[b1][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b0] * delta[b1][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][b1] * delta[b0][c1] * delta[d0][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d0] * delta[c1][d1] + delta[a0][c0] * delta[a1][b1] * delta[b0][d1] * delta[c1][d0] + delta[a0][c0] * delta[a1][c1] * delta[b0][b1] * delta[d0][d1] + delta[a0][c0] * delta[a1][c1] * delta[b0][d0] * delta[b1][d1] + delta[a0][c0] * delta[a1][c1] * delta[b0][d1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d0] * delta[b0][b1] * delta[c1][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][c1] * delta[b1][d1] + delta[a0][c0] * delta[a1][d0] * delta[b0][d1] * delta[b1][c1] + delta[a0][c0] * delta[a1][d1] * delta[b0][b1] * delta[c1][d0] + delta[a0][c0] * delta[a1][d1] * delta[b0][c1] * delta[b1][d0] + delta[a0][c0] * delta[a1][d1] * delta[b0][d0] * delta[b1][c1] + delta[a0][c1] * delta[a1][b0] * delta[b1][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b0] * delta[b1][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][b1] * delta[b0][c0] * delta[d0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d0] * delta[c0][d1] + delta[a0][c1] * delta[a1][b1] * delta[b0][d1] * delta[c0][d0] + delta[a0][c1] * delta[a1][c0] * delta[b0][b1] * delta[d0][d1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d0] * delta[b1][d1] + delta[a0][c1] * delta[a1][c0] * delta[b0][d1] * delta[b1][d0] + delta[a0][c1] * delta[a1][d0] * delta[b0][b1] * delta[c0][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][c0] * delta[b1][d1] + delta[a0][c1] * delta[a1][d0] * delta[b0][d1] * delta[b1][c0] + delta[a0][c1] * delta[a1][d1] * delta[b0][b1] * delta[c0][d0] + delta[a0][c1] * delta[a1][d1] * delta[b0][c0] * delta[b1][d0] + delta[a0][c1] * delta[a1][d1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d0] * delta[a1][b0] * delta[b1][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b0] * delta[b1][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c0] * delta[c1][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][c1] * delta[c0][d1] + delta[a0][d0] * delta[a1][b1] * delta[b0][d1] * delta[c0][c1] + delta[a0][d0] * delta[a1][c0] * delta[b0][b1] * delta[c1][d1] + delta[a0][d0] * delta[a1][c0] * delta[b0][c1] * delta[b1][d1] + delta[a0][d0] * delta[a1][c0] * delta[b0][d1] * delta[b1][c1] + delta[a0][d0] * delta[a1][c1] * delta[b0][b1] * delta[c0][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][c0] * delta[b1][d1] + delta[a0][d0] * delta[a1][c1] * delta[b0][d1] * delta[b1][c0] + delta[a0][d0] * delta[a1][d1] * delta[b0][b1] * delta[c0][c1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c0] * delta[b1][c1] + delta[a0][d0] * delta[a1][d1] * delta[b0][c1] * delta[b1][c0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b0] * delta[b1][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][b1] * delta[b0][c0] * delta[c1][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][c1] * delta[c0][d0] + delta[a0][d1] * delta[a1][b1] * delta[b0][d0] * delta[c0][c1] + delta[a0][d1] * delta[a1][c0] * delta[b0][b1] * delta[c1][d0] + delta[a0][d1] * delta[a1][c0] * delta[b0][c1] * delta[b1][d0] + delta[a0][d1] * delta[a1][c0] * delta[b0][d0] * delta[b1][c1] + delta[a0][d1] * delta[a1][c1] * delta[b0][b1] * delta[c0][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][c0] * delta[b1][d0] + delta[a0][d1] * delta[a1][c1] * delta[b0][d0] * delta[b1][c0] + delta[a0][d1] * delta[a1][d0] * delta[b0][b1] * delta[c0][c1] + delta[a0][d1] * delta[a1][d0] * delta[b0][c0] * delta[b1][c1] + delta[a0][d1] * delta[a1][d0] * delta[b0][c1] * delta[b1][c0])
                                )
                            );
						}
						else if constexpr(part == 73)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 74)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 75)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[4] * (

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
                            );
						}
						else if constexpr(part == 76)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 77)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 78)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 79)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 80)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
						else if constexpr(part == 81)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 82)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

                                + S1 * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PA_0 * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PB_1 * PA_1 * PQ[a0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_0 * PA_1 * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_0 * PA_1 * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )
                            );
						}
						else if constexpr(part == 83)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
                            );
						}
						else if constexpr(part == 84)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

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
						else if constexpr(part == 85)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[5] * (

                                + S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * QD_0 * QD_1 * QC_1 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * QD_0 * QD_1 * QC_0 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * QD_1 * QC_0 * QC_1 * (-1.0)
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d1] * QD_0 * QC_0 * QC_1 * (-1.0)
                                )
							);
                            
						}
						else if constexpr(part == 86)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

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
                            );
						}
						else if constexpr(part == 87)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

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
                            );
						}
						else if constexpr(part == 88)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

                                + S1 * S1 * S1 * S1 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PB_1 * PQ[a0] * PQ[a1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_0 * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_0 * PA_1 * PQ[a0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_0 * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PA_1 * PQ[a0] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_0 * PA_1 * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )
                            );
						}
						else if constexpr(part == 89)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

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
                            );
						}
						else if constexpr(part == 90)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

                                + S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * QD_0 * QD_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * QD_1 * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d1] * QD_0 * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * QD_1 * QC_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d1] * QD_0 * QC_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[d0] * PQ[d1] * QC_0 * QC_1
                                )
                            );
						}
						else if constexpr(part == 91)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[6] * (

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
                            );
						}
						else if constexpr(part == 92)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[7] * (

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
                            );
						}
						else if constexpr(part == 93)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[7] * (

                                + S1 * S1 * S1 * S1 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    
                                    + PB_0 * PQ[a0] * PQ[a1] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PB_1 * PQ[a0] * PQ[a1] * PQ[b0] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_0 * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                    + PA_1 * PQ[a0] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )
							);
                            
						}
						else if constexpr(part == 94)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[7] * (

                                + (-1.0) * S1 * S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * QD_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d1] * QD_0
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[d0] * PQ[d1] * QC_1
                                    + PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c1] * PQ[d0] * PQ[d1] * QC_0
                                )
                            );
						}
						else if constexpr(part == 95)
						{
                            return Lambda * S_ij_00 * S_kl_00 * F8_t[8] * (

                                S1 * S1 * S1 * S1 * S2 * S2 * S2 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * (
                                    PQ[a0] * PQ[a1] * PQ[b0] * PQ[b1] * PQ[c0] * PQ[c1] * PQ[d0] * PQ[d1]
                                )
                            );
						}
						else
						{
							static_assert(false);
						}
}

template<uint32_t part>
__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeFockDDDD(
                        double*         mat_K,
                        const DeviceStore<9>*         data,
                        const double*   mat_D_full_AO,
                        const uint32_t  pair_inds_count_for_K_dd,
                        const uint32_t  naos)
{

    __shared__ double   ERIs[TILE_SIZE_K];

    double K_ik = 0.0;
    const uint32_t ik = blockIdx.x;
    if (ik < pair_inds_count_for_K_dd )
    {
        const auto* storedData = data + (TILE_SIZE_K * ik);

        const uint32_t index = threadIdx.x;
        if (index < TILE_SIZE_K && storedData[index].valid_)
        {
            const auto& entry = storedData[index];
            ERIs[index] = computeExchangeFockDDDDInternal<part>(entry.boysFuncData_.data(), entry.PQ_.data(), entry.Lambda_, entry.S_ij_00_,
                            entry.S_kl_00_, entry.S1_, entry.S2_, entry.inv_S1_, entry.inv_S2_, entry.inv_S4_, entry.PA_0_,
                            entry.PA_1_, entry.PB_0_, entry.PB_1_, entry.QC_0_, entry.QC_1_, entry.QD_0_, entry.QD_1_,
                            entry.a0_, entry.a1_, entry.b0_, entry.b1_, entry.c0_, entry.c1_, entry.d0_, entry.d1_) *
                            mat_D_full_AO[entry.j_cgto_ * naos + entry.l_cgto_];
        }
        else
        {
            ERIs[index] = 0.0;
        }

        __syncthreads();

        if (threadIdx.x == 0)
        {
            for (uint32_t x = 0; x < TILE_SIZE_K; x++)
            {
               K_ik += ERIs[x];
            }
            mat_K[ik] += K_ik;
        }
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

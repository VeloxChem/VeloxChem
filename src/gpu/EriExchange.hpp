//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#ifndef EriExchange_hpp
#define EriExchange_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void computeExchangeFockSSSS(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_ss,
                                        const uint32_t* pair_inds_k_for_K_ss,
                                        const uint32_t  pair_inds_count_for_K_ss,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full,
                                        const double*   mat_Q_for_K,
                                        const uint32_t* density_inds_for_K,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSSSP(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_ss,
                                        const uint32_t* pair_inds_k_for_K_ss,
                                        const uint32_t  pair_inds_count_for_K_ss,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_ss,
                                        const double*   mat_Q_for_K_sp,
                                        const uint32_t* density_inds_for_K_ss,
                                        const uint32_t* density_inds_for_K_sp,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSPSS(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_ss,
                                        const uint32_t* pair_inds_k_for_K_ss,
                                        const uint32_t  pair_inds_count_for_K_ss,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_ss,
                                        const double*   mat_Q_for_K_sp,
                                        const uint32_t* density_inds_for_K_ss,
                                        const uint32_t* density_inds_for_K_sp,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSPSP(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_ss,
                                        const uint32_t* pair_inds_k_for_K_ss,
                                        const uint32_t  pair_inds_count_for_K_ss,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_sp,
                                        const uint32_t* density_inds_for_K_sp,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSSPS(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_sp,
                                        const uint32_t* pair_inds_k_for_K_sp,
                                        const uint32_t  pair_inds_count_for_K_sp,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_ss,
                                        const double*   mat_Q_for_K_ps,
                                        const uint32_t* density_inds_for_K_ss,
                                        const uint32_t* density_inds_for_K_ps,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSSPP(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_sp,
                                        const uint32_t* pair_inds_k_for_K_sp,
                                        const uint32_t  pair_inds_count_for_K_sp,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_ss,
                                        const double*   mat_Q_for_K_pp,
                                        const uint32_t* density_inds_for_K_ss,
                                        const uint32_t* density_inds_for_K_pp,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSPPS(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_sp,
                                        const uint32_t* pair_inds_k_for_K_sp,
                                        const uint32_t  pair_inds_count_for_K_sp,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_sp,
                                        const double*   mat_Q_for_K_ps,
                                        const uint32_t* density_inds_for_K_sp,
                                        const uint32_t* density_inds_for_K_ps,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void computeExchangeFockSPPP(double*         mat_K,
                                        const uint32_t* pair_inds_i_for_K_sp,
                                        const uint32_t* pair_inds_k_for_K_sp,
                                        const uint32_t  pair_inds_count_for_K_sp,
                                        const double*   s_prim_info,
                                        const uint32_t* s_prim_aoinds,
                                        const uint32_t  s_prim_count,
                                        const double*   p_prim_info,
                                        const uint32_t* p_prim_aoinds,
                                        const uint32_t  p_prim_count,
                                        const double    max_D,
                                        const double*   mat_D_full_AO,
                                        const double*   mat_Q_for_K_sp,
                                        const double*   mat_Q_for_K_pp,
                                        const uint32_t* density_inds_for_K_sp,
                                        const uint32_t* density_inds_for_K_pp,
                                        const uint32_t  naos,
                                        const double*   boys_func_table,
                                        const double*   boys_func_ft);

__global__ void
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
                        const double    max_D,
                        const double*   mat_D_full_AO,
                        const double*   mat_Q_for_K_ps,
                        const uint32_t* density_inds_for_K_ps,
                        const uint32_t  naos,
                        const double*   boys_func_table,
                        const double*   boys_func_ft);

__global__ void
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
                        const double    max_D,
                        const double*   mat_D_full_AO,
                        const double*   mat_Q_for_K_ps,
                        const double*   mat_Q_for_K_pp,
                        const uint32_t* density_inds_for_K_ps,
                        const uint32_t* density_inds_for_K_pp,
                        const uint32_t  naos,
                        const double*   boys_func_table,
                        const double*   boys_func_ft);

__global__ void
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
                        const double    max_D,
                        const double*   mat_D_full_AO,
                        const double*   mat_Q_for_K_pp,
                        const double*   mat_Q_for_K_ps,
                        const uint32_t* density_inds_for_K_pp,
                        const uint32_t* density_inds_for_K_ps,
                        const uint32_t  naos,
                        const double*   boys_func_table,
                        const double*   boys_func_ft);

__global__ void
computeExchangeFockPPPP(double*         mat_K,
                        const uint32_t* pair_inds_i_for_K_pp,
                        const uint32_t* pair_inds_k_for_K_pp,
                        const uint32_t  pair_inds_count_for_K_pp,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double    max_D,
                        const double*   mat_D_full_AO,
                        const double*   mat_Q_for_K_pp,
                        const uint32_t* density_inds_for_K_pp,
                        const uint32_t  naos,
                        const double*   boys_func_table,
                        const double*   boys_func_ft);

}  // namespace gpu

#endif

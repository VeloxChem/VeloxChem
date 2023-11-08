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

#ifndef EriCoulombFock_hpp
#define EriCoulombFock_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void computeCoulombFockSSSS(double*         mat_J,
          const double*   s_prim_info,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D,
          const double*   mat_Q,
          const uint32_t* first_inds,
          const uint32_t* second_inds,
          const uint32_t  ss_prim_pair_count,
          const double*   boys_func_table,
          const double*   boys_func_ft);

__global__ void computeCoulombFockSSSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   ss_mat_Q,
                       const double*   sp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockSSPP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   ss_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockPPSS(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   ss_mat_D,
                       const double*   ss_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockSPSS(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   ss_mat_D,
                       const double*   ss_mat_Q,
                       const double*   sp_mat_Q,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockSPSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   sp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockSPPP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   sp_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockPPSP(double*         mat_J,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   sp_mat_D,
                       const double*   sp_mat_Q,
                       const double*   pp_mat_Q,
                       const uint32_t* sp_first_inds,
                       const uint32_t* sp_second_inds,
                       const uint32_t  sp_prim_pair_count,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void computeCoulombFockPPPP(double*         mat_J,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double    max_D,
                       const double*   pp_mat_D,
                       const double*   pp_mat_Q,
                       const uint32_t* pp_first_inds,
                       const uint32_t* pp_second_inds,
                       const uint32_t  pp_prim_pair_count,
                       const double*   boys_func_table,
                       const double*   boys_func_ft);

__global__ void cudaFockK(double*         mat_K,
          const uint32_t* pari_inds_i,
          const uint32_t* pari_inds_k,
          const double*   s_prim_info,
          const uint32_t* s_prim_aoinds,
          const uint32_t  s_prim_count,
          const double    max_D,
          const double*   mat_D_full,
          const double*   mat_Q_for_K,
          const uint32_t* first_inds_for_K,
          const uint32_t* second_inds_for_K,
          const uint32_t  ss_prim_pair_count,
          const uint32_t  naos);

}  // namespace gpu

#endif

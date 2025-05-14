//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef OverlapGradient_hpp
#define OverlapGradient_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void
computeOverlapGradientSS(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   ss_mat_W_local,
                         const uint32_t* first_inds_local,
                         const uint32_t* second_inds_local,
                         const uint32_t  ss_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds);

__global__ void
computeOverlapGradientSP(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   ss_mat_W_local,
                         const uint32_t* sp_first_inds_local,
                         const uint32_t* sp_second_inds_local,
                         const uint32_t  sp_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds);

__global__ void
computeOverlapGradientSD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   s_prim_info,
                         const uint32_t  s_prim_count,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   sd_mat_W_local,
                         const uint32_t* sd_first_inds_local,
                         const uint32_t* sd_second_inds_local,
                         const uint32_t  sd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  p_prim_count);

__global__ void
computeOverlapGradientPP(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   pp_mat_W_local,
                         const uint32_t* pp_first_inds_local,
                         const uint32_t* pp_second_inds_local,
                         const uint32_t  pp_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count);

__global__ void
computeOverlapGradientPD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   p_prim_info,
                         const uint32_t  p_prim_count,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   pd_mat_W_local,
                         const uint32_t* pd_first_inds_local,
                         const uint32_t* pd_second_inds_local,
                         const uint32_t  pd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count);

__global__ void
computeOverlapGradientDD(double*         grad_x,
                         const uint32_t  grad_cart_ind,
                         const double*   d_prim_info,
                         const uint32_t  d_prim_count,
                         const double*   dd_mat_W_local,
                         const uint32_t* dd_first_inds_local,
                         const uint32_t* dd_second_inds_local,
                         const uint32_t  dd_prim_pair_count_local,
                         const uint32_t* prim_cart_ao_to_atom_inds,
                         const uint32_t  s_prim_count,
                         const uint32_t  p_prim_count);

}  // namespace gpu

#endif

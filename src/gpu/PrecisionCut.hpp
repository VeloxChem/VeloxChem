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

#ifndef PrecisionCut_hpp
#define PrecisionCut_hpp

#include <vector>
#include <cstdint>

#include "GpuRuntime.hpp"
#include "GpuWrapper.hpp"

// per ij-tile cut; single tile dimension for both ij and kl sides
std::vector<uint32_t>
build_cut_ij_tile(
    const std::vector<double>& Q_ij_local,
    const std::vector<double>& Q_kl,
    const std::vector<double>& D_kl,
    uint32_t ij_count_local,
    uint32_t kl_count,
    int tile_dim,
    double tau);

// per ij-tile cut with separate ij and kl tile dimensions
std::vector<uint32_t>
build_cut_ij_tile(
    const std::vector<double>& Q_ij_local,
    const std::vector<double>& Q_kl,
    const std::vector<double>& D_kl,
    uint32_t ij_count_local,
    uint32_t kl_count,
    int ij_tile_dim,
    int kl_tile_dim,
    double tau);

struct ExchangeCuts {
    std::vector<uint32_t> prec_cut_flat;    // optional host-side flat cut array
    std::vector<uint32_t> screen_cut_flat;  // optional host-side flat cut array
    std::vector<uint32_t> displ_cuts;       // [n_ik] offset into flat arrays
    // std::vector<uint32_t> cut_weights;   // unused: host-only m-tile weight per cut entry
    uint32_t total_cut_entries = 0;         // number of flat cut entries
};

ExchangeCuts
build_exchange_cut_layout(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_counts_AB,
    int      tile_dim_y);

namespace gpu {  // gpu namespace

void
build_exchange_cuts_device(
    uint32_t*       d_prec_cut_flat,
    uint32_t*       d_screen_cut_flat,
    const uint32_t* d_displ_cuts,
    const uint32_t* d_pair_inds_i,
    const uint32_t* d_pair_inds_k,
    const double*   d_Q_K_AB,
    const double*   d_Q_K_CD,
    const uint32_t* d_pair_displs_AB,
    const uint32_t* d_pair_displs_CD,
    const uint32_t* d_pair_counts_AB,
    const uint32_t* d_pair_counts_CD,
    uint32_t        n_ik,
    uint32_t        tile_dim_y,
    uint32_t        tile_dim_x,
    double          max_D,
    double          tau,
    double          eri_threshold,
    gpuStream_t     stream,
    unsigned long long* d_work_counts = nullptr);

}  // namespace gpu

#endif

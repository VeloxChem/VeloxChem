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


#include "SimdElectronRepulsionRecID.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_id_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_id_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2330, 2022, 168, nvalues);

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = nvalues;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto pi = mathconst::pi_value();

            const auto fj = 2.0 * a_norms[i] * b_norms[j] * pi * pi * std::sqrt(pi)
                            / (a_exps[i] * b_exps[j] * std::sqrt(p));

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 37, 0, 7, 8, 16, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 8, 9, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 9, 10, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 10, 11, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 11, 12, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 12, 13, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 13, 14, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 79, 0, 16, 19, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 89, 0, 19, 22, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 99, 0, 22, 25, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 109, 0, 25, 28, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 28, 31, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 129, 0, 37, 43, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 144, 0, 43, 49, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 159, 0, 49, 55, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 174, 0, 55, 61, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 189, 0, 61, 67, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 204, 0, 79, 89, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 225, 0, 89, 99, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 246, 0, 99, 109, 189, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 267, 0, 129, 144, 204, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 295, 0, 144, 159, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 323, 0, 159, 174, 246, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 351, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 354, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 357, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 360, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 363, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 366, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 369, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 378, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 387, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 396, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 405, 3, 13, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 414, 3, 14, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 423, 0, 3, 19, 378, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 441, 0, 3, 22, 387, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 459, 0, 3, 25, 396, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 477, 0, 3, 28, 405, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 495, 0, 3, 31, 414, 73, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 513, 0, 3, 49, 441, 89, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 543, 0, 3, 55, 459, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 573, 0, 3, 61, 477, 109, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 603, 0, 3, 67, 495, 119, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 633, 0, 3, 89, 543, 159, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 678, 0, 3, 99, 573, 174, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 723, 0, 3, 109, 603, 189, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 768, 0, 3, 159, 678, 225, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 831, 0, 3, 174, 723, 246, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 894, 0, 3, 225, 831, 323, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 978, 3, 9, 10, 354, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 984, 3, 10, 11, 357, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 990, 3, 11, 12, 360, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 996, 3, 12, 13, 363, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1002, 3, 13, 14, 366, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1008, 0, 3, 351, 978, 378, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1026, 0, 3, 354, 984, 387, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1044, 0, 3, 357, 990, 396, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1062, 0, 3, 360, 996, 405, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1080, 0, 3, 363, 1002, 414, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1098, 0, 3, 369, 1008, 37, 43, 423,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1134, 0, 3, 378, 1026, 43, 49, 441,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1170, 0, 3, 387, 1044, 49, 55, 459,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1206, 0, 3, 396, 1062, 55, 61, 477,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1242, 0, 3, 405, 1080, 61, 67, 495,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1278, 0, 3, 441, 1170, 79, 89, 543,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1338, 0, 3, 459, 1206, 89, 99, 573,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1398, 0, 3, 477, 1242, 99, 109, 603,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1458, 0, 3, 1098, 1134, 513, 1278, 129,
                                                 144, 633, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1548, 0, 3, 1134, 1170, 543, 1338, 144,
                                                 159, 678, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1638, 0, 3, 1170, 1206, 573, 1398, 159,
                                                 174, 723, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 1728, 0, 3, 1278, 1338, 678, 1638, 204,
                                                 225, 831, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 1854, 0, 3, 1458, 1548, 768, 1728, 267,
                                                 295, 894, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2022, 1854, 168, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 2190, 2022, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 2190, 5, nmax);
}

}  // namespace simdt2ceri

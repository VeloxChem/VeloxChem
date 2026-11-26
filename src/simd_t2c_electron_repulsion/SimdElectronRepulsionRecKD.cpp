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


#include "SimdElectronRepulsionRecKD.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kd_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3530, 3134, 216, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 85, 0, 16, 19, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 95, 0, 19, 22, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 22, 25, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 25, 28, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 28, 31, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 31, 34, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 34, 37, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 155, 0, 43, 49, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 170, 0, 49, 55, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 185, 0, 55, 61, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 200, 0, 61, 67, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 215, 0, 67, 73, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 230, 0, 85, 95, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 251, 0, 95, 105, 170, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 272, 0, 105, 115, 185, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 293, 0, 115, 125, 200, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 314, 0, 125, 135, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 335, 0, 155, 170, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 363, 0, 170, 185, 293, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 391, 0, 185, 200, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 419, 0, 230, 251, 335, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 455, 0, 251, 272, 363, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 491, 0, 272, 293, 391, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 527, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 530, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 533, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 536, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 539, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 542, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 545, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 554, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 563, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 572, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 581, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 590, 3, 14, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 599, 0, 3, 22, 545, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 617, 0, 3, 25, 554, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 635, 0, 3, 28, 563, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 653, 0, 3, 31, 572, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 671, 0, 3, 34, 581, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 689, 0, 3, 37, 590, 79, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 707, 0, 3, 49, 617, 105, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 737, 0, 3, 55, 635, 115, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 767, 0, 3, 61, 653, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 797, 0, 3, 67, 671, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 827, 0, 3, 73, 689, 145, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 857, 0, 3, 105, 737, 170, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 902, 0, 3, 115, 767, 185, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 947, 0, 3, 125, 797, 200, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 992, 0, 3, 135, 827, 215, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1037, 0, 3, 170, 902, 272, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1100, 0, 3, 185, 947, 293, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1163, 0, 3, 200, 992, 314, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1226, 0, 3, 272, 1100, 363, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1310, 0, 3, 293, 1163, 391, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1394, 0, 3, 363, 1310, 491, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1502, 3, 9, 10, 530, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1508, 3, 10, 11, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1514, 3, 11, 12, 536, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1520, 3, 12, 13, 539, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1526, 3, 13, 14, 542, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1532, 0, 3, 527, 1502, 554, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1550, 0, 3, 530, 1508, 563, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1568, 0, 3, 533, 1514, 572, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1586, 0, 3, 536, 1520, 581, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1604, 0, 3, 539, 1526, 590, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1622, 0, 3, 545, 1532, 43, 49, 617,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1658, 0, 3, 554, 1550, 49, 55, 635,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1694, 0, 3, 563, 1568, 55, 61, 653,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1730, 0, 3, 572, 1586, 61, 67, 671,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1766, 0, 3, 581, 1604, 67, 73, 689,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1802, 0, 3, 599, 1622, 85, 95, 707,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1862, 0, 3, 617, 1658, 95, 105, 737,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1922, 0, 3, 635, 1694, 105, 115, 767,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1982, 0, 3, 653, 1730, 115, 125, 797,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2042, 0, 3, 671, 1766, 125, 135, 827,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2102, 0, 3, 1622, 1658, 737, 1922, 155,
                                                 170, 902, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2192, 0, 3, 1658, 1694, 767, 1982, 170,
                                                 185, 947, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2282, 0, 3, 1694, 1730, 797, 2042, 185,
                                                 200, 992, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2372, 0, 3, 1802, 1862, 857, 2102, 230,
                                                 251, 1037, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2498, 0, 3, 1862, 1922, 902, 2192, 251,
                                                 272, 1100, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2624, 0, 3, 1922, 1982, 947, 2282, 272,
                                                 293, 1163, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 2750, 0, 3, 2102, 2192, 1100, 2624, 335,
                                                 363, 1310, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 2918, 0, 3, 2372, 2498, 1226, 2750, 419,
                                                 455, 1394, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3134, 2918, 216, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 3350, 3134, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 3350, 5, nmax);
}

}  // namespace simdt2ceri

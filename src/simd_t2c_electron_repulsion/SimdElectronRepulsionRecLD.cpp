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


#include "SimdElectronRepulsionRecLD.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ld_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ld_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5277, 4782, 270, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 13, 14, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 14, 15, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 15, 16, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 99, 0, 18, 21, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 109, 0, 21, 24, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 24, 27, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 27, 30, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 30, 33, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 33, 36, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 36, 39, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 169, 0, 45, 51, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 184, 0, 51, 57, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 199, 0, 57, 63, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 214, 0, 63, 69, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 229, 0, 69, 75, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 244, 0, 75, 81, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 259, 0, 81, 87, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 274, 0, 99, 109, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 295, 0, 109, 119, 214, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 316, 0, 119, 129, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 337, 0, 129, 139, 244, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 358, 0, 139, 149, 259, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 379, 0, 169, 184, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 407, 0, 184, 199, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 435, 0, 199, 214, 316, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 463, 0, 214, 229, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 491, 0, 229, 244, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 519, 0, 274, 295, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 555, 0, 295, 316, 463, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 591, 0, 316, 337, 491, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 627, 0, 379, 407, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 672, 0, 407, 435, 555, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 717, 0, 435, 463, 591, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 762, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 765, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 768, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 771, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 774, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 777, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 780, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 783, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 786, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 795, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 804, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 813, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 822, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 831, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 840, 3, 15, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 849, 3, 16, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 858, 0, 3, 21, 795, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 876, 0, 3, 24, 804, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 894, 0, 3, 27, 813, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 912, 0, 3, 30, 822, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 930, 0, 3, 33, 831, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 948, 0, 3, 36, 840, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 966, 0, 3, 39, 849, 93, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 984, 0, 3, 57, 876, 109, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1014, 0, 3, 63, 894, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1044, 0, 3, 69, 912, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1074, 0, 3, 75, 930, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1104, 0, 3, 81, 948, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1134, 0, 3, 87, 966, 159, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1164, 0, 3, 109, 1014, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1209, 0, 3, 119, 1044, 214, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1254, 0, 3, 129, 1074, 229, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1299, 0, 3, 139, 1104, 244, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1344, 0, 3, 149, 1134, 259, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1389, 0, 3, 199, 1209, 295, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1452, 0, 3, 214, 1254, 316, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1515, 0, 3, 229, 1299, 337, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1578, 0, 3, 244, 1344, 358, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1641, 0, 3, 295, 1452, 435, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1725, 0, 3, 316, 1515, 463, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1809, 0, 3, 337, 1578, 491, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1893, 0, 3, 435, 1725, 555, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2001, 0, 3, 463, 1809, 591, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 2109, 0, 3, 555, 2001, 717, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2244, 3, 9, 10, 765, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2250, 3, 10, 11, 768, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2256, 3, 11, 12, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2262, 3, 12, 13, 774, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2268, 3, 13, 14, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2274, 3, 14, 15, 780, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2280, 3, 15, 16, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2286, 0, 3, 762, 2244, 795, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2304, 0, 3, 765, 2250, 804, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2322, 0, 3, 768, 2256, 813, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2340, 0, 3, 771, 2262, 822, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2358, 0, 3, 774, 2268, 831, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2376, 0, 3, 777, 2274, 840, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2394, 0, 3, 780, 2280, 849, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2412, 0, 3, 786, 2286, 45, 51, 858,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2448, 0, 3, 795, 2304, 51, 57, 876,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2484, 0, 3, 804, 2322, 57, 63, 894,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2520, 0, 3, 813, 2340, 63, 69, 912,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2556, 0, 3, 822, 2358, 69, 75, 930,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2592, 0, 3, 831, 2376, 75, 81, 948,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2628, 0, 3, 840, 2394, 81, 87, 966,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2664, 0, 3, 876, 2484, 99, 109, 1014,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2724, 0, 3, 894, 2520, 109, 119, 1044,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2784, 0, 3, 912, 2556, 119, 129, 1074,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2844, 0, 3, 930, 2592, 129, 139, 1104,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2904, 0, 3, 948, 2628, 139, 149, 1134,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2964, 0, 3, 2412, 2448, 984, 2664, 169,
                                                 184, 1164, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3054, 0, 3, 2448, 2484, 1014, 2724, 184,
                                                 199, 1209, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3144, 0, 3, 2484, 2520, 1044, 2784, 199,
                                                 214, 1254, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3234, 0, 3, 2520, 2556, 1074, 2844, 214,
                                                 229, 1299, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3324, 0, 3, 2556, 2592, 1104, 2904, 229,
                                                 244, 1344, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3414, 0, 3, 2664, 2724, 1209, 3144, 274,
                                                 295, 1452, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3540, 0, 3, 2724, 2784, 1254, 3234, 295,
                                                 316, 1515, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3666, 0, 3, 2784, 2844, 1299, 3324, 316,
                                                 337, 1578, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3792, 0, 3, 2964, 3054, 1389, 3414, 379,
                                                 407, 1641, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3960, 0, 3, 3054, 3144, 1452, 3540, 407,
                                                 435, 1725, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4128, 0, 3, 3144, 3234, 1515, 3666, 435,
                                                 463, 1809, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 4296, 0, 3, 3414, 3540, 1725, 4128, 519,
                                                 555, 2001, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 4512, 0, 3, 3792, 3960, 1893, 4296, 627,
                                                 672, 2109, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4782, 4512, 270, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 5052, 4782, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 5052, 5, nmax);
}

}  // namespace simdt2ceri

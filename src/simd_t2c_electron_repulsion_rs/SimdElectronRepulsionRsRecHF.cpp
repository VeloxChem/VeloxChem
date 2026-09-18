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


#include "SimdElectronRepulsionRsRecHF.hpp"

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
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_hf_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_hf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5671, 5104, 420, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 15, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 7, 8, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 8, 9, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 9, 10, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 10, 11, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 11, 12, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 12, 13, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 16, 17, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 17, 18, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 18, 19, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 45, 48, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 48, 51, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 51, 54, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 54, 57, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 57, 60, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 238, 0, 66, 72, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 253, 0, 72, 78, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 268, 0, 78, 84, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 283, 0, 84, 90, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 298, 0, 102, 108, 198, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 313, 0, 108, 114, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 328, 0, 114, 120, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 343, 0, 120, 126, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 358, 0, 138, 148, 253, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 379, 0, 148, 158, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 400, 0, 158, 168, 283, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 421, 0, 188, 198, 313, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 442, 0, 198, 208, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 463, 0, 208, 218, 343, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 484, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 487, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 490, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 493, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 496, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 499, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 502, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 505, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 508, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 511, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 514, 3, 9, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 523, 3, 10, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 532, 3, 11, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 541, 3, 12, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 550, 3, 13, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 559, 3, 18, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 568, 3, 19, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 577, 3, 20, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 586, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 595, 3, 22, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 604, 0, 3, 27, 514, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 622, 0, 3, 30, 523, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 640, 0, 3, 33, 532, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 658, 0, 3, 36, 541, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 676, 0, 3, 39, 550, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 694, 0, 3, 48, 559, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 712, 0, 3, 51, 568, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 730, 0, 3, 54, 577, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 748, 0, 3, 57, 586, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 766, 0, 3, 60, 595, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 784, 0, 3, 66, 604, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 814, 0, 3, 72, 622, 148, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 844, 0, 3, 78, 640, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 874, 0, 3, 84, 658, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 904, 0, 3, 90, 676, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 934, 0, 3, 102, 694, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 964, 0, 3, 108, 712, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 994, 0, 3, 114, 730, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1024, 0, 3, 120, 748, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1054, 0, 3, 126, 766, 228, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1084, 0, 3, 148, 844, 253, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1129, 0, 3, 158, 874, 268, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1174, 0, 3, 168, 904, 283, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1219, 0, 3, 198, 994, 313, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1264, 0, 3, 208, 1024, 328, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1309, 0, 3, 218, 1054, 343, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1354, 0, 3, 238, 1084, 358, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1417, 0, 3, 253, 1129, 379, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1480, 0, 3, 268, 1174, 400, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1543, 0, 3, 298, 1219, 421, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1606, 0, 3, 313, 1264, 442, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1669, 0, 3, 328, 1309, 463, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1732, 3, 9, 10, 487, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1738, 3, 10, 11, 490, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1744, 3, 11, 12, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1750, 3, 12, 13, 496, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1756, 3, 18, 19, 502, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1762, 3, 19, 20, 505, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1768, 3, 20, 21, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1774, 3, 21, 22, 511, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1780, 0, 3, 484, 1732, 523, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1798, 0, 3, 487, 1738, 532, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1816, 0, 3, 490, 1744, 541, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1834, 0, 3, 493, 1750, 550, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1852, 0, 3, 499, 1756, 568, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1870, 0, 3, 502, 1762, 577, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1888, 0, 3, 505, 1768, 586, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1906, 0, 3, 508, 1774, 595, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1924, 0, 3, 514, 1780, 66, 72, 622,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1960, 0, 3, 523, 1798, 72, 78, 640,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1996, 0, 3, 532, 1816, 78, 84, 658,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2032, 0, 3, 541, 1834, 84, 90, 676,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2068, 0, 3, 559, 1852, 102, 108, 712,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2104, 0, 3, 568, 1870, 108, 114, 730,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2140, 0, 3, 577, 1888, 114, 120, 748,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2176, 0, 3, 586, 1906, 120, 126, 766,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2212, 0, 3, 622, 1960, 138, 148, 844,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2272, 0, 3, 640, 1996, 148, 158, 874,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2332, 0, 3, 658, 2032, 158, 168, 904,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2392, 0, 3, 712, 2104, 188, 198, 994,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2452, 0, 3, 730, 2140, 198, 208, 1024,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2512, 0, 3, 748, 2176, 208, 218, 1054,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2572, 0, 3, 1924, 1960, 844, 2272, 238,
                                                 253, 1129, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2662, 0, 3, 1960, 1996, 874, 2332, 253,
                                                 268, 1174, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2752, 0, 3, 2068, 2104, 994, 2452, 298,
                                                 313, 1264, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2842, 0, 3, 2104, 2140, 1024, 2512, 313,
                                                 328, 1309, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2932, 0, 3, 2212, 2272, 1129, 2662, 358,
                                                 379, 1480, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3058, 0, 3, 2392, 2452, 1264, 2842, 421,
                                                 442, 1669, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3184, 3, 484, 487, 1738, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3194, 3, 487, 490, 1744, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3204, 3, 490, 493, 1750, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3214, 3, 499, 502, 1762, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3224, 3, 502, 505, 1768, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3234, 3, 505, 508, 1774, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3244, 0, 3, 1732, 3184, 1798, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3274, 0, 3, 1738, 3194, 1816, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3304, 0, 3, 1744, 3204, 1834, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3334, 0, 3, 1756, 3214, 1870, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3364, 0, 3, 1762, 3224, 1888, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3394, 0, 3, 1768, 3234, 1906, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3424, 0, 3, 1780, 3244, 604, 622, 1960,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3484, 0, 3, 1798, 3274, 622, 640, 1996,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3544, 0, 3, 1816, 3304, 640, 658, 2032,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3604, 0, 3, 1852, 3334, 694, 712, 2104,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3664, 0, 3, 1870, 3364, 712, 730, 2140,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3724, 0, 3, 1888, 3394, 730, 748, 2176,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3784, 0, 3, 1924, 3424, 784, 814, 2212,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3884, 0, 3, 1960, 3484, 814, 844, 2272,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3984, 0, 3, 1996, 3544, 844, 874, 2332,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4084, 0, 3, 2068, 3604, 934, 964, 2392,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4184, 0, 3, 2104, 3664, 964, 994, 2452,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4284, 0, 3, 2140, 3724, 994, 1024, 2512,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4384, 0, 3, 3424, 3484, 2272, 3984,
                                                 1084, 1129, 2662, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4534, 0, 3, 3604, 3664, 2452, 4284,
                                                 1219, 1264, 2842, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 4684, 0, 3, 3784, 3884, 2572, 4384,
                                                 1354, 1417, 2932, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 4894, 0, 3, 4084, 4184, 2752, 4534,
                                                 1543, 1606, 3058, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5104, 4684, 420, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 5524, 5314, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 5524, 7, nmax);

    simdtrf::transform_f_inner(buffer, 5524, 5104, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 77 * nvalues, nvalues, buffer, 5524, 7, nmax);
}

}  // namespace simdt2ceri

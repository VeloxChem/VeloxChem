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


#include "SimdElectronRepulsionRsRecID.hpp"

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
compute_rs_id_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_id_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4514, 4038, 336, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 8, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 16, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 9, 10, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 10, 11, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 11, 12, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 12, 13, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 13, 14, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 17, 18, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 18, 19, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 19, 20, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 20, 21, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 21, 22, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 22, 23, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 23, 24, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 47, 50, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 50, 53, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 53, 56, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 56, 59, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 59, 62, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 252, 0, 68, 74, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 267, 0, 74, 80, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 282, 0, 80, 86, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 297, 0, 86, 92, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 312, 0, 92, 98, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 327, 0, 110, 116, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 342, 0, 116, 122, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 357, 0, 122, 128, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 372, 0, 128, 134, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 387, 0, 134, 140, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 402, 0, 152, 162, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 423, 0, 162, 172, 297, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 444, 0, 172, 182, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 465, 0, 202, 212, 357, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 486, 0, 212, 222, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 507, 0, 222, 232, 387, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 528, 0, 252, 267, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 556, 0, 267, 282, 423, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 584, 0, 282, 297, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 612, 0, 327, 342, 465, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 640, 0, 342, 357, 486, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 668, 0, 357, 372, 507, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 696, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 699, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 702, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 705, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 708, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 711, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 714, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 717, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 720, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 723, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 726, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 729, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 732, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 741, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 750, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 759, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 768, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 777, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 786, 3, 19, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 795, 3, 20, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 804, 3, 21, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 813, 3, 22, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 822, 3, 23, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 831, 3, 24, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 840, 0, 3, 29, 741, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 858, 0, 3, 32, 750, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 876, 0, 3, 35, 759, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 894, 0, 3, 38, 768, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 912, 0, 3, 41, 777, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 930, 0, 3, 50, 795, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 948, 0, 3, 53, 804, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 966, 0, 3, 56, 813, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 984, 0, 3, 59, 822, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1002, 0, 3, 62, 831, 146, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1020, 0, 3, 80, 858, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1050, 0, 3, 86, 876, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1080, 0, 3, 92, 894, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1110, 0, 3, 98, 912, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1140, 0, 3, 122, 948, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1170, 0, 3, 128, 966, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1200, 0, 3, 134, 984, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1230, 0, 3, 140, 1002, 242, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1260, 0, 3, 162, 1050, 282, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1305, 0, 3, 172, 1080, 297, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1350, 0, 3, 182, 1110, 312, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1395, 0, 3, 212, 1170, 357, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1440, 0, 3, 222, 1200, 372, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1485, 0, 3, 232, 1230, 387, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1530, 0, 3, 282, 1305, 423, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1593, 0, 3, 297, 1350, 444, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1656, 0, 3, 357, 1440, 486, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1719, 0, 3, 372, 1485, 507, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1782, 0, 3, 423, 1593, 584, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1866, 0, 3, 486, 1719, 668, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1950, 3, 9, 10, 699, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1956, 3, 10, 11, 702, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1962, 3, 11, 12, 705, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1968, 3, 12, 13, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1974, 3, 13, 14, 711, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1980, 3, 19, 20, 717, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1986, 3, 20, 21, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1992, 3, 21, 22, 723, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1998, 3, 22, 23, 726, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2004, 3, 23, 24, 729, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2010, 0, 3, 696, 1950, 741, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2028, 0, 3, 699, 1956, 750, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2046, 0, 3, 702, 1962, 759, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2064, 0, 3, 705, 1968, 768, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2082, 0, 3, 708, 1974, 777, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2100, 0, 3, 714, 1980, 795, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2118, 0, 3, 717, 1986, 804, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2136, 0, 3, 720, 1992, 813, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2154, 0, 3, 723, 1998, 822, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2172, 0, 3, 726, 2004, 831, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2190, 0, 3, 732, 2010, 68, 74, 840,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2226, 0, 3, 741, 2028, 74, 80, 858,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2262, 0, 3, 750, 2046, 80, 86, 876,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2298, 0, 3, 759, 2064, 86, 92, 894,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2334, 0, 3, 768, 2082, 92, 98, 912,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2370, 0, 3, 786, 2100, 110, 116, 930,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2406, 0, 3, 795, 2118, 116, 122, 948,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2442, 0, 3, 804, 2136, 122, 128, 966,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2478, 0, 3, 813, 2154, 128, 134, 984,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2514, 0, 3, 822, 2172, 134, 140, 1002,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2550, 0, 3, 858, 2262, 152, 162, 1050,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2610, 0, 3, 876, 2298, 162, 172, 1080,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2670, 0, 3, 894, 2334, 172, 182, 1110,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2730, 0, 3, 948, 2442, 202, 212, 1170,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2790, 0, 3, 966, 2478, 212, 222, 1200,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2850, 0, 3, 984, 2514, 222, 232, 1230,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2910, 0, 3, 2190, 2226, 1020, 2550, 252,
                                                 267, 1260, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3000, 0, 3, 2226, 2262, 1050, 2610, 267,
                                                 282, 1305, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3090, 0, 3, 2262, 2298, 1080, 2670, 282,
                                                 297, 1350, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3180, 0, 3, 2370, 2406, 1140, 2730, 327,
                                                 342, 1395, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3270, 0, 3, 2406, 2442, 1170, 2790, 342,
                                                 357, 1440, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3360, 0, 3, 2442, 2478, 1200, 2850, 357,
                                                 372, 1485, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3450, 0, 3, 2550, 2610, 1305, 3090, 402,
                                                 423, 1593, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3576, 0, 3, 2730, 2790, 1440, 3360, 465,
                                                 486, 1719, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3702, 0, 3, 2910, 3000, 1530, 3450, 528,
                                                 556, 1782, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 3870, 0, 3, 3180, 3270, 1656, 3576, 612,
                                                 640, 1866, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4038, 3702, 336, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 4374, 4206, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 4374, 5, nmax);

    simdtrf::transform_d_inner(buffer, 4374, 4038, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 65 * nvalues, nvalues, buffer, 4374, 5, nmax);
}

}  // namespace simdt2ceri

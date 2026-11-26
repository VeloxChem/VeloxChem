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


#include "SimdElectronRepulsionGeom10RsRecHF.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_hf_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_hf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 11471, 10064, 1260, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 16, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 7, 8, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 8, 9, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 9, 10, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 10, 11, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 11, 12, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 12, 13, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 13, 14, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 18, 19, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 19, 20, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 20, 21, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 21, 22, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 22, 23, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 23, 24, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 50, 53, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 53, 56, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 56, 59, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 59, 62, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 62, 65, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 65, 68, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 278, 0, 74, 80, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 293, 0, 80, 86, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 308, 0, 86, 92, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 323, 0, 92, 98, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 338, 0, 98, 104, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 353, 0, 116, 122, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 368, 0, 122, 128, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 383, 0, 128, 134, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 398, 0, 134, 140, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 413, 0, 140, 146, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 428, 0, 158, 168, 293, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 449, 0, 168, 178, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 470, 0, 178, 188, 323, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 491, 0, 188, 198, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 512, 0, 218, 228, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 533, 0, 228, 238, 383, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 554, 0, 238, 248, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 575, 0, 248, 258, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 596, 0, 278, 293, 449, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 624, 0, 293, 308, 470, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 652, 0, 308, 323, 491, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 680, 0, 353, 368, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 708, 0, 368, 383, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 736, 0, 383, 398, 575, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 764, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 767, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 770, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 773, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 776, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 779, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 782, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 785, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 788, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 791, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 794, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 797, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 800, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 803, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 806, 3, 8, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 815, 3, 9, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 824, 3, 10, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 833, 3, 11, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 842, 3, 12, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 851, 3, 13, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 860, 3, 14, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 869, 3, 18, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 878, 3, 19, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 887, 3, 20, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 896, 3, 21, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 905, 3, 22, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 914, 3, 23, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 923, 3, 24, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 932, 0, 3, 26, 806, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 950, 0, 3, 29, 815, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 968, 0, 3, 32, 824, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 986, 0, 3, 35, 833, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1004, 0, 3, 38, 842, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1022, 0, 3, 41, 851, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1040, 0, 3, 44, 860, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1058, 0, 3, 50, 869, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1076, 0, 3, 53, 878, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1094, 0, 3, 56, 887, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1112, 0, 3, 59, 896, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1130, 0, 3, 62, 905, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1148, 0, 3, 65, 914, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1166, 0, 3, 68, 923, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1184, 0, 3, 80, 968, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1214, 0, 3, 86, 986, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1244, 0, 3, 92, 1004, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1274, 0, 3, 98, 1022, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1304, 0, 3, 104, 1040, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1334, 0, 3, 122, 1094, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1364, 0, 3, 128, 1112, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1394, 0, 3, 134, 1130, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1424, 0, 3, 140, 1148, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1454, 0, 3, 146, 1166, 268, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1484, 0, 3, 158, 1184, 278, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1529, 0, 3, 168, 1214, 293, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1574, 0, 3, 178, 1244, 308, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1619, 0, 3, 188, 1274, 323, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1664, 0, 3, 198, 1304, 338, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1709, 0, 3, 218, 1334, 353, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1754, 0, 3, 228, 1364, 368, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1799, 0, 3, 238, 1394, 383, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1844, 0, 3, 248, 1424, 398, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1889, 0, 3, 258, 1454, 413, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1934, 0, 3, 293, 1574, 449, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1997, 0, 3, 308, 1619, 470, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2060, 0, 3, 323, 1664, 491, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2123, 0, 3, 368, 1799, 533, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2186, 0, 3, 383, 1844, 554, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2249, 0, 3, 398, 1889, 575, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2312, 0, 3, 428, 1934, 596, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2396, 0, 3, 449, 1997, 624, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2480, 0, 3, 470, 2060, 652, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2564, 0, 3, 512, 2123, 680, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2648, 0, 3, 533, 2186, 708, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2732, 0, 3, 554, 2249, 736, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2816, 3, 8, 9, 767, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2822, 3, 9, 10, 770, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2828, 3, 10, 11, 773, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2834, 3, 11, 12, 776, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2840, 3, 12, 13, 779, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2846, 3, 13, 14, 782, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2852, 3, 18, 19, 788, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2858, 3, 19, 20, 791, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2864, 3, 20, 21, 794, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2870, 3, 21, 22, 797, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2876, 3, 22, 23, 800, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2882, 3, 23, 24, 803, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2888, 0, 3, 764, 2816, 815, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2906, 0, 3, 767, 2822, 824, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2924, 0, 3, 770, 2828, 833, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2942, 0, 3, 773, 2834, 842, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2960, 0, 3, 776, 2840, 851, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2978, 0, 3, 779, 2846, 860, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2996, 0, 3, 785, 2852, 878, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3014, 0, 3, 788, 2858, 887, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3032, 0, 3, 791, 2864, 896, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3050, 0, 3, 794, 2870, 905, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3068, 0, 3, 797, 2876, 914, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3086, 0, 3, 800, 2882, 923, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3104, 0, 3, 815, 2906, 74, 80, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3140, 0, 3, 824, 2924, 80, 86, 986,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3176, 0, 3, 833, 2942, 86, 92, 1004,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3212, 0, 3, 842, 2960, 92, 98, 1022,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3248, 0, 3, 851, 2978, 98, 104, 1040,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3284, 0, 3, 878, 3014, 116, 122, 1094,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3320, 0, 3, 887, 3032, 122, 128, 1112,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3356, 0, 3, 896, 3050, 128, 134, 1130,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3392, 0, 3, 905, 3068, 134, 140, 1148,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3428, 0, 3, 914, 3086, 140, 146, 1166,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3464, 0, 3, 968, 3140, 158, 168, 1214,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3524, 0, 3, 986, 3176, 168, 178, 1244,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3584, 0, 3, 1004, 3212, 178, 188, 1274,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3644, 0, 3, 1022, 3248, 188, 198, 1304,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3704, 0, 3, 1094, 3320, 218, 228, 1364,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3764, 0, 3, 1112, 3356, 228, 238, 1394,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3824, 0, 3, 1130, 3392, 238, 248, 1424,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3884, 0, 3, 1148, 3428, 248, 258, 1454,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3944, 0, 3, 3104, 3140, 1214, 3524, 278,
                                                 293, 1574, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4034, 0, 3, 3140, 3176, 1244, 3584, 293,
                                                 308, 1619, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4124, 0, 3, 3176, 3212, 1274, 3644, 308,
                                                 323, 1664, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4214, 0, 3, 3284, 3320, 1364, 3764, 353,
                                                 368, 1799, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4304, 0, 3, 3320, 3356, 1394, 3824, 368,
                                                 383, 1844, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4394, 0, 3, 3356, 3392, 1424, 3884, 383,
                                                 398, 1889, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4484, 0, 3, 3464, 3524, 1574, 4034, 428,
                                                 449, 1997, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4610, 0, 3, 3524, 3584, 1619, 4124, 449,
                                                 470, 2060, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4736, 0, 3, 3704, 3764, 1799, 4304, 512,
                                                 533, 2186, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4862, 0, 3, 3764, 3824, 1844, 4394, 533,
                                                 554, 2249, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 4988, 0, 3, 3944, 4034, 1997, 4610, 596,
                                                 624, 2480, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5156, 0, 3, 4214, 4304, 2186, 4862, 680,
                                                 708, 2732, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5324, 3, 764, 767, 2822, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5334, 3, 767, 770, 2828, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5344, 3, 770, 773, 2834, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5354, 3, 773, 776, 2840, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5364, 3, 776, 779, 2846, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5374, 3, 785, 788, 2858, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5384, 3, 788, 791, 2864, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5394, 3, 791, 794, 2870, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5404, 3, 794, 797, 2876, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5414, 3, 797, 800, 2882, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 5424, 0, 3, 2816, 5324, 2906, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5454, 0, 3, 2822, 5334, 2924, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5484, 0, 3, 2828, 5344, 2942, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5514, 0, 3, 2834, 5354, 2960, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5544, 0, 3, 2840, 5364, 2978, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5574, 0, 3, 2852, 5374, 3014, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5604, 0, 3, 2858, 5384, 3032, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5634, 0, 3, 2864, 5394, 3050, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5664, 0, 3, 2870, 5404, 3068, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5694, 0, 3, 2876, 5414, 3086, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 5724, 0, 3, 2888, 5424, 932, 950, 3104,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5784, 0, 3, 2906, 5454, 950, 968, 3140,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5844, 0, 3, 2924, 5484, 968, 986, 3176,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5904, 0, 3, 2942, 5514, 986, 1004, 3212,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5964, 0, 3, 2960, 5544, 1004, 1022,
                                                 3248, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6024, 0, 3, 2996, 5574, 1058, 1076,
                                                 3284, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6084, 0, 3, 3014, 5604, 1076, 1094,
                                                 3320, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6144, 0, 3, 3032, 5634, 1094, 1112,
                                                 3356, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6204, 0, 3, 3050, 5664, 1112, 1130,
                                                 3392, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6264, 0, 3, 3068, 5694, 1130, 1148,
                                                 3428, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6324, 0, 3, 3140, 5844, 1184, 1214,
                                                 3524, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6424, 0, 3, 3176, 5904, 1214, 1244,
                                                 3584, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6524, 0, 3, 3212, 5964, 1244, 1274,
                                                 3644, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6624, 0, 3, 3320, 6144, 1334, 1364,
                                                 3764, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6724, 0, 3, 3356, 6204, 1364, 1394,
                                                 3824, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6824, 0, 3, 3392, 6264, 1394, 1424,
                                                 3884, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6924, 0, 3, 5724, 5784, 3464, 6324,
                                                 1484, 1529, 3944, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7074, 0, 3, 5784, 5844, 3524, 6424,
                                                 1529, 1574, 4034, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7224, 0, 3, 5844, 5904, 3584, 6524,
                                                 1574, 1619, 4124, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7374, 0, 3, 6024, 6084, 3704, 6624,
                                                 1709, 1754, 4214, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7524, 0, 3, 6084, 6144, 3764, 6724,
                                                 1754, 1799, 4304, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7674, 0, 3, 6144, 6204, 3824, 6824,
                                                 1799, 1844, 4394, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 7824, 0, 3, 6324, 6424, 4034, 7224,
                                                 1934, 1997, 4610, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8034, 0, 3, 6624, 6724, 4304, 7674,
                                                 2123, 2186, 4862, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8244, 0, 3, 6924, 7074, 4484, 7824,
                                                 2312, 2396, 4988, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 8524, 0, 3, 7374, 7524, 4736, 8034,
                                                 2564, 2648, 5156, ncols, alpha, beta, p);

            simdgeo::geom_h_x(buffer, 8804, 7374, 8524, 1, 10, ncols, alpha);

            simdgeo::geom_h_y(buffer, 9014, 7374, 8524, 1, 10, ncols, alpha);

            simdgeo::geom_h_z(buffer, 9224, 7374, 8524, 1, 10, ncols, alpha);

            simdgeo::geom_h_x(buffer, 9434, 6924, 8244, 1, 10, ncols, alpha);

            simdgeo::geom_h_y(buffer, 9644, 6924, 8244, 1, 10, ncols, alpha);

            simdgeo::geom_h_z(buffer, 9854, 6924, 8244, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 10064, 9434, 630, ncols);

            simdfunc::contract_primitives(buffer, 10694, 8804, 630, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 11324, 10694, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 11324, 7, nmax);

    simdtrf::transform_f_inner(buffer, 11324, 10904, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 77 * nvalues, nvalues, buffer, 11324, 7, nmax);

    simdtrf::transform_f_inner(buffer, 11324, 11114, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 154 * nvalues, nvalues, buffer, 11324, 7, nmax);

    simdtrf::transform_f_inner(buffer, 11324, 10064, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 231 * nvalues, nvalues, buffer, 11324, 7, nmax);

    simdtrf::transform_f_inner(buffer, 11324, 10274, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 308 * nvalues, nvalues, buffer, 11324, 7, nmax);

    simdtrf::transform_f_inner(buffer, 11324, 10484, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 385 * nvalues, nvalues, buffer, 11324, 7, nmax);
}

}  // namespace simdt2ceri

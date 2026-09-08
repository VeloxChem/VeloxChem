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


#include "SimdElectronRepulsionRecIK.hpp"

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
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ik_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ik_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 37506, 36078, 1008, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 7, 8, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 8, 9, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 9, 10, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 10, 11, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 11, 12, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 12, 13, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 13, 14, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 14, 15, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 15, 16, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 16, 17, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 17, 18, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 20, 23, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 23, 26, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 26, 29, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 29, 32, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 32, 35, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 35, 38, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 38, 41, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 41, 44, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 44, 47, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 47, 50, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 222, 0, 56, 62, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 237, 0, 62, 68, 142, ncols, alpha, beta,
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

            compute_prim_gs_electron_repulsion_0(buffer, 327, 0, 98, 104, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 342, 0, 104, 110, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 357, 0, 122, 132, 237, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 378, 0, 132, 142, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 399, 0, 142, 152, 267, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 420, 0, 152, 162, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 441, 0, 162, 172, 297, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 462, 0, 172, 182, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 483, 0, 182, 192, 327, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 504, 0, 192, 202, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 525, 0, 222, 237, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 553, 0, 237, 252, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 252, 267, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 267, 282, 441, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 282, 297, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 297, 312, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 312, 327, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 721, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 724, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 727, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 730, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 733, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 736, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 739, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 742, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 745, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 748, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 751, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 754, 3, 8, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 763, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 772, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 781, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 790, 3, 12, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 799, 3, 13, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 808, 3, 14, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 817, 3, 15, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 826, 3, 16, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 835, 3, 17, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 844, 3, 18, 53, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 853, 0, 3, 20, 754, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 871, 0, 3, 23, 763, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 889, 0, 3, 26, 772, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 907, 0, 3, 29, 781, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 925, 0, 3, 32, 790, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 943, 0, 3, 35, 799, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 961, 0, 3, 38, 808, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 979, 0, 3, 41, 817, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 997, 0, 3, 44, 826, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1015, 0, 3, 47, 835, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1033, 0, 3, 50, 844, 116, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1051, 0, 3, 62, 889, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1081, 0, 3, 68, 907, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1111, 0, 3, 74, 925, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1141, 0, 3, 80, 943, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1171, 0, 3, 86, 961, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1201, 0, 3, 92, 979, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1231, 0, 3, 98, 997, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1261, 0, 3, 104, 1015, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1291, 0, 3, 110, 1033, 212, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1321, 0, 3, 122, 1051, 222, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1366, 0, 3, 132, 1081, 237, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1411, 0, 3, 142, 1111, 252, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1456, 0, 3, 152, 1141, 267, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1501, 0, 3, 162, 1171, 282, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1546, 0, 3, 172, 1201, 297, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1591, 0, 3, 182, 1231, 312, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1636, 0, 3, 192, 1261, 327, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1681, 0, 3, 202, 1291, 342, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1726, 0, 3, 237, 1411, 378, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1789, 0, 3, 252, 1456, 399, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1852, 0, 3, 267, 1501, 420, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1915, 0, 3, 282, 1546, 441, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1978, 0, 3, 297, 1591, 462, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2041, 0, 3, 312, 1636, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2104, 0, 3, 327, 1681, 504, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2167, 0, 3, 357, 1726, 525, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2251, 0, 3, 378, 1789, 553, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2335, 0, 3, 399, 1852, 581, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2419, 0, 3, 420, 1915, 609, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2503, 0, 3, 441, 1978, 637, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2587, 0, 3, 462, 2041, 665, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2671, 0, 3, 483, 2104, 693, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2755, 3, 8, 9, 724, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2761, 3, 9, 10, 727, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2767, 3, 10, 11, 730, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2773, 3, 11, 12, 733, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2779, 3, 12, 13, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2785, 3, 13, 14, 739, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2791, 3, 14, 15, 742, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2797, 3, 15, 16, 745, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2803, 3, 16, 17, 748, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2809, 3, 17, 18, 751, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2815, 0, 3, 721, 2755, 763, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2833, 0, 3, 724, 2761, 772, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2851, 0, 3, 727, 2767, 781, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2869, 0, 3, 730, 2773, 790, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2887, 0, 3, 733, 2779, 799, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2905, 0, 3, 736, 2785, 808, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2923, 0, 3, 739, 2791, 817, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2941, 0, 3, 742, 2797, 826, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2959, 0, 3, 745, 2803, 835, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2977, 0, 3, 748, 2809, 844, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2995, 0, 3, 763, 2833, 56, 62, 889,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3031, 0, 3, 772, 2851, 62, 68, 907,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3067, 0, 3, 781, 2869, 68, 74, 925,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3103, 0, 3, 790, 2887, 74, 80, 943,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3139, 0, 3, 799, 2905, 80, 86, 961,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3175, 0, 3, 808, 2923, 86, 92, 979,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3211, 0, 3, 817, 2941, 92, 98, 997,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3247, 0, 3, 826, 2959, 98, 104, 1015,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3283, 0, 3, 835, 2977, 104, 110, 1033,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3319, 0, 3, 889, 3031, 122, 132, 1081,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3379, 0, 3, 907, 3067, 132, 142, 1111,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3439, 0, 3, 925, 3103, 142, 152, 1141,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3499, 0, 3, 943, 3139, 152, 162, 1171,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3559, 0, 3, 961, 3175, 162, 172, 1201,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3619, 0, 3, 979, 3211, 172, 182, 1231,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3679, 0, 3, 997, 3247, 182, 192, 1261,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3739, 0, 3, 1015, 3283, 192, 202, 1291,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3799, 0, 3, 2995, 3031, 1081, 3379, 222,
                                                 237, 1411, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3889, 0, 3, 3031, 3067, 1111, 3439, 237,
                                                 252, 1456, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3979, 0, 3, 3067, 3103, 1141, 3499, 252,
                                                 267, 1501, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4069, 0, 3, 3103, 3139, 1171, 3559, 267,
                                                 282, 1546, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4159, 0, 3, 3139, 3175, 1201, 3619, 282,
                                                 297, 1591, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4249, 0, 3, 3175, 3211, 1231, 3679, 297,
                                                 312, 1636, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4339, 0, 3, 3211, 3247, 1261, 3739, 312,
                                                 327, 1681, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4429, 0, 3, 3319, 3379, 1411, 3889, 357,
                                                 378, 1789, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4555, 0, 3, 3379, 3439, 1456, 3979, 378,
                                                 399, 1852, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4681, 0, 3, 3439, 3499, 1501, 4069, 399,
                                                 420, 1915, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4807, 0, 3, 3499, 3559, 1546, 4159, 420,
                                                 441, 1978, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 4933, 0, 3, 3559, 3619, 1591, 4249, 441,
                                                 462, 2041, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5059, 0, 3, 3619, 3679, 1636, 4339, 462,
                                                 483, 2104, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5185, 0, 3, 3799, 3889, 1789, 4555, 525,
                                                 553, 2335, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5353, 0, 3, 3889, 3979, 1852, 4681, 553,
                                                 581, 2419, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5521, 0, 3, 3979, 4069, 1915, 4807, 581,
                                                 609, 2503, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5689, 0, 3, 4069, 4159, 1978, 4933, 609,
                                                 637, 2587, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5857, 0, 3, 4159, 4249, 2041, 5059, 637,
                                                 665, 2671, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6025, 3, 721, 724, 2761, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6035, 3, 724, 727, 2767, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6045, 3, 727, 730, 2773, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6055, 3, 730, 733, 2779, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6065, 3, 733, 736, 2785, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6075, 3, 736, 739, 2791, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6085, 3, 739, 742, 2797, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6095, 3, 742, 745, 2803, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6105, 3, 745, 748, 2809, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 6115, 0, 3, 2755, 6025, 2833, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6145, 0, 3, 2761, 6035, 2851, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6175, 0, 3, 2767, 6045, 2869, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6205, 0, 3, 2773, 6055, 2887, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6235, 0, 3, 2779, 6065, 2905, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6265, 0, 3, 2785, 6075, 2923, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6295, 0, 3, 2791, 6085, 2941, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6325, 0, 3, 2797, 6095, 2959, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6355, 0, 3, 2803, 6105, 2977, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 6385, 0, 3, 2815, 6115, 853, 871, 2995,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6445, 0, 3, 2833, 6145, 871, 889, 3031,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6505, 0, 3, 2851, 6175, 889, 907, 3067,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6565, 0, 3, 2869, 6205, 907, 925, 3103,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6625, 0, 3, 2887, 6235, 925, 943, 3139,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6685, 0, 3, 2905, 6265, 943, 961, 3175,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6745, 0, 3, 2923, 6295, 961, 979, 3211,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6805, 0, 3, 2941, 6325, 979, 997, 3247,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6865, 0, 3, 2959, 6355, 997, 1015, 3283,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 6925, 0, 3, 3031, 6505, 1051, 1081,
                                                 3379, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7025, 0, 3, 3067, 6565, 1081, 1111,
                                                 3439, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7125, 0, 3, 3103, 6625, 1111, 1141,
                                                 3499, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7225, 0, 3, 3139, 6685, 1141, 1171,
                                                 3559, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7325, 0, 3, 3175, 6745, 1171, 1201,
                                                 3619, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7425, 0, 3, 3211, 6805, 1201, 1231,
                                                 3679, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7525, 0, 3, 3247, 6865, 1231, 1261,
                                                 3739, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7625, 0, 3, 6385, 6445, 3319, 6925,
                                                 1321, 1366, 3799, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7775, 0, 3, 6445, 6505, 3379, 7025,
                                                 1366, 1411, 3889, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 7925, 0, 3, 6505, 6565, 3439, 7125,
                                                 1411, 1456, 3979, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8075, 0, 3, 6565, 6625, 3499, 7225,
                                                 1456, 1501, 4069, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8225, 0, 3, 6625, 6685, 3559, 7325,
                                                 1501, 1546, 4159, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8375, 0, 3, 6685, 6745, 3619, 7425,
                                                 1546, 1591, 4249, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8525, 0, 3, 6745, 6805, 3679, 7525,
                                                 1591, 1636, 4339, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8675, 0, 3, 6925, 7025, 3889, 7925,
                                                 1726, 1789, 4555, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 8885, 0, 3, 7025, 7125, 3979, 8075,
                                                 1789, 1852, 4681, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9095, 0, 3, 7125, 7225, 4069, 8225,
                                                 1852, 1915, 4807, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9305, 0, 3, 7225, 7325, 4159, 8375,
                                                 1915, 1978, 4933, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9515, 0, 3, 7325, 7425, 4249, 8525,
                                                 1978, 2041, 5059, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 9725, 0, 3, 7625, 7775, 4429, 8675,
                                                 2167, 2251, 5185, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 10005, 0, 3, 7775, 7925, 4555, 8885,
                                                 2251, 2335, 5353, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 10285, 0, 3, 7925, 8075, 4681, 9095,
                                                 2335, 2419, 5521, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 10565, 0, 3, 8075, 8225, 4807, 9305,
                                                 2419, 2503, 5689, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 10845, 0, 3, 8225, 8375, 4933, 9515,
                                                 2503, 2587, 5857, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11125, 3, 2755, 2761, 6035, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11140, 3, 2761, 2767, 6045, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11155, 3, 2767, 2773, 6055, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11170, 3, 2773, 2779, 6065, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11185, 3, 2779, 2785, 6075, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11200, 3, 2785, 2791, 6085, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11215, 3, 2791, 2797, 6095, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11230, 3, 2797, 2803, 6105, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 11245, 0, 3, 6025, 11125, 6145, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11290, 0, 3, 6035, 11140, 6175, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11335, 0, 3, 6045, 11155, 6205, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11380, 0, 3, 6055, 11170, 6235, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11425, 0, 3, 6065, 11185, 6265, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11470, 0, 3, 6075, 11200, 6295, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11515, 0, 3, 6085, 11215, 6325, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11560, 0, 3, 6095, 11230, 6355, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 11605, 0, 3, 6145, 11290, 2995, 3031,
                                                 6505, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11695, 0, 3, 6175, 11335, 3031, 3067,
                                                 6565, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11785, 0, 3, 6205, 11380, 3067, 3103,
                                                 6625, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11875, 0, 3, 6235, 11425, 3103, 3139,
                                                 6685, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11965, 0, 3, 6265, 11470, 3139, 3175,
                                                 6745, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12055, 0, 3, 6295, 11515, 3175, 3211,
                                                 6805, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12145, 0, 3, 6325, 11560, 3211, 3247,
                                                 6865, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12235, 0, 3, 6505, 11695, 3319, 3379,
                                                 7025, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12385, 0, 3, 6565, 11785, 3379, 3439,
                                                 7125, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12535, 0, 3, 6625, 11875, 3439, 3499,
                                                 7225, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12685, 0, 3, 6685, 11965, 3499, 3559,
                                                 7325, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12835, 0, 3, 6745, 12055, 3559, 3619,
                                                 7425, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12985, 0, 3, 6805, 12145, 3619, 3679,
                                                 7525, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13135, 0, 3, 11605, 11695, 7025, 12385,
                                                 3799, 3889, 7925, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13360, 0, 3, 11695, 11785, 7125, 12535,
                                                 3889, 3979, 8075, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13585, 0, 3, 11785, 11875, 7225, 12685,
                                                 3979, 4069, 8225, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13810, 0, 3, 11875, 11965, 7325, 12835,
                                                 4069, 4159, 8375, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14035, 0, 3, 11965, 12055, 7425, 12985,
                                                 4159, 4249, 8525, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14260, 0, 3, 12235, 12385, 7925, 13360,
                                                 4429, 4555, 8885, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14575, 0, 3, 12385, 12535, 8075, 13585,
                                                 4555, 4681, 9095, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 14890, 0, 3, 12535, 12685, 8225, 13810,
                                                 4681, 4807, 9305, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 15205, 0, 3, 12685, 12835, 8375, 14035,
                                                 4807, 4933, 9515, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 15520, 0, 3, 13135, 13360, 8885, 14575,
                                                 5185, 5353, 10285, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 15940, 0, 3, 13360, 13585, 9095, 14890,
                                                 5353, 5521, 10565, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 16360, 0, 3, 13585, 13810, 9305, 15205,
                                                 5521, 5689, 10845, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16780, 3, 6025, 6035, 11140, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16801, 3, 6035, 6045, 11155, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16822, 3, 6045, 6055, 11170, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16843, 3, 6055, 6065, 11185, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16864, 3, 6065, 6075, 11200, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16885, 3, 6075, 6085, 11215, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 16906, 3, 6085, 6095, 11230, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 16927, 0, 3, 11125, 16780, 11290, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16990, 0, 3, 11140, 16801, 11335, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17053, 0, 3, 11155, 16822, 11380, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17116, 0, 3, 11170, 16843, 11425, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17179, 0, 3, 11185, 16864, 11470, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17242, 0, 3, 11200, 16885, 11515, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17305, 0, 3, 11215, 16906, 11560, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 17368, 0, 3, 11245, 16927, 6385, 6445,
                                                 11605, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17494, 0, 3, 11290, 16990, 6445, 6505,
                                                 11695, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17620, 0, 3, 11335, 17053, 6505, 6565,
                                                 11785, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17746, 0, 3, 11380, 17116, 6565, 6625,
                                                 11875, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17872, 0, 3, 11425, 17179, 6625, 6685,
                                                 11965, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17998, 0, 3, 11470, 17242, 6685, 6745,
                                                 12055, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 18124, 0, 3, 11515, 17305, 6745, 6805,
                                                 12145, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18250, 0, 3, 11695, 17620, 6925, 7025,
                                                 12385, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18460, 0, 3, 11785, 17746, 7025, 7125,
                                                 12535, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18670, 0, 3, 11875, 17872, 7125, 7225,
                                                 12685, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18880, 0, 3, 11965, 17998, 7225, 7325,
                                                 12835, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19090, 0, 3, 12055, 18124, 7325, 7425,
                                                 12985, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19300, 0, 3, 17368, 17494, 12235, 18250,
                                                 7625, 7775, 13135, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19615, 0, 3, 17494, 17620, 12385, 18460,
                                                 7775, 7925, 13360, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19930, 0, 3, 17620, 17746, 12535, 18670,
                                                 7925, 8075, 13585, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20245, 0, 3, 17746, 17872, 12685, 18880,
                                                 8075, 8225, 13810, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20560, 0, 3, 17872, 17998, 12835, 19090,
                                                 8225, 8375, 14035, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 20875, 0, 3, 18250, 18460, 13360, 19930,
                                                 8675, 8885, 14575, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 21316, 0, 3, 18460, 18670, 13585, 20245,
                                                 8885, 9095, 14890, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 21757, 0, 3, 18670, 18880, 13810, 20560,
                                                 9095, 9305, 15205, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 22198, 0, 3, 19300, 19615, 14260, 20875,
                                                 9725, 10005, 15520, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 22786, 0, 3, 19615, 19930, 14575, 21316,
                                                 10005, 10285, 15940, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 23374, 0, 3, 19930, 20245, 14890, 21757,
                                                 10285, 10565, 16360, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 23962, 3, 11125, 11140, 16801, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 23990, 3, 11140, 11155, 16822, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 24018, 3, 11155, 11170, 16843, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 24046, 3, 11170, 11185, 16864, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 24074, 3, 11185, 11200, 16885, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 24102, 3, 11200, 11215, 16906, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 24130, 0, 3, 16780, 23962, 16990, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 24214, 0, 3, 16801, 23990, 17053, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 24298, 0, 3, 16822, 24018, 17116, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 24382, 0, 3, 16843, 24046, 17179, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 24466, 0, 3, 16864, 24074, 17242, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 24550, 0, 3, 16885, 24102, 17305, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 24634, 0, 3, 16990, 24214, 11605, 11695,
                                                 17620, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 24802, 0, 3, 17053, 24298, 11695, 11785,
                                                 17746, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 24970, 0, 3, 17116, 24382, 11785, 11875,
                                                 17872, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 25138, 0, 3, 17179, 24466, 11875, 11965,
                                                 17998, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 25306, 0, 3, 17242, 24550, 11965, 12055,
                                                 18124, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25474, 0, 3, 17620, 24802, 12235, 12385,
                                                 18460, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25754, 0, 3, 17746, 24970, 12385, 12535,
                                                 18670, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 26034, 0, 3, 17872, 25138, 12535, 12685,
                                                 18880, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 26314, 0, 3, 17998, 25306, 12685, 12835,
                                                 19090, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 26594, 0, 3, 24634, 24802, 18460, 25754,
                                                 13135, 13360, 19930, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 27014, 0, 3, 24802, 24970, 18670, 26034,
                                                 13360, 13585, 20245, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 27434, 0, 3, 24970, 25138, 18880, 26314,
                                                 13585, 13810, 20560, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 27854, 0, 3, 25474, 25754, 19930, 27014,
                                                 14260, 14575, 21316, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 28442, 0, 3, 25754, 26034, 20245, 27434,
                                                 14575, 14890, 21757, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 29030, 0, 3, 26594, 27014, 21316, 28442,
                                                 15520, 15940, 23374, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 29814, 3, 16780, 16801, 23990, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 29850, 3, 16801, 16822, 24018, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 29886, 3, 16822, 16843, 24046, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 29922, 3, 16843, 16864, 24074, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 29958, 3, 16864, 16885, 24102, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 29994, 0, 3, 23962, 29814, 24214, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 30102, 0, 3, 23990, 29850, 24298, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 30210, 0, 3, 24018, 29886, 24382, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 30318, 0, 3, 24046, 29922, 24466, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 30426, 0, 3, 24074, 29958, 24550, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 30534, 0, 3, 24130, 29994, 17368, 17494,
                                                 24634, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 30750, 0, 3, 24214, 30102, 17494, 17620,
                                                 24802, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 30966, 0, 3, 24298, 30210, 17620, 17746,
                                                 24970, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 31182, 0, 3, 24382, 30318, 17746, 17872,
                                                 25138, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 31398, 0, 3, 24466, 30426, 17872, 17998,
                                                 25306, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 31614, 0, 3, 24802, 30966, 18250, 18460,
                                                 25754, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 31974, 0, 3, 24970, 31182, 18460, 18670,
                                                 26034, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 32334, 0, 3, 25138, 31398, 18670, 18880,
                                                 26314, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 32694, 0, 3, 30534, 30750, 25474, 31614,
                                                 19300, 19615, 26594, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 33234, 0, 3, 30750, 30966, 25754, 31974,
                                                 19615, 19930, 27014, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 33774, 0, 3, 30966, 31182, 26034, 32334,
                                                 19930, 20245, 27434, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 34314, 0, 3, 31614, 31974, 27014, 33774,
                                                 20875, 21316, 28442, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 35070, 0, 3, 32694, 33234, 27854, 34314,
                                                 22198, 22786, 29030, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 36078, 35070, 1008, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 37086, 36078, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 37086, 15, nmax);
}

}  // namespace simdt2ceri

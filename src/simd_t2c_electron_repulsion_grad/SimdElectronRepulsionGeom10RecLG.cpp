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


#include "SimdElectronRepulsionGeom10RecLG.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_lg_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_lg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 32004, 29574, 2025, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 9, 10, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 10, 11, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 11, 12, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 12, 13, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 13, 14, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 14, 15, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 15, 16, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 16, 17, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 119, 0, 17, 18, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 20, 23, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 23, 26, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 26, 29, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 29, 32, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 32, 35, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 35, 38, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 38, 41, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 41, 44, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 44, 47, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 47, 50, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 50, 53, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 235, 0, 59, 65, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 250, 0, 65, 71, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 265, 0, 71, 77, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 280, 0, 77, 83, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 295, 0, 83, 89, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 310, 0, 89, 95, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 325, 0, 95, 101, 205, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 340, 0, 101, 107, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 355, 0, 107, 113, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 370, 0, 125, 135, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 391, 0, 135, 145, 250, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 412, 0, 145, 155, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 433, 0, 155, 165, 280, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 454, 0, 165, 175, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 475, 0, 175, 185, 310, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 496, 0, 185, 195, 325, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 517, 0, 195, 205, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 538, 0, 205, 215, 355, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 559, 0, 235, 250, 412, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 587, 0, 250, 265, 433, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 615, 0, 265, 280, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 643, 0, 280, 295, 475, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 671, 0, 295, 310, 496, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 699, 0, 310, 325, 517, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 727, 0, 325, 340, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 755, 0, 370, 391, 559, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 791, 0, 391, 412, 587, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 827, 0, 412, 433, 615, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 863, 0, 433, 454, 643, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 899, 0, 454, 475, 671, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 935, 0, 475, 496, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 971, 0, 496, 517, 727, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1007, 0, 559, 587, 827, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1052, 0, 587, 615, 863, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1097, 0, 615, 643, 899, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1142, 0, 643, 671, 935, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1187, 0, 671, 699, 971, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1232, 0, 755, 791, 1007, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1287, 0, 791, 827, 1052, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1342, 0, 827, 863, 1097, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1397, 0, 863, 899, 1142, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1452, 0, 899, 935, 1187, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1507, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1510, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1513, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1516, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1519, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1522, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1525, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1528, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1531, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1534, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1537, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1546, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1555, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1564, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1573, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1582, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1591, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1600, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1609, 3, 17, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1618, 3, 18, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1627, 0, 3, 26, 1537, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1645, 0, 3, 29, 1546, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1663, 0, 3, 32, 1555, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1681, 0, 3, 35, 1564, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1699, 0, 3, 38, 1573, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1717, 0, 3, 41, 1582, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1735, 0, 3, 44, 1591, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1753, 0, 3, 47, 1600, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1771, 0, 3, 50, 1609, 113, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1789, 0, 3, 53, 1618, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1807, 0, 3, 65, 1645, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1837, 0, 3, 71, 1663, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1867, 0, 3, 77, 1681, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1897, 0, 3, 83, 1699, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1927, 0, 3, 89, 1717, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1957, 0, 3, 95, 1735, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1987, 0, 3, 101, 1753, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2017, 0, 3, 107, 1771, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2047, 0, 3, 113, 1789, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2077, 0, 3, 145, 1837, 250, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2122, 0, 3, 155, 1867, 265, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2167, 0, 3, 165, 1897, 280, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2212, 0, 3, 175, 1927, 295, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2257, 0, 3, 185, 1957, 310, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2302, 0, 3, 195, 1987, 325, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2347, 0, 3, 205, 2017, 340, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2392, 0, 3, 215, 2047, 355, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2437, 0, 3, 250, 2122, 412, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2500, 0, 3, 265, 2167, 433, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2563, 0, 3, 280, 2212, 454, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2626, 0, 3, 295, 2257, 475, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2689, 0, 3, 310, 2302, 496, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2752, 0, 3, 325, 2347, 517, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2815, 0, 3, 340, 2392, 538, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2878, 0, 3, 412, 2500, 587, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2962, 0, 3, 433, 2563, 615, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3046, 0, 3, 454, 2626, 643, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3130, 0, 3, 475, 2689, 671, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3214, 0, 3, 496, 2752, 699, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3298, 0, 3, 517, 2815, 727, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3382, 0, 3, 587, 2962, 827, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3490, 0, 3, 615, 3046, 863, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3598, 0, 3, 643, 3130, 899, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3706, 0, 3, 671, 3214, 935, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3814, 0, 3, 699, 3298, 971, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 3922, 0, 3, 827, 3490, 1052, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4057, 0, 3, 863, 3598, 1097, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4192, 0, 3, 899, 3706, 1142, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4327, 0, 3, 935, 3814, 1187, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 4462, 0, 3, 1052, 4057, 1342, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 4627, 0, 3, 1097, 4192, 1397, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 4792, 0, 3, 1142, 4327, 1452, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4957, 3, 9, 10, 1510, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4963, 3, 10, 11, 1513, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4969, 3, 11, 12, 1516, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4975, 3, 12, 13, 1519, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4981, 3, 13, 14, 1522, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4987, 3, 14, 15, 1525, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4993, 3, 15, 16, 1528, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4999, 3, 16, 17, 1531, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5005, 3, 17, 18, 1534, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 5011, 0, 3, 1507, 4957, 1546, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5029, 0, 3, 1510, 4963, 1555, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5047, 0, 3, 1513, 4969, 1564, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5065, 0, 3, 1516, 4975, 1573, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5083, 0, 3, 1519, 4981, 1582, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5101, 0, 3, 1522, 4987, 1591, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5119, 0, 3, 1525, 4993, 1600, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5137, 0, 3, 1528, 4999, 1609, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5155, 0, 3, 1531, 5005, 1618, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 5173, 0, 3, 1537, 5011, 59, 65, 1645,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5209, 0, 3, 1546, 5029, 65, 71, 1663,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5245, 0, 3, 1555, 5047, 71, 77, 1681,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5281, 0, 3, 1564, 5065, 77, 83, 1699,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5317, 0, 3, 1573, 5083, 83, 89, 1717,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5353, 0, 3, 1582, 5101, 89, 95, 1735,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5389, 0, 3, 1591, 5119, 95, 101, 1753,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5425, 0, 3, 1600, 5137, 101, 107, 1771,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5461, 0, 3, 1609, 5155, 107, 113, 1789,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5497, 0, 3, 1627, 5173, 125, 135, 1807,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5557, 0, 3, 1645, 5209, 135, 145, 1837,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5617, 0, 3, 1663, 5245, 145, 155, 1867,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5677, 0, 3, 1681, 5281, 155, 165, 1897,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5737, 0, 3, 1699, 5317, 165, 175, 1927,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5797, 0, 3, 1717, 5353, 175, 185, 1957,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5857, 0, 3, 1735, 5389, 185, 195, 1987,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5917, 0, 3, 1753, 5425, 195, 205, 2017,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5977, 0, 3, 1771, 5461, 205, 215, 2047,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6037, 0, 3, 5173, 5209, 1837, 5617, 235,
                                                 250, 2122, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6127, 0, 3, 5209, 5245, 1867, 5677, 250,
                                                 265, 2167, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6217, 0, 3, 5245, 5281, 1897, 5737, 265,
                                                 280, 2212, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6307, 0, 3, 5281, 5317, 1927, 5797, 280,
                                                 295, 2257, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6397, 0, 3, 5317, 5353, 1957, 5857, 295,
                                                 310, 2302, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6487, 0, 3, 5353, 5389, 1987, 5917, 310,
                                                 325, 2347, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6577, 0, 3, 5389, 5425, 2017, 5977, 325,
                                                 340, 2392, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6667, 0, 3, 5497, 5557, 2077, 6037, 370,
                                                 391, 2437, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6793, 0, 3, 5557, 5617, 2122, 6127, 391,
                                                 412, 2500, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6919, 0, 3, 5617, 5677, 2167, 6217, 412,
                                                 433, 2563, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7045, 0, 3, 5677, 5737, 2212, 6307, 433,
                                                 454, 2626, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7171, 0, 3, 5737, 5797, 2257, 6397, 454,
                                                 475, 2689, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7297, 0, 3, 5797, 5857, 2302, 6487, 475,
                                                 496, 2752, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7423, 0, 3, 5857, 5917, 2347, 6577, 496,
                                                 517, 2815, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7549, 0, 3, 6037, 6127, 2500, 6919, 559,
                                                 587, 2962, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7717, 0, 3, 6127, 6217, 2563, 7045, 587,
                                                 615, 3046, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7885, 0, 3, 6217, 6307, 2626, 7171, 615,
                                                 643, 3130, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8053, 0, 3, 6307, 6397, 2689, 7297, 643,
                                                 671, 3214, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8221, 0, 3, 6397, 6487, 2752, 7423, 671,
                                                 699, 3298, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8389, 0, 3, 6667, 6793, 2878, 7549, 755,
                                                 791, 3382, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8605, 0, 3, 6793, 6919, 2962, 7717, 791,
                                                 827, 3490, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8821, 0, 3, 6919, 7045, 3046, 7885, 827,
                                                 863, 3598, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9037, 0, 3, 7045, 7171, 3130, 8053, 863,
                                                 899, 3706, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9253, 0, 3, 7171, 7297, 3214, 8221, 899,
                                                 935, 3814, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 9469, 0, 3, 7549, 7717, 3490, 8821,
                                                 1007, 1052, 4057, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 9739, 0, 3, 7717, 7885, 3598, 9037,
                                                 1052, 1097, 4192, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10009, 0, 3, 7885, 8053, 3706, 9253,
                                                 1097, 1142, 4327, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 10279, 0, 3, 8389, 8605, 3922, 9469,
                                                 1232, 1287, 4462, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 10609, 0, 3, 8605, 8821, 4057, 9739,
                                                 1287, 1342, 4627, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 10939, 0, 3, 8821, 9037, 4192, 10009,
                                                 1342, 1397, 4792, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11269, 3, 1507, 1510, 4963, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11279, 3, 1510, 1513, 4969, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11289, 3, 1513, 1516, 4975, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11299, 3, 1516, 1519, 4981, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11309, 3, 1519, 1522, 4987, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11319, 3, 1522, 1525, 4993, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11329, 3, 1525, 1528, 4999, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11339, 3, 1528, 1531, 5005, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 11349, 0, 3, 4957, 11269, 5029, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11379, 0, 3, 4963, 11279, 5047, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11409, 0, 3, 4969, 11289, 5065, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11439, 0, 3, 4975, 11299, 5083, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11469, 0, 3, 4981, 11309, 5101, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11499, 0, 3, 4987, 11319, 5119, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11529, 0, 3, 4993, 11329, 5137, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11559, 0, 3, 4999, 11339, 5155, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 11589, 0, 3, 5011, 11349, 1627, 1645,
                                                 5209, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11649, 0, 3, 5029, 11379, 1645, 1663,
                                                 5245, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11709, 0, 3, 5047, 11409, 1663, 1681,
                                                 5281, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11769, 0, 3, 5065, 11439, 1681, 1699,
                                                 5317, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11829, 0, 3, 5083, 11469, 1699, 1717,
                                                 5353, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11889, 0, 3, 5101, 11499, 1717, 1735,
                                                 5389, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11949, 0, 3, 5119, 11529, 1735, 1753,
                                                 5425, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12009, 0, 3, 5137, 11559, 1753, 1771,
                                                 5461, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12069, 0, 3, 5209, 11649, 1807, 1837,
                                                 5617, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12169, 0, 3, 5245, 11709, 1837, 1867,
                                                 5677, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12269, 0, 3, 5281, 11769, 1867, 1897,
                                                 5737, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12369, 0, 3, 5317, 11829, 1897, 1927,
                                                 5797, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12469, 0, 3, 5353, 11889, 1927, 1957,
                                                 5857, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12569, 0, 3, 5389, 11949, 1957, 1987,
                                                 5917, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12669, 0, 3, 5425, 12009, 1987, 2017,
                                                 5977, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12769, 0, 3, 11589, 11649, 5617, 12169,
                                                 2077, 2122, 6127, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12919, 0, 3, 11649, 11709, 5677, 12269,
                                                 2122, 2167, 6217, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13069, 0, 3, 11709, 11769, 5737, 12369,
                                                 2167, 2212, 6307, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13219, 0, 3, 11769, 11829, 5797, 12469,
                                                 2212, 2257, 6397, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13369, 0, 3, 11829, 11889, 5857, 12569,
                                                 2257, 2302, 6487, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13519, 0, 3, 11889, 11949, 5917, 12669,
                                                 2302, 2347, 6577, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13669, 0, 3, 12069, 12169, 6127, 12919,
                                                 2437, 2500, 6919, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13879, 0, 3, 12169, 12269, 6217, 13069,
                                                 2500, 2563, 7045, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14089, 0, 3, 12269, 12369, 6307, 13219,
                                                 2563, 2626, 7171, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14299, 0, 3, 12369, 12469, 6397, 13369,
                                                 2626, 2689, 7297, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14509, 0, 3, 12469, 12569, 6487, 13519,
                                                 2689, 2752, 7423, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14719, 0, 3, 12769, 12919, 6919, 13879,
                                                 2878, 2962, 7717, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 14999, 0, 3, 12919, 13069, 7045, 14089,
                                                 2962, 3046, 7885, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15279, 0, 3, 13069, 13219, 7171, 14299,
                                                 3046, 3130, 8053, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15559, 0, 3, 13219, 13369, 7297, 14509,
                                                 3130, 3214, 8221, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 15839, 0, 3, 13669, 13879, 7717, 14999,
                                                 3382, 3490, 8821, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16199, 0, 3, 13879, 14089, 7885, 15279,
                                                 3490, 3598, 9037, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16559, 0, 3, 14089, 14299, 8053, 15559,
                                                 3598, 3706, 9253, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 16919, 0, 3, 14719, 14999, 8821, 16199,
                                                 3922, 4057, 9739, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 17369, 0, 3, 14999, 15279, 9037, 16559,
                                                 4057, 4192, 10009, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 17819, 0, 3, 15839, 16199, 9739, 17369,
                                                 4462, 4627, 10939, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18369, 3, 4957, 4963, 11279, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18384, 3, 4963, 4969, 11289, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18399, 3, 4969, 4975, 11299, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18414, 3, 4975, 4981, 11309, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18429, 3, 4981, 4987, 11319, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18444, 3, 4987, 4993, 11329, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 18459, 3, 4993, 4999, 11339, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 18474, 0, 3, 11269, 18369, 11379, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18519, 0, 3, 11279, 18384, 11409, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18564, 0, 3, 11289, 18399, 11439, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18609, 0, 3, 11299, 18414, 11469, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18654, 0, 3, 11309, 18429, 11499, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18699, 0, 3, 11319, 18444, 11529, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18744, 0, 3, 11329, 18459, 11559, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 18789, 0, 3, 11349, 18474, 5173, 5209,
                                                 11649, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18879, 0, 3, 11379, 18519, 5209, 5245,
                                                 11709, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18969, 0, 3, 11409, 18564, 5245, 5281,
                                                 11769, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19059, 0, 3, 11439, 18609, 5281, 5317,
                                                 11829, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19149, 0, 3, 11469, 18654, 5317, 5353,
                                                 11889, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19239, 0, 3, 11499, 18699, 5353, 5389,
                                                 11949, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 19329, 0, 3, 11529, 18744, 5389, 5425,
                                                 12009, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19419, 0, 3, 11589, 18789, 5497, 5557,
                                                 12069, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19569, 0, 3, 11649, 18879, 5557, 5617,
                                                 12169, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19719, 0, 3, 11709, 18969, 5617, 5677,
                                                 12269, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19869, 0, 3, 11769, 19059, 5677, 5737,
                                                 12369, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20019, 0, 3, 11829, 19149, 5737, 5797,
                                                 12469, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20169, 0, 3, 11889, 19239, 5797, 5857,
                                                 12569, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20319, 0, 3, 11949, 19329, 5857, 5917,
                                                 12669, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20469, 0, 3, 18789, 18879, 12169, 19719,
                                                 6037, 6127, 12919, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20694, 0, 3, 18879, 18969, 12269, 19869,
                                                 6127, 6217, 13069, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20919, 0, 3, 18969, 19059, 12369, 20019,
                                                 6217, 6307, 13219, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21144, 0, 3, 19059, 19149, 12469, 20169,
                                                 6307, 6397, 13369, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21369, 0, 3, 19149, 19239, 12569, 20319,
                                                 6397, 6487, 13519, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 21594, 0, 3, 19419, 19569, 12769, 20469,
                                                 6667, 6793, 13669, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 21909, 0, 3, 19569, 19719, 12919, 20694,
                                                 6793, 6919, 13879, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22224, 0, 3, 19719, 19869, 13069, 20919,
                                                 6919, 7045, 14089, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22539, 0, 3, 19869, 20019, 13219, 21144,
                                                 7045, 7171, 14299, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22854, 0, 3, 20019, 20169, 13369, 21369,
                                                 7171, 7297, 14509, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 23169, 0, 3, 20469, 20694, 13879, 22224,
                                                 7549, 7717, 14999, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 23589, 0, 3, 20694, 20919, 14089, 22539,
                                                 7717, 7885, 15279, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 24009, 0, 3, 20919, 21144, 14299, 22854,
                                                 7885, 8053, 15559, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 24429, 0, 3, 21594, 21909, 14719, 23169,
                                                 8389, 8605, 15839, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 24969, 0, 3, 21909, 22224, 14999, 23589,
                                                 8605, 8821, 16199, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 25509, 0, 3, 22224, 22539, 15279, 24009,
                                                 8821, 9037, 16559, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 26049, 0, 3, 23169, 23589, 16199, 25509,
                                                 9469, 9739, 17369, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 26724, 0, 3, 24429, 24969, 16919, 26049,
                                                 10279, 10609, 17819, ncols, alpha, beta, p);

            compute_prim_geom_10_lg_electron_repulsion_0(buffer, 27549, 24429, 26724, ncols,
                                                         alpha);

            compute_prim_geom_10_lg_electron_repulsion_1(buffer, 28224, 24429, 26724, ncols,
                                                         alpha);

            compute_prim_geom_10_lg_electron_repulsion_2(buffer, 28899, 24429, 26724, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 29574, 27549, 2025, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 31599, 29574, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 31599, 9, nmax);

    simdtrf::transform_g_inner(buffer, 31599, 30249, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 153 * nvalues, nvalues, buffer, 31599, 9, nmax);

    simdtrf::transform_g_inner(buffer, 31599, 30924, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 306 * nvalues, nvalues, buffer, 31599, 9, nmax);
}

}  // namespace simdt2ceri

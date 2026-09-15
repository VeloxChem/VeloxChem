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


#include "SimdElectronRepulsionGeom10RecLH.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_lh_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_lh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 53165, 49835, 2835, nvalues);

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
                                            10, 11, 12, 13, 14}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 20, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 16, 17, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 17, 18, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 18, 19, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 21, 24, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 39, 42, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 42, 45, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 45, 48, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 48, 51, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 51, 54, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 242, 0, 60, 66, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 257, 0, 66, 72, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 272, 0, 72, 78, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 287, 0, 78, 84, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 302, 0, 84, 90, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 317, 0, 90, 96, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 332, 0, 96, 102, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 347, 0, 102, 108, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 362, 0, 108, 114, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 377, 0, 114, 120, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 392, 0, 132, 142, 257, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 413, 0, 142, 152, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 434, 0, 152, 162, 287, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 455, 0, 162, 172, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 476, 0, 172, 182, 317, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 497, 0, 182, 192, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 518, 0, 192, 202, 347, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 539, 0, 202, 212, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 560, 0, 212, 222, 377, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 581, 0, 242, 257, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 609, 0, 257, 272, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 272, 287, 455, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 287, 302, 476, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 302, 317, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 721, 0, 317, 332, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 749, 0, 332, 347, 539, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 777, 0, 347, 362, 560, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 805, 0, 392, 413, 609, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 841, 0, 413, 434, 637, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 877, 0, 434, 455, 665, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 913, 0, 455, 476, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 949, 0, 476, 497, 721, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 985, 0, 497, 518, 749, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1021, 0, 518, 539, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1057, 0, 581, 609, 841, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1102, 0, 609, 637, 877, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1147, 0, 637, 665, 913, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1192, 0, 665, 693, 949, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1237, 0, 693, 721, 985, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1282, 0, 721, 749, 1021, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1327, 0, 805, 841, 1102, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1382, 0, 841, 877, 1147, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1437, 0, 877, 913, 1192, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1492, 0, 913, 949, 1237, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1547, 0, 949, 985, 1282, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1602, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1605, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1608, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1611, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1614, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1617, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1620, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1623, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1626, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1629, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1632, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1635, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1644, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1653, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1662, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1671, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1680, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1689, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1698, 3, 16, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1707, 3, 17, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1716, 3, 18, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1725, 3, 19, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1734, 0, 3, 24, 1635, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1752, 0, 3, 27, 1644, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1770, 0, 3, 30, 1653, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1788, 0, 3, 33, 1662, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1806, 0, 3, 36, 1671, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1824, 0, 3, 39, 1680, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1842, 0, 3, 42, 1689, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1860, 0, 3, 45, 1698, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1878, 0, 3, 48, 1707, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1896, 0, 3, 51, 1716, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1914, 0, 3, 54, 1725, 126, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1932, 0, 3, 60, 1734, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1962, 0, 3, 66, 1752, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1992, 0, 3, 72, 1770, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2022, 0, 3, 78, 1788, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2052, 0, 3, 84, 1806, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2082, 0, 3, 90, 1824, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2112, 0, 3, 96, 1842, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2142, 0, 3, 102, 1860, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2172, 0, 3, 108, 1878, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2202, 0, 3, 114, 1896, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2232, 0, 3, 120, 1914, 232, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2262, 0, 3, 142, 1992, 257, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2307, 0, 3, 152, 2022, 272, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2352, 0, 3, 162, 2052, 287, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2397, 0, 3, 172, 2082, 302, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2442, 0, 3, 182, 2112, 317, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2487, 0, 3, 192, 2142, 332, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2532, 0, 3, 202, 2172, 347, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2577, 0, 3, 212, 2202, 362, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2622, 0, 3, 222, 2232, 377, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2667, 0, 3, 242, 2262, 392, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2730, 0, 3, 257, 2307, 413, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2793, 0, 3, 272, 2352, 434, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2856, 0, 3, 287, 2397, 455, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2919, 0, 3, 302, 2442, 476, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2982, 0, 3, 317, 2487, 497, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3045, 0, 3, 332, 2532, 518, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3108, 0, 3, 347, 2577, 539, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3171, 0, 3, 362, 2622, 560, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3234, 0, 3, 413, 2793, 609, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3318, 0, 3, 434, 2856, 637, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3402, 0, 3, 455, 2919, 665, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3486, 0, 3, 476, 2982, 693, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3570, 0, 3, 497, 3045, 721, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3654, 0, 3, 518, 3108, 749, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3738, 0, 3, 539, 3171, 777, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3822, 0, 3, 581, 3234, 805, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3930, 0, 3, 609, 3318, 841, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4038, 0, 3, 637, 3402, 877, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4146, 0, 3, 665, 3486, 913, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4254, 0, 3, 693, 3570, 949, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4362, 0, 3, 721, 3654, 985, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4470, 0, 3, 749, 3738, 1021, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4578, 0, 3, 841, 4038, 1102, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4713, 0, 3, 877, 4146, 1147, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4848, 0, 3, 913, 4254, 1192, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4983, 0, 3, 949, 4362, 1237, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5118, 0, 3, 985, 4470, 1282, ncols, p);

            compute_prim_mp_electron_repulsion_0(buffer, 5253, 0, 3, 1057, 4578, 1327, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 5418, 0, 3, 1102, 4713, 1382, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 5583, 0, 3, 1147, 4848, 1437, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 5748, 0, 3, 1192, 4983, 1492, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 5913, 0, 3, 1237, 5118, 1547, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6078, 3, 9, 10, 1605, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6084, 3, 10, 11, 1608, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6090, 3, 11, 12, 1611, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6096, 3, 12, 13, 1614, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6102, 3, 13, 14, 1617, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6108, 3, 14, 15, 1620, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6114, 3, 15, 16, 1623, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6120, 3, 16, 17, 1626, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6126, 3, 17, 18, 1629, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6132, 3, 18, 19, 1632, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6138, 0, 3, 1602, 6078, 1644, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6156, 0, 3, 1605, 6084, 1653, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6174, 0, 3, 1608, 6090, 1662, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6192, 0, 3, 1611, 6096, 1671, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6210, 0, 3, 1614, 6102, 1680, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6228, 0, 3, 1617, 6108, 1689, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6246, 0, 3, 1620, 6114, 1698, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6264, 0, 3, 1623, 6120, 1707, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6282, 0, 3, 1626, 6126, 1716, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6300, 0, 3, 1629, 6132, 1725, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6318, 0, 3, 1635, 6138, 60, 66, 1752,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6354, 0, 3, 1644, 6156, 66, 72, 1770,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6390, 0, 3, 1653, 6174, 72, 78, 1788,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6426, 0, 3, 1662, 6192, 78, 84, 1806,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6462, 0, 3, 1671, 6210, 84, 90, 1824,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6498, 0, 3, 1680, 6228, 90, 96, 1842,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6534, 0, 3, 1689, 6246, 96, 102, 1860,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6570, 0, 3, 1698, 6264, 102, 108, 1878,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6606, 0, 3, 1707, 6282, 108, 114, 1896,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6642, 0, 3, 1716, 6300, 114, 120, 1914,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6678, 0, 3, 1752, 6354, 132, 142, 1992,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6738, 0, 3, 1770, 6390, 142, 152, 2022,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6798, 0, 3, 1788, 6426, 152, 162, 2052,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6858, 0, 3, 1806, 6462, 162, 172, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6918, 0, 3, 1824, 6498, 172, 182, 2112,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6978, 0, 3, 1842, 6534, 182, 192, 2142,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7038, 0, 3, 1860, 6570, 192, 202, 2172,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7098, 0, 3, 1878, 6606, 202, 212, 2202,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7158, 0, 3, 1896, 6642, 212, 222, 2232,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7218, 0, 3, 6318, 6354, 1992, 6738, 242,
                                                 257, 2307, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7308, 0, 3, 6354, 6390, 2022, 6798, 257,
                                                 272, 2352, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7398, 0, 3, 6390, 6426, 2052, 6858, 272,
                                                 287, 2397, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7488, 0, 3, 6426, 6462, 2082, 6918, 287,
                                                 302, 2442, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7578, 0, 3, 6462, 6498, 2112, 6978, 302,
                                                 317, 2487, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7668, 0, 3, 6498, 6534, 2142, 7038, 317,
                                                 332, 2532, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7758, 0, 3, 6534, 6570, 2172, 7098, 332,
                                                 347, 2577, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7848, 0, 3, 6570, 6606, 2202, 7158, 347,
                                                 362, 2622, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7938, 0, 3, 6678, 6738, 2307, 7308, 392,
                                                 413, 2793, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8064, 0, 3, 6738, 6798, 2352, 7398, 413,
                                                 434, 2856, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8190, 0, 3, 6798, 6858, 2397, 7488, 434,
                                                 455, 2919, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8316, 0, 3, 6858, 6918, 2442, 7578, 455,
                                                 476, 2982, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8442, 0, 3, 6918, 6978, 2487, 7668, 476,
                                                 497, 3045, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8568, 0, 3, 6978, 7038, 2532, 7758, 497,
                                                 518, 3108, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8694, 0, 3, 7038, 7098, 2577, 7848, 518,
                                                 539, 3171, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8820, 0, 3, 7218, 7308, 2793, 8064, 581,
                                                 609, 3318, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8988, 0, 3, 7308, 7398, 2856, 8190, 609,
                                                 637, 3402, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9156, 0, 3, 7398, 7488, 2919, 8316, 637,
                                                 665, 3486, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9324, 0, 3, 7488, 7578, 2982, 8442, 665,
                                                 693, 3570, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9492, 0, 3, 7578, 7668, 3045, 8568, 693,
                                                 721, 3654, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9660, 0, 3, 7668, 7758, 3108, 8694, 721,
                                                 749, 3738, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9828, 0, 3, 7938, 8064, 3318, 8988, 805,
                                                 841, 4038, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10044, 0, 3, 8064, 8190, 3402, 9156,
                                                 841, 877, 4146, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10260, 0, 3, 8190, 8316, 3486, 9324,
                                                 877, 913, 4254, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10476, 0, 3, 8316, 8442, 3570, 9492,
                                                 913, 949, 4362, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10692, 0, 3, 8442, 8568, 3654, 9660,
                                                 949, 985, 4470, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10908, 0, 3, 8820, 8988, 4038, 10044,
                                                 1057, 1102, 4713, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11178, 0, 3, 8988, 9156, 4146, 10260,
                                                 1102, 1147, 4848, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11448, 0, 3, 9156, 9324, 4254, 10476,
                                                 1147, 1192, 4983, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11718, 0, 3, 9324, 9492, 4362, 10692,
                                                 1192, 1237, 5118, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 11988, 0, 3, 9828, 10044, 4713, 11178,
                                                 1327, 1382, 5583, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 12318, 0, 3, 10044, 10260, 4848, 11448,
                                                 1382, 1437, 5748, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 12648, 0, 3, 10260, 10476, 4983, 11718,
                                                 1437, 1492, 5913, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12978, 3, 1602, 1605, 6084, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12988, 3, 1605, 1608, 6090, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12998, 3, 1608, 1611, 6096, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13008, 3, 1611, 1614, 6102, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13018, 3, 1614, 1617, 6108, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13028, 3, 1617, 1620, 6114, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13038, 3, 1620, 1623, 6120, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13048, 3, 1623, 1626, 6126, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 13058, 3, 1626, 1629, 6132, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 13068, 0, 3, 6078, 12978, 6156, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13098, 0, 3, 6084, 12988, 6174, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13128, 0, 3, 6090, 12998, 6192, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13158, 0, 3, 6096, 13008, 6210, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13188, 0, 3, 6102, 13018, 6228, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13218, 0, 3, 6108, 13028, 6246, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13248, 0, 3, 6114, 13038, 6264, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13278, 0, 3, 6120, 13048, 6282, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13308, 0, 3, 6126, 13058, 6300, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 13338, 0, 3, 6138, 13068, 1734, 1752,
                                                 6354, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13398, 0, 3, 6156, 13098, 1752, 1770,
                                                 6390, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13458, 0, 3, 6174, 13128, 1770, 1788,
                                                 6426, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13518, 0, 3, 6192, 13158, 1788, 1806,
                                                 6462, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13578, 0, 3, 6210, 13188, 1806, 1824,
                                                 6498, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13638, 0, 3, 6228, 13218, 1824, 1842,
                                                 6534, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13698, 0, 3, 6246, 13248, 1842, 1860,
                                                 6570, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13758, 0, 3, 6264, 13278, 1860, 1878,
                                                 6606, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13818, 0, 3, 6282, 13308, 1878, 1896,
                                                 6642, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13878, 0, 3, 6318, 13338, 1932, 1962,
                                                 6678, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13978, 0, 3, 6354, 13398, 1962, 1992,
                                                 6738, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14078, 0, 3, 6390, 13458, 1992, 2022,
                                                 6798, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14178, 0, 3, 6426, 13518, 2022, 2052,
                                                 6858, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14278, 0, 3, 6462, 13578, 2052, 2082,
                                                 6918, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14378, 0, 3, 6498, 13638, 2082, 2112,
                                                 6978, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14478, 0, 3, 6534, 13698, 2112, 2142,
                                                 7038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14578, 0, 3, 6570, 13758, 2142, 2172,
                                                 7098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14678, 0, 3, 6606, 13818, 2172, 2202,
                                                 7158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14778, 0, 3, 13338, 13398, 6738, 14078,
                                                 2262, 2307, 7308, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14928, 0, 3, 13398, 13458, 6798, 14178,
                                                 2307, 2352, 7398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15078, 0, 3, 13458, 13518, 6858, 14278,
                                                 2352, 2397, 7488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15228, 0, 3, 13518, 13578, 6918, 14378,
                                                 2397, 2442, 7578, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15378, 0, 3, 13578, 13638, 6978, 14478,
                                                 2442, 2487, 7668, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15528, 0, 3, 13638, 13698, 7038, 14578,
                                                 2487, 2532, 7758, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15678, 0, 3, 13698, 13758, 7098, 14678,
                                                 2532, 2577, 7848, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15828, 0, 3, 13878, 13978, 7218, 14778,
                                                 2667, 2730, 7938, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16038, 0, 3, 13978, 14078, 7308, 14928,
                                                 2730, 2793, 8064, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16248, 0, 3, 14078, 14178, 7398, 15078,
                                                 2793, 2856, 8190, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16458, 0, 3, 14178, 14278, 7488, 15228,
                                                 2856, 2919, 8316, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16668, 0, 3, 14278, 14378, 7578, 15378,
                                                 2919, 2982, 8442, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16878, 0, 3, 14378, 14478, 7668, 15528,
                                                 2982, 3045, 8568, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17088, 0, 3, 14478, 14578, 7758, 15678,
                                                 3045, 3108, 8694, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17298, 0, 3, 14778, 14928, 8064, 16248,
                                                 3234, 3318, 8988, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17578, 0, 3, 14928, 15078, 8190, 16458,
                                                 3318, 3402, 9156, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17858, 0, 3, 15078, 15228, 8316, 16668,
                                                 3402, 3486, 9324, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18138, 0, 3, 15228, 15378, 8442, 16878,
                                                 3486, 3570, 9492, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18418, 0, 3, 15378, 15528, 8568, 17088,
                                                 3570, 3654, 9660, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 18698, 0, 3, 15828, 16038, 8820, 17298,
                                                 3822, 3930, 9828, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19058, 0, 3, 16038, 16248, 8988, 17578,
                                                 3930, 4038, 10044, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19418, 0, 3, 16248, 16458, 9156, 17858,
                                                 4038, 4146, 10260, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19778, 0, 3, 16458, 16668, 9324, 18138,
                                                 4146, 4254, 10476, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 20138, 0, 3, 16668, 16878, 9492, 18418,
                                                 4254, 4362, 10692, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 20498, 0, 3, 17298, 17578, 10044, 19418,
                                                 4578, 4713, 11178, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 20948, 0, 3, 17578, 17858, 10260, 19778,
                                                 4713, 4848, 11448, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 21398, 0, 3, 17858, 18138, 10476, 20138,
                                                 4848, 4983, 11718, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 21848, 0, 3, 18698, 19058, 10908, 20498,
                                                 5253, 5418, 11988, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 22398, 0, 3, 19058, 19418, 11178, 20948,
                                                 5418, 5583, 12318, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 22948, 0, 3, 19418, 19778, 11448, 21398,
                                                 5583, 5748, 12648, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23498, 3, 6078, 6084, 12988, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23513, 3, 6084, 6090, 12998, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23528, 3, 6090, 6096, 13008, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23543, 3, 6096, 6102, 13018, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23558, 3, 6102, 6108, 13028, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23573, 3, 6108, 6114, 13038, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23588, 3, 6114, 6120, 13048, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23603, 3, 6120, 6126, 13058, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 23618, 0, 3, 12978, 23498, 13098, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23663, 0, 3, 12988, 23513, 13128, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23708, 0, 3, 12998, 23528, 13158, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23753, 0, 3, 13008, 23543, 13188, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23798, 0, 3, 13018, 23558, 13218, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23843, 0, 3, 13028, 23573, 13248, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23888, 0, 3, 13038, 23588, 13278, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23933, 0, 3, 13048, 23603, 13308, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 23978, 0, 3, 13068, 23618, 6318, 6354,
                                                 13398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24068, 0, 3, 13098, 23663, 6354, 6390,
                                                 13458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24158, 0, 3, 13128, 23708, 6390, 6426,
                                                 13518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24248, 0, 3, 13158, 23753, 6426, 6462,
                                                 13578, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24338, 0, 3, 13188, 23798, 6462, 6498,
                                                 13638, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24428, 0, 3, 13218, 23843, 6498, 6534,
                                                 13698, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24518, 0, 3, 13248, 23888, 6534, 6570,
                                                 13758, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24608, 0, 3, 13278, 23933, 6570, 6606,
                                                 13818, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24698, 0, 3, 13398, 24068, 6678, 6738,
                                                 14078, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24848, 0, 3, 13458, 24158, 6738, 6798,
                                                 14178, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24998, 0, 3, 13518, 24248, 6798, 6858,
                                                 14278, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25148, 0, 3, 13578, 24338, 6858, 6918,
                                                 14378, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25298, 0, 3, 13638, 24428, 6918, 6978,
                                                 14478, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25448, 0, 3, 13698, 24518, 6978, 7038,
                                                 14578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25598, 0, 3, 13758, 24608, 7038, 7098,
                                                 14678, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 25748, 0, 3, 23978, 24068, 14078, 24848,
                                                 7218, 7308, 14928, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 25973, 0, 3, 24068, 24158, 14178, 24998,
                                                 7308, 7398, 15078, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26198, 0, 3, 24158, 24248, 14278, 25148,
                                                 7398, 7488, 15228, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26423, 0, 3, 24248, 24338, 14378, 25298,
                                                 7488, 7578, 15378, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26648, 0, 3, 24338, 24428, 14478, 25448,
                                                 7578, 7668, 15528, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26873, 0, 3, 24428, 24518, 14578, 25598,
                                                 7668, 7758, 15678, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 27098, 0, 3, 24698, 24848, 14928, 25973,
                                                 7938, 8064, 16248, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 27413, 0, 3, 24848, 24998, 15078, 26198,
                                                 8064, 8190, 16458, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 27728, 0, 3, 24998, 25148, 15228, 26423,
                                                 8190, 8316, 16668, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28043, 0, 3, 25148, 25298, 15378, 26648,
                                                 8316, 8442, 16878, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28358, 0, 3, 25298, 25448, 15528, 26873,
                                                 8442, 8568, 17088, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 28673, 0, 3, 25748, 25973, 16248, 27413,
                                                 8820, 8988, 17578, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 29093, 0, 3, 25973, 26198, 16458, 27728,
                                                 8988, 9156, 17858, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 29513, 0, 3, 26198, 26423, 16668, 28043,
                                                 9156, 9324, 18138, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 29933, 0, 3, 26423, 26648, 16878, 28358,
                                                 9324, 9492, 18418, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 30353, 0, 3, 27098, 27413, 17578, 29093,
                                                 9828, 10044, 19418, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 30893, 0, 3, 27413, 27728, 17858, 29513,
                                                 10044, 10260, 19778, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 31433, 0, 3, 27728, 28043, 18138, 29933,
                                                 10260, 10476, 20138, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 31973, 0, 3, 28673, 29093, 19418, 30893,
                                                 10908, 11178, 20948, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 32648, 0, 3, 29093, 29513, 19778, 31433,
                                                 11178, 11448, 21398, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 33323, 0, 3, 30353, 30893, 20948, 32648,
                                                 11988, 12318, 22948, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34148, 3, 12978, 12988, 23513, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34169, 3, 12988, 12998, 23528, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34190, 3, 12998, 13008, 23543, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34211, 3, 13008, 13018, 23558, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34232, 3, 13018, 13028, 23573, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34253, 3, 13028, 13038, 23588, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 34274, 3, 13038, 13048, 23603, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 34295, 0, 3, 23498, 34148, 23663, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34358, 0, 3, 23513, 34169, 23708, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34421, 0, 3, 23528, 34190, 23753, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34484, 0, 3, 23543, 34211, 23798, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34547, 0, 3, 23558, 34232, 23843, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34610, 0, 3, 23573, 34253, 23888, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 34673, 0, 3, 23588, 34274, 23933, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 34736, 0, 3, 23618, 34295, 13338, 13398,
                                                 24068, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34862, 0, 3, 23663, 34358, 13398, 13458,
                                                 24158, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 34988, 0, 3, 23708, 34421, 13458, 13518,
                                                 24248, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35114, 0, 3, 23753, 34484, 13518, 13578,
                                                 24338, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35240, 0, 3, 23798, 34547, 13578, 13638,
                                                 24428, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35366, 0, 3, 23843, 34610, 13638, 13698,
                                                 24518, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 35492, 0, 3, 23888, 34673, 13698, 13758,
                                                 24608, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35618, 0, 3, 23978, 34736, 13878, 13978,
                                                 24698, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 35828, 0, 3, 24068, 34862, 13978, 14078,
                                                 24848, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36038, 0, 3, 24158, 34988, 14078, 14178,
                                                 24998, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36248, 0, 3, 24248, 35114, 14178, 14278,
                                                 25148, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36458, 0, 3, 24338, 35240, 14278, 14378,
                                                 25298, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36668, 0, 3, 24428, 35366, 14378, 14478,
                                                 25448, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 36878, 0, 3, 24518, 35492, 14478, 14578,
                                                 25598, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37088, 0, 3, 34736, 34862, 24848, 36038,
                                                 14778, 14928, 25973, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37403, 0, 3, 34862, 34988, 24998, 36248,
                                                 14928, 15078, 26198, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 37718, 0, 3, 34988, 35114, 25148, 36458,
                                                 15078, 15228, 26423, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38033, 0, 3, 35114, 35240, 25298, 36668,
                                                 15228, 15378, 26648, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 38348, 0, 3, 35240, 35366, 25448, 36878,
                                                 15378, 15528, 26873, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 38663, 0, 3, 35618, 35828, 25748, 37088,
                                                 15828, 16038, 27098, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39104, 0, 3, 35828, 36038, 25973, 37403,
                                                 16038, 16248, 27413, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39545, 0, 3, 36038, 36248, 26198, 37718,
                                                 16248, 16458, 27728, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 39986, 0, 3, 36248, 36458, 26423, 38033,
                                                 16458, 16668, 28043, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 40427, 0, 3, 36458, 36668, 26648, 38348,
                                                 16668, 16878, 28358, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 40868, 0, 3, 37088, 37403, 27413, 39545,
                                                 17298, 17578, 29093, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 41456, 0, 3, 37403, 37718, 27728, 39986,
                                                 17578, 17858, 29513, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 42044, 0, 3, 37718, 38033, 28043, 40427,
                                                 17858, 18138, 29933, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 42632, 0, 3, 38663, 39104, 28673, 40868,
                                                 18698, 19058, 30353, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 43388, 0, 3, 39104, 39545, 29093, 41456,
                                                 19058, 19418, 30893, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 44144, 0, 3, 39545, 39986, 29513, 42044,
                                                 19418, 19778, 31433, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 44900, 0, 3, 40868, 41456, 30893, 44144,
                                                 20498, 20948, 32648, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 45845, 0, 3, 42632, 43388, 31973, 44900,
                                                 21848, 22398, 33323, ncols, alpha, beta, p);

            compute_prim_geom_10_lh_electron_repulsion_0(buffer, 47000, 42632, 45845, ncols,
                                                         alpha);

            compute_prim_geom_10_lh_electron_repulsion_1(buffer, 47945, 42632, 45845, ncols,
                                                         alpha);

            compute_prim_geom_10_lh_electron_repulsion_2(buffer, 48890, 42632, 45845, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 49835, 47000, 2835, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 52670, 49835, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 52670, 11, nmax);

    simdtrf::transform_h_inner(buffer, 52670, 50780, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 187 * nvalues, nvalues, buffer, 52670, 11, nmax);

    simdtrf::transform_h_inner(buffer, 52670, 51725, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 374 * nvalues, nvalues, buffer, 52670, 11, nmax);
}

}  // namespace simdt2ceri

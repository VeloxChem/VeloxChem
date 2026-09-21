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


#include "SimdElectronRepulsionGeom10RecLK.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMI.hpp"
#include "SimdElectronRepulsionVrrRecMK.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
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
#include "SimdGeometryL1.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_lk_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_lk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 123947, 118412, 4860, nvalues);

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
                                            10, 11, 12, 13, 14, 15, 16}, ncols, fj, mu);

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

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 22, ncols);

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

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 14, 15, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 15, 16, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 16, 17, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 17, 18, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 18, 19, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 19, 20, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 20, 21, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 23, 26, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 44, 47, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 47, 50, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 50, 53, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 53, 56, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 262, 0, 56, 59, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 59, 62, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 282, 0, 68, 74, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 297, 0, 74, 80, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 312, 0, 80, 86, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 327, 0, 86, 92, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 342, 0, 92, 98, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 357, 0, 98, 104, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 372, 0, 104, 110, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 387, 0, 110, 116, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 402, 0, 116, 122, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 417, 0, 122, 128, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 432, 0, 128, 134, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 447, 0, 134, 140, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 462, 0, 152, 162, 297, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 483, 0, 162, 172, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 504, 0, 172, 182, 327, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 525, 0, 182, 192, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 546, 0, 192, 202, 357, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 567, 0, 202, 212, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 588, 0, 212, 222, 387, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 609, 0, 222, 232, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 630, 0, 232, 242, 417, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 651, 0, 242, 252, 432, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 672, 0, 252, 262, 447, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 282, 297, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 721, 0, 297, 312, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 749, 0, 312, 327, 525, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 777, 0, 327, 342, 546, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 805, 0, 342, 357, 567, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 833, 0, 357, 372, 588, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 861, 0, 372, 387, 609, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 889, 0, 387, 402, 630, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 917, 0, 402, 417, 651, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 945, 0, 417, 432, 672, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 973, 0, 462, 483, 721, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1009, 0, 483, 504, 749, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1045, 0, 504, 525, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1081, 0, 525, 546, 805, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1117, 0, 546, 567, 833, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1153, 0, 567, 588, 861, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1189, 0, 588, 609, 889, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1225, 0, 609, 630, 917, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1261, 0, 630, 651, 945, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1297, 0, 693, 721, 1009, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1342, 0, 721, 749, 1045, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1387, 0, 749, 777, 1081, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1432, 0, 777, 805, 1117, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1477, 0, 805, 833, 1153, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1522, 0, 833, 861, 1189, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1567, 0, 861, 889, 1225, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1612, 0, 889, 917, 1261, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1657, 0, 973, 1009, 1342, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1712, 0, 1009, 1045, 1387, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1767, 0, 1045, 1081, 1432, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1822, 0, 1081, 1117, 1477, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1877, 0, 1117, 1153, 1522, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1932, 0, 1153, 1189, 1567, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1987, 0, 1189, 1225, 1612, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2042, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2045, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2048, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2051, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2054, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2057, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2060, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2063, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2066, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2069, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2072, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2075, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2078, 3, 22, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2081, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2090, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2099, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2108, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2117, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2126, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2135, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2144, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2153, 3, 17, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2162, 3, 18, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2171, 3, 19, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2180, 3, 20, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2189, 3, 21, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2198, 0, 3, 26, 2081, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2216, 0, 3, 29, 2090, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2234, 0, 3, 32, 2099, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2252, 0, 3, 35, 2108, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2270, 0, 3, 38, 2117, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2288, 0, 3, 41, 2126, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2306, 0, 3, 44, 2135, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2324, 0, 3, 47, 2144, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2342, 0, 3, 50, 2153, 122, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2360, 0, 3, 53, 2162, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2378, 0, 3, 56, 2171, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2396, 0, 3, 59, 2180, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2414, 0, 3, 62, 2189, 146, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2432, 0, 3, 68, 2198, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2462, 0, 3, 74, 2216, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2492, 0, 3, 80, 2234, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2522, 0, 3, 86, 2252, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2552, 0, 3, 92, 2270, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2582, 0, 3, 98, 2288, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2612, 0, 3, 104, 2306, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2642, 0, 3, 110, 2324, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2672, 0, 3, 116, 2342, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2702, 0, 3, 122, 2360, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2732, 0, 3, 128, 2378, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2762, 0, 3, 134, 2396, 262, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2792, 0, 3, 140, 2414, 272, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2822, 0, 3, 162, 2492, 297, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2867, 0, 3, 172, 2522, 312, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2912, 0, 3, 182, 2552, 327, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2957, 0, 3, 192, 2582, 342, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3002, 0, 3, 202, 2612, 357, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3047, 0, 3, 212, 2642, 372, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3092, 0, 3, 222, 2672, 387, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3137, 0, 3, 232, 2702, 402, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3182, 0, 3, 242, 2732, 417, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3227, 0, 3, 252, 2762, 432, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3272, 0, 3, 262, 2792, 447, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3317, 0, 3, 282, 2822, 462, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3380, 0, 3, 297, 2867, 483, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3443, 0, 3, 312, 2912, 504, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3506, 0, 3, 327, 2957, 525, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3569, 0, 3, 342, 3002, 546, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3632, 0, 3, 357, 3047, 567, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3695, 0, 3, 372, 3092, 588, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3758, 0, 3, 387, 3137, 609, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3821, 0, 3, 402, 3182, 630, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3884, 0, 3, 417, 3227, 651, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3947, 0, 3, 432, 3272, 672, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4010, 0, 3, 483, 3443, 721, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4094, 0, 3, 504, 3506, 749, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4178, 0, 3, 525, 3569, 777, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4262, 0, 3, 546, 3632, 805, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4346, 0, 3, 567, 3695, 833, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4430, 0, 3, 588, 3758, 861, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4514, 0, 3, 609, 3821, 889, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4598, 0, 3, 630, 3884, 917, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4682, 0, 3, 651, 3947, 945, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4766, 0, 3, 693, 4010, 973, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4874, 0, 3, 721, 4094, 1009, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4982, 0, 3, 749, 4178, 1045, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5090, 0, 3, 777, 4262, 1081, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5198, 0, 3, 805, 4346, 1117, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5306, 0, 3, 833, 4430, 1153, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5414, 0, 3, 861, 4514, 1189, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5522, 0, 3, 889, 4598, 1225, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5630, 0, 3, 917, 4682, 1261, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5738, 0, 3, 1009, 4982, 1342, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5873, 0, 3, 1045, 5090, 1387, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6008, 0, 3, 1081, 5198, 1432, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6143, 0, 3, 1117, 5306, 1477, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6278, 0, 3, 1153, 5414, 1522, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6413, 0, 3, 1189, 5522, 1567, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6548, 0, 3, 1225, 5630, 1612, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6683, 0, 3, 1297, 5738, 1657, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6848, 0, 3, 1342, 5873, 1712, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7013, 0, 3, 1387, 6008, 1767, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7178, 0, 3, 1432, 6143, 1822, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7343, 0, 3, 1477, 6278, 1877, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7508, 0, 3, 1522, 6413, 1932, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7673, 0, 3, 1567, 6548, 1987, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 7838, 3, 9, 10, 2045, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7844, 3, 10, 11, 2048, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7850, 3, 11, 12, 2051, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7856, 3, 12, 13, 2054, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7862, 3, 13, 14, 2057, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7868, 3, 14, 15, 2060, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7874, 3, 15, 16, 2063, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7880, 3, 16, 17, 2066, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7886, 3, 17, 18, 2069, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7892, 3, 18, 19, 2072, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7898, 3, 19, 20, 2075, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 7904, 3, 20, 21, 2078, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 7910, 0, 3, 2042, 7838, 2090, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7928, 0, 3, 2045, 7844, 2099, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7946, 0, 3, 2048, 7850, 2108, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7964, 0, 3, 2051, 7856, 2117, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7982, 0, 3, 2054, 7862, 2126, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8000, 0, 3, 2057, 7868, 2135, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8018, 0, 3, 2060, 7874, 2144, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8036, 0, 3, 2063, 7880, 2153, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8054, 0, 3, 2066, 7886, 2162, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8072, 0, 3, 2069, 7892, 2171, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8090, 0, 3, 2072, 7898, 2180, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8108, 0, 3, 2075, 7904, 2189, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 8126, 0, 3, 2081, 7910, 68, 74, 2216,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8162, 0, 3, 2090, 7928, 74, 80, 2234,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8198, 0, 3, 2099, 7946, 80, 86, 2252,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8234, 0, 3, 2108, 7964, 86, 92, 2270,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8270, 0, 3, 2117, 7982, 92, 98, 2288,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8306, 0, 3, 2126, 8000, 98, 104, 2306,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8342, 0, 3, 2135, 8018, 104, 110, 2324,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8378, 0, 3, 2144, 8036, 110, 116, 2342,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8414, 0, 3, 2153, 8054, 116, 122, 2360,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8450, 0, 3, 2162, 8072, 122, 128, 2378,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8486, 0, 3, 2171, 8090, 128, 134, 2396,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8522, 0, 3, 2180, 8108, 134, 140, 2414,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8558, 0, 3, 2216, 8162, 152, 162, 2492,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8618, 0, 3, 2234, 8198, 162, 172, 2522,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8678, 0, 3, 2252, 8234, 172, 182, 2552,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8738, 0, 3, 2270, 8270, 182, 192, 2582,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8798, 0, 3, 2288, 8306, 192, 202, 2612,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8858, 0, 3, 2306, 8342, 202, 212, 2642,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8918, 0, 3, 2324, 8378, 212, 222, 2672,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8978, 0, 3, 2342, 8414, 222, 232, 2702,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9038, 0, 3, 2360, 8450, 232, 242, 2732,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9098, 0, 3, 2378, 8486, 242, 252, 2762,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9158, 0, 3, 2396, 8522, 252, 262, 2792,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9218, 0, 3, 8126, 8162, 2492, 8618, 282,
                                                 297, 2867, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9308, 0, 3, 8162, 8198, 2522, 8678, 297,
                                                 312, 2912, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9398, 0, 3, 8198, 8234, 2552, 8738, 312,
                                                 327, 2957, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9488, 0, 3, 8234, 8270, 2582, 8798, 327,
                                                 342, 3002, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9578, 0, 3, 8270, 8306, 2612, 8858, 342,
                                                 357, 3047, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9668, 0, 3, 8306, 8342, 2642, 8918, 357,
                                                 372, 3092, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9758, 0, 3, 8342, 8378, 2672, 8978, 372,
                                                 387, 3137, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9848, 0, 3, 8378, 8414, 2702, 9038, 387,
                                                 402, 3182, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9938, 0, 3, 8414, 8450, 2732, 9098, 402,
                                                 417, 3227, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10028, 0, 3, 8450, 8486, 2762, 9158,
                                                 417, 432, 3272, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10118, 0, 3, 8558, 8618, 2867, 9308,
                                                 462, 483, 3443, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10244, 0, 3, 8618, 8678, 2912, 9398,
                                                 483, 504, 3506, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10370, 0, 3, 8678, 8738, 2957, 9488,
                                                 504, 525, 3569, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10496, 0, 3, 8738, 8798, 3002, 9578,
                                                 525, 546, 3632, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10622, 0, 3, 8798, 8858, 3047, 9668,
                                                 546, 567, 3695, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10748, 0, 3, 8858, 8918, 3092, 9758,
                                                 567, 588, 3758, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10874, 0, 3, 8918, 8978, 3137, 9848,
                                                 588, 609, 3821, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11000, 0, 3, 8978, 9038, 3182, 9938,
                                                 609, 630, 3884, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11126, 0, 3, 9038, 9098, 3227, 10028,
                                                 630, 651, 3947, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11252, 0, 3, 9218, 9308, 3443, 10244,
                                                 693, 721, 4094, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11420, 0, 3, 9308, 9398, 3506, 10370,
                                                 721, 749, 4178, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11588, 0, 3, 9398, 9488, 3569, 10496,
                                                 749, 777, 4262, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11756, 0, 3, 9488, 9578, 3632, 10622,
                                                 777, 805, 4346, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11924, 0, 3, 9578, 9668, 3695, 10748,
                                                 805, 833, 4430, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12092, 0, 3, 9668, 9758, 3758, 10874,
                                                 833, 861, 4514, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12260, 0, 3, 9758, 9848, 3821, 11000,
                                                 861, 889, 4598, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12428, 0, 3, 9848, 9938, 3884, 11126,
                                                 889, 917, 4682, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12596, 0, 3, 10118, 10244, 4094, 11420,
                                                 973, 1009, 4982, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12812, 0, 3, 10244, 10370, 4178, 11588,
                                                 1009, 1045, 5090, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13028, 0, 3, 10370, 10496, 4262, 11756,
                                                 1045, 1081, 5198, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13244, 0, 3, 10496, 10622, 4346, 11924,
                                                 1081, 1117, 5306, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13460, 0, 3, 10622, 10748, 4430, 12092,
                                                 1117, 1153, 5414, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13676, 0, 3, 10748, 10874, 4514, 12260,
                                                 1153, 1189, 5522, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13892, 0, 3, 10874, 11000, 4598, 12428,
                                                 1189, 1225, 5630, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14108, 0, 3, 11252, 11420, 4982, 12812,
                                                 1297, 1342, 5873, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14378, 0, 3, 11420, 11588, 5090, 13028,
                                                 1342, 1387, 6008, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14648, 0, 3, 11588, 11756, 5198, 13244,
                                                 1387, 1432, 6143, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14918, 0, 3, 11756, 11924, 5306, 13460,
                                                 1432, 1477, 6278, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15188, 0, 3, 11924, 12092, 5414, 13676,
                                                 1477, 1522, 6413, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15458, 0, 3, 12092, 12260, 5522, 13892,
                                                 1522, 1567, 6548, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 15728, 0, 3, 12596, 12812, 5873, 14378,
                                                 1657, 1712, 7013, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 16058, 0, 3, 12812, 13028, 6008, 14648,
                                                 1712, 1767, 7178, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 16388, 0, 3, 13028, 13244, 6143, 14918,
                                                 1767, 1822, 7343, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 16718, 0, 3, 13244, 13460, 6278, 15188,
                                                 1822, 1877, 7508, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 17048, 0, 3, 13460, 13676, 6413, 15458,
                                                 1877, 1932, 7673, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17378, 3, 2042, 2045, 7844, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17388, 3, 2045, 2048, 7850, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17398, 3, 2048, 2051, 7856, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17408, 3, 2051, 2054, 7862, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17418, 3, 2054, 2057, 7868, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17428, 3, 2057, 2060, 7874, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17438, 3, 2060, 2063, 7880, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17448, 3, 2063, 2066, 7886, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17458, 3, 2066, 2069, 7892, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17468, 3, 2069, 2072, 7898, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 17478, 3, 2072, 2075, 7904, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 17488, 0, 3, 7838, 17378, 7928, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17518, 0, 3, 7844, 17388, 7946, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17548, 0, 3, 7850, 17398, 7964, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17578, 0, 3, 7856, 17408, 7982, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17608, 0, 3, 7862, 17418, 8000, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17638, 0, 3, 7868, 17428, 8018, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17668, 0, 3, 7874, 17438, 8036, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17698, 0, 3, 7880, 17448, 8054, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17728, 0, 3, 7886, 17458, 8072, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17758, 0, 3, 7892, 17468, 8090, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 17788, 0, 3, 7898, 17478, 8108, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 17818, 0, 3, 7910, 17488, 2198, 2216,
                                                 8162, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17878, 0, 3, 7928, 17518, 2216, 2234,
                                                 8198, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17938, 0, 3, 7946, 17548, 2234, 2252,
                                                 8234, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17998, 0, 3, 7964, 17578, 2252, 2270,
                                                 8270, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18058, 0, 3, 7982, 17608, 2270, 2288,
                                                 8306, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18118, 0, 3, 8000, 17638, 2288, 2306,
                                                 8342, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18178, 0, 3, 8018, 17668, 2306, 2324,
                                                 8378, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18238, 0, 3, 8036, 17698, 2324, 2342,
                                                 8414, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18298, 0, 3, 8054, 17728, 2342, 2360,
                                                 8450, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18358, 0, 3, 8072, 17758, 2360, 2378,
                                                 8486, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 18418, 0, 3, 8090, 17788, 2378, 2396,
                                                 8522, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18478, 0, 3, 8126, 17818, 2432, 2462,
                                                 8558, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18578, 0, 3, 8162, 17878, 2462, 2492,
                                                 8618, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18678, 0, 3, 8198, 17938, 2492, 2522,
                                                 8678, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18778, 0, 3, 8234, 17998, 2522, 2552,
                                                 8738, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18878, 0, 3, 8270, 18058, 2552, 2582,
                                                 8798, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18978, 0, 3, 8306, 18118, 2582, 2612,
                                                 8858, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19078, 0, 3, 8342, 18178, 2612, 2642,
                                                 8918, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19178, 0, 3, 8378, 18238, 2642, 2672,
                                                 8978, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19278, 0, 3, 8414, 18298, 2672, 2702,
                                                 9038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19378, 0, 3, 8450, 18358, 2702, 2732,
                                                 9098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 19478, 0, 3, 8486, 18418, 2732, 2762,
                                                 9158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19578, 0, 3, 17818, 17878, 8618, 18678,
                                                 2822, 2867, 9308, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19728, 0, 3, 17878, 17938, 8678, 18778,
                                                 2867, 2912, 9398, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19878, 0, 3, 17938, 17998, 8738, 18878,
                                                 2912, 2957, 9488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20028, 0, 3, 17998, 18058, 8798, 18978,
                                                 2957, 3002, 9578, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20178, 0, 3, 18058, 18118, 8858, 19078,
                                                 3002, 3047, 9668, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20328, 0, 3, 18118, 18178, 8918, 19178,
                                                 3047, 3092, 9758, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20478, 0, 3, 18178, 18238, 8978, 19278,
                                                 3092, 3137, 9848, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20628, 0, 3, 18238, 18298, 9038, 19378,
                                                 3137, 3182, 9938, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20778, 0, 3, 18298, 18358, 9098, 19478,
                                                 3182, 3227, 10028, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20928, 0, 3, 18478, 18578, 9218, 19578,
                                                 3317, 3380, 10118, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21138, 0, 3, 18578, 18678, 9308, 19728,
                                                 3380, 3443, 10244, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21348, 0, 3, 18678, 18778, 9398, 19878,
                                                 3443, 3506, 10370, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21558, 0, 3, 18778, 18878, 9488, 20028,
                                                 3506, 3569, 10496, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21768, 0, 3, 18878, 18978, 9578, 20178,
                                                 3569, 3632, 10622, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21978, 0, 3, 18978, 19078, 9668, 20328,
                                                 3632, 3695, 10748, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22188, 0, 3, 19078, 19178, 9758, 20478,
                                                 3695, 3758, 10874, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22398, 0, 3, 19178, 19278, 9848, 20628,
                                                 3758, 3821, 11000, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22608, 0, 3, 19278, 19378, 9938, 20778,
                                                 3821, 3884, 11126, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22818, 0, 3, 19578, 19728, 10244, 21348,
                                                 4010, 4094, 11420, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23098, 0, 3, 19728, 19878, 10370, 21558,
                                                 4094, 4178, 11588, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23378, 0, 3, 19878, 20028, 10496, 21768,
                                                 4178, 4262, 11756, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23658, 0, 3, 20028, 20178, 10622, 21978,
                                                 4262, 4346, 11924, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23938, 0, 3, 20178, 20328, 10748, 22188,
                                                 4346, 4430, 12092, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24218, 0, 3, 20328, 20478, 10874, 22398,
                                                 4430, 4514, 12260, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24498, 0, 3, 20478, 20628, 11000, 22608,
                                                 4514, 4598, 12428, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 24778, 0, 3, 20928, 21138, 11252, 22818,
                                                 4766, 4874, 12596, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25138, 0, 3, 21138, 21348, 11420, 23098,
                                                 4874, 4982, 12812, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25498, 0, 3, 21348, 21558, 11588, 23378,
                                                 4982, 5090, 13028, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25858, 0, 3, 21558, 21768, 11756, 23658,
                                                 5090, 5198, 13244, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26218, 0, 3, 21768, 21978, 11924, 23938,
                                                 5198, 5306, 13460, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26578, 0, 3, 21978, 22188, 12092, 24218,
                                                 5306, 5414, 13676, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26938, 0, 3, 22188, 22398, 12260, 24498,
                                                 5414, 5522, 13892, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 27298, 0, 3, 22818, 23098, 12812, 25498,
                                                 5738, 5873, 14378, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 27748, 0, 3, 23098, 23378, 13028, 25858,
                                                 5873, 6008, 14648, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 28198, 0, 3, 23378, 23658, 13244, 26218,
                                                 6008, 6143, 14918, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 28648, 0, 3, 23658, 23938, 13460, 26578,
                                                 6143, 6278, 15188, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 29098, 0, 3, 23938, 24218, 13676, 26938,
                                                 6278, 6413, 15458, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 29548, 0, 3, 24778, 25138, 14108, 27298,
                                                 6683, 6848, 15728, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 30098, 0, 3, 25138, 25498, 14378, 27748,
                                                 6848, 7013, 16058, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 30648, 0, 3, 25498, 25858, 14648, 28198,
                                                 7013, 7178, 16388, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 31198, 0, 3, 25858, 26218, 14918, 28648,
                                                 7178, 7343, 16718, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 31748, 0, 3, 26218, 26578, 15188, 29098,
                                                 7343, 7508, 17048, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32298, 3, 7838, 7844, 17388, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32313, 3, 7844, 7850, 17398, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32328, 3, 7850, 7856, 17408, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32343, 3, 7856, 7862, 17418, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32358, 3, 7862, 7868, 17428, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32373, 3, 7868, 7874, 17438, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32388, 3, 7874, 7880, 17448, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32403, 3, 7880, 7886, 17458, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32418, 3, 7886, 7892, 17468, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 32433, 3, 7892, 7898, 17478, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 32448, 0, 3, 17378, 32298, 17518, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32493, 0, 3, 17388, 32313, 17548, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32538, 0, 3, 17398, 32328, 17578, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32583, 0, 3, 17408, 32343, 17608, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32628, 0, 3, 17418, 32358, 17638, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32673, 0, 3, 17428, 32373, 17668, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32718, 0, 3, 17438, 32388, 17698, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32763, 0, 3, 17448, 32403, 17728, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32808, 0, 3, 17458, 32418, 17758, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 32853, 0, 3, 17468, 32433, 17788, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 32898, 0, 3, 17488, 32448, 8126, 8162,
                                                 17878, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 32988, 0, 3, 17518, 32493, 8162, 8198,
                                                 17938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33078, 0, 3, 17548, 32538, 8198, 8234,
                                                 17998, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33168, 0, 3, 17578, 32583, 8234, 8270,
                                                 18058, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33258, 0, 3, 17608, 32628, 8270, 8306,
                                                 18118, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33348, 0, 3, 17638, 32673, 8306, 8342,
                                                 18178, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33438, 0, 3, 17668, 32718, 8342, 8378,
                                                 18238, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33528, 0, 3, 17698, 32763, 8378, 8414,
                                                 18298, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33618, 0, 3, 17728, 32808, 8414, 8450,
                                                 18358, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 33708, 0, 3, 17758, 32853, 8450, 8486,
                                                 18418, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 33798, 0, 3, 17878, 32988, 8558, 8618,
                                                 18678, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 33948, 0, 3, 17938, 33078, 8618, 8678,
                                                 18778, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34098, 0, 3, 17998, 33168, 8678, 8738,
                                                 18878, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34248, 0, 3, 18058, 33258, 8738, 8798,
                                                 18978, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34398, 0, 3, 18118, 33348, 8798, 8858,
                                                 19078, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34548, 0, 3, 18178, 33438, 8858, 8918,
                                                 19178, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34698, 0, 3, 18238, 33528, 8918, 8978,
                                                 19278, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34848, 0, 3, 18298, 33618, 8978, 9038,
                                                 19378, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 34998, 0, 3, 18358, 33708, 9038, 9098,
                                                 19478, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 35148, 0, 3, 32898, 32988, 18678, 33948,
                                                 9218, 9308, 19728, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 35373, 0, 3, 32988, 33078, 18778, 34098,
                                                 9308, 9398, 19878, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 35598, 0, 3, 33078, 33168, 18878, 34248,
                                                 9398, 9488, 20028, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 35823, 0, 3, 33168, 33258, 18978, 34398,
                                                 9488, 9578, 20178, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 36048, 0, 3, 33258, 33348, 19078, 34548,
                                                 9578, 9668, 20328, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 36273, 0, 3, 33348, 33438, 19178, 34698,
                                                 9668, 9758, 20478, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 36498, 0, 3, 33438, 33528, 19278, 34848,
                                                 9758, 9848, 20628, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 36723, 0, 3, 33528, 33618, 19378, 34998,
                                                 9848, 9938, 20778, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 36948, 0, 3, 33798, 33948, 19728, 35373,
                                                 10118, 10244, 21348, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 37263, 0, 3, 33948, 34098, 19878, 35598,
                                                 10244, 10370, 21558, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 37578, 0, 3, 34098, 34248, 20028, 35823,
                                                 10370, 10496, 21768, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 37893, 0, 3, 34248, 34398, 20178, 36048,
                                                 10496, 10622, 21978, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 38208, 0, 3, 34398, 34548, 20328, 36273,
                                                 10622, 10748, 22188, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 38523, 0, 3, 34548, 34698, 20478, 36498,
                                                 10748, 10874, 22398, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 38838, 0, 3, 34698, 34848, 20628, 36723,
                                                 10874, 11000, 22608, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39153, 0, 3, 35148, 35373, 21348, 37263,
                                                 11252, 11420, 23098, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39573, 0, 3, 35373, 35598, 21558, 37578,
                                                 11420, 11588, 23378, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39993, 0, 3, 35598, 35823, 21768, 37893,
                                                 11588, 11756, 23658, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 40413, 0, 3, 35823, 36048, 21978, 38208,
                                                 11756, 11924, 23938, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 40833, 0, 3, 36048, 36273, 22188, 38523,
                                                 11924, 12092, 24218, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 41253, 0, 3, 36273, 36498, 22398, 38838,
                                                 12092, 12260, 24498, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 41673, 0, 3, 36948, 37263, 23098, 39573,
                                                 12596, 12812, 25498, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 42213, 0, 3, 37263, 37578, 23378, 39993,
                                                 12812, 13028, 25858, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 42753, 0, 3, 37578, 37893, 23658, 40413,
                                                 13028, 13244, 26218, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 43293, 0, 3, 37893, 38208, 23938, 40833,
                                                 13244, 13460, 26578, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 43833, 0, 3, 38208, 38523, 24218, 41253,
                                                 13460, 13676, 26938, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 44373, 0, 3, 39153, 39573, 25498, 42213,
                                                 14108, 14378, 27748, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 45048, 0, 3, 39573, 39993, 25858, 42753,
                                                 14378, 14648, 28198, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 45723, 0, 3, 39993, 40413, 26218, 43293,
                                                 14648, 14918, 28648, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 46398, 0, 3, 40413, 40833, 26578, 43833,
                                                 14918, 15188, 29098, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 47073, 0, 3, 41673, 42213, 27748, 45048,
                                                 15728, 16058, 30648, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 47898, 0, 3, 42213, 42753, 28198, 45723,
                                                 16058, 16388, 31198, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 48723, 0, 3, 42753, 43293, 28648, 46398,
                                                 16388, 16718, 31748, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49548, 3, 17378, 17388, 32313, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49569, 3, 17388, 17398, 32328, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49590, 3, 17398, 17408, 32343, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49611, 3, 17408, 17418, 32358, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49632, 3, 17418, 17428, 32373, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49653, 3, 17428, 17438, 32388, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49674, 3, 17438, 17448, 32403, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49695, 3, 17448, 17458, 32418, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 49716, 3, 17458, 17468, 32433, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 49737, 0, 3, 32298, 49548, 32493, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49800, 0, 3, 32313, 49569, 32538, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49863, 0, 3, 32328, 49590, 32583, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49926, 0, 3, 32343, 49611, 32628, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49989, 0, 3, 32358, 49632, 32673, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50052, 0, 3, 32373, 49653, 32718, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50115, 0, 3, 32388, 49674, 32763, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50178, 0, 3, 32403, 49695, 32808, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 50241, 0, 3, 32418, 49716, 32853, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 50304, 0, 3, 32448, 49737, 17818, 17878,
                                                 32988, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50430, 0, 3, 32493, 49800, 17878, 17938,
                                                 33078, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50556, 0, 3, 32538, 49863, 17938, 17998,
                                                 33168, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50682, 0, 3, 32583, 49926, 17998, 18058,
                                                 33258, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50808, 0, 3, 32628, 49989, 18058, 18118,
                                                 33348, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50934, 0, 3, 32673, 50052, 18118, 18178,
                                                 33438, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51060, 0, 3, 32718, 50115, 18178, 18238,
                                                 33528, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51186, 0, 3, 32763, 50178, 18238, 18298,
                                                 33618, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51312, 0, 3, 32808, 50241, 18298, 18358,
                                                 33708, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51438, 0, 3, 32898, 50304, 18478, 18578,
                                                 33798, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51648, 0, 3, 32988, 50430, 18578, 18678,
                                                 33948, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51858, 0, 3, 33078, 50556, 18678, 18778,
                                                 34098, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52068, 0, 3, 33168, 50682, 18778, 18878,
                                                 34248, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52278, 0, 3, 33258, 50808, 18878, 18978,
                                                 34398, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52488, 0, 3, 33348, 50934, 18978, 19078,
                                                 34548, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52698, 0, 3, 33438, 51060, 19078, 19178,
                                                 34698, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52908, 0, 3, 33528, 51186, 19178, 19278,
                                                 34848, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53118, 0, 3, 33618, 51312, 19278, 19378,
                                                 34998, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 53328, 0, 3, 50304, 50430, 33948, 51858,
                                                 19578, 19728, 35373, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 53643, 0, 3, 50430, 50556, 34098, 52068,
                                                 19728, 19878, 35598, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 53958, 0, 3, 50556, 50682, 34248, 52278,
                                                 19878, 20028, 35823, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54273, 0, 3, 50682, 50808, 34398, 52488,
                                                 20028, 20178, 36048, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54588, 0, 3, 50808, 50934, 34548, 52698,
                                                 20178, 20328, 36273, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54903, 0, 3, 50934, 51060, 34698, 52908,
                                                 20328, 20478, 36498, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55218, 0, 3, 51060, 51186, 34848, 53118,
                                                 20478, 20628, 36723, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 55533, 0, 3, 51438, 51648, 35148, 53328,
                                                 20928, 21138, 36948, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 55974, 0, 3, 51648, 51858, 35373, 53643,
                                                 21138, 21348, 37263, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 56415, 0, 3, 51858, 52068, 35598, 53958,
                                                 21348, 21558, 37578, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 56856, 0, 3, 52068, 52278, 35823, 54273,
                                                 21558, 21768, 37893, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 57297, 0, 3, 52278, 52488, 36048, 54588,
                                                 21768, 21978, 38208, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 57738, 0, 3, 52488, 52698, 36273, 54903,
                                                 21978, 22188, 38523, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 58179, 0, 3, 52698, 52908, 36498, 55218,
                                                 22188, 22398, 38838, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 58620, 0, 3, 53328, 53643, 37263, 56415,
                                                 22818, 23098, 39573, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 59208, 0, 3, 53643, 53958, 37578, 56856,
                                                 23098, 23378, 39993, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 59796, 0, 3, 53958, 54273, 37893, 57297,
                                                 23378, 23658, 40413, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 60384, 0, 3, 54273, 54588, 38208, 57738,
                                                 23658, 23938, 40833, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 60972, 0, 3, 54588, 54903, 38523, 58179,
                                                 23938, 24218, 41253, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 61560, 0, 3, 55533, 55974, 39153, 58620,
                                                 24778, 25138, 41673, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 62316, 0, 3, 55974, 56415, 39573, 59208,
                                                 25138, 25498, 42213, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 63072, 0, 3, 56415, 56856, 39993, 59796,
                                                 25498, 25858, 42753, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 63828, 0, 3, 56856, 57297, 40413, 60384,
                                                 25858, 26218, 43293, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 64584, 0, 3, 57297, 57738, 40833, 60972,
                                                 26218, 26578, 43833, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 65340, 0, 3, 58620, 59208, 42213, 63072,
                                                 27298, 27748, 45048, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 66285, 0, 3, 59208, 59796, 42753, 63828,
                                                 27748, 28198, 45723, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 67230, 0, 3, 59796, 60384, 43293, 64584,
                                                 28198, 28648, 46398, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 68175, 0, 3, 61560, 62316, 44373, 65340,
                                                 29548, 30098, 47073, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 69330, 0, 3, 62316, 63072, 45048, 66285,
                                                 30098, 30648, 47898, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 70485, 0, 3, 63072, 63828, 45723, 67230,
                                                 30648, 31198, 48723, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71640, 3, 32298, 32313, 49569, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71668, 3, 32313, 32328, 49590, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71696, 3, 32328, 32343, 49611, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71724, 3, 32343, 32358, 49632, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71752, 3, 32358, 32373, 49653, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71780, 3, 32373, 32388, 49674, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71808, 3, 32388, 32403, 49695, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 71836, 3, 32403, 32418, 49716, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 71864, 0, 3, 49548, 71640, 49800, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 71948, 0, 3, 49569, 71668, 49863, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72032, 0, 3, 49590, 71696, 49926, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72116, 0, 3, 49611, 71724, 49989, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72200, 0, 3, 49632, 71752, 50052, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72284, 0, 3, 49653, 71780, 50115, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72368, 0, 3, 49674, 71808, 50178, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 72452, 0, 3, 49695, 71836, 50241, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 72536, 0, 3, 49737, 71864, 32898, 32988,
                                                 50430, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 72704, 0, 3, 49800, 71948, 32988, 33078,
                                                 50556, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 72872, 0, 3, 49863, 72032, 33078, 33168,
                                                 50682, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73040, 0, 3, 49926, 72116, 33168, 33258,
                                                 50808, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73208, 0, 3, 49989, 72200, 33258, 33348,
                                                 50934, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73376, 0, 3, 50052, 72284, 33348, 33438,
                                                 51060, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73544, 0, 3, 50115, 72368, 33438, 33528,
                                                 51186, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 73712, 0, 3, 50178, 72452, 33528, 33618,
                                                 51312, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 73880, 0, 3, 50430, 72704, 33798, 33948,
                                                 51858, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 74160, 0, 3, 50556, 72872, 33948, 34098,
                                                 52068, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 74440, 0, 3, 50682, 73040, 34098, 34248,
                                                 52278, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 74720, 0, 3, 50808, 73208, 34248, 34398,
                                                 52488, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75000, 0, 3, 50934, 73376, 34398, 34548,
                                                 52698, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75280, 0, 3, 51060, 73544, 34548, 34698,
                                                 52908, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 75560, 0, 3, 51186, 73712, 34698, 34848,
                                                 53118, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 75840, 0, 3, 72536, 72704, 51858, 74160,
                                                 35148, 35373, 53643, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 76260, 0, 3, 72704, 72872, 52068, 74440,
                                                 35373, 35598, 53958, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 76680, 0, 3, 72872, 73040, 52278, 74720,
                                                 35598, 35823, 54273, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 77100, 0, 3, 73040, 73208, 52488, 75000,
                                                 35823, 36048, 54588, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 77520, 0, 3, 73208, 73376, 52698, 75280,
                                                 36048, 36273, 54903, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 77940, 0, 3, 73376, 73544, 52908, 75560,
                                                 36273, 36498, 55218, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 78360, 0, 3, 73880, 74160, 53643, 76260,
                                                 36948, 37263, 56415, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 78948, 0, 3, 74160, 74440, 53958, 76680,
                                                 37263, 37578, 56856, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 79536, 0, 3, 74440, 74720, 54273, 77100,
                                                 37578, 37893, 57297, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 80124, 0, 3, 74720, 75000, 54588, 77520,
                                                 37893, 38208, 57738, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 80712, 0, 3, 75000, 75280, 54903, 77940,
                                                 38208, 38523, 58179, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 81300, 0, 3, 75840, 76260, 56415, 78948,
                                                 39153, 39573, 59208, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 82084, 0, 3, 76260, 76680, 56856, 79536,
                                                 39573, 39993, 59796, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 82868, 0, 3, 76680, 77100, 57297, 80124,
                                                 39993, 40413, 60384, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 83652, 0, 3, 77100, 77520, 57738, 80712,
                                                 40413, 40833, 60972, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 84436, 0, 3, 78360, 78948, 59208, 82084,
                                                 41673, 42213, 63072, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 85444, 0, 3, 78948, 79536, 59796, 82868,
                                                 42213, 42753, 63828, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 86452, 0, 3, 79536, 80124, 60384, 83652,
                                                 42753, 43293, 64584, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 87460, 0, 3, 81300, 82084, 63072, 85444,
                                                 44373, 45048, 66285, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 88720, 0, 3, 82084, 82868, 63828, 86452,
                                                 45048, 45723, 67230, ncols, alpha, beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 89980, 0, 3, 84436, 85444, 66285, 88720,
                                                 47073, 47898, 70485, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91520, 3, 49548, 49569, 71668, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91556, 3, 49569, 49590, 71696, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91592, 3, 49590, 49611, 71724, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91628, 3, 49611, 49632, 71752, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91664, 3, 49632, 49653, 71780, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91700, 3, 49653, 49674, 71808, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 91736, 3, 49674, 49695, 71836, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 91772, 0, 3, 71640, 91520, 71948, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 91880, 0, 3, 71668, 91556, 72032, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 91988, 0, 3, 71696, 91592, 72116, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 92096, 0, 3, 71724, 91628, 72200, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 92204, 0, 3, 71752, 91664, 72284, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 92312, 0, 3, 71780, 91700, 72368, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 92420, 0, 3, 71808, 91736, 72452, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 92528, 0, 3, 71864, 91772, 50304, 50430,
                                                 72704, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92744, 0, 3, 71948, 91880, 50430, 50556,
                                                 72872, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 92960, 0, 3, 72032, 91988, 50556, 50682,
                                                 73040, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 93176, 0, 3, 72116, 92096, 50682, 50808,
                                                 73208, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 93392, 0, 3, 72200, 92204, 50808, 50934,
                                                 73376, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 93608, 0, 3, 72284, 92312, 50934, 51060,
                                                 73544, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 93824, 0, 3, 72368, 92420, 51060, 51186,
                                                 73712, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94040, 0, 3, 72536, 92528, 51438, 51648,
                                                 73880, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94400, 0, 3, 72704, 92744, 51648, 51858,
                                                 74160, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 94760, 0, 3, 72872, 92960, 51858, 52068,
                                                 74440, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95120, 0, 3, 73040, 93176, 52068, 52278,
                                                 74720, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95480, 0, 3, 73208, 93392, 52278, 52488,
                                                 75000, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 95840, 0, 3, 73376, 93608, 52488, 52698,
                                                 75280, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 96200, 0, 3, 73544, 93824, 52698, 52908,
                                                 75560, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 96560, 0, 3, 92528, 92744, 74160, 94760,
                                                 53328, 53643, 76260, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 97100, 0, 3, 92744, 92960, 74440, 95120,
                                                 53643, 53958, 76680, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 97640, 0, 3, 92960, 93176, 74720, 95480,
                                                 53958, 54273, 77100, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 98180, 0, 3, 93176, 93392, 75000, 95840,
                                                 54273, 54588, 77520, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 98720, 0, 3, 93392, 93608, 75280, 96200,
                                                 54588, 54903, 77940, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 99260, 0, 3, 94040, 94400, 75840, 96560,
                                                 55533, 55974, 78360, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 100016, 0, 3, 94400, 94760, 76260,
                                                 97100, 55974, 56415, 78948, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 100772, 0, 3, 94760, 95120, 76680,
                                                 97640, 56415, 56856, 79536, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 101528, 0, 3, 95120, 95480, 77100,
                                                 98180, 56856, 57297, 80124, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 102284, 0, 3, 95480, 95840, 77520,
                                                 98720, 57297, 57738, 80712, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 103040, 0, 3, 96560, 97100, 78948,
                                                 100772, 58620, 59208, 82084, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 104048, 0, 3, 97100, 97640, 79536,
                                                 101528, 59208, 59796, 82868, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 105056, 0, 3, 97640, 98180, 80124,
                                                 102284, 59796, 60384, 83652, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 106064, 0, 3, 99260, 100016, 81300,
                                                 103040, 61560, 62316, 84436, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 107360, 0, 3, 100016, 100772, 82084,
                                                 104048, 62316, 63072, 85444, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 108656, 0, 3, 100772, 101528, 82868,
                                                 105056, 63072, 63828, 86452, ncols, alpha, beta,
                                                 p);

            compute_prim_lk_electron_repulsion_0(buffer, 109952, 0, 3, 103040, 104048, 85444,
                                                 108656, 65340, 66285, 88720, ncols, alpha, beta,
                                                 p);

            compute_prim_mk_electron_repulsion_0(buffer, 111572, 0, 3, 106064, 107360, 87460,
                                                 109952, 68175, 69330, 89980, ncols, alpha, beta,
                                                 p);

            simdgeo::geom_l_x(buffer, 113552, 106064, 111572, 1, 36, ncols, alpha);

            simdgeo::geom_l_y(buffer, 115172, 106064, 111572, 1, 36, ncols, alpha);

            simdgeo::geom_l_z(buffer, 116792, 106064, 111572, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 118412, 113552, 4860, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 123272, 118412, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 123272, 15, nmax);

    simdtrf::transform_k_inner(buffer, 123272, 120032, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 255 * nvalues, nvalues, buffer, 123272, 15, nmax);

    simdtrf::transform_k_inner(buffer, 123272, 121652, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 510 * nvalues, nvalues, buffer, 123272, 15, nmax);
}

}  // namespace simdt2ceri

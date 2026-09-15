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


#include "SimdElectronRepulsionGeom10RecLI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMI.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_li_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_li_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 83140, 78775, 3780, nvalues);

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
                                            10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 7, 8, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 8, 9, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 9, 10, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 85, 0, 10, 11, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 91, 0, 11, 12, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 97, 0, 12, 13, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 103, 0, 13, 14, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 109, 0, 14, 15, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 115, 0, 15, 16, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 121, 0, 16, 17, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 127, 0, 17, 18, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 133, 0, 18, 19, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 139, 0, 19, 20, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 22, 25, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 25, 28, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 28, 31, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 31, 34, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 34, 37, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 37, 40, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 40, 43, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 43, 46, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 46, 49, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 235, 0, 49, 52, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 245, 0, 52, 55, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 255, 0, 55, 58, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 265, 0, 58, 61, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 275, 0, 67, 73, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 290, 0, 73, 79, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 305, 0, 79, 85, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 320, 0, 85, 91, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 335, 0, 91, 97, 205, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 350, 0, 97, 103, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 365, 0, 103, 109, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 380, 0, 109, 115, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 395, 0, 115, 121, 245, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 410, 0, 121, 127, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 425, 0, 127, 133, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 440, 0, 145, 155, 275, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 461, 0, 155, 165, 290, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 482, 0, 165, 175, 305, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 503, 0, 175, 185, 320, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 524, 0, 185, 195, 335, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 545, 0, 195, 205, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 566, 0, 205, 215, 365, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 587, 0, 215, 225, 380, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 608, 0, 225, 235, 395, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 629, 0, 235, 245, 410, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 650, 0, 245, 255, 425, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 671, 0, 275, 290, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 699, 0, 290, 305, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 727, 0, 305, 320, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 755, 0, 320, 335, 545, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 783, 0, 335, 350, 566, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 811, 0, 350, 365, 587, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 839, 0, 365, 380, 608, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 867, 0, 380, 395, 629, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 895, 0, 395, 410, 650, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 923, 0, 440, 461, 671, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 959, 0, 461, 482, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 995, 0, 482, 503, 727, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1031, 0, 503, 524, 755, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1067, 0, 524, 545, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1103, 0, 545, 566, 811, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1139, 0, 566, 587, 839, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1175, 0, 587, 608, 867, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1211, 0, 608, 629, 895, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1247, 0, 671, 699, 995, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1292, 0, 699, 727, 1031, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1337, 0, 727, 755, 1067, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1382, 0, 755, 783, 1103, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1427, 0, 783, 811, 1139, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1472, 0, 811, 839, 1175, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1517, 0, 839, 867, 1211, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1562, 0, 923, 959, 1247, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1617, 0, 959, 995, 1292, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1672, 0, 995, 1031, 1337, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1727, 0, 1031, 1067, 1382, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1782, 0, 1067, 1103, 1427, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1837, 0, 1103, 1139, 1472, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1892, 0, 1139, 1175, 1517, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1947, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1950, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1953, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1956, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1959, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1962, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1965, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1968, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1971, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1974, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1977, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1980, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1983, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1992, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2001, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2010, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2019, 3, 13, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2028, 3, 14, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2037, 3, 15, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2046, 3, 16, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2055, 3, 17, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2064, 3, 18, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2073, 3, 19, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2082, 3, 20, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2091, 0, 3, 28, 1983, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2109, 0, 3, 31, 1992, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2127, 0, 3, 34, 2001, 85, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2145, 0, 3, 37, 2010, 91, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2163, 0, 3, 40, 2019, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2181, 0, 3, 43, 2028, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2199, 0, 3, 46, 2037, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2217, 0, 3, 49, 2046, 115, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2235, 0, 3, 52, 2055, 121, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2253, 0, 3, 55, 2064, 127, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2271, 0, 3, 58, 2073, 133, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2289, 0, 3, 61, 2082, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2307, 0, 3, 73, 2109, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2337, 0, 3, 79, 2127, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2367, 0, 3, 85, 2145, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2397, 0, 3, 91, 2163, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2427, 0, 3, 97, 2181, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2457, 0, 3, 103, 2199, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2487, 0, 3, 109, 2217, 225, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2517, 0, 3, 115, 2235, 235, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2547, 0, 3, 121, 2253, 245, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2577, 0, 3, 127, 2271, 255, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2607, 0, 3, 133, 2289, 265, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2637, 0, 3, 165, 2337, 290, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2682, 0, 3, 175, 2367, 305, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2727, 0, 3, 185, 2397, 320, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2772, 0, 3, 195, 2427, 335, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2817, 0, 3, 205, 2457, 350, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2862, 0, 3, 215, 2487, 365, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2907, 0, 3, 225, 2517, 380, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2952, 0, 3, 235, 2547, 395, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2997, 0, 3, 245, 2577, 410, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3042, 0, 3, 255, 2607, 425, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3087, 0, 3, 290, 2682, 482, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3150, 0, 3, 305, 2727, 503, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3213, 0, 3, 320, 2772, 524, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3276, 0, 3, 335, 2817, 545, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3339, 0, 3, 350, 2862, 566, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3402, 0, 3, 365, 2907, 587, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3465, 0, 3, 380, 2952, 608, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3528, 0, 3, 395, 2997, 629, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3591, 0, 3, 410, 3042, 650, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3654, 0, 3, 482, 3150, 699, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3738, 0, 3, 503, 3213, 727, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3822, 0, 3, 524, 3276, 755, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3906, 0, 3, 545, 3339, 783, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3990, 0, 3, 566, 3402, 811, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4074, 0, 3, 587, 3465, 839, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4158, 0, 3, 608, 3528, 867, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4242, 0, 3, 629, 3591, 895, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4326, 0, 3, 699, 3738, 995, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4434, 0, 3, 727, 3822, 1031, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4542, 0, 3, 755, 3906, 1067, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4650, 0, 3, 783, 3990, 1103, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4758, 0, 3, 811, 4074, 1139, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4866, 0, 3, 839, 4158, 1175, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4974, 0, 3, 867, 4242, 1211, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5082, 0, 3, 995, 4434, 1292, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5217, 0, 3, 1031, 4542, 1337, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5352, 0, 3, 1067, 4650, 1382, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5487, 0, 3, 1103, 4758, 1427, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5622, 0, 3, 1139, 4866, 1472, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5757, 0, 3, 1175, 4974, 1517, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 5892, 0, 3, 1292, 5217, 1672, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6057, 0, 3, 1337, 5352, 1727, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6222, 0, 3, 1382, 5487, 1782, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6387, 0, 3, 1427, 5622, 1837, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6552, 0, 3, 1472, 5757, 1892, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6717, 3, 9, 10, 1950, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6723, 3, 10, 11, 1953, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6729, 3, 11, 12, 1956, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6735, 3, 12, 13, 1959, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6741, 3, 13, 14, 1962, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6747, 3, 14, 15, 1965, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6753, 3, 15, 16, 1968, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6759, 3, 16, 17, 1971, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6765, 3, 17, 18, 1974, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6771, 3, 18, 19, 1977, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6777, 3, 19, 20, 1980, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6783, 0, 3, 1947, 6717, 1992, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6801, 0, 3, 1950, 6723, 2001, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6819, 0, 3, 1953, 6729, 2010, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6837, 0, 3, 1956, 6735, 2019, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6855, 0, 3, 1959, 6741, 2028, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6873, 0, 3, 1962, 6747, 2037, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6891, 0, 3, 1965, 6753, 2046, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6909, 0, 3, 1968, 6759, 2055, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6927, 0, 3, 1971, 6765, 2064, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6945, 0, 3, 1974, 6771, 2073, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6963, 0, 3, 1977, 6777, 2082, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6981, 0, 3, 1983, 6783, 67, 73, 2109,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7017, 0, 3, 1992, 6801, 73, 79, 2127,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7053, 0, 3, 2001, 6819, 79, 85, 2145,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7089, 0, 3, 2010, 6837, 85, 91, 2163,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7125, 0, 3, 2019, 6855, 91, 97, 2181,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7161, 0, 3, 2028, 6873, 97, 103, 2199,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7197, 0, 3, 2037, 6891, 103, 109, 2217,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7233, 0, 3, 2046, 6909, 109, 115, 2235,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7269, 0, 3, 2055, 6927, 115, 121, 2253,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7305, 0, 3, 2064, 6945, 121, 127, 2271,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7341, 0, 3, 2073, 6963, 127, 133, 2289,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7377, 0, 3, 2091, 6981, 145, 155, 2307,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7437, 0, 3, 2109, 7017, 155, 165, 2337,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7497, 0, 3, 2127, 7053, 165, 175, 2367,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7557, 0, 3, 2145, 7089, 175, 185, 2397,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7617, 0, 3, 2163, 7125, 185, 195, 2427,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7677, 0, 3, 2181, 7161, 195, 205, 2457,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7737, 0, 3, 2199, 7197, 205, 215, 2487,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7797, 0, 3, 2217, 7233, 215, 225, 2517,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7857, 0, 3, 2235, 7269, 225, 235, 2547,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7917, 0, 3, 2253, 7305, 235, 245, 2577,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7977, 0, 3, 2271, 7341, 245, 255, 2607,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8037, 0, 3, 6981, 7017, 2337, 7497, 275,
                                                 290, 2682, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8127, 0, 3, 7017, 7053, 2367, 7557, 290,
                                                 305, 2727, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8217, 0, 3, 7053, 7089, 2397, 7617, 305,
                                                 320, 2772, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8307, 0, 3, 7089, 7125, 2427, 7677, 320,
                                                 335, 2817, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8397, 0, 3, 7125, 7161, 2457, 7737, 335,
                                                 350, 2862, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8487, 0, 3, 7161, 7197, 2487, 7797, 350,
                                                 365, 2907, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8577, 0, 3, 7197, 7233, 2517, 7857, 365,
                                                 380, 2952, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8667, 0, 3, 7233, 7269, 2547, 7917, 380,
                                                 395, 2997, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8757, 0, 3, 7269, 7305, 2577, 7977, 395,
                                                 410, 3042, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8847, 0, 3, 7377, 7437, 2637, 8037, 440,
                                                 461, 3087, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8973, 0, 3, 7437, 7497, 2682, 8127, 461,
                                                 482, 3150, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9099, 0, 3, 7497, 7557, 2727, 8217, 482,
                                                 503, 3213, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9225, 0, 3, 7557, 7617, 2772, 8307, 503,
                                                 524, 3276, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9351, 0, 3, 7617, 7677, 2817, 8397, 524,
                                                 545, 3339, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9477, 0, 3, 7677, 7737, 2862, 8487, 545,
                                                 566, 3402, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9603, 0, 3, 7737, 7797, 2907, 8577, 566,
                                                 587, 3465, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9729, 0, 3, 7797, 7857, 2952, 8667, 587,
                                                 608, 3528, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9855, 0, 3, 7857, 7917, 2997, 8757, 608,
                                                 629, 3591, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9981, 0, 3, 8037, 8127, 3150, 9099, 671,
                                                 699, 3738, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10149, 0, 3, 8127, 8217, 3213, 9225,
                                                 699, 727, 3822, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10317, 0, 3, 8217, 8307, 3276, 9351,
                                                 727, 755, 3906, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10485, 0, 3, 8307, 8397, 3339, 9477,
                                                 755, 783, 3990, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10653, 0, 3, 8397, 8487, 3402, 9603,
                                                 783, 811, 4074, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10821, 0, 3, 8487, 8577, 3465, 9729,
                                                 811, 839, 4158, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10989, 0, 3, 8577, 8667, 3528, 9855,
                                                 839, 867, 4242, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11157, 0, 3, 8847, 8973, 3654, 9981,
                                                 923, 959, 4326, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11373, 0, 3, 8973, 9099, 3738, 10149,
                                                 959, 995, 4434, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11589, 0, 3, 9099, 9225, 3822, 10317,
                                                 995, 1031, 4542, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11805, 0, 3, 9225, 9351, 3906, 10485,
                                                 1031, 1067, 4650, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12021, 0, 3, 9351, 9477, 3990, 10653,
                                                 1067, 1103, 4758, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12237, 0, 3, 9477, 9603, 4074, 10821,
                                                 1103, 1139, 4866, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12453, 0, 3, 9603, 9729, 4158, 10989,
                                                 1139, 1175, 4974, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12669, 0, 3, 9981, 10149, 4434, 11589,
                                                 1247, 1292, 5217, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12939, 0, 3, 10149, 10317, 4542, 11805,
                                                 1292, 1337, 5352, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13209, 0, 3, 10317, 10485, 4650, 12021,
                                                 1337, 1382, 5487, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13479, 0, 3, 10485, 10653, 4758, 12237,
                                                 1382, 1427, 5622, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13749, 0, 3, 10653, 10821, 4866, 12453,
                                                 1427, 1472, 5757, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 14019, 0, 3, 11157, 11373, 5082, 12669,
                                                 1562, 1617, 5892, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 14349, 0, 3, 11373, 11589, 5217, 12939,
                                                 1617, 1672, 6057, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 14679, 0, 3, 11589, 11805, 5352, 13209,
                                                 1672, 1727, 6222, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 15009, 0, 3, 11805, 12021, 5487, 13479,
                                                 1727, 1782, 6387, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 15339, 0, 3, 12021, 12237, 5622, 13749,
                                                 1782, 1837, 6552, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15669, 3, 1947, 1950, 6723, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15679, 3, 1950, 1953, 6729, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15689, 3, 1953, 1956, 6735, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15699, 3, 1956, 1959, 6741, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15709, 3, 1959, 1962, 6747, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15719, 3, 1962, 1965, 6753, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15729, 3, 1965, 1968, 6759, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15739, 3, 1968, 1971, 6765, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15749, 3, 1971, 1974, 6771, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15759, 3, 1974, 1977, 6777, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 15769, 0, 3, 6717, 15669, 6801, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15799, 0, 3, 6723, 15679, 6819, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15829, 0, 3, 6729, 15689, 6837, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15859, 0, 3, 6735, 15699, 6855, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15889, 0, 3, 6741, 15709, 6873, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15919, 0, 3, 6747, 15719, 6891, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15949, 0, 3, 6753, 15729, 6909, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15979, 0, 3, 6759, 15739, 6927, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16009, 0, 3, 6765, 15749, 6945, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16039, 0, 3, 6771, 15759, 6963, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 16069, 0, 3, 6783, 15769, 2091, 2109,
                                                 7017, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16129, 0, 3, 6801, 15799, 2109, 2127,
                                                 7053, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16189, 0, 3, 6819, 15829, 2127, 2145,
                                                 7089, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16249, 0, 3, 6837, 15859, 2145, 2163,
                                                 7125, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16309, 0, 3, 6855, 15889, 2163, 2181,
                                                 7161, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16369, 0, 3, 6873, 15919, 2181, 2199,
                                                 7197, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16429, 0, 3, 6891, 15949, 2199, 2217,
                                                 7233, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16489, 0, 3, 6909, 15979, 2217, 2235,
                                                 7269, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16549, 0, 3, 6927, 16009, 2235, 2253,
                                                 7305, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16609, 0, 3, 6945, 16039, 2253, 2271,
                                                 7341, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16669, 0, 3, 7017, 16129, 2307, 2337,
                                                 7497, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16769, 0, 3, 7053, 16189, 2337, 2367,
                                                 7557, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16869, 0, 3, 7089, 16249, 2367, 2397,
                                                 7617, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16969, 0, 3, 7125, 16309, 2397, 2427,
                                                 7677, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17069, 0, 3, 7161, 16369, 2427, 2457,
                                                 7737, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17169, 0, 3, 7197, 16429, 2457, 2487,
                                                 7797, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17269, 0, 3, 7233, 16489, 2487, 2517,
                                                 7857, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17369, 0, 3, 7269, 16549, 2517, 2547,
                                                 7917, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17469, 0, 3, 7305, 16609, 2547, 2577,
                                                 7977, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17569, 0, 3, 16069, 16129, 7497, 16769,
                                                 2637, 2682, 8127, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17719, 0, 3, 16129, 16189, 7557, 16869,
                                                 2682, 2727, 8217, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17869, 0, 3, 16189, 16249, 7617, 16969,
                                                 2727, 2772, 8307, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18019, 0, 3, 16249, 16309, 7677, 17069,
                                                 2772, 2817, 8397, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18169, 0, 3, 16309, 16369, 7737, 17169,
                                                 2817, 2862, 8487, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18319, 0, 3, 16369, 16429, 7797, 17269,
                                                 2862, 2907, 8577, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18469, 0, 3, 16429, 16489, 7857, 17369,
                                                 2907, 2952, 8667, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18619, 0, 3, 16489, 16549, 7917, 17469,
                                                 2952, 2997, 8757, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18769, 0, 3, 16669, 16769, 8127, 17719,
                                                 3087, 3150, 9099, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18979, 0, 3, 16769, 16869, 8217, 17869,
                                                 3150, 3213, 9225, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19189, 0, 3, 16869, 16969, 8307, 18019,
                                                 3213, 3276, 9351, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19399, 0, 3, 16969, 17069, 8397, 18169,
                                                 3276, 3339, 9477, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19609, 0, 3, 17069, 17169, 8487, 18319,
                                                 3339, 3402, 9603, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19819, 0, 3, 17169, 17269, 8577, 18469,
                                                 3402, 3465, 9729, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20029, 0, 3, 17269, 17369, 8667, 18619,
                                                 3465, 3528, 9855, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20239, 0, 3, 17569, 17719, 9099, 18979,
                                                 3654, 3738, 10149, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20519, 0, 3, 17719, 17869, 9225, 19189,
                                                 3738, 3822, 10317, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20799, 0, 3, 17869, 18019, 9351, 19399,
                                                 3822, 3906, 10485, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21079, 0, 3, 18019, 18169, 9477, 19609,
                                                 3906, 3990, 10653, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21359, 0, 3, 18169, 18319, 9603, 19819,
                                                 3990, 4074, 10821, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21639, 0, 3, 18319, 18469, 9729, 20029,
                                                 4074, 4158, 10989, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 21919, 0, 3, 18769, 18979, 10149, 20519,
                                                 4326, 4434, 11589, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22279, 0, 3, 18979, 19189, 10317, 20799,
                                                 4434, 4542, 11805, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22639, 0, 3, 19189, 19399, 10485, 21079,
                                                 4542, 4650, 12021, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22999, 0, 3, 19399, 19609, 10653, 21359,
                                                 4650, 4758, 12237, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 23359, 0, 3, 19609, 19819, 10821, 21639,
                                                 4758, 4866, 12453, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 23719, 0, 3, 20239, 20519, 11589, 22279,
                                                 5082, 5217, 12939, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 24169, 0, 3, 20519, 20799, 11805, 22639,
                                                 5217, 5352, 13209, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 24619, 0, 3, 20799, 21079, 12021, 22999,
                                                 5352, 5487, 13479, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 25069, 0, 3, 21079, 21359, 12237, 23359,
                                                 5487, 5622, 13749, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 25519, 0, 3, 21919, 22279, 12939, 24169,
                                                 5892, 6057, 14679, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 26069, 0, 3, 22279, 22639, 13209, 24619,
                                                 6057, 6222, 15009, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 26619, 0, 3, 22639, 22999, 13479, 25069,
                                                 6222, 6387, 15339, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27169, 3, 6717, 6723, 15679, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27184, 3, 6723, 6729, 15689, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27199, 3, 6729, 6735, 15699, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27214, 3, 6735, 6741, 15709, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27229, 3, 6741, 6747, 15719, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27244, 3, 6747, 6753, 15729, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27259, 3, 6753, 6759, 15739, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27274, 3, 6759, 6765, 15749, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27289, 3, 6765, 6771, 15759, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 27304, 0, 3, 15669, 27169, 15799, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27349, 0, 3, 15679, 27184, 15829, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27394, 0, 3, 15689, 27199, 15859, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27439, 0, 3, 15699, 27214, 15889, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27484, 0, 3, 15709, 27229, 15919, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27529, 0, 3, 15719, 27244, 15949, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27574, 0, 3, 15729, 27259, 15979, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27619, 0, 3, 15739, 27274, 16009, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27664, 0, 3, 15749, 27289, 16039, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 27709, 0, 3, 15769, 27304, 6981, 7017,
                                                 16129, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27799, 0, 3, 15799, 27349, 7017, 7053,
                                                 16189, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27889, 0, 3, 15829, 27394, 7053, 7089,
                                                 16249, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27979, 0, 3, 15859, 27439, 7089, 7125,
                                                 16309, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28069, 0, 3, 15889, 27484, 7125, 7161,
                                                 16369, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28159, 0, 3, 15919, 27529, 7161, 7197,
                                                 16429, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28249, 0, 3, 15949, 27574, 7197, 7233,
                                                 16489, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28339, 0, 3, 15979, 27619, 7233, 7269,
                                                 16549, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28429, 0, 3, 16009, 27664, 7269, 7305,
                                                 16609, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28519, 0, 3, 16069, 27709, 7377, 7437,
                                                 16669, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28669, 0, 3, 16129, 27799, 7437, 7497,
                                                 16769, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28819, 0, 3, 16189, 27889, 7497, 7557,
                                                 16869, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28969, 0, 3, 16249, 27979, 7557, 7617,
                                                 16969, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29119, 0, 3, 16309, 28069, 7617, 7677,
                                                 17069, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29269, 0, 3, 16369, 28159, 7677, 7737,
                                                 17169, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29419, 0, 3, 16429, 28249, 7737, 7797,
                                                 17269, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29569, 0, 3, 16489, 28339, 7797, 7857,
                                                 17369, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29719, 0, 3, 16549, 28429, 7857, 7917,
                                                 17469, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29869, 0, 3, 27709, 27799, 16769, 28819,
                                                 8037, 8127, 17719, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30094, 0, 3, 27799, 27889, 16869, 28969,
                                                 8127, 8217, 17869, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30319, 0, 3, 27889, 27979, 16969, 29119,
                                                 8217, 8307, 18019, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30544, 0, 3, 27979, 28069, 17069, 29269,
                                                 8307, 8397, 18169, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30769, 0, 3, 28069, 28159, 17169, 29419,
                                                 8397, 8487, 18319, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30994, 0, 3, 28159, 28249, 17269, 29569,
                                                 8487, 8577, 18469, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31219, 0, 3, 28249, 28339, 17369, 29719,
                                                 8577, 8667, 18619, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31444, 0, 3, 28519, 28669, 17569, 29869,
                                                 8847, 8973, 18769, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31759, 0, 3, 28669, 28819, 17719, 30094,
                                                 8973, 9099, 18979, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32074, 0, 3, 28819, 28969, 17869, 30319,
                                                 9099, 9225, 19189, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32389, 0, 3, 28969, 29119, 18019, 30544,
                                                 9225, 9351, 19399, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32704, 0, 3, 29119, 29269, 18169, 30769,
                                                 9351, 9477, 19609, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33019, 0, 3, 29269, 29419, 18319, 30994,
                                                 9477, 9603, 19819, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33334, 0, 3, 29419, 29569, 18469, 31219,
                                                 9603, 9729, 20029, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 33649, 0, 3, 29869, 30094, 18979, 32074,
                                                 9981, 10149, 20519, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34069, 0, 3, 30094, 30319, 19189, 32389,
                                                 10149, 10317, 20799, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34489, 0, 3, 30319, 30544, 19399, 32704,
                                                 10317, 10485, 21079, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34909, 0, 3, 30544, 30769, 19609, 33019,
                                                 10485, 10653, 21359, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35329, 0, 3, 30769, 30994, 19819, 33334,
                                                 10653, 10821, 21639, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 35749, 0, 3, 31444, 31759, 20239, 33649,
                                                 11157, 11373, 21919, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36289, 0, 3, 31759, 32074, 20519, 34069,
                                                 11373, 11589, 22279, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36829, 0, 3, 32074, 32389, 20799, 34489,
                                                 11589, 11805, 22639, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 37369, 0, 3, 32389, 32704, 21079, 34909,
                                                 11805, 12021, 22999, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 37909, 0, 3, 32704, 33019, 21359, 35329,
                                                 12021, 12237, 23359, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 38449, 0, 3, 33649, 34069, 22279, 36829,
                                                 12669, 12939, 24169, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 39124, 0, 3, 34069, 34489, 22639, 37369,
                                                 12939, 13209, 24619, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 39799, 0, 3, 34489, 34909, 22999, 37909,
                                                 13209, 13479, 25069, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 40474, 0, 3, 35749, 36289, 23719, 38449,
                                                 14019, 14349, 25519, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 41299, 0, 3, 36289, 36829, 24169, 39124,
                                                 14349, 14679, 26069, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 42124, 0, 3, 36829, 37369, 24619, 39799,
                                                 14679, 15009, 26619, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42949, 3, 15669, 15679, 27184, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42970, 3, 15679, 15689, 27199, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42991, 3, 15689, 15699, 27214, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 43012, 3, 15699, 15709, 27229, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 43033, 3, 15709, 15719, 27244, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 43054, 3, 15719, 15729, 27259, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 43075, 3, 15729, 15739, 27274, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 43096, 3, 15739, 15749, 27289, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 43117, 0, 3, 27169, 42949, 27349, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43180, 0, 3, 27184, 42970, 27394, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43243, 0, 3, 27199, 42991, 27439, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43306, 0, 3, 27214, 43012, 27484, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43369, 0, 3, 27229, 43033, 27529, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43432, 0, 3, 27244, 43054, 27574, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43495, 0, 3, 27259, 43075, 27619, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43558, 0, 3, 27274, 43096, 27664, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 43621, 0, 3, 27304, 43117, 16069, 16129,
                                                 27799, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43747, 0, 3, 27349, 43180, 16129, 16189,
                                                 27889, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43873, 0, 3, 27394, 43243, 16189, 16249,
                                                 27979, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43999, 0, 3, 27439, 43306, 16249, 16309,
                                                 28069, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44125, 0, 3, 27484, 43369, 16309, 16369,
                                                 28159, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44251, 0, 3, 27529, 43432, 16369, 16429,
                                                 28249, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44377, 0, 3, 27574, 43495, 16429, 16489,
                                                 28339, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44503, 0, 3, 27619, 43558, 16489, 16549,
                                                 28429, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44629, 0, 3, 27799, 43747, 16669, 16769,
                                                 28819, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44839, 0, 3, 27889, 43873, 16769, 16869,
                                                 28969, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45049, 0, 3, 27979, 43999, 16869, 16969,
                                                 29119, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45259, 0, 3, 28069, 44125, 16969, 17069,
                                                 29269, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45469, 0, 3, 28159, 44251, 17069, 17169,
                                                 29419, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45679, 0, 3, 28249, 44377, 17169, 17269,
                                                 29569, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45889, 0, 3, 28339, 44503, 17269, 17369,
                                                 29719, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46099, 0, 3, 43621, 43747, 28819, 44839,
                                                 17569, 17719, 30094, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46414, 0, 3, 43747, 43873, 28969, 45049,
                                                 17719, 17869, 30319, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46729, 0, 3, 43873, 43999, 29119, 45259,
                                                 17869, 18019, 30544, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47044, 0, 3, 43999, 44125, 29269, 45469,
                                                 18019, 18169, 30769, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47359, 0, 3, 44125, 44251, 29419, 45679,
                                                 18169, 18319, 30994, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47674, 0, 3, 44251, 44377, 29569, 45889,
                                                 18319, 18469, 31219, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 47989, 0, 3, 44629, 44839, 30094, 46414,
                                                 18769, 18979, 32074, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 48430, 0, 3, 44839, 45049, 30319, 46729,
                                                 18979, 19189, 32389, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 48871, 0, 3, 45049, 45259, 30544, 47044,
                                                 19189, 19399, 32704, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49312, 0, 3, 45259, 45469, 30769, 47359,
                                                 19399, 19609, 33019, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49753, 0, 3, 45469, 45679, 30994, 47674,
                                                 19609, 19819, 33334, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 50194, 0, 3, 46099, 46414, 32074, 48430,
                                                 20239, 20519, 34069, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 50782, 0, 3, 46414, 46729, 32389, 48871,
                                                 20519, 20799, 34489, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 51370, 0, 3, 46729, 47044, 32704, 49312,
                                                 20799, 21079, 34909, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 51958, 0, 3, 47044, 47359, 33019, 49753,
                                                 21079, 21359, 35329, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 52546, 0, 3, 47989, 48430, 34069, 50782,
                                                 21919, 22279, 36829, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 53302, 0, 3, 48430, 48871, 34489, 51370,
                                                 22279, 22639, 37369, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 54058, 0, 3, 48871, 49312, 34909, 51958,
                                                 22639, 22999, 37909, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 54814, 0, 3, 50194, 50782, 36829, 53302,
                                                 23719, 24169, 39124, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 55759, 0, 3, 50782, 51370, 37369, 54058,
                                                 24169, 24619, 39799, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 56704, 0, 3, 52546, 53302, 39124, 55759,
                                                 25519, 26069, 42124, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57859, 3, 27169, 27184, 42970, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57887, 3, 27184, 27199, 42991, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57915, 3, 27199, 27214, 43012, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57943, 3, 27214, 27229, 43033, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57971, 3, 27229, 27244, 43054, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 57999, 3, 27244, 27259, 43075, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 58027, 3, 27259, 27274, 43096, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 58055, 0, 3, 42949, 57859, 43180, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58139, 0, 3, 42970, 57887, 43243, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58223, 0, 3, 42991, 57915, 43306, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58307, 0, 3, 43012, 57943, 43369, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58391, 0, 3, 43033, 57971, 43432, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58475, 0, 3, 43054, 57999, 43495, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 58559, 0, 3, 43075, 58027, 43558, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 58643, 0, 3, 43117, 58055, 27709, 27799,
                                                 43747, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58811, 0, 3, 43180, 58139, 27799, 27889,
                                                 43873, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58979, 0, 3, 43243, 58223, 27889, 27979,
                                                 43999, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59147, 0, 3, 43306, 58307, 27979, 28069,
                                                 44125, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59315, 0, 3, 43369, 58391, 28069, 28159,
                                                 44251, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59483, 0, 3, 43432, 58475, 28159, 28249,
                                                 44377, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59651, 0, 3, 43495, 58559, 28249, 28339,
                                                 44503, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 59819, 0, 3, 43621, 58643, 28519, 28669,
                                                 44629, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60099, 0, 3, 43747, 58811, 28669, 28819,
                                                 44839, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60379, 0, 3, 43873, 58979, 28819, 28969,
                                                 45049, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60659, 0, 3, 43999, 59147, 28969, 29119,
                                                 45259, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60939, 0, 3, 44125, 59315, 29119, 29269,
                                                 45469, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61219, 0, 3, 44251, 59483, 29269, 29419,
                                                 45679, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61499, 0, 3, 44377, 59651, 29419, 29569,
                                                 45889, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 61779, 0, 3, 58643, 58811, 44839, 60379,
                                                 29869, 30094, 46414, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 62199, 0, 3, 58811, 58979, 45049, 60659,
                                                 30094, 30319, 46729, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 62619, 0, 3, 58979, 59147, 45259, 60939,
                                                 30319, 30544, 47044, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63039, 0, 3, 59147, 59315, 45469, 61219,
                                                 30544, 30769, 47359, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63459, 0, 3, 59315, 59483, 45679, 61499,
                                                 30769, 30994, 47674, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 63879, 0, 3, 59819, 60099, 46099, 61779,
                                                 31444, 31759, 47989, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 64467, 0, 3, 60099, 60379, 46414, 62199,
                                                 31759, 32074, 48430, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 65055, 0, 3, 60379, 60659, 46729, 62619,
                                                 32074, 32389, 48871, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 65643, 0, 3, 60659, 60939, 47044, 63039,
                                                 32389, 32704, 49312, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 66231, 0, 3, 60939, 61219, 47359, 63459,
                                                 32704, 33019, 49753, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 66819, 0, 3, 61779, 62199, 48430, 65055,
                                                 33649, 34069, 50782, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 67603, 0, 3, 62199, 62619, 48871, 65643,
                                                 34069, 34489, 51370, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 68387, 0, 3, 62619, 63039, 49312, 66231,
                                                 34489, 34909, 51958, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 69171, 0, 3, 63879, 64467, 50194, 66819,
                                                 35749, 36289, 52546, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 70179, 0, 3, 64467, 65055, 50782, 67603,
                                                 36289, 36829, 53302, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 71187, 0, 3, 65055, 65643, 51370, 68387,
                                                 36829, 37369, 54058, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 72195, 0, 3, 66819, 67603, 53302, 71187,
                                                 38449, 39124, 55759, ncols, alpha, beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 73455, 0, 3, 69171, 70179, 54814, 72195,
                                                 40474, 41299, 56704, ncols, alpha, beta, p);

            compute_prim_geom_10_li_electron_repulsion_0(buffer, 74995, 69171, 73455, ncols,
                                                         alpha);

            compute_prim_geom_10_li_electron_repulsion_1(buffer, 76255, 69171, 73455, ncols,
                                                         alpha);

            compute_prim_geom_10_li_electron_repulsion_2(buffer, 77515, 69171, 73455, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 78775, 74995, 3780, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 82555, 78775, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 82555, 13, nmax);

    simdtrf::transform_i_inner(buffer, 82555, 80035, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 221 * nvalues, nvalues, buffer, 82555, 13, nmax);

    simdtrf::transform_i_inner(buffer, 82555, 81295, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 442 * nvalues, nvalues, buffer, 82555, 13, nmax);
}

}  // namespace simdt2ceri

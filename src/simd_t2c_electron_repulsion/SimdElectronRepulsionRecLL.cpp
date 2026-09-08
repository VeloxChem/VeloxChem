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


#include "SimdElectronRepulsionRecLL.hpp"

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
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFL.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGL.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHL.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIL.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKL.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLL.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ll_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ll_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 120464, 117674, 2025, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 16, ncols, fj, mu);

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

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 105, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 111, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 117, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 123, 0, 16, 17, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 129, 0, 17, 18, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 135, 0, 18, 19, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 141, 0, 19, 20, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 147, 0, 20, 21, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 153, 0, 21, 22, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 24, 27, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 27, 30, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 30, 33, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 33, 36, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 36, 39, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 209, 0, 39, 42, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 219, 0, 42, 45, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 229, 0, 45, 48, 123, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 239, 0, 48, 51, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 249, 0, 51, 54, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 259, 0, 54, 57, 141, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 269, 0, 57, 60, 147, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 279, 0, 60, 63, 153, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 289, 0, 69, 75, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 304, 0, 75, 81, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 319, 0, 81, 87, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 334, 0, 87, 93, 189, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 349, 0, 93, 99, 199, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 364, 0, 99, 105, 209, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 379, 0, 105, 111, 219, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 394, 0, 111, 117, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 409, 0, 117, 123, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 424, 0, 123, 129, 249, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 439, 0, 129, 135, 259, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 454, 0, 135, 141, 269, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 469, 0, 141, 147, 279, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 484, 0, 159, 169, 319, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 505, 0, 169, 179, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 526, 0, 179, 189, 349, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 547, 0, 189, 199, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 568, 0, 199, 209, 379, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 589, 0, 209, 219, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 610, 0, 219, 229, 409, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 631, 0, 229, 239, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 652, 0, 239, 249, 439, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 673, 0, 249, 259, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 694, 0, 259, 269, 469, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 715, 0, 289, 304, 484, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 743, 0, 304, 319, 505, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 771, 0, 319, 334, 526, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 799, 0, 334, 349, 547, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 827, 0, 349, 364, 568, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 855, 0, 364, 379, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 883, 0, 379, 394, 610, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 911, 0, 394, 409, 631, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 939, 0, 409, 424, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 967, 0, 424, 439, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 995, 0, 439, 454, 694, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1023, 0, 484, 505, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1059, 0, 505, 526, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1095, 0, 526, 547, 827, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1131, 0, 547, 568, 855, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1167, 0, 568, 589, 883, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1203, 0, 589, 610, 911, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1239, 0, 610, 631, 939, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1275, 0, 631, 652, 967, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1311, 0, 652, 673, 995, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1347, 0, 715, 743, 1023, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1392, 0, 743, 771, 1059, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1437, 0, 771, 799, 1095, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1482, 0, 799, 827, 1131, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1527, 0, 827, 855, 1167, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1572, 0, 855, 883, 1203, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1617, 0, 883, 911, 1239, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1662, 0, 911, 939, 1275, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1707, 0, 939, 967, 1311, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1752, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1755, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1758, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1761, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1764, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1767, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1770, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1773, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1776, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1779, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1782, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1785, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1788, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1791, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1794, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1803, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1812, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1821, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1830, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1839, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1848, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1857, 3, 16, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1866, 3, 17, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1875, 3, 18, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1884, 3, 19, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1893, 3, 20, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1902, 3, 21, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1911, 3, 22, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1920, 0, 3, 27, 1803, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1938, 0, 3, 30, 1812, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1956, 0, 3, 33, 1821, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1974, 0, 3, 36, 1830, 99, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1992, 0, 3, 39, 1839, 105, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2010, 0, 3, 42, 1848, 111, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2028, 0, 3, 45, 1857, 117, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2046, 0, 3, 48, 1866, 123, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2064, 0, 3, 51, 1875, 129, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2082, 0, 3, 54, 1884, 135, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2100, 0, 3, 57, 1893, 141, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2118, 0, 3, 60, 1902, 147, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2136, 0, 3, 63, 1911, 153, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2154, 0, 3, 81, 1938, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2184, 0, 3, 87, 1956, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2214, 0, 3, 93, 1974, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2244, 0, 3, 99, 1992, 199, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2274, 0, 3, 105, 2010, 209, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2304, 0, 3, 111, 2028, 219, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2334, 0, 3, 117, 2046, 229, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2364, 0, 3, 123, 2064, 239, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2394, 0, 3, 129, 2082, 249, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2424, 0, 3, 135, 2100, 259, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2454, 0, 3, 141, 2118, 269, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2484, 0, 3, 147, 2136, 279, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2514, 0, 3, 169, 2184, 319, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2559, 0, 3, 179, 2214, 334, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2604, 0, 3, 189, 2244, 349, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2649, 0, 3, 199, 2274, 364, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2694, 0, 3, 209, 2304, 379, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2739, 0, 3, 219, 2334, 394, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2784, 0, 3, 229, 2364, 409, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2829, 0, 3, 239, 2394, 424, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2874, 0, 3, 249, 2424, 439, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2919, 0, 3, 259, 2454, 454, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2964, 0, 3, 269, 2484, 469, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3009, 0, 3, 319, 2559, 505, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3072, 0, 3, 334, 2604, 526, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3135, 0, 3, 349, 2649, 547, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3198, 0, 3, 364, 2694, 568, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3261, 0, 3, 379, 2739, 589, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3324, 0, 3, 394, 2784, 610, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3387, 0, 3, 409, 2829, 631, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3450, 0, 3, 424, 2874, 652, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3513, 0, 3, 439, 2919, 673, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3576, 0, 3, 454, 2964, 694, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3639, 0, 3, 505, 3072, 771, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3723, 0, 3, 526, 3135, 799, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3807, 0, 3, 547, 3198, 827, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3891, 0, 3, 568, 3261, 855, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3975, 0, 3, 589, 3324, 883, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4059, 0, 3, 610, 3387, 911, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4143, 0, 3, 631, 3450, 939, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4227, 0, 3, 652, 3513, 967, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4311, 0, 3, 673, 3576, 995, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4395, 0, 3, 771, 3723, 1059, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4503, 0, 3, 799, 3807, 1095, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4611, 0, 3, 827, 3891, 1131, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4719, 0, 3, 855, 3975, 1167, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4827, 0, 3, 883, 4059, 1203, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4935, 0, 3, 911, 4143, 1239, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5043, 0, 3, 939, 4227, 1275, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5151, 0, 3, 967, 4311, 1311, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5259, 0, 3, 1059, 4503, 1437, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5394, 0, 3, 1095, 4611, 1482, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5529, 0, 3, 1131, 4719, 1527, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5664, 0, 3, 1167, 4827, 1572, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5799, 0, 3, 1203, 4935, 1617, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5934, 0, 3, 1239, 5043, 1662, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6069, 0, 3, 1275, 5151, 1707, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6204, 3, 9, 10, 1755, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6210, 3, 10, 11, 1758, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6216, 3, 11, 12, 1761, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6222, 3, 12, 13, 1764, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6228, 3, 13, 14, 1767, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6234, 3, 14, 15, 1770, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6240, 3, 15, 16, 1773, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6246, 3, 16, 17, 1776, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6252, 3, 17, 18, 1779, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6258, 3, 18, 19, 1782, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6264, 3, 19, 20, 1785, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6270, 3, 20, 21, 1788, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6276, 3, 21, 22, 1791, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6282, 0, 3, 1752, 6204, 1803, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6300, 0, 3, 1755, 6210, 1812, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6318, 0, 3, 1758, 6216, 1821, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6336, 0, 3, 1761, 6222, 1830, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6354, 0, 3, 1764, 6228, 1839, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6372, 0, 3, 1767, 6234, 1848, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6390, 0, 3, 1770, 6240, 1857, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6408, 0, 3, 1773, 6246, 1866, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6426, 0, 3, 1776, 6252, 1875, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6444, 0, 3, 1779, 6258, 1884, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6462, 0, 3, 1782, 6264, 1893, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6480, 0, 3, 1785, 6270, 1902, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6498, 0, 3, 1788, 6276, 1911, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6516, 0, 3, 1794, 6282, 69, 75, 1920,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6552, 0, 3, 1803, 6300, 75, 81, 1938,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6588, 0, 3, 1812, 6318, 81, 87, 1956,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6624, 0, 3, 1821, 6336, 87, 93, 1974,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6660, 0, 3, 1830, 6354, 93, 99, 1992,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6696, 0, 3, 1839, 6372, 99, 105, 2010,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6732, 0, 3, 1848, 6390, 105, 111, 2028,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6768, 0, 3, 1857, 6408, 111, 117, 2046,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6804, 0, 3, 1866, 6426, 117, 123, 2064,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6840, 0, 3, 1875, 6444, 123, 129, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6876, 0, 3, 1884, 6462, 129, 135, 2100,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6912, 0, 3, 1893, 6480, 135, 141, 2118,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6948, 0, 3, 1902, 6498, 141, 147, 2136,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6984, 0, 3, 1938, 6588, 159, 169, 2184,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7044, 0, 3, 1956, 6624, 169, 179, 2214,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7104, 0, 3, 1974, 6660, 179, 189, 2244,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7164, 0, 3, 1992, 6696, 189, 199, 2274,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7224, 0, 3, 2010, 6732, 199, 209, 2304,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7284, 0, 3, 2028, 6768, 209, 219, 2334,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7344, 0, 3, 2046, 6804, 219, 229, 2364,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7404, 0, 3, 2064, 6840, 229, 239, 2394,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7464, 0, 3, 2082, 6876, 239, 249, 2424,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7524, 0, 3, 2100, 6912, 249, 259, 2454,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7584, 0, 3, 2118, 6948, 259, 269, 2484,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7644, 0, 3, 6516, 6552, 2154, 6984, 289,
                                                 304, 2514, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7734, 0, 3, 6552, 6588, 2184, 7044, 304,
                                                 319, 2559, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7824, 0, 3, 6588, 6624, 2214, 7104, 319,
                                                 334, 2604, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7914, 0, 3, 6624, 6660, 2244, 7164, 334,
                                                 349, 2649, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8004, 0, 3, 6660, 6696, 2274, 7224, 349,
                                                 364, 2694, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8094, 0, 3, 6696, 6732, 2304, 7284, 364,
                                                 379, 2739, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8184, 0, 3, 6732, 6768, 2334, 7344, 379,
                                                 394, 2784, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8274, 0, 3, 6768, 6804, 2364, 7404, 394,
                                                 409, 2829, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8364, 0, 3, 6804, 6840, 2394, 7464, 409,
                                                 424, 2874, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8454, 0, 3, 6840, 6876, 2424, 7524, 424,
                                                 439, 2919, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8544, 0, 3, 6876, 6912, 2454, 7584, 439,
                                                 454, 2964, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8634, 0, 3, 6984, 7044, 2559, 7824, 484,
                                                 505, 3072, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8760, 0, 3, 7044, 7104, 2604, 7914, 505,
                                                 526, 3135, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8886, 0, 3, 7104, 7164, 2649, 8004, 526,
                                                 547, 3198, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9012, 0, 3, 7164, 7224, 2694, 8094, 547,
                                                 568, 3261, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9138, 0, 3, 7224, 7284, 2739, 8184, 568,
                                                 589, 3324, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9264, 0, 3, 7284, 7344, 2784, 8274, 589,
                                                 610, 3387, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9390, 0, 3, 7344, 7404, 2829, 8364, 610,
                                                 631, 3450, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9516, 0, 3, 7404, 7464, 2874, 8454, 631,
                                                 652, 3513, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9642, 0, 3, 7464, 7524, 2919, 8544, 652,
                                                 673, 3576, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9768, 0, 3, 7644, 7734, 3009, 8634, 715,
                                                 743, 3639, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9936, 0, 3, 7734, 7824, 3072, 8760, 743,
                                                 771, 3723, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10104, 0, 3, 7824, 7914, 3135, 8886,
                                                 771, 799, 3807, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10272, 0, 3, 7914, 8004, 3198, 9012,
                                                 799, 827, 3891, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10440, 0, 3, 8004, 8094, 3261, 9138,
                                                 827, 855, 3975, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10608, 0, 3, 8094, 8184, 3324, 9264,
                                                 855, 883, 4059, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10776, 0, 3, 8184, 8274, 3387, 9390,
                                                 883, 911, 4143, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10944, 0, 3, 8274, 8364, 3450, 9516,
                                                 911, 939, 4227, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11112, 0, 3, 8364, 8454, 3513, 9642,
                                                 939, 967, 4311, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11280, 0, 3, 8634, 8760, 3723, 10104,
                                                 1023, 1059, 4503, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11496, 0, 3, 8760, 8886, 3807, 10272,
                                                 1059, 1095, 4611, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11712, 0, 3, 8886, 9012, 3891, 10440,
                                                 1095, 1131, 4719, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11928, 0, 3, 9012, 9138, 3975, 10608,
                                                 1131, 1167, 4827, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12144, 0, 3, 9138, 9264, 4059, 10776,
                                                 1167, 1203, 4935, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12360, 0, 3, 9264, 9390, 4143, 10944,
                                                 1203, 1239, 5043, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12576, 0, 3, 9390, 9516, 4227, 11112,
                                                 1239, 1275, 5151, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12792, 0, 3, 9768, 9936, 4395, 11280,
                                                 1347, 1392, 5259, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13062, 0, 3, 9936, 10104, 4503, 11496,
                                                 1392, 1437, 5394, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13332, 0, 3, 10104, 10272, 4611, 11712,
                                                 1437, 1482, 5529, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13602, 0, 3, 10272, 10440, 4719, 11928,
                                                 1482, 1527, 5664, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 13872, 0, 3, 10440, 10608, 4827, 12144,
                                                 1527, 1572, 5799, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14142, 0, 3, 10608, 10776, 4935, 12360,
                                                 1572, 1617, 5934, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 14412, 0, 3, 10776, 10944, 5043, 12576,
                                                 1617, 1662, 6069, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14682, 3, 1752, 1755, 6210, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14692, 3, 1755, 1758, 6216, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14702, 3, 1758, 1761, 6222, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14712, 3, 1761, 1764, 6228, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14722, 3, 1764, 1767, 6234, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14732, 3, 1767, 1770, 6240, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14742, 3, 1770, 1773, 6246, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14752, 3, 1773, 1776, 6252, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14762, 3, 1776, 1779, 6258, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14772, 3, 1779, 1782, 6264, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14782, 3, 1782, 1785, 6270, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14792, 3, 1785, 1788, 6276, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 14802, 0, 3, 6204, 14682, 6300, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14832, 0, 3, 6210, 14692, 6318, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14862, 0, 3, 6216, 14702, 6336, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14892, 0, 3, 6222, 14712, 6354, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14922, 0, 3, 6228, 14722, 6372, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14952, 0, 3, 6234, 14732, 6390, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14982, 0, 3, 6240, 14742, 6408, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15012, 0, 3, 6246, 14752, 6426, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15042, 0, 3, 6252, 14762, 6444, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15072, 0, 3, 6258, 14772, 6462, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15102, 0, 3, 6264, 14782, 6480, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 15132, 0, 3, 6270, 14792, 6498, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 15162, 0, 3, 6300, 14832, 1920, 1938,
                                                 6588, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15222, 0, 3, 6318, 14862, 1938, 1956,
                                                 6624, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15282, 0, 3, 6336, 14892, 1956, 1974,
                                                 6660, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15342, 0, 3, 6354, 14922, 1974, 1992,
                                                 6696, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15402, 0, 3, 6372, 14952, 1992, 2010,
                                                 6732, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15462, 0, 3, 6390, 14982, 2010, 2028,
                                                 6768, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15522, 0, 3, 6408, 15012, 2028, 2046,
                                                 6804, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15582, 0, 3, 6426, 15042, 2046, 2064,
                                                 6840, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15642, 0, 3, 6444, 15072, 2064, 2082,
                                                 6876, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15702, 0, 3, 6462, 15102, 2082, 2100,
                                                 6912, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15762, 0, 3, 6480, 15132, 2100, 2118,
                                                 6948, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15822, 0, 3, 6588, 15222, 2154, 2184,
                                                 7044, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15922, 0, 3, 6624, 15282, 2184, 2214,
                                                 7104, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16022, 0, 3, 6660, 15342, 2214, 2244,
                                                 7164, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16122, 0, 3, 6696, 15402, 2244, 2274,
                                                 7224, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16222, 0, 3, 6732, 15462, 2274, 2304,
                                                 7284, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16322, 0, 3, 6768, 15522, 2304, 2334,
                                                 7344, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16422, 0, 3, 6804, 15582, 2334, 2364,
                                                 7404, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16522, 0, 3, 6840, 15642, 2364, 2394,
                                                 7464, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16622, 0, 3, 6876, 15702, 2394, 2424,
                                                 7524, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16722, 0, 3, 6912, 15762, 2424, 2454,
                                                 7584, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16822, 0, 3, 15162, 15222, 7044, 15922,
                                                 2514, 2559, 7824, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 16972, 0, 3, 15222, 15282, 7104, 16022,
                                                 2559, 2604, 7914, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17122, 0, 3, 15282, 15342, 7164, 16122,
                                                 2604, 2649, 8004, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17272, 0, 3, 15342, 15402, 7224, 16222,
                                                 2649, 2694, 8094, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17422, 0, 3, 15402, 15462, 7284, 16322,
                                                 2694, 2739, 8184, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17572, 0, 3, 15462, 15522, 7344, 16422,
                                                 2739, 2784, 8274, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17722, 0, 3, 15522, 15582, 7404, 16522,
                                                 2784, 2829, 8364, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17872, 0, 3, 15582, 15642, 7464, 16622,
                                                 2829, 2874, 8454, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18022, 0, 3, 15642, 15702, 7524, 16722,
                                                 2874, 2919, 8544, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18172, 0, 3, 15822, 15922, 7824, 16972,
                                                 3009, 3072, 8760, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18382, 0, 3, 15922, 16022, 7914, 17122,
                                                 3072, 3135, 8886, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18592, 0, 3, 16022, 16122, 8004, 17272,
                                                 3135, 3198, 9012, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 18802, 0, 3, 16122, 16222, 8094, 17422,
                                                 3198, 3261, 9138, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19012, 0, 3, 16222, 16322, 8184, 17572,
                                                 3261, 3324, 9264, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19222, 0, 3, 16322, 16422, 8274, 17722,
                                                 3324, 3387, 9390, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19432, 0, 3, 16422, 16522, 8364, 17872,
                                                 3387, 3450, 9516, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19642, 0, 3, 16522, 16622, 8454, 18022,
                                                 3450, 3513, 9642, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19852, 0, 3, 16822, 16972, 8760, 18382,
                                                 3639, 3723, 10104, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20132, 0, 3, 16972, 17122, 8886, 18592,
                                                 3723, 3807, 10272, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20412, 0, 3, 17122, 17272, 9012, 18802,
                                                 3807, 3891, 10440, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20692, 0, 3, 17272, 17422, 9138, 19012,
                                                 3891, 3975, 10608, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 20972, 0, 3, 17422, 17572, 9264, 19222,
                                                 3975, 4059, 10776, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21252, 0, 3, 17572, 17722, 9390, 19432,
                                                 4059, 4143, 10944, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 21532, 0, 3, 17722, 17872, 9516, 19642,
                                                 4143, 4227, 11112, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 21812, 0, 3, 18172, 18382, 10104, 20132,
                                                 4395, 4503, 11496, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22172, 0, 3, 18382, 18592, 10272, 20412,
                                                 4503, 4611, 11712, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22532, 0, 3, 18592, 18802, 10440, 20692,
                                                 4611, 4719, 11928, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 22892, 0, 3, 18802, 19012, 10608, 20972,
                                                 4719, 4827, 12144, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 23252, 0, 3, 19012, 19222, 10776, 21252,
                                                 4827, 4935, 12360, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 23612, 0, 3, 19222, 19432, 10944, 21532,
                                                 4935, 5043, 12576, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 23972, 0, 3, 19852, 20132, 11496, 22172,
                                                 5259, 5394, 13332, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 24422, 0, 3, 20132, 20412, 11712, 22532,
                                                 5394, 5529, 13602, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 24872, 0, 3, 20412, 20692, 11928, 22892,
                                                 5529, 5664, 13872, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 25322, 0, 3, 20692, 20972, 12144, 23252,
                                                 5664, 5799, 14142, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 25772, 0, 3, 20972, 21252, 12360, 23612,
                                                 5799, 5934, 14412, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26222, 3, 6204, 6210, 14692, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26237, 3, 6210, 6216, 14702, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26252, 3, 6216, 6222, 14712, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26267, 3, 6222, 6228, 14722, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26282, 3, 6228, 6234, 14732, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26297, 3, 6234, 6240, 14742, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26312, 3, 6240, 6246, 14752, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26327, 3, 6246, 6252, 14762, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26342, 3, 6252, 6258, 14772, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26357, 3, 6258, 6264, 14782, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 26372, 3, 6264, 6270, 14792, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 26387, 0, 3, 14682, 26222, 14832, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26432, 0, 3, 14692, 26237, 14862, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26477, 0, 3, 14702, 26252, 14892, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26522, 0, 3, 14712, 26267, 14922, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26567, 0, 3, 14722, 26282, 14952, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26612, 0, 3, 14732, 26297, 14982, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26657, 0, 3, 14742, 26312, 15012, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26702, 0, 3, 14752, 26327, 15042, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26747, 0, 3, 14762, 26342, 15072, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26792, 0, 3, 14772, 26357, 15102, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26837, 0, 3, 14782, 26372, 15132, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 26882, 0, 3, 14802, 26387, 6516, 6552,
                                                 15162, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26972, 0, 3, 14832, 26432, 6552, 6588,
                                                 15222, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27062, 0, 3, 14862, 26477, 6588, 6624,
                                                 15282, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27152, 0, 3, 14892, 26522, 6624, 6660,
                                                 15342, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27242, 0, 3, 14922, 26567, 6660, 6696,
                                                 15402, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27332, 0, 3, 14952, 26612, 6696, 6732,
                                                 15462, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27422, 0, 3, 14982, 26657, 6732, 6768,
                                                 15522, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27512, 0, 3, 15012, 26702, 6768, 6804,
                                                 15582, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27602, 0, 3, 15042, 26747, 6804, 6840,
                                                 15642, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27692, 0, 3, 15072, 26792, 6840, 6876,
                                                 15702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27782, 0, 3, 15102, 26837, 6876, 6912,
                                                 15762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 27872, 0, 3, 15222, 27062, 6984, 7044,
                                                 15922, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28022, 0, 3, 15282, 27152, 7044, 7104,
                                                 16022, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28172, 0, 3, 15342, 27242, 7104, 7164,
                                                 16122, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28322, 0, 3, 15402, 27332, 7164, 7224,
                                                 16222, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28472, 0, 3, 15462, 27422, 7224, 7284,
                                                 16322, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28622, 0, 3, 15522, 27512, 7284, 7344,
                                                 16422, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28772, 0, 3, 15582, 27602, 7344, 7404,
                                                 16522, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28922, 0, 3, 15642, 27692, 7404, 7464,
                                                 16622, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29072, 0, 3, 15702, 27782, 7464, 7524,
                                                 16722, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29222, 0, 3, 26882, 26972, 15822, 27872,
                                                 7644, 7734, 16822, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29447, 0, 3, 26972, 27062, 15922, 28022,
                                                 7734, 7824, 16972, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29672, 0, 3, 27062, 27152, 16022, 28172,
                                                 7824, 7914, 17122, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29897, 0, 3, 27152, 27242, 16122, 28322,
                                                 7914, 8004, 17272, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30122, 0, 3, 27242, 27332, 16222, 28472,
                                                 8004, 8094, 17422, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30347, 0, 3, 27332, 27422, 16322, 28622,
                                                 8094, 8184, 17572, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30572, 0, 3, 27422, 27512, 16422, 28772,
                                                 8184, 8274, 17722, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30797, 0, 3, 27512, 27602, 16522, 28922,
                                                 8274, 8364, 17872, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31022, 0, 3, 27602, 27692, 16622, 29072,
                                                 8364, 8454, 18022, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31247, 0, 3, 27872, 28022, 16972, 29672,
                                                 8634, 8760, 18382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31562, 0, 3, 28022, 28172, 17122, 29897,
                                                 8760, 8886, 18592, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 31877, 0, 3, 28172, 28322, 17272, 30122,
                                                 8886, 9012, 18802, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32192, 0, 3, 28322, 28472, 17422, 30347,
                                                 9012, 9138, 19012, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32507, 0, 3, 28472, 28622, 17572, 30572,
                                                 9138, 9264, 19222, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32822, 0, 3, 28622, 28772, 17722, 30797,
                                                 9264, 9390, 19432, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33137, 0, 3, 28772, 28922, 17872, 31022,
                                                 9390, 9516, 19642, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 33452, 0, 3, 29222, 29447, 18172, 31247,
                                                 9768, 9936, 19852, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 33872, 0, 3, 29447, 29672, 18382, 31562,
                                                 9936, 10104, 20132, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34292, 0, 3, 29672, 29897, 18592, 31877,
                                                 10104, 10272, 20412, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 34712, 0, 3, 29897, 30122, 18802, 32192,
                                                 10272, 10440, 20692, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35132, 0, 3, 30122, 30347, 19012, 32507,
                                                 10440, 10608, 20972, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35552, 0, 3, 30347, 30572, 19222, 32822,
                                                 10608, 10776, 21252, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 35972, 0, 3, 30572, 30797, 19432, 33137,
                                                 10776, 10944, 21532, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36392, 0, 3, 31247, 31562, 20132, 34292,
                                                 11280, 11496, 22172, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 36932, 0, 3, 31562, 31877, 20412, 34712,
                                                 11496, 11712, 22532, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 37472, 0, 3, 31877, 32192, 20692, 35132,
                                                 11712, 11928, 22892, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 38012, 0, 3, 32192, 32507, 20972, 35552,
                                                 11928, 12144, 23252, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 38552, 0, 3, 32507, 32822, 21252, 35972,
                                                 12144, 12360, 23612, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 39092, 0, 3, 33452, 33872, 21812, 36392,
                                                 12792, 13062, 23972, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 39767, 0, 3, 33872, 34292, 22172, 36932,
                                                 13062, 13332, 24422, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 40442, 0, 3, 34292, 34712, 22532, 37472,
                                                 13332, 13602, 24872, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 41117, 0, 3, 34712, 35132, 22892, 38012,
                                                 13602, 13872, 25322, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 41792, 0, 3, 35132, 35552, 23252, 38552,
                                                 13872, 14142, 25772, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42467, 3, 14682, 14692, 26237, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42488, 3, 14692, 14702, 26252, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42509, 3, 14702, 14712, 26267, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42530, 3, 14712, 14722, 26282, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42551, 3, 14722, 14732, 26297, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42572, 3, 14732, 14742, 26312, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42593, 3, 14742, 14752, 26327, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42614, 3, 14752, 14762, 26342, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42635, 3, 14762, 14772, 26357, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42656, 3, 14772, 14782, 26372, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 42677, 0, 3, 26222, 42467, 26432, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42740, 0, 3, 26237, 42488, 26477, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42803, 0, 3, 26252, 42509, 26522, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42866, 0, 3, 26267, 42530, 26567, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42929, 0, 3, 26282, 42551, 26612, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42992, 0, 3, 26297, 42572, 26657, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43055, 0, 3, 26312, 42593, 26702, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43118, 0, 3, 26327, 42614, 26747, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43181, 0, 3, 26342, 42635, 26792, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43244, 0, 3, 26357, 42656, 26837, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 43307, 0, 3, 26432, 42740, 15162, 15222,
                                                 27062, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43433, 0, 3, 26477, 42803, 15222, 15282,
                                                 27152, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43559, 0, 3, 26522, 42866, 15282, 15342,
                                                 27242, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43685, 0, 3, 26567, 42929, 15342, 15402,
                                                 27332, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43811, 0, 3, 26612, 42992, 15402, 15462,
                                                 27422, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43937, 0, 3, 26657, 43055, 15462, 15522,
                                                 27512, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44063, 0, 3, 26702, 43118, 15522, 15582,
                                                 27602, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44189, 0, 3, 26747, 43181, 15582, 15642,
                                                 27692, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44315, 0, 3, 26792, 43244, 15642, 15702,
                                                 27782, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44441, 0, 3, 27062, 43433, 15822, 15922,
                                                 28022, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44651, 0, 3, 27152, 43559, 15922, 16022,
                                                 28172, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44861, 0, 3, 27242, 43685, 16022, 16122,
                                                 28322, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45071, 0, 3, 27332, 43811, 16122, 16222,
                                                 28472, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45281, 0, 3, 27422, 43937, 16222, 16322,
                                                 28622, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45491, 0, 3, 27512, 44063, 16322, 16422,
                                                 28772, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45701, 0, 3, 27602, 44189, 16422, 16522,
                                                 28922, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45911, 0, 3, 27692, 44315, 16522, 16622,
                                                 29072, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46121, 0, 3, 43307, 43433, 28022, 44651,
                                                 16822, 16972, 29672, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46436, 0, 3, 43433, 43559, 28172, 44861,
                                                 16972, 17122, 29897, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46751, 0, 3, 43559, 43685, 28322, 45071,
                                                 17122, 17272, 30122, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47066, 0, 3, 43685, 43811, 28472, 45281,
                                                 17272, 17422, 30347, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47381, 0, 3, 43811, 43937, 28622, 45491,
                                                 17422, 17572, 30572, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47696, 0, 3, 43937, 44063, 28772, 45701,
                                                 17572, 17722, 30797, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48011, 0, 3, 44063, 44189, 28922, 45911,
                                                 17722, 17872, 31022, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 48326, 0, 3, 44441, 44651, 29672, 46436,
                                                 18172, 18382, 31562, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 48767, 0, 3, 44651, 44861, 29897, 46751,
                                                 18382, 18592, 31877, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49208, 0, 3, 44861, 45071, 30122, 47066,
                                                 18592, 18802, 32192, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49649, 0, 3, 45071, 45281, 30347, 47381,
                                                 18802, 19012, 32507, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50090, 0, 3, 45281, 45491, 30572, 47696,
                                                 19012, 19222, 32822, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50531, 0, 3, 45491, 45701, 30797, 48011,
                                                 19222, 19432, 33137, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 50972, 0, 3, 46121, 46436, 31562, 48767,
                                                 19852, 20132, 34292, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 51560, 0, 3, 46436, 46751, 31877, 49208,
                                                 20132, 20412, 34712, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 52148, 0, 3, 46751, 47066, 32192, 49649,
                                                 20412, 20692, 35132, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 52736, 0, 3, 47066, 47381, 32507, 50090,
                                                 20692, 20972, 35552, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 53324, 0, 3, 47381, 47696, 32822, 50531,
                                                 20972, 21252, 35972, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 53912, 0, 3, 48326, 48767, 34292, 51560,
                                                 21812, 22172, 36932, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 54668, 0, 3, 48767, 49208, 34712, 52148,
                                                 22172, 22532, 37472, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 55424, 0, 3, 49208, 49649, 35132, 52736,
                                                 22532, 22892, 38012, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 56180, 0, 3, 49649, 50090, 35552, 53324,
                                                 22892, 23252, 38552, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 56936, 0, 3, 50972, 51560, 36932, 54668,
                                                 23972, 24422, 40442, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 57881, 0, 3, 51560, 52148, 37472, 55424,
                                                 24422, 24872, 41117, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 58826, 0, 3, 52148, 52736, 38012, 56180,
                                                 24872, 25322, 41792, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59771, 3, 26222, 26237, 42488, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59799, 3, 26237, 26252, 42509, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59827, 3, 26252, 26267, 42530, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59855, 3, 26267, 26282, 42551, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59883, 3, 26282, 26297, 42572, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59911, 3, 26297, 26312, 42593, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59939, 3, 26312, 26327, 42614, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59967, 3, 26327, 26342, 42635, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 59995, 3, 26342, 26357, 42656, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 60023, 0, 3, 42467, 59771, 42740, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60107, 0, 3, 42488, 59799, 42803, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60191, 0, 3, 42509, 59827, 42866, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60275, 0, 3, 42530, 59855, 42929, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60359, 0, 3, 42551, 59883, 42992, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60443, 0, 3, 42572, 59911, 43055, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60527, 0, 3, 42593, 59939, 43118, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60611, 0, 3, 42614, 59967, 43181, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 60695, 0, 3, 42635, 59995, 43244, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 60779, 0, 3, 42677, 60023, 26882, 26972,
                                                 43307, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 60947, 0, 3, 42740, 60107, 26972, 27062,
                                                 43433, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61115, 0, 3, 42803, 60191, 27062, 27152,
                                                 43559, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61283, 0, 3, 42866, 60275, 27152, 27242,
                                                 43685, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61451, 0, 3, 42929, 60359, 27242, 27332,
                                                 43811, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61619, 0, 3, 42992, 60443, 27332, 27422,
                                                 43937, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61787, 0, 3, 43055, 60527, 27422, 27512,
                                                 44063, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 61955, 0, 3, 43118, 60611, 27512, 27602,
                                                 44189, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 62123, 0, 3, 43181, 60695, 27602, 27692,
                                                 44315, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62291, 0, 3, 43433, 61115, 27872, 28022,
                                                 44651, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62571, 0, 3, 43559, 61283, 28022, 28172,
                                                 44861, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62851, 0, 3, 43685, 61451, 28172, 28322,
                                                 45071, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 63131, 0, 3, 43811, 61619, 28322, 28472,
                                                 45281, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 63411, 0, 3, 43937, 61787, 28472, 28622,
                                                 45491, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 63691, 0, 3, 44063, 61955, 28622, 28772,
                                                 45701, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 63971, 0, 3, 44189, 62123, 28772, 28922,
                                                 45911, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 64251, 0, 3, 60779, 60947, 44441, 62291,
                                                 29222, 29447, 46121, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 64671, 0, 3, 60947, 61115, 44651, 62571,
                                                 29447, 29672, 46436, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 65091, 0, 3, 61115, 61283, 44861, 62851,
                                                 29672, 29897, 46751, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 65511, 0, 3, 61283, 61451, 45071, 63131,
                                                 29897, 30122, 47066, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 65931, 0, 3, 61451, 61619, 45281, 63411,
                                                 30122, 30347, 47381, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 66351, 0, 3, 61619, 61787, 45491, 63691,
                                                 30347, 30572, 47696, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 66771, 0, 3, 61787, 61955, 45701, 63971,
                                                 30572, 30797, 48011, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 67191, 0, 3, 62291, 62571, 46436, 65091,
                                                 31247, 31562, 48767, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 67779, 0, 3, 62571, 62851, 46751, 65511,
                                                 31562, 31877, 49208, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 68367, 0, 3, 62851, 63131, 47066, 65931,
                                                 31877, 32192, 49649, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 68955, 0, 3, 63131, 63411, 47381, 66351,
                                                 32192, 32507, 50090, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 69543, 0, 3, 63411, 63691, 47696, 66771,
                                                 32507, 32822, 50531, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 70131, 0, 3, 64251, 64671, 48326, 67191,
                                                 33452, 33872, 50972, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 70915, 0, 3, 64671, 65091, 48767, 67779,
                                                 33872, 34292, 51560, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 71699, 0, 3, 65091, 65511, 49208, 68367,
                                                 34292, 34712, 52148, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 72483, 0, 3, 65511, 65931, 49649, 68955,
                                                 34712, 35132, 52736, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 73267, 0, 3, 65931, 66351, 50090, 69543,
                                                 35132, 35552, 53324, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 74051, 0, 3, 67191, 67779, 51560, 71699,
                                                 36392, 36932, 54668, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 75059, 0, 3, 67779, 68367, 52148, 72483,
                                                 36932, 37472, 55424, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 76067, 0, 3, 68367, 68955, 52736, 73267,
                                                 37472, 38012, 56180, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 77075, 0, 3, 70131, 70915, 53912, 74051,
                                                 39092, 39767, 56936, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 78335, 0, 3, 70915, 71699, 54668, 75059,
                                                 39767, 40442, 57881, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 79595, 0, 3, 71699, 72483, 55424, 76067,
                                                 40442, 41117, 58826, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 80855, 3, 42467, 42488, 59799, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 80891, 3, 42488, 42509, 59827, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 80927, 3, 42509, 42530, 59855, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 80963, 3, 42530, 42551, 59883, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 80999, 3, 42551, 42572, 59911, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 81035, 3, 42572, 42593, 59939, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 81071, 3, 42593, 42614, 59967, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 81107, 3, 42614, 42635, 59995, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 81143, 0, 3, 59771, 80855, 60107, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81251, 0, 3, 59799, 80891, 60191, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81359, 0, 3, 59827, 80927, 60275, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81467, 0, 3, 59855, 80963, 60359, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81575, 0, 3, 59883, 80999, 60443, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81683, 0, 3, 59911, 81035, 60527, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81791, 0, 3, 59939, 81071, 60611, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 81899, 0, 3, 59967, 81107, 60695, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 82007, 0, 3, 60107, 81251, 43307, 43433,
                                                 61115, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 82223, 0, 3, 60191, 81359, 43433, 43559,
                                                 61283, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 82439, 0, 3, 60275, 81467, 43559, 43685,
                                                 61451, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 82655, 0, 3, 60359, 81575, 43685, 43811,
                                                 61619, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 82871, 0, 3, 60443, 81683, 43811, 43937,
                                                 61787, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 83087, 0, 3, 60527, 81791, 43937, 44063,
                                                 61955, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 83303, 0, 3, 60611, 81899, 44063, 44189,
                                                 62123, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 83519, 0, 3, 61115, 82223, 44441, 44651,
                                                 62571, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 83879, 0, 3, 61283, 82439, 44651, 44861,
                                                 62851, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 84239, 0, 3, 61451, 82655, 44861, 45071,
                                                 63131, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 84599, 0, 3, 61619, 82871, 45071, 45281,
                                                 63411, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 84959, 0, 3, 61787, 83087, 45281, 45491,
                                                 63691, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 85319, 0, 3, 61955, 83303, 45491, 45701,
                                                 63971, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 85679, 0, 3, 82007, 82223, 62571, 83879,
                                                 46121, 46436, 65091, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 86219, 0, 3, 82223, 82439, 62851, 84239,
                                                 46436, 46751, 65511, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 86759, 0, 3, 82439, 82655, 63131, 84599,
                                                 46751, 47066, 65931, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 87299, 0, 3, 82655, 82871, 63411, 84959,
                                                 47066, 47381, 66351, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 87839, 0, 3, 82871, 83087, 63691, 85319,
                                                 47381, 47696, 66771, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 88379, 0, 3, 83519, 83879, 65091, 86219,
                                                 48326, 48767, 67779, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 89135, 0, 3, 83879, 84239, 65511, 86759,
                                                 48767, 49208, 68367, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 89891, 0, 3, 84239, 84599, 65931, 87299,
                                                 49208, 49649, 68955, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 90647, 0, 3, 84599, 84959, 66351, 87839,
                                                 49649, 50090, 69543, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 91403, 0, 3, 85679, 86219, 67779, 89135,
                                                 50972, 51560, 71699, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 92411, 0, 3, 86219, 86759, 68367, 89891,
                                                 51560, 52148, 72483, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 93419, 0, 3, 86759, 87299, 68955, 90647,
                                                 52148, 52736, 73267, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 94427, 0, 3, 88379, 89135, 71699, 92411,
                                                 53912, 54668, 75059, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 95723, 0, 3, 89135, 89891, 72483, 93419,
                                                 54668, 55424, 76067, ncols, alpha, beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 97019, 0, 3, 91403, 92411, 75059, 95723,
                                                 56936, 57881, 79595, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98639, 3, 59771, 59799, 80891, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98684, 3, 59799, 59827, 80927, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98729, 3, 59827, 59855, 80963, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98774, 3, 59855, 59883, 80999, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98819, 3, 59883, 59911, 81035, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98864, 3, 59911, 59939, 81071, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 98909, 3, 59939, 59967, 81107, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 98954, 0, 3, 80855, 98639, 81251, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99089, 0, 3, 80891, 98684, 81359, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99224, 0, 3, 80927, 98729, 81467, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99359, 0, 3, 80963, 98774, 81575, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99494, 0, 3, 80999, 98819, 81683, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99629, 0, 3, 81035, 98864, 81791, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 99764, 0, 3, 81071, 98909, 81899, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 99899, 0, 3, 81143, 98954, 60779, 60947,
                                                 82007, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 100169, 0, 3, 81251, 99089, 60947,
                                                 61115, 82223, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 100439, 0, 3, 81359, 99224, 61115,
                                                 61283, 82439, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 100709, 0, 3, 81467, 99359, 61283,
                                                 61451, 82655, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 100979, 0, 3, 81575, 99494, 61451,
                                                 61619, 82871, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 101249, 0, 3, 81683, 99629, 61619,
                                                 61787, 83087, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 101519, 0, 3, 81791, 99764, 61787,
                                                 61955, 83303, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 101789, 0, 3, 82223, 100439, 62291,
                                                 62571, 83879, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 102239, 0, 3, 82439, 100709, 62571,
                                                 62851, 84239, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 102689, 0, 3, 82655, 100979, 62851,
                                                 63131, 84599, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 103139, 0, 3, 82871, 101249, 63131,
                                                 63411, 84959, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 103589, 0, 3, 83087, 101519, 63411,
                                                 63691, 85319, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 104039, 0, 3, 99899, 100169, 83519,
                                                 101789, 64251, 64671, 85679, ncols, alpha, beta,
                                                 p);

            compute_prim_gl_electron_repulsion_0(buffer, 104714, 0, 3, 100169, 100439, 83879,
                                                 102239, 64671, 65091, 86219, ncols, alpha, beta,
                                                 p);

            compute_prim_gl_electron_repulsion_0(buffer, 105389, 0, 3, 100439, 100709, 84239,
                                                 102689, 65091, 65511, 86759, ncols, alpha, beta,
                                                 p);

            compute_prim_gl_electron_repulsion_0(buffer, 106064, 0, 3, 100709, 100979, 84599,
                                                 103139, 65511, 65931, 87299, ncols, alpha, beta,
                                                 p);

            compute_prim_gl_electron_repulsion_0(buffer, 106739, 0, 3, 100979, 101249, 84959,
                                                 103589, 65931, 66351, 87839, ncols, alpha, beta,
                                                 p);

            compute_prim_hl_electron_repulsion_0(buffer, 107414, 0, 3, 101789, 102239, 86219,
                                                 105389, 67191, 67779, 89135, ncols, alpha, beta,
                                                 p);

            compute_prim_hl_electron_repulsion_0(buffer, 108359, 0, 3, 102239, 102689, 86759,
                                                 106064, 67779, 68367, 89891, ncols, alpha, beta,
                                                 p);

            compute_prim_hl_electron_repulsion_0(buffer, 109304, 0, 3, 102689, 103139, 87299,
                                                 106739, 68367, 68955, 90647, ncols, alpha, beta,
                                                 p);

            compute_prim_il_electron_repulsion_0(buffer, 110249, 0, 3, 104039, 104714, 88379,
                                                 107414, 70131, 70915, 91403, ncols, alpha, beta,
                                                 p);

            compute_prim_il_electron_repulsion_0(buffer, 111509, 0, 3, 104714, 105389, 89135,
                                                 108359, 70915, 71699, 92411, ncols, alpha, beta,
                                                 p);

            compute_prim_il_electron_repulsion_0(buffer, 112769, 0, 3, 105389, 106064, 89891,
                                                 109304, 71699, 72483, 93419, ncols, alpha, beta,
                                                 p);

            compute_prim_kl_electron_repulsion_0(buffer, 114029, 0, 3, 107414, 108359, 92411,
                                                 112769, 74051, 75059, 95723, ncols, alpha, beta,
                                                 p);

            compute_prim_ll_electron_repulsion_0(buffer, 115649, 0, 3, 110249, 111509, 94427,
                                                 114029, 77075, 78335, 97019, ncols, alpha, beta,
                                                 p);

            simdfunc::contract_primitives(buffer, 117674, 115649, 2025, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 119699, 117674, 45, nmax);

    simdtrf::transform_l_outer_tri(values, nvalues, buffer, 119699, nmax);
}

}  // namespace simdt2ceri

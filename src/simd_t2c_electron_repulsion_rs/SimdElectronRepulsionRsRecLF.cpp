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


#include "SimdElectronRepulsionRsRecLF.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_lf_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_lf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 21779, 20564, 900, nvalues);

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
                                                9, 10, 11}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 18, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 7, 8, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 8, 9, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 9, 10, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 10, 11, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 11, 12, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 12, 13, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 13, 14, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 14, 15, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 15, 16, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 19, 20, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 20, 21, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 21, 22, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 22, 23, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 23, 24, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 24, 25, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 25, 26, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 26, 27, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 27, 28, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 51, 54, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 60, 63, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 63, 66, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 66, 69, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 69, 72, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 72, 75, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 75, 78, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 78, 81, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 81, 84, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 358, 0, 90, 96, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 373, 0, 96, 102, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 388, 0, 102, 108, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 403, 0, 108, 114, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 418, 0, 114, 120, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 433, 0, 120, 126, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 448, 0, 126, 132, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 463, 0, 144, 150, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 478, 0, 150, 156, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 493, 0, 156, 162, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 508, 0, 162, 168, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 523, 0, 168, 174, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 538, 0, 174, 180, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 553, 0, 180, 186, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 568, 0, 198, 208, 373, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 589, 0, 208, 218, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 610, 0, 218, 228, 403, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 631, 0, 228, 238, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 652, 0, 238, 248, 433, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 673, 0, 248, 258, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 694, 0, 278, 288, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 715, 0, 288, 298, 493, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 736, 0, 298, 308, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 757, 0, 308, 318, 523, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 778, 0, 318, 328, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 799, 0, 328, 338, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 820, 0, 358, 373, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 848, 0, 373, 388, 610, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 876, 0, 388, 403, 631, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 904, 0, 403, 418, 652, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 932, 0, 418, 433, 673, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 960, 0, 463, 478, 715, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 988, 0, 478, 493, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1016, 0, 493, 508, 757, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1044, 0, 508, 523, 778, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1072, 0, 523, 538, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1100, 0, 568, 589, 848, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1136, 0, 589, 610, 876, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1172, 0, 610, 631, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1208, 0, 631, 652, 932, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1244, 0, 694, 715, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1280, 0, 715, 736, 1016, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1316, 0, 736, 757, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1352, 0, 757, 778, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1388, 0, 820, 848, 1136, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1433, 0, 848, 876, 1172, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1478, 0, 876, 904, 1208, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1523, 0, 960, 988, 1280, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1568, 0, 988, 1016, 1316, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1613, 0, 1016, 1044, 1352, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1658, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1661, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1664, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1667, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1670, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1673, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1676, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1679, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1682, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1685, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1688, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1691, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1694, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1697, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1700, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1703, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1706, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1709, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1712, 3, 8, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1721, 3, 9, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1730, 3, 10, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1739, 3, 11, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1748, 3, 12, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1757, 3, 13, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1766, 3, 14, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1775, 3, 15, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1784, 3, 16, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1793, 3, 20, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1802, 3, 21, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1811, 3, 22, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1820, 3, 23, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1829, 3, 24, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1838, 3, 25, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1847, 3, 26, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1856, 3, 27, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1865, 3, 28, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1874, 0, 3, 30, 1712, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1892, 0, 3, 33, 1721, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1910, 0, 3, 36, 1730, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1928, 0, 3, 39, 1739, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1946, 0, 3, 42, 1748, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1964, 0, 3, 45, 1757, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1982, 0, 3, 48, 1766, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2000, 0, 3, 51, 1775, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2018, 0, 3, 54, 1784, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2036, 0, 3, 60, 1793, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2054, 0, 3, 63, 1802, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2072, 0, 3, 66, 1811, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2090, 0, 3, 69, 1820, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2108, 0, 3, 72, 1829, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2126, 0, 3, 75, 1838, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2144, 0, 3, 78, 1847, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2162, 0, 3, 81, 1856, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2180, 0, 3, 84, 1865, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2198, 0, 3, 96, 1910, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2228, 0, 3, 102, 1928, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2258, 0, 3, 108, 1946, 228, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2288, 0, 3, 114, 1964, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2318, 0, 3, 120, 1982, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2348, 0, 3, 126, 2000, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2378, 0, 3, 132, 2018, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2408, 0, 3, 150, 2072, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2438, 0, 3, 156, 2090, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2468, 0, 3, 162, 2108, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2498, 0, 3, 168, 2126, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2528, 0, 3, 174, 2144, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2558, 0, 3, 180, 2162, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2588, 0, 3, 186, 2180, 348, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2618, 0, 3, 198, 2198, 358, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2663, 0, 3, 208, 2228, 373, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2708, 0, 3, 218, 2258, 388, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2753, 0, 3, 228, 2288, 403, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2798, 0, 3, 238, 2318, 418, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2843, 0, 3, 248, 2348, 433, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2888, 0, 3, 258, 2378, 448, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2933, 0, 3, 278, 2408, 463, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2978, 0, 3, 288, 2438, 478, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3023, 0, 3, 298, 2468, 493, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3068, 0, 3, 308, 2498, 508, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3113, 0, 3, 318, 2528, 523, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3158, 0, 3, 328, 2558, 538, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3203, 0, 3, 338, 2588, 553, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3248, 0, 3, 373, 2708, 589, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3311, 0, 3, 388, 2753, 610, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3374, 0, 3, 403, 2798, 631, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3437, 0, 3, 418, 2843, 652, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3500, 0, 3, 433, 2888, 673, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3563, 0, 3, 478, 3023, 715, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3626, 0, 3, 493, 3068, 736, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3689, 0, 3, 508, 3113, 757, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3752, 0, 3, 523, 3158, 778, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3815, 0, 3, 538, 3203, 799, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3878, 0, 3, 568, 3248, 820, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3962, 0, 3, 589, 3311, 848, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4046, 0, 3, 610, 3374, 876, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4130, 0, 3, 631, 3437, 904, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4214, 0, 3, 652, 3500, 932, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4298, 0, 3, 694, 3563, 960, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4382, 0, 3, 715, 3626, 988, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4466, 0, 3, 736, 3689, 1016, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4550, 0, 3, 757, 3752, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4634, 0, 3, 778, 3815, 1072, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4718, 0, 3, 848, 4046, 1136, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4826, 0, 3, 876, 4130, 1172, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4934, 0, 3, 904, 4214, 1208, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5042, 0, 3, 988, 4466, 1280, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5150, 0, 3, 1016, 4550, 1316, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5258, 0, 3, 1044, 4634, 1352, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5366, 0, 3, 1100, 4718, 1388, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5501, 0, 3, 1136, 4826, 1433, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5636, 0, 3, 1172, 4934, 1478, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5771, 0, 3, 1244, 5042, 1523, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5906, 0, 3, 1280, 5150, 1568, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6041, 0, 3, 1316, 5258, 1613, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6176, 3, 8, 9, 1661, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6182, 3, 9, 10, 1664, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6188, 3, 10, 11, 1667, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6194, 3, 11, 12, 1670, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6200, 3, 12, 13, 1673, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6206, 3, 13, 14, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6212, 3, 14, 15, 1679, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6218, 3, 15, 16, 1682, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6224, 3, 20, 21, 1688, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6230, 3, 21, 22, 1691, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6236, 3, 22, 23, 1694, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6242, 3, 23, 24, 1697, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6248, 3, 24, 25, 1700, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6254, 3, 25, 26, 1703, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6260, 3, 26, 27, 1706, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6266, 3, 27, 28, 1709, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6272, 0, 3, 1658, 6176, 1721, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6290, 0, 3, 1661, 6182, 1730, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6308, 0, 3, 1664, 6188, 1739, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6326, 0, 3, 1667, 6194, 1748, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6344, 0, 3, 1670, 6200, 1757, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6362, 0, 3, 1673, 6206, 1766, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6380, 0, 3, 1676, 6212, 1775, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6398, 0, 3, 1679, 6218, 1784, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6416, 0, 3, 1685, 6224, 1802, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6434, 0, 3, 1688, 6230, 1811, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6452, 0, 3, 1691, 6236, 1820, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6470, 0, 3, 1694, 6242, 1829, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6488, 0, 3, 1697, 6248, 1838, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6506, 0, 3, 1700, 6254, 1847, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6524, 0, 3, 1703, 6260, 1856, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6542, 0, 3, 1706, 6266, 1865, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6560, 0, 3, 1721, 6290, 90, 96, 1910,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6596, 0, 3, 1730, 6308, 96, 102, 1928,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6632, 0, 3, 1739, 6326, 102, 108, 1946,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6668, 0, 3, 1748, 6344, 108, 114, 1964,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6704, 0, 3, 1757, 6362, 114, 120, 1982,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6740, 0, 3, 1766, 6380, 120, 126, 2000,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6776, 0, 3, 1775, 6398, 126, 132, 2018,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6812, 0, 3, 1802, 6434, 144, 150, 2072,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6848, 0, 3, 1811, 6452, 150, 156, 2090,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6884, 0, 3, 1820, 6470, 156, 162, 2108,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6920, 0, 3, 1829, 6488, 162, 168, 2126,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6956, 0, 3, 1838, 6506, 168, 174, 2144,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6992, 0, 3, 1847, 6524, 174, 180, 2162,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7028, 0, 3, 1856, 6542, 180, 186, 2180,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7064, 0, 3, 1910, 6596, 198, 208, 2228,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7124, 0, 3, 1928, 6632, 208, 218, 2258,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7184, 0, 3, 1946, 6668, 218, 228, 2288,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7244, 0, 3, 1964, 6704, 228, 238, 2318,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7304, 0, 3, 1982, 6740, 238, 248, 2348,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7364, 0, 3, 2000, 6776, 248, 258, 2378,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7424, 0, 3, 2072, 6848, 278, 288, 2438,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7484, 0, 3, 2090, 6884, 288, 298, 2468,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7544, 0, 3, 2108, 6920, 298, 308, 2498,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7604, 0, 3, 2126, 6956, 308, 318, 2528,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7664, 0, 3, 2144, 6992, 318, 328, 2558,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7724, 0, 3, 2162, 7028, 328, 338, 2588,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7784, 0, 3, 6560, 6596, 2228, 7124, 358,
                                                 373, 2708, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7874, 0, 3, 6596, 6632, 2258, 7184, 373,
                                                 388, 2753, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7964, 0, 3, 6632, 6668, 2288, 7244, 388,
                                                 403, 2798, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8054, 0, 3, 6668, 6704, 2318, 7304, 403,
                                                 418, 2843, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8144, 0, 3, 6704, 6740, 2348, 7364, 418,
                                                 433, 2888, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8234, 0, 3, 6812, 6848, 2438, 7484, 463,
                                                 478, 3023, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8324, 0, 3, 6848, 6884, 2468, 7544, 478,
                                                 493, 3068, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8414, 0, 3, 6884, 6920, 2498, 7604, 493,
                                                 508, 3113, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8504, 0, 3, 6920, 6956, 2528, 7664, 508,
                                                 523, 3158, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8594, 0, 3, 6956, 6992, 2558, 7724, 523,
                                                 538, 3203, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8684, 0, 3, 7064, 7124, 2708, 7874, 568,
                                                 589, 3311, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8810, 0, 3, 7124, 7184, 2753, 7964, 589,
                                                 610, 3374, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8936, 0, 3, 7184, 7244, 2798, 8054, 610,
                                                 631, 3437, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9062, 0, 3, 7244, 7304, 2843, 8144, 631,
                                                 652, 3500, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9188, 0, 3, 7424, 7484, 3023, 8324, 694,
                                                 715, 3626, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9314, 0, 3, 7484, 7544, 3068, 8414, 715,
                                                 736, 3689, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9440, 0, 3, 7544, 7604, 3113, 8504, 736,
                                                 757, 3752, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9566, 0, 3, 7604, 7664, 3158, 8594, 757,
                                                 778, 3815, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9692, 0, 3, 7784, 7874, 3311, 8810, 820,
                                                 848, 4046, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9860, 0, 3, 7874, 7964, 3374, 8936, 848,
                                                 876, 4130, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10028, 0, 3, 7964, 8054, 3437, 9062,
                                                 876, 904, 4214, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10196, 0, 3, 8234, 8324, 3626, 9314,
                                                 960, 988, 4466, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10364, 0, 3, 8324, 8414, 3689, 9440,
                                                 988, 1016, 4550, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10532, 0, 3, 8414, 8504, 3752, 9566,
                                                 1016, 1044, 4634, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10700, 0, 3, 8684, 8810, 4046, 9860,
                                                 1100, 1136, 4826, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10916, 0, 3, 8810, 8936, 4130, 10028,
                                                 1136, 1172, 4934, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11132, 0, 3, 9188, 9314, 4466, 10364,
                                                 1244, 1280, 5150, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11348, 0, 3, 9314, 9440, 4550, 10532,
                                                 1280, 1316, 5258, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11564, 0, 3, 9692, 9860, 4826, 10916,
                                                 1388, 1433, 5636, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11834, 0, 3, 10196, 10364, 5150, 11348,
                                                 1523, 1568, 6041, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12104, 3, 1658, 1661, 6182, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12114, 3, 1661, 1664, 6188, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12124, 3, 1664, 1667, 6194, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12134, 3, 1667, 1670, 6200, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12144, 3, 1670, 1673, 6206, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12154, 3, 1673, 1676, 6212, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12164, 3, 1676, 1679, 6218, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12174, 3, 1685, 1688, 6230, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12184, 3, 1688, 1691, 6236, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12194, 3, 1691, 1694, 6242, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12204, 3, 1694, 1697, 6248, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12214, 3, 1697, 1700, 6254, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12224, 3, 1700, 1703, 6260, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12234, 3, 1703, 1706, 6266, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 12244, 0, 3, 6176, 12104, 6290, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12274, 0, 3, 6182, 12114, 6308, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12304, 0, 3, 6188, 12124, 6326, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12334, 0, 3, 6194, 12134, 6344, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12364, 0, 3, 6200, 12144, 6362, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12394, 0, 3, 6206, 12154, 6380, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12424, 0, 3, 6212, 12164, 6398, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12454, 0, 3, 6224, 12174, 6434, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12484, 0, 3, 6230, 12184, 6452, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12514, 0, 3, 6236, 12194, 6470, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12544, 0, 3, 6242, 12204, 6488, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12574, 0, 3, 6248, 12214, 6506, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12604, 0, 3, 6254, 12224, 6524, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12634, 0, 3, 6260, 12234, 6542, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 12664, 0, 3, 6272, 12244, 1874, 1892,
                                                 6560, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12724, 0, 3, 6290, 12274, 1892, 1910,
                                                 6596, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12784, 0, 3, 6308, 12304, 1910, 1928,
                                                 6632, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12844, 0, 3, 6326, 12334, 1928, 1946,
                                                 6668, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12904, 0, 3, 6344, 12364, 1946, 1964,
                                                 6704, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12964, 0, 3, 6362, 12394, 1964, 1982,
                                                 6740, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13024, 0, 3, 6380, 12424, 1982, 2000,
                                                 6776, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13084, 0, 3, 6416, 12454, 2036, 2054,
                                                 6812, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13144, 0, 3, 6434, 12484, 2054, 2072,
                                                 6848, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13204, 0, 3, 6452, 12514, 2072, 2090,
                                                 6884, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13264, 0, 3, 6470, 12544, 2090, 2108,
                                                 6920, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13324, 0, 3, 6488, 12574, 2108, 2126,
                                                 6956, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13384, 0, 3, 6506, 12604, 2126, 2144,
                                                 6992, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13444, 0, 3, 6524, 12634, 2144, 2162,
                                                 7028, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13504, 0, 3, 6596, 12784, 2198, 2228,
                                                 7124, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13604, 0, 3, 6632, 12844, 2228, 2258,
                                                 7184, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13704, 0, 3, 6668, 12904, 2258, 2288,
                                                 7244, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13804, 0, 3, 6704, 12964, 2288, 2318,
                                                 7304, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13904, 0, 3, 6740, 13024, 2318, 2348,
                                                 7364, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14004, 0, 3, 6848, 13204, 2408, 2438,
                                                 7484, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14104, 0, 3, 6884, 13264, 2438, 2468,
                                                 7544, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14204, 0, 3, 6920, 13324, 2468, 2498,
                                                 7604, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14304, 0, 3, 6956, 13384, 2498, 2528,
                                                 7664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14404, 0, 3, 6992, 13444, 2528, 2558,
                                                 7724, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14504, 0, 3, 12664, 12724, 7064, 13504,
                                                 2618, 2663, 7784, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14654, 0, 3, 12724, 12784, 7124, 13604,
                                                 2663, 2708, 7874, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14804, 0, 3, 12784, 12844, 7184, 13704,
                                                 2708, 2753, 7964, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14954, 0, 3, 12844, 12904, 7244, 13804,
                                                 2753, 2798, 8054, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15104, 0, 3, 12904, 12964, 7304, 13904,
                                                 2798, 2843, 8144, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15254, 0, 3, 13084, 13144, 7424, 14004,
                                                 2933, 2978, 8234, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15404, 0, 3, 13144, 13204, 7484, 14104,
                                                 2978, 3023, 8324, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15554, 0, 3, 13204, 13264, 7544, 14204,
                                                 3023, 3068, 8414, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15704, 0, 3, 13264, 13324, 7604, 14304,
                                                 3068, 3113, 8504, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15854, 0, 3, 13324, 13384, 7664, 14404,
                                                 3113, 3158, 8594, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16004, 0, 3, 13504, 13604, 7874, 14804,
                                                 3248, 3311, 8810, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16214, 0, 3, 13604, 13704, 7964, 14954,
                                                 3311, 3374, 8936, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16424, 0, 3, 13704, 13804, 8054, 15104,
                                                 3374, 3437, 9062, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16634, 0, 3, 14004, 14104, 8324, 15554,
                                                 3563, 3626, 9314, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16844, 0, 3, 14104, 14204, 8414, 15704,
                                                 3626, 3689, 9440, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17054, 0, 3, 14204, 14304, 8504, 15854,
                                                 3689, 3752, 9566, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17264, 0, 3, 14504, 14654, 8684, 16004,
                                                 3878, 3962, 9692, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17544, 0, 3, 14654, 14804, 8810, 16214,
                                                 3962, 4046, 9860, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17824, 0, 3, 14804, 14954, 8936, 16424,
                                                 4046, 4130, 10028, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18104, 0, 3, 15254, 15404, 9188, 16634,
                                                 4298, 4382, 10196, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18384, 0, 3, 15404, 15554, 9314, 16844,
                                                 4382, 4466, 10364, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18664, 0, 3, 15554, 15704, 9440, 17054,
                                                 4466, 4550, 10532, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 18944, 0, 3, 16004, 16214, 9860, 17824,
                                                 4718, 4826, 10916, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19304, 0, 3, 16634, 16844, 10364, 18664,
                                                 5042, 5150, 11348, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 19664, 0, 3, 17264, 17544, 10700, 18944,
                                                 5366, 5501, 11564, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 20114, 0, 3, 18104, 18384, 11132, 19304,
                                                 5771, 5906, 11834, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 20564, 19664, 900, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 21464, 21014, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 21464, 7, nmax);

    simdtrf::transform_f_inner(buffer, 21464, 20564, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 119 * nvalues, nvalues, buffer, 21464, 7, nmax);
}

}  // namespace simdt2ceri

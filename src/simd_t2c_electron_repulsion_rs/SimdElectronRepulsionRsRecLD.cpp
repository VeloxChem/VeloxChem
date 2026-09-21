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


#include "SimdElectronRepulsionRsRecLD.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_ld_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ld_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 10323, 9558, 540, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 18, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 7, 8, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 8, 9, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 9, 10, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 10, 11, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 11, 12, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 12, 13, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 13, 14, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 14, 15, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 15, 16, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 22, 23, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 23, 24, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 24, 25, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 25, 26, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 26, 27, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 27, 28, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 262, 0, 57, 60, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 60, 63, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 282, 0, 63, 66, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 292, 0, 66, 69, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 302, 0, 69, 72, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 72, 75, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 75, 78, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 332, 0, 84, 90, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 347, 0, 90, 96, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 362, 0, 96, 102, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 377, 0, 102, 108, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 392, 0, 108, 114, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 407, 0, 114, 120, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 422, 0, 120, 126, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 437, 0, 138, 144, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 452, 0, 144, 150, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 467, 0, 150, 156, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 482, 0, 156, 162, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 497, 0, 162, 168, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 512, 0, 168, 174, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 527, 0, 174, 180, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 542, 0, 192, 202, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 563, 0, 202, 212, 377, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 584, 0, 212, 222, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 605, 0, 222, 232, 407, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 626, 0, 232, 242, 422, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 647, 0, 262, 272, 467, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 668, 0, 272, 282, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 689, 0, 282, 292, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 710, 0, 292, 302, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 731, 0, 302, 312, 527, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 752, 0, 332, 347, 542, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 780, 0, 347, 362, 563, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 808, 0, 362, 377, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 836, 0, 377, 392, 605, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 864, 0, 392, 407, 626, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 892, 0, 437, 452, 647, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 920, 0, 452, 467, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 948, 0, 467, 482, 689, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 976, 0, 482, 497, 710, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1004, 0, 497, 512, 731, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1032, 0, 542, 563, 808, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1068, 0, 563, 584, 836, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1104, 0, 584, 605, 864, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1140, 0, 647, 668, 948, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1176, 0, 668, 689, 976, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1212, 0, 689, 710, 1004, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1248, 0, 752, 780, 1032, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1293, 0, 780, 808, 1068, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1338, 0, 808, 836, 1104, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1383, 0, 892, 920, 1140, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1428, 0, 920, 948, 1176, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1473, 0, 948, 976, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1518, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1521, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1524, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1527, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1530, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1533, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1536, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1539, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1542, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1545, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1548, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1551, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1554, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1557, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1560, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1563, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1566, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1575, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1584, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1593, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1602, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1611, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1620, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1629, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1638, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1647, 3, 22, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1656, 3, 23, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1665, 3, 24, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1674, 3, 25, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1683, 3, 26, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1692, 3, 27, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1701, 3, 28, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1710, 0, 3, 33, 1575, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1728, 0, 3, 36, 1584, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1746, 0, 3, 39, 1593, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1764, 0, 3, 42, 1602, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1782, 0, 3, 45, 1611, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1800, 0, 3, 48, 1620, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1818, 0, 3, 51, 1629, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1836, 0, 3, 60, 1647, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1854, 0, 3, 63, 1656, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1872, 0, 3, 66, 1665, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1890, 0, 3, 69, 1674, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1908, 0, 3, 72, 1683, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1926, 0, 3, 75, 1692, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1944, 0, 3, 78, 1701, 186, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1962, 0, 3, 96, 1728, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1992, 0, 3, 102, 1746, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2022, 0, 3, 108, 1764, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2052, 0, 3, 114, 1782, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2082, 0, 3, 120, 1800, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2112, 0, 3, 126, 1818, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2142, 0, 3, 150, 1854, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2172, 0, 3, 156, 1872, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2202, 0, 3, 162, 1890, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2232, 0, 3, 168, 1908, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2262, 0, 3, 174, 1926, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2292, 0, 3, 180, 1944, 322, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2322, 0, 3, 202, 1992, 362, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2367, 0, 3, 212, 2022, 377, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2412, 0, 3, 222, 2052, 392, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2457, 0, 3, 232, 2082, 407, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2502, 0, 3, 242, 2112, 422, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2547, 0, 3, 272, 2172, 467, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2592, 0, 3, 282, 2202, 482, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2637, 0, 3, 292, 2232, 497, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2682, 0, 3, 302, 2262, 512, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2727, 0, 3, 312, 2292, 527, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2772, 0, 3, 362, 2367, 563, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2835, 0, 3, 377, 2412, 584, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2898, 0, 3, 392, 2457, 605, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2961, 0, 3, 407, 2502, 626, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3024, 0, 3, 467, 2592, 668, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3087, 0, 3, 482, 2637, 689, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3150, 0, 3, 497, 2682, 710, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3213, 0, 3, 512, 2727, 731, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3276, 0, 3, 563, 2835, 808, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3360, 0, 3, 584, 2898, 836, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3444, 0, 3, 605, 2961, 864, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3528, 0, 3, 668, 3087, 948, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3612, 0, 3, 689, 3150, 976, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3696, 0, 3, 710, 3213, 1004, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3780, 0, 3, 808, 3360, 1068, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3888, 0, 3, 836, 3444, 1104, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3996, 0, 3, 948, 3612, 1176, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4104, 0, 3, 976, 3696, 1212, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4212, 0, 3, 1068, 3888, 1338, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 4347, 0, 3, 1176, 4104, 1473, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4482, 3, 9, 10, 1521, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4488, 3, 10, 11, 1524, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4494, 3, 11, 12, 1527, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4500, 3, 12, 13, 1530, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4506, 3, 13, 14, 1533, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4512, 3, 14, 15, 1536, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4518, 3, 15, 16, 1539, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4524, 3, 21, 22, 1545, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4530, 3, 22, 23, 1548, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4536, 3, 23, 24, 1551, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4542, 3, 24, 25, 1554, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4548, 3, 25, 26, 1557, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4554, 3, 26, 27, 1560, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4560, 3, 27, 28, 1563, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4566, 0, 3, 1518, 4482, 1575, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4584, 0, 3, 1521, 4488, 1584, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4602, 0, 3, 1524, 4494, 1593, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4620, 0, 3, 1527, 4500, 1602, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4638, 0, 3, 1530, 4506, 1611, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4656, 0, 3, 1533, 4512, 1620, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4674, 0, 3, 1536, 4518, 1629, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4692, 0, 3, 1542, 4524, 1647, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4710, 0, 3, 1545, 4530, 1656, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4728, 0, 3, 1548, 4536, 1665, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4746, 0, 3, 1551, 4542, 1674, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4764, 0, 3, 1554, 4548, 1683, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4782, 0, 3, 1557, 4554, 1692, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4800, 0, 3, 1560, 4560, 1701, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4818, 0, 3, 1566, 4566, 84, 90, 1710,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4854, 0, 3, 1575, 4584, 90, 96, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4890, 0, 3, 1584, 4602, 96, 102, 1746,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4926, 0, 3, 1593, 4620, 102, 108, 1764,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4962, 0, 3, 1602, 4638, 108, 114, 1782,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4998, 0, 3, 1611, 4656, 114, 120, 1800,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5034, 0, 3, 1620, 4674, 120, 126, 1818,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5070, 0, 3, 1638, 4692, 138, 144, 1836,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5106, 0, 3, 1647, 4710, 144, 150, 1854,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5142, 0, 3, 1656, 4728, 150, 156, 1872,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5178, 0, 3, 1665, 4746, 156, 162, 1890,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5214, 0, 3, 1674, 4764, 162, 168, 1908,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5250, 0, 3, 1683, 4782, 168, 174, 1926,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5286, 0, 3, 1692, 4800, 174, 180, 1944,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5322, 0, 3, 1728, 4890, 192, 202, 1992,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5382, 0, 3, 1746, 4926, 202, 212, 2022,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5442, 0, 3, 1764, 4962, 212, 222, 2052,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5502, 0, 3, 1782, 4998, 222, 232, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5562, 0, 3, 1800, 5034, 232, 242, 2112,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5622, 0, 3, 1854, 5142, 262, 272, 2172,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5682, 0, 3, 1872, 5178, 272, 282, 2202,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5742, 0, 3, 1890, 5214, 282, 292, 2232,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5802, 0, 3, 1908, 5250, 292, 302, 2262,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5862, 0, 3, 1926, 5286, 302, 312, 2292,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5922, 0, 3, 4818, 4854, 1962, 5322, 332,
                                                 347, 2322, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6012, 0, 3, 4854, 4890, 1992, 5382, 347,
                                                 362, 2367, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6102, 0, 3, 4890, 4926, 2022, 5442, 362,
                                                 377, 2412, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6192, 0, 3, 4926, 4962, 2052, 5502, 377,
                                                 392, 2457, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6282, 0, 3, 4962, 4998, 2082, 5562, 392,
                                                 407, 2502, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6372, 0, 3, 5070, 5106, 2142, 5622, 437,
                                                 452, 2547, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6462, 0, 3, 5106, 5142, 2172, 5682, 452,
                                                 467, 2592, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6552, 0, 3, 5142, 5178, 2202, 5742, 467,
                                                 482, 2637, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6642, 0, 3, 5178, 5214, 2232, 5802, 482,
                                                 497, 2682, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6732, 0, 3, 5214, 5250, 2262, 5862, 497,
                                                 512, 2727, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6822, 0, 3, 5322, 5382, 2367, 6102, 542,
                                                 563, 2835, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6948, 0, 3, 5382, 5442, 2412, 6192, 563,
                                                 584, 2898, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7074, 0, 3, 5442, 5502, 2457, 6282, 584,
                                                 605, 2961, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7200, 0, 3, 5622, 5682, 2592, 6552, 647,
                                                 668, 3087, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7326, 0, 3, 5682, 5742, 2637, 6642, 668,
                                                 689, 3150, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7452, 0, 3, 5742, 5802, 2682, 6732, 689,
                                                 710, 3213, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7578, 0, 3, 5922, 6012, 2772, 6822, 752,
                                                 780, 3276, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7746, 0, 3, 6012, 6102, 2835, 6948, 780,
                                                 808, 3360, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7914, 0, 3, 6102, 6192, 2898, 7074, 808,
                                                 836, 3444, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8082, 0, 3, 6372, 6462, 3024, 7200, 892,
                                                 920, 3528, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8250, 0, 3, 6462, 6552, 3087, 7326, 920,
                                                 948, 3612, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8418, 0, 3, 6552, 6642, 3150, 7452, 948,
                                                 976, 3696, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8586, 0, 3, 6822, 6948, 3360, 7914,
                                                 1032, 1068, 3888, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8802, 0, 3, 7200, 7326, 3612, 8418,
                                                 1140, 1176, 4104, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 9018, 0, 3, 7578, 7746, 3780, 8586,
                                                 1248, 1293, 4212, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 9288, 0, 3, 8082, 8250, 3996, 8802,
                                                 1383, 1428, 4347, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9558, 9018, 540, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 10098, 9828, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 10098, 5, nmax);

    simdtrf::transform_d_inner(buffer, 10098, 9558, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 85 * nvalues, nvalues, buffer, 10098, 5, nmax);
}

}  // namespace simdt2ceri

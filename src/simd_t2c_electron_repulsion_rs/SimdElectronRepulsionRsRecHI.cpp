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


#include "SimdElectronRepulsionRsRecHI.hpp"

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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_hi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 29605, 28156, 1176, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 7, 8, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 8, 9, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 9, 10, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 10, 11, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 11, 12, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 12, 13, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 13, 14, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 14, 15, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 15, 16, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 19, 20, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 20, 21, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 21, 22, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 22, 23, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 23, 24, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 24, 25, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 25, 26, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 26, 27, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 27, 28, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 51, 54, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 54, 57, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 63, 66, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 304, 0, 66, 69, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 314, 0, 69, 72, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 72, 75, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 75, 78, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 78, 81, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 81, 84, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 84, 87, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 87, 90, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 384, 0, 96, 102, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 399, 0, 102, 108, 234, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 414, 0, 108, 114, 244, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 429, 0, 114, 120, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 444, 0, 120, 126, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 459, 0, 126, 132, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 474, 0, 132, 138, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 489, 0, 150, 156, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 504, 0, 156, 162, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 519, 0, 162, 168, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 534, 0, 168, 174, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 549, 0, 174, 180, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 564, 0, 180, 186, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 579, 0, 186, 192, 374, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 594, 0, 204, 214, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 615, 0, 214, 224, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 636, 0, 224, 234, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 657, 0, 234, 244, 429, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 678, 0, 244, 254, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 699, 0, 254, 264, 459, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 720, 0, 264, 274, 474, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 741, 0, 294, 304, 489, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 762, 0, 304, 314, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 783, 0, 314, 324, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 804, 0, 324, 334, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 825, 0, 334, 344, 549, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 846, 0, 344, 354, 564, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 867, 0, 354, 364, 579, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 888, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 891, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 894, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 897, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 900, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 903, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 906, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 909, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 912, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 915, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 918, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 921, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 924, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 927, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 930, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 933, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 936, 3, 9, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 945, 3, 10, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 954, 3, 11, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 963, 3, 12, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 972, 3, 13, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 981, 3, 14, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 990, 3, 15, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 999, 3, 16, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1008, 3, 21, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1017, 3, 22, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1026, 3, 23, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1035, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1044, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1053, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1062, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1071, 3, 28, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1080, 0, 3, 36, 936, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1098, 0, 3, 39, 945, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1116, 0, 3, 42, 954, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1134, 0, 3, 45, 963, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1152, 0, 3, 48, 972, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1170, 0, 3, 51, 981, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1188, 0, 3, 54, 990, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1206, 0, 3, 57, 999, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1224, 0, 3, 69, 1008, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1242, 0, 3, 72, 1017, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1260, 0, 3, 75, 1026, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1278, 0, 3, 78, 1035, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1296, 0, 3, 81, 1044, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1314, 0, 3, 84, 1053, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1332, 0, 3, 87, 1062, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1350, 0, 3, 90, 1071, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1368, 0, 3, 102, 1098, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1398, 0, 3, 108, 1116, 234, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1428, 0, 3, 114, 1134, 244, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1458, 0, 3, 120, 1152, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1488, 0, 3, 126, 1170, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1518, 0, 3, 132, 1188, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1548, 0, 3, 138, 1206, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1578, 0, 3, 156, 1242, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1608, 0, 3, 162, 1260, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1638, 0, 3, 168, 1278, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1668, 0, 3, 174, 1296, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1698, 0, 3, 180, 1314, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1728, 0, 3, 186, 1332, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1758, 0, 3, 192, 1350, 374, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1788, 0, 3, 224, 1398, 399, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1833, 0, 3, 234, 1428, 414, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1878, 0, 3, 244, 1458, 429, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1923, 0, 3, 254, 1488, 444, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1968, 0, 3, 264, 1518, 459, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2013, 0, 3, 274, 1548, 474, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2058, 0, 3, 314, 1608, 504, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2103, 0, 3, 324, 1638, 519, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2148, 0, 3, 334, 1668, 534, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2193, 0, 3, 344, 1698, 549, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2238, 0, 3, 354, 1728, 564, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2283, 0, 3, 364, 1758, 579, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2328, 0, 3, 399, 1833, 636, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2391, 0, 3, 414, 1878, 657, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2454, 0, 3, 429, 1923, 678, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2517, 0, 3, 444, 1968, 699, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2580, 0, 3, 459, 2013, 720, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2643, 0, 3, 504, 2103, 783, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2706, 0, 3, 519, 2148, 804, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2769, 0, 3, 534, 2193, 825, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2832, 0, 3, 549, 2238, 846, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2895, 0, 3, 564, 2283, 867, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2958, 3, 9, 10, 891, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2964, 3, 10, 11, 894, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2970, 3, 11, 12, 897, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2976, 3, 12, 13, 900, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2982, 3, 13, 14, 903, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2988, 3, 14, 15, 906, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2994, 3, 15, 16, 909, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3000, 3, 21, 22, 915, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3006, 3, 22, 23, 918, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3012, 3, 23, 24, 921, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3018, 3, 24, 25, 924, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3024, 3, 25, 26, 927, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3030, 3, 26, 27, 930, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3036, 3, 27, 28, 933, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3042, 0, 3, 888, 2958, 945, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3060, 0, 3, 891, 2964, 954, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3078, 0, 3, 894, 2970, 963, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3096, 0, 3, 897, 2976, 972, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3114, 0, 3, 900, 2982, 981, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3132, 0, 3, 903, 2988, 990, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3150, 0, 3, 906, 2994, 999, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3168, 0, 3, 912, 3000, 1017, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3186, 0, 3, 915, 3006, 1026, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3204, 0, 3, 918, 3012, 1035, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3222, 0, 3, 921, 3018, 1044, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3240, 0, 3, 924, 3024, 1053, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3258, 0, 3, 927, 3030, 1062, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3276, 0, 3, 930, 3036, 1071, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3294, 0, 3, 936, 3042, 96, 102, 1098,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3330, 0, 3, 945, 3060, 102, 108, 1116,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3366, 0, 3, 954, 3078, 108, 114, 1134,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3402, 0, 3, 963, 3096, 114, 120, 1152,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3438, 0, 3, 972, 3114, 120, 126, 1170,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3474, 0, 3, 981, 3132, 126, 132, 1188,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3510, 0, 3, 990, 3150, 132, 138, 1206,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3546, 0, 3, 1008, 3168, 150, 156, 1242,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3582, 0, 3, 1017, 3186, 156, 162, 1260,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3618, 0, 3, 1026, 3204, 162, 168, 1278,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3654, 0, 3, 1035, 3222, 168, 174, 1296,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3690, 0, 3, 1044, 3240, 174, 180, 1314,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3726, 0, 3, 1053, 3258, 180, 186, 1332,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3762, 0, 3, 1062, 3276, 186, 192, 1350,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3798, 0, 3, 1080, 3294, 204, 214, 1368,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3858, 0, 3, 1098, 3330, 214, 224, 1398,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3918, 0, 3, 1116, 3366, 224, 234, 1428,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3978, 0, 3, 1134, 3402, 234, 244, 1458,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4038, 0, 3, 1152, 3438, 244, 254, 1488,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4098, 0, 3, 1170, 3474, 254, 264, 1518,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4158, 0, 3, 1188, 3510, 264, 274, 1548,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4218, 0, 3, 1224, 3546, 294, 304, 1578,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4278, 0, 3, 1242, 3582, 304, 314, 1608,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4338, 0, 3, 1260, 3618, 314, 324, 1638,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4398, 0, 3, 1278, 3654, 324, 334, 1668,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4458, 0, 3, 1296, 3690, 334, 344, 1698,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4518, 0, 3, 1314, 3726, 344, 354, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4578, 0, 3, 1332, 3762, 354, 364, 1758,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4638, 0, 3, 3294, 3330, 1398, 3918, 384,
                                                 399, 1833, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4728, 0, 3, 3330, 3366, 1428, 3978, 399,
                                                 414, 1878, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4818, 0, 3, 3366, 3402, 1458, 4038, 414,
                                                 429, 1923, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4908, 0, 3, 3402, 3438, 1488, 4098, 429,
                                                 444, 1968, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4998, 0, 3, 3438, 3474, 1518, 4158, 444,
                                                 459, 2013, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5088, 0, 3, 3546, 3582, 1608, 4338, 489,
                                                 504, 2103, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5178, 0, 3, 3582, 3618, 1638, 4398, 504,
                                                 519, 2148, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5268, 0, 3, 3618, 3654, 1668, 4458, 519,
                                                 534, 2193, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5358, 0, 3, 3654, 3690, 1698, 4518, 534,
                                                 549, 2238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5448, 0, 3, 3690, 3726, 1728, 4578, 549,
                                                 564, 2283, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5538, 0, 3, 3798, 3858, 1788, 4638, 594,
                                                 615, 2328, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5664, 0, 3, 3858, 3918, 1833, 4728, 615,
                                                 636, 2391, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5790, 0, 3, 3918, 3978, 1878, 4818, 636,
                                                 657, 2454, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5916, 0, 3, 3978, 4038, 1923, 4908, 657,
                                                 678, 2517, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6042, 0, 3, 4038, 4098, 1968, 4998, 678,
                                                 699, 2580, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6168, 0, 3, 4218, 4278, 2058, 5088, 741,
                                                 762, 2643, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6294, 0, 3, 4278, 4338, 2103, 5178, 762,
                                                 783, 2706, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6420, 0, 3, 4338, 4398, 2148, 5268, 783,
                                                 804, 2769, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6546, 0, 3, 4398, 4458, 2193, 5358, 804,
                                                 825, 2832, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6672, 0, 3, 4458, 4518, 2238, 5448, 825,
                                                 846, 2895, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6798, 3, 888, 891, 2964, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6808, 3, 891, 894, 2970, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6818, 3, 894, 897, 2976, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6828, 3, 897, 900, 2982, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6838, 3, 900, 903, 2988, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6848, 3, 903, 906, 2994, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6858, 3, 912, 915, 3006, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6868, 3, 915, 918, 3012, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6878, 3, 918, 921, 3018, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6888, 3, 921, 924, 3024, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6898, 3, 924, 927, 3030, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 6908, 3, 927, 930, 3036, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 6918, 0, 3, 2958, 6798, 3060, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6948, 0, 3, 2964, 6808, 3078, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6978, 0, 3, 2970, 6818, 3096, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7008, 0, 3, 2976, 6828, 3114, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7038, 0, 3, 2982, 6838, 3132, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7068, 0, 3, 2988, 6848, 3150, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7098, 0, 3, 3000, 6858, 3186, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7128, 0, 3, 3006, 6868, 3204, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7158, 0, 3, 3012, 6878, 3222, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7188, 0, 3, 3018, 6888, 3240, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7218, 0, 3, 3024, 6898, 3258, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7248, 0, 3, 3030, 6908, 3276, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 7278, 0, 3, 3042, 6918, 1080, 1098,
                                                 3330, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7338, 0, 3, 3060, 6948, 1098, 1116,
                                                 3366, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7398, 0, 3, 3078, 6978, 1116, 1134,
                                                 3402, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7458, 0, 3, 3096, 7008, 1134, 1152,
                                                 3438, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7518, 0, 3, 3114, 7038, 1152, 1170,
                                                 3474, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7578, 0, 3, 3132, 7068, 1170, 1188,
                                                 3510, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7638, 0, 3, 3168, 7098, 1224, 1242,
                                                 3582, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7698, 0, 3, 3186, 7128, 1242, 1260,
                                                 3618, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7758, 0, 3, 3204, 7158, 1260, 1278,
                                                 3654, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7818, 0, 3, 3222, 7188, 1278, 1296,
                                                 3690, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7878, 0, 3, 3240, 7218, 1296, 1314,
                                                 3726, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7938, 0, 3, 3258, 7248, 1314, 1332,
                                                 3762, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7998, 0, 3, 3330, 7338, 1368, 1398,
                                                 3918, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8098, 0, 3, 3366, 7398, 1398, 1428,
                                                 3978, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8198, 0, 3, 3402, 7458, 1428, 1458,
                                                 4038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8298, 0, 3, 3438, 7518, 1458, 1488,
                                                 4098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8398, 0, 3, 3474, 7578, 1488, 1518,
                                                 4158, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8498, 0, 3, 3582, 7698, 1578, 1608,
                                                 4338, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8598, 0, 3, 3618, 7758, 1608, 1638,
                                                 4398, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8698, 0, 3, 3654, 7818, 1638, 1668,
                                                 4458, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8798, 0, 3, 3690, 7878, 1668, 1698,
                                                 4518, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8898, 0, 3, 3726, 7938, 1698, 1728,
                                                 4578, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8998, 0, 3, 7278, 7338, 3918, 8098,
                                                 1788, 1833, 4728, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9148, 0, 3, 7338, 7398, 3978, 8198,
                                                 1833, 1878, 4818, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9298, 0, 3, 7398, 7458, 4038, 8298,
                                                 1878, 1923, 4908, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9448, 0, 3, 7458, 7518, 4098, 8398,
                                                 1923, 1968, 4998, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9598, 0, 3, 7638, 7698, 4338, 8598,
                                                 2058, 2103, 5178, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9748, 0, 3, 7698, 7758, 4398, 8698,
                                                 2103, 2148, 5268, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9898, 0, 3, 7758, 7818, 4458, 8798,
                                                 2148, 2193, 5358, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10048, 0, 3, 7818, 7878, 4518, 8898,
                                                 2193, 2238, 5448, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10198, 0, 3, 7998, 8098, 4728, 9148,
                                                 2328, 2391, 5790, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10408, 0, 3, 8098, 8198, 4818, 9298,
                                                 2391, 2454, 5916, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10618, 0, 3, 8198, 8298, 4908, 9448,
                                                 2454, 2517, 6042, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10828, 0, 3, 8498, 8598, 5178, 9748,
                                                 2643, 2706, 6420, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11038, 0, 3, 8598, 8698, 5268, 9898,
                                                 2706, 2769, 6546, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11248, 0, 3, 8698, 8798, 5358, 10048,
                                                 2769, 2832, 6672, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11458, 3, 2958, 2964, 6808, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11473, 3, 2964, 2970, 6818, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11488, 3, 2970, 2976, 6828, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11503, 3, 2976, 2982, 6838, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11518, 3, 2982, 2988, 6848, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11533, 3, 3000, 3006, 6868, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11548, 3, 3006, 3012, 6878, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11563, 3, 3012, 3018, 6888, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11578, 3, 3018, 3024, 6898, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 11593, 3, 3024, 3030, 6908, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 11608, 0, 3, 6798, 11458, 6948, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11653, 0, 3, 6808, 11473, 6978, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11698, 0, 3, 6818, 11488, 7008, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11743, 0, 3, 6828, 11503, 7038, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11788, 0, 3, 6838, 11518, 7068, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11833, 0, 3, 6858, 11533, 7128, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11878, 0, 3, 6868, 11548, 7158, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11923, 0, 3, 6878, 11563, 7188, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 11968, 0, 3, 6888, 11578, 7218, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12013, 0, 3, 6898, 11593, 7248, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 12058, 0, 3, 6918, 11608, 3294, 3330,
                                                 7338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12148, 0, 3, 6948, 11653, 3330, 3366,
                                                 7398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12238, 0, 3, 6978, 11698, 3366, 3402,
                                                 7458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12328, 0, 3, 7008, 11743, 3402, 3438,
                                                 7518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12418, 0, 3, 7038, 11788, 3438, 3474,
                                                 7578, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12508, 0, 3, 7098, 11833, 3546, 3582,
                                                 7698, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12598, 0, 3, 7128, 11878, 3582, 3618,
                                                 7758, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12688, 0, 3, 7158, 11923, 3618, 3654,
                                                 7818, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12778, 0, 3, 7188, 11968, 3654, 3690,
                                                 7878, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 12868, 0, 3, 7218, 12013, 3690, 3726,
                                                 7938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12958, 0, 3, 7278, 12058, 3798, 3858,
                                                 7998, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13108, 0, 3, 7338, 12148, 3858, 3918,
                                                 8098, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13258, 0, 3, 7398, 12238, 3918, 3978,
                                                 8198, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13408, 0, 3, 7458, 12328, 3978, 4038,
                                                 8298, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13558, 0, 3, 7518, 12418, 4038, 4098,
                                                 8398, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13708, 0, 3, 7638, 12508, 4218, 4278,
                                                 8498, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13858, 0, 3, 7698, 12598, 4278, 4338,
                                                 8598, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14008, 0, 3, 7758, 12688, 4338, 4398,
                                                 8698, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14158, 0, 3, 7818, 12778, 4398, 4458,
                                                 8798, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14308, 0, 3, 7878, 12868, 4458, 4518,
                                                 8898, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14458, 0, 3, 12058, 12148, 8098, 13258,
                                                 4638, 4728, 9148, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14683, 0, 3, 12148, 12238, 8198, 13408,
                                                 4728, 4818, 9298, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14908, 0, 3, 12238, 12328, 8298, 13558,
                                                 4818, 4908, 9448, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15133, 0, 3, 12508, 12598, 8598, 14008,
                                                 5088, 5178, 9748, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15358, 0, 3, 12598, 12688, 8698, 14158,
                                                 5178, 5268, 9898, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15583, 0, 3, 12688, 12778, 8798, 14308,
                                                 5268, 5358, 10048, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 15808, 0, 3, 12958, 13108, 8998, 14458,
                                                 5538, 5664, 10198, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16123, 0, 3, 13108, 13258, 9148, 14683,
                                                 5664, 5790, 10408, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16438, 0, 3, 13258, 13408, 9298, 14908,
                                                 5790, 5916, 10618, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16753, 0, 3, 13708, 13858, 9598, 15133,
                                                 6168, 6294, 10828, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17068, 0, 3, 13858, 14008, 9748, 15358,
                                                 6294, 6420, 11038, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17383, 0, 3, 14008, 14158, 9898, 15583,
                                                 6420, 6546, 11248, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17698, 3, 6798, 6808, 11473, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17719, 3, 6808, 6818, 11488, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17740, 3, 6818, 6828, 11503, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17761, 3, 6828, 6838, 11518, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17782, 3, 6858, 6868, 11548, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17803, 3, 6868, 6878, 11563, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17824, 3, 6878, 6888, 11578, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 17845, 3, 6888, 6898, 11593, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 17866, 0, 3, 11458, 17698, 11653, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17929, 0, 3, 11473, 17719, 11698, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 17992, 0, 3, 11488, 17740, 11743, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 18055, 0, 3, 11503, 17761, 11788, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 18118, 0, 3, 11533, 17782, 11878, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 18181, 0, 3, 11548, 17803, 11923, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 18244, 0, 3, 11563, 17824, 11968, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 18307, 0, 3, 11578, 17845, 12013, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 18370, 0, 3, 11608, 17866, 7278, 7338,
                                                 12148, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 18496, 0, 3, 11653, 17929, 7338, 7398,
                                                 12238, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 18622, 0, 3, 11698, 17992, 7398, 7458,
                                                 12328, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 18748, 0, 3, 11743, 18055, 7458, 7518,
                                                 12418, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 18874, 0, 3, 11833, 18118, 7638, 7698,
                                                 12598, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 19000, 0, 3, 11878, 18181, 7698, 7758,
                                                 12688, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 19126, 0, 3, 11923, 18244, 7758, 7818,
                                                 12778, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 19252, 0, 3, 11968, 18307, 7818, 7878,
                                                 12868, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19378, 0, 3, 12148, 18496, 7998, 8098,
                                                 13258, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19588, 0, 3, 12238, 18622, 8098, 8198,
                                                 13408, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19798, 0, 3, 12328, 18748, 8198, 8298,
                                                 13558, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 20008, 0, 3, 12598, 19000, 8498, 8598,
                                                 14008, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 20218, 0, 3, 12688, 19126, 8598, 8698,
                                                 14158, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 20428, 0, 3, 12778, 19252, 8698, 8798,
                                                 14308, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20638, 0, 3, 18370, 18496, 13258, 19588,
                                                 8998, 9148, 14683, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20953, 0, 3, 18496, 18622, 13408, 19798,
                                                 9148, 9298, 14908, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 21268, 0, 3, 18874, 19000, 14008, 20218,
                                                 9598, 9748, 15358, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 21583, 0, 3, 19000, 19126, 14158, 20428,
                                                 9748, 9898, 15583, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 21898, 0, 3, 19378, 19588, 14683, 20953,
                                                 10198, 10408, 16438, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 22339, 0, 3, 20008, 20218, 15358, 21583,
                                                 10828, 11038, 17383, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22780, 3, 11458, 11473, 17719, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22808, 3, 11473, 11488, 17740, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22836, 3, 11488, 11503, 17761, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22864, 3, 11533, 11548, 17803, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22892, 3, 11548, 11563, 17824, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 22920, 3, 11563, 11578, 17845, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 22948, 0, 3, 17698, 22780, 17929, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 23032, 0, 3, 17719, 22808, 17992, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 23116, 0, 3, 17740, 22836, 18055, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 23200, 0, 3, 17782, 22864, 18181, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 23284, 0, 3, 17803, 22892, 18244, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 23368, 0, 3, 17824, 22920, 18307, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 23452, 0, 3, 17866, 22948, 12058, 12148,
                                                 18496, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23620, 0, 3, 17929, 23032, 12148, 12238,
                                                 18622, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23788, 0, 3, 17992, 23116, 12238, 12328,
                                                 18748, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23956, 0, 3, 18118, 23200, 12508, 12598,
                                                 19000, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 24124, 0, 3, 18181, 23284, 12598, 12688,
                                                 19126, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 24292, 0, 3, 18244, 23368, 12688, 12778,
                                                 19252, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 24460, 0, 3, 18370, 23452, 12958, 13108,
                                                 19378, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 24740, 0, 3, 18496, 23620, 13108, 13258,
                                                 19588, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25020, 0, 3, 18622, 23788, 13258, 13408,
                                                 19798, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25300, 0, 3, 18874, 23956, 13708, 13858,
                                                 20008, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25580, 0, 3, 19000, 24124, 13858, 14008,
                                                 20218, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25860, 0, 3, 19126, 24292, 14008, 14158,
                                                 20428, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 26140, 0, 3, 23452, 23620, 19588, 25020,
                                                 14458, 14683, 20953, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 26560, 0, 3, 23956, 24124, 20218, 25860,
                                                 15133, 15358, 21583, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 26980, 0, 3, 24460, 24740, 20638, 26140,
                                                 15808, 16123, 21898, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 27568, 0, 3, 25300, 25580, 21268, 26560,
                                                 16753, 17068, 22339, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 28156, 26980, 1176, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 29332, 28744, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 29332, 13, nmax);

    simdtrf::transform_i_inner(buffer, 29332, 28156, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 143 * nvalues, nvalues, buffer, 29332, 13, nmax);
}

}  // namespace simdt2ceri

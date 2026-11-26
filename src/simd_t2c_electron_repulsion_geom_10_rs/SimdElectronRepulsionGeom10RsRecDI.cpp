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


#include "SimdElectronRepulsionGeom10RsRecDI.hpp"

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
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_di_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_di_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 10408, 9322, 1008, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 25, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 7, 8, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 8, 9, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 9, 10, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 98, 0, 10, 11, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 104, 0, 11, 12, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 110, 0, 12, 13, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 13, 14, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 17, 18, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 18, 19, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 19, 20, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 20, 21, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 21, 22, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 22, 23, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 26, 29, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 29, 32, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 32, 35, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 35, 38, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 38, 41, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 41, 44, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 44, 47, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 53, 56, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 56, 59, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 59, 62, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 62, 65, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 65, 68, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 68, 71, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 71, 74, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 304, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 307, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 310, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 313, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 316, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 319, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 322, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 325, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 328, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 331, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 334, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 337, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 340, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 343, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 346, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 349, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 352, 3, 9, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 361, 3, 10, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 370, 3, 11, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 379, 3, 12, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 388, 3, 13, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 397, 3, 14, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 406, 3, 19, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 415, 3, 20, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 424, 3, 21, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 433, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 442, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 451, 3, 24, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 460, 0, 3, 32, 352, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 478, 0, 3, 35, 361, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 496, 0, 3, 38, 370, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 514, 0, 3, 41, 379, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 532, 0, 3, 44, 388, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 550, 0, 3, 47, 397, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 568, 0, 3, 59, 406, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 586, 0, 3, 62, 415, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 604, 0, 3, 65, 424, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 622, 0, 3, 68, 433, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 640, 0, 3, 71, 442, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 658, 0, 3, 74, 451, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 676, 0, 3, 86, 478, 184, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 706, 0, 3, 92, 496, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 736, 0, 3, 98, 514, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 766, 0, 3, 104, 532, 214, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 796, 0, 3, 110, 550, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 826, 0, 3, 128, 586, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 856, 0, 3, 134, 604, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 886, 0, 3, 140, 622, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 916, 0, 3, 146, 640, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 946, 0, 3, 152, 658, 294, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 976, 3, 7, 8, 307, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 982, 3, 8, 9, 310, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 988, 3, 9, 10, 313, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 994, 3, 10, 11, 316, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1000, 3, 11, 12, 319, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1006, 3, 12, 13, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1012, 3, 13, 14, 325, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1018, 3, 17, 18, 331, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1024, 3, 18, 19, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1030, 3, 19, 20, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1036, 3, 20, 21, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1042, 3, 21, 22, 343, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1048, 3, 22, 23, 346, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1054, 3, 23, 24, 349, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1060, 0, 3, 310, 988, 361, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1078, 0, 3, 313, 994, 370, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1096, 0, 3, 316, 1000, 379, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1114, 0, 3, 319, 1006, 388, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1132, 0, 3, 322, 1012, 397, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1150, 0, 3, 334, 1030, 415, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1168, 0, 3, 337, 1036, 424, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1186, 0, 3, 340, 1042, 433, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1204, 0, 3, 343, 1048, 442, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1222, 0, 3, 346, 1054, 451, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1240, 0, 3, 352, 1060, 80, 86, 478,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1276, 0, 3, 361, 1078, 86, 92, 496,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1312, 0, 3, 370, 1096, 92, 98, 514,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1348, 0, 3, 379, 1114, 98, 104, 532,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1384, 0, 3, 388, 1132, 104, 110, 550,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1420, 0, 3, 406, 1150, 122, 128, 586,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1456, 0, 3, 415, 1168, 128, 134, 604,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1492, 0, 3, 424, 1186, 134, 140, 622,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1528, 0, 3, 433, 1204, 140, 146, 640,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1564, 0, 3, 442, 1222, 146, 152, 658,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1600, 0, 3, 460, 1240, 164, 174, 676,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1660, 0, 3, 478, 1276, 174, 184, 706,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1720, 0, 3, 496, 1312, 184, 194, 736,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1780, 0, 3, 514, 1348, 194, 204, 766,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1840, 0, 3, 532, 1384, 204, 214, 796,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1900, 0, 3, 568, 1420, 234, 244, 826,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1960, 0, 3, 586, 1456, 244, 254, 856,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2020, 0, 3, 604, 1492, 254, 264, 886,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2080, 0, 3, 622, 1528, 264, 274, 916,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2140, 0, 3, 640, 1564, 274, 284, 946,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2200, 3, 304, 307, 982, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2210, 3, 307, 310, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2220, 3, 310, 313, 994, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2230, 3, 313, 316, 1000, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2240, 3, 316, 319, 1006, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2250, 3, 319, 322, 1012, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2260, 3, 328, 331, 1024, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2270, 3, 331, 334, 1030, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2280, 3, 334, 337, 1036, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2290, 3, 337, 340, 1042, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2300, 3, 340, 343, 1048, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2310, 3, 343, 346, 1054, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2320, 0, 3, 988, 2220, 1078, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2350, 0, 3, 994, 2230, 1096, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2380, 0, 3, 1000, 2240, 1114, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2410, 0, 3, 1006, 2250, 1132, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2440, 0, 3, 1030, 2280, 1168, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2470, 0, 3, 1036, 2290, 1186, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2500, 0, 3, 1042, 2300, 1204, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2530, 0, 3, 1048, 2310, 1222, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2560, 0, 3, 1060, 2320, 460, 478, 1276,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2620, 0, 3, 1078, 2350, 478, 496, 1312,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2680, 0, 3, 1096, 2380, 496, 514, 1348,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2740, 0, 3, 1114, 2410, 514, 532, 1384,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2800, 0, 3, 1150, 2440, 568, 586, 1456,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2860, 0, 3, 1168, 2470, 586, 604, 1492,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2920, 0, 3, 1186, 2500, 604, 622, 1528,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2980, 0, 3, 1204, 2530, 622, 640, 1564,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3040, 0, 3, 1276, 2620, 676, 706, 1720,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3140, 0, 3, 1312, 2680, 706, 736, 1780,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3240, 0, 3, 1348, 2740, 736, 766, 1840,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3340, 0, 3, 1456, 2860, 826, 856, 2020,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3440, 0, 3, 1492, 2920, 856, 886, 2080,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3540, 0, 3, 1528, 2980, 886, 916, 2140,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3640, 3, 976, 982, 2210, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3655, 3, 982, 988, 2220, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3670, 3, 988, 994, 2230, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3685, 3, 994, 1000, 2240, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3700, 3, 1000, 1006, 2250, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3715, 3, 1018, 1024, 2270, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3730, 3, 1024, 1030, 2280, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3745, 3, 1030, 1036, 2290, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3760, 3, 1036, 1042, 2300, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3775, 3, 1042, 1048, 2310, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 3790, 0, 3, 2220, 3670, 2350, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3835, 0, 3, 2230, 3685, 2380, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3880, 0, 3, 2240, 3700, 2410, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3925, 0, 3, 2280, 3745, 2470, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3970, 0, 3, 2290, 3760, 2500, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4015, 0, 3, 2300, 3775, 2530, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 4060, 0, 3, 2320, 3790, 1240, 1276,
                                                 2620, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4150, 0, 3, 2350, 3835, 1276, 1312,
                                                 2680, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4240, 0, 3, 2380, 3880, 1312, 1348,
                                                 2740, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4330, 0, 3, 2440, 3925, 1420, 1456,
                                                 2860, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4420, 0, 3, 2470, 3970, 1456, 1492,
                                                 2920, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4510, 0, 3, 2500, 4015, 1492, 1528,
                                                 2980, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4600, 0, 3, 2560, 4060, 1600, 1660,
                                                 3040, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4750, 0, 3, 2620, 4150, 1660, 1720,
                                                 3140, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4900, 0, 3, 2680, 4240, 1720, 1780,
                                                 3240, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5050, 0, 3, 2800, 4330, 1900, 1960,
                                                 3340, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5200, 0, 3, 2860, 4420, 1960, 2020,
                                                 3440, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5350, 0, 3, 2920, 4510, 2020, 2080,
                                                 3540, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5500, 3, 2200, 2210, 3655, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5521, 3, 2210, 2220, 3670, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5542, 3, 2220, 2230, 3685, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5563, 3, 2230, 2240, 3700, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5584, 3, 2260, 2270, 3730, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5605, 3, 2270, 2280, 3745, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5626, 3, 2280, 2290, 3760, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5647, 3, 2290, 2300, 3775, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 5668, 0, 3, 3655, 5521, 3790, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5731, 0, 3, 3670, 5542, 3835, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5794, 0, 3, 3685, 5563, 3880, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5857, 0, 3, 3730, 5605, 3925, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5920, 0, 3, 3745, 5626, 3970, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5983, 0, 3, 3760, 5647, 4015, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 6046, 0, 3, 3790, 5731, 2560, 2620,
                                                 4150, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6172, 0, 3, 3835, 5794, 2620, 2680,
                                                 4240, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6298, 0, 3, 3925, 5920, 2800, 2860,
                                                 4420, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6424, 0, 3, 3970, 5983, 2860, 2920,
                                                 4510, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6550, 0, 3, 4150, 6172, 3040, 3140,
                                                 4900, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6760, 0, 3, 4420, 6424, 3340, 3440,
                                                 5350, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6970, 3, 3640, 3655, 5521, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6998, 3, 3670, 3685, 5563, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 7026, 3, 3715, 3730, 5605, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 7054, 3, 3745, 3760, 5647, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 7082, 0, 3, 5500, 6970, 5668, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 7166, 0, 3, 5542, 6998, 5794, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 7250, 0, 3, 5584, 7026, 5857, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 7334, 0, 3, 5626, 7054, 5983, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 7418, 0, 3, 5731, 7166, 4060, 4150,
                                                 6172, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7586, 0, 3, 5920, 7334, 4330, 4420,
                                                 6424, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 7754, 0, 3, 6046, 7418, 4600, 4750,
                                                 6550, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 8034, 0, 3, 6298, 7586, 5050, 5200,
                                                 6760, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 8314, 7250, 8034, 1, 28, ncols, alpha);

            simdgeo::geom_d_y(buffer, 8482, 7250, 8034, 1, 28, ncols, alpha);

            simdgeo::geom_d_z(buffer, 8650, 7250, 8034, 1, 28, ncols, alpha);

            simdgeo::geom_d_x(buffer, 8818, 7082, 7754, 1, 28, ncols, alpha);

            simdgeo::geom_d_y(buffer, 8986, 7082, 7754, 1, 28, ncols, alpha);

            simdgeo::geom_d_z(buffer, 9154, 7082, 7754, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 9322, 8818, 504, ncols);

            simdfunc::contract_primitives(buffer, 9826, 8314, 504, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 10330, 9826, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 10330, 13, nmax);

    simdtrf::transform_i_inner(buffer, 10330, 9994, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 65 * nvalues, nvalues, buffer, 10330, 13, nmax);

    simdtrf::transform_i_inner(buffer, 10330, 10162, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 130 * nvalues, nvalues, buffer, 10330, 13, nmax);

    simdtrf::transform_i_inner(buffer, 10330, 9322, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 195 * nvalues, nvalues, buffer, 10330, 13, nmax);

    simdtrf::transform_i_inner(buffer, 10330, 9490, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 260 * nvalues, nvalues, buffer, 10330, 13, nmax);

    simdtrf::transform_i_inner(buffer, 10330, 9658, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 325 * nvalues, nvalues, buffer, 10330, 13, nmax);
}

}  // namespace simdt2ceri

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


#include "SimdElectronRepulsionGeom10RsRecFI.hpp"

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
#include "SimdGeometryF1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_fi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_fi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 20134, 18324, 1680, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 542, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 545, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 548, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 551, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 554, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 557, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 560, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 563, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 566, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 569, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 572, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 575, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 578, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 581, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 584, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 587, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 590, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 599, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 608, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 617, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 626, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 635, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 644, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 653, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 662, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 671, 3, 22, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 680, 3, 23, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 689, 3, 24, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 698, 3, 25, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 707, 3, 26, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 716, 3, 27, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 725, 3, 28, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 734, 0, 3, 33, 599, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 752, 0, 3, 36, 608, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 770, 0, 3, 39, 617, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 788, 0, 3, 42, 626, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 806, 0, 3, 45, 635, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 824, 0, 3, 48, 644, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 842, 0, 3, 51, 653, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 860, 0, 3, 60, 671, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 878, 0, 3, 63, 680, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 896, 0, 3, 66, 689, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 914, 0, 3, 69, 698, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 932, 0, 3, 72, 707, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 950, 0, 3, 75, 716, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 968, 0, 3, 78, 725, 186, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 986, 0, 3, 96, 752, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1016, 0, 3, 102, 770, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1046, 0, 3, 108, 788, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1076, 0, 3, 114, 806, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1106, 0, 3, 120, 824, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1136, 0, 3, 126, 842, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1166, 0, 3, 150, 878, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1196, 0, 3, 156, 896, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1226, 0, 3, 162, 914, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1256, 0, 3, 168, 932, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1286, 0, 3, 174, 950, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1316, 0, 3, 180, 968, 322, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1346, 0, 3, 202, 1016, 362, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1391, 0, 3, 212, 1046, 377, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1436, 0, 3, 222, 1076, 392, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1481, 0, 3, 232, 1106, 407, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1526, 0, 3, 242, 1136, 422, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1571, 0, 3, 272, 1196, 467, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1616, 0, 3, 282, 1226, 482, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1661, 0, 3, 292, 1256, 497, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1706, 0, 3, 302, 1286, 512, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1751, 0, 3, 312, 1316, 527, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1796, 3, 9, 10, 545, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1802, 3, 10, 11, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1808, 3, 11, 12, 551, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1814, 3, 12, 13, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1820, 3, 13, 14, 557, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1826, 3, 14, 15, 560, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1832, 3, 15, 16, 563, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1838, 3, 21, 22, 569, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1844, 3, 22, 23, 572, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1850, 3, 23, 24, 575, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1856, 3, 24, 25, 578, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1862, 3, 25, 26, 581, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1868, 3, 26, 27, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1874, 3, 27, 28, 587, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1880, 0, 3, 542, 1796, 599, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1898, 0, 3, 545, 1802, 608, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1916, 0, 3, 548, 1808, 617, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1934, 0, 3, 551, 1814, 626, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1952, 0, 3, 554, 1820, 635, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1970, 0, 3, 557, 1826, 644, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1988, 0, 3, 560, 1832, 653, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2006, 0, 3, 566, 1838, 671, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2024, 0, 3, 569, 1844, 680, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2042, 0, 3, 572, 1850, 689, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2060, 0, 3, 575, 1856, 698, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2078, 0, 3, 578, 1862, 707, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2096, 0, 3, 581, 1868, 716, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2114, 0, 3, 584, 1874, 725, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2132, 0, 3, 590, 1880, 84, 90, 734,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2168, 0, 3, 599, 1898, 90, 96, 752,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2204, 0, 3, 608, 1916, 96, 102, 770,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2240, 0, 3, 617, 1934, 102, 108, 788,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2276, 0, 3, 626, 1952, 108, 114, 806,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2312, 0, 3, 635, 1970, 114, 120, 824,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2348, 0, 3, 644, 1988, 120, 126, 842,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2384, 0, 3, 662, 2006, 138, 144, 860,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2420, 0, 3, 671, 2024, 144, 150, 878,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2456, 0, 3, 680, 2042, 150, 156, 896,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2492, 0, 3, 689, 2060, 156, 162, 914,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2528, 0, 3, 698, 2078, 162, 168, 932,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2564, 0, 3, 707, 2096, 168, 174, 950,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2600, 0, 3, 716, 2114, 174, 180, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2636, 0, 3, 752, 2204, 192, 202, 1016,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2696, 0, 3, 770, 2240, 202, 212, 1046,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2756, 0, 3, 788, 2276, 212, 222, 1076,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2816, 0, 3, 806, 2312, 222, 232, 1106,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2876, 0, 3, 824, 2348, 232, 242, 1136,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2936, 0, 3, 878, 2456, 262, 272, 1196,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2996, 0, 3, 896, 2492, 272, 282, 1226,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3056, 0, 3, 914, 2528, 282, 292, 1256,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3116, 0, 3, 932, 2564, 292, 302, 1286,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3176, 0, 3, 950, 2600, 302, 312, 1316,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3236, 0, 3, 2132, 2168, 986, 2636, 332,
                                                 347, 1346, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3326, 0, 3, 2168, 2204, 1016, 2696, 347,
                                                 362, 1391, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3416, 0, 3, 2204, 2240, 1046, 2756, 362,
                                                 377, 1436, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3506, 0, 3, 2240, 2276, 1076, 2816, 377,
                                                 392, 1481, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3596, 0, 3, 2276, 2312, 1106, 2876, 392,
                                                 407, 1526, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3686, 0, 3, 2384, 2420, 1166, 2936, 437,
                                                 452, 1571, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3776, 0, 3, 2420, 2456, 1196, 2996, 452,
                                                 467, 1616, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3866, 0, 3, 2456, 2492, 1226, 3056, 467,
                                                 482, 1661, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 3956, 0, 3, 2492, 2528, 1256, 3116, 482,
                                                 497, 1706, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4046, 0, 3, 2528, 2564, 1286, 3176, 497,
                                                 512, 1751, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4136, 3, 542, 545, 1802, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4146, 3, 545, 548, 1808, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4156, 3, 548, 551, 1814, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4166, 3, 551, 554, 1820, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4176, 3, 554, 557, 1826, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4186, 3, 557, 560, 1832, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4196, 3, 566, 569, 1844, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4206, 3, 569, 572, 1850, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4216, 3, 572, 575, 1856, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4226, 3, 575, 578, 1862, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4236, 3, 578, 581, 1868, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 4246, 3, 581, 584, 1874, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 4256, 0, 3, 1796, 4136, 1898, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4286, 0, 3, 1802, 4146, 1916, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4316, 0, 3, 1808, 4156, 1934, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4346, 0, 3, 1814, 4166, 1952, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4376, 0, 3, 1820, 4176, 1970, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4406, 0, 3, 1826, 4186, 1988, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4436, 0, 3, 1838, 4196, 2024, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4466, 0, 3, 1844, 4206, 2042, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4496, 0, 3, 1850, 4216, 2060, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4526, 0, 3, 1856, 4226, 2078, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4556, 0, 3, 1862, 4236, 2096, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 4586, 0, 3, 1868, 4246, 2114, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 4616, 0, 3, 1898, 4286, 734, 752, 2204,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4676, 0, 3, 1916, 4316, 752, 770, 2240,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4736, 0, 3, 1934, 4346, 770, 788, 2276,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4796, 0, 3, 1952, 4376, 788, 806, 2312,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4856, 0, 3, 1970, 4406, 806, 824, 2348,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4916, 0, 3, 2024, 4466, 860, 878, 2456,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 4976, 0, 3, 2042, 4496, 878, 896, 2492,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5036, 0, 3, 2060, 4526, 896, 914, 2528,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5096, 0, 3, 2078, 4556, 914, 932, 2564,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 5156, 0, 3, 2096, 4586, 932, 950, 2600,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5216, 0, 3, 2204, 4676, 986, 1016, 2696,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5316, 0, 3, 2240, 4736, 1016, 1046,
                                                 2756, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5416, 0, 3, 2276, 4796, 1046, 1076,
                                                 2816, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5516, 0, 3, 2312, 4856, 1076, 1106,
                                                 2876, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5616, 0, 3, 2456, 4976, 1166, 1196,
                                                 2996, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5716, 0, 3, 2492, 5036, 1196, 1226,
                                                 3056, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5816, 0, 3, 2528, 5096, 1226, 1256,
                                                 3116, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 5916, 0, 3, 2564, 5156, 1256, 1286,
                                                 3176, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6016, 0, 3, 4616, 4676, 2696, 5316,
                                                 1346, 1391, 3416, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6166, 0, 3, 4676, 4736, 2756, 5416,
                                                 1391, 1436, 3506, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6316, 0, 3, 4736, 4796, 2816, 5516,
                                                 1436, 1481, 3596, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6466, 0, 3, 4916, 4976, 2996, 5716,
                                                 1571, 1616, 3866, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6616, 0, 3, 4976, 5036, 3056, 5816,
                                                 1616, 1661, 3956, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 6766, 0, 3, 5036, 5096, 3116, 5916,
                                                 1661, 1706, 4046, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6916, 3, 1796, 1802, 4146, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6931, 3, 1802, 1808, 4156, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6946, 3, 1808, 1814, 4166, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6961, 3, 1814, 1820, 4176, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6976, 3, 1820, 1826, 4186, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 6991, 3, 1838, 1844, 4206, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7006, 3, 1844, 1850, 4216, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7021, 3, 1850, 1856, 4226, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7036, 3, 1856, 1862, 4236, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 7051, 3, 1862, 1868, 4246, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 7066, 0, 3, 4136, 6916, 4286, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7111, 0, 3, 4146, 6931, 4316, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7156, 0, 3, 4156, 6946, 4346, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7201, 0, 3, 4166, 6961, 4376, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7246, 0, 3, 4176, 6976, 4406, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7291, 0, 3, 4196, 6991, 4466, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7336, 0, 3, 4206, 7006, 4496, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7381, 0, 3, 4216, 7021, 4526, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7426, 0, 3, 4226, 7036, 4556, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 7471, 0, 3, 4236, 7051, 4586, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 7516, 0, 3, 4256, 7066, 2132, 2168,
                                                 4616, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7606, 0, 3, 4286, 7111, 2168, 2204,
                                                 4676, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7696, 0, 3, 4316, 7156, 2204, 2240,
                                                 4736, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7786, 0, 3, 4346, 7201, 2240, 2276,
                                                 4796, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7876, 0, 3, 4376, 7246, 2276, 2312,
                                                 4856, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 7966, 0, 3, 4436, 7291, 2384, 2420,
                                                 4916, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8056, 0, 3, 4466, 7336, 2420, 2456,
                                                 4976, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8146, 0, 3, 4496, 7381, 2456, 2492,
                                                 5036, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8236, 0, 3, 4526, 7426, 2492, 2528,
                                                 5096, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 8326, 0, 3, 4556, 7471, 2528, 2564,
                                                 5156, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8416, 0, 3, 4676, 7696, 2636, 2696,
                                                 5316, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8566, 0, 3, 4736, 7786, 2696, 2756,
                                                 5416, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8716, 0, 3, 4796, 7876, 2756, 2816,
                                                 5516, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 8866, 0, 3, 4976, 8146, 2936, 2996,
                                                 5716, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9016, 0, 3, 5036, 8236, 2996, 3056,
                                                 5816, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 9166, 0, 3, 5096, 8326, 3056, 3116,
                                                 5916, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9316, 0, 3, 7516, 7606, 5216, 8416,
                                                 3236, 3326, 6016, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9541, 0, 3, 7606, 7696, 5316, 8566,
                                                 3326, 3416, 6166, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9766, 0, 3, 7696, 7786, 5416, 8716,
                                                 3416, 3506, 6316, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 9991, 0, 3, 7966, 8056, 5616, 8866,
                                                 3686, 3776, 6466, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10216, 0, 3, 8056, 8146, 5716, 9016,
                                                 3776, 3866, 6616, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 10441, 0, 3, 8146, 8236, 5816, 9166,
                                                 3866, 3956, 6766, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10666, 3, 4136, 4146, 6931, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10687, 3, 4146, 4156, 6946, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10708, 3, 4156, 4166, 6961, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10729, 3, 4166, 4176, 6976, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10750, 3, 4196, 4206, 7006, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10771, 3, 4206, 4216, 7021, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10792, 3, 4216, 4226, 7036, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 10813, 3, 4226, 4236, 7051, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 10834, 0, 3, 6916, 10666, 7111, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10897, 0, 3, 6931, 10687, 7156, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 10960, 0, 3, 6946, 10708, 7201, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11023, 0, 3, 6961, 10729, 7246, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11086, 0, 3, 6991, 10750, 7336, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11149, 0, 3, 7006, 10771, 7381, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11212, 0, 3, 7021, 10792, 7426, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 11275, 0, 3, 7036, 10813, 7471, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 11338, 0, 3, 7111, 10897, 4616, 4676,
                                                 7696, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11464, 0, 3, 7156, 10960, 4676, 4736,
                                                 7786, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11590, 0, 3, 7201, 11023, 4736, 4796,
                                                 7876, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11716, 0, 3, 7336, 11149, 4916, 4976,
                                                 8146, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11842, 0, 3, 7381, 11212, 4976, 5036,
                                                 8236, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 11968, 0, 3, 7426, 11275, 5036, 5096,
                                                 8326, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12094, 0, 3, 7696, 11464, 5216, 5316,
                                                 8566, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12304, 0, 3, 7786, 11590, 5316, 5416,
                                                 8716, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12514, 0, 3, 8146, 11842, 5616, 5716,
                                                 9016, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 12724, 0, 3, 8236, 11968, 5716, 5816,
                                                 9166, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 12934, 0, 3, 11338, 11464, 8566, 12304,
                                                 6016, 6166, 9766, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 13249, 0, 3, 11716, 11842, 9016, 12724,
                                                 6466, 6616, 10441, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13564, 3, 6916, 6931, 10687, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13592, 3, 6931, 6946, 10708, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13620, 3, 6946, 6961, 10729, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13648, 3, 6991, 7006, 10771, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13676, 3, 7006, 7021, 10792, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 13704, 3, 7021, 7036, 10813, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 13732, 0, 3, 10666, 13564, 10897, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 13816, 0, 3, 10687, 13592, 10960, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 13900, 0, 3, 10708, 13620, 11023, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 13984, 0, 3, 10750, 13648, 11149, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 14068, 0, 3, 10771, 13676, 11212, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 14152, 0, 3, 10792, 13704, 11275, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 14236, 0, 3, 10834, 13732, 7516, 7606,
                                                 11338, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 14404, 0, 3, 10897, 13816, 7606, 7696,
                                                 11464, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 14572, 0, 3, 10960, 13900, 7696, 7786,
                                                 11590, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 14740, 0, 3, 11086, 13984, 7966, 8056,
                                                 11716, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 14908, 0, 3, 11149, 14068, 8056, 8146,
                                                 11842, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 15076, 0, 3, 11212, 14152, 8146, 8236,
                                                 11968, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 15244, 0, 3, 11464, 14572, 8416, 8566,
                                                 12304, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 15524, 0, 3, 11842, 15076, 8866, 9016,
                                                 12724, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 15804, 0, 3, 14236, 14404, 12094, 15244,
                                                 9316, 9541, 12934, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 16224, 0, 3, 14740, 14908, 12514, 15524,
                                                 9991, 10216, 13249, ncols, alpha, beta, p);

            simdgeo::geom_f_x(buffer, 16644, 14740, 16224, 1, 28, ncols, alpha);

            simdgeo::geom_f_y(buffer, 16924, 14740, 16224, 1, 28, ncols, alpha);

            simdgeo::geom_f_z(buffer, 17204, 14740, 16224, 1, 28, ncols, alpha);

            simdgeo::geom_f_x(buffer, 17484, 14236, 15804, 1, 28, ncols, alpha);

            simdgeo::geom_f_y(buffer, 17764, 14236, 15804, 1, 28, ncols, alpha);

            simdgeo::geom_f_z(buffer, 18044, 14236, 15804, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 18324, 17484, 840, ncols);

            simdfunc::contract_primitives(buffer, 19164, 16644, 840, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 20004, 19164, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 20004, 13, nmax);

    simdtrf::transform_i_inner(buffer, 20004, 19444, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 91 * nvalues, nvalues, buffer, 20004, 13, nmax);

    simdtrf::transform_i_inner(buffer, 20004, 19724, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 182 * nvalues, nvalues, buffer, 20004, 13, nmax);

    simdtrf::transform_i_inner(buffer, 20004, 18324, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 273 * nvalues, nvalues, buffer, 20004, 13, nmax);

    simdtrf::transform_i_inner(buffer, 20004, 18604, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 364 * nvalues, nvalues, buffer, 20004, 13, nmax);

    simdtrf::transform_i_inner(buffer, 20004, 18884, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 455 * nvalues, nvalues, buffer, 20004, 13, nmax);
}

}  // namespace simdt2ceri

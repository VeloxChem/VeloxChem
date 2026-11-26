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


#include "SimdElectronRepulsionGeom10RsRecDH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
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
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_dh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 6590, 5768, 756, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 15, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 7, 8, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 8, 9, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 9, 10, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 10, 11, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 11, 12, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 12, 13, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 16, 17, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 17, 18, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 18, 19, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 45, 48, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 48, 51, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 51, 54, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 54, 57, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 57, 60, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 238, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 241, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 244, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 247, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 250, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 253, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 256, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 259, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 262, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 265, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 268, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 271, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 274, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 277, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 280, 3, 9, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 289, 3, 10, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 298, 3, 11, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 307, 3, 12, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 316, 3, 13, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 325, 3, 18, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 334, 3, 19, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 343, 3, 20, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 352, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 361, 3, 22, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 370, 0, 3, 27, 280, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 388, 0, 3, 30, 289, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 406, 0, 3, 33, 298, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 424, 0, 3, 36, 307, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 442, 0, 3, 39, 316, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 460, 0, 3, 48, 325, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 478, 0, 3, 51, 334, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 496, 0, 3, 54, 343, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 514, 0, 3, 57, 352, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 532, 0, 3, 60, 361, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 550, 0, 3, 66, 370, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 580, 0, 3, 72, 388, 148, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 610, 0, 3, 78, 406, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 640, 0, 3, 84, 424, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 670, 0, 3, 90, 442, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 700, 0, 3, 102, 460, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 730, 0, 3, 108, 478, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 760, 0, 3, 114, 496, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 790, 0, 3, 120, 514, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 820, 0, 3, 126, 532, 228, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 850, 3, 7, 8, 241, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 856, 3, 8, 9, 244, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 862, 3, 9, 10, 247, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 868, 3, 10, 11, 250, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 874, 3, 11, 12, 253, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 880, 3, 12, 13, 256, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 886, 3, 16, 17, 262, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 892, 3, 17, 18, 265, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 898, 3, 18, 19, 268, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 904, 3, 19, 20, 271, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 910, 3, 20, 21, 274, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 916, 3, 21, 22, 277, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 922, 0, 3, 244, 862, 289, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 940, 0, 3, 247, 868, 298, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 958, 0, 3, 250, 874, 307, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 976, 0, 3, 253, 880, 316, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 994, 0, 3, 265, 898, 334, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1012, 0, 3, 268, 904, 343, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1030, 0, 3, 271, 910, 352, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1048, 0, 3, 274, 916, 361, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1066, 0, 3, 280, 922, 66, 72, 388,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1102, 0, 3, 289, 940, 72, 78, 406,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1138, 0, 3, 298, 958, 78, 84, 424,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1174, 0, 3, 307, 976, 84, 90, 442,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1210, 0, 3, 325, 994, 102, 108, 478,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1246, 0, 3, 334, 1012, 108, 114, 496,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1282, 0, 3, 343, 1030, 114, 120, 514,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1318, 0, 3, 352, 1048, 120, 126, 532,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1354, 0, 3, 388, 1102, 138, 148, 610,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1414, 0, 3, 406, 1138, 148, 158, 640,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1474, 0, 3, 424, 1174, 158, 168, 670,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1534, 0, 3, 478, 1246, 188, 198, 760,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1594, 0, 3, 496, 1282, 198, 208, 790,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1654, 0, 3, 514, 1318, 208, 218, 820,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1714, 3, 238, 241, 856, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1724, 3, 241, 244, 862, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1734, 3, 244, 247, 868, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1744, 3, 247, 250, 874, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1754, 3, 250, 253, 880, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1764, 3, 259, 262, 892, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1774, 3, 262, 265, 898, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1784, 3, 265, 268, 904, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1794, 3, 268, 271, 910, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1804, 3, 271, 274, 916, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1814, 0, 3, 862, 1734, 940, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1844, 0, 3, 868, 1744, 958, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1874, 0, 3, 874, 1754, 976, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1904, 0, 3, 898, 1784, 1012, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1934, 0, 3, 904, 1794, 1030, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1964, 0, 3, 910, 1804, 1048, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1994, 0, 3, 922, 1814, 370, 388, 1102,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2054, 0, 3, 940, 1844, 388, 406, 1138,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2114, 0, 3, 958, 1874, 406, 424, 1174,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2174, 0, 3, 994, 1904, 460, 478, 1246,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2234, 0, 3, 1012, 1934, 478, 496, 1282,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2294, 0, 3, 1030, 1964, 496, 514, 1318,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2354, 0, 3, 1066, 1994, 550, 580, 1354,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2454, 0, 3, 1102, 2054, 580, 610, 1414,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2554, 0, 3, 1138, 2114, 610, 640, 1474,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2654, 0, 3, 1210, 2174, 700, 730, 1534,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2754, 0, 3, 1246, 2234, 730, 760, 1594,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2854, 0, 3, 1282, 2294, 760, 790, 1654,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2954, 3, 850, 856, 1724, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2969, 3, 856, 862, 1734, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2984, 3, 862, 868, 1744, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2999, 3, 868, 874, 1754, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3014, 3, 886, 892, 1774, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3029, 3, 892, 898, 1784, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3044, 3, 898, 904, 1794, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3059, 3, 904, 910, 1804, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 3074, 0, 3, 1724, 2969, 1814, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3119, 0, 3, 1734, 2984, 1844, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3164, 0, 3, 1744, 2999, 1874, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3209, 0, 3, 1774, 3029, 1904, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3254, 0, 3, 1784, 3044, 1934, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3299, 0, 3, 1794, 3059, 1964, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3344, 0, 3, 1814, 3119, 1066, 1102,
                                                 2054, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3434, 0, 3, 1844, 3164, 1102, 1138,
                                                 2114, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3524, 0, 3, 1904, 3254, 1210, 1246,
                                                 2234, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3614, 0, 3, 1934, 3299, 1246, 1282,
                                                 2294, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3704, 0, 3, 2054, 3434, 1354, 1414,
                                                 2554, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3854, 0, 3, 2234, 3614, 1534, 1594,
                                                 2854, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4004, 3, 1714, 1724, 2969, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4025, 3, 1734, 1744, 2999, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4046, 3, 1764, 1774, 3029, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4067, 3, 1784, 1794, 3059, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 4088, 0, 3, 2954, 4004, 3074, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4151, 0, 3, 2984, 4025, 3164, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4214, 0, 3, 3014, 4046, 3209, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4277, 0, 3, 3044, 4067, 3299, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4340, 0, 3, 3119, 4151, 1994, 2054,
                                                 3434, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4466, 0, 3, 3254, 4277, 2174, 2234,
                                                 3614, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4592, 0, 3, 3344, 4340, 2354, 2454,
                                                 3704, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4802, 0, 3, 3524, 4466, 2654, 2754,
                                                 3854, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 5012, 4214, 4802, 1, 21, ncols, alpha);

            simdgeo::geom_d_y(buffer, 5138, 4214, 4802, 1, 21, ncols, alpha);

            simdgeo::geom_d_z(buffer, 5264, 4214, 4802, 1, 21, ncols, alpha);

            simdgeo::geom_d_x(buffer, 5390, 4088, 4592, 1, 21, ncols, alpha);

            simdgeo::geom_d_y(buffer, 5516, 4088, 4592, 1, 21, ncols, alpha);

            simdgeo::geom_d_z(buffer, 5642, 4088, 4592, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 5768, 5390, 378, ncols);

            simdfunc::contract_primitives(buffer, 6146, 5012, 378, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 6524, 6146, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 6524, 11, nmax);

    simdtrf::transform_h_inner(buffer, 6524, 6272, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 55 * nvalues, nvalues, buffer, 6524, 11, nmax);

    simdtrf::transform_h_inner(buffer, 6524, 6398, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 110 * nvalues, nvalues, buffer, 6524, 11, nmax);

    simdtrf::transform_h_inner(buffer, 6524, 5768, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 165 * nvalues, nvalues, buffer, 6524, 11, nmax);

    simdtrf::transform_h_inner(buffer, 6524, 5894, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 220 * nvalues, nvalues, buffer, 6524, 11, nmax);

    simdtrf::transform_h_inner(buffer, 6524, 6020, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 275 * nvalues, nvalues, buffer, 6524, 11, nmax);
}

}  // namespace simdt2ceri

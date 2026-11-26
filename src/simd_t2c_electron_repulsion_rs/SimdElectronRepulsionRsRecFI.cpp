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


#include "SimdElectronRepulsionRsRecFI.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_fi_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 8434, 7744, 560, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 304, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 307, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 310, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 313, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 316, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 319, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 322, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 325, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 328, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 331, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 334, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 337, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 340, 3, 9, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 349, 3, 10, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 358, 3, 11, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 367, 3, 12, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 376, 3, 13, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 385, 3, 14, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 394, 3, 19, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 403, 3, 20, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 412, 3, 21, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 421, 3, 22, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 430, 3, 23, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 439, 3, 24, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 448, 0, 3, 32, 340, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 466, 0, 3, 35, 349, 92, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 484, 0, 3, 38, 358, 98, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 502, 0, 3, 41, 367, 104, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 520, 0, 3, 44, 376, 110, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 538, 0, 3, 47, 385, 116, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 556, 0, 3, 59, 394, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 574, 0, 3, 62, 403, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 592, 0, 3, 65, 412, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 610, 0, 3, 68, 421, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 628, 0, 3, 71, 430, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 646, 0, 3, 74, 439, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 664, 0, 3, 86, 466, 184, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 694, 0, 3, 92, 484, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 724, 0, 3, 98, 502, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 754, 0, 3, 104, 520, 214, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 784, 0, 3, 110, 538, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 814, 0, 3, 128, 574, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 844, 0, 3, 134, 592, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 874, 0, 3, 140, 610, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 904, 0, 3, 146, 628, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 934, 0, 3, 152, 646, 294, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 964, 3, 9, 10, 307, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 970, 3, 10, 11, 310, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 976, 3, 11, 12, 313, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 982, 3, 12, 13, 316, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 988, 3, 13, 14, 319, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 994, 3, 19, 20, 325, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1000, 3, 20, 21, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1006, 3, 21, 22, 331, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1012, 3, 22, 23, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1018, 3, 23, 24, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1024, 0, 3, 304, 964, 349, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1042, 0, 3, 307, 970, 358, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1060, 0, 3, 310, 976, 367, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1078, 0, 3, 313, 982, 376, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1096, 0, 3, 316, 988, 385, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1114, 0, 3, 322, 994, 403, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1132, 0, 3, 325, 1000, 412, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1150, 0, 3, 328, 1006, 421, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1168, 0, 3, 331, 1012, 430, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1186, 0, 3, 334, 1018, 439, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1204, 0, 3, 340, 1024, 80, 86, 466,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1240, 0, 3, 349, 1042, 86, 92, 484,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1276, 0, 3, 358, 1060, 92, 98, 502,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1312, 0, 3, 367, 1078, 98, 104, 520,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1348, 0, 3, 376, 1096, 104, 110, 538,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1384, 0, 3, 394, 1114, 122, 128, 574,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1420, 0, 3, 403, 1132, 128, 134, 592,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1456, 0, 3, 412, 1150, 134, 140, 610,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1492, 0, 3, 421, 1168, 140, 146, 628,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1528, 0, 3, 430, 1186, 146, 152, 646,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1564, 0, 3, 448, 1204, 164, 174, 664,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1624, 0, 3, 466, 1240, 174, 184, 694,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1684, 0, 3, 484, 1276, 184, 194, 724,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1744, 0, 3, 502, 1312, 194, 204, 754,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1804, 0, 3, 520, 1348, 204, 214, 784,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1864, 0, 3, 556, 1384, 234, 244, 814,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1924, 0, 3, 574, 1420, 244, 254, 844,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1984, 0, 3, 592, 1456, 254, 264, 874,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2044, 0, 3, 610, 1492, 264, 274, 904,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2104, 0, 3, 628, 1528, 274, 284, 934,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2164, 3, 304, 307, 970, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2174, 3, 307, 310, 976, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2184, 3, 310, 313, 982, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2194, 3, 313, 316, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2204, 3, 322, 325, 1000, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2214, 3, 325, 328, 1006, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2224, 3, 328, 331, 1012, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2234, 3, 331, 334, 1018, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2244, 0, 3, 964, 2164, 1042, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2274, 0, 3, 970, 2174, 1060, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2304, 0, 3, 976, 2184, 1078, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2334, 0, 3, 982, 2194, 1096, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2364, 0, 3, 994, 2204, 1132, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2394, 0, 3, 1000, 2214, 1150, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2424, 0, 3, 1006, 2224, 1168, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2454, 0, 3, 1012, 2234, 1186, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2484, 0, 3, 1024, 2244, 448, 466, 1240,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2544, 0, 3, 1042, 2274, 466, 484, 1276,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2604, 0, 3, 1060, 2304, 484, 502, 1312,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2664, 0, 3, 1078, 2334, 502, 520, 1348,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2724, 0, 3, 1114, 2364, 556, 574, 1420,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2784, 0, 3, 1132, 2394, 574, 592, 1456,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2844, 0, 3, 1150, 2424, 592, 610, 1492,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2904, 0, 3, 1168, 2454, 610, 628, 1528,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2964, 0, 3, 1240, 2544, 664, 694, 1684,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3064, 0, 3, 1276, 2604, 694, 724, 1744,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3164, 0, 3, 1312, 2664, 724, 754, 1804,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3264, 0, 3, 1420, 2784, 814, 844, 1984,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3364, 0, 3, 1456, 2844, 844, 874, 2044,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3464, 0, 3, 1492, 2904, 874, 904, 2104,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3564, 3, 964, 970, 2174, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3579, 3, 970, 976, 2184, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3594, 3, 976, 982, 2194, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3609, 3, 994, 1000, 2214, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3624, 3, 1000, 1006, 2224, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3639, 3, 1006, 1012, 2234, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 3654, 0, 3, 2164, 3564, 2274, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3699, 0, 3, 2174, 3579, 2304, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3744, 0, 3, 2184, 3594, 2334, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3789, 0, 3, 2204, 3609, 2394, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3834, 0, 3, 2214, 3624, 2424, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3879, 0, 3, 2224, 3639, 2454, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3924, 0, 3, 2244, 3654, 1204, 1240,
                                                 2544, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4014, 0, 3, 2274, 3699, 1240, 1276,
                                                 2604, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4104, 0, 3, 2304, 3744, 1276, 1312,
                                                 2664, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4194, 0, 3, 2364, 3789, 1384, 1420,
                                                 2784, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4284, 0, 3, 2394, 3834, 1420, 1456,
                                                 2844, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4374, 0, 3, 2424, 3879, 1456, 1492,
                                                 2904, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4464, 0, 3, 2484, 3924, 1564, 1624,
                                                 2964, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4614, 0, 3, 2544, 4014, 1624, 1684,
                                                 3064, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4764, 0, 3, 2604, 4104, 1684, 1744,
                                                 3164, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4914, 0, 3, 2724, 4194, 1864, 1924,
                                                 3264, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5064, 0, 3, 2784, 4284, 1924, 1984,
                                                 3364, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5214, 0, 3, 2844, 4374, 1984, 2044,
                                                 3464, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5364, 3, 2164, 2174, 3579, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5385, 3, 2174, 2184, 3594, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5406, 3, 2204, 2214, 3624, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5427, 3, 2214, 2224, 3639, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 5448, 0, 3, 3564, 5364, 3699, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5511, 0, 3, 3579, 5385, 3744, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5574, 0, 3, 3609, 5406, 3834, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5637, 0, 3, 3624, 5427, 3879, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 5700, 0, 3, 3654, 5448, 2484, 2544,
                                                 4014, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5826, 0, 3, 3699, 5511, 2544, 2604,
                                                 4104, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5952, 0, 3, 3789, 5574, 2724, 2784,
                                                 4284, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6078, 0, 3, 3834, 5637, 2784, 2844,
                                                 4374, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6204, 0, 3, 4014, 5826, 2964, 3064,
                                                 4764, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6414, 0, 3, 4284, 6078, 3264, 3364,
                                                 5214, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6624, 3, 3564, 3579, 5385, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6652, 3, 3609, 3624, 5427, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 6680, 0, 3, 5364, 6624, 5511, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6764, 0, 3, 5406, 6652, 5637, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 6848, 0, 3, 5448, 6680, 3924, 4014,
                                                 5826, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7016, 0, 3, 5574, 6764, 4194, 4284,
                                                 6078, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 7184, 0, 3, 5700, 6848, 4464, 4614,
                                                 6204, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 7464, 0, 3, 5952, 7016, 4914, 5064,
                                                 6414, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 7744, 7184, 560, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 8304, 8024, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 8304, 13, nmax);

    simdtrf::transform_i_inner(buffer, 8304, 7744, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 91 * nvalues, nvalues, buffer, 8304, 13, nmax);
}

}  // namespace simdt2ceri

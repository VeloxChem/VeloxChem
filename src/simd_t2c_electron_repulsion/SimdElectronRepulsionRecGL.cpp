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


#include "SimdElectronRepulsionRecGL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

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
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gl_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(20094, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 12, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 15, 16, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 16, 17, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 17, 18, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 20, 23, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 23, 26, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 26, 29, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 29, 32, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 32, 35, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 35, 38, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 38, 41, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 41, 44, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 44, 47, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 209, 0, 53, 59, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 224, 0, 59, 65, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 239, 0, 65, 71, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 254, 0, 71, 77, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 269, 0, 77, 83, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 284, 0, 83, 89, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 299, 0, 89, 95, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 314, 0, 95, 101, 189, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 329, 0, 101, 107, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 344, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 347, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 350, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 353, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 356, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 359, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 362, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 365, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 368, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 371, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 374, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 383, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 392, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 401, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 410, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 419, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 428, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 437, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 446, 3, 17, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 455, 3, 18, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 464, 0, 3, 23, 383, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 482, 0, 3, 26, 392, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 500, 0, 3, 29, 401, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 518, 0, 3, 32, 410, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 536, 0, 3, 35, 419, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 554, 0, 3, 38, 428, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 572, 0, 3, 41, 437, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 590, 0, 3, 44, 446, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 608, 0, 3, 47, 455, 113, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 626, 0, 3, 65, 482, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 656, 0, 3, 71, 500, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 686, 0, 3, 77, 518, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 716, 0, 3, 83, 536, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 746, 0, 3, 89, 554, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 776, 0, 3, 95, 572, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 806, 0, 3, 101, 590, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 836, 0, 3, 107, 608, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 866, 0, 3, 129, 656, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 911, 0, 3, 139, 686, 254, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 956, 0, 3, 149, 716, 269, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1001, 0, 3, 159, 746, 284, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1046, 0, 3, 169, 776, 299, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1091, 0, 3, 179, 806, 314, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1136, 0, 3, 189, 836, 329, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1181, 3, 9, 10, 347, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1187, 3, 10, 11, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1193, 3, 11, 12, 353, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1199, 3, 12, 13, 356, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1205, 3, 13, 14, 359, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1211, 3, 14, 15, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1217, 3, 15, 16, 365, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1223, 3, 16, 17, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1229, 3, 17, 18, 371, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1235, 0, 3, 344, 1181, 383, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1253, 0, 3, 347, 1187, 392, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1271, 0, 3, 350, 1193, 401, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1289, 0, 3, 353, 1199, 410, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1307, 0, 3, 356, 1205, 419, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1325, 0, 3, 359, 1211, 428, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1343, 0, 3, 362, 1217, 437, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1361, 0, 3, 365, 1223, 446, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1379, 0, 3, 368, 1229, 455, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1397, 0, 3, 374, 1235, 53, 59, 464,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1433, 0, 3, 383, 1253, 59, 65, 482,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1469, 0, 3, 392, 1271, 65, 71, 500,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1505, 0, 3, 401, 1289, 71, 77, 518,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1541, 0, 3, 410, 1307, 77, 83, 536,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1577, 0, 3, 419, 1325, 83, 89, 554,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1613, 0, 3, 428, 1343, 89, 95, 572,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1649, 0, 3, 437, 1361, 95, 101, 590,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1685, 0, 3, 446, 1379, 101, 107, 608,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1721, 0, 3, 482, 1469, 119, 129, 656,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1781, 0, 3, 500, 1505, 129, 139, 686,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1841, 0, 3, 518, 1541, 139, 149, 716,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1901, 0, 3, 536, 1577, 149, 159, 746,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1961, 0, 3, 554, 1613, 159, 169, 776,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2021, 0, 3, 572, 1649, 169, 179, 806,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2081, 0, 3, 590, 1685, 179, 189, 836,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2141, 0, 3, 1397, 1433, 626, 1721, 209,
                                                 224, 866, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2231, 0, 3, 1433, 1469, 656, 1781, 224,
                                                 239, 911, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2321, 0, 3, 1469, 1505, 686, 1841, 239,
                                                 254, 956, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2411, 0, 3, 1505, 1541, 716, 1901, 254,
                                                 269, 1001, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2501, 0, 3, 1541, 1577, 746, 1961, 269,
                                                 284, 1046, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2591, 0, 3, 1577, 1613, 776, 2021, 284,
                                                 299, 1091, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2681, 0, 3, 1613, 1649, 806, 2081, 299,
                                                 314, 1136, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2771, 3, 344, 347, 1187, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2781, 3, 347, 350, 1193, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2791, 3, 350, 353, 1199, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2801, 3, 353, 356, 1205, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2811, 3, 356, 359, 1211, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2821, 3, 359, 362, 1217, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2831, 3, 362, 365, 1223, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2841, 3, 365, 368, 1229, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2851, 0, 3, 1181, 2771, 1253, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2881, 0, 3, 1187, 2781, 1271, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2911, 0, 3, 1193, 2791, 1289, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2941, 0, 3, 1199, 2801, 1307, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2971, 0, 3, 1205, 2811, 1325, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3001, 0, 3, 1211, 2821, 1343, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3031, 0, 3, 1217, 2831, 1361, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3061, 0, 3, 1223, 2841, 1379, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3091, 0, 3, 1253, 2881, 464, 482, 1469,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3151, 0, 3, 1271, 2911, 482, 500, 1505,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3211, 0, 3, 1289, 2941, 500, 518, 1541,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3271, 0, 3, 1307, 2971, 518, 536, 1577,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3331, 0, 3, 1325, 3001, 536, 554, 1613,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3391, 0, 3, 1343, 3031, 554, 572, 1649,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3451, 0, 3, 1361, 3061, 572, 590, 1685,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3511, 0, 3, 1469, 3151, 626, 656, 1781,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3611, 0, 3, 1505, 3211, 656, 686, 1841,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3711, 0, 3, 1541, 3271, 686, 716, 1901,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3811, 0, 3, 1577, 3331, 716, 746, 1961,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3911, 0, 3, 1613, 3391, 746, 776, 2021,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4011, 0, 3, 1649, 3451, 776, 806, 2081,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4111, 0, 3, 3091, 3151, 1781, 3611, 866,
                                                 911, 2321, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4261, 0, 3, 3151, 3211, 1841, 3711, 911,
                                                 956, 2411, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4411, 0, 3, 3211, 3271, 1901, 3811, 956,
                                                 1001, 2501, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4561, 0, 3, 3271, 3331, 1961, 3911,
                                                 1001, 1046, 2591, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4711, 0, 3, 3331, 3391, 2021, 4011,
                                                 1046, 1091, 2681, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4861, 3, 1181, 1187, 2781, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4876, 3, 1187, 1193, 2791, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4891, 3, 1193, 1199, 2801, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4906, 3, 1199, 1205, 2811, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4921, 3, 1205, 1211, 2821, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4936, 3, 1211, 1217, 2831, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4951, 3, 1217, 1223, 2841, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4966, 0, 3, 2771, 4861, 2881, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5011, 0, 3, 2781, 4876, 2911, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5056, 0, 3, 2791, 4891, 2941, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5101, 0, 3, 2801, 4906, 2971, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5146, 0, 3, 2811, 4921, 3001, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5191, 0, 3, 2821, 4936, 3031, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5236, 0, 3, 2831, 4951, 3061, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 5281, 0, 3, 2851, 4966, 1397, 1433,
                                                 3091, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5371, 0, 3, 2881, 5011, 1433, 1469,
                                                 3151, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5461, 0, 3, 2911, 5056, 1469, 1505,
                                                 3211, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5551, 0, 3, 2941, 5101, 1505, 1541,
                                                 3271, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5641, 0, 3, 2971, 5146, 1541, 1577,
                                                 3331, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5731, 0, 3, 3001, 5191, 1577, 1613,
                                                 3391, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5821, 0, 3, 3031, 5236, 1613, 1649,
                                                 3451, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5911, 0, 3, 3151, 5461, 1721, 1781,
                                                 3611, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6061, 0, 3, 3211, 5551, 1781, 1841,
                                                 3711, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6211, 0, 3, 3271, 5641, 1841, 1901,
                                                 3811, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6361, 0, 3, 3331, 5731, 1901, 1961,
                                                 3911, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6511, 0, 3, 3391, 5821, 1961, 2021,
                                                 4011, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 6661, 0, 3, 5281, 5371, 3511, 5911,
                                                 2141, 2231, 4111, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 6886, 0, 3, 5371, 5461, 3611, 6061,
                                                 2231, 2321, 4261, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7111, 0, 3, 5461, 5551, 3711, 6211,
                                                 2321, 2411, 4411, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7336, 0, 3, 5551, 5641, 3811, 6361,
                                                 2411, 2501, 4561, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7561, 0, 3, 5641, 5731, 3911, 6511,
                                                 2501, 2591, 4711, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7786, 3, 2771, 2781, 4876, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7807, 3, 2781, 2791, 4891, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7828, 3, 2791, 2801, 4906, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7849, 3, 2801, 2811, 4921, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7870, 3, 2811, 2821, 4936, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 7891, 3, 2821, 2831, 4951, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 7912, 0, 3, 4861, 7786, 5011, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 7975, 0, 3, 4876, 7807, 5056, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8038, 0, 3, 4891, 7828, 5101, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8101, 0, 3, 4906, 7849, 5146, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8164, 0, 3, 4921, 7870, 5191, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8227, 0, 3, 4936, 7891, 5236, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 8290, 0, 3, 5011, 7975, 3091, 3151,
                                                 5461, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8416, 0, 3, 5056, 8038, 3151, 3211,
                                                 5551, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8542, 0, 3, 5101, 8101, 3211, 3271,
                                                 5641, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8668, 0, 3, 5146, 8164, 3271, 3331,
                                                 5731, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 8794, 0, 3, 5191, 8227, 3331, 3391,
                                                 5821, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 8920, 0, 3, 5461, 8416, 3511, 3611,
                                                 6061, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9130, 0, 3, 5551, 8542, 3611, 3711,
                                                 6211, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9340, 0, 3, 5641, 8668, 3711, 3811,
                                                 6361, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9550, 0, 3, 5731, 8794, 3811, 3911,
                                                 6511, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 9760, 0, 3, 8290, 8416, 6061, 9130,
                                                 4111, 4261, 7111, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 10075, 0, 3, 8416, 8542, 6211, 9340,
                                                 4261, 4411, 7336, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 10390, 0, 3, 8542, 8668, 6361, 9550,
                                                 4411, 4561, 7561, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 10705, 3, 4861, 4876, 7807, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 10733, 3, 4876, 4891, 7828, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 10761, 3, 4891, 4906, 7849, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 10789, 3, 4906, 4921, 7870, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 10817, 3, 4921, 4936, 7891, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 10845, 0, 3, 7786, 10705, 7975, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 10929, 0, 3, 7807, 10733, 8038, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11013, 0, 3, 7828, 10761, 8101, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11097, 0, 3, 7849, 10789, 8164, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11181, 0, 3, 7870, 10817, 8227, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 11265, 0, 3, 7912, 10845, 5281, 5371,
                                                 8290, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11433, 0, 3, 7975, 10929, 5371, 5461,
                                                 8416, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11601, 0, 3, 8038, 11013, 5461, 5551,
                                                 8542, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11769, 0, 3, 8101, 11097, 5551, 5641,
                                                 8668, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11937, 0, 3, 8164, 11181, 5641, 5731,
                                                 8794, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12105, 0, 3, 8416, 11601, 5911, 6061,
                                                 9130, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12385, 0, 3, 8542, 11769, 6061, 6211,
                                                 9340, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12665, 0, 3, 8668, 11937, 6211, 6361,
                                                 9550, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 12945, 0, 3, 11265, 11433, 8920, 12105,
                                                 6661, 6886, 9760, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 13365, 0, 3, 11433, 11601, 9130, 12385,
                                                 6886, 7111, 10075, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 13785, 0, 3, 11601, 11769, 9340, 12665,
                                                 7111, 7336, 10390, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14205, 3, 7786, 7807, 10733, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14241, 3, 7807, 7828, 10761, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14277, 3, 7828, 7849, 10789, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 14313, 3, 7849, 7870, 10817, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 14349, 0, 3, 10705, 14205, 10929, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14457, 0, 3, 10733, 14241, 11013, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14565, 0, 3, 10761, 14277, 11097, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 14673, 0, 3, 10789, 14313, 11181, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 14781, 0, 3, 10929, 14457, 8290, 8416,
                                                 11601, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 14997, 0, 3, 11013, 14565, 8416, 8542,
                                                 11769, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 15213, 0, 3, 11097, 14673, 8542, 8668,
                                                 11937, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 15429, 0, 3, 11601, 14997, 8920, 9130,
                                                 12385, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 15789, 0, 3, 11769, 15213, 9130, 9340,
                                                 12665, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 16149, 0, 3, 14781, 14997, 12385, 15789,
                                                 9760, 10075, 13785, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 16689, 3, 10705, 10733, 14241, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 16734, 3, 10733, 10761, 14277, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 16779, 3, 10761, 10789, 14313, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 16824, 0, 3, 14205, 16689, 14457, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 16959, 0, 3, 14241, 16734, 14565, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 17094, 0, 3, 14277, 16779, 14673, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 17229, 0, 3, 14349, 16824, 11265, 11433,
                                                 14781, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 17499, 0, 3, 14457, 16959, 11433, 11601,
                                                 14997, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 17769, 0, 3, 14565, 17094, 11601, 11769,
                                                 15213, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 18039, 0, 3, 14997, 17769, 12105, 12385,
                                                 15789, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 18489, 0, 3, 17229, 17499, 15429, 18039,
                                                 12945, 13365, 16149, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 19164, 18489, 675, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 19839, 19164, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 19839, 17, nmax);
}

}  // namespace simdt2ceri

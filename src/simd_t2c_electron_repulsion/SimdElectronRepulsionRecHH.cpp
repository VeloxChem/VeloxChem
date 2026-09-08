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


#include "SimdElectronRepulsionRecHH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformHH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hh_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2947, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 47, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 50, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 53, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 12, 13, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 13, 14, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 14, 15, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 68, 0, 17, 20, 47, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 74, 0, 20, 23, 50, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 80, 0, 23, 26, 53, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 86, 0, 26, 29, 56, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 92, 0, 29, 32, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 98, 0, 32, 35, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 104, 0, 35, 38, 65, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 110, 0, 44, 47, 74, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 116, 0, 47, 50, 80, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 125, 0, 50, 53, 86, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 134, 0, 53, 56, 92, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 143, 0, 56, 59, 98, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 152, 0, 59, 62, 104, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 161, 0, 68, 74, 116, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 170, 0, 74, 80, 125, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 179, 0, 80, 86, 134, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 188, 0, 86, 92, 143, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 197, 0, 92, 98, 152, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 206, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 209, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 212, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 215, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 218, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 221, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 224, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 227, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 230, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 233, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 236, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 239, 3, 14, 38, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 242, 3, 20, 47, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 245, 3, 23, 50, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 248, 3, 26, 53, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 257, 3, 29, 56, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 266, 3, 32, 59, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 275, 3, 35, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 284, 3, 38, 65, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 287, 3, 44, 68, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 290, 3, 47, 74, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 293, 3, 50, 80, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 302, 3, 53, 86, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 311, 3, 56, 92, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 320, 3, 59, 98, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 329, 3, 62, 104, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 338, 3, 74, 116, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 350, 3, 80, 125, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 362, 3, 86, 134, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 374, 3, 92, 143, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 386, 3, 98, 152, ncols, p);

            compute_prim_hp_electron_repulsion_2(buffer, 398, 3, 110, 161, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 413, 3, 116, 170, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 428, 3, 125, 179, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 443, 3, 134, 188, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 458, 3, 143, 197, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 473, 3, 9, 10, 209, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 476, 3, 10, 11, 212, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 479, 3, 11, 12, 215, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 482, 3, 12, 13, 218, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 485, 3, 13, 14, 221, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 488, 0, 206, 473, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 491, 0, 209, 476, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 494, 0, 212, 479, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 497, 0, 215, 482, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 500, 0, 218, 485, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 503, 3, 224, 44, 47, 245, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 506, 3, 227, 47, 50, 248, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 509, 3, 230, 50, 53, 257, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 518, 3, 233, 53, 56, 266, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 527, 3, 236, 56, 59, 275, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 536, 3, 239, 59, 62, 284, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 539, 0, 3, 245, 506, 68, 74, 293, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 548, 0, 3, 248, 509, 74, 80, 302, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_12(buffer, 563, 0, 3, 257, 518, 80, 86, 311, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 581, 0, 3, 266, 527, 86, 92, 320, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 596, 0, 3, 275, 536, 92, 98, 329, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_14(buffer, 611, 0, 3, 503, 506, 293, 548, 110, 116, 350, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_15(buffer, 632, 0, 3, 506, 509, 302, 563, 116, 125, 362, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_16(buffer, 656, 0, 3, 509, 518, 311, 581, 125, 134, 374, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 683, 0, 3, 518, 527, 320, 596, 134, 143, 386, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_5(buffer, 707, 0, 3, 539, 548, 350, 632, 161, 170, 428, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_6(buffer, 734, 0, 3, 548, 563, 362, 656, 170, 179, 443, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_7(buffer, 761, 0, 3, 563, 581, 374, 683, 179, 188, 458, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 788, 3, 206, 209, 476, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 791, 3, 209, 212, 479, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 794, 3, 212, 215, 482, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 797, 3, 215, 218, 485, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 800, 0, 473, 788, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 803, 0, 476, 791, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 806, 0, 479, 794, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 809, 0, 482, 797, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 812, 3, 488, 242, 245, 506, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 815, 3, 491, 245, 248, 509, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_22(buffer, 818, 0, 3, 494, 806, 248, 257, 518, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 839, 3, 497, 257, 266, 527, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 857, 3, 500, 266, 275, 536, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 860, 0, 3, 503, 812, 287, 290, 539, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_20(buffer, 869, 0, 3, 506, 815, 290, 293, 548, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_21(buffer, 878, 0, 3, 509, 818, 293, 302, 563, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_22(buffer, 920, 0, 3, 518, 839, 302, 311, 581, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_19(buffer, 956, 0, 3, 527, 857, 311, 320, 596, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_11(buffer, 980, 0, 3, 812, 815, 548, 878, 338, 350, 632, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_12(buffer, 1022, 0, 3, 815, 818, 563, 920, 350, 362, 656, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_13(buffer, 1085, 0, 3, 818, 839, 581, 956, 362, 374, 683, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_2(buffer, 1133, 0, 3, 860, 869, 611, 980, 398, 413, 707, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_3(buffer, 1187, 0, 3, 869, 878, 632, 1022, 413, 428, 734, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_4(buffer, 1241, 0, 3, 878, 920, 656, 1085, 428, 443, 761, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1319, 3, 473, 476, 791, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1322, 3, 476, 479, 794, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 1325, 3, 479, 482, 797, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1328, 0, 788, 1319, ncols, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1331, 0, 791, 1322, ncols, p);

            compute_prim_pg_electron_repulsion_14(buffer, 1334, 0, 794, 1325, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1337, 3, 800, 503, 506, 815, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_19(buffer, 1340, 3, 803, 506, 509, 818, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_20(buffer, 1343, 3, 806, 509, 518, 839, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_21(buffer, 1361, 3, 809, 518, 527, 857, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_13(buffer, 1364, 0, 3, 815, 1340, 539, 548, 878, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_14(buffer, 1373, 0, 3, 818, 1343, 548, 563, 920, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 1432, 0, 3, 839, 1361, 563, 581, 956, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_6(buffer, 1465, 0, 3, 1337, 1340, 878, 1373, 611, 632, 1022, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_7(buffer, 1607, 0, 3, 1340, 1343, 920, 1432, 632, 656, 1085, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_1(buffer, 1693, 0, 3, 1364, 1373, 1022, 1607, 707, 734, 1241, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_5(buffer, 1897, 3, 1328, 812, 815, 1340, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_11(buffer, 1900, 3, 1331, 815, 818, 1343, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_12(buffer, 1903, 3, 1334, 818, 839, 1361, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_5(buffer, 1906, 0, 3, 1337, 1897, 860, 869, 1364, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_6(buffer, 1915, 0, 3, 1340, 1900, 869, 878, 1373, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_7(buffer, 1924, 0, 3, 1343, 1903, 878, 920, 1432, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_2(buffer, 1957, 0, 3, 1897, 1900, 1373, 1924, 980, 1022, 1607, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 2065, 0, 3, 1906, 1915, 1465, 1957, 1133, 1187, 1693, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2506, 2065, 441, ncols);
        }
    }

    simdtrf::transform_hh_tri(values, nvalues, buffer, 2506, nmax);
}

}  // namespace simdt2ceri

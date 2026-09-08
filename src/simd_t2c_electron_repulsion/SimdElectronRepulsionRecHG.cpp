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


#include "SimdElectronRepulsionRecHG.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5436, 4932, 315, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 85, 0, 16, 19, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 95, 0, 19, 22, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 22, 25, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 25, 28, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 28, 31, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 31, 34, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 34, 37, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 155, 0, 43, 49, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 170, 0, 49, 55, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 185, 0, 55, 61, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 200, 0, 61, 67, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 215, 0, 67, 73, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 230, 0, 85, 95, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 251, 0, 95, 105, 170, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 272, 0, 105, 115, 185, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 293, 0, 115, 125, 200, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 314, 0, 125, 135, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 335, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 338, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 341, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 344, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 347, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 350, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 353, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 362, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 371, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 380, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 389, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 398, 3, 14, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 407, 0, 3, 22, 353, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 425, 0, 3, 25, 362, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 443, 0, 3, 28, 371, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 461, 0, 3, 31, 380, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 479, 0, 3, 34, 389, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 497, 0, 3, 37, 398, 79, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 515, 0, 3, 49, 425, 105, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 545, 0, 3, 55, 443, 115, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 575, 0, 3, 61, 461, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 605, 0, 3, 67, 479, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 635, 0, 3, 73, 497, 145, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 665, 0, 3, 105, 545, 170, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 710, 0, 3, 115, 575, 185, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 755, 0, 3, 125, 605, 200, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 800, 0, 3, 135, 635, 215, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 845, 0, 3, 170, 710, 272, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 908, 0, 3, 185, 755, 293, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 971, 0, 3, 200, 800, 314, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1034, 3, 9, 10, 338, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1040, 3, 10, 11, 341, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1046, 3, 11, 12, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1052, 3, 12, 13, 347, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1058, 3, 13, 14, 350, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1064, 0, 3, 335, 1034, 362, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1082, 0, 3, 338, 1040, 371, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1100, 0, 3, 341, 1046, 380, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1118, 0, 3, 344, 1052, 389, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1136, 0, 3, 347, 1058, 398, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1154, 0, 3, 353, 1064, 43, 49, 425,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1190, 0, 3, 362, 1082, 49, 55, 443,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1226, 0, 3, 371, 1100, 55, 61, 461,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1262, 0, 3, 380, 1118, 61, 67, 479,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1298, 0, 3, 389, 1136, 67, 73, 497,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1334, 0, 3, 407, 1154, 85, 95, 515,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1394, 0, 3, 425, 1190, 95, 105, 545,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1454, 0, 3, 443, 1226, 105, 115, 575,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1514, 0, 3, 461, 1262, 115, 125, 605,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1574, 0, 3, 479, 1298, 125, 135, 635,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1634, 0, 3, 1154, 1190, 545, 1454, 155,
                                                 170, 710, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1724, 0, 3, 1190, 1226, 575, 1514, 170,
                                                 185, 755, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1814, 0, 3, 1226, 1262, 605, 1574, 185,
                                                 200, 800, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 1904, 0, 3, 1334, 1394, 665, 1634, 230,
                                                 251, 845, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2030, 0, 3, 1394, 1454, 710, 1724, 251,
                                                 272, 908, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2156, 0, 3, 1454, 1514, 755, 1814, 272,
                                                 293, 971, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2282, 3, 335, 338, 1040, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2292, 3, 338, 341, 1046, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2302, 3, 341, 344, 1052, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2312, 3, 344, 347, 1058, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2322, 0, 3, 1034, 2282, 1082, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2352, 0, 3, 1040, 2292, 1100, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2382, 0, 3, 1046, 2302, 1118, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2412, 0, 3, 1052, 2312, 1136, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2442, 0, 3, 1064, 2322, 407, 425, 1190,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2502, 0, 3, 1082, 2352, 425, 443, 1226,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2562, 0, 3, 1100, 2382, 443, 461, 1262,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2622, 0, 3, 1118, 2412, 461, 479, 1298,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2682, 0, 3, 1190, 2502, 515, 545, 1454,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2782, 0, 3, 1226, 2562, 545, 575, 1514,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2882, 0, 3, 1262, 2622, 575, 605, 1574,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2982, 0, 3, 2442, 2502, 1454, 2782, 665,
                                                 710, 1724, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3132, 0, 3, 2502, 2562, 1514, 2882, 710,
                                                 755, 1814, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 3282, 0, 3, 2682, 2782, 1724, 3132, 845,
                                                 908, 2156, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3492, 3, 1034, 1040, 2292, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3507, 3, 1040, 1046, 2302, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3522, 3, 1046, 1052, 2312, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 3537, 0, 3, 2282, 3492, 2352, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3582, 0, 3, 2292, 3507, 2382, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3627, 0, 3, 2302, 3522, 2412, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3672, 0, 3, 2322, 3537, 1154, 1190,
                                                 2502, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3762, 0, 3, 2352, 3582, 1190, 1226,
                                                 2562, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3852, 0, 3, 2382, 3627, 1226, 1262,
                                                 2622, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3942, 0, 3, 2442, 3672, 1334, 1394,
                                                 2682, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4092, 0, 3, 2502, 3762, 1394, 1454,
                                                 2782, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4242, 0, 3, 2562, 3852, 1454, 1514,
                                                 2882, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 4392, 0, 3, 3672, 3762, 2782, 4242,
                                                 1634, 1724, 3132, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 4617, 0, 3, 3942, 4092, 2982, 4392,
                                                 1904, 2030, 3282, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4932, 4617, 315, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 5247, 4932, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 5247, 9, nmax);
}

}  // namespace simdt2ceri

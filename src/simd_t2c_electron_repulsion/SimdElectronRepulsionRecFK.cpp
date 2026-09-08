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


#include "SimdElectronRepulsionRecFK.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(6622, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 12, 13, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 13, 14, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 14, 15, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 17, 20, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 20, 23, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 23, 26, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 26, 29, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 29, 32, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 32, 35, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 35, 38, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 162, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 165, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 168, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 171, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 174, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 177, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 180, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 183, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 192, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 201, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 210, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 219, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 228, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 237, 3, 15, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 246, 0, 3, 20, 183, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 264, 0, 3, 23, 192, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 282, 0, 3, 26, 201, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 300, 0, 3, 29, 210, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 318, 0, 3, 32, 219, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 336, 0, 3, 35, 228, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 354, 0, 3, 38, 237, 86, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 372, 0, 3, 44, 246, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 402, 0, 3, 50, 264, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 432, 0, 3, 56, 282, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 462, 0, 3, 62, 300, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 492, 0, 3, 68, 318, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 522, 0, 3, 74, 336, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 552, 0, 3, 80, 354, 152, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 582, 3, 9, 10, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 588, 3, 10, 11, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 594, 3, 11, 12, 171, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 600, 3, 12, 13, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 606, 3, 13, 14, 177, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 612, 3, 14, 15, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 618, 0, 3, 162, 582, 192, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 636, 0, 3, 165, 588, 201, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 654, 0, 3, 168, 594, 210, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 672, 0, 3, 171, 600, 219, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 690, 0, 3, 174, 606, 228, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 708, 0, 3, 177, 612, 237, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 726, 0, 3, 183, 618, 44, 50, 264, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 762, 0, 3, 192, 636, 50, 56, 282, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 798, 0, 3, 201, 654, 56, 62, 300, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 834, 0, 3, 210, 672, 62, 68, 318, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 870, 0, 3, 219, 690, 68, 74, 336, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 906, 0, 3, 228, 708, 74, 80, 354, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 942, 0, 3, 264, 762, 92, 102, 432,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1002, 0, 3, 282, 798, 102, 112, 462,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1062, 0, 3, 300, 834, 112, 122, 492,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1122, 0, 3, 318, 870, 122, 132, 522,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1182, 0, 3, 336, 906, 132, 142, 552,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1242, 3, 162, 165, 588, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1252, 3, 165, 168, 594, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1262, 3, 168, 171, 600, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1272, 3, 171, 174, 606, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1282, 3, 174, 177, 612, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1292, 0, 3, 582, 1242, 636, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1322, 0, 3, 588, 1252, 654, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1352, 0, 3, 594, 1262, 672, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1382, 0, 3, 600, 1272, 690, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1412, 0, 3, 606, 1282, 708, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1442, 0, 3, 618, 1292, 246, 264, 762,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1502, 0, 3, 636, 1322, 264, 282, 798,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1562, 0, 3, 654, 1352, 282, 300, 834,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1622, 0, 3, 672, 1382, 300, 318, 870,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1682, 0, 3, 690, 1412, 318, 336, 906,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1742, 0, 3, 726, 1442, 372, 402, 942,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1842, 0, 3, 762, 1502, 402, 432, 1002,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1942, 0, 3, 798, 1562, 432, 462, 1062,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2042, 0, 3, 834, 1622, 462, 492, 1122,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2142, 0, 3, 870, 1682, 492, 522, 1182,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2242, 3, 582, 588, 1252, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2257, 3, 588, 594, 1262, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2272, 3, 594, 600, 1272, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2287, 3, 600, 606, 1282, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2302, 0, 3, 1242, 2242, 1322, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2347, 0, 3, 1252, 2257, 1352, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2392, 0, 3, 1262, 2272, 1382, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2437, 0, 3, 1272, 2287, 1412, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2482, 0, 3, 1292, 2302, 726, 762, 1502,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2572, 0, 3, 1322, 2347, 762, 798, 1562,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2662, 0, 3, 1352, 2392, 798, 834, 1622,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2752, 0, 3, 1382, 2437, 834, 870, 1682,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2842, 0, 3, 1502, 2572, 942, 1002, 1942,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2992, 0, 3, 1562, 2662, 1002, 1062,
                                                 2042, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3142, 0, 3, 1622, 2752, 1062, 1122,
                                                 2142, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3292, 3, 1242, 1252, 2257, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3313, 3, 1252, 1262, 2272, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3334, 3, 1262, 1272, 2287, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3355, 0, 3, 2242, 3292, 2347, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3418, 0, 3, 2257, 3313, 2392, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3481, 0, 3, 2272, 3334, 2437, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3544, 0, 3, 2302, 3355, 1442, 1502,
                                                 2572, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3670, 0, 3, 2347, 3418, 1502, 1562,
                                                 2662, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3796, 0, 3, 2392, 3481, 1562, 1622,
                                                 2752, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 3922, 0, 3, 2482, 3544, 1742, 1842,
                                                 2842, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4132, 0, 3, 2572, 3670, 1842, 1942,
                                                 2992, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4342, 0, 3, 2662, 3796, 1942, 2042,
                                                 3142, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4552, 3, 2242, 2257, 3313, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4580, 3, 2257, 2272, 3334, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 4608, 0, 3, 3292, 4552, 3418, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4692, 0, 3, 3313, 4580, 3481, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 4776, 0, 3, 3355, 4608, 2482, 2572,
                                                 3670, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 4944, 0, 3, 3418, 4692, 2572, 2662,
                                                 3796, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 5112, 0, 3, 3670, 4944, 2842, 2992,
                                                 4342, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5392, 3, 3292, 3313, 4580, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 5428, 0, 3, 4552, 5392, 4692, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 5536, 0, 3, 4608, 5428, 3544, 3670,
                                                 4944, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 5752, 0, 3, 4776, 5536, 3922, 4132,
                                                 5112, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 6112, 5752, 360, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 6472, 6112, 10, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 6472, 15, nmax);
}

}  // namespace simdt2ceri

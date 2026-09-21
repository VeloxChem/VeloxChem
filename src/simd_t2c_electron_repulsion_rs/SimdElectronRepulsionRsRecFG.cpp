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


#include "SimdElectronRepulsionRsRecFG.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_fg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2954, 2564, 300, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 7, 8, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 8, 9, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 9, 10, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 10, 11, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 11, 12, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 15, 16, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 16, 17, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 17, 18, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 124, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 134, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 144, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 154, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 34, 37, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 43, 46, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 46, 49, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 49, 52, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 52, 55, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 55, 58, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 224, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 227, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 230, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 233, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 236, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 239, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 242, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 245, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 248, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 257, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 266, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 275, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 284, 3, 17, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 293, 3, 18, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 302, 3, 19, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 311, 3, 20, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 320, 0, 3, 28, 248, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 338, 0, 3, 31, 257, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 356, 0, 3, 34, 266, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 374, 0, 3, 37, 275, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 392, 0, 3, 49, 284, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 410, 0, 3, 52, 293, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 428, 0, 3, 55, 302, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 446, 0, 3, 58, 311, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 464, 0, 3, 70, 338, 144, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 494, 0, 3, 76, 356, 154, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 524, 0, 3, 82, 374, 164, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 554, 0, 3, 100, 410, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 584, 0, 3, 106, 428, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 614, 0, 3, 112, 446, 214, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 644, 3, 9, 10, 227, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 650, 3, 10, 11, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 656, 3, 11, 12, 233, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 662, 3, 17, 18, 239, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 668, 3, 18, 19, 242, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 674, 3, 19, 20, 245, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 680, 0, 3, 224, 644, 257, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 698, 0, 3, 227, 650, 266, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 716, 0, 3, 230, 656, 275, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 734, 0, 3, 236, 662, 293, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 752, 0, 3, 239, 668, 302, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 770, 0, 3, 242, 674, 311, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 788, 0, 3, 248, 680, 64, 70, 338, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 824, 0, 3, 257, 698, 70, 76, 356, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 860, 0, 3, 266, 716, 76, 82, 374, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 896, 0, 3, 284, 734, 94, 100, 410,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 932, 0, 3, 293, 752, 100, 106, 428,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 968, 0, 3, 302, 770, 106, 112, 446,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1004, 0, 3, 320, 788, 124, 134, 464,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1064, 0, 3, 338, 824, 134, 144, 494,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1124, 0, 3, 356, 860, 144, 154, 524,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1184, 0, 3, 392, 896, 174, 184, 554,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1244, 0, 3, 410, 932, 184, 194, 584,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1304, 0, 3, 428, 968, 194, 204, 614,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1364, 3, 224, 227, 650, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1374, 3, 227, 230, 656, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1384, 3, 236, 239, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1394, 3, 239, 242, 674, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1404, 0, 3, 644, 1364, 698, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1434, 0, 3, 650, 1374, 716, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1464, 0, 3, 662, 1384, 752, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1494, 0, 3, 668, 1394, 770, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1524, 0, 3, 680, 1404, 320, 338, 824,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1584, 0, 3, 698, 1434, 338, 356, 860,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1644, 0, 3, 734, 1464, 392, 410, 932,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1704, 0, 3, 752, 1494, 410, 428, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1764, 0, 3, 824, 1584, 464, 494, 1124,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1864, 0, 3, 932, 1704, 554, 584, 1304,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1964, 3, 644, 650, 1374, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1979, 3, 662, 668, 1394, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1994, 0, 3, 1364, 1964, 1434, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2039, 0, 3, 1384, 1979, 1494, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2084, 0, 3, 1404, 1994, 788, 824, 1584,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2174, 0, 3, 1464, 2039, 896, 932, 1704,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2264, 0, 3, 1524, 2084, 1004, 1064,
                                                 1764, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2414, 0, 3, 1644, 2174, 1184, 1244,
                                                 1864, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2564, 2264, 300, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 2864, 2714, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 2864, 9, nmax);

    simdtrf::transform_g_inner(buffer, 2864, 2564, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 63 * nvalues, nvalues, buffer, 2864, 9, nmax);
}

}  // namespace simdt2ceri

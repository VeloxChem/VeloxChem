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


#include "SimdElectronRepulsionRecGI.hpp"

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
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gi_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 8940, 8325, 420, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 13, 14, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 14, 15, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 15, 16, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 99, 0, 18, 21, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 109, 0, 21, 24, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 24, 27, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 129, 0, 27, 30, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 30, 33, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 33, 36, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 36, 39, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 169, 0, 45, 51, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 184, 0, 51, 57, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 199, 0, 57, 63, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 214, 0, 63, 69, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 229, 0, 69, 75, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 244, 0, 75, 81, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 259, 0, 81, 87, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 274, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 277, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 280, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 283, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 286, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 289, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 292, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 295, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 298, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 307, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 316, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 325, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 334, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 343, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 352, 3, 15, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 361, 3, 16, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 370, 0, 3, 21, 307, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 388, 0, 3, 24, 316, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 406, 0, 3, 27, 325, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 424, 0, 3, 30, 334, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 442, 0, 3, 33, 343, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 460, 0, 3, 36, 352, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 478, 0, 3, 39, 361, 93, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 496, 0, 3, 57, 388, 109, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 526, 0, 3, 63, 406, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 556, 0, 3, 69, 424, 129, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 586, 0, 3, 75, 442, 139, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 616, 0, 3, 81, 460, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 646, 0, 3, 87, 478, 159, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 676, 0, 3, 109, 526, 199, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 721, 0, 3, 119, 556, 214, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 766, 0, 3, 129, 586, 229, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 811, 0, 3, 139, 616, 244, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 856, 0, 3, 149, 646, 259, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 901, 3, 9, 10, 277, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 907, 3, 10, 11, 280, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 913, 3, 11, 12, 283, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 919, 3, 12, 13, 286, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 925, 3, 13, 14, 289, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 931, 3, 14, 15, 292, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 937, 3, 15, 16, 295, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 943, 0, 3, 274, 901, 307, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 961, 0, 3, 277, 907, 316, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 979, 0, 3, 280, 913, 325, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 997, 0, 3, 283, 919, 334, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1015, 0, 3, 286, 925, 343, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1033, 0, 3, 289, 931, 352, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1051, 0, 3, 292, 937, 361, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1069, 0, 3, 298, 943, 45, 51, 370,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1105, 0, 3, 307, 961, 51, 57, 388,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1141, 0, 3, 316, 979, 57, 63, 406,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1177, 0, 3, 325, 997, 63, 69, 424,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1213, 0, 3, 334, 1015, 69, 75, 442,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1249, 0, 3, 343, 1033, 75, 81, 460,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1285, 0, 3, 352, 1051, 81, 87, 478,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1321, 0, 3, 388, 1141, 99, 109, 526,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1381, 0, 3, 406, 1177, 109, 119, 556,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1441, 0, 3, 424, 1213, 119, 129, 586,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1501, 0, 3, 442, 1249, 129, 139, 616,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1561, 0, 3, 460, 1285, 139, 149, 646,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1621, 0, 3, 1069, 1105, 496, 1321, 169,
                                                 184, 676, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1711, 0, 3, 1105, 1141, 526, 1381, 184,
                                                 199, 721, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1801, 0, 3, 1141, 1177, 556, 1441, 199,
                                                 214, 766, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1891, 0, 3, 1177, 1213, 586, 1501, 214,
                                                 229, 811, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1981, 0, 3, 1213, 1249, 616, 1561, 229,
                                                 244, 856, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2071, 3, 274, 277, 907, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2081, 3, 277, 280, 913, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2091, 3, 280, 283, 919, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2101, 3, 283, 286, 925, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2111, 3, 286, 289, 931, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2121, 3, 289, 292, 937, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2131, 0, 3, 901, 2071, 961, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2161, 0, 3, 907, 2081, 979, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2191, 0, 3, 913, 2091, 997, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2221, 0, 3, 919, 2101, 1015, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2251, 0, 3, 925, 2111, 1033, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2281, 0, 3, 931, 2121, 1051, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 2311, 0, 3, 961, 2161, 370, 388, 1141,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2371, 0, 3, 979, 2191, 388, 406, 1177,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2431, 0, 3, 997, 2221, 406, 424, 1213,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2491, 0, 3, 1015, 2251, 424, 442, 1249,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2551, 0, 3, 1033, 2281, 442, 460, 1285,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2611, 0, 3, 1141, 2371, 496, 526, 1381,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2711, 0, 3, 1177, 2431, 526, 556, 1441,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2811, 0, 3, 1213, 2491, 556, 586, 1501,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2911, 0, 3, 1249, 2551, 586, 616, 1561,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3011, 0, 3, 2311, 2371, 1381, 2711, 676,
                                                 721, 1801, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3161, 0, 3, 2371, 2431, 1441, 2811, 721,
                                                 766, 1891, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3311, 0, 3, 2431, 2491, 1501, 2911, 766,
                                                 811, 1981, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3461, 3, 901, 907, 2081, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3476, 3, 907, 913, 2091, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3491, 3, 913, 919, 2101, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3506, 3, 919, 925, 2111, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 3521, 3, 925, 931, 2121, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 3536, 0, 3, 2071, 3461, 2161, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3581, 0, 3, 2081, 3476, 2191, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3626, 0, 3, 2091, 3491, 2221, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3671, 0, 3, 2101, 3506, 2251, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3716, 0, 3, 2111, 3521, 2281, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3761, 0, 3, 2131, 3536, 1069, 1105,
                                                 2311, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3851, 0, 3, 2161, 3581, 1105, 1141,
                                                 2371, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3941, 0, 3, 2191, 3626, 1141, 1177,
                                                 2431, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4031, 0, 3, 2221, 3671, 1177, 1213,
                                                 2491, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4121, 0, 3, 2251, 3716, 1213, 1249,
                                                 2551, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4211, 0, 3, 2371, 3941, 1321, 1381,
                                                 2711, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4361, 0, 3, 2431, 4031, 1381, 1441,
                                                 2811, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4511, 0, 3, 2491, 4121, 1441, 1501,
                                                 2911, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 4661, 0, 3, 3761, 3851, 2611, 4211,
                                                 1621, 1711, 3011, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 4886, 0, 3, 3851, 3941, 2711, 4361,
                                                 1711, 1801, 3161, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5111, 0, 3, 3941, 4031, 2811, 4511,
                                                 1801, 1891, 3311, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5336, 3, 2071, 2081, 3476, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5357, 3, 2081, 2091, 3491, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5378, 3, 2091, 2101, 3506, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 5399, 3, 2101, 2111, 3521, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 5420, 0, 3, 3461, 5336, 3581, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5483, 0, 3, 3476, 5357, 3626, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5546, 0, 3, 3491, 5378, 3671, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 5609, 0, 3, 3506, 5399, 3716, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 5672, 0, 3, 3581, 5483, 2311, 2371,
                                                 3941, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5798, 0, 3, 3626, 5546, 2371, 2431,
                                                 4031, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5924, 0, 3, 3671, 5609, 2431, 2491,
                                                 4121, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6050, 0, 3, 3941, 5798, 2611, 2711,
                                                 4361, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 6260, 0, 3, 4031, 5924, 2711, 2811,
                                                 4511, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 6470, 0, 3, 5672, 5798, 4361, 6260,
                                                 3011, 3161, 5111, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6785, 3, 3461, 3476, 5357, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6813, 3, 3476, 3491, 5378, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 6841, 3, 3491, 3506, 5399, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 6869, 0, 3, 5336, 6785, 5483, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6953, 0, 3, 5357, 6813, 5546, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 7037, 0, 3, 5378, 6841, 5609, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 7121, 0, 3, 5420, 6869, 3761, 3851,
                                                 5672, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7289, 0, 3, 5483, 6953, 3851, 3941,
                                                 5798, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7457, 0, 3, 5546, 7037, 3941, 4031,
                                                 5924, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 7625, 0, 3, 5798, 7457, 4211, 4361,
                                                 6260, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 7905, 0, 3, 7121, 7289, 6050, 7625,
                                                 4661, 4886, 6470, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8325, 7905, 420, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 8745, 8325, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 8745, 13, nmax);
}

}  // namespace simdt2ceri

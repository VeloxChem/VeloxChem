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


#include "SimdElectronRepulsionGeom10RecDK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecDK.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_dk_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_dk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 7892, 7154, 648, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 162, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 165, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 168, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 171, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 174, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 177, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 180, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 183, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 186, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 189, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 198, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 207, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 216, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 225, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 234, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 243, 3, 15, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 252, 0, 3, 20, 189, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 270, 0, 3, 23, 198, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 288, 0, 3, 26, 207, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 306, 0, 3, 29, 216, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 324, 0, 3, 32, 225, 74, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 342, 0, 3, 35, 234, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 360, 0, 3, 38, 243, 86, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 378, 0, 3, 44, 252, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 408, 0, 3, 50, 270, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 438, 0, 3, 56, 288, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 468, 0, 3, 62, 306, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 498, 0, 3, 68, 324, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 528, 0, 3, 74, 342, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 558, 0, 3, 80, 360, 152, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 588, 3, 7, 8, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 594, 3, 8, 9, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 600, 3, 9, 10, 171, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 606, 3, 10, 11, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 612, 3, 11, 12, 177, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 618, 3, 12, 13, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 624, 3, 13, 14, 183, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 630, 3, 14, 15, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 636, 0, 3, 168, 600, 198, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 654, 0, 3, 171, 606, 207, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 672, 0, 3, 174, 612, 216, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 690, 0, 3, 177, 618, 225, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 708, 0, 3, 180, 624, 234, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 726, 0, 3, 183, 630, 243, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 744, 0, 3, 189, 636, 44, 50, 270, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 780, 0, 3, 198, 654, 50, 56, 288, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 816, 0, 3, 207, 672, 56, 62, 306, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 852, 0, 3, 216, 690, 62, 68, 324, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 888, 0, 3, 225, 708, 68, 74, 342, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 924, 0, 3, 234, 726, 74, 80, 360, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 960, 0, 3, 270, 780, 92, 102, 438,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1020, 0, 3, 288, 816, 102, 112, 468,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1080, 0, 3, 306, 852, 112, 122, 498,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1140, 0, 3, 324, 888, 122, 132, 528,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1200, 0, 3, 342, 924, 132, 142, 558,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1260, 3, 162, 165, 594, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1270, 3, 165, 168, 600, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1280, 3, 168, 171, 606, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1290, 3, 171, 174, 612, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1300, 3, 174, 177, 618, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1310, 3, 177, 180, 624, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1320, 3, 180, 183, 630, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1330, 0, 3, 600, 1280, 654, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1360, 0, 3, 606, 1290, 672, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1390, 0, 3, 612, 1300, 690, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1420, 0, 3, 618, 1310, 708, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1450, 0, 3, 624, 1320, 726, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1480, 0, 3, 636, 1330, 252, 270, 780,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1540, 0, 3, 654, 1360, 270, 288, 816,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1600, 0, 3, 672, 1390, 288, 306, 852,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1660, 0, 3, 690, 1420, 306, 324, 888,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1720, 0, 3, 708, 1450, 324, 342, 924,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1780, 0, 3, 744, 1480, 378, 408, 960,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1880, 0, 3, 780, 1540, 408, 438, 1020,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1980, 0, 3, 816, 1600, 438, 468, 1080,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2080, 0, 3, 852, 1660, 468, 498, 1140,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2180, 0, 3, 888, 1720, 498, 528, 1200,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2280, 3, 588, 594, 1270, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2295, 3, 594, 600, 1280, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2310, 3, 600, 606, 1290, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2325, 3, 606, 612, 1300, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2340, 3, 612, 618, 1310, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2355, 3, 618, 624, 1320, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2370, 0, 3, 1280, 2310, 1360, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2415, 0, 3, 1290, 2325, 1390, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2460, 0, 3, 1300, 2340, 1420, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2505, 0, 3, 1310, 2355, 1450, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2550, 0, 3, 1330, 2370, 744, 780, 1540,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2640, 0, 3, 1360, 2415, 780, 816, 1600,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2730, 0, 3, 1390, 2460, 816, 852, 1660,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2820, 0, 3, 1420, 2505, 852, 888, 1720,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2910, 0, 3, 1540, 2640, 960, 1020, 1980,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3060, 0, 3, 1600, 2730, 1020, 1080,
                                                 2080, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3210, 0, 3, 1660, 2820, 1080, 1140,
                                                 2180, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3360, 3, 1260, 1270, 2295, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3381, 3, 1270, 1280, 2310, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3402, 3, 1280, 1290, 2325, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3423, 3, 1290, 1300, 2340, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3444, 3, 1300, 1310, 2355, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3465, 0, 3, 2310, 3402, 2415, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3528, 0, 3, 2325, 3423, 2460, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3591, 0, 3, 2340, 3444, 2505, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3654, 0, 3, 2370, 3465, 1480, 1540,
                                                 2640, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3780, 0, 3, 2415, 3528, 1540, 1600,
                                                 2730, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3906, 0, 3, 2460, 3591, 1600, 1660,
                                                 2820, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4032, 0, 3, 2550, 3654, 1780, 1880,
                                                 2910, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4242, 0, 3, 2640, 3780, 1880, 1980,
                                                 3060, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4452, 0, 3, 2730, 3906, 1980, 2080,
                                                 3210, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4662, 3, 2280, 2295, 3381, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4690, 3, 2295, 2310, 3402, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4718, 3, 2310, 2325, 3423, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 4746, 3, 2325, 2340, 3444, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 4774, 0, 3, 3381, 4690, 3465, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4858, 0, 3, 3402, 4718, 3528, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 4942, 0, 3, 3423, 4746, 3591, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 5026, 0, 3, 3465, 4858, 2550, 2640,
                                                 3780, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 5194, 0, 3, 3528, 4942, 2640, 2730,
                                                 3906, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 5362, 0, 3, 3780, 5194, 2910, 3060,
                                                 4452, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5642, 3, 3360, 3381, 4690, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 5678, 3, 3402, 3423, 4746, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 5714, 0, 3, 4662, 5642, 4774, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 5822, 0, 3, 4718, 5678, 4942, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 5930, 0, 3, 4858, 5822, 3654, 3780,
                                                 5194, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 6146, 0, 3, 5026, 5930, 4032, 4242,
                                                 5362, ncols, alpha, beta, p);

            compute_prim_geom_10_dk_electron_repulsion_0(buffer, 6506, 5714, 6146, ncols,
                                                         alpha);

            compute_prim_geom_10_dk_electron_repulsion_1(buffer, 6722, 5714, 6146, ncols,
                                                         alpha);

            compute_prim_geom_10_dk_electron_repulsion_2(buffer, 6938, 5714, 6146, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 7154, 6506, 648, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 7802, 7154, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 7802, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7802, 7370, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 75 * nvalues, nvalues, buffer, 7802, 15, nmax);

    simdtrf::transform_k_inner(buffer, 7802, 7586, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 150 * nvalues, nvalues, buffer, 7802, 15, nmax);
}

}  // namespace simdt2ceri

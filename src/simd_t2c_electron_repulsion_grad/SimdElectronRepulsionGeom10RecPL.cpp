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


#include "SimdElectronRepulsionGeom10RecPL.hpp"

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
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
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
#include "SimdGeometryP1.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_pl_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_pl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5613, 5157, 405, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 99, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 102, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 105, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 108, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 111, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 114, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 117, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 120, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 123, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 126, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 135, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 144, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 153, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 162, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 171, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 180, 3, 15, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 189, 3, 16, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 198, 0, 3, 21, 135, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 216, 0, 3, 24, 144, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 234, 0, 3, 27, 153, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 252, 0, 3, 30, 162, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 270, 0, 3, 33, 171, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 288, 0, 3, 36, 180, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 306, 0, 3, 39, 189, 93, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 324, 3, 7, 8, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 330, 3, 8, 9, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 336, 3, 9, 10, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 342, 3, 10, 11, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 348, 3, 11, 12, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 354, 3, 12, 13, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 360, 3, 13, 14, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 366, 3, 14, 15, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 372, 3, 15, 16, 123, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 378, 0, 3, 102, 336, 135, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 396, 0, 3, 105, 342, 144, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 414, 0, 3, 108, 348, 153, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 432, 0, 3, 111, 354, 162, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 450, 0, 3, 114, 360, 171, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 468, 0, 3, 117, 366, 180, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 486, 0, 3, 120, 372, 189, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 504, 0, 3, 126, 378, 45, 51, 198, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 540, 0, 3, 135, 396, 51, 57, 216, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 576, 0, 3, 144, 414, 57, 63, 234, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 612, 0, 3, 153, 432, 63, 69, 252, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 648, 0, 3, 162, 450, 69, 75, 270, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 684, 0, 3, 171, 468, 75, 81, 288, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 720, 0, 3, 180, 486, 81, 87, 306, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 756, 3, 99, 102, 336, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 766, 3, 102, 105, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 776, 3, 105, 108, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 786, 3, 108, 111, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 796, 3, 111, 114, 360, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 806, 3, 114, 117, 366, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 816, 3, 117, 120, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 826, 0, 3, 336, 766, 396, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 856, 0, 3, 342, 776, 414, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 886, 0, 3, 348, 786, 432, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 916, 0, 3, 354, 796, 450, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 946, 0, 3, 360, 806, 468, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 976, 0, 3, 366, 816, 486, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1006, 0, 3, 396, 856, 198, 216, 576,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1066, 0, 3, 414, 886, 216, 234, 612,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1126, 0, 3, 432, 916, 234, 252, 648,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1186, 0, 3, 450, 946, 252, 270, 684,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1246, 0, 3, 468, 976, 270, 288, 720,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1306, 3, 324, 330, 756, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1321, 3, 330, 336, 766, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1336, 3, 336, 342, 776, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1351, 3, 342, 348, 786, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1366, 3, 348, 354, 796, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1381, 3, 354, 360, 806, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1396, 3, 360, 366, 816, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1411, 0, 3, 766, 1336, 856, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1456, 0, 3, 776, 1351, 886, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1501, 0, 3, 786, 1366, 916, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1546, 0, 3, 796, 1381, 946, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1591, 0, 3, 806, 1396, 976, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1636, 0, 3, 826, 1411, 504, 540, 1006,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1726, 0, 3, 856, 1456, 540, 576, 1066,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1816, 0, 3, 886, 1501, 576, 612, 1126,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1906, 0, 3, 916, 1546, 612, 648, 1186,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1996, 0, 3, 946, 1591, 648, 684, 1246,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2086, 3, 756, 766, 1336, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2107, 3, 766, 776, 1351, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2128, 3, 776, 786, 1366, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2149, 3, 786, 796, 1381, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2170, 3, 796, 806, 1396, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2191, 0, 3, 1336, 2107, 1456, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2254, 0, 3, 1351, 2128, 1501, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2317, 0, 3, 1366, 2149, 1546, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2380, 0, 3, 1381, 2170, 1591, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2443, 0, 3, 1456, 2254, 1006, 1066,
                                                 1816, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2569, 0, 3, 1501, 2317, 1066, 1126,
                                                 1906, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2695, 0, 3, 1546, 2380, 1126, 1186,
                                                 1996, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2821, 3, 1306, 1321, 2086, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2849, 3, 1321, 1336, 2107, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2877, 3, 1336, 1351, 2128, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2905, 3, 1351, 1366, 2149, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2933, 3, 1366, 1381, 2170, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 2961, 0, 3, 2107, 2877, 2254, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 3045, 0, 3, 2128, 2905, 2317, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 3129, 0, 3, 2149, 2933, 2380, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3213, 0, 3, 2191, 2961, 1636, 1726,
                                                 2443, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3381, 0, 3, 2254, 3045, 1726, 1816,
                                                 2569, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3549, 0, 3, 2317, 3129, 1816, 1906,
                                                 2695, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 3717, 3, 2086, 2107, 2877, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 3753, 3, 2107, 2128, 2905, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 3789, 3, 2128, 2149, 2933, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 3825, 0, 3, 2877, 3753, 3045, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 3933, 0, 3, 2905, 3789, 3129, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 4041, 0, 3, 3045, 3933, 2443, 2569,
                                                 3549, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 4257, 3, 2821, 2849, 3717, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 4302, 3, 2877, 2905, 3789, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 4347, 0, 3, 3753, 4302, 3933, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 4482, 0, 3, 3825, 4347, 3213, 3381,
                                                 4041, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 4752, 4257, 4482, 1, 45, ncols, alpha);

            simdgeo::geom_p_y(buffer, 4887, 4257, 4482, 1, 45, ncols, alpha);

            simdgeo::geom_p_z(buffer, 5022, 4257, 4482, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 5157, 4752, 405, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 5562, 5157, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 5562, 17, nmax);

    simdtrf::transform_l_inner(buffer, 5562, 5292, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 51 * nvalues, nvalues, buffer, 5562, 17, nmax);

    simdtrf::transform_l_inner(buffer, 5562, 5427, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 102 * nvalues, nvalues, buffer, 5562, 17, nmax);
}

}  // namespace simdt2ceri

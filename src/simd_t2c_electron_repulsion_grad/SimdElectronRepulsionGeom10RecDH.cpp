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


#include "SimdElectronRepulsionGeom10RecDH.hpp"

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
compute_geom_10_dh_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_dh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3331, 2887, 378, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 15, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 36, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 42, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 72, 0, 15, 18, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 82, 0, 18, 21, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 21, 24, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 24, 27, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 27, 30, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 122, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 125, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 128, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 131, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 134, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 137, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 140, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 143, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 152, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 161, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 170, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 179, 3, 13, 33, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 188, 0, 3, 18, 143, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 206, 0, 3, 21, 152, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 224, 0, 3, 24, 161, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 242, 0, 3, 27, 170, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 260, 0, 3, 30, 179, 66, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 278, 0, 3, 36, 188, 72, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 308, 0, 3, 42, 206, 82, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 338, 0, 3, 48, 224, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 368, 0, 3, 54, 242, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 398, 0, 3, 60, 260, 112, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 428, 3, 7, 8, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 434, 3, 8, 9, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 440, 3, 9, 10, 131, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 446, 3, 10, 11, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 452, 3, 11, 12, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 458, 3, 12, 13, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 464, 0, 3, 128, 440, 152, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 482, 0, 3, 131, 446, 161, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 500, 0, 3, 134, 452, 170, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 518, 0, 3, 137, 458, 179, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 536, 0, 3, 143, 464, 36, 42, 206, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 572, 0, 3, 152, 482, 42, 48, 224, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 608, 0, 3, 161, 500, 48, 54, 242, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 644, 0, 3, 170, 518, 54, 60, 260, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 680, 0, 3, 206, 572, 72, 82, 338, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 740, 0, 3, 224, 608, 82, 92, 368, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 800, 0, 3, 242, 644, 92, 102, 398,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 860, 3, 122, 125, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 870, 3, 125, 128, 440, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 880, 3, 128, 131, 446, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 890, 3, 131, 134, 452, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 900, 3, 134, 137, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 910, 0, 3, 440, 880, 482, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 940, 0, 3, 446, 890, 500, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 970, 0, 3, 452, 900, 518, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1000, 0, 3, 464, 910, 188, 206, 572,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1060, 0, 3, 482, 940, 206, 224, 608,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1120, 0, 3, 500, 970, 224, 242, 644,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1180, 0, 3, 536, 1000, 278, 308, 680,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1280, 0, 3, 572, 1060, 308, 338, 740,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1380, 0, 3, 608, 1120, 338, 368, 800,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1480, 3, 428, 434, 870, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1495, 3, 434, 440, 880, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1510, 3, 440, 446, 890, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1525, 3, 446, 452, 900, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1540, 0, 3, 870, 1495, 910, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1585, 0, 3, 880, 1510, 940, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1630, 0, 3, 890, 1525, 970, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1675, 0, 3, 910, 1585, 536, 572, 1060,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1765, 0, 3, 940, 1630, 572, 608, 1120,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 1855, 0, 3, 1060, 1765, 680, 740, 1380,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2005, 3, 860, 870, 1495, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2026, 3, 880, 890, 1525, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2047, 0, 3, 1480, 2005, 1540, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2110, 0, 3, 1510, 2026, 1630, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2173, 0, 3, 1585, 2110, 1000, 1060,
                                                 1765, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 2299, 0, 3, 1675, 2173, 1180, 1280,
                                                 1855, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 2509, 2047, 2299, 1, 21, ncols, alpha);

            simdgeo::geom_d_y(buffer, 2635, 2047, 2299, 1, 21, ncols, alpha);

            simdgeo::geom_d_z(buffer, 2761, 2047, 2299, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 2887, 2509, 378, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 3265, 2887, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 3265, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3265, 3013, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 55 * nvalues, nvalues, buffer, 3265, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3265, 3139, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 110 * nvalues, nvalues, buffer, 3265, 11, nmax);
}

}  // namespace simdt2ceri

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


#include "SimdElectronRepulsionRecFH.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2632, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 122, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 125, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 128, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 131, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 134, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 137, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 146, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 155, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 164, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 173, 3, 13, 33, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 182, 0, 3, 18, 137, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 200, 0, 3, 21, 146, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 218, 0, 3, 24, 155, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 236, 0, 3, 27, 164, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 254, 0, 3, 30, 173, 66, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 272, 0, 3, 36, 182, 72, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 302, 0, 3, 42, 200, 82, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 332, 0, 3, 48, 218, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 362, 0, 3, 54, 236, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 392, 0, 3, 60, 254, 112, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 422, 3, 9, 10, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 428, 3, 10, 11, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 434, 3, 11, 12, 131, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 440, 3, 12, 13, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 446, 0, 3, 122, 422, 146, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 464, 0, 3, 125, 428, 155, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 482, 0, 3, 128, 434, 164, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 500, 0, 3, 131, 440, 173, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 518, 0, 3, 137, 446, 36, 42, 200, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 554, 0, 3, 146, 464, 42, 48, 218, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 590, 0, 3, 155, 482, 48, 54, 236, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 626, 0, 3, 164, 500, 54, 60, 254, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 662, 0, 3, 200, 554, 72, 82, 332, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 722, 0, 3, 218, 590, 82, 92, 362, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 782, 0, 3, 236, 626, 92, 102, 392,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 842, 3, 122, 125, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 852, 3, 125, 128, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 862, 3, 128, 131, 440, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 872, 0, 3, 422, 842, 464, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 902, 0, 3, 428, 852, 482, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 932, 0, 3, 434, 862, 500, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 962, 0, 3, 446, 872, 182, 200, 554,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1022, 0, 3, 464, 902, 200, 218, 590,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1082, 0, 3, 482, 932, 218, 236, 626,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1142, 0, 3, 518, 962, 272, 302, 662,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1242, 0, 3, 554, 1022, 302, 332, 722,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1342, 0, 3, 590, 1082, 332, 362, 782,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1442, 3, 422, 428, 852, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1457, 3, 428, 434, 862, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1472, 0, 3, 842, 1442, 902, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1517, 0, 3, 852, 1457, 932, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1562, 0, 3, 872, 1472, 518, 554, 1022,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1652, 0, 3, 902, 1517, 554, 590, 1082,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 1742, 0, 3, 1022, 1652, 662, 722, 1342,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1892, 3, 842, 852, 1457, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1913, 0, 3, 1442, 1892, 1517, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 1976, 0, 3, 1472, 1913, 962, 1022, 1652,
                                                 ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 2102, 0, 3, 1562, 1976, 1142, 1242,
                                                 1742, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2312, 2102, 210, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 2522, 2312, 10, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 2522, 11, nmax);
}

}  // namespace simdt2ceri

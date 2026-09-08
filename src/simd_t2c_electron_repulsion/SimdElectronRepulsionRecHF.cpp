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


#include "SimdElectronRepulsionRecHF.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hf_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2912, nvalues);

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

            compute_prim_gs_electron_repulsion_0(buffer, 122, 0, 36, 42, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 137, 0, 42, 48, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 152, 0, 48, 54, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 167, 0, 54, 60, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 182, 0, 72, 82, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 203, 0, 82, 92, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 224, 0, 92, 102, 167, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 245, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 248, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 251, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 254, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 257, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 260, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 269, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 278, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 287, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 296, 3, 13, 33, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 305, 0, 3, 18, 260, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 323, 0, 3, 21, 269, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 341, 0, 3, 24, 278, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 359, 0, 3, 27, 287, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 377, 0, 3, 30, 296, 66, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 395, 0, 3, 36, 305, 72, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 425, 0, 3, 42, 323, 82, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 455, 0, 3, 48, 341, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 485, 0, 3, 54, 359, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 515, 0, 3, 60, 377, 112, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 545, 0, 3, 82, 455, 137, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 590, 0, 3, 92, 485, 152, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 635, 0, 3, 102, 515, 167, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 680, 0, 3, 122, 545, 182, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 743, 0, 3, 137, 590, 203, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 806, 0, 3, 152, 635, 224, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 869, 3, 9, 10, 248, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 875, 3, 10, 11, 251, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 881, 3, 11, 12, 254, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 887, 3, 12, 13, 257, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 893, 0, 3, 245, 869, 269, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 911, 0, 3, 248, 875, 278, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 929, 0, 3, 251, 881, 287, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 947, 0, 3, 254, 887, 296, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 965, 0, 3, 260, 893, 36, 42, 323, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1001, 0, 3, 269, 911, 42, 48, 341,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1037, 0, 3, 278, 929, 48, 54, 359,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1073, 0, 3, 287, 947, 54, 60, 377,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1109, 0, 3, 323, 1001, 72, 82, 455,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1169, 0, 3, 341, 1037, 82, 92, 485,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1229, 0, 3, 359, 1073, 92, 102, 515,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1289, 0, 3, 965, 1001, 455, 1169, 122,
                                                 137, 590, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1379, 0, 3, 1001, 1037, 485, 1229, 137,
                                                 152, 635, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 1469, 0, 3, 1109, 1169, 590, 1379, 182,
                                                 203, 806, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1595, 3, 245, 248, 875, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1605, 3, 248, 251, 881, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1615, 3, 251, 254, 887, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1625, 0, 3, 869, 1595, 911, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1655, 0, 3, 875, 1605, 929, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1685, 0, 3, 881, 1615, 947, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1715, 0, 3, 893, 1625, 305, 323, 1001,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1775, 0, 3, 911, 1655, 323, 341, 1037,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1835, 0, 3, 929, 1685, 341, 359, 1073,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1895, 0, 3, 965, 1715, 395, 425, 1109,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1995, 0, 3, 1001, 1775, 425, 455, 1169,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2095, 0, 3, 1037, 1835, 455, 485, 1229,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2195, 0, 3, 1715, 1775, 1169, 2095, 545,
                                                 590, 1379, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 2345, 0, 3, 1895, 1995, 1289, 2195, 680,
                                                 743, 1469, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2555, 2345, 210, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 2765, 2555, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 2765, 7, nmax);
}

}  // namespace simdt2ceri

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


#include "SimdElectronRepulsionRecIF.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_if_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_if_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4881, 4405, 280, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 8, 9, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 9, 10, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 10, 11, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 11, 12, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 12, 13, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 13, 14, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 82, 0, 16, 19, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 19, 22, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 22, 25, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 25, 28, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 28, 31, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 31, 34, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 142, 0, 40, 46, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 157, 0, 46, 52, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 172, 0, 52, 58, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 187, 0, 58, 64, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 202, 0, 64, 70, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 217, 0, 82, 92, 157, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 238, 0, 92, 102, 172, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 259, 0, 102, 112, 187, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 280, 0, 112, 122, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 301, 0, 142, 157, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 329, 0, 157, 172, 259, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 357, 0, 172, 187, 280, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 385, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 388, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 391, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 394, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 397, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 400, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 403, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 406, 3, 8, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 415, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 424, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 433, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 442, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 451, 3, 13, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 460, 3, 14, 37, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 469, 0, 3, 16, 406, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 487, 0, 3, 19, 415, 46, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 505, 0, 3, 22, 424, 52, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 523, 0, 3, 25, 433, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 541, 0, 3, 28, 442, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 559, 0, 3, 31, 451, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 577, 0, 3, 34, 460, 76, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 595, 0, 3, 46, 505, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 625, 0, 3, 52, 523, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 655, 0, 3, 58, 541, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 685, 0, 3, 64, 559, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 715, 0, 3, 70, 577, 132, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 745, 0, 3, 82, 595, 142, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 790, 0, 3, 92, 625, 157, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 835, 0, 3, 102, 655, 172, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 880, 0, 3, 112, 685, 187, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 925, 0, 3, 122, 715, 202, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 970, 0, 3, 157, 835, 238, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1033, 0, 3, 172, 880, 259, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1096, 0, 3, 187, 925, 280, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1159, 0, 3, 217, 970, 301, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1243, 0, 3, 238, 1033, 329, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1327, 0, 3, 259, 1096, 357, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1411, 3, 8, 9, 388, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1417, 3, 9, 10, 391, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1423, 3, 10, 11, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1429, 3, 11, 12, 397, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1435, 3, 12, 13, 400, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1441, 3, 13, 14, 403, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1447, 0, 3, 385, 1411, 415, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1465, 0, 3, 388, 1417, 424, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1483, 0, 3, 391, 1423, 433, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1501, 0, 3, 394, 1429, 442, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1519, 0, 3, 397, 1435, 451, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1537, 0, 3, 400, 1441, 460, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1555, 0, 3, 415, 1465, 40, 46, 505,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1591, 0, 3, 424, 1483, 46, 52, 523,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1627, 0, 3, 433, 1501, 52, 58, 541,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1663, 0, 3, 442, 1519, 58, 64, 559,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1699, 0, 3, 451, 1537, 64, 70, 577,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1735, 0, 3, 505, 1591, 82, 92, 625,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1795, 0, 3, 523, 1627, 92, 102, 655,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1855, 0, 3, 541, 1663, 102, 112, 685,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1915, 0, 3, 559, 1699, 112, 122, 715,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1975, 0, 3, 1555, 1591, 625, 1795, 142,
                                                 157, 835, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2065, 0, 3, 1591, 1627, 655, 1855, 157,
                                                 172, 880, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2155, 0, 3, 1627, 1663, 685, 1915, 172,
                                                 187, 925, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2245, 0, 3, 1735, 1795, 835, 2065, 217,
                                                 238, 1033, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2371, 0, 3, 1795, 1855, 880, 2155, 238,
                                                 259, 1096, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 2497, 0, 3, 1975, 2065, 1033, 2371, 301,
                                                 329, 1327, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2665, 3, 385, 388, 1417, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2675, 3, 388, 391, 1423, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2685, 3, 391, 394, 1429, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2695, 3, 394, 397, 1435, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2705, 3, 397, 400, 1441, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2715, 0, 3, 1411, 2665, 1465, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2745, 0, 3, 1417, 2675, 1483, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2775, 0, 3, 1423, 2685, 1501, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2805, 0, 3, 1429, 2695, 1519, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2835, 0, 3, 1435, 2705, 1537, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2865, 0, 3, 1447, 2715, 469, 487, 1555,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2925, 0, 3, 1465, 2745, 487, 505, 1591,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2985, 0, 3, 1483, 2775, 505, 523, 1627,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3045, 0, 3, 1501, 2805, 523, 541, 1663,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3105, 0, 3, 1519, 2835, 541, 559, 1699,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3165, 0, 3, 1591, 2985, 595, 625, 1795,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3265, 0, 3, 1627, 3045, 625, 655, 1855,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3365, 0, 3, 1663, 3105, 655, 685, 1915,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3465, 0, 3, 2865, 2925, 1735, 3165, 745,
                                                 790, 1975, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3615, 0, 3, 2925, 2985, 1795, 3265, 790,
                                                 835, 2065, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3765, 0, 3, 2985, 3045, 1855, 3365, 835,
                                                 880, 2155, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 3915, 0, 3, 3165, 3265, 2065, 3765, 970,
                                                 1033, 2371, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 4125, 0, 3, 3465, 3615, 2245, 3915,
                                                 1159, 1243, 2497, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4405, 4125, 280, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 4685, 4405, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 4685, 7, nmax);
}

}  // namespace simdt2ceri

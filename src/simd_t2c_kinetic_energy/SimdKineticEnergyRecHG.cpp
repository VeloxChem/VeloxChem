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


#include "SimdKineticEnergyRecHG.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdKineticEnergyVrrRecDD.hpp"
#include "SimdKineticEnergyVrrRecDF.hpp"
#include "SimdKineticEnergyVrrRecDG.hpp"
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFG.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGG.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHG.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPG.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSG.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDG.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFG.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGG.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHG.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPG.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hg_kinetic_energy(double               *values,
                          const size_t          nvalues,
                          const CBasisFunction &bra,
                          const CBasisFunction &ket,
                          const CSimdMatrix    &coordinates,
                          const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hg_kinetic_energy: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as every integral is a sum over them
    // and the error of a sum is bounded by the number of its terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound, threshold / static_cast<double>(nprims));

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 4430);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 99 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    errors::assertMsgCritical(dimensions.size() == nprim_a * nprim_b,
                              std::string("Dimensions do not match the pairs of primitives"));

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto fpi = mathconst::pi_value() / p;

            const auto fovl = a_norms[i] * b_norms[j] * fpi * std::sqrt(fpi);

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdovl::compute_prim_ss_overlap(buffer, coordinates, 6, ncols, fovl, mu);

            simdovl::compute_prim_sp_overlap_0(buffer, 7, 3, 6, ncols);

            simdovl::compute_prim_sd_overlap_0(buffer, 10, 3, 6, 7, ncols, p);

            simdovl::compute_prim_sf_overlap_0(buffer, 16, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_0(buffer, 26, 3, 10, 16, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 41, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 42, 3, 7, 41, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_0(buffer, 45, 3, 6, 10, 41, 42, ncols, alpha, beta,
                                             p);

            compute_prim_sf_kinetic_energy_0(buffer, 51, 3, 7, 16, 42, 45, ncols, alpha, beta,
                                             p);

            compute_prim_sg_kinetic_energy_0(buffer, 61, 3, 10, 26, 45, 51, ncols, alpha, beta,
                                             p);

            simdovl::compute_prim_ps_overlap_0(buffer, 76, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 79, 3, 6, 76, ncols, p);

            simdovl::compute_prim_pd_overlap_0(buffer, 88, 0, 3, 7, 10, 79, ncols, p);

            simdovl::compute_prim_pf_overlap_0(buffer, 106, 0, 3, 10, 16, 79, 88, ncols, p);

            simdovl::compute_prim_pg_overlap_0(buffer, 136, 0, 3, 16, 26, 88, 106, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 181, 0, 41, 76, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 184, 3, 41, 79, 181, ncols, alpha, beta,
                                             p);

            compute_prim_pd_kinetic_energy_0(buffer, 193, 0, 3, 42, 45, 88, 184, ncols, alpha,
                                             beta, p);

            compute_prim_pf_kinetic_energy_0(buffer, 211, 0, 3, 45, 51, 106, 193, ncols, alpha,
                                             beta, p);

            compute_prim_pg_kinetic_energy_0(buffer, 241, 0, 3, 51, 61, 136, 211, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_ds_overlap_0(buffer, 286, 0, 6, 76, ncols, p);

            simdovl::compute_prim_dp_overlap_0(buffer, 292, 0, 3, 76, 79, 286, ncols, p);

            simdovl::compute_prim_dd_overlap_0(buffer, 310, 0, 3, 79, 88, 286, 292, ncols, p);

            simdovl::compute_prim_df_overlap_0(buffer, 346, 0, 3, 88, 106, 292, 310, ncols, p);

            simdovl::compute_prim_dg_overlap_0(buffer, 406, 0, 3, 106, 136, 310, 346, ncols, p);

            compute_prim_ds_kinetic_energy_0(buffer, 496, 0, 6, 41, 181, 286, ncols, alpha, beta,
                                             p);

            compute_prim_dp_kinetic_energy_0(buffer, 502, 0, 3, 181, 184, 292, 496, ncols, alpha,
                                             beta, p);

            compute_prim_dd_kinetic_energy_0(buffer, 520, 0, 3, 184, 193, 286, 310, 496, 502,
                                             ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_0(buffer, 556, 0, 3, 193, 211, 292, 346, 502, 520,
                                             ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_0(buffer, 616, 0, 3, 211, 241, 310, 406, 520, 556,
                                             ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_0(buffer, 706, 0, 76, 286, ncols, p);

            simdovl::compute_prim_fp_overlap_0(buffer, 716, 0, 3, 286, 292, 706, ncols, p);

            simdovl::compute_prim_fd_overlap_0(buffer, 746, 0, 3, 88, 292, 310, 706, 716, ncols,
                                               p);

            simdovl::compute_prim_ff_overlap_0(buffer, 806, 0, 3, 106, 310, 346, 716, 746, ncols,
                                               p);

            simdovl::compute_prim_fg_overlap_0(buffer, 906, 0, 3, 136, 346, 406, 746, 806, ncols,
                                               p);

            compute_prim_fs_kinetic_energy_0(buffer, 1056, 0, 76, 181, 496, 706, ncols, alpha,
                                             beta, p);

            compute_prim_fp_kinetic_energy_0(buffer, 1066, 0, 3, 496, 502, 716, 1056, ncols,
                                             alpha, beta, p);

            compute_prim_fd_kinetic_energy_0(buffer, 1096, 0, 3, 88, 193, 502, 520, 706, 746,
                                             1056, 1066, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_0(buffer, 1156, 0, 3, 106, 211, 520, 556, 716, 806,
                                             1066, 1096, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_0(buffer, 1256, 0, 3, 136, 241, 556, 616, 746, 906,
                                             1096, 1156, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_0(buffer, 1406, 0, 286, 706, ncols, p);

            simdovl::compute_prim_gp_overlap_0(buffer, 1421, 0, 3, 706, 716, 1406, ncols, p);

            simdovl::compute_prim_gd_overlap_0(buffer, 1466, 0, 3, 310, 716, 746, 1406, 1421,
                                               ncols, p);

            simdovl::compute_prim_gf_overlap_0(buffer, 1556, 0, 3, 346, 746, 806, 1421, 1466,
                                               ncols, p);

            simdovl::compute_prim_gg_overlap_0(buffer, 1706, 0, 3, 406, 806, 906, 1466, 1556,
                                               ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 1931, 0, 286, 496, 1056, 1406, ncols, alpha,
                                             beta, p);

            compute_prim_gp_kinetic_energy_0(buffer, 1946, 0, 3, 1056, 1066, 1421, 1931, ncols,
                                             alpha, beta, p);

            compute_prim_gd_kinetic_energy_0(buffer, 1991, 0, 3, 310, 520, 1066, 1096, 1406,
                                             1466, 1931, 1946, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_0(buffer, 2081, 0, 3, 346, 556, 1096, 1156, 1421,
                                             1556, 1946, 1991, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_0(buffer, 2231, 0, 3, 406, 616, 1156, 1256, 1466,
                                             1706, 1991, 2081, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_0(buffer, 2456, 0, 706, 1406, ncols, p);

            simdovl::compute_prim_hp_overlap_0(buffer, 2477, 0, 3, 1406, 1421, 2456, ncols, p);

            simdovl::compute_prim_hd_overlap_0(buffer, 2540, 0, 3, 746, 1421, 1466, 2456, 2477,
                                               ncols, p);

            simdovl::compute_prim_hf_overlap_0(buffer, 2666, 0, 3, 806, 1466, 1556, 2477, 2540,
                                               ncols, p);

            simdovl::compute_prim_hg_overlap_0(buffer, 2876, 0, 3, 906, 1556, 1706, 2540, 2666,
                                               ncols, p);

            compute_prim_hs_kinetic_energy_0(buffer, 3191, 0, 706, 1056, 1931, 2456, ncols,
                                             alpha, beta, p);

            compute_prim_hp_kinetic_energy_0(buffer, 3212, 0, 3, 1931, 1946, 2477, 3191, ncols,
                                             alpha, beta, p);

            compute_prim_hd_kinetic_energy_0(buffer, 3275, 0, 3, 746, 1096, 1946, 1991, 2456,
                                             2540, 3191, 3212, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_0(buffer, 3401, 0, 3, 806, 1156, 1991, 2081, 2477,
                                             2666, 3212, 3275, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_0(buffer, 3611, 0, 3, 906, 1256, 2081, 2231, 2540,
                                             2876, 3275, 3401, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3926, 3611, 315, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 4241, 3926, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 4241, 9, nmax);

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

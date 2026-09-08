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


#include "SimdKineticEnergyRecIG.hpp"

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
#include "SimdKineticEnergyVrrRecID.hpp"
#include "SimdKineticEnergyVrrRecIF.hpp"
#include "SimdKineticEnergyVrrRecIG.hpp"
#include "SimdKineticEnergyVrrRecIP.hpp"
#include "SimdKineticEnergyVrrRecIS.hpp"
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
#include "SimdOverlapVrrRecID.hpp"
#include "SimdOverlapVrrRecIF.hpp"
#include "SimdOverlapVrrRecIG.hpp"
#include "SimdOverlapVrrRecIP.hpp"
#include "SimdOverlapVrrRecIS.hpp"
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
#include "SimdTransformIG.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ig_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ig_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 3179);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

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

            simdovl::compute_prim_sd_overlap_1(buffer, 10, 3, 6, 7, ncols, p);

            simdovl::compute_prim_sf_overlap_8(buffer, 13, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_9(buffer, 18, 3, 10, 13, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 21, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 22, 3, 7, 21, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 25, 3, 6, 10, 21, 22, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_6(buffer, 28, 3, 7, 13, 22, 25, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_6(buffer, 33, 3, 10, 18, 25, 28, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 36, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 39, 3, 6, 36, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 42, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_15(buffer, 45, 0, 3, 10, 13, 42, ncols, p);

            simdovl::compute_prim_pg_overlap_12(buffer, 49, 0, 3, 13, 18, 42, 45, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 52, 0, 21, 36, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 55, 3, 21, 39, 52, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 58, 0, 22, 25, 42, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_15(buffer, 61, 0, 3, 25, 28, 45, 58, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_9(buffer, 65, 0, 28, 33, 49, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 68, 0, 6, 36, ncols, p);

            simdovl::compute_prim_dp_overlap_11(buffer, 71, 3, 36, 68, ncols, p);

            simdovl::compute_prim_dd_overlap_17(buffer, 77, 0, 3, 39, 42, 68, 71, ncols, p);

            simdovl::compute_prim_df_overlap_17(buffer, 86, 0, 3, 42, 45, 71, 77, ncols, p);

            simdovl::compute_prim_dg_overlap_12(buffer, 96, 0, 3, 45, 49, 77, 86, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 107, 0, 6, 21, 52, 68, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_12(buffer, 110, 3, 52, 71, 107, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_17(buffer, 116, 0, 3, 55, 58, 68, 77, 107, 110, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_17(buffer, 125, 0, 3, 58, 61, 71, 86, 110, 116, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_12(buffer, 135, 0, 3, 61, 65, 77, 96, 116, 125, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_3(buffer, 146, 0, 36, 68, ncols, p);

            simdovl::compute_prim_fp_overlap_14(buffer, 152, 3, 68, 146, ncols, p);

            simdovl::compute_prim_fd_overlap_14(buffer, 160, 0, 3, 42, 71, 77, 146, 152, ncols, p);

            simdovl::compute_prim_ff_overlap_13(buffer, 180, 0, 3, 45, 77, 86, 152, 160, ncols, p);

            simdovl::compute_prim_fg_overlap_9(buffer, 202, 0, 3, 49, 86, 96, 160, 180, ncols, p);

            compute_prim_fs_kinetic_energy_4(buffer, 223, 0, 36, 52, 107, 146, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_14(buffer, 229, 3, 107, 152, 223, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_16(buffer, 237, 0, 3, 42, 58, 110, 116, 146, 160, 223, 229, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_13(buffer, 257, 0, 3, 45, 61, 116, 125, 152, 180, 229, 237, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_9(buffer, 279, 0, 3, 49, 65, 125, 135, 160, 202, 237, 257, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_2(buffer, 300, 0, 68, 146, ncols, p);

            simdovl::compute_prim_gp_overlap_11(buffer, 311, 3, 146, 300, ncols, p);

            simdovl::compute_prim_gd_overlap_10(buffer, 322, 0, 3, 77, 152, 160, 300, 311, ncols, p);

            simdovl::compute_prim_gf_overlap_9(buffer, 358, 0, 3, 86, 160, 180, 311, 322, ncols, p);

            simdovl::compute_prim_gg_overlap_6(buffer, 411, 0, 3, 96, 180, 202, 322, 358, ncols, p);

            compute_prim_gs_kinetic_energy_3(buffer, 452, 0, 68, 107, 223, 300, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_11(buffer, 463, 3, 223, 311, 452, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_12(buffer, 474, 0, 3, 77, 116, 229, 237, 300, 322, 452, 463, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_9(buffer, 510, 0, 3, 86, 125, 237, 257, 311, 358, 463, 474, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_6(buffer, 563, 0, 3, 96, 135, 257, 279, 322, 411, 474, 510, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_10(buffer, 604, 0, 146, 300, ncols, p);

            simdovl::compute_prim_hp_overlap_7(buffer, 618, 0, 3, 300, 311, 604, ncols, p);

            simdovl::compute_prim_hd_overlap_6(buffer, 640, 0, 3, 160, 311, 322, 604, 618, ncols, p);

            simdovl::compute_prim_hf_overlap_5(buffer, 700, 0, 3, 180, 322, 358, 618, 640, ncols, p);

            simdovl::compute_prim_hg_overlap_3(buffer, 822, 0, 3, 202, 358, 411, 640, 700, ncols, p);

            compute_prim_hs_kinetic_energy_7(buffer, 935, 0, 146, 223, 452, 604, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_7(buffer, 949, 0, 3, 452, 463, 618, 935, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_7(buffer, 970, 0, 3, 160, 237, 463, 474, 604, 640, 935, 949, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_5(buffer, 1030, 0, 3, 180, 257, 474, 510, 618, 700, 949, 970, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_3(buffer, 1152, 0, 3, 202, 279, 510, 563, 640, 822, 970, 1030, ncols, alpha, beta, p);

            simdovl::compute_prim_is_overlap_6(buffer, 1265, 0, 300, 604, ncols, p);

            simdovl::compute_prim_ip_overlap_3(buffer, 1280, 0, 3, 604, 618, 1265, ncols, p);

            simdovl::compute_prim_id_overlap_2(buffer, 1311, 0, 3, 322, 618, 640, 1265, 1280, ncols, p);

            simdovl::compute_prim_if_overlap_1(buffer, 1398, 0, 3, 358, 640, 700, 1280, 1311, ncols, p);

            simdovl::compute_prim_ig_overlap_0(buffer, 1598, 0, 3, 411, 700, 822, 1311, 1398, ncols, p);

            compute_prim_is_kinetic_energy_3(buffer, 2018, 0, 300, 452, 935, 1265, ncols, alpha, beta, p);

            compute_prim_ip_kinetic_energy_3(buffer, 2032, 3, 935, 1280, 2018, ncols, alpha, beta, p);

            compute_prim_id_kinetic_energy_2(buffer, 2061, 0, 3, 322, 474, 949, 970, 1265, 1311, 2018, 2032, ncols, alpha, beta, p);

            compute_prim_if_kinetic_energy_1(buffer, 2145, 0, 3, 358, 510, 970, 1030, 1280, 1398, 2032, 2061, ncols, alpha, beta, p);

            compute_prim_ig_kinetic_energy_0(buffer, 2339, 0, 3, 411, 563, 1030, 1152, 1311, 1598, 2061, 2145, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2759, 2339, 420, ncols);
        }
    }

    simdtrf::transform_ig(values, nvalues, buffer, 2759, nmax);

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

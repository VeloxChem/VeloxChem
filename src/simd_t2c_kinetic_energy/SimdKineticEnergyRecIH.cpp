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


#include "SimdKineticEnergyRecIH.hpp"

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
#include "SimdKineticEnergyVrrRecDH.hpp"
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFG.hpp"
#include "SimdKineticEnergyVrrRecFH.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGG.hpp"
#include "SimdKineticEnergyVrrRecGH.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHG.hpp"
#include "SimdKineticEnergyVrrRecHH.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecID.hpp"
#include "SimdKineticEnergyVrrRecIF.hpp"
#include "SimdKineticEnergyVrrRecIG.hpp"
#include "SimdKineticEnergyVrrRecIH.hpp"
#include "SimdKineticEnergyVrrRecIP.hpp"
#include "SimdKineticEnergyVrrRecIS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPG.hpp"
#include "SimdKineticEnergyVrrRecPH.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSG.hpp"
#include "SimdKineticEnergyVrrRecSH.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDG.hpp"
#include "SimdOverlapVrrRecDH.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFG.hpp"
#include "SimdOverlapVrrRecFH.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGG.hpp"
#include "SimdOverlapVrrRecGH.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHG.hpp"
#include "SimdOverlapVrrRecHH.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecID.hpp"
#include "SimdOverlapVrrRecIF.hpp"
#include "SimdOverlapVrrRecIG.hpp"
#include "SimdOverlapVrrRecIH.hpp"
#include "SimdOverlapVrrRecIP.hpp"
#include "SimdOverlapVrrRecIS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPG.hpp"
#include "SimdOverlapVrrRecPH.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSH.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformIH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ih_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ih_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 4816);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

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

            simdovl::compute_prim_sg_overlap_11(buffer, 18, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_12(buffer, 24, 3, 13, 18, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 28, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 29, 3, 7, 28, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 32, 3, 6, 10, 28, 29, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_6(buffer, 35, 3, 7, 13, 29, 32, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_8(buffer, 40, 3, 10, 18, 32, 35, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_9(buffer, 46, 3, 13, 24, 35, 40, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 50, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 53, 3, 6, 50, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 56, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_15(buffer, 59, 0, 3, 10, 13, 56, ncols, p);

            simdovl::compute_prim_pg_overlap_13(buffer, 63, 0, 3, 13, 18, 56, 59, ncols, p);

            simdovl::compute_prim_ph_overlap_8(buffer, 68, 0, 3, 18, 24, 59, 63, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 72, 0, 28, 50, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 75, 3, 28, 53, 72, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 78, 0, 29, 32, 56, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_15(buffer, 81, 0, 3, 32, 35, 59, 78, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_12(buffer, 85, 0, 3, 35, 40, 63, 81, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_8(buffer, 90, 0, 40, 46, 68, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 94, 0, 6, 50, ncols, p);

            simdovl::compute_prim_dp_overlap_11(buffer, 97, 3, 50, 94, ncols, p);

            simdovl::compute_prim_dd_overlap_17(buffer, 103, 0, 3, 53, 56, 94, 97, ncols, p);

            simdovl::compute_prim_df_overlap_18(buffer, 112, 0, 3, 56, 59, 97, 103, ncols, p);

            simdovl::compute_prim_dg_overlap_13(buffer, 126, 0, 3, 59, 63, 103, 112, ncols, p);

            simdovl::compute_prim_dh_overlap_8(buffer, 142, 0, 3, 63, 68, 112, 126, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 157, 0, 6, 28, 72, 94, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_12(buffer, 160, 3, 72, 97, 157, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_17(buffer, 166, 0, 3, 75, 78, 94, 103, 157, 160, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_18(buffer, 175, 0, 3, 78, 81, 97, 112, 160, 166, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_13(buffer, 189, 0, 3, 81, 85, 103, 126, 166, 175, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_8(buffer, 205, 0, 3, 85, 90, 112, 142, 175, 189, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_3(buffer, 220, 0, 50, 94, ncols, p);

            simdovl::compute_prim_fp_overlap_14(buffer, 226, 3, 94, 220, ncols, p);

            simdovl::compute_prim_fd_overlap_15(buffer, 234, 0, 3, 56, 97, 103, 220, 226, ncols, p);

            simdovl::compute_prim_ff_overlap_14(buffer, 251, 0, 3, 59, 103, 112, 226, 234, ncols, p);

            simdovl::compute_prim_fg_overlap_10(buffer, 282, 0, 3, 63, 112, 126, 234, 251, ncols, p);

            simdovl::compute_prim_fh_overlap_6(buffer, 317, 0, 3, 68, 126, 142, 251, 282, ncols, p);

            compute_prim_fs_kinetic_energy_4(buffer, 345, 0, 50, 72, 157, 220, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_14(buffer, 351, 3, 157, 226, 345, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_17(buffer, 359, 0, 3, 56, 78, 160, 166, 220, 234, 345, 351, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_14(buffer, 376, 0, 3, 59, 81, 166, 175, 226, 251, 351, 359, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_10(buffer, 407, 0, 3, 63, 85, 175, 189, 234, 282, 359, 376, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_6(buffer, 442, 0, 3, 68, 90, 189, 205, 251, 317, 376, 407, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_2(buffer, 470, 0, 94, 220, ncols, p);

            simdovl::compute_prim_gp_overlap_11(buffer, 481, 3, 220, 470, ncols, p);

            simdovl::compute_prim_gd_overlap_11(buffer, 492, 0, 3, 103, 226, 234, 470, 481, ncols, p);

            simdovl::compute_prim_gf_overlap_10(buffer, 520, 0, 3, 112, 234, 251, 481, 492, ncols, p);

            simdovl::compute_prim_gg_overlap_7(buffer, 576, 0, 3, 126, 251, 282, 492, 520, ncols, p);

            simdovl::compute_prim_gh_overlap_4(buffer, 655, 0, 3, 142, 282, 317, 520, 576, ncols, p);

            compute_prim_gs_kinetic_energy_3(buffer, 713, 0, 94, 157, 345, 470, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_11(buffer, 724, 3, 345, 481, 713, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_13(buffer, 735, 0, 3, 103, 166, 351, 359, 470, 492, 713, 724, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_10(buffer, 763, 0, 3, 112, 175, 359, 376, 481, 520, 724, 735, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_7(buffer, 819, 0, 3, 126, 189, 376, 407, 492, 576, 735, 763, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_4(buffer, 898, 0, 3, 142, 205, 407, 442, 520, 655, 763, 819, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_10(buffer, 956, 0, 220, 470, ncols, p);

            simdovl::compute_prim_hp_overlap_7(buffer, 970, 0, 3, 470, 481, 956, ncols, p);

            simdovl::compute_prim_hd_overlap_7(buffer, 992, 0, 3, 234, 481, 492, 956, 970, ncols, p);

            simdovl::compute_prim_hf_overlap_6(buffer, 1037, 0, 3, 251, 492, 520, 970, 992, ncols, p);

            simdovl::compute_prim_hg_overlap_4(buffer, 1133, 0, 3, 282, 520, 576, 992, 1037, ncols, p);

            simdovl::compute_prim_hh_overlap_2(buffer, 1318, 0, 3, 317, 576, 655, 1037, 1133, ncols, p);

            compute_prim_hs_kinetic_energy_7(buffer, 1474, 0, 220, 345, 713, 956, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_7(buffer, 1488, 0, 3, 713, 724, 970, 1474, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_8(buffer, 1509, 0, 3, 234, 359, 724, 735, 956, 992, 1474, 1488, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_6(buffer, 1552, 0, 3, 251, 376, 735, 763, 970, 1037, 1488, 1509, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_4(buffer, 1647, 0, 3, 282, 407, 763, 819, 992, 1133, 1509, 1552, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_2(buffer, 1832, 0, 3, 317, 442, 819, 898, 1037, 1318, 1552, 1647, ncols, alpha, beta, p);

            simdovl::compute_prim_is_overlap_6(buffer, 1988, 0, 470, 956, ncols, p);

            simdovl::compute_prim_ip_overlap_3(buffer, 2003, 0, 3, 956, 970, 1988, ncols, p);

            simdovl::compute_prim_id_overlap_3(buffer, 2034, 0, 3, 492, 970, 992, 1988, 2003, ncols, p);

            simdovl::compute_prim_if_overlap_2(buffer, 2099, 0, 3, 520, 992, 1037, 2003, 2034, ncols, p);

            simdovl::compute_prim_ig_overlap_1(buffer, 2241, 0, 3, 576, 1037, 1133, 2034, 2099, ncols, p);

            simdovl::compute_prim_ih_overlap_0(buffer, 2539, 0, 3, 655, 1133, 1318, 2099, 2241, ncols, p);

            compute_prim_is_kinetic_energy_3(buffer, 3127, 0, 470, 713, 1474, 1988, ncols, alpha, beta, p);

            compute_prim_ip_kinetic_energy_3(buffer, 3141, 3, 1474, 2003, 3127, ncols, alpha, beta, p);

            compute_prim_id_kinetic_energy_3(buffer, 3170, 0, 3, 492, 735, 1488, 1509, 1988, 2034, 3127, 3141, ncols, alpha, beta, p);

            compute_prim_if_kinetic_energy_2(buffer, 3226, 0, 3, 520, 763, 1509, 1552, 2003, 2099, 3141, 3170, ncols, alpha, beta, p);

            compute_prim_ig_kinetic_energy_1(buffer, 3356, 0, 3, 576, 819, 1552, 1647, 2034, 2241, 3170, 3226, ncols, alpha, beta, p);

            compute_prim_ih_kinetic_energy_0(buffer, 3640, 0, 3, 655, 898, 1647, 1832, 2099, 2539, 3226, 3356, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4228, 3640, 588, ncols);
        }
    }

    simdtrf::transform_ih(values, nvalues, buffer, 4228, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

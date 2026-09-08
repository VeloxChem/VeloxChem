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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ih_kinetic_energy(double               *values,
                          const size_t          nvalues,
                          const CBasisFunction &bra,
                          const CBasisFunction &ket,
                          const CSimdMatrix    &coordinates,
                          CSimdMatrix          &buffer,
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10310, 9414, 588, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

        return;
    }

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

            simdovl::compute_prim_sh_overlap_0(buffer, 41, 3, 16, 26, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 62, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 63, 3, 7, 62, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_0(buffer, 66, 3, 6, 10, 62, 63, ncols, alpha, beta,
                                             p);

            compute_prim_sf_kinetic_energy_0(buffer, 72, 3, 7, 16, 63, 66, ncols, alpha, beta,
                                             p);

            compute_prim_sg_kinetic_energy_0(buffer, 82, 3, 10, 26, 66, 72, ncols, alpha, beta,
                                             p);

            compute_prim_sh_kinetic_energy_0(buffer, 97, 3, 16, 41, 72, 82, ncols, alpha, beta,
                                             p);

            simdovl::compute_prim_ps_overlap_0(buffer, 118, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 121, 3, 6, 118, ncols, p);

            simdovl::compute_prim_pd_overlap_0(buffer, 130, 0, 3, 7, 10, 121, ncols, p);

            simdovl::compute_prim_pf_overlap_0(buffer, 148, 0, 3, 10, 16, 121, 130, ncols, p);

            simdovl::compute_prim_pg_overlap_0(buffer, 178, 0, 3, 16, 26, 130, 148, ncols, p);

            simdovl::compute_prim_ph_overlap_0(buffer, 223, 0, 3, 26, 41, 148, 178, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 286, 0, 62, 118, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 289, 3, 62, 121, 286, ncols, alpha, beta,
                                             p);

            compute_prim_pd_kinetic_energy_0(buffer, 298, 0, 3, 63, 66, 130, 289, ncols, alpha,
                                             beta, p);

            compute_prim_pf_kinetic_energy_0(buffer, 316, 0, 3, 66, 72, 148, 298, ncols, alpha,
                                             beta, p);

            compute_prim_pg_kinetic_energy_0(buffer, 346, 0, 3, 72, 82, 178, 316, ncols, alpha,
                                             beta, p);

            compute_prim_ph_kinetic_energy_0(buffer, 391, 0, 3, 82, 97, 223, 346, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_ds_overlap_0(buffer, 454, 0, 6, 118, ncols, p);

            simdovl::compute_prim_dp_overlap_0(buffer, 460, 0, 3, 118, 121, 454, ncols, p);

            simdovl::compute_prim_dd_overlap_0(buffer, 478, 0, 3, 121, 130, 454, 460, ncols, p);

            simdovl::compute_prim_df_overlap_0(buffer, 514, 0, 3, 130, 148, 460, 478, ncols, p);

            simdovl::compute_prim_dg_overlap_0(buffer, 574, 0, 3, 148, 178, 478, 514, ncols, p);

            simdovl::compute_prim_dh_overlap_0(buffer, 664, 0, 3, 178, 223, 514, 574, ncols, p);

            compute_prim_ds_kinetic_energy_0(buffer, 790, 0, 6, 62, 286, 454, ncols, alpha, beta,
                                             p);

            compute_prim_dp_kinetic_energy_0(buffer, 796, 0, 3, 286, 289, 460, 790, ncols, alpha,
                                             beta, p);

            compute_prim_dd_kinetic_energy_0(buffer, 814, 0, 3, 289, 298, 454, 478, 790, 796,
                                             ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_0(buffer, 850, 0, 3, 298, 316, 460, 514, 796, 814,
                                             ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_0(buffer, 910, 0, 3, 316, 346, 478, 574, 814, 850,
                                             ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_0(buffer, 1000, 0, 3, 346, 391, 514, 664, 850, 910,
                                             ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_0(buffer, 1126, 0, 118, 454, ncols, p);

            simdovl::compute_prim_fp_overlap_0(buffer, 1136, 0, 3, 454, 460, 1126, ncols, p);

            simdovl::compute_prim_fd_overlap_0(buffer, 1166, 0, 3, 130, 460, 478, 1126, 1136,
                                               ncols, p);

            simdovl::compute_prim_ff_overlap_0(buffer, 1226, 0, 3, 148, 478, 514, 1136, 1166,
                                               ncols, p);

            simdovl::compute_prim_fg_overlap_0(buffer, 1326, 0, 3, 178, 514, 574, 1166, 1226,
                                               ncols, p);

            simdovl::compute_prim_fh_overlap_0(buffer, 1476, 0, 3, 223, 574, 664, 1226, 1326,
                                               ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 1686, 0, 118, 286, 790, 1126, ncols, alpha,
                                             beta, p);

            compute_prim_fp_kinetic_energy_0(buffer, 1696, 0, 3, 790, 796, 1136, 1686, ncols,
                                             alpha, beta, p);

            compute_prim_fd_kinetic_energy_0(buffer, 1726, 0, 3, 130, 298, 796, 814, 1126, 1166,
                                             1686, 1696, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_0(buffer, 1786, 0, 3, 148, 316, 814, 850, 1136, 1226,
                                             1696, 1726, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_0(buffer, 1886, 0, 3, 178, 346, 850, 910, 1166, 1326,
                                             1726, 1786, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_0(buffer, 2036, 0, 3, 223, 391, 910, 1000, 1226, 1476,
                                             1786, 1886, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_0(buffer, 2246, 0, 454, 1126, ncols, p);

            simdovl::compute_prim_gp_overlap_0(buffer, 2261, 0, 3, 1126, 1136, 2246, ncols, p);

            simdovl::compute_prim_gd_overlap_0(buffer, 2306, 0, 3, 478, 1136, 1166, 2246, 2261,
                                               ncols, p);

            simdovl::compute_prim_gf_overlap_0(buffer, 2396, 0, 3, 514, 1166, 1226, 2261, 2306,
                                               ncols, p);

            simdovl::compute_prim_gg_overlap_0(buffer, 2546, 0, 3, 574, 1226, 1326, 2306, 2396,
                                               ncols, p);

            simdovl::compute_prim_gh_overlap_0(buffer, 2771, 0, 3, 664, 1326, 1476, 2396, 2546,
                                               ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 3086, 0, 454, 790, 1686, 2246, ncols, alpha,
                                             beta, p);

            compute_prim_gp_kinetic_energy_0(buffer, 3101, 0, 3, 1686, 1696, 2261, 3086, ncols,
                                             alpha, beta, p);

            compute_prim_gd_kinetic_energy_0(buffer, 3146, 0, 3, 478, 814, 1696, 1726, 2246,
                                             2306, 3086, 3101, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_0(buffer, 3236, 0, 3, 514, 850, 1726, 1786, 2261,
                                             2396, 3101, 3146, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_0(buffer, 3386, 0, 3, 574, 910, 1786, 1886, 2306,
                                             2546, 3146, 3236, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_0(buffer, 3611, 0, 3, 664, 1000, 1886, 2036, 2396,
                                             2771, 3236, 3386, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_0(buffer, 3926, 0, 1126, 2246, ncols, p);

            simdovl::compute_prim_hp_overlap_0(buffer, 3947, 0, 3, 2246, 2261, 3926, ncols, p);

            simdovl::compute_prim_hd_overlap_0(buffer, 4010, 0, 3, 1166, 2261, 2306, 3926, 3947,
                                               ncols, p);

            simdovl::compute_prim_hf_overlap_0(buffer, 4136, 0, 3, 1226, 2306, 2396, 3947, 4010,
                                               ncols, p);

            simdovl::compute_prim_hg_overlap_0(buffer, 4346, 0, 3, 1326, 2396, 2546, 4010, 4136,
                                               ncols, p);

            simdovl::compute_prim_hh_overlap_0(buffer, 4661, 0, 3, 1476, 2546, 2771, 4136, 4346,
                                               ncols, p);

            compute_prim_hs_kinetic_energy_0(buffer, 5102, 0, 1126, 1686, 3086, 3926, ncols,
                                             alpha, beta, p);

            compute_prim_hp_kinetic_energy_0(buffer, 5123, 0, 3, 3086, 3101, 3947, 5102, ncols,
                                             alpha, beta, p);

            compute_prim_hd_kinetic_energy_0(buffer, 5186, 0, 3, 1166, 1726, 3101, 3146, 3926,
                                             4010, 5102, 5123, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_0(buffer, 5312, 0, 3, 1226, 1786, 3146, 3236, 3947,
                                             4136, 5123, 5186, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_0(buffer, 5522, 0, 3, 1326, 1886, 3236, 3386, 4010,
                                             4346, 5186, 5312, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_0(buffer, 5837, 0, 3, 1476, 2036, 3386, 3611, 4136,
                                             4661, 5312, 5522, ncols, alpha, beta, p);

            simdovl::compute_prim_is_overlap_0(buffer, 6278, 0, 2246, 3926, ncols, p);

            simdovl::compute_prim_ip_overlap_0(buffer, 6306, 0, 3, 3926, 3947, 6278, ncols, p);

            simdovl::compute_prim_id_overlap_0(buffer, 6390, 0, 3, 2306, 3947, 4010, 6278, 6306,
                                               ncols, p);

            simdovl::compute_prim_if_overlap_0(buffer, 6558, 0, 3, 2396, 4010, 4136, 6306, 6390,
                                               ncols, p);

            simdovl::compute_prim_ig_overlap_0(buffer, 6838, 0, 3, 2546, 4136, 4346, 6390, 6558,
                                               ncols, p);

            simdovl::compute_prim_ih_overlap_0(buffer, 7258, 0, 3, 2771, 4346, 4661, 6558, 6838,
                                               ncols, p);

            compute_prim_is_kinetic_energy_0(buffer, 7846, 0, 2246, 3086, 5102, 6278, ncols,
                                             alpha, beta, p);

            compute_prim_ip_kinetic_energy_0(buffer, 7874, 0, 3, 5102, 5123, 6306, 7846, ncols,
                                             alpha, beta, p);

            compute_prim_id_kinetic_energy_0(buffer, 7958, 0, 3, 2306, 3146, 5123, 5186, 6278,
                                             6390, 7846, 7874, ncols, alpha, beta, p);

            compute_prim_if_kinetic_energy_0(buffer, 8126, 0, 3, 2396, 3236, 5186, 5312, 6306,
                                             6558, 7874, 7958, ncols, alpha, beta, p);

            compute_prim_ig_kinetic_energy_0(buffer, 8406, 0, 3, 2546, 3386, 5312, 5522, 6390,
                                             6838, 7958, 8126, ncols, alpha, beta, p);

            compute_prim_ih_kinetic_energy_0(buffer, 8826, 0, 3, 2771, 3611, 5522, 5837, 6558,
                                             7258, 8126, 8406, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9414, 8826, 588, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 10002, 9414, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 10002, 11, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

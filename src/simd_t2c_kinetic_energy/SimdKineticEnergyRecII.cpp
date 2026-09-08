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


#include "SimdKineticEnergyRecII.hpp"

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
#include "SimdKineticEnergyVrrRecDI.hpp"
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFG.hpp"
#include "SimdKineticEnergyVrrRecFH.hpp"
#include "SimdKineticEnergyVrrRecFI.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGG.hpp"
#include "SimdKineticEnergyVrrRecGH.hpp"
#include "SimdKineticEnergyVrrRecGI.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHG.hpp"
#include "SimdKineticEnergyVrrRecHH.hpp"
#include "SimdKineticEnergyVrrRecHI.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecID.hpp"
#include "SimdKineticEnergyVrrRecIF.hpp"
#include "SimdKineticEnergyVrrRecIG.hpp"
#include "SimdKineticEnergyVrrRecIH.hpp"
#include "SimdKineticEnergyVrrRecII.hpp"
#include "SimdKineticEnergyVrrRecIP.hpp"
#include "SimdKineticEnergyVrrRecIS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPG.hpp"
#include "SimdKineticEnergyVrrRecPH.hpp"
#include "SimdKineticEnergyVrrRecPI.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSG.hpp"
#include "SimdKineticEnergyVrrRecSH.hpp"
#include "SimdKineticEnergyVrrRecSI.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDG.hpp"
#include "SimdOverlapVrrRecDH.hpp"
#include "SimdOverlapVrrRecDI.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFG.hpp"
#include "SimdOverlapVrrRecFH.hpp"
#include "SimdOverlapVrrRecFI.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGG.hpp"
#include "SimdOverlapVrrRecGH.hpp"
#include "SimdOverlapVrrRecGI.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHG.hpp"
#include "SimdOverlapVrrRecHH.hpp"
#include "SimdOverlapVrrRecHI.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecID.hpp"
#include "SimdOverlapVrrRecIF.hpp"
#include "SimdOverlapVrrRecIG.hpp"
#include "SimdOverlapVrrRecIH.hpp"
#include "SimdOverlapVrrRecII.hpp"
#include "SimdOverlapVrrRecIP.hpp"
#include "SimdOverlapVrrRecIS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPG.hpp"
#include "SimdOverlapVrrRecPH.hpp"
#include "SimdOverlapVrrRecPI.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSH.hpp"
#include "SimdOverlapVrrRecSI.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformII.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ii_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ii_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 6902);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 169 * nvalues, 0.0);

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

            simdovl::compute_prim_sh_overlap_13(buffer, 24, 3, 13, 18, ncols, p);

            simdovl::compute_prim_si_overlap_7(buffer, 31, 3, 18, 24, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 36, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 37, 3, 7, 36, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 40, 3, 6, 10, 36, 37, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_6(buffer, 43, 3, 7, 13, 37, 40, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_8(buffer, 48, 3, 10, 18, 40, 43, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_10(buffer, 54, 3, 13, 24, 43, 48, ncols, alpha, beta, p);

            compute_prim_si_kinetic_energy_4(buffer, 61, 3, 18, 31, 48, 54, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 66, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 69, 3, 6, 66, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 72, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_15(buffer, 75, 0, 3, 10, 13, 72, ncols, p);

            simdovl::compute_prim_pg_overlap_13(buffer, 79, 0, 3, 13, 18, 72, 75, ncols, p);

            simdovl::compute_prim_ph_overlap_9(buffer, 84, 0, 3, 18, 24, 75, 79, ncols, p);

            simdovl::compute_prim_pi_overlap_4(buffer, 90, 0, 3, 24, 31, 79, 84, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 95, 0, 36, 66, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 98, 3, 36, 69, 95, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 101, 0, 37, 40, 72, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_16(buffer, 104, 0, 40, 43, 75, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_13(buffer, 107, 0, 3, 43, 48, 79, 104, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_9(buffer, 112, 0, 3, 48, 54, 84, 107, ncols, alpha, beta, p);

            compute_prim_pi_kinetic_energy_4(buffer, 118, 0, 54, 61, 90, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 123, 0, 6, 66, ncols, p);

            simdovl::compute_prim_dp_overlap_11(buffer, 126, 3, 66, 123, ncols, p);

            simdovl::compute_prim_dd_overlap_17(buffer, 132, 0, 3, 69, 72, 123, 126, ncols, p);

            simdovl::compute_prim_df_overlap_18(buffer, 141, 0, 3, 72, 75, 126, 132, ncols, p);

            simdovl::compute_prim_dg_overlap_14(buffer, 155, 0, 3, 75, 79, 132, 141, ncols, p);

            simdovl::compute_prim_dh_overlap_9(buffer, 174, 0, 3, 79, 84, 141, 155, ncols, p);

            simdovl::compute_prim_di_overlap_4(buffer, 196, 0, 3, 84, 90, 155, 174, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 215, 0, 6, 36, 95, 123, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_12(buffer, 218, 3, 95, 126, 215, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_17(buffer, 224, 0, 3, 98, 101, 123, 132, 215, 218, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_19(buffer, 233, 0, 3, 101, 104, 126, 141, 218, 224, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_14(buffer, 247, 0, 3, 104, 107, 132, 155, 224, 233, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_9(buffer, 266, 0, 3, 107, 112, 141, 174, 233, 247, ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_4(buffer, 288, 0, 3, 112, 118, 155, 196, 247, 266, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_3(buffer, 307, 0, 66, 123, ncols, p);

            simdovl::compute_prim_fp_overlap_14(buffer, 313, 3, 123, 307, ncols, p);

            simdovl::compute_prim_fd_overlap_15(buffer, 321, 0, 3, 72, 126, 132, 307, 313, ncols, p);

            simdovl::compute_prim_ff_overlap_15(buffer, 338, 0, 3, 75, 132, 141, 313, 321, ncols, p);

            simdovl::compute_prim_fg_overlap_11(buffer, 366, 0, 3, 79, 141, 155, 321, 338, ncols, p);

            simdovl::compute_prim_fh_overlap_7(buffer, 407, 0, 3, 84, 155, 174, 338, 366, ncols, p);

            simdovl::compute_prim_fi_overlap_3(buffer, 455, 0, 3, 90, 174, 196, 366, 407, ncols, p);

            compute_prim_fs_kinetic_energy_4(buffer, 490, 0, 66, 95, 215, 307, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_14(buffer, 496, 3, 215, 313, 490, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_17(buffer, 504, 0, 3, 72, 101, 218, 224, 307, 321, 490, 496, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_15(buffer, 521, 0, 3, 75, 104, 224, 233, 313, 338, 496, 504, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_11(buffer, 549, 0, 3, 79, 107, 233, 247, 321, 366, 504, 521, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_7(buffer, 590, 0, 3, 84, 112, 247, 266, 338, 407, 521, 549, ncols, alpha, beta, p);

            compute_prim_fi_kinetic_energy_3(buffer, 638, 0, 3, 90, 118, 266, 288, 366, 455, 549, 590, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_2(buffer, 673, 0, 123, 307, ncols, p);

            simdovl::compute_prim_gp_overlap_11(buffer, 684, 3, 307, 673, ncols, p);

            simdovl::compute_prim_gd_overlap_11(buffer, 695, 0, 3, 132, 313, 321, 673, 684, ncols, p);

            simdovl::compute_prim_gf_overlap_11(buffer, 723, 0, 3, 141, 321, 338, 684, 695, ncols, p);

            simdovl::compute_prim_gg_overlap_8(buffer, 771, 0, 3, 155, 338, 366, 695, 723, ncols, p);

            simdovl::compute_prim_gh_overlap_5(buffer, 846, 0, 3, 174, 366, 407, 723, 771, ncols, p);

            simdovl::compute_prim_gi_overlap_2(buffer, 957, 0, 3, 196, 407, 455, 771, 846, ncols, p);

            compute_prim_gs_kinetic_energy_3(buffer, 1035, 0, 123, 215, 490, 673, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_11(buffer, 1046, 3, 490, 684, 1035, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_14(buffer, 1057, 0, 3, 132, 224, 496, 504, 673, 695, 1035, 1046, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_11(buffer, 1084, 0, 3, 141, 233, 504, 521, 684, 723, 1046, 1057, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_8(buffer, 1131, 0, 3, 155, 247, 521, 549, 695, 771, 1057, 1084, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_5(buffer, 1206, 0, 3, 174, 266, 549, 590, 723, 846, 1084, 1131, ncols, alpha, beta, p);

            compute_prim_gi_kinetic_energy_2(buffer, 1317, 0, 3, 196, 288, 590, 638, 771, 957, 1131, 1206, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_10(buffer, 1395, 0, 307, 673, ncols, p);

            simdovl::compute_prim_hp_overlap_7(buffer, 1409, 0, 3, 673, 684, 1395, ncols, p);

            simdovl::compute_prim_hd_overlap_7(buffer, 1431, 0, 3, 321, 684, 695, 1395, 1409, ncols, p);

            simdovl::compute_prim_hf_overlap_7(buffer, 1476, 0, 3, 338, 695, 723, 1409, 1431, ncols, p);

            simdovl::compute_prim_hg_overlap_5(buffer, 1557, 0, 3, 366, 723, 771, 1431, 1476, ncols, p);

            simdovl::compute_prim_hh_overlap_3(buffer, 1693, 0, 3, 407, 771, 846, 1476, 1557, ncols, p);

            simdovl::compute_prim_hi_overlap_1(buffer, 1955, 0, 3, 455, 846, 957, 1557, 1693, ncols, p);

            compute_prim_hs_kinetic_energy_7(buffer, 2166, 0, 307, 490, 1035, 1395, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_7(buffer, 2180, 0, 3, 1035, 1046, 1409, 2166, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_9(buffer, 2201, 0, 3, 321, 504, 1046, 1057, 1395, 1431, 2166, 2180, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_7(buffer, 2243, 0, 3, 338, 521, 1057, 1084, 1409, 1476, 2180, 2201, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_5(buffer, 2320, 0, 3, 366, 549, 1084, 1131, 1431, 1557, 2201, 2243, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_3(buffer, 2454, 0, 3, 407, 590, 1131, 1206, 1476, 1693, 2243, 2320, ncols, alpha, beta, p);

            compute_prim_hi_kinetic_energy_1(buffer, 2716, 0, 3, 455, 638, 1206, 1317, 1557, 1955, 2320, 2454, ncols, alpha, beta, p);

            simdovl::compute_prim_is_overlap_6(buffer, 2927, 0, 673, 1395, ncols, p);

            simdovl::compute_prim_ip_overlap_3(buffer, 2942, 0, 3, 1395, 1409, 2927, ncols, p);

            simdovl::compute_prim_id_overlap_3(buffer, 2973, 0, 3, 695, 1409, 1431, 2927, 2942, ncols, p);

            simdovl::compute_prim_if_overlap_3(buffer, 3038, 0, 3, 723, 1431, 1476, 2942, 2973, ncols, p);

            simdovl::compute_prim_ig_overlap_2(buffer, 3158, 0, 3, 771, 1476, 1557, 2973, 3038, ncols, p);

            simdovl::compute_prim_ih_overlap_1(buffer, 3366, 0, 3, 846, 1557, 1693, 3038, 3158, ncols, p);

            simdovl::compute_prim_ii_overlap_0(buffer, 3781, 0, 3, 957, 1693, 1955, 3158, 3366, ncols, p);

            compute_prim_is_kinetic_energy_3(buffer, 4565, 0, 673, 1035, 2166, 2927, ncols, alpha, beta, p);

            compute_prim_ip_kinetic_energy_3(buffer, 4579, 3, 2166, 2942, 4565, ncols, alpha, beta, p);

            compute_prim_id_kinetic_energy_4(buffer, 4608, 0, 3, 695, 1057, 2180, 2201, 2927, 2973, 4565, 4579, ncols, alpha, beta, p);

            compute_prim_if_kinetic_energy_3(buffer, 4664, 0, 3, 723, 1084, 2201, 2243, 2942, 3038, 4579, 4608, ncols, alpha, beta, p);

            compute_prim_ig_kinetic_energy_2(buffer, 4766, 0, 3, 771, 1131, 2243, 2320, 2973, 3158, 4608, 4664, ncols, alpha, beta, p);

            compute_prim_ih_kinetic_energy_1(buffer, 4947, 0, 3, 846, 1206, 2320, 2454, 3038, 3366, 4664, 4766, ncols, alpha, beta, p);

            compute_prim_ii_kinetic_energy_0(buffer, 5334, 0, 3, 957, 1317, 2454, 2716, 3158, 3781, 4766, 4947, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 6118, 5334, 784, ncols);
        }
    }

    simdtrf::transform_ii_tri(values, nvalues, buffer, 6118, nmax);

    for (size_t m = 0; m < 169; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

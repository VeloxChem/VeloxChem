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


#include "SimdKineticEnergyRecHH.hpp"

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
#include "SimdTransformHH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hh_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hh_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 3461);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 121 * nvalues, 0.0);

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

            simdovl::compute_prim_dp_overlap_10(buffer, 97, 3, 50, 94, ncols, p);

            simdovl::compute_prim_dd_overlap_15(buffer, 105, 0, 3, 53, 56, 94, 97, ncols, p);

            simdovl::compute_prim_df_overlap_14(buffer, 116, 0, 3, 56, 59, 97, 105, ncols, p);

            simdovl::compute_prim_dg_overlap_10(buffer, 137, 0, 3, 59, 63, 105, 116, ncols, p);

            simdovl::compute_prim_dh_overlap_6(buffer, 157, 0, 3, 63, 68, 116, 137, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 172, 0, 6, 28, 72, 94, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_11(buffer, 175, 3, 72, 97, 172, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_15(buffer, 183, 0, 3, 75, 78, 94, 105, 172, 175, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_14(buffer, 194, 0, 3, 78, 81, 97, 116, 175, 183, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_10(buffer, 215, 0, 3, 81, 85, 105, 137, 183, 194, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_6(buffer, 235, 0, 3, 85, 90, 116, 157, 194, 215, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 250, 0, 50, 94, ncols, p);

            simdovl::compute_prim_fp_overlap_11(buffer, 258, 3, 94, 250, ncols, p);

            simdovl::compute_prim_fd_overlap_11(buffer, 267, 0, 3, 56, 97, 105, 250, 258, ncols, p);

            simdovl::compute_prim_ff_overlap_10(buffer, 287, 0, 3, 59, 105, 116, 258, 267, ncols, p);

            simdovl::compute_prim_fg_overlap_7(buffer, 324, 0, 3, 63, 116, 137, 267, 287, ncols, p);

            simdovl::compute_prim_fh_overlap_4(buffer, 380, 0, 3, 68, 137, 157, 287, 324, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 418, 0, 50, 72, 172, 250, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_11(buffer, 426, 3, 172, 258, 418, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_12(buffer, 435, 0, 3, 56, 78, 175, 183, 250, 267, 418, 426, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_10(buffer, 455, 0, 3, 59, 81, 183, 194, 258, 287, 426, 435, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_7(buffer, 492, 0, 3, 63, 85, 194, 215, 267, 324, 435, 455, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_4(buffer, 548, 0, 3, 68, 90, 215, 235, 287, 380, 455, 492, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_7(buffer, 586, 0, 94, 250, ncols, p);

            simdovl::compute_prim_gp_overlap_7(buffer, 597, 0, 3, 250, 258, 586, ncols, p);

            simdovl::compute_prim_gd_overlap_7(buffer, 614, 0, 3, 105, 258, 267, 586, 597, ncols, p);

            simdovl::compute_prim_gf_overlap_6(buffer, 649, 0, 3, 116, 267, 287, 597, 614, ncols, p);

            simdovl::compute_prim_gg_overlap_4(buffer, 717, 0, 3, 137, 287, 324, 614, 649, ncols, p);

            simdovl::compute_prim_gh_overlap_2(buffer, 847, 0, 3, 157, 324, 380, 649, 717, ncols, p);

            compute_prim_gs_kinetic_energy_4(buffer, 961, 0, 94, 172, 418, 586, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_7(buffer, 972, 0, 3, 418, 426, 597, 961, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_8(buffer, 988, 0, 3, 105, 183, 426, 435, 586, 614, 961, 972, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_6(buffer, 1021, 0, 3, 116, 194, 435, 455, 597, 649, 972, 988, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_4(buffer, 1088, 0, 3, 137, 215, 455, 492, 614, 717, 988, 1021, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_2(buffer, 1218, 0, 3, 157, 235, 492, 548, 649, 847, 1021, 1088, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_7(buffer, 1332, 0, 250, 586, ncols, p);

            simdovl::compute_prim_hp_overlap_3(buffer, 1344, 0, 3, 586, 597, 1332, ncols, p);

            simdovl::compute_prim_hd_overlap_3(buffer, 1369, 0, 3, 267, 597, 614, 1332, 1344, ncols, p);

            simdovl::compute_prim_hf_overlap_2(buffer, 1421, 0, 3, 287, 614, 649, 1344, 1369, ncols, p);

            simdovl::compute_prim_hg_overlap_1(buffer, 1530, 0, 3, 324, 649, 717, 1369, 1421, ncols, p);

            simdovl::compute_prim_hh_overlap_0(buffer, 1754, 0, 3, 380, 717, 847, 1421, 1530, ncols, p);

            compute_prim_hs_kinetic_energy_3(buffer, 2195, 0, 250, 418, 961, 1332, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_3(buffer, 2206, 3, 961, 1344, 2195, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_3(buffer, 2229, 0, 3, 267, 435, 972, 988, 1332, 1369, 2195, 2206, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_2(buffer, 2272, 0, 3, 287, 455, 988, 1021, 1344, 1421, 2206, 2229, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_1(buffer, 2369, 0, 3, 324, 492, 1021, 1088, 1369, 1530, 2229, 2272, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_0(buffer, 2579, 0, 3, 380, 548, 1088, 1218, 1421, 1754, 2272, 2369, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3020, 2579, 441, ncols);
        }
    }

    simdtrf::transform_hh_tri(values, nvalues, buffer, 3020, nmax);

    for (size_t m = 0; m < 121; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

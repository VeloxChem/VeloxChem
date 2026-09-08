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


#include "SimdKineticEnergyRecGH.hpp"

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
#include "SimdTransformGH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_gh_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gh_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 2363);

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

            simdovl::compute_prim_sd_overlap_1(buffer, 10, 3, 6, 7, ncols, p);

            simdovl::compute_prim_sf_overlap_5(buffer, 13, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_10(buffer, 22, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_10(buffer, 30, 3, 13, 22, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 34, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 35, 3, 7, 34, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 38, 3, 6, 10, 34, 35, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 41, 3, 7, 13, 35, 38, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_7(buffer, 50, 3, 10, 22, 38, 41, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_7(buffer, 58, 3, 13, 30, 41, 50, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 62, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 65, 3, 6, 62, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 68, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_13(buffer, 71, 0, 3, 10, 13, 68, ncols, p);

            simdovl::compute_prim_pg_overlap_10(buffer, 78, 0, 3, 13, 22, 68, 71, ncols, p);

            simdovl::compute_prim_ph_overlap_6(buffer, 85, 0, 3, 22, 30, 71, 78, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 89, 0, 34, 62, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 92, 3, 34, 65, 89, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 95, 0, 35, 38, 68, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_13(buffer, 98, 0, 3, 38, 41, 71, 95, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_10(buffer, 105, 0, 3, 41, 50, 78, 98, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_6(buffer, 112, 0, 50, 58, 85, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 116, 0, 6, 62, ncols, p);

            simdovl::compute_prim_dp_overlap_7(buffer, 119, 3, 62, 116, ncols, p);

            simdovl::compute_prim_dd_overlap_11(buffer, 128, 0, 3, 65, 68, 116, 119, ncols, p);

            simdovl::compute_prim_df_overlap_10(buffer, 141, 0, 3, 68, 71, 119, 128, ncols, p);

            simdovl::compute_prim_dg_overlap_7(buffer, 165, 0, 3, 71, 78, 128, 141, ncols, p);

            simdovl::compute_prim_dh_overlap_4(buffer, 204, 0, 3, 78, 85, 141, 165, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 229, 0, 6, 34, 89, 116, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_8(buffer, 232, 3, 89, 119, 229, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_11(buffer, 241, 0, 3, 92, 95, 116, 128, 229, 232, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_10(buffer, 254, 0, 3, 95, 98, 119, 141, 232, 241, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_7(buffer, 278, 0, 3, 98, 105, 128, 165, 241, 254, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_4(buffer, 317, 0, 3, 105, 112, 141, 204, 254, 278, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 342, 0, 62, 116, ncols, p);

            simdovl::compute_prim_fp_overlap_7(buffer, 350, 0, 3, 116, 119, 342, ncols, p);

            simdovl::compute_prim_fd_overlap_7(buffer, 362, 0, 3, 68, 119, 128, 342, 350, ncols, p);

            simdovl::compute_prim_ff_overlap_6(buffer, 387, 0, 3, 71, 128, 141, 350, 362, ncols, p);

            simdovl::compute_prim_fg_overlap_4(buffer, 431, 0, 3, 78, 141, 165, 362, 387, ncols, p);

            simdovl::compute_prim_fh_overlap_2(buffer, 515, 0, 3, 85, 165, 204, 387, 431, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 596, 0, 62, 89, 229, 342, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_7(buffer, 604, 0, 3, 229, 232, 350, 596, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_7(buffer, 615, 0, 3, 68, 95, 232, 241, 342, 362, 596, 604, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_6(buffer, 638, 0, 3, 71, 98, 241, 254, 350, 387, 604, 615, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_4(buffer, 681, 0, 3, 78, 105, 254, 278, 362, 431, 615, 638, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_2(buffer, 765, 0, 3, 85, 112, 278, 317, 387, 515, 638, 681, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_6(buffer, 846, 0, 116, 342, ncols, p);

            simdovl::compute_prim_gp_overlap_3(buffer, 855, 0, 3, 342, 350, 846, ncols, p);

            simdovl::compute_prim_gd_overlap_3(buffer, 874, 0, 3, 128, 350, 362, 846, 855, ncols, p);

            simdovl::compute_prim_gf_overlap_2(buffer, 913, 0, 3, 141, 362, 387, 855, 874, ncols, p);

            simdovl::compute_prim_gg_overlap_1(buffer, 992, 0, 3, 165, 387, 431, 874, 913, ncols, p);

            simdovl::compute_prim_gh_overlap_0(buffer, 1151, 0, 3, 204, 431, 515, 913, 992, ncols, p);

            compute_prim_gs_kinetic_energy_2(buffer, 1466, 0, 116, 229, 596, 846, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_3(buffer, 1474, 3, 596, 855, 1466, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_3(buffer, 1491, 0, 3, 128, 241, 604, 615, 846, 874, 1466, 1474, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_2(buffer, 1521, 0, 3, 141, 254, 615, 638, 855, 913, 1474, 1491, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_1(buffer, 1588, 0, 3, 165, 278, 638, 681, 874, 992, 1491, 1521, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_0(buffer, 1733, 0, 3, 204, 317, 681, 765, 913, 1151, 1521, 1588, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2048, 1733, 315, ncols);
        }
    }

    simdtrf::transform_gh(values, nvalues, buffer, 2048, nmax);

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

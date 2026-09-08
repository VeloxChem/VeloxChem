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


#include "SimdKineticEnergyRecFH.hpp"

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
#include "SimdTransformFH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_fh_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fh_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1489);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 77 * nvalues, 0.0);

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

            simdovl::compute_prim_sg_overlap_6(buffer, 22, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_8(buffer, 34, 3, 13, 22, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 43, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 44, 3, 7, 43, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 47, 3, 6, 10, 43, 44, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 50, 3, 7, 13, 44, 47, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_3(buffer, 59, 3, 10, 22, 47, 50, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_5(buffer, 71, 3, 13, 34, 50, 59, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 80, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 83, 3, 6, 80, ncols, p);

            simdovl::compute_prim_pd_overlap_9(buffer, 86, 0, 3, 7, 10, 83, ncols, p);

            simdovl::compute_prim_pf_overlap_10(buffer, 91, 0, 3, 10, 13, 83, 86, ncols, p);

            simdovl::compute_prim_pg_overlap_7(buffer, 101, 0, 3, 13, 22, 86, 91, ncols, p);

            simdovl::compute_prim_ph_overlap_4(buffer, 117, 0, 3, 22, 34, 91, 101, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 126, 0, 43, 80, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 129, 3, 43, 83, 126, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_7(buffer, 132, 0, 44, 47, 86, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_9(buffer, 136, 0, 3, 47, 50, 91, 132, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_7(buffer, 143, 0, 3, 50, 59, 101, 136, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_4(buffer, 159, 0, 59, 71, 117, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 168, 0, 6, 80, ncols, p);

            simdovl::compute_prim_dp_overlap_3(buffer, 171, 0, 3, 80, 83, 168, ncols, p);

            simdovl::compute_prim_dd_overlap_7(buffer, 181, 0, 3, 83, 86, 168, 171, ncols, p);

            simdovl::compute_prim_df_overlap_6(buffer, 198, 0, 3, 86, 91, 171, 181, ncols, p);

            simdovl::compute_prim_dg_overlap_4(buffer, 226, 0, 3, 91, 101, 181, 198, ncols, p);

            simdovl::compute_prim_dh_overlap_2(buffer, 278, 0, 3, 101, 117, 198, 226, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 334, 0, 6, 43, 126, 168, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_3(buffer, 337, 3, 126, 171, 334, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_7(buffer, 346, 0, 3, 129, 132, 168, 181, 334, 337, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_6(buffer, 361, 0, 3, 132, 136, 171, 198, 337, 346, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_4(buffer, 389, 0, 3, 136, 143, 181, 226, 346, 361, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_2(buffer, 441, 0, 3, 143, 159, 198, 278, 361, 389, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_6(buffer, 497, 0, 80, 168, ncols, p);

            simdovl::compute_prim_fp_overlap_3(buffer, 503, 0, 3, 168, 171, 497, ncols, p);

            simdovl::compute_prim_fd_overlap_3(buffer, 516, 0, 3, 86, 171, 181, 497, 503, ncols, p);

            simdovl::compute_prim_ff_overlap_2(buffer, 542, 0, 3, 91, 181, 198, 503, 516, ncols, p);

            simdovl::compute_prim_fg_overlap_1(buffer, 594, 0, 3, 101, 198, 226, 516, 542, ncols, p);

            simdovl::compute_prim_fh_overlap_0(buffer, 697, 0, 3, 117, 226, 278, 542, 594, ncols, p);

            compute_prim_fs_kinetic_energy_2(buffer, 907, 0, 80, 126, 334, 497, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_3(buffer, 912, 3, 334, 503, 907, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_3(buffer, 923, 0, 3, 86, 132, 337, 346, 497, 516, 907, 912, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_2(buffer, 940, 0, 3, 91, 136, 346, 361, 503, 542, 912, 923, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_1(buffer, 980, 0, 3, 101, 143, 361, 389, 516, 594, 923, 940, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_0(buffer, 1069, 0, 3, 117, 159, 389, 441, 542, 697, 940, 980, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1279, 1069, 210, ncols);
        }
    }

    simdtrf::transform_fh(values, nvalues, buffer, 1279, nmax);

    for (size_t m = 0; m < 77; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

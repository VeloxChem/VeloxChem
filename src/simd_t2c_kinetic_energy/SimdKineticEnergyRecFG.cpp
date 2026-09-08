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


#include "SimdKineticEnergyRecFG.hpp"

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
#include "SimdTransformFG.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_fg_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fg_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 984);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 63 * nvalues, 0.0);

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

            simdovl::compute_prim_sg_overlap_8(buffer, 22, 3, 10, 13, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 28, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 29, 3, 7, 28, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 32, 3, 6, 10, 28, 29, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 35, 3, 7, 13, 29, 32, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_5(buffer, 44, 3, 10, 22, 32, 35, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 50, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 53, 3, 6, 50, ncols, p);

            simdovl::compute_prim_pd_overlap_9(buffer, 56, 0, 3, 7, 10, 53, ncols, p);

            simdovl::compute_prim_pf_overlap_9(buffer, 61, 0, 3, 10, 13, 53, 56, ncols, p);

            simdovl::compute_prim_pg_overlap_6(buffer, 72, 0, 3, 13, 22, 56, 61, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 78, 0, 28, 50, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 81, 3, 28, 53, 78, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_7(buffer, 84, 0, 29, 32, 56, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_8(buffer, 88, 0, 3, 32, 35, 61, 84, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_6(buffer, 99, 0, 35, 44, 72, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 105, 0, 6, 50, ncols, p);

            simdovl::compute_prim_dp_overlap_3(buffer, 108, 0, 3, 50, 53, 105, ncols, p);

            simdovl::compute_prim_dd_overlap_6(buffer, 118, 0, 3, 53, 56, 105, 108, ncols, p);

            simdovl::compute_prim_df_overlap_5(buffer, 136, 0, 3, 56, 61, 108, 118, ncols, p);

            simdovl::compute_prim_dg_overlap_3(buffer, 169, 0, 3, 61, 72, 118, 136, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 209, 0, 6, 28, 78, 105, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_3(buffer, 212, 3, 78, 108, 209, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_6(buffer, 221, 0, 3, 81, 84, 105, 118, 209, 212, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_5(buffer, 239, 0, 3, 84, 88, 108, 136, 212, 221, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_3(buffer, 272, 0, 3, 88, 99, 118, 169, 221, 239, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_6(buffer, 312, 0, 50, 105, ncols, p);

            simdovl::compute_prim_fp_overlap_3(buffer, 318, 0, 3, 105, 108, 312, ncols, p);

            simdovl::compute_prim_fd_overlap_2(buffer, 331, 0, 3, 56, 108, 118, 312, 318, ncols, p);

            simdovl::compute_prim_ff_overlap_1(buffer, 361, 0, 3, 61, 118, 136, 318, 331, ncols, p);

            simdovl::compute_prim_fg_overlap_0(buffer, 429, 0, 3, 72, 136, 169, 331, 361, ncols, p);

            compute_prim_fs_kinetic_energy_2(buffer, 579, 0, 50, 78, 209, 312, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_3(buffer, 584, 3, 209, 318, 579, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_2(buffer, 595, 0, 3, 56, 84, 212, 221, 312, 331, 579, 584, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_1(buffer, 622, 0, 3, 61, 88, 221, 239, 318, 361, 584, 595, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_0(buffer, 684, 0, 3, 72, 99, 239, 272, 331, 429, 595, 622, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 834, 684, 150, ncols);
        }
    }

    simdtrf::transform_fg(values, nvalues, buffer, 834, nmax);

    for (size_t m = 0; m < 63; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

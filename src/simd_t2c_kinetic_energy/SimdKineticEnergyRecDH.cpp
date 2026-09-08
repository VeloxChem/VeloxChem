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


#include "SimdKineticEnergyRecDH.hpp"

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
#include "SimdTransformDH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_dh_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dh_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 841);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 55 * nvalues, 0.0);

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

            simdovl::compute_prim_sh_overlap_6(buffer, 34, 3, 13, 22, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 47, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 48, 3, 7, 47, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 51, 3, 6, 10, 47, 48, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 54, 3, 7, 13, 48, 51, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_3(buffer, 63, 3, 10, 22, 51, 54, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_3(buffer, 75, 3, 13, 34, 54, 63, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 88, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 91, 3, 6, 88, ncols, p);

            simdovl::compute_prim_pd_overlap_6(buffer, 94, 0, 3, 7, 10, 91, ncols, p);

            simdovl::compute_prim_pf_overlap_6(buffer, 101, 0, 3, 10, 13, 91, 94, ncols, p);

            simdovl::compute_prim_pg_overlap_4(buffer, 116, 0, 3, 13, 22, 94, 101, ncols, p);

            simdovl::compute_prim_ph_overlap_2(buffer, 141, 0, 3, 22, 34, 101, 116, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 168, 0, 47, 88, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 171, 3, 47, 91, 168, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_4(buffer, 174, 0, 48, 51, 94, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_5(buffer, 178, 0, 3, 51, 54, 101, 174, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_4(buffer, 189, 0, 3, 54, 63, 116, 178, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_2(buffer, 209, 0, 3, 63, 75, 141, 189, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 236, 0, 6, 88, ncols, p);

            simdovl::compute_prim_dp_overlap_3(buffer, 239, 0, 3, 88, 91, 236, ncols, p);

            simdovl::compute_prim_dd_overlap_3(buffer, 249, 0, 3, 91, 94, 236, 239, ncols, p);

            simdovl::compute_prim_df_overlap_2(buffer, 265, 0, 3, 94, 101, 239, 249, ncols, p);

            simdovl::compute_prim_dg_overlap_1(buffer, 297, 0, 3, 101, 116, 249, 265, ncols, p);

            simdovl::compute_prim_dh_overlap_0(buffer, 356, 0, 3, 116, 141, 265, 297, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 482, 0, 6, 47, 168, 236, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_3(buffer, 485, 3, 168, 239, 482, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_3(buffer, 494, 3, 171, 236, 249, 482, 485, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_2(buffer, 506, 0, 3, 174, 178, 239, 265, 485, 494, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_1(buffer, 534, 0, 3, 178, 189, 249, 297, 494, 506, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_0(buffer, 589, 0, 3, 189, 209, 265, 356, 506, 534, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 715, 589, 126, ncols);
        }
    }

    simdtrf::transform_dh(values, nvalues, buffer, 715, nmax);

    for (size_t m = 0; m < 55; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

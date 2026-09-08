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


#include "SimdKineticEnergyRecDI.hpp"

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
#include "SimdTransformDI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_di_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_di_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1207);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 65 * nvalues, 0.0);

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

            simdovl::compute_prim_sh_overlap_7(buffer, 34, 3, 13, 22, ncols, p);

            simdovl::compute_prim_si_overlap_4(buffer, 52, 3, 22, 34, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 70, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 71, 3, 7, 70, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 74, 3, 6, 10, 70, 71, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 77, 3, 7, 13, 71, 74, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_3(buffer, 86, 3, 10, 22, 74, 77, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_4(buffer, 98, 3, 13, 34, 77, 86, ncols, alpha, beta, p);

            compute_prim_si_kinetic_energy_1(buffer, 116, 3, 22, 52, 86, 98, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 134, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 137, 3, 6, 134, ncols, p);

            simdovl::compute_prim_pd_overlap_6(buffer, 140, 0, 3, 7, 10, 137, ncols, p);

            simdovl::compute_prim_pf_overlap_7(buffer, 147, 0, 3, 10, 13, 137, 140, ncols, p);

            simdovl::compute_prim_pg_overlap_5(buffer, 160, 0, 3, 13, 22, 140, 147, ncols, p);

            simdovl::compute_prim_ph_overlap_3(buffer, 182, 0, 3, 22, 34, 147, 160, ncols, p);

            simdovl::compute_prim_pi_overlap_1(buffer, 218, 0, 3, 34, 52, 160, 182, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 255, 0, 70, 134, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 258, 3, 70, 137, 255, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_4(buffer, 261, 0, 71, 74, 140, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_6(buffer, 265, 0, 74, 77, 147, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_5(buffer, 269, 0, 3, 77, 86, 160, 265, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_3(buffer, 284, 0, 3, 86, 98, 182, 269, ncols, alpha, beta, p);

            compute_prim_pi_kinetic_energy_1(buffer, 314, 0, 3, 98, 116, 218, 284, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 351, 0, 6, 134, ncols, p);

            simdovl::compute_prim_dp_overlap_3(buffer, 354, 0, 3, 134, 137, 351, ncols, p);

            simdovl::compute_prim_dd_overlap_3(buffer, 364, 0, 3, 137, 140, 351, 354, ncols, p);

            simdovl::compute_prim_df_overlap_3(buffer, 380, 0, 3, 140, 147, 354, 364, ncols, p);

            simdovl::compute_prim_dg_overlap_2(buffer, 409, 0, 3, 147, 160, 364, 380, ncols, p);

            simdovl::compute_prim_dh_overlap_1(buffer, 458, 0, 3, 160, 182, 380, 409, ncols, p);

            simdovl::compute_prim_di_overlap_0(buffer, 543, 0, 3, 182, 218, 409, 458, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 711, 0, 6, 70, 255, 351, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_3(buffer, 714, 3, 255, 354, 711, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_3(buffer, 723, 3, 258, 351, 364, 711, 714, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_3(buffer, 735, 3, 261, 354, 380, 714, 723, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_2(buffer, 757, 0, 3, 265, 269, 364, 409, 723, 735, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_1(buffer, 796, 0, 3, 269, 284, 380, 458, 735, 757, ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_0(buffer, 871, 0, 3, 284, 314, 409, 543, 757, 796, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1039, 871, 168, ncols);
        }
    }

    simdtrf::transform_di(values, nvalues, buffer, 1039, nmax);

    for (size_t m = 0; m < 65; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

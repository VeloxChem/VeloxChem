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


#include "SimdKineticEnergyRecIP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecIP.hpp"
#include "SimdKineticEnergyVrrRecIS.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecIP.hpp"
#include "SimdOverlapVrrRecIS.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ip_kinetic_energy(double               *values,
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
            false, std::string("compute_ip_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 840, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 39 * nvalues, 0.0);

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

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 7, 6, ncols, mu);

            simdovl::compute_prim_ps_overlap_0(buffer, 8, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 11, 3, 6, 8, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 20, 0, 7, 8, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 23, 3, 7, 11, 20, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_0(buffer, 32, 0, 6, 8, ncols, p);

            simdovl::compute_prim_dp_overlap_0(buffer, 38, 0, 3, 8, 11, 32, ncols, p);

            compute_prim_ds_kinetic_energy_0(buffer, 56, 0, 6, 7, 20, 32, ncols, alpha, beta,
                                             p);

            compute_prim_dp_kinetic_energy_0(buffer, 62, 0, 3, 20, 23, 38, 56, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_fs_overlap_0(buffer, 80, 0, 8, 32, ncols, p);

            simdovl::compute_prim_fp_overlap_0(buffer, 90, 0, 3, 32, 38, 80, ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 120, 0, 8, 20, 56, 80, ncols, alpha, beta,
                                             p);

            compute_prim_fp_kinetic_energy_0(buffer, 130, 0, 3, 56, 62, 90, 120, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_gs_overlap_0(buffer, 160, 0, 32, 80, ncols, p);

            simdovl::compute_prim_gp_overlap_0(buffer, 175, 0, 3, 80, 90, 160, ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 220, 0, 32, 56, 120, 160, ncols, alpha,
                                             beta, p);

            compute_prim_gp_kinetic_energy_0(buffer, 235, 0, 3, 120, 130, 175, 220, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_hs_overlap_0(buffer, 280, 0, 80, 160, ncols, p);

            simdovl::compute_prim_hp_overlap_0(buffer, 301, 0, 3, 160, 175, 280, ncols, p);

            compute_prim_hs_kinetic_energy_0(buffer, 364, 0, 80, 120, 220, 280, ncols, alpha,
                                             beta, p);

            compute_prim_hp_kinetic_energy_0(buffer, 385, 0, 3, 220, 235, 301, 364, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_is_overlap_0(buffer, 448, 0, 160, 280, ncols, p);

            simdovl::compute_prim_ip_overlap_0(buffer, 476, 0, 3, 280, 301, 448, ncols, p);

            compute_prim_is_kinetic_energy_0(buffer, 560, 0, 160, 220, 364, 448, ncols, alpha,
                                             beta, p);

            compute_prim_ip_kinetic_energy_0(buffer, 588, 0, 3, 364, 385, 476, 560, ncols, alpha,
                                             beta, p);

            simdfunc::contract_primitives(buffer, 672, 588, 84, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 756, 672, 28, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 756, 3, nmax);

    for (size_t m = 0; m < 39; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

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


#include "SimdKineticEnergyRecPF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_pf_kinetic_energy(double               *values,
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
            false, std::string("compute_pf_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 217, 166, 30, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 21 * nvalues, 0.0);

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

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 26, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 27, 3, 7, 26, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_0(buffer, 30, 3, 6, 10, 26, 27, ncols, alpha, beta,
                                             p);

            compute_prim_sf_kinetic_energy_0(buffer, 36, 3, 7, 16, 27, 30, ncols, alpha, beta,
                                             p);

            simdovl::compute_prim_ps_overlap_0(buffer, 46, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 49, 3, 6, 46, ncols, p);

            simdovl::compute_prim_pd_overlap_0(buffer, 58, 0, 3, 7, 10, 49, ncols, p);

            simdovl::compute_prim_pf_overlap_0(buffer, 76, 0, 3, 10, 16, 49, 58, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 106, 0, 26, 46, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 109, 3, 26, 49, 106, ncols, alpha, beta,
                                             p);

            compute_prim_pd_kinetic_energy_0(buffer, 118, 0, 3, 27, 30, 58, 109, ncols, alpha,
                                             beta, p);

            compute_prim_pf_kinetic_energy_0(buffer, 136, 0, 3, 30, 36, 76, 118, ncols, alpha,
                                             beta, p);

            simdfunc::contract_primitives(buffer, 166, 136, 30, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 196, 166, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 196, 7, nmax);

    for (size_t m = 0; m < 21; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin

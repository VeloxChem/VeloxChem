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



#include "SimdTwoCenterElectronRepulsionRecDS.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ds_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &harmonics,
                              const CSimdMatrix    &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionRecDS.compute_ds_electron_repulsion: Basis functions must be of angular momenta two and zero"));
    }

    if (harmonics.number_of_rows() != 5)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionRecDS.compute_ds_electron_repulsion: Harmonics must have 5 rows"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionRecDS.compute_ds_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the squared distances of the atom pairs are carried by the
    // coordinates, so that they are formed once for the whole block instead of
    // once for every combination of basis functions.

    const auto *ab_2 = coordinates.data(6);

    // NOTE: the pairs of primitives are not screened and every one of them
    // reaches every atom pair, as the Coulomb operator decays as the inverse of
    // the interatomic distance, so every row spans the atom pairs of the block.

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 4);

    // set up the arguments of the Boys function of each pair of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto fmu = aexp * b_exps[j] / (aexp + b_exps[j]);

            auto *bargs = boys.data(0, i * nprim_b + j);

#pragma omp simd aligned(bargs, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                bargs[k] = fmu * ab_2[k];
            }
        }
    }

    // NOTE: the Boys function of every pair of primitives is computed by one
    // call, which fills the orders zero to two of every row. The integrals need
    // the order two alone, and the lower orders are formed on the way to it by
    // the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the first row accumulates the factor which the angular components
    // share, and the remaining rows hold the integrals of the angular components.

    auto buffer = CSimdMatrix(6, nvalues);

    auto *prim = buffer.data(0);

    std::fill(prim, prim + nvalues, 0.0);

    constexpr auto fpi = mathconst::pi_value();

    // NOTE: the two-center repulsion of a pair of s type primitives is two pi to
    // the five halves over the exponents times the square root of their sum,
    // times the Boys function of the order the angular momenta ask for.

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    // accumulate the integrals of each pair of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto fexp = aexp + b_exps[j];

            auto ffact = fcoul * anorm * b_norms[j] / (aexp * b_exps[j] * std::sqrt(fexp));

            // NOTE: the ratio is raised to the angular momentum by repeated
            // multiplication, as the angular momentum is small.

            const auto frat = -b_exps[j] / fexp;

            ffact *= frat;

            ffact *= frat;

            const auto *bvals = boys.data(3, i * nprim_b + j);

            // NOTE: the row of the buffer and the row of the Boys function start at
            // a cache line boundary, so the loop is vectorized with aligned loads
            // and stores.

#pragma omp simd aligned(prim, bvals : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                prim[k] += ffact * bvals[k];
            }
        }
    }

    // NOTE: the integral of an angular component is the accumulated factor
    // times the solid harmonic of that component.

    const auto *ph_m2 = harmonics.data(0);
    const auto *ph_m1 = harmonics.data(1);
    const auto *ph_0 = harmonics.data(2);
    const auto *ph_p1 = harmonics.data(3);
    const auto *ph_p2 = harmonics.data(4);

    auto *pb_m2 = buffer.data(1);
    auto *pb_m1 = buffer.data(2);
    auto *pb_0 = buffer.data(3);
    auto *pb_p1 = buffer.data(4);
    auto *pb_p2 = buffer.data(5);

#pragma omp simd aligned(prim, ph_m2, ph_m1, ph_0, ph_p1, ph_p2, pb_m2, pb_m1, pb_0, pb_p1, pb_p2 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto fval = prim[k];

        pb_m2[k] = fval * ph_m2[k];
        pb_m1[k] = fval * ph_m1[k];
        pb_0[k] = fval * ph_0[k];
        pb_p1[k] = fval * ph_p1[k];
        pb_p2[k] = fval * ph_p2[k];
    }

    // NOTE: the values of an angular component are stored as one row of nvalues
    // columns.

    for (size_t m = 0; m < 5; m++)
    {
        const auto *pb = buffer.data(m + 1);

        std::copy(pb, pb + nvalues, values + m * nvalues);
    }
}

}  // namespace simdt2ceri

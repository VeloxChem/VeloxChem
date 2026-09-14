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


#include "SimdNuclearPotentialRecSS.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ss_nuclear_potential(double                    *values,
                             const size_t               nvalues,
                             const CBasisFunction      &bra,
                             const CBasisFunction      &ket,
                             const CSimdMatrix         &coordinates,
                             const std::vector<double> &charges,
                             const std::vector<double> &points,
                             CSimdMatrix               &buffer,
                             const double               threshold) -> void
{
    errors::assertMsgCritical(
        nvalues <= coordinates.number_of_columns(),
        std::string("compute_ss_nuclear_potential: Number of values exceeds number of atom pairs"));

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ss_nuclear_potential: Expecting three coordinates for each charge"));

    if (nvalues == 0) return;

    // NOTE: the values are zeroed before anything writes them, the sums over the
    // primitives and over the charges accumulating into them. The atom pairs no pair
    // of primitives reaches keep the zeros set here.

    std::fill(values, values + 1 * nvalues, 0.0);

    if (charges.empty()) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the integrals
    // divided by their number and by the number of charges, as every integral is a
    // sum over both and the error of a sum is bounded by the number of its terms.

    const auto terms = static_cast<double>(nprims * charges.size());

    const auto dimensions = simdfunc::make_column_dimensions(bra,
                                                             ket,
                                                             nvalues,
                                                             coordinates,
                                                             screenfunc::two_center_nuclear_potential_primitive_bound,
                                                             threshold / terms);

    // NOTE: five rows -- the Gaussian product centre, the argument of the Boys
    // function and its value -- and none of them accumulated into, as every one is
    // written before it is read within a pair of primitives.

    simdfunc::prepare_buffer(buffer, 5, 0, 0, dimensions);

    errors::assertMsgCritical(dimensions.size() == nprims,
                              std::string("Dimensions do not match the pairs of primitives"));

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);

    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    const auto *ab2 = coordinates.data(9);

    // NOTE: the Gaussian product centre and the vector from it to a charge are held
    // in the buffer of the block rather than formed into locals, so that the loops
    // below vectorise over the atom pairs as the kernels of the other operators do.

    auto *p_x = buffer.data(0);
    auto *p_y = buffer.data(1);
    auto *p_z = buffer.data(2);

    // NOTE: the Boys function reads its argument from one row of the buffer and
    // writes its values to the rows which follow, so the argument goes in row three
    // and the value of order zero comes back in row four.

    const size_t argument_row = 3;

    auto *pc_2 = buffer.data(argument_row);

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto fa = a_exps[i] / p;

            const auto fb = b_exps[j] / p;

            // NOTE: two pi over p and not the pi to the three halves of the overlap:
            // the Coulomb operator is integrated over the third coordinate, which
            // replaces one of the Gaussian integrals with the Boys function.

            // NOTE: the integral is positive, the charge of the electron not being
            // carried here. This is the convention of CNuclearPotentialDriver, whose
            // matrix this one has to agree with, and the factor of minus one is
            // applied by the caller -- oneeints.compute_nuclear_potential_integrals
            // does it for the plain driver.

            const auto fnpot = 2.0 * mathconst::pi_value() / p * a_norms[i] * b_norms[j];

#pragma omp simd aligned(a_x, a_y, a_z, b_x, b_y, b_z : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                p_x[k] = fa * a_x[k] + fb * b_x[k];
                p_y[k] = fa * a_y[k] + fb * b_y[k];
                p_z[k] = fa * a_z[k] + fb * b_z[k];
            }

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto c_x = points[3 * ic + 0];
                const auto c_y = points[3 * ic + 1];
                const auto c_z = points[3 * ic + 2];

#pragma omp simd
                for (size_t k = 0; k < ncols; k++)
                {
                    const auto dx = p_x[k] - c_x;
                    const auto dy = p_y[k] - c_y;
                    const auto dz = p_z[k] - c_z;

                    pc_2[k] = p * (dx * dx + dy * dy + dz * dz);
                }

                // NOTE: the Boys function is evaluated for the whole row of atom
                // pairs at once rather than one pair at a time. Only order zero is
                // wanted here; a kernel of higher angular momenta asks for the
                // orders its recursion climbs.

                simdfunc::compute_boys_values(buffer, argument_row, 0, ncols);

                const auto *boys = buffer.data(argument_row + 1);

                const auto factor = fnpot * charges[ic];

#pragma omp simd aligned(ab2 : simd::cache_line_size())
                for (size_t k = 0; k < ncols; k++)
                {
                    values[k] += factor * std::exp(-mu * ab2[k]) * boys[k];
                }
            }
        }
    }
}

}  // namespace simdnpot

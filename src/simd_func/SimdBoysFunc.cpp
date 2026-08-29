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


#include "SimdBoysFunc.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "SimdBoysTable.hpp"

namespace simdfunc {  // simdfunc namespace

auto
compute_boys_function(CSimdVariableMatrix &matrix) -> void
{
    errors::assertMsgCritical(matrix.number_of_blocks() >= 2,
                              std::string("SimdBoysFunc.compute_boys_function: Matrix must have at least two blocks"));

    const auto order = static_cast<int>(matrix.number_of_blocks()) - 2;

    errors::assertMsgCritical(order <= max_boys_order(),
                              std::string("SimdBoysFunc.compute_boys_function: Order of Boys function is out of range"));

    const auto *table = boys_table() + static_cast<size_t>(order) * boys_grid_points() * boys_coefficients();

    const auto *factors = boys_factors();

    const auto fpi = 0.5 * std::sqrt(mathconst::pi_value());

    std::vector<double *> rows(static_cast<size_t>(order) + 1, nullptr);

    for (size_t irow = 0; irow < matrix.number_of_rows_in_block(); irow++)
    {
        const auto ncols = matrix.number_of_columns(irow);

        const auto *args = matrix.data(0, irow);

        for (int j = 0; j <= order; j++)
        {
            rows[static_cast<size_t>(j)] = matrix.data(static_cast<size_t>(j) + 1, irow);
        }

        for (size_t i = 0; i < ncols; i++)
        {
            const auto fa = args[i];

            const auto pnt = (fa > 1.0e5) ? 1000000 : static_cast<int>(10.0 * fa + 0.5);

            if (pnt < boys_grid_points())
            {
                // NOTE: the highest order is taken from the Taylor expansion on
                // the grid of arguments, and the lower orders follow from the
                // downward recursion, which is stable for small arguments.

                const auto *coefs = table + static_cast<size_t>(pnt) * boys_coefficients();

                const auto w = fa - 0.1 * pnt;

                const auto w2 = w * w;

                const auto w4 = w2 * w2;

                rows[static_cast<size_t>(order)][i] = coefs[0] + coefs[1] * w + coefs[2] * w2 + coefs[3] * w2 * w + coefs[4] * w4 +
                                                      coefs[5] * w4 * w + coefs[6] * w4 * w2 + coefs[7] * w4 * w2 * w;

                const auto f2a = fa + fa;

                const auto fx = std::exp(-fa);

                for (int j = order - 1; j >= 0; j--)
                {
                    rows[static_cast<size_t>(j)][i] =
                        factors[j] * (f2a * rows[static_cast<size_t>(j) + 1][i] + fx);
                }
            }
            else
            {
                // NOTE: the lowest order is taken from the asymptotic expansion,
                // and the higher orders follow from the upward recursion, which
                // is stable for large arguments.

                const auto fia = 1.0 / fa;

                const auto fia2 = fia * fia;

                const auto f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                // NOTE: the exponential term is kept for every argument, as it
                // underflows to zero on its own beyond an argument of about seven
                // hundred. Dropping it earlier costs accuracy for the highest
                // orders, whose values approach it well before that.

                const auto fx = std::exp(-fa);

                rows[0][i] = fpi * std::sqrt(fia) - f * fx;

                const auto rterm = 0.5 * fia * fx;

                // NOTE: the factor of the recursion is recomputed at every order
                // rather than accumulated, as accumulating it compounds the
                // rounding of the additions into every value which follows and
                // costs an order of magnitude of accuracy at the highest orders.

                for (int j = 1; j <= order; j++)
                {
                    const auto pf = (static_cast<double>(j) - 0.5) * fia;

                    rows[static_cast<size_t>(j)][i] = pf * rows[static_cast<size_t>(j) - 1][i] - rterm;
                }
            }
        }
    }
}

}  // namespace simdfunc

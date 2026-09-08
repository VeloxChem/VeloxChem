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

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysTable.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Writes the argument of Boys function, the scaled squared distance of the
/// atom pair.
static auto
_make_arguments(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target, const size_t ncols, const double mu)
    -> void
{
    auto *args = buffer.data(target);

    const auto *ab_2 = coordinates.data(9);

#pragma omp simd aligned(args, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        args[k] = mu * ab_2[k];
    }
}

/// @brief Computes the values of Boys function of order zero to order for one
/// argument, into the leading order + 1 entries of values.
/// @note The highest order is taken from the Taylor expansion on the grid and the
/// lower orders follow from the downward recursion, which is stable for small
/// arguments; above the grid the lowest order is taken from the asymptotic expansion
/// and the higher orders follow from the upward recursion, which is stable there.
static inline auto
_boys_ladder(double *values, const double fa, const int order, const double *table, const double *factors, const double fpi)
    -> void
{
    const auto pnt = (fa > 1.0e5) ? 1000000 : static_cast<int>(10.0 * fa + 0.5);

    if (pnt < boys_grid_points())
    {
        const auto *coefs = table + static_cast<size_t>(pnt) * boys_coefficients();

        const auto w = fa - 0.1 * pnt;

        const auto w2 = w * w;

        const auto w4 = w2 * w2;

        values[order] = coefs[0] + coefs[1] * w + coefs[2] * w2 + coefs[3] * w2 * w + coefs[4] * w4 + coefs[5] * w4 * w +
                        coefs[6] * w4 * w2 + coefs[7] * w4 * w2 * w;

        const auto f2a = fa + fa;

        const auto fx = std::exp(-fa);

        for (int j = order - 1; j >= 0; j--)
        {
            values[j] = factors[j] * (f2a * values[j + 1] + fx);
        }
    }
    else
    {
        const auto fia = 1.0 / fa;

        const auto fia2 = fia * fia;

        const auto f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

        // NOTE: the exponential term is kept for every argument, as it underflows to
        // zero on its own beyond an argument of about seven hundred. Dropping it
        // earlier costs accuracy for the highest orders, whose values approach it
        // well before that.

        const auto fx = std::exp(-fa);

        values[0] = fpi * std::sqrt(fia) - f * fx;

        const auto rterm = 0.5 * fia * fx;

        // NOTE: the factor of the recursion is recomputed at every order rather than
        // accumulated, as accumulating it compounds the rounding of the additions
        // into every value which follows and costs an order of magnitude of accuracy
        // at the highest orders.

        for (int j = 1; j <= order; j++)
        {
            const auto pf = (static_cast<double>(j) - 0.5) * fia;

            values[j] = pf * values[j - 1] - rterm;
        }
    }
}

auto
compute_full_boys_function(CSimdMatrix       &buffer,
                           const CSimdMatrix &coordinates,
                           const size_t       target,
                           const size_t       order,
                           const size_t       ncols,
                           const double       fj,
                           const double       mu) -> void
{
    errors::assertMsgCritical(static_cast<int>(order) <= max_boys_order(),
                              std::string("SimdBoysFunc.compute_full_boys_function: Order of Boys function is out of range"));

    errors::assertMsgCritical(target + order + 1 < buffer.number_of_rows(),
                              std::string("SimdBoysFunc.compute_full_boys_function: Buffer has too few rows"));

    _make_arguments(buffer, coordinates, target, ncols, mu);

    const auto iorder = static_cast<int>(order);

    const auto *table = boys_table() + static_cast<size_t>(iorder) * boys_grid_points() * boys_coefficients();

    const auto *factors = boys_factors();

    const auto fpi = 0.5 * std::sqrt(mathconst::pi_value());

    const auto *args = buffer.data(target);

    // NOTE: every order is wanted here, so the values are formed straight into the
    // rows of the buffer and the ladder is the rows themselves.

    std::vector<double *> rows(order + 1, nullptr);

    for (size_t j = 0; j <= order; j++)
    {
        rows[j] = buffer.data(target + 1 + j);
    }

    std::vector<double> ladder(order + 1, 0.0);

    for (size_t i = 0; i < ncols; i++)
    {
        _boys_ladder(ladder.data(), args[i], iorder, table, factors, fpi);

        for (size_t j = 0; j <= order; j++)
        {
            rows[j][i] = fj * ladder[j];
        }
    }
}

auto
compute_boys_function(CSimdMatrix                        &buffer,
                      const CSimdMatrix                  &coordinates,
                      const size_t                        target,
                      const std::initializer_list<size_t> orders,
                      const size_t                        ncols,
                      const double                        fj,
                      const double                        mu) -> void
{
    errors::assertMsgCritical(orders.size() > 0, std::string("SimdBoysFunc.compute_boys_function: No order was requested"));

    const auto top = *std::ranges::max_element(orders);

    errors::assertMsgCritical(static_cast<int>(top) <= max_boys_order(),
                              std::string("SimdBoysFunc.compute_boys_function: Order of Boys function is out of range"));

    errors::assertMsgCritical(target + orders.size() < buffer.number_of_rows(),
                              std::string("SimdBoysFunc.compute_boys_function: Buffer has too few rows"));

    _make_arguments(buffer, coordinates, target, ncols, mu);

    const auto itop = static_cast<int>(top);

    const auto *table = boys_table() + static_cast<size_t>(itop) * boys_grid_points() * boys_coefficients();

    const auto *factors = boys_factors();

    const auto fpi = 0.5 * std::sqrt(mathconst::pi_value());

    const auto *args = buffer.data(target);

    // NOTE: the ladder is built to the highest requested order in a scratch of its
    // own and only the requested orders are written out, as the rows of the buffer
    // hold those alone and the orders below them have nowhere to go.

    std::vector<double> ladder(top + 1, 0.0);

    std::vector<double *> rows;

    rows.reserve(orders.size());

    for (size_t j = 0; j < orders.size(); j++)
    {
        rows.push_back(buffer.data(target + 1 + j));
    }

    for (size_t i = 0; i < ncols; i++)
    {
        _boys_ladder(ladder.data(), args[i], itop, table, factors, fpi);

        size_t irow = 0;

        for (const auto j : orders)
        {
            rows[irow][i] = fj * ladder[j];

            irow++;
        }
    }
}

}  // namespace simdfunc

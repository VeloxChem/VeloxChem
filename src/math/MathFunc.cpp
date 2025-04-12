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

#include "MathFunc.hpp"

#include <cmath>
#include <cstdlib>

#include "MathConst.hpp"

namespace mathfunc {  // mathfunc namespace

auto
countSignificantElements(const std::vector<int64_t>& mask) -> int64_t
{
    int64_t nelems = 0;

    for (auto mvalue : mask)
    {
        if (mvalue == 1) nelems++;
    }

    return nelems;
}

auto
zero(std::vector<double>& values) -> void
{
    const auto ndim = values.size();

    auto ptr_values = values.data();

#pragma omp simd
    for (size_t i = 0; i < ndim; i++)
    {
        ptr_values[i] = 0.0;
    }
}

auto
quadChebyshevOfKindTwo(double* coordinates, double* weights, const int64_t nPoints) -> void
{
    // prefactor

    auto fstep = mathconst::getPiValue() / (static_cast<double>(nPoints) + 1.0);

    // loop over grid points

    for (int64_t i = 1; i < nPoints + 1; i++)
    {
        auto farg = static_cast<double>(i) * fstep;

        coordinates[i - 1] = std::cos(farg);

        auto warg = std::sin(farg);

        weights[i - 1] = fstep * warg * warg;
    }
}

auto
batch_sizes(const int64_t nElements, const int64_t nodes) -> std::vector<int64_t>
{
    int64_t ave = nElements / nodes;

    int64_t rem = nElements % nodes;

    std::vector<int64_t> counts;

    for (int64_t p = 0; p < nodes; p++)
    {
        counts.push_back((p < rem) ? (ave + 1) : ave);
    }

    return counts;
}

auto
batch_offsets(const int64_t nElements, const int64_t nodes) -> std::vector<int64_t>
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    std::vector<int64_t> displs;

    int64_t index = 0;

    for (int64_t p = 0; p < nodes; p++)
    {
        displs.push_back(index);

        index += counts[p];
    }

    return displs;
}

auto
batch_size(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    return counts[rank];
}

auto
batch_offset(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t
{
    auto displs = mathfunc::batch_offsets(nElements, nodes);

    return displs[rank];
}

}  // namespace mathfunc

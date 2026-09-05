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


#ifndef SimdStorage_hpp
#define SimdStorage_hpp

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Stores the integrals of the angular components accumulated in a buffer
/// into the values of a combination of basis functions.
/// @param values The values of the combination of basis functions in the values
/// block of the sparsity pattern.
/// @param nvalues The number of values of one angular component, i.e. the number
/// of atom pairs surviving the screening of the combination of basis functions.
/// @param buffer The buffer holding the integrals of the angular components, one
/// row each, spanning the atom pairs reached by the pairs of primitives.
/// @param first The index of the row of the buffer holding the first angular
/// component, which follows the rows the buffer uses for its accumulators.
/// @param ncomps The number of angular components of the combination of basis
/// functions.
/// @note The values of an angular component span nvalues columns and the
/// components follow one another, with the components of the bra side running
/// slowest. The buffer is shorter than that whenever the screening of the pairs
/// of primitives reaches fewer atom pairs than the screening of the basis
/// functions did, so the atom pairs beyond its reach are set to zero here.
/// @note An empty buffer means no pair of primitives reached any atom pair, and
/// every value of every component is then zero.
inline auto
store_components(double *values, const size_t nvalues, const CSimdMatrix &buffer, const size_t first, const size_t ncomps) -> void
{
    const auto nmax = buffer.number_of_columns();

    if (nmax == 0)
    {
        std::fill(values, values + ncomps * nvalues, 0.0);

        return;
    }

    errors::assertMsgCritical(first + ncomps <= buffer.number_of_rows(),
                              std::string("SimdStorage.store_components: Angular components exceed the rows of the buffer"));

    errors::assertMsgCritical(nmax <= nvalues,
                              std::string("SimdStorage.store_components: Buffer is wider than the values of a component"));

    for (size_t m = 0; m < ncomps; m++)
    {
        const auto *comp = buffer.data(first + m);

        auto *out = values + m * nvalues;

        std::copy(comp, comp + nmax, out);

        std::fill(out + nmax, out + nvalues, 0.0);
    }
}

}  // namespace simdfunc

#endif /* SimdStorage_hpp */

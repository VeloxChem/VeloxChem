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


#include "SimdCoordinates.hpp"

#include <algorithm>
#include <ranges>

#include "ErrorHandler.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_coordinates(const CAtomBasisPairSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix
{
    const auto coords = molecule.coordinates();

    const auto &bra_atoms = sparsity.bra_atoms();

    const auto &ket_atoms = sparsity.ket_atoms();

    errors::assertMsgCritical(
        std::ranges::all_of(bra_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }) &&
            std::ranges::all_of(ket_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }),
        std::string("SimdCoordinates.make_coordinates: Atomic index out of range of molecule"));

    const auto npairs = sparsity.number_of_pairs();

    auto matrix = CSimdMatrix(6, npairs);

    // NOTE: the coordinates of an axis are stored contiguously over the atom
    // pairs, so that a recursion loads them with aligned SIMD instructions.

    std::ranges::for_each(std::views::iota(size_t{0}, npairs), [&](const auto i) {
        const auto r_a = coords[bra_atoms[i]].coordinates();

        const auto r_b = coords[ket_atoms[i]].coordinates();

        matrix.data(0)[i] = r_a[0];

        matrix.data(1)[i] = r_a[1];

        matrix.data(2)[i] = r_a[2];

        matrix.data(3)[i] = r_b[0];

        matrix.data(4)[i] = r_b[1];

        matrix.data(5)[i] = r_b[2];
    });

    return matrix;
}

}  // namespace simdfunc

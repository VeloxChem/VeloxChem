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

#include <vector>

#include "ErrorHandler.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the coordinates of the atom pairs given by their atomic indices.
/// @param bra_atoms The atomic indices of the atom pairs on bra side.
/// @param ket_atoms The atomic indices of the atom pairs on ket side.
/// @param npairs The number of atom pairs.
/// @param molecule The molecule to take the atomic coordinates from.
/// @return The matrix of coordinates with seven rows and npairs columns.
static auto
_make_coordinates(const std::vector<int> &bra_atoms, const std::vector<int> &ket_atoms, const size_t npairs,
                  const CMolecule &molecule) -> CSimdMatrix
{
    const auto &coords = molecule.coordinates();

    const auto natoms = static_cast<int>(coords.size());

    const auto in_range = [&](const std::vector<int> &atoms) {
        for (size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i] >= natoms) return false;
        }

        return true;
    };

    errors::assertMsgCritical(in_range(bra_atoms) && in_range(ket_atoms),
                              std::string("SimdCoordinates.make_coordinates: Atomic index out of range of molecule"));

    auto matrix = CSimdMatrix(7, npairs);

    // NOTE: the coordinates of an axis are stored contiguously over the atom
    // pairs, so that a recursion loads them with aligned SIMD instructions. The
    // pointers to the rows are taken once, as the accessor of a row is bounds
    // checked and would otherwise be called for every coordinate of every pair.

    auto *a_x = matrix.data(0);
    auto *a_y = matrix.data(1);
    auto *a_z = matrix.data(2);
    auto *b_x = matrix.data(3);
    auto *b_y = matrix.data(4);
    auto *b_z = matrix.data(5);

    auto *ab_2 = matrix.data(6);

    const auto *rxyz = coords.data();

    const auto *bra = bra_atoms.data();

    const auto *ket = ket_atoms.data();

    for (size_t i = 0; i < npairs; i++)
    {
        const auto r_a = rxyz[bra[i]].coordinates();

        const auto r_b = rxyz[ket[i]].coordinates();

        a_x[i] = r_a[0];

        a_y[i] = r_a[1];

        a_z[i] = r_a[2];

        b_x[i] = r_b[0];

        b_y[i] = r_b[1];

        b_z[i] = r_b[2];

        // NOTE: the squared distance of an atom pair is stored alongside the
        // coordinates, as the recursions and the screening would otherwise
        // recompute it for every combination of basis functions of the block.

        const auto rx = r_a[0] - r_b[0];

        const auto ry = r_a[1] - r_b[1];

        const auto rz = r_a[2] - r_b[2];

        ab_2[i] = rx * rx + ry * ry + rz * rz;
    }

    return matrix;
}

auto
make_coordinates(const CAtomBasisPairSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix
{
    return _make_coordinates(sparsity.bra_atoms(), sparsity.ket_atoms(), sparsity.number_of_pairs(), molecule);
}

auto
make_coordinates(const CAtomBasisPairGroup &group, const CMolecule &molecule) -> CSimdMatrix
{
    return _make_coordinates(group.bra_atoms(), group.ket_atoms(), group.number_of_pairs(), molecule);
}

auto
make_coordinates(const CAtomBasisTripleSparsity &sparsity, const CMolecule &molecule) -> CSimdMatrix
{
    return _make_coordinates(sparsity.a_atoms(), sparsity.b_atoms(), sparsity.number_of_pairs(), molecule);
}

}  // namespace simdfunc

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


#include "SimdOverlapDriver.hpp"

#include <algorithm>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "ScreeningFunc.hpp"
#include "SimdOverlapFunc.hpp"

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix
{
    // NOTE: the overlap operator is spherically symmetric about an atom, so the
    // diagonal atom pair blocks hold one value for each pair of basis functions
    // with the same angular momentum.

    auto matrix = CSparseMatrix(molecule, basis, screener::overlap, _threshold, mat_t::symmetric, diagstor::scalar);

    matrix.allocate();

    matrix.zero();

    _compute_pair_blocks(matrix, molecule, basis, basis);

    _compute_diagonal_blocks(matrix, molecule, basis, basis);

    return matrix;
}

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
    -> CSparseMatrix
{
    auto matrix = CSparseMatrix(molecule, bra_basis, ket_basis, screener::overlap, _threshold, mat_t::general, diagstor::scalar);

    matrix.allocate();

    matrix.zero();

    _compute_pair_blocks(matrix, molecule, bra_basis, ket_basis);

    _compute_diagonal_blocks(matrix, molecule, bra_basis, ket_basis);

    return matrix;
}

auto
CSimdOverlapDriver::_compute_pair_blocks(CSparseMatrix         &matrix,
                                         const CMolecule       &molecule,
                                         const CMolecularBasis &bra_basis,
                                         const CMolecularBasis &ket_basis) const -> void
{
    // TODO: for each off-diagonal block of the matrix
    //         create the coordinates of its atom pairs
    //         for each combination of basis functions on bra and ket sides
    //           determine the number of surviving atom pairs of each pair of
    //             primitives
    //           set up the primitive factors of the pairs of primitives
    //           compute the primitive integrals with the overlap recursion
    //           contract the primitive integrals
    //           transform the contracted integrals to the spherical basis
    //           add the integrals to the values of the block

    errors::assertMsgCritical(false, std::string("SimdOverlapDriver: Computation of off-diagonal blocks is not implemented"));
}

/// @brief Gets the angular momentum and the index within it of each basis
/// function of an atom basis.
/// @param basis The atom basis to index the basis functions of.
/// @return The vector of angular momenta and indices within them.
static auto
_index_functions(const CAtomBasis &basis) -> std::vector<std::pair<int, size_t>>
{
    std::vector<std::pair<int, size_t>> indices;

    std::vector<size_t> counts(static_cast<size_t>(basis.max_angular_momentum() + 1), 0);

    std::ranges::for_each(basis.functions(), [&](const auto &bfn) {
        const auto lval = bfn.get_angular_momentum();

        indices.push_back({lval, counts[lval]});

        counts[lval]++;
    });

    return indices;
}

auto
CSimdOverlapDriver::_compute_diagonal_blocks(CSparseMatrix         &matrix,
                                             const CMolecule       &molecule,
                                             const CMolecularBasis &bra_basis,
                                             const CMolecularBasis &ket_basis) const -> void
{
    // NOTE: the overlap of two basis functions on the same atom does not depend
    // on the position of the atom, so a single value is stored for each pair of
    // basis functions with the same angular momentum and the molecule is not
    // needed here.

    std::ranges::for_each(std::views::iota(size_t{0}, matrix.number_of_diagonal_blocks()), [&](const auto iblk) {
        const auto &block = matrix.diagonal_block(iblk);

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        const auto a_indices = _index_functions(a_basis);

        const auto b_indices = _index_functions(b_basis);

        auto *values = matrix.diagonal_values(iblk);

        std::ranges::for_each(std::views::iota(size_t{0}, a_indices.size()), [&](const auto i) {
            std::ranges::for_each(std::views::iota(size_t{0}, b_indices.size()), [&](const auto j) {
                // NOTE: only the stored combinations are computed, as the values
                // of the reverse order share their storage.

                if (block.is_triangular() && (i > j)) return;

                const auto [la, ia] = a_indices[i];

                const auto [lb, jb] = b_indices[j];

                if (block.number_of_elements(la, ia, lb, jb) == 0) return;

                values[block.element_offset(la, ia, lb, jb)] =
                    simdovl::one_center_overlap(a_basis.functions()[i], b_basis.functions()[j]);
            });
        });
    });
}

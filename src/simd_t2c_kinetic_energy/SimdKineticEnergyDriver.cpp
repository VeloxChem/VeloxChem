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


#include "SimdKineticEnergyDriver.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdHarmonics.hpp"
#include "SimdKineticEnergyFunc.hpp"

auto
CSimdKineticEnergyDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix
{
    // NOTE: the kernels of the kinetic energy are not written yet, so the driver
    // refuses here rather than where they would be called. The kernels are called
    // from inside a parallel region, which an exception must not leave, and a
    // critical error there would end the interpreter of a caller who merely asked
    // for a matrix which cannot be computed yet.

    throw std::runtime_error("SimdKineticEnergyDriver.compute: Kernels of the kinetic energy are not implemented");

    // NOTE: the kinetic energy operator is spherically symmetric about an atom,
    // so the diagonal atom pair blocks hold one value for each pair of basis
    // functions with the same angular momentum, as the overlap does.

    auto matrix = CSparseMatrix(molecule, basis, screener::kinetic_energy, _threshold, mat_t::symmetric, diagstor::scalar);

    // NOTE: the values blocks are not set to zero after they are allocated, as
    // every value of every block is written below. A combination of basis
    // functions reaching no atom pair holds no values, and a kernel writes the
    // integrals of the atom pairs it reaches and zeros of the remaining ones, so
    // no value is left with the undefined content of the allocation.

    matrix.allocate();

    _compute_pair_blocks(matrix, molecule, basis, basis);

    _compute_diagonal_blocks(matrix, molecule, basis, basis);

    return matrix;
}

auto
CSimdKineticEnergyDriver::compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
    -> CSparseMatrix
{
    // NOTE: the kernels of the kinetic energy are not written yet, so the driver
    // refuses here rather than where they would be called. The kernels are called
    // from inside a parallel region, which an exception must not leave, and a
    // critical error there would end the interpreter of a caller who merely asked
    // for a matrix which cannot be computed yet.

    throw std::runtime_error("SimdKineticEnergyDriver.compute: Kernels of the kinetic energy are not implemented");

    auto matrix = CSparseMatrix(molecule, bra_basis, ket_basis, screener::kinetic_energy, _threshold, mat_t::general, diagstor::scalar);

    // NOTE: the values blocks are not set to zero after they are allocated, as
    // every value of every block is written below. A combination of basis
    // functions reaching no atom pair holds no values, and a kernel writes the
    // integrals of the atom pairs it reaches and zeros of the remaining ones, so
    // no value is left with the undefined content of the allocation.

    matrix.allocate();

    _compute_pair_blocks(matrix, molecule, bra_basis, ket_basis);

    _compute_diagonal_blocks(matrix, molecule, bra_basis, ket_basis);

    return matrix;
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

    const auto &functions = basis.functions();

    for (size_t i = 0; i < functions.size(); i++)
    {
        const auto lval = functions[i].get_angular_momentum();

        indices.push_back({lval, counts[lval]});

        counts[lval]++;
    }

    return indices;
}

auto
CSimdKineticEnergyDriver::_compute_pair_blocks(CSparseMatrix         &matrix,
                                         const CMolecule       &molecule,
                                         const CMolecularBasis &bra_basis,
                                         const CMolecularBasis &ket_basis) const -> void
{
    // NOTE: the blocks are independent, as each of them forms its own coordinates
    // and solid harmonics and writes the values of its own combinations of basis
    // functions, which no other block addresses. Dynamic scheduling is used as
    // the blocks hold a comparable number of atom pairs but differ in the number
    // of the combinations of basis functions and in the cost of their kernels.

    const auto nblocks = static_cast<int>(matrix.number_of_pair_blocks());

    // NOTE: the blocks are visited from the most costly to the least, so that a
    // costly block is taken while there is still work to fill the other threads
    // with. The threads draw two or three blocks each, so a costly block drawn
    // last is finished alone and sets the time of the whole loop.

    // NOTE: the cost of a block is the number of its atom pairs times the cost of
    // its combinations of basis functions, each weighted by the sum of the
    // angular momenta it carries, as the recursions of the integrals grow with
    // that sum. The weight is an estimate and orders the blocks, it is not used
    // for anything else.

    std::vector<size_t> order(static_cast<size_t>(nblocks));

    std::vector<double> costs(static_cast<size_t>(nblocks), 0.0);

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = matrix.pair_block(static_cast<size_t>(iblk));

        const auto a_indices = _index_functions(bra_basis.basis_set(block.bra_index()));

        const auto b_indices = _index_functions(ket_basis.basis_set(block.ket_index()));

        double weight = 0.0;

        for (const auto &[la, ia] : a_indices)
        {
            for (const auto &[lb, jb] : b_indices)
            {
                const auto lsum = static_cast<double>(la + lb + 1);

                weight += lsum * lsum;
            }
        }

        order[static_cast<size_t>(iblk)] = static_cast<size_t>(iblk);

        costs[static_cast<size_t>(iblk)] = static_cast<double>(block.number_of_pairs()) * weight;
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto jblk = order[static_cast<size_t>(iblk)];

        const auto &block = matrix.pair_block(jblk);

        if (block.number_of_pairs() == 0) continue;

        // NOTE: the coordinates of the atom pairs are created once for the whole
        // block, as all combinations of basis functions of the block share them.

        const auto coordinates = simdfunc::make_coordinates(block, molecule);

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        // NOTE: the solid harmonics of the vectors between the atoms are created
        // once for the whole block, as all combinations of basis functions of the
        // block share them. They reach the sum of the angular momenta of the atom
        // bases of the block, as that is the highest angular momentum the
        // recursions of the integrals reach, and are empty when both atom bases
        // carry S type functions only.

        // NOTE: the recursions of the kinetic energy reach two angular momenta
        // above the sum of the angular momenta of the atom bases, as the operator
        // differentiates twice, so the harmonics are created to that sum plus two.

        const auto lmax = a_basis.max_angular_momentum() + b_basis.max_angular_momentum() + 2;

        const auto harmonics = simdfunc::make_solid_harmonics(coordinates, lmax);

        const auto a_indices = _index_functions(a_basis);

        const auto b_indices = _index_functions(b_basis);

        // NOTE: the atom bases of an off-diagonal block sit on different atoms,
        // so all combinations of basis functions are computed and none of them
        // shares its storage with the reverse order.

        // NOTE: the combinations of basis functions are independent, as each of
        // them writes its own values and reads the coordinates and the harmonics
        // of the block without changing them.

        for (size_t i = 0; i < a_indices.size(); i++)
        {
            for (size_t j = 0; j < b_indices.size(); j++)
            {
                const auto [la, ia] = a_indices[i];

                const auto [lb, jb] = b_indices[j];

                const auto nvalues = block.number_of_pairs(la, ia, lb, jb);

                if (nvalues == 0) continue;

                simdkin::compute_kinetic_energy(matrix.pair_values(jblk, la, ia, lb, jb),
                                         nvalues,
                                         a_basis.functions()[i],
                                         b_basis.functions()[j],
                                         harmonics,
                                         coordinates,
                                         _threshold);
            }
        }

    }
}

auto
CSimdKineticEnergyDriver::_compute_diagonal_blocks(CSparseMatrix         &matrix,
                                             const CMolecule       &molecule,
                                             const CMolecularBasis &bra_basis,
                                             const CMolecularBasis &ket_basis) const -> void
{
    // NOTE: the kinetic energy of two basis functions on the same atom does not
    // depend on the position of the atom, so a single value is stored for each
    // pair of basis functions with the same angular momentum and the molecule is
    // not needed here.

    for (size_t iblk = 0; iblk < matrix.number_of_diagonal_blocks(); iblk++)
    {
        const auto &block = matrix.diagonal_block(iblk);

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        const auto a_indices = _index_functions(a_basis);

        const auto b_indices = _index_functions(b_basis);

        auto *values = matrix.diagonal_values(iblk);

        for (size_t i = 0; i < a_indices.size(); i++)
        {
            for (size_t j = 0; j < b_indices.size(); j++)
            {
                // NOTE: only the stored combinations are computed, as the values
                // of the reverse order share their storage.

                if (block.is_triangular() && (i > j)) continue;

                const auto [la, ia] = a_indices[i];

                const auto [lb, jb] = b_indices[j];

                if (block.number_of_elements(la, ia, lb, jb) == 0) continue;

                values[block.element_offset(la, ia, lb, jb)] =
                    simdkin::one_center_kinetic_energy(a_basis.functions()[i], b_basis.functions()[j]);
            }
        }
    }
}

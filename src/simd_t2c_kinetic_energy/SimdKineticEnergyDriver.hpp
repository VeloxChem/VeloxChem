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



#ifndef SimdKineticEnergyDriver_hpp
#define SimdKineticEnergyDriver_hpp

#include <algorithm>
#include <array>
#include <cstddef>
#include <ranges>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdKineticEnergyBufferRows.hpp"
#include "SimdKineticEnergyFunc.hpp"
#include "SimdT2CDistributor.hpp"
#include "SparseMatrix.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

/// @brief Class CSimdKineticEnergyDriver computes the two-center kinetic energy
/// integrals of a molecular basis and hands them to a distributor, using a sparsity
/// pattern to skip the atom pairs and the combinations of basis functions whose
/// integrals are below the threshold.
///
/// @note The kinetic energy operator is symmetric in the two sides, so the driver
/// takes a single molecular basis and describes a symmetric quantity. There is no
/// form which takes a pair of bases.
class CSimdKineticEnergyDriver
{
   public:
    /// @brief The constructor with screening threshold and target block size.
    /// @param threshold The screening threshold of the integrals.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdKineticEnergyDriver(const double threshold = 1.0e-14, const size_t block_size = 0)

        : _threshold(threshold)

        , _block_size(block_size)
    {
    }

    /// @brief Gets screening threshold of the integrals.
    /// @return The screening threshold.
    auto
    get_threshold() const -> double
    {
        return _threshold;
    }

    /// @brief Gets target number of atom pairs of a block.
    /// @return The target number of atom pairs, zero if it is chosen automatically.
    auto
    get_block_size() const -> size_t
    {
        return _block_size;
    }

    /// @brief Creates the sparsity pattern the driver computes in.
    /// @param molecule The molecule to compute the interatomic distances from.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The sparsity pattern.
    auto
    make_pattern(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparsityPattern
    {
        return sparsity::make_pattern(
            molecule, basis, basis, screener::kinetic_energy, _threshold, mat_t::symmetric, diagstor::scalar, _block_size);
    }

    /// @brief Computes the kinetic energy integrals of a sparsity pattern and hands
    /// them to a distributor.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param distributor The distributor to hand the integrals to.
    /// @note This is the form which does not name a container of the values. The
    /// overloads which return a sparse matrix are written in terms of it.
    /// @note The pattern carries the threshold its atom pairs were screened with and
    /// the integrals are screened with that one, so the threshold of the driver is
    /// the one a pattern is formed with and not the one a pattern is computed with.
    template <class D>
    auto
    compute(const CSparsityPattern &pattern, const CMolecule &molecule, const CMolecularBasis &basis, D &distributor) const
        -> void
    {
        // NOTE: the basis functions of an atom basis are indexed once here rather
        // than once per block, as the index depends on the atom basis alone and the
        // blocks of one pair of atom bases are many.

        const auto indices = denseidx::index_functions(basis);

        sparsity::check_pattern(pattern, indices, indices);

        _compute_pair_blocks(pattern, molecule, basis, indices, distributor);

        _compute_diagonal_blocks(pattern, basis, indices, distributor);
    }

    /// @brief Computes the kinetic energy matrix of a sparsity pattern.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the kinetic energy matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse kinetic energy matrix of the pattern.
    /// @note This is the form which takes a pattern the caller already holds, so
    /// that a pattern is formed once and the matrices of several operators of the
    /// same basis and threshold are computed in it.
    auto compute_matrix(const CSparsityPattern &pattern, const CMolecule &molecule, const CMolecularBasis &basis) const
        -> CSparseMatrix;

    /// @brief Computes the kinetic energy matrix of a molecular basis.
    /// @param molecule The molecule to compute the kinetic energy matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse kinetic energy matrix.
    auto compute_matrix(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix;

   private:
    /// @brief Computes the integrals of the off-diagonal blocks.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param indices The index of the basis functions of the molecular basis.
    /// @param distributor The distributor to hand the integrals to.
    template <class D>
    auto _compute_pair_blocks(const CSparsityPattern              &pattern,
                              const CMolecule                     &molecule,
                              const CMolecularBasis               &basis,
                              const denseidx::TBasisFunctionIndex &indices,
                              D                                   &distributor) const -> void;

    /// @brief Computes the integrals of the diagonal blocks.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param indices The index of the basis functions of the molecular basis.
    /// @param distributor The distributor to hand the integrals to.
    /// @note The kinetic energy of two basis functions on the same atom does not
    /// depend on the position of the atom, so the molecule is not needed here.
    template <class D>
    auto _compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                  const CMolecularBasis               &basis,
                                  const denseidx::TBasisFunctionIndex &indices,
                                  D                                   &distributor) const -> void;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block, zero to choose it from
    /// the number of the threads and the number of the atom pairs.
    size_t _block_size;
};

template <class D>
auto
CSimdKineticEnergyDriver::_compute_pair_blocks(const CSparsityPattern              &pattern,
                                               const CMolecule                     &molecule,
                                               const CMolecularBasis               &basis,
                                               const denseidx::TBasisFunctionIndex &indices,
                                               D                                   &distributor) const -> void
{
    // NOTE: the blocks are independent, as each of them forms its own coordinates
    // and hands the distributor the values of its own combinations of basis
    // functions, which no other block addresses. Dynamic scheduling is used as the
    // blocks hold a comparable number of atom pairs but differ in the number of the
    // combinations of basis functions and in the cost of their kernels.

    const auto nblocks = static_cast<int>(pattern.number_of_pair_blocks());

    // NOTE: the blocks are visited from the most costly to the least, so that a
    // costly block is taken while there is still work to fill the other threads
    // with. The weight is the number of atom pairs surviving the screening of each
    // combination of basis functions, by the spherical components it carries and by
    // the pairs of primitives it sums over. It orders the blocks and is used for
    // nothing else.

    std::vector<size_t> order(static_cast<size_t>(nblocks));

    std::vector<double> costs(static_cast<size_t>(nblocks), 0.0);

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    // NOTE: the cost of a block is read off the block and not recomputed here, and
    // the shape of the arena is gathered in the same pass. See the overlap driver.

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(static_cast<size_t>(iblk));

        order[static_cast<size_t>(iblk)] = static_cast<size_t>(iblk);

        costs[static_cast<size_t>(iblk)] = block.weight();

        arena_rows = std::max(arena_rows,
                              simdkin::number_of_buffer_rows(basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                             basis.basis_set(block.ket_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });


    // NOTE: the arena is formed once per thread and not once per block. Its
    // largest shape serves every block, and holding it over the whole loop costs
    // no more memory than a block at a time did, as every thread held one of them
    // at once in any case. What it saves is the allocation and the page faults of
    // a mapping this large, which the cache of blocks is too small to hold back.

#pragma omp parallel if (nblocks > 1)
    {
        auto arena = CSimdMatrix(arena_rows, arena_cols);

        // NOTE: the combinations work over a view of the arena and not over the
        // arena itself, so that each takes the shape of the atom pairs it reaches.
        // See the overlap driver for why.

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            const auto jblk = order[static_cast<size_t>(iblk)];

            const auto &block = pattern.pair_block(jblk);

            // NOTE: the coordinates of the atom pairs are created once for the whole
            // block, as all combinations of basis functions of the block share them.

            const auto coordinates = simdfunc::make_coordinates(block, molecule);

            const auto &a_basis = basis.basis_set(block.bra_index());

            const auto &b_basis = basis.basis_set(block.ket_index());

            const auto &a_index = indices[block.bra_index()];

            const auto &b_index = indices[block.ket_index()];

            // NOTE: the atom bases of an off-diagonal block sit on different atoms, so
            // all combinations of basis functions are computed and none of them shares
            // the storage of its values with the reverse order.

            for (size_t i = 0; i < a_index.size(); i++)
            {
                for (size_t j = 0; j < b_index.size(); j++)
                {
                    const auto [la, ia] = a_index[i];

                    const auto [lb, jb] = b_index[j];

                    const auto nvalues = block.number_of_pairs(la, ia, lb, jb);

                    if (nvalues == 0) continue;

                    const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

                    auto *values = distributor.target(block, jblk, la, ia, lb, jb, nvalues, ncomps);

                    // NOTE: the pairs of primitives are screened with the threshold the
                    // atom pairs of the pattern were screened with, so that the two
                    // screenings of a computation cannot disagree when the pattern comes
                    // from the caller rather than from this driver.

                    simdkin::compute_kinetic_energy(values,
                                                    nvalues,
                                                    a_basis.functions()[i],
                                                    b_basis.functions()[j],
                                                    coordinates,
                                                    buffer,
                                                    pattern.get_threshold());

                    distributor.commit(block, jblk, la, ia, lb, jb, nvalues, ncomps);
                }
            }
        }
    }
}

template <class D>
auto
CSimdKineticEnergyDriver::_compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                                   const CMolecularBasis               &basis,
                                                   const denseidx::TBasisFunctionIndex &indices,
                                                   D                                   &distributor) const -> void
{
    // NOTE: the kinetic energy operator is spherically symmetric about an atom, so a
    // single value is stored for each pair of basis functions with the same angular
    // momentum and the position of the atom does not enter.

    // NOTE: the loop is serial. The diagonal blocks are as many as the unique atom
    // bases of the molecule, which is set by the variety of its elements and not by
    // its size, and each combination is one scalar over the pairs of primitives.

    for (size_t iblk = 0; iblk < pattern.number_of_diagonal_blocks(); iblk++)
    {
        const auto &block = pattern.diagonal_block(iblk);

        const auto &a_basis = basis.basis_set(block.bra_index());

        const auto &b_basis = basis.basis_set(block.ket_index());

        const auto &a_index = indices[block.bra_index()];

        const auto &b_index = indices[block.ket_index()];

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                // NOTE: only the stored combinations are computed, as the values of
                // the reverse order share their storage.

                if (block.is_triangular() && (i > j)) continue;

                const auto [la, ia] = a_index[i];

                const auto [lb, jb] = b_index[j];

                if (block.number_of_elements(la, ia, lb, jb) == 0) continue;

                auto *values = distributor.diagonal_target(block, iblk, la, ia, lb, jb);

                *values = simdkin::one_center_kinetic_energy(a_basis.functions()[i], b_basis.functions()[j]);

                distributor.diagonal_commit(block, iblk, la, ia, lb, jb);
            }
        }
    }
}

#endif /* SimdKineticEnergyDriver_hpp */

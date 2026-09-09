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


#ifndef SimdOverlapDriver_hpp
#define SimdOverlapDriver_hpp

#include <algorithm>
#include <cstddef>
#include <map>
#include <ranges>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdOverlapBufferRows.hpp"
#include "SimdOverlapFunc.hpp"
#include "SimdT2CDistributor.hpp"
#include "SparseMatrix.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

/// @brief Class CSimdOverlapDriver computes the two-center overlap integrals of
/// a molecular basis and stores them in a sparse matrix, using the sparsity
/// patterns of the atom basis pair groups to skip the atom pairs and the
/// combinations of basis functions whose integrals are below the threshold.
class CSimdOverlapDriver
{
    /// @brief One combination of basis functions of one block, which is the unit of
    /// work the threads draw on.
    struct TPairTask
    {
        /// @brief The index of the block among the off-diagonal blocks.
        size_t iblock;

        /// @brief The index of the basis function on bra side within its atom basis.
        size_t i;

        /// @brief The index of the basis function on ket side within its atom basis.
        size_t j;

        /// @brief The number of atom pairs the combination reaches.
        size_t nvalues;
    };

   public:
    /// @brief The constructor with screening threshold and target block size.
    /// @param threshold The screening threshold of the integrals.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdOverlapDriver(const double threshold = 1.0e-14, const size_t block_size = 0)

        : _threshold(threshold)

        , _block_size(block_size)
    {
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
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param mat_type The type of the quantity to describe.
    /// @return The sparsity pattern.
    auto
    make_pattern(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis, const mat_t mat_type)
        const -> CSparsityPattern
    {
        return sparsity::make_pattern(
            molecule, bra_basis, ket_basis, screener::overlap, _threshold, mat_type, diagstor::scalar, _block_size);
    }

    /// @brief Computes the overlap integrals of a sparsity pattern and hands them to
    /// a distributor.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param distributor The distributor to hand the integrals to.
    /// @note This is the form which does not name a container of the values. The
    /// overloads which return a sparse matrix are written in terms of it.
    /// @note The pattern carries the threshold its atom pairs were screened with and
    /// the integrals are screened with that one, so the threshold of the driver is
    /// the one a pattern is formed with and not the one a pattern is computed with.
    template <class D>
    auto
    compute(const CSparsityPattern &pattern,
            const CMolecule        &molecule,
            const CMolecularBasis  &bra_basis,
            const CMolecularBasis  &ket_basis,
            D                      &distributor) const -> void
    {
        // NOTE: the basis functions of an atom basis are indexed once here rather
        // than once per block, as the index depends on the atom basis alone and the
        // blocks of one pair of atom bases are many.

        const auto a_indices = denseidx::index_functions(bra_basis);

        const auto b_indices = denseidx::index_functions(ket_basis);

        sparsity::check_pattern(pattern, a_indices, b_indices);

        _compute_pair_blocks(pattern, molecule, bra_basis, ket_basis, a_indices, b_indices, distributor);

        _compute_diagonal_blocks(pattern, bra_basis, ket_basis, a_indices, b_indices, distributor);
    }

    /// @brief Computes the overlap matrix of a sparsity pattern.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the overlap matrix of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @return The sparse overlap matrix of the pattern.
    /// @note This is the form which takes a pattern the caller already holds, so
    /// that a pattern is formed once and the matrices of several operators of the
    /// same basis and threshold are computed in it. The overloads below form the
    /// pattern themselves and are written in terms of it.
    auto compute_matrix(const CSparsityPattern &pattern,
                        const CMolecule        &molecule,
                        const CMolecularBasis  &bra_basis,
                        const CMolecularBasis  &ket_basis) const -> CSparseMatrix;

    /// @brief Computes the overlap matrix of a molecular basis.
    /// @param molecule The molecule to compute the overlap matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse overlap matrix.
    auto compute_matrix(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix;

    /// @brief Computes the overlap matrix of a pair of molecular bases.
    /// @param molecule The molecule to compute the overlap matrix of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @return The general sparse overlap matrix.
    auto compute_matrix(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
        -> CSparseMatrix;

    /// @brief Gets screening threshold of the integrals.
    /// @return The screening threshold.
    auto
    get_threshold() const -> double
    {
        return _threshold;
    }

   private:
    /// @brief Computes the integrals of the off-diagonal blocks.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param a_indices The index of the basis functions of the basis on bra side.
    /// @param b_indices The index of the basis functions of the basis on ket side.
    /// @param distributor The distributor to hand the integrals to.
    template <class D>
    auto _compute_pair_blocks(const CSparsityPattern              &pattern,
                              const CMolecule                     &molecule,
                              const CMolecularBasis               &bra_basis,
                              const CMolecularBasis               &ket_basis,
                              const denseidx::TBasisFunctionIndex &a_indices,
                              const denseidx::TBasisFunctionIndex &b_indices,
                              D                                   &distributor) const -> void;

    /// @brief Computes the integrals of the diagonal blocks.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param a_indices The index of the basis functions of the basis on bra side.
    /// @param b_indices The index of the basis functions of the basis on ket side.
    /// @param distributor The distributor to hand the integrals to.
    /// @note The overlap of two basis functions on the same atom does not depend on
    /// the position of the atom, so the molecule is not needed here.
    template <class D>
    auto _compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                  const CMolecularBasis               &bra_basis,
                                  const CMolecularBasis               &ket_basis,
                                  const denseidx::TBasisFunctionIndex &a_indices,
                                  const denseidx::TBasisFunctionIndex &b_indices,
                                  D                                   &distributor) const -> void;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block, zero to choose it from
    /// the number of the threads and the number of the atom pairs.
    size_t _block_size;
};

template <class D>
auto
CSimdOverlapDriver::_compute_pair_blocks(const CSparsityPattern              &pattern,
                                         const CMolecule                     &molecule,
                                         const CMolecularBasis               &bra_basis,
                                         const CMolecularBasis               &ket_basis,
                                         const denseidx::TBasisFunctionIndex &a_indices,
                                         const denseidx::TBasisFunctionIndex &b_indices,
                                         D                                   &distributor) const -> void
{
    // NOTE: the unit of work is one combination of basis functions of one block and
    // not one block. A block is what the sparsity and the storage are described in,
    // and it carries a fixed cost of some microseconds, so it cannot be made small
    // enough to feed a large machine: an ordinary molecule holds tens of blocks
    // whatever the number of the threads, as the target size of a block is bounded
    // from below. The combinations of a block are tens to hundreds and cost nothing
    // to enumerate, and each of them writes its own values and reads the coordinates
    // of its block without changing them, so they are independent of one another.

    const auto nblocks = static_cast<size_t>(pattern.number_of_pair_blocks());

    if (nblocks == 0) return;

    // NOTE: the blocks are visited from the most costly to the least, and within a
    // block the combinations likewise, so that a costly task is taken while there is
    // still work to fill the other threads with. A costly task drawn last is
    // finished alone and sets the time of the whole loop.

    std::vector<size_t> order(nblocks);

    std::vector<double> costs(nblocks, 0.0);

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(iblk);

        order[iblk] = iblk;

        costs[iblk] = block.weight();

        arena_rows = std::max(arena_rows,
                              simdovl::number_of_buffer_rows(bra_basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                             ket_basis.basis_set(block.ket_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });

    // NOTE: the order of the combinations by cost depends on the pair of atom bases
    // and not on the block, as the components and the primitives of a combination are
    // properties of its two basis functions. It is therefore formed once for every
    // pair of atom bases the blocks draw on, and not once per block. Sorting the
    // tasks themselves would cost more than the ordering saves, as they are hundreds
    // of thousands of them for a large fitting set.

    std::map<std::pair<int, int>, std::vector<std::pair<size_t, size_t>>> orders;

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(iblk);

        const auto key = std::pair<int, int>{block.bra_index(), block.ket_index()};

        if (orders.contains(key)) continue;

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        const auto &a_index = a_indices[block.bra_index()];

        const auto &b_index = b_indices[block.ket_index()];

        std::vector<std::pair<size_t, size_t>> combinations;

        std::vector<double> weights;

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                const auto [la, ia] = a_index[i];

                const auto [lb, jb] = b_index[j];

                const auto ncomps = static_cast<double>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

                const auto nprims = static_cast<double>(a_basis.functions()[i].exponents().size() *
                                                        b_basis.functions()[j].exponents().size());

                combinations.emplace_back(i, j);

                weights.push_back(ncomps * nprims);
            }
        }

        std::vector<size_t> perm(combinations.size());

        std::ranges::copy(std::views::iota(size_t{0}, combinations.size()), perm.begin());

        std::ranges::sort(perm, [&](const size_t a, const size_t b) { return weights[a] > weights[b]; });

        std::vector<std::pair<size_t, size_t>> sorted;

        sorted.reserve(perm.size());

        for (const auto k : perm) sorted.push_back(combinations[k]);

        orders.emplace(key, std::move(sorted));
    }

    // NOTE: a combination which reaches no atom pair is left out here rather than
    // skipped in the loop, so that every task the threads draw carries work.

    std::vector<TPairTask> tasks;

    for (const auto iblk : order)
    {
        const auto &block = pattern.pair_block(iblk);

        const auto &a_index = a_indices[block.bra_index()];

        const auto &b_index = b_indices[block.ket_index()];

        for (const auto [i, j] : orders.at({block.bra_index(), block.ket_index()}))
        {
            const auto [la, ia] = a_index[i];

            const auto [lb, jb] = b_index[j];

            if (const auto nvalues = block.number_of_pairs(la, ia, lb, jb); nvalues > 0)
            {
                tasks.push_back({iblk, i, j, nvalues});
            }
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return;

    // NOTE: the coordinates of a block are formed before the tasks are drawn and not
    // inside a task, as the combinations of a block are computed by different threads
    // and all of them read the same coordinates.

    std::vector<CSimdMatrix> coordinates(nblocks);

    const auto nblk = static_cast<int>(nblocks);

#pragma omp parallel if (ntasks > 1)
    {
#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblk; iblk++)
        {
            coordinates[static_cast<size_t>(iblk)] = simdfunc::make_coordinates(pattern.pair_block(static_cast<size_t>(iblk)), molecule);
        }

        // NOTE: the arena is formed once per thread and spans the largest combination
        // any block carries. The combinations work over a view of it, so that each
        // takes the shape of the atom pairs it reaches.

        auto arena = CSimdMatrix(arena_rows, arena_cols);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int itask = 0; itask < ntasks; itask++)
        {
            const auto &task = tasks[static_cast<size_t>(itask)];

            const auto &block = pattern.pair_block(task.iblock);

            const auto &a_basis = bra_basis.basis_set(block.bra_index());

            const auto &b_basis = ket_basis.basis_set(block.ket_index());

            const auto [la, ia] = a_indices[block.bra_index()][task.i];

            const auto [lb, jb] = b_indices[block.ket_index()][task.j];

            const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

            auto *values = distributor.target(block, task.iblock, la, ia, lb, jb, task.nvalues, ncomps);

            // NOTE: the pairs of primitives are screened with the threshold the atom
            // pairs of the pattern were screened with, so that the two screenings of a
            // computation cannot disagree when the pattern comes from the caller
            // rather than from this driver.

            simdovl::compute_overlap(values,
                                     task.nvalues,
                                     a_basis.functions()[task.i],
                                     b_basis.functions()[task.j],
                                     coordinates[task.iblock],
                                     buffer,
                                     pattern.get_threshold());

            distributor.commit(block, task.iblock, la, ia, lb, jb, task.nvalues, ncomps);
        }
    }
}

template <class D>
auto
CSimdOverlapDriver::_compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                             const CMolecularBasis               &bra_basis,
                                             const CMolecularBasis               &ket_basis,
                                             const denseidx::TBasisFunctionIndex &a_indices,
                                             const denseidx::TBasisFunctionIndex &b_indices,
                                             D                                   &distributor) const -> void
{
    // NOTE: the overlap operator is spherically symmetric about an atom, so a single
    // value is stored for each pair of basis functions with the same angular
    // momentum and the position of the atom does not enter.

    for (size_t iblk = 0; iblk < pattern.number_of_diagonal_blocks(); iblk++)
    {
        const auto &block = pattern.diagonal_block(iblk);

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        const auto &a_index = a_indices[block.bra_index()];

        const auto &b_index = b_indices[block.ket_index()];

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

                *values = simdovl::one_center_overlap(a_basis.functions()[i], b_basis.functions()[j]);

                distributor.diagonal_commit(block, iblk, la, ia, lb, jb);
            }
        }
    }
}

#endif /* SimdOverlapDriver_hpp */

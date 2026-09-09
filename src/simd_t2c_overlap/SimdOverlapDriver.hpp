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
    // NOTE: the blocks are independent, as each of them forms its own coordinates
    // and hands the distributor the values of its own combinations of basis
    // functions, which no other block addresses. Dynamic scheduling is used as the
    // blocks hold a comparable number of atom pairs but differ in the number of the
    // combinations of basis functions and in the cost of their kernels.

    const auto nblocks = static_cast<int>(pattern.number_of_pair_blocks());

    // NOTE: the blocks are visited from the most costly to the least, so that a
    // costly block is taken while there is still work to fill the other threads
    // with. The threads draw two or three blocks each, so a costly block drawn last
    // is finished alone and sets the time of the whole loop.

    // NOTE: the cost of a block is the number of atom pairs surviving the screening
    // of each of its combinations of basis functions, weighted by the number of the
    // spherical components the combination carries and by the number of the pairs of
    // primitives it sums over. The weight is an estimate and orders the blocks, it is
    // not used for anything else.

    std::vector<size_t> order(static_cast<size_t>(nblocks));

    std::vector<double> costs(static_cast<size_t>(nblocks), 0.0);

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    // NOTE: the cost of a block is read off the block and not recomputed here. It is
    // accumulated where the sparsity of the block is described, which is inside a
    // parallel region, so this pass is a read of one number per block rather than a
    // walk over every combination of basis functions of every block. That walk was
    // serial and its cost grew with the number of threads, as the number of blocks
    // does.

    // NOTE: the shape of the arena is gathered in the same pass. It spans the largest
    // combination any block carries, which is the one of the highest angular momenta
    // of the two atom bases of the block, as the rows a combination needs do not
    // decrease with either momentum.

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(static_cast<size_t>(iblk));

        order[static_cast<size_t>(iblk)] = static_cast<size_t>(iblk);

        costs[static_cast<size_t>(iblk)] = block.weight();

        arena_rows = std::max(arena_rows,
                              simdovl::number_of_buffer_rows(bra_basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                             ket_basis.basis_set(block.ket_index()).max_angular_momentum()));

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
        // A view stretched to the pairs of the whole block would leave every row of
        // a combination which reaches fewer of them a page away from the next, and
        // the zeroing alone then costs more than the allocations this saves.

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            const auto jblk = order[static_cast<size_t>(iblk)];

            const auto &block = pattern.pair_block(jblk);

            // NOTE: the coordinates of the atom pairs are created once for the whole
            // block, as all combinations of basis functions of the block share them.

            const auto coordinates = simdfunc::make_coordinates(block, molecule);

            const auto &a_basis = bra_basis.basis_set(block.bra_index());

            const auto &b_basis = ket_basis.basis_set(block.ket_index());

            const auto &a_index = a_indices[block.bra_index()];

            const auto &b_index = b_indices[block.ket_index()];

            // NOTE: the atom bases of an off-diagonal block sit on different atoms, so
            // all combinations of basis functions are computed and none of them shares
            // the storage of its values with the reverse order.

            // NOTE: the combinations of basis functions are independent, as each of them
            // writes its own values and reads the coordinates of the block without
            // changing them.

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

                    simdovl::compute_overlap(values,
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

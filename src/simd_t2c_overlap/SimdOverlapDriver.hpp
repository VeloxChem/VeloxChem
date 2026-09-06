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
    template <class D>
    auto
    compute(const CSparsityPattern &pattern,
            const CMolecule        &molecule,
            const CMolecularBasis  &bra_basis,
            const CMolecularBasis  &ket_basis,
            D                      &distributor) const -> void
    {
        _compute_pair_blocks(pattern, molecule, bra_basis, ket_basis, distributor);

        _compute_diagonal_blocks(pattern, bra_basis, ket_basis, distributor);
    }

    /// @brief Computes the overlap matrix of a molecular basis.
    /// @param molecule The molecule to compute the overlap matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse overlap matrix.
    auto compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix;

    /// @brief Computes the overlap matrix of a pair of molecular bases.
    /// @param molecule The molecule to compute the overlap matrix of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @return The general sparse overlap matrix.
    auto compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
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
    /// @param distributor The distributor to hand the integrals to.
    template <class D>
    auto _compute_pair_blocks(const CSparsityPattern &pattern,
                              const CMolecule        &molecule,
                              const CMolecularBasis  &bra_basis,
                              const CMolecularBasis  &ket_basis,
                              D                      &distributor) const -> void;

    /// @brief Computes the integrals of the diagonal blocks.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param distributor The distributor to hand the integrals to.
    /// @note The overlap of two basis functions on the same atom does not depend on
    /// the position of the atom, so the molecule is not needed here.
    template <class D>
    auto _compute_diagonal_blocks(const CSparsityPattern &pattern,
                                  const CMolecularBasis  &bra_basis,
                                  const CMolecularBasis  &ket_basis,
                                  D                      &distributor) const -> void;

    /// @brief Indexes the basis functions of every unique atom basis of a molecular
    /// basis by their angular momentum and their order within it.
    /// @param basis The molecular basis to index the atom bases of.
    /// @return The vector of indices, one entry per unique atom basis.
    static auto
    _index_functions(const CMolecularBasis &basis) -> std::vector<std::vector<std::pair<int, size_t>>>
    {
        std::vector<std::vector<std::pair<int, size_t>>> indices;

        for (const auto &atom_basis : basis.basis_sets())
        {
            indices.push_back(denseidx::index_functions(atom_basis));
        }

        return indices;
    }

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block, zero to choose it from
    /// the number of the threads and the number of the atom pairs.
    size_t _block_size;
};

template <class D>
auto
CSimdOverlapDriver::_compute_pair_blocks(const CSparsityPattern &pattern,
                                         const CMolecule        &molecule,
                                         const CMolecularBasis  &bra_basis,
                                         const CMolecularBasis  &ket_basis,
                                         D                      &distributor) const -> void
{
    // NOTE: the blocks are independent, as each of them forms its own coordinates
    // and hands the distributor the values of its own combinations of basis
    // functions, which no other block addresses. Dynamic scheduling is used as the
    // blocks hold a comparable number of atom pairs but differ in the number of the
    // combinations of basis functions and in the cost of their kernels.

    const auto nblocks = static_cast<int>(pattern.number_of_pair_blocks());

    // NOTE: the basis functions of an atom basis are indexed once per atom basis
    // rather than once per block, as the index depends on the atom basis alone and
    // the blocks of one pair of atom bases are many.

    const auto a_indices = _index_functions(bra_basis);

    const auto b_indices = _index_functions(ket_basis);

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

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(static_cast<size_t>(iblk));

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        double weight = 0.0;

        for (size_t i = 0; i < a_indices[block.bra_index()].size(); i++)
        {
            for (size_t j = 0; j < b_indices[block.ket_index()].size(); j++)
            {
                const auto [la, ia] = a_indices[block.bra_index()][i];

                const auto [lb, jb] = b_indices[block.ket_index()][j];

                const auto ncomps = static_cast<double>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

                const auto nprims = static_cast<double>(a_basis.functions()[i].exponents().size() *
                                                        b_basis.functions()[j].exponents().size());

                weight += static_cast<double>(block.number_of_pairs(la, ia, lb, jb)) * ncomps * nprims;
            }
        }

        order[static_cast<size_t>(iblk)] = static_cast<size_t>(iblk);

        costs[static_cast<size_t>(iblk)] = weight;
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
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

                simdovl::compute_overlap(values, nvalues, a_basis.functions()[i], b_basis.functions()[j], coordinates, _threshold);

                distributor.commit(block, jblk, la, ia, lb, jb, nvalues, ncomps);
            }
        }
    }
}

template <class D>
auto
CSimdOverlapDriver::_compute_diagonal_blocks(const CSparsityPattern &pattern,
                                             const CMolecularBasis  &bra_basis,
                                             const CMolecularBasis  &ket_basis,
                                             D                      &distributor) const -> void
{
    // NOTE: the overlap operator is spherically symmetric about an atom, so a single
    // value is stored for each pair of basis functions with the same angular
    // momentum and the position of the atom does not enter.

    const auto a_indices = _index_functions(bra_basis);

    const auto b_indices = _index_functions(ket_basis);

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

                *distributor.diagonal_target(block, iblk, la, ia, lb, jb) =
                    simdovl::one_center_overlap(a_basis.functions()[i], b_basis.functions()[j]);
            }
        }
    }
}

#endif /* SimdOverlapDriver_hpp */

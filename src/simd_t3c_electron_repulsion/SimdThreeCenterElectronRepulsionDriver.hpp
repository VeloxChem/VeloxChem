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



#ifndef SimdThreeCenterElectronRepulsionDriver_hpp
#define SimdThreeCenterElectronRepulsionDriver_hpp

#include <array>
#include <cstddef>
#include <map>
#include <ranges>
#include <vector>

#include "AtomBasisTripleSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdMatrix.hpp"
#include "SimdT3CDistributor.hpp"
#include "SimdThreeCenterElectronRepulsionBufferRows.hpp"
#include "SimdThreeCenterElectronRepulsionFunc.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"
#include "TensorComponents.hpp"

/// @brief Class CSimdThreeCenterElectronRepulsionDriver computes the three-center
/// electron repulsion integrals of a molecular basis and an auxiliary molecular
/// basis and hands them to a distributor.
///
/// @note The threshold is an argument of compute rather than a property of the
/// driver, as the tensor of one molecule is built at several thresholds and the
/// bound is what decides its storage.
/// @note Only the charge distribution of the atom pair on a and b sides screens. The
/// Coulomb operator decays as the inverse of the distance to the atom on c side, so
/// no atom on c side falls below a threshold and every one of them is computed.
/// @note The atoms on c side may be given, in which case the tensor holds the part
/// of the whole which those atoms carry. This is how a tensor too large for the
/// memory is formed in batches: the caller divides the atoms and the parts are
/// computed one after another.
class CSimdThreeCenterElectronRepulsionDriver
{
    /// @brief One combination of basis functions of one block, which is the unit of
    /// work the threads draw on.
    struct TTripleTask
    {
        /// @brief The index of the block among the blocks of the pattern.
        size_t iblock;

        /// @brief The index of the basis function on a side within its atom basis.
        size_t i;

        /// @brief The index of the basis function on b side within its atom basis.
        size_t j;

        /// @brief The index of the basis function on c side within its atom basis.
        size_t k;

        /// @brief The number of atom pairs the combination reaches.
        size_t npairs;
    };

   public:
    /// @brief The default constructor.
    CSimdThreeCenterElectronRepulsionDriver() = default;

    /// @brief Creates the sparsity pattern the driver computes in, for all atoms on
    /// c side.
    /// @param molecule The molecule to compute the interatomic distances from.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param threshold The screening threshold.
    /// @return The sparsity pattern.
    auto
    make_pattern(const CMolecule &molecule, const CMolecularBasis &basis, const CMolecularBasis &aux_basis,
                 const double threshold) const -> CTripleSparsityPattern
    {
        return sparsity::make_triple_pattern(
            molecule, basis, aux_basis, screenfunc::three_center_electron_repulsion_bound, threshold, mat_t::symmetric);
    }

    /// @brief Creates the sparsity pattern the driver computes in, for the given
    /// atoms on c side.
    /// @param atoms The atoms on c side, as their indices in the molecule.
    auto
    make_pattern(const CMolecule &molecule, const CMolecularBasis &basis, const CMolecularBasis &aux_basis,
                 const double threshold, const std::vector<int> &atoms) const -> CTripleSparsityPattern
    {
        return sparsity::make_triple_pattern(molecule,
                                             basis,
                                             aux_basis,
                                             screenfunc::three_center_electron_repulsion_bound,
                                             threshold,
                                             mat_t::symmetric,
                                             atoms);
    }

    /// @brief Computes the integrals of a sparsity pattern and hands them to a
    /// distributor.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to take the atomic coordinates from.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param distributor The distributor to hand the integrals to.
    /// @note This is the form which does not name a container of the values. The
    /// overloads which return a sparse tensor are written in terms of it.
    template <class D>
    auto compute(const CTripleSparsityPattern &pattern,
                 const CMolecule              &molecule,
                 const CMolecularBasis        &basis,
                 const CMolecularBasis        &aux_basis,
                 D                            &distributor) const -> void;

    /// @brief Computes the three-center electron repulsion integrals of a molecular
    /// basis and an auxiliary molecular basis, for all atoms on c side.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param threshold The screening threshold.
    /// @return The sparse tensor of the integrals.
    auto compute(const CMolecule       &molecule,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const double           threshold) const -> CSparseTensor;

    /// @brief Computes the three-center electron repulsion integrals of a molecular
    /// basis and an auxiliary molecular basis, for the given atoms on c side.
    /// @param molecule The molecule to compute the integrals of.
    /// @param basis The molecular basis on a and b sides.
    /// @param aux_basis The auxiliary molecular basis on c side.
    /// @param threshold The screening threshold.
    /// @param atoms The atoms on c side to compute the integrals of, as their indices
    /// in the molecule and without repetition.
    /// @return The sparse tensor of the integrals.
    auto compute(const CMolecule        &molecule,
                 const CMolecularBasis  &basis,
                 const CMolecularBasis  &aux_basis,
                 const double            threshold,
                 const std::vector<int> &atoms) const -> CSparseTensor;

   private:
    /// @brief Creates the coordinates of the atoms on c side of a block.
    /// @param block The sparsity pattern of the block.
    /// @param molecule The molecule to take the atomic coordinates from.
    /// @return The matrix of three rows and as many columns as there are atoms on c
    /// side, holding their coordinates in atomic units.
    /// @note The atoms on c side are a separate and shorter dimension than the atom
    /// pairs, so they carry their own matrix rather than rows of the coordinates of
    /// the pairs.
    static auto _make_c_coordinates(const CAtomBasisTripleSparsity &block, const CMolecule &molecule) -> CSimdMatrix;
};

template <class D>
auto
CSimdThreeCenterElectronRepulsionDriver::compute(const CTripleSparsityPattern &pattern,
                                                 const CMolecule              &molecule,
                                                 const CMolecularBasis        &basis,
                                                 const CMolecularBasis        &aux_basis,
                                                 D                            &distributor) const -> void
{
    // NOTE: the basis functions of an atom basis are indexed once here rather than
    // once per block, as the index depends on the atom basis alone and the blocks of
    // one triple of atom bases are many.

    const auto indices = denseidx::index_functions(basis);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    sparsity::check_triple_pattern(pattern, indices, aux_indices);

    const auto nblocks = static_cast<size_t>(pattern.number_of_blocks());

    if (nblocks == 0) return;

    // NOTE: the unit of work is one combination of basis functions of one block and
    // not one block. A block is what the sparsity and the storage are described in
    // and it carries a fixed cost, so it cannot be made small enough to feed a large
    // machine, while the combinations of a block are the product of the basis
    // functions of its three atom bases and are many. Each of them writes its own
    // values and reads the coordinates of its block without changing them, so they
    // are independent of one another. See the overlap driver.

    // NOTE: the blocks are drawn from the largest to the smallest, and within a block
    // the combinations from the most costly to the least, so that a costly task is
    // taken while there is still work to fill the other threads with. The cost of a
    // block is estimated from its atom pairs and the atoms it carries on c side,
    // which needs no walk over its combinations.

    // NOTE: the arena spans the largest combination any block carries, which is the
    // one of the highest angular momenta of its three atom bases, as the rows a
    // combination needs do not decrease with any of them. It is formed once per
    // thread and the combinations work over a view of it. See the overlap driver.

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.block(iblk);

        arena_rows = std::max(arena_rows,
                              simdt3ceri::number_of_buffer_rows(basis.basis_set(block.a_index()).max_angular_momentum(),
                                                                basis.basis_set(block.b_index()).max_angular_momentum(),
                                                                aux_basis.basis_set(block.c_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    std::vector<size_t> border(nblocks);

    std::ranges::copy(std::views::iota(size_t{0}, nblocks), border.begin());

    std::ranges::sort(border, [&](const size_t a, const size_t b) {
        return pattern.block(a).number_of_pairs() * pattern.block(a).number_of_c_atoms() >
               pattern.block(b).number_of_pairs() * pattern.block(b).number_of_c_atoms();
    });

    // NOTE: the order of the combinations by cost belongs to the triple of atom bases
    // and not to the block, as the components and the primitives of a combination are
    // properties of its three basis functions. It is therefore formed once for every
    // triple the blocks draw on, and not once per block.

    std::map<std::array<int, 3>, std::vector<std::array<size_t, 3>>> orders;

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.block(iblk);

        const auto key = std::array<int, 3>{block.a_index(), block.b_index(), block.c_index()};

        if (orders.contains(key)) continue;

        const auto &a_basis = basis.basis_set(block.a_index());

        const auto &b_basis = basis.basis_set(block.b_index());

        const auto &c_basis = aux_basis.basis_set(block.c_index());

        const auto &a_index = indices[static_cast<size_t>(block.a_index())];

        const auto &b_index = indices[static_cast<size_t>(block.b_index())];

        const auto &c_index = aux_indices[static_cast<size_t>(block.c_index())];

        std::vector<std::array<size_t, 3>> combinations;

        std::vector<double> weights;

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                for (size_t k = 0; k < c_index.size(); k++)
                {
                    const auto [la, ia] = a_index[i];

                    const auto [lb, jb] = b_index[j];

                    const auto [lc, kc] = c_index[k];

                    const auto ncomps =
                        static_cast<double>(tensor::number_of_spherical_components(std::array<int, 3>{la, lb, lc}));

                    const auto nprims = static_cast<double>(a_basis.functions()[i].exponents().size() *
                                                            b_basis.functions()[j].exponents().size() *
                                                            c_basis.functions()[k].exponents().size());

                    combinations.push_back({i, j, k});

                    weights.push_back(ncomps * nprims);
                }
            }
        }

        std::vector<size_t> perm(combinations.size());

        std::ranges::copy(std::views::iota(size_t{0}, combinations.size()), perm.begin());

        std::ranges::sort(perm, [&](const size_t a, const size_t b) { return weights[a] > weights[b]; });

        std::vector<std::array<size_t, 3>> sorted;

        sorted.reserve(perm.size());

        for (const auto m : perm) sorted.push_back(combinations[m]);

        orders.emplace(key, std::move(sorted));
    }

    // NOTE: a combination which reaches no atom pair is left out here rather than
    // skipped in the loop, so that every task the threads draw carries work.

    std::vector<TTripleTask> tasks;

    for (const auto iblk : border)
    {
        const auto &block = pattern.block(iblk);

        if ((block.number_of_pairs() == 0) || (block.number_of_c_atoms() == 0)) continue;

        const auto &a_index = indices[static_cast<size_t>(block.a_index())];

        const auto &b_index = indices[static_cast<size_t>(block.b_index())];

        const auto &c_index = aux_indices[static_cast<size_t>(block.c_index())];

        for (const auto [i, j, k] : orders.at({block.a_index(), block.b_index(), block.c_index()}))
        {
            const auto [la, ia] = a_index[i];

            const auto [lb, jb] = b_index[j];

            const auto [lc, kc] = c_index[k];

            if (const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc); npairs > 0)
            {
                tasks.push_back({iblk, i, j, k, npairs});
            }
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return;

    // NOTE: the coordinates of a block are formed before the tasks are drawn and not
    // inside a task, as the combinations of a block are computed by different threads
    // and all of them read the same coordinates.

    std::vector<CSimdMatrix> coordinates(nblocks);

    std::vector<CSimdMatrix> c_coordinates(nblocks);

    const auto nblk = static_cast<int>(nblocks);

#pragma omp parallel if (ntasks > 1)
    {
#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblk; iblk++)
        {
            const auto &block = pattern.block(static_cast<size_t>(iblk));

            if ((block.number_of_pairs() == 0) || (block.number_of_c_atoms() == 0)) continue;

            coordinates[static_cast<size_t>(iblk)] = simdfunc::make_coordinates(block, molecule);

            c_coordinates[static_cast<size_t>(iblk)] = _make_c_coordinates(block, molecule);
        }

        auto arena = CSimdMatrix(arena_rows, arena_cols);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int itask = 0; itask < ntasks; itask++)
        {
            const auto &task = tasks[static_cast<size_t>(itask)];

            const auto &block = pattern.block(task.iblock);

            const auto natoms = block.number_of_c_atoms();

            const auto &a_basis = basis.basis_set(block.a_index());

            const auto &b_basis = basis.basis_set(block.b_index());

            const auto &c_basis = aux_basis.basis_set(block.c_index());

            const auto [la, ia] = indices[static_cast<size_t>(block.a_index())][task.i];

            const auto [lb, jb] = indices[static_cast<size_t>(block.b_index())][task.j];

            const auto [lc, kc] = aux_indices[static_cast<size_t>(block.c_index())][task.k];

            const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 3>{la, lb, lc}));

            auto *values = distributor.target(block, task.iblock, la, ia, lb, jb, lc, kc, task.npairs, natoms, ncomps);

            simdt3ceri::compute_electron_repulsion(values,
                                                   task.npairs,
                                                   natoms,
                                                   a_basis.functions()[task.i],
                                                   b_basis.functions()[task.j],
                                                   c_basis.functions()[task.k],
                                                   coordinates[task.iblock],
                                                   c_coordinates[task.iblock],
                                                   buffer,
                                                   pattern.get_threshold());

            distributor.commit(block, task.iblock, la, ia, lb, jb, lc, kc, task.npairs, natoms, ncomps);
        }
    }
}

#endif /* SimdThreeCenterElectronRepulsionDriver_hpp */

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



#include "SimdTwoCenterElectronRepulsionDriver.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <ranges>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdMatrix.hpp"
#include "SimdTwoCenterElectronRepulsionBufferRows.hpp"
#include "SimdTwoCenterElectronRepulsionFunc.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

auto
CSimdTwoCenterElectronRepulsionDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CPackedMatrix
{
    // NOTE: the matrix is symmetric and is stored as its lower triangle. Its
    // constructor makes the one allocation and zeroes it with the threads, so the
    // atom pairs below write into a matrix which is already zero and the atoms which
    // no block reaches keep the zeros set there.

    auto matrix = CPackedMatrix(basis, mat_t::symmetric);

    // NOTE: only the geometry half of the pattern is formed. The blocks divide the
    // atom pairs for the threads and order them by interatomic distance; nothing is
    // described and nothing is screened, as no atom pair of the Coulomb operator
    // falls below a threshold.

    auto groups = basis.basis_pair_groups();

    const auto nblock_pairs = (_block_size == 0)
                                  ? CAtomBasisPairGroup::make_block_size(
                                        groups, sparsity::blocks_per_thread, sparsity::min_block_size, max_block_size)
                                  : _block_size;

    auto blocks = (nblock_pairs == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, nblock_pairs);

    CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

    _compute_pair_blocks(matrix, molecule, basis, blocks);

    _compute_diagonal_blocks(matrix, basis);

    return matrix;
}

auto
CSimdTwoCenterElectronRepulsionDriver::_compute_pair_blocks(CPackedMatrix                          &matrix,
                                                            const CMolecule                        &molecule,
                                                            const CMolecularBasis                  &basis,
                                                            const std::vector<CAtomBasisPairGroup> &blocks) const -> void
{
    // NOTE: the basis functions of an atom basis are indexed once here rather than
    // once per block, as the index depends on the atom basis alone and the blocks of
    // one pair of atom bases are many.

    const auto indices = denseidx::index_functions(basis);

    // NOTE: the dense index of an atomic orbital is looked up rather than
    // recomputed, as the innermost loop below runs over the atom pairs and would
    // otherwise scan the basis for every value.

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = strides.size();

    const auto nblocks = static_cast<int>(blocks.size());

    // NOTE: an element of the matrix belongs to a single atom pair, and an atom pair
    // belongs to a single block, so the blocks write to disjoint elements and need no
    // synchronization. The unit of work is one combination of basis functions of one
    // block and not one block: a block carries a fixed cost, so it cannot be made
    // small enough to feed a large machine, while the combinations of a block are
    // tens to hundreds and cost nothing to enumerate. See the overlap driver.

    // NOTE: the arena spans the largest combination of basis functions any block
    // carries, and every combination works over a view of it shaped to the atom pairs
    // it reaches. See the overlap driver for why it is shaped this way.

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = blocks[static_cast<size_t>(iblk)];

        arena_rows = std::max(arena_rows,
                              simdt2ceri::number_of_buffer_rows(basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                                basis.basis_set(block.ket_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    // NOTE: the order of the combinations by cost belongs to the pair of atom bases
    // and not to the block, so it is formed once per pair of atom bases rather than
    // once per block. Nothing is screened here, so a combination of a block always
    // reaches every atom pair of it and the cost of a task is the atom pairs of its
    // block times what its combination computes over them.

    std::map<std::pair<int, int>, std::vector<std::pair<size_t, size_t>>> orders;

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = blocks[static_cast<size_t>(iblk)];

        const auto key = std::pair<int, int>{block.bra_index(), block.ket_index()};

        if (orders.contains(key)) continue;

        const auto &a_basis = basis.basis_set(block.bra_index());

        const auto &b_basis = basis.basis_set(block.ket_index());

        const auto &a_index = indices[static_cast<size_t>(block.bra_index())];

        const auto &b_index = indices[static_cast<size_t>(block.ket_index())];

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

    // NOTE: the blocks are visited from the most costly to the least, and within a
    // block the combinations likewise, so that a costly task is drawn while there is
    // still work to fill the other threads with.

    std::vector<size_t> border(static_cast<size_t>(nblocks));

    std::ranges::copy(std::views::iota(size_t{0}, static_cast<size_t>(nblocks)), border.begin());

    std::ranges::sort(border, [&](const size_t a, const size_t b) { return blocks[a].number_of_pairs() > blocks[b].number_of_pairs(); });

    std::vector<TPairTask> tasks;

    for (const auto iblk : border)
    {
        const auto &block = blocks[iblk];

        if (block.number_of_pairs() == 0) continue;

        for (const auto [i, j] : orders.at({block.bra_index(), block.ket_index()}))
        {
            tasks.push_back({iblk, i, j});
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return;

    // NOTE: the coordinates of a block are formed before the tasks are drawn, as the
    // combinations of a block are computed by different threads and all of them read
    // the same coordinates.

    std::vector<CSimdMatrix> coordinates(static_cast<size_t>(nblocks));

#pragma omp parallel if (ntasks > 1)
    {
#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            coordinates[static_cast<size_t>(iblk)] = simdfunc::make_coordinates(blocks[static_cast<size_t>(iblk)], molecule);
        }

        auto arena = CSimdMatrix(arena_rows, arena_cols);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int itask = 0; itask < ntasks; itask++)
        {
            const auto &task = tasks[static_cast<size_t>(itask)];

            const auto &block = blocks[task.iblock];

            const auto npairs = block.number_of_pairs();

            const auto &bra_atoms = block.bra_atoms();

            const auto &ket_atoms = block.ket_atoms();

            const auto &a_basis = basis.basis_set(block.bra_index());

            const auto &b_basis = basis.basis_set(block.ket_index());

            const auto [la, ia] = indices[static_cast<size_t>(block.bra_index())][task.i];

            const auto [lb, jb] = indices[static_cast<size_t>(block.ket_index())][task.j];

            const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

            // NOTE: the values of a combination are held in a scratch of one row per
            // pair of angular components, as the elements of the matrix a combination
            // reaches are not contiguous.

            // NOTE: the scratch is contiguous and not a CSimdMatrix. A kernel addresses
            // the row of a component as values + m * nvalues, while the rows of a
            // CSimdMatrix are padded to a cache line, so every component past the first
            // would land in the padding of the row before it.

            std::vector<double> scratch(ncomps * npairs, 0.0);

            simdt2ceri::compute_electron_repulsion(
                scratch.data(), npairs, a_basis.functions()[task.i], b_basis.functions()[task.j], coordinates[task.iblock], buffer);

            const auto a_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

            const auto b_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

            for (size_t ma = 0; ma < a_ncomps; ma++)
            {
                for (size_t mb = 0; mb < b_ncomps; mb++)
                {
                    const auto *cell = scratch.data() + (ma * b_ncomps + mb) * npairs;

                    for (size_t k = 0; k < npairs; k++)
                    {
                        const auto row = starts[static_cast<size_t>(bra_atoms[k]) * nmoms + la] + ia + ma * strides[la];

                        const auto col = starts[static_cast<size_t>(ket_atoms[k]) * nmoms + lb] + jb + mb * strides[lb];

                        matrix.data()[matrix.index(row, col)] = cell[k];
                    }
                }
            }
        }
    }
}

auto
CSimdTwoCenterElectronRepulsionDriver::_compute_diagonal_blocks(CPackedMatrix &matrix, const CMolecularBasis &basis) const -> void
{
    // NOTE: the Coulomb operator is spherically symmetric about an atom, so the
    // integral of two basis functions on the same atom is diagonal in the angular
    // components, is the same for every component and does not depend on the position
    // of the atom. One value therefore serves every atom carrying that atom basis and
    // every component of the combination.

    // NOTE: the loop is serial. The atoms are as many as the molecule holds and each
    // combination is one scalar over the pairs of primitives, which is a fraction of
    // the atom pairs above.

    const auto indices = denseidx::index_functions(basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = strides.size();

    const auto atoms = basis.basis_sets_indices();

    for (size_t iatom = 0; iatom < atoms.size(); iatom++)
    {
        const auto &atom_basis = basis.basis_set(atoms[iatom]);

        const auto &index = indices[static_cast<size_t>(atoms[iatom])];

        for (size_t i = 0; i < index.size(); i++)
        {
            for (size_t j = 0; j < index.size(); j++)
            {
                const auto [la, ia] = index[i];

                const auto [lb, jb] = index[j];

                if (la != lb) continue;

                const auto fval = simdt2ceri::one_center_electron_repulsion(atom_basis.functions()[i], atom_basis.functions()[j]);

                const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

                for (size_t m = 0; m < ncomps; m++)
                {
                    const auto row = starts[iatom * nmoms + la] + ia + m * strides[la];

                    const auto col = starts[iatom * nmoms + lb] + jb + m * strides[lb];

                    // NOTE: the matrix is stored as its lower triangle, so only the
                    // combinations which land there are written and the reverse order
                    // is read off them.

                    if (row >= col) matrix.data()[matrix.index(row, col)] = fval;
                }
            }
        }
    }
}

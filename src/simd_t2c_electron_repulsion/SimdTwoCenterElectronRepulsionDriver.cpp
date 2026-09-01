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
#include <cstdint>
#include <ranges>
#include <utility>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "Matrix.hpp"
#include "SimdCoordinates.hpp"
#include "SimdHarmonics.hpp"
#include "SimdTwoCenterElectronRepulsionFunc.hpp"
#include "TensorComponents.hpp"

/// @brief Gets the number of angular components of an angular momentum.
/// @param angular_momentum The angular momentum.
/// @return The number of spherical components.
static auto
_number_of_components(const int angular_momentum) -> size_t
{
    return static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{angular_momentum}));
}

auto
CSimdTwoCenterElectronRepulsionDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CPackedMatrix
{
    // NOTE: the matrix is stored as its lower triangle, so an atom pair and the
    // reverse of it share their elements and only one of the two is computed, as
    // the atom basis pair groups hold only one of the two.

    auto matrix = CPackedMatrix(basis, mat_t::symmetric);

    auto groups = basis.basis_pair_groups();

    // NOTE: the atom basis pair groups are as many as the pairs of the unique
    // atom bases, so their number is set by the variety of the elements of the
    // molecule and not by its size, and the largest of them holds a third of the
    // atom pairs. Dividing them into blocks of a target number of atom pairs
    // makes the number of the blocks follow the size of the molecule instead, so
    // the work divides for any number of threads.

    const auto block_size = CAtomBasisPairGroup::make_block_size(groups, blocks_per_thread, min_block_size, max_block_size);

    // NOTE: the atom pairs are not ordered by interatomic distance, unlike those
    // of a sparse matrix. The ordering is there for the bisection of the
    // screening to read the surviving atom pairs off a leading subrange, and
    // nothing is screened here.

    const auto blocks = (block_size == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, block_size);

    _compute_pair_blocks(matrix, molecule, basis, blocks);

    _compute_diagonal_blocks(matrix, basis, blocks);

    return matrix;
}

auto
CSimdTwoCenterElectronRepulsionDriver::_compute_pair_blocks(CPackedMatrix                          &matrix,
                                                            const CMolecule                        &molecule,
                                                            const CMolecularBasis                  &basis,
                                                            const std::vector<CAtomBasisPairGroup> &blocks) const -> void
{
    // NOTE: the blocks are independent, as each of them forms its own coordinates
    // and solid harmonics and writes the elements of its own atom pairs, which no
    // other block addresses. Dynamic scheduling is used as the blocks hold a
    // comparable number of atom pairs but differ in the number of the
    // combinations of basis functions and in the cost of their kernels.

    const auto nblocks = static_cast<int>(blocks.size());

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
        const auto &block = blocks[static_cast<size_t>(iblk)];

        const auto a_indices = denseidx::index_functions(block.bra_basis());

        const auto b_indices = denseidx::index_functions(block.ket_basis());

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

    // NOTE: the dense index of an atomic orbital is looked up instead of
    // recomputed, as the innermost loops below run over the atom pairs and would
    // otherwise scan the basis for every element.

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = strides.size();

    auto *values = matrix.data();

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = blocks[order[static_cast<size_t>(iblk)]];

        const auto npairs = block.number_of_pairs();

        if (npairs == 0) continue;

        // NOTE: the coordinates of the atom pairs are created once for the whole
        // block, as all combinations of basis functions of the block share them.

        const auto coordinates = simdfunc::make_coordinates(block, molecule);

        const auto &a_basis = block.bra_basis();

        const auto &b_basis = block.ket_basis();

        // NOTE: the solid harmonics of the vectors between the atoms are created
        // once for the whole block, as all combinations of basis functions of the
        // block share them. They reach the sum of the angular momenta of the atom
        // bases of the block, as that is the highest angular momentum the
        // recursions of the integrals reach, and are empty when both atom bases
        // carry S type functions only.

        const auto lmax = a_basis.max_angular_momentum() + b_basis.max_angular_momentum();

        const auto harmonics = simdfunc::make_solid_harmonics(coordinates, lmax);

        const auto a_indices = denseidx::index_functions(a_basis);

        const auto b_indices = denseidx::index_functions(b_basis);

        // NOTE: one combination of basis functions is computed at a time into a
        // buffer of the block and added to the matrix while it is still in the
        // cache. The buffer holds the largest combination of the block, so it is
        // allocated once and the combinations reuse it.

        const auto ncomps_max = _number_of_components(a_basis.max_angular_momentum()) * _number_of_components(b_basis.max_angular_momentum());

        std::vector<double> buffer(ncomps_max * npairs, 0.0);

        const auto &bra_atoms = block.bra_atoms();

        const auto &ket_atoms = block.ket_atoms();

        // NOTE: the atom bases of a block sit on different atoms, so all
        // combinations of basis functions are computed and none of them shares
        // its elements with the reverse order.

        for (size_t i = 0; i < a_indices.size(); i++)
        {
            for (size_t j = 0; j < b_indices.size(); j++)
            {
                const auto [la, ia] = a_indices[i];

                const auto [lb, jb] = b_indices[j];

                simdt2ceri::compute_electron_repulsion(
                    buffer.data(), npairs, a_basis.functions()[i], b_basis.functions()[j], harmonics, coordinates);

                const auto bra_ncomps = _number_of_components(la);

                const auto ket_ncomps = _number_of_components(lb);

                // NOTE: the values of a pair of angular components are stored
                // contiguously over the atom pairs, with the component on bra
                // side as the slowest running index.

                for (size_t ma = 0; ma < bra_ncomps; ma++)
                {
                    for (size_t mb = 0; mb < ket_ncomps; mb++)
                    {
                        const auto *prim = buffer.data() + (ma * ket_ncomps + mb) * npairs;

                        for (size_t k = 0; k < npairs; k++)
                        {
                            const auto row = starts[static_cast<size_t>(bra_atoms[k]) * nmoms + la] + ia + ma * strides[la];

                            const auto col = starts[static_cast<size_t>(ket_atoms[k]) * nmoms + lb] + jb + mb * strides[lb];

                            // NOTE: which of the two indices is the larger one is
                            // decided element by element. The atomic orbitals are
                            // ordered by ascending angular momentum first and by
                            // atom within it, so the side of the diagonal an
                            // element falls on follows from the angular momenta
                            // of the combination as well as from the atom pair.

                            values[matrix.index(row, col)] = prim[k];
                        }
                    }
                }
            }
        }
    }
}

auto
CSimdTwoCenterElectronRepulsionDriver::_compute_diagonal_blocks(CPackedMatrix                          &matrix,
                                                                const CMolecularBasis                  &basis,
                                                                const std::vector<CAtomBasisPairGroup> &blocks) const -> void
{
    // NOTE: the repulsion of two basis functions on the same atom does not depend
    // on the position of the atom and is diagonal in the angular components, so
    // one value is computed for each pair of basis functions of the same angular
    // momentum and written to the angular components of every atom. The molecule
    // is not needed here.

    // NOTE: the elements of a diagonal atom pair which the integral does not
    // reach, i.e. those of two angular momenta which differ and those of two
    // angular components which differ, are left at the zero the matrix was
    // constructed with.

    // NOTE: this is not divided over the threads. The atoms are as many as the
    // molecule holds and the value of a combination of basis functions is
    // computed once for all of them, which is a thousandth of the work of the
    // atom pairs even for a small molecule.

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = strides.size();

    auto *values = matrix.data();

    for (const auto &block : blocks)
    {
        const auto &atoms = block.diagonal_atoms();

        if (atoms.empty()) continue;

        const auto &a_basis = block.bra_basis();

        const auto &b_basis = block.ket_basis();

        const auto a_indices = denseidx::index_functions(a_basis);

        const auto b_indices = denseidx::index_functions(b_basis);

        for (size_t i = 0; i < a_indices.size(); i++)
        {
            for (size_t j = 0; j < b_indices.size(); j++)
            {
                // NOTE: the two orders of a combination of basis functions of a
                // symmetric block share their elements, so only one is computed.

                if (block.is_symmetric() && (i > j)) continue;

                const auto [la, ia] = a_indices[i];

                const auto [lb, jb] = b_indices[j];

                if (la != lb) continue;

                const auto fval = simdt2ceri::one_center_electron_repulsion(a_basis.functions()[i], b_basis.functions()[j]);

                const auto ncomps = _number_of_components(la);

                for (const auto atom : atoms)
                {
                    const auto katom = static_cast<size_t>(atom);

                    for (size_t m = 0; m < ncomps; m++)
                    {
                        const auto row = starts[katom * nmoms + la] + ia + m * strides[la];

                        const auto col = starts[katom * nmoms + lb] + jb + m * strides[lb];

                        values[matrix.index(row, col)] = fval;
                    }
                }
            }
        }
    }
}

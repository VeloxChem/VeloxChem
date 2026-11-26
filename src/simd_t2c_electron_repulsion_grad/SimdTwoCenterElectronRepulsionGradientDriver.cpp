//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdTwoCenterElectronRepulsionGradientDriver.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include <omp.h>

#include "OpenMPFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdTwoCenterElectronRepulsionGeom10BufferRows.hpp"
#include "SimdTwoCenterElectronRepulsionGeom10Func.hpp"
#include "SimdMatrix.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

auto
CSimdTwoCenterElectronRepulsionGradientDriver::compute(const CMolecule        &molecule,
                                                       const CMolecularBasis  &basis,
                                                       const CPackedMatrix    &omega,
                                                       const std::vector<int> &atoms) const -> CPackedMatrix
{
    const auto natoms = molecule.number_of_atoms();

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdTwoCenterElectronRepulsionGradientDriver: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdTwoCenterElectronRepulsionGradientDriver: An atom is named more than once"));

    errors::assertMsgCritical(
        omega.get_type() == mat_t::symmetric,
        std::string("SimdTwoCenterElectronRepulsionGradientDriver: Omega is expected to be symmetric"));

    errors::assertMsgCritical(
        omega.number_of_rows() == basis.dimensions_of_basis(),
        std::string("SimdTwoCenterElectronRepulsionGradientDriver: Omega is not of the auxiliary basis"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    if (atoms.empty()) return gradient;

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms) wanted[static_cast<size_t>(iatom)] = true;

    // NOTE: only the geometry half of the pattern is formed, as the driver of the
    // integrals themselves does. Nothing is described and nothing is screened: no
    // atom pair of the Coulomb operator falls below a threshold, and its
    // derivative falls off one power faster without becoming sparse.

    auto all = basis.basis_pair_groups();

    // NOTE: the pairs which touch none of the atoms asked for are dropped here,
    // before the blocks are formed. A kernel computes every atom pair of a block
    // in one call, so a pair dropped later would have been computed anyway; and
    // the blocks are what the threads draw on, so dropping first is also what
    // gives them the work which is actually wanted.

    auto groups = std::vector<CAtomBasisPairGroup>();

    groups.reserve(all.size());

    for (const auto &group : all)
    {
        auto kept = group.select_atoms(wanted);

        if (kept.number_of_pairs() > 0) groups.push_back(std::move(kept));
    }

    if (groups.empty()) return gradient;

    const auto nblock_pairs = (_block_size == 0)
                                  ? CAtomBasisPairGroup::make_block_size(
                                        groups, sparsity::blocks_per_thread, sparsity::min_block_size, max_block_size)
                                  : _block_size;

    auto blocks = (nblock_pairs == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, nblock_pairs);

    CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

    _compute_pair_blocks(gradient, molecule, basis, omega, blocks, wanted);

    // NOTE: no diagonal blocks. The two centers of a pair of basis functions on
    // one atom move together, so the integral does not change and its derivative
    // is zero. The driver of the integrals themselves has such a pass and this
    // one has nothing to put in it.

    return gradient;
}

auto
CSimdTwoCenterElectronRepulsionGradientDriver::compute(const CMolecule       &molecule,
                                                       const CMolecularBasis &basis,
                                                       const CPackedMatrix   &omega) const -> CPackedMatrix
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute(molecule, basis, omega, atoms);
}

auto
CSimdTwoCenterElectronRepulsionGradientDriver::_compute_pair_blocks(CPackedMatrix                          &gradient,
                                                                    const CMolecule                        &molecule,
                                                                    const CMolecularBasis                  &basis,
                                                                    const CPackedMatrix                    &omega,
                                                                    const std::vector<CAtomBasisPairGroup> &blocks,
                                                                    const std::vector<bool>                &wanted) const
    -> void
{
    const auto indices = denseidx::index_functions(basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto nmoms = strides.size();

    const auto nblocks = static_cast<int>(blocks.size());

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = blocks[static_cast<size_t>(iblk)];

        arena_rows = std::max(arena_rows,
                              simdt2ceri::number_of_geom_10_buffer_rows(basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                                    basis.basis_set(block.ket_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    std::vector<TPairTask> tasks;

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = blocks[static_cast<size_t>(iblk)];

        if (block.number_of_pairs() == 0) continue;

        const auto &a_index = indices[static_cast<size_t>(block.bra_index())];

        const auto &b_index = indices[static_cast<size_t>(block.ket_index())];

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                tasks.push_back({static_cast<size_t>(iblk), i, j});
            }
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return;

    std::vector<CSimdMatrix> coordinates(static_cast<size_t>(nblocks));

    const auto natoms = gradient.number_of_rows();

    // NOTE: two atom pairs share an atom whenever they share a center, so the
    // threads do not write to disjoint rows the way the driver of the integrals
    // does. Each holds its own gradient and they are summed at the end, which is
    // three doubles an atom a thread.

    const auto nthreads = omp::get_number_of_threads();

    auto partials = std::vector<std::vector<double>>(static_cast<size_t>(nthreads),
                                                     std::vector<double>(natoms * 3, 0.0));

#pragma omp parallel if (ntasks > 1)
    {
#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblocks; iblk++)
        {
            coordinates[static_cast<size_t>(iblk)] = simdfunc::make_coordinates(blocks[static_cast<size_t>(iblk)], molecule);
        }

        auto arena = CSimdMatrix(arena_rows, arena_cols);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

        auto &partial = partials[static_cast<size_t>(omp_get_thread_num())];

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

            // NOTE: three components of the derivative for every pair of angular
            // components, laid out as the kernels of the integrals lay out one.

            std::vector<double> scratch(3 * ncomps * npairs, 0.0);

            simdt2cerigrad::compute_electron_repulsion_geom_10(
                scratch.data(), npairs, a_basis.functions()[task.i], b_basis.functions()[task.j],
                coordinates[task.iblock], buffer);

            const auto a_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

            const auto b_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

            for (size_t k = 0; k < npairs; k++)
            {
                const auto iatom = static_cast<size_t>(bra_atoms[k]);

                const auto jatom = static_cast<size_t>(ket_atoms[k]);

                std::array<double, 3> taken{0.0, 0.0, 0.0};

                for (size_t ma = 0; ma < a_ncomps; ma++)
                {
                    const auto row = starts[iatom * nmoms + la] + ia + ma * strides[la];

                    for (size_t mb = 0; mb < b_ncomps; mb++)
                    {
                        const auto col = starts[jatom * nmoms + lb] + jb + mb * strides[lb];

                        // NOTE: the pattern holds the strict upper triangle of the
                        // atom pairs, so the pair with its two atoms the other way
                        // round is not swept. Omega and the integral are both
                        // symmetric and its derivative with respect to this atom is
                        // the same either way, so the pair stands for two terms of
                        // the sum and is counted twice.

                        const auto weight = 2.0 * omega.at(row, col);

                        for (size_t c = 0; c < 3; c++)
                        {
                            taken[c] += weight * scratch[(c * ncomps + ma * b_ncomps + mb) * npairs + k];
                        }
                    }
                }

                // NOTE: the derivative of the atom on ket side is the negative of
                // the one on bra side, which is why the pair is swept once.

                for (size_t c = 0; c < 3; c++)
                {
                    if (wanted[iatom]) partial[iatom * 3 + c] += taken[c];

                    if (wanted[jatom]) partial[jatom * 3 + c] -= taken[c];
                }
            }
        }
    }

    for (const auto &partial : partials)
    {
        for (size_t iatom = 0; iatom < natoms; iatom++)
        {
            for (size_t c = 0; c < 3; c++)
            {
                gradient.data()[gradient.index(iatom, c)] += partial[iatom * 3 + c];
            }
        }
    }
}

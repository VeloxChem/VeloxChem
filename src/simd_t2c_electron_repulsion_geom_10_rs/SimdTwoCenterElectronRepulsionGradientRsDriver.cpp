//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdTwoCenterElectronRepulsionGradientRsDriver.hpp"

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
#include "SimdTwoCenterElectronRepulsionGeom10RsBufferRows.hpp"
#include "SimdTwoCenterElectronRepulsionGeom10RsFunc.hpp"
#include "SimdMatrix.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

auto
CSimdTwoCenterElectronRepulsionGradientRsDriver::compute(const CMolecule        &molecule,
                                                       const CMolecularBasis  &basis,
                                                       const CPackedMatrix    &omega_coulomb,
                                                       const CPackedMatrix    &omega_attenuated,
                                                       const double            omega,
                                                       const std::vector<int> &atoms) const
    -> std::pair<CPackedMatrix, CPackedMatrix>
{
    errors::assertMsgCritical(
        omega > 0.0,
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: The range separation parameter must be positive"));

    const auto natoms = molecule.number_of_atoms();

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical(
            (iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
            std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: Atom is not an atom of the molecule"));
    }

    errors::assertMsgCritical(
        std::set<int>(atoms.begin(), atoms.end()).size() == atoms.size(),
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: An atom is named more than once"));

    errors::assertMsgCritical(
        omega_coulomb.get_type() == mat_t::symmetric,
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: Omega is expected to be symmetric"));

    errors::assertMsgCritical(
        omega_coulomb.number_of_rows() == basis.dimensions_of_basis(),
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: Omega is not of the auxiliary basis"));

    errors::assertMsgCritical(
        omega_attenuated.get_type() == mat_t::symmetric,
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: The attenuated Omega is expected to be "
                    "symmetric"));

    errors::assertMsgCritical(
        omega_attenuated.number_of_rows() == basis.dimensions_of_basis(),
        std::string("SimdTwoCenterElectronRepulsionGradientRsDriver: The attenuated Omega is not of the auxiliary "
                    "basis"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    auto gradient_attenuated = gradient;

    if (atoms.empty()) return {std::move(gradient), std::move(gradient_attenuated)};

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

    if (groups.empty()) return {std::move(gradient), std::move(gradient_attenuated)};

    const auto nblock_pairs = (_block_size == 0)
                                  ? CAtomBasisPairGroup::make_block_size(
                                        groups, sparsity::blocks_per_thread, sparsity::min_block_size, max_block_size)
                                  : _block_size;

    auto blocks = (nblock_pairs == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, nblock_pairs);

    CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

    _compute_pair_blocks(gradient, gradient_attenuated, molecule, basis, omega_coulomb, omega_attenuated,
                         omega, blocks, wanted);

    // NOTE: no diagonal blocks. The two centers of a pair of basis functions on
    // one atom move together, so the integral does not change and its derivative
    // is zero. The driver of the integrals themselves has such a pass and this
    // one has nothing to put in it.

    return {std::move(gradient), std::move(gradient_attenuated)};
}

auto
CSimdTwoCenterElectronRepulsionGradientRsDriver::compute(const CMolecule       &molecule,
                                                       const CMolecularBasis &basis,
                                                       const CPackedMatrix   &omega_coulomb,
                                                       const CPackedMatrix   &omega_attenuated,
                                                       const double           omega) const
    -> std::pair<CPackedMatrix, CPackedMatrix>
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute(molecule, basis, omega_coulomb, omega_attenuated, omega, atoms);
}

auto
CSimdTwoCenterElectronRepulsionGradientRsDriver::_compute_pair_blocks(
    CPackedMatrix                          &gradient_coulomb,
    CPackedMatrix                          &gradient_attenuated,
    const CMolecule                        &molecule,
    const CMolecularBasis                  &basis,
    const CPackedMatrix                    &omega_coulomb,
    const CPackedMatrix                    &omega_attenuated,
    const double                            omega,
    const std::vector<CAtomBasisPairGroup> &blocks,
    const std::vector<bool>                &wanted) const -> void
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
                              simdt2ceri::number_of_geom_10_rs_buffer_rows(basis.basis_set(block.bra_index()).max_angular_momentum(),
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

    const auto natoms = gradient_coulomb.number_of_rows();

    // NOTE: two atom pairs share an atom whenever they share a center, so the
    // threads do not write to disjoint rows the way the driver of the integrals
    // does. Each holds its own gradient and they are summed at the end, which is
    // three doubles an atom a thread.

    const auto nthreads = omp::get_number_of_threads();

    // NOTE: a partial for each operator as well as for each thread. The two are
    // accumulated in one sweep because they come out of one kernel call; keeping
    // them apart costs three doubles an atom a thread more than the unattenuated
    // driver and is what lets the caller weight them differently.

    auto partials = std::vector<std::vector<double>>(static_cast<size_t>(nthreads),
                                                     std::vector<double>(natoms * 3, 0.0));

    auto partials_attenuated = partials;

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

        auto &partial_attenuated = partials_attenuated[static_cast<size_t>(omp_get_thread_num())];

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

            // NOTE: **six** blocks and not three: three Cartesian components of
            // each of the two operators, laid out one operator after the other.

            std::vector<double> scratch(6 * ncomps * npairs, 0.0);

            simdt2cerigrad::compute_rs_electron_repulsion_geom_10(
                scratch.data(), npairs, a_basis.functions()[task.i], b_basis.functions()[task.j],
                coordinates[task.iblock], buffer, omega);

            const auto a_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

            const auto b_ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

            for (size_t k = 0; k < npairs; k++)
            {
                const auto iatom = static_cast<size_t>(bra_atoms[k]);

                const auto jatom = static_cast<size_t>(ket_atoms[k]);

                std::array<double, 3> taken{0.0, 0.0, 0.0};

                std::array<double, 3> taken_attenuated{0.0, 0.0, 0.0};

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

                        const auto weight = 2.0 * omega_coulomb.at(row, col);

                        const auto weight_attenuated = 2.0 * omega_attenuated.at(row, col);

                        // NOTE: each operator is contracted with its own weighted
                        // density. They are fitted in different metrics and the two
                        // matrices are not interchangeable.

                        for (size_t c = 0; c < 3; c++)
                        {
                            const auto element = ma * b_ncomps + mb;

                            taken[c] += weight *
                                        scratch[((_coulomb_block + c) * ncomps + element) * npairs + k];

                            taken_attenuated[c] +=
                                weight_attenuated *
                                scratch[((_attenuated_block + c) * ncomps + element) * npairs + k];
                        }
                    }
                }

                // NOTE: the derivative of the atom on ket side is the negative of
                // the one on bra side, which is why the pair is swept once.

                for (size_t c = 0; c < 3; c++)
                {
                    if (wanted[iatom])
                    {
                        partial[iatom * 3 + c] += taken[c];

                        partial_attenuated[iatom * 3 + c] += taken_attenuated[c];
                    }

                    if (wanted[jatom])
                    {
                        partial[jatom * 3 + c] -= taken[c];

                        partial_attenuated[jatom * 3 + c] -= taken_attenuated[c];
                    }
                }
            }
        }
    }

    for (size_t ithread = 0; ithread < partials.size(); ithread++)
    {
        const auto &partial = partials[ithread];

        const auto &attenuated = partials_attenuated[ithread];

        for (size_t iatom = 0; iatom < natoms; iatom++)
        {
            for (size_t c = 0; c < 3; c++)
            {
                gradient_coulomb.data()[gradient_coulomb.index(iatom, c)] += partial[iatom * 3 + c];

                gradient_attenuated.data()[gradient_attenuated.index(iatom, c)] += attenuated[iatom * 3 + c];
            }
        }
    }
}

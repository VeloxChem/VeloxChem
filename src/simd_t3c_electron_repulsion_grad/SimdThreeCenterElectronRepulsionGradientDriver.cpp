//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"

#include <algorithm>
#include <array>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "SimdCoordinates.hpp"
#include "SimdMatrix.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100Func.hpp"
#include "TensorComponents.hpp"

namespace {

/// @brief One combination of basis functions of one block, the unit of work.
struct TTripleTask
{
    size_t iblock;

    size_t i;

    size_t j;

    size_t k;

    size_t npairs;
};

/// @brief The coordinates of the atoms on the auxiliary side of a block.
auto
_make_c_coordinates(const CAtomBasisTripleSparsity &block, const CMolecule &molecule) -> CSimdMatrix
{
    const auto &coords = molecule.coordinates();

    const auto &atoms = block.c_atoms();

    auto matrix = CSimdMatrix(3, atoms.size());

    for (size_t i = 0; i < atoms.size(); i++)
    {
        const auto r = coords[static_cast<size_t>(atoms[i])].coordinates();

        for (size_t c = 0; c < 3; c++) matrix.data(c)[i] = r[c];
    }

    return matrix;
}

}  // namespace

auto
CSimdThreeCenterElectronRepulsionGradientDriver::compute(const CTripleSparsityPattern &pattern,
                                                         const CMolecule              &molecule,
                                                         const CMolecularBasis        &basis,
                                                         const CMolecularBasis        &aux_basis) const -> CSparseTensor
{
    const auto indices = denseidx::index_functions(basis);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    sparsity::check_triple_pattern(pattern, indices, aux_indices);

    auto tensor = CSparseTensor(pattern);

    // NOTE: six values an element, the three directions of each of the two atoms
    // on bra side. The components are set before the values are allocated, as the
    // size of a block follows them.

    tensor.set_number_of_components(simdt3cerigrad::number_of_components);

    tensor.allocate();

    tensor.zero();

    const auto nblocks = static_cast<size_t>(pattern.number_of_blocks());

    if (nblocks == 0) return tensor;

    // NOTE: the arena spans the largest combination any block carries. The
    // derivative raises the momentum of the differentiated center by one and
    // writes six components, so it is larger than the integral's and may force a
    // smaller block than the energy driver uses. That is a measurement to make
    // once the kernels exist.

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.block(iblk);

        if ((block.number_of_pairs() == 0) || (block.number_of_c_atoms() == 0)) continue;

        arena_rows = std::max(arena_rows,
                              simdt3cerigrad::number_of_buffer_rows(basis.basis_set(block.a_index()).max_angular_momentum(),
                                                                    basis.basis_set(block.b_index()).max_angular_momentum(),
                                                                    aux_basis.basis_set(block.c_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs() * block.number_of_c_atoms());
    }

    std::vector<TTripleTask> tasks;

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.block(iblk);

        if ((block.number_of_pairs() == 0) || (block.number_of_c_atoms() == 0)) continue;

        const auto &a_index = indices[static_cast<size_t>(block.a_index())];

        const auto &b_index = indices[static_cast<size_t>(block.b_index())];

        const auto &c_index = aux_indices[static_cast<size_t>(block.c_index())];

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                for (size_t k = 0; k < c_index.size(); k++)
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
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return tensor;

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

            // NOTE: the offsets of the tensor scale by the components, so this is
            // the start of the six runs of the combination.

            auto *values = tensor.values(task.iblock, la, ia, lb, jb, lc, kc);

            simdt3cerigrad::compute_electron_repulsion_geom_100(values,
                                                               task.npairs,
                                                               natoms,
                                                               a_basis.functions()[task.i],
                                                               b_basis.functions()[task.j],
                                                               c_basis.functions()[task.k],
                                                               coordinates[task.iblock],
                                                               c_coordinates[task.iblock],
                                                               buffer,
                                                               pattern.get_threshold());
        }
    }

    return tensor;
}

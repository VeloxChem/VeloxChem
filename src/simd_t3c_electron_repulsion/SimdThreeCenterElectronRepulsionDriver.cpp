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



#include "SimdThreeCenterElectronRepulsionDriver.hpp"

#include <algorithm>
#include <cstdint>
#include <ranges>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdHarmonics.hpp"
#include "SimdThreeCenterElectronRepulsionFunc.hpp"

auto
CSimdThreeCenterElectronRepulsionDriver::compute(const CMolecule       &molecule,
                                                 const CMolecularBasis &basis,
                                                 const CMolecularBasis &aux_basis,
                                                 const double           threshold) const -> CSparseTensor
{
    auto tensor = CSparseTensor(molecule, basis, aux_basis, screenfunc::three_center_electron_repulsion_bound, threshold, mat_t::symmetric);

    // NOTE: the values blocks are not set to zero after they are allocated, as
    // every value of every block is written below. A combination of basis
    // functions which no atom pair survives holds no values, and a kernel writes
    // the integrals of the atom pairs it reaches, so no value is left with the
    // undefined content of the allocation.

    tensor.allocate();

    _compute_blocks(tensor, molecule, basis, aux_basis);

    return tensor;
}

auto
CSimdThreeCenterElectronRepulsionDriver::compute(const CMolecule        &molecule,
                                                 const CMolecularBasis  &basis,
                                                 const CMolecularBasis  &aux_basis,
                                                 const double            threshold,
                                                 const std::vector<int> &atoms) const -> CSparseTensor
{
    auto tensor =
        CSparseTensor(molecule, basis, aux_basis, screenfunc::three_center_electron_repulsion_bound, threshold, mat_t::symmetric, atoms);

    tensor.allocate();

    _compute_blocks(tensor, molecule, basis, aux_basis);

    return tensor;
}

auto
CSimdThreeCenterElectronRepulsionDriver::_compute_blocks(CSparseTensor         &tensor,
                                                         const CMolecule       &molecule,
                                                         const CMolecularBasis &basis,
                                                         const CMolecularBasis &aux_basis) const -> void
{
    // NOTE: the blocks are independent, as each of them forms its own
    // coordinates and solid harmonics and writes the values of its own block,
    // which no other block addresses. Dynamic scheduling is used as the blocks
    // differ both in the atom pairs they hold and in the cost of their kernels.

    const auto nblocks = static_cast<int>(tensor.number_of_blocks());

    // NOTE: the blocks are visited from the most costly to the least, so that a
    // costly block is taken while there is still work to fill the other threads
    // with. A costly block drawn last is finished alone and sets the time of the
    // whole loop.

    // NOTE: the cost of a block is the atom pairs it holds times the atoms on c
    // side times the cost of its combinations of basis functions, each weighted
    // by the sum of the three angular momenta, as the recursions of the
    // integrals grow with that sum. The weight is an estimate and orders the
    // blocks, it is not used for anything else.

    std::vector<size_t> order(static_cast<size_t>(nblocks));

    std::vector<double> costs(static_cast<size_t>(nblocks), 0.0);

    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = tensor.block(static_cast<size_t>(iblk));

        const auto a_indices = denseidx::index_functions(basis.basis_set(block.a_index()));

        const auto b_indices = denseidx::index_functions(basis.basis_set(block.b_index()));

        const auto c_indices = denseidx::index_functions(aux_basis.basis_set(block.c_index()));

        double weight = 0.0;

        for (const auto &[la, ia] : a_indices)
        {
            for (const auto &[lb, jb] : b_indices)
            {
                for (const auto &[lc, kc] : c_indices)
                {
                    const auto lsum = static_cast<double>(la + lb + lc + 1);

                    weight += lsum * lsum;
                }
            }
        }

        order[static_cast<size_t>(iblk)] = static_cast<size_t>(iblk);

        costs[static_cast<size_t>(iblk)] =
            static_cast<double>(block.number_of_pairs()) * static_cast<double>(block.number_of_c_atoms()) * weight;
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto jblk = order[static_cast<size_t>(iblk)];

        const auto &block = tensor.block(jblk);

        if (block.number_of_pairs() == 0) continue;

        // NOTE: the coordinates of the atom pairs and of the atoms on c side are
        // created once for the whole block, as all combinations of basis
        // functions of the block share them.

        const auto coordinates = simdfunc::make_coordinates(block, molecule);

        const auto c_coordinates = _make_c_coordinates(block, molecule);

        const auto &a_basis = basis.basis_set(block.a_index());

        const auto &b_basis = basis.basis_set(block.b_index());

        const auto &c_basis = aux_basis.basis_set(block.c_index());

        // NOTE: the solid harmonics of the vectors between the atoms of the atom
        // pairs are created once for the whole block. They reach the sum of the
        // angular momenta of the atom bases on the a and b sides, as that is the
        // highest angular momentum the recursions over the atom pairs reach.

        const auto lmax = a_basis.max_angular_momentum() + b_basis.max_angular_momentum();

        const auto harmonics = simdfunc::make_solid_harmonics(coordinates, lmax);

        const auto a_indices = denseidx::index_functions(a_basis);

        const auto b_indices = denseidx::index_functions(b_basis);

        const auto c_indices = denseidx::index_functions(c_basis);

        const auto natoms = block.number_of_c_atoms();

        // NOTE: the combinations of basis functions are independent, as each of
        // them writes its own values and reads the coordinates and the harmonics
        // of the block without changing them. All combinations are computed, the
        // reverse order of the a and b sides included, as the sparsity pattern
        // gives each of them its own values.

        for (size_t i = 0; i < a_indices.size(); i++)
        {
            for (size_t j = 0; j < b_indices.size(); j++)
            {
                for (size_t k = 0; k < c_indices.size(); k++)
                {
                    const auto [la, ia] = a_indices[i];

                    const auto [lb, jb] = b_indices[j];

                    const auto [lc, kc] = c_indices[k];

                    // NOTE: the atom pairs which survive the screening differ
                    // from one combination of basis functions to another, as the
                    // bound carries the angular momenta and the exponents of the
                    // three basis functions.

                    const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                    if (npairs == 0) continue;

                    simdt3ceri::compute_electron_repulsion(tensor.values(jblk, la, ia, lb, jb, lc, kc),
                                                           npairs,
                                                           natoms,
                                                           a_basis.functions()[i],
                                                           b_basis.functions()[j],
                                                           c_basis.functions()[k],
                                                           harmonics,
                                                           coordinates,
                                                           c_coordinates);
                }
            }
        }
    }
}

auto
CSimdThreeCenterElectronRepulsionDriver::_make_c_coordinates(const CAtomBasisTripleSparsity &block, const CMolecule &molecule) const
    -> CSimdMatrix
{
    // NOTE: the atoms on c side are a separate and shorter dimension than the
    // atom pairs, so their coordinates are held as three rows of as many columns
    // as there are atoms, unlike the seven rows of the atom pairs.

    const auto &atoms = block.c_atoms();

    const auto natoms = atoms.size();

    const auto &coords = molecule.coordinates();

    errors::assertMsgCritical(std::ranges::none_of(atoms, [&](const auto atom) { return atom >= static_cast<int>(coords.size()); }),
                              std::string("SimdThreeCenterElectronRepulsionDriver: Atomic index out of range of molecule"));

    auto matrix = CSimdMatrix(3, natoms);

    // NOTE: the coordinates of an axis are stored contiguously over the atoms on
    // c side, as they are over the atom pairs, so that a recursion loads them
    // with aligned SIMD instructions.

    auto *c_x = matrix.data(0);

    auto *c_y = matrix.data(1);

    auto *c_z = matrix.data(2);

    const auto *rxyz = coords.data();

    for (size_t i = 0; i < natoms; i++)
    {
        const auto r_c = rxyz[atoms[i]].coordinates();

        c_x[i] = r_c[0];

        c_y[i] = r_c[1];

        c_z[i] = r_c[2];
    }

    return matrix;
}

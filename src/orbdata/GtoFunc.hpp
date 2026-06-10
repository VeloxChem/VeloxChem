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

#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <algorithm>
#include <vector>

#include "AtomBasis.hpp"
#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "ScreenedBasisFunctionPair.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis functions blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule) -> std::vector<CGtoBlock>;

/// @brief Creates vector of basis functions blocks for selected atoms.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param atoms The vector of atoms to select.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms) -> std::vector<CGtoBlock>;

/**
 Gets number of atomic orbitals from vector of contracted GTOs blocks.
 @param gto_blocks the vector of contracted GTOs blocks.
 @return the number of atomic orbitals.
 */
auto getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int;

/// @brief Gets Cartesian coordinates (in au) of atoms associated with a
/// specific atom basis in molecule.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param atom_basis The atom basis to select atoms for; matched by basis name
/// and chemical element identifier.
/// @return The vector of atom Cartesian coordinates.
auto atom_coordinates(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomBasis &atom_basis) -> std::vector<TPoint<double>>;

/// @brief Gets indices of atoms associated with a specific atom basis in
/// molecular basis. The atoms are listed in the same order as the coordinates
/// returned by atom_coordinates for the same atom basis.
/// @param basis The molecular basis.
/// @param atom_basis The atom basis to select atoms for; matched by basis name
/// and chemical element identifier.
/// @return The vector of atom indices.
auto atom_indices(const CMolecularBasis &basis, const CAtomBasis &atom_basis) -> std::vector<int>;

/// @brief Creates vector of screened basis function pairs for a molecular basis
/// and molecule. The unique symmetric pairs of atom bases are enumerated; for
/// each atom basis pair a screened basis function pair is built for every basis
/// function pair combination (a triangular loop over basis functions, with the
/// symmetric single-function constructor on the diagonal, when the two atom
/// bases are identical; a full loop otherwise). The per combination
/// constructions are independent and run in parallel with OpenMP; pairs with no
/// surviving atom pair are dropped.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param make_keep_pair The predicate factory; called as
/// make_keep_pair(bra_function, ket_function) and returning a center predicate
/// keep_pair(bra_center, ket_center) -> bool selecting important atom pairs. A
/// fresh predicate is built for each basis function pair combination, so the
/// exponent dependent data can be precomputed once per combination. The factory
/// must be safe to call concurrently from multiple threads.
/// @return The vector of screened basis function pairs.
template <typename PredicateFactory>
auto
make_screened_basis_function_pairs(const CMolecularBasis &basis, const CMolecule &molecule, PredicateFactory &&make_keep_pair)
    -> std::vector<CScreenedBasisFunctionPair>
{
    const auto basis_sets = basis.basis_sets();

    const auto nsets = basis_sets.size();

    // precompute per unique atom basis: atom coordinates, atom indices and basis functions

    std::vector<std::vector<TPoint<double>>> coordinates(nsets);

    std::vector<std::vector<int>> atomic_indices(nsets);

    std::vector<std::vector<CBasisFunction>> functions(nsets);

    for (size_t k = 0; k < nsets; k++)
    {
        coordinates[k] = atom_coordinates(basis, molecule, basis_sets[k]);

        atomic_indices[k] = atom_indices(basis, basis_sets[k]);

        functions[k] = basis_sets[k].basis_functions();
    }

    // build flat list of basis function pair tasks over unique symmetric atom basis pairs

    struct CScreenedPairTask
    {
        size_t bra_slot;
        size_t bra_func;
        size_t ket_slot;
        size_t ket_func;
    };

    std::vector<CScreenedPairTask> tasks;

    for (size_t bi = 0; bi < nsets; bi++)
    {
        const auto nbra = functions[bi].size();

        for (size_t bj = bi; bj < nsets; bj++)
        {
            if (bi == bj)
            {
                for (size_t fi = 0; fi < nbra; fi++)
                    for (size_t fj = fi; fj < nbra; fj++) tasks.push_back({bi, fi, bi, fj});
            }
            else
            {
                const auto nket = functions[bj].size();

                for (size_t fi = 0; fi < nbra; fi++)
                    for (size_t fj = 0; fj < nket; fj++) tasks.push_back({bi, fi, bj, fj});
            }
        }
    }

    const auto ntasks = tasks.size();

    // schedule larger atom-pair work first for better dynamic load balance

    std::vector<size_t> schedule(ntasks);

    for (size_t t = 0; t < ntasks; t++) schedule[t] = t;

    std::ranges::sort(schedule, [&](const size_t a, const size_t b) {
        return coordinates[tasks[a].bra_slot].size() * coordinates[tasks[a].ket_slot].size() >
               coordinates[tasks[b].bra_slot].size() * coordinates[tasks[b].ket_slot].size();
    });

    // construct screened basis function pairs in parallel into per-task slots

    std::vector<CScreenedBasisFunctionPair> results(ntasks);

    const auto nsched = static_cast<long>(ntasks);

#pragma omp parallel for schedule(dynamic)
    for (long s = 0; s < nsched; s++)
    {
        const auto  t   = schedule[static_cast<size_t>(s)];
        const auto &task = tasks[t];

        const auto &bra_function = functions[task.bra_slot][task.bra_func];
        const auto &ket_function = functions[task.ket_slot][task.ket_func];

        auto keep_pair = make_keep_pair(bra_function, ket_function);

        if ((task.bra_slot == task.ket_slot) && (task.bra_func == task.ket_func))
        {
            // identical function on the identical atom set: triangular atom pairs

            results[t] = CScreenedBasisFunctionPair(
                bra_function, static_cast<int>(task.bra_func), coordinates[task.bra_slot], atomic_indices[task.bra_slot], keep_pair);
        }
        else
        {
            // distinct functions and/or atom sets: full atom product

            results[t] = CScreenedBasisFunctionPair(bra_function,
                                                    static_cast<int>(task.bra_func),
                                                    coordinates[task.bra_slot],
                                                    atomic_indices[task.bra_slot],
                                                    ket_function,
                                                    static_cast<int>(task.ket_func),
                                                    coordinates[task.ket_slot],
                                                    atomic_indices[task.ket_slot],
                                                    keep_pair);
        }
    }

    // keep only pairs with at least one surviving atom pair, in task order

    std::vector<CScreenedBasisFunctionPair> screened_pairs;

    screened_pairs.reserve(ntasks);

    for (size_t t = 0; t < ntasks; t++)
    {
        if (results[t].number_of_pairs() > 0) screened_pairs.push_back(std::move(results[t]));
    }

    return screened_pairs;
}

}  // namespace gtofunc

#endif /* GtoFunc_hpp */

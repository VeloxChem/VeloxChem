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
#include <vector>

#include "AtomBasisTripleSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdMatrix.hpp"
#include "SimdT3CDistributor.hpp"
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

    const auto nblocks = static_cast<int>(pattern.number_of_blocks());

    // NOTE: the blocks are independent, as each of them forms its own coordinates and
    // hands the distributor the values of its own combinations of basis functions,
    // which no other block addresses. Dynamic scheduling is used as the blocks differ
    // in the number of the combinations and in the atoms they carry on c side.

#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
    for (int iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.block(static_cast<size_t>(iblk));

        const auto natoms = block.number_of_c_atoms();

        if ((block.number_of_pairs() == 0) || (natoms == 0)) continue;

        // NOTE: the coordinates of the atom pairs and of the atoms on c side are
        // created once for the whole block, as all combinations of basis functions of
        // the block share them.

        const auto coordinates = simdfunc::make_coordinates(block, molecule);

        const auto c_coordinates = _make_c_coordinates(block, molecule);

        const auto &a_basis = basis.basis_set(block.a_index());

        const auto &b_basis = basis.basis_set(block.b_index());

        const auto &c_basis = aux_basis.basis_set(block.c_index());

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

                    const auto npairs = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                    if (npairs == 0) continue;

                    const auto ncomps =
                        static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 3>{la, lb, lc}));

                    auto *values = distributor.target(
                        block, static_cast<size_t>(iblk), la, ia, lb, jb, lc, kc, npairs, natoms, ncomps);

                    simdt3ceri::compute_electron_repulsion(values,
                                                           npairs,
                                                           natoms,
                                                           a_basis.functions()[i],
                                                           b_basis.functions()[j],
                                                           c_basis.functions()[k],
                                                           coordinates,
                                                           c_coordinates);

                    distributor.commit(
                        block, static_cast<size_t>(iblk), la, ia, lb, jb, lc, kc, npairs, natoms, ncomps);
                }
            }
        }
    }
}

#endif /* SimdThreeCenterElectronRepulsionDriver_hpp */

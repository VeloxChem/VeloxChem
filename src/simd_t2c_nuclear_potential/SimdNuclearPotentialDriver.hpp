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


#ifndef SimdNuclearPotentialDriver_hpp
#define SimdNuclearPotentialDriver_hpp

#include <algorithm>
#include <cstddef>
#include <map>
#include <ranges>
#include <vector>

#include "DenseIndexFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"
#include "SimdCoordinates.hpp"
#include "SimdNuclearPotentialBufferRows.hpp"
#include "SimdNuclearPotentialFunc.hpp"
#include "SimdT2CDistributor.hpp"
#include "SparseMatrix.hpp"
#include "SparsityPattern.hpp"
#include "TensorComponents.hpp"

/// @brief Class CSimdNuclearPotentialDriver computes the two-center nuclear
/// potential integrals of a molecular basis and a set of point charges, and stores
/// them in a sparse matrix, using the sparsity patterns of the atom basis pair
/// groups to skip the atom pairs and the combinations of basis functions whose
/// integrals are below the threshold.
///
/// @note The operator is a sum over point charges of the Coulomb attraction of an
/// electron to each of them, so the integrals carry the charges and their positions
/// beside the basis. Two forms are offered: one which takes them, for a calculation
/// whose charges are not the nuclei -- an embedding, a set of external charges, a
/// subset of the molecule -- and one which takes them from the molecule.
///
/// @note The screening bound of the operator depends on the two basis functions and
/// the distance between their atoms and not on where the charges are, so a sparsity
/// pattern serves any set of charges and is formed once for a molecule and a basis.
///
/// @note Unlike the overlap and the kinetic energy, there is no closed form for two
/// basis functions on the same atom: the operator is centred on the charges rather
/// than on the atom, so the diagonal blocks carry the same kernels as the others and
/// are not a cheap case.
class CSimdNuclearPotentialDriver
{
    /// @brief One combination of basis functions of one block, which is the unit of
    /// work the threads draw on.
    struct TPairTask
    {
        /// @brief The index of the block among the off-diagonal blocks.
        size_t iblock;

        /// @brief The index of the basis function on bra side within its atom basis.
        size_t i;

        /// @brief The index of the basis function on ket side within its atom basis.
        size_t j;

        /// @brief The number of atom pairs the combination reaches.
        size_t nvalues;
    };

   public:
    /// @brief The constructor with screening threshold and target block size.
    /// @param threshold The screening threshold of the integrals.
    /// @param block_size The target number of atom pairs of a block, or zero to
    /// choose it from the number of the threads and the number of the atom pairs.
    explicit CSimdNuclearPotentialDriver(const double threshold = 1.0e-14, const size_t block_size = 0)

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

    /// @brief Gets screening threshold of the integrals.
    /// @return The screening threshold.
    auto
    get_threshold() const -> double
    {
        return _threshold;
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
        return sparsity::make_pattern(molecule,
                                      bra_basis,
                                      ket_basis,
                                      screener::nuclear_potential,
                                      _threshold,
                                      mat_type,
                                      diagstor::full,
                                      _block_size,
                                      _min_block_pairs);
    }

    /// @brief Computes the nuclear potential integrals of a sparsity pattern and
    /// hands them to a distributor.
    /// @param pattern The sparsity pattern to compute the integrals of.
    /// @param molecule The molecule to compute the integrals of.
    /// @param bra_basis The molecular basis on bra side.
    /// @param ket_basis The molecular basis on ket side.
    /// @param charges The magnitude of each point charge.
    /// @param points The position of each point charge, as a flat array of three
    /// values per charge.
    /// @param distributor The distributor to hand the integrals to.
    /// @note This is the form which does not name a container of the values. The
    /// overloads which return a sparse matrix are written in terms of it.
    template <class D>
    auto
    compute(const CSparsityPattern    &pattern,
            const CMolecule           &molecule,
            const CMolecularBasis     &bra_basis,
            const CMolecularBasis     &ket_basis,
            const std::vector<double> &charges,
            const std::vector<double> &points,
            D                         &distributor) const -> void
    {
        errors::assertMsgCritical(points.size() == 3 * charges.size(),
                                  std::string("SimdNuclearPotentialDriver.compute: Expecting three coordinates for each charge"));

        const auto a_indices = denseidx::index_functions(bra_basis);

        const auto b_indices = denseidx::index_functions(ket_basis);

        sparsity::check_pattern(pattern, a_indices, b_indices);

        _compute_pair_blocks(pattern, molecule, bra_basis, ket_basis, charges, points, a_indices, b_indices, distributor);

        _compute_diagonal_blocks(pattern, molecule, bra_basis, ket_basis, charges, points, a_indices, b_indices, distributor);
    }

    /// @brief Computes the nuclear potential matrix of a sparsity pattern.
    auto compute_matrix(const CSparsityPattern    &pattern,
                        const CMolecule           &molecule,
                        const CMolecularBasis     &bra_basis,
                        const CMolecularBasis     &ket_basis,
                        const std::vector<double> &charges,
                        const std::vector<double> &points) const -> CSparseMatrix;

    /// @brief Computes the nuclear potential matrix of a molecular basis and a set
    /// of point charges.
    /// @param molecule The molecule to compute the matrix of.
    /// @param basis The molecular basis on bra and ket sides.
    /// @param charges The magnitude of each point charge.
    /// @param points The position of each point charge, as a flat array of three
    /// values per charge, in bohr.
    /// @return The symmetric sparse nuclear potential matrix.
    auto compute_matrix(const CMolecule           &molecule,
                        const CMolecularBasis     &basis,
                        const std::vector<double> &charges,
                        const std::vector<double> &points) const -> CSparseMatrix;

    /// @brief Computes the nuclear potential matrix of a molecular basis and the
    /// nuclei of the molecule.
    /// @param molecule The molecule to compute the matrix of, whose charges and
    /// coordinates are the point charges.
    /// @param basis The molecular basis on bra and ket sides.
    /// @return The symmetric sparse nuclear potential matrix.
    auto compute_matrix(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix;

   private:
    /// @brief Gets the charges and their positions from the nuclei of a molecule.
    /// @param molecule The molecule.
    /// @return The charges and the flat array of three coordinates each, in bohr.
    static auto nuclei_of(const CMolecule &molecule) -> std::pair<std::vector<double>, std::vector<double>>;

    /// @brief Computes the integrals of the off-diagonal blocks.
    template <class D>
    auto _compute_pair_blocks(const CSparsityPattern              &pattern,
                              const CMolecule                     &molecule,
                              const CMolecularBasis               &bra_basis,
                              const CMolecularBasis               &ket_basis,
                              const std::vector<double>           &charges,
                              const std::vector<double>           &points,
                              const denseidx::TBasisFunctionIndex &a_indices,
                              const denseidx::TBasisFunctionIndex &b_indices,
                              D                                   &distributor) const -> void;

    /// @brief Computes the integrals of the diagonal blocks.
    /// @note The operator does not vanish for two basis functions on the same atom
    /// and has no closed form there, so these blocks are computed with the same
    /// kernels as the others and the position of the atom does enter.
    template <class D>
    auto _compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                  const CMolecule                     &molecule,
                                  const CMolecularBasis               &bra_basis,
                                  const CMolecularBasis               &ket_basis,
                                  const std::vector<double>           &charges,
                                  const std::vector<double>           &points,
                                  const denseidx::TBasisFunctionIndex &a_indices,
                                  const denseidx::TBasisFunctionIndex &b_indices,
                                  D                                   &distributor) const -> void;

    /// @brief The screening threshold of the integrals.
    double _threshold;

    /// @brief The target number of atom pairs of a block, zero to choose it from
    /// the number of the threads and the number of the atom pairs.
    size_t _block_size;

    /// @brief The smallest target number of atom pairs of a block this driver lets
    /// the chosen block size fall to, against the 256 of sparsity::min_block_size.
    /// @note Measured on 2026-09-15 by sweeping the block size at one, four and
    /// fourteen threads. Every case has the same shape: the time rises steeply below
    /// 512 atom pairs and is flat from there to 32768, so the cost is a fixed one per
    /// block and not the buffer, which runs from 0.1 to 56 MB over that range with no
    /// structure in the timings. A nuclear potential kernel evaluates every pair
    /// against every charge, which is what makes its work per block large enough for
    /// the floor of the overlap to bind too early: 256 is what the chosen size falls
    /// to for anything under about 150 atoms, and it costs those 4 to 6 per cent.
    /// The molecules above that choose a block size far above either floor and are
    /// untouched by this.
    static constexpr size_t _min_block_pairs = 512;
};

template <class D>
auto
CSimdNuclearPotentialDriver::_compute_pair_blocks(const CSparsityPattern              &pattern,
                                                   const CMolecule                     &molecule,
                                                   const CMolecularBasis               &bra_basis,
                                                   const CMolecularBasis               &ket_basis,
                                                   const std::vector<double>           &charges,
                                                   const std::vector<double>           &points,
                                                   const denseidx::TBasisFunctionIndex &a_indices,
                                                   const denseidx::TBasisFunctionIndex &b_indices,
                                                   D                                   &distributor) const -> void
{
    // NOTE: the unit of work is one combination of basis functions of one block and
    // not one block, for the reason the overlap driver gives: a block carries a fixed
    // cost of some microseconds and an ordinary molecule holds tens of them whatever
    // the number of the threads, while the combinations are tens to hundreds and are
    // independent of one another.

    const auto nblocks = static_cast<size_t>(pattern.number_of_pair_blocks());

    if (nblocks == 0) return;

    std::vector<size_t> order(nblocks);

    std::vector<double> costs(nblocks, 0.0);

    auto arena_rows = size_t{0};

    auto arena_cols = size_t{0};

    for (size_t iblk = 0; iblk < nblocks; iblk++)
    {
        const auto &block = pattern.pair_block(iblk);

        order[iblk] = iblk;

        costs[iblk] = block.weight();

        arena_rows = std::max(arena_rows,
                              simdnpot::number_of_buffer_rows(bra_basis.basis_set(block.bra_index()).max_angular_momentum(),
                                                              ket_basis.basis_set(block.ket_index()).max_angular_momentum()));

        arena_cols = std::max(arena_cols, block.number_of_pairs());
    }

    std::ranges::sort(order, [&](const size_t a, const size_t b) { return costs[a] > costs[b]; });

    std::vector<TPairTask> tasks;

    for (const auto iblk : order)
    {
        const auto &block = pattern.pair_block(iblk);

        const auto &a_index = a_indices[block.bra_index()];

        const auto &b_index = b_indices[block.ket_index()];

        for (size_t i = 0; i < a_index.size(); i++)
        {
            for (size_t j = 0; j < b_index.size(); j++)
            {
                const auto [la, ia] = a_index[i];

                const auto [lb, jb] = b_index[j];

                if (const auto nvalues = block.number_of_pairs(la, ia, lb, jb); nvalues > 0)
                {
                    tasks.push_back({iblk, i, j, nvalues});
                }
            }
        }
    }

    const auto ntasks = static_cast<int>(tasks.size());

    if (ntasks == 0) return;

    std::vector<CSimdMatrix> coordinates(nblocks);

    const auto nblk = static_cast<int>(nblocks);

#pragma omp parallel if (ntasks > 1)
    {
#pragma omp for schedule(dynamic)
        for (int iblk = 0; iblk < nblk; iblk++)
        {
            coordinates[static_cast<size_t>(iblk)] =
                simdfunc::make_coordinates(pattern.pair_block(static_cast<size_t>(iblk)), molecule);
        }

        auto arena = CSimdMatrix(arena_rows, arena_cols);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

#pragma omp for schedule(dynamic)
        for (int itask = 0; itask < ntasks; itask++)
        {
            const auto &task = tasks[static_cast<size_t>(itask)];

            const auto &block = pattern.pair_block(task.iblock);

            const auto &a_basis = bra_basis.basis_set(block.bra_index());

            const auto &b_basis = ket_basis.basis_set(block.ket_index());

            const auto [la, ia] = a_indices[block.bra_index()][task.i];

            const auto [lb, jb] = b_indices[block.ket_index()][task.j];

            const auto ncomps = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 2>{la, lb}));

            auto *values = distributor.target(block, task.iblock, la, ia, lb, jb, task.nvalues, ncomps);

            simdnpot::compute_nuclear_potential(values,
                                                task.nvalues,
                                                a_basis.functions()[task.i],
                                                b_basis.functions()[task.j],
                                                coordinates[task.iblock],
                                                charges,
                                                points,
                                                buffer,
                                                pattern.get_threshold());

            distributor.commit(block, task.iblock, la, ia, lb, jb, task.nvalues, ncomps);
        }
    }
}

template <class D>
auto
CSimdNuclearPotentialDriver::_compute_diagonal_blocks(const CSparsityPattern              &pattern,
                                                       const CMolecule                     &molecule,
                                                       const CMolecularBasis               &bra_basis,
                                                       const CMolecularBasis               &ket_basis,
                                                       const std::vector<double>           &charges,
                                                       const std::vector<double>           &points,
                                                       const denseidx::TBasisFunctionIndex &a_indices,
                                                       const denseidx::TBasisFunctionIndex &b_indices,
                                                       D                                   &distributor) const -> void
{
    // NOTE: the atoms of a diagonal block stand for atom pairs whose two atoms are
    // the same one, so the coordinates are built here rather than by make_coordinates,
    // which takes the pairs of an off-diagonal block. Bra and ket hold the same
    // position, the vector between them is zero and so is its square.

    const auto &coords = molecule.coordinates("au");

    for (size_t iblk = 0; iblk < pattern.number_of_diagonal_blocks(); iblk++)
    {
        const auto &block = pattern.diagonal_block(iblk);

        errors::assertMsgCritical(block.get_storage() == diagstor::full,
                                  std::string("SimdNuclearPotentialDriver: The diagonal blocks of the nuclear potential "
                                              "hold every angular component and cannot be stored as one value"));

        const auto &atoms = block.atoms();

        const auto natoms = atoms.size();

        if (natoms == 0) continue;

        auto positions = CSimdMatrix(10, natoms);

        for (size_t k = 0; k < natoms; k++)
        {
            const auto r = coords[static_cast<size_t>(atoms[k])].coordinates();

            for (size_t axis = 0; axis < 3; axis++)
            {
                positions.data(axis)[k]     = r[axis];
                positions.data(axis + 3)[k] = r[axis];
                positions.data(axis + 6)[k] = 0.0;
            }

            positions.data(9)[k] = 0.0;
        }

        const auto &a_basis = bra_basis.basis_set(block.bra_index());

        const auto &b_basis = ket_basis.basis_set(block.ket_index());

        const auto &a_index = a_indices[block.bra_index()];

        const auto &b_index = b_indices[block.ket_index()];

        const auto rows = simdnpot::number_of_buffer_rows(a_basis.max_angular_momentum(), b_basis.max_angular_momentum());

        auto arena = CSimdMatrix(rows, natoms);

        auto buffer = CSimdMatrix(arena.data(), arena.capacity());

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

                auto *values = distributor.diagonal_target(block, iblk, la, ia, lb, jb);

                simdnpot::compute_nuclear_potential(values,
                                                    natoms,
                                                    a_basis.functions()[i],
                                                    b_basis.functions()[j],
                                                    positions,
                                                    charges,
                                                    points,
                                                    buffer,
                                                    pattern.get_threshold());

                distributor.diagonal_commit(block, iblk, la, ia, lb, jb);
            }
        }
    }
}

#endif /* SimdNuclearPotentialDriver_hpp */

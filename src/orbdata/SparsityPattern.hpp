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


#ifndef SparsityPattern_hpp
#define SparsityPattern_hpp

#include <cstddef>
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairGroup.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "ScreeningFunc.hpp"

/// @brief Class CSparsityPattern holds the sparsity patterns of the atom pair
/// blocks a two-center quantity is computed in, without the values themselves.
/// @note The pattern is what an integrals driver needs to know before it computes
/// anything: which atom pairs each block holds, in which order, and how many of
/// them survive the screening of each combination of basis functions. It is
/// separated from the containers of the values, so that one pattern serves a
/// sparse matrix, a contraction with a density, or any other consumer, and so
/// that the choice of the number of atom pairs of a block is made by the caller.
class CSparsityPattern
{
   public:
    /// @brief The default constructor.
    CSparsityPattern()

        : _pair_blocks{}

        , _diagonal_blocks{}

        , _type(mat_t::general)

        , _threshold(0.0)
    {
    }

    /// @brief The constructor with sparsity patterns of the blocks and matrix type.
    /// @param pair_blocks The sparsity patterns of the off-diagonal blocks.
    /// @param diagonal_blocks The sparsity patterns of the diagonal blocks.
    /// @param mat_type The type of quantity the pattern describes.
    /// @param threshold The screening threshold the blocks were described under.
    CSparsityPattern(std::vector<CAtomBasisPairSparsity>     pair_blocks,
                     std::vector<CAtomBasisDiagonalSparsity> diagonal_blocks,
                     const mat_t                             mat_type,
                     const double                            threshold)

        : _pair_blocks(std::move(pair_blocks))

        , _diagonal_blocks(std::move(diagonal_blocks))

        , _type(mat_type)

        , _threshold(threshold)
    {
    }

    /// @brief Gets number of off-diagonal blocks.
    /// @return The number of off-diagonal blocks.
    auto
    number_of_pair_blocks() const -> size_t
    {
        return _pair_blocks.size();
    }

    /// @brief Gets number of diagonal blocks.
    /// @return The number of diagonal blocks.
    auto
    number_of_diagonal_blocks() const -> size_t
    {
        return _diagonal_blocks.size();
    }

    /// @brief Gets sparsity pattern of the off-diagonal block with specific index.
    /// @param index The index of off-diagonal block.
    /// @return The constant reference to the sparsity pattern.
    auto
    pair_block(const size_t index) const -> const CAtomBasisPairSparsity &
    {
        errors::assertMsgCritical(index < _pair_blocks.size(),
                                  std::string("SparsityPattern.pair_block: Index of block is out of range"));

        return _pair_blocks[index];
    }

    /// @brief Gets sparsity pattern of the diagonal block with specific index.
    /// @param index The index of diagonal block.
    /// @return The constant reference to the sparsity pattern.
    auto
    diagonal_block(const size_t index) const -> const CAtomBasisDiagonalSparsity &
    {
        errors::assertMsgCritical(index < _diagonal_blocks.size(),
                                  std::string("SparsityPattern.diagonal_block: Index of block is out of range"));

        return _diagonal_blocks[index];
    }

    /// @brief Gets sparsity patterns of all off-diagonal blocks.
    /// @return The constant reference to the vector of sparsity patterns.
    auto
    pair_blocks() const -> const std::vector<CAtomBasisPairSparsity> &
    {
        return _pair_blocks;
    }

    /// @brief Gets sparsity patterns of all diagonal blocks.
    /// @return The constant reference to the vector of sparsity patterns.
    auto
    diagonal_blocks() const -> const std::vector<CAtomBasisDiagonalSparsity> &
    {
        return _diagonal_blocks;
    }

    /// @brief Gets type of the quantity the pattern describes.
    /// @return The type of quantity.
    auto
    get_type() const -> mat_t
    {
        return _type;
    }

    /// @brief Gets screening threshold the blocks were described under.
    /// @return The screening threshold.
    /// @note This is the threshold which dropped the atom pairs of the blocks, so
    /// it is the threshold a driver screens the pairs of primitives with. Taking it
    /// from the pattern rather than from the driver keeps the two screenings of a
    /// computation from disagreeing when a pattern is shared or is formed by another
    /// driver.
    auto
    get_threshold() const -> double
    {
        return _threshold;
    }

   private:
    /// @brief The sparsity patterns of the off-diagonal blocks.
    std::vector<CAtomBasisPairSparsity> _pair_blocks;

    /// @brief The sparsity patterns of the diagonal blocks.
    std::vector<CAtomBasisDiagonalSparsity> _diagonal_blocks;

    /// @brief The type of the quantity the pattern describes.
    mat_t _type;

    /// @brief The screening threshold the blocks were described under.
    double _threshold;
};

namespace sparsity {  // sparsity namespace

/// @brief The number of blocks per thread aimed at when the target number of atom
/// pairs of a block is chosen. The blocks are a few per thread, so that dynamic
/// scheduling has enough of them to even out the ones which differ in cost, and no
/// more, as a block carries a fixed cost and the blocks contend for the memory.
/// Measured on fourteen threads, where two per thread is five percent better than
/// four and four is twice as good as sixteen.
inline constexpr size_t blocks_per_thread = 2;

/// @brief The smallest target number of atom pairs of a block chosen. A block
/// carries a fixed cost which does not shrink with the atom pairs it holds, chiefly
/// the bisection of the screening over the pairs of primitives, so a molecule too
/// small to fill the threads is divided into fewer blocks rather than into blocks
/// whose fixed cost outweighs their work.
inline constexpr size_t min_block_size = 2048;

/// @brief Checks that a sparsity pattern describes a pair of molecular bases.
/// @param pattern The sparsity pattern to check.
/// @param a_indices The index of the basis functions of the basis on bra side.
/// @param b_indices The index of the basis functions of the basis on ket side.
/// @note A pattern carries neither the molecule nor the bases it was formed from,
/// so a driver handed a pattern of other bases would address an atom basis which is
/// not there, or write the values of a number of combinations of basis functions
/// which is not the number the pattern reserved. Both are caught here, before
/// anything is computed. A pattern of another molecule is caught by the coordinates,
/// which assert that every atom of a block is an atom of the molecule.
inline auto
check_pattern(const CSparsityPattern              &pattern,
              const denseidx::TBasisFunctionIndex &a_indices,
              const denseidx::TBasisFunctionIndex &b_indices) -> void
{
    const auto check_block = [&](const auto &block) {
        const auto ibra = static_cast<size_t>(block.bra_index());

        const auto iket = static_cast<size_t>(block.ket_index());

        errors::assertMsgCritical((block.bra_index() >= 0) && (ibra < a_indices.size()) && (block.ket_index() >= 0) &&
                                      (iket < b_indices.size()),
                                  std::string("SparsityPattern.check_pattern: Block addresses an atom basis which is not in the basis"));

        errors::assertMsgCritical((block.number_of_bra_basis_functions() == a_indices[ibra].size()) &&
                                      (block.number_of_ket_basis_functions() == b_indices[iket].size()),
                                  std::string("SparsityPattern.check_pattern: Block was described for other molecular bases"));
    };

    for (const auto &block : pattern.pair_blocks())
    {
        check_block(block);
    }

    for (const auto &block : pattern.diagonal_blocks())
    {
        check_block(block);
    }
}

/// @brief Divides the atom basis pair groups into the blocks the integrals are
/// computed in, ordered by interatomic distance.
/// @param molecule The molecule to compute interatomic distances from.
/// @param groups The atom basis pair groups to divide.
/// @param block_size The target number of atom pairs of a block, or zero to choose
/// it from the number of the threads and the number of the atom pairs.
/// @return The vector of blocks.
/// @note This is the half of the pattern which depends on the geometry and the
/// bases alone, and not on the operator or the threshold. It is kept apart so that
/// the screening below can be repeated on it, and so that it can be shared between
/// operators without changing any caller when that becomes worth doing.
inline auto
make_blocks(const CMolecule &molecule, std::vector<CAtomBasisPairGroup> &groups, const size_t block_size)
    -> std::vector<CAtomBasisPairGroup>
{
    // NOTE: the atom basis pair groups are as many as the pairs of the unique atom
    // bases, so their number is set by the variety of the elements of the molecule
    // and not by its size, and the largest of them holds a third of the atom pairs.
    // Dividing them into blocks of a target number of atom pairs makes the number of
    // the blocks follow the size of the molecule instead, so the work of every stage
    // below divides for any number of threads.

    const auto nblock_pairs =
        (block_size == 0) ? CAtomBasisPairGroup::make_block_size(groups, blocks_per_thread, min_block_size) : block_size;

    auto blocks = (nblock_pairs == 0) ? std::move(groups) : CAtomBasisPairGroup::divide(groups, nblock_pairs);

    // NOTE: the atom pairs of all the blocks are ordered by interatomic distance
    // before the sparsity patterns are described, as the patterns are read off the
    // leading atom pairs which survive the screening. A block holds a subrange of the
    // atom pairs of its group and is ordered within itself, which is all the bisection
    // of the screening needs, as the screening keeps an atom pair or drops it on its
    // own distance.

    CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

    return blocks;
}

/// @brief Describes the sparsity of the blocks under an integral bound.
/// @param blocks The blocks, ordered by interatomic distance.
/// @param bound The integral bound, evaluated as bound(bra_function, ket_function, distance).
/// @param threshold The screening threshold.
/// @param storage The storage layout of the diagonal blocks.
/// @param mat_type The type of quantity the pattern describes.
/// @return The sparsity pattern.
template <typename B>
inline auto
describe(const std::vector<CAtomBasisPairGroup> &blocks,
         const B                                &bound,
         const double                            threshold,
         const diagstor                          storage,
         const mat_t                             mat_type) -> CSparsityPattern
{
    // NOTE: the patterns are held in a vector indexed by the block, so that they are
    // described in any order and added in the order of the blocks. The layout of the
    // values of a consumer therefore does not depend on the scheduling.

    const auto nblocks = static_cast<int>(blocks.size());

    std::vector<std::optional<CAtomBasisPairSparsity>> pair_blocks(blocks.size());

    std::vector<std::optional<CAtomBasisDiagonalSparsity>> diagonal_blocks(blocks.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nblocks; i++)
    {
        pair_blocks[i].emplace(blocks[i], bound, threshold);

        diagonal_blocks[i].emplace(blocks[i], storage);
    }

    // NOTE: an empty pattern is left out, as it carries no values of its own.

    std::vector<CAtomBasisPairSparsity> pairs;

    std::vector<CAtomBasisDiagonalSparsity> diagonals;

    for (int i = 0; i < nblocks; i++)
    {
        if (pair_blocks[i]->number_of_pairs() > 0) pairs.push_back(std::move(*pair_blocks[i]));

        if (diagonal_blocks[i]->number_of_atoms() > 0) diagonals.push_back(std::move(*diagonal_blocks[i]));
    }

    return CSparsityPattern(std::move(pairs), std::move(diagonals), mat_type, threshold);
}

/// @brief Creates the sparsity pattern of a two-center quantity.
/// @param molecule The molecule to compute interatomic distances from.
/// @param bra_basis The molecular basis on bra side.
/// @param ket_basis The molecular basis on ket side.
/// @param bound The integral bound, evaluated as bound(bra_function, ket_function, distance).
/// @param threshold The screening threshold.
/// @param mat_type The type of quantity the pattern describes.
/// @param storage The storage layout of the diagonal blocks.
/// @param block_size The target number of atom pairs of a block, or zero to choose
/// it from the number of the threads and the number of the atom pairs.
/// @return The sparsity pattern.
template <typename B>
inline auto
make_pattern(const CMolecule       &molecule,
             const CMolecularBasis &bra_basis,
             const CMolecularBasis &ket_basis,
             const B               &bound,
             const double           threshold,
             const mat_t            mat_type,
             const diagstor         storage,
             const size_t           block_size = 0) -> CSparsityPattern
{
    // NOTE: a symmetric or antisymmetric quantity needs the upper triangle of the
    // atom basis pair groups only, while a general one needs their full direct
    // product, which the two molecular bases factory delivers even when handed the
    // same molecular basis twice.

    auto groups = ((mat_type != mat_t::general) && (&bra_basis == &ket_basis)) ? bra_basis.basis_pair_groups()
                                                                              : bra_basis.basis_pair_groups(ket_basis);

    auto blocks = make_blocks(molecule, groups, block_size);

    return describe(blocks, bound, threshold, storage, mat_type);
}

/// @brief Creates the sparsity pattern of a two-center quantity from a named bound.
/// @param bound The named integral bound to screen atom pairs with.
/// @return The sparsity pattern.
inline auto
make_pattern(const CMolecule       &molecule,
             const CMolecularBasis &bra_basis,
             const CMolecularBasis &ket_basis,
             const screener         bound,
             const double           threshold,
             const mat_t            mat_type,
             const diagstor         storage,
             const size_t           block_size = 0) -> CSparsityPattern
{
    if (bound == screener::overlap)
    {
        return make_pattern(molecule, bra_basis, ket_basis, screenfunc::two_center_overlap_bound, threshold, mat_type, storage, block_size);
    }

    if (bound == screener::kinetic_energy)
    {
        return make_pattern(
            molecule, bra_basis, ket_basis, screenfunc::two_center_kinetic_energy_bound, threshold, mat_type, storage, block_size);
    }

    if (bound == screener::nuclear_potential)
    {
        return make_pattern(
            molecule, bra_basis, ket_basis, screenfunc::two_center_nuclear_potential_bound, threshold, mat_type, storage, block_size);
    }

    errors::assertMsgCritical(false, std::string("SparsityPattern.make_pattern: Integral bound is not a two-center bound"));

    return CSparsityPattern();
}

}  // namespace sparsity

#endif /* SparsityPattern_hpp */

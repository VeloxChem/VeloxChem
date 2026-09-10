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



#ifndef TripleSparsityPattern_hpp
#define TripleSparsityPattern_hpp

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "AtomBasisGroup.hpp"
#include "AtomBasisPairGroup.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningFunc.hpp"

/// @brief Class CTripleSparsityPattern holds the sparsity patterns of the blocks a
/// three-center quantity is computed in, without the values themselves.
///
/// @note The pattern is what an integrals driver needs to know before it computes
/// anything: which atom pairs each block holds on a and b sides, which atoms it
/// holds on c side, and how many of the atom pairs survive the screening of each
/// combination of basis functions. It is separated from the containers of the
/// values, so that one pattern serves a sparse tensor, a contraction with fitting
/// coefficients, or any other consumer.
class CTripleSparsityPattern
{
   public:
    /// @brief The default constructor.
    CTripleSparsityPattern()

        : _blocks{}

        , _type(mat_t::general)

        , _threshold(0.0)
    {
    }

    /// @brief The constructor with sparsity patterns of the blocks, matrix type and
    /// screening threshold.
    /// @param blocks The sparsity patterns of the blocks.
    /// @param mat_type The type of quantity the pattern describes.
    /// @param threshold The screening threshold the blocks were described under.
    CTripleSparsityPattern(std::vector<CAtomBasisTripleSparsity> blocks, const mat_t mat_type, const double threshold)

        : _blocks(std::move(blocks))

        , _type(mat_type)

        , _threshold(threshold)
    {
    }

    /// @brief Gets sparsity patterns of all blocks.
    /// @return The constant reference to the vector of sparsity patterns.
    auto
    blocks() const -> const std::vector<CAtomBasisTripleSparsity> &
    {
        return _blocks;
    }

    /// @brief Gets sparsity pattern of the block with specific index.
    /// @param index The index of block.
    /// @return The constant reference to the sparsity pattern.
    auto
    block(const size_t index) const -> const CAtomBasisTripleSparsity &
    {
        errors::assertMsgCritical(index < _blocks.size(),
                                  std::string("TripleSparsityPattern.block: Index of block is out of range"));

        return _blocks[index];
    }

    /// @brief Gets number of blocks.
    /// @return The number of blocks.
    auto
    number_of_blocks() const -> size_t
    {
        return _blocks.size();
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
    /// @note This is the threshold which dropped the atom pairs of the blocks, so it
    /// is the threshold a driver screens with. Taking it from the pattern rather than
    /// from the caller keeps the two screenings of a computation from disagreeing
    /// when a pattern is shared.
    auto
    get_threshold() const -> double
    {
        return _threshold;
    }

   private:
    /// @brief The sparsity patterns of the blocks.
    std::vector<CAtomBasisTripleSparsity> _blocks;

    /// @brief The type of the quantity the pattern describes.
    mat_t _type;

    /// @brief The screening threshold the blocks were described under.
    double _threshold;
};

namespace sparsity {  // sparsity namespace

/// @brief The target number of atom pairs of a block of a three-center quantity when
/// the caller does not choose one.
/// @note A caller which knows the operator should choose the size itself, from the
/// rows its largest combination needs and a budget for the buffer. This is what the
/// callers which cannot do that fall back on. See the driver of the three-center
/// Coulomb quantity for the rule.
inline constexpr size_t triple_default_block_size = 32;

/// @brief Selects the atom basis groups on c side which carry the given atoms.
/// @param molecule The molecule the atoms belong to.
/// @param aux_basis The auxiliary molecular basis on c side.
/// @param aux_atoms The atoms on c side to select, as their indices in the molecule
/// and without repetition.
/// @return The vector of atom basis groups holding those atoms alone.
/// @note This is how a tensor too large for the memory is formed in batches: the
/// caller divides the atoms and the parts are described one after another.
inline auto
select_aux_groups(const CMolecule &molecule, const CMolecularBasis &aux_basis, const std::vector<int> &aux_atoms)
    -> std::vector<CAtomBasisGroup>
{
    errors::assertMsgCritical(!aux_atoms.empty(),
                              std::string("TripleSparsityPattern.select_aux_groups: The atoms on c side must not be empty"));

    const auto natoms = molecule.number_of_atoms();

    std::ranges::for_each(aux_atoms, [&](const auto atom) {
        errors::assertMsgCritical(
            (atom >= 0) && (atom < natoms),
            std::string("TripleSparsityPattern.select_aux_groups: Index of atom on c side is out of range"));
    });

    std::vector<CAtomBasisGroup> groups;

    std::ranges::for_each(aux_basis.basis_groups(), [&](const auto &group) {
        const std::unordered_set<int> atoms(group.atoms().begin(), group.atoms().end());

        std::vector<int> selected;

        selected.reserve(aux_atoms.size());

        std::ranges::copy_if(aux_atoms, std::back_inserter(selected), [&](const auto atom) { return atoms.contains(atom); });

        if (!selected.empty()) groups.push_back(CAtomBasisGroup(group.basis(), selected, group.index()));
    });

    return groups;
}

/// @brief Divides the atom basis pair groups into the blocks a three-center quantity
/// is computed in, ordered by interatomic distance.
/// @param molecule The molecule to compute interatomic distances from.
/// @param groups The atom basis pair groups on a and b sides to divide.
/// @return The vector of blocks.
/// @note This is the half of the pattern which depends on the geometry and the bases
/// alone, and not on the operator or the threshold.
inline auto
make_triple_blocks(const CMolecule &molecule, std::vector<CAtomBasisPairGroup> &groups, const size_t block_size)
    -> std::vector<CAtomBasisPairGroup>
{
    // NOTE: the atom basis pair groups are as many as the pairs of the unique atom
    // bases, so their number is set by the variety of the elements of the molecule
    // and not by its size. Dividing them into blocks of a target number of atom pairs
    // makes the number of the blocks follow the size of the molecule instead, so the
    // work below divides for any number of threads. Batching the c side does not do
    // this, as the cost of a block follows the atom pairs on the a and b sides.

    // NOTE: the size does not follow the number of the threads. The parallelism of a
    // three-center quantity comes from the combinations of basis functions of a block
    // and not from the blocks, and the size is bounded by the buffer a combination
    // needs rather than by the machine.

    const auto nblock_pairs = (block_size == 0) ? triple_default_block_size : block_size;

    auto blocks = CAtomBasisPairGroup::divide(groups, nblock_pairs);

    // NOTE: the atom pairs of all the blocks are ordered by interatomic distance
    // before the sparsity patterns are described, as the patterns are read off the
    // leading atom pairs which survive the screening. A block holds a subrange of the
    // atom pairs of its group and is ordered within itself, which is all the
    // bisection of the screening needs.

    CAtomBasisPairGroup::sort_by_distance(blocks, molecule);

    return blocks;
}

/// @brief Describes the sparsity of the blocks of a three-center quantity under an
/// integral bound.
/// @param blocks The blocks on a and b sides, ordered by interatomic distance.
/// @param aux_groups The atom basis groups on c side.
/// @param bound The integral bound.
/// @param threshold The screening threshold.
/// @param mat_type The type of quantity the pattern describes.
/// @return The sparsity pattern.
template <typename B>
inline auto
describe_triples(const std::vector<CAtomBasisPairGroup> &blocks,
                 const std::vector<CAtomBasisGroup>     &aux_groups,
                 const B                                &bound,
                 const double                            threshold,
                 const mat_t                             mat_type) -> CTripleSparsityPattern
{
    // NOTE: the patterns are held in a vector indexed by the block and the atom basis
    // group on c side, so that they are described in any order and added in the order
    // of that index. The layout of the values blocks of a consumer therefore does not
    // depend on the scheduling.

    const auto naux = aux_groups.size();

    const auto npatterns = static_cast<int>(blocks.size() * naux);

    std::vector<std::optional<CAtomBasisTripleSparsity>> patterns(blocks.size() * naux);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < npatterns; i++)
    {
        const auto index = static_cast<size_t>(i);

        patterns[index].emplace(blocks[index / naux], aux_groups[index % naux], bound, threshold);
    }

    // NOTE: an empty pattern is left out, as it carries no values block of its own.

    std::vector<CAtomBasisTripleSparsity> kept;

    for (int i = 0; i < npatterns; i++)
    {
        auto &pattern = patterns[static_cast<size_t>(i)];

        if (pattern->number_of_pairs() > 0) kept.push_back(std::move(*pattern));
    }

    return CTripleSparsityPattern(std::move(kept), mat_type, threshold);
}

/// @brief Creates the sparsity pattern of a three-center quantity, for the given
/// atoms on c side.
/// @param molecule The molecule to compute interatomic distances from.
/// @param basis The molecular basis on a and b sides.
/// @param aux_basis The auxiliary molecular basis on c side.
/// @param bound The integral bound.
/// @param threshold The screening threshold.
/// @param mat_type The type of quantity the pattern describes.
/// @param aux_atoms The atoms on c side to describe.
/// @return The sparsity pattern.
template <typename B>
inline auto
make_triple_pattern(const CMolecule        &molecule,
                    const CMolecularBasis  &basis,
                    const CMolecularBasis  &aux_basis,
                    const B                &bound,
                    const double            threshold,
                    const mat_t             mat_type,
                    const std::vector<int> &aux_atoms,
                    const size_t            block_size = 0) -> CTripleSparsityPattern
{
    auto groups = (mat_type == mat_t::general) ? basis.basis_pair_groups(basis) : basis.basis_pair_groups();

    auto blocks = make_triple_blocks(molecule, groups, block_size);

    return describe_triples(blocks, select_aux_groups(molecule, aux_basis, aux_atoms), bound, threshold, mat_type);
}

/// @brief Creates the sparsity pattern of a three-center quantity, for all atoms on
/// c side.
template <typename B>
inline auto
make_triple_pattern(const CMolecule       &molecule,
                    const CMolecularBasis &basis,
                    const CMolecularBasis &aux_basis,
                    const B               &bound,
                    const double           threshold,
                    const mat_t            mat_type,
                    const size_t           block_size = 0) -> CTripleSparsityPattern
{
    auto groups = (mat_type == mat_t::general) ? basis.basis_pair_groups(basis) : basis.basis_pair_groups();

    auto blocks = make_triple_blocks(molecule, groups, block_size);

    return describe_triples(blocks, aux_basis.basis_groups(), bound, threshold, mat_type);
}

/// @brief Checks that a sparsity pattern describes a molecular basis and an
/// auxiliary molecular basis.
/// @param pattern The sparsity pattern to check.
/// @param indices The index of the basis functions of the basis on a and b sides.
/// @param aux_indices The index of the basis functions of the auxiliary basis.
/// @note A pattern carries neither the molecule nor the bases it was formed from, so
/// a driver handed a pattern of other bases would address an atom basis which is not
/// there, or write a number of values which is not the number the pattern reserved.
inline auto
check_triple_pattern(const CTripleSparsityPattern        &pattern,
                     const denseidx::TBasisFunctionIndex &indices,
                     const denseidx::TBasisFunctionIndex &aux_indices) -> void
{
    for (const auto &block : pattern.blocks())
    {
        const auto ia = static_cast<size_t>(block.a_index());

        const auto ib = static_cast<size_t>(block.b_index());

        const auto ic = static_cast<size_t>(block.c_index());

        errors::assertMsgCritical(
            (block.a_index() >= 0) && (ia < indices.size()) && (block.b_index() >= 0) && (ib < indices.size()) &&
                (block.c_index() >= 0) && (ic < aux_indices.size()),
            std::string("TripleSparsityPattern.check_triple_pattern: Block addresses an atom basis which is not in the basis"));

        errors::assertMsgCritical(
            (block.number_of_a_basis_functions() == indices[ia].size()) &&
                (block.number_of_b_basis_functions() == indices[ib].size()) &&
                (block.number_of_c_basis_functions() == aux_indices[ic].size()),
            std::string("TripleSparsityPattern.check_triple_pattern: Block was described for other molecular bases"));
    }
}

}  // namespace sparsity

#endif /* TripleSparsityPattern_hpp */

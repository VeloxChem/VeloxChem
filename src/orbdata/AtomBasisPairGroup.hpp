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


#ifndef AtomBasisPairGroup_hpp
#define AtomBasisPairGroup_hpp

#include <algorithm>
#include <functional>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

#include "AtomBasis.hpp"
#include "AtomBasisGroup.hpp"

class CMolecule;

/// @brief Defines supported orderings of atom pairs:
/// pairord::none     - the atom pairs are in their construction order
/// pairord::distance - the atom pairs are ordered by ascending interatomic distance
enum class pairord
{
    none,
    distance
};

/// @brief Class CAtomBasisPairGroup associates a pair of atom bases with the
/// atom pairs formed by the atoms sharing them. Off-diagonal atom pairs are
/// stored as parallel bra and ket vectors, while atoms of diagonal atom pairs
/// are stored separately. Symmetry is established by construction, not derived
/// from the atom bases: the single argument constructor builds a symmetric
/// group, the two arguments constructor a non-symmetric one. The atom bases are
/// not owned by the group: they must outlive the group, and modifying the
/// molecular bases they originate from invalidates the group. Atom pairs are
/// created unordered and are reordered by sort_by_distance(), which also fills
/// in the interatomic distances.
class CAtomBasisPairGroup
{
   public:
    /// @brief The constructor with atom basis group shared by bra and ket sides.
    /// @param group The atom basis group on bra and ket sides.
    explicit CAtomBasisPairGroup(const CAtomBasisGroup &group)

        : _bra_basis(group.basis())

        , _ket_basis(group.basis())

        , _bra_atoms{}

        , _ket_atoms{}

        , _diagonal_atoms(group.atoms())

        , _distances{}

        , _bra_index(group.index())

        , _ket_index(group.index())

        , _symmetric(true)

        , _order(pairord::none)
    {
        const auto npairs = group.number_of_atoms() * (group.number_of_atoms() - 1) / 2;

        _bra_atoms.reserve(npairs);

        _ket_atoms.reserve(npairs);

        // NOTE: diagonal atom pairs are excluded here, as their atoms are all
        // atoms in group. Only the strict upper triangle is retained.

        std::ranges::for_each(group.atoms(), [&](const int i) {
            std::ranges::for_each(group.atoms(), [&](const int j) {
                if (i >= j) return;

                _bra_atoms.push_back(i);

                _ket_atoms.push_back(j);
            });
        });
    }

    /// @brief The constructor with bra and ket atom basis groups.
    /// @param bra The atom basis group on bra side.
    /// @param ket The atom basis group on ket side.
    CAtomBasisPairGroup(const CAtomBasisGroup &bra, const CAtomBasisGroup &ket)

        : _bra_basis(bra.basis())

        , _ket_basis(ket.basis())

        , _bra_atoms{}

        , _ket_atoms{}

        , _diagonal_atoms{}

        , _distances{}

        , _bra_index(bra.index())

        , _ket_index(ket.index())

        , _symmetric(false)

        , _order(pairord::none)
    {
        _bra_atoms.reserve(bra.number_of_atoms() * ket.number_of_atoms());

        _ket_atoms.reserve(bra.number_of_atoms() * ket.number_of_atoms());

        // NOTE: diagonal atom pairs are excluded here and collected below. The
        // full direct product of bra and ket atoms is retained, even if bra and
        // ket atom bases happen to be equal.

        std::ranges::for_each(bra.atoms(), [&](const int i) {
            std::ranges::for_each(ket.atoms(), [&](const int j) {
                if (i == j) return;

                _bra_atoms.push_back(i);

                _ket_atoms.push_back(j);
            });
        });

        // NOTE: atoms of diagonal atom pairs are the atoms carrying both atom
        // bases. Within a single molecular basis this is empty, but bra and ket
        // may come from different molecular bases, where the atom sets do
        // overlap. This assumes ascending atomic indices in both groups, as
        // provided by CMolecularBasis::basis_groups().

        std::ranges::set_intersection(bra.atoms(), ket.atoms(), std::back_inserter(_diagonal_atoms));
    }

    /// @brief Orders atom pairs of the groups by ascending interatomic distance
    /// and fills in the interatomic distances. Atoms of diagonal atom pairs are
    /// not affected, as their interatomic distance is zero.
    /// @param groups The atom basis pair groups to order.
    /// @param molecule The molecule providing the atomic coordinates.
    /// @note All the groups are ordered by one call, as they are ordered
    /// independently of each other and one group is thus ordered by one thread.
    static auto sort_by_distance(std::vector<CAtomBasisPairGroup> &groups, const CMolecule &molecule) -> void;

    /// @brief Partitions the atom pairs and the atoms of the diagonal atom pairs
    /// of the group into groups of similar size.
    /// @param nblocks The number of groups to partition the group into.
    /// @return The vector of partitioned groups.
    /// @note The atom pairs and the atoms of the diagonal atom pairs are divided
    /// independently, as their numbers are unrelated, and the remainder of a
    /// division is spread one over the leading groups, so the groups differ in
    /// size by at most one atom pair and one atom. A group left with no atom pair
    /// and no atom is dropped, so fewer groups than requested are returned when
    /// the group holds fewer atom pairs and atoms than that.
    /// @note The ordering of the atom pairs is carried over, as a leading
    /// subrange of atom pairs ordered by ascending interatomic distance is
    /// ordered by ascending interatomic distance as well.
    auto partition(const size_t nblocks) const -> std::vector<CAtomBasisPairGroup>;

    /// @brief Creates the atom basis pair groups of the atom basis groups of one
    /// molecular basis, i.e. the symmetric group of every atom basis group and
    /// the group of every pair of them.
    /// @param groups The atom basis groups of the molecular basis.
    /// @return The vector of atom basis pair groups, in the order of the groups.
    /// @note The atom pairs are formed by the threads. An atom of a molecular
    /// basis carries one atom basis, so the atom basis groups hold no atom in
    /// common: the atom pairs of a symmetric group are the strict upper triangle
    /// of its atoms and the atom pairs of a pair of groups are all their
    /// combinations, and the atom pair with a given index within a group follows
    /// from the index alone. The atom pairs are formed in the order the
    /// constructors form them, so that the groups are the same however they were
    /// made.
    static auto make_pair_groups(const std::vector<CAtomBasisGroup> &groups) -> std::vector<CAtomBasisPairGroup>;

    /// @brief The number of atom pairs one thread forms at a time.
    static constexpr size_t _pairs_per_chunk = 4096;

    /// @brief Gets atom basis on bra side.
    /// @return The constant reference to atom basis.
    auto
    bra_basis() const -> const CAtomBasis &
    {
        return _bra_basis.get();
    }

    /// @brief Gets atom basis on ket side.
    /// @return The constant reference to atom basis.
    auto
    ket_basis() const -> const CAtomBasis &
    {
        return _ket_basis.get();
    }

    /// @brief Gets index of atom basis on bra side among the unique atom bases
    /// of the molecular basis it originates from.
    /// @return The index of atom basis.
    auto
    bra_index() const -> int
    {
        return _bra_index;
    }

    /// @brief Gets index of atom basis on ket side among the unique atom bases
    /// of the molecular basis it originates from.
    /// @return The index of atom basis.
    auto
    ket_index() const -> int
    {
        return _ket_index;
    }

    /// @brief Gets maximum angular momentum of basis functions on bra side.
    /// @return The maximum angular momentum, -1 if atom basis is empty.
    auto
    bra_max_angular_momentum() const -> int
    {
        return _bra_basis.get().max_angular_momentum();
    }

    /// @brief Gets maximum angular momentum of basis functions on ket side.
    /// @return The maximum angular momentum, -1 if atom basis is empty.
    auto
    ket_max_angular_momentum() const -> int
    {
        return _ket_basis.get().max_angular_momentum();
    }

    /// @brief Gets number of basis functions with specific angular momentum on bra side.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of basis functions.
    auto
    bra_number_of_basis_functions(const int angular_momentum) const -> size_t
    {
        return _bra_basis.get().number_of_basis_functions(angular_momentum);
    }

    /// @brief Gets number of basis functions with specific angular momentum on ket side.
    /// @param angular_momentum The angular momentum of basis functions.
    /// @return The number of basis functions.
    auto
    ket_number_of_basis_functions(const int angular_momentum) const -> size_t
    {
        return _ket_basis.get().number_of_basis_functions(angular_momentum);
    }

    /// @brief Gets basis function with specific angular momentum and index on bra side.
    /// @param angular_momentum The angular momentum of basis function.
    /// @param index The index of basis function with given angular momentum.
    /// @return The constant reference to basis function.
    auto
    bra_basis_function(const int angular_momentum, const size_t index) const -> const CBasisFunction &
    {
        return _bra_basis.get().basis_function(angular_momentum, index);
    }

    /// @brief Gets basis function with specific angular momentum and index on ket side.
    /// @param angular_momentum The angular momentum of basis function.
    /// @param index The index of basis function with given angular momentum.
    /// @return The constant reference to basis function.
    auto
    ket_basis_function(const int angular_momentum, const size_t index) const -> const CBasisFunction &
    {
        return _ket_basis.get().basis_function(angular_momentum, index);
    }

    /// @brief Checks if bra and ket sides are interchangeable, i.e. atom pairs
    /// are stored as strict upper triangle of the atom pairs matrix.
    /// @return True if pair group is symmetric, False otherwise.
    auto
    is_symmetric() const -> bool
    {
        return _symmetric;
    }

    /// @brief Gets atomic indices of off-diagonal atom pairs on bra side.
    /// @return The constant reference to vector of atomic indices.
    auto
    bra_atoms() const -> const std::vector<int> &
    {
        return _bra_atoms;
    }

    /// @brief Gets atomic indices of off-diagonal atom pairs on ket side.
    /// @return The constant reference to vector of atomic indices.
    auto
    ket_atoms() const -> const std::vector<int> &
    {
        return _ket_atoms;
    }

    /// @brief Gets atomic indices of atoms in diagonal atom pairs.
    /// @return The constant reference to vector of atomic indices.
    auto
    diagonal_atoms() const -> const std::vector<int> &
    {
        return _diagonal_atoms;
    }

    /// @brief Gets interatomic distances of off-diagonal atom pairs.
    /// @return The constant reference to vector of interatomic distances, empty
    /// if atom pairs are not ordered by interatomic distance.
    auto
    distances() const -> const std::vector<double> &
    {
        return _distances;
    }

    /// @brief Gets pointer to interatomic distances of off-diagonal atom pairs.
    /// @return The constant pointer to interatomic distances. The pointed to
    /// range is empty, and must not be dereferenced, if atom pairs are not
    /// ordered by interatomic distance.
    auto
    distances_data() const -> const double *
    {
        return _distances.data();
    }

    /// @brief Gets ordering of atom pairs in group.
    /// @return The ordering of atom pairs.
    auto
    ordering() const -> pairord
    {
        return _order;
    }

    /// @brief Gets number of off-diagonal atom pairs in group.
    /// @return The number of off-diagonal atom pairs.
    auto
    number_of_pairs() const -> size_t
    {
        return _bra_atoms.size();
    }

    /// @brief Gets number of atoms in diagonal atom pairs.
    /// @return The number of atoms in diagonal atom pairs.
    auto
    number_of_diagonal_atoms() const -> size_t
    {
        return _diagonal_atoms.size();
    }

   private:
    /// @brief The constructor with the atom bases of a group and the number of its
    /// atom pairs, which are sized but not formed.
    /// @param bra The atom basis group on bra side.
    /// @param ket The atom basis group on ket side.
    /// @param npairs The number of off-diagonal atom pairs of the group.
    /// @param symmetric The flag indicating that bra and ket sides are interchangeable.
    /// @note The atom pairs are left with the content of the allocation, so this
    /// constructor is private and the caller which uses it fills them.
    CAtomBasisPairGroup(const CAtomBasisGroup &bra, const CAtomBasisGroup &ket, const size_t npairs, const bool symmetric)

        : _bra_basis(bra.basis())

        , _ket_basis(ket.basis())

        , _bra_atoms(npairs)

        , _ket_atoms(npairs)

        , _diagonal_atoms{}

        , _distances{}

        , _bra_index(bra.index())

        , _ket_index(ket.index())

        , _symmetric(symmetric)

        , _order(pairord::none)
    {
        if (symmetric)
        {
            _diagonal_atoms = bra.atoms();
        }
        else
        {
            std::ranges::set_intersection(bra.atoms(), ket.atoms(), std::back_inserter(_diagonal_atoms));
        }
    }

    /// @brief The constructor with the atom bases and the atom pairs of a group.
    /// @param bra_basis The atom basis on bra side.
    /// @param ket_basis The atom basis on ket side.
    /// @param bra_atoms The atomic indices of off-diagonal atom pairs on bra side.
    /// @param ket_atoms The atomic indices of off-diagonal atom pairs on ket side.
    /// @param diagonal_atoms The atomic indices of atoms in diagonal atom pairs.
    /// @param distances The interatomic distances of off-diagonal atom pairs.
    /// @param bra_index The index of atom basis on bra side.
    /// @param ket_index The index of atom basis on ket side.
    /// @param symmetric The flag indicating that bra and ket sides are interchangeable.
    /// @param order The ordering of the atom pairs.
    /// @note This constructor takes the state of a group as it is, rather than
    /// forming the atom pairs from the atom basis groups, so it is private: the
    /// ordering it is given is not checked against the atom pairs it is given.
    CAtomBasisPairGroup(const CAtomBasis   &bra_basis,
                        const CAtomBasis   &ket_basis,
                        std::vector<int>    bra_atoms,
                        std::vector<int>    ket_atoms,
                        std::vector<int>    diagonal_atoms,
                        std::vector<double> distances,
                        const int           bra_index,
                        const int           ket_index,
                        const bool          symmetric,
                        const pairord       order)

        : _bra_basis(bra_basis)

        , _ket_basis(ket_basis)

        , _bra_atoms(std::move(bra_atoms))

        , _ket_atoms(std::move(ket_atoms))

        , _diagonal_atoms(std::move(diagonal_atoms))

        , _distances(std::move(distances))

        , _bra_index(bra_index)

        , _ket_index(ket_index)

        , _symmetric(symmetric)

        , _order(order)
    {
    }

    /// @brief The atom basis on bra side.
    std::reference_wrapper<const CAtomBasis> _bra_basis;

    /// @brief The atom basis on ket side.
    std::reference_wrapper<const CAtomBasis> _ket_basis;

    /// @brief The atomic indices of off-diagonal atom pairs on bra side.
    std::vector<int> _bra_atoms;

    /// @brief The atomic indices of off-diagonal atom pairs on ket side.
    std::vector<int> _ket_atoms;

    /// @brief The atomic indices of atoms in diagonal atom pairs.
    std::vector<int> _diagonal_atoms;

    /// @brief The interatomic distances of off-diagonal atom pairs.
    std::vector<double> _distances;

    /// @brief The index of atom basis on bra side among the unique atom bases
    /// of the molecular basis.
    int _bra_index;

    /// @brief The index of atom basis on ket side among the unique atom bases
    /// of the molecular basis.
    int _ket_index;

    /// @brief The flag indicating that bra and ket sides are interchangeable.
    bool _symmetric;

    /// @brief The ordering of atom pairs in group.
    pairord _order;
};

#endif /* AtomBasisPairGroup_hpp */

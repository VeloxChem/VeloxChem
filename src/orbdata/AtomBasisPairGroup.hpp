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

    /// @brief Orders atom pairs by ascending interatomic distance and fills in
    /// the interatomic distances. Atoms of diagonal atom pairs are not affected,
    /// as their interatomic distance is zero.
    /// @param molecule The molecule providing the atomic coordinates.
    auto sort_by_distance(const CMolecule &molecule) -> void;

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

    /// @brief The flag indicating that bra and ket sides are interchangeable.
    bool _symmetric;

    /// @brief The ordering of atom pairs in group.
    pairord _order;
};

#endif /* AtomBasisPairGroup_hpp */

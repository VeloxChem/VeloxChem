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

/// @brief Class CAtomBasisPairGroup associates a pair of atom bases with the
/// atom pairs formed by the atoms sharing them. Off-diagonal atom pairs are
/// stored as parallel bra and ket vectors, while atoms of diagonal atom pairs
/// are stored separately. The atom bases are not owned by the group: they must
/// outlive the group, and modifying the molecular bases they originate from
/// invalidates the group.
class CAtomBasisPairGroup
{
   public:
    /// @brief The constructor with bra and ket atom basis groups.
    /// @param bra The atom basis group on bra side.
    /// @param ket The atom basis group on ket side.
    CAtomBasisPairGroup(const CAtomBasisGroup &bra, const CAtomBasisGroup &ket)

        : _bra_basis(bra.basis())

        , _ket_basis(ket.basis())

        , _bra_atoms{}

        , _ket_atoms{}

        , _diagonal_atoms{}

        , _symmetric((&bra.basis() == &ket.basis()) || (bra.basis() == ket.basis()))
    {
        _bra_atoms.reserve(bra.number_of_atoms() * ket.number_of_atoms());

        _ket_atoms.reserve(bra.number_of_atoms() * ket.number_of_atoms());

        // NOTE: diagonal atom pairs are excluded here and collected below. For
        // symmetric pair groups only the strict upper triangle is retained.

        std::ranges::for_each(bra.atoms(), [&](const int i) {
            std::ranges::for_each(ket.atoms(), [&](const int j) {
                if (i == j) return;

                if (_symmetric && (i > j)) return;

                _bra_atoms.push_back(i);

                _ket_atoms.push_back(j);
            });
        });

        // NOTE: atoms of diagonal atom pairs are the atoms carrying both atom
        // bases. Within a single molecular basis this is empty unless the pair
        // group is symmetric, but bra and ket may come from different molecular
        // bases, where the atom sets do overlap.

        std::ranges::set_intersection(bra.atoms(), ket.atoms(), std::back_inserter(_diagonal_atoms));
    }

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

    /// @brief Checks if bra and ket sides share the same atom basis.
    /// @return True if atom bases are the same, False otherwise.
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

    /// @brief The flag indicating that bra and ket atom bases are the same.
    bool _symmetric;
};

#endif /* AtomBasisPairGroup_hpp */

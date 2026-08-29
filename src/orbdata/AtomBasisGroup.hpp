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


#ifndef AtomBasisGroup_hpp
#define AtomBasisGroup_hpp

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

#include "AtomBasis.hpp"

/// @brief Class CAtomBasisGroup associates an atom basis with the atoms which
/// share it, and with its index among the unique atom bases of the molecular
/// basis it originates from. The atom basis is not owned by the group: it must
/// outlive the group, and modifying the molecular basis it originates from
/// invalidates the group.
class CAtomBasisGroup
{
   public:
    /// @brief The constructor with atom basis, its atoms and its index.
    /// @param basis The atom basis shared by the atoms.
    /// @param atoms The vector of atomic indices using the atom basis.
    /// @param index The index of atom basis among the unique atom bases of the
    /// molecular basis.
    CAtomBasisGroup(const CAtomBasis &basis, const std::vector<int> &atoms, const int index)

        : _basis(basis)

        , _atoms(atoms)

        , _index(index)
    {
    }

    /// @brief Slices atom basis group by selecting a range of its atoms.
    /// @param offset The index of the first selected atom in group.
    /// @param natoms The number of selected atoms.
    /// @return The atom basis group with the selected atoms.
    /// @note The selected range is clamped to the atoms of the group, so an
    /// offset beyond the last atom gives an empty atom basis group and the last
    /// slice of a fixed stride is short rather than out of range. The sliced
    /// group refers to the same atom basis and keeps its index, as it describes
    /// the same atom basis of the same molecular basis.
    auto
    slice(const size_t offset, const size_t natoms) const -> CAtomBasisGroup
    {
        const auto first = std::min(offset, _atoms.size());

        // NOTE: the count is clamped against the atoms left after the offset, so
        // that a large number of atoms cannot overflow into an inverted range.

        const auto count = std::min(natoms, _atoms.size() - first);

        return CAtomBasisGroup(_basis.get(), std::vector<int>(_atoms.begin() + first, _atoms.begin() + first + count), _index);
    }

    /// @brief Gets atom basis shared by the atoms in group.
    /// @return The constant reference to atom basis.
    auto
    basis() const -> const CAtomBasis &
    {
        return _basis.get();
    }

    /// @brief Gets atomic indices of atoms in group.
    /// @return The constant reference to vector of atomic indices.
    auto
    atoms() const -> const std::vector<int> &
    {
        return _atoms;
    }

    /// @brief Gets index of atom basis among the unique atom bases of the
    /// molecular basis the group originates from.
    /// @return The index of atom basis.
    auto
    index() const -> int
    {
        return _index;
    }

    /// @brief Gets number of atoms in group.
    /// @return The number of atoms in group.
    auto
    number_of_atoms() const -> size_t
    {
        return _atoms.size();
    }

   private:
    /// @brief The atom basis shared by the atoms in group.
    std::reference_wrapper<const CAtomBasis> _basis;

    /// @brief The vector of atomic indices using the atom basis.
    std::vector<int> _atoms;

    /// @brief The index of atom basis among the unique atom bases of the
    /// molecular basis.
    int _index;
};

#endif /* AtomBasisGroup_hpp */

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

#include <functional>
#include <vector>

#include "AtomBasis.hpp"

/// @brief Class CAtomBasisGroup associates an atom basis with the atoms which
/// share it. The atom basis is not owned by the group: it must outlive the
/// group, and modifying the molecular basis it originates from invalidates the
/// group.
class CAtomBasisGroup
{
   public:
    /// @brief The constructor with atom basis and its atoms.
    /// @param basis The atom basis shared by the atoms.
    /// @param atoms The vector of atomic indices using the atom basis.
    CAtomBasisGroup(const CAtomBasis &basis, const std::vector<int> &atoms)

        : _basis(basis)

        , _atoms(atoms)
    {
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
};

#endif /* AtomBasisGroup_hpp */

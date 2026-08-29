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


#ifndef AtomBasisDiagonalSparsity_hpp
#define AtomBasisDiagonalSparsity_hpp

#include <cstddef>
#include <vector>

#include "AtomBasisPairGroup.hpp"

/// @brief Class CAtomBasisDiagonalSparsity stores the atoms of the diagonal atom
/// pairs of an atom basis pair group. Diagonal atom pairs are never screened
/// out, as bra and ket basis functions are centered on the same atom, so no
/// counts of surviving atom pairs are needed and all combinations of basis
/// functions are retained. The off-diagonal atom pairs of the group are
/// described by CAtomBasisPairSparsity, whose storage requirements differ.
/// @note The atom bases are referred to by their index among the unique atom
/// bases of the molecular basis they originate from, so the sparsity pattern is
/// a self contained value. Within a single molecular basis a diagonal atom pair
/// carries the same atom basis on both sides, and the indices are equal. For a
/// pair of molecular bases the same atom carries a different atom basis on each
/// side, and the bra index refers to the bra molecular basis while the ket index
/// refers to the ket molecular basis. Resolving an index against the wrong
/// molecular basis addresses the wrong atom basis instead of failing.
class CAtomBasisDiagonalSparsity
{
   public:
    /// @brief The constructor with atom basis pair group.
    /// @param group The atom basis pair group to describe.
    explicit CAtomBasisDiagonalSparsity(const CAtomBasisPairGroup &group)

        : _bra_index(group.bra_index())

        , _ket_index(group.ket_index())

        , _atoms(group.diagonal_atoms())
    {
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

    /// @brief Gets atomic indices of atoms in diagonal atom pairs.
    /// @return The constant reference to vector of atomic indices.
    auto
    atoms() const -> const std::vector<int> &
    {
        return _atoms;
    }

    /// @brief Gets number of atoms in diagonal atom pairs.
    /// @return The number of atoms in diagonal atom pairs.
    auto
    number_of_atoms() const -> size_t
    {
        return _atoms.size();
    }

   private:
    /// @brief The index of atom basis on bra side among the unique atom bases
    /// of the molecular basis.
    int _bra_index;

    /// @brief The index of atom basis on ket side among the unique atom bases
    /// of the molecular basis.
    int _ket_index;

    /// @brief The atomic indices of atoms in diagonal atom pairs.
    std::vector<int> _atoms;
};

#endif /* AtomBasisDiagonalSparsity_hpp */

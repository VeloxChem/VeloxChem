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

/// @brief Defines supported storage layouts of the diagonal atom pair blocks:
/// diagstor::scalar - one value for each pair of basis functions with the same
///                    angular momentum, and none for the remaining pairs. The
///                    operator is spherically symmetric about the atom, so the
///                    block is diagonal in the angular components and its value
///                    does not depend on the position of the atom. Used by the
///                    overlap, kinetic energy and Coulomb metric integrals.
/// diagstor::full   - all angular components for each pair of basis functions,
///                    for each atom. The operator is not spherically symmetric
///                    about the atom, so different angular momenta couple and
///                    the values depend on the position of the atom. Used by the
///                    nuclear potential integrals.
enum class diagstor
{
    scalar,
    full
};

/// @brief Class CAtomBasisDiagonalSparsity stores the atoms of the diagonal atom
/// pairs of an atom basis pair group. Diagonal atom pairs are never screened
/// out, as bra and ket basis functions are centered on the same atom, so no
/// counts of surviving atom pairs are needed and all combinations of basis
/// functions are retained. The block is symmetric, so only the combinations of
/// basis functions with the index on bra side not exceeding the index on ket
/// side are stored. As basis functions are kept sorted by ascending angular
/// momentum, this retains the triangle of the combinations with equal angular
/// momenta and one of the two orderings of the combinations with different
/// angular momenta. This applies only when bra and ket sides carry the same atom
/// basis. A diagonal atom pair spanning two molecular bases carries a different
/// atom basis on each side, so the block is not symmetric and all combinations
/// of basis functions are stored. The off-diagonal atom pairs of the group are
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
    /// @brief The constructor with atom basis pair group and storage layout.
    /// @param group The atom basis pair group to describe.
    /// @param storage The storage layout of the diagonal atom pair blocks.
    CAtomBasisDiagonalSparsity(const CAtomBasisPairGroup &group, const diagstor storage);

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

    /// @brief Gets storage layout of the diagonal atom pair blocks.
    /// @return The storage layout.
    auto
    get_storage() const -> diagstor
    {
        return _storage;
    }

    /// @brief Gets number of values required to store the integrals of specific
    /// combination of basis functions.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The number of values.
    auto number_of_elements(const int    bra_angular_momentum,
                            const size_t bra_index,
                            const int    ket_angular_momentum,
                            const size_t ket_index) const -> size_t;

    /// @brief Gets offset of the values of specific combination of basis
    /// functions in the values block.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The offset of the values.
    auto element_offset(const int    bra_angular_momentum,
                        const size_t bra_index,
                        const int    ket_angular_momentum,
                        const size_t ket_index) const -> size_t;

    /// @brief Gets number of values required to store the integrals of all
    /// combinations of basis functions, i.e. the size of the values block.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        return _value_offsets.back();
    }

    /// @brief Checks if only the triangle of the combinations of basis functions
    /// is stored, which is the case when bra and ket sides carry the same atom
    /// basis.
    /// @return True if the combinations are stored as a triangle, False if all
    /// combinations are stored.
    auto
    is_triangular() const -> bool
    {
        return _triangular;
    }

    /// @brief Gets number of basis functions on bra side.
    /// @return The number of basis functions.
    auto
    number_of_bra_basis_functions() const -> size_t
    {
        return _bra_offsets.back();
    }

    /// @brief Gets number of basis functions on ket side.
    /// @return The number of basis functions.
    auto
    number_of_ket_basis_functions() const -> size_t
    {
        return _ket_offsets.back();
    }

   private:
    /// @brief Computes the flat index of a combination of basis functions, which
    /// is swapped into the stored triangle if given in the reverse order and only
    /// the triangle is stored.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The flat index of the combination of basis functions.
    auto _cell_index(const int    bra_angular_momentum,
                     const size_t bra_index,
                     const int    ket_angular_momentum,
                     const size_t ket_index) const -> size_t;

    /// @brief The index of atom basis on bra side among the unique atom bases
    /// of the molecular basis.
    int _bra_index;

    /// @brief The index of atom basis on ket side among the unique atom bases
    /// of the molecular basis.
    int _ket_index;

    /// @brief The atomic indices of atoms in diagonal atom pairs.
    std::vector<int> _atoms;

    /// @brief The storage layout of the diagonal atom pair blocks.
    diagstor _storage;

    /// @brief The flag indicating that only the triangle of the combinations of
    /// basis functions is stored.
    bool _triangular;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on bra side, with the total number of basis functions as last element.
    std::vector<size_t> _bra_offsets;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on ket side, with the total number of basis functions as last element.
    std::vector<size_t> _ket_offsets;

    /// @brief The offsets of the values of each stored combination of basis
    /// functions in the values block, with the total number of values as last
    /// element. The combinations are ordered row wise.
    std::vector<size_t> _value_offsets;
};

#endif /* AtomBasisDiagonalSparsity_hpp */

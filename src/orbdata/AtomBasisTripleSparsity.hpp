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


#ifndef AtomBasisTripleSparsity_hpp
#define AtomBasisTripleSparsity_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <vector>

#include "AtomBasis.hpp"
#include "AtomBasisGroup.hpp"
#include "AtomBasisPairGroup.hpp"
#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"

/// @brief Class CAtomBasisTripleSparsity stores the sparsity pattern of a
/// partially sparse block of three-center integrals, where the atom pairs on the
/// a and b sides are sparse while the atoms on the c side are not. It holds the
/// atom pairs of the a and b sides, the atoms of the c side, and the number of
/// atom pairs surviving screening for each combination of basis functions on the
/// a, b and c sides. The c side enters through the contraction and exponents of
/// its basis functions only, not through its distances, so its atoms are never
/// screened out.
/// @note The atom pairs of the diagonal atom pairs on the a and b sides are
/// placed at the front of the atom pair list, as their interatomic distance is
/// zero and they therefore precede every screened atom pair. Every combination
/// of basis functions thus retains a leading subrange of the atom pair list, and
/// the atom pairs retained by no combination are dropped.
/// @note The atom bases are referred to by their index among the unique atom
/// bases of the molecular basis they originate from, so the sparsity pattern is
/// a self contained value. The indices do not record which molecular basis they
/// refer to, and resolving an index against the wrong molecular basis addresses
/// the wrong atom basis instead of failing.
class CAtomBasisTripleSparsity
{
   public:
    /// @brief The constructor with atom basis pair group and atom basis group.
    /// All atom pairs are retained for every combination of basis functions,
    /// i.e. the sparsity pattern is dense until screening is applied.
    /// @param pair_group The atom basis pair group on a and b sides.
    /// @param group The atom basis group on c side.
    CAtomBasisTripleSparsity(const CAtomBasisPairGroup &pair_group, const CAtomBasisGroup &group);

    /// @brief The constructor with atom basis pair group, atom basis group,
    /// integral bound and screening threshold. Only the atom pairs whose
    /// integral bound exceeds the threshold are retained for each combination of
    /// basis functions.
    /// @param pair_group The atom basis pair group on a and b sides, ordered by
    /// ascending interatomic distance.
    /// @param group The atom basis group on c side.
    /// @param bound The integral bound, evaluated as bound(a_function,
    /// b_function, c_function, distance), decreasing monotonically with the
    /// distance between the atoms on the a and b sides.
    /// @param threshold The screening threshold.
    template <typename B>
    CAtomBasisTripleSparsity(const CAtomBasisPairGroup &pair_group, const CAtomBasisGroup &group, const B &bound, const double threshold)

        : _a_index(pair_group.bra_index())

        , _b_index(pair_group.ket_index())

        , _c_index(group.index())

        , _a_atoms(pair_group.diagonal_atoms())

        , _b_atoms(pair_group.diagonal_atoms())

        , _c_atoms(group.atoms())

        , _counts{}

        , _a_offsets(pair_group.bra_basis().basis_function_offsets())

        , _b_offsets(pair_group.ket_basis().basis_function_offsets())

        , _c_offsets(group.basis().basis_function_offsets())
    {
        errors::assertMsgCritical(pair_group.ordering() == pairord::distance,
                                  std::string("AtomBasisTripleSparsity: Atom pairs are not ordered by interatomic distance"));

        const auto ndiag = _a_atoms.size();

        const auto &a_functions = pair_group.bra_basis().functions();

        const auto &b_functions = pair_group.ket_basis().functions();

        const auto &c_functions = group.basis().functions();

        _counts.reserve(a_functions.size() * b_functions.size() * c_functions.size());

        // NOTE: the diagonal atom pairs are at zero interatomic distance, so they
        // survive any threshold and are counted in addition to the bisected
        // number of surviving off-diagonal atom pairs.

        std::ranges::for_each(a_functions, [&](const auto &a_function) {
            std::ranges::for_each(b_functions, [&](const auto &b_function) {
                std::ranges::for_each(c_functions, [&](const auto &c_function) {
                    _counts.push_back(ndiag + screenfunc::number_of_leading(pair_group.distances(), [&](const double r) {
                                          return bound(a_function, b_function, c_function, r) > threshold;
                                      }));
                });
            });
        });

        const auto nkept = _counts.empty() ? ndiag : std::ranges::max(_counts);

        _a_atoms.insert(_a_atoms.end(), pair_group.bra_atoms().begin(), pair_group.bra_atoms().begin() + (nkept - ndiag));

        _b_atoms.insert(_b_atoms.end(), pair_group.ket_atoms().begin(), pair_group.ket_atoms().begin() + (nkept - ndiag));

        _a_atoms.shrink_to_fit();

        _b_atoms.shrink_to_fit();
    }

    /// @brief Gets index of atom basis on a side among the unique atom bases of
    /// the molecular basis it originates from.
    /// @return The index of atom basis.
    auto
    a_index() const -> int
    {
        return _a_index;
    }

    /// @brief Gets index of atom basis on b side among the unique atom bases of
    /// the molecular basis it originates from.
    /// @return The index of atom basis.
    auto
    b_index() const -> int
    {
        return _b_index;
    }

    /// @brief Gets index of atom basis on c side among the unique atom bases of
    /// the molecular basis it originates from.
    /// @return The index of atom basis.
    auto
    c_index() const -> int
    {
        return _c_index;
    }

    /// @brief Gets atomic indices of atom pairs on a side.
    /// @return The constant reference to vector of atomic indices.
    auto
    a_atoms() const -> const std::vector<int> &
    {
        return _a_atoms;
    }

    /// @brief Gets atomic indices of atom pairs on b side.
    /// @return The constant reference to vector of atomic indices.
    auto
    b_atoms() const -> const std::vector<int> &
    {
        return _b_atoms;
    }

    /// @brief Gets atomic indices of atoms on c side.
    /// @return The constant reference to vector of atomic indices.
    auto
    c_atoms() const -> const std::vector<int> &
    {
        return _c_atoms;
    }

    /// @brief Gets number of atom pairs in sparsity pattern, i.e. the atom pairs
    /// retained by at least one combination of basis functions. This is the
    /// largest number of surviving atom pairs over all combinations of basis
    /// functions.
    /// @return The number of atom pairs.
    auto
    number_of_pairs() const -> size_t
    {
        return _a_atoms.size();
    }

    /// @brief Gets number of surviving atom pairs for specific combination of
    /// basis functions on a, b and c sides.
    /// @param a_angular_momentum The angular momentum of basis function on a side.
    /// @param a_index The index of basis function on a side.
    /// @param b_angular_momentum The angular momentum of basis function on b side.
    /// @param b_index The index of basis function on b side.
    /// @param c_angular_momentum The angular momentum of basis function on c side.
    /// @param c_index The index of basis function on c side.
    /// @return The number of surviving atom pairs.
    auto number_of_pairs(const int    a_angular_momentum,
                         const size_t a_index,
                         const int    b_angular_momentum,
                         const size_t b_index,
                         const int    c_angular_momentum,
                         const size_t c_index) const -> size_t;

    /// @brief Gets number of atoms on c side.
    /// @return The number of atoms.
    auto
    number_of_c_atoms() const -> size_t
    {
        return _c_atoms.size();
    }

    /// @brief Gets number of basis functions on a side.
    /// @return The number of basis functions.
    auto
    number_of_a_basis_functions() const -> size_t
    {
        return _a_offsets.back();
    }

    /// @brief Gets number of basis functions on b side.
    /// @return The number of basis functions.
    auto
    number_of_b_basis_functions() const -> size_t
    {
        return _b_offsets.back();
    }

    /// @brief Gets number of basis functions on c side.
    /// @return The number of basis functions.
    auto
    number_of_c_basis_functions() const -> size_t
    {
        return _c_offsets.back();
    }

   private:
    /// @brief The index of atom basis on a side.
    int _a_index;

    /// @brief The index of atom basis on b side.
    int _b_index;

    /// @brief The index of atom basis on c side.
    int _c_index;

    /// @brief The atomic indices of atom pairs on a side.
    std::vector<int> _a_atoms;

    /// @brief The atomic indices of atom pairs on b side.
    std::vector<int> _b_atoms;

    /// @brief The atomic indices of atoms on c side.
    std::vector<int> _c_atoms;

    /// @brief The number of surviving atom pairs for each combination of basis
    /// functions on a, b and c sides, stored with the a side as the slowest and
    /// the c side as the fastest running index.
    std::vector<size_t> _counts;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on a side, with the total number of basis functions as last element.
    std::vector<size_t> _a_offsets;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on b side, with the total number of basis functions as last element.
    std::vector<size_t> _b_offsets;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on c side, with the total number of basis functions as last element.
    std::vector<size_t> _c_offsets;
};

#endif /* AtomBasisTripleSparsity_hpp */

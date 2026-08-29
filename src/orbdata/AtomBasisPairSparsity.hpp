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


#ifndef AtomBasisPairSparsity_hpp
#define AtomBasisPairSparsity_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <vector>

#include "AtomBasis.hpp"
#include "AtomBasisPairGroup.hpp"
#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"

/// @brief Class CAtomBasisPairSparsity stores the sparsity pattern of the
/// off-diagonal atom pairs of an atom basis pair group, i.e. the atom pairs
/// themselves and the number of atom pairs surviving screening for each
/// combination of basis functions on bra and ket sides. Atom pairs of diagonal
/// atom pairs are not described here, as their storage requirements differ. The
/// atom pairs retained by no combination of basis functions are dropped, so the
/// atom pair list is as long as the largest number of surviving atom pairs. The
/// atom bases are referred to by their index among the unique atom bases of the
/// molecular basis they originate from, so the sparsity pattern is a self
/// contained value which outlives the atom bases it was constructed from.
/// @note The indices do not record which molecular basis they refer to. For a
/// sparsity pattern of a single molecular basis both indices refer to it, while
/// for a pair of molecular bases the bra index refers to the bra molecular basis
/// and the ket index to the ket molecular basis. Resolving an index against the
/// wrong molecular basis addresses the wrong atom basis instead of failing.
class CAtomBasisPairSparsity
{
   public:
    /// @brief The constructor with atom basis pair group. All atom pairs of the
    /// group are retained for every combination of basis functions, i.e. the
    /// sparsity pattern is dense until screening is applied.
    /// @param group The atom basis pair group to describe.
    explicit CAtomBasisPairSparsity(const CAtomBasisPairGroup &group);

    /// @brief The constructor with atom basis pair group, integral bound and
    /// screening threshold. Only the atom pairs whose integral bound exceeds
    /// the threshold are retained for each combination of basis functions.
    /// @param group The atom basis pair group to describe, ordered by ascending
    /// interatomic distance.
    /// @param bound The integral bound, evaluated as bound(bra_function,
    /// ket_function, distance), decreasing monotonically with distance.
    /// @param threshold The screening threshold.
    template <typename B>
    CAtomBasisPairSparsity(const CAtomBasisPairGroup &group, const B &bound, const double threshold)

        : _bra_index(group.bra_index())

        , _ket_index(group.ket_index())

        , _bra_atoms(group.bra_atoms())

        , _ket_atoms(group.ket_atoms())

        , _counts{}

        , _bra_offsets(group.bra_basis().basis_function_offsets())

        , _ket_offsets(group.ket_basis().basis_function_offsets())
    {
        // NOTE: bisection over the atom pairs requires them to be ordered by
        // ascending interatomic distance, as the integral bound decreases
        // monotonically with distance and surviving atom pairs thus form a
        // leading subrange of the atom pair list.

        errors::assertMsgCritical(group.ordering() == pairord::distance,
                                  std::string("AtomBasisPairSparsity: Atom pairs are not ordered by interatomic distance"));

        const auto &bra_functions = group.bra_basis().functions();

        const auto &ket_functions = group.ket_basis().functions();

        _counts.reserve(bra_functions.size() * ket_functions.size());

        std::ranges::for_each(bra_functions, [&](const auto &bra_function) {
            std::ranges::for_each(ket_functions, [&](const auto &ket_function) {
                _counts.push_back(screenfunc::number_of_leading(
                    group.distances(), [&](const double r) { return bound(bra_function, ket_function, r) > threshold; }));
            });
        });

        // NOTE: every combination of basis functions retains a leading subrange
        // of the atom pairs, so the atom pairs beyond the largest count are used
        // by no combination and are dropped to compress the sparsity pattern.

        const auto nkept = _counts.empty() ? size_t{0} : std::ranges::max(_counts);

        _bra_atoms.resize(nkept);

        _ket_atoms.resize(nkept);

        _bra_atoms.shrink_to_fit();

        _ket_atoms.shrink_to_fit();
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

    /// @brief Gets atomic indices of atom pairs on bra side.
    /// @return The constant reference to vector of atomic indices.
    auto
    bra_atoms() const -> const std::vector<int> &
    {
        return _bra_atoms;
    }

    /// @brief Gets atomic indices of atom pairs on ket side.
    /// @return The constant reference to vector of atomic indices.
    auto
    ket_atoms() const -> const std::vector<int> &
    {
        return _ket_atoms;
    }

    /// @brief Gets number of atom pairs in sparsity pattern, i.e. the atom pairs
    /// retained by at least one combination of basis functions. This is the
    /// largest number of surviving atom pairs over all combinations of basis
    /// functions, and zero if the sparsity pattern is screened out entirely.
    /// @return The number of atom pairs.
    auto
    number_of_pairs() const -> size_t
    {
        return _bra_atoms.size();
    }

    /// @brief Gets number of surviving atom pairs for specific combination of
    /// basis functions on bra and ket sides.
    /// @param bra_angular_momentum The angular momentum of basis function on bra side.
    /// @param bra_index The index of basis function on bra side.
    /// @param ket_angular_momentum The angular momentum of basis function on ket side.
    /// @param ket_index The index of basis function on ket side.
    /// @return The number of surviving atom pairs.
    auto number_of_pairs(const int    bra_angular_momentum,
                         const size_t bra_index,
                         const int    ket_angular_momentum,
                         const size_t ket_index) const -> size_t;

    /// @brief Gets number of basis functions on bra side.
    /// @return The number of basis functions.
    auto
    number_of_bra_basis_functions() const -> size_t
    {
        // NOTE: the last offset is the total number of basis functions.

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
    /// @brief The index of atom basis on bra side among the unique atom bases
    /// of the molecular basis.
    int _bra_index;

    /// @brief The index of atom basis on ket side among the unique atom bases
    /// of the molecular basis.
    int _ket_index;

    /// @brief The atomic indices of atom pairs on bra side.
    std::vector<int> _bra_atoms;

    /// @brief The atomic indices of atom pairs on ket side.
    std::vector<int> _ket_atoms;

    /// @brief The number of surviving atom pairs for each combination of basis
    /// functions on bra and ket sides, stored row wise with basis functions on
    /// bra side as rows.
    std::vector<size_t> _counts;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on bra side, with the total number of basis functions as last element.
    std::vector<size_t> _bra_offsets;

    /// @brief The offsets of the first basis function of each angular momentum
    /// on ket side, with the total number of basis functions as last element.
    std::vector<size_t> _ket_offsets;
};

#endif /* AtomBasisPairSparsity_hpp */

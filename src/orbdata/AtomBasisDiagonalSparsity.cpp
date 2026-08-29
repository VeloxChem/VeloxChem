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


#include "AtomBasisDiagonalSparsity.hpp"

#include <algorithm>
#include <ranges>

#include "ErrorHandler.hpp"
#include "TensorComponents.hpp"

CAtomBasisDiagonalSparsity::CAtomBasisDiagonalSparsity(const CAtomBasisPairGroup &group, const diagstor storage)

    : _bra_index(group.bra_index())

    , _ket_index(group.ket_index())

    , _atoms(group.diagonal_atoms())

    , _storage(storage)

    , _triangular((&group.bra_basis() == &group.ket_basis()) || (group.bra_basis() == group.ket_basis()))

    , _bra_offsets(group.bra_basis().basis_function_offsets())

    , _ket_offsets(group.ket_basis().basis_function_offsets())

    , _value_offsets{}
{
    // NOTE: the angular momentum of a basis function follows from the offsets,
    // as the basis functions are kept sorted by ascending angular momentum.

    const auto momenta = [](const std::vector<size_t> &offsets) {
        std::vector<int> moments(offsets.back(), 0);

        std::ranges::for_each(std::views::iota(size_t{0}, offsets.size() - 1), [&](const auto i) {
            std::ranges::fill(moments.begin() + offsets[i], moments.begin() + offsets[i + 1], static_cast<int>(i));
        });

        return moments;
    };

    const auto bra_moments = momenta(_bra_offsets);

    const auto ket_moments = momenta(_ket_offsets);

    const auto nbra = bra_moments.size();

    const auto nket = ket_moments.size();

    // NOTE: only the triangle of the combinations is stored when bra and ket
    // sides carry the same atom basis, as the block is then symmetric. For a
    // diagonal atom pair spanning two molecular bases the block is not symmetric
    // and all combinations are stored.

    _value_offsets.reserve((_triangular ? nbra * (nbra + 1) / 2 : nbra * nket) + 1);

    _value_offsets.push_back(0);

    std::ranges::for_each(std::views::iota(size_t{0}, nbra), [&](const auto i) {
        std::ranges::for_each(std::views::iota(_triangular ? i : size_t{0}, nket), [&](const auto j) {
            const auto ncomps = tensor::number_of_spherical_components(std::array<int, 2>{bra_moments[i], ket_moments[j]});

            const auto nvals = (_storage == diagstor::scalar) ? ((bra_moments[i] == ket_moments[j]) ? size_t{1} : size_t{0})
                                                              : _atoms.size() * static_cast<size_t>(ncomps);

            _value_offsets.push_back(_value_offsets.back() + nvals);
        });
    });
}

auto
CAtomBasisDiagonalSparsity::_cell_index(const int    bra_angular_momentum,
                                        const size_t bra_index,
                                        const int    ket_angular_momentum,
                                        const size_t ket_index) const -> size_t
{
    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum + 2 <= static_cast<int>(_bra_offsets.size())) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum + 2 <= static_cast<int>(_ket_offsets.size())),
                              std::string("AtomBasisDiagonalSparsity: Angular momentum is out of range"));

    auto bra_row = _bra_offsets[bra_angular_momentum] + bra_index;

    auto ket_col = _ket_offsets[ket_angular_momentum] + ket_index;

    errors::assertMsgCritical((bra_row < _bra_offsets[bra_angular_momentum + 1]) && (ket_col < _ket_offsets[ket_angular_momentum + 1]),
                              std::string("AtomBasisDiagonalSparsity: Index of basis function is out of range"));

    if (!_triangular) return bra_row * number_of_ket_basis_functions() + ket_col;

    // NOTE: the combination is swapped into the stored triangle if given in the
    // reverse order, as the diagonal atom pair block is then symmetric.

    if (bra_row > ket_col) std::swap(bra_row, ket_col);

    const auto nbfs = number_of_bra_basis_functions();

    return bra_row * nbfs - bra_row * (bra_row + 1) / 2 + ket_col;
}

auto
CAtomBasisDiagonalSparsity::number_of_elements(const int    bra_angular_momentum,
                                               const size_t bra_index,
                                               const int    ket_angular_momentum,
                                               const size_t ket_index) const -> size_t
{
    const auto index = _cell_index(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index);

    return _value_offsets[index + 1] - _value_offsets[index];
}

auto
CAtomBasisDiagonalSparsity::element_offset(const int    bra_angular_momentum,
                                           const size_t bra_index,
                                           const int    ket_angular_momentum,
                                           const size_t ket_index) const -> size_t
{
    return _value_offsets[_cell_index(bra_angular_momentum, bra_index, ket_angular_momentum, ket_index)];
}

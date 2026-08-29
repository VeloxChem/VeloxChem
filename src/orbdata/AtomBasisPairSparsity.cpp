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


#include "AtomBasisPairSparsity.hpp"

#include <algorithm>
#include <ranges>

#include "ErrorHandler.hpp"

auto
CAtomBasisPairSparsity::_make_offsets(const CAtomBasis &basis) -> std::vector<size_t>
{
    std::vector<size_t> offsets(static_cast<size_t>(basis.max_angular_momentum() + 2), 0);

    std::ranges::for_each(std::views::iota(0, basis.max_angular_momentum() + 1),
                          [&](const int i) { offsets[i + 1] = offsets[i] + basis.number_of_basis_functions(i); });

    return offsets;
}

CAtomBasisPairSparsity::CAtomBasisPairSparsity(const CAtomBasisPairGroup &group)

    : _bra_index(group.bra_index())

    , _ket_index(group.ket_index())

    , _bra_atoms(group.bra_atoms())

    , _ket_atoms(group.ket_atoms())

    , _counts{}

    , _bra_offsets(_make_offsets(group.bra_basis()))

    , _ket_offsets(_make_offsets(group.ket_basis()))
{
    // NOTE: all atom pairs of the group survive for every combination of basis
    // functions until screening is applied, which can only lower the counts.

    _counts.assign(number_of_bra_basis_functions() * number_of_ket_basis_functions(), _bra_atoms.size());
}

auto
CAtomBasisPairSparsity::number_of_pairs(const int    bra_angular_momentum,
                                        const size_t bra_index,
                                        const int    ket_angular_momentum,
                                        const size_t ket_index) const -> size_t
{
    errors::assertMsgCritical((bra_angular_momentum >= 0) && (bra_angular_momentum + 2 <= static_cast<int>(_bra_offsets.size())) &&
                                  (ket_angular_momentum >= 0) && (ket_angular_momentum + 2 <= static_cast<int>(_ket_offsets.size())),
                              std::string("AtomBasisPairSparsity.number_of_pairs: Angular momentum is out of range"));

    const auto bra_row = _bra_offsets[bra_angular_momentum] + bra_index;

    const auto ket_col = _ket_offsets[ket_angular_momentum] + ket_index;

    errors::assertMsgCritical((bra_row < _bra_offsets[bra_angular_momentum + 1]) && (ket_col < _ket_offsets[ket_angular_momentum + 1]),
                              std::string("AtomBasisPairSparsity.number_of_pairs: Index of basis function is out of range"));

    return _counts[bra_row * number_of_ket_basis_functions() + ket_col];
}

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


#include "AtomBasisTripleSparsity.hpp"

#include <algorithm>
#include <ranges>

#include "ErrorHandler.hpp"

CAtomBasisTripleSparsity::CAtomBasisTripleSparsity(const CAtomBasisPairGroup &pair_group, const CAtomBasisGroup &group)

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
    // NOTE: the diagonal atom pairs precede the off-diagonal atom pairs, as their
    // interatomic distance is zero. All atom pairs survive for every combination
    // of basis functions until screening is applied.

    _a_atoms.insert(_a_atoms.end(), pair_group.bra_atoms().begin(), pair_group.bra_atoms().end());

    _b_atoms.insert(_b_atoms.end(), pair_group.ket_atoms().begin(), pair_group.ket_atoms().end());

    _counts.assign(number_of_a_basis_functions() * number_of_b_basis_functions() * number_of_c_basis_functions(), _a_atoms.size());
}

auto
CAtomBasisTripleSparsity::number_of_pairs(const int    a_angular_momentum,
                                          const size_t a_index,
                                          const int    b_angular_momentum,
                                          const size_t b_index,
                                          const int    c_angular_momentum,
                                          const size_t c_index) const -> size_t
{
    errors::assertMsgCritical((a_angular_momentum >= 0) && (a_angular_momentum + 2 <= static_cast<int>(_a_offsets.size())) &&
                                  (b_angular_momentum >= 0) && (b_angular_momentum + 2 <= static_cast<int>(_b_offsets.size())) &&
                                  (c_angular_momentum >= 0) && (c_angular_momentum + 2 <= static_cast<int>(_c_offsets.size())),
                              std::string("AtomBasisTripleSparsity.number_of_pairs: Angular momentum is out of range"));

    const auto arow = _a_offsets[a_angular_momentum] + a_index;

    const auto brow = _b_offsets[b_angular_momentum] + b_index;

    const auto ccol = _c_offsets[c_angular_momentum] + c_index;

    errors::assertMsgCritical((arow < _a_offsets[a_angular_momentum + 1]) && (brow < _b_offsets[b_angular_momentum + 1]) &&
                                  (ccol < _c_offsets[c_angular_momentum + 1]),
                              std::string("AtomBasisTripleSparsity.number_of_pairs: Index of basis function is out of range"));

    return _counts[(arow * number_of_b_basis_functions() + brow) * number_of_c_basis_functions() + ccol];
}

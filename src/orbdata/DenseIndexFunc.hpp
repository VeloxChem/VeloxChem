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



#ifndef DenseIndexFunc_hpp
#define DenseIndexFunc_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "AtomBasis.hpp"
#include "MolecularBasis.hpp"

namespace denseidx {  // denseidx namespace

/// @brief Gets the angular momentum and the index within it of each basis
/// function of an atom basis.
/// @param basis The atom basis to index the basis functions of.
/// @return The vector of angular momenta and indices within them.
inline auto
index_functions(const CAtomBasis &basis) -> std::vector<std::pair<int, size_t>>
{
    std::vector<std::pair<int, size_t>> indices;

    std::vector<size_t> counts(static_cast<size_t>(basis.max_angular_momentum() + 1), 0);

    const auto &functions = basis.functions();

    for (size_t i = 0; i < functions.size(); i++)
    {
        const auto lval = functions[i].get_angular_momentum();

        indices.push_back({lval, counts[lval]});

        counts[lval]++;
    }

    return indices;
}

/// @brief Creates the dense index of the first angular component of the first
/// basis function of each angular momentum of each atom.
/// @param basis The molecular basis to index the basis functions of.
/// @return The vector of dense indices, with the atom as the slowest running
/// index and the angular momentum as the fastest.
inline auto
make_dense_starts(const CMolecularBasis &basis) -> std::vector<size_t>
{
    const auto indices = basis.basis_sets_indices();

    const auto nmoms = static_cast<size_t>(basis.max_angular_momentum() + 1);

    std::vector<size_t> starts(indices.size() * nmoms, 0);

    // NOTE: the basis functions of an angular momentum are numbered from the
    // dimension of the basis functions of the lower angular momenta, as the
    // dense matrix is ordered by ascending angular momentum.

    std::vector<size_t> counts(nmoms, 0);

    for (size_t l = 0; l < nmoms; l++)
    {
        counts[l] = basis.dimensions_of_basis(static_cast<int>(l));
    }

    for (size_t i = 0; i < indices.size(); i++)
    {
        const auto &abas = basis.basis_set(indices[i]);

        for (size_t l = 0; l < nmoms; l++)
        {
            starts[i * nmoms + l] = counts[l];

            counts[l] += abas.number_of_basis_functions(static_cast<int>(l));
        }
    }

    return starts;
}

/// @brief Creates the distance between the dense indices of two consecutive
/// angular components of a basis function of each angular momentum.
/// @param basis The molecular basis to index the basis functions of.
/// @return The vector of distances.
inline auto
make_dense_strides(const CMolecularBasis &basis) -> std::vector<size_t>
{
    const auto nmoms = static_cast<size_t>(basis.max_angular_momentum() + 1);

    std::vector<size_t> strides(nmoms, 0);

    for (size_t l = 0; l < nmoms; l++)
    {
        strides[l] = basis.number_of_basis_functions(static_cast<int>(l));
    }

    return strides;
}

}  // namespace denseidx

#endif /* DenseIndexFunc_hpp */

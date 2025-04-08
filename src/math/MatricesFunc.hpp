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

#ifndef MatricesFunc_hpp
#define MatricesFunc_hpp

#include <algorithm>
#include <array>
#include <ranges>
#include <string>

#include "CustomViews.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "MatrixFunc.hpp"
#include "MolecularBasis.hpp"
#include "TensorLabels.hpp"

namespace matfunc {

/// @brief Creates matrices for given tensor orders.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @param basis The molecular basis to create matrices.
/// @param mtype  The matrix type of created matrices.
/// @return The matrices.
template <size_t N>
auto
make_matrices(const std::array<int, N>& orders, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices
{
    CMatrices matrices;

    const auto matrix = matfunc::make_matrix(basis, mtype);

    auto labels = tensor::cartesian_labels(orders[0]);

    std::ranges::for_each(std::views::iota(size_t{1}, N), [&](const auto i) {
        std::vector<std::string> new_labels;
        auto                     clabels = tensor::cartesian_labels(orders[i]);
        std::ranges::for_each(views::rectangular(labels.size(), clabels.size()),
                              [&](const auto& index) { new_labels.push_back(labels[index.first] + "_" + clabels[index.second]); });
        labels = new_labels;
    });

    std::ranges::for_each(labels, [&](const auto& label) { matrices.add(matrix, label); });

    return matrices;
}

/// @brief Creates matrices for given tensor orders.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @param bra_basis The molecular basis to create matrix (rows data).
/// @param ket_basis The molecular basis to create matrix (columns data).
/// @return The matrices.
template <size_t N>
auto
make_matrices(const std::array<int, N>& orders, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices
{
    CMatrices matrices;

    const auto matrix = matfunc::make_matrix(bra_basis, ket_basis);

    auto labels = tensor::cartesian_labels(orders[0]);

    std::ranges::for_each(std::views::iota(size_t{1}, N), [&](const auto i) {
        std::vector<std::string> new_labels;
        auto                     clabels = tensor::cartesian_labels(orders[i]);
        std::ranges::for_each(views::rectangular(labels.size(), clabels.size()),
                              [&](const auto& index) { new_labels.push_back(labels[index.first] + "_" + clabels[index.second]); });
        labels = new_labels;
    });

    std::ranges::for_each(labels, [&](const auto& label) { matrices.add(matrix, label); });

    return matrices;
}

}  // namespace matfunc

#endif /* MatricesFunc_hpp */

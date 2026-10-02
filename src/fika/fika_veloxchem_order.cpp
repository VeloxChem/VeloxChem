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

#include "fika_veloxchem_order.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace fika {

auto veloxchem_order(const MolecularBasis& basis) -> std::vector<std::size_t> {
  int max_l = -1;
  for (const MolecularShell& shell : basis.shells()) {
    max_l = std::max(max_l, shell.angular_momentum);
  }
  std::vector<std::size_t> order(basis.function_count());
  std::size_t next = 0;
  for (int l = 0; l <= max_l; ++l) {
    for (int m = -l; m <= l; ++m) {
      // shells() runs over the atoms in order and each atom's shells in basis order.
      for (const MolecularShell& shell : basis.shells()) {
        if (shell.angular_momentum != l) {
          continue;
        }
        for (std::size_t k = 0; k < shell.contraction_count; ++k) {
          order[shell.function_offset + k * static_cast<std::size_t>(2 * l + 1) +
                static_cast<std::size_t>(m + l)] = next++;
        }
      }
    }
  }
  return order;
}

namespace {

void check_size(const DenseMatrix& matrix, const MolecularBasis& basis, const std::string& name) {
  const std::size_t n = basis.function_count();
  if (matrix.rows() != n || matrix.columns() != n) {
    throw std::invalid_argument("fika::" + name + ": matrix is " + std::to_string(matrix.rows()) +
                                " x " + std::to_string(matrix.columns()) + " for " +
                                std::to_string(n) + " basis functions");
  }
}

/// result(i, j) = matrix(order[i], order[j]), of the same symmetry.
auto permuted(const DenseMatrix& matrix, const std::vector<std::size_t>& order) -> DenseMatrix {
  const std::size_t n = order.size();
  DenseMatrix result(n, n, matrix.symmetry());
  const bool full = matrix.symmetry() == MatrixSymmetry::general;
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t end = full ? n : i + 1;
    for (std::size_t j = 0; j < end; ++j) {
      result.set(i, j, matrix(order[i], order[j]));
    }
  }
  return result;
}

}  // namespace

auto veloxchem_to_fika(const DenseMatrix& matrix, const MolecularBasis& basis) -> DenseMatrix {
  check_size(matrix, basis, "veloxchem_to_fika");
  return permuted(matrix, veloxchem_order(basis));
}

auto fika_to_veloxchem(const DenseMatrix& matrix, const MolecularBasis& basis) -> DenseMatrix {
  check_size(matrix, basis, "fika_to_veloxchem");
  const auto order = veloxchem_order(basis);
  std::vector<std::size_t> inverse(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    inverse[order[i]] = i;
  }
  return permuted(matrix, inverse);
}

}  // namespace fika

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

#include "fika_dense_matrix.hpp"

#include <stdexcept>
#include <string>

namespace fika {

namespace {

auto stored_count(std::size_t rows, std::size_t columns, MatrixSymmetry symmetry) -> std::size_t {
  switch (symmetry) {
    case MatrixSymmetry::general:
      return rows * columns;
    case MatrixSymmetry::symmetric:
      return rows * (rows + 1) / 2;
    case MatrixSymmetry::antisymmetric:
      return rows == 0 ? 0 : rows * (rows - 1) / 2;
  }
  return 0;
}

}  // namespace

DenseMatrix::DenseMatrix(std::size_t rows, std::size_t columns, MatrixSymmetry symmetry)
    : rows_(rows), columns_(columns), symmetry_(symmetry) {
  if (symmetry != MatrixSymmetry::general && rows != columns) {
    throw std::invalid_argument("fika::DenseMatrix: a " + std::to_string(rows) + " x " +
                                std::to_string(columns) +
                                " matrix cannot be symmetric or antisymmetric");
  }
  values_.assign(stored_count(rows, columns, symmetry), 0.0);
}

auto DenseMatrix::to_full() const -> std::vector<double> {
  if (symmetry_ == MatrixSymmetry::general) {
    return values_;
  }
  std::vector<double> full(rows_ * columns_);
  const double sign = symmetry_ == MatrixSymmetry::symmetric ? 1.0 : -1.0;
  const std::size_t n = rows_;
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t end = symmetry_ == MatrixSymmetry::symmetric ? i + 1 : i;  // j < end
    for (std::size_t j = 0; j < end; ++j) {
      full[i * n + j] = values_[index];
      full[j * n + i] = sign * values_[index];
      ++index;
    }
  }
  return full;
}

}  // namespace fika

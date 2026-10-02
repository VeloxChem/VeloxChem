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

#ifndef fika_dense_matrix_hpp
#define fika_dense_matrix_hpp

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

#include "fika_matrix_symmetry.hpp"

namespace fika {

/// Dense matrix of doubles, stored by symmetry:
///  - general (rows x columns): full, row-major, element (i, j) at i columns + j;
///  - symmetric (n x n): the lower triangle j <= i row by row, (i, j) at i (i + 1) / 2 + j,
///    n (n + 1) / 2 values (LAPACK packed storage with UPLO = 'U');
///  - antisymmetric (n x n): the strictly lower triangle j < i row by row, (i, j) at
///    i (i - 1) / 2 + j, n (n - 1) / 2 values (the diagonal is zero).
/// Elements of the other triangle read and write through the stored one with sign + (symmetric)
/// or - (antisymmetric).
class DenseMatrix {
 public:
  /// Zero matrix; symmetric and antisymmetric matrices must be square. Throws
  /// std::invalid_argument otherwise.
  DenseMatrix(std::size_t rows, std::size_t columns,
              MatrixSymmetry symmetry = MatrixSymmetry::general);

  /// Zero n x n matrix.
  DenseMatrix(std::size_t n, MatrixSymmetry symmetry) : DenseMatrix(n, n, symmetry) {}

  auto rows() const noexcept -> std::size_t { return rows_; }
  auto columns() const noexcept -> std::size_t { return columns_; }
  auto symmetry() const noexcept -> MatrixSymmetry { return symmetry_; }

  /// Number of stored values.
  auto size() const noexcept -> std::size_t { return values_.size(); }

  auto values() const noexcept -> std::span<const double> { return values_; }
  auto values() noexcept -> std::span<double> { return values_; }

  /// Element (i, j), from either triangle (0 on the diagonal of an antisymmetric matrix).
  auto operator()(std::size_t i, std::size_t j) const noexcept -> double {
    assert(i < rows_ && j < columns_);
    switch (symmetry_) {
      case MatrixSymmetry::general:
        return values_[i * columns_ + j];
      case MatrixSymmetry::symmetric:
        return i >= j ? values_[i * (i + 1) / 2 + j] : values_[j * (j + 1) / 2 + i];
      case MatrixSymmetry::antisymmetric:
        return i > j ? values_[i * (i - 1) / 2 + j] : i < j ? -values_[j * (j - 1) / 2 + i] : 0.0;
    }
    return 0.0;
  }

  /// Sets element (i, j) (and its mirror); the antisymmetric diagonal only takes 0.
  void set(std::size_t i, std::size_t j, double value) noexcept {
    if (double* stored = element(i, j, value)) {
      *stored = value;
    }
  }

  /// Adds to element (i, j) (and its mirror); the antisymmetric diagonal only takes 0.
  void add(std::size_t i, std::size_t j, double value) noexcept {
    if (double* stored = element(i, j, value)) {
      *stored += value;
    }
  }

  /// Full row-major matrix (rows x columns), mirrored triangle filled in.
  auto to_full() const -> std::vector<double>;

 private:
  /// Stored value of (i, j), with `value` negated when (i, j) is the upper triangle of an
  /// antisymmetric matrix; nullptr for its diagonal.
  auto element(std::size_t i, std::size_t j, double& value) noexcept -> double* {
    assert(i < rows_ && j < columns_);
    switch (symmetry_) {
      case MatrixSymmetry::general:
        return &values_[i * columns_ + j];
      case MatrixSymmetry::symmetric:
        return i >= j ? &values_[i * (i + 1) / 2 + j] : &values_[j * (j + 1) / 2 + i];
      case MatrixSymmetry::antisymmetric:
        if (i == j) {
          assert(value == 0.0);
          return nullptr;
        }
        if (i < j) {
          value = -value;
          return &values_[j * (j - 1) / 2 + i];
        }
        return &values_[i * (i - 1) / 2 + j];
    }
    return nullptr;
  }

  std::size_t rows_;
  std::size_t columns_;
  MatrixSymmetry symmetry_;
  std::vector<double> values_;
};

}  // namespace fika

#endif  // fika_dense_matrix_hpp

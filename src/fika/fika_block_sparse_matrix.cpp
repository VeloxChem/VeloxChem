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

#include "fika_block_sparse_matrix.hpp"

#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::BlockSparseMatrix: " + reason);
}

auto components(const MolecularShell& shell) -> std::size_t {
  return static_cast<std::size_t>(2 * shell.angular_momentum + 1);
}

auto is_mirrored(MatrixSymmetry symmetry) -> bool {
  return symmetry != MatrixSymmetry::general;
}

auto function_count(std::span<const MolecularShell> shells) -> std::size_t {
  return shells.empty() ? 0 : shells.back().function_offset + shells.back().function_count();
}

/// Number of index pairs i, j < n stored for a diagonal block: all (general), i <= j (symmetric)
/// or i < j (antisymmetric).
auto diagonal_pair_count(std::size_t n, MatrixSymmetry symmetry) -> std::size_t {
  switch (symmetry) {
    case MatrixSymmetry::general:
      return n * n;
    case MatrixSymmetry::symmetric:
      return n * (n + 1) / 2;
    case MatrixSymmetry::antisymmetric:
      return n * (n - 1) / 2;
  }
  return 0;
}

/// First stored column of row i in a diagonal block of the given symmetry.
auto first_column(std::size_t i, MatrixSymmetry symmetry) -> std::size_t {
  switch (symmetry) {
    case MatrixSymmetry::general:
      return 0;
    case MatrixSymmetry::symmetric:
      return i;
    case MatrixSymmetry::antisymmetric:
      return i + 1;
  }
  return 0;
}

auto square_shape(std::span<const MolecularShell> shells, MatrixSymmetry symmetry,
                  DiagonalFormat diagonal_format) -> detail::MatrixShape {
  return {{shells.begin(), shells.end()}, {}, true, symmetry, diagonal_format};
}

auto rectangular_shape(std::span<const MolecularShell> bra_shells,
                       std::span<const MolecularShell> ket_shells) -> detail::MatrixShape {
  return {{bra_shells.begin(), bra_shells.end()},
          {ket_shells.begin(), ket_shells.end()},
          false,
          MatrixSymmetry::general,
          DiagonalFormat::full};
}

// Dense element count from which the conversions below run in parallel.
constexpr std::size_t parallel_dense_elements = std::size_t{1} << 16;

/// Calls element(value_index, row, column) for every dense element a stored value represents,
/// in its stored (upper) position; scalar diagonal values are visited once per m. Every call has
/// its own value and (row, column), so bra shells run in parallel for large patterns.
template <typename Element>
void for_each_stored_element(const SparsityPattern& pattern, Element&& element) {
  const auto bra_shells = pattern.bra_shells();
  const auto ket_shells = pattern.ket_shells();
  const MatrixSymmetry symmetry = pattern.symmetry();
  const bool parallel = pattern.row_count() * pattern.column_count() >= parallel_dense_elements;
#pragma omp parallel for schedule(dynamic, 8) if (parallel)
  for (std::size_t bra = 0; bra < bra_shells.size(); ++bra) {
    const MolecularShell& a = bra_shells[bra];
    for (std::size_t block = pattern.row_offsets()[bra]; block < pattern.row_offsets()[bra + 1];
         ++block) {
      const MolecularShell& b = ket_shells[pattern.block_ket(block)];
      const BlockKind kind = block_kind(a, b, bra == pattern.block_ket(block), pattern.square());
      const bool scalar = pattern.diagonal_format() == DiagonalFormat::scalar;
      std::size_t value = pattern.block_value_offset(block);

      if (kind != BlockKind::dense && scalar) {
        const std::size_t n = components(a);
        const bool same_shell = kind == BlockKind::same_shell;
        for (std::size_t k_a = 0; k_a < a.contraction_count; ++k_a) {
          for (std::size_t k_b = same_shell ? first_column(k_a, symmetry) : 0;
               k_b < b.contraction_count; ++k_b) {
            for (std::size_t m = 0; m < n; ++m) {
              element(value, a.function_offset + k_a * n + m, b.function_offset + k_b * n + m);
            }
            ++value;
          }
        }
      } else if (kind == BlockKind::same_shell) {
        const std::size_t n = a.function_count();
        for (std::size_t i = 0; i < n; ++i) {
          for (std::size_t j = first_column(i, symmetry); j < n; ++j) {
            element(value++, a.function_offset + i, a.function_offset + j);
          }
        }
      } else {
        for (std::size_t i = 0; i < a.function_count(); ++i) {
          for (std::size_t j = 0; j < b.function_count(); ++j) {
            element(value++, a.function_offset + i, b.function_offset + j);
          }
        }
      }
    }
  }
}

}  // namespace

auto block_kind(const MolecularShell& bra, const MolecularShell& ket, bool same_shell,
                bool square) noexcept -> BlockKind {
  if (!square) {
    return BlockKind::dense;
  }
  if (same_shell) {
    return BlockKind::same_shell;
  }
  return bra.atom == ket.atom && bra.angular_momentum == ket.angular_momentum
             ? BlockKind::same_centre
             : BlockKind::dense;
}

auto block_value_count(const MolecularShell& bra, const MolecularShell& ket, BlockKind kind,
                       MatrixSymmetry symmetry, DiagonalFormat format) noexcept -> std::size_t {
  const bool scalar = format == DiagonalFormat::scalar;
  switch (kind) {
    case BlockKind::dense:
      break;
    case BlockKind::same_shell:
      return diagonal_pair_count(scalar ? bra.contraction_count : bra.function_count(), symmetry);
    case BlockKind::same_centre:
      if (scalar) {
        return bra.contraction_count * ket.contraction_count;
      }
      break;
  }
  return bra.function_count() * ket.function_count();
}

auto SparsityPattern::row_count() const noexcept -> std::size_t {
  return function_count(shape_.bra_shells);
}

auto SparsityPattern::column_count() const noexcept -> std::size_t {
  return function_count(shape_.kets());
}

auto SparsityPattern::from_rows(std::span<const MolecularShell> shells, MatrixSymmetry symmetry,
                                DiagonalFormat diagonal_format,
                                std::vector<std::size_t> row_offsets,
                                std::vector<std::uint32_t> block_kets)
    -> std::shared_ptr<const SparsityPattern> {
  return from_rows(square_shape(shells, symmetry, diagonal_format), std::move(row_offsets),
                   std::move(block_kets));
}

auto SparsityPattern::from_rows(std::span<const MolecularShell> bra_shells,
                                std::span<const MolecularShell> ket_shells,
                                std::vector<std::size_t> row_offsets,
                                std::vector<std::uint32_t> block_kets)
    -> std::shared_ptr<const SparsityPattern> {
  return from_rows(rectangular_shape(bra_shells, ket_shells), std::move(row_offsets),
                   std::move(block_kets));
}

auto SparsityPattern::from_rows(detail::MatrixShape shape, std::vector<std::size_t> row_offsets,
                                std::vector<std::uint32_t> block_kets)
    -> std::shared_ptr<const SparsityPattern> {
  const std::size_t rows = shape.bra_shells.size();
  if (row_offsets.size() != rows + 1 || row_offsets.front() != 0 ||
      row_offsets.back() != block_kets.size() || !std::ranges::is_sorted(row_offsets)) {
    fail("row offsets must increase from 0 to the block count, one per bra shell and one more");
  }
  // Validate each row and store each block's value count, then turn the counts into offsets
  // row after row.
  auto pattern = std::shared_ptr<SparsityPattern>(new SparsityPattern());
  std::vector<std::size_t>& offsets = pattern->block_value_offsets_;
  offsets.resize(block_kets.size() + 1);
  const bool mirrored = is_mirrored(shape.symmetry);
  std::vector<std::size_t> row_values(rows + 1, 0);
  bool invalid = false;
#pragma omp parallel for schedule(dynamic, 64) reduction(|| : invalid)
  for (std::size_t bra = 0; bra < rows; ++bra) {
    std::size_t total = 0;
    for (std::size_t block = row_offsets[bra]; block < row_offsets[bra + 1]; ++block) {
      const std::size_t ket = block_kets[block];
      if (ket >= shape.kets().size() || (mirrored && ket < bra) ||
          (block > row_offsets[bra] && ket <= block_kets[block - 1])) {
        invalid = true;
        break;
      }
      const MolecularShell& a = shape.bra_shells[bra];
      const MolecularShell& b = shape.kets()[ket];
      offsets[block] = fika::block_value_count(a, b, block_kind(a, b, bra == ket, shape.square),
                                               shape.symmetry, shape.diagonal_format);
      total += offsets[block];
    }
    row_values[bra + 1] = total;
  }
  if (invalid) {
    fail(mirrored ? "ket shells must increase within a row, exist and not precede the bra shell"
                  : "ket shells must increase within a row and exist");
  }
  for (std::size_t bra = 0; bra < rows; ++bra) {
    row_values[bra + 1] += row_values[bra];
  }
  offsets.back() = row_values[rows];
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t bra = 0; bra < rows; ++bra) {
    std::size_t offset = row_values[bra];
    for (std::size_t block = row_offsets[bra]; block < row_offsets[bra + 1]; ++block) {
      offset += std::exchange(offsets[block], offset);
    }
  }
  pattern->shape_ = std::move(shape);
  pattern->row_offsets_ = std::move(row_offsets);
  pattern->block_ket_ = std::move(block_kets);
  return pattern;
}

auto SparsityPattern::find_block(std::size_t bra_start, std::size_t ket_start) const
    -> std::optional<std::size_t> {
  const auto shell_index = [](std::span<const MolecularShell> shells, std::size_t start) {
    const auto shell =
        std::ranges::lower_bound(shells, start, {}, &MolecularShell::function_offset);
    if (shell == shells.end() || shell->function_offset != start) {
      throw std::invalid_argument("fika::SparsityPattern: no shell starts at basis function " +
                                  std::to_string(start));
    }
    return static_cast<std::size_t>(shell - shells.begin());
  };
  const std::size_t bra = shell_index(shape_.bra_shells, bra_start);
  const std::size_t ket = shell_index(shape_.kets(), ket_start);
  const auto first = block_ket_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[bra]);
  const auto last = block_ket_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[bra + 1]);
  const auto block = std::lower_bound(first, last, ket);
  if (block == last || *block != ket) {
    return std::nullopt;
  }
  return static_cast<std::size_t>(block - block_ket_.begin());
}

BlockSparseMatrix::BlockSparseMatrix(std::shared_ptr<const SparsityPattern> pattern)
    : pattern_(std::move(pattern)) {
  // Allocate without initializing, then zero in parallel: pages are first touched by all
  // threads instead of one (large matrices spend most of their setup here otherwise).
  values_.resize(pattern_->value_count());
  double* values = values_.data();
  const std::size_t count = values_.size();
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = 0.0;
  }
}

auto BlockSparseMatrix::from_dense(std::shared_ptr<const SparsityPattern> pattern,
                                   std::span<const double> dense) -> BlockSparseMatrix {
  const std::size_t rows = pattern->row_count();
  const std::size_t columns = pattern->column_count();
  if (dense.size() != rows * columns) {
    fail("dense matrix of " + std::to_string(dense.size()) + " elements, expected " +
         std::to_string(rows * columns));
  }
  BlockSparseMatrix matrix(std::move(pattern));
  for_each_stored_element(matrix.pattern(),
                          [&](std::size_t value, std::size_t row, std::size_t column) {
                            matrix.values_[value] = dense[row * columns + column];
                          });
  return matrix;
}

auto BlockSparseMatrix::to_dense() const -> std::vector<double> {
  const std::size_t columns = pattern_->column_count();
  std::vector<double> dense(pattern_->row_count() * columns, 0.0);
  const bool mirrored = is_mirrored(pattern_->symmetry());
  const double mirror_sign = pattern_->symmetry() == MatrixSymmetry::antisymmetric ? -1.0 : 1.0;
  for_each_stored_element(*pattern_, [&](std::size_t value, std::size_t row, std::size_t column) {
    dense[row * columns + column] = values_[value];
    if (mirrored && row != column) {
      dense[column * columns + row] = mirror_sign * values_[value];
    }
  });
  return dense;
}

auto BlockSparseMatrix::to_dense_matrix() const -> DenseMatrix {
  DenseMatrix dense(pattern_->row_count(), pattern_->column_count(), pattern_->symmetry());
  for_each_stored_element(*pattern_, [&](std::size_t value, std::size_t row, std::size_t column) {
    dense.set(row, column, values_[value]);
  });
  return dense;
}

}  // namespace fika

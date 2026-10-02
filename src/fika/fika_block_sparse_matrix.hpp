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

#ifndef fika_block_sparse_matrix_hpp
#define fika_block_sparse_matrix_hpp

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "fika_molecular_basis.hpp"
#include "fika_dense_matrix.hpp"
#include "fika_matrix_symmetry.hpp"

namespace fika {

/// Storage of same-centre blocks (shells of equal angular momentum on one atom, square matrices):
/// the full block, or one value per contracted pair for integrals diagonal in m (value times
/// identity).
enum class DiagonalFormat { full, scalar };

/// How a block is stored: dense, a shell with itself, or two different shells of equal angular
/// momentum on the same atom (square matrices only; rectangular matrices are always dense).
enum class BlockKind { dense, same_shell, same_centre };

class BlockSparseMatrix;

namespace detail {

/// Allocator that default-initializes: resizing a vector of doubles leaves the new elements
/// uninitialized, so large value arrays can be first touched (and zeroed) in parallel.
template <typename T>
struct DefaultInitAllocator : std::allocator<T> {
  template <typename U>
  struct rebind {
    using other = DefaultInitAllocator<U>;
  };

  using std::allocator<T>::allocator;

  template <typename U>
  void construct(U* p) noexcept(std::is_nothrow_default_constructible_v<U>) {
    ::new (static_cast<void*>(p)) U;
  }
  template <typename U, typename... Args>
  void construct(U* p, Args&&... args) {
    ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
  }
};

/// Value storage of block-sparse matrices.
using ValueVector = std::vector<double, DefaultInitAllocator<double>>;

/// Shell lists and storage rules shared by patterns and builders.
struct MatrixShape {
  std::vector<MolecularShell> bra_shells;
  std::vector<MolecularShell> ket_shells;
  bool square = true;  // ket shells are the bra shells (diagonal blocks exist)
  MatrixSymmetry symmetry = MatrixSymmetry::general;
  DiagonalFormat diagonal_format = DiagonalFormat::full;
};

/// Blocks appended in (bra, ket) order: bra and ket shell per block and its values.
struct BlockBuffer {
  std::vector<std::uint32_t> bras;
  std::vector<std::uint32_t> kets;
  std::vector<std::size_t> value_ends;
  std::vector<double> values;
};

/// Joins buffered blocks, already in (bra, ket) shell order across the buffers, into a matrix;
/// blocks without values are dropped and buffers are released as they are consumed.
auto assemble(MatrixShape shape, std::vector<BlockBuffer*> buffers) -> BlockSparseMatrix;

}  // namespace detail

/// Block-sparse structure over pairs of molecular shells.
///
/// A block couples a bra shell and a ket shell (keyed by their first basis-function indices) and
/// stores their whole sub-matrix. A shell has N = n_c (2l + 1) functions, ordered by contracted
/// function k, then m = -l..l. Square matrices use one shell list for bra and ket; rectangular ones
/// (e.g. orbital x fitting basis) have separate lists and are always general and dense.
/// Symmetric and antisymmetric matrices store only blocks with bra shell <= ket shell.
///
/// Values per block, row-major:
///  - dense blocks (different atoms, different l, rectangular matrices, and same-centre blocks
///    with DiagonalFormat::full): N_bra x N_ket;
///  - a shell with itself, DiagonalFormat::full: N x N (general), the upper triangle i <= j
///    (symmetric) or i < j (antisymmetric) over the shell's functions;
///  - same-centre blocks with DiagonalFormat::scalar: one value per contracted pair (k_a, k_b)
///    standing for value x identity in m; for a shell with itself all pairs (general), k_a <= k_b
///    (symmetric) or k_a < k_b (antisymmetric), for two different shells all pairs.
class SparsityPattern {
 public:
  /// Square pattern over `shells` from its rows: the blocks of bra shell s have the ket shells
  /// block_kets[row_offsets[s] .. row_offsets[s + 1]], strictly increasing (and >= s for
  /// symmetric and antisymmetric matrices). Value offsets follow from the storage rules. Throws
  /// std::invalid_argument for rows that break these rules.
  static auto from_rows(std::span<const MolecularShell> shells, MatrixSymmetry symmetry,
                        DiagonalFormat diagonal_format, std::vector<std::size_t> row_offsets,
                        std::vector<std::uint32_t> block_kets)
      -> std::shared_ptr<const SparsityPattern>;

  /// Rectangular (general) pattern between `bra_shells` and `ket_shells` from its rows.
  static auto from_rows(std::span<const MolecularShell> bra_shells,
                        std::span<const MolecularShell> ket_shells,
                        std::vector<std::size_t> row_offsets, std::vector<std::uint32_t> block_kets)
      -> std::shared_ptr<const SparsityPattern>;

  auto symmetry() const noexcept -> MatrixSymmetry { return shape_.symmetry; }
  auto diagonal_format() const noexcept -> DiagonalFormat { return shape_.diagonal_format; }
  auto square() const noexcept -> bool { return shape_.square; }
  auto bra_shells() const noexcept -> std::span<const MolecularShell> { return shape_.bra_shells; }
  auto ket_shells() const noexcept -> std::span<const MolecularShell> { return shape_.ket_shells; }

  /// Dimensions of the dense matrix (bra and ket basis functions).
  auto row_count() const noexcept -> std::size_t;
  auto column_count() const noexcept -> std::size_t;

  auto block_count() const noexcept -> std::size_t { return block_ket_.size(); }
  auto value_count() const noexcept -> std::size_t { return block_value_offsets_.back(); }

  /// Blocks of bra shell s are row_offsets()[s] .. row_offsets()[s + 1].
  auto row_offsets() const noexcept -> std::span<const std::size_t> { return row_offsets_; }

  auto block_ket(std::size_t block) const noexcept -> std::size_t { return block_ket_[block]; }

  /// Ket shells of all blocks, row after row.
  auto block_kets() const noexcept -> std::span<const std::uint32_t> { return block_ket_; }

  /// Offset of a block's values in the value array, and their number.
  auto block_value_offset(std::size_t block) const noexcept -> std::size_t {
    return block_value_offsets_[block];
  }
  auto block_value_count(std::size_t block) const noexcept -> std::size_t {
    return block_value_offsets_[block + 1] - block_value_offsets_[block];
  }

  /// Block between the bra shell starting at basis function bra_start and the ket shell starting
  /// at ket_start, if stored. Throws std::invalid_argument if either is not a shell start.
  auto find_block(std::size_t bra_start, std::size_t ket_start) const -> std::optional<std::size_t>;

 private:
  friend auto detail::assemble(detail::MatrixShape, std::vector<detail::BlockBuffer*>)
      -> BlockSparseMatrix;
  SparsityPattern() = default;

  static auto from_rows(detail::MatrixShape shape, std::vector<std::size_t> row_offsets,
                        std::vector<std::uint32_t> block_kets)
      -> std::shared_ptr<const SparsityPattern>;

  detail::MatrixShape shape_;
  std::vector<std::size_t> row_offsets_;          // bra shells + 1
  std::vector<std::uint32_t> block_ket_;          // blocks
  std::vector<std::size_t> block_value_offsets_;  // blocks + 1
};

/// Storage kind of the block between bra shell `bra` and ket shell `ket` (indices into the shell
/// lists) of a square or rectangular matrix.
auto block_kind(const MolecularShell& bra, const MolecularShell& ket, bool same_shell,
                bool square) noexcept -> BlockKind;

/// Number of values stored for a block of the given kind.
auto block_value_count(const MolecularShell& bra, const MolecularShell& ket, BlockKind kind,
                       MatrixSymmetry symmetry, DiagonalFormat format) noexcept -> std::size_t;

/// Values on a (possibly shared) sparsity pattern.
class BlockSparseMatrix {
 public:
  /// Zero matrix on an existing pattern (zeroed in parallel).
  explicit BlockSparseMatrix(std::shared_ptr<const SparsityPattern> pattern);

  /// Matrix on `pattern` with values taken from a dense row-major matrix; entries outside the
  /// pattern (and the mirrored triangle of symmetric/antisymmetric matrices) are ignored.
  static auto from_dense(std::shared_ptr<const SparsityPattern> pattern,
                         std::span<const double> dense) -> BlockSparseMatrix;

  auto pattern() const noexcept -> const SparsityPattern& { return *pattern_; }
  auto shared_pattern() const noexcept -> const std::shared_ptr<const SparsityPattern>& {
    return pattern_;
  }

  auto values() const noexcept -> std::span<const double> { return values_; }
  auto values() noexcept -> std::span<double> { return values_; }

  /// Values of `block`.
  auto block_values(std::size_t block) const noexcept -> std::span<const double> {
    return std::span(values_).subspan(pattern_->block_value_offset(block),
                                      pattern_->block_value_count(block));
  }
  auto block_values(std::size_t block) noexcept -> std::span<double> {
    return std::span(values_).subspan(pattern_->block_value_offset(block),
                                      pattern_->block_value_count(block));
  }

  /// Dense row-major matrix (row_count x column_count), with mirrored blocks filled in for
  /// symmetric (+) and antisymmetric (-) matrices and scalar diagonal values expanded.
  auto to_dense() const -> std::vector<double>;

  /// Dense matrix of the same symmetry (packed triangle for symmetric and antisymmetric
  /// matrices), scalar diagonal values expanded.
  auto to_dense_matrix() const -> DenseMatrix;

 private:
  friend auto detail::assemble(detail::MatrixShape, std::vector<detail::BlockBuffer*>)
      -> BlockSparseMatrix;
  BlockSparseMatrix(std::shared_ptr<const SparsityPattern> pattern, detail::ValueVector values)
      : pattern_(std::move(pattern)), values_(std::move(values)) {}

  std::shared_ptr<const SparsityPattern> pattern_;
  detail::ValueVector values_;
};

/// Builds pattern and values in one pass with one buffer per bra shell (row), so different rows
/// may be filled concurrently; a row must be filled by one thread at a time, kets increasing.
///
///   std::span<double> v = builder.row(a).add_block(b);
class BlockSparseBuilder {
 public:
  class Row {
   public:
    /// Starts the block with ket shell `ket_shell` (increasing within a row) and returns its
    /// values to fill: zero-initialized, valid until the next add_block on this row.
    auto add_block(std::size_t ket_shell) -> std::span<double>;

   private:
    friend class BlockSparseBuilder;
    Row(const detail::MatrixShape& shape, std::size_t bra) : shape_(&shape), bra_(bra) {}
    const detail::MatrixShape* shape_;
    std::size_t bra_;
    detail::BlockBuffer buffer_;
  };

  /// Square matrix over `shells`.
  BlockSparseBuilder(std::span<const MolecularShell> shells, MatrixSymmetry symmetry,
                     DiagonalFormat diagonal_format);

  /// Rectangular (general) matrix between `bra_shells` and `ket_shells`.
  BlockSparseBuilder(std::span<const MolecularShell> bra_shells,
                     std::span<const MolecularShell> ket_shells);

  // Rows refer to the builder's shape, so it stays in place.
  BlockSparseBuilder(const BlockSparseBuilder&) = delete;
  auto operator=(const BlockSparseBuilder&) -> BlockSparseBuilder& = delete;

  auto row(std::size_t bra_shell) -> Row& {
    assert(bra_shell < rows_.size());
    return rows_[bra_shell];
  }

  auto finish() && -> BlockSparseMatrix;

 private:
  detail::MatrixShape shape_;
  std::vector<Row> rows_;
};

}  // namespace fika

#endif  // fika_block_sparse_matrix_hpp

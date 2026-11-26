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

#ifndef fika_uniform_m2l_hpp
#define fika_uniform_m2l_hpp

// Internal: multipole-to-local translations of a uniform octree as batched matrix products.
//
// In a uniform octree where cells interact once separated by s cells, a cell's interaction list
// holds the cells whose centres differ from its own by h v, h the cell edge and v an integer
// offset with components in -(2s + 1)..2s + 1 and largest magnitude above s: 316 offsets for
// s = 1, 1206 for s = 2. Since I_lm(h v) = I_lm(v) / h^(l+1),
// the M2L operator of an offset at any level is a fixed matrix between diagonal scalings, and
// every translation with the same offset and level is one row of a matrix product. The operators
// act on expansions as real vectors (real parts, then imaginary parts, of the m >= 0
// coefficients) and are stored for the offsets with non-negative components (56 for s = 1, 189
// for s = 2); the others
// follow by reflecting axes, which maps an expansion X_lm to (-1)^m conj(X_lm) (x), conj(X_lm)
// (y) or (-1)^(l+m) X_lm (z).

#include <array>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace fika::detail {

using CellOffset = std::array<int, 3>;

/// Scratch storage of UniformM2L::apply (per thread).
struct UniformM2LWorkspace {
  std::vector<double> packed;
  std::vector<double> product;
  std::vector<double> in_factor;
  std::vector<double> out_factor;
};

class UniformM2L {
 public:
  /// Operators of expansion order `order` (1..max_expansion_order) for separation 1 or 2; they
  /// take 56 or 189 (2P)^2 doubles, P = (order + 1)(order + 2) / 2 (42 or 141 MB at order 16).
  /// Throws std::invalid_argument otherwise.
  explicit UniformM2L(int order, int separation = 1);

  auto order() const noexcept -> int { return order_; }
  auto separation() const noexcept -> int { return separation_; }

  /// The interaction-list offsets (target centre - source centre) / edge in a fixed order.
  auto offsets() const noexcept -> std::span<const CellOffset> { return offsets_; }

  /// The offsets of separation 1 or 2 in the order of offsets().
  static auto interaction_offsets(int separation) -> std::vector<CellOffset>;

  /// Bytes taken by the operators.
  auto operator_bytes() const noexcept -> std::size_t {
    return kernels_.size() * 4 * size_ * size_ * sizeof(double);
  }

  /// Whether `offset` is in the interaction list.
  auto interacts(const CellOffset& offset) const noexcept -> bool;

  /// locals[i] += M2L of multipoles[i] (about the source centre) to the target centre, for n
  /// translations with target - source = edge * offset (an interaction-list offset); expansions
  /// of order order() stored contiguously, one after the other. Equal to multipole_to_local up
  /// to rounding. Throws std::invalid_argument for another offset or wrong sizes.
  void apply(const CellOffset& offset, double edge,
             std::span<const std::complex<double>> multipoles,
             std::span<std::complex<double>> locals, std::size_t n,
             UniformM2LWorkspace& workspace) const;

 private:
  int order_;
  int separation_;
  std::size_t size_;  // P
  std::vector<CellOffset> offsets_;
  std::vector<int> kernel_index_;  // canonical offset (x, y, z in 0..2s + 1) -> kernel, or -1
  std::vector<std::vector<double>> kernels_;  // (2P) x (2P), row-major, rows = outputs
};

}  // namespace fika::detail

#endif  // fika_uniform_m2l_hpp

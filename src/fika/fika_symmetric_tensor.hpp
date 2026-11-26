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

#ifndef fika_symmetric_tensor_hpp
#define fika_symmetric_tensor_hpp

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>

namespace fika {

/// Symmetric Cartesian tensor of rank n (a permanent multipole of an MM site), storing its
/// (n + 1)(n + 2) / 2 unique components in lexicographic order of the sorted axis indices
/// (0 = x, 1 = y, 2 = z):
///   rank 1 (dipole):     x, y, z
///   rank 2 (quadrupole): xx, xy, xz, yy, yz, zz
///   rank 3 (octupole):   xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
/// The components are stored as given: conventions (primitive or traceless moments, prefactors)
/// belong to the drivers that use them.
template <int Rank>
struct SymmetricTensor {
  static_assert(Rank >= 1, "a symmetric tensor has rank >= 1");

  static constexpr std::size_t size = static_cast<std::size_t>((Rank + 1) * (Rank + 2) / 2);

  std::array<double, size> components{};

  /// Position in `components` of the component with the given axes, in any order.
  static constexpr auto index(std::array<int, Rank> axes) noexcept -> std::size_t {
    for (std::size_t i = 1; i < axes.size(); ++i) {  // insertion sort
      for (std::size_t j = i; j > 0 && axes[j - 1] > axes[j]; --j) {
        const int swap = axes[j - 1];
        axes[j - 1] = axes[j];
        axes[j] = swap;
      }
    }
    // Count the non-decreasing tuples before `axes` in lexicographic order.
    std::array<int, Rank> tuple{};
    std::size_t position = 0;
    while (tuple != axes) {
      std::size_t i = tuple.size();
      while (tuple[i - 1] == 2) {
        --i;
      }
      const int value = tuple[i - 1] + 1;
      for (std::size_t j = i - 1; j < tuple.size(); ++j) {
        tuple[j] = value;
      }
      ++position;
    }
    return position;
  }

  /// Component with the given axes (0..2) in any order, e.g. q(2, 0) == q(0, 2).
  template <std::integral... Axes>
    requires(sizeof...(Axes) == Rank)
  constexpr auto operator()(Axes... axes) const noexcept -> double {
    assert(((axes >= 0 && axes <= 2) && ...));
    return components[index({static_cast<int>(axes)...})];
  }

  template <std::integral... Axes>
    requires(sizeof...(Axes) == Rank)
  constexpr auto operator()(Axes... axes) noexcept -> double& {
    assert(((axes >= 0 && axes <= 2) && ...));
    return components[index({static_cast<int>(axes)...})];
  }

  friend constexpr auto operator==(const SymmetricTensor&, const SymmetricTensor&)
      -> bool = default;
};

using Dipole = SymmetricTensor<1>;
using Quadrupole = SymmetricTensor<2>;
using Octupole = SymmetricTensor<3>;

}  // namespace fika

#endif  // fika_symmetric_tensor_hpp

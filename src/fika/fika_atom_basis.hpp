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

#ifndef fika_atom_basis_hpp
#define fika_atom_basis_hpp

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_element.hpp"

namespace fika {

/// Basis set of one element: shells ordered by angular momentum (all s, then all p, ...), shells
/// of equal angular momentum in their input order; gaps in angular momentum are allowed.
class AtomBasis {
 public:
  /// Shells may come in any order and several per angular momentum; they are sorted by angular
  /// momentum, keeping the input order within each. The label is stored lower-cased. Throws
  /// std::invalid_argument for an empty label or no shells.
  AtomBasis(Element element, std::string_view label, std::vector<BasisShell> shells);

  auto element() const noexcept -> Element { return element_; }

  /// Standard basis-set label, lower-cased, e.g. "cc-pvdz".
  auto label() const noexcept -> std::string_view { return label_; }

  /// Shells in increasing angular momentum.
  auto shells() const noexcept -> std::span<const BasisShell> { return shells_; }

  /// Highest angular momentum in the basis.
  auto max_angular_momentum() const noexcept -> int { return angular_momentum(shells_.back()); }

  /// Consecutive shells with angular momentum l (empty if the basis has none).
  auto shells_with_angular_momentum(int l) const noexcept -> std::span<const BasisShell>;

  /// Index into shells() of the first shell with angular momentum l (shells().size() if none).
  auto first_shell_with_angular_momentum(int l) const noexcept -> std::size_t;

  /// Number of basis functions: sum over shells of (2l + 1) times contracted functions.
  auto function_count() const noexcept -> std::size_t { return shell_offsets_.back(); }

  /// Offset of each shell's first function within the atom (shells().size() + 1 entries); a
  /// shell's functions are ordered by contracted function, then m = -l..l.
  auto shell_offsets() const noexcept -> std::span<const std::size_t> { return shell_offsets_; }

 private:
  Element element_;
  std::string label_;
  std::vector<BasisShell> shells_;
  std::vector<std::size_t> shell_offsets_;
};

}  // namespace fika

#endif  // fika_atom_basis_hpp

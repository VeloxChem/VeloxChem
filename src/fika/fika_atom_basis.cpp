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

#include "fika_atom_basis.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace fika {

namespace {

constexpr auto to_lower(char c) noexcept -> char {
  return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
}

}  // namespace

AtomBasis::AtomBasis(Element element, std::string_view label, std::vector<BasisShell> shells)
    : element_(element) {
  if (label.empty()) {
    throw std::invalid_argument("fika::AtomBasis: empty basis-set label");
  }
  if (shells.empty()) {
    throw std::invalid_argument("fika::AtomBasis: no shells for " + std::string(element.label()));
  }
  std::ranges::transform(label, std::back_inserter(label_), to_lower);

  // Stable, so shells of equal l keep their input order.
  std::ranges::stable_sort(shells, {},
                           [](const BasisShell& shell) { return angular_momentum(shell); });
  shells_ = std::move(shells);

  shell_offsets_.reserve(shells_.size() + 1);
  shell_offsets_.push_back(0);
  for (const BasisShell& shell : shells_) {
    shell_offsets_.push_back(shell_offsets_.back() + fika::function_count(shell));
  }
}

auto AtomBasis::first_shell_with_angular_momentum(int l) const noexcept -> std::size_t {
  const auto shell = std::ranges::lower_bound(
      shells_, l, {}, [](const BasisShell& s) { return angular_momentum(s); });
  return static_cast<std::size_t>(shell - shells_.begin());
}

auto AtomBasis::shells_with_angular_momentum(int l) const noexcept -> std::span<const BasisShell> {
  const auto [first, last] = std::ranges::equal_range(
      shells_, l, {}, [](const BasisShell& s) { return angular_momentum(s); });
  return {first, last};
}

}  // namespace fika

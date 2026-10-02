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

#include "fika_atom_pair_group.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fika {

namespace {

void check_block_size(std::size_t block_size) {
  if (block_size == 0) {
    throw std::invalid_argument("fika::AtomPairGroupFactory: block size must be positive");
  }
}

/// First position of row p in the upper triangle (including the diagonal) of an n x n matrix.
auto triangle_row_start(std::size_t n, std::size_t p) -> std::size_t {
  return p * n - (p == 0 ? 0 : p * (p - 1) / 2);
}

/// Row of position t in the row-major upper triangle (including the diagonal) of an n x n matrix.
auto triangle_row(std::size_t n, std::size_t t) -> std::size_t {
  // Largest p with p * n - p * (p - 1) / 2 <= t, from the quadratic; then correct rounding.
  const double b = 2.0 * static_cast<double>(n) + 1.0;
  const double root = std::sqrt(b * b - 8.0 * static_cast<double>(t));
  auto p = static_cast<std::size_t>(std::max(0.0, std::floor((b - root) / 2.0)));
  p = std::min(p, n - 1);
  while (p > 0 && triangle_row_start(n, p) > t) {
    --p;
  }
  while (p + 1 < n && triangle_row_start(n, p + 1) <= t) {
    ++p;
  }
  return p;
}

}  // namespace

auto automatic_block_size(std::size_t pair_count, std::size_t threads, PairCost cost)
    -> std::size_t {
  const bool light = cost == PairCost::light;
  const std::size_t blocks_per_thread = light ? 16 : 64;
  const std::size_t smallest = light ? 256 : 32;
  const std::size_t largest = light ? 1024 : 256;
  return std::clamp(pair_count / (blocks_per_thread * std::max<std::size_t>(threads, 1)), smallest,
                    largest);
}

AtomPairGroupFactory::AtomPairGroupFactory(const MolecularBasis& basis, std::size_t block_size)
    : symmetric_(true), block_size_(block_size) {
  check_block_size(block_size);
  for (const AtomBasisPair pair : basis.unique_basis_pairs()) {
    add_group(pair, basis.atoms_with_basis(pair.bra), basis.atoms_with_basis(pair.ket),
              pair.bra == pair.ket);
  }
}

AtomPairGroupFactory::AtomPairGroupFactory(const MolecularBasis& bra, const MolecularBasis& ket,
                                           std::size_t block_size)
    : symmetric_(false), block_size_(block_size) {
  check_block_size(block_size);
  if (bra.atom_count() != ket.atom_count()) {
    throw std::invalid_argument("fika::AtomPairGroupFactory: bra basis has " +
                                std::to_string(bra.atom_count()) + " atoms, ket basis " +
                                std::to_string(ket.atom_count()));
  }
  for (const AtomBasisPair pair : bra.unique_basis_pairs(ket)) {
    add_group(pair, bra.atoms_with_basis(pair.bra), ket.atoms_with_basis(pair.ket), false);
  }
}

void AtomPairGroupFactory::add_group(AtomBasisPair basis_pair,
                                     std::span<const std::size_t> bra_atoms,
                                     std::span<const std::size_t> ket_atoms, bool triangular) {
  const std::size_t n = bra_atoms.size();
  const std::size_t size = triangular ? n * (n + 1) / 2 : n * ket_atoms.size();
  groups_.push_back({basis_pair, bra_atoms, ket_atoms, triangular, size});
  group_block_offsets_.push_back(group_block_offsets_.back() +
                                 (size + block_size_ - 1) / block_size_);
}

auto AtomPairGroupFactory::next(AtomPairGroup& group) -> bool {
  const std::size_t index = next_block_.fetch_add(1, std::memory_order_relaxed);
  if (index >= block_count()) {
    return false;
  }
  block(index, group);
  return true;
}

void AtomPairGroupFactory::block(std::size_t index, AtomPairGroup& group) const {
  if (index >= block_count()) {
    throw std::out_of_range("fika::AtomPairGroupFactory: block " + std::to_string(index) + " of " +
                            std::to_string(block_count()));
  }
  const auto offset = std::ranges::upper_bound(group_block_offsets_, index) - 1;
  const Group& source = groups_[static_cast<std::size_t>(offset - group_block_offsets_.begin())];
  const std::size_t begin = (index - *offset) * block_size_;
  const std::size_t end = std::min(begin + block_size_, source.size);

  group.index = index;
  group.basis_pair = source.basis_pair;
  group.symmetric = symmetric_;
  group.off_diagonal_pairs.clear();
  group.diagonal_atoms.clear();

  const auto add = [&group](std::size_t bra, std::size_t ket) {
    if (bra == ket) {
      group.diagonal_atoms.push_back(bra);
    } else {
      group.off_diagonal_pairs.push_back({bra, ket});
    }
  };

  if (source.triangular) {
    const std::size_t n = source.bra_atoms.size();
    std::size_t p = triangle_row(n, begin);
    std::size_t q = p + (begin - triangle_row_start(n, p));
    for (std::size_t t = begin; t < end; ++t) {
      add(source.bra_atoms[p], source.bra_atoms[q]);
      if (++q == n) {
        ++p;
        q = p;
      }
    }
  } else {
    const std::size_t n_ket = source.ket_atoms.size();
    for (std::size_t t = begin; t < end; ++t) {
      add(source.bra_atoms[t / n_ket], source.ket_atoms[t % n_ket]);
    }
  }
}

}  // namespace fika

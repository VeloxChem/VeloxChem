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

#include "fika_molecular_basis.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::MolecularBasis: " + reason);
}

}  // namespace

MolecularBasis::MolecularBasis(std::vector<AtomBasis> unique_bases,
                               std::vector<std::uint16_t> atom_basis_indices)
    : unique_bases_(std::move(unique_bases)), atom_basis_indices_(std::move(atom_basis_indices)) {
  if (unique_bases_.size() > std::size_t{std::numeric_limits<std::uint16_t>::max()} + 1) {
    fail("too many unique bases");
  }
  for (std::size_t i = 0; i < unique_bases_.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (unique_bases_[i].element() == unique_bases_[j].element() &&
          unique_bases_[i].label() == unique_bases_[j].label()) {
        fail("duplicate basis " + std::string(unique_bases_[i].label()) + " for " +
             std::string(unique_bases_[i].element().label()));
      }
    }
  }
  for (std::size_t atom = 0; atom < atom_basis_indices_.size(); ++atom) {
    if (atom_basis_indices_[atom] >= unique_bases_.size()) {
      fail("atom " + std::to_string(atom) + " refers to basis " +
           std::to_string(atom_basis_indices_[atom]) + " of " +
           std::to_string(unique_bases_.size()));
    }
  }

  atom_offsets_.reserve(atom_basis_indices_.size() + 1);
  atom_offsets_.push_back(0);
  for (const std::uint16_t index : atom_basis_indices_) {
    atom_offsets_.push_back(atom_offsets_.back() + unique_bases_[index].function_count());
  }

  index_atoms_by_basis();

  for (std::size_t atom = 0; atom < atom_basis_indices_.size(); ++atom) {
    const AtomBasis& basis = atom_basis(atom);
    for (std::size_t shell = 0; shell < basis.shells().size(); ++shell) {
      const BasisShell& basis_shell = basis.shells()[shell];
      shells_.push_back({atom, shell, angular_momentum(basis_shell), contraction_count(basis_shell),
                         atom_offsets_[atom] + basis.shell_offsets()[shell], kind(basis_shell)});
    }
  }
}

void MolecularBasis::index_atoms_by_basis() {
  // Counting sort of atoms by basis index keeps atoms in increasing order within each basis.
  basis_atom_offsets_.assign(unique_bases_.size() + 1, 0);
  for (const std::uint16_t index : atom_basis_indices_) {
    ++basis_atom_offsets_[index + 1];
  }
  for (std::size_t basis = 0; basis < unique_bases_.size(); ++basis) {
    basis_atom_offsets_[basis + 1] += basis_atom_offsets_[basis];
  }
  basis_atoms_.resize(atom_basis_indices_.size());
  std::vector<std::size_t> next(basis_atom_offsets_.begin(), basis_atom_offsets_.end() - 1);
  for (std::size_t atom = 0; atom < atom_basis_indices_.size(); ++atom) {
    basis_atoms_[next[atom_basis_indices_[atom]]++] = atom;
  }

  for (std::size_t bra = 0; bra < unique_bases_.size(); ++bra) {
    for (std::size_t ket = bra; ket < unique_bases_.size(); ++ket) {
      if (!atoms_with_basis(bra).empty() && !atoms_with_basis(ket).empty()) {
        unique_basis_pairs_.push_back(
            {static_cast<std::uint16_t>(bra), static_cast<std::uint16_t>(ket)});
      }
    }
  }
}

auto MolecularBasis::unique_basis_pairs(const MolecularBasis& ket_basis) const
    -> std::vector<AtomBasisPair> {
  std::vector<AtomBasisPair> pairs;
  for (std::size_t bra = 0; bra < unique_bases_.size(); ++bra) {
    if (atoms_with_basis(bra).empty()) {
      continue;
    }
    for (std::size_t ket = 0; ket < ket_basis.unique_bases().size(); ++ket) {
      if (!ket_basis.atoms_with_basis(ket).empty()) {
        pairs.push_back({static_cast<std::uint16_t>(bra), static_cast<std::uint16_t>(ket)});
      }
    }
  }
  return pairs;
}

auto MolecularBasis::shell_at(std::size_t function_offset) const -> std::size_t {
  const auto shell =
      std::ranges::lower_bound(shells_, function_offset, {}, &MolecularShell::function_offset);
  if (shell == shells_.end() || shell->function_offset != function_offset) {
    fail("no shell starts at basis function " + std::to_string(function_offset));
  }
  return static_cast<std::size_t>(shell - shells_.begin());
}

auto MolecularBasis::function_location(std::size_t index) const -> BasisFunctionLocation {
  if (index >= function_count()) {
    throw std::out_of_range("fika::MolecularBasis: basis function " + std::to_string(index) +
                            " of " + std::to_string(function_count()));
  }
  // Last atom whose offset is <= index (atoms without functions are skipped).
  const auto atom_end = std::ranges::upper_bound(atom_offsets_, index);
  const auto atom = static_cast<std::size_t>(atom_end - atom_offsets_.begin()) - 1;
  const AtomBasis& basis = atom_basis(atom);

  const std::size_t local = index - atom_offsets_[atom];
  const auto offsets = basis.shell_offsets();
  const auto shell_end = std::ranges::upper_bound(offsets, local);
  const auto shell = static_cast<std::size_t>(shell_end - offsets.begin()) - 1;

  const int l = angular_momentum(basis.shells()[shell]);
  const auto components = static_cast<std::size_t>(2 * l + 1);
  const std::size_t within = local - offsets[shell];
  return {atom, shell, within / components, static_cast<int>(within % components) - l};
}

auto MolecularBasis::for_elements(std::span<const Element> elements,
                                  std::span<const AtomBasis> bases) -> MolecularBasis {
  constexpr std::uint16_t unassigned = std::numeric_limits<std::uint16_t>::max();
  std::array<std::uint16_t, 119> index_by_atomic_number;
  index_by_atomic_number.fill(unassigned);

  std::vector<AtomBasis> unique_bases;
  std::vector<std::uint16_t> indices;
  indices.reserve(elements.size());
  for (const Element element : elements) {
    auto& index = index_by_atomic_number[static_cast<std::size_t>(element.atomic_number())];
    if (index == unassigned) {
      const auto matches = std::ranges::count(bases, element, &AtomBasis::element);
      if (matches != 1) {
        fail(std::string(matches == 0 ? "no basis" : "several bases") + " for " +
             std::string(element.label()));
      }
      index = static_cast<std::uint16_t>(unique_bases.size());
      unique_bases.push_back(*std::ranges::find(bases, element, &AtomBasis::element));
    }
    indices.push_back(index);
  }
  return MolecularBasis(std::move(unique_bases), std::move(indices));
}

}  // namespace fika

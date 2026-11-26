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

#ifndef fika_molecular_basis_hpp
#define fika_molecular_basis_hpp

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "fika_atom_basis.hpp"
#include "fika_element.hpp"

namespace fika {

/// Position of a basis function in the molecular basis ordering.
struct BasisFunctionLocation {
  std::size_t atom;
  std::size_t shell;        // index into atom_basis(atom).shells()
  std::size_t contraction;  // contracted function within the shell
  int m;                    // -l..l

  friend auto operator==(const BasisFunctionLocation&, const BasisFunctionLocation&)
      -> bool = default;
};

/// One shell of one atom in the molecular basis.
struct MolecularShell {
  std::size_t atom;
  std::size_t shell;  // index into atom_basis(atom).shells()
  int angular_momentum;
  std::size_t contraction_count;
  std::size_t function_offset;  // index of the shell's first basis function
  ShellKind kind;

  /// Number of basis functions: (2l + 1) times contracted functions.
  auto function_count() const noexcept -> std::size_t {
    return static_cast<std::size_t>(2 * angular_momentum + 1) * contraction_count;
  }
};

/// Pair of unique atom-basis indices (into MolecularBasis::unique_bases()).
struct AtomBasisPair {
  std::uint16_t bra;
  std::uint16_t ket;

  friend auto operator==(const AtomBasisPair&, const AtomBasisPair&) -> bool = default;
};

/// Basis of a molecule: unique atom bases and, for every atom, the index of its basis in
/// unique_bases(). An atom may use a basis of a different element.
///
/// Basis functions are ordered by atom, then shell (increasing l), then contracted function,
/// then component m = -l..l.
class MolecularBasis {
 public:
  /// Throws std::invalid_argument if an index is out of range or two unique bases share
  /// the same (element, label).
  MolecularBasis(std::vector<AtomBasis> unique_bases,
                 std::vector<std::uint16_t> atom_basis_indices);

  /// Assigns every atom the basis of its element from `bases`; only bases actually used are kept,
  /// in order of first use. Throws std::invalid_argument if an element has no basis or several.
  static auto for_elements(std::span<const Element> elements, std::span<const AtomBasis> bases)
      -> MolecularBasis;

  auto atom_count() const noexcept -> std::size_t { return atom_basis_indices_.size(); }

  auto unique_bases() const noexcept -> std::span<const AtomBasis> { return unique_bases_; }

  auto atom_basis_indices() const noexcept -> std::span<const std::uint16_t> {
    return atom_basis_indices_;
  }

  /// Total number of basis functions.
  auto function_count() const noexcept -> std::size_t { return atom_offsets_.back(); }

  /// Index of each atom's first basis function (atom_count() + 1 entries).
  auto atom_offsets() const noexcept -> std::span<const std::size_t> { return atom_offsets_; }

  /// Index of component m (-l..l) of contracted function k in shell `shell` of atom `atom`.
  auto function_index(std::size_t atom, std::size_t shell, std::size_t k, int m) const noexcept
      -> std::size_t {
    const AtomBasis& basis = atom_basis(atom);
    assert(shell < basis.shells().size());
    const int l = fika::angular_momentum(basis.shells()[shell]);
    assert(k < fika::contraction_count(basis.shells()[shell]));
    assert(m >= -l && m <= l);
    return atom_offsets_[atom] + basis.shell_offsets()[shell] +
           k * static_cast<std::size_t>(2 * l + 1) + static_cast<std::size_t>(m + l);
  }

  /// Pairs of unique bases that are each used by at least one atom, bra <= ket, ordered by
  /// (bra, ket).
  auto unique_basis_pairs() const noexcept -> std::span<const AtomBasisPair> {
    return unique_basis_pairs_;
  }

  /// Pairs (bra, ket) of a unique basis of this molecular basis (bra) and one of `ket_basis`
  /// (ket), both used by at least one atom: the direct product, ordered by (bra, ket).
  auto unique_basis_pairs(const MolecularBasis& ket_basis) const -> std::vector<AtomBasisPair>;

  /// Atoms using unique basis `basis`, in increasing order (empty for an unused basis).
  auto atoms_with_basis(std::size_t basis) const noexcept -> std::span<const std::size_t> {
    assert(basis < unique_bases_.size());
    return std::span(basis_atoms_)
        .subspan(basis_atom_offsets_[basis],
                 basis_atom_offsets_[basis + 1] - basis_atom_offsets_[basis]);
  }

  /// All shells of all atoms in basis-function order.
  auto shells() const noexcept -> std::span<const MolecularShell> { return shells_; }

  /// Index into shells() of the shell whose first basis function is `function_offset`; throws
  /// std::invalid_argument if no shell starts there.
  auto shell_at(std::size_t function_offset) const -> std::size_t;

  /// Inverse of function_index; throws std::out_of_range for index >= function_count().
  auto function_location(std::size_t index) const -> BasisFunctionLocation;

  /// Basis of atom `atom`, i.e. unique_bases()[atom_basis_indices()[atom]].
  auto atom_basis(std::size_t atom) const noexcept -> const AtomBasis& {
    assert(atom < atom_count());
    return unique_bases_[atom_basis_indices_[atom]];
  }

 private:
  void index_atoms_by_basis();

  std::vector<AtomBasis> unique_bases_;
  std::vector<std::uint16_t> atom_basis_indices_;
  std::vector<std::size_t> atom_offsets_;
  std::vector<MolecularShell> shells_;
  std::vector<std::size_t> basis_atom_offsets_;  // unique bases + 1
  std::vector<std::size_t> basis_atoms_;         // atoms grouped by unique basis
  std::vector<AtomBasisPair> unique_basis_pairs_;
};

}  // namespace fika

#endif  // fika_molecular_basis_hpp

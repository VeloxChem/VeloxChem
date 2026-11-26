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

#ifndef fika_atom_pair_group_hpp
#define fika_atom_pair_group_hpp

#include <atomic>
#include <cstddef>
#include <span>
#include <vector>

#include "fika_molecular_basis.hpp"

namespace fika {

/// Two atoms by index: bra side and ket side.
struct AtomPair {
  std::size_t bra;
  std::size_t ket;

  friend auto operator==(const AtomPair&, const AtomPair&) -> bool = default;
};

/// One block of atom pairs sharing a pair of unique atom bases, filled by AtomPairGroupFactory.
struct AtomPairGroup {
  std::size_t index = 0;  // block index, independent of thread count and request order
  AtomBasisPair basis_pair{};
  bool symmetric = false;                    // bra and ket from the same molecular basis
  std::vector<AtomPair> off_diagonal_pairs;  // different atoms, in pair-space order
  std::vector<std::size_t> diagonal_atoms;   // atoms paired with themselves (R_AB = 0)
};

/// Cost of one atom pair: light when bookkeeping dominates (overlap, kinetic), heavy when the
/// kernels do (point-source potentials).
enum class PairCost { light, heavy };

/// Atom pairs per block for `pair_count` pairs on `threads` threads: light pairs aim at 16 blocks
/// per thread within [256, 1024], heavy pairs at 64 within [32, 256] (from a scan over molecules
/// of 24 to 1231 atoms; larger blocks amortize bookkeeping, more blocks balance the load).
auto automatic_block_size(std::size_t pair_count, std::size_t threads, PairCost cost)
    -> std::size_t;

/// Lazily hands out blocks of at most `block_size` atom pairs; next() may be called from any
/// thread. A block never mixes pairs of unique bases.
///
/// Each unique-basis pair (i, j) spans a pair space: for symmetric factories with i == j the
/// upper triangle A <= B over the atoms with basis i, otherwise all (A, B) with A using basis i
/// on the bra side and B basis j on the ket side. Blocks are consecutive ranges of these spaces,
/// numbered group after group; positions with A == B become diagonal atoms.
///
/// The molecular bases must outlive the factory.
class AtomPairGroupFactory {
 public:
  /// Symmetric: bra and ket are `basis`; groups follow basis.unique_basis_pairs().
  AtomPairGroupFactory(const MolecularBasis& basis, std::size_t block_size);

  /// Non-symmetric: bra and ket describe the same molecule (equal atom counts); groups follow
  /// bra.unique_basis_pairs(ket).
  AtomPairGroupFactory(const MolecularBasis& bra, const MolecularBasis& ket,
                       std::size_t block_size);

  AtomPairGroupFactory(const AtomPairGroupFactory&) = delete;
  auto operator=(const AtomPairGroupFactory&) -> AtomPairGroupFactory& = delete;

  auto symmetric() const noexcept -> bool { return symmetric_; }
  auto block_size() const noexcept -> std::size_t { return block_size_; }
  auto block_count() const noexcept -> std::size_t { return group_block_offsets_.back(); }

  /// Fills `group` with the next unclaimed block, reusing its storage; false when none remain.
  auto next(AtomPairGroup& group) -> bool;

  /// Fills `group` with block `index` (< block_count()), independently of next().
  void block(std::size_t index, AtomPairGroup& group) const;

  /// Makes all blocks available to next() again; not while other threads call next().
  void reset() noexcept { next_block_.store(0, std::memory_order_relaxed); }

 private:
  struct Group {
    AtomBasisPair basis_pair;
    std::span<const std::size_t> bra_atoms;
    std::span<const std::size_t> ket_atoms;
    bool triangular;
    std::size_t size;  // positions in the pair space
  };

  void add_group(AtomBasisPair basis_pair, std::span<const std::size_t> bra_atoms,
                 std::span<const std::size_t> ket_atoms, bool triangular);

  bool symmetric_;
  std::size_t block_size_;
  std::vector<Group> groups_;
  std::vector<std::size_t> group_block_offsets_{0};  // groups + 1
  std::atomic<std::size_t> next_block_{0};
};

}  // namespace fika

#endif  // fika_atom_pair_group_hpp

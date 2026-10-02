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

#include "fika_two_centre_driver.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <omp.h>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "fika_atom_pair_group.hpp"

namespace fika::detail {

namespace {

/// Index of each atom's first shell in basis.shells() (atoms + 1 entries).
auto atom_shell_offsets(const MolecularBasis& basis) -> std::vector<std::size_t> {
  std::vector<std::size_t> offsets(basis.atom_count() + 1, 0);
  for (const MolecularShell& shell : basis.shells()) {
    ++offsets[shell.atom + 1];
  }
  for (std::size_t atom = 0; atom < basis.atom_count(); ++atom) {
    offsets[atom + 1] += offsets[atom];
  }
  return offsets;
}

/// Same-centre block between two shells of equal angular momentum on one atom.
struct SameCentreBlock {
  std::size_t bra_shell;  // index into the bra atom basis shells
  std::size_t ket_shell;  // index into the ket atom basis shells
  std::size_t bra_contractions;
  std::size_t ket_contractions;
  std::vector<double> values;  // v(k_a, k_b), row-major
};

/// All same-centre blocks between two atom bases (shell pairs of equal l); with `upper`, bra and
/// ket are the same basis and only pairs with bra shell <= ket shell are kept.
auto same_centre_blocks(const TwoCentreOperator& op, const AtomBasis& bra, const AtomBasis& ket,
                        bool upper) -> std::vector<SameCentreBlock> {
  std::vector<SameCentreBlock> blocks;
  for (std::size_t s_a = 0; s_a < bra.shells().size(); ++s_a) {
    const BasisShell& a = bra.shells()[s_a];
    const int l = angular_momentum(a);
    const std::size_t first = ket.first_shell_with_angular_momentum(l);
    const std::size_t count = ket.shells_with_angular_momentum(l).size();
    for (std::size_t s_b = std::max(first, upper ? s_a : first); s_b < first + count; ++s_b) {
      const BasisShell& b = ket.shells()[s_b];
      blocks.push_back(
          {s_a, s_b, contraction_count(a), contraction_count(b), op.same_centre(a, b)});
    }
  }
  return blocks;
}

/// Same-centre blocks of the current group's basis pair, recomputed when the pair changes.
class SameCentreCache {
 public:
  auto blocks(const TwoCentreOperator& op, const AtomBasis& bra, const AtomBasis& ket,
              AtomBasisPair pair, bool upper) -> const std::vector<SameCentreBlock>& {
    if (!valid_ || !(pair == pair_)) {
      blocks_ = same_centre_blocks(op, bra, ket, upper);
      pair_ = pair;
      valid_ = true;
    }
    return blocks_;
  }

 private:
  bool valid_ = false;
  AtomBasisPair pair_{};
  std::vector<SameCentreBlock> blocks_;
};

/// Shell pair of a pair of unique atom bases: screening cutoff and kernel.
struct ShellPairCutoff {
  double cutoff;  // R^2 beyond which the shell pair is insignificant
  int order;      // l_a + l_b: highest solid-harmonic order the shell pair needs
  KernelShape shape;
  std::unique_ptr<ShellPairKernel> kernel;
};

/// All shell pairs (a, b) of a pair of unique atom bases, bra shell-major: shell pair (a, b) is
/// shell_pairs[a * ket_shell_count + b].
struct BasisPairShells {
  std::size_t ket_shell_count = 0;
  std::vector<ShellPairCutoff> shell_pairs;
  double max_cutoff = -1.0;  // largest cutoff of the shell pairs
};

/// Shell-pair data of every pair of unique bra and ket bases, index bra * (unique ket bases) +
/// ket. Pattern and kernels read the same cutoffs, so they agree on every screened block.
auto basis_pair_shells(const TwoCentreOperator& op, const MolecularBasis& bra,
                       const MolecularBasis& ket, double threshold)
    -> std::vector<BasisPairShells> {
  const std::size_t ket_count = ket.unique_bases().size();
  std::vector<BasisPairShells> tables(bra.unique_bases().size() * ket_count);
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t index = 0; index < tables.size(); ++index) {
    const AtomBasis& bra_basis = bra.unique_bases()[index / ket_count];
    const AtomBasis& ket_basis = ket.unique_bases()[index % ket_count];
    BasisPairShells& table = tables[index];
    table.ket_shell_count = ket_basis.shells().size();
    for (const BasisShell& bra_shell : bra_basis.shells()) {
      for (const BasisShell& ket_shell : ket_basis.shells()) {
        const double cutoff = op.bound(bra_shell, ket_shell).cutoff_distance_squared(threshold);
        const int order = angular_momentum(bra_shell) + angular_momentum(ket_shell);
        auto kernel = op.kernel(bra_shell, ket_shell);
        const KernelShape shape = kernel->shape();
        table.shell_pairs.push_back({cutoff, order, shape, std::move(kernel)});
        table.max_cutoff = std::max(table.max_cutoff, cutoff);
      }
    }
  }
  return tables;
}

/// Atom pairs of one group, kept from the first pass (pattern) to the second (values).
struct GroupPairs {
  AtomBasisPair basis_pair{};
  SortedAtomPairs sorted;                   // pairs within the largest cutoff, by distance
  std::vector<std::size_t> diagonal_atoms;  // atoms paired with themselves
  std::vector<std::size_t> entries;         // neighbour-list entry of each sorted pair
};

/// Storage orientation of an atom pair: in symmetric matrices a pair whose bra atom follows its
/// ket atom is stored transposed, in the rows of the ket atom's shells.
auto is_transposed(const AtomPair& pair, bool symmetric) -> bool {
  return symmetric && pair.bra > pair.ket;
}

/// Atom pair seen from the atom whose shells are its storage rows.
struct Neighbour {
  std::uint32_t atom;            // atom of the stored blocks' ket shells
  std::uint32_t group;           // index into the groups
  std::uint32_t pair;            // index into the group's sorted pairs, or same_atom
  bool transposed;               // the row atom is the pair's ket atom
  double distance_squared;       // R_AB^2 of the pair (copied for the pattern loop)
  const BasisPairShells* table;  // shell-pair table of the group's basis pair
};

constexpr std::uint32_t same_atom = std::numeric_limits<std::uint32_t>::max();

/// For each atom, the atom pairs stored in its rows (itself included), by increasing ket atom,
/// and, filled with the pattern, where each entry's blocks start in each row of the atom.
struct NeighbourLists {
  std::vector<std::size_t> offsets;  // atoms + 1
  std::vector<Neighbour> entries;
  /// Row start of atom A's entry e (offsets[A] <= e < offsets[A + 1]) in local row shell r:
  /// row_starts[start_offsets[A] + (e - offsets[A]) * (shells of A) + r], the position of the
  /// entry's first block within the row.
  std::vector<std::size_t> start_offsets;  // atoms + 1
  std::vector<std::uint32_t> row_starts;

  auto row_start_index(std::size_t atom, std::size_t entry, std::size_t shells, std::size_t r) const
      -> std::size_t {
    return start_offsets[atom] + (entry - offsets[atom]) * shells + r;
  }
};

/// Neighbour lists of the groups' atom pairs; records each sorted pair's entry in its group.
auto neighbour_lists(std::span<GroupPairs> groups, std::size_t atom_count, bool symmetric,
                     const std::function<const BasisPairShells*(AtomBasisPair)>& table_of)
    -> NeighbourLists {
  // Count per row atom, place with atomic cursors, then sort each list: every (row atom, ket
  // atom) occurs once, so the result does not depend on the placement order.
  std::vector<std::atomic<std::size_t>> cursors(atom_count + 1);
  const auto for_each_entry = [&](auto&& f) {
#pragma omp parallel for schedule(dynamic, 16)
    for (std::size_t g = 0; g < groups.size(); ++g) {
      const GroupPairs& group = groups[g];
      const BasisPairShells* table = table_of(group.basis_pair);
      for (std::size_t i = 0; i < group.sorted.pairs.size(); ++i) {
        const AtomPair& pair = group.sorted.pairs[i];
        const bool transposed = is_transposed(pair, symmetric);
        f(transposed ? pair.ket : pair.bra,
          Neighbour{static_cast<std::uint32_t>(transposed ? pair.bra : pair.ket),
                    static_cast<std::uint32_t>(g), static_cast<std::uint32_t>(i), transposed,
                    group.sorted.distances_squared[i], table});
      }
      for (const std::size_t atom : group.diagonal_atoms) {
        f(atom, Neighbour{static_cast<std::uint32_t>(atom), static_cast<std::uint32_t>(g),
                          same_atom, false, 0.0, table});
      }
    }
  };
  for_each_entry([&](std::size_t row_atom, const Neighbour&) {
    cursors[row_atom + 1].fetch_add(1, std::memory_order_relaxed);
  });
  NeighbourLists lists;
  lists.offsets.assign(atom_count + 1, 0);
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    lists.offsets[atom + 1] = lists.offsets[atom] + cursors[atom + 1].load();
    cursors[atom].store(lists.offsets[atom]);
  }
  lists.entries.resize(lists.offsets[atom_count]);
  for (GroupPairs& group : groups) {
    group.entries.resize(group.sorted.pairs.size());
  }
  for_each_entry([&](std::size_t row_atom, const Neighbour& neighbour) {
    lists.entries[cursors[row_atom].fetch_add(1, std::memory_order_relaxed)] = neighbour;
  });
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    std::sort(lists.entries.begin() + static_cast<std::ptrdiff_t>(lists.offsets[atom]),
              lists.entries.begin() + static_cast<std::ptrdiff_t>(lists.offsets[atom + 1]),
              [](const Neighbour& x, const Neighbour& y) { return x.atom < y.atom; });
    for (std::size_t e = lists.offsets[atom]; e < lists.offsets[atom + 1]; ++e) {
      const Neighbour& neighbour = lists.entries[e];
      if (neighbour.pair != same_atom) {
        groups[neighbour.group].entries[neighbour.pair] = e;  // each pair has one entry
      }
    }
  }
  return lists;
}

/// Operator, shells, shell offsets per atom, groups and shell-pair data of one computation.
struct Context {
  const TwoCentreOperator& op;
  const MolecularBasis& bra;
  const MolecularBasis& ket;
  bool symmetric;
  std::vector<std::size_t> bra_first_shell;  // atoms + 1
  std::vector<std::size_t> ket_first_shell;  // atoms + 1
  std::vector<BasisPairShells> tables;
  std::vector<GroupPairs> groups;

  auto table(AtomBasisPair pair) const -> const BasisPairShells& {
    return tables[pair.bra * ket.unique_bases().size() + pair.ket];
  }
};

/// Sparsity pattern of the matrix, built row by row in (bra, ket) order. Bra shell r (local
/// index) of a row atom stores ket shell c of a neighbour if both are on one atom with equal l
/// (and c >= r when symmetric), or if R_AB^2 <= cutoff of the shell pair: the rule of
/// significant_pair_count, which the kernels use.
/// Also records the neighbours' row starts.
auto matrix_pattern(const Context& context, NeighbourLists& neighbours)
    -> std::shared_ptr<const SparsityPattern> {
  const auto bra_shells = context.bra.shells();
  const auto ket_shells = context.ket.shells();
  const std::size_t atom_count = context.bra.atom_count();
  const bool diagonal = context.op.same_centre_blocks_diagonal();
  // Calls start(entry) before, and f(ket shell) for, the stored blocks of each neighbour entry
  // of local bra shell r of row atom `atom`.
  const auto for_each_block = [&](std::size_t atom, std::size_t r, auto&& start, auto&& f) {
    const std::size_t row = context.bra_first_shell[atom] + r;
    for (std::size_t e = neighbours.offsets[atom]; e < neighbours.offsets[atom + 1]; ++e) {
      const Neighbour& neighbour = neighbours.entries[e];
      start(e);
      const std::size_t first = context.ket_first_shell[neighbour.atom];
      const std::size_t count = context.ket_first_shell[neighbour.atom + 1] - first;
      if (neighbour.pair == same_atom) {
        for (std::size_t c = context.symmetric ? r : 0; c < count; ++c) {
          if (!diagonal ||
              ket_shells[first + c].angular_momentum == bra_shells[row].angular_momentum) {
            f(first + c);
          }
        }
        continue;
      }
      const BasisPairShells& table = *neighbour.table;
      const double distance_squared = neighbour.distance_squared;
      for (std::size_t c = 0; c < count; ++c) {
        const std::size_t s =
            neighbour.transposed ? c * table.ket_shell_count + r : r * table.ket_shell_count + c;
        if (distance_squared <= table.shell_pairs[s].cutoff) {
          f(first + c);
        }
      }
    }
  };

  std::vector<std::size_t> row_offsets(bra_shells.size() + 1, 0);
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    const std::size_t first = context.bra_first_shell[atom];
    for (std::size_t r = 0; r < context.bra_first_shell[atom + 1] - first; ++r) {
      std::size_t count = 0;
      for_each_block(atom, r, [](std::size_t) {}, [&](std::size_t) { ++count; });
      row_offsets[first + r + 1] = count;
    }
  }
  for (std::size_t row = 0; row < bra_shells.size(); ++row) {
    row_offsets[row + 1] += row_offsets[row];
  }
  neighbours.start_offsets.assign(atom_count + 1, 0);
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    neighbours.start_offsets[atom + 1] =
        neighbours.start_offsets[atom] +
        (neighbours.offsets[atom + 1] - neighbours.offsets[atom]) *
            (context.bra_first_shell[atom + 1] - context.bra_first_shell[atom]);
  }
  neighbours.row_starts.resize(neighbours.start_offsets.back());
  std::vector<std::uint32_t> block_kets(row_offsets.back());
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    const std::size_t first = context.bra_first_shell[atom];
    const std::size_t shells = context.bra_first_shell[atom + 1] - first;
    for (std::size_t r = 0; r < shells; ++r) {
      const std::size_t row_first = row_offsets[first + r];
      std::size_t block = row_first;
      for_each_block(
          atom, r,
          [&](std::size_t entry) {
            neighbours.row_starts[neighbours.row_start_index(atom, entry, shells, r)] =
                static_cast<std::uint32_t>(block - row_first);
          },
          [&](std::size_t ket_shell) {
            block_kets[block++] = static_cast<std::uint32_t>(ket_shell);
          });
    }
  }
  if (context.symmetric) {
    return SparsityPattern::from_rows(bra_shells, MatrixSymmetry::symmetric,
                                      diagonal ? DiagonalFormat::scalar : DiagonalFormat::full,
                                      std::move(row_offsets), std::move(block_kets));
  }
  return SparsityPattern::from_rows(bra_shells, ket_shells, std::move(row_offsets),
                                    std::move(block_kets));
}

constexpr std::size_t no_block = std::numeric_limits<std::size_t>::max();

/// Value offset of block (row shell, ket shell), or no_block if the pattern lacks it.
auto block_value_offset(const SparsityPattern& pattern, std::size_t row, std::size_t ket)
    -> std::size_t {
  const auto kets = pattern.block_kets();
  const auto first = kets.begin() + static_cast<std::ptrdiff_t>(pattern.row_offsets()[row]);
  const auto last = kets.begin() + static_cast<std::ptrdiff_t>(pattern.row_offsets()[row + 1]);
  const auto block = std::lower_bound(first, last, ket);
  return block != last && *block == ket
             ? pattern.block_value_offset(static_cast<std::size_t>(block - kets.begin()))
             : no_block;
}

/// Per-thread storage for off-diagonal blocks, reused across groups.
struct OffDiagonalWorkspace {
  std::vector<std::size_t> counts;        // significant atom pairs per shell pair
  std::vector<int> orders;                // l_a + l_b per shell pair
  std::vector<std::size_t> order_counts;  // atom pairs needing each solid-harmonic order
  SolidHarmonics harmonics;               // R, R^2 and S_{L,M}(R_AB) of the sorted pairs
  SolidHarmonics centre_harmonics;        // at R_AB = 0 (same-atom blocks by the kernels)
  std::vector<Point3D<double>> origins;   // R_AB = 0 per diagonal atom
  std::vector<AtomPair> diagonal_pairs;   // (A, A) per diagonal atom
  KernelWorkspace kernel;
  std::vector<double> values;        // N_A N_B (2l+1)(2l'+1) x count shell-pair values
  std::vector<std::size_t> rows;     // value row per block entry of the current shell pair
  std::vector<std::size_t> offsets;  // block value offset per (shell pair, sorted atom pair)
};

/// Fills offsets[s * n + i] with the value offset of the block of shell pair s and sorted atom
/// pair i (n sorted pairs) for the significant blocks (i < counts[s]): per atom pair and row
/// shell, a walk over the ket atom's blocks from the row start recorded with the pattern. Returns
/// false unless the walks find exactly the significant blocks.
auto find_block_offsets(const Context& context, const GroupPairs& group,
                        const NeighbourLists& neighbours, const SparsityPattern& pattern,
                        std::span<const std::size_t> counts, std::vector<std::size_t>& offsets)
    -> bool {
  const BasisPairShells& table = context.table(group.basis_pair);
  const std::size_t n = group.sorted.pairs.size();
  offsets.resize(table.shell_pairs.size() * n);
  const auto kets = pattern.block_kets();
  std::size_t found = 0;
  bool significant = true;
  for (std::size_t i = 0; i < n; ++i) {
    const AtomPair& pair = group.sorted.pairs[i];
    const bool transposed = is_transposed(pair, context.symmetric);
    const std::size_t row_atom = transposed ? pair.ket : pair.bra;
    const std::size_t ket_atom = transposed ? pair.bra : pair.ket;
    const std::size_t ket_first = context.ket_first_shell[ket_atom];
    const std::size_t ket_end = context.ket_first_shell[ket_atom + 1];
    const std::size_t row_first = context.bra_first_shell[row_atom];
    const std::size_t shells = context.bra_first_shell[row_atom + 1] - row_first;
    for (std::size_t r = 0; r < shells; ++r) {
      const std::size_t row = row_first + r;
      const std::size_t last = pattern.row_offsets()[row + 1];
      std::size_t block =
          pattern.row_offsets()[row] +
          neighbours.row_starts[neighbours.row_start_index(row_atom, group.entries[i], shells, r)];
      for (; block < last && kets[block] < ket_end; ++block) {
        const std::size_t c = kets[block] - ket_first;
        const std::size_t s =
            transposed ? c * table.ket_shell_count + r : r * table.ket_shell_count + c;
        significant = significant && i < counts[s];
        offsets[s * n + i] = pattern.block_value_offset(block);
        ++found;
      }
    }
  }
  std::size_t expected = 0;
  for (const std::size_t count : counts) {
    expected += count;
  }
  return significant && found == expected;
}

/// Row of the kernel values for bra function (I, m) and ket function (J, m') of a shell pair,
/// row-major over the block's functions: rows ordered by contracted pair (I, J), then component
/// (m, m'); mirrored shapes store only m <= m'.
void value_rows(const KernelShape& shape, std::vector<std::size_t>& rows) {
  const auto bra_components = static_cast<std::size_t>(2 * shape.l + 1);
  const auto ket_components = static_cast<std::size_t>(2 * shape.l_prime + 1);
  const std::size_t components = bra_components * ket_components;
  const std::size_t bra_functions = shape.bra_contractions * bra_components;
  const std::size_t ket_functions = shape.ket_contractions * ket_components;
  rows.resize(bra_functions * ket_functions);
  for (std::size_t bra = 0; bra < bra_functions; ++bra) {
    for (std::size_t ket = 0; ket < ket_functions; ++ket) {
      const std::size_t m = bra % bra_components;
      const std::size_t m_prime = ket % ket_components;
      const std::size_t component = shape.mirrored && m_prime < m ? m_prime * ket_components + m
                                                                  : m * ket_components + m_prime;
      rows[bra * ket_functions + ket] =
          ((bra / bra_components) * shape.ket_contractions + ket / ket_components) * components +
          component;
    }
  }
}

/// Writes the blocks of one shell pair for the first `count` sorted atom pairs into `matrix`
/// at `offsets` (one per atom pair) from the kernel values (see value_rows). Block entries are
/// ordered (I, m) x (J, m'); transposed pairs are stored as (J, m') x (I, m).
void write_blocks(const KernelShape& shape, std::size_t count, const SortedAtomPairs& sorted,
                  bool symmetric, std::span<const double> values,
                  std::span<const std::size_t> offsets, std::vector<std::size_t>& rows,
                  std::span<double> matrix) {
  const std::size_t bra_functions =
      shape.bra_contractions * static_cast<std::size_t>(2 * shape.l + 1);
  const std::size_t ket_functions =
      shape.ket_contractions * static_cast<std::size_t>(2 * shape.l_prime + 1);
  value_rows(shape, rows);
  for (std::size_t i = 0; i < count; ++i) {
    double* block = matrix.data() + offsets[i];
    if (is_transposed(sorted.pairs[i], symmetric)) {
      for (std::size_t ket = 0; ket < ket_functions; ++ket) {
        for (std::size_t bra = 0; bra < bra_functions; ++bra) {
          block[ket * bra_functions + bra] = values[rows[bra * ket_functions + ket] * count + i];
        }
      }
    } else {
      for (std::size_t bra = 0; bra < bra_functions; ++bra) {
        for (std::size_t ket = 0; ket < ket_functions; ++ket) {
          block[bra * ket_functions + ket] = values[rows[bra * ket_functions + ket] * count + i];
        }
      }
    }
  }
}

/// Kernel values of a shell pair for n atom pairs, in the workspace.
auto kernel_values(const ShellPairCutoff& shell_pair, const SolidHarmonics& harmonics,
                   std::span<const AtomPair> pairs, std::size_t n, OffDiagonalWorkspace& workspace)
    -> std::span<const double> {
  const KernelShape& shape = shell_pair.shape;
  const std::size_t size = shape.bra_contractions * shape.ket_contractions *
                           static_cast<std::size_t>((2 * shape.l + 1) * (2 * shape.l_prime + 1)) *
                           n;
  if (workspace.values.size() < size) {
    workspace.values.resize(size);
  }
  const std::span<double> values(workspace.values.data(), size);
  shell_pair.kernel->compute(harmonics, pairs.first(n), n, workspace.kernel, values);
  return values;
}

/// Off-diagonal blocks (A != B) of one group: counts of significant atom pairs per shell pair,
/// solid harmonics (each order for the pairs that need it), block offsets, then kernel and
/// blocks per shell pair. Returns the number of blocks written; sets `missing` if the pattern
/// lacks a block.
auto add_off_diagonal_blocks(const Context& context, const GroupPairs& group,
                             const NeighbourLists& neighbours, const SparsityPattern& pattern,
                             OffDiagonalWorkspace& workspace, std::span<double> matrix,
                             bool& missing) -> std::size_t {
  const SortedAtomPairs& sorted = group.sorted;
  if (sorted.pairs.empty()) {
    return 0;
  }
  const BasisPairShells& table = context.table(group.basis_pair);
  workspace.counts.clear();
  workspace.orders.clear();
  std::size_t max_count = 0;
  for (const ShellPairCutoff& shell_pair : table.shell_pairs) {
    workspace.counts.push_back(significant_pair_count(shell_pair.cutoff, sorted.distances_squared));
    workspace.orders.push_back(shell_pair.order);
    max_count = std::max(max_count, workspace.counts.back());
  }
  if (max_count == 0) {
    return 0;
  }
  counts_per_order(workspace.counts, workspace.orders, workspace.order_counts);
  workspace.harmonics.compute(sorted.separations, workspace.order_counts);
  if (!find_block_offsets(context, group, neighbours, pattern, workspace.counts,
                          workspace.offsets)) {
    missing = true;
    return 0;
  }

  const std::size_t n = sorted.pairs.size();
  std::size_t written = 0;
  for (std::size_t s = 0; s < table.shell_pairs.size(); ++s) {
    const std::size_t count = workspace.counts[s];
    if (count == 0) {
      continue;
    }
    const ShellPairCutoff& shell_pair = table.shell_pairs[s];
    const auto values =
        kernel_values(shell_pair, workspace.harmonics, sorted.pairs, count, workspace);
    write_blocks(shell_pair.shape, count, sorted, context.symmetric, values,
                 std::span(workspace.offsets).subspan(s * n, count), workspace.rows, matrix);
    written += count;
  }
  return written;
}

/// Same-atom blocks of one group: in symmetric matrices one scalar per contracted pair (k_a <=
/// k_b for a shell with itself), otherwise dense blocks over (k_a, m_a) x (k_b, m_b) with s(k_a,
/// k_b) where m_a == m_b. Returns the number of blocks written; sets `missing` if the pattern
/// lacks a block.
auto add_same_centre_blocks(const Context& context, const GroupPairs& group,
                            const SparsityPattern& pattern, SameCentreCache& cache,
                            std::span<double> matrix, bool& missing) -> std::size_t {
  if (group.diagonal_atoms.empty()) {
    return 0;
  }
  const AtomBasis& bra_basis = context.bra.unique_bases()[group.basis_pair.bra];
  const AtomBasis& ket_basis = context.ket.unique_bases()[group.basis_pair.ket];
  std::size_t written = 0;
  for (const SameCentreBlock& same :
       cache.blocks(context.op, bra_basis, ket_basis, group.basis_pair, context.symmetric)) {
    const bool same_shell = same.bra_shell == same.ket_shell;
    const auto components =
        static_cast<std::size_t>(2 * angular_momentum(bra_basis.shells()[same.bra_shell]) + 1);
    const std::size_t columns = same.ket_contractions * components;
    for (const std::size_t atom : group.diagonal_atoms) {
      const std::size_t offset =
          block_value_offset(pattern, context.bra_first_shell[atom] + same.bra_shell,
                             context.ket_first_shell[atom] + same.ket_shell);
      if (offset == no_block) {
        missing = true;
        continue;
      }
      double* values = matrix.data() + offset;
      ++written;
      if (context.symmetric) {
        std::size_t value = 0;
        for (std::size_t k_a = 0; k_a < same.bra_contractions; ++k_a) {
          for (std::size_t k_b = same_shell ? k_a : 0; k_b < same.ket_contractions; ++k_b) {
            values[value++] = same.values[k_a * same.ket_contractions + k_b];
          }
        }
        continue;
      }
      std::fill(values, values + same.bra_contractions * components * columns, 0.0);
      for (std::size_t k_a = 0; k_a < same.bra_contractions; ++k_a) {
        for (std::size_t k_b = 0; k_b < same.ket_contractions; ++k_b) {
          for (std::size_t m = 0; m < components; ++m) {
            values[(k_a * components + m) * columns + k_b * components + m] =
                same.values[k_a * same.ket_contractions + k_b];
          }
        }
      }
    }
  }
  return written;
}

/// Same-atom blocks of one group computed by the kernels at R_AB = 0 (operators whose same-atom
/// blocks are not diagonal in m): every shell pair of the atom, bra shell <= ket shell when
/// symmetric, never screened. A shell with itself stores the upper triangle over its functions
/// in symmetric matrices. Returns the number of blocks written; sets `missing` if the pattern
/// lacks a block.
auto add_same_atom_kernel_blocks(const Context& context, const GroupPairs& group,
                                 const SparsityPattern& pattern, OffDiagonalWorkspace& workspace,
                                 std::span<double> matrix, bool& missing) -> std::size_t {
  const std::size_t n = group.diagonal_atoms.size();
  if (n == 0) {
    return 0;
  }
  const BasisPairShells& table = context.table(group.basis_pair);
  workspace.origins.assign(n, Point3D<double>{0.0, 0.0, 0.0});
  workspace.diagonal_pairs.clear();
  for (const std::size_t atom : group.diagonal_atoms) {
    workspace.diagonal_pairs.push_back({atom, atom});
  }
  int max_order = 0;
  for (const ShellPairCutoff& shell_pair : table.shell_pairs) {
    max_order = std::max(max_order, shell_pair.order);
  }
  workspace.order_counts.assign(static_cast<std::size_t>(max_order) + 1, n);
  workspace.centre_harmonics.compute(workspace.origins, workspace.order_counts);

  const std::size_t bra_shells = table.shell_pairs.size() / table.ket_shell_count;
  std::size_t written = 0;
  for (std::size_t a = 0; a < bra_shells; ++a) {
    for (std::size_t b = context.symmetric ? a : 0; b < table.ket_shell_count; ++b) {
      const ShellPairCutoff& shell_pair = table.shell_pairs[a * table.ket_shell_count + b];
      const KernelShape& shape = shell_pair.shape;
      const auto values = kernel_values(shell_pair, workspace.centre_harmonics,
                                        workspace.diagonal_pairs, n, workspace);
      value_rows(shape, workspace.rows);
      const std::size_t bra_functions =
          shape.bra_contractions * static_cast<std::size_t>(2 * shape.l + 1);
      const std::size_t ket_functions =
          shape.ket_contractions * static_cast<std::size_t>(2 * shape.l_prime + 1);
      const bool triangle = context.symmetric && a == b;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t atom = group.diagonal_atoms[i];
        const std::size_t offset = block_value_offset(pattern, context.bra_first_shell[atom] + a,
                                                      context.ket_first_shell[atom] + b);
        if (offset == no_block) {
          missing = true;
          continue;
        }
        double* block = matrix.data() + offset;
        std::size_t value = 0;
        for (std::size_t bra = 0; bra < bra_functions; ++bra) {
          for (std::size_t ket = triangle ? bra : 0; ket < ket_functions; ++ket) {
            block[value++] = values[workspace.rows[bra * ket_functions + ket] * n + i];
          }
        }
        ++written;
      }
    }
  }
  return written;
}

/// Kernel of an uncontracted or segmented shell pair.
class SegmentedKernel final : public ShellPairKernel {
 public:
  explicit SegmentedKernel(SegmentedShellPair pair) : pair_(std::move(pair)) {}

  auto shape() const noexcept -> KernelShape override {
    return {pair_.l, pair_.l_prime, 1, 1, pair_.l == pair_.l_prime};
  }

  void compute(const SolidHarmonics& harmonics, std::span<const AtomPair>, std::size_t n,
               KernelWorkspace& workspace, std::span<double> values) const override {
    segmented_shell_pair_values(pair_, harmonics, n, workspace, values);
  }

 private:
  SegmentedShellPair pair_;
};

/// Kernel of a shell pair with a general contraction.
class GeneralKernel final : public ShellPairKernel {
 public:
  explicit GeneralKernel(GeneralShellPair pair) : pair_(std::move(pair)) {}

  auto shape() const noexcept -> KernelShape override {
    return {pair_.l, pair_.l_prime, pair_.bra_contractions, pair_.ket_contractions,
            pair_.l == pair_.l_prime};
  }

  void compute(const SolidHarmonics& harmonics, std::span<const AtomPair>, std::size_t n,
               KernelWorkspace& workspace, std::span<double> values) const override {
    general_shell_pair_values(pair_, harmonics, n, workspace, values);
  }

 private:
  GeneralShellPair pair_;
};

}  // namespace

auto make_integral_kernel(const BasisShell& bra, const BasisShell& ket, TwoCentreIntegral integral)
    -> std::unique_ptr<ShellPairKernel> {
  if (kind(bra) == ShellKind::general || kind(ket) == ShellKind::general) {
    return std::make_unique<GeneralKernel>(make_general_shell_pair(bra, ket, integral));
  }
  return std::make_unique<SegmentedKernel>(make_segmented_shell_pair(bra, ket, integral));
}

auto compute_two_centre(const TwoCentreOperator& op, const Molecule<double>& molecule,
                        const MolecularBasis& bra, const MolecularBasis& ket, bool symmetric,
                        double threshold, std::size_t block_size, PairCost cost)
    -> BlockSparseMatrix {
  // Two passes over the atom-pair groups: 1) sort each group's atom pairs by distance and keep
  // them, 2) build the pattern row by row from the resulting neighbour lists, 3) compute the
  // blocks and write them in place.
  if (block_size == 0) {
    const std::size_t atoms = molecule.size();
    const std::size_t pairs = symmetric ? atoms * (atoms + 1) / 2 : atoms * atoms;
    block_size = automatic_block_size(pairs, static_cast<std::size_t>(omp_get_max_threads()), cost);
  }
  const auto factory = symmetric ? std::make_unique<AtomPairGroupFactory>(bra, block_size)
                                 : std::make_unique<AtomPairGroupFactory>(bra, ket, block_size);
  Context context{op,
                  bra,
                  ket,
                  symmetric,
                  atom_shell_offsets(bra),
                  atom_shell_offsets(ket),
                  basis_pair_shells(op, bra, ket, threshold),
                  std::vector<GroupPairs>(factory->block_count())};

  // 1) Atom pairs within each group's largest cutoff, by distance.
#pragma omp parallel
  {
    AtomPairGroup group;
    std::vector<SortedAtomPairs::SortKey> keys;  // sort scratch, reused across groups
    while (factory->next(group)) {
      GroupPairs& stored = context.groups[group.index];
      stored.basis_pair = group.basis_pair;
      stored.diagonal_atoms = group.diagonal_atoms;
      if (!group.off_diagonal_pairs.empty()) {
        stored.sorted.keys.swap(keys);
        sort_by_distance(group.off_diagonal_pairs, molecule.coordinates(), stored.sorted,
                         context.table(group.basis_pair).max_cutoff);
        keys.swap(stored.sorted.keys);  // the group keeps no scratch
      }
    }
  }

  // 2) Pattern.
  auto neighbours = neighbour_lists(context.groups, bra.atom_count(), symmetric,
                                    [&](AtomBasisPair pair) { return &context.table(pair); });
  const auto pattern = matrix_pattern(context, neighbours);
  BlockSparseMatrix matrix(pattern);
  const std::span<double> values = matrix.values();

  // 3) Blocks, written in place.
  bool missing = false;
  std::size_t written = 0;
#pragma omp parallel reduction(|| : missing) reduction(+ : written)
  {
    SameCentreCache cache;
    OffDiagonalWorkspace workspace;
    const bool diagonal = op.same_centre_blocks_diagonal();
#pragma omp for schedule(dynamic, 1)
    for (std::size_t g = 0; g < context.groups.size(); ++g) {
      const GroupPairs& group = context.groups[g];
      written += diagonal ? add_same_centre_blocks(context, group, *pattern, cache, values, missing)
                          : add_same_atom_kernel_blocks(context, group, *pattern, workspace, values,
                                                        missing);
      written +=
          add_off_diagonal_blocks(context, group, neighbours, *pattern, workspace, values, missing);
    }
  }
  if (missing || written != pattern->block_count()) {
    throw std::logic_error("fika: sparsity pattern and computed blocks disagree");
  }
  return matrix;
}

}  // namespace fika::detail

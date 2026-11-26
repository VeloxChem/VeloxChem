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

#ifndef fika_two_centre_driver_hpp
#define fika_two_centre_driver_hpp

// Internal: driver shared by the two-centre one-electron integrals (overlap, kinetic energy, ...):
// shell pairs screened by an operator's bound, the sparsity pattern built first from the
// screened atom pairs, then the operator's kernels write the blocks in place.

#include <cstddef>
#include <memory>
#include <span>
#include <string_view>
#include <variant>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_molecular_basis.hpp"
#include "fika_atom_pair_group.hpp"
#include "fika_overlap_screening.hpp"
#include "fika_two_centre_kernels.hpp"
#include "fika_block_sparse_matrix.hpp"
#include "fika_solid_harmonics.hpp"
#include "fika_molecule.hpp"

namespace fika::detail {

/// Shape of the kernel values of a shell pair: angular momenta, contracted functions, and
/// whether blocks with l == l' are symmetric in m, m' (only m <= m' rows are then filled).
struct KernelShape {
  int l;
  int l_prime;
  std::size_t bra_contractions;
  std::size_t ket_contractions;
  bool mirrored;
};

/// Geometry-independent data of one shell pair (bra shell on A, ket shell on B) and its kernel.
class ShellPairKernel {
 public:
  virtual ~ShellPairKernel() = default;

  virtual auto shape() const noexcept -> KernelShape = 0;

  /// Integrals for the n atom pairs `pairs` (bra atom A, ket atom B) whose separations
  /// R_AB = A - B are those of `harmonics` (which covers the orders l + l'; R_AB = 0 for
  /// same-atom blocks): N_A N_B (2l + 1)(2l' + 1) rows of n entries, row
  /// ((I N_B + J)(2l + 1) + m + l)(2l' + 1) + m' + l' for contracted functions I, J and component
  /// (m, m'). For mirrored shapes (l == l') only rows with m <= m' are read.
  virtual void compute(const SolidHarmonics& harmonics, std::span<const AtomPair> pairs,
                       std::size_t n, KernelWorkspace& workspace,
                       std::span<double> values) const = 0;
};

/// A two-centre one-electron operator whose same-centre blocks are diagonal in m (rotationally
/// invariant operators such as the overlap and the kinetic energy).
class TwoCentreOperator {
 public:
  virtual ~TwoCentreOperator() = default;

  /// Upper bound of the shell pair's integrals as a function of R_AB^2 (screening).
  virtual auto bound(const BasisShell& bra, const BasisShell& ket) const -> ShellPairBound = 0;

  /// Kernel of a shell pair on different atoms.
  virtual auto kernel(const BasisShell& bra, const BasisShell& ket) const
      -> std::unique_ptr<ShellPairKernel> = 0;

  /// Whether same-atom blocks are diagonal in m and vanish between different l (true for
  /// rotationally invariant operators). If so they are computed by same_centre() and stored as
  /// scalars (symmetric matrices); otherwise every same-atom shell pair is stored in full and
  /// computed by the kernels at R_AB = 0.
  virtual auto same_centre_blocks_diagonal() const -> bool { return true; }

  /// Same-centre values v(k_a, k_b) of two shells of equal angular momentum on one atom,
  /// row-major (contracted functions of bra) x (of ket): the block is v times identity in m.
  /// Used only if same_centre_blocks_diagonal().
  virtual auto same_centre(const BasisShell& bra, const BasisShell& ket) const
      -> std::vector<double> = 0;
};

/// Kernel of `integral` for a shell pair on different atoms: early contraction for uncontracted
/// and segmented shells, matrix products if either shell is a general contraction.
auto make_integral_kernel(const BasisShell& bra, const BasisShell& ket, TwoCentreIntegral integral)
    -> std::unique_ptr<ShellPairKernel>;

/// Contracted values v(k_a, k_b) = sum_ij c_i(k_a) d_j(k_b) primitive(alpha_i, beta_j) of two
/// shells, row-major.
template <typename F>
auto contract_primitives(const BasisShell& bra, const BasisShell& ket, F&& primitive)
    -> std::vector<double> {
  return std::visit(
      [&](const auto& a, const auto& b) {
        const auto exponents_a = a.exponents();
        const auto exponents_b = b.exponents();
        std::vector<double> primitives(exponents_a.size() * exponents_b.size());
        for (std::size_t i = 0; i < exponents_a.size(); ++i) {
          for (std::size_t j = 0; j < exponents_b.size(); ++j) {
            primitives[i * exponents_b.size() + j] = primitive(exponents_a[i], exponents_b[j]);
          }
        }
        const std::size_t contractions_a = a.contraction_count();
        const std::size_t contractions_b = b.contraction_count();
        std::vector<double> contracted(contractions_a * contractions_b);
        for (std::size_t k_a = 0; k_a < contractions_a; ++k_a) {
          for (std::size_t k_b = 0; k_b < contractions_b; ++k_b) {
            double sum = 0.0;
            for_each_coefficient(a, k_a, [&](std::size_t i, double c_a) {
              const double* row = &primitives[i * exponents_b.size()];
              double inner = 0.0;
              for_each_coefficient(b, k_b,
                                   [&](std::size_t j, double c_b) { inner += c_b * row[j]; });
              sum += c_a * inner;
            });
            contracted[k_a * contractions_b + k_b] = sum;
          }
        }
        return contracted;
      },
      bra, ket);
}

/// Throws std::invalid_argument, prefixed with `caller`, for a negative or non-finite threshold or
/// a basis whose atom count differs from the molecule's.
void check_two_centre_input(std::string_view caller, const Molecule<double>& molecule,
                            const MolecularBasis& basis, double threshold);

/// Matrix of `op` between `bra` and `ket` (bases of `molecule`, already validated): symmetric
/// with scalar (or, for operators without diagonal same-atom blocks, full) same-centre blocks
/// when `symmetric` (bra and ket are the same basis), otherwise
/// rectangular. Off-diagonal blocks are screened with `threshold`; same-atom blocks never.
/// Work is handed out in blocks of `block_size` atom pairs; 0 picks automatic_block_size for
/// pairs of `cost` and the OpenMP thread count.
auto compute_two_centre(const TwoCentreOperator& op, const Molecule<double>& molecule,
                        const MolecularBasis& bra, const MolecularBasis& ket, bool symmetric,
                        double threshold, std::size_t block_size, PairCost cost)
    -> BlockSparseMatrix;

}  // namespace fika::detail

#endif  // fika_two_centre_driver_hpp

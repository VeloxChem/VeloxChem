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

#ifndef fika_overlap_screening_hpp
#define fika_overlap_screening_hpp

// Internal: screening of off-diagonal (different-atom) two-centre overlap blocks.

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_atom_pair_group.hpp"

namespace fika::detail {

/// Atom pairs sorted by increasing R_AB^2, with R_AB = R_A - R_B (A bra, B ket).
struct SortedAtomPairs {
  std::vector<AtomPair> pairs;
  std::vector<Point3D<double>> separations;  // R_AB
  std::vector<double> distances_squared;

  struct SortKey {
    double distance_squared;
    std::size_t index;  // position in the input: makes the order total and deterministic
  };
  std::vector<SortKey> keys;  // scratch
};

/// Sorts the pairs with R_AB^2 <= max_distance_squared by increasing distance (equal distances
/// keep their input order) into `sorted`, reusing its storage; farther pairs are dropped.
void sort_by_distance(std::span<const AtomPair> pairs, std::span<const Point3D<double>> coordinates,
                      SortedAtomPairs& sorted,
                      double max_distance_squared = std::numeric_limits<double>::infinity());

/// Simplified bound of one shell pair as a function of R = |A - B|,
///   B(R) = prefactor max(1, R^L) (c0 + c1 R^2) exp(-mu R^2),  L = l_a + l_b.
/// Overlap (c0 = 1, c1 = 0):
///   prefactor = max sum|c| max sum|d| (pi/p)^(3/2), mu = ab/p, p = a + b,
/// with a and b the smallest exponents of the two shells and the largest sums of |coefficient|
/// over their contracted functions. Kinetic energy: the same prefactor and mu times
/// mu'(3 + 2 mu' R^2), mu' = a'b'/(a' + b') from the largest exponents a', b'. Nuclear
/// attraction: see nuclear_attraction_shell_pair_bound. Only the decay
/// with R matters: it ranks atom pairs, it is not a rigorous bound.
struct ShellPairBound {
  double prefactor;
  double exponent;         // mu
  int angular_momentum;    // L
  double constant = 1.0;   // c0
  double quadratic = 0.0;  // c1

  auto operator()(double distance_squared) const noexcept -> double;

  /// Largest R^2 at which the bound still reaches `threshold`: beyond it the bound stays below
  /// the threshold. Negative if the bound never reaches it; infinite for threshold <= 0. The
  /// bound may rise with R before it decays (diffuse shells, R^L and R^2 factors) and may lie
  /// below the integrals of close pairs; screening by this cutoff only drops pairs in the
  /// decaying tail, where it covers them (measured within the threshold up to 0.1).
  auto cutoff_distance_squared(double threshold) const -> double;
};

/// Overlap bound of the shell pair (a, b).
auto shell_pair_bound(const BasisShell& a, const BasisShell& b) -> ShellPairBound;

/// Kinetic-energy bound of the shell pair (a, b).
auto kinetic_shell_pair_bound(const BasisShell& a, const BasisShell& b) -> ShellPairBound;

/// Nuclear-attraction bound of the shell pair (a, b) for point charges q_C of total magnitude
/// `charge_sum` = sum_C |q_C|, independent of their positions:
///   charge_sum max sum|c| max sum|d| (L + 1) (2 pi / p) max(1, R^L) exp(-mu R^2),
/// the unit-charge s-type integral (2 pi / p) exp(-mu R^2) F_0(p R_PC^2) with F_0 <= 1 (smallest
/// exponents, as for the overlap).
auto nuclear_attraction_shell_pair_bound(const BasisShell& a, const BasisShell& b,
                                         double charge_sum) -> ShellPairBound;

/// Upper bound on sqrt(T) F_1(T) over T >= 0 (the maximum 0.1896518... is reached at T ~ 0.937).
inline constexpr double dipole_boys_bound = 0.18966;

/// Bound of the shell pair (a, b) on the potential of point dipoles mu_C of total magnitude
/// `dipole_sum` = sum_C |mu_C|, independent of their positions:
///   dipole_sum max sum|c| max sum|d| (L + 1) (2 pi / p) 2 sqrt(p) s_1 max(1, R^L) exp(-mu R^2),
/// from the unit-dipole s-type integral |mu . grad_C (2 pi / p) exp(-mu R^2) F_0(p U^2)|
/// = (2 pi / p) exp(-mu R^2) 2 sqrt(p) sqrt(T) F_1(T), T = p U^2, with sqrt(T) F_1(T) <= s_1
/// (smallest exponents, as for the overlap).
auto dipole_potential_shell_pair_bound(const BasisShell& a, const BasisShell& b, double dipole_sum)
    -> ShellPairBound;

/// Upper bound on T F_2(T) over T >= 0 (the maximum 0.1086573... is reached at T ~ 1.58).
inline constexpr double quadrupole_boys_bound = 0.10866;

/// ||Theta||_2, the largest absolute eigenvalue of the traceless part Theta = Q - tr(Q)/3 of a
/// symmetric quadrupole (closed-form eigenvalues of a symmetric 3 x 3 matrix).
auto traceless_quadrupole_norm(const Quadrupole& quadrupole) -> double;

/// Bound of the shell pair (a, b) on the potential of point quadrupoles of total norm
/// `quadrupole_sum` = sum_C ||Theta_C||_2 (traceless parts), independent of their positions
/// (quadrupole notes, Eq. 26):
///   quadrupole_sum max sum|c| max sum|d| (L + 1) 4 pi s_2 max(1, R^L) exp(-mu R^2),
/// from the unit s-type integral 1/2 Theta : grad_C grad_C (2 pi / p) exp(-mu R^2) F_0(p U^2)
/// = 4 pi exp(-mu R^2) p F_2(p U^2) U . Theta . U with |U . Theta . U| <= ||Theta||_2 U^2 and
/// T F_2(T) <= s_2 (smallest exponents, as for the overlap).
auto quadrupole_potential_shell_pair_bound(const BasisShell& a, const BasisShell& b,
                                           double quadrupole_sum) -> ShellPairBound;

/// Number of leading pairs (sorted by increasing distance) within the cutoff distance, i.e. up
/// to the last pair whose bound reaches the threshold (binary search on R^2).
auto significant_pair_count(double cutoff_distance_squared,
                            std::span<const double> distances_squared) -> std::size_t;

/// Separations needing each order of the solid harmonics: result[L] is the largest count of the
/// shell pairs with l_a + l_b >= L (non-increasing in L), for L = 0..max(orders).
/// counts[i] and orders[i] belong to shell pair i.
void counts_per_order(std::span<const std::size_t> counts, std::span<const int> orders,
                      std::vector<std::size_t>& result);

}  // namespace fika::detail

#endif  // fika_overlap_screening_hpp

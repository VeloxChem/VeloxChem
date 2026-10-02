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

#ifndef fika_two_centre_kernels_hpp
#define fika_two_centre_kernels_hpp

// Internal: kernels of two-centre one-electron integrals (overlap and kinetic energy) over shell
// pairs on different atoms: early contraction for uncontracted and segmented shells (overlap
// notes Sections 4.8 and 6.3, kinetic notes Section 5) and matrix products for shell pairs
// involving general contractions (overlap notes Section 7, kinetic notes Section 6).
//
// Both integrals share one form. With p = a + b, mu = ab/p, n = l + l', k_max = min(l, l') and
// the primitive factors
//   w_ab,j = (-1)^l b^l a^l' p^-n (pi/p)^(3/2) mu_ab^-j exp(-mu_ab R^2),
// contracted to W_j (early: sum_ab c_a d_b w_ab,j; general: c^T w_j d), the shell-pair block is
//   w0 S_lm(R) S_l'm'(R) + sum_{L < n} T_L Phi^(L)_mm'(R),  Phi^(L) = sum_M C^{LM}_{lm,l'm'} S_LM,
// with T_L = sum_{t=1..deg} r_t R^(2(deg - t)) W(t) (Horner in R^2). Accumulator rows index the
// powers j = first_power.. of mu^-j; row t is W_{t + first_power}.
//  - Overlap: j = 0..k_max, w0 = W_0, deg = k, r_t = (-1)^t D^(t)_kL.
//  - Kinetic: j = -2..k_max - 1, w0 = -2 R^2 W_-2 + 2 (n + 3/2) W_-1, deg = k + 1,
//    r_1 = 2 D^(1)_kL, r_t = -2 (-1)^t D^(t)_(k+1)L (t >= 2).

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_solid_harmonics.hpp"
#include "fika_multipole_expansion.hpp"

namespace fika::detail {

enum class TwoCentreIntegral { overlap, kinetic };

/// Operator-dependent coefficients shared by both kernel kinds.
struct IntegralForm {
  TwoCentreIntegral integral = TwoCentreIntegral::overlap;
  int first_power = 0;          // mu^-j powers j = first_power .. first_power + power_count - 1
  std::size_t power_count = 0;  // accumulator rows
  double w0_quadratic = 0.0;    // w0 = w0_quadratic R^2 W(0) + w0_constant W(w0_constant_row)
  double w0_constant = 1.0;
  std::size_t w0_constant_row = 0;
  std::vector<int> radial_offsets;          // per L < n: start in radial_coefficients (+1)
  std::vector<double> radial_coefficients;  // r_t, t = 1..deg, per L
};

/// Geometry-independent data of a shell pair (bra shell l on A, ket shell l' on B), both
/// uncontracted or segmented: mu_ab and the factors q_ab,j = c_a d_b w_ab,j / exp(-mu_ab R^2)
/// (power_count rows of primitive_pairs).
struct SegmentedShellPair {
  int l = 0;
  int l_prime = 0;
  std::size_t primitive_pairs = 0;
  std::vector<double> mu;
  std::vector<double> factors;
  IntegralForm form;
};

/// Builds the shell-pair data; both shells must be uncontracted or segmented.
auto make_segmented_shell_pair(const BasisShell& bra, const BasisShell& ket,
                               TwoCentreIntegral integral) -> SegmentedShellPair;

/// Reusable scratch storage of the integral kernels.
struct KernelWorkspace {
  std::vector<double> exponentials;  // primitive pairs x n
  std::vector<double> accumulators;  // W: power_count x n; general: power_count N_A N_B x n
  std::vector<double> radial;        // T_L: couplings x n
  std::vector<double> leading;       // w0 over n (kinetic)
  std::vector<double> weighted;      // general: primitive factors w_ab,j over n
  std::vector<double> half;          // general: half-transformed factors
  std::vector<double> kernels;       // general: Psi_j, power_count (2l + 1)(2l' + 1) x n
  std::vector<double> powers;        // general: R^(2p), p = 1.., over n
  // Nuclear attraction.
  SolidHarmonics charge_harmonics;                  // S_LM(W) of near charges
  std::vector<Point3D<double>> charge_separations;  // W = C - A of near charges
  std::vector<double> boys;                         // F_0..F_n
  std::vector<double> field;                        // Y^(L)_M,ab accumulators
  std::vector<double> contracted;                   // contracted field
  std::vector<double> charge_weights;               // far weights q/U^(2L+1), 1/U^2 per charge
  std::vector<double> far_field;                    // Phi_LM of the far charges
  std::vector<std::uint32_t> near_charges;          // charges with p U^2 < T_far
  std::vector<double> translation_bra;   // G^(lambda)_m,mu(R) over n (factorized assembly)
  std::vector<double> translation_ket;   // H^(lambda')_m',mu'(R) over n
  std::vector<double> gaunt_contracted;  // Z^(lambda lambda')_mu,mu' over n
  std::vector<double> bra_transformed;   // W^(lambda')_m,mu' over n
  ExpansionWorkspace expansion;          // far-field evaluation (L2P)
  std::vector<double> gathered_charges;  // near charges of a far-field leaf
  std::vector<Point3D<double>> gathered_coordinates;
  // Dipoles.
  SolidHarmonics dipole_harmonics;                  // S_LM(U) of the dipoles, L <= n + 1
  std::vector<Point3D<double>> dipole_separations;  // U = D - P
  std::vector<double> dipole_weights;               // far weights mu_(m) / U^(2L+3), 1/U^2
  std::vector<std::uint32_t> near_dipoles;          // dipoles with p U^2 < T_far(n + 1)
  std::vector<double> dipole_couplings;             // A+-_LM of one dipole, and Psi
  std::vector<Dipole> gathered_dipoles;             // near dipoles of a far-field leaf
  std::vector<Point3D<double>> gathered_dipole_coordinates;
  std::vector<Point3D<double>> gathered_quadrupole_coordinates;
  // Quadrupoles.
  SolidHarmonics quadrupole_harmonics;                  // S_LM(U), L <= n + 2
  std::vector<Point3D<double>> quadrupole_separations;  // U = E - P
  std::vector<std::array<double, 5>> thetas;            // theta_mu per gathered quadrupole
  std::vector<double> quadrupole_weights;               // far weights theta_mu / U^(2L+5), 1/U^2
  std::vector<std::uint32_t> near_quadrupoles;          // p U^2 < T_far(n + 2)
  std::vector<double> quadrupole_couplings;             // B^(L+2), B^L, B^(L-2) of one site; Psi
};

/// Integrals of a segmented shell pair for the first n separations of `harmonics` (which covers
/// the orders l + l'). `values` holds (2l + 1)(2l' + 1) rows of n entries, row (m + l)(2l' + 1) +
/// (m' + l') for component (m, m'). For l == l' only rows with m <= m' are written (the block is
/// symmetric in m, m'). A single primitive pair takes a simplified path.
void segmented_shell_pair_values(const SegmentedShellPair& pair, const SolidHarmonics& harmonics,
                                 std::size_t n, KernelWorkspace& workspace,
                                 std::span<double> values);

/// Geometry-independent data of a shell pair (bra shell l on A with K_A primitives and N_A
/// contracted functions, ket shell l' on B with K_B and N_B) of which at least one is a general
/// contraction; uncontracted and segmented shells take part with N = 1: mu_ab (a-major), the
/// factors w_ab,j / exp(-mu_ab R^2) without coefficients (power_count rows of K_A K_B), and the
/// dense coefficient matrices c (K_A x N_A) and d (K_B x N_B) that contract them to
/// W_j = c^T w_j d.
struct GeneralShellPair {
  int l = 0;
  int l_prime = 0;
  std::size_t bra_primitives = 0;    // K_A
  std::size_t ket_primitives = 0;    // K_B
  std::size_t bra_contractions = 0;  // N_A
  std::size_t ket_contractions = 0;  // N_B
  bool bra_first = true;             // contract the bra primitives first (the cheaper order)
  std::vector<double> mu;
  std::vector<double> factors;
  std::vector<double> bra_coefficients;  // c^T: N_A x K_A, row-major
  std::vector<double> ket_coefficients;  // d^T: N_B x K_B, row-major
  IntegralForm form;
};

/// Builds the shell-pair data of any two shells (intended for pairs with a general contraction).
auto make_general_shell_pair(const BasisShell& bra, const BasisShell& ket,
                             TwoCentreIntegral integral) -> GeneralShellPair;

/// Integrals of a shell pair with general contractions for the first n separations of
/// `harmonics`: sum_t W(t)_IJ Psi_t,mm' with the kernels Psi_t collecting, per accumulator row,
/// the w0 and T_L terms of the common form. `values` holds N_A N_B (2l + 1)(2l' + 1) rows of n
/// entries, row ((I N_B + J)(2l + 1) + m + l)(2l' + 1) + m' + l' for contracted functions I, J
/// and component (m, m'). For l == l' only rows with m <= m' are written.
void general_shell_pair_values(const GeneralShellPair& pair, const SolidHarmonics& harmonics,
                               std::size_t n, KernelWorkspace& workspace, std::span<double> values);

}  // namespace fika::detail

#endif  // fika_two_centre_kernels_hpp

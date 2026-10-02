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

#ifndef fika_mm_induced_dipoles_hpp
#define fika_mm_induced_dipoles_hpp

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_summation.hpp"
#include "fika_classical_system.hpp"
#include "fika_dipole_interaction.hpp"
#include "fika_electric_field.hpp"
#include "fika_polarizable_sites.hpp"

namespace fika {

struct InducedDipoleOptions {
  /// Convergence: largest site residual |F_i - (B mu)_i| (field a.u.), B = alpha^-1 - T.
  double tolerance = 1e-8;
  std::size_t max_iterations = 100;
  /// Starting dipoles (one per site, e.g. of the previous SCF iteration); empty: mu = alpha F.
  std::vector<Point3D<double>> initial_guess;
  /// Rescale the initial guess by the s minimizing the CG functional Q(s mu_g), Q(mu) =
  /// 1/2 mu . B mu - mu . F, i.e. s = (mu_g . F) / (mu_g . B mu_g), and start from alpha F instead
  /// when s is not positive or the scaled guess has the higher Q (false: use the guess as given).
  bool scale_initial_guess = true;
  /// induced_dipoles only: direct sums, fast multipole methods, or automatic (the permanent field
  /// from multipole_field_sites sites, the dipole coupling from multipole_coupling_sites).
  ChargeSummation summation = ChargeSummation::automatic;
  /// induced_dipoles only: field error target of the fast multipole methods; 0: tolerance / 10.
  double field_accuracy = 0.0;
  /// induced_dipoles only: include the field of the permanent charges (false: the external
  /// contributions alone, e.g. of a perturbed QM density in response theory).
  bool permanent_field = true;
};

/// Polarizable sites from which induced_dipoles uses the fast multipole methods automatically
/// (measured crossovers of full solves on Ahlstrom water droplets, 14 threads: the FMM coupling
/// breaks even near 80000 sites and wins 1.2x at 160000).
inline constexpr std::size_t multipole_field_sites = 160000;
inline constexpr std::size_t multipole_coupling_sites = 120000;

struct InducedDipoles {
  std::vector<Point3D<double>> dipoles;  // one per site (a.u.)
  /// induced_dipoles only: the field F the dipoles respond to (permanent charges and external
  /// contributions), one per site.
  std::vector<Point3D<double>> field;
  std::size_t iterations = 0;  // operator applications in the CG loop
  double residual = 0.0;       // largest site residual of the returned dipoles
  double rms_residual = 0.0;   // root mean square site residual
  bool converged = false;
  /// Whether the CG started from the initial guess, and the factor it was scaled by (1 when used
  /// as given; 0 when no guess was given or it was replaced by alpha F).
  bool guess_used = false;
  double guess_scale = 0.0;
  /// induced_dipoles only: how the permanent field and the coupling were summed (direct or
  /// multipole) and the final FMM expansion orders (0 for direct sums).
  ChargeSummation field_summation = ChargeSummation::direct;
  ChargeSummation coupling_summation = ChargeSummation::direct;
  int field_order = 0;
  int coupling_order = 0;
};

/// Solves (alpha^-1 - T) mu = F for the induced dipoles by conjugate gradients with the
/// block-Jacobi preconditioner alpha. Converged when the largest site residual, recomputed from
/// the dipoles (not only recursively updated), is at most the tolerance. Independent of the thread
/// count. Throws std::invalid_argument for inconsistent sizes, a nonpositive tolerance, an initial
/// guess that is not finite, or a polarizability that is not positive definite, and
/// std::runtime_error when B turns out not to be positive definite (p . B p <= 0 for a search
/// direction or a nonzero guess: a polarization catastrophe, e.g. of undamped close sites).
auto solve_induced_dipoles(const PolarizableSites& sites, std::span<const Point3D<double>> field,
                           const DipoleInteraction& interaction,
                           const InducedDipoleOptions& options = {}) -> InducedDipoles;

/// Induced dipoles of the polarizable region of `system` in the field of its permanent charges,
/// summed directly (PermanentChargeField, DirectDipoleInteraction) or through fast multipole
/// methods (FmmChargeField, FmmDipoleInteraction with dipole scale max alpha * max |F|) as
/// options.summation selects, plus the `external` contributions (e.g. of a QM region) added in
/// order; the dipoles are in the order of polarizable_sites(system). The external contributions
/// see the sites of polarizable_sites(system) and should meet options.field_accuracy themselves.
auto induced_dipoles(const ClassicalSystem& system, std::optional<TholeDamping> damping,
                     const InducedDipoleOptions& options = {},
                     std::span<const FieldContribution* const> external = {}) -> InducedDipoles;

/// Induction energy -1/2 sum_i mu_i . F_i (Hartree).
auto induction_energy(std::span<const Point3D<double>> dipoles,
                      std::span<const Point3D<double>> field) -> double;

}  // namespace fika

#endif  // fika_mm_induced_dipoles_hpp

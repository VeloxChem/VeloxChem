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

#ifndef fika_multipole_expansion_hpp
#define fika_multipole_expansion_hpp

// Internal: multipole and local expansions of the potential of point charges (the operators of
// a fast multipole method). Expansions use the scaled complex solid harmonics (Helgaker,
// Molecular Electronic-Structure Theory, 9.13)
//   R_lm(r) = C_lm(r) / sqrt((l - m)! (l + m)!),  I_lm(r) = sqrt((l - m)! (l + m)!) C_lm(r) /
//   r^(2l + 1),
// with C_lm the complex Racah-normalized solid harmonics, which obey
//   1 / |r1 - r2| = sum_lm conj(R_lm(r1)) I_lm(r2)                        (|r1| < |r2|),
//   R_lm(a + b) = sum_jk R_jk(a) R_(l-j),(m-k)(b),
//   I_lm(a + b) = sum_jk (-1)^j conj(R_jk(a)) I_(l+j),(m+k)(b)            (|a| < |b|).
// About a centre O, the multipole expansion of charges q at s is M_lm = sum q conj(R_lm(s - O))
// and the potential phi(x) = sum_lm M_lm I_lm(x - O); about a centre Q, the local expansion is
// phi(x) = sum_lm L_lm conj(R_lm(x - Q)) with L_lm = sum q I_lm(s - Q) for direct charges. The
// charges are real, so X_(l,-m) = (-1)^m conj(X_lm) for every expansion X and only m >= 0 is
// stored, at index l (l + 1) / 2 + m.
//
// An expansion of order p keeps l <= p; multipole-to-local keeps total degree j + l <= p, so a
// field tensor of rank Lambda evaluated from it holds every term of degree <= p - Lambda in the
// combined displacement u = (s - O) - (x - Q) (see field_tensor_error_bound).

#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"

namespace fika::detail {

/// Highest supported expansion order.
inline constexpr int max_expansion_order = 64;

/// Coefficients of an expansion of order p (m >= 0).
constexpr auto expansion_size(int order) noexcept -> std::size_t {
  return static_cast<std::size_t>((order + 1) * (order + 2) / 2);
}

constexpr auto expansion_index(int l, int m) noexcept -> std::size_t {
  return static_cast<std::size_t>(l * (l + 1) / 2 + m);
}

/// R_lm(r), l <= order, m >= 0.
void regular_harmonics(const Point3D<double>& r, int order, std::span<std::complex<double>> values);

/// I_lm(r), l <= order, m >= 0 (r != 0).
void irregular_harmonics(const Point3D<double>& r, int order,
                         std::span<std::complex<double>> values);

/// Reusable scratch storage of the expansion operators.
struct ExpansionWorkspace {
  std::vector<std::complex<double>> harmonics;
  std::vector<std::complex<double>> shifted;
  std::vector<double> real;  // expansions with all components m = -l..l
  std::vector<double> imaginary;
  std::vector<double> kernel_real;
  std::vector<double> kernel_imaginary;
};

/// P2M: adds the charges to the multipole expansion (order p) about `centre`.
void add_charges_to_multipole(std::span<const double> charges,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order,
                              std::span<std::complex<double>> multipole,
                              ExpansionWorkspace& workspace);

/// Harmonics R_jk(from - to) of a multipole translation (all components m), for reuse by many
/// translations with the same shift.
struct MultipoleShift {
  int order = 0;
  std::vector<double> real;
  std::vector<double> imaginary;
};

auto make_multipole_shift(const Point3D<double>& from, const Point3D<double>& to, int order)
    -> MultipoleShift;

/// P2M of point quadrupoles (primitive moments): adds 1/2 Theta : grad_s grad_s of each charge's
/// multipole (the trace of Q drops out).
void add_quadrupoles_to_multipole(std::span<const Quadrupole> quadrupoles,
                                  std::span<const Point3D<double>> coordinates,
                                  const Point3D<double>& centre, int order,
                                  std::span<std::complex<double>> multipole,
                                  ExpansionWorkspace& workspace);

/// M2M with a precomputed shift (order p of the shift).
void translate_multipole(std::span<const std::complex<double>> source, const MultipoleShift& shift,
                         std::span<std::complex<double>> target, ExpansionWorkspace& workspace);

/// P2M of point dipoles: adds mu . grad_s of each charge's multipole (the exact derivative of
/// add_charges_to_multipole with respect to the source position).
void add_dipoles_to_multipole(std::span<const Dipole> dipoles,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order,
                              std::span<std::complex<double>> multipole,
                              ExpansionWorkspace& workspace);

/// P2M of a point multipole of rank k at `position` given by real moments Q_LM (L <= k, rows
/// L^2 + M + L, fika's real solid harmonics S): its potential sum_LM Q_LM S_LM(x - position) /
/// |x - position|^(2L + 1) (for charges Q_LM = sum q S_LM(s - position)). As complex coefficients
/// about the position, M_l0 = Q_l0 / l! and M_lm = (-1)^m (Q_l,m - i Q_l,-m) / (sqrt 2
/// sqrt((l - m)! (l + m)!)) for m > 0, translated to the expansion (order p) about `centre` over
/// the source degrees l <= k alone (exact up to order p).
void add_real_multipole_to_multipole(std::span<const double> moments, int rank,
                                     const Point3D<double>& position, const Point3D<double>& centre,
                                     int order, std::span<std::complex<double>> multipole,
                                     ExpansionWorkspace& workspace);

/// M2M: adds the multipole expansion about `from` to the one about `to` (both order p; exact).
void translate_multipole(std::span<const std::complex<double>> source, const Point3D<double>& from,
                         const Point3D<double>& to, int order,
                         std::span<std::complex<double>> target, ExpansionWorkspace& workspace);

/// M2L: adds the multipole expansion about `from` to the local expansion (order p, total
/// degree <= p) about `to`.
void multipole_to_local(std::span<const std::complex<double>> multipole,
                        const Point3D<double>& from, const Point3D<double>& to, int order,
                        std::span<std::complex<double>> local, ExpansionWorkspace& workspace);

/// P2L: adds the charges directly to the local expansion (order p) about `centre`.
void add_charges_to_local(std::span<const double> charges,
                          std::span<const Point3D<double>> coordinates,
                          const Point3D<double>& centre, int order,
                          std::span<std::complex<double>> local, ExpansionWorkspace& workspace);

/// P2L of point dipoles: adds mu . grad_s of each charge's local expansion (order p; uses
/// irregular harmonics of order p + 1 <= max_expansion_order).
void add_dipoles_to_local(std::span<const Dipole> dipoles,
                          std::span<const Point3D<double>> coordinates,
                          const Point3D<double>& centre, int order,
                          std::span<std::complex<double>> local, ExpansionWorkspace& workspace);

/// P2L of point quadrupoles: adds 1/2 Theta : grad_s grad_s of each charge's local expansion
/// (order p; uses irregular harmonics of order p + 2 <= max_expansion_order).
void add_quadrupoles_to_local(std::span<const Quadrupole> quadrupoles,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order,
                              std::span<std::complex<double>> local, ExpansionWorkspace& workspace);

/// L2L: adds the local expansion about `from` to the one about `to` (both order p; exact).
void translate_local(std::span<const std::complex<double>> source, const Point3D<double>& from,
                     const Point3D<double>& to, int order, std::span<std::complex<double>> target,
                     ExpansionWorkspace& workspace);

/// L2P: the field tensor Phi_Lambda,M(x) = sum q S_Lambda,M(s - x) / |s - x|^(2 Lambda + 1)
/// (fika's real solid harmonics S), Lambda <= rank <= order, from the local expansion (order p)
/// about `centre`; `phi` holds (rank + 1)^2 entries, row Lambda^2 + M + Lambda.
void local_field_tensor(std::span<const std::complex<double>> local, const Point3D<double>& centre,
                        int order, const Point3D<double>& point, int rank, std::span<double> phi,
                        ExpansionWorkspace& workspace);

/// Bound on the error of Phi_Lambda (Euclidean norm over M) evaluated through multipole
/// expansions (order p) of charges within `source_radius` of O, translated to local expansions
/// about Q with |Q - O| = distance, at points within `target_radius` of Q:
///   sum|q| C(p + 1, Lambda) theta^(p - Lambda + 1) / ((1 - theta)^(Lambda + 1) d^(Lambda + 1)),
/// theta = (source_radius + target_radius) / d < 1. It sums the truncated terms of the
/// expansion of Phi_Lambda in u, whose degree-N part is at most C(N + Lambda, N) |u|^N /
/// d^(N + Lambda + 1) (the Gegenbauer bound).
auto field_tensor_error_bound(int order, int rank, double charge_sum, double source_radius,
                              double target_radius, double distance) -> double;

/// The same bound for point dipoles of total magnitude sum |mu| within `source_radius` of O:
/// the charge error is harmonic in the source position, so the gradient estimate for harmonic
/// functions gives |mu| (3 / delta) times the charge bound for radius source_radius + delta,
/// minimized over delta = gap / 2^k (k = 1..10), gap = d - source_radius - target_radius > 0.
auto dipole_field_tensor_error_bound(int order, int rank, double dipole_sum, double source_radius,
                                     double target_radius, double distance) -> double;

/// Charge-bound tolerance per unit dipole at ball radius delta: delta / 3 (the dipole bound is the
/// charge bound for radius source_radius + delta divided by it). Shared with FarFieldExpansion.
inline auto dipole_tolerance_scale(double delta) noexcept -> double {
  return delta / 3.0;
}

/// Charge-bound tolerance per unit quadrupole at ball radius delta: delta^2 3 sqrt 3 / 20.
inline auto quadrupole_tolerance_scale(double delta) noexcept -> double {
  return delta * delta * 3.0 * std::sqrt(3.0) / 20.0;
}

/// The same bound for point quadrupoles of total norm sum ||Theta||_2 within `source_radius` of
/// O: the charge error is harmonic in the source position, and its second derivatives obey
/// |1/2 Theta : grad grad e| <= (20 / (3 sqrt 3)) ||Theta||_2 sup |e| / delta^2 over a ball of
/// radius delta (degree-2 spherical-harmonic projection), so the bound is (20 / (3 sqrt 3)) /
/// delta^2 times the charge bound for radius source_radius + delta, minimized over
/// delta = gap / 2^k (k = 1..10), gap = d - source_radius - target_radius > 0.
auto quadrupole_field_tensor_error_bound(int order, int rank, double quadrupole_sum,
                                         double source_radius, double target_radius,
                                         double distance) -> double;

/// The same bound for point multipoles given by real moments (add_real_multipole_to_multipole)
/// within `source_radius` of O, with moment_sums[k] = sum over sources of ||Q_k||_2 (Euclidean
/// norm over M of the rank-k moments), k <= 2 max_angular_momentum. The charge error e(s) is
/// harmonic in the source position and a rank-k source picks (-1)^k sum_M Q_kM b_M from its
/// degree-k part h_k(y) = sum_M b_M S_kM(y) about the source (S_kM(grad) S_kM'(y) =
/// (2k - 1)!! delta_MM'). On the sphere |y| = delta, |h_k| <= L_k sup|e| (Funk-Hecke projection,
/// L_k = (2k + 1) / 2 integral_-1^1 |P_k|) and ||b||_2 <= sqrt(2k + 1) max|h_k| / delta^k, so the
/// rank-k sources add ||Q_k|| sqrt(2k + 1) L_k / delta^k times the charge bound for radius
/// source_radius + delta, minimized over delta = gap / 2^j (j = 1..10), gap = d - source_radius -
/// target_radius > 0; rank 0 is the charge bound.
auto real_multipole_field_tensor_error_bound(int order, int rank,
                                             std::span<const double> moment_sums,
                                             double source_radius, double target_radius,
                                             double distance) -> double;

}  // namespace fika::detail

#endif  // fika_multipole_expansion_hpp

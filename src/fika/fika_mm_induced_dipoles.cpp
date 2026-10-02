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

#include "fika_mm_induced_dipoles.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "fika_compensated_sum.hpp"
#include "fika_electric_field.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::solve_induced_dipoles: " + reason);
}

using Vector = std::vector<Point3D<double>>;

/// Symmetric 3 x 3 matrix as its 6 packed components (xx, xy, xz, yy, yz, zz).
using Symmetric = std::array<double, 6>;

auto times(const Symmetric& m, const Point3D<double>& v) -> Point3D<double> {
  return {m[0] * v.x + m[1] * v.y + m[2] * v.z, m[1] * v.x + m[3] * v.y + m[4] * v.z,
          m[2] * v.x + m[4] * v.y + m[5] * v.z};
}

/// Inverse of a positive definite polarizability; throws if it is not.
auto inverse(const Polarizability& alpha, std::size_t site) -> Symmetric {
  const auto& a = alpha.components;
  const double m00 = a[3] * a[5] - a[4] * a[4];
  const double m01 = a[2] * a[4] - a[1] * a[5];
  const double m02 = a[1] * a[4] - a[2] * a[3];
  const double m11 = a[0] * a[5] - a[2] * a[2];
  const double m12 = a[1] * a[2] - a[0] * a[4];
  const double m22 = a[0] * a[3] - a[1] * a[1];
  const double det = a[0] * m00 + a[1] * m01 + a[2] * m02;
  if (!(a[0] > 0.0 && m22 > 0.0 && det > 0.0) || !std::isfinite(det)) {
    fail("polarizability of site " + std::to_string(site) + " is not positive definite");
  }
  return {m00 / det, m01 / det, m02 / det, m11 / det, m12 / det, m22 / det};
}

auto dot(const Vector& a, const Vector& b) -> double {
  detail::NeumaierSum sum;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum.add(a[i].x * b[i].x + a[i].y * b[i].y + a[i].z * b[i].z);
  }
  return sum.value();
}

struct ResidualNorms {
  double largest = 0.0;
  double rms = 0.0;
};

auto norms(const Vector& r) -> ResidualNorms {
  detail::NeumaierSum squares;
  double largest = 0.0;
  for (const auto& v : r) {
    const double s = v.x * v.x + v.y * v.y + v.z * v.z;
    squares.add(s);
    largest = std::max(largest, s);
  }
  return {std::sqrt(largest), r.empty() ? 0.0 : std::sqrt(squares.value() / double(r.size()))};
}

}  // namespace

auto solve_induced_dipoles(const PolarizableSites& sites, std::span<const Point3D<double>> field,
                           const DipoleInteraction& interaction,
                           const InducedDipoleOptions& options) -> InducedDipoles {
  const std::size_t n = sites.positions.size();
  if (sites.polarizabilities.size() != n || field.size() != n) {
    fail("the field and polarizabilities need one entry per site");
  }
  if (!options.initial_guess.empty() && options.initial_guess.size() != n) {
    fail("the initial guess needs one entry per site");
  }
  for (const auto& mu : options.initial_guess) {
    if (!std::isfinite(mu.x) || !std::isfinite(mu.y) || !std::isfinite(mu.z)) {
      fail("the initial guess must be finite");
    }
  }
  if (!(options.tolerance > 0.0) || !std::isfinite(options.tolerance)) {
    fail("the tolerance must be positive and finite");
  }
  std::vector<Symmetric> alpha(n), alpha_inverse(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& c = sites.polarizabilities[i].components;
    alpha[i] = {c[0], c[1], c[2], c[3], c[4], c[5]};
    alpha_inverse[i] = inverse(sites.polarizabilities[i], i);
  }
  const Vector f(field.begin(), field.end());

  // q = B p = alpha^-1 p - T p.
  Vector t(n);
  const auto apply_b = [&](const Vector& p, Vector& q) {
    interaction.apply(p, t);
    for (std::size_t i = 0; i < n; ++i) {
      const auto a = times(alpha_inverse[i], p[i]);
      q[i] = {a.x - t[i].x, a.y - t[i].y, a.z - t[i].z};
    }
  };
  const auto precondition = [&](const Vector& r, Vector& z) {
    for (std::size_t i = 0; i < n; ++i) {
      z[i] = times(alpha[i], r[i]);
    }
  };
  // True residual r = F - B x.
  Vector q(n);
  const auto residual = [&](const Vector& x, Vector& r) {
    apply_b(x, q);
    for (std::size_t i = 0; i < n; ++i) {
      r[i] = {f[i].x - q[i].x, f[i].y - q[i].y, f[i].z - q[i].z};
    }
  };

  const auto not_positive_definite = [] {
    throw std::runtime_error(
        "fika::solve_induced_dipoles: B = alpha^-1 - T is not positive definite (polarization "
        "catastrophe; damping may be needed)");
  };

  InducedDipoles result;
  Vector& x = result.dipoles;
  x.resize(n);
  precondition(f, x);  // alpha F
  if (!options.initial_guess.empty() && !options.scale_initial_guess) {
    x = options.initial_guess;
    result.guess_used = true;
    result.guess_scale = 1.0;
  } else if (!options.initial_guess.empty()) {
    // Q(s g) is smallest at s = (g . F) / (g . B g), where Q = -1/2 s (g . F); keep it when below
    // Q(alpha F).
    const Vector& g = options.initial_guess;
    Vector bg(n);
    apply_b(g, bg);
    const double gbg = dot(g, bg);
    const bool zero = std::ranges::all_of(
        g, [](const Point3D<double>& mu) { return mu.x == 0.0 && mu.y == 0.0 && mu.z == 0.0; });
    if (!zero && !(gbg > 0.0)) {
      not_positive_definite();
    }
    if (!zero) {
      const double gf = dot(g, f);
      const double scale = gf / gbg;
      Vector bx(n);
      apply_b(x, bx);
      const double q_default = 0.5 * dot(x, bx) - dot(x, f);
      if (scale > 0.0 && std::isfinite(scale) && -0.5 * scale * gf < q_default) {
        for (std::size_t i = 0; i < n; ++i) {
          x[i] = {scale * g[i].x, scale * g[i].y, scale * g[i].z};
        }
        result.guess_used = true;
        result.guess_scale = scale;
      }
    }
  }
  Vector r(n), z(n), p(n);
  residual(x, r);
  ResidualNorms norm = norms(r);
  while (norm.largest > options.tolerance && result.iterations < options.max_iterations) {
    // (Re)start from the true residual.
    precondition(r, z);
    p = z;
    double rz = dot(r, z);
    while (result.iterations < options.max_iterations) {
      apply_b(p, q);
      ++result.iterations;
      const double curvature = dot(p, q);
      if (!(curvature > 0.0)) {
        not_positive_definite();
      }
      const double step = rz / curvature;
      for (std::size_t i = 0; i < n; ++i) {
        x[i] = {x[i].x + step * p[i].x, x[i].y + step * p[i].y, x[i].z + step * p[i].z};
        r[i] = {r[i].x - step * q[i].x, r[i].y - step * q[i].y, r[i].z - step * q[i].z};
      }
      if (norms(r).largest <= options.tolerance) {
        break;
      }
      precondition(r, z);
      const double rz_next = dot(r, z);
      const double beta = rz_next / rz;
      rz = rz_next;
      for (std::size_t i = 0; i < n; ++i) {
        p[i] = {z[i].x + beta * p[i].x, z[i].y + beta * p[i].y, z[i].z + beta * p[i].z};
      }
    }
    residual(x, r);  // the recursive residual may drift
    norm = norms(r);
  }
  result.residual = norm.largest;
  result.rms_residual = norm.rms;
  result.converged = norm.largest <= options.tolerance;
  return result;
}

auto induced_dipoles(const ClassicalSystem& system, std::optional<TholeDamping> damping,
                     const InducedDipoleOptions& options,
                     std::span<const FieldContribution* const> external) -> InducedDipoles {
  const PolarizableSites sites = polarizable_sites(system);
  const std::size_t n = sites.positions.size();
  const double accuracy =
      options.field_accuracy > 0.0 ? options.field_accuracy : 0.1 * options.tolerance;
  const auto multipole = [&](std::size_t crossover) {
    return options.summation == ChargeSummation::multipole ||
           (options.summation == ChargeSummation::automatic && n >= crossover);
  };

  std::vector<Point3D<double>> field(n, Point3D<double>{0.0, 0.0, 0.0});
  ChargeSummation field_summation = ChargeSummation::direct;
  int field_order = 0;
  if (options.permanent_field && multipole(multipole_field_sites)) {
    const FmmChargeField source(system, {.absolute_accuracy = accuracy});
    source.add_field(sites, field);
    field_summation = ChargeSummation::multipole;
    field_order = source.order();
  } else if (options.permanent_field) {
    PermanentChargeField(system).add_field(sites, field);
  }
  for (const FieldContribution* contribution : external) {
    if (contribution == nullptr) {
      throw std::invalid_argument("fika::induced_dipoles: null external field contribution");
    }
    contribution->add_field(sites, field);
  }

  // Dipole scale of the coupling: the starting dipoles alpha F, or the initial guess if larger.
  double largest_field = 0.0, largest_alpha = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    largest_field = std::max(largest_field, std::hypot(field[i].x, field[i].y, field[i].z));
    for (const double a : sites.polarizabilities[i].components) {
      largest_alpha = std::max(largest_alpha, std::abs(a));
    }
  }
  double largest_guess = 0.0;
  for (const auto& mu : options.initial_guess) {
    largest_guess = std::max(largest_guess, std::hypot(mu.x, mu.y, mu.z));
  }
  const double dipole_scale =
      std::max(std::sqrt(3.0) * largest_alpha * largest_field, largest_guess);
  InducedDipoles result;
  if (multipole(multipole_coupling_sites) && dipole_scale > 0.0) {
    const FmmDipoleInteraction coupling(
        sites, damping, {.absolute_accuracy = accuracy, .dipole_scale = dipole_scale});
    result = solve_induced_dipoles(sites, field, coupling, options);
    result.coupling_summation = ChargeSummation::multipole;
    result.coupling_order = coupling.order();
  } else {
    result = solve_induced_dipoles(sites, field, DirectDipoleInteraction(sites, damping), options);
  }
  result.field_summation = field_summation;
  result.field_order = field_order;
  result.field = std::move(field);
  return result;
}

auto induction_energy(std::span<const Point3D<double>> dipoles,
                      std::span<const Point3D<double>> field) -> double {
  if (dipoles.size() != field.size()) {
    throw std::invalid_argument("fika::induction_energy: dipole and field counts differ");
  }
  detail::NeumaierSum sum;
  for (std::size_t i = 0; i < dipoles.size(); ++i) {
    sum.add(dipoles[i].x * field[i].x + dipoles[i].y * field[i].y + dipoles[i].z * field[i].z);
  }
  return -0.5 * sum.value();
}

}  // namespace fika

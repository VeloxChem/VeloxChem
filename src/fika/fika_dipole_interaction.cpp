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

#include "fika_dipole_interaction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace fika {

namespace {

/// Largest 1 - lambda of the damping kept (relative change of the undamped tensor).
constexpr double damping_tolerance = 1e-15;

/// v beyond which (1 + v + v^2/2 + v^3/6) e^-v (>= 1 - lambda_3) is below damping_tolerance.
auto damping_cutoff() -> double {
  double low = 0.0, high = 100.0;
  for (int i = 0; i < 200; ++i) {
    const double v = 0.5 * (low + high);
    ((1.0 + v + v * v / 2.0 + v * v * v / 6.0) * std::exp(-v) > damping_tolerance ? low : high) = v;
  }
  return high;
}

/// y += c5 3 r (r . mu) / r^5 - c3 mu / r^3 for r = (dx, dy, dz).
inline void add_pair(double dx, double dy, double dz, double c3, double c5, const double* mu,
                     double* y) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double inverse3 = 1.0 / (r2 * std::sqrt(r2));
  const double projection = 3.0 * c5 * (dx * mu[0] + dy * mu[1] + dz * mu[2]) * inverse3 / r2;
  y[0] += projection * dx - c3 * inverse3 * mu[0];
  y[1] += projection * dy - c3 * inverse3 * mu[1];
  y[2] += projection * dz - c3 * inverse3 * mu[2];
}

}  // namespace

detail::DipoleCorrections::DipoleCorrections(const PolarizableSites& sites,
                                             std::optional<TholeDamping> damping,
                                             const char* caller)
    : damping_(damping) {
  const std::size_t n = sites.positions.size();
  if (sites.polarizabilities.size() != n || sites.owners.size() != n) {
    throw std::invalid_argument(std::string(caller) + ": inconsistent site arrays");
  }
  if (damping && !(damping->a > 0.0 && std::isfinite(damping->a))) {
    throw std::invalid_argument(std::string(caller) + ": damping parameter must be positive");
  }
  x_.resize(n);
  y_.resize(n);
  z_.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Point3D<double>& p = sites.positions[i];
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) {
      throw std::invalid_argument(std::string(caller) + ": site " + std::to_string(i) +
                                  " has a non-finite position");
    }
    x_[i] = p.x;
    y_[i] = p.y;
    z_[i] = p.z;
  }

  // Residue groups: sites by increasing owner, in site order within an owner (polarizable_sites
  // emits them that way already).
  owner_sites_.resize(n);
  std::iota(owner_sites_.begin(), owner_sites_.end(), std::uint32_t{0});
  if (!std::ranges::is_sorted(sites.owners)) {
    std::ranges::stable_sort(owner_sites_, {}, [&](std::uint32_t i) { return sites.owners[i]; });
  }
  owner_group_.resize(n);
  owner_offsets_.push_back(0);
  for (std::size_t k = 0; k < n; ++k) {
    if (k > 0 && sites.owners[owner_sites_[k]] != sites.owners[owner_sites_[k - 1]]) {
      owner_offsets_.push_back(k);
    }
    owner_group_[owner_sites_[k]] = owner_offsets_.size() - 1;
  }
  if (n > 0) {
    owner_offsets_.push_back(n);
  }

  if (!damping || n == 0) {
    return;
  }
  double largest = 0.0;
  scales_.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& alpha = sites.polarizabilities[i];
    const double mean = (alpha(0, 0) + alpha(1, 1) + alpha(2, 2)) / 3.0;
    if (!(mean > 0.0) || !std::isfinite(mean)) {
      throw std::invalid_argument(std::string(caller) + ": site " + std::to_string(i) +
                                  " needs a positive polarizability for damping");
    }
    scales_[i] = std::pow(mean, 1.0 / 6.0);
    largest = std::max(largest, scales_[i]);
  }
  // v = a r / (s_i s_j) <= damping_cutoff() for r <= cutoff.
  cutoff_ = damping_cutoff() * largest * largest / damping->a;

  Point3D<double> high{-std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};
  origin_ = {-high.x, -high.y, -high.z};
  for (std::size_t i = 0; i < n; ++i) {
    origin_ = {std::min(origin_.x, x_[i]), std::min(origin_.y, y_[i]), std::min(origin_.z, z_[i])};
    high = {std::max(high.x, x_[i]), std::max(high.y, y_[i]), std::max(high.z, z_[i])};
  }
  const double extent[3] = {high.x - origin_.x, high.y - origin_.y, high.z - origin_.z};
  // Cells of edge cutoff_, widened when the sites are so spread out (or the cutoff so small) that
  // the grid would exceed max_cells: neighbouring cells still hold every pair within the cutoff.
  const double max_cells = std::max(double{1 << 20}, 8.0 * static_cast<double>(n));
  cell_edge_ = cutoff_;
  const auto cell_count = [&] {
    return (std::floor(extent[0] / cell_edge_) + 1.0) * (std::floor(extent[1] / cell_edge_) + 1.0) *
           (std::floor(extent[2] / cell_edge_) + 1.0);
  };
  while (cell_count() > max_cells) {
    cell_edge_ *= 1.01 * std::cbrt(cell_count() / max_cells);
  }
  for (int d = 0; d < 3; ++d) {
    cells_[d] = static_cast<std::size_t>(extent[d] / cell_edge_) + 1;
  }
  const auto cell_of = [&](std::size_t i) {
    const std::size_t cx =
        std::min(cells_[0] - 1, static_cast<std::size_t>((x_[i] - origin_.x) / cell_edge_));
    const std::size_t cy =
        std::min(cells_[1] - 1, static_cast<std::size_t>((y_[i] - origin_.y) / cell_edge_));
    const std::size_t cz =
        std::min(cells_[2] - 1, static_cast<std::size_t>((z_[i] - origin_.z) / cell_edge_));
    return (cx * cells_[1] + cy) * cells_[2] + cz;
  };
  cell_offsets_.assign(cells_[0] * cells_[1] * cells_[2] + 1, 0);
  for (std::size_t i = 0; i < n; ++i) {
    ++cell_offsets_[cell_of(i) + 1];
  }
  for (std::size_t c = 0; c + 1 < cell_offsets_.size(); ++c) {
    cell_offsets_[c + 1] += cell_offsets_[c];
  }
  cell_sites_.resize(n);
  std::vector<std::size_t> next(cell_offsets_.begin(), cell_offsets_.end() - 1);
  for (std::size_t i = 0; i < n; ++i) {
    cell_sites_[next[cell_of(i)]++] = static_cast<std::uint32_t>(i);
  }
}

DirectDipoleInteraction::DirectDipoleInteraction(const PolarizableSites& sites,
                                                 std::optional<TholeDamping> damping)
    : corrections_(sites, damping, "fika::DirectDipoleInteraction") {}

void DirectDipoleInteraction::apply(std::span<const Point3D<double>> mu,
                                    std::span<Point3D<double>> y) const {
  const std::size_t n = corrections_.size();
  if (mu.size() != n || y.size() != n) {
    throw std::invalid_argument("fika::DirectDipoleInteraction: vectors need one entry per site");
  }
  std::vector<double> mx(n), my(n), mz(n), fx(n), fy(n), fz(n);
  for (std::size_t j = 0; j < n; ++j) {
    mx[j] = mu[j].x;
    my[j] = mu[j].y;
    mz[j] = mu[j].z;
  }
  const double* px = corrections_.x();
  const double* py = corrections_.y();
  const double* pz = corrections_.z();
  const double* qx = mx.data();
  const double* qy = my.data();
  const double* qz = mz.data();
  double* ox = fx.data();
  double* oy = fy.data();
  double* oz = fz.data();

  // Undamped T mu: each site sums the field of every other site itself, j in increasing order, so
  // sites run in parallel without shared writes or synchronization and every site receives its
  // contributions in a fixed order (independent of the thread count). This does each pair twice
  // (T_ij = T_ji is not used), but the inner loop is a plain reduction that vectorizes.
  const auto row = [&](std::size_t i, std::size_t j0, std::size_t j1, double& sx, double& sy,
                       double& sz) {
    const double xi = px[i], yi = py[i], zi = pz[i];
    for (std::size_t j = j0; j < j1; ++j) {
      const double dx = xi - px[j], dy = yi - py[j], dz = zi - pz[j];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double inverse3 = 1.0 / (r2 * std::sqrt(r2));
      const double inverse5 = 3.0 * inverse3 / r2;
      const double pj = (dx * qx[j] + dy * qy[j] + dz * qz[j]) * inverse5;
      sx += pj * dx - inverse3 * qx[j];
      sy += pj * dy - inverse3 * qy[j];
      sz += pj * dz - inverse3 * qz[j];
    }
  };
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t i = 0; i < n; ++i) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    row(i, 0, i, sx, sy, sz);
    row(i, i + 1, n, sx, sy, sz);
    ox[i] = sx;
    oy[i] = sy;
    oz[i] = sz;
  }

  corrections_.apply(qx, qy, qz, ox, oy, oz, y);
}

void detail::DipoleCorrections::apply(const double* qx, const double* qy, const double* qz,
                                      const double* ox, const double* oy, const double* oz,
                                      std::span<Point3D<double>> y) const {
  const std::size_t n = x_.size();
  const double* px = x_.data();
  const double* py = y_.data();
  const double* pz = z_.data();
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t i = 0; i < n; ++i) {
    const double xi = px[i], yi = py[i], zi = pz[i];
    double field[3] = {ox[i], oy[i], oz[i]};
    // Same-residue pairs removed.
    const std::size_t group = owner_group_[i];
    for (std::size_t k = owner_offsets_[group]; k < owner_offsets_[group + 1]; ++k) {
      const std::size_t j = owner_sites_[k];
      if (j != i) {
        const double m[3] = {-qx[j], -qy[j], -qz[j]};
        add_pair(xi - px[j], yi - py[j], zi - pz[j], 1.0, 1.0, m, field);
      }
    }

    // Damping corrections (lambda - 1) for close pairs of different residues.
    if (damping_) {
      const double a = damping_->a;
      const double cutoff2 = cutoff_ * cutoff_;
      const auto cell = [&](double value, double low, std::size_t count) {
        return std::min(count - 1, static_cast<std::size_t>((value - low) / cell_edge_));
      };
      const std::size_t c[3] = {cell(xi, origin_.x, cells_[0]), cell(yi, origin_.y, cells_[1]),
                                cell(zi, origin_.z, cells_[2])};
      for (std::size_t cx = c[0] > 0 ? c[0] - 1 : 0; cx <= std::min(c[0] + 1, cells_[0] - 1);
           ++cx) {
        for (std::size_t cy = c[1] > 0 ? c[1] - 1 : 0; cy <= std::min(c[1] + 1, cells_[1] - 1);
             ++cy) {
          for (std::size_t cz = c[2] > 0 ? c[2] - 1 : 0; cz <= std::min(c[2] + 1, cells_[2] - 1);
               ++cz) {
            const std::size_t index = (cx * cells_[1] + cy) * cells_[2] + cz;
            for (std::size_t k = cell_offsets_[index]; k < cell_offsets_[index + 1]; ++k) {
              const std::size_t j = cell_sites_[k];
              if (j == i || owner_group_[j] == group) {
                continue;
              }
              const double dx = xi - px[j], dy = yi - py[j], dz = zi - pz[j];
              const double r2 = dx * dx + dy * dy + dz * dz;
              if (r2 >= cutoff2) {
                continue;
              }
              const double v = a * std::sqrt(r2) / (scales_[i] * scales_[j]);
              const double e = std::exp(-v);
              const double p2 = 1.0 + v + v * v / 2.0;
              // lambda3 - 1 and lambda5 - 1 directly (no cancellation against 1).
              const double m[3] = {qx[j], qy[j], qz[j]};
              add_pair(dx, dy, dz, -p2 * e, -(p2 + v * v * v / 6.0) * e, m, field);
            }
          }
        }
      }
    }
    y[i] = {field[0], field[1], field[2]};
  }
}

FmmDipoleInteraction::FmmDipoleInteraction(const PolarizableSites& sites,
                                           std::optional<TholeDamping> damping,
                                           const FmmDipoleOptions& options)
    : corrections_(sites, damping, "fika::FmmDipoleInteraction"),
      positions_(sites.positions),
      options_(options) {
  if (!(options.dipole_scale > 0.0) || !std::isfinite(options.dipole_scale) ||
      !(options.absolute_accuracy > 0.0) || !std::isfinite(options.absolute_accuracy)) {
    throw std::invalid_argument(
        "fika::FmmDipoleInteraction: dipole scale and accuracy must be positive");
  }
  build(options.order);
}

void FmmDipoleInteraction::build(int order) const {
  const std::vector<Point3D<double>> none;
  fmm_ = std::make_unique<detail::VolumeFmm>(
      none, positions_, positions_,
      detail::VolumeFmmOptions{.absolute_accuracy = options_.absolute_accuracy,
                               .order = order,
                               .dipole_scale = options_.dipole_scale});
}

void FmmDipoleInteraction::apply(std::span<const Point3D<double>> mu,
                                 std::span<Point3D<double>> y) const {
  const std::size_t n = corrections_.size();
  if (mu.size() != n || y.size() != n) {
    throw std::invalid_argument("fika::FmmDipoleInteraction: vectors need one entry per site");
  }
  std::vector<Dipole> dipoles(n);
  std::vector<double> mx(n), my(n), mz(n);
  for (std::size_t j = 0; j < n; ++j) {
    dipoles[j] = {{mu[j].x, mu[j].y, mu[j].z}};
    mx[j] = mu[j].x;
    my[j] = mu[j].y;
    mz[j] = mu[j].z;
  }
  const std::vector<double> no_charges;
  std::vector<Point3D<double>> undamped(n);
  if (!checked_) {
    // First product (the starting dipoles, of the scale the order was chosen for): raise the
    // order while the sampled error exceeds a tenth of the accuracy.
    while (true) {
      detail::VolumeFmmReport report;
      fmm_->field(no_charges, dipoles, undamped, &report);
      if (report.sampled_error <= 0.1 * options_.absolute_accuracy ||
          fmm_->order() + 2 > detail::largest_automatic_fmm_order) {
        break;
      }
      build(fmm_->order() + 2);
      ++retries_;
    }
    checked_ = true;
  } else {
    fmm_->field(no_charges, dipoles, undamped);
  }
  std::vector<double> ux(n), uy(n), uz(n);
  for (std::size_t i = 0; i < n; ++i) {
    ux[i] = undamped[i].x;
    uy[i] = undamped[i].y;
    uz[i] = undamped[i].z;
  }
  corrections_.apply(mx.data(), my.data(), mz.data(), ux.data(), uy.data(), uz.data(), y);
}

}  // namespace fika

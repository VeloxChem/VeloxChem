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
#include <map>
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
    x_[i] = sites.positions[i].x;
    y_[i] = sites.positions[i].y;
    z_[i] = sites.positions[i].z;
  }

  // Residue groups.
  std::map<std::size_t, std::vector<std::uint32_t>> groups;
  for (std::size_t i = 0; i < n; ++i) {
    groups[sites.owners[i]].push_back(static_cast<std::uint32_t>(i));
  }
  owner_group_.resize(n);
  owner_offsets_.push_back(0);
  for (const auto& [owner, members] : groups) {
    for (const std::uint32_t i : members) {
      owner_group_[i] = owner_offsets_.size() - 1;
    }
    owner_sites_.insert(owner_sites_.end(), members.begin(), members.end());
    owner_offsets_.push_back(owner_sites_.size());
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
  for (int d = 0; d < 3; ++d) {
    cells_[d] = static_cast<std::size_t>(extent[d] / cutoff_) + 1;
  }
  const auto cell_of = [&](std::size_t i) {
    const std::size_t cx =
        std::min(cells_[0] - 1, static_cast<std::size_t>((x_[i] - origin_.x) / cutoff_));
    const std::size_t cy =
        std::min(cells_[1] - 1, static_cast<std::size_t>((y_[i] - origin_.y) / cutoff_));
    const std::size_t cz =
        std::min(cells_[2] - 1, static_cast<std::size_t>((z_[i] - origin_.z) / cutoff_));
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

namespace {

/// Blocks of the symmetric pair sweep: about 112 (so a round has ~56 independent block pairs),
/// 16..512 sites each; depends only on the site count.
auto sweep_blocks(std::size_t n) -> std::size_t {
  const std::size_t size = std::clamp<std::size_t>((n + 111) / 112, 16, 512);
  return (n + size - 1) / size;
}

}  // namespace

void DirectDipoleInteraction::apply(std::span<const Point3D<double>> mu,
                                    std::span<Point3D<double>> y) const {
  const std::size_t n = corrections_.size();
  if (mu.size() != n || y.size() != n) {
    throw std::invalid_argument("fika::DirectDipoleInteraction: vectors need one entry per site");
  }
  std::vector<double> mx(n), my(n), mz(n), fx(n, 0.0), fy(n, 0.0), fz(n, 0.0);
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

  // Undamped T mu over all pairs, each pair once: T_ij = T_ji (r r^T is even in r), so a pair adds
  // T mu_j to site i and T mu_i to site j. Blocks of sites are paired in a round-robin tournament:
  // the block pairs of a round are disjoint, so they run in parallel and every site receives its
  // contributions in a fixed order (independent of the thread count).
  const std::size_t blocks = sweep_blocks(n);
  const std::size_t block_size = blocks > 0 ? (n + blocks - 1) / blocks : 0;
  const auto first = [&](std::size_t b) { return std::min(n, b * block_size); };
  // Sites [i0, i1) with [j0, j1); within one block (diagonal) only j > i.
  const auto sweep = [&](std::size_t i0, std::size_t i1, std::size_t j0, std::size_t j1,
                         bool diagonal, double* tx, double* ty, double* tz) {
    for (std::size_t i = i0; i < i1; ++i) {
      const double xi = px[i], yi = py[i], zi = pz[i];
      const double ux = qx[i], uy = qy[i], uz = qz[i];
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (std::size_t j = diagonal ? i + 1 : j0; j < j1; ++j) {
        const double dx = xi - px[j], dy = yi - py[j], dz = zi - pz[j];
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double inverse3 = 1.0 / (r2 * std::sqrt(r2));
        const double inverse5 = 3.0 * inverse3 / r2;
        const double pj = (dx * qx[j] + dy * qy[j] + dz * qz[j]) * inverse5;
        const double pi = (dx * ux + dy * uy + dz * uz) * inverse5;
        sx += pj * dx - inverse3 * qx[j];
        sy += pj * dy - inverse3 * qy[j];
        sz += pj * dz - inverse3 * qz[j];
        tx[j - j0] += pi * dx - inverse3 * ux;
        ty[j - j0] += pi * dy - inverse3 * uy;
        tz[j - j0] += pi * dz - inverse3 * uz;
      }
      ox[i] += sx;
      oy[i] += sy;
      oz[i] += sz;
    }
  };
  const std::size_t players = blocks + (blocks % 2);  // an odd count gets an idle player
#pragma omp parallel
  {
    std::vector<double> tx(block_size), ty(block_size), tz(block_size);
    const auto add_to = [&](std::size_t j0, std::size_t j1) {
      for (std::size_t j = j0; j < j1; ++j) {
        ox[j] += tx[j - j0];
        oy[j] += ty[j - j0];
        oz[j] += tz[j - j0];
      }
    };
#pragma omp for schedule(dynamic, 1)
    for (std::size_t b = 0; b < blocks; ++b) {
      const std::size_t i0 = first(b), i1 = first(b + 1);
      std::fill_n(tx.begin(), i1 - i0, 0.0);
      std::fill_n(ty.begin(), i1 - i0, 0.0);
      std::fill_n(tz.begin(), i1 - i0, 0.0);
      sweep(i0, i1, i0, i1, true, tx.data(), ty.data(), tz.data());
      add_to(i0, i1);
    }
    for (std::size_t round = 0; round + 1 < players; ++round) {
#pragma omp for schedule(dynamic, 1)
      for (std::size_t k = 0; k < players / 2; ++k) {
        const std::size_t a = k == 0 ? players - 1 : (round + k) % (players - 1);
        const std::size_t b = (round + players - 1 - k) % (players - 1);
        if (a >= blocks || b >= blocks) {
          continue;  // the idle player
        }
        const std::size_t i0 = first(std::min(a, b)), i1 = first(std::min(a, b) + 1);
        const std::size_t j0 = first(std::max(a, b)), j1 = first(std::max(a, b) + 1);
        std::fill_n(tx.begin(), j1 - j0, 0.0);
        std::fill_n(ty.begin(), j1 - j0, 0.0);
        std::fill_n(tz.begin(), j1 - j0, 0.0);
        sweep(i0, i1, j0, j1, false, tx.data(), ty.data(), tz.data());
        add_to(j0, j1);
      }
    }
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
        return std::min(count - 1, static_cast<std::size_t>((value - low) / cutoff_));
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
              const double lambda3 = 1.0 - p2 * e;
              const double lambda5 = 1.0 - (p2 + v * v * v / 6.0) * e;
              const double m[3] = {qx[j], qy[j], qz[j]};
              add_pair(dx, dy, dz, lambda3 - 1.0, lambda5 - 1.0, m, field);
            }
          }
        }
      }
    }
    y[i] = {field[0], field[1], field[2]};
  }
}

namespace {

/// Largest automatic order of the volume FMM (see volume_fmm.hpp).
constexpr int largest_fmm_order = 26;

}  // namespace

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
          fmm_->order() + 2 > largest_fmm_order) {
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

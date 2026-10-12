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

#include "fika_electric_field.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "fika_volume_fmm.hpp"

namespace fika {

namespace {

/// Adds sum_{c in [begin, end)} q_c S(|s - C|) (s - C) / |s - C|^3 to `field`, with the inner
/// switching factor S = 1 when `switching` is unset. The weight check also covers a charge at
/// `s` (r = 0), which the switch auto-excludes (S(0) = 0) like the FMM near field excludes a
/// source at the target point itself.
void add_charge_field(std::span<const double> charges, std::span<const Point3D<double>> positions,
                      std::size_t begin, std::size_t end, const Point3D<double>& s,
                      const std::optional<InnerSwitch>& switching, Point3D<double>& field) {
  double fx = 0.0, fy = 0.0, fz = 0.0;
  for (std::size_t c = begin; c < end; ++c) {
    const double dx = s.x - positions[c].x;
    const double dy = s.y - positions[c].y;
    const double dz = s.z - positions[c].z;
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double r = std::sqrt(r2);
    const double weight = switching ? inner_switch(r, *switching) : 1.0;
    if (weight == 0.0) {
      continue;
    }
    const double scale = weight * charges[c] / (r2 * r);
    fx += scale * dx;
    fy += scale * dy;
    fz += scale * dz;
  }
  field.x += fx;
  field.y += fy;
  field.z += fz;
}

}  // namespace

PermanentChargeField::PermanentChargeField(const ClassicalSystem& system,
                                           std::optional<InnerSwitch> switching)
    : charges_(classical_charges(system)), switching_(switching) {
  if (switching_) {
    check_inner_switch(*switching_, "fika::PermanentChargeField");
  }
}

void PermanentChargeField::add_field(const PolarizableSites& sites,
                                     std::span<Point3D<double>> field) const {
  const std::size_t n = sites.positions.size();
  if (field.size() != n) {
    throw std::invalid_argument("fika::PermanentChargeField: field has " +
                                std::to_string(field.size()) + " entries for " + std::to_string(n) +
                                " sites");
  }
  if (sites.owners.size() != n) {
    throw std::invalid_argument("fika::PermanentChargeField: " + std::to_string(sites.owners.size()) +
                                " owners for " + std::to_string(n) + " sites");
  }
  const std::size_t residues = charges_.residue_offsets.size() - 1;
  for (const std::size_t owner : sites.owners) {
    if (owner >= residues) {
      throw std::invalid_argument("fika::PermanentChargeField: site of residue " +
                                  std::to_string(owner) + " outside the classical system");
    }
  }
  const std::span<const double> q = charges_.charges;
  const std::span<const Point3D<double>> positions = charges_.coordinates;
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t i = 0; i < n; ++i) {
    // Charges outside the site's own residue: [0, begin) and [end, all).
    const std::size_t begin = charges_.residue_offsets[sites.owners[i]];
    const std::size_t end = charges_.residue_offsets[sites.owners[i] + 1];
    add_charge_field(q, positions, 0, begin, sites.positions[i], switching_, field[i]);
    add_charge_field(q, positions, end, q.size(), sites.positions[i], switching_, field[i]);
  }
}

FmmChargeField::FmmChargeField(const ClassicalSystem& system, const FmmFieldOptions& options,
                               std::optional<InnerSwitch> switching)
    : charges_(classical_charges(system)), options_(options), switching_(switching) {
  if (!(options.absolute_accuracy > 0.0) || !std::isfinite(options.absolute_accuracy)) {
    throw std::invalid_argument("fika::FmmChargeField: the accuracy must be positive");
  }
  if (!switching_) {
    return;
  }
  check_inner_switch(*switching_, "fika::FmmChargeField");
  const std::size_t nc = charges_.coordinates.size();
  if (nc == 0) {
    return;
  }
  charge_residue_.resize(nc);
  for (std::size_t r = 0; r + 1 < charges_.residue_offsets.size(); ++r) {
    for (std::size_t c = charges_.residue_offsets[r]; c < charges_.residue_offsets[r + 1]; ++c) {
      charge_residue_[c] = static_cast<std::uint32_t>(r);
    }
  }
}

void FmmChargeField::add_field(const PolarizableSites& sites,
                               std::span<Point3D<double>> field) const {
  const std::size_t n = sites.positions.size();
  if (field.size() != n) {
    throw std::invalid_argument("fika::FmmChargeField: field has " + std::to_string(field.size()) +
                                " entries for " + std::to_string(n) + " sites");
  }
  if (sites.owners.size() != n) {
    throw std::invalid_argument("fika::FmmChargeField: " + std::to_string(sites.owners.size()) +
                                " owners for " + std::to_string(n) + " sites");
  }
  const std::size_t residues = charges_.residue_offsets.size() - 1;
  for (const std::size_t owner : sites.owners) {
    if (owner >= residues) {
      throw std::invalid_argument("fika::FmmChargeField: site of residue " + std::to_string(owner) +
                                  " outside the classical system");
    }
  }
  double charge_scale = 0.0;
  for (const double q : charges_.charges) {
    charge_scale = std::max(charge_scale, std::abs(q));
  }
  std::vector<Point3D<double>> fmm_field(n);
  bool direct = false;
  if (charge_scale > 0.0 && n > 0) {
    const std::vector<Point3D<double>> no_dipole_positions;
    const std::vector<Dipole> no_dipoles;
    int order = options_.order;
    retries_ = 0;
    while (true) {
      detail::VolumeFmmOptions fmm_options{.absolute_accuracy = options_.absolute_accuracy,
                                           .order = order,
                                           .charge_scale = charge_scale};
      if (switching_) {
        // Keep every switched pair in the P2P lists: the near-field weight then sees all of
        // them, without cancellation between a huge bare term and its correction.
        fmm_options.minimum_near_radius = switching_->r_off;
      }
      const detail::VolumeFmm fmm(charges_.coordinates, no_dipole_positions, sites.positions,
                                  fmm_options);
      const bool weighted_near = switching_ && switching_->r_off <= fmm.near_radius();
      if (switching_ && !weighted_near) {
        // A system too small for the tree: sum the switched field directly instead.
        direct = true;
        break;
      }
      detail::VolumeFmmReport report;
      // Same-residue pairs keep the bare field in the near sums and are removed exactly below;
      // the weight applies S to the different-residue near-field pairs.
      if (switching_) {
        const detail::NearFieldWeight weight =
            [this, &sites](std::size_t target, std::size_t source, double r2) -> double {
          if (charge_residue_[source] == sites.owners[target]) {
            return 1.0;
          }
          return inner_switch(std::sqrt(r2), *switching_);
        };
        fmm.field(charges_.charges, no_dipoles, fmm_field, &report, &weight);
      } else {
        fmm.field(charges_.charges, no_dipoles, fmm_field, &report);
      }
      order_ = fmm.order();
      if (report.sampled_error <= 0.1 * options_.absolute_accuracy ||
          order_ + 2 > detail::largest_automatic_fmm_order) {
        break;
      }
      order = order_ + 2;
      ++retries_;
    }
  }
  const std::span<const double> q = charges_.charges;
  const std::span<const Point3D<double>> positions = charges_.coordinates;
  if (direct) {
    // The tree cannot put every switched pair in the near field: sum the switched field directly.
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t begin = charges_.residue_offsets[sites.owners[i]];
      const std::size_t end = charges_.residue_offsets[sites.owners[i] + 1];
      Point3D<double> e = field[i];
      add_charge_field(q, positions, 0, begin, sites.positions[i], switching_, e);
      add_charge_field(q, positions, end, q.size(), sites.positions[i], switching_, e);
      field[i] = e;
    }
    return;
  }
  // Charges of the site's own residue (other than one at the site itself) removed.
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t begin = charges_.residue_offsets[sites.owners[i]];
    const std::size_t end = charges_.residue_offsets[sites.owners[i] + 1];
    Point3D<double> e = fmm_field[i];
    const Point3D<double>& s = sites.positions[i];
    for (std::size_t c = begin; c < end; ++c) {
      const double dx = s.x - positions[c].x, dy = s.y - positions[c].y, dz = s.z - positions[c].z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 == 0.0) {
        continue;
      }
      const double scale = q[c] / (r2 * std::sqrt(r2));
      e = {e.x - scale * dx, e.y - scale * dy, e.z - scale * dz};
    }
    field[i] = {field[i].x + e.x, field[i].y + e.y, field[i].z + e.z};
  }
}

}  // namespace fika

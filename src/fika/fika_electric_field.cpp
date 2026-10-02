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

#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "fika_volume_fmm.hpp"

namespace fika {

namespace {

/// Adds sum_{c in [begin, end)} q_c (s - C) / |s - C|^3 to `field`.
void add_charge_field(std::span<const double> charges, std::span<const Point3D<double>> positions,
                      std::size_t begin, std::size_t end, const Point3D<double>& s,
                      Point3D<double>& field) {
  double fx = 0.0, fy = 0.0, fz = 0.0;
  for (std::size_t c = begin; c < end; ++c) {
    const double dx = s.x - positions[c].x;
    const double dy = s.y - positions[c].y;
    const double dz = s.z - positions[c].z;
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double scale = charges[c] / (r2 * std::sqrt(r2));
    fx += scale * dx;
    fy += scale * dy;
    fz += scale * dz;
  }
  field.x += fx;
  field.y += fy;
  field.z += fz;
}

}  // namespace

PermanentChargeField::PermanentChargeField(const ClassicalSystem& system)
    : charges_(classical_charges(system)) {}

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
    add_charge_field(q, positions, 0, begin, sites.positions[i], field[i]);
    add_charge_field(q, positions, end, q.size(), sites.positions[i], field[i]);
  }
}

FmmChargeField::FmmChargeField(const ClassicalSystem& system, const FmmFieldOptions& options)
    : charges_(classical_charges(system)), options_(options) {
  if (!(options.absolute_accuracy > 0.0) || !std::isfinite(options.absolute_accuracy)) {
    throw std::invalid_argument("fika::FmmChargeField: the accuracy must be positive");
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
  if (charge_scale > 0.0 && n > 0) {
    const std::vector<Point3D<double>> no_dipole_positions;
    const std::vector<Dipole> no_dipoles;
    int order = options_.order;
    retries_ = 0;
    while (true) {
      const detail::VolumeFmm fmm(
          charges_.coordinates, no_dipole_positions, sites.positions,
          detail::VolumeFmmOptions{.absolute_accuracy = options_.absolute_accuracy,
                                   .order = order,
                                   .charge_scale = charge_scale});
      detail::VolumeFmmReport report;
      fmm.field(charges_.charges, no_dipoles, fmm_field, &report);
      order_ = fmm.order();
      if (report.sampled_error <= 0.1 * options_.absolute_accuracy ||
          order_ + 2 > detail::largest_automatic_fmm_order) {
        break;
      }
      order = order_ + 2;
      ++retries_;
    }
  }
  // Charges of the site's own residue (other than one at the site itself) removed.
  const std::span<const double> q = charges_.charges;
  const std::span<const Point3D<double>> positions = charges_.coordinates;
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

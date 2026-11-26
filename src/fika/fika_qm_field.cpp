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

#include "fika_qm_field.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "fika_veloxchem_order.hpp"

namespace fika {

QmNuclearField::QmNuclearField(const Molecule<double>& molecule)
    : coordinates_(molecule.coordinates().begin(), molecule.coordinates().end()) {
  charges_.reserve(molecule.size());
  for (const Element& element : molecule.elements()) {
    charges_.push_back(static_cast<double>(element.atomic_number()));
  }
  for (const Point3D<double>& r : coordinates_) {
    if (!std::isfinite(r.x) || !std::isfinite(r.y) || !std::isfinite(r.z)) {
      throw std::invalid_argument("fika::QmNuclearField: nuclear coordinates must be finite");
    }
  }
}

void QmNuclearField::add_field(const PolarizableSites& sites,
                               std::span<Point3D<double>> field) const {
  const std::size_t n = sites.positions.size();
  if (field.size() != n) {
    throw std::invalid_argument("fika::QmNuclearField: field has " + std::to_string(field.size()) +
                                " entries for " + std::to_string(n) + " sites");
  }
  bool on_nucleus = false;
#pragma omp parallel for schedule(static) reduction(|| : on_nucleus)
  for (std::size_t i = 0; i < n; ++i) {
    const Point3D<double>& s = sites.positions[i];
    double fx = 0.0, fy = 0.0, fz = 0.0;
    for (std::size_t a = 0; a < charges_.size(); ++a) {
      const double dx = s.x - coordinates_[a].x;
      const double dy = s.y - coordinates_[a].y;
      const double dz = s.z - coordinates_[a].z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 == 0.0) {
        on_nucleus = true;
        continue;
      }
      const double scale = charges_[a] / (r2 * std::sqrt(r2));
      fx += scale * dx;
      fy += scale * dy;
      fz += scale * dz;
    }
    field[i].x += fx;
    field[i].y += fy;
    field[i].z += fz;
  }
  if (on_nucleus) {
    throw std::domain_error("fika::QmNuclearField: a polarizable site sits on a nucleus");
  }
}

namespace {

auto fika_density(const MolecularBasis& basis, const DenseMatrix& density) -> DenseMatrix {
  const DenseMatrix permuted = veloxchem_to_fika(density, basis);
  if (permuted.symmetry() == MatrixSymmetry::symmetric) {
    return permuted;
  }
  // Only the symmetric part (D + D^T) / 2 creates a field (the integrals are symmetric).
  const std::size_t n = permuted.rows();
  DenseMatrix symmetric(n, MatrixSymmetry::symmetric);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      symmetric.set(i, j, 0.5 * (permuted(i, j) + permuted(j, i)));
    }
  }
  return symmetric;
}

}  // namespace

QmElectronicField::QmElectronicField(const Molecule<double>& molecule, const MolecularBasis& basis,
                                     const DenseMatrix& density, double accuracy,
                                     ChargeSummation summation)
    : sources_(detail::density_multipoles(molecule, basis, fika_density(basis, density),
                                          0.5 * accuracy)),
      accuracy_(accuracy),
      summation_(summation) {}

void QmElectronicField::add_field(const PolarizableSites& sites,
                                  std::span<Point3D<double>> field) const {
  if (field.size() != sites.positions.size()) {
    throw std::invalid_argument("fika::QmElectronicField: field has " +
                                std::to_string(field.size()) + " entries for " +
                                std::to_string(sites.positions.size()) + " sites");
  }
  report_ = {};
  const bool multipole =
      summation_ == ChargeSummation::multipole ||
      (summation_ == ChargeSummation::automatic && sites.positions.size() >= multipole_site_count);
  if (multipole) {
    report_ = detail::add_density_field_fmm(sources_, sites.positions, field,
                                            {.accuracy = 0.5 * accuracy_});
  } else {
    detail::add_density_field(sources_, sites.positions, field);
  }
}

}  // namespace fika

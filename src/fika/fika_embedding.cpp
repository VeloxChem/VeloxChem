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

#include "fika_embedding.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_veloxchem_order.hpp"
#include "fika_compensated_sum.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_dipole_potential_driver.hpp"
#include "fika_nuclear_attraction_driver.hpp"
#include "fika_point_sources.hpp"
#include "fika_polarizable_sites.hpp"
#include "fika_qmmm_induced_dipoles.hpp"
#include "fika_qm_field.hpp"

namespace fika {

namespace {

/// sum_ij A_ij B_ij of a density (any symmetry) and a symmetric matrix, compensated.
auto trace_product(const DenseMatrix& a, const DenseMatrix& b) -> double {
  const std::size_t n = a.rows();
  detail::NeumaierSum sum;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      sum.add(a(i, j) * b(i, j));
    }
  }
  return sum.value();
}

}  // namespace

auto QmmmEmbedding::fock() const -> DenseMatrix {
  if (permanent_fock.rows() == 0) {
    return induced_fock;
  }
  DenseMatrix total = permanent_fock;
  const auto addend = induced_fock.values();
  auto values = total.values();
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += addend[i];
  }
  return total;
}

auto induced_dipole_fock(const Molecule<double>& molecule, const MolecularBasis& basis,
                         std::span<const Point3D<double>> positions,
                         std::span<const Point3D<double>> dipoles, double threshold,
                         ChargeSummation summation) -> DenseMatrix {
  if (positions.size() != dipoles.size()) {
    throw std::invalid_argument("fika::induced_dipole_fock: " + std::to_string(dipoles.size()) +
                                " dipoles at " + std::to_string(positions.size()) + " positions");
  }
  const auto finite = [](const Point3D<double>& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  const auto zero = [](const Point3D<double>& p) {
    return p.x == 0.0 && p.y == 0.0 && p.z == 0.0;
  };
  // All-zero dipoles (e.g. in response to a vanishing density) contribute nothing: skip the
  // integrals, after the checks the driver would make.
  if (std::ranges::all_of(dipoles, zero)) {
    if (!std::isfinite(threshold) || threshold < 0.0) {
      throw std::invalid_argument("fika::induced_dipole_fock: threshold must be finite and "
                                  "non-negative");
    }
    if (!std::ranges::all_of(positions, finite)) {
      throw std::invalid_argument("fika::induced_dipole_fock: positions must be finite");
    }
    return DenseMatrix(basis.function_count(), MatrixSymmetry::symmetric);
  }
  std::vector<Dipole> negated(dipoles.size());
  for (std::size_t s = 0; s < dipoles.size(); ++s) {
    negated[s] = Dipole{{-dipoles[s].x, -dipoles[s].y, -dipoles[s].z}};
  }
  return fika_to_veloxchem(DipolePotentialDriver(0, summation)
                               .compute(molecule, basis, negated, positions, threshold)
                               .to_dense_matrix(),
                           basis);
}

QmmmEmbeddingDriver::QmmmEmbeddingDriver(const Molecule<double>& molecule,
                                         const MolecularBasis& basis,
                                         const ClassicalSystem& system,
                                         std::optional<TholeDamping> damping,
                                         const QmmmEmbeddingOptions& options)
    : molecule_(molecule),
      basis_(basis),
      damping_(damping),
      options_(options),
      sites_(polarizable_sites(system)),
      system_(system) {
  options_.induced.initial_guess.clear();
  const InducedDipoleOptions& induced = options_.induced;
  if (induced.summation == ChargeSummation::direct ||
      (induced.summation == ChargeSummation::automatic &&
       sites_.positions.size() < multipole_coupling_sites)) {
    direct_.emplace(sites_, damping_);
  }
}

auto QmmmEmbeddingDriver::compute(const DenseMatrix& density, QmmmSources sources,
                                  std::span<const Point3D<double>> initial_guess)
    -> QmmmEmbedding {
  const bool ground = sources == QmmmSources::all;
  if (ground && !permanent_) {
    // Field in the order induced_dipoles sums it: MM permanent charges, then the QM nuclei.
    Permanent permanent;
    permanent.field.assign(sites_.positions.size(), Point3D<double>{0.0, 0.0, 0.0});
    if (options_.induced.permanent_field) {
      permanent.field_summation =
          detail::add_permanent_field(*system_, sites_, options_.induced, permanent.field);
    }
    QmNuclearField(molecule_).add_field(sites_, permanent.field);
    permanent.fock =
        fika_to_veloxchem(NuclearAttractionDriver(0, options_.fock_summation)
                              .compute(molecule_, basis_, *system_, options_.fock_threshold)
                              .to_dense_matrix(),
                          basis_);
    permanent.nuclear_energy = nuclear_classical_energy(molecule_, *system_);
    permanent_ = std::move(permanent);
    system_.reset();
  }

  InducedDipoleOptions induced = options_.induced;
  induced.initial_guess.assign(initial_guess.begin(), initial_guess.end());
  if (!ground) {
    induced.permanent_field = false;
  }
  std::vector<Point3D<double>> field =
      ground ? permanent_->field
             : std::vector<Point3D<double>>(sites_.positions.size(), Point3D<double>{});
  const QmElectronicField electrons(molecule_, basis_, density, field_accuracy(induced),
                                    induced.summation);
  electrons.add_field(sites_, field);

  QmmmEmbedding result;
  result.induced = detail::solve_with_selected_coupling(sites_, std::move(field), damping_, induced,
                                                        direct_ ? &*direct_ : nullptr);
  if (ground) {
    result.induced.field_summation = permanent_->field_summation.summation;
    result.induced.field_order = permanent_->field_summation.order;
  }
  result.electron_field_order = electrons.report().order;
  result.electron_field_summation =
      result.electron_field_order > 0 ? ChargeSummation::multipole : ChargeSummation::direct;
  if (!result.induced.converged && !options_.allow_unconverged) {
    throw std::runtime_error("fika::qmmm_embedding: induced dipoles not converged after " +
                             std::to_string(result.induced.iterations) +
                             " iterations (residual " + std::to_string(result.induced.residual) +
                             ")");
  }
  result.induced_fock =
      induced_dipole_fock(molecule_, basis_, sites_.positions, result.induced.dipoles,
                          options_.fock_threshold, options_.fock_summation);
  if (ground) {
    result.permanent_fock = permanent_->fock;
    result.electron_permanent_energy = trace_product(density, result.permanent_fock);
    result.nuclear_permanent_energy = permanent_->nuclear_energy;
    result.polarization_energy = induction_energy(result.induced.dipoles, result.induced.field);
  }
  return result;
}

auto qmmm_embedding(const Molecule<double>& molecule, const MolecularBasis& basis,
                    const DenseMatrix& density, const ClassicalSystem& system,
                    std::optional<TholeDamping> damping, const QmmmEmbeddingOptions& options)
    -> QmmmEmbedding {
  return QmmmEmbeddingDriver(molecule, basis, system, damping, options)
      .compute(density, options.sources, options.induced.initial_guess);
}

}  // namespace fika

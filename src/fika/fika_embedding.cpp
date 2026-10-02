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

#include <array>
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
  std::vector<Dipole> negated(dipoles.size());
  for (std::size_t s = 0; s < dipoles.size(); ++s) {
    negated[s] = Dipole{{-dipoles[s].x, -dipoles[s].y, -dipoles[s].z}};
  }
  return fika_to_veloxchem(DipolePotentialDriver(0, summation)
                               .compute(molecule, basis, negated, positions, threshold)
                               .to_dense_matrix(),
                           basis);
}

auto qmmm_embedding(const Molecule<double>& molecule, const MolecularBasis& basis,
                    const DenseMatrix& density, const ClassicalSystem& system,
                    std::optional<TholeDamping> damping, const QmmmEmbeddingOptions& options)
    -> QmmmEmbedding {
  QmmmEmbedding result;
  const auto sites = polarizable_sites(system);
  if (options.sources == QmmmSources::all) {
    auto induced = qmmm_induced_dipoles(molecule, basis, density, system, damping, options.induced);
    result.induced = std::move(induced.induced);
    result.electron_field_summation = induced.electron_field_summation;
    result.electron_field_order = induced.electron_field_order;
  } else {
    const InducedDipoleOptions& induced = options.induced;
    const double accuracy =
        induced.field_accuracy > 0.0 ? induced.field_accuracy : 0.1 * induced.tolerance;
    const QmElectronicField electrons(molecule, basis, density, accuracy, induced.summation);
    const std::array<const FieldContribution*, 1> external{&electrons};
    InducedDipoleOptions electrons_only = induced;
    electrons_only.permanent_field = false;
    result.induced = induced_dipoles(system, damping, electrons_only, external);
    result.electron_field_order = electrons.report().order;
    result.electron_field_summation =
        result.electron_field_order > 0 ? ChargeSummation::multipole : ChargeSummation::direct;
  }
  result.induced_fock =
      induced_dipole_fock(molecule, basis, sites.positions, result.induced.dipoles,
                          options.fock_threshold, options.fock_summation);
  if (options.sources == QmmmSources::all) {
    result.permanent_fock =
        fika_to_veloxchem(NuclearAttractionDriver(0, options.fock_summation)
                              .compute(molecule, basis, system, options.fock_threshold)
                              .to_dense_matrix(),
                          basis);
    result.electron_permanent_energy = trace_product(density, result.permanent_fock);
    result.nuclear_permanent_energy = nuclear_classical_energy(molecule, system);
    result.polarization_energy = induction_energy(result.induced.dipoles, result.induced.field);
  }
  return result;
}

}  // namespace fika

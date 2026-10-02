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

#include "fika_nuclear_attraction_driver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_nuclear_attraction_kernels.hpp"
#include "fika_overlap_screening.hpp"
#include "fika_two_centre_driver.hpp"
#include "fika_point_sources.hpp"
#include "fika_far_field_expansion.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::NuclearAttractionDriver: " + reason);
}

/// Kernel of a shell pair: (A|V|A) blocks with ordering A, (A|V|B) blocks with Scheme II.
class NuclearAttractionKernel final : public detail::ShellPairKernel {
 public:
  NuclearAttractionKernel(const BasisShell& bra, const BasisShell& ket,
                          const detail::SourceFields& fields,
                          const detail::FarFieldExpansion* far_field, double charge_sum,
                          double threshold)
      : same_atom_(detail::make_same_atom_shell_pair(bra, ket)),
        two_atom_(detail::make_two_atom_shell_pair(bra, ket, {.charges = charge_sum}, threshold)),
        fields_(&fields),
        far_field_(far_field) {}

  // The blocks are not symmetric in m, m' for l == l' (the charges break the symmetry).
  auto shape() const noexcept -> detail::KernelShape override {
    return {same_atom_.l, same_atom_.l_prime, same_atom_.bra_contractions,
            same_atom_.ket_contractions, false};
  }

  void compute(const SolidHarmonics& harmonics, std::span<const AtomPair> pairs, std::size_t n,
               detail::KernelWorkspace& workspace, std::span<double> values) const override {
    if (n > 0 && pairs[0].bra == pairs[0].ket) {
      detail::same_atom_potential_values(same_atom_, *fields_, pairs, n, workspace, values);
    } else {
      detail::two_atom_potential_values(two_atom_, *fields_, harmonics, pairs, n, workspace, values,
                                        far_field_);
    }
  }

 private:
  detail::SameAtomShellPair same_atom_;
  detail::TwoAtomShellPair two_atom_;
  const detail::SourceFields* fields_;
  const detail::FarFieldExpansion* far_field_;
};

class NuclearAttractionOperator final : public detail::TwoCentreOperator {
 public:
  NuclearAttractionOperator(const Molecule<double>& molecule, const MolecularBasis& basis,
                            std::span<const double> charges,
                            std::span<const Point3D<double>> coordinates, double threshold,
                            bool multipole)
      : charge_sum_(absolute_sum(charges)),
        threshold_(threshold),
        far_field_(
            multipole
                ? std::make_unique<detail::FarFieldExpansion>(detail::make_source_far_field(
                      molecule, basis,
                      detail::PointSources{.charges = charges, .charge_coordinates = coordinates},
                      threshold))
                : nullptr),
        fields_(molecule, basis,
                detail::PointSources{.charges = charges, .charge_coordinates = coordinates},
                far_field_.get()) {}

  auto bound(const BasisShell& bra, const BasisShell& ket) const
      -> detail::ShellPairBound override {
    return detail::nuclear_attraction_shell_pair_bound(bra, ket, charge_sum_);
  }

  auto kernel(const BasisShell& bra, const BasisShell& ket) const
      -> std::unique_ptr<detail::ShellPairKernel> override {
    return std::make_unique<NuclearAttractionKernel>(bra, ket, fields_, far_field_.get(),
                                                     charge_sum_, threshold_);
  }

  auto same_centre_blocks_diagonal() const -> bool override { return false; }

  auto same_centre(const BasisShell&, const BasisShell&) const -> std::vector<double> override {
    throw std::logic_error("fika: nuclear-attraction same-atom blocks come from the kernels");
  }

 private:
  static auto absolute_sum(std::span<const double> charges) -> double {
    double sum = 0.0;
    for (const double charge : charges) {
      sum += std::abs(charge);
    }
    return sum;
  }

  // Built in this order: the charge fields use the far-field expansion.
  double charge_sum_;
  double threshold_;
  std::unique_ptr<detail::FarFieldExpansion> far_field_;
  detail::SourceFields fields_;
};

}  // namespace

NuclearAttractionDriver::NuclearAttractionDriver(std::size_t block_size, ChargeSummation summation)
    : block_size_(block_size), summation_(summation) {}

auto NuclearAttractionDriver::compute(const Molecule<double>& molecule, const MolecularBasis& basis,
                                      std::span<const double> charges,
                                      std::span<const Point3D<double>> charge_coordinates,
                                      double threshold) const -> BlockSparseMatrix {
  if (!std::isfinite(threshold) || threshold < 0.0) {
    fail("threshold must be finite and non-negative");
  }
  if (molecule.size() != basis.atom_count()) {
    fail("molecule has " + std::to_string(molecule.size()) + " atoms, basis " +
         std::to_string(basis.atom_count()));
  }
  if (charges.size() != charge_coordinates.size()) {
    fail(std::to_string(charges.size()) + " charges but " +
         std::to_string(charge_coordinates.size()) + " charge coordinates");
  }
  const auto finite = [](const Point3D<double>& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  if (!std::ranges::all_of(charges, [](double q) { return std::isfinite(q); }) ||
      !std::ranges::all_of(charge_coordinates, finite)) {
    fail("charges and their coordinates must be finite");
  }
  if (summation_ == ChargeSummation::multipole && threshold == 0.0) {
    fail("the multipole charge summation needs a positive threshold");
  }
  const bool multipole = summation_ == ChargeSummation::multipole ||
                         (summation_ == ChargeSummation::automatic && threshold > 0.0 &&
                          charges.size() >= multipole_charge_count);
  return detail::compute_two_centre(
      NuclearAttractionOperator(molecule, basis, charges, charge_coordinates, threshold, multipole),
      molecule, basis, basis, true, threshold, block_size_, PairCost::heavy);
}

auto NuclearAttractionDriver::compute(const Molecule<double>& molecule, const MolecularBasis& basis,
                                      double threshold) const -> BlockSparseMatrix {
  std::vector<double> charges;
  charges.reserve(molecule.size());
  for (const Element& element : molecule.elements()) {
    charges.push_back(-static_cast<double>(element.atomic_number()));
  }
  return compute(molecule, basis, charges, molecule.coordinates(), threshold);
}

auto NuclearAttractionDriver::compute(const Molecule<double>& molecule, const MolecularBasis& basis,
                                      const ClassicalSystem& system, double threshold) const
    -> BlockSparseMatrix {
  PointCharges sources = classical_charges(system);
  for (double& charge : sources.charges) {
    charge = -charge;
  }
  return compute(molecule, basis, sources.charges, sources.coordinates, threshold);
}

}  // namespace fika

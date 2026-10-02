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

#include "fika_dipole_potential_driver.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_nuclear_attraction_kernels.hpp"
#include "fika_overlap_screening.hpp"
#include "fika_two_centre_driver.hpp"
#include "fika_far_field_expansion.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::DipolePotentialDriver: " + reason);
}

/// Kernel of a shell pair: (A|D|A) blocks with ordering A, (A|D|B) blocks with Scheme II (the
/// point-source kernels with dipoles only).
class DipolePotentialKernel final : public detail::ShellPairKernel {
 public:
  DipolePotentialKernel(const BasisShell& bra, const BasisShell& ket,
                        const detail::SourceFields& fields,
                        const detail::FarFieldExpansion* far_field, double dipole_sum,
                        double threshold)
      : same_atom_(detail::make_same_atom_shell_pair(bra, ket)),
        two_atom_(detail::make_two_atom_shell_pair(bra, ket, {.dipoles = dipole_sum}, threshold)),
        fields_(&fields),
        far_field_(far_field) {}

  // The blocks are not symmetric in m, m' for l == l' (the dipoles break the symmetry).
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

class DipolePotentialOperator final : public detail::TwoCentreOperator {
 public:
  DipolePotentialOperator(const Molecule<double>& molecule, const MolecularBasis& basis,
                          std::span<const Dipole> dipoles,
                          std::span<const Point3D<double>> coordinates, double threshold,
                          bool multipole)
      : threshold_(threshold),
        dipole_sum_(absolute_sum(dipoles)),
        far_field_(
            multipole
                ? std::make_unique<detail::FarFieldExpansion>(detail::make_source_far_field(
                      molecule, basis,
                      detail::PointSources{.dipoles = dipoles, .dipole_coordinates = coordinates},
                      threshold))
                : nullptr),
        fields_(molecule, basis,
                detail::PointSources{.dipoles = dipoles, .dipole_coordinates = coordinates},
                far_field_.get()) {}

  auto bound(const BasisShell& bra, const BasisShell& ket) const
      -> detail::ShellPairBound override {
    return detail::dipole_potential_shell_pair_bound(bra, ket, dipole_sum_);
  }

  auto kernel(const BasisShell& bra, const BasisShell& ket) const
      -> std::unique_ptr<detail::ShellPairKernel> override {
    return std::make_unique<DipolePotentialKernel>(bra, ket, fields_, far_field_.get(), dipole_sum_,
                                                   threshold_);
  }

  auto same_centre_blocks_diagonal() const -> bool override { return false; }

  auto same_centre(const BasisShell&, const BasisShell&) const -> std::vector<double> override {
    throw std::logic_error("fika: dipole-potential same-atom blocks come from the kernels");
  }

 private:
  static auto absolute_sum(std::span<const Dipole> dipoles) -> double {
    double sum = 0.0;
    for (const Dipole& dipole : dipoles) {
      sum += std::hypot(dipole(0), dipole(1), dipole(2));
    }
    return sum;
  }

  // Built in this order: the dipole fields use the far-field expansion.
  double threshold_;
  double dipole_sum_;
  std::unique_ptr<detail::FarFieldExpansion> far_field_;
  detail::SourceFields fields_;
};

}  // namespace

DipolePotentialDriver::DipolePotentialDriver(std::size_t block_size, ChargeSummation summation)
    : block_size_(block_size), summation_(summation) {}

auto DipolePotentialDriver::compute(const Molecule<double>& molecule, const MolecularBasis& basis,
                                    std::span<const Dipole> dipoles,
                                    std::span<const Point3D<double>> dipole_coordinates,
                                    double threshold) const -> BlockSparseMatrix {
  if (!std::isfinite(threshold) || threshold < 0.0) {
    fail("threshold must be finite and non-negative");
  }
  if (molecule.size() != basis.atom_count()) {
    fail("molecule has " + std::to_string(molecule.size()) + " atoms, basis " +
         std::to_string(basis.atom_count()));
  }
  if (dipoles.size() != dipole_coordinates.size()) {
    fail(std::to_string(dipoles.size()) + " dipoles but " +
         std::to_string(dipole_coordinates.size()) + " dipole coordinates");
  }
  const auto finite_point = [](const Point3D<double>& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  const auto finite_dipole = [](const Dipole& d) {
    return std::ranges::all_of(d.components, [](double v) { return std::isfinite(v); });
  };
  if (!std::ranges::all_of(dipoles, finite_dipole) ||
      !std::ranges::all_of(dipole_coordinates, finite_point)) {
    fail("dipoles and their coordinates must be finite");
  }
  if (summation_ == ChargeSummation::multipole && threshold == 0.0) {
    fail("the multipole summation needs a positive threshold");
  }
  const bool multipole = summation_ == ChargeSummation::multipole ||
                         (summation_ == ChargeSummation::automatic && threshold > 0.0 &&
                          dipoles.size() >= multipole_dipole_count);
  return detail::compute_two_centre(
      DipolePotentialOperator(molecule, basis, dipoles, dipole_coordinates, threshold, multipole),
      molecule, basis, basis, true, threshold, block_size_, PairCost::heavy);
}

}  // namespace fika

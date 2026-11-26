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

#include "fika_point_source_operator.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_nuclear_attraction_kernels.hpp"
#include "fika_overlap_screening.hpp"
#include "fika_two_centre_driver.hpp"

namespace fika::detail {

namespace {

/// Kernel of a shell pair: (A|V|A) blocks with ordering A, (A|V|B) blocks with Scheme II.
class PointSourceKernel final : public ShellPairKernel {
 public:
  PointSourceKernel(const BasisShell& bra, const BasisShell& ket, const SourceFields& fields,
                    const FarFieldExpansion* far_field, const SourceSums& sums, double threshold)
      : same_atom_(make_same_atom_shell_pair(bra, ket)),
        two_atom_(make_two_atom_shell_pair(bra, ket, sums, threshold)),
        fields_(&fields),
        far_field_(far_field) {}

  // The blocks are not symmetric in m, m' for l == l' (the sources break the symmetry).
  auto shape() const noexcept -> KernelShape override {
    return {same_atom_.l, same_atom_.l_prime, same_atom_.bra_contractions,
            same_atom_.ket_contractions, false};
  }

  void compute(const SolidHarmonics& harmonics, std::span<const AtomPair> pairs, std::size_t n,
               KernelWorkspace& workspace, std::span<double> values) const override {
    if (n > 0 && pairs[0].bra == pairs[0].ket) {
      same_atom_potential_values(same_atom_, *fields_, pairs, n, workspace, values);
    } else {
      two_atom_potential_values(two_atom_, *fields_, harmonics, pairs, n, workspace, values,
                                far_field_);
    }
  }

 private:
  SameAtomShellPair same_atom_;
  TwoAtomShellPair two_atom_;
  const SourceFields* fields_;
  const FarFieldExpansion* far_field_;
};

class PointSourceOperator final : public TwoCentreOperator {
 public:
  PointSourceOperator(const Molecule<double>& molecule, const MolecularBasis& basis,
                      const PointSources& sources, PointSourceKind kind, double threshold,
                      bool multipole)
      : kind_(kind),
        sums_(norm_sums(sources)),
        threshold_(threshold),
        far_field_(multipole ? std::make_unique<FarFieldExpansion>(
                                   make_source_far_field(molecule, basis, sources, threshold))
                             : nullptr),
        fields_(molecule, basis, sources, far_field_.get()) {}

  auto bound(const BasisShell& bra, const BasisShell& ket) const -> ShellPairBound override {
    switch (kind_) {
      case PointSourceKind::charges:
        return nuclear_attraction_shell_pair_bound(bra, ket, sums_.charges);
      case PointSourceKind::dipoles:
        return dipole_potential_shell_pair_bound(bra, ket, sums_.dipoles);
      case PointSourceKind::quadrupoles:
        return quadrupole_potential_shell_pair_bound(bra, ket, sums_.quadrupoles);
    }
    throw std::logic_error("fika: unknown point-source kind");
  }

  auto kernel(const BasisShell& bra, const BasisShell& ket) const
      -> std::unique_ptr<ShellPairKernel> override {
    return std::make_unique<PointSourceKernel>(bra, ket, fields_, far_field_.get(), sums_,
                                               threshold_);
  }

  auto same_centre_blocks_diagonal() const -> bool override { return false; }

  auto same_centre(const BasisShell&, const BasisShell&) const -> std::vector<double> override {
    throw std::logic_error("fika: point-source same-atom blocks come from the kernels");
  }

 private:
  /// sum |q|, sum |mu| and sum ||Theta|| (traceless norm) of the sources.
  static auto norm_sums(const PointSources& sources) -> SourceSums {
    SourceSums sums;
    for (const double charge : sources.charges) {
      sums.charges += std::abs(charge);
    }
    for (const Dipole& dipole : sources.dipoles) {
      sums.dipoles += std::hypot(dipole(0), dipole(1), dipole(2));
    }
    for (const Quadrupole& quadrupole : sources.quadrupoles) {
      sums.quadrupoles += traceless_quadrupole_norm(quadrupole);
    }
    return sums;
  }

  // Built in this order: the source fields use the far-field expansion.
  PointSourceKind kind_;
  SourceSums sums_;
  double threshold_;
  std::unique_ptr<FarFieldExpansion> far_field_;
  SourceFields fields_;
};

auto finite(const Point3D<double>& p) -> bool {
  return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

template <typename Source>
auto finite(const Source& source) -> bool {
  return std::ranges::all_of(source.components, [](double v) { return std::isfinite(v); });
}

auto finite(double charge) -> bool {
  return std::isfinite(charge);
}

/// Checks one kind's sources and coordinates; returns their count.
template <typename Source>
auto checked_count(std::span<const Source> values, std::span<const Point3D<double>> coordinates,
                   const std::string& noun, const auto& fail) -> std::size_t {
  if (values.size() != coordinates.size()) {
    fail(std::to_string(values.size()) + " " + noun + "s but " +
         std::to_string(coordinates.size()) + " " + noun + " coordinates");
  }
  if (!std::ranges::all_of(values, [](const Source& v) { return finite(v); }) ||
      !std::ranges::all_of(coordinates, [](const Point3D<double>& p) { return finite(p); })) {
    fail(noun + "s and their coordinates must be finite");
  }
  return values.size();
}

}  // namespace

auto point_source_matrix(const Molecule<double>& molecule, const MolecularBasis& basis,
                         const PointSources& sources, PointSourceKind kind, double threshold,
                         const PointSourceSettings& settings) -> BlockSparseMatrix {
  const auto fail = [&](const std::string& reason) {
    throw std::invalid_argument(std::string(settings.caller) + ": " + reason);
  };
  check_two_centre_input(settings.caller, molecule, basis, threshold);
  std::size_t count = 0;
  switch (kind) {
    case PointSourceKind::charges:
      count = checked_count(sources.charges, sources.charge_coordinates, "charge", fail);
      break;
    case PointSourceKind::dipoles:
      count = checked_count(sources.dipoles, sources.dipole_coordinates, "dipole", fail);
      break;
    case PointSourceKind::quadrupoles:
      count =
          checked_count(sources.quadrupoles, sources.quadrupole_coordinates, "quadrupole", fail);
      break;
  }
  if (settings.summation == ChargeSummation::multipole && threshold == 0.0) {
    fail("the multipole summation needs a positive threshold");
  }
  const bool multipole = settings.summation == ChargeSummation::multipole ||
                         (settings.summation == ChargeSummation::automatic && threshold > 0.0 &&
                          count >= settings.multipole_count);
  return compute_two_centre(
      PointSourceOperator(molecule, basis, sources, kind, threshold, multipole), molecule, basis,
      basis, true, threshold, settings.block_size, PairCost::heavy);
}

}  // namespace fika::detail

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

#ifndef fika_electric_field_hpp
#define fika_electric_field_hpp

#include <span>

#include "fika_point3d.hpp"
#include "fika_classical_system.hpp"
#include "fika_point_sources.hpp"
#include "fika_polarizable_sites.hpp"

namespace fika {

/// A source of the electric field at the polarizable sites (the right-hand side of the induced
/// dipole equations); the total field is the sum of the contributions.
class FieldContribution {
 public:
  virtual ~FieldContribution() = default;

  /// Adds this source's electric field (a.u.) at every site to `field` (one entry per site).
  virtual void add_field(const PolarizableSites& sites, std::span<Point3D<double>> field) const = 0;
};

/// Field of the permanent charges of a classical system (both regions), E(s) = sum_C q_C (s - C)
/// / |s - C|^3, without the charges of each site's own residue. Summed directly per site (in a
/// fixed order, so independent of the thread count).
class PermanentChargeField final : public FieldContribution {
 public:
  /// Throws as classical_charges.
  explicit PermanentChargeField(const ClassicalSystem& system);

  /// The sites must come from the same classical system (their owners index its residues).
  void add_field(const PolarizableSites& sites, std::span<Point3D<double>> field) const override;

 private:
  PointCharges charges_;
};

/// Options of the fast multipole permanent field.
struct FmmFieldOptions {
  double absolute_accuracy = 1e-9;  // field error target (a.u.)
  int order = 0;                    // starting expansion order; 0: from the accuracy
};

/// The field of PermanentChargeField through a fast multipole method (detail::VolumeFmm,
/// separation 2): the FMM sums every charge except one at the site itself, and the other charges
/// of the site's own residue are then subtracted directly. The evaluation checks the error on a
/// sample of sites; above accuracy / 10 the order is raised by 2 (up to 26) and the field
/// recomputed. add_field is therefore not safe to call concurrently.
class FmmChargeField final : public FieldContribution {
 public:
  /// Throws as classical_charges, or std::invalid_argument for a nonpositive accuracy.
  explicit FmmChargeField(const ClassicalSystem& system, const FmmFieldOptions& options = {});

  void add_field(const PolarizableSites& sites, std::span<Point3D<double>> field) const override;

  /// Expansion order of the last evaluation and its order increases.
  auto order() const noexcept -> int { return order_; }
  auto retries() const noexcept -> int { return retries_; }

 private:
  PointCharges charges_;
  FmmFieldOptions options_;
  mutable int order_ = 0;
  mutable int retries_ = 0;
};

}  // namespace fika

#endif  // fika_electric_field_hpp

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

#ifndef fika_dipole_interaction_hpp
#define fika_dipole_interaction_hpp

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_polarizable_sites.hpp"
#include "fika_volume_fmm.hpp"

namespace fika {

/// Thole exponential damping of the dipole-dipole interaction (van Duijnen and Swart, J. Phys.
/// Chem. A 102, 2399 (1998)): with v = a r / (alpha_i alpha_j)^(1/6) (mean polarizabilities),
///   lambda_3 = 1 - (1 + v + v^2/2) e^-v,  lambda_5 = 1 - (1 + v + v^2/2 + v^3/6) e^-v.
struct TholeDamping {
  double a = 2.1304;
};

/// The coupling of induced dipoles: y = T mu, the field at every site of the dipoles at the other
/// sites (sites of the same residue excluded).
class DipoleInteraction {
 public:
  virtual ~DipoleInteraction() = default;

  /// Sets y (one entry per site) to T mu.
  virtual void apply(std::span<const Point3D<double>> mu, std::span<Point3D<double>> y) const = 0;
};

/// T mu summed directly over all pairs (O(N^2)), with T_ij mu = lambda_5 3 r (r . mu) / r^5 -
/// lambda_3 mu / r^3, r = s_i - s_j, lambda = 1 without damping. The sum over pairs is undamped
/// and evaluates each pair once for both sites (T_ij = T_ji); the pairs where damping changes
/// lambda by more than 1e-15 (found with a cell grid) are then corrected, and pairs of sites of
/// the same residue removed. Independent of the thread count.
namespace detail {

/// Corrections of an undamped sum of T mu over all pairs of different sites: pairs of sites of the
/// same residue removed, and for close pairs of different residues the Thole damping added as
/// lambda - 1 terms (where it changes lambda by more than 1e-15, found with a cell grid). Shared
/// by the direct and fast multipole dipole couplings.
class DipoleCorrections {
 public:
  /// Throws std::invalid_argument (prefixed with `caller`) for inconsistent site arrays, a
  /// nonpositive damping parameter, or a site whose mean polarizability is not positive (with
  /// damping).
  DipoleCorrections(const PolarizableSites& sites, std::optional<TholeDamping> damping,
                    const char* caller);

  /// y[i] = undamped[i] + the corrections at site i for dipoles mu (structure of arrays);
  /// parallel over sites, independent of the thread count.
  void apply(const double* mu_x, const double* mu_y, const double* mu_z, const double* undamped_x,
             const double* undamped_y, const double* undamped_z,
             std::span<Point3D<double>> y) const;

  auto size() const noexcept -> std::size_t { return x_.size(); }
  auto x() const noexcept -> const double* { return x_.data(); }
  auto y() const noexcept -> const double* { return y_.data(); }
  auto z() const noexcept -> const double* { return z_.data(); }

 private:
  std::vector<double> x_, y_, z_;
  std::vector<double> scales_;  // mean polarizability^(1/6), with damping
  std::optional<TholeDamping> damping_;
  double cutoff_ = 0.0;  // damping negligible beyond (bohr)
  // Sites grouped by residue: the sites of site i's residue are
  // owner_sites_[owner_offsets_[g] .. owner_offsets_[g + 1]), g = owner_group_[i].
  std::vector<std::size_t> owner_offsets_;
  std::vector<std::uint32_t> owner_sites_;
  std::vector<std::size_t> owner_group_;
  // Cell grid of edge cutoff_ (with damping): sites of cell c are
  // cell_sites_[cell_offsets_[c] .. cell_offsets_[c + 1]).
  Point3D<double> origin_{};
  std::size_t cells_[3] = {0, 0, 0};
  std::vector<std::size_t> cell_offsets_;
  std::vector<std::uint32_t> cell_sites_;
};

}  // namespace detail

class DirectDipoleInteraction final : public DipoleInteraction {
 public:
  /// Throws std::invalid_argument for a nonpositive damping parameter or a site whose mean
  /// polarizability is not positive (with damping).
  DirectDipoleInteraction(const PolarizableSites& sites, std::optional<TholeDamping> damping);

  void apply(std::span<const Point3D<double>> mu, std::span<Point3D<double>> y) const override;

 private:
  detail::DipoleCorrections corrections_;
};

/// Options of the fast multipole dipole coupling.
struct FmmDipoleOptions {
  double absolute_accuracy = 1e-9;  // field error target (a.u.) for dipoles of dipole_scale
  double dipole_scale = 0.0;        // largest |mu| expected (e.g. max alpha |F|); required
  int order = 0;                    // starting expansion order; 0: from the accuracy
};

/// T mu through a fast multipole method (detail::VolumeFmm, separation 2) for the undamped sum
/// over all pairs, then the same corrections as DirectDipoleInteraction (same-residue pairs,
/// Thole damping). The FMM error is linear in mu, so products of smaller vectors (CG search
/// directions) are proportionally more accurate. The first apply checks the error on a sample of
/// sites; above accuracy / 10 the order is raised by 2 (up to 26) and the product recomputed.
/// apply is therefore not safe to call concurrently.
class FmmDipoleInteraction final : public DipoleInteraction {
 public:
  /// Throws std::invalid_argument as DirectDipoleInteraction, or for a nonpositive dipole scale
  /// or accuracy.
  FmmDipoleInteraction(const PolarizableSites& sites, std::optional<TholeDamping> damping,
                       const FmmDipoleOptions& options);

  void apply(std::span<const Point3D<double>> mu, std::span<Point3D<double>> y) const override;

  /// Expansion order in use, and the number of order increases after the sampled check.
  auto order() const noexcept -> int { return fmm_->order(); }
  auto retries() const noexcept -> int { return retries_; }

 private:
  void build(int order) const;

  detail::DipoleCorrections corrections_;
  std::vector<Point3D<double>> positions_;
  FmmDipoleOptions options_;
  mutable std::unique_ptr<detail::VolumeFmm> fmm_;
  mutable bool checked_ = false;
  mutable int retries_ = 0;
};

}  // namespace fika

#endif  // fika_dipole_interaction_hpp

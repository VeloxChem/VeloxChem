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

#ifndef fika_qm_field_hpp
#define fika_qm_field_hpp

#include <span>
#include <vector>

#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_summation.hpp"
#include "fika_dense_matrix.hpp"
#include "fika_electric_field.hpp"
#include "fika_polarizable_sites.hpp"
#include "fika_molecule.hpp"
#include "fika_density_fmm.hpp"
#include "fika_density_multipoles.hpp"

namespace fika {

/// Field of the nuclei of the QM region at the polarizable sites,
/// E(s) = sum_A Z_A (s - R_A) / |s - R_A|^3, summed directly per site over the atoms in order
/// (independent of the thread count).
class QmNuclearField final : public FieldContribution {
 public:
  /// Throws std::invalid_argument if a coordinate of `molecule` is not finite.
  explicit QmNuclearField(const Molecule<double>& molecule);

  /// Throws std::invalid_argument if `field` has not one entry per site, and std::domain_error if
  /// a site sits on a nucleus.
  void add_field(const PolarizableSites& sites, std::span<Point3D<double>> field) const override;

 private:
  std::vector<double> charges_;
  std::vector<Point3D<double>> coordinates_;
};

/// Field of the electrons of the QM region at the polarizable sites,
///   E(s) = -grad_s phi(s),  phi(s) = -sum_ab D_ab <a| 1 / |r - s| |b>,
/// for the total density D (alpha + beta): equivalently E_k(s) = Tr(D V(e_k at s)) with V the
/// dipole-potential integrals. The density is held as one source per primitive pair
/// (detail::DensityMultipoles): sites outside a source's penetration radius see its point
/// multipoles, sites inside it the exact near-field kernels. The sources are summed directly per
/// site, or through a fast multipole method (detail::add_density_field_fmm, with the order from
/// its error bound); half of the accuracy goes to dropping primitive pairs and, with the multipole
/// summation, half to the expansions. Either way each site sums in a fixed order (independent of
/// the thread count).
class QmElectronicField final : public FieldContribution {
 public:
  /// Site count from which ChargeSummation::automatic uses the multipole summation (measured
  /// break-even near 9300 sites: osimertinib/def2-SVP, accuracy 1e-9, the water oxygens nearest
  /// the molecule, 14 threads; the FMM takes 1.28x the direct time at 8000, 0.85x at 10000 and
  /// 0.57-0.66x at 14000-20000).
  static constexpr std::size_t multipole_site_count = 10000;

  /// `density`: over the functions of `basis` in VeloxChem's order (veloxchem_order), of any
  /// symmetry; only its symmetric part (D + D^T) / 2 creates a field (e.g. of a perturbed
  /// density), so an antisymmetric density gives none. The field error stays below `accuracy`
  /// (a.u., Euclidean norm per site). Throws std::invalid_argument for a molecule and basis of
  /// different atom counts, a density of the wrong size, or an accuracy that is not positive and
  /// finite.
  QmElectronicField(const Molecule<double>& molecule, const MolecularBasis& basis,
                    const DenseMatrix& density, double accuracy = 1e-10,
                    ChargeSummation summation = ChargeSummation::automatic);

  /// Throws std::invalid_argument if `field` has not one entry per site, and std::runtime_error if
  /// the multipole summation cannot reach the accuracy (detail::add_density_field_fmm).
  void add_field(const PolarizableSites& sites, std::span<Point3D<double>> field) const override;

  /// Number of primitive-pair sources kept after screening.
  auto source_count() const noexcept -> std::size_t { return sources_.size(); }

  /// Report of the last multipole summation (default when it summed directly); add_field is
  /// therefore not safe to call concurrently.
  auto report() const noexcept -> const detail::DensityFmmReport& { return report_; }

 private:
  detail::DensityMultipoles sources_;
  double accuracy_;
  ChargeSummation summation_;
  mutable detail::DensityFmmReport report_;
};

}  // namespace fika

#endif  // fika_qm_field_hpp

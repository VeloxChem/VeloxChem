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

#ifndef fika_density_multipoles_hpp
#define fika_density_multipoles_hpp

// Internal: the electron density of a QM region as point multipoles, one per primitive pair.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_dense_matrix.hpp"
#include "fika_solid_harmonics.hpp"
#include "fika_molecule.hpp"

namespace fika::detail {

/// The density rho(r) = sum_ab D_ab chi_a(r) chi_b(r) (D symmetric, fika's function order) as one
/// source per primitive pair (shells l, l', exponent p, product centre P): the field-row weights
///   W_row = sum D_(Im),(Jm') c_aI d_bJ W_(m m'),row      (primitive_pair_field_weights)
/// summed over the components and contracted functions of both shells (and over the mirrored
/// block of a pair of different shells), and from them the multipole moments
///   Q_LM = sum_kappa W_(L kappa),M Gamma(kappa + L + 3/2) / (2L + 1) p^-(kappa + L + 3/2)
///        = integral rho_source S_LM(r - P) d^3r,  L <= l + l'.
/// Outside its penetration radius sqrt(T_far(l + l' + 1) / p) a source's field equals that of its
/// point multipoles to rounding; inside, the weights give it exactly through the near-dipole
/// kernels of the dipole-potential integrals.
///
/// Primitive pairs whose contribution to the field cannot exceed their share of `accuracy` are
/// dropped: a shell-pair block with entry sum ||D_block||_1 gets the potential-element threshold
/// accuracy / (blocks ||D_block||_1) (doubled count for mirrored blocks) of the dipole-potential
/// screening (unit dipole), so the dropped pairs change no field component by more than
/// `accuracy` in total, anywhere.
struct DensityMultipoles {
  std::vector<Point3D<double>> centres;
  std::vector<double> exponents;            // p = alpha + beta
  std::vector<double> penetration_squared;  // squared penetration radius
  std::vector<int> bra_l;
  std::vector<int> ket_l;
  std::vector<std::size_t> weight_offsets;  // sources + 1, into weights
  std::vector<double> weights;              // per source field_row_count(l, l') rows
  std::vector<std::size_t> offsets;         // sources + 1, into moments
  std::vector<double> moments;              // per source (l + l' + 1)^2 rows L^2 + M + L
  std::size_t primitive_pairs = 0;          // primitive pairs considered (before screening)

  auto size() const noexcept -> std::size_t { return centres.size(); }

  auto rank(std::size_t source) const noexcept -> int { return bra_l[source] + ket_l[source]; }

  auto moments_of(std::size_t source) const -> std::span<const double> {
    return std::span(moments).subspan(offsets[source], offsets[source + 1] - offsets[source]);
  }

  auto weights_of(std::size_t source) const -> std::span<const double> {
    return std::span(weights).subspan(weight_offsets[source],
                                      weight_offsets[source + 1] - weight_offsets[source]);
  }
};

/// Throws std::invalid_argument if the molecule and basis have different atom counts, the density
/// is not a symmetric function_count() x function_count() matrix, or the accuracy is not positive
/// and finite.
auto density_multipoles(const Molecule<double>& molecule, const MolecularBasis& basis,
                        const DenseMatrix& density, double accuracy) -> DensityMultipoles;

/// Scratch storage of add_density_field_block.
struct DensityFieldWorkspace {
  std::vector<Point3D<double>> separations;
  std::vector<double> boys;
  SolidHarmonics harmonics;
};

/// Adds the exact field of the sources `which` (indices, summed in this order) at `sites` to `sums`
/// (3 per site, Racah components m = -1, 0, 1: y, z, x).
void add_density_field_block(const DensityMultipoles& sources, std::span<const std::uint32_t> which,
                             std::span<const Point3D<double>> sites, std::span<double> sums,
                             DensityFieldWorkspace& workspace);

/// Adds the exact field of the electrons (charge -1) of `sources` at `sites`: far sources through
/// their multipoles (as add_far_density_field), near ones through their field-row weights and the
/// near-dipole kernels. Sources are summed directly, in order for each site (in parallel over
/// blocks of sites, so independent of the thread count). Throws std::invalid_argument if `field`
/// has not one entry per site.
void add_density_field(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                       std::span<Point3D<double>> field);

/// Adds the field of the electrons (charge -1) of `sources` at `sites`,
///   E(C) = grad_C sum_LM Q_LM S_LM(C - P) / |C - P|^(2L + 1),
/// summed directly over the sources in order. Throws std::domain_error if a site is within a
/// source's penetration radius and std::invalid_argument if `field` has not one entry per site.
void add_far_density_field(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                           std::span<Point3D<double>> field);

}  // namespace fika::detail

#endif  // fika_density_multipoles_hpp

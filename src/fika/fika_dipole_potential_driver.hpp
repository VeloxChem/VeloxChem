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

#ifndef fika_dipole_potential_driver_hpp
#define fika_dipole_potential_driver_hpp

#include <cstddef>
#include <span>

#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_nuclear_attraction_driver.hpp"
#include "fika_block_sparse_matrix.hpp"
#include "fika_molecule.hpp"

namespace fika {

/// Two-centre integrals of the electrostatic potential of permanent point dipoles (of an MM
/// region):
///   V_ab = sum_C <a| mu_C . (r - C) / |r - C|^3 |b>
/// over contracted real solid-harmonic Gaussians, for dipoles mu_C at C (the same sign
/// convention as NuclearAttractionDriver's sum_C q_C <a|1/|r - C||b>: the potential of the
/// sources; pass -mu for the energy of an electron).
///
/// The dipoles break rotational symmetry, so blocks of equal atoms hold every shell pair of the
/// atom (DiagonalFormat::full) and are always stored; the threshold screens only blocks of
/// different atoms, with a bound independent of the dipole positions (scaled by
/// sum_C |mu_C|).
///
/// Blocks of equal atoms use ordering A (dipole notes, Section 3), blocks of different atoms
/// Scheme II (Section 6); far dipoles enter through their field tensor. With
/// ChargeSummation::multipole the dipoles beyond the penetration radius of a product centre come
/// from a fast multipole expansion (detail::FarFieldExpansion) whose tolerances keep its error in
/// every element of V below the threshold (in addition to the screening).
class DipolePotentialDriver {
 public:
  /// Dipole count from which ChargeSummation::automatic uses the multipole expansion.
  static constexpr std::size_t multipole_dipole_count = 600;

  /// `block_size`: atom pairs per work block handed to a thread; 0 picks it from the system size
  /// and thread count (automatic_block_size).
  explicit DipolePotentialDriver(std::size_t block_size = 0,
                                 ChargeSummation summation = ChargeSummation::automatic);

  /// Symmetric matrix V of `basis` for the dipoles `dipoles` at `dipole_coordinates` (bohr,
  /// atomic units). Throws std::invalid_argument if the molecule and basis have different
  /// numbers of atoms, the dipole and coordinate counts differ, a dipole component or coordinate
  /// is not finite, the threshold is negative or not finite, or ChargeSummation::multipole is
  /// requested with a zero threshold.
  auto compute(const Molecule<double>& molecule, const MolecularBasis& basis,
               std::span<const Dipole> dipoles, std::span<const Point3D<double>> dipole_coordinates,
               double threshold) const -> BlockSparseMatrix;

 private:
  std::size_t block_size_;
  ChargeSummation summation_;
};

}  // namespace fika

#endif  // fika_dipole_potential_driver_hpp

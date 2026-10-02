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

#ifndef fika_nuclear_attraction_driver_hpp
#define fika_nuclear_attraction_driver_hpp

#include <cstddef>
#include <span>

#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_summation.hpp"
#include "fika_block_sparse_matrix.hpp"
#include "fika_classical_system.hpp"
#include "fika_molecule.hpp"

namespace fika {

// ChargeSummation (core/summation.hpp): the (A|V|B) kernels sum the point charges directly per
// primitive pair, through a fast multipole far-field expansion (near charges still summed
// directly), or automatically (the expansion for many charges and a positive threshold).

/// Two-centre nuclear-attraction (point-charge potential) integrals
///   V_ab = sum_C q_C <a| 1 / |r - C| |b>
/// over contracted real solid-harmonic Gaussians, for point charges q_C at C.
///
/// Work is distributed over OpenMP threads by atom-pair blocks (AtomPairGroupFactory). The
/// charges break rotational symmetry, so blocks of equal atoms hold every shell pair of the atom
/// (DiagonalFormat::full) and are always stored; the threshold screens only blocks of different
/// atoms, with a bound independent of the charge positions (scaled by sum_C |q_C|).
///
/// Blocks of equal atoms accumulate the charges at the primitive level (charges far from the
/// atom through its multipole field); blocks of different atoms accumulate the charge field at
/// each primitive product centre (Scheme II), far charges through the asymptotic Coulomb kernel.
///
/// With ChargeSummation::multipole the charges beyond the penetration radius of a product centre
/// come from a fast multipole expansion (detail::FarFieldExpansion) whose field-tensor tolerances
/// keep its error in every element of V below the threshold (in addition to the screening).
class NuclearAttractionDriver {
 public:
  /// Charge count from which ChargeSummation::automatic uses the multipole expansion.
  static constexpr std::size_t multipole_charge_count = 1000;

  /// `block_size`: atom pairs per work block handed to a thread; 0 picks it from the system size
  /// and thread count (automatic_block_size).
  explicit NuclearAttractionDriver(std::size_t block_size = 0,
                                   ChargeSummation summation = ChargeSummation::automatic);

  /// Symmetric matrix V of `basis` for the point charges `charges` at `charge_coordinates`
  /// (bohr). Throws std::invalid_argument if the molecule and basis have different numbers of
  /// atoms, the charge and coordinate counts differ, a charge or coordinate is not finite, the
  /// threshold is negative or not finite, or ChargeSummation::multipole is requested with a zero
  /// threshold.
  auto compute(const Molecule<double>& molecule, const MolecularBasis& basis,
               std::span<const double> charges, std::span<const Point3D<double>> charge_coordinates,
               double threshold) const -> BlockSparseMatrix;

  /// Electron-nucleus attraction: the charges are q_A = -Z_A at the atoms of `molecule`.
  auto compute(const Molecule<double>& molecule, const MolecularBasis& basis,
               double threshold) const -> BlockSparseMatrix;

  /// Electron-MM attraction: the charges are q = -q_C at the atoms of both regions of `system`,
  /// with q_C the permanent charges of their force fields (classical_charges), so the matrix adds
  /// to the core Hamiltonian as is. Throws as the charge overload, and as classical_charges.
  auto compute(const Molecule<double>& molecule, const MolecularBasis& basis,
               const ClassicalSystem& system, double threshold) const -> BlockSparseMatrix;

 private:
  std::size_t block_size_;
  ChargeSummation summation_;
};

}  // namespace fika

#endif  // fika_nuclear_attraction_driver_hpp

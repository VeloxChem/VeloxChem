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

#ifndef fika_embedding_hpp
#define fika_embedding_hpp

#include <optional>
#include <span>

#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_summation.hpp"
#include "fika_dense_matrix.hpp"
#include "fika_classical_system.hpp"
#include "fika_dipole_interaction.hpp"
#include "fika_mm_induced_dipoles.hpp"
#include "fika_molecule.hpp"

namespace fika {

/// Which sources polarize the MM region: all of them (the ground state: MM permanent charges, QM
/// nuclei and electrons), or the electrons of the density alone (response theory: a perturbed
/// density D^1 induces mu^1 with B mu^1 = F_e(D^1)).
enum class QmmmSources { all, electrons_only };

struct QmmmEmbeddingOptions {
  InducedDipoleOptions induced;  // solver, field accuracy and summation of the fields
  QmmmSources sources = QmmmSources::all;
  double fock_threshold = 1e-12;  // screening threshold of the Fock integrals
  ChargeSummation fock_summation = ChargeSummation::automatic;
};

/// QM/MM embedding of a QM density in a classical system: the induced dipoles and the energy and
/// Fock contributions, by component. Fock matrices are symmetric, over the AO basis in VeloxChem's
/// order; the permanent ones are empty (0 x 0) and the permanent energies zero for
/// QmmmSources::electrons_only.
struct QmmmEmbedding {
  InducedDipoles induced;  // dipoles (order of polarizable_sites), field, solver statistics
  ChargeSummation electron_field_summation = ChargeSummation::direct;
  int electron_field_order = 0;  // final FMM order of the electron field (0: direct)

  double electron_permanent_energy = 0.0;                    // Tr(D V_es)
  double nuclear_permanent_energy = 0.0;                     // sum_A sum_C Z_A q_C / |R_A - C|
  double polarization_energy = 0.0;                          // -1/2 sum_s mu_s . F_s
  DenseMatrix permanent_fock{0, MatrixSymmetry::symmetric};  // V_es = -sum_C q_C <a|1/|r-C||b>
  DenseMatrix induced_fock{0, MatrixSymmetry::symmetric};    // -sum_s mu_s . <a|(r-s)/|r-s|^3|b>

  /// Sum of the energy components (Hartree).
  auto energy() const noexcept -> double {
    return electron_permanent_energy + nuclear_permanent_energy + polarization_energy;
  }

  /// permanent_fock + induced_fock (induced_fock alone when there is no permanent part).
  auto fock() const -> DenseMatrix;
};

/// Embedding of the density `density` (alpha + beta, over `basis` in VeloxChem's order; only its
/// symmetric part acts) of `molecule` in `system`. With QmmmSources::all the energy
///   E(D) = Tr(D V_es) + E_nuc-MM - 1/2 sum_s mu_s(D) . F_s(D),
/// F the field of the MM permanent charges, the QM nuclei and the electrons (qmmm_induced_dipoles),
/// is stationary in the dipoles, so dE/dD = V_es + V_ind with V_ind = induced_dipole_fock(mu).
/// With QmmmSources::electrons_only the dipoles respond to the electron field alone and only
/// induced_fock is formed. The MM-MM energies (independent of D) are not included. Throws as
/// qmmm_induced_dipoles and the integral drivers.
auto qmmm_embedding(const Molecule<double>& molecule, const MolecularBasis& basis,
                    const DenseMatrix& density, const ClassicalSystem& system,
                    std::optional<TholeDamping> damping, const QmmmEmbeddingOptions& options = {})
    -> QmmmEmbedding;

/// Fock contribution of point dipoles `dipoles` at `positions` (bohr, a.u.) on the electrons,
///   V_ab = -sum_s <a| mu_s . (r - s) / |r - s|^3 |b>
/// (DipolePotentialDriver with -mu), symmetric over `basis` in VeloxChem's order.
auto induced_dipole_fock(const Molecule<double>& molecule, const MolecularBasis& basis,
                         std::span<const Point3D<double>> positions,
                         std::span<const Point3D<double>> dipoles, double threshold = 1e-12,
                         ChargeSummation summation = ChargeSummation::automatic) -> DenseMatrix;

}  // namespace fika

#endif  // fika_embedding_hpp

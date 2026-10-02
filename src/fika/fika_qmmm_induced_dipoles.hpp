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

#ifndef fika_qmmm_induced_dipoles_hpp
#define fika_qmmm_induced_dipoles_hpp

#include <optional>

#include "fika_molecular_basis.hpp"
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

/// Induced dipoles of a QM/MM system and how the field of the QM electrons was summed.
struct QmmmInducedDipoles {
  InducedDipoles induced;  // dipoles, the field they respond to, and the solver statistics
  ChargeSummation electron_field_summation = ChargeSummation::direct;
  int electron_field_order = 0;  // final FMM order (0 for the direct sum)
};

/// Induced dipoles of the polarizable region of `system` in the field of the permanent charges of
/// both regions and of the QM region: its nuclei (QmNuclearField) and electrons
/// (QmElectronicField, total density `density` over `basis` in VeloxChem's AO order). The electron
/// field meets options.field_accuracy (0: tolerance / 10), as the MM field does, and is summed as
/// options.summation selects (automatic: its own multipole crossover). With
/// QmmmSources::electrons_only the dipoles respond to the electrons alone (no permanent charges or
/// nuclei). Throws as induced_dipoles, QmNuclearField and QmElectronicField.
auto qmmm_induced_dipoles(const Molecule<double>& molecule, const MolecularBasis& basis,
                          const DenseMatrix& density, const ClassicalSystem& system,
                          std::optional<TholeDamping> damping,
                          const InducedDipoleOptions& options = {},
                          QmmmSources sources = QmmmSources::all) -> QmmmInducedDipoles;

}  // namespace fika

#endif  // fika_qmmm_induced_dipoles_hpp

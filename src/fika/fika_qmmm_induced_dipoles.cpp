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

#include "fika_qmmm_induced_dipoles.hpp"

#include <array>

#include "fika_qm_field.hpp"

namespace fika {

auto qmmm_induced_dipoles(const Molecule<double>& molecule, const MolecularBasis& basis,
                          const DenseMatrix& density, const ClassicalSystem& system,
                          std::optional<TholeDamping> damping, const InducedDipoleOptions& options)
    -> QmmmInducedDipoles {
  const double accuracy =
      options.field_accuracy > 0.0 ? options.field_accuracy : 0.1 * options.tolerance;
  const QmNuclearField nuclei(molecule);
  const QmElectronicField electrons(molecule, basis, density, accuracy, options.summation);
  const std::array<const FieldContribution*, 2> external{&nuclei, &electrons};
  QmmmInducedDipoles result;
  result.induced = induced_dipoles(system, damping, options, external);
  result.electron_field_order = electrons.report().order;
  result.electron_field_summation =
      result.electron_field_order > 0 ? ChargeSummation::multipole : ChargeSummation::direct;
  return result;
}

}  // namespace fika

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

#ifndef fika_point_sources_hpp
#define fika_point_sources_hpp

#include <cstddef>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_classical_system.hpp"
#include "fika_molecule.hpp"

namespace fika {

/// Point charges and their positions (bohr). From a classical system, the charges of its k-th
/// residue (polarizable region first, then the nonpolarizable one) are
/// [residue_offsets[k], residue_offsets[k + 1]).
struct PointCharges {
  std::vector<double> charges;
  std::vector<Point3D<double>> coordinates;
  std::vector<std::size_t> residue_offsets;
};

/// Permanent charges of every atom of a classical system, from the force field of each residue
/// (loaded from the force-field library once per force field and residue name): the polarizable
/// region first, then the nonpolarizable one, residues and atoms in order (zero charges kept).
/// Throws std::invalid_argument if a residue's atom count differs from its force field's, and
/// std::runtime_error if a force field is not found.
auto classical_charges(const ClassicalSystem& system) -> PointCharges;

/// Interaction energy (Hartree) of the nuclei of `molecule` with the permanent charges of both
/// regions of `system`: sum_A sum_C Z_A q_C / |R_A - C| (nuclear_point_charge_energy). Throws as
/// classical_charges and nuclear_point_charge_energy.
auto nuclear_classical_energy(const Molecule<double>& molecule, const ClassicalSystem& system)
    -> double;

}  // namespace fika

#endif  // fika_point_sources_hpp

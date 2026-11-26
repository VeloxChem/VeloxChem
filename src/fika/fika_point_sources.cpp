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

#include "fika_point_sources.hpp"

#include "fika_force_field_cache.hpp"
#include "fika_residue.hpp"
#include "fika_nuclear_repulsion.hpp"

namespace fika {

auto classical_charges(const ClassicalSystem& system) -> PointCharges {
  const detail::ForceFieldCache cache(system, "fika::classical_charges");
  PointCharges result;
  result.residue_offsets.push_back(0);
  for (const ClassicalRegion* region :
       {&system.polarizable_region(), &system.nonpolarizable_region()}) {
    for (const Residue& residue : region->residues()) {
      const auto charges = cache.of(residue).charges();
      const auto coordinates = residue.molecule().coordinates();
      result.charges.insert(result.charges.end(), charges.begin(), charges.end());
      result.coordinates.insert(result.coordinates.end(), coordinates.begin(), coordinates.end());
      result.residue_offsets.push_back(result.charges.size());
    }
  }
  return result;
}

auto nuclear_classical_energy(const Molecule<double>& molecule, const ClassicalSystem& system)
    -> double {
  const PointCharges sources = classical_charges(system);
  return nuclear_point_charge_energy(molecule, sources.charges, sources.coordinates);
}

}  // namespace fika

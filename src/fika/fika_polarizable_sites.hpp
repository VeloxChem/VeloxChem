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

#ifndef fika_polarizable_sites_hpp
#define fika_polarizable_sites_hpp

#include <cstddef>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_classical_system.hpp"
#include "fika_force_field.hpp"

namespace fika {

/// Polarizable sites (bohr, polarizabilities in a.u.). owners[i] is the residue of site i as a
/// classical-system residue index: polarizable region first, then the nonpolarizable one (the
/// order of classical_charges).
struct PolarizableSites {
  std::vector<Point3D<double>> positions;
  std::vector<Polarizability> polarizabilities;  // 3 x 3 (isotropic alpha as alpha I)
  std::vector<std::size_t> owners;
};

/// The atoms of the polarizable region with a nonzero polarizability in their force field, residue
/// by residue in atom order. Throws std::invalid_argument for anisotropic polarizabilities (they
/// need the force field mapped onto the residue geometry, not supported yet) or a residue whose
/// atom count differs from its force field's, and std::runtime_error if a force field is not
/// found.
auto polarizable_sites(const ClassicalSystem& system) -> PolarizableSites;

}  // namespace fika

#endif  // fika_polarizable_sites_hpp

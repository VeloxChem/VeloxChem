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

#ifndef fika_nuclear_repulsion_hpp
#define fika_nuclear_repulsion_hpp

#include <span>

#include "fika_point3d.hpp"
#include "fika_molecule.hpp"

namespace fika {

// Defined for float and double in nuclear_repulsion.cpp. Energies are accumulated and returned in
// double for both; coincident nuclei throw std::domain_error.

/// Nuclear repulsion energy (Hartree) of a molecule: sum over atom pairs of Z_i Z_j / r_ij.
template <real_scalar T>
auto nuclear_repulsion_energy(const Molecule<T>& molecule) -> double;

/// Intermolecular nuclear repulsion energy (Hartree): sum of Z_i Z_j / r_ij, i in a, j in b.
template <real_scalar T>
auto nuclear_repulsion_energy(const Molecule<T>& a, const Molecule<T>& b) -> double;

/// Interaction energy (Hartree) of the nuclei of a molecule with point charges:
/// sum_A sum_C Z_A q_C / |R_A - C|, charges q_C at `coordinates` (bohr). The charges are cut into
/// fixed blocks summed in parallel, each block and the blocks in order with compensated
/// (Neumaier) summation, so the result is independent of the thread count. Throws
/// std::invalid_argument if the charge and coordinate counts differ or a charge or coordinate is
/// not finite, and std::domain_error if a charge sits on a nucleus.
auto nuclear_point_charge_energy(const Molecule<double>& molecule, std::span<const double> charges,
                                 std::span<const Point3D<double>> coordinates) -> double;

}  // namespace fika

#endif  // fika_nuclear_repulsion_hpp

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

#ifndef fika_veloxchem_hpp
#define fika_veloxchem_hpp

#include "AtomBasis.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "fika_molecular_basis.hpp"
#include "fika_molecule.hpp"

namespace fika {

/// Molecule of VeloxChem's `molecule`: elements from the atom identifiers, coordinates in bohr.
/// Throws std::invalid_argument for an atom without nuclear charge (a ghost atom, identifier 0),
/// which fika does not support.
auto from_veloxchem(const CMolecule& molecule) -> Molecule<double>;

/// Atom basis of VeloxChem's `basis`: one shell per basis function (VeloxChem has no general
/// contractions), in VeloxChem's order within each angular momentum, with VeloxChem's
/// normalization factors as effective coefficients (not renormalized). Throws
/// std::invalid_argument for an identifier without nuclear charge or an effective core potential.
auto from_veloxchem(const CAtomBasis& basis) -> AtomBasis;

/// Molecular basis of VeloxChem's `basis` for `molecule`: its unique atom bases and their
/// assignment to atoms. With veloxchem_order the basis functions keep VeloxChem's AO order.
/// Throws std::invalid_argument if the basis does not match the molecule (atom count or
/// elements) or as the atom-basis overload.
auto from_veloxchem(const CMolecularBasis& basis, const CMolecule& molecule) -> MolecularBasis;

}  // namespace fika

#endif  // fika_veloxchem_hpp

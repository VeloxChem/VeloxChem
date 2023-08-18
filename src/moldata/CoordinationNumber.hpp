//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#ifndef CoordinationNumber_hpp
#define CoordinationNumber_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

namespace coordnum {  // coordnum namespace

/**
 Creates covalent radii (reference: dftd4 v2.4.0).

 @return a vector of covalent radii with nuclear charge as index.
 */
std::vector<double> getCovalentRadius();

/**
 Creates atomic coordination numbers (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @return a vector of atomic coordination numbers for a molecule.
 */
std::vector<double> getCoordinationNumber(const CMolecule& molecule);

/**
 Creates atomic coordination numbers (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param dcndr the derivative matrix of dimension (3N,N).
 @return a vector of atomic coordination numbers for a molecule.
 */
std::vector<double> getCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcndr);

/**
 Creates Pauling electronegativity (reference: dftd4 v2.4.0).

 @return a vector of Pauling electronegativity with nuclear charge as index.
 */
std::vector<double> getPaulingElectronegativity();

/**
 Creates atomic covalent coordination numbers (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param dcovcndr the derivative matrix of dimension (3N,N).
 @return a vector of atomic covalent coordination numbers for a molecule.
 */
std::vector<double> getCovalentCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcovcndr);

}  // namespace coordnum

#endif /* CoordinationNumber_hpp */

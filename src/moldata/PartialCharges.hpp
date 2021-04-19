//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#ifndef PartialCharges_hpp
#define PartialCharges_hpp

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

namespace parchg {  // parchg namespace

/**
 Creates atomic partial charges (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param netcharge net charge of the molecule.
 @return a vector of atomic partial charges for a molecule.
 */
std::vector<double> getPartialCharges(const CMolecule& molecule, const double netcharge);

/**
 Creates atomic partial charges (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param netcharge net charge of the molecule.
 @param dqdr the derivative matrix of dimension (3N,N+1).
 @return a vector of atomic partial charges for a molecule.
 */
std::vector<double> getPartialCharges(const CMolecule& molecule, const double netcharge, CDenseMatrix& dqdr);

}  // namespace parchg

#endif /* PartialCharges_hpp */

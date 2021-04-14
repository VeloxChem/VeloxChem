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

#ifndef ExportOrbData_hpp
#define ExportOrbData_hpp

#include <pybind11/pybind11.h>

class CMolecule;
class CMolecularBasis;

namespace py = pybind11;

namespace vlx_orbdata {  // vlx_orbdata namespace

/**
 Gets number of atomic orbitals.

 @param molecule the molecule.
 @param basis the AO basis set.
 @return the number of atomic orbitals.
 */
int32_t get_number_of_atomic_orbitals(const CMolecule& molecule, const CMolecularBasis& basis);

/**
 Exports classes/functions in src/orbdata to python.
 */
void export_orbdata(py::module& m);

}  // namespace vlx_orbdata

#endif /* ExportOrbData_hpp */

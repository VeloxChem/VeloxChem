//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef XCMolecularHessianForGGA_hpp
#define XCMolecularHessianForGGA_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCFunctional.hpp"

namespace xchessgga {  // xchessgga namespace

/**
 Integrates GGA exchnage-correlation functional contribution to molecular
 Hessian.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular Hessian.
 */
auto integrateExcHessianForGGA(const CMolecule&                  molecule,
                               const CMolecularBasis&            basis,
                               const std::vector<const double*>& gsDensityPointers,
                               const CMolecularGrid&             molecularGrid,
                               const double                      screeningThresholdForGTOValues,
                               const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates GGA exchnage-correlation functional contribution to molecular
 gradient of Vxc matrix element.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @param atomIdx the index of the atom with respect to which gradient is
 computed.
 @return the Vxc gradient.
 */
auto integrateVxcFockGradientForGGA(const CMolecule&                  molecule,
                                    const CMolecularBasis&            basis,
                                    const std::vector<const double*>& gsDensityPointers,
                                    const CMolecularGrid&             molecularGrid,
                                    const double                      screeningThresholdForGTOValues,
                                    const CXCFunctional&              xcFunctional,
                                    const int32_t                     atomIdx) -> std::vector<CDenseMatrix>;

}  // namespace xchessgga

#endif /* XCMolecularHessianForGGA_hpp */

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

#ifndef XCMolecularHessianForMGGA_hpp
#define XCMolecularHessianForMGGA_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCFunctional.hpp"

namespace xchessmgga {  // xchessmgga namespace

/**
 Integrates meta-GGA exchnage-correlation functional contribution to
 closed-shell molecular Hessian.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular Hessian.
 */
auto integrateExcHessianForMetaGgaClosedShell(const CMolecule&                  molecule,
                                              const CMolecularBasis&            basis,
                                              const std::vector<const double*>& gsDensityPointers,
                                              const CMolecularGrid&             molecularGrid,
                                              const double                      screeningThresholdForGTOValues,
                                              const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates meta-GGA exchnage-correlation functional contribution to molecular
 gradient of Vxc matrix element.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @param atomIdxVec the indices of the atoms with respect to which gradient is
 computed.
 @return the Vxc gradient.
 */
auto integrateVxcFockGradientForMetaGGA(const CMolecule&                  molecule,
                                        const CMolecularBasis&            basis,
                                        const std::vector<const double*>& gsDensityPointers,
                                        const CMolecularGrid&             molecularGrid,
                                        const double                      screeningThresholdForGTOValues,
                                        const CXCFunctional&              xcFunctional,
                                        const std::vector<int>&           atomIdxVec) -> std::vector<CDenseMatrix>;

}  // namespace xchessmgga

#endif /* XCMolecularHessianForMGGA_hpp */

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

#ifndef XCMolecularGradientForGGA_hpp
#define XCMolecularGradientForGGA_hpp

#include <array>
#include <list>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GridBox.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcgradgga {  // xcgradgga namespace

/**
 Integrates first-order GGA exchnage-correlation functional contribution to
 closed-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityMatrix the perturbed AO density matrix.
 @param gsDensityMatrix the ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateVxcGradientForGgaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates first-order GGA exchnage-correlation functional contribution to
 open-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed AO density matrix.
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateVxcGradientForGgaOpenShell(const CMolecule&                  molecule,
                                         const CMolecularBasis&            basis,
                                         const std::vector<const double*>& rwDensityPointers,
                                         const std::vector<const double*>& gsDensityPointers,
                                         const CMolecularGrid&             molecularGrid,
                                         const double                      screeningThresholdForGTOValues,
                                         const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates second-order GGA exchnage-correlation functional contribution
 to closed-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointersOne the pointers to perturbed AO density matrix.
 @param rwDensityPointersTwo the pointers to perturbed AO density matrix (to be
        contracted with GTO gradient).
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateFxcGradientForGgaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates third-order GGA exchnage-correlation functional contribution to
 closed-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointersOne the pointers to perturbed AO density matrix.
 @param rwDensityPointersTwo the pointers to perturbed AO density matrix.
 @param gsDensityPointers the pointers to ground state AO density matrix (to be
        contracted with GTO gradient).
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateKxcGradientForGgaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CDenseMatrix;

}  // namespace xcgradgga

#endif /* XCMolecularGradientForGGA_hpp */

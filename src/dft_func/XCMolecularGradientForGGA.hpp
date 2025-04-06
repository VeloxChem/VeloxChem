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
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

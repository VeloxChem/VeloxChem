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

#ifndef XCMolecularGradientForLDA_hpp
#define XCMolecularGradientForLDA_hpp

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

namespace xcgradlda {  // xcgradlda namespace

/**
 Integrates first-order LDA exchnage-correlation functional contribution to
 closed-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityMatrix the perturbed AO density matrix (to be contracted
        with GTO gradient).
 @param gsDensityMatrix the ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateVxcGradientForLdaClosedShell(const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&   molecularGrid,
                                           const double            screeningThresholdForGTOValues,
                                           const CXCFunctional&    xcFunctional) -> CDenseMatrix;

/**
 Integrates first-order LDA exchnage-correlation functional contribution to
 open-shell molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed AO density matrix (to be contracted
        with GTO gradient).
 @param gsDensityPointers the pointers to ground state AO density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @return the molecular gradient.
 */
auto integrateVxcGradientForLdaOpenShell(const CMolecule&                  molecule,
                                         const CMolecularBasis&            basis,
                                         const std::vector<const double*>& rwDensityPointers,
                                         const std::vector<const double*>& gsDensityPointers,
                                         const CMolecularGrid&             molecularGrid,
                                         const double                      screeningThresholdForGTOValues,
                                         const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates second-order LDA exchnage-correlation functional contribution
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
auto integrateFxcGradientForLdaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CDenseMatrix;

/**
 Integrates third-order LDA exchnage-correlation functional contribution to
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
auto integrateKxcGradientForLdaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointersOne,
                                           const std::vector<const double*>& rwDensityPointersTwo,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CDenseMatrix;

}  // namespace xcgradlda

#endif /* XCMolecularGradientForLDA_hpp */

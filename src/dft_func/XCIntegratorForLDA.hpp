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

#ifndef XCIntegratorForLDA_hpp
#define XCIntegratorForLDA_hpp

#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridCubic.hpp"
#include "DensityGridQuad.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcintlda {  // xcintlda namespace

/**
 Integrates first-order LDA exchange-correlation functional contribution to
 closed-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForLdaClosedShell(const CMolecule&                  molecule,
                                       const CMolecularBasis&            basis,
                                       const std::vector<const double*>& gsDensityPointers,
                                       const CMolecularGrid&             molecularGrid,
                                       const double                      screeningThresholdForGTOValues,
                                       const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix;

/**
 Integrates first-order LDA exchange-correlation functional contribution to
 open-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForLdaOpenShell(const CMolecule&                  molecule,
                                     const CMolecularBasis&            basis,
                                     const std::vector<const double*>& gsDensityPointers,
                                     const CMolecularGrid&             molecularGrid,
                                     const double                      screeningThresholdForGTOValues,
                                     const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix;

/**
 Integrates second-order LDA exchange-correlation functional contribution
 to closed-shell AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed density matrices.
 @param gsDensityPointers the pointers to ground state density matrices.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 */
auto integrateFxcFockForLdaClosedShell(const std::vector<double*>&       aoFockPointers,
                                       const CMolecule&                  molecule,
                                       const CMolecularBasis&            basis,
                                       const std::vector<const double*>& rwDensityPointers,
                                       const std::vector<const double*>& gsDensityPointers,
                                       const CMolecularGrid&             molecularGrid,
                                       const double                      screeningThresholdForGTOValues,
                                       const CXCFunctional&              xcFunctional) -> void;

/**
 Integrates third-order LDA exchange-correlation functional contribution
 to closed-shell AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed one-time transformed densities.
 @param rw2DensityPointers the pointers to two-time transformed densities.
 @param gsDensityPointers the pointers to ground state density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @param quadMode a string that specifies which densities should be combined.
 */
auto integrateKxcFockForLdaClosedShell(const std::vector<double*>& aoFockPointers,
                                       const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const std::vector<const double*>& rwDensityPointers,
                                       const std::vector<const double*>& rw2DensityPointers,
                                       const std::vector<const double*>& gsDensityPointers,
                                       const CMolecularGrid&   molecularGrid,
                                       const double            screeningThresholdForGTOValues,
                                       const CXCFunctional&    xcFunctional,
                                       const std::string&      quadMode) -> void;

/**
 Integrates fourth-order LDA exchnage-correlation functional contribution
 to closed-shell AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed one-time transformed densities.
 @param rw2DensityPointers the pointers to two-time transformed densities.
 @param rw3DensityPointers the pointers to three-time transformed densities.
 @param gsDensityPointers the pointers to ground state density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @param cubeMode a string that specifies which densities should be combined.
 */
auto integrateKxcLxcFockForLdaClosedShell(const std::vector<double*>& aoFockPointers,
                                          const CMolecule&            molecule,
                                          const CMolecularBasis&      basis,
                                          const std::vector<const double*>& rwDensityPointers,
                                          const std::vector<const double*>& rw2DensityPointers,
                                          const std::vector<const double*>& rw3DensityPointers,
                                          const std::vector<const double*>& gsDensityPointers,
                                          const CMolecularGrid&       molecularGrid,
                                          const double                screeningThresholdForGTOValues,
                                          const CXCFunctional&        xcFunctional,
                                          const std::string&          cubeMode) -> void;

/**
 Integrates LDA contribution to (third-order) closed-shell Kxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param v3rho3 the 3rd-order functional derivative wrt density.
 @param rwDensityGridPointers the pointers to the products of one and two-time
        transformed densities on grid points.
 @param rw2DensityGridPointers the pointers to the two-time transformed
        densities on grid points.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialKxcFockForLdaClosedShell(const CXCFunctional&     xcFunctional,
                                              const double*            weights,
                                              const CDenseMatrix&      gtoValues,
                                              const double*            v2rho2,
                                              const double*            v3rho3,
                                              const std::vector<const double*>& rwDensityGridPointers,
                                              const std::vector<const double*>& rw2DensityGridPointers,
                                              CMultiTimer&             timer) -> CDenseMatrix;

/**
 Integrates LDA contribution to (fourth-order) closed-shell Lxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param v3rho3 the 3rd-order functional derivative wrt density.
 @param v4rho4 the 4rd-order functional derivative wrt density.
 @param rwDensityGridCubic the products of one and two-time transformed densities on grid points.
 @param rw3DensityGrid the three-time transformed densities on grid points.
 @param iFock the index of the AO Fock matrix.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialLxcFockForLdaClosedShell(const CXCFunctional&     xcFunctional,
                                              const double*            weights,
                                              const CDenseMatrix&      gtoValues,
                                              const double*            v2rho2,
                                              const double*            v3rho3,
                                              const double*            v4rho4,
                                              const CDensityGridCubic& rwDensityGridCubic,
                                              const CDensityGrid&      rw3DensityGrid,
                                              const int                iFock,
                                              CMultiTimer&             timer) -> CDenseMatrix;

}  // namespace xcintlda

#endif /* XCIntegratorForLDA_hpp */

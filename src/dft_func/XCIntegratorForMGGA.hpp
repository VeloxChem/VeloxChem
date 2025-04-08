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

#ifndef XCIntegratorForMGGA_hpp
#define XCIntegratorForMGGA_hpp

#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridCubic.hpp"
#include "DensityGridQuad.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcintmgga {  // xcintmgga namespace

/**
 Integrates first-order meta-GGA exchange-correlation functional
 contribution to closed-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForMetaGgaClosedShell(const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix;

/**
 Integrates first-order meta-GGA exchange-correlation functional
 contribution to open-shell AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForMetaGgaOpenShell(const CMolecule&                  molecule,
                                         const CMolecularBasis&            basis,
                                         const std::vector<const double*>& gsDensityPointers,
                                         const CMolecularGrid&             molecularGrid,
                                         const double                      screeningThresholdForGTOValues,
                                         const CXCFunctional&              xcFunctional) -> CAOKohnShamMatrix;

/**
 Integrates second-order meta-GGA exchange-correlation functional
 contribution to AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed density matrices.
 @param gsDensityPointers the pointers to ground state density matrices.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 */
auto integrateFxcFockForMetaGgaClosedShell(const std::vector<double*>&       aoFockPointers,
                                           const CMolecule&                  molecule,
                                           const CMolecularBasis&            basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&             molecularGrid,
                                           const double                      screeningThresholdForGTOValues,
                                           const CXCFunctional&              xcFunctional) -> void;

/**
 Integrates third-order meta-GGA exchange-correlation functional
 contribution to AO Fock matrix.

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
auto integrateKxcFockForMetaGgaClosedShell(const std::vector<double*>& aoFockPointers,
                                           const CMolecule&            molecule,
                                           const CMolecularBasis&      basis,
                                           const std::vector<const double*>& rwDensityPointers,
                                           const std::vector<const double*>& rw2DensityPointers,
                                           const std::vector<const double*>& gsDensityPointers,
                                           const CMolecularGrid&       molecularGrid,
                                           const double                screeningThresholdForGTOValues,
                                           const CXCFunctional&        xcFunctional,
                                           const std::string&          quadMode) -> void;

/**
 Integrates fourth-order meta-GGA exchange-correlation functional
 contribution to AO Fock matrix.

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
auto integrateKxcLxcFockForMetaGgaClosedShell(const std::vector<double*>& aoFockPointers,
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
 Integrates meta-GGA contribution to (third-order) Kxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the density gradient.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma ,
 @param v2rholapl ,
 @param v2rhotau ,
 @param v2sigma2 ,
 @param v2sigmalapl ,
 @param v2sigmatau ,
 @param v2lapl2 ,
 @param v2lapltau ,
 @param v2tau2 ,
 @param v3rho3 ,
 @param v3rho2sigma ,
 @param v3rho2lapl ,
 @param v3rho2tau ,
 @param v3rhosigma2 ,
 @param v3rhosigmalapl ,
 @param v3rhosigmatau ,
 @param v3rholapl2 ,
 @param v3rholapltau ,
 @param v3rhotau2 ,
 @param v3sigma3 ,
 @param v3sigma2lapl ,
 @param v3sigma2tau ,
 @param v3sigmalapl2 ,
 @param v3sigmalapltau ,
 @param v3sigmatau2 ,
 @param v3lapl3 ,
 @param v3lapl2tau ,
 @param v3lapltau2 ,
 @param v3tau3 ,
 @param rwDensityGridPointers the pointers to the products of one-time
        transformed densities on grid points.
 @param rw2DensityGridPointers the pointers to the two-time transformed
        densities on grid points.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialKxcFockForMetaGgaClosedShell(const CXCFunctional&     xcFunctional,
                                                  const double*            weights,
                                                  const CDenseMatrix&      gtoValues,
                                                  const CDenseMatrix&      gtoValuesX,
                                                  const CDenseMatrix&      gtoValuesY,
                                                  const CDenseMatrix&      gtoValuesZ,
                                                  const double*            rhograd,
                                                  const double*            vsigma,
                                                  const double*            v2rho2,
                                                  const double*            v2rhosigma,
                                                  const double*            v2rholapl,
                                                  const double*            v2rhotau,
                                                  const double*            v2sigma2,
                                                  const double*            v2sigmalapl,
                                                  const double*            v2sigmatau,
                                                  const double*            v2lapl2,
                                                  const double*            v2lapltau,
                                                  const double*            v2tau2,
                                                  const double*            v3rho3,
                                                  const double*            v3rho2sigma,
                                                  const double*            v3rho2lapl,
                                                  const double*            v3rho2tau,
                                                  const double*            v3rhosigma2,
                                                  const double*            v3rhosigmalapl,
                                                  const double*            v3rhosigmatau,
                                                  const double*            v3rholapl2,
                                                  const double*            v3rholapltau,
                                                  const double*            v3rhotau2,
                                                  const double*            v3sigma3,
                                                  const double*            v3sigma2lapl,
                                                  const double*            v3sigma2tau,
                                                  const double*            v3sigmalapl2,
                                                  const double*            v3sigmalapltau,
                                                  const double*            v3sigmatau2,
                                                  const double*            v3lapl3,
                                                  const double*            v3lapl2tau,
                                                  const double*            v3lapltau2,
                                                  const double*            v3tau3,
                                                  const std::vector<const double*>& rwDensityGridPointers,
                                                  const std::vector<const double*>& rw2DensityGridPointers,
                                                  CMultiTimer&             timer) -> CDenseMatrix;

/**
 Integrates meta-GGA contribution to (fourth-order) Lxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the density gradient.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 ,
 @param v2rhosigma ,
 @param v2rholapl ,
 @param v2rhotau ,
 @param v2sigma2 ,
 @param v2sigmalapl ,
 @param v2sigmatau ,
 @param v2lapl2 ,
 @param v2lapltau ,
 @param v2tau2 ,
 @param v3rho3 ,
 @param v3rho2sigma ,
 @param v3rho2lapl ,
 @param v3rho2tau ,
 @param v3rhosigma2 ,
 @param v3rhosigmalapl ,
 @param v3rhosigmatau ,
 @param v3rholapl2 ,
 @param v3rholapltau ,
 @param v3rhotau2 ,
 @param v3sigma3 ,
 @param v3sigma2lapl ,
 @param v3sigma2tau ,
 @param v3sigmalapl2 ,
 @param v3sigmalapltau ,
 @param v3sigmatau2 ,
 @param v3lapl3 ,
 @param v3lapl2tau ,
 @param v3lapltau2 ,
 @param v3tau3 ,
 @param v4rho4 ,
 @param v4rho3sigma ,
 @param v4rho3lapl ,
 @param v4rho3tau ,
 @param v4rho2sigma2 ,
 @param v4rho2sigmalapl ,
 @param v4rho2sigmatau ,
 @param v4rho2lapl2 ,
 @param v4rho2lapltau ,
 @param v4rho2tau2 ,
 @param v4rhosigma3 ,
 @param v4rhosigma2lapl ,
 @param v4rhosigma2tau ,
 @param v4rhosigmalapl2 ,
 @param v4rhosigmalapltau ,
 @param v4rhosigmatau2 ,
 @param v4rholapl3 ,
 @param v4rholapl2tau ,
 @param v4rholapltau2 ,
 @param v4rhotau3 ,
 @param v4sigma4 ,
 @param v4sigma3lapl ,
 @param v4sigma3tau ,
 @param v4sigma2lapl2 ,
 @param v4sigma2lapltau ,
 @param v4sigma2tau2 ,
 @param v4sigmalapl3 ,
 @param v4sigmalapl2tau ,
 @param v4sigmalapltau2 ,
 @param v4sigmatau3 ,
 @param v4lapl4 ,
 @param v4lapl3tau ,
 @param v4lapl2tau2 ,
 @param v4lapltau3 ,
 @param v4tau4 ,
 @param rwDensityGridQuad the products of one-time transformed densities on grid points.
 @param rw3DensityGrid the three-time transformed densities on grid points.
 @param iFock the index of the AO Fock matrix.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialLxcFockForMGGA(const CXCFunctional&     xcFunctional,
                                    const double*            weights,
                                    const CDenseMatrix&      gtoValues,
                                    const CDenseMatrix&      gtoValuesX,
                                    const CDenseMatrix&      gtoValuesY,
                                    const CDenseMatrix&      gtoValuesZ,
                                    const double*            rhograd,
                                    const double*            vsigma,
                                    const double*            v2rho2,
                                    const double*            v2rhosigma,
                                    const double*            v2rholapl,
                                    const double*            v2rhotau,
                                    const double*            v2sigma2,
                                    const double*            v2sigmalapl,
                                    const double*            v2sigmatau,
                                    const double*            v2lapl2,
                                    const double*            v2lapltau,
                                    const double*            v2tau2,
                                    const double*            v3rho3,
                                    const double*            v3rho2sigma,
                                    const double*            v3rho2lapl,
                                    const double*            v3rho2tau,
                                    const double*            v3rhosigma2,
                                    const double*            v3rhosigmalapl,
                                    const double*            v3rhosigmatau,
                                    const double*            v3rholapl2,
                                    const double*            v3rholapltau,
                                    const double*            v3rhotau2,
                                    const double*            v3sigma3,
                                    const double*            v3sigma2lapl,
                                    const double*            v3sigma2tau,
                                    const double*            v3sigmalapl2,
                                    const double*            v3sigmalapltau,
                                    const double*            v3sigmatau2,
                                    const double*            v3lapl3,
                                    const double*            v3lapl2tau,
                                    const double*            v3lapltau2,
                                    const double*            v3tau3,
                                    const double*            v4rho4,
                                    const double*            v4rho3sigma,
                                    const double*            v4rho3lapl,
                                    const double*            v4rho3tau,
                                    const double*            v4rho2sigma2,
                                    const double*            v4rho2sigmalapl,
                                    const double*            v4rho2sigmatau,
                                    const double*            v4rho2lapl2,
                                    const double*            v4rho2lapltau,
                                    const double*            v4rho2tau2,
                                    const double*            v4rhosigma3,
                                    const double*            v4rhosigma2lapl,
                                    const double*            v4rhosigma2tau,
                                    const double*            v4rhosigmalapl2,
                                    const double*            v4rhosigmalapltau,
                                    const double*            v4rhosigmatau2,
                                    const double*            v4rholapl3,
                                    const double*            v4rholapl2tau,
                                    const double*            v4rholapltau2,
                                    const double*            v4rhotau3,
                                    const double*            v4sigma4,
                                    const double*            v4sigma3lapl,
                                    const double*            v4sigma3tau,
                                    const double*            v4sigma2lapl2,
                                    const double*            v4sigma2lapltau,
                                    const double*            v4sigma2tau2,
                                    const double*            v4sigmalapl3,
                                    const double*            v4sigmalapl2tau,
                                    const double*            v4sigmalapltau2,
                                    const double*            v4sigmatau3,
                                    const double*            v4lapl4,
                                    const double*            v4lapl3tau,
                                    const double*            v4lapl2tau2,
                                    const double*            v4lapltau3,
                                    const double*            v4tau4,
                                    const CDensityGridCubic& rwDensityGridCubic,
                                    const CDensityGrid&      rw3DensityGrid,
                                    const int                iFock,
                                    CMultiTimer&             timer) -> CDenseMatrix;

}  // namespace xcintmgga

#endif /* XCIntegratorForMGGA_hpp */

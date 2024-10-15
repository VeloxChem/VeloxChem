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
 contribution to AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @param flag the flag for closed/open shell.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForMGGA(const CMolecule&                  molecule,
                             const CMolecularBasis&            basis,
                             const std::vector<const double*>& gsDensityPointers,
                             const CMolecularGrid&             molecularGrid,
                             const double                      screeningThresholdForGTOValues,
                             const CXCFunctional&              xcFunctional,
                             const std::string&                flag = std::string("closedshell")) -> CAOKohnShamMatrix;

/**
 Integrates meta-GGA contribution to AO Kohn-Sham matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the gradient density.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param vlapl the 1st-order functional derivative wrt laplacian.
 @param vtau the 1st-order functional derivative wrt tau.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialVxcFockForMGGA(const double*       weights,
                                    const CDenseMatrix& gtoValues,
                                    const CDenseMatrix& gtoValuesX,
                                    const CDenseMatrix& gtoValuesY,
                                    const CDenseMatrix& gtoValuesZ,
                                    const double*       rhograd,
                                    const double*       vrho,
                                    const double*       vsigma,
                                    const double*       vlapl,
                                    const double*       vtau,
                                    CMultiTimer&        timer) -> CDenseMatrix;

/**
 Integrates meta-GGA contribution to AO Kohn-Sham matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the gradient density.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param vlapl the 1st-order functional derivative wrt laplacian.
 @param vtau the 1st-order functional derivative wrt tau.
 @param timer the timer.
 @return the alpha and beta contribution as a list of CDenseMatrix objects.
 */
auto integratePartialVxcFockForMGGAOpenShell(const double*       weights,
                                             const CDenseMatrix& gtoValues,
                                             const CDenseMatrix& gtoValuesX,
                                             const CDenseMatrix& gtoValuesY,
                                             const CDenseMatrix& gtoValuesZ,
                                             const double*       rhograd,
                                             const double*       vrho,
                                             const double*       vsigma,
                                             const double*       vlapl,
                                             const double*       vtau,
                                             CMultiTimer&        timer) -> std::vector<CDenseMatrix>;


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
auto integrateFxcFockForMGGA(const std::vector<double*>&       aoFockPointers,
                             const CMolecule&                  molecule,
                             const CMolecularBasis&            basis,
                             const std::vector<const double*>& rwDensityPointers,
                             const std::vector<const double*>& gsDensityPointers,
                             const CMolecularGrid&             molecularGrid,
                             const double                      screeningThresholdForGTOValues,
                             const CXCFunctional&              xcFunctional) -> void;

/**
 Integrates meta-GGA contribution to (second-order) Fxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhow the pointer to perturbed density.
 @param rhograd the pointer to density gradient.
 @param rhowgrad the pointer to perturbed density gradient.
 @param tauw ,
 @param laplw ,
 @param vrho ,
 @param vsigma ,
 @param vlapl ,
 @param vtau ,
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param v2lapl2 ,
 @param v2tau2 ,
 @param v2rholapl ,
 @param v2rhotau ,
 @param v2lapltau ,
 @param v2rhosigma the 2nd-order functional derivative wrt density and
        density gradient.
 @param v2sigmalapl ,
 @param v2sigmatau ,
 @param v2sigma2 the 2nd-order functional derivative wrt density gradient.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialFxcFockForMGGA(const CXCFunctional& xcFunctional,
                                    const double*        weights,
                                    const CDenseMatrix&  gtoValues,
                                    const CDenseMatrix&  gtoValuesX,
                                    const CDenseMatrix&  gtoValuesY,
                                    const CDenseMatrix&  gtoValuesZ,
                                    const double*        rhow,
                                    const double*        rhograd,
                                    const double*        rhowgrad,
                                    const double*        tauw,
                                    const double*        laplw,
                                    const double*        vrho,
                                    const double*        vsigma,
                                    const double*        vlapl,
                                    const double*        vtau,
                                    const double*        v2rho2,
                                    const double*        v2lapl2,
                                    const double*        v2tau2,
                                    const double*        v2rholapl,
                                    const double*        v2rhotau,
                                    const double*        v2lapltau,
                                    const double*        v2rhosigma,
                                    const double*        v2sigmalapl,
                                    const double*        v2sigmatau,
                                    const double*        v2sigma2,
                                    CMultiTimer&         timer) -> CDenseMatrix;

/**
 Integrates third-order meta-GGA exchange-correlation functional
 contribution to AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityMatrix the perturbed one-time transformed densities.
 @param rw2DensityMatrix the two-time transformed densities.
 @param gsDensityMatrix the ground state density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @param quadMode a string that specifies which densities should be combined.
 */
auto integrateKxcFockForMGGA(const std::vector<double*>& aoFockPointers,
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
 @param rwDensityGridQuad the products of one-time transformed densities on grid points.
 @param rw2DensityMatrix the two-time transformed densities on grid points.
 @param iFock the index of the AO Fock matrix.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialKxcFockForMGGA(const CXCFunctional&    xcFunctional,
                                    const double*           weights,
                                    const CDenseMatrix&     gtoValues,
                                    const CDenseMatrix&     gtoValuesX,
                                    const CDenseMatrix&     gtoValuesY,
                                    const CDenseMatrix&     gtoValuesZ,
                                    const double*           rhograd,
                                    const double*           vsigma,
                                    const double*           v2rho2,
                                    const double*           v2rhosigma,
                                    const double*           v2rholapl,
                                    const double*           v2rhotau,
                                    const double*           v2sigma2,
                                    const double*           v2sigmalapl,
                                    const double*           v2sigmatau,
                                    const double*           v2lapl2,
                                    const double*           v2lapltau,
                                    const double*           v2tau2,
                                    const double*           v3rho3,
                                    const double*           v3rho2sigma,
                                    const double*           v3rho2lapl,
                                    const double*           v3rho2tau,
                                    const double*           v3rhosigma2,
                                    const double*           v3rhosigmalapl,
                                    const double*           v3rhosigmatau,
                                    const double*           v3rholapl2,
                                    const double*           v3rholapltau,
                                    const double*           v3rhotau2,
                                    const double*           v3sigma3,
                                    const double*           v3sigma2lapl,
                                    const double*           v3sigma2tau,
                                    const double*           v3sigmalapl2,
                                    const double*           v3sigmalapltau,
                                    const double*           v3sigmatau2,
                                    const double*           v3lapl3,
                                    const double*           v3lapl2tau,
                                    const double*           v3lapltau2,
                                    const double*           v3tau3,
                                    const CDensityGridQuad& rwDensityGridQuad,
                                    const CDensityGrid&     rw2DensityGrid,
                                    const int               iFock,
                                    CMultiTimer&            timer) -> CDenseMatrix;

/**
 Integrates fourth-order meta-GGA exchange-correlation functional
 contribution to AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityMatrix the perturbed one-time transformed densities.
 @param rw2DensityMatrix the two-time transformed densities.
 @param rw3DensityMatrix the three-time transformed densities.
 @param gsDensityMatrix the ground state density matrix.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 @param cubeMode a string that specifies which densities should be combined.
 */
auto integrateKxcLxcFockForMGGA(const std::vector<double*>& aoFockPointers,
                                const CMolecule&            molecule,
                                const CMolecularBasis&      basis,
                                const CAODensityMatrix&     rwDensityMatrix,
                                const CAODensityMatrix&     rw2DensityMatrix,
                                const CAODensityMatrix&     rw3DensityMatrix,
                                const CAODensityMatrix&     gsDensityMatrix,
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
 @param rwDensityGridQuad the products of one-time transformed densities on grid points.
 @param rw2DensityMatrix the two-time transformed densities on grid points.
 @param iFock the index of the AO Fock matrix.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialKxcFockForMGGA2(const CXCFunctional&     xcFunctional,
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
                                     const CDensityGridCubic& rwDensityGridCubic,
                                     const CDensityGrid&      rw2DensityGrid,
                                     const int                iFock,
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
 @param rw3DensityMatrix the three-time transformed densities on grid points.
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

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

#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
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

}  // namespace xcintmgga

#endif /* XCIntegratorForMGGA_hpp */

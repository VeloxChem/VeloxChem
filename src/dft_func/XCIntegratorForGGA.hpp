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

#ifndef XCIntegratorForGGA_hpp
#define XCIntegratorForGGA_hpp

#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcintgga {  // xcintgga namespace

/**
 Integrates first-order GGA exchange-correlation functional contribution to
 AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @param flag the flag for closed/open shell.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForGGA(const CMolecule&                  molecule,
                            const CMolecularBasis&            basis,
                            const std::vector<const double*>& gsDensityPointers,
                            const CMolecularGrid&             molecularGrid,
                            const double                      screeningThresholdForGTOValues,
                            const CXCFunctional&              xcFunctional,
                            const std::string&                flag = std::string("closedshell")) -> CAOKohnShamMatrix;

/**
 Integrates GGA contribution to (first-order) Vxc matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the gradient density.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialVxcFockForGGA(const double*       weights,
                                   const CDenseMatrix& gtoValues,
                                   const CDenseMatrix& gtoValuesX,
                                   const CDenseMatrix& gtoValuesY,
                                   const CDenseMatrix& gtoValuesZ,
                                   const double*       rhograd,
                                   const double*       vrho,
                                   const double*       vsigma,
                                   CMultiTimer&        timer) -> CDenseMatrix;

/**
 Integrates GGA contribution to AO Kohn-Sham matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhograd the gradient density.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param timer the timer.
 @return the alpha and beta contribution as a list of CDenseMatrix objects.
 */
auto integratePartialVxcFockForGGAOpenShell(const double*       weights,
                                            const CDenseMatrix& gtoValues,
                                            const CDenseMatrix& gtoValuesX,
                                            const CDenseMatrix& gtoValuesY,
                                            const CDenseMatrix& gtoValuesZ,
                                            const double*       rhograd,
                                            const double*       vrho,
                                            const double*       vsigma,
                                            CMultiTimer&        timer) -> std::vector<CDenseMatrix>;

/**
 Integrates second-order GGA exchange-correlation functional contribution
 to AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param rwDensityPointers the pointers to perturbed density matrices.
 @param gsDensityPointers the pointers to ground state density matrices.
 @param molecularGrid the molecular grid.
 @param screeningThresholdForGTOValues the screening threshold for GTO values.
 @param xcFunctional the exchange-correlation functional.
 */
auto integrateFxcFockForGGA(const std::vector<double*>&       aoFockPointers,
                            const CMolecule&                  molecule,
                            const CMolecularBasis&            basis,
                            const std::vector<const double*>& rwDensityPointers,
                            const std::vector<const double*>& gsDensityPointers,
                            const CMolecularGrid&             molecularGrid,
                            const double                      screeningThresholdForGTOValues,
                            const CXCFunctional&              xcFunctional) -> void;

/**
 Integrates GGA contribution to (second-order) Fxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param rhow the pointer to perturbed density.
 @param rhograd the pointer to density gradient.
 @param rhowgrad the pointer to perturbed density gradient.
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param v2rhosigma the 2nd-order functional derivative wrt density and
        density gradient.
 @param v2sigma2 the 2nd-order functional derivative wrt density gradient.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialFxcFockForGGA(const CXCFunctional& xcFunctional,
                                   const double*        weights,
                                   const CDenseMatrix&  gtoValues,
                                   const CDenseMatrix&  gtoValuesX,
                                   const CDenseMatrix&  gtoValuesY,
                                   const CDenseMatrix&  gtoValuesZ,
                                   const double*        rhow,
                                   const double*        rhograd,
                                   const double*        rhowgrad,
                                   const double*        vsigma,
                                   const double*        v2rho2,
                                   const double*        v2rhosigma,
                                   const double*        v2sigma2,
                                   CMultiTimer&         timer) -> CDenseMatrix;

}  // namespace xcintgga

#endif /* XCIntegratorForGGA_hpp */

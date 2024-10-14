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

#ifndef XCIntegratorForLDA_hpp
#define XCIntegratorForLDA_hpp

#include <string>

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridQuad.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

namespace xcintlda {  // xcintlda namespace

/**
 Integrates first-order LDA exchange-correlation functional contribution to
 AO Kohn-Sham matrix.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param gsDensityPointers the pointers to AO density matrices.
 @param molecularGrid the molecular grid.
 @param xcFunctional the exchange-correlation functional.
 @param flag the flag for closed/open shell.
 @return the AO Kohn-Sham matrix.
 */
auto integrateVxcFockForLDA(const CMolecule&                  molecule,
                            const CMolecularBasis&            basis,
                            const std::vector<const double*>& gsDensityPointers,
                            const CMolecularGrid&             molecularGrid,
                            const double                      screeningThresholdForGTOValues,
                            const CXCFunctional&              xcFunctional,
                            const std::string&                flag = std::string("closedshell")) -> CAOKohnShamMatrix;

/**
 Integrates LDA contribution to (first-order) Vxc matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param vrho the 1st-order functional derivative wrt density.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialVxcFockForLDA(const double* weights, const CDenseMatrix& gtoValues, const double* vrho, CMultiTimer& timer) -> CDenseMatrix;

/**
 Integrates LDA contribution to (first-order) Vxc matrix.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param vrho the 1st-order functional derivative wrt density.
 @param timer the timer.
 @return the alpha and beta contribution as a list of CDenseMatrix objects.
 */
auto integratePartialVxcFockForLDAOpenShell(const double*       weights,
                                            const CDenseMatrix& gtoValues,
                                            const double*       vrho,
                                            CMultiTimer&        timer) -> std::vector<CDenseMatrix>;

/**
 Integrates second-order LDA exchange-correlation functional contribution
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
auto integrateFxcFockForLDA(const std::vector<double*>&       aoFockPointers,
                            const CMolecule&                  molecule,
                            const CMolecularBasis&            basis,
                            const std::vector<const double*>& rwDensityPointers,
                            const std::vector<const double*>& gsDensityPointers,
                            const CMolecularGrid&             molecularGrid,
                            const double                      screeningThresholdForGTOValues,
                            const CXCFunctional&              xcFunctional) -> void;

/**
 Integrates LDA contribution to (second-order) Fxc matrix.

 @param xcFunctional the exchange-correlation functional.
 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param rhow the pointer to perturbed density.
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param timer the timer.
 @return the contribution as a CDenseMatrix object.
 */
auto integratePartialFxcFockForLDA(const CXCFunctional& xcFunctional,
                                   const double*        weights,
                                   const CDenseMatrix&  gtoValues,
                                   const double*        rhow,
                                   const double*        v2rho2,
                                   CMultiTimer&         timer) -> CDenseMatrix;

auto
integrateKxcFockForLDA(const std::vector<double*>& aoFockPointers,
                       const CMolecule&        molecule,
                       const CMolecularBasis&  basis,
                       const CAODensityMatrix& rwDensityMatrix,
                       const CAODensityMatrix& rw2DensityMatrix,
                       const CAODensityMatrix& gsDensityMatrix,
                       const CMolecularGrid&   molecularGrid,
                       const double            screeningThresholdForGTOValues,
                       const CXCFunctional&    xcFunctional,
                       const std::string&      quadMode) -> void;

auto
integratePartialKxcFockForLDA(const CXCFunctional&    xcFunctional,
                              const double*           weights,
                              const CDenseMatrix&     gtoValues,
                              const double*           v2rho2,
                              const double*           v3rho3,
                              const CDensityGridQuad& rwDensityGridQuad,
                              const CDensityGrid&     rw2DensityGrid,
                              const int           iFock,
                              CMultiTimer&            timer) -> CDenseMatrix;

}  // namespace xcintlda

#endif /* XCIntegratorForLDA_hpp */

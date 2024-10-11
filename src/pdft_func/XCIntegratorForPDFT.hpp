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

#ifndef XCIntegratorForPDFT_hpp
#define XCIntegratorForPDFT_hpp

#include <array>
#include <string>

#include "AOKohnShamMatrix.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCPairDensityFunctional.hpp"

namespace xcintpdft {  // xcintpdft namespace

/**
 Integrates first-order LDA pair-density functional contribution to
 AO Kohn-Sham matrix and MO "Q-matrix".

 @param aoFockMatrix the AO Fock matrix.
 @param tensorWxc the MO Two-body energy gradient term.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param densityMatrixPointer the pointer to AO density matrix.
 @param twoBodyDensityMatrix the MO two-body active density matrix.
 @param activeMOs the active molecular orbitals.
 @param molecularGrid the molecular grid.
 @param fvxc the exchange-correlation functional.
 */
void integrateVxcPDFTForLDA(CAOKohnShamMatrix&              aoFockMatrix,
                            CDense4DTensor&                 tensorWxc,
                            const CMolecule&                molecule,
                            const CMolecularBasis&          basis,
                            const double*                   densityMatrixPointer,
                            const CDenseMatrix&             twoBodyDensityMatrix,
                            const CDenseMatrix&             activeMOs,
                            const CMolecularGrid&           molecularGrid,
                            const double                    screeningThresholdForGTOValues,
                            const CXCPairDensityFunctional& xcFunctional);

/**
 Integrates first-order GGA pair-density functional contribution to
 AO Kohn-Sham matrix and MO "Q-matrix".

 @param aoFockMatrix the AO Fock matrix.
 @param tensorWxc the MO Two-body energy gradient term.
 @param molecule the molecule.
 @param basis the molecular basis.
 @param densityMatrixPointer the pointer to AO density matrix.
 @param twoBodyDensityMatrix the MO two-body active density matrix.
 @param activeMOs the active molecular orbitals.
 @param molecularGrid the molecular grid.
 @param fvxc the exchange-correlation functional.
 */
void integrateVxcPDFTForGGA(CAOKohnShamMatrix&              aoFockMatrix,
                            CDense4DTensor&                 tensorWxc,
                            const CMolecule&                molecule,
                            const CMolecularBasis&          basis,
                            const double*                   densityMatrixPointer,
                            const CDenseMatrix&             twoBodyDensityMatrix,
                            const CDenseMatrix&             activeMOs,
                            const CMolecularGrid&           molecularGrid,
                            const double                    screeningThresholdForGTOValues,
                            const CXCPairDensityFunctional& xcFunctional);

/**
 Integrates PLDA contribution to (first-order) Wxc tensor.

 @param weights the weights of grid points.
 @param gtoValues the GTO values on grid points.
 @param activeMOs the active molecular orbitals.
 @param vrho the 1st-order functional derivative wrt density.
 @param timer the timer.
 @return the contribution as a CDense4DTensor object.
 */
CDenseMatrix integratePartialWxcFockForPLDA(const double*       weights,
                                            const CDenseMatrix& gtoValues,
                                            const CDenseMatrix& activeMOs,
                                            const double*       vrho,
                                            CMultiTimer&        timer);

}  // namespace xcintpdft

#endif /* XCIntegratorForPDFT_hpp */

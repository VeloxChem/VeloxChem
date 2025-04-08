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
 @param rs_omega the range-separation parameter.
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
                            const CXCPairDensityFunctional& xcFunctional,
                            const double                    rs_omega);

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
 @param rs_omega the range-separation parameter.
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
                            const CXCPairDensityFunctional& xcFunctional,
                            const double                    rs_omega);

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

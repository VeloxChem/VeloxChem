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

#ifndef XCMolecularGradientForPDFT_hpp
#define XCMolecularGradientForPDFT_hpp

#include "XCPairDensityFunctional.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"

namespace xcgradpdft {  // xcgradpdft namespace

/**
 Integrates first-order LDA exchange-correlation functional contribution to
 PDFT molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param densityMatrix the AO density matrix object.
 @param twoBodyDensityMatrix the MO two-body active density matrix.
 @param activeMOs the active molecular orbitals.
 @param molecularGrid the molecular grid.
 @param fvxc exchange-correlation functional.
 @param rs_omega range-separation parameter.
 @return the molecular gradient.
 */
 CDenseMatrix integrateVxcPDFTGradientForLDA(const CMolecule&                molecule,
                                             const CMolecularBasis&          basis,
                                             const double*                   densityMatrixPointer,
                                             const CDenseMatrix&             twoBodyDensityMatrix,
                                             const CDenseMatrix&             activeMOs,
                                             const CMolecularGrid&           molecularGrid,
                                             const double                    screeningThresholdForGTOValues,
                                             const CXCPairDensityFunctional& xcFunctional,
                                             const double                    rs_omega);
/**
 Integrates first-order GGA exchange-correlation functional contribution to
 PDFT molecular gradient.

 @param molecule the molecule.
 @param basis the molecular basis.
 @param densityMatrix the AO density matrix object.
 @param twoBodyDensityMatrix the MO two-body active density matrix.
 @param activeMOs the active molecular orbitals.
 @param molecularGrid the molecular grid.
 @param fvxc exchange-correlation functional.
 @param rs_omega range-separation parameter.
 @return the molecular gradient.
 */
 CDenseMatrix integrateVxcPDFTGradientForGGA(const CMolecule&                molecule,
                                             const CMolecularBasis&          basis,
                                             const double*                   densityMatrixPointer,
                                             const CDenseMatrix&             twoBodyDensityMatrix,
                                             const CDenseMatrix&             activeMOs,
                                             const CMolecularGrid&           molecularGrid,
                                             const double                    screeningThresholdForGTOValues,
                                             const CXCPairDensityFunctional& xcFunctional,
                                             const double                    rs_omega);

}   // namespace xcgradpdft

#endif /* XCMolecularGradientForPDFT_hpp */

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

#ifndef PairDensityGridGenerator_hpp
#define PairDensityGridGenerator_hpp

#include "DenseMatrix.hpp"
#include "MultiTimer.hpp"

namespace pairdengridgen {  // pairdengridgen namespace

/**
 Generates density and on-top pair-density for LDA.

 @param rho the pointer to density.
 @param gtoValues the GTO values on grid points.
 @param densityMatrix the total density matrix.
 @param activeMOs the MO coefficients.
 @param twoBodyDensityMatrix the MO two-body density matrix.
 */
void serialGeneratePairDensityForLDA(double*               rho,
                                     const CDenseMatrix&   gtoValues,
                                     const CDenseMatrix&   densityMatrix,
                                     const CDenseMatrix&   activeMOs,
                                     const CDenseMatrix&   twoBodyDensityMatrix);

/**
 Generates density and on-top pair-density for GGA.

 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO X derivative values on grid points.
 @param gtoValuesY the GTO Y derivative values on grid points.
 @param gtoValuesZ the GTO Z derivative values on grid points.
 @param densityMatrix the total density matrix.
 @param activeMOs the MO coefficients.
 @param twoBodyDensityMatrix the MO two-body density matrix.
 @param timer the timer.
 */

void generatePairDensityForGGA(double*               rho,
                               double*               rhograd,
                               double*               sigma,
                               const CDenseMatrix&   gtoValues,
                               const CDenseMatrix&   gtoValuesX,
                               const CDenseMatrix&   gtoValuesY,
                               const CDenseMatrix&   gtoValuesZ,
                               const CDenseMatrix&   densityMatrix,
                               const CDenseMatrix&   activeMOs,
                               const CDenseMatrix&   twoBodyDensityMatrix,
                               CMultiTimer&          timer);

}  // namespace pairdengridgen

#endif /* PairDensityGridGenerator_hpp */

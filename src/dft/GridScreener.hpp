//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#ifndef GridScreener_hpp
#define GridScreener_hpp

#include <cstdint>
#include <vector>

#include "DensityGrid.hpp"

namespace gridscreen {  // gridscreen namespace

/**
 Screens Vxc Fock for LDA.

 @param rho the density.
 @param exc the functional value.
 @param vrho the 1st-order functional derivative wrt density.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenVxcFockForLDA(double*               rho,
                         double*               exc,
                         double*               vrho,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Fxc Fock for LDA.

 @param rho the density.
 @param v2rho2 the 2nd-order functional derivative wrt density.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenFxcFockForLDA(double*               rho,
                         double*               v2rho2,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Kxc Fock for LDA.

 @param rho the density.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v3rho3 the 3rd-order functional derivative wrt rho.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenKxcFockForLDA(double*               rho,
                         double*               v2rho2,
                         double*               v3rho3,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Lxc Fock for LDA.

 @param rho the density.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v3rho3 the 3rd-order functional derivative wrt rho.
 @param v4rho4 the 4rd-order functional derivative wrt rho.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenLxcFockForLDA(double*               rho,
                         double*               v2rho2,
                         double*               v3rho3,
                         double*               v4rho4,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Vxc Fock for PLDA.

 @param rho the density.
 @param exc the functional value.
 @param vrho the 1st-order functional derivative wrt density.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenVxcFockForPLDA(double*               rho,
                          double*               exc,
                          double*               vrho,
                          const int32_t         npoints,
                          const double          densityThreshold);

/**
 Screens Vxc Fock for GGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param exc the functional value.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenVxcFockForGGA(double*               rho,
                         double*               sigma,
                         double*               exc,
                         double*               vrho,
                         double*               vsigma,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Fxc Fock for GGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenFxcFockForGGA(double*               rho,
                         double*               sigma,
                         double*               vrho,
                         double*               vsigma,
                         double*               v2rho2,
                         double*               v2rhosigma,
                         double*               v2sigma2,
                         const int32_t         npoints,
                         const double          densityThreshold);


/**
 Screens Fxc Fock for MGGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma
 
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenFxcFockForMGGA(double*                rho, 
                          double*                sigma, 
                          double*                lapl, 
                          double*                tau, 
                          double*                v2rho2, 
                          double*                v2rhosigma,
                          double*                v2rholapl,
                          double*                v2rhotau,
                          double*                v2sigma2,
                          double*                v2sigmalapl,
                          double*                v2sigmatau,
                          double*                v2lapl2,
                          double*                v2lapltau,
                          double*                v2tau2,
                          const int32_t         npoints,
                          const double          densityThreshold);  

/**
 Screens Kxc Fock for MGGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma
 
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenKxcFockForMGGA(double*               rho, 
                          double*               sigma, 
                          double*               lapl, 
                          double*               tau,
                          double*               v3rho3,
                          double*               v3rho2sigma,
                          double*               v3rho2lapl,
                          double*               v3rho2tau,
                          double*               v3rhosigma2,
                          double*               v3rhosigmalapl,
                          double*               v3rhosigmatau,
                          double*               v3rholapl2,
                          double*               v3rholapltau,
                          double*               v3rhotau2,
                          double*               v3sigma3,
                          double*               v3sigma2lapl,
                          double*               v3sigma2tau,
                          double*               v3sigmalapl2,
                          double*               v3sigmalapltau,
                          double*               v3sigmatau2,
                          double*               v3lapl3,
                          double*               v3lapl2tau,
                          double*               v3lapltau2,
                          double*               v3tau3,      
                          const int32_t         npoints,
                          const double          densityThreshold);

/**
 Screens Lxc Fock for MGGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma
 
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */

void screenLxcFockForMGGA(double*             rho, 
                          double*             sigma,
                          double*             lapl, 
                          double*             tau,
                          double*             v4rho4,
                          double*             v4rho3sigma,
                          double*             v4rho3lapl,
                          double*             v4rho3tau,
                          double*             v4rho2sigma2,
                          double*             v4rho2sigmalapl,
                          double*             v4rho2sigmatau,
                          double*             v4rho2lapl2,
                          double*             v4rho2lapltau,
                          double*             v4rho2tau2,
                          double*             v4rhosigma3,
                          double*             v4rhosigma2lapl,
                          double*             v4rhosigma2tau,
                          double*             v4rhosigmalapl2,
                          double*             v4rhosigmalapltau,
                          double*             v4rhosigmatau2,
                          double*             v4rholapl3,
                          double*             v4rholapl2tau,
                          double*             v4rholapltau2,
                          double*             v4rhotau3,
                          double*             v4sigma4,
                          double*             v4sigma3lapl,
                          double*             v4sigma3tau,
                          double*             v4sigma2lapl2,
                          double*             v4sigma2lapltau,
                          double*             v4sigma2tau2,
                          double*             v4sigmalapl3,
                          double*             v4sigmalapl2tau,
                          double*             v4sigmalapltau2,
                          double*             v4sigmatau3,
                          double*             v4lapl4,
                          double*             v4lapl3tau,
                          double*             v4lapl2tau2,
                          double*             v4lapltau3,
                          double*             v4tau4,
                          const int32_t         npoints,
                          const double          densityThreshold);




/**
 Screens Kxc Fock for GGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma.
 @param v3rho2 the 3rd-order functional derivative wrt rho.
 @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
 @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
 @param v3sigma2 the 3rd-order functional derivative wrt sigma.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenKxcFockForGGA(double*               rho,
                         double*               sigma,
                         double*               vrho,
                         double*               vsigma,
                         double*               v2rho2,
                         double*               v2rhosigma,
                         double*               v2sigma2,
                         double*               v3rho3,
                         double*               v3rho2sigma,
                         double*               v3rhosigma2,
                         double*               v3sigma3,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens Lxc Fock for GGA.

 @param rho the density.
 @param sigma the dot product of density gradient.
 @param vrho the 1st-order functional derivative wrt rho.
 @param vsigma the 1st-order functional derivative wrt sigma.
 @param v2rho2 the 2nd-order functional derivative wrt rho.
 @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
 @param v2sigma2 the 2nd-order functional derivative wrt sigma.
 @param v3rho2 the 3rd-order functional derivative wrt rho.
 @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
 @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
 @param v3sigma2 the 3rd-order functional derivative wrt sigma.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 */
void screenLxcFockForGGA(double*               rho,
                         double*               sigma,
                         double*               vrho,
                         double*               vsigma,
                         double*               v2rho2,
                         double*               v2rhosigma,
                         double*               v2sigma2,
                         double*               v3rho3,
                         double*               v3rho2sigma,
                         double*               v3rhosigma2,
                         double*               v3sigma3,
                         double*               v4rho4,
                         double*               v4rho3sigma,
                         double*               v4rho2sigma2,
                         double*               v4rhosigma3,
                         double*               v4sigma4,
                         const int32_t         npoints,
                         const double          densityThreshold);

/**
 Screens LDA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param rho the pointer to density.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 @return the number of grid points after screening.
 */
int32_t screenDensityForLDA(std::vector<int32_t>& gridPointInds,
                            double*               rho,
                            const int32_t         npoints,
                            const double          densityThreshold);

/**
 Screens GGA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param rho the pointer to density.
 @param rhograd the pointer to density gradient.
 @param sigma the pointer to dot product of density gradient.
 @param npoints the number of grid points.
 @param densityThreshold the threshold for density grid screening.
 @return the number of grid points after screening.
 */
int32_t screenDensityForGGA(std::vector<int32_t>& gridPointInds,
                            double*               rho,
                            double*               rhograd,
                            double*               sigma,
                            const int32_t         npoints,
                            const double          densityThreshold);

/**
 Screens LDA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param destDensityGrid the density grid to store the result (destination).
 @param srcDensityGrid the density grid to provide the grid points (source).
 @param densityThreshold the threshold for density grid screening.
 */
void screenDensityGridForLDA(std::vector<int32_t>& gridPointInds,
                             CDensityGrid&         destDensityGrid,
                             const CDensityGrid&   srcDensityGrid,
                             const double          densityThreshold);

/**
 Screens GGA density grid to get rid of invalid grid points.

 @param gridPointInds mapping between grid points before and after screening
 @param destDensityGrid the density grid to store the result (destination).
 @param srcDensityGrid the density grid to provide the grid points (source).
 @param densityThreshold the threshold for density grid screening.
 */
void screenDensityGridForGGA(std::vector<int32_t>& gridPointInds,
                             CDensityGrid&         destDensityGrid,
                             const CDensityGrid&   srcDensityGrid,
                             const double          densityThreshold);

/**
 Copies weights.

 @param screenedWeights pointer to the screened weights.
 @param gridBlockPosition the starting position of the grid box.
 @param weights pointer to the original weights of all grid points.
 @param nScreenedGridPoints the number of grid points after screening.
 */
void
copyWeights(double*                     screenedWeights,
            const int32_t               gridBlockPosition,
            const double*               weights,
            const int32_t               nGridPoints);

/**
 Screens weights.

 @param screenedWeights pointer to the screened weights.
 @param gridBlockPosition the starting position of the grid box.
 @param weights pointer to the original weights of all grid points.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenWeights(double*                     screenedWeights,
                   const int32_t               gridBlockPosition,
                   const double*               weights,
                   const std::vector<int32_t>& gridPointInds,
                   const int32_t               nScreenedGridPoints);

/**
 Screens GTO matrix for LDA.

 @param screenedGtoValues the matrix containing screened GTO values.
 @param originalGtoValues the matrix containing original GTO values.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenGtoMatrixForLDA(CDenseMatrix&               screenedGtoValues,
                           const CDenseMatrix&         originalGtoValues,
                           const std::vector<int32_t>& gridPointInds,
                           const int32_t               nScreenedGridPoints);

/**
 Screens GTO matrix for GGA.

 @param screenedGtoValues the matrix containing screened GTO values.
 @param screenedGtoValuesX the matrix containing screened GTO X derivatives.
 @param screenedGtoValuesY the matrix containing screened GTO Y derivatives.
 @param screenedGtoValuesZ the matrix containing screened GTO Z derivatives.
 @param originalGtoValues the matrix containing original GTO values.
 @param originalGtoValuesX the matrix containing original GTO X derivatives.
 @param originalGtoValuesY the matrix containing original GTO Y derivatives.
 @param originalGtoValuesZ the matrix containing original GTO Z derivatives.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenGtoMatrixForGGA(CDenseMatrix&               screenedGtoValues,
                           CDenseMatrix&               screenedGtoValuesX,
                           CDenseMatrix&               screenedGtoValuesY,
                           CDenseMatrix&               screenedGtoValuesZ,
                           const CDenseMatrix&         originalGtoValues,
                           const CDenseMatrix&         originalGtoValuesX,
                           const CDenseMatrix&         originalGtoValuesY,
                           const CDenseMatrix&         originalGtoValuesZ,
                           const std::vector<int32_t>& gridPointInds,
                           const int32_t               nScreenedGridPoints);

/**
 Screens GTO matrix for meta-GGA.

 @param screenedGtoValues the matrix containing screened GTO values.
 @param screenedGtoValuesX the matrix containing screened GTO X derivatives.
 @param screenedGtoValuesY the matrix containing screened GTO Y derivatives.
 @param screenedGtoValuesZ the matrix containing screened GTO Z derivatives.
 @param screenedGtoValuesXX the matrix containing screened GTO XX derivatives.
 @param screenedGtoValuesXY the matrix containing screened GTO XY derivatives.
 @param screenedGtoValuesXZ the matrix containing screened GTO XZ derivatives.
 @param screenedGtoValuesYY the matrix containing screened GTO YY derivatives.
 @param screenedGtoValuesYZ the matrix containing screened GTO YZ derivatives.
 @param screenedGtoValuesZZ the matrix containing screened GTO ZZ derivatives.
 @param originalGtoValues the matrix containing original GTO values.
 @param originalGtoValuesX the matrix containing original GTO X derivatives.
 @param originalGtoValuesY the matrix containing original GTO Y derivatives.
 @param originalGtoValuesZ the matrix containing original GTO Z derivatives.
 @param originalGtoValuesXX the matrix containing original GTO XX derivatives.
 @param originalGtoValuesXY the matrix containing original GTO XY derivatives.
 @param originalGtoValuesXZ the matrix containing original GTO XZ derivatives.
 @param originalGtoValuesYY the matrix containing original GTO YY derivatives.
 @param originalGtoValuesYZ the matrix containing original GTO YZ derivatives.
 @param originalGtoValuesZZ the matrix containing original GTO ZZ derivatives.
 @param gridPointInds mapping between grid points before and after screening
 @param nScreenedGridPoints the number of grid points after screening.
 */
void screenGtoMatrixForMetaGGA(CDenseMatrix&               screenedGtoValues,
                               CDenseMatrix&               screenedGtoValuesX,
                               CDenseMatrix&               screenedGtoValuesY,
                               CDenseMatrix&               screenedGtoValuesZ,
                               CDenseMatrix&               screenedGtoValuesXX,
                               CDenseMatrix&               screenedGtoValuesXY,
                               CDenseMatrix&               screenedGtoValuesXZ,
                               CDenseMatrix&               screenedGtoValuesYY,
                               CDenseMatrix&               screenedGtoValuesYZ,
                               CDenseMatrix&               screenedGtoValuesZZ,
                               const CDenseMatrix&         originalGtoValues,
                               const CDenseMatrix&         originalGtoValuesX,
                               const CDenseMatrix&         originalGtoValuesY,
                               const CDenseMatrix&         originalGtoValuesZ,
                               const CDenseMatrix&         originalGtoValuesXX,
                               const CDenseMatrix&         originalGtoValuesXY,
                               const CDenseMatrix&         originalGtoValuesXZ,
                               const CDenseMatrix&         originalGtoValuesYY,
                               const CDenseMatrix&         originalGtoValuesYZ,
                               const CDenseMatrix&         originalGtoValuesZZ,
                               const std::vector<int32_t>& gridPointInds,
                               const int32_t               nScreenedGridPoints);

}  // namespace gridscreen

#endif /* GridScreener_hpp */

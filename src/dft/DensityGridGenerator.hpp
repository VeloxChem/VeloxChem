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

#ifndef DensityGridGenerator_hpp
#define DensityGridGenerator_hpp

#include "AODensityMatrix.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "MultiTimer.hpp"
#include "XCFuncType.hpp"

namespace dengridgen {  // dengridgen namespace

/**
 Generates density grid for LDA.

 @param npoints the number of grid points.
 @param gtoValues the GTO values on grid points.
 @param densityMatrix the density matrix.
 @param xcFunType the type of exchange-correlation functional.
 @param timer the timer.
 @return the density grid.
 */
CDensityGrid generateDensityGridForLDA(const int32_t       npoints,
                                       const CDenseMatrix& gtoValues,
                                       const CDenseMatrix& densityMatrix,
                                       const xcfun         xcFunType,
                                       CMultiTimer&        timer);

CDensityGrid generateDensityGridForLDA(const int32_t           npoints,
                                       const CDenseMatrix&     gtoValues,
                                       const CAODensityMatrix& densityMatrix,
                                       const xcfun             xcFunType,
                                       CMultiTimer&            timer);

/**
 Generates density grid for GGA.

 @param npoints the number of grid points.
 @param gtoValues the GTO values on grid points.
 @param gtoValuesX the GTO gradient X values on grid points.
 @param gtoValuesY the GTO gradient Y values on grid points.
 @param gtoValuesZ the GTO gradient Z values on grid points.
 @param densityMatrix the density matrix.
 @param xcFunType the type of exchange-correlation functional.
 @param timer the timer.
 @return the density grid.
 */
CDensityGrid generateDensityGridForGGA(const int32_t       npoints,
                                       const CDenseMatrix& gtoValues,
                                       const CDenseMatrix& gtoValuesX,
                                       const CDenseMatrix& gtoValuesY,
                                       const CDenseMatrix& gtoValuesZ,
                                       const CDenseMatrix& densityMatrix,
                                       const xcfun         xcFunType,
                                       CMultiTimer&        timer);

}  // namespace dengridgen

#endif /* DensityGridGenerator_hpp */

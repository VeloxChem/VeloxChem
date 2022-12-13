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

#ifndef DensityGridCubic_hpp
#define DensityGridCubic_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "DensityGrid.hpp"
#include "DensityGridType.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

/**
 Class CDensityGridCubic class implements density grid.

 @author K. Ahmadzadeh
 */
class CDensityGridCubic
{
    /**
     The type of density grid.
     */
    dengrid _gridType;

    /**
     The number of associated density matrices.
     */
    int32_t _nDensityMatrices;

    /**
     The density variables values at grid points.
     */
    CMemBlock2D<double> _densityValues;

   public:
    /**
     Creates an empty density grid object.
     */
    CDensityGridCubic();

    /**
     Creates a density grid object.

     @param densityValues the 2D memory block object with density values data.
     @param gridType the type of density grid.
     */
    CDensityGridCubic(const CMemBlock2D<double>& densityValues, const dengrid gridType);

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGridCubic(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

    /**
     Creates a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic(const CDensityGridCubic& source);

    /**
     Creates a density grid object by by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic(CDensityGridCubic&& source) noexcept;

    /**
     Destroys a density grid object.
     */
    ~CDensityGridCubic();

    /**
     Assigns a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic& operator=(const CDensityGridCubic& source);

    /**
     Assigns a density grid object by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGridCubic& operator=(CDensityGridCubic&& source) noexcept;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are equal, false otherwise.
     */
    bool operator==(const CDensityGridCubic& other) const;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are not equal, false otherwise.
     */
    bool operator!=(const CDensityGridCubic& other) const;

    /**
     Initialize density values at grid point to zero.
     */
    void zero();

    /**
     Gets number of grid points in density grid object.

     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;

    /**
     Gets number of densitu matrices associated with density grid.

     @return the number of density matrices.
     */
    int32_t getNumberOfDensityMatrices() const;

    /**
     Gets type of density grid type object.

     @return the type of density grid.
     */
    dengrid getDensityGridType() const;

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrD(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrD(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDY(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDXX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDXX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDXY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDXY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDXZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDXZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDYX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDYX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDYY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDYY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDYZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDYZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDZX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDZX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDZY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDZY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* rBrCrDZZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* rBrCrDZZ(const int32_t iDensityMatrix);


    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDXZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDXZZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDYZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDYZZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* rBrCrDZZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* rBrCrDZZZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product with X gradient of  one-time transformed density.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    const double* RhoBCRhoD(const int32_t iDensityMatrix) const;

    /**
     Gets  pointer to one-time transformed density product with X gradient of  one-time transformed density.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    double* RhoBCRhoD(const int32_t iDensityMatrix);


    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDZZ(const int32_t iDensityMatrix);



    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* RhoBCRhoDZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* RhoBCRhoDZ(const int32_t iDensityMatrix);

    /**
     Generates products of one-time transformed densities to be used for quadratic response.

     @param densityGridab the screened density grid.
     @param molecularGridab the screened molecular grid.
     @param rwDensityGrid  one-time transformed densities evaluated on the grid
     @param densityThreshold the screening threshold for density values.
     @param numdens the number of densities.
     @param quadMode a string that specifies which densities should be combined.
     */
    void DensityProd(const CDensityGrid& rwDensityGrid,
                     const CDensityGrid& rw2DensityGrid,
                     const xcfun         xcFuncType,
                     const int32_t             numdens,
                     const std::string&  quadMode);

    /**
     Converts density grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGridCubic& source);
};

#endif /* DensityGridQuad_hpp */

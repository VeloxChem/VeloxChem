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
    const double* pi(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* pi(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* gam2(const int32_t iDensityMatrix);


    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* gam2(const int32_t iDensityMatrix) const;

/**
         Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        double* gam2X(const int32_t iDensityMatrix);


        /**
         Gets constant pointer to one-time transformed density product values.

         @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        const double* gam2X(const int32_t iDensityMatrix) const;
        
        /**
         Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        double* gam2Y(const int32_t iDensityMatrix);


        /**
         Gets constant pointer to one-time transformed density product values.

         @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        const double* gam2Y(const int32_t iDensityMatrix) const;
        
        /**
         Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        double* gam2Z(const int32_t iDensityMatrix);

        /**
         Gets constant pointer to one-time transformed density product values.

         @param iDensityMatrix the index of density matrix.
         @return the pointer to alpha density values.
         */
        const double* gam2Z(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2XX(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2XX(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2XY(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2XY(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2XZ(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2XZ(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2YX(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2YX(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2YY(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2YY(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2YZ(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2YZ(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2ZX(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2ZX(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2ZY(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2ZY(const int32_t iDensityMatrix) const;
        
        /**
            Gets constant pointer to one-time transformed density product values.

        @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        double* gam2ZZ(const int32_t iDensityMatrix);


        /**
            Gets constant pointer to one-time transformed density product values.

            @param iDensityMatrix the index of density matrix.
            @return the pointer to alpha density values.
            */
        const double* gam2ZZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piY(const int32_t iDensityMatrix);


    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piXX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piXX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piXY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piXY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piXZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piXZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piYX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piYX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piYY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piYY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piYZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piYZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piZX(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piZX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piZY(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piZY(const int32_t iDensityMatrix);

        /**
     Gets constant pointer to one-time transformed density product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* piZZ(const int32_t iDensityMatrix) const;

    /**
     Gets constant pointer to one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* piZZ(const int32_t iDensityMatrix);


    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piXZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piXZZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piYZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piYZZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    const double* piZZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values. 
    */
    double* piZZZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to one-time transformed density product with X gradient of  one-time transformed density.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    const double* gam(const int32_t iDensityMatrix) const;

    /**
     Gets  pointer to one-time transformed density product with X gradient of  one-time transformed density.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    double* gam(const int32_t iDensityMatrix);


    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamXX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamXX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamXY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamXY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamXZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamXZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamYX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamYX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamYY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamYY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamYZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamYZ(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamZX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamZX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamZY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamZY(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamZZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamZZ(const int32_t iDensityMatrix);



    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamX(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamX(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamY(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamY(const int32_t iDensityMatrix);


    double prod3_r(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i);
    
    double prod3_i(double B_r, double B_i, double C_r, double C_i, double D_r, double D_i);

    double prod2_r(double B_r, double B_i, double C_r, double C_i);

    double prod2_i(double B_r, double B_i, double C_r, double C_i);

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    const double* gamZ(const int32_t iDensityMatrix) const;

    /**
    Gets constant pointer to products of one-time and two-time transformed density product values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to alpha density values.
    */
    double* gamZ(const int32_t iDensityMatrix);

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

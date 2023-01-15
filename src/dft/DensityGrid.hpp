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

#ifndef DensityGrid_hpp
#define DensityGrid_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "DensityGridType.hpp"
#include "MemBlock2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"

/**
 Class CDensityGrid class implements density grid.

 @author Z. Rinkevicius
 */
class CDensityGrid
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
    CDensityGrid();

    /**
     Creates a density grid object.

     @param densityValues the 2D memory block object with density values data.
     @param gridType the type of density grid.
     */
    CDensityGrid(const CMemBlock2D<double>& densityValues, const dengrid gridType);

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGrid(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

    /**
     Creates a density grid object.
     
     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param nComponents the number of components per single AO density matrix.
     @param gridType the type of density grid.
     */
    CDensityGrid(const int32_t nGridPoints,
                 const int32_t nDensityMatrices,
                 const int32_t nComponents,
                 const dengrid gridType);
    
    /**
     Creates a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGrid(const CDensityGrid& source);

    /**
     Creates a density grid object by by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGrid(CDensityGrid&& source) noexcept;

    /**
     Destroys a density grid object.
     */
    ~CDensityGrid();

    /**
     Assigns a density grid object by copying other density grid object.

     @param source the density grid object.
     */
    CDensityGrid& operator=(const CDensityGrid& source);

    /**
     Assigns a density grid object by moving other density grid object.

     @param source the density grid object.
     */
    CDensityGrid& operator=(CDensityGrid&& source) noexcept;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are equal, false otherwise.
     */
    bool operator==(const CDensityGrid& other) const;

    /**
     Compares density grid object with other density grid object.

     @param other the density grid object.
     @return true if density grid objects are not equal, false otherwise.
     */
    bool operator!=(const CDensityGrid& other) const;

    /**
     Initialize density values at grid point to zero.
     */
    void zero();

    /**
     Reduces number of grid points by slicing off all grid points beoynd specific number of grid points.

     @param nGridPoints the number of grid points.
     */
    void slice(const int32_t nGridPoints);

    /**
     Updates beta densites (density, it's gradient components) by assigning alpha density values.
     */
    void updateBetaDensities();

    /**
     Computes density norms.
     */
    void computeDensityNorms();

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
     Gets constant pointer to alpha density values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    const double* alphaDensity(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* alphaDensity(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to beta density values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density values.
     */
    const double* betaDensity(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density values.
     */
    double* betaDensity(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient norm values.
     */
    const double* alphaDensityGradient(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient norm values.
     */
    double* alphaDensityGradient(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    const double* betaDensityGradient(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    double* betaDensityGradient(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha * beta density gradients product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha * beta density gradients product values.
     */
    const double* mixedDensityGradient(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha * beta density gradients product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha * beta density gradients product values.
     */
    double* mixedDensityGradient(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    const double* alphaDensityGradientX(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    double* alphaDensityGradientX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Y values.
     */
    const double* alphaDensityGradientY(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Y values.
     */
    double* alphaDensityGradientY(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Z values.
     */
    const double* alphaDensityGradientZ(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Z values.
     */
    double* alphaDensityGradientZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient X values.
     */
    const double* betaDensityGradientX(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient X values.
     */
    double* betaDensityGradientX(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Y values.
     */
    const double* betaDensityGradientY(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Y values.
     */
    double* betaDensityGradientY(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Z values.
     */
    const double* betaDensityGradientZ(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Z values.
     */
    double* betaDensityGradientZ(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to alpha density tau norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    const double* alphaDensitytau(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density tau norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    double* alphaDensitytau(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to beta density tau norm values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to beta density gradient norm values.
    */
    const double* betaDensitytau(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density tau norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    double* betaDensitytau(const int32_t iDensityMatrix);


    /**
     Gets constant pointer to alpha density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha laplacian values.
     */
    const double* alphaDensitylapl(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to alpha density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to  alpha laplacian values.
     */
    double* alphaDensitylapl(const int32_t iDensityMatrix);

    /**
    Gets constant pointer to beta density laplacian  values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to beta laplacian values.
    */
    const double* betaDensitylapl(const int32_t iDensityMatrix) const;

    /**
     Gets pointer to beta density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta laplacian values.
     */
    double* betaDensitylapl(const int32_t iDensityMatrix);

    /**
     Gets constant pointer to specific component of density grid.
     
     @param iComponent the index of density grid component.
     @return the constant pointer to density grid component.
     */
    const double* getComponent(const int32_t iComponent) const;
    
    /**
     Gets pointer to specific component of density grid.
     
     @param iComponent the index of density grid component.
     @return the constant pointer to density grid component.
     */
    double* getComponent(const int32_t iComponent);
    
    /**
     Converts density grid object to text and insert it into output text
     stream.

     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGrid& source);
};

#endif /* DensityGrid_hpp */

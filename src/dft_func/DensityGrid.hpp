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

#ifndef DensityGrid_hpp
#define DensityGrid_hpp

#include <cstdint>
#include <vector>

#include "DensityGridData2D.hpp"
#include "MolecularGrid.hpp"
#include "XCFunctionalType.hpp"

enum class dengrid
{
    ab, lima, limb, undefined
};

/**
 Class CDensityGrid class implements density grid.
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
    int _nDensityMatrices;

    /**
     The density variables values at grid points.
     */
    CDensityGridData2D _densityValues;

   public:
    /**
     Creates an empty density grid object.
     */
    CDensityGrid();

    /**
     Creates a density grid object.

     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGrid(const int nGridPoints, const int nDensityMatrices, const xcfun xcFuncType, const dengrid gridType);

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
     Gets number of grid points in density grid object.

     @return the number of grid points.
     */
    int getNumberOfGridPoints() const;

    /**
     Gets number of densitu matrices associated with density grid.

     @return the number of density matrices.
     */
    int getNumberOfDensityMatrices() const;

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
    const double* alphaDensity(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density values.
     */
    double* alphaDensity(const int iDensityMatrix);

    /**
     Gets constant pointer to beta density values.

    @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density values.
     */
    const double* betaDensity(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density values.
     */
    double* betaDensity(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient norm values.
     */
    const double* alphaDensityGradient(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient norm values.
     */
    double* alphaDensityGradient(const int iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    const double* betaDensityGradient(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient norm values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient norm values.
     */
    double* betaDensityGradient(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha * beta density gradients product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha * beta density gradients product values.
     */
    const double* mixedDensityGradient(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha * beta density gradients product values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha * beta density gradients product values.
     */
    double* mixedDensityGradient(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    const double* alphaDensityGradientX(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient X values.
     */
    double* alphaDensityGradientX(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Y values.
     */
    const double* alphaDensityGradientY(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Y values.
     */
    double* alphaDensityGradientY(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Z values.
     */
    const double* alphaDensityGradientZ(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha density gradient Z values.
     */
    double* alphaDensityGradientZ(const int iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient X values.
     */
    const double* betaDensityGradientX(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient X values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient X values.
     */
    double* betaDensityGradientX(const int iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Y values.
     */
    const double* betaDensityGradientY(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient Y values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Y values.
     */
    double* betaDensityGradientY(const int iDensityMatrix);

    /**
     Gets constant pointer to beta density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Z values.
     */
    const double* betaDensityGradientZ(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density gradient Z values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta density gradient Z values.
     */
    double* betaDensityGradientZ(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density tau values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha tau values.
     */
    const double* alphaDensitytau(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density tau values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha tau values.
     */
    double* alphaDensitytau(const int iDensityMatrix);

    /**
    Gets constant pointer to beta density tau values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to beta tau values.
    */
    const double* betaDensitytau(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density tauvalues.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta tau values.
     */
    double* betaDensitytau(const int iDensityMatrix);

    /**
     Gets constant pointer to alpha density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to alpha laplacian values.
     */
    const double* alphaDensitylapl(const int iDensityMatrix) const;

    /**
     Gets pointer to alpha density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to  alpha laplacian values.
     */
    double* alphaDensitylapl(const int iDensityMatrix);

    /**
    Gets constant pointer to beta density laplacian  values.

    @param iDensityMatrix the index of density matrix.
    @return the pointer to beta laplacian values.
    */
    const double* betaDensitylapl(const int iDensityMatrix) const;

    /**
     Gets pointer to beta density laplacian values.

     @param iDensityMatrix the index of density matrix.
     @return the pointer to beta laplacian values.
     */
    double* betaDensitylapl(const int iDensityMatrix);

    /**
     Gets constant pointer to specific component of density grid.
     
     @param iComponent the index of density grid component.
     @return the constant pointer to density grid component.
     */
    const double* getComponent(const int iComponent) const;
    
    /**
     Gets pointer to specific component of density grid.
     
     @param iComponent the index of density grid component.
     @return the constant pointer to density grid component.
     */
    double* getComponent(const int iComponent);
};

#endif /* DensityGrid_hpp */

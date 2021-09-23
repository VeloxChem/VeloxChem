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

#include "MemBlock2D.hpp"
#include "DensityGridType.hpp"
#include "XCFuncType.hpp"
#include "MolecularGrid.hpp"

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
    
    /**
     Checks if density grid point is within range of allowe density values for LDA case.

     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForLda(const double alphaDensity,
                                 const double betaDensity,
                                 const double densityThreshold) const;
    
    
    
    /**
     Checks if density grid point is within range of allowe density values for LDA case.

     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForLdaUnrestricted(const double alphaDensity,
                                       const double betaDensity,
                                       const double densityThreshold) const;

   /**
     Checks if density grid point is within range of allowe density values for LDA case.
     Only valid if beta density is greater than the threshhold but not the alpha density.

     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForLdaA(const double alphaDensity,
                                  const double betaDensity,
                                  const double densityThreshold) const;

    /**
     Checks if density grid point is within range of allowe density values for LDA case.
     Only valid if alpha density is greater than the threshhold but not the beta density.

     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForLdaB(const double alphaDensity,
                                 const double betaDensity,
                                 const double densityThreshold) const;


    /**
     Checks if density grid point is within range of allowe density values for GGA case.
     
     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param alphaDensityGradient the value of alpha density gradient at grid point.
     @param betaDensityGradient the value of beta density gradient at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForGga(const double alphaDensity,
                                 const double betaDensity,
                                 const double alphaDensityGradient,
                                 const double betaDensityGradient,
                                 const double densityThreshold) const;


    /**
     Checks if density grid point is within range of allowe density values for GGA case.
     
     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param alphaDensityGradient the value of alpha density gradient at grid point.
     @param betaDensityGradient the value of beta density gradient at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForGgaUnrestrictedAB(const double alphaDensity,
                                               const double betaDensity,
                                               const double alphaDensityGradient,
                                               const double betaDensityGradient,
                                               const double densityThreshold) const;


    /**
     Checks if density grid point is within range of allowe density values for GGA case.
     
     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param alphaDensityGradient the value of alpha density gradient at grid point.
     @param betaDensityGradient the value of beta density gradient at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForGgaUnrestrictedA(const double alphaDensity,
                                              const double betaDensity,
                                              const double alphaDensityGradient,
                                              const double betaDensityGradient,
                                              const double densityThreshold) const;


    /**
     Checks if density grid point is within range of allowe density values for GGA case.
     
     @param alphaDensity the value of alpha density at grid point.
     @param betaDensity the value of beta density at grid point.
     @param alphaDensityGradient the value of alpha density gradient at grid point.
     @param betaDensityGradient the value of beta density gradient at grid point.
     @param densityThreshold the threshold for density variable.
     @return true if valid density grid point, false otherwise.
     */
    bool _isValidGridPointForGgaUnrestrictedB(const double alphaDensity,
                                              const double betaDensity,
                                              const double alphaDensityGradient,
                                              const double betaDensityGradient,
                                              const double densityThreshold) const;
    
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
    CDensityGrid(const CMemBlock2D<double>& densityValues,
                 const dengrid              gridType);
    
    /**
     Creates a density grid object.
     
     @param nGridPoints the number of grid points with density values.
     @param nDensityMatrices the number of AO density matrices.
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGrid(const int32_t nGridPoints,
                 const int32_t nDensityMatrices,
                 const xcfun   xcFuncType,
                 const dengrid gridType);
    
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
     Generates pair of screened molecular and density grids by removing grid points with specific density/density
     gradient values bellow given threshold. NOTE: This method is exclusive to dengrid::ab type.
     
     @param densityGridab the screened density grid.
     @param molecularGridab the screened molecular grid.
     @param iDensityMatrix the index of density matrix.
     @param densityThreshold the screening threshold for density values.
     @param xcFuncType the type of exchange-correlation functional.
     */
    void getScreenedGridsPair(      CDensityGrid&   densityGridab,
                                    CMolecularGrid& molecularGridab,
                              const int32_t         iDensityMatrix,
                              const double          densityThreshold,
                              const xcfun           xcFuncType) const;


    /**
     Generates pair of screened molecular and density grids by removing grid points with specific density/density
     gradient values bellow given threshold. Extends to other types of densities than dengrid::ab

     @param densityGridab the screened density grid, a and b type.
     @param densityGrida the screened density grid lima .
     @param densityGridb the screened density grid limb.
     @param molecularGridab the screened molecular grid.
     @param iDensityMatrix the index of density matrix.
     @param densityThreshold the screening threshold for density values.
     @param xcFuncType the type of exchange-correlation functional.
     */
    void getScreenedGridPairUnrestricted(       CDensityGrid&   densityGridab,
                                                CDensityGrid&   densityGrida,
                                                CDensityGrid&   densityGridb,
                                                CMolecularGrid& molecularGridab,
                                                CMolecularGrid& molecularGrida,
                                                CMolecularGrid& molecularGridb,
                                          const int32_t         iDensityMatrix,
                                          const double          densityThreshold,
                                          const xcfun           xcFuncType) const;
    
    /**
     Generates screened molecular grid from given molecular grid by removing grid points with specific density/density
     gradient values bellow given threshold. NOTE: This method is exclusive to dengrid::ab type.

     @param molecularGrids the molecular grid.
     @param densityThreshold the density/density gradient screening threshold.
     @param iDensityMatrix the index of density matrix.
     @param xcFuncType the type of exchange-correlation functional.
     @return the screened molecular grid.
     */
    CMolecularGrid getScreenedGrid(      CMolecularGrid& molecularGrids,
                                   const int32_t         iDensityMatrix,
                                   const double          densityThreshold,
                                   const xcfun           xcFuncType) const;

    /**
    Generates screened molecular grid from given molecular grid by removing grid points with specific density/density
    gradient values bellow given threshold.

    @param molecularGridsAB the molecular grid for ab.
    @param molecularGridsA the molecular grid for lima.
    @param molecularGridsB the molecular grid for limb.
    @param densityThreshold the density/density gradient screening threshold.
    @param iDensityMatrix the index of density matrix.
    @param xcFuncType the type of exchange-correlation functional.
    @return the screened molecular grid.
    */
    void getScreenedGridUnrestricted(      CMolecularGrid& molecularGridsAB,
                                           CMolecularGrid& molecularGridsA,
                                           CMolecularGrid& molecularGridsB,
                                     const int32_t         iDensityMatrix,
                                     const double          densityThreshold,
                                     const xcfun           xcFuncType) const;
    
    /**
     Converts density grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGrid& source);
};

#endif /* DensityGrid_hpp */

//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
     Generates pair of screened molecular and density grids by removing grid points with specific density/density
     gradient values bellow given threshold. NOTE: This method is exclusive to dengrid::ab type.
     
     @param densityGrid the screened density grid.
     @param molecularGrid the screened molecular grid.
     @param iDensityMatrix the index of density matrix.
     @param densityThreshold the screening threshold for density values.
     @param xcFuncType the type of exchange-correlation functional.
     */
    void getScreenedGridsPair(      CDensityGrid&   densityGrid,
                                    CMolecularGrid& molecularGrid,
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
     Converts density grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGrid& source);
};

#endif /* DensityGrid_hpp */

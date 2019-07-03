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

#include "MemBlock2D.hpp"
#include "DensityGridType.hpp"
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
     @param xcFuncType the type of exchange-correlation functional.
     @param gridType the type of density grid.
     */
    CDensityGrid(const int32_t nGridPoints, const xcfun xcFuncType, const dengrid gridType);
    
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
     Gets number of grid points in density grid object.
     
     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;
    
    /**
     Gets constant pointer to alpha density values.
     
     @return the pointer to alpha density values .
     */
    const double* alphaDensity() const;
    
    /**
     Gets pointer to alpha density values.
     
     @return the pointer to alpha density values .
     */
    double* alphaDensity();
    
    /**
     Gets constant pointer to beta density values.
     
     @return the pointer to beta density values .
     */
    const double* betaDensity() const;
    
    /**
     Gets pointer to beta density values.
     
     @return the pointer to beta density values .
     */
    double* betaDensity();
    
    /**
     Gets constant pointer to alpha density gradient norm values.
     
     @return the pointer to alpha density gradient norm values .
     */
    const double* alphaDensityGradient() const;
    
    /**
     Gets pointer to alpha density gradient norm values.
     
     @return the pointer to alpha density gradient norm values .
     */
    double* alphaDensityGradient();
    
    /**
     Gets constant pointer to beta density gradient norm values.
     
     @return the pointer to beta density gradient norm values .
     */
    const double* betaDensityGradient() const;
    
    /**
     Gets pointer to beta density gradient norm values.
     
     @return the pointer to beta density gradient norm values .
     */
    double* betaDensityGradient();
    
    /**
     Gets constant pointer to alpha * beta density gradients product values.
     
     @return the pointer to alpha * beta density gradients product values .
     */
    const double* mixedDensityGradient() const;
    
    /**
     Gets pointer to alpha * beta density gradients product values.
     
     @return the pointer to alpha * beta density gradients product  values .
     */
    double* mixedDensityGradient();
    
    /**
     Converts density grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the density grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CDensityGrid& source);
};


#endif /* DensityGrid_hpp */

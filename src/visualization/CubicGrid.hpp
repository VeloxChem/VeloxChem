//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CubicGrid_hpp
#define CubicGrid_hpp

#include <vector>
#include <cstdint>

#include "MemBlock.hpp"

/**
 Class CCubicGrid contains cubic grid points.
 
 @author X. Li
 */
class CCubicGrid
{
    /**
     Origin of the cubic grid points.
     */
    std::vector<double> _origin;

    /**
     Step size in three dimensions.
     */
    std::vector<double> _stepSize;

    /**
     Number of points in three dimensions.
     */
    std::vector<int32_t> _numPoints;

    /**
     Values at the grid points.
     */
    CMemBlock<double> _values;

public:

    /**
     Creates an empty cubic grid object.
     */
    CCubicGrid();
    
    /**
     Creates a cubic grid object from origin, step size, and
     number of points.
     */
    CCubicGrid(const std::vector<double>&  origin,
               const std::vector<double>&  stepSize,
               const std::vector<int32_t>& numPoints);
    
    /**
     Creates a cubic grid object by copying other cubic grid object.
     
     @param source the cubic grid object.
     */
    CCubicGrid(const CCubicGrid& source);
    
    /**
     Creates a cubic grid object by moving other cubic grid object.
     
     @param source the cubic grid object.
     */
    CCubicGrid(CCubicGrid&& source) noexcept;
    
    /**
     Destroys a cubic grid object.
     */
    ~CCubicGrid();
    
    /**
     Assigns a cubic grid object by copying other cubic grid object.
     
     @param source the cubic grid object.
     */
    CCubicGrid& operator=(const CCubicGrid& source);
    
    /**
     Assigns a cubic grid object by moving other cubic grid object.
     
     @param source the cubic grid object.
     */
    CCubicGrid& operator=(CCubicGrid&& source) noexcept;
    
    /**
     Compares cubic grid object with other cubic grid object.
     
     @param other the cubic grid object.
     @return true if cubic grid objects are equal, false otherwise.
     */
    bool operator==(const CCubicGrid& other) const;
    
    /**
     Compares cubic grid object with other cubic grid object.
     
     @param other the cubic grid object.
     @return true if cubic grid objects are not equal, false otherwise.
     */
    bool operator!=(const CCubicGrid& other) const;
    
    /**
     Gets X coordinate of the origin.
     
     @return X coordinate of the origin.
     */
    double originX() const;
    
    /**
     Gets Y coordinate of the origin.
     
     @return Y coordinate of the origin.
     */
    double originY() const;
    
    /**
     Gets Z coordinate of the origin.
     
     @return Z coordinate of the origin.
     */
    double originZ() const;
    
    /**
     Gets step size in X direction.
     
     @return step size in X direction.
     */
    double stepSizeX() const;

    /**
     Gets step size in Y direction.
     
     @return step size in Y direction.
     */
    double stepSizeY() const;

    /**
     Gets step size in Z direction.
     
     @return step size in Z direction.
     */
    double stepSizeZ() const;

    /**
     Gets number of points in X direction.
     
     @return number of points in X direction.
     */
    int32_t numPointsX() const;

    /**
     Gets number of points in Y direction.
     
     @return number of points in Y direction.
     */
    int32_t numPointsY() const;

    /**
     Gets number of points in Z direction.
     
     @return number of points in Z direction.
     */
    int32_t numPointsZ() const;
     
    /**
     Gets constant pointer to first element of cubic grid values.
     
     @return the constant pointer to first element of cubic grid values.
     */
    const double* values() const;
     
    /**
     Gets pointer to first element of cubic grid values.
     
     @return the pointer to first element of cubic grid values.
     */
    double* values();
};

#endif /* CubicGrid_hpp */

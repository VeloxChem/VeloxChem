//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef XCGradientGrid_hpp
#define XCGradientGrid_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock2D.hpp"
#include "DensityGridType.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"
#include "DensityGrid.hpp"

/**
 Class CXCGradientGrid class implements exchange-correlation grid.
 
 @author Z. Rinkevicius
 */
class CXCGradientGrid
{
    /**
     The type of desnity grid associated with exchange-correlation gradient grid.
     */
    dengrid _densityGridType;
    
    /**
    The type of exchange-correlation gradient grid.
    */
    xcfun _xcGridType;
    
    /**
     The exchange-correlation functional and it's first derrivative values at grid points.
     */
    CMemBlock2D<double> _xcValues;
    
public:
    /**
     Creates an empty exchange-correlation gradient grid object.
     */
    CXCGradientGrid();
    
    /**
     Creates a exchange-correlation gradient grid object.
     
     @param xcValues the 2D memory block object with exchange-correlation functional and it's first derrivatives
            values data.
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation gradient grid.
     */
    CXCGradientGrid(const CMemBlock2D<double>& xcValues,
                    const dengrid              densityGridType,
                    const xcfun                xcGridType);
    
    /**
     Creates a exchange-correlation gradient grid object.
     
     @param nGridPoints the number of grid points in exchange-correlation gradient grid. 
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation gradient grid.
     */
    CXCGradientGrid(const int32_t nGridPoints,
                    const dengrid densityGridType,
                    const xcfun   xcGridType);

    /**
     Creates a exchange-correlation gradient grid object by copying other exchange-correlation gradient grid object.
     
     @param source the exchange-correlation gradient grid object.
     */
    CXCGradientGrid(const CXCGradientGrid& source);
    
    /**
     Creates a exchange-correlation gradient grid object by moving other exchange-correlation gradient grid object.
     
     @param source the exchange-correlation gradient grid object.
     */
    CXCGradientGrid(CXCGradientGrid&& source) noexcept;
    
    /**
     Destroys a exchange-correlation gradient grid object.
     */
    ~CXCGradientGrid();
    
    /**
     Assigns a exchange-correlation gradient grid object by copying other exchange-correlation gradient grid object.
     
     @param source the exchange-correlation gradient grid object.
     */
    CXCGradientGrid& operator=(const CXCGradientGrid& source);
    
    /**
     Assigns a exchange-correlation gradient grid object by moving other exchange-correlation gradient grid object.
     
     @param source the exchange-correlation gradient grid object.
     */
    CXCGradientGrid& operator=(CXCGradientGrid&& source) noexcept;
    
    /**
     Compares exchange-correlation gradient grid object with other exchange-correlation gradient grid object.
     
     @param other the exchange-correlation gradient grid object.
     @return true if exchange-correlation gradient grid objects are equal, false otherwise.
     */
    bool operator==(const CXCGradientGrid& other) const;
    
    /**
     Compares exchange-correlation gradient grid object with other exchange-correlation gradient grid object.
     
     @param other the exchange-correlation gradient grid object.
     @return true if exchange-correlation gradient grid objects are not equal, false otherwise.
     */
    bool operator!=(const CXCGradientGrid& other) const;
    
    /**
     Initialize density values at grid point to zero.
     */
    void zero();
    
    /**
     Gets number of grid points in exchange-correlation gradient grid object.
     
     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;
    
    /**
     Gets constant pointer to exchange-correlation functional values.
     
     @return the constant pointer to exchage-correlation functional values.
     */
    const double* xcFunctionalValues() const;
    
    /**
     Gets pointer to exchange-correlation functional values.
     
     @return the pointer to exchage-correlation functional values.
     */
    double* xcFunctionalValues();
    
    /**
     Gets constant pointer to exchange-correlation functional gradient values.
     
     @param gradientType the type of exchange-correlation functional gradient.
     @return the constant pointer to exchage-correlation functional gradient values.
     */
    const double* xcGradientValues(const xcvars gradientType) const;
    
    /**
    Gets pointer to exchange-correlation functional gradient values.
    
    @param gradientType the type of exchange-correlation functional gradient.
    @return the pointer to exchage-correlation functional gradient values.
    */
    double* xcGradientValues(const xcvars gradientType);
    
    /**
     Converts exchange-correlation gradient grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the exchange-correlation gradient grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CXCGradientGrid& source);
};

#endif /* XCGradientGrid_hpp */

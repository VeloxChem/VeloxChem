//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef XCCubicHessianGrid_hpp
#define XCCubicHessianGrid_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock2D.hpp"
#include "DensityGridType.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"
#include "DensityGrid.hpp"

/**
 Class CXCCubicHessianGrid class implements exchange-correlation functional hessian grid.
 
 @author Z. Rinkevicius
 */
class CXCCubicHessianGrid
{
    /**
     The type of desnity grid associated with exchange-correlation hessian grid.
     */
    dengrid _densityGridType;
    
    /**
     The type of exchange-correlation hessian grid.
     */
    xcfun _xcGridType;
    
    /**
     The exchange-correlation functional hessian values at grid points.
     */
    CMemBlock2D<double> _xcValues;
    
public:
    /**
     Creates an empty exchange-correlation hessian grid object.
     */
    CXCCubicHessianGrid();
    
    /**
     Creates a exchange-correlation hessian grid object.
     
     @param xcValues the 2D memory block object with exchange-correlation functional hessian values data.
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation hessian grid.
     */
    CXCCubicHessianGrid(const CMemBlock2D<double>& xcValues,
                        const dengrid              densityGridType,
                        const xcfun                xcGridType);
    
    /**
     Creates a exchange-correlation hessian grid object.
     
     @param nGridPoints the number of grid points in exchange-correlation hessian grid.
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation hessian grid.
     */
    CXCCubicHessianGrid(const int32_t nGridPoints,
                        const dengrid densityGridType,
                        const xcfun   xcGridType);
    
    /**
     Creates a exchange-correlation hessian grid object by copying other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCCubicHessianGrid(const CXCCubicHessianGrid& source);
    
    /**
     Creates a exchange-correlation hessian grid object by moving other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCCubicHessianGrid(CXCCubicHessianGrid&& source) noexcept;
    
    /**
     Destroys a exchange-correlation hessian grid object.
     */
    ~CXCCubicHessianGrid();
    
    /**
     Assigns a exchange-correlation hessian grid object by copying other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCCubicHessianGrid& operator=(const CXCCubicHessianGrid& source);
    
    /**
     Assigns a exchange-correlation hessian grid object by moving other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCCubicHessianGrid& operator=(CXCCubicHessianGrid&& source) noexcept;
    
    /**
     Compares exchange-correlation hessian grid object with other exchange-correlation hessian grid object.
     
     @param other the exchange-correlation hessian grid object.
     @return true if exchange-correlation hessian grid objects are equal, false otherwise.
     */
    bool operator==(const CXCCubicHessianGrid& other) const;
    
    /**
     Compares exchange-correlation hessian grid object with other exchange-correlation hessian grid object.
     
     @param other the exchange-correlation hessian grid object.
     @return true if exchange-correlation hessian grid objects are not equal, false otherwise.
     */
    bool operator!=(const CXCCubicHessianGrid& other) const;
    
    /**
     Initialize hessian values at grid point to zero.
     */
    void zero();
    
    /**
     Gets number of grid points in exchange-correlation hessian grid object.
     
     @return the number of grid points.
     */
    int32_t getNumberOfGridPoints() const;
    
    /**
     Gets constant pointer to exchange-correlation functional hessian values.
     
     @param iVariableType the type of first differentiation variable in hessian.
     @param jVariableType the type of second differentiation variable in hessian.
     @param jVariableType the type of third differentiation variable in hessian.
     @return the constant pointer to exchage-correlation functional hessian values.
     */
    const double* xcCubicHessianValues(const xcvars iVariableType,
                                       const xcvars jVariableType,
                                       const xcvars kVariableType) const;
    
    /**
     Gets pointer to exchange-correlation functional hessian values.
     
     @param iVariableType the type of first differentiation variable in hessian.
     @param jVariableType the type of second differentiation variable in hessian.
     @param jVariableType the type of third differentiation variable in hessian.
     @return the pointer to exchage-correlation functional hessian values.
     */
    double* xcCubicHessianValues(const xcvars iVariableType,
                                 const xcvars jVariableType,
                                 const xcvars kVariableType);
    
    /**
     Converts exchange-correlation hessian grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the exchange-correlation hessian grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CXCCubicHessianGrid& source);
};

#endif /* XCHessianGrid_hpp */

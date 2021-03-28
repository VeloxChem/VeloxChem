//
//                           VELOXCHEM 1.0-RC
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

#ifndef XCHessianGrid_hpp
#define XCHessianGrid_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock2D.hpp"
#include "DensityGridType.hpp"
#include "XCFuncType.hpp"
#include "XCVarsType.hpp"
#include "DensityGrid.hpp"

/**
 Class CXCHessianGrid class implements exchange-correlation functional hessian grid.
 
 @author Z. Rinkevicius
 */
class CXCHessianGrid
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
    CXCHessianGrid();
    
    /**
     Creates a exchange-correlation hessian grid object.
     
     @param xcValues the 2D memory block object with exchange-correlation functional hessian values data.
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation hessian grid.
     */
    CXCHessianGrid(const CMemBlock2D<double>& xcValues,
                   const dengrid              densityGridType,
                   const xcfun                xcGridType);
    
    /**
     Creates a exchange-correlation hessian grid object.
     
     @param nGridPoints the number of grid points in exchange-correlation hessian grid.
     @param densityGridType the type of density grid.
     @param xcGridType the type of echange-correlation hessian grid.
     */
    CXCHessianGrid(const int32_t nGridPoints,
                   const dengrid densityGridType,
                   const xcfun   xcGridType);
    
    /**
     Creates a exchange-correlation hessian grid object by copying other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCHessianGrid(const CXCHessianGrid& source);
    
    /**
     Creates a exchange-correlation hessian grid object by moving other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCHessianGrid(CXCHessianGrid&& source) noexcept;
    
    /**
     Destroys a exchange-correlation hessian grid object.
     */
    ~CXCHessianGrid();
    
    /**
     Assigns a exchange-correlation hessian grid object by copying other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCHessianGrid& operator=(const CXCHessianGrid& source);
    
    /**
     Assigns a exchange-correlation hessian grid object by moving other exchange-correlation hessian grid object.
     
     @param source the exchange-correlation hessian grid object.
     */
    CXCHessianGrid& operator=(CXCHessianGrid&& source) noexcept;
    
    /**
     Compares exchange-correlation hessian grid object with other exchange-correlation hessian grid object.
     
     @param other the exchange-correlation hessian grid object.
     @return true if exchange-correlation hessian grid objects are equal, false otherwise.
     */
    bool operator==(const CXCHessianGrid& other) const;
    
    /**
     Compares exchange-correlation hessian grid object with other exchange-correlation hessian grid object.
     
     @param other the exchange-correlation hessian grid object.
     @return true if exchange-correlation hessian grid objects are not equal, false otherwise.
     */
    bool operator!=(const CXCHessianGrid& other) const;
    
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
     @return the constant pointer to exchage-correlation functional hessian values.
     */
    const double* xcHessianValues(const xcvars iVariableType,
                                  const xcvars jVariableType) const;
    
    /**
     Gets pointer to exchange-correlation functional hessian values.
     
     @param iVariableType the type of first differentiation variable in hessian.
     @param jVariableType the type of second differentiation variable in hessian.
     @return the pointer to exchage-correlation functional hessian values.
     */
    double* xcHessianValues(const xcvars iVariableType,
                            const xcvars jVariableType);
    
    /**
     Converts exchange-correlation hessian grid object to text and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the exchange-correlation hessian grid.
     */
    friend std::ostream& operator<<(std::ostream& output, const CXCHessianGrid& source);
};

#endif /* XCHessianGrid_hpp */

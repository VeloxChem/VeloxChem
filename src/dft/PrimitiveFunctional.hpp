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

#ifndef PrimitiveFunctional_hpp
#define PrimitiveFunctional_hpp

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>

#include "XCFuncType.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"
#include "DensityGrid.hpp"

using def_vxc_func_typ = void(CXCGradientGrid&, const double factor, const CDensityGrid&);

using def_vxc2_func_typ = void(CXCHessianGrid&, const double factor, const CDensityGrid&);

// using def_vxc2_func_typ = void(CXCGradientGrid&, const CDensityGrid&);

/**
 Class CPrimitiveFunctional stores information about primitive exchange-correlation functional
 and provides methods to perform actions with stored exchange-correlation functional.
 
 @author Z. Rinkevicius
 */
class CPrimitiveFunctional
{
    /**
     The label of primitive exchange-correlation functional.
     */
    std::string _label;
    
    /**
     The type of primitive exchange-correlation functional.
     */
    xcfun _xcFuncType;
    
    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc_func_typ> _abFirstOrderFunction;
    
    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::lima type.
     */
    std::function<def_vxc_func_typ> _aFirstOrderFunction;
    
    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::limb type.
     */
    std::function<def_vxc_func_typ> _bFirstOrderFunction;
    
    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc2_func_typ> _abSecondOrderFunction;
    
    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::lima type.
     */
    std::function<def_vxc2_func_typ> _aSecondOrderFunction;
    
    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::limb type.
     */
    std::function<def_vxc2_func_typ> _bSecondOrderFunction;
    
public:
    /**
     Creates an empty primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional();
    
    /**
     Creates a primitive exchange-correlation functional object.
     
     @param label the label of primitive exchange-correlation functional.
     @param xcFuncType the type of primitive exchange-correlation functional.
     @param abFirstOrderFunction the first-order derivative function (dengrid::ab).
     @param aFirstOrderFunction the first-order derivative function (dengrid::lima).
     @param bFirstOrderFunction the first-order derivative function (dengrid::limb).
     @param abSecondOrderFunction the second-order derivative function (dengrid::ab).
     @param aSecondOrderFunction the second-order derivative function (dengrid::lima).
     @param bSecondOrderFunction the second-order derivative function (dengrid::limb).
     */
    CPrimitiveFunctional(const std::string&                      label,
                         const xcfun                             xcFuncType,
                         const std::function<def_vxc_func_typ>&  abFirstOrderFunction,
                         const std::function<def_vxc_func_typ>&  aFirstOrderFunction,
                         const std::function<def_vxc_func_typ>&  bFirstOrderFunction,
                         const std::function<def_vxc2_func_typ>& abSecondOrderFunction,
                         const std::function<def_vxc2_func_typ>& aSecondOrderFunction,
                         const std::function<def_vxc2_func_typ>& bSecondOrderFunction);
    
    /**
     Creates a primitive exchange-correlation functional object by copying other primitive exchange-correlation functional
     object.
     
     @param source the recursion term object.
     */
    CPrimitiveFunctional(const CPrimitiveFunctional& source);
    
    /**
     Creates a primitive exchange-correlation functional object by moving other primitive exchange-correlation functional
     object.
     
     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional(CPrimitiveFunctional&& source) noexcept;
    
    /**
     Destroys a primitive exchange-correlation functional object.
     */
    ~CPrimitiveFunctional();
    
    /**
     Assigns a primitive exchange-correlation functional object by copying other primitive exchange-correlation functional
     object.
     
     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional& operator=(const CPrimitiveFunctional& source);
    
    /**
     Assigns a primitive exchange-correlation functional object by moving other primitive exchange-correlation functional
     object.
     
     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional& operator=(CPrimitiveFunctional&& source) noexcept;
    
    /**
     Compares primitive exchange-correlation functional object with other primitive exchange-correlation functional object.
     
     @param other the primitive exchange-correlation functional object.
     @return true if primitive exchange-correlation functional objects are equal, false otherwise.
     */
    bool operator==(const CPrimitiveFunctional& other) const;
    
    /**
     Compares primitive exchange-correlation functional object with other primitive exchange-correlation functional object.
     
     @param other the primitive exchange-correlation functional object.
     @return true if primitive exchange-correlation functional objects are not equal, false otherwise.
     */
    bool operator!=(const CPrimitiveFunctional& other) const;
    
    /**
     Gets label of primitive exchange-correlation functional.
     
     @return the label.
     */
    std::string getLabel() const;
    
    /**
     Gets type of primitive exchange-correlation functional.

     @return the type of exchange-correlation functional.
     */
    xcfun getFunctionalType() const;
    
    /**
     Computes first derivative of exchange-correlation functional for given density grid.

     @param xcGradientGrid the exchange-correlation gradient grid object.
     @param factor the scaling factor of functional contribution.
     @param densityGrid the density grid object.
     */
    void compute(      CXCGradientGrid& xcGradientGrid,
                 const double           factor,
                 const CDensityGrid&    densityGrid) const;
    
    /**
     Computes second derivative of exchange-correlation functional for given density grid.
     
     @param xcHessianGrid the exchange-correlation hessian grid object.
     @param factor the scaling factor of functional contribution.
     @param densityGrid the density grid object.
     */
    void compute(      CXCHessianGrid& xcHessianGrid,
                 const double          factor,
                 const CDensityGrid&   densityGrid) const;
    
    /**
     Converts primitive exchange-correlation functional object to text format and insert it into output text stream.
     
     @param output the output text stream.
     @param source the primitive exchange-correlation functional object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CPrimitiveFunctional& source);
};

#endif /* PrimitiveFunctional_hpp */

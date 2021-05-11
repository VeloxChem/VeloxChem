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

#ifndef XCFunctional_hpp
#define XCFunctional_hpp

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>

#include "XCFuncType.hpp"
#include "PrimitiveFunctional.hpp"

/**
 Class CXCFunctional stores information about exchange-correlation functional and provides methods
 to perform actions with stored exchange-correlation functional.
 
 @author Z. Rinkevicius
 */
class CXCFunctional
{
    /**
     The label of exchange-correlation functional.
     */
    std::string _label;
    
    /**
     The type of exchange-correlation functional.
     */
    xcfun _xcFuncType;
    
    
    /**
     The fraction of exact Hatree-Fock exchange in functional.
     */
    double _fractionOfExactExchange;
    
    /**
     The vector of primitive exchange-correlation functionals.
     */
    std::vector<CPrimitiveFunctional> _primitiveFunctionals;
    
    /**
    The vector of weights of primitive exchange-correlation functionals.
    */
    std::vector<double> _weightsOfPrimitiveFunctionals;
    
public:
    /**
     Creates an empty exchange-correlation functional object.
     */
    CXCFunctional();
    
    /**
     Creates a exchange-correlation functional object.
     
     @param label the label of primitive exchange-correlation functional.
     @param xcFuncType the type of primitive exchange-correlation functional.
     @param fractionOfExactExchange the fraction of exact Hatree-Fock exchange.
     @param weightsOfPrimitiveFunctionals the vector of weights of primitive functionals.
     */
    CXCFunctional(const std::string&                       label,
                  const xcfun                              xcFuncType,
                  const double                             fractionOfExactExchange,
                  const std::vector<CPrimitiveFunctional>& primitiveFunctionals,
                  const std::vector<double>&               weightsOfPrimitiveFunctionals);
    
    /**
     Creates a exchange-correlation functional object by copying other exchange-correlation functional object.
     
     @param source the recursion term object.
     */
    CXCFunctional(const CXCFunctional& source);
    
    /**
     Creates a exchange-correlation functional object by moving other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional(CXCFunctional&& source) noexcept;
    
    /**
     Destroys a exchange-correlation functional object.
     */
    ~CXCFunctional();
    
    /**
     Assigns a exchange-correlation functional object by copying other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional& operator=(const CXCFunctional& source);
    
    /**
     Assigns a exchange-correlation functional object by moving other exchange-correlation functional object.
     
     @param source the exchange-correlation functional object.
     */
    CXCFunctional& operator=(CXCFunctional&& source) noexcept;
    
    /**
     Compares exchange-correlation functional object with other exchange-correlation functional object.
     
     @param other the exchange-correlation functional object.
     @return true if exchange-correlation functional objects are equal, false otherwise.
     */
    bool operator==(const CXCFunctional& other) const;
    
    /**
     Compares exchange-correlation functional object with other exchange-correlation functional object.
     
     @param other the exchange-correlation functional object.
     @return true if exchange-correlation functional objects are not equal, false otherwise.
     */
    bool operator!=(const CXCFunctional& other) const;
    
    /**
     Gets label of exchange-correlation functional.
     
     @return the label.
     */
    std::string getLabel() const;
    
    /**
     Gets type of exchange-correlation functional.
     
     @return the type of exchange-correlation functional.
     */
    xcfun getFunctionalType() const;
    
    /**
     Gets fraction of exact Hatree-Fock exchange in exchange-correlation functional.

     @return the fraction of exact Hatree-Fock exchange.
     */
    double getFractionOfExactExchange() const;
    
    /**
     Determines if exchange-correlation functional is of hybrid type i.e. non-zero fraction of exact Hatree-Fock
     exchange.

     @return true if hybrid exchange-correlation functional, false otherwise.
     */
    bool isHybridFunctional() const;
    
    /**
     Determines if exchange-correlation function is undefined.

     @return true if functional is undefined, false otherwise.
     */
    bool isUndefined() const;
    
    /**
     Computes first derivative of exchange-correlation functional for given density grid.
     
     @param xcGradientGrid the exchange-correlation gradient grid object.
     @param densityGrid the density grid object.
     */
    void compute(      CXCGradientGrid& xcGradientGrid,
                 const CDensityGrid&    densityGrid) const;
    
    /**
     Computes second derivative of exchange-correlation functional for given density grid.
     
     @param xcHessianGrid the exchange-correlation hessian grid object.
     @param densityGrid the density grid object.
     */
    void compute(      CXCHessianGrid& xcHessianGrid,
                 const CDensityGrid&   densityGrid) const;
    
    /**
     Converts exchange-correlation functional object to text format and insert it into output text stream.
     
     @param output the output text stream.
     @param source the exchange-correlation functional object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CXCFunctional& source);
};

#endif /* XCFunctional_hpp */

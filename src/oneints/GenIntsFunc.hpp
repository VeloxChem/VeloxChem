//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef GenIntsFunc_hpp
#define GenIntsFunc_hpp

#include "RecursionBlock.hpp"
#include "RecursionTerm.hpp"
#include "RecursionMap.hpp"
#include "RecursionFunctionsList.hpp"

namespace gintsfunc { // gintsfunc namespace
    
    /**
     Generates recursion map object for given recursion term object.

     @param recursionTerm the recursion term object.
     @param angularForm the angular form of recursion map object.
     @param recursionFunctionsList the list of recursion functions.
     @return the recursion map object.
     */
    CRecursionMap genRecursionMap(const CRecursionTerm&          recursionTerm,
                                  const recblock                 angularForm, 
                                  const CRecursionFunctionsList& recursionFunctionsList);
    
    /**
     Generates recursion term corresponding to one-electron integral of given type.

     @param labelOfOperator the label of integrand operator.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @param ordderOfOperator the order of integral.
     @return the recursion term object. 
     */
    CRecursionTerm genIntegral(const std::string& labelOfOperator,
                               const int32_t      braAngularMomentum,
                               const int32_t      ketAngularMomentum,
                               const int32_t      ordderOfOperator);
    
} // gintsfunc namespace

#endif /* GenIntsFunc_hpp */

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
    
} // gintsfunc namespace

#endif /* GenIntsFunc_hpp */

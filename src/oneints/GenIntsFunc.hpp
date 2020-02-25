//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef GenIntsFunc_hpp
#define GenIntsFunc_hpp

#include "RecursionBlock.hpp"
#include "RecursionFunctionsList.hpp"
#include "RecursionMap.hpp"
#include "RecursionTerm.hpp"

namespace gintsfunc {  // gintsfunc namespace

/**
 Generates recursion map object for given recursion term object.

 @param recursionTerm the recursion term object.
 @param angularForm the angular form of recursion map object.
 @param maxRepeatUnits the maximum number of repeated units in recursion
        term subspace.
 @param recursionFunctionsList the list of recursion functions.
 @return the recursion map object.
 */
CRecursionMap genRecursionMap(const CRecursionTerm&          recursionTerm,
                              const recblock                 angularForm,
                              const int32_t                  maxRepeatUnits,
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
    
/**
 Generates recursion term corresponding to two electron repulsion integral.

 @param angularMomentumA the angular momentum of center A.
 @param angularMomentumB the angular momentum of center B.
 @param angularMomentumC the angular momentum of center C.
 @param angularMomentumD the angular momentum of center D.
 @return the recursion term object.
 */
CRecursionTerm genElectronRepulsionIntegral(const int32_t      angularMomentumA,
                                            const int32_t      angularMomentumB,
                                            const int32_t      angularMomentumC,
                                            const int32_t      angularMomentumD);
    

/**
 Generates recursion term corresponding to two electron repulsion integral.

 @param angularMomentumB the angular momentum of center B.
 @param angularMomentumC the angular momentum of center C.
 @param angularMomentumD the angular momentum of center D.
 @return the recursion term object.
 */
CRecursionTerm genElectronRepulsionIntegral(const int32_t      angularMomentumB,
                                            const int32_t      angularMomentumC,
                                            const int32_t      angularMomentumD);

}  // namespace gintsfunc

#endif /* GenIntsFunc_hpp */

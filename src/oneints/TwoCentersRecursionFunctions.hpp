//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef TwoCentersRecursionFunctions_hpp
#define TwoCentersRecursionFunctions_hpp

#include <cstdint>
#include <vector>

#include "RecursionTerm.hpp"

namespace t2crecfunc {  // t2crecfunc namespace

/**
 Applies Obara-Saika overlap recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForOverlap(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika kinetic energy recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForKineticEnergy(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika nuclear potential recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForNuclearPotential(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika electric dipole recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForElectricDipole(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika linear momentum recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForLinearMomentum(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika angular momentum recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForAngularMomentum(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika electric field recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForElectricField(const CRecursionTerm& recursionTerm);

/**
 Applies Obara-Saika electric field gradient recursion to recursion term object.

 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
 */
std::vector<CRecursionTerm> obRecursionForElectricFieldGradient(const CRecursionTerm& recursionTerm);
    
/**
 Applies Obara-Saika vertical electron repulsion recursion to recursion term object.
     
 @param recursionTerm the recursion term object.
 @return the vector of recursion term objects.
*/
std::vector<CRecursionTerm> obRecursionForElectronRepulsion(const CRecursionTerm& recursionTerm);

}  // namespace t2crecfunc

#endif /* TwoCentersRecursionFunctions_hpp */

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

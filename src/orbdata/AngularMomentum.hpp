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

#ifndef AngularMomentum_hpp
#define AngularMomentum_hpp

#include <cstdint>
#include <string>


namespace angmom {  // angmom namespace

/**
 Determines number of spherical components for given angular momentum.

 @param angularMomentum the angular momentum.
 @return the number of spherical components.
*/
int32_t to_SphericalComponents(const int32_t angularMomentum);

/**
 Determines number of spherical components for given pair of angular momentums.

  @param angularMomentumA the first angular momentum.
  @param angularMomentumB the second angular momentum.
  @return the number of spherical components.
*/
int32_t to_SphericalComponents(const int32_t angularMomentumA, const int32_t angularMomentumB);

/**
 Determines number of Cartesian components for given angular momentum.

 @param angularMomentum the angular momentum.
 @return the number of Cartesian momentum.
 */
int32_t to_CartesianComponents(const int32_t angularMomentum);

/**
 Determines number of Cartesian components for given pair of angular momentums.

 @param angularMomentumA the first angular momentum.
 @param angularMomentumB the second angular momentum.
 @return the number of Cartesian momentum.
 */
int32_t to_CartesianComponents(const int32_t angularMomentumA, const int32_t angularMomentumB);

/**
 Gets string representation of spherical angular momentum component.

 @param angularMomentum the angular momentum.
 @param sphericalComponent the spherical component of angular momentum.
 @return the string of angular momentum component.
 */
std::string getStringOfAngularMomentum(const int32_t angularMomentum, const int32_t sphericalComponent);

}  // namespace angmom

#endif /* AngularMomentum_hpp */

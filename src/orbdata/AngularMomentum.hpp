//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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

//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AngularMomentum_hpp
#define AngularMomentum_hpp

#include <cstdint>

namespace angmom { // angmom namespace

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
int32_t to_SphericalComponents(const int32_t angularMomentumA,
                               const int32_t angularMomentumB);
    
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
int32_t to_CartesianComponents(const int32_t angularMomentumA,
                               const int32_t angularMomentumB);

} // angmom namespace

#endif /* AngularMomentum_hpp */

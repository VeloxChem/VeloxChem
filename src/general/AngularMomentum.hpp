//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
 Determines number of Cartesian components for given angular momentum.

 @param angularMomentum the angular momentum.
 @return the number of Cartesian momentum.
 */
int32_t to_CartesianComponents(const int32_t angularMomentum);

} // angmom namespace

#endif /* AngularMomentum_hpp */

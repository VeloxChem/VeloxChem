//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentum.hpp"

namespace angmom { // angmom namespace

int32_t
to_SphericalComponents(const int32_t angularMomentum)
{
    return 2 * angularMomentum + 1;
}
    
int32_t
to_SphericalComponents(const int32_t angularMomentumA,
                       const int32_t angularMomentumB)
{
    return (2 * angularMomentumA + 1) * (2 * angularMomentumB + 1) ;
}

int32_t
to_CartesianComponents(const int32_t angularMomentum)
{
    return (angularMomentum + 1) * (angularMomentum + 2) / 2;
}
    
int32_t
to_CartesianComponents(const int32_t angularMomentumA,
                       const int32_t angularMomentumB)
{
    auto ndim = (angularMomentumA + 1) * (angularMomentumA + 2) / 2;
    
    ndim *= (angularMomentumB + 1) * (angularMomentumB + 2) / 2;
    
    return ndim;
}

} // angmom namespace

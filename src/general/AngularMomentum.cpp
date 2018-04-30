//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AngularMomentum.hpp"

namespace angmom { // angmom namespace

int32_t
to_SphericalComponents(const int32_t angularMomentum)
{
    return 2 * angularMomentum + 1;
}

int32_t
to_CartesianComponents(const int32_t angularMomentum)
{
    return (angularMomentum + 1) * (angularMomentum + 2) / 2;
}

} // angmom namespace

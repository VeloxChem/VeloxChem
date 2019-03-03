//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentum.hpp"

#include <array>

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

std::string
getStringOfAngularMomentum(const int32_t angularMomentum,
                           const int32_t sphericalComponent)
    
{
    if (angularMomentum == 0) return std::string("s  ");
        
    if (angularMomentum == 1)
    {
        std::array<std::string, 3> clist({"p-1", "p0 ", "p+1"});
            
        return clist[sphericalComponent];
    }
        
    if (angularMomentum == 2)
    {
        std::array<std::string, 5> clist({"d-2", "d-1", "d0 ", "d+1", "d+2"});
            
        return clist[sphericalComponent];
    }
        
    if (angularMomentum == 3)
    {
        std::array<std::string, 7> clist({"f-3", "f-2", "f-1", "f0 ", "f+1",
                                          "f+2", "f+3"});
            
        return clist[sphericalComponent];
    }
        
    if (angularMomentum == 4)
    {
        std::array<std::string, 9> clist({"g-4", "g-3", "g-2", "g-1", "g0 ",
                                          "g+1", "g+2", "g+3", "g+4"});
            
        return clist[sphericalComponent];
    }
        
    if (angularMomentum == 5)
    {
        std::array<std::string, 11> clist({"h-5", "h-4", "h-3", "h-2", "h-1",
                                           "h0 ", "h+1", "h+2", "h+3", "h+4",
                                           "h+5"});
            
        return clist[sphericalComponent];
    }
        
    if (angularMomentum == 6)
    {
        std::array<std::string, 13> clist({"i-6", "i-5", "i-4", "i-3", "i-2",
                                           "i-1", "i0 ", "i+1", "i+2", "i+3",
                                           "i+4", "i+5", "i+6"});
            
        return clist[sphericalComponent];
    }
        
    return std::string();
}

} // angmom namespace

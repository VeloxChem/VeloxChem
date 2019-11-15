//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef VdwRadii_hpp
#define VdwRadii_hpp

#include <vector>

namespace vdwradii {  // vdwradii namespace

/**
 Creates VDW radii.

 @return a vector of VDW radii with nuclear charge as index.
 */
std::vector<double> buildVdwRadii();

}  // namespace vdwradii

#endif /* VdwRadii_hpp */

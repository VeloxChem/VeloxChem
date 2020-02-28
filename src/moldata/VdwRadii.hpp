//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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

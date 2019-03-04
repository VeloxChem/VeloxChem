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

#include "Molecule.hpp"

namespace vdwradii { // vdwradii namespace

/**
 Creates VDW radii.
 
 @return a vector with nuclear charge as index.
 */
const std::vector<double> buildVdwRadii();

/**
 Creates VDW radii for a molecule.
 
 @return a vector indexed by atom.
 */
const std::vector<double> getRadii(const CMolecule& mol);

} // vdwradii namespace

#endif /* VdwRadii_hpp */

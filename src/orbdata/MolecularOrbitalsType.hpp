//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularOrbitalsType_hpp
#define MolecularOrbitalsType_hpp

#include <string>

/**
 Enumerate class morb:
 
 Defines supported molecular orbital types:
 morb::rest   - the spin restricted molecular orbitals.
 morb::unrest - the spin unrestricted molecular orbitals.
 */
enum class morb
{
    rest,
    unrest
};

/**
 Converts enumerate class value to it's string label.
 
 @param molecularOrbitals the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const morb molecularOrbitals)
{
    if (molecularOrbitals == morb::rest)
    {
        return std::string("Spin Restricted Molecular Orbitals");
    }
   
    if (molecularOrbitals == morb::unrest)
    {
        return std::string("Spin Unrestricted Molecular Orbitals");
    }
    
    return std::string("UNKNOWN");
}

#endif /* MolecularOrbitalsType_hpp */

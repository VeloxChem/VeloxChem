//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef MolecularOrbitalsType_hpp
#define MolecularOrbitalsType_hpp

#include <string>

/**
 Enumerate class molorb:

 Defines supported molecular orbital types:
 morb::rest   - the spin restricted molecular orbitals.
 morb::unrest - the spin unrestricted molecular orbitals.
 */
enum class molorb
{
    rest,
    unrest
};

/**
 Converts key value of molecular orbitals type to integer number.

 @param motyp the enumerate class value.
 @return the integer number.
 */
inline int32_t
to_int(const molorb motyp)
{
    return static_cast<int32_t>(motyp);
}

/**
 Converts integer key value to molecular orbitals type.

 @param keyValue the integer key value.
 @return the molecular orbital matrix type.
 */
inline molorb
to_molorb(const int32_t keyValue)
{
    if (keyValue == to_int(molorb::rest)) return molorb::rest;

    if (keyValue == to_int(molorb::unrest)) return molorb::unrest;

    return molorb::rest;
}

/**
 Converts enumerate class value to it's string label.

 @param molecularOrbitals the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const molorb molecularOrbitals)
{
    if (molecularOrbitals == molorb::rest)
    {
        return std::string("Spin Restricted Molecular Orbitals");
    }

    if (molecularOrbitals == molorb::unrest)
    {
        return std::string("Spin Unrestricted Molecular Orbitals");
    }

    return std::string("UNKNOWN");
}

#endif /* MolecularOrbitalsType_hpp */

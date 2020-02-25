//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef TwoIntsDistType_hpp
#define TwoIntsDistType_hpp

#include <string>

/**
 Enumerate class dist2e:
 
 Defines supported two electron integrals distribution keys:
 dist2e::batch   - the batch with natural order of data
 dist2e::fock    - the Fock matrix
 dist2e::qvalues - the Q values vectors for bra and ket sides
 */
enum class dist2e
{
    batch,
    fock,
    qvalues 
};

/**
 Converts enumerate class value to it's string label.
 
 @param distPattern the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const dist2e distPattern)
{
    if (distPattern == dist2e::batch)
    {
        return std::string("Raw Integrals Batch");
    }
    
    if (distPattern == dist2e::fock)
    {
        return std::string("Fock Matrix");
    }
    
    if (distPattern == dist2e::qvalues)
    {
        return std::string("Q Values Vector");
    }
    
    return std::string("UNKNOWN");
}

#endif /* TwoIntsDistType_hpp */

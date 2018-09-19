//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TwoIntsDistType_hpp
#define TwoIntsDistType_hpp

#include <string>

/**
 Enumerate class dist2e:
 
 Defines supported two electron integrals distribution keys:
 dist2e::batch   - the batch with natural order of data
 dist2e::rfock   - the restricted Fock matrix for Hatree-Fock method
 dist2e::ufock   - the unrestricted Fock matrix for Hatree-Fock method
 dist2e::qvalues - the Q values vectors for bra and ket sides
 */
enum class dist2e
{
    batch,
    rfock,
    ufock,
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
    
    if (distPattern == dist2e::rfock)
    {
        return std::string("Restricted Fock Matrix");
    }
    
    if (distPattern == dist2e::ufock)
    {
        return std::string("Unrestricted Fock Matrix");
    }
    
    if (distPattern == dist2e::qvalues)
    {
        return std::string("Q Values Vector");
    }
    
    return std::string("UNKNOWN");
}

#endif /* TwoIntsDistType_hpp */

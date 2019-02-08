//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef SpinBlock_hpp
#define SpinBlock_hpp

#include <cstdint>

/**
 Enumerate class szblock:

 Defines all allowed key values for spin block:
 szblock::aa - the alpha-alpha spin block
 szblock::bb - the beta-beta spin block
 szblock::ab - the alpha-beta spin block
 szblock::ba - the beta-alpha spin block
 */

enum class szblock
{
    aa,
    bb,
    ab,
    ba
};

/**
 Converts key value of spin block to integer number.
 
 @param szBlockKey the key value of spin block.
 @return the integer number.
 */
inline int32_t to_int(const szblock szBlockKey)
{
    return static_cast<int32_t>(szBlockKey);
}

/**
 Converts enumerate class value to it's string label.
 
 @param szBlockKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const szblock szBlockKey)
{
    if (szBlockKey == szblock::aa)
    {
        return std::string("Spin Block: Alpha-Alpha");
    }
    
    if (szBlockKey == szblock::bb)
    {
        return std::string("Spin Block: Beta-Beta");
    }
    
    if (szBlockKey == szblock::ab)
    {
        return std::string("Spin Block: Alpha-Beta");
    }
    
    if (szBlockKey == szblock::ba)
    {
        return std::string("Spin Block: Beta-Alpha");
    }
    
    return std::string("UNKNOWN");
}

#endif /* SpinBlock_hpp */

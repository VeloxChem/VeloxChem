//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

#ifndef JobType_hpp
#define JobType_hpp

#include <cstdint>

/**
 Enumerate class job:

 Defines all allowed key values for derrived class of base job class:
 job::sp_energy    - the single point energy calculation
 job::opt_geometry - the geometry optimzation

 */

enum class job : int32_t
{
    sp_energy, opt_geometry
};

/**
 Converts key value for derrived job class to integer number.

 @param jobkey the key value of derrived class.
 @return the integer number.
 */
inline int32_t to_int(const job jobkey)
{
    return static_cast<int32_t>(jobkey);
}

#endif /* JobType_hpp */

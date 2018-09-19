//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ExecMode_hpp
#define ExecMode_hpp

#include <cstdint>

/**
 Enumerate class execmode:

 Defines all allowed key values for jobs execution mode:
 execmode::cpu     - the CPU jobs execution mode
 execmode::cpu_gpu - the hybrid CPU/GPU jobs execution mode
 */

enum class execmode : int32_t
{
    cpu,
    cpu_gpu
};

/**
 Converts key value of job execution mode to integer number.
 
 @param execmodekey the key value of job execution mode.
 @return the integer number.
 */
inline int32_t to_int(const execmode execmodekey)
{
    return static_cast<int32_t>(execmodekey);
}

#endif /* ExecMode_hpp */

//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExecMode_hpp
#define ExecMode_hpp

#include <cstdint>
#include <string>

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

 @param execModeKey the key value of job execution mode.
 @return the integer number.
 */
inline int32_t
to_int(const execmode execModeKey)
{
    return static_cast<int32_t>(execModeKey);
}

/**
 Converts enumerate class value to it's string label.

 @param execModeKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const execmode execModeKey)
{
    if (execModeKey == execmode::cpu)
    {
        return std::string("Execution Mode: CPU");
    }

    if (execModeKey == execmode::cpu_gpu)
    {
        return std::string("Execution Mode: CPU-GPU");
    }

    return std::string("UNKNOWN");
}

#endif /* ExecMode_hpp */

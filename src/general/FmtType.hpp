//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef FmtType_hpp
#define FmtType_hpp

/**
 Enumerate class fmt:

 Defines supported formatting keys:
 fmt::left   - the left alignment of output line
 fmt::center - the center alignment of output line
 fmt::right  - the right alignment of output line
 */
enum class fmt
{
    left,
    center,
    right
};

/**
 Converts enumerate class value to it's string label.
 
 @param formatKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const fmt formatKey)
{
    if (formatKey == fmt::left)
    {
        return std::string("Format Key: Left");
    }
    
    if (formatKey == fmt::center)
    {
        return std::string("Format Key: Center");
    }
    
    if (formatKey == fmt::right)
    {
        return std::string("Format Key: Right");
    }
    
    return std::string("UNKNOWN");
}

#endif /* FmtType_hpp */

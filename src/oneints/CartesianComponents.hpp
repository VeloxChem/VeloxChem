//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef CartesianComponents_hpp
#define CartesianComponents_hpp

/**
 Enumerate class cartesians:

 Defines all allowed key values for Cartesian directions:
 cartesians::X - x direction
 cartesians::Y - Y direction
 cartesians::Z - Z direction
 */

enum class cartesians : int32_t
{
    X = 0,
    Y = 1,
    Z = 2
};

#endif /* CartesianComponents_hpp */

//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef VecIndexes_hpp
#define VecIndexes_hpp

#include <vector>

#include "FourIndexes.hpp"
#include "ThreeIndexes.hpp"
#include "TwoIndexes.hpp"

/**
 Defines alias to STL vector of two indexes objects.
 */
using CVecTwoIndexes = std::vector<CTwoIndexes>;

/**
 Defines alias to STL vector of three indexes objects.
 */
using CVecThreeIndexes = std::vector<CThreeIndexes>;

/**
 Defines alias to STL vector of four indexes objects.
 */
using CVecFourIndexes = std::vector<CFourIndexes>;

#endif /* VecIndexes_hpp */

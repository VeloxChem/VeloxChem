//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef VecIndexes_hpp
#define VecIndexes_hpp

#include <vector>

#include "TwoIndexes.hpp"
#include "ThreeIndexes.hpp"
#include "FourIndexes.hpp"

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

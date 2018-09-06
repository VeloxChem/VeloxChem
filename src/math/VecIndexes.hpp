//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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

//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef VecMemBlocks_hpp
#define VecMemBlocks_hpp

#include <vector>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"

/**
 Defines alias to STL vector of memory block objects.
 */
template <class T>
using CVecMemBlock = std::vector<CMemBlock<T>>;

/**
 Defines alias to STL vector of 2D memory block objects.
 */
template <class T>
using CVecMemBlock2D = std::vector<CMemBlock2D<T>>;

#endif /* VecMemBlocks_hpp */

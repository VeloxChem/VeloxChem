//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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

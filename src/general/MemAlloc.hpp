//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MemAlloc_hpp
#define MemAlloc_hpp

#include <cstdlib>
#include <cstdint>

namespace mem // mem namespace
{

/**
 Allocates memory block aligned to VLX_ALIGN-bit boundary.

 @param size the size of memory block.
 @return the pointer to memory block.
 */
void* malloc(const size_t size);

/**
 Deallocates memory block pointed to by pointer.

 @param pointer the pointer to memory block.
 */
void free(void* pointer);

} // mem namespace

#endif /* MemAlloc_hpp */

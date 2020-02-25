//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef MemAlloc_hpp
#define MemAlloc_hpp

#include <cstdint>
#include <cstdlib>

namespace mem  // mem namespace
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

}  // namespace mem

#endif /* MemAlloc_hpp */

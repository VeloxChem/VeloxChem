//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "MemAlloc.hpp"
#include "ErrorHandler.hpp"

#ifdef ENABLE_MKL
#include "mkl.h"
#else
#include <cstdlib>
#endif

namespace mem { // mem namespace

void*
malloc(const size_t size)
{
    #ifdef ENABLE_MKL

    return MKL_malloc(size, VLX_ALIGN);

    #else

    void* ptr = nullptr;

    int ierr = ::posix_memalign(&ptr, VLX_ALIGN, size);

    errors::assertMsgCritical(ierr == 0, "malloc: posix_memalign failed!");

    return ptr;

    #endif
}

void
free(void* pointer)
{
    #ifdef ENABLE_MKL

    if (pointer != nullptr) MKL_free(pointer);

    #else

    if (pointer != nullptr) ::free(pointer);

    #endif
}

} // mem namespace

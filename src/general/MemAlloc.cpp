//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "MemAlloc.hpp"
#include "ErrorHandler.hpp"

#ifdef ENABLE_MKL
#include "mkl.h"
#else
#include <cstdlib>
#endif

namespace mem {  // mem namespace

void*
malloc(const size_t size)
{
#ifdef ENABLE_MKL

    return MKL_malloc(size, VLX_ALIGN);

#else

    void* ptr = nullptr;

    int ierr = ::posix_memalign(&ptr, VLX_ALIGN, size);

    errors::assertMsgCritical(ierr == 0, "malloc: posix_memalign failed");

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

}  // namespace mem

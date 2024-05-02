//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef GpuSafeChecks_hpp
#define GpuSafeChecks_hpp

#include <iostream>

#define hipSafe(e)                                                                                                       \
    {                                                                                                                    \
        hipError_t err = (e);                                                                                            \
        if (err != hipSuccess)                                                                                           \
        {                                                                                                                \
            std::cerr << "HIP error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl;  \
            std::exit(EXIT_FAILURE);                                                                                     \
        }                                                                                                                \
    }

#define hipblasSafe(e)                                                                            \
    {                                                                                             \
        hipblasStatus_t err = (e);                                                                \
        if (err != HIPBLAS_STATUS_SUCCESS)                                                        \
        {                                                                                         \
            std::cerr << "hipBLAS error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                              \
        }                                                                                         \
    }

#define magmaSafe(e)                                                                                                   \
    {                                                                                                                  \
        magma_int_t err = (e);                                                                                         \
        if (err != MAGMA_SUCCESS) {                                                                                    \
            std::cerr << "MAGMA error in " << __FILE__ << ":" << __LINE__ << ": " << magma_strerror(err) << std::endl; \
            std::exit(EXIT_FAILURE);                                                                                   \
        }                                                                                                              \
    }

/*
#define hipsolverSafe(e)                                                                            \
    {                                                                                               \
        hipsolverStatus_t err = (e);                                                                \
        if (err != HIPSOLVER_STATUS_SUCCESS)                                                        \
        {                                                                                           \
            std::cerr << "hipSolver error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                                \
        }                                                                                           \
    }
*/

#endif /* GpuSafeChecks_hpp */

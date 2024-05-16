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

#define cudaSafe(e)                                                                                                       \
    {                                                                                                                    \
        cudaError_t err = (e);                                                                                            \
        if (err != cudaSuccess)                                                                                           \
        {                                                                                                                \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl;  \
            std::exit(EXIT_FAILURE);                                                                                     \
        }                                                                                                                \
    }

#define cublasSafe(e)                                                                            \
    {                                                                                             \
        cublasStatus_t err = (e);                                                                \
        if (err != CUBLAS_STATUS_SUCCESS)                                                        \
        {                                                                                         \
            std::cerr << "cudaBLAS error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                              \
        }                                                                                         \
    }

#define cusolverSafe(e)                                                                            \
    {                                                                                               \
        cusolverStatus_t err = (e);                                                                \
        if (err != CUSOLVER_STATUS_SUCCESS)                                                        \
        {                                                                                           \
            std::cerr << "cudaSolver error in " << __FILE__ << ":" << __LINE__ << ": " << std::endl; \
            std::exit(EXIT_FAILURE);                                                                \
        }                                                                                           \
    }

#endif /* GpuSafeChecks_hpp */

//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#ifndef GtoValuesGgaGPU_hpp
#define GtoValuesGgaGPU_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void gtoValuesGgaRecS(double*        gto_values,
                                 double*        gto_values_x,
                                 double*        gto_values_y,
                                 double*        gto_values_z,
                                 const uint32_t row_offset,
                                 const double*  gto_info,
                                 const double*  grid_x,
                                 const double*  grid_y,
                                 const double*  grid_z,
                                 const uint32_t grid_offset,
                                 const uint32_t nrows,
                                 const uint32_t npgtos,
                                 const uint32_t ncols);

__global__ void gtoValuesGgaRecP(double*        gto_values_p3,
                                 double*        gto_values_p3_x,
                                 double*        gto_values_p3_y,
                                 double*        gto_values_p3_z,
                                 const uint32_t row_offset,
                                 const double*  gto_info,
                                 const double*  grid_x,
                                 const double*  grid_y,
                                 const double*  grid_z,
                                 const uint32_t grid_offset,
                                 const uint32_t nrows,
                                 const uint32_t npgtos,
                                 const uint32_t ncols);

__global__ void gtoValuesGgaRecD(double*        gto_values_d5,
                                 double*        gto_values_d5_x,
                                 double*        gto_values_d5_y,
                                 double*        gto_values_d5_z,
                                 const uint32_t row_offset,
                                 const double   f2_3,
                                 const double*  gto_info,
                                 const double*  grid_x,
                                 const double*  grid_y,
                                 const double*  grid_z,
                                 const uint32_t grid_offset,
                                 const uint32_t nrows,
                                 const uint32_t npgtos,
                                 const uint32_t ncols);

}  // namespace gpu

#endif

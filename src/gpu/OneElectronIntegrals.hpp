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

#ifndef OneElectronIntegrals_hpp
#define OneElectronIntegrals_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void
computeOverlapAndKineticEnergySS(double*         mat_S,
                                 double*         mat_T,
                                 const double*   s_prim_info,
                                 const uint32_t  s_prim_count,
                                 const uint32_t* first_inds_local,
                                 const uint32_t* second_inds_local,
                                 const uint32_t  ss_prim_pair_count_local);

__global__ void
computeOverlapAndKineticEnergySP(double*         mat_S,
                                 double*         mat_T,
                                 const double*   s_prim_info,
                                 const uint32_t  s_prim_count,
                                 const double*   p_prim_info,
                                 const uint32_t  p_prim_count,
                                 const uint32_t* sp_first_inds_local,
                                 const uint32_t* sp_second_inds_local,
                                 const uint32_t  sp_prim_pair_count_local);

__global__ void
computeOverlapAndKineticEnergyPP(double*         mat_S,
                                 double*         mat_T,
                                 const double*   p_prim_info,
                                 const uint32_t  p_prim_count,
                                 const uint32_t* pp_first_inds_local,
                                 const uint32_t* pp_second_inds_local,
                                 const uint32_t  pp_prim_pair_count_local);

}  // namespace gpu

#endif

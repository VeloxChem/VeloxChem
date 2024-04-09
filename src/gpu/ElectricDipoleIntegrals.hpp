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

#ifndef ElectricDipole_hpp
#define ElectricDipole_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void
computeElectricDipoleSS(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const uint32_t* first_inds_local,
                        const uint32_t* second_inds_local,
                        const uint32_t  ss_prim_pair_count_local);

__global__ void
computeElectricDipoleSP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const uint32_t* sp_first_inds_local,
                        const uint32_t* sp_second_inds_local,
                        const uint32_t  sp_prim_pair_count_local);

__global__ void
computeElectricDipoleSD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   s_prim_info,
                        const uint32_t  s_prim_count,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* sd_first_inds_local,
                        const uint32_t* sd_second_inds_local,
                        const uint32_t  sd_prim_pair_count_local);

__global__ void
computeElectricDipolePP(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const uint32_t* pp_first_inds_local,
                        const uint32_t* pp_second_inds_local,
                        const uint32_t  pp_prim_pair_count_local);

__global__ void
computeElectricDipolePD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   p_prim_info,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* pd_first_inds_local,
                        const uint32_t* pd_second_inds_local,
                        const uint32_t  pd_prim_pair_count_local);

__global__ void
computeElectricDipoleDD(double*         mat_mu_X,
                        double*         mat_mu_Y,
                        double*         mat_mu_Z,
                        const double    origin_X,
                        const double    origin_Y,
                        const double    origin_Z,
                        const double*   d_prim_info,
                        const uint32_t  d_prim_count,
                        const uint32_t* dd_first_inds_local,
                        const uint32_t* dd_second_inds_local,
                        const uint32_t  dd_prim_pair_count_local);

}  // namespace gpu

#endif

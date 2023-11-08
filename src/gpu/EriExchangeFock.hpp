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

#ifndef EriExchangeFock_hpp
#define EriExchangeFock_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void cudaFockK(double*         mat_K,
                          const uint32_t* pari_inds_i,
                          const uint32_t* pari_inds_k,
                          const double*   s_prim_info,
                          const uint32_t* s_prim_aoinds,
                          const uint32_t  s_prim_count,
                          const double    max_D,
                          const double*   mat_D_full,
                          const double*   mat_Q_for_K,
                          const uint32_t* first_inds_for_K,
                          const uint32_t* second_inds_for_K,
                          const uint32_t  ss_prim_pair_count,
                          const uint32_t  naos);

}  // namespace gpu

#endif

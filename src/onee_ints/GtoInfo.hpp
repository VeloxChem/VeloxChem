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

#ifndef GtoInfo_hpp
#define GtoInfo_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"

namespace gtoinfo {

auto updatePrimitiveInfoForS(double* s_prim_info, int* s_prim_aoinds, const int s_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void;

auto updatePrimitiveInfoForP(double* p_prim_info, int* p_prim_aoinds, const int p_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void;

auto updatePrimitiveInfoForD(double* d_prim_info, int* d_prim_aoinds, const int d_prim_count, const std::vector<CGtoBlock>& gto_blocks) -> void;

auto updatePrimitiveInfoForF(double* d_prim_info, int* d_prim_aoinds, const int d_prim_count, const std::vector<CGtoBlock>& gto_blocks, const int ncgtos_d) -> void;

}  // namespace gtoinfo

#endif

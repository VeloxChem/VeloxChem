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

#ifndef BoysFuncTable_hpp
#define BoysFuncTable_hpp

#include <array>
#include <cstdint>
#include <vector>

namespace boysfunc {

auto getBoysFuncFactors() -> std::vector<double>;

auto getFullBoysFuncTable() -> std::vector<double>;

auto getBoysFuncTable(const int64_t N) -> std::array<std::array<double, 7>, 121>;

auto getBoysFunction(const double fa, const uint32_t N, const double* bf_table, const double* ft) -> std::vector<double>;

}  // namespace boysfunc

#endif /* BoysFuncTable_hpp */

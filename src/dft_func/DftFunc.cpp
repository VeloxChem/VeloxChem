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

#include "DftFunc.hpp"

#include <algorithm>
#include <ranges>

namespace gtoval {  // gtoval namespace

auto
distribute(CSubMatrix* matrix, const std::vector<double>& values, const size_t irow) -> void
{
    if (const size_t ncols = values.size(); ncols > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, ncols), [&](const auto i) {
            matrix->at({irow, i}) += values[i]; 
        });
    }
}

auto
distribute(CSubMatrix* matrix, const std::vector<double>& values, const double factor, const size_t irow) -> void
{
    if (const size_t ncols = values.size(); ncols > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, ncols), [&](const auto i) {
            matrix->at({irow, i}) += factor * values[i];
        });
    }
}

}  // namespace gtoval

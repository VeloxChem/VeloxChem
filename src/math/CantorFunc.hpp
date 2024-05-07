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

#ifndef CantorFunc_hpp
#define CantorFunc_hpp

#include <cmath>

#include "T2Pair.hpp"

namespace mathfunc {  // mathfunc namespace

/**
 Computes Cantor pairing index for two integer non-negative numbers.

 @return the Cantor index.
 */
inline auto
getCantorIndex(const T2Pair& pair) -> int64_t
{
    const auto x = pair.first;

    const auto y = pair.second;

    return (x + y + 1) * (x + y) / 2 + y;
}

/**
 Reduces Cantor pairing index to Cantor pair i.e. two integer non-negative numbers.

 @return the Cantor pair.
 */
inline auto
getCantorPair(const int64_t index) -> T2Pair
{
    const auto w = static_cast<int64_t>(std::floor(0.5 * (std::sqrt(8 * index + 1) - 1)));

    const auto y = index - (w * w + w) / 2;

    return {w - y, y};
}

}  // namespace mathfunc

#endif /* CantorFunc_hpp */

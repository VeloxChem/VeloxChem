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

#ifndef MatrixIndex_hpp
#define MatrixIndex_hpp

#include <cstdint>

namespace mathfunc {  // mathfunc namespace

/**
 Gets upper triangular matrix index (C++ indexing scheme).

 @param i the index of row in matrix.
 @param j the index of collumn in matrix..
 @return the upper triangular matrix index.
 */
inline auto
uplo_index(const int64_t i, const int64_t j) -> int64_t
{
    return i + j * (j + 1) / 2;
}

}  // namespace mathfunc

#endif /* MatrixIndex_hpp */
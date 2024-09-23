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

#ifndef DftFunc_hpp
#define DftFunc_hpp

#include <cstddef>
#include <vector>

#include "SubMatrix.hpp"

namespace gtoval {  // gtoval namespace

/// @brief Distributes basis function values into given submatrix.
/// @param matrix The submatrix to distribute basis function values.
/// @param values The vector of basis function values to distribute.
/// @param irow The index of row to distribute basis function values.
auto distribute(CSubMatrix* matrix, const std::vector<double>& values, const size_t irow) -> void;

// @brief Distributes basis function values into given submatrix.
/// @param matrix The submatrix to distribute basis function values.
/// @param values The vector of basis function values to distribute.
/// @param factor The scaling factor of values.
/// @param irow The index of row to distribute basis function values.
auto distribute(CSubMatrix* matrix, const std::vector<double>& values, const double factor, const size_t irow) -> void;

}  // namespace gtoval

#endif /* DftFunc_hpp */

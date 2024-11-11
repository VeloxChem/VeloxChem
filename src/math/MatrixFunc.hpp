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

#ifndef MatrixFunc_hpp
#define MatrixFunc_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {

/// @brief Creates a matrix.
/// @param basis The molecular basis to create matrix.
/// @param mtype  The matrix type of created matrix.
/// @return The matrix.
auto make_matrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix;

/// @brief Creates a matrix.
/// @param bra_basis The molecular basis to create matrix (rows data).
/// @param ket_basis The molecular basis to create matrix (columns data).
/// @return The matrix.
auto make_matrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix;

/// @brief Creates matrix.
/// @param label The exchange-correlation functional type label.
/// @param nrows The number of rows in submatrix.
/// @param ncols The number of columns in submatrix.
/// @return The matrix.
auto make_matrix(const std::string& label, const size_t nrows, const size_t ncols) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */

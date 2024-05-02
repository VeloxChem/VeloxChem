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

#include <cstdint>
#include <string>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {  // matfunc namespace

/**
 Creates matrix.

 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrix.
 */
auto makeMatrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix;

/**
 Creates matrix.

 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrix.
 */
auto makeMatrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix;

/**
 Creates matrix.

 @param label the exchange-correlation functional type label.
 @param nrows the number of rows in submatrix.
 @param ncols the number of columns in submatrix.
 @return the matrix.
 */
auto makeMatrix(const std::string& label, const int64_t nrows, const int64_t ncols) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */

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

#ifndef PrimitiveOverlapDD_YY_T
#define PrimitiveOverlapDD_YY_T

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <D_YY|1|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveOverlapDD_YY_T(TDoubleArray&       buffer_xx,
                                 TDoubleArray&       buffer_xy,
                                 TDoubleArray&       buffer_xz,
                                 TDoubleArray&       buffer_yy,
                                 TDoubleArray&       buffer_yz,
                                 TDoubleArray&       buffer_zz,
                                 const double        bra_exp,
                                 const double        bra_norm,
                                 const TPoint3D&     bra_coord,
                                 const TDoubleArray& ket_exps,
                                 const TDoubleArray& ket_norms,
                                 const TDoubleArray& ket_coords_x,
                                 const TDoubleArray& ket_coords_y,
                                 const TDoubleArray& ket_coords_z,
                                 const int64_t       ket_dim) -> void;

}  // namespace ovlrec

#endif /* PrimitiveOverlapDD_YY_T */

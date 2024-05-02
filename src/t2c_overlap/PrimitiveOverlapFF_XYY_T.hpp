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

#ifndef PrimitiveOverlapFF_XYY_T
#define PrimitiveOverlapFF_XYY_T

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <F_XYY|1|F>  integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
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
auto compPrimitiveOverlapFF_XYY_T(TDoubleArray&       buffer_xxx,
                                  TDoubleArray&       buffer_xxy,
                                  TDoubleArray&       buffer_xxz,
                                  TDoubleArray&       buffer_xyy,
                                  TDoubleArray&       buffer_xyz,
                                  TDoubleArray&       buffer_xzz,
                                  TDoubleArray&       buffer_yyy,
                                  TDoubleArray&       buffer_yyz,
                                  TDoubleArray&       buffer_yzz,
                                  TDoubleArray&       buffer_zzz,
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

#endif /* PrimitiveOverlapFF_XYY_T */

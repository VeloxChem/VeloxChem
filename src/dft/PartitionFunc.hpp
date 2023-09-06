//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#ifndef PartitionFunc_hpp
#define PartitionFunc_hpp

#include <cstdint>

#include "DenseMatrix.hpp"
#include "Point.hpp"

namespace partfunc {  // partfunc namespace

/**
 Applies SSF partitioning scheme to selecte grid weights.
 Reference: R. E. Stratmann, G. E. Scuseria, and M. J. Frisch, Chem. Phys.
 Lett., 213 (257), 1996.

 @param rawGridPoints the raw grid points.
 @param minDistance the distance to closest neighbouring atom.
 @param gridOffset the atom grid points offset in raw grid points.
 @param nGridPoints the number of grid points.
 @param atomCoordinates the Cartesian coordinates of atoms.
 @param nAtoms the number of atoms.
 @param idAtomic the index of atom.
 */
auto ssf(CDenseMatrix*   rawGridPoints,
         const double    minDistance,
         const int64_t   gridOffset,
         const int64_t   nGridPoints,
         const TPoint3D* atomCoordinates,
         const int64_t   nAtoms,
         const int64_t   idAtomic) -> void;

/**
 Computes polynomial weight function in eliptical coordinates. See Eq. 14
 in R. E. Stratmann, G. E. Scuseria, and M. J. Frisch, Chem. Phys. Lett.,
 213 (257), 1996.

 @param eRadius the eliptical coordinate.
 @return the polynomial weight.
 */
inline auto zeta(const double eRadius) -> double;

}  // namespace partfunc

#endif /* PartitionFunc_hpp */

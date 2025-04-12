//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

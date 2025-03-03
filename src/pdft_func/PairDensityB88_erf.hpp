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

#ifndef PairDensityB88_erf_hpp
#define PairDensityB88_erf_hpp

#include <cstdint>

namespace pdftb88_erf {  // pdftb88_erf namespace

/**
 Computes Exc and Vxc Fock for short-ranged version of pair-density Becke88 exchange.

 @param np the number of grid points.
 @param rho the density.
 @param exc the functional value.
 @param vrho the 1st-order functional derivative wrt density.
 @param mu the range-separation parameter.
 */
void compute_exc_vxc(const int np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma, const double mu);

}  // namespace pdftb88_erf

#endif /* PairDensityB88_erf_hpp */

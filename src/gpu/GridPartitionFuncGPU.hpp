//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#ifndef GridPartitionFuncGPU_hpp
#define GridPartitionFuncGPU_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "Point.hpp"

namespace gpu {

auto applyGridPartitionFunc(CDenseMatrix*                rawGridPoints,
                            const std::vector<uint32_t>& atomIdsOfGridPoints,
                            const std::vector<double>&   atomMinDistances,
                            const int64_t                nGridPoints,
                            const TPoint3D*              atomCoordinates,
                            const int64_t                nAtoms,
                            const int64_t                numGpusPerNode) -> void;

}  // namespace gpu

#endif

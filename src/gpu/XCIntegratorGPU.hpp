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

#ifndef XCIntegratorGPU_hpp
#define XCIntegratorGPU_hpp

#include "DenseMatrix.hpp"
#include "GtoBlock.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"

namespace gpu {

auto computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> CDenseMatrix;

auto getGtoValuesForLda(const CGtoBlock&            gto_block,
                        const std::vector<double>&  grid_coords_x,
                        const std::vector<double>&  grid_coords_y,
                        const std::vector<double>&  grid_coords_z,
                        const std::vector<int64_t>& gtos_mask) -> CMatrix;

auto getLdaValuesRecS(const CGtoBlock&            gto_block,
                      const std::vector<double>&  grid_coords_x,
                      const std::vector<double>&  grid_coords_y,
                      const std::vector<double>&  grid_coords_z,
                      const std::vector<int64_t>& gtos_mask) -> CMatrix;

}  // namespace gpu

#endif

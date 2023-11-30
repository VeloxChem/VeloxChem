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

#ifndef FockDriverGPU_hpp
#define FockDriverGPU_hpp

#include <vector>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ScreeningData.hpp"

namespace gpu {

auto computeOneElectronIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>;

auto computeFockOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix, CScreeningData& screening) -> CDenseMatrix;

auto computeMatrixMultiplication(double* C, const double* A, const double* B, const int64_t nrows_A, const int64_t ncols_A, const int64_t ncols_B) -> void;

auto diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void;

}  // namespace gpu

#endif

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

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "GtoBlock.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCFunctional.hpp"

namespace gpu {

auto computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) -> CDenseMatrix;

auto computeGtoValuesAndDerivativesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid)
    -> std::vector<CDenseMatrix>;

auto integrateVxcFock(const CMolecule&        molecule,
                      const CMolecularBasis&  basis,
                      const CAODensityMatrix& densityMatrix,
                      const CMolecularGrid&   molecularGrid,
                      const std::string&      xcFuncLabel,
                      const int64_t           numGpusPerNode) -> CAOKohnShamMatrix;

auto integrateFxcFock(CDenseMatrix&           aoFockMatrix,
                      const CMolecule&        molecule,
                      const CMolecularBasis&  basis,
                      const CAODensityMatrix& rwDensityMatrix,
                      const CAODensityMatrix& gsDensityMatrix,
                      const CMolecularGrid&   molecularGrid,
                      const std::string&      xcFuncLabel,
                      const int64_t           numGpusPerNode) -> void;

}  // namespace gpu

#endif

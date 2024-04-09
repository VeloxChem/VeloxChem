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

auto computeQMatrixOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> CDenseMatrix;

auto computeOneElectronIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>;

auto computeElectricDipoleIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const std::vector<double>& origin, const CScreeningData& screening) -> std::vector<CDenseMatrix>;

auto computeLinearMomentumIntegralsOnGPU(const CMolecule& molecule, const CMolecularBasis& basis, const CScreeningData& screening) -> std::vector<CDenseMatrix>;

auto transformDensity(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> CDenseMatrix;

auto computeFockOnGPU(const              CMolecule& molecule,
                      const              CMolecularBasis& basis,
                      const              CAODensityMatrix& densityMatrix,
                      const double       prefac_coulomb,
                      const double       frac_exact_exchange,
                      const double       omega,
                      const std::string& flag_K,
                      const double       eri_threshold,
                      const double       prelink_threshold,
                      CScreeningData&    screening) -> CDenseMatrix;

auto computeDotProduct(const double* A, const double* B, const int64_t size) -> double;

auto computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& pointers, const int64_t size) -> void;

auto computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                        const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void;

auto transformMatrix(double* transformed_F, const double* X, const double* F,
                     const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void;

auto computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                                 const int64_t m_inp, const int64_t k_inp, const int64_t n_inp) -> void;

auto diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void;

auto diagonalizeMatrixMultiGPU(double* A, double* W, const int64_t nrows_A, const int64_t num_gpus_per_node) -> void;

}  // namespace gpu

#endif

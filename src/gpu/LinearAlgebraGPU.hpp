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

#ifndef LinearAlgebraGPU_hpp
#define LinearAlgebraGPU_hpp

#include <string>
#include <vector>

namespace gpu {

auto computeDotProduct(const double* A, const double* B, const int64_t size) -> double;

auto computeWeightedSum(double* weighted_data, const std::vector<double>& weights, const std::vector<const double*>& pointers, const int64_t size) -> void;

auto computeErrorVector(double* errvec, const double* X, const double* F, const double* D, const double* S,
                        const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void;

auto transformMatrix(double* transformed_F, const double* X, const double* F,
                     const int64_t nmo_inp, const int64_t nao_inp, const std::string& trans_X) -> void;

auto computeMatrixMultiplication(double* C, const double* A, const double* B, const std::string& trans_A, const std::string& trans_B,
                                 const int64_t m_inp, const int64_t k_inp, const int64_t n_inp) -> void;

auto diagonalizeMatrix(double* A, double* D, const int64_t nrows_A) -> void;

}  // namespace gpu

#endif

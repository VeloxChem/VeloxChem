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

#ifndef ElectricFieldRec_hpp
#define ElectricFieldRec_hpp

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace onee {  // onee namespace

auto computeElectricFieldRecSS(const double F1_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m) -> double;

auto computeElectricFieldRecSP(const double F2_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    b0,
                               const double PB_0) -> double;

auto computeElectricFieldRecSD(const double F3_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    b0,
                               const int    b1,
                               const double PB_0,
                               const double PB_1) -> double;

auto computeElectricFieldRecSF(const double F4_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    b0,
                               const int    b1,
                               const int    b2,
                               const double PB_0,
                               const double PB_1,
                               const double PB_2) -> double;

auto computeElectricFieldRecPP(const double F3_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const double PA_0,
                               const int    b0,
                               const double PB_0) -> double;

auto computeElectricFieldRecPD(const double F4_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const double PA_0,
                               const int    b0,
                               const int    b1,
                               const double PB_0,
                               const double PB_1) -> double;

auto computeElectricFieldRecPF(const double F5_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const double PA_0,
                               const int    b0,
                               const int    b1,
                               const int    b2,
                               const double PB_0,
                               const double PB_1,
                               const double PB_2) -> double;

auto computeElectricFieldRecDD(const double F5_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const int    a1,
                               const double PA_0,
                               const double PA_1,
                               const int    b0,
                               const int    b1,
                               const double PB_0,
                               const double PB_1) -> double;

auto computeElectricFieldRecDF(const double F6_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const int    a1,
                               const double PA_0,
                               const double PA_1,
                               const int    b0,
                               const int    b1,
                               const int    b2,
                               const double PB_0,
                               const double PB_1,
                               const double PB_2) -> double;

auto computeElectricFieldRecFF(const double F7_t[],
                               const double a_i,
                               const double a_j,
                               const double PC[],
                               const int    m,
                               const double delta[][3],
                               const int    a0,
                               const int    a1,
                               const int    a2,
                               const double PA_0,
                               const double PA_1,
                               const double PA_2,
                               const int    b0,
                               const int    b1,
                               const int    b2,
                               const double PB_0,
                               const double PB_1,
                               const double PB_2) -> double;

}  // namespace onee

#endif /* ElectricFieldRec_hpp */

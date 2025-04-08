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

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

#include "ElectricFieldRec.hpp"

#include <omp.h>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"

#define PAD_SIZE 8

#define MATH_CONST_PI 3.14159265358979323846

#define MATH_CONST_TWO_OVER_SQRT_PI 1.12837916709551255856

namespace onee {  // onee namespace

// J. Chem. Phys. 84, 3963-3974 (1986)

auto
computeElectricFieldRecSS(const double F1_t[],
                          const double a_i,
                          const double a_j,
                          const double PC[],
                          const int    m) -> double
{
    return (

        F1_t[1] * (

            2.0 * (a_i + a_j) * (
                PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecSP(const double F2_t[],
                          const double a_i,
                          const double a_j,
                          const double PC[],
                          const int    m,
                          const double delta[][3],
                          const int    b0,
                          const double PB_0) -> double
{
    return (

        F2_t[1] * (

            2.0 * (a_i + a_j) * (
                PB_0 * PC[m]
            )

            + (
                delta[b0][m]
            )

        )

        + F2_t[2] * (

            (-2.0) * (a_i + a_j) * (
                PC[b0] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecSD(const double F3_t[],
                          const double a_i,
                          const double a_j,
                          const double PC[],
                          const int    m,
                          const double delta[][3],
                          const int    b0,
                          const int    b1,
                          const double PB_0,
                          const double PB_1) -> double
{
    return (

        F3_t[1] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PC[m])
            )

            + 2.0 * (a_i + a_j) * (
                PB_0 * PB_1 * PC[m]
            )

            + (
                delta[b1][m] * (PB_0)
                + delta[b0][m] * (PB_1)
            )

        )

        + F3_t[2] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PC[m])
            )

            + (-2.0) * (a_i + a_j) * (
                PB_0 * PC[b1] * PC[m]
                + PB_1 * PC[b0] * PC[m]
            )

            + (-1.0) * (
                delta[b1][m] * (PC[b0])
                + delta[b0][m] * (PC[b1])
            )

        )

        + F3_t[3] * (

            2.0 * (a_i + a_j) * (
                PC[b0] * PC[b1] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecSF(const double F4_t[],
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
                          const double PB_2) -> double
{
    return (

        F4_t[1] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PB_0 * PC[m])
                + delta[b0][b2] * (PB_1 * PC[m])
                + delta[b0][b1] * (PB_2 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2])
            )

            + 2.0 * (a_i + a_j) * (
                PB_0 * PB_1 * PB_2 * PC[m]
            )

            + (
                delta[b2][m] * (PB_0 * PB_1)
                + delta[b1][m] * (PB_0 * PB_2)
                + delta[b0][m] * (PB_1 * PB_2)
            )

        )

        + F4_t[2] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PB_0 * PC[m] + PC[b0] * PC[m])
                + delta[b0][b2] * (PB_1 * PC[m] + PC[b1] * PC[m])
                + delta[b0][b1] * (PB_2 * PC[m] + PC[b2] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2])
            )

            + (-2.0) * (a_i + a_j) * (
                PB_0 * PB_1 * PC[b2] * PC[m]
                + PB_0 * PB_2 * PC[b1] * PC[m]
                + PB_1 * PB_2 * PC[b0] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                + delta[b1][m] * (PB_0 * PC[b2] + PB_2 * PC[b0])
                + delta[b0][m] * (PB_1 * PC[b2] + PB_2 * PC[b1])
            )

        )

        + F4_t[3] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PC[b0] * PC[m])
                + delta[b0][b2] * (PC[b1] * PC[m])
                + delta[b0][b1] * (PC[b2] * PC[m])
            )

            + 2.0 * (a_i + a_j) * (
                PB_0 * PC[b1] * PC[b2] * PC[m]
                + PB_1 * PC[b0] * PC[b2] * PC[m]
                + PB_2 * PC[b0] * PC[b1] * PC[m]
            )

            + (
                delta[b2][m] * (PC[b0] * PC[b1])
                + delta[b1][m] * (PC[b0] * PC[b2])
                + delta[b0][m] * (PC[b1] * PC[b2])
            )

        )

        + F4_t[4] * (

            (-2.0) * (a_i + a_j) * (
                PC[b0] * PC[b1] * PC[b2] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecPP(const double F3_t[],
                          const double a_i,
                          const double a_j,
                          const double PC[],
                          const int    m,
                          const double delta[][3],
                          const int    a0,
                          const double PA_0,
                          const int    b0,
                          const double PB_0) -> double
{
    return (

        F3_t[1] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[a0][b0] * (PC[m])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PB_0 * PC[m]
            )

            + (
                delta[b0][m] * (PA_0)
                + delta[a0][m] * (PB_0)
            )

        )

        + F3_t[2] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[a0][b0] * (PC[m])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PC[b0] * PC[m]
                + PB_0 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b0][m] * (PC[a0])
                + delta[a0][m] * (PC[b0])
            )

        )

        + F3_t[3] * (

            2.0 * (a_i + a_j) * (
                PC[a0] * PC[b0] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecPD(const double F4_t[],
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
                          const double PB_1) -> double
{
    return (

        F4_t[1] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PA_0 * PC[m])
                + delta[a0][b1] * (PB_0 * PC[m])
                + delta[a0][b0] * (PB_1 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PB_0 * PB_1 * PC[m]
            )

            + (
                delta[b1][m] * (PA_0 * PB_0)
                + delta[b0][m] * (PA_0 * PB_1)
                + delta[a0][m] * (PB_0 * PB_1)
            )

        )

        + F4_t[2] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PA_0 * PC[m] + PC[a0] * PC[m])
                + delta[a0][b1] * (PB_0 * PC[m] + PC[b0] * PC[m])
                + delta[a0][b0] * (PB_1 * PC[m] + PC[b1] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PB_0 * PC[b1] * PC[m]
                + PA_0 * PB_1 * PC[b0] * PC[m]
                + PB_0 * PB_1 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b1][m] * (PA_0 * PC[b0] + PB_0 * PC[a0])
                + delta[b0][m] * (PA_0 * PC[b1] + PB_1 * PC[a0])
                + delta[a0][m] * (PB_0 * PC[b1] + PB_1 * PC[b0])
            )

        )

        + F4_t[3] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PC[a0] * PC[m])
                + delta[a0][b1] * (PC[b0] * PC[m])
                + delta[a0][b0] * (PC[b1] * PC[m])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PC[b0] * PC[b1] * PC[m]
                + PB_0 * PC[a0] * PC[b1] * PC[m]
                + PB_1 * PC[a0] * PC[b0] * PC[m]
            )

            + (
                delta[b1][m] * (PC[a0] * PC[b0])
                + delta[b0][m] * (PC[a0] * PC[b1])
                + delta[a0][m] * (PC[b0] * PC[b1])
            )

        )

        + F4_t[4] * (

            (-2.0) * (a_i + a_j) * (
                PC[a0] * PC[b0] * PC[b1] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecPF(const double F5_t[],
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
                          const double PB_2) -> double
{
    return (

        F5_t[1] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PC[m])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PB_0 * PC[m])
                + delta[b0][b2] * (PA_0 * PB_1 * PC[m])
                + delta[b0][b1] * (PA_0 * PB_2 * PC[m])
                + delta[a0][b2] * (PB_0 * PB_1 * PC[m])
                + delta[a0][b1] * (PB_0 * PB_2 * PC[m])
                + delta[a0][b0] * (PB_1 * PB_2 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0)
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PB_0)
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PB_1)
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_2)
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PB_0 * PB_1 * PB_2 * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PB_0 * PB_1)
                + delta[b1][m] * (PA_0 * PB_0 * PB_2)
                + delta[b0][m] * (PA_0 * PB_1 * PB_2)
                + delta[a0][m] * (PB_0 * PB_1 * PB_2)
            )

        )

        + F5_t[2] * (

            (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PC[m])
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                + delta[b0][b2] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                + delta[b0][b1] * (PA_0 * PB_2 * PC[m] + PA_0 * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[m])
                + delta[a0][b2] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                + delta[a0][b1] * (PB_0 * PB_2 * PC[m] + PB_0 * PC[b2] * PC[m] + PB_2 * PC[b0] * PC[m])
                + delta[a0][b0] * (PB_1 * PB_2 * PC[m] + PB_1 * PC[b2] * PC[m] + PB_2 * PC[b1] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 + PC[a0])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PB_0 + PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PB_1 + PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_2 + PC[b2])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PB_0 * PB_1 * PC[b2] * PC[m]
                + PA_0 * PB_0 * PB_2 * PC[b1] * PC[m]
                + PA_0 * PB_1 * PB_2 * PC[b0] * PC[m]
                + PB_0 * PB_1 * PB_2 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                + delta[b1][m] * (PA_0 * PB_0 * PC[b2] + PA_0 * PB_2 * PC[b0] + PB_0 * PB_2 * PC[a0])
                + delta[b0][m] * (PA_0 * PB_1 * PC[b2] + PA_0 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[a0])
                + delta[a0][m] * (PB_0 * PB_1 * PC[b2] + PB_0 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[b0])
            )

        )

        + F5_t[3] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PC[m])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                + delta[b0][b2] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                + delta[b0][b1] * (PA_0 * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[m] + PC[a0] * PC[b2] * PC[m])
                + delta[a0][b2] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PB_0 * PC[b2] * PC[m] + PB_2 * PC[b0] * PC[m] + PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PB_1 * PC[b2] * PC[m] + PB_2 * PC[b1] * PC[m] + PC[b1] * PC[b2] * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PC[a0])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[b2])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PB_0 * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PB_1 * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PB_2 * PC[b0] * PC[b1] * PC[m]
                + PB_0 * PB_1 * PC[a0] * PC[b2] * PC[m]
                + PB_0 * PB_2 * PC[a0] * PC[b1] * PC[m]
                + PB_1 * PB_2 * PC[a0] * PC[b0] * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                + delta[b1][m] * (PA_0 * PC[b0] * PC[b2] + PB_0 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b0])
                + delta[b0][m] * (PA_0 * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b1])
                + delta[a0][m] * (PB_0 * PC[b1] * PC[b2] + PB_1 * PC[b0] * PC[b2] + PB_2 * PC[b0] * PC[b1])
            )

        )

        + F5_t[4] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PC[a0] * PC[b0] * PC[m])
                + delta[b0][b2] * (PC[a0] * PC[b1] * PC[m])
                + delta[b0][b1] * (PC[a0] * PC[b2] * PC[m])
                + delta[a0][b2] * (PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PC[b1] * PC[b2] * PC[m])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PB_0 * PC[a0] * PC[b1] * PC[b2] * PC[m]
                + PB_1 * PC[a0] * PC[b0] * PC[b2] * PC[m]
                + PB_2 * PC[a0] * PC[b0] * PC[b1] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PC[a0] * PC[b0] * PC[b1])
                + delta[b1][m] * (PC[a0] * PC[b0] * PC[b2])
                + delta[b0][m] * (PC[a0] * PC[b1] * PC[b2])
                + delta[a0][m] * (PC[b0] * PC[b1] * PC[b2])
            )

        )

        + F5_t[5] * (

            2.0 * (a_i + a_j) * (
                PC[a0] * PC[b0] * PC[b1] * PC[b2] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecDD(const double F5_t[],
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
                          const double PB_1) -> double
{
    return (

        F5_t[1] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[m])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PA_0 * PA_1 * PC[m])
                + delta[a1][b1] * (PA_0 * PB_0 * PC[m])
                + delta[a1][b0] * (PA_0 * PB_1 * PC[m])
                + delta[a0][b1] * (PA_1 * PB_0 * PC[m])
                + delta[a0][b0] * (PA_1 * PB_1 * PC[m])
                + delta[a0][a1] * (PB_0 * PB_1 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0)
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1)
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0)
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1)
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PB_0 * PB_1 * PC[m]
            )

            + (
                delta[b1][m] * (PA_0 * PA_1 * PB_0)
                + delta[b0][m] * (PA_0 * PA_1 * PB_1)
                + delta[a1][m] * (PA_0 * PB_0 * PB_1)
                + delta[a0][m] * (PA_1 * PB_0 * PB_1)
            )

        )

        + F5_t[2] * (

            (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[m])
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                + delta[a1][b1] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                + delta[a1][b0] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                + delta[a0][b1] * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                + delta[a0][b0] * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                + delta[a0][a1] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 + PC[a0])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 + PC[a1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 + PC[b0])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 + PC[b1])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PA_1 * PB_0 * PC[b1] * PC[m]
                + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m]
                + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m]
                + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b1][m] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                + delta[b0][m] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                + delta[a1][m] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                + delta[a0][m] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
            )

        )

        + F5_t[3] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[m])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m] + PC[a0] * PC[a1] * PC[m])
                + delta[a1][b1] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                + delta[a1][b0] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m] + PC[a1] * PC[b0] * PC[m])
                + delta[a0][b0] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m] + PC[a1] * PC[b1] * PC[m])
                + delta[a0][a1] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m]
                + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m]
                + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m]
                + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m]
                + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m]
                + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m]
            )

            + (
                delta[b1][m] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
            )

        )

        + F5_t[4] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b0][b1] * (PC[a0] * PC[a1] * PC[m])
                + delta[a1][b1] * (PC[a0] * PC[b0] * PC[m])
                + delta[a1][b0] * (PC[a0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PC[a1] * PC[b0] * PC[m])
                + delta[a0][b0] * (PC[a1] * PC[b1] * PC[m])
                + delta[a0][a1] * (PC[b0] * PC[b1] * PC[m])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m]
                + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m]
                + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m]
                + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m]
            )

            + (-1.0) * (
                delta[b1][m] * (PC[a0] * PC[a1] * PC[b0])
                + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1])
                + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1])
                + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1])
            )

        )

        + F5_t[5] * (

            2.0 * (a_i + a_j) * (
                PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecDF(const double F6_t[],
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
                          const double PB_2) -> double
{
    return (

        F6_t[1] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PC[m])
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PC[m])
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PB_0 * PC[m])
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PB_1 * PC[m])
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_2 * PC[m])
            )

            + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PB_0 * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PB_1 * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PB_2 * PC[m])
                + delta[a1][b2] * (PA_0 * PB_0 * PB_1 * PC[m])
                + delta[a1][b1] * (PA_0 * PB_0 * PB_2 * PC[m])
                + delta[a1][b0] * (PA_0 * PB_1 * PB_2 * PC[m])
                + delta[a0][b2] * (PA_1 * PB_0 * PB_1 * PC[m])
                + delta[a0][b1] * (PA_1 * PB_0 * PB_2 * PC[m])
                + delta[a0][b0] * (PA_1 * PB_1 * PB_2 * PC[m])
                + delta[a0][a1] * (PB_0 * PB_1 * PB_2 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1)
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PB_0)
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PB_1)
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_2)
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PB_0)
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PB_1)
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_2)
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PB_0 * PB_1)
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_2)
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_2)
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PB_0 * PB_1 * PB_2 * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PA_1 * PB_0 * PB_1)
                + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PB_2)
                + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PB_2)
                + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PB_2)
                + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PB_2)
            )

        )

        + F6_t[2] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PC[m] * (-2.0) + PC[a0] * PC[m] * (-1.0))
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PC[m] * (-2.0) + PC[a1] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PB_0 * PC[m] * (-2.0) + PC[b0] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PB_1 * PC[m] * (-2.0) + PC[b1] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_2 * PC[m] * (-2.0) + PC[b2] * PC[m] * (-1.0))
            )

            + (-0.5) / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1])
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PB_0 * PC[m] + PA_0 * PA_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PB_1 * PC[m] + PA_0 * PA_1 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PB_2 * PC[m] + PA_0 * PA_1 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[m])
                + delta[a1][b2] * (PA_0 * PB_0 * PB_1 * PC[m] + PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                + delta[a1][b1] * (PA_0 * PB_0 * PB_2 * PC[m] + PA_0 * PB_0 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b0] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[m])
                + delta[a1][b0] * (PA_0 * PB_1 * PB_2 * PC[m] + PA_0 * PB_1 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[m])
                + delta[a0][b2] * (PA_1 * PB_0 * PB_1 * PC[m] + PA_1 * PB_0 * PC[b1] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[m])
                + delta[a0][b1] * (PA_1 * PB_0 * PB_2 * PC[m] + PA_1 * PB_0 * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b0] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[m])
                + delta[a0][b0] * (PA_1 * PB_1 * PB_2 * PC[m] + PA_1 * PB_1 * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[m])
                + delta[a0][a1] * (PB_0 * PB_1 * PB_2 * PC[m] + PB_0 * PB_1 * PC[b2] * PC[m] + PB_0 * PB_2 * PC[b1] * PC[m] + PB_1 * PB_2 * PC[b0] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1 + PA_0 * PC[a1] + PA_1 * PC[a0])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_2 + PA_0 * PC[b2] + PB_2 * PC[a0])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PB_0 + PA_1 * PC[b0] + PB_0 * PC[a1])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PB_1 + PA_1 * PC[b1] + PB_1 * PC[a1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_2 + PA_1 * PC[b2] + PB_2 * PC[a1])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_2 + PB_0 * PC[b2] + PB_2 * PC[b0])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_2 + PB_1 * PC[b2] + PB_2 * PC[b1])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PA_1 * PB_0 * PB_1 * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_0 * PB_2 * PC[b1] * PC[m]
                + PA_0 * PA_1 * PB_1 * PB_2 * PC[b0] * PC[m]
                + PA_0 * PB_0 * PB_1 * PB_2 * PC[a1] * PC[m]
                + PA_1 * PB_0 * PB_1 * PB_2 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PA_0 * PA_1 * PB_0 * PC[b1] + PA_0 * PA_1 * PB_1 * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] + PA_1 * PB_0 * PB_1 * PC[a0])
                + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PC[b2] + PA_0 * PA_1 * PB_2 * PC[b0] + PA_0 * PB_0 * PB_2 * PC[a1] + PA_1 * PB_0 * PB_2 * PC[a0])
                + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PC[b2] + PA_0 * PA_1 * PB_2 * PC[b1] + PA_0 * PB_1 * PB_2 * PC[a1] + PA_1 * PB_1 * PB_2 * PC[a0])
                + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PC[b2] + PA_0 * PB_0 * PB_2 * PC[b1] + PA_0 * PB_1 * PB_2 * PC[b0] + PB_0 * PB_1 * PB_2 * PC[a0])
                + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PC[b2] + PA_1 * PB_0 * PB_2 * PC[b1] + PA_1 * PB_1 * PB_2 * PC[b0] + PB_0 * PB_1 * PB_2 * PC[a1])
            )

        )

        + F6_t[3] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PC[m] + PC[a0] * PC[m] * 2.0)
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PC[m] + PC[a1] * PC[m] * 2.0)
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PB_0 * PC[m] + PC[b0] * PC[m] * 2.0)
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PB_1 * PC[m] + PC[b1] * PC[m] * 2.0)
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_2 * PC[m] + PC[b2] * PC[m] * 2.0)
            )

            + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[m] + PA_0 * PC[a1] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[m] + PA_1 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[m])
                + delta[a1][b2] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                + delta[a1][b1] * (PA_0 * PB_0 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b0] * PC[m])
                + delta[a1][b0] * (PA_0 * PB_1 * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b1] * PC[m] + PA_0 * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[m] + PB_1 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b1] * PC[m])
                + delta[a0][b2] * (PA_1 * PB_0 * PC[b1] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[m])
                + delta[a0][b1] * (PA_1 * PB_0 * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b0] * PC[m] + PA_1 * PC[b0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[m] + PB_0 * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b0] * PC[m])
                + delta[a0][b0] * (PA_1 * PB_1 * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b1] * PC[m] + PA_1 * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[m] + PB_1 * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b1] * PC[m])
                + delta[a0][a1] * (PB_0 * PB_1 * PC[b2] * PC[m] + PB_0 * PB_2 * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[b0] * PC[m] + PB_1 * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[b0] * PC[b1] * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PC[a1] + PA_1 * PC[a0] + PC[a0] * PC[a1])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[b2] + PB_2 * PC[a0] + PC[a0] * PC[b2])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PC[b0] + PB_0 * PC[a1] + PC[a1] * PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PC[b1] + PB_1 * PC[a1] + PC[a1] * PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[b2] + PB_2 * PC[a1] + PC[a1] * PC[b2])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PC[b2] + PB_2 * PC[b0] + PC[b0] * PC[b2])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PC[b2] + PB_2 * PC[b1] + PC[b1] * PC[b2])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PB_0 * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_1 * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_2 * PC[b0] * PC[b1] * PC[m]
                + PA_0 * PB_0 * PB_1 * PC[a1] * PC[b2] * PC[m]
                + PA_0 * PB_0 * PB_2 * PC[a1] * PC[b1] * PC[m]
                + PA_0 * PB_1 * PB_2 * PC[a1] * PC[b0] * PC[m]
                + PA_1 * PB_0 * PB_1 * PC[a0] * PC[b2] * PC[m]
                + PA_1 * PB_0 * PB_2 * PC[a0] * PC[b1] * PC[m]
                + PA_1 * PB_1 * PB_2 * PC[a0] * PC[b0] * PC[m]
                + PB_0 * PB_1 * PB_2 * PC[a0] * PC[a1] * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[b1] + PA_0 * PB_1 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] * PC[b1] + PA_1 * PB_1 * PC[a0] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1])
                + delta[b1][m] * (PA_0 * PA_1 * PC[b0] * PC[b2] + PA_0 * PB_0 * PC[a1] * PC[b2] + PA_0 * PB_2 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] * PC[b2] + PA_1 * PB_2 * PC[a0] * PC[b0] + PB_0 * PB_2 * PC[a0] * PC[a1])
                + delta[b0][m] * (PA_0 * PA_1 * PC[b1] * PC[b2] + PA_0 * PB_1 * PC[a1] * PC[b2] + PA_0 * PB_2 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] * PC[b2] + PA_1 * PB_2 * PC[a0] * PC[b1] + PB_1 * PB_2 * PC[a0] * PC[a1])
                + delta[a1][m] * (PA_0 * PB_0 * PC[b1] * PC[b2] + PA_0 * PB_1 * PC[b0] * PC[b2] + PA_0 * PB_2 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[b2] + PB_0 * PB_2 * PC[a0] * PC[b1] + PB_1 * PB_2 * PC[a0] * PC[b0])
                + delta[a0][m] * (PA_1 * PB_0 * PC[b1] * PC[b2] + PA_1 * PB_1 * PC[b0] * PC[b2] + PA_1 * PB_2 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] * PC[b2] + PB_0 * PB_2 * PC[a1] * PC[b1] + PB_1 * PB_2 * PC[a1] * PC[b0])
            )

        )

        + F6_t[4] * (

            (-0.5) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PC[a0] * PC[m])
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PC[a1] * PC[m])
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PC[b0] * PC[m])
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PC[b1] * PC[m])
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[b2] * PC[m])
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PC[a1] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[m])
                + delta[b0][b2] * (PA_0 * PC[a1] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[m])
                + delta[b0][b1] * (PA_0 * PC[a1] * PC[b2] * PC[m] + PA_1 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[m] + PC[a0] * PC[a1] * PC[b2] * PC[m])
                + delta[a1][b2] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[m])
                + delta[a1][b1] * (PA_0 * PC[b0] * PC[b2] * PC[m] + PB_0 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b2] * PC[m])
                + delta[a1][b0] * (PA_0 * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b1] * PC[m] + PC[a0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][b2] * (PA_1 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PA_1 * PC[b0] * PC[b2] * PC[m] + PB_0 * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b0] * PC[m] + PC[a1] * PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PA_1 * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b1] * PC[m] + PC[a1] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a1] * (PB_0 * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[b0] * PC[b1] * PC[m] + PC[b0] * PC[b1] * PC[b2] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PC[a0] * PC[a1])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PC[a0] * PC[b0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PC[a0] * PC[b1])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[b2])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PC[a1] * PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PC[a1] * PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[b2])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PC[b0] * PC[b1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0] * PC[b2])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1] * PC[b2])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PA_1 * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PB_2 * PC[a1] * PC[b0] * PC[b1] * PC[m]
                + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[b2] * PC[m]
                + PA_1 * PB_2 * PC[a0] * PC[b0] * PC[b1] * PC[m]
                + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[b2] * PC[m]
                + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[b1] * PC[m]
                + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] + PA_1 * PC[a0] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[b1] + PB_1 * PC[a0] * PC[a1] * PC[b0])
                + delta[b1][m] * (PA_0 * PC[a1] * PC[b0] * PC[b2] + PA_1 * PC[a0] * PC[b0] * PC[b2] + PB_0 * PC[a0] * PC[a1] * PC[b2] + PB_2 * PC[a0] * PC[a1] * PC[b0])
                + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] * PC[b2] + PA_1 * PC[a0] * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[a1] * PC[b2] + PB_2 * PC[a0] * PC[a1] * PC[b1])
                + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] * PC[b2] + PB_0 * PC[a0] * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[b0] * PC[b2] + PB_2 * PC[a0] * PC[b0] * PC[b1])
                + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] * PC[b2] + PB_0 * PC[a1] * PC[b1] * PC[b2] + PB_1 * PC[a1] * PC[b0] * PC[b2] + PB_2 * PC[a1] * PC[b0] * PC[b1])
            )

        )

        + F6_t[5] * (

            1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PC[a0] * PC[a1] * PC[b0] * PC[m])
                + delta[b0][b2] * (PC[a0] * PC[a1] * PC[b1] * PC[m])
                + delta[b0][b1] * (PC[a0] * PC[a1] * PC[b2] * PC[m])
                + delta[a1][b2] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                + delta[a1][b1] * (PC[a0] * PC[b0] * PC[b2] * PC[m])
                + delta[a1][b0] * (PC[a0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][b2] * (PC[a1] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PC[a1] * PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PC[a1] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a1] * (PC[b0] * PC[b1] * PC[b2] * PC[m])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[b2] * PC[m]
                + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[b2] * PC[m]
                + PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
            )

            + (
                delta[b2][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1])
                + delta[b1][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b2])
                + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1] * PC[b2])
                + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1] * PC[b2])
                + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1] * PC[b2])
            )

        )

        + F6_t[6] * (

            (-2.0) * (a_i + a_j) * (
                PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m]
            )

        )

    );
}

auto
computeElectricFieldRecFF(const double F7_t[],
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
                          const double PB_2) -> double
{
    return (

        F7_t[1] * (

            0.25 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0]) * (PC[m])
            )

            + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PA_0 * PA_1 * PC[m])
                + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PA_2 * PC[m])
                + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PA_0 * PB_0 * PC[m])
                + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PA_0 * PB_1 * PC[m])
                + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PA_0 * PB_2 * PC[m])
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PA_2 * PC[m])
                + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PA_1 * PB_0 * PC[m])
                + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PA_1 * PB_1 * PC[m])
                + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PA_1 * PB_2 * PC[m])
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PA_2 * PB_0 * PC[m])
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PA_2 * PB_1 * PC[m])
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_2 * PB_2 * PC[m])
                + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PB_0 * PB_1 * PC[m])
                + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PB_0 * PB_2 * PC[m])
                + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PB_1 * PB_2 * PC[m])
            )

            + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (PA_0)
                + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (PA_1)
                + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (PA_2)
                + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (PB_0)
                + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (PB_1)
                + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (PB_2)
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PA_2 * PB_0 * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PA_2 * PB_1 * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PA_2 * PB_2 * PC[m])
                + delta[a2][b2] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                + delta[a2][b1] * (PA_0 * PA_1 * PB_0 * PB_2 * PC[m])
                + delta[a2][b0] * (PA_0 * PA_1 * PB_1 * PB_2 * PC[m])
                + delta[a1][b2] * (PA_0 * PA_2 * PB_0 * PB_1 * PC[m])
                + delta[a1][b1] * (PA_0 * PA_2 * PB_0 * PB_2 * PC[m])
                + delta[a1][b0] * (PA_0 * PA_2 * PB_1 * PB_2 * PC[m])
                + delta[a1][a2] * (PA_0 * PB_0 * PB_1 * PB_2 * PC[m])
                + delta[a0][b2] * (PA_1 * PA_2 * PB_0 * PB_1 * PC[m])
                + delta[a0][b1] * (PA_1 * PA_2 * PB_0 * PB_2 * PC[m])
                + delta[a0][b0] * (PA_1 * PA_2 * PB_1 * PB_2 * PC[m])
                + delta[a0][a2] * (PA_1 * PB_0 * PB_1 * PB_2 * PC[m])
                + delta[a0][a1] * (PA_2 * PB_0 * PB_1 * PB_2 * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1 * PA_2)
                + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PA_0 * PA_1 * PB_0)
                + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PA_0 * PA_1 * PB_1)
                + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PA_0 * PA_1 * PB_2)
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PA_2 * PB_0)
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PA_2 * PB_1)
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_2 * PB_2)
                + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PA_0 * PB_0 * PB_1)
                + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PA_0 * PB_0 * PB_2)
                + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PA_0 * PB_1 * PB_2)
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PA_2 * PB_0)
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PA_2 * PB_1)
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_2 * PB_2)
                + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PA_1 * PB_0 * PB_1)
                + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PA_1 * PB_0 * PB_2)
                + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PA_1 * PB_1 * PB_2)
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PA_2 * PB_0 * PB_1)
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_2 * PB_0 * PB_2)
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_2 * PB_1 * PB_2)
                + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PB_0 * PB_1 * PB_2)
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PA_2 * PB_0 * PB_1 * PB_2 * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_1)
                + delta[b1][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_2)
                + delta[b0][m] * (PA_0 * PA_1 * PA_2 * PB_1 * PB_2)
                + delta[a2][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_2)
                + delta[a1][m] * (PA_0 * PA_2 * PB_0 * PB_1 * PB_2)
                + delta[a0][m] * (PA_1 * PA_2 * PB_0 * PB_1 * PB_2)
            )

        )

        + F7_t[2] * (

            (-0.75) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0]) * (PC[m])
            )

            + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PA_0 * PA_1 * PC[m] * (-2.0) + PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0))
                + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PA_2 * PC[m] * (-2.0) + PA_0 * PC[a2] * PC[m] * (-1.0) + PA_2 * PC[a0] * PC[m] * (-1.0))
                + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PA_0 * PB_0 * PC[m] * (-2.0) + PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0))
                + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PA_0 * PB_1 * PC[m] * (-2.0) + PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0))
                + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PA_0 * PB_2 * PC[m] * (-2.0) + PA_0 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a0] * PC[m] * (-1.0))
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PA_2 * PC[m] * (-2.0) + PA_1 * PC[a2] * PC[m] * (-1.0) + PA_2 * PC[a1] * PC[m] * (-1.0))
                + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PA_1 * PB_0 * PC[m] * (-2.0) + PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0))
                + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PA_1 * PB_1 * PC[m] * (-2.0) + PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0))
                + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PA_1 * PB_2 * PC[m] * (-2.0) + PA_1 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a1] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PA_2 * PB_0 * PC[m] * (-2.0) + PA_2 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a2] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PA_2 * PB_1 * PC[m] * (-2.0) + PA_2 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a2] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_2 * PB_2 * PC[m] * (-2.0) + PA_2 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a2] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PB_0 * PB_1 * PC[m] * (-2.0) + PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PB_0 * PB_2 * PC[m] * (-2.0) + PB_0 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[b0] * PC[m] * (-1.0))
                + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PB_1 * PB_2 * PC[m] * (-2.0) + PB_1 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[b1] * PC[m] * (-1.0))
            )

            + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (PA_2 * (-2.0) + PC[a2] * (-1.0))
                + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (PB_2 * (-2.0) + PC[b2] * (-1.0))
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PA_2 * PB_0 * PC[m] + PA_0 * PA_1 * PA_2 * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[a2] * PC[m] + PA_0 * PA_2 * PB_0 * PC[a1] * PC[m] + PA_1 * PA_2 * PB_0 * PC[a0] * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PA_2 * PB_1 * PC[m] + PA_0 * PA_1 * PA_2 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[a2] * PC[m] + PA_0 * PA_2 * PB_1 * PC[a1] * PC[m] + PA_1 * PA_2 * PB_1 * PC[a0] * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PA_2 * PB_2 * PC[m] + PA_0 * PA_1 * PA_2 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[a2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[a1] * PC[m] + PA_1 * PA_2 * PB_2 * PC[a0] * PC[m])
                + delta[a2][b2] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                + delta[a2][b1] * (PA_0 * PA_1 * PB_0 * PB_2 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_2 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_2 * PC[a0] * PC[m])
                + delta[a2][b0] * (PA_0 * PA_1 * PB_1 * PB_2 * PC[m] + PA_0 * PA_1 * PB_1 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[b1] * PC[m] + PA_0 * PB_1 * PB_2 * PC[a1] * PC[m] + PA_1 * PB_1 * PB_2 * PC[a0] * PC[m])
                + delta[a1][b2] * (PA_0 * PA_2 * PB_0 * PB_1 * PC[m] + PA_0 * PA_2 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_2 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a2] * PC[m] + PA_2 * PB_0 * PB_1 * PC[a0] * PC[m])
                + delta[a1][b1] * (PA_0 * PA_2 * PB_0 * PB_2 * PC[m] + PA_0 * PA_2 * PB_0 * PC[b2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_2 * PC[a2] * PC[m] + PA_2 * PB_0 * PB_2 * PC[a0] * PC[m])
                + delta[a1][b0] * (PA_0 * PA_2 * PB_1 * PB_2 * PC[m] + PA_0 * PA_2 * PB_1 * PC[b2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[b1] * PC[m] + PA_0 * PB_1 * PB_2 * PC[a2] * PC[m] + PA_2 * PB_1 * PB_2 * PC[a0] * PC[m])
                + delta[a1][a2] * (PA_0 * PB_0 * PB_1 * PB_2 * PC[m] + PA_0 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_0 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_0 * PB_1 * PB_2 * PC[b0] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a0] * PC[m])
                + delta[a0][b2] * (PA_1 * PA_2 * PB_0 * PB_1 * PC[m] + PA_1 * PA_2 * PB_0 * PC[b1] * PC[m] + PA_1 * PA_2 * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a2] * PC[m] + PA_2 * PB_0 * PB_1 * PC[a1] * PC[m])
                + delta[a0][b1] * (PA_1 * PA_2 * PB_0 * PB_2 * PC[m] + PA_1 * PA_2 * PB_0 * PC[b2] * PC[m] + PA_1 * PA_2 * PB_2 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_2 * PC[a2] * PC[m] + PA_2 * PB_0 * PB_2 * PC[a1] * PC[m])
                + delta[a0][b0] * (PA_1 * PA_2 * PB_1 * PB_2 * PC[m] + PA_1 * PA_2 * PB_1 * PC[b2] * PC[m] + PA_1 * PA_2 * PB_2 * PC[b1] * PC[m] + PA_1 * PB_1 * PB_2 * PC[a2] * PC[m] + PA_2 * PB_1 * PB_2 * PC[a1] * PC[m])
                + delta[a0][a2] * (PA_1 * PB_0 * PB_1 * PB_2 * PC[m] + PA_1 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_1 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_1 * PB_1 * PB_2 * PC[b0] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a1] * PC[m])
                + delta[a0][a1] * (PA_2 * PB_0 * PB_1 * PB_2 * PC[m] + PA_2 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_2 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_2 * PB_1 * PB_2 * PC[b0] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a2] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1 * PA_2 + PA_0 * PA_1 * PC[a2] + PA_0 * PA_2 * PC[a1] + PA_1 * PA_2 * PC[a0])
                + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PA_0 * PA_1 * PB_2 + PA_0 * PA_1 * PC[b2] + PA_0 * PB_2 * PC[a1] + PA_1 * PB_2 * PC[a0])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PA_2 * PB_0 + PA_0 * PA_2 * PC[b0] + PA_0 * PB_0 * PC[a2] + PA_2 * PB_0 * PC[a0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PA_2 * PB_1 + PA_0 * PA_2 * PC[b1] + PA_0 * PB_1 * PC[a2] + PA_2 * PB_1 * PC[a0])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_2 * PB_2 + PA_0 * PA_2 * PC[b2] + PA_0 * PB_2 * PC[a2] + PA_2 * PB_2 * PC[a0])
                + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PA_0 * PB_0 * PB_2 + PA_0 * PB_0 * PC[b2] + PA_0 * PB_2 * PC[b0] + PB_0 * PB_2 * PC[a0])
                + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PA_0 * PB_1 * PB_2 + PA_0 * PB_1 * PC[b2] + PA_0 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[a0])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PA_2 * PB_0 + PA_1 * PA_2 * PC[b0] + PA_1 * PB_0 * PC[a2] + PA_2 * PB_0 * PC[a1])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PA_2 * PB_1 + PA_1 * PA_2 * PC[b1] + PA_1 * PB_1 * PC[a2] + PA_2 * PB_1 * PC[a1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_2 * PB_2 + PA_1 * PA_2 * PC[b2] + PA_1 * PB_2 * PC[a2] + PA_2 * PB_2 * PC[a1])
                + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PA_1 * PB_0 * PB_2 + PA_1 * PB_0 * PC[b2] + PA_1 * PB_2 * PC[b0] + PB_0 * PB_2 * PC[a1])
                + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PA_1 * PB_1 * PB_2 + PA_1 * PB_1 * PC[b2] + PA_1 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[a1])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PA_2 * PB_0 * PB_1 + PA_2 * PB_0 * PC[b1] + PA_2 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a2])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_2 * PB_0 * PB_2 + PA_2 * PB_0 * PC[b2] + PA_2 * PB_2 * PC[b0] + PB_0 * PB_2 * PC[a2])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_2 * PB_1 * PB_2 + PA_2 * PB_1 * PC[b2] + PA_2 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[a2])
                + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PB_0 * PB_1 * PB_2 + PB_0 * PB_1 * PC[b2] + PB_0 * PB_2 * PC[b1] + PB_1 * PB_2 * PC[b0])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PA_1 * PA_2 * PB_0 * PB_1 * PC[b2] * PC[m]
                + PA_0 * PA_1 * PA_2 * PB_0 * PB_2 * PC[b1] * PC[m]
                + PA_0 * PA_1 * PA_2 * PB_1 * PB_2 * PC[b0] * PC[m]
                + PA_0 * PA_1 * PB_0 * PB_1 * PB_2 * PC[a2] * PC[m]
                + PA_0 * PA_2 * PB_0 * PB_1 * PB_2 * PC[a1] * PC[m]
                + PA_1 * PA_2 * PB_0 * PB_1 * PB_2 * PC[a0] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PC[b1] + PA_0 * PA_1 * PA_2 * PB_1 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_1 * PC[a2] + PA_0 * PA_2 * PB_0 * PB_1 * PC[a1] + PA_1 * PA_2 * PB_0 * PB_1 * PC[a0])
                + delta[b1][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PC[b2] + PA_0 * PA_1 * PA_2 * PB_2 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_2 * PC[a2] + PA_0 * PA_2 * PB_0 * PB_2 * PC[a1] + PA_1 * PA_2 * PB_0 * PB_2 * PC[a0])
                + delta[b0][m] * (PA_0 * PA_1 * PA_2 * PB_1 * PC[b2] + PA_0 * PA_1 * PA_2 * PB_2 * PC[b1] + PA_0 * PA_1 * PB_1 * PB_2 * PC[a2] + PA_0 * PA_2 * PB_1 * PB_2 * PC[a1] + PA_1 * PA_2 * PB_1 * PB_2 * PC[a0])
                + delta[a2][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[b2] + PA_0 * PA_1 * PB_0 * PB_2 * PC[b1] + PA_0 * PA_1 * PB_1 * PB_2 * PC[b0] + PA_0 * PB_0 * PB_1 * PB_2 * PC[a1] + PA_1 * PB_0 * PB_1 * PB_2 * PC[a0])
                + delta[a1][m] * (PA_0 * PA_2 * PB_0 * PB_1 * PC[b2] + PA_0 * PA_2 * PB_0 * PB_2 * PC[b1] + PA_0 * PA_2 * PB_1 * PB_2 * PC[b0] + PA_0 * PB_0 * PB_1 * PB_2 * PC[a2] + PA_2 * PB_0 * PB_1 * PB_2 * PC[a0])
                + delta[a0][m] * (PA_1 * PA_2 * PB_0 * PB_1 * PC[b2] + PA_1 * PA_2 * PB_0 * PB_2 * PC[b1] + PA_1 * PA_2 * PB_1 * PB_2 * PC[b0] + PA_1 * PB_0 * PB_1 * PB_2 * PC[a2] + PA_2 * PB_0 * PB_1 * PB_2 * PC[a1])
            )

        )

        + F7_t[3] * (

            0.75 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0]) * (PC[m])
            )

            + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a1] * PC[m])
                + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PA_2 * PC[m] + PA_0 * PC[a2] * PC[m] * 2.0 + PA_2 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a2] * PC[m])
                + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b0] * PC[m])
                + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b1] * PC[m])
                + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PA_0 * PB_2 * PC[m] + PA_0 * PC[b2] * PC[m] * 2.0 + PB_2 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b2] * PC[m])
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PA_2 * PC[m] + PA_1 * PC[a2] * PC[m] * 2.0 + PA_2 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[a2] * PC[m])
                + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b0] * PC[m])
                + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b1] * PC[m])
                + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PA_1 * PB_2 * PC[m] + PA_1 * PC[b2] * PC[m] * 2.0 + PB_2 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PA_2 * PB_0 * PC[m] + PA_2 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a2] * PC[m] * 2.0 + PC[a2] * PC[b0] * PC[m])
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PA_2 * PB_1 * PC[m] + PA_2 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a2] * PC[m] * 2.0 + PC[a2] * PC[b1] * PC[m])
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_2 * PB_2 * PC[m] + PA_2 * PC[b2] * PC[m] * 2.0 + PB_2 * PC[a2] * PC[m] * 2.0 + PC[a2] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b1] * PC[m])
                + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PB_0 * PB_2 * PC[m] + PB_0 * PC[b2] * PC[m] * 2.0 + PB_2 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PB_1 * PB_2 * PC[m] + PB_1 * PC[b2] * PC[m] * 2.0 + PB_2 * PC[b1] * PC[m] * 2.0 + PC[b1] * PC[b2] * PC[m])
            )

            + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (PA_0 + PC[a0] * 2.0)
                + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (PA_1 + PC[a1] * 2.0)
                + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (PA_2 + PC[a2] * 2.0)
                + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (PB_0 + PC[b0] * 2.0)
                + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (PB_1 + PC[b1] * 2.0)
                + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (PB_2 + PC[b2] * 2.0)
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PA_2 * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[a2] * PC[m] + PA_0 * PA_1 * PC[a2] * PC[b0] * PC[m] + PA_0 * PA_2 * PB_0 * PC[a1] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[a2] * PC[m] + PA_1 * PA_2 * PB_0 * PC[a0] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[a2] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[a1] * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PA_2 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[a2] * PC[m] + PA_0 * PA_1 * PC[a2] * PC[b1] * PC[m] + PA_0 * PA_2 * PB_1 * PC[a1] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[a2] * PC[m] + PA_1 * PA_2 * PB_1 * PC[a0] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[a2] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[a1] * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PA_2 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[a2] * PC[m] + PA_0 * PA_1 * PC[a2] * PC[b2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[a1] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[a2] * PC[m] + PA_1 * PA_2 * PB_2 * PC[a0] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[a2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[a1] * PC[m])
                + delta[a2][b2] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                + delta[a2][b1] * (PA_0 * PA_1 * PB_0 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_0 * PB_2 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_2 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[m])
                + delta[a2][b0] * (PA_0 * PA_1 * PB_1 * PC[b2] * PC[m] + PA_0 * PA_1 * PB_2 * PC[b1] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PB_2 * PC[a1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PB_2 * PC[a0] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[m])
                + delta[a1][b2] * (PA_0 * PA_2 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_2 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_2 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a2] * PC[m] + PA_0 * PB_0 * PC[a2] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a2] * PC[b0] * PC[m] + PA_2 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a2] * PC[m])
                + delta[a1][b1] * (PA_0 * PA_2 * PB_0 * PC[b2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[b0] * PC[m] + PA_0 * PA_2 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_0 * PB_2 * PC[a2] * PC[m] + PA_0 * PB_0 * PC[a2] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a2] * PC[b0] * PC[m] + PA_2 * PB_0 * PB_2 * PC[a0] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[a2] * PC[m])
                + delta[a1][b0] * (PA_0 * PA_2 * PB_1 * PC[b2] * PC[m] + PA_0 * PA_2 * PB_2 * PC[b1] * PC[m] + PA_0 * PA_2 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PB_2 * PC[a2] * PC[m] + PA_0 * PB_1 * PC[a2] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a2] * PC[b1] * PC[m] + PA_2 * PB_1 * PB_2 * PC[a0] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[a2] * PC[m])
                + delta[a1][a2] * (PA_0 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_0 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PB_2 * PC[b0] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[b0] * PC[m])
                + delta[a0][b2] * (PA_1 * PA_2 * PB_0 * PC[b1] * PC[m] + PA_1 * PA_2 * PB_1 * PC[b0] * PC[m] + PA_1 * PA_2 * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a2] * PC[m] + PA_1 * PB_0 * PC[a2] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a2] * PC[b0] * PC[m] + PA_2 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_2 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_2 * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[a2] * PC[m])
                + delta[a0][b1] * (PA_1 * PA_2 * PB_0 * PC[b2] * PC[m] + PA_1 * PA_2 * PB_2 * PC[b0] * PC[m] + PA_1 * PA_2 * PC[b0] * PC[b2] * PC[m] + PA_1 * PB_0 * PB_2 * PC[a2] * PC[m] + PA_1 * PB_0 * PC[a2] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a2] * PC[b0] * PC[m] + PA_2 * PB_0 * PB_2 * PC[a1] * PC[m] + PA_2 * PB_0 * PC[a1] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[a2] * PC[m])
                + delta[a0][b0] * (PA_1 * PA_2 * PB_1 * PC[b2] * PC[m] + PA_1 * PA_2 * PB_2 * PC[b1] * PC[m] + PA_1 * PA_2 * PC[b1] * PC[b2] * PC[m] + PA_1 * PB_1 * PB_2 * PC[a2] * PC[m] + PA_1 * PB_1 * PC[a2] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a2] * PC[b1] * PC[m] + PA_2 * PB_1 * PB_2 * PC[a1] * PC[m] + PA_2 * PB_1 * PC[a1] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[a2] * PC[m])
                + delta[a0][a2] * (PA_1 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_1 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_1 * PB_1 * PB_2 * PC[b0] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[b0] * PC[m])
                + delta[a0][a1] * (PA_2 * PB_0 * PB_1 * PC[b2] * PC[m] + PA_2 * PB_0 * PB_2 * PC[b1] * PC[m] + PA_2 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_2 * PB_1 * PB_2 * PC[b0] * PC[m] + PA_2 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_2 * PC[a2] * PC[m] + PB_0 * PB_1 * PC[a2] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a2] * PC[b1] * PC[m] + PB_1 * PB_2 * PC[a2] * PC[b0] * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1 * PC[a2] + PA_0 * PA_2 * PC[a1] + PA_0 * PC[a1] * PC[a2] + PA_1 * PA_2 * PC[a0] + PA_1 * PC[a0] * PC[a2] + PA_2 * PC[a0] * PC[a1])
                + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PA_0 * PA_1 * PC[b2] + PA_0 * PB_2 * PC[a1] + PA_0 * PC[a1] * PC[b2] + PA_1 * PB_2 * PC[a0] + PA_1 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[a1])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PA_2 * PC[b0] + PA_0 * PB_0 * PC[a2] + PA_0 * PC[a2] * PC[b0] + PA_2 * PB_0 * PC[a0] + PA_2 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a2])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PA_2 * PC[b1] + PA_0 * PB_1 * PC[a2] + PA_0 * PC[a2] * PC[b1] + PA_2 * PB_1 * PC[a0] + PA_2 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a2])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_2 * PC[b2] + PA_0 * PB_2 * PC[a2] + PA_0 * PC[a2] * PC[b2] + PA_2 * PB_2 * PC[a0] + PA_2 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[a2])
                + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PA_0 * PB_0 * PC[b2] + PA_0 * PB_2 * PC[b0] + PA_0 * PC[b0] * PC[b2] + PB_0 * PB_2 * PC[a0] + PB_0 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b0])
                + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PA_0 * PB_1 * PC[b2] + PA_0 * PB_2 * PC[b1] + PA_0 * PC[b1] * PC[b2] + PB_1 * PB_2 * PC[a0] + PB_1 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b1])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PA_2 * PC[b0] + PA_1 * PB_0 * PC[a2] + PA_1 * PC[a2] * PC[b0] + PA_2 * PB_0 * PC[a1] + PA_2 * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[a2])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PA_2 * PC[b1] + PA_1 * PB_1 * PC[a2] + PA_1 * PC[a2] * PC[b1] + PA_2 * PB_1 * PC[a1] + PA_2 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[a2])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_2 * PC[b2] + PA_1 * PB_2 * PC[a2] + PA_1 * PC[a2] * PC[b2] + PA_2 * PB_2 * PC[a1] + PA_2 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[a2])
                + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PA_1 * PB_0 * PC[b2] + PA_1 * PB_2 * PC[b0] + PA_1 * PC[b0] * PC[b2] + PB_0 * PB_2 * PC[a1] + PB_0 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[b0])
                + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PA_1 * PB_1 * PC[b2] + PA_1 * PB_2 * PC[b1] + PA_1 * PC[b1] * PC[b2] + PB_1 * PB_2 * PC[a1] + PB_1 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[b1])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PA_2 * PB_0 * PC[b1] + PA_2 * PB_1 * PC[b0] + PA_2 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a2] + PB_0 * PC[a2] * PC[b1] + PB_1 * PC[a2] * PC[b0])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_2 * PB_0 * PC[b2] + PA_2 * PB_2 * PC[b0] + PA_2 * PC[b0] * PC[b2] + PB_0 * PB_2 * PC[a2] + PB_0 * PC[a2] * PC[b2] + PB_2 * PC[a2] * PC[b0])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_2 * PB_1 * PC[b2] + PA_2 * PB_2 * PC[b1] + PA_2 * PC[b1] * PC[b2] + PB_1 * PB_2 * PC[a2] + PB_1 * PC[a2] * PC[b2] + PB_2 * PC[a2] * PC[b1])
                + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PB_0 * PB_1 * PC[b2] + PB_0 * PB_2 * PC[b1] + PB_0 * PC[b1] * PC[b2] + PB_1 * PB_2 * PC[b0] + PB_1 * PC[b0] * PC[b2] + PB_2 * PC[b0] * PC[b1])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PA_2 * PB_0 * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PA_2 * PB_1 * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PA_2 * PB_2 * PC[b0] * PC[b1] * PC[m]
                + PA_0 * PA_1 * PB_0 * PB_1 * PC[a2] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_0 * PB_2 * PC[a2] * PC[b1] * PC[m]
                + PA_0 * PA_1 * PB_1 * PB_2 * PC[a2] * PC[b0] * PC[m]
                + PA_0 * PA_2 * PB_0 * PB_1 * PC[a1] * PC[b2] * PC[m]
                + PA_0 * PA_2 * PB_0 * PB_2 * PC[a1] * PC[b1] * PC[m]
                + PA_0 * PA_2 * PB_1 * PB_2 * PC[a1] * PC[b0] * PC[m]
                + PA_0 * PB_0 * PB_1 * PB_2 * PC[a1] * PC[a2] * PC[m]
                + PA_1 * PA_2 * PB_0 * PB_1 * PC[a0] * PC[b2] * PC[m]
                + PA_1 * PA_2 * PB_0 * PB_2 * PC[a0] * PC[b1] * PC[m]
                + PA_1 * PA_2 * PB_1 * PB_2 * PC[a0] * PC[b0] * PC[m]
                + PA_1 * PB_0 * PB_1 * PB_2 * PC[a0] * PC[a2] * PC[m]
                + PA_2 * PB_0 * PB_1 * PB_2 * PC[a0] * PC[a1] * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PA_1 * PA_2 * PC[b0] * PC[b1] + PA_0 * PA_1 * PB_0 * PC[a2] * PC[b1] + PA_0 * PA_1 * PB_1 * PC[a2] * PC[b0] + PA_0 * PA_2 * PB_0 * PC[a1] * PC[b1] + PA_0 * PA_2 * PB_1 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[a2] + PA_1 * PA_2 * PB_0 * PC[a0] * PC[b1] + PA_1 * PA_2 * PB_1 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[a2] + PA_2 * PB_0 * PB_1 * PC[a0] * PC[a1])
                + delta[b1][m] * (PA_0 * PA_1 * PA_2 * PC[b0] * PC[b2] + PA_0 * PA_1 * PB_0 * PC[a2] * PC[b2] + PA_0 * PA_1 * PB_2 * PC[a2] * PC[b0] + PA_0 * PA_2 * PB_0 * PC[a1] * PC[b2] + PA_0 * PA_2 * PB_2 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_2 * PC[a1] * PC[a2] + PA_1 * PA_2 * PB_0 * PC[a0] * PC[b2] + PA_1 * PA_2 * PB_2 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_2 * PC[a0] * PC[a2] + PA_2 * PB_0 * PB_2 * PC[a0] * PC[a1])
                + delta[b0][m] * (PA_0 * PA_1 * PA_2 * PC[b1] * PC[b2] + PA_0 * PA_1 * PB_1 * PC[a2] * PC[b2] + PA_0 * PA_1 * PB_2 * PC[a2] * PC[b1] + PA_0 * PA_2 * PB_1 * PC[a1] * PC[b2] + PA_0 * PA_2 * PB_2 * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_2 * PC[a1] * PC[a2] + PA_1 * PA_2 * PB_1 * PC[a0] * PC[b2] + PA_1 * PA_2 * PB_2 * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_2 * PC[a0] * PC[a2] + PA_2 * PB_1 * PB_2 * PC[a0] * PC[a1])
                + delta[a2][m] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[b2] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[b2] + PA_0 * PA_1 * PB_2 * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[b2] + PA_0 * PB_0 * PB_2 * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_2 * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[b2] + PA_1 * PB_0 * PB_2 * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_2 * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_2 * PC[a0] * PC[a1])
                + delta[a1][m] * (PA_0 * PA_2 * PB_0 * PC[b1] * PC[b2] + PA_0 * PA_2 * PB_1 * PC[b0] * PC[b2] + PA_0 * PA_2 * PB_2 * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a2] * PC[b2] + PA_0 * PB_0 * PB_2 * PC[a2] * PC[b1] + PA_0 * PB_1 * PB_2 * PC[a2] * PC[b0] + PA_2 * PB_0 * PB_1 * PC[a0] * PC[b2] + PA_2 * PB_0 * PB_2 * PC[a0] * PC[b1] + PA_2 * PB_1 * PB_2 * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_2 * PC[a0] * PC[a2])
                + delta[a0][m] * (PA_1 * PA_2 * PB_0 * PC[b1] * PC[b2] + PA_1 * PA_2 * PB_1 * PC[b0] * PC[b2] + PA_1 * PA_2 * PB_2 * PC[b0] * PC[b1] + PA_1 * PB_0 * PB_1 * PC[a2] * PC[b2] + PA_1 * PB_0 * PB_2 * PC[a2] * PC[b1] + PA_1 * PB_1 * PB_2 * PC[a2] * PC[b0] + PA_2 * PB_0 * PB_1 * PC[a1] * PC[b2] + PA_2 * PB_0 * PB_2 * PC[a1] * PC[b1] + PA_2 * PB_1 * PB_2 * PC[a1] * PC[b0] + PB_0 * PB_1 * PB_2 * PC[a1] * PC[a2])
            )

        )

        + F7_t[4] * (

            (-0.25) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0]) * (PC[m])
            )

            + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PC[a2] * PC[m] * (-1.0) + PA_2 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a2] * PC[m] * (-2.0))
                + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-2.0))
                + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-2.0))
                + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PA_0 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b2] * PC[m] * (-2.0))
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PC[a2] * PC[m] * (-1.0) + PA_2 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[a2] * PC[m] * (-2.0))
                + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-2.0))
                + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-2.0))
                + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PA_1 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b2] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PA_2 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a2] * PC[m] * (-1.0) + PC[a2] * PC[b0] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PA_2 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a2] * PC[m] * (-1.0) + PC[a2] * PC[b1] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_2 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[a2] * PC[m] * (-1.0) + PC[a2] * PC[b2] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PB_0 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b2] * PC[m] * (-2.0))
                + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PB_1 * PC[b2] * PC[m] * (-1.0) + PB_2 * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[b2] * PC[m] * (-2.0))
            )

            + (-0.25) / ( (a_i + a_j) * (a_i + a_j) ) * (
                (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (PC[a0])
                + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (PC[a1])
                + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (PC[a2])
                + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (PC[b0])
                + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (PC[b1])
                + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (PC[b2])
            )

            + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PA_1 * PC[a2] * PC[b0] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[a2] * PC[m] + PA_0 * PC[a1] * PC[a2] * PC[b0] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[a2] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b0] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[a1] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[a2] * PC[m])
                + delta[b0][b2] * (PA_0 * PA_1 * PC[a2] * PC[b1] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[a2] * PC[m] + PA_0 * PC[a1] * PC[a2] * PC[b1] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[a2] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b1] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[a1] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[m])
                + delta[b0][b1] * (PA_0 * PA_1 * PC[a2] * PC[b2] * PC[m] + PA_0 * PA_2 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[a2] * PC[m] + PA_0 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PA_1 * PA_2 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[a2] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[a1] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[m])
                + delta[a2][b2] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                + delta[a2][b1] * (PA_0 * PA_1 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                + delta[a2][b0] * (PA_0 * PA_1 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a1] * PC[b1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[b1] * PC[m])
                + delta[a1][b2] * (PA_0 * PA_2 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a2] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a2] * PC[b0] * PC[m] + PA_0 * PC[a2] * PC[b0] * PC[b1] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_2 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a2] * PC[m] + PB_0 * PC[a0] * PC[a2] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a2] * PC[b0] * PC[m])
                + delta[a1][b1] * (PA_0 * PA_2 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_0 * PC[a2] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a2] * PC[b0] * PC[m] + PA_0 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PA_2 * PB_0 * PC[a0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[b0] * PC[m] + PA_2 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[a2] * PC[m] + PB_0 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a2] * PC[b0] * PC[m])
                + delta[a1][b0] * (PA_0 * PA_2 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PC[a2] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[a2] * PC[b1] * PC[m] + PA_0 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PA_2 * PB_1 * PC[a0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a0] * PC[b1] * PC[m] + PA_2 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[a2] * PC[m] + PB_1 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a2] * PC[b1] * PC[m])
                + delta[a1][a2] * (PA_0 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_0 * PB_2 * PC[b0] * PC[b1] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a0] * PC[b0] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][b2] * (PA_1 * PA_2 * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a2] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a2] * PC[b0] * PC[m] + PA_1 * PC[a2] * PC[b0] * PC[b1] * PC[m] + PA_2 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_2 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_2 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[a2] * PC[m] + PB_0 * PC[a1] * PC[a2] * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[a2] * PC[b0] * PC[m])
                + delta[a0][b1] * (PA_1 * PA_2 * PC[b0] * PC[b2] * PC[m] + PA_1 * PB_0 * PC[a2] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a2] * PC[b0] * PC[m] + PA_1 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PA_2 * PB_0 * PC[a1] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a1] * PC[b0] * PC[m] + PA_2 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[a2] * PC[m] + PB_0 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[a2] * PC[b0] * PC[m])
                + delta[a0][b0] * (PA_1 * PA_2 * PC[b1] * PC[b2] * PC[m] + PA_1 * PB_1 * PC[a2] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[a2] * PC[b1] * PC[m] + PA_1 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PA_2 * PB_1 * PC[a1] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[a1] * PC[b1] * PC[m] + PA_2 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[a2] * PC[m] + PB_1 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[a2] * PC[b1] * PC[m])
                + delta[a0][a2] * (PA_1 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_1 * PB_2 * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a1] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a1] * PC[b0] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][a1] * (PA_2 * PB_0 * PC[b1] * PC[b2] * PC[m] + PA_2 * PB_1 * PC[b0] * PC[b2] * PC[m] + PA_2 * PB_2 * PC[b0] * PC[b1] * PC[m] + PA_2 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PB_1 * PC[a2] * PC[b2] * PC[m] + PB_0 * PB_2 * PC[a2] * PC[b1] * PC[m] + PB_0 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PB_1 * PB_2 * PC[a2] * PC[b0] * PC[m] + PB_1 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a2] * PC[b0] * PC[b1] * PC[m])
            )

            + (-0.5) / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PC[a1] * PC[a2] + PA_1 * PC[a0] * PC[a2] + PA_2 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[a2])
                + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PA_0 * PC[a1] * PC[b2] + PA_1 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b2])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PC[a2] * PC[b0] + PA_2 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a2] + PC[a0] * PC[a2] * PC[b0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PC[a2] * PC[b1] + PA_2 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a2] + PC[a0] * PC[a2] * PC[b1])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[a2] * PC[b2] + PA_2 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[a2] + PC[a0] * PC[a2] * PC[b2])
                + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PA_0 * PC[b0] * PC[b2] + PB_0 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b2])
                + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PA_0 * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[b2] + PB_2 * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[b2])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PC[a2] * PC[b0] + PA_2 * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[a2] + PC[a1] * PC[a2] * PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PC[a2] * PC[b1] + PA_2 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[a2] + PC[a1] * PC[a2] * PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[a2] * PC[b2] + PA_2 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[a2] + PC[a1] * PC[a2] * PC[b2])
                + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PA_1 * PC[b0] * PC[b2] + PB_0 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b2])
                + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PA_1 * PC[b1] * PC[b2] + PB_1 * PC[a1] * PC[b2] + PB_2 * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[b2])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PA_2 * PC[b0] * PC[b1] + PB_0 * PC[a2] * PC[b1] + PB_1 * PC[a2] * PC[b0] + PC[a2] * PC[b0] * PC[b1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_2 * PC[b0] * PC[b2] + PB_0 * PC[a2] * PC[b2] + PB_2 * PC[a2] * PC[b0] + PC[a2] * PC[b0] * PC[b2])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_2 * PC[b1] * PC[b2] + PB_1 * PC[a2] * PC[b2] + PB_2 * PC[a2] * PC[b1] + PC[a2] * PC[b1] * PC[b2])
                + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PB_0 * PC[b1] * PC[b2] + PB_1 * PC[b0] * PC[b2] + PB_2 * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[b2])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PA_1 * PA_2 * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_0 * PC[a2] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_1 * PC[a2] * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PA_1 * PB_2 * PC[a2] * PC[b0] * PC[b1] * PC[m]
                + PA_0 * PA_2 * PB_0 * PC[a1] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_2 * PB_1 * PC[a1] * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PA_2 * PB_2 * PC[a1] * PC[b0] * PC[b1] * PC[m]
                + PA_0 * PB_0 * PB_1 * PC[a1] * PC[a2] * PC[b2] * PC[m]
                + PA_0 * PB_0 * PB_2 * PC[a1] * PC[a2] * PC[b1] * PC[m]
                + PA_0 * PB_1 * PB_2 * PC[a1] * PC[a2] * PC[b0] * PC[m]
                + PA_1 * PA_2 * PB_0 * PC[a0] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PA_2 * PB_1 * PC[a0] * PC[b0] * PC[b2] * PC[m]
                + PA_1 * PA_2 * PB_2 * PC[a0] * PC[b0] * PC[b1] * PC[m]
                + PA_1 * PB_0 * PB_1 * PC[a0] * PC[a2] * PC[b2] * PC[m]
                + PA_1 * PB_0 * PB_2 * PC[a0] * PC[a2] * PC[b1] * PC[m]
                + PA_1 * PB_1 * PB_2 * PC[a0] * PC[a2] * PC[b0] * PC[m]
                + PA_2 * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[b2] * PC[m]
                + PA_2 * PB_0 * PB_2 * PC[a0] * PC[a1] * PC[b1] * PC[m]
                + PA_2 * PB_1 * PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[m]
                + PB_0 * PB_1 * PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PA_0 * PA_1 * PC[a2] * PC[b0] * PC[b1] + PA_0 * PA_2 * PC[a1] * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[a2] * PC[b1] + PA_0 * PB_1 * PC[a1] * PC[a2] * PC[b0] + PA_1 * PA_2 * PC[a0] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[a2] * PC[b1] + PA_1 * PB_1 * PC[a0] * PC[a2] * PC[b0] + PA_2 * PB_0 * PC[a0] * PC[a1] * PC[b1] + PA_2 * PB_1 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[a2])
                + delta[b1][m] * (PA_0 * PA_1 * PC[a2] * PC[b0] * PC[b2] + PA_0 * PA_2 * PC[a1] * PC[b0] * PC[b2] + PA_0 * PB_0 * PC[a1] * PC[a2] * PC[b2] + PA_0 * PB_2 * PC[a1] * PC[a2] * PC[b0] + PA_1 * PA_2 * PC[a0] * PC[b0] * PC[b2] + PA_1 * PB_0 * PC[a0] * PC[a2] * PC[b2] + PA_1 * PB_2 * PC[a0] * PC[a2] * PC[b0] + PA_2 * PB_0 * PC[a0] * PC[a1] * PC[b2] + PA_2 * PB_2 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[a2])
                + delta[b0][m] * (PA_0 * PA_1 * PC[a2] * PC[b1] * PC[b2] + PA_0 * PA_2 * PC[a1] * PC[b1] * PC[b2] + PA_0 * PB_1 * PC[a1] * PC[a2] * PC[b2] + PA_0 * PB_2 * PC[a1] * PC[a2] * PC[b1] + PA_1 * PA_2 * PC[a0] * PC[b1] * PC[b2] + PA_1 * PB_1 * PC[a0] * PC[a2] * PC[b2] + PA_1 * PB_2 * PC[a0] * PC[a2] * PC[b1] + PA_2 * PB_1 * PC[a0] * PC[a1] * PC[b2] + PA_2 * PB_2 * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[a2])
                + delta[a2][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[b2] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[b2] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[b2] + PA_0 * PB_2 * PC[a1] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[b2] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[b2] + PA_1 * PB_2 * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[b2] + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[b0])
                + delta[a1][m] * (PA_0 * PA_2 * PC[b0] * PC[b1] * PC[b2] + PA_0 * PB_0 * PC[a2] * PC[b1] * PC[b2] + PA_0 * PB_1 * PC[a2] * PC[b0] * PC[b2] + PA_0 * PB_2 * PC[a2] * PC[b0] * PC[b1] + PA_2 * PB_0 * PC[a0] * PC[b1] * PC[b2] + PA_2 * PB_1 * PC[a0] * PC[b0] * PC[b2] + PA_2 * PB_2 * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a2] * PC[b2] + PB_0 * PB_2 * PC[a0] * PC[a2] * PC[b1] + PB_1 * PB_2 * PC[a0] * PC[a2] * PC[b0])
                + delta[a0][m] * (PA_1 * PA_2 * PC[b0] * PC[b1] * PC[b2] + PA_1 * PB_0 * PC[a2] * PC[b1] * PC[b2] + PA_1 * PB_1 * PC[a2] * PC[b0] * PC[b2] + PA_1 * PB_2 * PC[a2] * PC[b0] * PC[b1] + PA_2 * PB_0 * PC[a1] * PC[b1] * PC[b2] + PA_2 * PB_1 * PC[a1] * PC[b0] * PC[b2] + PA_2 * PB_2 * PC[a1] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] * PC[a2] * PC[b2] + PB_0 * PB_2 * PC[a1] * PC[a2] * PC[b1] + PB_1 * PB_2 * PC[a1] * PC[a2] * PC[b0])
            )

        )

        + F7_t[5] * (

            0.5 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * (
                (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PC[a0] * PC[a1] * PC[m])
                + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PC[a0] * PC[a2] * PC[m])
                + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PC[a0] * PC[b0] * PC[m])
                + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PC[a0] * PC[b1] * PC[m])
                + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PC[a0] * PC[b2] * PC[m])
                + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PC[a1] * PC[a2] * PC[m])
                + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PC[a1] * PC[b0] * PC[m])
                + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PC[a1] * PC[b1] * PC[m])
                + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PC[a1] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PC[a2] * PC[b0] * PC[m])
                + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PC[a2] * PC[b1] * PC[m])
                + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[a2] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PC[b0] * PC[b1] * PC[m])
                + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PC[b0] * PC[b2] * PC[m])
                + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PC[b1] * PC[b2] * PC[m])
            )

            + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PA_0 * PC[a1] * PC[a2] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b0] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[a2] * PC[m] + PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[m])
                + delta[b0][b2] * (PA_0 * PC[a1] * PC[a2] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b1] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[m] + PC[a0] * PC[a1] * PC[a2] * PC[b1] * PC[m])
                + delta[b0][b1] * (PA_0 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PA_1 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PA_2 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[m] + PC[a0] * PC[a1] * PC[a2] * PC[b2] * PC[m])
                + delta[a2][b2] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                + delta[a2][b1] * (PA_0 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b2] * PC[m])
                + delta[a2][b0] * (PA_0 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[b2] * PC[m])
                + delta[a1][b2] * (PA_0 * PC[a2] * PC[b0] * PC[b1] * PC[m] + PA_2 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a2] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a2] * PC[b0] * PC[m] + PC[a0] * PC[a2] * PC[b0] * PC[b1] * PC[m])
                + delta[a1][b1] * (PA_0 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PA_2 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_0 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a2] * PC[b0] * PC[m] + PC[a0] * PC[a2] * PC[b0] * PC[b2] * PC[m])
                + delta[a1][b0] * (PA_0 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PA_2 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a0] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[a2] * PC[b1] * PC[m] + PC[a0] * PC[a2] * PC[b1] * PC[b2] * PC[m])
                + delta[a1][a2] * (PA_0 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][b2] * (PA_1 * PC[a2] * PC[b0] * PC[b1] * PC[m] + PA_2 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[a2] * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[a2] * PC[b0] * PC[m] + PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PA_1 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PA_2 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PB_0 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[a2] * PC[b0] * PC[m] + PC[a1] * PC[a2] * PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PA_1 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PA_2 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a1] * PC[a2] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[a2] * PC[b1] * PC[m] + PC[a1] * PC[a2] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a2] * (PA_1 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a1] * (PA_2 * PC[b0] * PC[b1] * PC[b2] * PC[m] + PB_0 * PC[a2] * PC[b1] * PC[b2] * PC[m] + PB_1 * PC[a2] * PC[b0] * PC[b2] * PC[m] + PB_2 * PC[a2] * PC[b0] * PC[b1] * PC[m] + PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m])
            )

            + 0.5 / (a_i + a_j) * (
                (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PC[a0] * PC[a1] * PC[a2])
                + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PC[a0] * PC[a1] * PC[b0])
                + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PC[a0] * PC[a1] * PC[b1])
                + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PC[a0] * PC[a1] * PC[b2])
                + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PC[a0] * PC[a2] * PC[b0])
                + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PC[a0] * PC[a2] * PC[b1])
                + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[a2] * PC[b2])
                + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PC[a0] * PC[b0] * PC[b1])
                + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PC[a0] * PC[b0] * PC[b2])
                + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PC[a0] * PC[b1] * PC[b2])
                + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PC[a1] * PC[a2] * PC[b0])
                + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PC[a1] * PC[a2] * PC[b1])
                + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[a2] * PC[b2])
                + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PC[a1] * PC[b0] * PC[b1])
                + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PC[a1] * PC[b0] * PC[b2])
                + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PC[a1] * PC[b1] * PC[b2])
                + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PC[a2] * PC[b0] * PC[b1])
                + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[a2] * PC[b0] * PC[b2])
                + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[a2] * PC[b1] * PC[b2])
                + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PC[b0] * PC[b1] * PC[b2])
            )

            + 2.0 * (a_i + a_j) * (
                PA_0 * PA_1 * PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PA_2 * PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PB_0 * PC[a1] * PC[a2] * PC[b1] * PC[b2] * PC[m]
                + PA_0 * PB_1 * PC[a1] * PC[a2] * PC[b0] * PC[b2] * PC[m]
                + PA_0 * PB_2 * PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[m]
                + PA_1 * PA_2 * PC[a0] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PB_0 * PC[a0] * PC[a2] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PB_1 * PC[a0] * PC[a2] * PC[b0] * PC[b2] * PC[m]
                + PA_1 * PB_2 * PC[a0] * PC[a2] * PC[b0] * PC[b1] * PC[m]
                + PA_2 * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[b2] * PC[m]
                + PA_2 * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[b2] * PC[m]
                + PA_2 * PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
                + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[b2] * PC[m]
                + PB_0 * PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[b1] * PC[m]
                + PB_1 * PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[m]
            )

            + (
                delta[b2][m] * (PA_0 * PC[a1] * PC[a2] * PC[b0] * PC[b1] + PA_1 * PC[a0] * PC[a2] * PC[b0] * PC[b1] + PA_2 * PC[a0] * PC[a1] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[a2] * PC[b1] + PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[b0])
                + delta[b1][m] * (PA_0 * PC[a1] * PC[a2] * PC[b0] * PC[b2] + PA_1 * PC[a0] * PC[a2] * PC[b0] * PC[b2] + PA_2 * PC[a0] * PC[a1] * PC[b0] * PC[b2] + PB_0 * PC[a0] * PC[a1] * PC[a2] * PC[b2] + PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[b0])
                + delta[b0][m] * (PA_0 * PC[a1] * PC[a2] * PC[b1] * PC[b2] + PA_1 * PC[a0] * PC[a2] * PC[b1] * PC[b2] + PA_2 * PC[a0] * PC[a1] * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[b2] + PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[b1])
                + delta[a2][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[b2] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[b2] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[b2] + PB_2 * PC[a0] * PC[a1] * PC[b0] * PC[b1])
                + delta[a1][m] * (PA_0 * PC[a2] * PC[b0] * PC[b1] * PC[b2] + PA_2 * PC[a0] * PC[b0] * PC[b1] * PC[b2] + PB_0 * PC[a0] * PC[a2] * PC[b1] * PC[b2] + PB_1 * PC[a0] * PC[a2] * PC[b0] * PC[b2] + PB_2 * PC[a0] * PC[a2] * PC[b0] * PC[b1])
                + delta[a0][m] * (PA_1 * PC[a2] * PC[b0] * PC[b1] * PC[b2] + PA_2 * PC[a1] * PC[b0] * PC[b1] * PC[b2] + PB_0 * PC[a1] * PC[a2] * PC[b1] * PC[b2] + PB_1 * PC[a1] * PC[a2] * PC[b0] * PC[b2] + PB_2 * PC[a1] * PC[a2] * PC[b0] * PC[b1])
            )

        )

        + F7_t[6] * (

            (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                delta[b1][b2] * (PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[m])
                + delta[b0][b2] * (PC[a0] * PC[a1] * PC[a2] * PC[b1] * PC[m])
                + delta[b0][b1] * (PC[a0] * PC[a1] * PC[a2] * PC[b2] * PC[m])
                + delta[a2][b2] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                + delta[a2][b1] * (PC[a0] * PC[a1] * PC[b0] * PC[b2] * PC[m])
                + delta[a2][b0] * (PC[a0] * PC[a1] * PC[b1] * PC[b2] * PC[m])
                + delta[a1][b2] * (PC[a0] * PC[a2] * PC[b0] * PC[b1] * PC[m])
                + delta[a1][b1] * (PC[a0] * PC[a2] * PC[b0] * PC[b2] * PC[m])
                + delta[a1][b0] * (PC[a0] * PC[a2] * PC[b1] * PC[b2] * PC[m])
                + delta[a1][a2] * (PC[a0] * PC[b0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][b2] * (PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[m])
                + delta[a0][b1] * (PC[a1] * PC[a2] * PC[b0] * PC[b2] * PC[m])
                + delta[a0][b0] * (PC[a1] * PC[a2] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a2] * (PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m])
                + delta[a0][a1] * (PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m])
            )

            + (-2.0) * (a_i + a_j) * (
                PA_0 * PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_1 * PC[a0] * PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PA_2 * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[b2] * PC[m]
                + PB_0 * PC[a0] * PC[a1] * PC[a2] * PC[b1] * PC[b2] * PC[m]
                + PB_1 * PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[b2] * PC[m]
                + PB_2 * PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[m]
            )

            + (-1.0) * (
                delta[b2][m] * (PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[b1])
                + delta[b1][m] * (PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[b2])
                + delta[b0][m] * (PC[a0] * PC[a1] * PC[a2] * PC[b1] * PC[b2])
                + delta[a2][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[b2])
                + delta[a1][m] * (PC[a0] * PC[a2] * PC[b0] * PC[b1] * PC[b2])
                + delta[a0][m] * (PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[b2])
            )

        )

        + F7_t[7] * (

            2.0 * (a_i + a_j) * (
                PC[a0] * PC[a1] * PC[a2] * PC[b0] * PC[b1] * PC[b2] * PC[m]
            )

        )

    );
}

}  // namespace onee

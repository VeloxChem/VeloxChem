//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include "ScreeningData.hpp"

#include <omp.h>

#include <cmath>
#include <string>
#include <iostream>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

#define MATH_CONST_PI 3.14159265358979323846

CScreeningData::CScreeningData(const CMolecule& molecule, const CMolecularBasis& basis, const int64_t num_gpus_per_node)
{
    _num_gpus_per_node = num_gpus_per_node;

    _computeQMatrices(molecule, basis);

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    _sortQ(s_prim_count, p_prim_count, d_prim_count);

    form_mat_Q_and_density_inds_for_K(s_prim_count, p_prim_count, d_prim_count);
}

auto
CScreeningData::_computeQMatrices(const CMolecule& molecule, const CMolecularBasis& basis) -> void
{
    // Boys function

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    // GTO blocks

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // S block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    // P block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    // D block

    std::vector<double>   d_prim_info(5 * d_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    // GTO block pairs

    _Q_matrix_ss = CDenseMatrix(s_prim_count, s_prim_count);
    _Q_matrix_sp = CDenseMatrix(s_prim_count, p_prim_count * 3);
    _Q_matrix_sd = CDenseMatrix(s_prim_count, d_prim_count * 6);
    _Q_matrix_pp = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);
    _Q_matrix_pd = CDenseMatrix(p_prim_count * 3, d_prim_count * 6);
    _Q_matrix_dd = CDenseMatrix(d_prim_count * 6, d_prim_count * 6);

    // TODO distribute computation of Q matrices

    const double delta[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    const int64_t d_cart_ind[6][2] = {{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

    // S-S and S-P block pairs

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto a_j = s_prim_info[j + s_prim_count * 0];
            const auto c_j = s_prim_info[j + s_prim_count * 1];
            const auto x_j = s_prim_info[j + s_prim_count * 2];
            const auto y_j = s_prim_info[j + s_prim_count * 3];
            const auto z_j = s_prim_info[j + s_prim_count * 4];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto Lambda = std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto F0_t = boysfunc::getBoysFunction(0.0, 0, boys_func_table.data(), boys_func_ft.data());

            const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * F0_t[0]);

            _Q_matrix_ss.row(i)[j] = sqrt_eri_ijij;

            if (i != j) _Q_matrix_ss.row(j)[i] = sqrt_eri_ijij;
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            const auto a_j = p_prim_info[j + p_prim_count * 0];
            const auto c_j = p_prim_info[j + p_prim_count * 1];
            const auto x_j = p_prim_info[j + p_prim_count * 2];
            const auto y_j = p_prim_info[j + p_prim_count * 3];
            const auto z_j = p_prim_info[j + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto Lambda = std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

            const auto F1_t = boysfunc::getBoysFunction(0.0, 1, boys_func_table.data(), boys_func_ft.data());

            for (int64_t j_cart = 0; j_cart < 3; j_cart++)
            {
                // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
                // J. Chem. Phys. 84, 3963-3974 (1986)

                const auto b0 = j_cart;

                const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

                const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * (

                            F1_t[0] * PB_0 * PB_0 + F1_t[1] * 0.25 / a_ij

                            ));

                // TODO: think about the ordering of cartesian components

                _Q_matrix_sp.row(i)[j * 3 + j_cart] = sqrt_eri_ijij;
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            const auto a_j = d_prim_info[j + d_prim_count * 0];
            const auto c_j = d_prim_info[j + d_prim_count * 1];
            const auto x_j = d_prim_info[j + d_prim_count * 2];
            const auto y_j = d_prim_info[j + d_prim_count * 3];
            const auto z_j = d_prim_info[j + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_i + a_j) / (MATH_CONST_PI * (a_i + a_j + a_i + a_j)));

            const auto F2_t = boysfunc::getBoysFunction(0.0, 2, boys_func_table.data(), boys_func_ft.data());

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S1 = a_i + a_j;

            const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

            for (int64_t j_cart = 0; j_cart < 6; j_cart++)
            {
                const auto b0 = d_cart_ind[j_cart][0];
                const auto b1 = d_cart_ind[j_cart][1];

                const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
                const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

                const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * (

                            F2_t[0] * (

                                0.25 / ( S1 * S1 ) * (
                                    delta[b0][b1] * delta[b0][b1]
                                )

                                + 0.5 / S1 * (
                                    delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                                )

                                + (
                                    PB_0 * PB_0 * PB_1 * PB_1
                                )

                            )

                            + F2_t[1] * (

                                0.125 / ( S1 * S1 ) * (
                                    delta[b0][b1] * delta[b0][b1] * ((-2.0))
                                )

                                + 0.25 / S1 * (
                                    + (PB_0 * PB_0)
                                    + (PB_1 * PB_1)
                                )

                            )

                            + F2_t[2] * (

                                0.0625 / ( S1 * S1 ) * (
                                    1.0
                                    + delta[b0][b1] * delta[b0][b1] * 2.0
                                )

                            )

                            ));

                _Q_matrix_sd.row(i)[j * 6 + j_cart] = sqrt_eri_ijij;
            }
        }
    }

    // P-P gto block pair

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        const auto a_i = p_prim_info[i + p_prim_count * 0];
        const auto c_i = p_prim_info[i + p_prim_count * 1];
        const auto x_i = p_prim_info[i + p_prim_count * 2];
        const auto y_i = p_prim_info[i + p_prim_count * 3];
        const auto z_i = p_prim_info[i + p_prim_count * 4];

        for (int64_t j = i; j < p_prim_count; j++)
        {
            const auto a_j = p_prim_info[j + p_prim_count * 0];
            const auto c_j = p_prim_info[j + p_prim_count * 1];
            const auto x_j = p_prim_info[j + p_prim_count * 2];
            const auto y_j = p_prim_info[j + p_prim_count * 3];
            const auto z_j = p_prim_info[j + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_i + a_j) / (MATH_CONST_PI * (a_i + a_j + a_i + a_j)));

            const auto F2_t = boysfunc::getBoysFunction(0.0, 2, boys_func_table.data(), boys_func_ft.data());

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S1 = a_i + a_j;

            const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto a0 = i_cart;

                const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto b0 = j_cart;

                    const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];

                    const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * (

                                F2_t[0] * (

                                    0.5 / S1 * (
                                        delta[a0][b0] * (PB_0 * PA_0 * 2.0)
                                    )

                                    + (

                                        + PB_0 * PB_0 * PA_0 * PA_0
                                    )

                                    + 0.25 / ( S1 * S1 ) * (
                                        delta[a0][b0] * delta[a0][b0]
                                    )

                                )

                                + F2_t[1] * (

                                    0.125 / ( S1 * S1 ) * (
                                        delta[a0][b0] * delta[a0][b0] * ((-2.0))
                                    )

                                    + 0.25 / S1 * (
                                        (PA_0 * PA_0)
                                        + (PB_0 * PB_0)
                                    )

                                )

                                + F2_t[2] * (

                                    0.0625 / ( S1 * S1 ) * (
                                        1.0
                                        + delta[a0][b0] * delta[a0][b0] * 2.0
                                    )

                                )

                                ));

                    // TODO: think about the ordering of cartesian components

                    _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart] = sqrt_eri_ijij;

                    if (i * 3 + i_cart != j * 3 + j_cart) _Q_matrix_pp.row(j * 3 + j_cart)[i * 3 + i_cart] = sqrt_eri_ijij;
                }
            }
        }

        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            const auto a_j = d_prim_info[j + d_prim_count * 0];
            const auto c_j = d_prim_info[j + d_prim_count * 1];
            const auto x_j = d_prim_info[j + d_prim_count * 2];
            const auto y_j = d_prim_info[j + d_prim_count * 3];
            const auto z_j = d_prim_info[j + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto Lambda = std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            const auto F3_t = boysfunc::getBoysFunction(0.0, 3, boys_func_table.data(), boys_func_ft.data());

            const auto S1 = a_ij;

            const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto a0 = i_cart;

                const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto b0 = d_cart_ind[j_cart][0];
                    const auto b1 = d_cart_ind[j_cart][1];

                    const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
                    const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

                    const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * (

                                F3_t[0] * (

                                    0.25 / ( S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                                    )

                                    + 0.5 / S1 * (
                                        delta[b0][b1] * (PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                                        + delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                                        + delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                                    )

                                    + (

                                        + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0
                                    )

                                )

                                + F3_t[1] * (

                                    0.125 / ( S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][b0] * (PB_0 * PA_0 * 2.0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * (-2.0))
                                        + delta[a0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                                        + delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a0][b1] * delta[a0][b0] * (PB_0 * PB_1)
                                    )

                                    + 0.25 / S1 * (
                                        (PB_0 * PB_0 * PA_0 * PA_0)
                                        + (PB_1 * PB_1 * PA_0 * PA_0)
                                        + (PB_0 * PB_0 * PB_1 * PB_1)
                                    )

                                    + 0.0625 / ( S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1]
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 4.0
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1]
                                        + delta[a0][b0] * delta[a0][b0]
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                                        + delta[a0][b1] * delta[a0][b1]
                                    )

                                )

                                + F3_t[2] * (

                                    0.03125 / ( S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ((-8.0))
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b0] * delta[a0][b0] * ((-2.0))
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b1] * delta[a0][b1] * ((-2.0))
                                    )

                                    + 0.0625 / ( S1 * S1 ) * (
                                        (PA_0 * PA_0)
                                        + delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0 * (-1.0))
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0)
                                        + delta[a0][b1] * (PB_1 * PA_0)
                                        + delta[a0][b1] * (PB_1 * PA_0 * (-1.0))
                                        + (PB_0 * PB_0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * 2.0)
                                        + (PB_1 * PB_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * 2.0)
                                    )

                                )

                                + F3_t[3] * (

                                    0.015625 / ( S1 * S1 * S1 ) * (
                                        1.0
                                        + delta[b0][b1] * delta[b0][b1] * 2.0
                                        + delta[a0][b0] * delta[a0][b0] * 2.0
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 5.0
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                                        + delta[a0][b1] * delta[a0][b1]
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                                        + delta[a0][b1] * delta[a0][b1]
                                    )

                                )

                                ));

                    _Q_matrix_pd.row(i * 3 + i_cart)[j * 6 + j_cart] = sqrt_eri_ijij;
                }
            }
        }
    }

    // D-D gto block pair

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        const auto a_i = d_prim_info[i + d_prim_count * 0];
        const auto c_i = d_prim_info[i + d_prim_count * 1];
        const auto x_i = d_prim_info[i + d_prim_count * 2];
        const auto y_i = d_prim_info[i + d_prim_count * 3];
        const auto z_i = d_prim_info[i + d_prim_count * 4];

        for (int64_t j = i; j < d_prim_count; j++)
        {
            const auto a_j = d_prim_info[j + d_prim_count * 0];
            const auto c_j = d_prim_info[j + d_prim_count * 1];
            const auto x_j = d_prim_info[j + d_prim_count * 2];
            const auto y_j = d_prim_info[j + d_prim_count * 3];
            const auto z_j = d_prim_info[j + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_i + a_j) / (MATH_CONST_PI * (a_i + a_j + a_i + a_j)));

            const auto F4_t = boysfunc::getBoysFunction(0.0, 4, boys_func_table.data(), boys_func_ft.data());

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            const auto S1 = a_i + a_j;

            const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                const auto a0 = d_cart_ind[i_cart][0];
                const auto a1 = d_cart_ind[i_cart][1];

                const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
                const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto b0 = d_cart_ind[j_cart][0];
                    const auto b1 = d_cart_ind[j_cart][1];

                    const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
                    const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];

                    const auto sqrt_eri_ijij = std::sqrt(Lambda * S_ij_00 * S_ij_00 * (

                                F4_t[0] * (

                                    0.125 / ( S1 * S1 * S1 ) * (
                                        delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * (PA_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                                    )

                                    + 0.25 / ( S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * PA_1 * PA_1)
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[b0][b1] * (PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * PA_0 * PA_0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                                        + delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                                        + delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * PA_0 * PA_1 * 4.0)
                                        + delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                                        + delta[a0][a1] * delta[a1][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * PA_1 * PA_1)
                                        + delta[a0][a1] * delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1)
                                    )

                                    + 0.5 / S1 * (
                                        delta[b0][b1] * (PB_0 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a1][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                                    )

                                    + (

                                        + PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1
                                    )

                                    + 0.0625 / ( S1 * S1 * S1 * S1 ) * (
                                        delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1]
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 2.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 2.0
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1]
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 2.0
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1]
                                    )

                                )

                                + F4_t[1] * (

                                    0.03125 / ( S1 * S1 * S1 * S1 ) * (
                                        delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * ((-4.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * ((-8.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * ((-8.0))
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ((-4.0))
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * ((-8.0))
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * ((-4.0))
                                    )

                                    + 0.0625 / ( S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0)
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_0 * 4.0)
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * (PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b0] * (PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * (PA_0 * PA_0)
                                        + delta[a1][b1] * delta[a1][b1] * (PA_0 * PA_0)
                                        + delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_1 * (-4.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_1 * (-2.0))
                                        + delta[a0][b1] * delta[a1][b0] * delta[b0][b1] * (PA_0 * PA_1)
                                        + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * (PA_0 * PA_1 * (-2.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * (PA_0 * PA_1)
                                        + delta[a0][b0] * delta[a1][b0] * (PA_0 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * (PA_0 * PA_1)
                                        + delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * (PA_0 * PA_1)
                                        + delta[a0][b1] * delta[a1][b1] * (PA_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0 * (-2.0))
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (PB_0 * PA_0 * (-4.0))
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b1] * (PB_0 * PA_0)
                                        + delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][a1] * delta[a1][b0] * (PB_0 * PA_0 * 2.0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PA_0 * (-2.0))
                                        + delta[a1][b1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PA_0)
                                        + delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * (PB_1 * PA_0 * (-2.0))
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * (PB_1 * PA_0)
                                        + delta[a0][a1] * delta[a1][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a0][b0] * delta[a1][b1] * (PB_1 * PA_0)
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PA_0)
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_0 * (-5.0))
                                        + delta[b0][b1] * delta[b0][b1] * (PA_1 * PA_1)
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (PA_1 * PA_1 * 4.0)
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * (PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * (PA_1 * PA_1)
                                        + delta[a0][b1] * delta[a0][b1] * (PA_1 * PA_1)
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_1 * (-2.0))
                                        + delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PA_1 * (-2.0))
                                        + delta[a0][b1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PA_1)
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PA_1)
                                        + delta[a0][a1] * delta[a0][b0] * (PB_0 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_1)
                                        + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * (PB_0 * PA_1 * (-4.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b1] * (PB_0 * PA_1)
                                        + delta[a1][b0] * delta[b0][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * (PB_1 * PA_1)
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_1 * (-2.0))
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PA_1)
                                        + delta[a0][a1] * delta[a0][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * (PB_1 * PA_1)
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * (PB_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * (PB_1 * PA_1 * (-5.0))
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * (-4.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-2.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-1.0))
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0)
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * 4.0)
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0)
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                                        + delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_1)
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * (PB_0 * PB_1)
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1)
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1)
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * 4.0)
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1)
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                                    )

                                    + 0.125 / ( S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * PA_1 * PA_1 * (-2.0))
                                        + delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * (PB_0 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0 * PA_0 * PA_1 * (-2.0))
                                        + delta[a1][b1] * (PB_1 * PA_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * PA_0 * PA_0 * (-2.0))
                                        + delta[b0][b1] * (PB_0 * PB_1 * PA_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * PA_0 * PA_0 * (-1.0))
                                        + delta[a1][b1] * delta[a1][b0] * (PB_0 * PB_1 * PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * PA_0 * PA_0 * (-2.0))
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * (PB_0 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * PA_1 * PA_1 * (-2.0))
                                        + delta[a0][b1] * (PB_1 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * PA_0 * PA_1 * (-1.0))
                                        + delta[a0][a1] * (PB_0 * PB_0 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0 * PA_0 * PA_1)
                                        + delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * PA_0 * PA_1 * (-2.0))
                                        + delta[a0][b1] * delta[a1][b0] * (PB_0 * PB_1 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b1] * (PB_0 * PB_1 * PA_0 * PA_1 * (-2.0))
                                        + delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1 * PA_0 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * PA_0 * PA_1 * (-1.0))
                                        + delta[a0][a1] * (PB_1 * PB_1 * PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1 * PA_0 * PA_1)
                                        + delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * 2.0)
                                        + delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * PA_1 * PA_1 * (-2.0))
                                        + delta[b0][b1] * (PB_0 * PB_1 * PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * PA_1 * PA_1 * (-1.0))
                                        + delta[a0][b1] * delta[a0][b0] * (PB_0 * PB_1 * PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * PA_1 * PA_1 * (-2.0))
                                        + delta[a1][b1] * (PB_0 * PB_0 * PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_1 * (-1.0))
                                        + delta[a0][b1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PA_1)
                                        + delta[a1][b0] * (PB_0 * PB_1 * PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][a1] * (PB_0 * PB_1 * PB_1 * PA_1)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1 * (-2.0))
                                    )

                                    + 0.25 / S1 * (
                                        (PB_0 * PB_0 * PA_0 * PA_0 * PA_1 * PA_1)
                                        + (PB_1 * PB_1 * PA_0 * PA_0 * PA_1 * PA_1)
                                        + (PB_0 * PB_0 * PB_1 * PB_1 * PA_0 * PA_0)
                                        + (PB_0 * PB_0 * PB_1 * PB_1 * PA_1 * PA_1)
                                    )

                                )

                                + F4_t[2] * (

                                    0.03125 / ( S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_0 * (-8.0))
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a1][b0] * delta[a1][b0] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a1][b1] * delta[a1][b1] * (PA_0 * PA_0 * (-2.0))
                                        + delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_1 * 4.0)
                                        + delta[a0][b0] * delta[a1][b0] * (PA_0 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_1)
                                        + delta[a0][b1] * delta[a1][b1] * (PA_0 * PA_1)
                                        + delta[a0][a1] * (PA_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b0] * (PA_0 * PA_1)
                                        + delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * (PA_0 * PA_1)
                                        + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * (PA_0 * PA_1 * (-1.0))
                                        + delta[a1][b1] * delta[a0][b1] * (PA_0 * PA_1)
                                        + delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * (PA_0 * PA_1 * (-1.0))
                                        + delta[a0][b1] * delta[a1][b1] * (PA_0 * PA_1 * (-2.0))
                                        + delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_0 * (-1.0))
                                        + delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0)
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0)
                                        + delta[a0][b0] * (PB_0 * PA_0 * 2.0)
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * (-2.0))
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0)
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PA_0 * 2.0)
                                        + delta[a0][a1] * delta[a1][b1] * (PB_1 * PA_0 * (-1.0))
                                        + delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_0)
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * (PB_1 * PA_0 * (-2.0))
                                        + delta[a0][a1] * delta[a1][b1] * (PB_1 * PA_0)
                                        + delta[a0][b0] * delta[b0][b1] * (PB_1 * PA_0 * (-1.0))
                                        + delta[a0][b1] * (PB_1 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_0 * 6.0)
                                        + delta[b0][b1] * delta[b0][b1] * (PA_1 * PA_1 * (-2.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (PA_1 * PA_1 * (-8.0))
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * (PA_1 * PA_1 * (-2.0))
                                        + delta[a0][b0] * delta[a0][b0] * (PA_1 * PA_1 * (-2.0))
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * (PA_1 * PA_1 * (-2.0))
                                        + delta[a0][b1] * delta[a0][b1] * (PA_1 * PA_1 * (-2.0))
                                        + delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * (PB_0 * PA_1 * 2.0)
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_1 * (-2.0))
                                        + delta[a0][a1] * delta[a0][b0] * (PB_0 * PA_1 * (-1.0))
                                        + delta[a0][a1] * delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_1)
                                        + delta[a0][a1] * delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][a1] * (PB_0 * PA_1)
                                        + delta[a0][b1] * delta[a1][b0] * delta[a0][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * (PB_0 * PA_1 * 2.0)
                                        + delta[a1][b1] * (PB_1 * PA_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b1] * (PB_1 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * (PB_1 * PA_1 * 6.0)
                                        + delta[a0][b1] * delta[a0][a1] * (PB_1 * PA_1)
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * (PB_1 * PA_1 * (-2.0))
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * (-8.0))
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * (-2.0))
                                        + delta[a1][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-1.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1)
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * (PB_0 * PB_1 * 4.0)
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1 * (-1.0))
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1 * (-2.0))
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * (-2.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * (-8.0))
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1 * (-2.0))
                                        + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1 * (-2.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * (-2.0))
                                        + delta[b0][b1] * (PB_0 * PB_1 * 2.0)
                                        + delta[a1][b1] * delta[a1][b0] * (PB_0 * PB_1)
                                        + delta[a0][b0] * delta[a0][b1] * (PB_0 * PB_1)
                                        + delta[a0][b1] * delta[a0][b0] * (PB_0 * PB_1)
                                    )

                                    + 0.0625 / ( S1 * S1 ) * (
                                        (PA_0 * PA_0 * PA_1 * PA_1)
                                        + delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * PA_1 * PA_1 * 2.0)
                                        + delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_0 * PA_0 * PA_1 * (-1.0))
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_0 * PA_0 * PA_1)
                                        + delta[a1][b1] * (PB_1 * PA_0 * PA_0 * PA_1)
                                        + delta[a1][b1] * (PB_1 * PA_0 * PA_0 * PA_1 * (-1.0))
                                        + (PB_0 * PB_0 * PA_0 * PA_0)
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * PA_0 * PA_0 * 2.0)
                                        + (PB_1 * PB_1 * PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * PA_0 * PA_0 * 2.0)
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0 * PA_1 * PA_1 * (-1.0))
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0 * PA_1 * PA_1)
                                        + delta[a0][b1] * (PB_1 * PA_0 * PA_1 * PA_1)
                                        + delta[a0][b1] * (PB_1 * PA_0 * PA_1 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * PA_0 * PA_1)
                                        + delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1 * PA_0 * PA_1 * (-1.0))
                                        + delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_0)
                                        + delta[a0][b1] * (PB_0 * PB_0 * PB_1 * PA_0 * (-1.0))
                                        + delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_0)
                                        + delta[a0][b0] * (PB_0 * PB_1 * PB_1 * PA_0 * (-1.0))
                                        + (PB_0 * PB_0 * PA_1 * PA_1)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0 * PA_1 * PA_1 * 2.0)
                                        + (PB_1 * PB_1 * PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1 * PA_1 * PA_1 * 2.0)
                                        + (PB_0 * PB_0 * PB_1 * PB_1)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * PB_1 * PB_1 * 2.0)
                                    )

                                    + 0.015625 / ( S1 * S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1]
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * 4.0
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1]
                                        + delta[a1][b0] * delta[a1][b0]
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1]
                                        + delta[a1][b1] * delta[a1][b1]
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * 8.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 19.0
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * 4.0
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 17.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * 2.0
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * 4.0
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1]
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0]
                                        + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * 2.0
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1]
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * 7.0
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 15.0
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * 3.0
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1]
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1]
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * 3.0
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1]
                                        + delta[a0][a1] * delta[a0][a1]
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0]
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1]
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * 2.0
                                        + delta[a0][b0] * delta[a0][b0]
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 3.0
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                                        + delta[a0][b1] * delta[a0][b1]
                                        + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * 3.0
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * 7.0
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1]
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * 3.0
                                    )

                                )

                                + F4_t[3] * (

                                    0.0078125 / ( S1 * S1 * S1 * S1 ) * (
                                        delta[b0][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * ((-8.0))
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * ((-2.0))
                                        + delta[a1][b0] * delta[a1][b0] * ((-2.0))
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a1][b1] * delta[a1][b1] * ((-2.0))
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * ((-8.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * ((-22.0))
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * ((-4.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * ((-8.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * ((-18.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * ((-4.0))
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * ((-8.0))
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * ((-2.0))
                                        + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * ((-4.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ((-2.0))
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * ((-6.0))
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * ((-14.0))
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * ((-6.0))
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * ((-2.0))
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * ((-6.0))
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * ((-2.0))
                                        + delta[a0][a1] * delta[a0][a1] * ((-2.0))
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * ((-2.0))
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * ((-2.0))
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * ((-4.0))
                                        + delta[a0][b0] * delta[a0][b0] * ((-2.0))
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * ((-6.0))
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * ((-2.0))
                                        + delta[a0][b1] * delta[a0][b1] * ((-2.0))
                                        + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * ((-6.0))
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * ((-6.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * ((-2.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * ((-6.0))
                                    )

                                    + 0.015625 / ( S1 * S1 * S1 ) * (
                                        (PA_0 * PA_0)
                                        + delta[b0][b1] * delta[b0][b1] * (PA_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[a1][b0] * (PA_0 * PA_0 * 2.0)
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * (PA_0 * PA_0 * 5.0)
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * (PA_0 * PA_0 * 2.0)
                                        + delta[a1][b1] * delta[a1][b1] * (PA_0 * PA_0)
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1] * (PA_0 * PA_0)
                                        + delta[a1][b1] * delta[a1][b1] * (PA_0 * PA_0)
                                        + delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * (PA_0 * PA_1 * (-1.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * (PA_0 * PA_1)
                                        + delta[a0][b0] * (PB_0 * PA_0)
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * (PB_0 * PA_0)
                                        + delta[a0][b0] * (PB_0 * PA_0 * (-1.0))
                                        + delta[b0][b1] * delta[a0][b1] * (PB_0 * PA_0)
                                        + delta[a0][b1] * delta[b0][b1] * (PB_0 * PA_0 * (-1.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PA_0 * (-1.0))
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b0] * (PB_1 * PA_0)
                                        + delta[a0][b1] * (PB_1 * PA_0)
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * (PB_1 * PA_0 * 2.0)
                                        + delta[a0][b1] * (PB_1 * PA_0 * (-1.0))
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_0 * (-2.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PA_0 * (-1.0))
                                        + (PA_1 * PA_1)
                                        + delta[b0][b1] * delta[b0][b1] * (PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b0] * (PA_1 * PA_1 * 2.0)
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * (PA_1 * PA_1 * 5.0)
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * (PA_1 * PA_1 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PA_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1] * (PA_1 * PA_1)
                                        + delta[a0][b1] * delta[a0][b1] * (PA_1 * PA_1)
                                        + delta[a1][b1] * delta[b0][b1] * (PB_0 * PA_1 * (-1.0))
                                        + delta[b0][b1] * delta[a1][b1] * (PB_0 * PA_1)
                                        + delta[a0][b0] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PA_1)
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b1] * (PB_0 * PA_1 * (-1.0))
                                        + delta[a1][b1] * (PB_1 * PA_1)
                                        + delta[a1][b1] * (PB_1 * PA_1 * (-1.0))
                                        + delta[a0][b0] * delta[a1][b1] * delta[a0][b0] * (PB_1 * PA_1)
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * (PB_1 * PA_1 * (-2.0))
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * (PB_1 * PA_1)
                                        + delta[a0][b1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PA_1)
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PA_1 * (-1.0))
                                        + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * (PB_1 * PA_1 * (-1.0))
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * (PB_1 * PA_1)
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * (PB_0 * PB_1)
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * (PB_0 * PB_1 * (-1.0))
                                        + (PB_0 * PB_0)
                                        + delta[a1][b1] * delta[a1][b1] * (PB_0 * PB_0 * 2.0)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_0 * PB_0 * 2.0)
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * (PB_0 * PB_0 * 2.0)
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * (PB_0 * PB_0 * 4.0)
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * (PB_0 * PB_0 * 2.0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                                        + delta[a0][b1] * delta[a0][b1] * (PB_0 * PB_0)
                                        + (PB_1 * PB_1)
                                        + delta[a1][b0] * delta[a1][b0] * (PB_1 * PB_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][a1] * (PB_1 * PB_1 * 2.0)
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * (PB_1 * PB_1 * 2.0)
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * (PB_1 * PB_1 * 4.0)
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * (PB_1 * PB_1 * 2.0)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                                        + delta[a0][b0] * delta[a0][b0] * (PB_1 * PB_1)
                                    )

                                )

                                + F4_t[4] * (

                                    0.00390625 / ( S1 * S1 * S1 * S1 ) * (
                                        1.0
                                        + delta[b0][b1] * delta[b0][b1] * 2.0
                                        + delta[a1][b0] * delta[a1][b0] * 2.0
                                        + delta[a1][b0] * delta[b0][b1] * delta[a1][b1] * 5.0
                                        + delta[a1][b1] * delta[a1][b0] * delta[b0][b1] * 2.0
                                        + delta[a1][b1] * delta[a1][b1]
                                        + delta[a1][b0] * delta[a1][b1] * delta[b0][b1]
                                        + delta[a1][b1] * delta[a1][b1]
                                        + delta[a0][a1] * delta[a0][a1] * 2.0
                                        + delta[a0][a1] * delta[a0][a1] * delta[b0][b1] * delta[b0][b1] * 4.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b0] * 2.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * delta[a0][b1] * 10.0
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b0] * delta[b0][b1] * 3.0
                                        + delta[a0][a1] * delta[a1][b1] * delta[a0][b1] * 2.0
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b0] * 4.0
                                        + delta[a0][a1] * delta[a0][b0] * delta[b0][b1] * delta[a1][b1] * 7.0
                                        + delta[a0][a1] * delta[a1][b0] * delta[a0][b1] * delta[b0][b1] * 2.0
                                        + delta[a0][a1] * delta[a0][b1] * delta[a1][b1] * 4.0
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b1] * delta[b0][b1]
                                        + delta[a0][b0] * delta[a0][a1] * delta[a1][b0] * 2.0
                                        + delta[a0][b0] * delta[a0][a1] * delta[b0][b1] * delta[a1][b1] * 3.0
                                        + delta[a0][b0] * delta[a0][b0]
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 2.0
                                        + delta[a0][b0] * delta[a1][b1] * delta[a0][b0] * delta[a1][b1]
                                        + delta[a0][b0] * delta[a1][b1] * delta[a1][b0] * delta[a0][b1]
                                        + delta[a0][b0] * delta[a0][b0] * delta[a1][b1] * delta[a1][b1] * 2.0
                                        + delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * delta[a1][b1] * 5.0
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b0] * delta[b0][b1] * 4.0
                                        + delta[a0][b1] * delta[a0][a1] * delta[a1][b1] * 2.0
                                        + delta[a0][b1] * delta[a0][b0] * delta[b0][b1] * 2.0
                                        + delta[a0][b1] * delta[a0][b1]
                                        + delta[a0][b1] * delta[a1][b0] * delta[a0][b0] * delta[a1][b1]
                                        + delta[a0][b1] * delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * 2.0
                                        + delta[a0][b1] * delta[a0][b0] * delta[a1][b0] * delta[a1][b1] * 3.0
                                        + delta[a0][a1] * delta[a0][b0] * delta[a1][b1] * delta[b0][b1] * 2.0
                                        + delta[a0][b0] * delta[a0][b0]
                                        + delta[a0][b0] * delta[b0][b1] * delta[a0][b1] * 3.0
                                        + delta[a0][b0] * delta[a0][b1] * delta[b0][b1]
                                        + delta[a0][b1] * delta[a0][b1]
                                        + delta[a1][b0] * delta[a0][b0] * delta[a0][b1] * delta[a1][b1] * 3.0
                                        + delta[a1][b0] * delta[a1][b0] * delta[a0][b1] * delta[a0][b1] * 2.0
                                        + delta[a1][b1] * delta[a0][b0] * delta[a0][b0] * delta[a1][b1]
                                        + delta[a1][b1] * delta[a0][b0] * delta[a1][b0] * delta[a0][b1] * 3.0
                                    )

                                )

                                ));

                    _Q_matrix_dd.row(i * 6 + i_cart)[j * 6 + j_cart] = sqrt_eri_ijij;

                    if (i * 6 + i_cart != j * 6 + j_cart) _Q_matrix_dd.row(j * 6 + j_cart)[i * 6 + i_cart] = sqrt_eri_ijij;
                }
            }
        }
    }
}

auto
CScreeningData::getNumGpusPerNode() const -> const int64_t
{
    return _num_gpus_per_node;
}

auto
CScreeningData::getQMatrixSS() const -> const CDenseMatrix&
{
    return _Q_matrix_ss;
}

auto
CScreeningData::getQMatrixSP() const -> const CDenseMatrix&
{
    return _Q_matrix_sp;
}

auto
CScreeningData::getQMatrixSD() const -> const CDenseMatrix&
{
    return _Q_matrix_sd;
}

auto
CScreeningData::getQMatrixPP() const -> const CDenseMatrix&
{
    return _Q_matrix_pp;
}

auto
CScreeningData::getQMatrixPD() const -> const CDenseMatrix&
{
    return _Q_matrix_pd;
}

auto
CScreeningData::getQMatrixDD() const -> const CDenseMatrix&
{
    return _Q_matrix_dd;
}

auto
CScreeningData::_sortQ(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) -> void
{
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_ss_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_sp_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_sd_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_pp_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_pd_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_dd_mat_Q;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            sorted_ss_mat_Q.push_back(std::make_tuple(Q_ij, i, j));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + s];
                sorted_sp_mat_Q.push_back(std::make_tuple(Q_ij, i, j * 3 + s));
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t s = 0; s < 6; s++)
            {
                const auto Q_ij = _Q_matrix_sd.row(i)[j * 6 + s];
                sorted_sd_mat_Q.push_back(std::make_tuple(Q_ij, i, j * 6 + s));
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    sorted_pp_mat_Q.push_back(std::make_tuple(Q_ij, i * 3 + i_cart, j * 3 + j_cart));
                }
            }
        }

        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pd.row(i * 3 + i_cart)[j * 6 + j_cart];
                    sorted_pd_mat_Q.push_back(std::make_tuple(Q_ij, i * 3 + i_cart, j * 6 + j_cart));
                }
            }
        }
    }

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_dd.row(i * 6 + i_cart)[j * 6 + j_cart];
                    sorted_dd_mat_Q.push_back(std::make_tuple(Q_ij, i * 6 + i_cart, j * 6 + j_cart));
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q.begin(), sorted_ss_mat_Q.end());
    std::sort(sorted_sp_mat_Q.begin(), sorted_sp_mat_Q.end());
    std::sort(sorted_sd_mat_Q.begin(), sorted_sd_mat_Q.end());
    std::sort(sorted_pp_mat_Q.begin(), sorted_pp_mat_Q.end());
    std::sort(sorted_pd_mat_Q.begin(), sorted_pd_mat_Q.end());
    std::sort(sorted_dd_mat_Q.begin(), sorted_dd_mat_Q.end());

    const auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    const auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    const auto sd_prim_pair_count = s_prim_count * d_prim_count * 6;
    const auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;
    const auto pd_prim_pair_count = p_prim_count * 3 * d_prim_count * 6;
    const auto dd_prim_pair_count = d_prim_count * 6 * (d_prim_count * 6 + 1) / 2;

    // form local vectors

    _ss_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _ss_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _ss_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _sp_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sp_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sp_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _sd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _pp_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pp_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pp_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _pd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _dd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _dd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _dd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
        auto gpu_count = nnodes * _num_gpus_per_node;

        auto ss_batch_size = mathfunc::batch_size(ss_prim_pair_count, gpu_rank, gpu_count);

        _ss_first_inds_local[gpu_id]  = std::vector<uint32_t>(ss_batch_size);
        _ss_second_inds_local[gpu_id] = std::vector<uint32_t>(ss_batch_size);
        _ss_mat_Q_local[gpu_id]       = std::vector<double>(ss_batch_size);

        auto sp_batch_size = mathfunc::batch_size(sp_prim_pair_count, gpu_rank, gpu_count);

        _sp_first_inds_local[gpu_id]  = std::vector<uint32_t>(sp_batch_size);
        _sp_second_inds_local[gpu_id] = std::vector<uint32_t>(sp_batch_size);
        _sp_mat_Q_local[gpu_id]       = std::vector<double>(sp_batch_size);

        auto sd_batch_size = mathfunc::batch_size(sd_prim_pair_count, gpu_rank, gpu_count);

        _sd_first_inds_local[gpu_id]  = std::vector<uint32_t>(sd_batch_size);
        _sd_second_inds_local[gpu_id] = std::vector<uint32_t>(sd_batch_size);
        _sd_mat_Q_local[gpu_id]       = std::vector<double>(sd_batch_size);

        auto pp_batch_size = mathfunc::batch_size(pp_prim_pair_count, gpu_rank, gpu_count);

        _pp_first_inds_local[gpu_id]  = std::vector<uint32_t>(pp_batch_size);
        _pp_second_inds_local[gpu_id] = std::vector<uint32_t>(pp_batch_size);
        _pp_mat_Q_local[gpu_id]       = std::vector<double>(pp_batch_size);

        auto pd_batch_size = mathfunc::batch_size(pd_prim_pair_count, gpu_rank, gpu_count);

        _pd_first_inds_local[gpu_id]  = std::vector<uint32_t>(pd_batch_size);
        _pd_second_inds_local[gpu_id] = std::vector<uint32_t>(pd_batch_size);
        _pd_mat_Q_local[gpu_id]       = std::vector<double>(pd_batch_size);

        auto dd_batch_size = mathfunc::batch_size(dd_prim_pair_count, gpu_rank, gpu_count);

        _dd_first_inds_local[gpu_id]  = std::vector<uint32_t>(dd_batch_size);
        _dd_second_inds_local[gpu_id] = std::vector<uint32_t>(dd_batch_size);
        _dd_mat_Q_local[gpu_id]       = std::vector<double>(dd_batch_size);
    }

    auto nthreads = omp_get_max_threads();
    auto num_threads_per_gpu = nthreads / _num_gpus_per_node;

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        if (thread_id % num_threads_per_gpu == 0)
        {
            auto gpu_id = thread_id / num_threads_per_gpu;
            auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
            auto gpu_count = nnodes * _num_gpus_per_node;

            for (int64_t ij = gpu_rank, idx = 0; ij < ss_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_ss_mat_Q[ss_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _ss_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _ss_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _ss_mat_Q_local[gpu_id][idx]       = Q_ij;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < sp_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_sp_mat_Q[sp_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _sp_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _sp_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _sp_mat_Q_local[gpu_id][idx]       = Q_ij;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < sd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_sd_mat_Q[sd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _sd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _sd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _sd_mat_Q_local[gpu_id][idx]       = Q_ij;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < pp_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_pp_mat_Q[pp_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _pp_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _pp_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _pp_mat_Q_local[gpu_id][idx]       = Q_ij;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < pd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_pd_mat_Q[pd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _pd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _pd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _pd_mat_Q_local[gpu_id][idx]       = Q_ij;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < dd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_dd_mat_Q[dd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _dd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _dd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _dd_mat_Q_local[gpu_id][idx]       = Q_ij;
            }
        }
    }
}

auto
CScreeningData::sortQD(const int64_t s_prim_count,
                       const int64_t p_prim_count,
                       const int64_t d_prim_count,
                       const std::vector<uint32_t>& s_prim_aoinds,
                       const std::vector<uint32_t>& p_prim_aoinds,
                       const std::vector<uint32_t>& d_prim_aoinds,
                       const int64_t naos,
                       const double* dens_ptr) -> void
{
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_ss_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_sp_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_sd_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_pp_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_pd_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_dd_mat_Q_D;

    auto max_Q_ss_ptr = std::max_element(_Q_matrix_ss.values(), _Q_matrix_ss.values() + _Q_matrix_ss.getNumberOfElements());
    double max_Q = *max_Q_ss_ptr;

    if (p_prim_count > 0)
    {
        auto max_Q_sp_ptr = std::max_element(_Q_matrix_sp.values(), _Q_matrix_sp.values() + _Q_matrix_sp.getNumberOfElements());
        auto max_Q_pp_ptr = std::max_element(_Q_matrix_pp.values(), _Q_matrix_pp.values() + _Q_matrix_pp.getNumberOfElements());
        max_Q = std::max({max_Q, *max_Q_sp_ptr, *max_Q_pp_ptr});
    }

    if (d_prim_count > 0)
    {
        auto max_Q_sd_ptr = std::max_element(_Q_matrix_sd.values(), _Q_matrix_sd.values() + _Q_matrix_sd.getNumberOfElements());
        auto max_Q_pd_ptr = std::max_element(_Q_matrix_pd.values(), _Q_matrix_pd.values() + _Q_matrix_pd.getNumberOfElements());
        auto max_Q_dd_ptr = std::max_element(_Q_matrix_dd.values(), _Q_matrix_dd.values() + _Q_matrix_dd.getNumberOfElements());
        max_Q = std::max({max_Q, *max_Q_sd_ptr, *max_Q_pd_ptr, *max_Q_dd_ptr});
    }

    // TODO use J threshold from arguments
    const double thresh_QD = 1.0e-13 / max_Q;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

            if (QD_abs_ij > thresh_QD) sorted_ss_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j, Q_ij, D_ij));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 3; j_cart++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + j_cart];
                const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                if (QD_abs_ij > thresh_QD) sorted_sp_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j * 3 + j_cart, Q_ij, D_ij));
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 6; j_cart++)
            {
                const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                const auto Q_ij = _Q_matrix_sd.row(i)[j * 6 + j_cart];
                const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                if (QD_abs_ij > thresh_QD) sorted_sd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j * 6 + j_cart, Q_ij, D_ij));
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if (QD_abs_ij > thresh_QD) sorted_pp_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 3 + i_cart, j * 3 + j_cart, Q_ij, D_ij));
                }
            }
        }
    
        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_pd.row(i * 3 + i_cart)[j * 6 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if (QD_abs_ij > thresh_QD) sorted_pd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 3 + i_cart, j * 6 + j_cart, Q_ij, D_ij));
                }
            }
        }    
    }

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                const auto i_cgto = d_prim_aoinds[i + d_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_dd.row(i * 6 + i_cart)[j * 6 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if (QD_abs_ij > thresh_QD) sorted_dd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 6 + i_cart, j * 6 + j_cart, Q_ij, D_ij));
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q_D.begin(), sorted_ss_mat_Q_D.end());
    std::sort(sorted_sp_mat_Q_D.begin(), sorted_sp_mat_Q_D.end());
    std::sort(sorted_sd_mat_Q_D.begin(), sorted_sd_mat_Q_D.end());
    std::sort(sorted_pp_mat_Q_D.begin(), sorted_pp_mat_Q_D.end());
    std::sort(sorted_pd_mat_Q_D.begin(), sorted_pd_mat_Q_D.end());
    std::sort(sorted_dd_mat_Q_D.begin(), sorted_dd_mat_Q_D.end());

    const auto ss_prim_pair_count = static_cast<int64_t>(sorted_ss_mat_Q_D.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sorted_sp_mat_Q_D.size());
    const auto sd_prim_pair_count = static_cast<int64_t>(sorted_sd_mat_Q_D.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(sorted_pp_mat_Q_D.size());
    const auto pd_prim_pair_count = static_cast<int64_t>(sorted_pd_mat_Q_D.size());
    const auto dd_prim_pair_count = static_cast<int64_t>(sorted_dd_mat_Q_D.size());

    // const auto ss_prim_pair_total = s_prim_count * (s_prim_count + 1) / 2;
    // const auto sp_prim_pair_total = s_prim_count * p_prim_count * 3;
    // const auto pp_prim_pair_total = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    // const auto prim_pair_count = ss_prim_pair_count + sp_prim_pair_count + pp_prim_pair_count;
    // const auto prim_pair_total = ss_prim_pair_total + sp_prim_pair_total + pp_prim_pair_total;

    // std::cout << "\nJ screening total ket prim. pairs: " << prim_pair_count << "\n";
    // std::cout << "J screening valid ket prim. pairs: " << prim_pair_total << "\n";
    // std::cout << "J screening ratio on the ket side: " << 1.0 - static_cast<double>(prim_pair_count) / prim_pair_total << "\n";

    _ss_max_D = 0.0;
    _sp_max_D = 0.0;
    _sd_max_D = 0.0;
    _pp_max_D = 0.0;
    _pd_max_D = 0.0;
    _dd_max_D = 0.0;

    _ss_mat_Q       = std::vector<double>  (ss_prim_pair_count);
    _ss_mat_D       = std::vector<double>  (ss_prim_pair_count);
    _ss_first_inds  = std::vector<uint32_t>(ss_prim_pair_count);
    _ss_second_inds = std::vector<uint32_t>(ss_prim_pair_count);

    for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto& vals = sorted_ss_mat_Q_D[ss_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _ss_mat_Q[ij] = Q_ij;
        _ss_mat_D[ij] = D_ij;

        _ss_first_inds[ij]  = static_cast<uint32_t>(i);
        _ss_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _ss_max_D) _ss_max_D = std::fabs(D_ij);
    }

    _sp_mat_Q       = std::vector<double>  (sp_prim_pair_count);
    _sp_mat_D       = std::vector<double>  (sp_prim_pair_count);
    _sp_first_inds  = std::vector<uint32_t>(sp_prim_pair_count);
    _sp_second_inds = std::vector<uint32_t>(sp_prim_pair_count);

    for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto& vals = sorted_sp_mat_Q_D[sp_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _sp_mat_Q[ij] = Q_ij;
        _sp_mat_D[ij] = D_ij;

        _sp_first_inds[ij]  = static_cast<uint32_t>(i);
        _sp_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _sp_max_D) _sp_max_D = std::fabs(D_ij);
    }

    _sd_mat_Q       = std::vector<double>  (sd_prim_pair_count);
    _sd_mat_D       = std::vector<double>  (sd_prim_pair_count);
    _sd_first_inds  = std::vector<uint32_t>(sd_prim_pair_count);
    _sd_second_inds = std::vector<uint32_t>(sd_prim_pair_count);

    for (int64_t ij = 0; ij < sd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_sd_mat_Q_D[sd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _sd_mat_Q[ij] = Q_ij;
        _sd_mat_D[ij] = D_ij;

        _sd_first_inds[ij]  = static_cast<uint32_t>(i);
        _sd_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _sd_max_D) _sd_max_D = std::fabs(D_ij);
    }

    _pp_mat_Q       = std::vector<double>  (pp_prim_pair_count);
    _pp_mat_D       = std::vector<double>  (pp_prim_pair_count);
    _pp_first_inds  = std::vector<uint32_t>(pp_prim_pair_count);
    _pp_second_inds = std::vector<uint32_t>(pp_prim_pair_count);

    for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto& vals = sorted_pp_mat_Q_D[pp_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _pp_mat_Q[ij] = Q_ij;
        _pp_mat_D[ij] = D_ij;

        _pp_first_inds[ij]  = static_cast<uint32_t>(i);
        _pp_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _pp_max_D) _pp_max_D = std::fabs(D_ij);
    }

    _pd_mat_Q       = std::vector<double>  (pd_prim_pair_count);
    _pd_mat_D       = std::vector<double>  (pd_prim_pair_count);
    _pd_first_inds  = std::vector<uint32_t>(pd_prim_pair_count);
    _pd_second_inds = std::vector<uint32_t>(pd_prim_pair_count);

    for (int64_t ij = 0; ij < pd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_pd_mat_Q_D[pd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _pd_mat_Q[ij] = Q_ij;
        _pd_mat_D[ij] = D_ij;

        _pd_first_inds[ij]  = static_cast<uint32_t>(i);
        _pd_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _pd_max_D) _pd_max_D = std::fabs(D_ij);
    }

    _dd_mat_Q       = std::vector<double>  (dd_prim_pair_count);
    _dd_mat_D       = std::vector<double>  (dd_prim_pair_count);
    _dd_first_inds  = std::vector<uint32_t>(dd_prim_pair_count);
    _dd_second_inds = std::vector<uint32_t>(dd_prim_pair_count);

    for (int64_t ij = 0; ij < dd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_dd_mat_Q_D[dd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _dd_mat_Q[ij] = Q_ij;
        _dd_mat_D[ij] = D_ij;

        _dd_first_inds[ij]  = static_cast<uint32_t>(i);
        _dd_second_inds[ij] = static_cast<uint32_t>(j);

        if (std::fabs(D_ij) > _dd_max_D) _dd_max_D = std::fabs(D_ij);
    }
}

auto CScreeningData::get_ss_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _ss_first_inds_local[gpu_id]; }
auto CScreeningData::get_sp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sp_first_inds_local[gpu_id]; }
auto CScreeningData::get_sd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sd_first_inds_local[gpu_id]; }
auto CScreeningData::get_pp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pp_first_inds_local[gpu_id]; }
auto CScreeningData::get_pd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pd_first_inds_local[gpu_id]; }
auto CScreeningData::get_dd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _dd_first_inds_local[gpu_id]; }

auto CScreeningData::get_ss_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _ss_second_inds_local[gpu_id]; }
auto CScreeningData::get_sp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sp_second_inds_local[gpu_id]; }
auto CScreeningData::get_sd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sd_second_inds_local[gpu_id]; }
auto CScreeningData::get_pp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pp_second_inds_local[gpu_id]; }
auto CScreeningData::get_pd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pd_second_inds_local[gpu_id]; }
auto CScreeningData::get_dd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _dd_second_inds_local[gpu_id]; }

auto CScreeningData::get_ss_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _ss_mat_Q_local[gpu_id]; }
auto CScreeningData::get_sp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sp_mat_Q_local[gpu_id]; }
auto CScreeningData::get_sd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sd_mat_Q_local[gpu_id]; }
auto CScreeningData::get_pp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pp_mat_Q_local[gpu_id]; }
auto CScreeningData::get_pd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pd_mat_Q_local[gpu_id]; }
auto CScreeningData::get_dd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _dd_mat_Q_local[gpu_id]; }

auto CScreeningData::get_ss_first_inds() const -> const std::vector<uint32_t>& { return _ss_first_inds; }
auto CScreeningData::get_sp_first_inds() const -> const std::vector<uint32_t>& { return _sp_first_inds; }
auto CScreeningData::get_sd_first_inds() const -> const std::vector<uint32_t>& { return _sd_first_inds; }
auto CScreeningData::get_pp_first_inds() const -> const std::vector<uint32_t>& { return _pp_first_inds; }
auto CScreeningData::get_pd_first_inds() const -> const std::vector<uint32_t>& { return _pd_first_inds; }
auto CScreeningData::get_dd_first_inds() const -> const std::vector<uint32_t>& { return _dd_first_inds; }

auto CScreeningData::get_ss_second_inds() const -> const std::vector<uint32_t>& { return _ss_second_inds; }
auto CScreeningData::get_sp_second_inds() const -> const std::vector<uint32_t>& { return _sp_second_inds; }
auto CScreeningData::get_sd_second_inds() const -> const std::vector<uint32_t>& { return _sd_second_inds; }
auto CScreeningData::get_pp_second_inds() const -> const std::vector<uint32_t>& { return _pp_second_inds; }
auto CScreeningData::get_pd_second_inds() const -> const std::vector<uint32_t>& { return _pd_second_inds; }
auto CScreeningData::get_dd_second_inds() const -> const std::vector<uint32_t>& { return _dd_second_inds; }

auto CScreeningData::get_ss_mat_Q() const -> const std::vector<double>& { return _ss_mat_Q; }
auto CScreeningData::get_sp_mat_Q() const -> const std::vector<double>& { return _sp_mat_Q; }
auto CScreeningData::get_sd_mat_Q() const -> const std::vector<double>& { return _sd_mat_Q; }
auto CScreeningData::get_pp_mat_Q() const -> const std::vector<double>& { return _pp_mat_Q; }
auto CScreeningData::get_pd_mat_Q() const -> const std::vector<double>& { return _pd_mat_Q; }
auto CScreeningData::get_dd_mat_Q() const -> const std::vector<double>& { return _dd_mat_Q; }

auto CScreeningData::get_ss_mat_D() const -> const std::vector<double>& { return _ss_mat_D; }
auto CScreeningData::get_sp_mat_D() const -> const std::vector<double>& { return _sp_mat_D; }
auto CScreeningData::get_sd_mat_D() const -> const std::vector<double>& { return _sd_mat_D; }
auto CScreeningData::get_pp_mat_D() const -> const std::vector<double>& { return _pp_mat_D; }
auto CScreeningData::get_pd_mat_D() const -> const std::vector<double>& { return _pd_mat_D; }
auto CScreeningData::get_dd_mat_D() const -> const std::vector<double>& { return _dd_mat_D; }

auto CScreeningData::get_ss_max_D() const -> double { return _ss_max_D; }
auto CScreeningData::get_sp_max_D() const -> double { return _sp_max_D; }
auto CScreeningData::get_sd_max_D() const -> double { return _sd_max_D; }
auto CScreeningData::get_pp_max_D() const -> double { return _pp_max_D; }
auto CScreeningData::get_pd_max_D() const -> double { return _pd_max_D; }
auto CScreeningData::get_dd_max_D() const -> double { return _dd_max_D; }

auto CScreeningData::get_mat_Q_for_K_ss() const -> const std::vector<double>& { return _mat_Q_for_K_ss; };
auto CScreeningData::get_mat_Q_for_K_sp() const -> const std::vector<double>& { return _mat_Q_for_K_sp; };
auto CScreeningData::get_mat_Q_for_K_ps() const -> const std::vector<double>& { return _mat_Q_for_K_ps; };
auto CScreeningData::get_mat_Q_for_K_sd() const -> const std::vector<double>& { return _mat_Q_for_K_sd; };
auto CScreeningData::get_mat_Q_for_K_ds() const -> const std::vector<double>& { return _mat_Q_for_K_ds; };
auto CScreeningData::get_mat_Q_for_K_pp() const -> const std::vector<double>& { return _mat_Q_for_K_pp; };
auto CScreeningData::get_mat_Q_for_K_pd() const -> const std::vector<double>& { return _mat_Q_for_K_pd; };
auto CScreeningData::get_mat_Q_for_K_dp() const -> const std::vector<double>& { return _mat_Q_for_K_dp; };
auto CScreeningData::get_mat_Q_for_K_dd() const -> const std::vector<double>& { return _mat_Q_for_K_dd; };

auto CScreeningData::get_density_inds_for_K_ss() const -> const std::vector<uint32_t>& { return _density_inds_for_K_ss; }
auto CScreeningData::get_density_inds_for_K_sp() const -> const std::vector<uint32_t>& { return _density_inds_for_K_sp; }
auto CScreeningData::get_density_inds_for_K_ps() const -> const std::vector<uint32_t>& { return _density_inds_for_K_ps; }
auto CScreeningData::get_density_inds_for_K_sd() const -> const std::vector<uint32_t>& { return _density_inds_for_K_sd; }
auto CScreeningData::get_density_inds_for_K_ds() const -> const std::vector<uint32_t>& { return _density_inds_for_K_ds; }
auto CScreeningData::get_density_inds_for_K_pp() const -> const std::vector<uint32_t>& { return _density_inds_for_K_pp; }
auto CScreeningData::get_density_inds_for_K_pd() const -> const std::vector<uint32_t>& { return _density_inds_for_K_pd; }
auto CScreeningData::get_density_inds_for_K_dp() const -> const std::vector<uint32_t>& { return _density_inds_for_K_dp; }
auto CScreeningData::get_density_inds_for_K_dd() const -> const std::vector<uint32_t>& { return _density_inds_for_K_dd; }

auto CScreeningData::get_local_pair_inds_i_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_ss[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_ss[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_sp[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_sp[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_sd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_sd[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_pp[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_pp[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_pd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_pd[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_dd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_dd[gpu_id]; }

auto CScreeningData::get_mat_Q_full(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) const -> CDenseMatrix
{
    const auto cart_naos = s_prim_count + p_prim_count * 3 + d_prim_count * 6;

    CDenseMatrix mat_Q_full(cart_naos, cart_naos);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            mat_Q_full.row(i)[j] = _Q_matrix_ss.row(i)[j];
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            auto j_full = s_prim_count + j;

            mat_Q_full.row(i)[j_full] = _Q_matrix_sp.row(i)[j];
            mat_Q_full.row(j_full)[i] = _Q_matrix_sp.row(i)[j];
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_Q_full.row(i)[j_full] = _Q_matrix_sd.row(i)[j];
            mat_Q_full.row(j_full)[i] = _Q_matrix_sd.row(i)[j];
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        auto i_full = s_prim_count + i;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            auto j_full = s_prim_count + j;

            mat_Q_full.row(i_full)[j_full] = _Q_matrix_pp.row(i)[j];
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_Q_full.row(i_full)[j_full] = _Q_matrix_pd.row(i)[j];
            mat_Q_full.row(j_full)[i_full] = _Q_matrix_pd.row(i)[j];
        }
    }

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        auto i_full = s_prim_count + p_prim_count * 3 + i;

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_Q_full.row(i_full)[j_full] = _Q_matrix_dd.row(i)[j];
        }
    }

    return mat_Q_full;
}

auto CScreeningData::get_mat_D_abs_full(const int64_t s_prim_count,
                                        const int64_t p_prim_count,
                                        const int64_t d_prim_count,
                                        const std::vector<uint32_t>& s_prim_aoinds,
                                        const std::vector<uint32_t>& p_prim_aoinds,
                                        const std::vector<uint32_t>& d_prim_aoinds,
                                        const int64_t naos,
                                        const double* dens_ptr) const -> CDenseMatrix
{
    const auto cart_naos = s_prim_count + p_prim_count * 3 + d_prim_count * 6;

    CDenseMatrix mat_D_abs_full(cart_naos, cart_naos);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            mat_D_abs_full.row(i)[j] = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + j;

            mat_D_abs_full.row(i)[j_full] = D_ij;
            mat_D_abs_full.row(j_full)[i] = D_ij;
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_D_abs_full.row(i)[j_full] = D_ij;
            mat_D_abs_full.row(j_full)[i] = D_ij;
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];

        auto i_full = s_prim_count + i;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            auto j_full = s_prim_count + j;

            mat_D_abs_full.row(i_full)[j_full] = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_D_abs_full.row(i_full)[j_full] = D_ij;
            mat_D_abs_full.row(j_full)[i_full] = D_ij;
        }
    }

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];

        auto i_full = s_prim_count + p_prim_count * 3 + i;

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            mat_D_abs_full.row(i_full)[j_full] = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);
        }
    }

    return mat_D_abs_full;
}

auto CScreeningData::form_mat_Q_and_density_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) -> void
{
    // Q_ss and D_ss for K

    _mat_Q_for_K_ss = std::vector<double> (s_prim_count * s_prim_count);
    _density_inds_for_K_ss = std::vector<uint32_t> (s_prim_count * s_prim_count);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_ss.row(i);

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_ss[i * s_prim_count + j]        = q_val;
            _density_inds_for_K_ss[i * s_prim_count + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_sp and D_sp for K

    _mat_Q_for_K_sp = std::vector<double> (s_prim_count * p_prim_count * 3);
    _density_inds_for_K_sp = std::vector<uint32_t> (s_prim_count * p_prim_count * 3);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_sp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_sp[i * p_prim_count * 3 + j]        = q_val;
            _density_inds_for_K_sp[i * p_prim_count * 3 + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_ps and D_ps for K

    _mat_Q_for_K_ps = std::vector<double> (p_prim_count * 3 * s_prim_count);
    _density_inds_for_K_ps = std::vector<uint32_t> (p_prim_count * 3 * s_prim_count);

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_sp.row(j)[i], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_ps[i * s_prim_count + j]        = q_val;
            _density_inds_for_K_ps[i * s_prim_count + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_sd and D_sd for K

    _mat_Q_for_K_sd = std::vector<double> (s_prim_count * d_prim_count * 6);
    _density_inds_for_K_sd = std::vector<uint32_t> (s_prim_count * d_prim_count * 6);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_sd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_sd[i * d_prim_count * 6 + j]        = q_val;
            _density_inds_for_K_sd[i * d_prim_count * 6 + j] = static_cast<uint32_t>(j_idx);
        }
    
    }

    // Q_ds and D_ds for K

    _mat_Q_for_K_ds = std::vector<double> (d_prim_count * 6 * s_prim_count);
    _density_inds_for_K_ds = std::vector<uint32_t> (d_prim_count * 6 * s_prim_count);

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_sd.row(j)[i], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_ds[i * s_prim_count + j]        = q_val;
            _density_inds_for_K_ds[i * s_prim_count + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_pp and D_pp for K

    _mat_Q_for_K_pp = std::vector<double> (p_prim_count * 3 * p_prim_count * 3);
    _density_inds_for_K_pp = std::vector<uint32_t> (p_prim_count * 3 * p_prim_count * 3);

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_pp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_pp[i * p_prim_count * 3 + j]        = q_val;
            _density_inds_for_K_pp[i * p_prim_count * 3 + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_pd and D_pd for K

    _mat_Q_for_K_pd = std::vector<double> (p_prim_count * 3 * d_prim_count * 6);
    _density_inds_for_K_pd = std::vector<uint32_t> (p_prim_count * 3 * d_prim_count * 6);

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_pd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_pd[i * d_prim_count * 6 + j]        = q_val;
            _density_inds_for_K_pd[i * d_prim_count * 6 + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_dp and D_dp for K

    _mat_Q_for_K_dp = std::vector<double> (d_prim_count * 6 * p_prim_count * 3);
    _density_inds_for_K_dp = std::vector<uint32_t> (d_prim_count * 6 * p_prim_count * 3);

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_pd.row(j)[i], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_dp[i * p_prim_count * 3 + j]        = q_val;
            _density_inds_for_K_dp[i * p_prim_count * 3 + j] = static_cast<uint32_t>(j_idx);
        }
    }

    // Q_dd and D_dd for K

    _mat_Q_for_K_dd = std::vector<double> (d_prim_count * 6 * d_prim_count * 6);
    _density_inds_for_K_dd = std::vector<uint32_t> (d_prim_count * 6 * d_prim_count * 6);

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_dd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            Q_vec_sorted.push_back(std::make_tuple(Q_i[j], i, j));
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            // auto i_idx = std::get<1>(q_ij);
            auto j_idx = std::get<2>(q_ij);

            _mat_Q_for_K_dd[i * d_prim_count * 6 + j]        = q_val;
            _density_inds_for_K_dd[i * d_prim_count * 6 + j] = static_cast<uint32_t>(j_idx);
        }
    }
}

auto CScreeningData::form_pair_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count, const CDenseMatrix& Q_prime, const double Q_prime_thresh) -> void
{
    // TODO think about sorting by Q_prime bound

    // ss, sp, sd blocks

    std::vector<uint32_t> pair_inds_i_for_K_ss;
    std::vector<uint32_t> pair_inds_k_for_K_ss;

    std::vector<uint32_t> pair_inds_i_for_K_sp;
    std::vector<uint32_t> pair_inds_k_for_K_sp;

    std::vector<uint32_t> pair_inds_i_for_K_sd;
    std::vector<uint32_t> pair_inds_k_for_K_sd;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = i; k < s_prim_count; k++)
        {
            if (std::fabs(Q_prime.row(i)[k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_ss.push_back(i);
                pair_inds_k_for_K_ss.push_back(k);
            }
        }

        for (int64_t k = 0; k < p_prim_count * 3; k++)
        {
            if (std::fabs(Q_prime.row(i)[s_prim_count + k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_sp.push_back(i);
                pair_inds_k_for_K_sp.push_back(k);
            }
        }

        for (int64_t k = 0; k < d_prim_count * 6; k++)
        {
            if (std::fabs(Q_prime.row(i)[s_prim_count + p_prim_count * 3 + k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_sd.push_back(i);
                pair_inds_k_for_K_sd.push_back(k);
            }
        }
    }

    // pp, pd blocks

    std::vector<uint32_t> pair_inds_i_for_K_pp;
    std::vector<uint32_t> pair_inds_k_for_K_pp;

    std::vector<uint32_t> pair_inds_i_for_K_pd;
    std::vector<uint32_t> pair_inds_k_for_K_pd;

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t k = i; k < p_prim_count * 3; k++)
        {
            if (std::fabs(Q_prime.row(s_prim_count + i)[s_prim_count + k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_pp.push_back(i);
                pair_inds_k_for_K_pp.push_back(k);
            }
        }

        for (int64_t k = 0; k < d_prim_count * 6; k++)
        {
            if (std::fabs(Q_prime.row(s_prim_count + i)[s_prim_count + p_prim_count * 3 + k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_pd.push_back(i);
                pair_inds_k_for_K_pd.push_back(k);
            }
        }
    }

    // dd block

    std::vector<uint32_t> pair_inds_i_for_K_dd;
    std::vector<uint32_t> pair_inds_k_for_K_dd;

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        for (int64_t k = i; k < d_prim_count * 6; k++)
        {
            if (std::fabs(Q_prime.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + k]) > Q_prime_thresh)
            {
                pair_inds_i_for_K_dd.push_back(i);
                pair_inds_k_for_K_dd.push_back(k);
            }
        }
    }

    const auto ss_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_ss.size());
    const auto sp_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_sp.size());
    const auto sd_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_sd.size());
    const auto pp_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_pp.size());
    const auto pd_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_pd.size());
    const auto dd_pair_count_for_K = static_cast<int64_t>(pair_inds_i_for_K_dd.size());

    // form local vectors

    _local_pair_inds_i_for_K_ss = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_ss = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_sp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_sp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_sd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_sd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_pp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_pp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_pd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_pd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_dd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_dd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    // TODO use communicator from arguments
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
        auto gpu_count = nnodes * _num_gpus_per_node;

        auto ss_batch_size = mathfunc::batch_size(ss_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_ss[gpu_id] = std::vector<uint32_t>(ss_batch_size);
        _local_pair_inds_k_for_K_ss[gpu_id] = std::vector<uint32_t>(ss_batch_size);

        auto sp_batch_size = mathfunc::batch_size(sp_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_sp[gpu_id] = std::vector<uint32_t>(sp_batch_size);
        _local_pair_inds_k_for_K_sp[gpu_id] = std::vector<uint32_t>(sp_batch_size);

        auto sd_batch_size = mathfunc::batch_size(sd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_sd[gpu_id] = std::vector<uint32_t>(sd_batch_size);
        _local_pair_inds_k_for_K_sd[gpu_id] = std::vector<uint32_t>(sd_batch_size);

        auto pp_batch_size = mathfunc::batch_size(pp_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_pp[gpu_id] = std::vector<uint32_t>(pp_batch_size);
        _local_pair_inds_k_for_K_pp[gpu_id] = std::vector<uint32_t>(pp_batch_size);

        auto pd_batch_size = mathfunc::batch_size(pd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_pd[gpu_id] = std::vector<uint32_t>(pd_batch_size);
        _local_pair_inds_k_for_K_pd[gpu_id] = std::vector<uint32_t>(pd_batch_size);

        auto dd_batch_size = mathfunc::batch_size(dd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_dd[gpu_id] = std::vector<uint32_t>(dd_batch_size);
        _local_pair_inds_k_for_K_dd[gpu_id] = std::vector<uint32_t>(dd_batch_size);
    }

    auto nthreads = omp_get_max_threads();
    auto num_threads_per_gpu = nthreads / _num_gpus_per_node;

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        if (thread_id % num_threads_per_gpu == 0)
        {
            auto gpu_id = thread_id / num_threads_per_gpu;
            auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
            auto gpu_count = nnodes * _num_gpus_per_node;

            for (int64_t ik = gpu_rank, idx = 0; ik < ss_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_ss[gpu_id][idx] = pair_inds_i_for_K_ss[ik];
                _local_pair_inds_k_for_K_ss[gpu_id][idx] = pair_inds_k_for_K_ss[ik];
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < sp_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_sp[gpu_id][idx] = pair_inds_i_for_K_sp[ik];
                _local_pair_inds_k_for_K_sp[gpu_id][idx] = pair_inds_k_for_K_sp[ik];
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < sd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_sd[gpu_id][idx] = pair_inds_i_for_K_sd[ik];
                _local_pair_inds_k_for_K_sd[gpu_id][idx] = pair_inds_k_for_K_sd[ik];
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < pp_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_pp[gpu_id][idx] = pair_inds_i_for_K_pp[ik];
                _local_pair_inds_k_for_K_pp[gpu_id][idx] = pair_inds_k_for_K_pp[ik];
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < pd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_pd[gpu_id][idx] = pair_inds_i_for_K_pd[ik];
                _local_pair_inds_k_for_K_pd[gpu_id][idx] = pair_inds_k_for_K_pd[ik];
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < dd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                _local_pair_inds_i_for_K_dd[gpu_id][idx] = pair_inds_i_for_K_dd[ik];
                _local_pair_inds_k_for_K_dd[gpu_id][idx] = pair_inds_k_for_K_dd[ik];
            }
        }
    }
}

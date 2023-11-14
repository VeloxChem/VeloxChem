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

#include "EriScreener.hpp"

#include <cmath>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"

#define MATH_CONST_PI 3.14159265358979323846

namespace eriscreen {  // eriscreen namespace

auto
computeQMatrices(const CMolecule& molecule, const CMolecularBasis& basis) -> std::vector<CDenseMatrix>
{
    std::vector<CDenseMatrix> Q_matrices(3);

    // Boys function

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();

    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    // GTO blocks

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
    }

    // S block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    // P block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    // GTO block pairs

    Q_matrices[0] = CDenseMatrix(s_prim_count, s_prim_count);
    Q_matrices[1] = CDenseMatrix(s_prim_count, p_prim_count * 3);
    Q_matrices[2] = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);

    // TODO distribute computation of Q matrices

    // S-S and S-P block pairs

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++, ij++)
        {
            const auto a_j = s_prim_info[j + s_prim_count * 0];
            const auto c_j = s_prim_info[j + s_prim_count * 1];
            const auto x_j = s_prim_info[j + s_prim_count * 2];
            const auto y_j = s_prim_info[j + s_prim_count * 3];
            const auto z_j = s_prim_info[j + s_prim_count * 4];

            // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto F0_t = boysfunc::getBoysFunction(0.0, 0, boys_func_table.data(), boys_func_ft.data());

            const auto eri_ijij = F0_t[0] * S_ij * S_ij * std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            const auto sqrt_ijij = std::sqrt(eri_ijij);

            Q_matrices[0].row(i)[j] = sqrt_ijij;

            if (i != j) Q_matrices[0].row(j)[i] = sqrt_ijij;
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            const auto a_j = p_prim_info[j + p_prim_count * 0];
            const auto c_j = p_prim_info[j + p_prim_count * 1];
            const auto x_j = p_prim_info[j + p_prim_count * 2];
            const auto y_j = p_prim_info[j + p_prim_count * 3];
            const auto z_j = p_prim_info[j + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto a_ij = a_i + a_j;

            const auto S_ij_0 = c_i * c_j * std::pow(MATH_CONST_PI / a_ij, 1.5) * std::exp(-a_i * a_j / a_ij * r2_ij);

            const auto Lambda = std::sqrt(4.0 * a_ij * a_ij / (a_ij + a_ij) / MATH_CONST_PI);

            for (int64_t s = 0; s < 3; s++)
            {
                // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)

                // p-1: py
                // p_0: pz
                // p+1: px

                double S_ij_1;

                if (s == 0) S_ij_1 = -(a_i / a_ij) * (y_j - y_i) * S_ij_0;
                if (s == 1) S_ij_1 = -(a_i / a_ij) * (z_j - z_i) * S_ij_0;
                if (s == 2) S_ij_1 = -(a_i / a_ij) * (x_j - x_i) * S_ij_0;

                const auto F1_t = boysfunc::getBoysFunction(0.0, 1, boys_func_table.data(), boys_func_ft.data());

                const double eri_ijij = Lambda * (S_ij_1 * S_ij_1 * F1_t[0] + S_ij_0 * S_ij_0 * F1_t[1] / (2.0 * (a_ij + a_ij)));

                const auto sqrt_ijij = std::sqrt(eri_ijij);

                // TODO: think about the ordering of cartesian components

                Q_matrices[1].row(i)[j * 3 + s] = sqrt_ijij;
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

            const auto Lambda = std::sqrt(4.0 * (a_i + a_j) * (a_i + a_j) / (MATH_CONST_PI * (a_i + a_j + a_i + a_j)));

            const auto F2_t = boysfunc::getBoysFunction(0.0, 2, boys_func_table.data(), boys_func_ft.data());

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            // p-1: py
            // p_0: pz
            // p+1: px

            const double rij[3] = {y_j - y_i, z_j - z_i, x_j - x_i};

            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto S_ij_10 = (a_j / (a_i + a_j)) * rij[i_cart] * S_ij_00;
                    const auto S_ij_01 = -(a_i / (a_i + a_j)) * rij[j_cart] * S_ij_00;

                    const auto S_ij_11 = (

                            (i_cart == j_cart ? 1.0 : 0.0) * (0.5 / (a_i + a_j)) +

                            (a_i * a_j / ((a_i + a_j) * (a_i + a_j))) * (-rij[i_cart]) * rij[j_cart]

                            ) * S_ij_00;

                    const double eri_ijij = Lambda * (

                            S_ij_11 * S_ij_11 * F2_t[0]

                            + S_ij_11 * S_ij_00 * (a_i + a_j) / (a_i + a_j + a_i + a_j) * ( -(i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F2_t[1] )

                            + S_ij_10 * S_ij_10 / (a_i + a_j + a_i + a_j) * ( 0.5 * F2_t[1] )

                            + S_ij_10 * S_ij_01 / (a_i + a_j + a_i + a_j) * ( (j_cart == i_cart ? 1.0 : 0.0) * 0.5 * F2_t[1] )

                            + S_ij_01 * S_ij_10 / (a_i + a_j + a_i + a_j) * ( (i_cart == j_cart ? 1.0 : 0.0) * 0.5 * F2_t[1] )

                            + S_ij_01 * S_ij_01 / (a_i + a_j + a_i + a_j) * ( 0.5 * F2_t[1] )

                            + S_ij_00 * S_ij_11 * (a_i + a_j) / (a_i + a_j + a_i + a_j) * ( -(i_cart == j_cart ? 1.0 : 0.0) * 0.5 / (a_i + a_j) * F2_t[1] )

                            + S_ij_00 * S_ij_00 / ((a_i + a_j + a_i + a_j) * (a_i + a_j + a_i + a_j)) * (

                                  ((i_cart == j_cart ? 1.0 : 0.0) * (i_cart == j_cart ? 1.0 : 0.0) + 1.0 +
                                   (i_cart == j_cart ? 1.0 : 0.0) * (j_cart == i_cart ? 1.0 : 0.0)) * 0.25 * F2_t[2]

                            ));

                    const auto sqrt_ijij = std::sqrt(eri_ijij);

                    // TODO: think about the ordering of cartesian components

                    Q_matrices[2].row(i * 3 + i_cart)[j * 3 + j_cart] = sqrt_ijij;

                    if (i * 3 + i_cart != j * 3 + j_cart) Q_matrices[2].row(j * 3 + j_cart)[i * 3 + i_cart] = sqrt_ijij;
                }
            }
        }
    }

    return Q_matrices;
}

auto
computeDMatrices(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> std::vector<CDenseMatrix>
{
    std::vector<CDenseMatrix> D_matrices(3);

    // GTO blocks

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::string errnaos("eriscreen::computeDMatrices: Inconsistent number of AOs");

    errors::assertMsgCritical((naos == densityMatrix.getNumberOfRows(0)) && (naos == densityMatrix.getNumberOfColumns(0)), errnaos);

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
    }

    // S block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    // P block

    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    // GTO block pairs

    D_matrices[0] = CDenseMatrix(s_prim_count, s_prim_count);
    D_matrices[1] = CDenseMatrix(s_prim_count, p_prim_count * 3);
    D_matrices[2] = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);

    auto dens_ptr = densityMatrix.alphaDensity(0);

    // S-S and S-P gto block pairs

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        const auto i_cgto = s_prim_aoinds[i];

        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++, ij++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

            D_matrices[0].row(i)[j] = D_ij;

            if (i != j) D_matrices[0].row(j)[i] = D_ij;
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * s];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                // TODO: think about the ordering of cartesian components

                D_matrices[1].row(i)[j * 3 + s] = D_ij;
            }
        }
    }

    // P-P gto block pair

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                    // TODO: think about the ordering of cartesian components

                    D_matrices[2].row(i * 3 + i_cart)[j * 3 + j_cart] = D_ij;

                    if (i * 3 + i_cart != j * 3 + j_cart) D_matrices[2].row(j * 3 + j_cart)[i * 3 + i_cart] = D_ij;
                }
            }
        }
    }

    return D_matrices;
}

}  // namespace eriscreen

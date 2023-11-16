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

#include <cmath>
#include <string>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

#define MATH_CONST_PI 3.14159265358979323846

CScreeningData::CScreeningData(const CMolecule& molecule, const CMolecularBasis& basis)
{
    _computeQMatrices(molecule, basis);

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

    _sortQ(s_prim_count, p_prim_count);
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

    _Q_matrix_ss = CDenseMatrix(s_prim_count, s_prim_count);
    _Q_matrix_sp = CDenseMatrix(s_prim_count, p_prim_count * 3);
    _Q_matrix_pp = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);

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

            _Q_matrix_ss.row(i)[j] = sqrt_ijij;

            if (i != j) _Q_matrix_ss.row(j)[i] = sqrt_ijij;
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

                _Q_matrix_sp.row(i)[j * 3 + s] = sqrt_ijij;
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

                    _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart] = sqrt_ijij;

                    if (i * 3 + i_cart != j * 3 + j_cart) _Q_matrix_pp.row(j * 3 + j_cart)[i * 3 + i_cart] = sqrt_ijij;
                }
            }
        }
    }
}

auto
CScreeningData::updateDensityMatrix(const CMolecule& molecule, const CMolecularBasis& basis, const CAODensityMatrix& densityMatrix) -> void
{
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

    _D_matrix_ss = CDenseMatrix(s_prim_count, s_prim_count);
    _D_matrix_sp = CDenseMatrix(s_prim_count, p_prim_count * 3);
    _D_matrix_pp = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);

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

            _D_matrix_ss.row(i)[j] = D_ij;

            if (i != j) _D_matrix_ss.row(j)[i] = D_ij;
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * s];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                // TODO: think about the ordering of cartesian components

                _D_matrix_sp.row(i)[j * 3 + s] = D_ij;
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

                    _D_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart] = D_ij;

                    if (i * 3 + i_cart != j * 3 + j_cart) _D_matrix_pp.row(j * 3 + j_cart)[i * 3 + i_cart] = D_ij;
                }
            }
        }
    }
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
CScreeningData::getQMatrixPP() const -> const CDenseMatrix&
{
    return _Q_matrix_pp;
}

auto
CScreeningData::getDMatrixSS() const -> const CDenseMatrix&
{
    return _D_matrix_ss;
}

auto
CScreeningData::getDMatrixSP() const -> const CDenseMatrix&
{
    return _D_matrix_sp;
}

auto
CScreeningData::getDMatrixPP() const -> const CDenseMatrix&
{
    return _D_matrix_pp;
}

auto
CScreeningData::_sortQ(const int64_t s_prim_count, const int64_t p_prim_count) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    std::vector<std::tuple<double, int64_t, int64_t>> sorted_ss_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_sp_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_pp_mat_Q;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++, ij++)
        {
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            sorted_ss_mat_Q.push_back(std::make_tuple(Q_ij, i, j));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + s];
                sorted_sp_mat_Q.push_back(std::make_tuple(Q_ij, i, j * 3 + s));
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
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    sorted_pp_mat_Q.push_back(std::make_tuple(Q_ij, i * 3 + i_cart, j * 3 + j_cart));
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q.begin(), sorted_ss_mat_Q.end());
    std::sort(sorted_sp_mat_Q.begin(), sorted_sp_mat_Q.end());
    std::sort(sorted_pp_mat_Q.begin(), sorted_pp_mat_Q.end());

    const auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    const auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    const auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    _ss_first_inds_local  = std::vector<uint32_t>();
    _ss_second_inds_local = std::vector<uint32_t>();
    _ss_mat_Q_local       = std::vector<double>();

    for (int64_t ij = rank; ij < ss_prim_pair_count; ij+=nnodes)
    {
        const auto& vals = sorted_ss_mat_Q[ss_prim_pair_count - 1 - ij];

        auto Q_ij = std::get<0>(vals);
        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);

        _ss_first_inds_local.push_back(static_cast<uint32_t>(i));
        _ss_second_inds_local.push_back(static_cast<uint32_t>(j));
        _ss_mat_Q_local.push_back(Q_ij);
    }

    _sp_first_inds_local  = std::vector<uint32_t>();
    _sp_second_inds_local = std::vector<uint32_t>();
    _sp_mat_Q_local       = std::vector<double>();

    for (int64_t ij = rank; ij < sp_prim_pair_count; ij+=nnodes)
    {
        const auto& vals = sorted_sp_mat_Q[sp_prim_pair_count - 1 - ij];

        auto Q_ij = std::get<0>(vals);
        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);

        _sp_first_inds_local.push_back(static_cast<uint32_t>(i));
        _sp_second_inds_local.push_back(static_cast<uint32_t>(j));
        _sp_mat_Q_local.push_back(Q_ij);
    }

    _pp_first_inds_local  = std::vector<uint32_t>();
    _pp_second_inds_local = std::vector<uint32_t>();
    _pp_mat_Q_local       = std::vector<double>();

    for (int64_t ij = rank; ij < pp_prim_pair_count; ij+=nnodes)
    {
        const auto& vals = sorted_pp_mat_Q[pp_prim_pair_count - 1 - ij];

        auto Q_ij = std::get<0>(vals);
        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);

        _pp_first_inds_local.push_back(static_cast<uint32_t>(i));
        _pp_second_inds_local.push_back(static_cast<uint32_t>(j));
        _pp_mat_Q_local.push_back(Q_ij);
    }
}

auto
CScreeningData::sortQD(const int64_t s_prim_count, const int64_t p_prim_count) -> void
{
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_ss_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_sp_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_pp_mat_Q_D;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0, ij = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++, ij++)
        {
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            const auto D_ij = _D_matrix_ss.row(i)[j];
            sorted_ss_mat_Q_D.push_back(std::make_tuple(Q_ij * std::fabs(D_ij), i, j, Q_ij, D_ij));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++, ij++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + s];
                const auto D_ij = _D_matrix_sp.row(i)[j * 3 + s];
                sorted_sp_mat_Q_D.push_back(std::make_tuple(Q_ij * std::fabs(D_ij), i, j * 3 + s, Q_ij, D_ij));
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
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    const auto D_ij = _D_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    sorted_pp_mat_Q_D.push_back(std::make_tuple(Q_ij * std::fabs(D_ij), i * 3 + i_cart, j * 3 + j_cart, Q_ij, D_ij));
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q_D.begin(), sorted_ss_mat_Q_D.end());
    std::sort(sorted_sp_mat_Q_D.begin(), sorted_sp_mat_Q_D.end());
    std::sort(sorted_pp_mat_Q_D.begin(), sorted_pp_mat_Q_D.end());

    const auto ss_prim_pair_count = s_prim_count * (s_prim_count + 1) / 2;
    const auto sp_prim_pair_count = s_prim_count * p_prim_count * 3;
    const auto pp_prim_pair_count = p_prim_count * 3 * (p_prim_count * 3 + 1) / 2;

    _ss_max_D = 0.0;
    _sp_max_D = 0.0;
    _pp_max_D = 0.0;

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
}

auto CScreeningData::get_ss_first_inds_local() const -> const std::vector<uint32_t>& { return _ss_first_inds_local; }
auto CScreeningData::get_sp_first_inds_local() const -> const std::vector<uint32_t>& { return _sp_first_inds_local; }
auto CScreeningData::get_pp_first_inds_local() const -> const std::vector<uint32_t>& { return _pp_first_inds_local; }

auto CScreeningData::get_ss_second_inds_local() const -> const std::vector<uint32_t>& { return _ss_second_inds_local; }
auto CScreeningData::get_sp_second_inds_local() const -> const std::vector<uint32_t>& { return _sp_second_inds_local; }
auto CScreeningData::get_pp_second_inds_local() const -> const std::vector<uint32_t>& { return _pp_second_inds_local; }

auto CScreeningData::get_ss_mat_Q_local() const -> const std::vector<double>& { return _ss_mat_Q_local; }
auto CScreeningData::get_sp_mat_Q_local() const -> const std::vector<double>& { return _sp_mat_Q_local; }
auto CScreeningData::get_pp_mat_Q_local() const -> const std::vector<double>& { return _pp_mat_Q_local; }

auto CScreeningData::get_ss_first_inds() const -> const std::vector<uint32_t>& { return _ss_first_inds; }
auto CScreeningData::get_sp_first_inds() const -> const std::vector<uint32_t>& { return _sp_first_inds; }
auto CScreeningData::get_pp_first_inds() const -> const std::vector<uint32_t>& { return _pp_first_inds; }

auto CScreeningData::get_ss_second_inds() const -> const std::vector<uint32_t>& { return _ss_second_inds; }
auto CScreeningData::get_sp_second_inds() const -> const std::vector<uint32_t>& { return _sp_second_inds; }
auto CScreeningData::get_pp_second_inds() const -> const std::vector<uint32_t>& { return _pp_second_inds; }

auto CScreeningData::get_ss_mat_Q() const -> const std::vector<double>& { return _ss_mat_Q; }
auto CScreeningData::get_sp_mat_Q() const -> const std::vector<double>& { return _sp_mat_Q; }
auto CScreeningData::get_pp_mat_Q() const -> const std::vector<double>& { return _pp_mat_Q; }

auto CScreeningData::get_ss_mat_D() const -> const std::vector<double>& { return _ss_mat_D; }
auto CScreeningData::get_sp_mat_D() const -> const std::vector<double>& { return _sp_mat_D; }
auto CScreeningData::get_pp_mat_D() const -> const std::vector<double>& { return _pp_mat_D; }

auto CScreeningData::get_ss_max_D() const -> double { return _ss_max_D; }
auto CScreeningData::get_sp_max_D() const -> double { return _sp_max_D; }
auto CScreeningData::get_pp_max_D() const -> double { return _pp_max_D; }

auto CScreeningData::get_mat_Q_for_K_ss() const -> const std::vector<double>& { return _mat_Q_for_K_ss; };
auto CScreeningData::get_mat_Q_for_K_sp() const -> const std::vector<double>& { return _mat_Q_for_K_sp; };
auto CScreeningData::get_mat_Q_for_K_ps() const -> const std::vector<double>& { return _mat_Q_for_K_ps; };
auto CScreeningData::get_mat_Q_for_K_pp() const -> const std::vector<double>& { return _mat_Q_for_K_pp; };

auto CScreeningData::get_density_inds_for_K_ss() const -> const std::vector<uint32_t>& { return _density_inds_for_K_ss; }
auto CScreeningData::get_density_inds_for_K_sp() const -> const std::vector<uint32_t>& { return _density_inds_for_K_sp; }
auto CScreeningData::get_density_inds_for_K_ps() const -> const std::vector<uint32_t>& { return _density_inds_for_K_ps; }
auto CScreeningData::get_density_inds_for_K_pp() const -> const std::vector<uint32_t>& { return _density_inds_for_K_pp; }

auto CScreeningData::get_pair_inds_i_for_K_ss() const -> const std::vector<uint32_t>& { return _pair_inds_i_for_K_ss; };
auto CScreeningData::get_pair_inds_k_for_K_ss() const -> const std::vector<uint32_t>& { return _pair_inds_k_for_K_ss; };

auto CScreeningData::get_pair_inds_i_for_K_sp() const -> const std::vector<uint32_t>& { return _pair_inds_i_for_K_sp; };
auto CScreeningData::get_pair_inds_k_for_K_sp() const -> const std::vector<uint32_t>& { return _pair_inds_k_for_K_sp; };

auto CScreeningData::get_pair_inds_i_for_K_pp() const -> const std::vector<uint32_t>& { return _pair_inds_i_for_K_pp; };
auto CScreeningData::get_pair_inds_k_for_K_pp() const -> const std::vector<uint32_t>& { return _pair_inds_k_for_K_pp; };

auto CScreeningData::get_mat_Q_full(const int64_t s_prim_count, const int64_t p_prim_count) const -> CDenseMatrix
{
    CDenseMatrix mat_Q_full(s_prim_count + p_prim_count * 3, s_prim_count + p_prim_count * 3);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            mat_Q_full.row(i)[j] = _Q_matrix_ss.row(i)[j];
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_Q_full.row(i)[s_prim_count + j] = _Q_matrix_sp.row(i)[j];
            mat_Q_full.row(s_prim_count + j)[i] = _Q_matrix_sp.row(i)[j];
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_Q_full.row(s_prim_count + i)[s_prim_count + j] = _Q_matrix_pp.row(i)[j];
        }
    }

    return mat_Q_full;
}

auto CScreeningData::get_mat_D_abs_full(const int64_t s_prim_count, const int64_t p_prim_count) const -> CDenseMatrix
{
    CDenseMatrix mat_D_abs_full(s_prim_count + p_prim_count * 3, s_prim_count + p_prim_count * 3);

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            mat_D_abs_full.row(i)[j] = std::fabs(_D_matrix_ss.row(i)[j]);
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_D_abs_full.row(i)[s_prim_count + j] = std::fabs(_D_matrix_sp.row(i)[j]);
            mat_D_abs_full.row(s_prim_count + j)[i] = std::fabs(_D_matrix_sp.row(i)[j]);
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            mat_D_abs_full.row(s_prim_count + i)[s_prim_count + j] = std::fabs(_D_matrix_pp.row(i)[j]);
        }
    }

    return mat_D_abs_full;
}

auto CScreeningData::form_mat_Q_and_density_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count) -> void
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
}

auto CScreeningData::form_pair_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count, const CDenseMatrix& Q_prime, const double Q_prime_thresh) -> void
{
    auto rank = mpi::rank(MPI_COMM_WORLD);
    auto nnodes = mpi::nodes(MPI_COMM_WORLD);

    // ss block

    _pair_inds_i_for_K_ss = std::vector<uint32_t>();
    _pair_inds_k_for_K_ss = std::vector<uint32_t>();

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = i; k < s_prim_count; k++)
        {
            if (std::fabs(Q_prime.row(i)[k]) > Q_prime_thresh)
            {
                _pair_inds_i_for_K_ss.push_back(i);
                _pair_inds_k_for_K_ss.push_back(k);
            }
        }
    }

    std::vector<uint32_t> local_pair_inds_i_for_K_ss;
    std::vector<uint32_t> local_pair_inds_k_for_K_ss;

    for (int64_t ik = rank; ik < static_cast<int64_t>(_pair_inds_i_for_K_ss.size()); ik+=nnodes)
    {
        local_pair_inds_i_for_K_ss.push_back(_pair_inds_i_for_K_ss[ik]);
        local_pair_inds_k_for_K_ss.push_back(_pair_inds_k_for_K_ss[ik]);
    }

    _pair_inds_i_for_K_ss = local_pair_inds_i_for_K_ss;
    _pair_inds_k_for_K_ss = local_pair_inds_k_for_K_ss;

    // sp block

    _pair_inds_i_for_K_sp = std::vector<uint32_t>();
    _pair_inds_k_for_K_sp = std::vector<uint32_t>();

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = 0; k < p_prim_count * 3; k++)
        {
            if (std::fabs(Q_prime.row(i)[s_prim_count + k]) > Q_prime_thresh)
            {
                _pair_inds_i_for_K_sp.push_back(i);
                _pair_inds_k_for_K_sp.push_back(k);
            }
        }
    }

    std::vector<uint32_t> local_pair_inds_i_for_K_sp;
    std::vector<uint32_t> local_pair_inds_k_for_K_sp;

    for (int64_t ik = rank; ik < static_cast<int64_t>(_pair_inds_i_for_K_sp.size()); ik+=nnodes)
    {
        local_pair_inds_i_for_K_sp.push_back(_pair_inds_i_for_K_sp[ik]);
        local_pair_inds_k_for_K_sp.push_back(_pair_inds_k_for_K_sp[ik]);
    }

    _pair_inds_i_for_K_sp = local_pair_inds_i_for_K_sp;
    _pair_inds_k_for_K_sp = local_pair_inds_k_for_K_sp;

    // pp block

    _pair_inds_i_for_K_pp = std::vector<uint32_t>();
    _pair_inds_k_for_K_pp = std::vector<uint32_t>();

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t k = i; k < p_prim_count * 3; k++)
        {
            if (std::fabs(Q_prime.row(s_prim_count + i)[s_prim_count + k]) > Q_prime_thresh)
            {
                _pair_inds_i_for_K_pp.push_back(i);
                _pair_inds_k_for_K_pp.push_back(k);
            }
        }
    }

    std::vector<uint32_t> local_pair_inds_i_for_K_pp;
    std::vector<uint32_t> local_pair_inds_k_for_K_pp;

    for (int64_t ik = rank; ik < static_cast<int64_t>(_pair_inds_i_for_K_pp.size()); ik+=nnodes)
    {
        local_pair_inds_i_for_K_pp.push_back(_pair_inds_i_for_K_pp[ik]);
        local_pair_inds_k_for_K_pp.push_back(_pair_inds_k_for_K_pp[ik]);
    }

    _pair_inds_i_for_K_pp = local_pair_inds_i_for_K_pp;
    _pair_inds_k_for_K_pp = local_pair_inds_k_for_K_pp;
}

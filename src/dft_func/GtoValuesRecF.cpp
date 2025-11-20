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

#include "GtoValuesRecF.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "SphericalMomentum.hpp"

#define ANGULAR_MOMENTUM_F 3

namespace gtoval {  // gtoval namespace

auto
get_lda_values_rec_f(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>& gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_F * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_F>(isph);
    }

    // set up GTO values storage

    const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1));

    const size_t ncols = grid_coords_x.size();

    auto gto_values = matfunc::make_matrix("LDA", 7 * nrows, ncols);

    auto submat = gto_values.sub_matrix({0, 0});

    submat->zero();

    // set up GTOs data

    const auto gto_exps = gto_block.exponents();

    const auto gto_norms = gto_block.normalization_factors();

    const auto gto_coords = gto_block.coordinates();

    // set up grid data

    auto g_x = grid_coords_x.data();

    auto g_y = grid_coords_y.data();

    auto g_z = grid_coords_z.data();

    // set GTOs block dimensions

    const auto ncgtos = gto_block.number_of_basis_functions();

    const auto npgtos = gto_block.number_of_primitives();

    // compute GTO values for F type GTOs on grid

    std::vector<double> buffer_xxx(ncols);

    std::vector<double> buffer_xxy(ncols);

    std::vector<double> buffer_xxz(ncols);

    std::vector<double> buffer_xyy(ncols);

    std::vector<double> buffer_xyz(ncols);

    std::vector<double> buffer_xzz(ncols);

    std::vector<double> buffer_yyy(ncols);

    std::vector<double> buffer_yyz(ncols);

    std::vector<double> buffer_yzz(ncols);

    std::vector<double> buffer_zzz(ncols);

    std::vector<std::vector<double>*> buffer_refs({
        &buffer_xxx, &buffer_xxy, &buffer_xxz, &buffer_xyy, &buffer_xyz, &buffer_xzz,
        &buffer_yyy, &buffer_yyz, &buffer_yzz, &buffer_zzz});

    auto ptr_buffer_xxx = buffer_xxx.data();

    auto ptr_buffer_xxy = buffer_xxy.data();

    auto ptr_buffer_xxz = buffer_xxz.data();

    auto ptr_buffer_xyy = buffer_xyy.data();

    auto ptr_buffer_xyz = buffer_xyz.data();

    auto ptr_buffer_xzz = buffer_xzz.data();

    auto ptr_buffer_yyy = buffer_yyy.data();

    auto ptr_buffer_yyz = buffer_yyz.data();

    auto ptr_buffer_yzz = buffer_yzz.data();

    auto ptr_buffer_zzz = buffer_zzz.data();

    int irow = 0;

    for (int i = 0; i < ncgtos; i++)
    {
        if (gtos_mask[i] == 1)
        {
            // set up GTO coordinates

            const auto xyz = gto_coords[i].coordinates();
        
            const auto r_x = xyz[0];

            const auto r_y = xyz[1];

            const auto r_z = xyz[2];

            // compute GTO values on grid

            std::ranges::fill(buffer_xxx, 0.0);

            std::ranges::fill(buffer_xxy, 0.0);

            std::ranges::fill(buffer_xxz, 0.0);

            std::ranges::fill(buffer_xyy, 0.0);

            std::ranges::fill(buffer_xyz, 0.0);

            std::ranges::fill(buffer_xzz, 0.0);

            std::ranges::fill(buffer_yyy, 0.0);

            std::ranges::fill(buffer_yyz, 0.0);

            std::ranges::fill(buffer_yzz, 0.0);

            std::ranges::fill(buffer_zzz, 0.0);

            for (int j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int k = 0; k < ncols; k++)
                {
                    const auto gr_x = g_x[k] - r_x;

                    const auto gr_y = g_y[k] - r_y;

                    const auto gr_z = g_z[k] - r_z;

                    const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    ptr_buffer_xxx[k] += gr_x * gr_x * gr_x * fss;

                    ptr_buffer_xxy[k] += gr_x * gr_x * gr_y * fss;

                    ptr_buffer_xxz[k] += gr_x * gr_x * gr_z * fss;

                    ptr_buffer_xyy[k] += gr_x * gr_y * gr_y * fss;

                    ptr_buffer_xyz[k] += gr_x * gr_y * gr_z * fss;

                    ptr_buffer_xzz[k] += gr_x * gr_z * gr_z * fss;

                    ptr_buffer_yyy[k] += gr_y * gr_y * gr_y * fss;

                    ptr_buffer_yyz[k] += gr_y * gr_y * gr_z * fss;

                    ptr_buffer_yzz[k] += gr_y * gr_z * gr_z * fss;

                    ptr_buffer_zzz[k] += gr_z * gr_z * gr_z * fss;
                }
            }

            // distribute GTO values into submatrix

            for (int isph = 0; isph < nsph; isph++)
            {
                for (int icomp = 0; icomp < static_cast<int>(sph_cart_coefs[isph].size()); icomp++)
                {
                    auto icart = sph_cart_coefs[isph][icomp].first;
                    auto fcart = sph_cart_coefs[isph][icomp].second;

                    gtoval::distribute(submat, *(buffer_refs[icart]), fcart, isph * nrows + irow);
                }
            }

            irow++;
        }
    }

    return gto_values;
}

auto
get_gga_values_rec_f(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_F * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_F>(isph);
    }

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", 7 * nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x = gto_values.sub_matrix({1, 0});
        
        auto submat_y = gto_values.sub_matrix({1, 1});
        
        auto submat_z = gto_values.sub_matrix({1, 2});
        
        // compute GTO values for F type GTOs on grid
        
        std::vector<double> buffer_0_xxx(ncols);
        std::vector<double> buffer_0_xxy(ncols);
        std::vector<double> buffer_0_xxz(ncols);
        std::vector<double> buffer_0_xyy(ncols);
        std::vector<double> buffer_0_xyz(ncols);
        std::vector<double> buffer_0_xzz(ncols);
        std::vector<double> buffer_0_yyy(ncols);
        std::vector<double> buffer_0_yyz(ncols);
        std::vector<double> buffer_0_yzz(ncols);
        std::vector<double> buffer_0_zzz(ncols);

        std::vector<double> buffer_x_xxx(ncols);
        std::vector<double> buffer_x_xxy(ncols);
        std::vector<double> buffer_x_xxz(ncols);
        std::vector<double> buffer_x_xyy(ncols);
        std::vector<double> buffer_x_xyz(ncols);
        std::vector<double> buffer_x_xzz(ncols);
        std::vector<double> buffer_x_yyy(ncols);
        std::vector<double> buffer_x_yyz(ncols);
        std::vector<double> buffer_x_yzz(ncols);
        std::vector<double> buffer_x_zzz(ncols);

        std::vector<double> buffer_y_xxx(ncols);
        std::vector<double> buffer_y_xxy(ncols);
        std::vector<double> buffer_y_xxz(ncols);
        std::vector<double> buffer_y_xyy(ncols);
        std::vector<double> buffer_y_xyz(ncols);
        std::vector<double> buffer_y_xzz(ncols);
        std::vector<double> buffer_y_yyy(ncols);
        std::vector<double> buffer_y_yyz(ncols);
        std::vector<double> buffer_y_yzz(ncols);
        std::vector<double> buffer_y_zzz(ncols);

        std::vector<double> buffer_z_xxx(ncols);
        std::vector<double> buffer_z_xxy(ncols);
        std::vector<double> buffer_z_xxz(ncols);
        std::vector<double> buffer_z_xyy(ncols);
        std::vector<double> buffer_z_xyz(ncols);
        std::vector<double> buffer_z_xzz(ncols);
        std::vector<double> buffer_z_yyy(ncols);
        std::vector<double> buffer_z_yyz(ncols);
        std::vector<double> buffer_z_yzz(ncols);
        std::vector<double> buffer_z_zzz(ncols);

        std::vector<std::vector<double>*> buffer_0_refs({
            &buffer_0_xxx, &buffer_0_xxy, &buffer_0_xxz, &buffer_0_xyy, &buffer_0_xyz, &buffer_0_xzz,
            &buffer_0_yyy, &buffer_0_yyz, &buffer_0_yzz, &buffer_0_zzz});

        std::vector<std::vector<double>*> buffer_x_refs({
            &buffer_x_xxx, &buffer_x_xxy, &buffer_x_xxz, &buffer_x_xyy, &buffer_x_xyz, &buffer_x_xzz,
            &buffer_x_yyy, &buffer_x_yyz, &buffer_x_yzz, &buffer_x_zzz});
        std::vector<std::vector<double>*> buffer_y_refs({
            &buffer_y_xxx, &buffer_y_xxy, &buffer_y_xxz, &buffer_y_xyy, &buffer_y_xyz, &buffer_y_xzz,
            &buffer_y_yyy, &buffer_y_yyz, &buffer_y_yzz, &buffer_y_zzz});
        std::vector<std::vector<double>*> buffer_z_refs({
            &buffer_z_xxx, &buffer_z_xxy, &buffer_z_xxz, &buffer_z_xyy, &buffer_z_xyz, &buffer_z_xzz,
            &buffer_z_yyy, &buffer_z_yyz, &buffer_z_yzz, &buffer_z_zzz});

        auto ptr_buffer_0_xxx = buffer_0_xxx.data();
        auto ptr_buffer_0_xxy = buffer_0_xxy.data();
        auto ptr_buffer_0_xxz = buffer_0_xxz.data();
        auto ptr_buffer_0_xyy = buffer_0_xyy.data();
        auto ptr_buffer_0_xyz = buffer_0_xyz.data();
        auto ptr_buffer_0_xzz = buffer_0_xzz.data();
        auto ptr_buffer_0_yyy = buffer_0_yyy.data();
        auto ptr_buffer_0_yyz = buffer_0_yyz.data();
        auto ptr_buffer_0_yzz = buffer_0_yzz.data();
        auto ptr_buffer_0_zzz = buffer_0_zzz.data();

        auto ptr_buffer_x_xxx = buffer_x_xxx.data();
        auto ptr_buffer_x_xxy = buffer_x_xxy.data();
        auto ptr_buffer_x_xxz = buffer_x_xxz.data();
        auto ptr_buffer_x_xyy = buffer_x_xyy.data();
        auto ptr_buffer_x_xyz = buffer_x_xyz.data();
        auto ptr_buffer_x_xzz = buffer_x_xzz.data();
        auto ptr_buffer_x_yyy = buffer_x_yyy.data();
        auto ptr_buffer_x_yyz = buffer_x_yyz.data();
        auto ptr_buffer_x_yzz = buffer_x_yzz.data();
        auto ptr_buffer_x_zzz = buffer_x_zzz.data();

        auto ptr_buffer_y_xxx = buffer_y_xxx.data();
        auto ptr_buffer_y_xxy = buffer_y_xxy.data();
        auto ptr_buffer_y_xxz = buffer_y_xxz.data();
        auto ptr_buffer_y_xyy = buffer_y_xyy.data();
        auto ptr_buffer_y_xyz = buffer_y_xyz.data();
        auto ptr_buffer_y_xzz = buffer_y_xzz.data();
        auto ptr_buffer_y_yyy = buffer_y_yyy.data();
        auto ptr_buffer_y_yyz = buffer_y_yyz.data();
        auto ptr_buffer_y_yzz = buffer_y_yzz.data();
        auto ptr_buffer_y_zzz = buffer_y_zzz.data();

        auto ptr_buffer_z_xxx = buffer_z_xxx.data();
        auto ptr_buffer_z_xxy = buffer_z_xxy.data();
        auto ptr_buffer_z_xxz = buffer_z_xxz.data();
        auto ptr_buffer_z_xyy = buffer_z_xyy.data();
        auto ptr_buffer_z_xyz = buffer_z_xyz.data();
        auto ptr_buffer_z_xzz = buffer_z_xzz.data();
        auto ptr_buffer_z_yyy = buffer_z_yyy.data();
        auto ptr_buffer_z_yyz = buffer_z_yyz.data();
        auto ptr_buffer_z_yzz = buffer_z_yzz.data();
        auto ptr_buffer_z_zzz = buffer_z_zzz.data();

        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];

                std::ranges::fill(buffer_0_xxx, 0.0);
                std::ranges::fill(buffer_0_xxy, 0.0);
                std::ranges::fill(buffer_0_xxz, 0.0);
                std::ranges::fill(buffer_0_xyy, 0.0);
                std::ranges::fill(buffer_0_xyz, 0.0);
                std::ranges::fill(buffer_0_xzz, 0.0);
                std::ranges::fill(buffer_0_yyy, 0.0);
                std::ranges::fill(buffer_0_yyz, 0.0);
                std::ranges::fill(buffer_0_yzz, 0.0);
                std::ranges::fill(buffer_0_zzz, 0.0);

                std::ranges::fill(buffer_x_xxx, 0.0);
                std::ranges::fill(buffer_x_xxy, 0.0);
                std::ranges::fill(buffer_x_xxz, 0.0);
                std::ranges::fill(buffer_x_xyy, 0.0);
                std::ranges::fill(buffer_x_xyz, 0.0);
                std::ranges::fill(buffer_x_xzz, 0.0);
                std::ranges::fill(buffer_x_yyy, 0.0);
                std::ranges::fill(buffer_x_yyz, 0.0);
                std::ranges::fill(buffer_x_yzz, 0.0);
                std::ranges::fill(buffer_x_zzz, 0.0);

                std::ranges::fill(buffer_y_xxx, 0.0);
                std::ranges::fill(buffer_y_xxy, 0.0);
                std::ranges::fill(buffer_y_xxz, 0.0);
                std::ranges::fill(buffer_y_xyy, 0.0);
                std::ranges::fill(buffer_y_xyz, 0.0);
                std::ranges::fill(buffer_y_xzz, 0.0);
                std::ranges::fill(buffer_y_yyy, 0.0);
                std::ranges::fill(buffer_y_yyz, 0.0);
                std::ranges::fill(buffer_y_yzz, 0.0);
                std::ranges::fill(buffer_y_zzz, 0.0);

                std::ranges::fill(buffer_z_xxx, 0.0);
                std::ranges::fill(buffer_z_xxy, 0.0);
                std::ranges::fill(buffer_z_xxz, 0.0);
                std::ranges::fill(buffer_z_xyy, 0.0);
                std::ranges::fill(buffer_z_xyz, 0.0);
                std::ranges::fill(buffer_z_xzz, 0.0);
                std::ranges::fill(buffer_z_yyy, 0.0);
                std::ranges::fill(buffer_z_yyz, 0.0);
                std::ranges::fill(buffer_z_yzz, 0.0);
                std::ranges::fill(buffer_z_zzz, 0.0);

                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
                        
                        const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const auto fg_0 = -2.0 * fexp;

                        const auto fg_1 = f00 * fg_0;

                        ptr_buffer_0_xxx[k] += f00 * gr_x * gr_x * gr_x;
                        ptr_buffer_x_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_x + 3.0 * f00 * gr_x * gr_x;
                        ptr_buffer_y_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xxy[k] += f00 * gr_x * gr_x * gr_y;
                        ptr_buffer_x_xxy[k] += fg_1 * gr_x * gr_x * gr_x * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_y_xxy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_x * gr_x;
                        ptr_buffer_z_xxy[k] += fg_1 * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xxz[k] += f00 * gr_x * gr_x * gr_z;
                        ptr_buffer_x_xxz[k] += fg_1 * gr_x * gr_x * gr_x * gr_z + 2.0 * f00 * gr_x * gr_z;
                        ptr_buffer_y_xxz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xxz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_x * gr_x;

                        ptr_buffer_0_xyy[k] += f00 * gr_x * gr_y * gr_y;
                        ptr_buffer_x_xyy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_y * gr_y;
                        ptr_buffer_y_xyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_z_xyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_0_xyz[k] += f00 * gr_x * gr_y * gr_z;
                        ptr_buffer_x_xyz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z + f00 * gr_y * gr_z;
                        ptr_buffer_y_xyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z + f00 * gr_x * gr_z;
                        ptr_buffer_z_xyz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z + f00 * gr_x * gr_y;

                        ptr_buffer_0_xzz[k] += f00 * gr_x * gr_z * gr_z;
                        ptr_buffer_x_xzz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_y_xzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_z_xzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z + 2.0 * f00 * gr_x * gr_z;

                        ptr_buffer_0_yyy[k] += f00 * gr_y * gr_y * gr_y;
                        ptr_buffer_x_yyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_y_yyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_y + 3.0 * f00 * gr_y * gr_y;
                        ptr_buffer_z_yyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yyz[k] += f00 * gr_y * gr_y * gr_z;
                        ptr_buffer_x_yyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_y_yyz[k] += fg_1 * gr_y * gr_y * gr_y * gr_z + 2.0 * f00 * gr_y * gr_z;
                        ptr_buffer_z_yyz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_y * gr_y;

                        ptr_buffer_0_yzz[k] += f00 * gr_y * gr_z * gr_z;
                        ptr_buffer_x_yzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_y_yzz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_z_yzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z + 2.0 * f00 * gr_y * gr_z;

                        ptr_buffer_0_zzz[k] += f00 * gr_z * gr_z * gr_z;
                        ptr_buffer_x_zzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_y_zzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_z_zzz[k] += fg_1 * gr_z * gr_z * gr_z * gr_z + 3.0 * f00 * gr_z * gr_z;
                     }
                }
        
                // distribute GTO values into submatrix
        
                for (int isph = 0; isph < nsph; isph++)
                {
                    for (int icomp = 0; icomp < static_cast<int>(sph_cart_coefs[isph].size()); icomp++)
                    {
                        auto icart = sph_cart_coefs[isph][icomp].first;
                        auto fcart = sph_cart_coefs[isph][icomp].second;

                        gtoval::distribute(submat_0, *(buffer_0_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_x, *(buffer_x_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_y, *(buffer_y_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_z, *(buffer_z_refs[icart]), fcart, isph * nrows + irow);
                    }
                }

                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

auto
get_mgga_values_rec_f(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_F * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_F>(isph);
    }

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("MGGA", 7 * nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x = gto_values.sub_matrix({1, 0});
        auto submat_y = gto_values.sub_matrix({1, 1});
        auto submat_z = gto_values.sub_matrix({1, 2});

        auto submat_xx = gto_values.sub_matrix({2, 0});
        auto submat_xy = gto_values.sub_matrix({2, 1});
        auto submat_xz = gto_values.sub_matrix({2, 2});
        auto submat_yy = gto_values.sub_matrix({2, 3});
        auto submat_yz = gto_values.sub_matrix({2, 4});
        auto submat_zz = gto_values.sub_matrix({2, 5});
        
        // compute GTO values for F type GTOs on grid
        
        std::vector<double> buffer_0_xxx(ncols);
        std::vector<double> buffer_0_xxy(ncols);
        std::vector<double> buffer_0_xxz(ncols);
        std::vector<double> buffer_0_xyy(ncols);
        std::vector<double> buffer_0_xyz(ncols);
        std::vector<double> buffer_0_xzz(ncols);
        std::vector<double> buffer_0_yyy(ncols);
        std::vector<double> buffer_0_yyz(ncols);
        std::vector<double> buffer_0_yzz(ncols);
        std::vector<double> buffer_0_zzz(ncols);

        std::vector<double> buffer_x_xxx(ncols);
        std::vector<double> buffer_x_xxy(ncols);
        std::vector<double> buffer_x_xxz(ncols);
        std::vector<double> buffer_x_xyy(ncols);
        std::vector<double> buffer_x_xyz(ncols);
        std::vector<double> buffer_x_xzz(ncols);
        std::vector<double> buffer_x_yyy(ncols);
        std::vector<double> buffer_x_yyz(ncols);
        std::vector<double> buffer_x_yzz(ncols);
        std::vector<double> buffer_x_zzz(ncols);

        std::vector<double> buffer_y_xxx(ncols);
        std::vector<double> buffer_y_xxy(ncols);
        std::vector<double> buffer_y_xxz(ncols);
        std::vector<double> buffer_y_xyy(ncols);
        std::vector<double> buffer_y_xyz(ncols);
        std::vector<double> buffer_y_xzz(ncols);
        std::vector<double> buffer_y_yyy(ncols);
        std::vector<double> buffer_y_yyz(ncols);
        std::vector<double> buffer_y_yzz(ncols);
        std::vector<double> buffer_y_zzz(ncols);

        std::vector<double> buffer_z_xxx(ncols);
        std::vector<double> buffer_z_xxy(ncols);
        std::vector<double> buffer_z_xxz(ncols);
        std::vector<double> buffer_z_xyy(ncols);
        std::vector<double> buffer_z_xyz(ncols);
        std::vector<double> buffer_z_xzz(ncols);
        std::vector<double> buffer_z_yyy(ncols);
        std::vector<double> buffer_z_yyz(ncols);
        std::vector<double> buffer_z_yzz(ncols);
        std::vector<double> buffer_z_zzz(ncols);

        std::vector<double> buffer_xx_xxx(ncols);
        std::vector<double> buffer_xx_xxy(ncols);
        std::vector<double> buffer_xx_xxz(ncols);
        std::vector<double> buffer_xx_xyy(ncols);
        std::vector<double> buffer_xx_xyz(ncols);
        std::vector<double> buffer_xx_xzz(ncols);
        std::vector<double> buffer_xx_yyy(ncols);
        std::vector<double> buffer_xx_yyz(ncols);
        std::vector<double> buffer_xx_yzz(ncols);
        std::vector<double> buffer_xx_zzz(ncols);

        std::vector<double> buffer_xy_xxx(ncols);
        std::vector<double> buffer_xy_xxy(ncols);
        std::vector<double> buffer_xy_xxz(ncols);
        std::vector<double> buffer_xy_xyy(ncols);
        std::vector<double> buffer_xy_xyz(ncols);
        std::vector<double> buffer_xy_xzz(ncols);
        std::vector<double> buffer_xy_yyy(ncols);
        std::vector<double> buffer_xy_yyz(ncols);
        std::vector<double> buffer_xy_yzz(ncols);
        std::vector<double> buffer_xy_zzz(ncols);

        std::vector<double> buffer_xz_xxx(ncols);
        std::vector<double> buffer_xz_xxy(ncols);
        std::vector<double> buffer_xz_xxz(ncols);
        std::vector<double> buffer_xz_xyy(ncols);
        std::vector<double> buffer_xz_xyz(ncols);
        std::vector<double> buffer_xz_xzz(ncols);
        std::vector<double> buffer_xz_yyy(ncols);
        std::vector<double> buffer_xz_yyz(ncols);
        std::vector<double> buffer_xz_yzz(ncols);
        std::vector<double> buffer_xz_zzz(ncols);

        std::vector<double> buffer_yy_xxx(ncols);
        std::vector<double> buffer_yy_xxy(ncols);
        std::vector<double> buffer_yy_xxz(ncols);
        std::vector<double> buffer_yy_xyy(ncols);
        std::vector<double> buffer_yy_xyz(ncols);
        std::vector<double> buffer_yy_xzz(ncols);
        std::vector<double> buffer_yy_yyy(ncols);
        std::vector<double> buffer_yy_yyz(ncols);
        std::vector<double> buffer_yy_yzz(ncols);
        std::vector<double> buffer_yy_zzz(ncols);

        std::vector<double> buffer_yz_xxx(ncols);
        std::vector<double> buffer_yz_xxy(ncols);
        std::vector<double> buffer_yz_xxz(ncols);
        std::vector<double> buffer_yz_xyy(ncols);
        std::vector<double> buffer_yz_xyz(ncols);
        std::vector<double> buffer_yz_xzz(ncols);
        std::vector<double> buffer_yz_yyy(ncols);
        std::vector<double> buffer_yz_yyz(ncols);
        std::vector<double> buffer_yz_yzz(ncols);
        std::vector<double> buffer_yz_zzz(ncols);

        std::vector<double> buffer_zz_xxx(ncols);
        std::vector<double> buffer_zz_xxy(ncols);
        std::vector<double> buffer_zz_xxz(ncols);
        std::vector<double> buffer_zz_xyy(ncols);
        std::vector<double> buffer_zz_xyz(ncols);
        std::vector<double> buffer_zz_xzz(ncols);
        std::vector<double> buffer_zz_yyy(ncols);
        std::vector<double> buffer_zz_yyz(ncols);
        std::vector<double> buffer_zz_yzz(ncols);
        std::vector<double> buffer_zz_zzz(ncols);

        std::vector<std::vector<double>*> buffer_0_refs({
            &buffer_0_xxx, &buffer_0_xxy, &buffer_0_xxz, &buffer_0_xyy, &buffer_0_xyz, &buffer_0_xzz,
            &buffer_0_yyy, &buffer_0_yyz, &buffer_0_yzz, &buffer_0_zzz});

        std::vector<std::vector<double>*> buffer_x_refs({
            &buffer_x_xxx, &buffer_x_xxy, &buffer_x_xxz, &buffer_x_xyy, &buffer_x_xyz, &buffer_x_xzz,
            &buffer_x_yyy, &buffer_x_yyz, &buffer_x_yzz, &buffer_x_zzz});
        std::vector<std::vector<double>*> buffer_y_refs({
            &buffer_y_xxx, &buffer_y_xxy, &buffer_y_xxz, &buffer_y_xyy, &buffer_y_xyz, &buffer_y_xzz,
            &buffer_y_yyy, &buffer_y_yyz, &buffer_y_yzz, &buffer_y_zzz});
        std::vector<std::vector<double>*> buffer_z_refs({
            &buffer_z_xxx, &buffer_z_xxy, &buffer_z_xxz, &buffer_z_xyy, &buffer_z_xyz, &buffer_z_xzz,
            &buffer_z_yyy, &buffer_z_yyz, &buffer_z_yzz, &buffer_z_zzz});

        std::vector<std::vector<double>*> buffer_xx_refs({
            &buffer_xx_xxx, &buffer_xx_xxy, &buffer_xx_xxz, &buffer_xx_xyy, &buffer_xx_xyz, &buffer_xx_xzz,
            &buffer_xx_yyy, &buffer_xx_yyz, &buffer_xx_yzz, &buffer_xx_zzz});
        std::vector<std::vector<double>*> buffer_xy_refs({
            &buffer_xy_xxx, &buffer_xy_xxy, &buffer_xy_xxz, &buffer_xy_xyy, &buffer_xy_xyz, &buffer_xy_xzz,
            &buffer_xy_yyy, &buffer_xy_yyz, &buffer_xy_yzz, &buffer_xy_zzz});
        std::vector<std::vector<double>*> buffer_xz_refs({
            &buffer_xz_xxx, &buffer_xz_xxy, &buffer_xz_xxz, &buffer_xz_xyy, &buffer_xz_xyz, &buffer_xz_xzz,
            &buffer_xz_yyy, &buffer_xz_yyz, &buffer_xz_yzz, &buffer_xz_zzz});
        std::vector<std::vector<double>*> buffer_yy_refs({
            &buffer_yy_xxx, &buffer_yy_xxy, &buffer_yy_xxz, &buffer_yy_xyy, &buffer_yy_xyz, &buffer_yy_xzz,
            &buffer_yy_yyy, &buffer_yy_yyz, &buffer_yy_yzz, &buffer_yy_zzz});
        std::vector<std::vector<double>*> buffer_yz_refs({
            &buffer_yz_xxx, &buffer_yz_xxy, &buffer_yz_xxz, &buffer_yz_xyy, &buffer_yz_xyz, &buffer_yz_xzz,
            &buffer_yz_yyy, &buffer_yz_yyz, &buffer_yz_yzz, &buffer_yz_zzz});
        std::vector<std::vector<double>*> buffer_zz_refs({
            &buffer_zz_xxx, &buffer_zz_xxy, &buffer_zz_xxz, &buffer_zz_xyy, &buffer_zz_xyz, &buffer_zz_xzz,
            &buffer_zz_yyy, &buffer_zz_yyz, &buffer_zz_yzz, &buffer_zz_zzz});

        auto ptr_buffer_0_xxx = buffer_0_xxx.data();
        auto ptr_buffer_0_xxy = buffer_0_xxy.data();
        auto ptr_buffer_0_xxz = buffer_0_xxz.data();
        auto ptr_buffer_0_xyy = buffer_0_xyy.data();
        auto ptr_buffer_0_xyz = buffer_0_xyz.data();
        auto ptr_buffer_0_xzz = buffer_0_xzz.data();
        auto ptr_buffer_0_yyy = buffer_0_yyy.data();
        auto ptr_buffer_0_yyz = buffer_0_yyz.data();
        auto ptr_buffer_0_yzz = buffer_0_yzz.data();
        auto ptr_buffer_0_zzz = buffer_0_zzz.data();

        auto ptr_buffer_x_xxx = buffer_x_xxx.data();
        auto ptr_buffer_x_xxy = buffer_x_xxy.data();
        auto ptr_buffer_x_xxz = buffer_x_xxz.data();
        auto ptr_buffer_x_xyy = buffer_x_xyy.data();
        auto ptr_buffer_x_xyz = buffer_x_xyz.data();
        auto ptr_buffer_x_xzz = buffer_x_xzz.data();
        auto ptr_buffer_x_yyy = buffer_x_yyy.data();
        auto ptr_buffer_x_yyz = buffer_x_yyz.data();
        auto ptr_buffer_x_yzz = buffer_x_yzz.data();
        auto ptr_buffer_x_zzz = buffer_x_zzz.data();

        auto ptr_buffer_y_xxx = buffer_y_xxx.data();
        auto ptr_buffer_y_xxy = buffer_y_xxy.data();
        auto ptr_buffer_y_xxz = buffer_y_xxz.data();
        auto ptr_buffer_y_xyy = buffer_y_xyy.data();
        auto ptr_buffer_y_xyz = buffer_y_xyz.data();
        auto ptr_buffer_y_xzz = buffer_y_xzz.data();
        auto ptr_buffer_y_yyy = buffer_y_yyy.data();
        auto ptr_buffer_y_yyz = buffer_y_yyz.data();
        auto ptr_buffer_y_yzz = buffer_y_yzz.data();
        auto ptr_buffer_y_zzz = buffer_y_zzz.data();

        auto ptr_buffer_z_xxx = buffer_z_xxx.data();
        auto ptr_buffer_z_xxy = buffer_z_xxy.data();
        auto ptr_buffer_z_xxz = buffer_z_xxz.data();
        auto ptr_buffer_z_xyy = buffer_z_xyy.data();
        auto ptr_buffer_z_xyz = buffer_z_xyz.data();
        auto ptr_buffer_z_xzz = buffer_z_xzz.data();
        auto ptr_buffer_z_yyy = buffer_z_yyy.data();
        auto ptr_buffer_z_yyz = buffer_z_yyz.data();
        auto ptr_buffer_z_yzz = buffer_z_yzz.data();
        auto ptr_buffer_z_zzz = buffer_z_zzz.data();

        auto ptr_buffer_xx_xxx = buffer_xx_xxx.data();
        auto ptr_buffer_xx_xxy = buffer_xx_xxy.data();
        auto ptr_buffer_xx_xxz = buffer_xx_xxz.data();
        auto ptr_buffer_xx_xyy = buffer_xx_xyy.data();
        auto ptr_buffer_xx_xyz = buffer_xx_xyz.data();
        auto ptr_buffer_xx_xzz = buffer_xx_xzz.data();
        auto ptr_buffer_xx_yyy = buffer_xx_yyy.data();
        auto ptr_buffer_xx_yyz = buffer_xx_yyz.data();
        auto ptr_buffer_xx_yzz = buffer_xx_yzz.data();
        auto ptr_buffer_xx_zzz = buffer_xx_zzz.data();

        auto ptr_buffer_xy_xxx = buffer_xy_xxx.data();
        auto ptr_buffer_xy_xxy = buffer_xy_xxy.data();
        auto ptr_buffer_xy_xxz = buffer_xy_xxz.data();
        auto ptr_buffer_xy_xyy = buffer_xy_xyy.data();
        auto ptr_buffer_xy_xyz = buffer_xy_xyz.data();
        auto ptr_buffer_xy_xzz = buffer_xy_xzz.data();
        auto ptr_buffer_xy_yyy = buffer_xy_yyy.data();
        auto ptr_buffer_xy_yyz = buffer_xy_yyz.data();
        auto ptr_buffer_xy_yzz = buffer_xy_yzz.data();
        auto ptr_buffer_xy_zzz = buffer_xy_zzz.data();

        auto ptr_buffer_xz_xxx = buffer_xz_xxx.data();
        auto ptr_buffer_xz_xxy = buffer_xz_xxy.data();
        auto ptr_buffer_xz_xxz = buffer_xz_xxz.data();
        auto ptr_buffer_xz_xyy = buffer_xz_xyy.data();
        auto ptr_buffer_xz_xyz = buffer_xz_xyz.data();
        auto ptr_buffer_xz_xzz = buffer_xz_xzz.data();
        auto ptr_buffer_xz_yyy = buffer_xz_yyy.data();
        auto ptr_buffer_xz_yyz = buffer_xz_yyz.data();
        auto ptr_buffer_xz_yzz = buffer_xz_yzz.data();
        auto ptr_buffer_xz_zzz = buffer_xz_zzz.data();

        auto ptr_buffer_yy_xxx = buffer_yy_xxx.data();
        auto ptr_buffer_yy_xxy = buffer_yy_xxy.data();
        auto ptr_buffer_yy_xxz = buffer_yy_xxz.data();
        auto ptr_buffer_yy_xyy = buffer_yy_xyy.data();
        auto ptr_buffer_yy_xyz = buffer_yy_xyz.data();
        auto ptr_buffer_yy_xzz = buffer_yy_xzz.data();
        auto ptr_buffer_yy_yyy = buffer_yy_yyy.data();
        auto ptr_buffer_yy_yyz = buffer_yy_yyz.data();
        auto ptr_buffer_yy_yzz = buffer_yy_yzz.data();
        auto ptr_buffer_yy_zzz = buffer_yy_zzz.data();

        auto ptr_buffer_yz_xxx = buffer_yz_xxx.data();
        auto ptr_buffer_yz_xxy = buffer_yz_xxy.data();
        auto ptr_buffer_yz_xxz = buffer_yz_xxz.data();
        auto ptr_buffer_yz_xyy = buffer_yz_xyy.data();
        auto ptr_buffer_yz_xyz = buffer_yz_xyz.data();
        auto ptr_buffer_yz_xzz = buffer_yz_xzz.data();
        auto ptr_buffer_yz_yyy = buffer_yz_yyy.data();
        auto ptr_buffer_yz_yyz = buffer_yz_yyz.data();
        auto ptr_buffer_yz_yzz = buffer_yz_yzz.data();
        auto ptr_buffer_yz_zzz = buffer_yz_zzz.data();

        auto ptr_buffer_zz_xxx = buffer_zz_xxx.data();
        auto ptr_buffer_zz_xxy = buffer_zz_xxy.data();
        auto ptr_buffer_zz_xxz = buffer_zz_xxz.data();
        auto ptr_buffer_zz_xyy = buffer_zz_xyy.data();
        auto ptr_buffer_zz_xyz = buffer_zz_xyz.data();
        auto ptr_buffer_zz_xzz = buffer_zz_xzz.data();
        auto ptr_buffer_zz_yyy = buffer_zz_yyy.data();
        auto ptr_buffer_zz_yyz = buffer_zz_yyz.data();
        auto ptr_buffer_zz_yzz = buffer_zz_yzz.data();
        auto ptr_buffer_zz_zzz = buffer_zz_zzz.data();

        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];

                std::ranges::fill(buffer_0_xxx, 0.0);
                std::ranges::fill(buffer_0_xxy, 0.0);
                std::ranges::fill(buffer_0_xxz, 0.0);
                std::ranges::fill(buffer_0_xyy, 0.0);
                std::ranges::fill(buffer_0_xyz, 0.0);
                std::ranges::fill(buffer_0_xzz, 0.0);
                std::ranges::fill(buffer_0_yyy, 0.0);
                std::ranges::fill(buffer_0_yyz, 0.0);
                std::ranges::fill(buffer_0_yzz, 0.0);
                std::ranges::fill(buffer_0_zzz, 0.0);

                std::ranges::fill(buffer_x_xxx, 0.0);
                std::ranges::fill(buffer_x_xxy, 0.0);
                std::ranges::fill(buffer_x_xxz, 0.0);
                std::ranges::fill(buffer_x_xyy, 0.0);
                std::ranges::fill(buffer_x_xyz, 0.0);
                std::ranges::fill(buffer_x_xzz, 0.0);
                std::ranges::fill(buffer_x_yyy, 0.0);
                std::ranges::fill(buffer_x_yyz, 0.0);
                std::ranges::fill(buffer_x_yzz, 0.0);
                std::ranges::fill(buffer_x_zzz, 0.0);

                std::ranges::fill(buffer_y_xxx, 0.0);
                std::ranges::fill(buffer_y_xxy, 0.0);
                std::ranges::fill(buffer_y_xxz, 0.0);
                std::ranges::fill(buffer_y_xyy, 0.0);
                std::ranges::fill(buffer_y_xyz, 0.0);
                std::ranges::fill(buffer_y_xzz, 0.0);
                std::ranges::fill(buffer_y_yyy, 0.0);
                std::ranges::fill(buffer_y_yyz, 0.0);
                std::ranges::fill(buffer_y_yzz, 0.0);
                std::ranges::fill(buffer_y_zzz, 0.0);

                std::ranges::fill(buffer_z_xxx, 0.0);
                std::ranges::fill(buffer_z_xxy, 0.0);
                std::ranges::fill(buffer_z_xxz, 0.0);
                std::ranges::fill(buffer_z_xyy, 0.0);
                std::ranges::fill(buffer_z_xyz, 0.0);
                std::ranges::fill(buffer_z_xzz, 0.0);
                std::ranges::fill(buffer_z_yyy, 0.0);
                std::ranges::fill(buffer_z_yyz, 0.0);
                std::ranges::fill(buffer_z_yzz, 0.0);
                std::ranges::fill(buffer_z_zzz, 0.0);

                std::ranges::fill(buffer_xx_xxx, 0.0);
                std::ranges::fill(buffer_xx_xxy, 0.0);
                std::ranges::fill(buffer_xx_xxz, 0.0);
                std::ranges::fill(buffer_xx_xyy, 0.0);
                std::ranges::fill(buffer_xx_xyz, 0.0);
                std::ranges::fill(buffer_xx_xzz, 0.0);
                std::ranges::fill(buffer_xx_yyy, 0.0);
                std::ranges::fill(buffer_xx_yyz, 0.0);
                std::ranges::fill(buffer_xx_yzz, 0.0);
                std::ranges::fill(buffer_xx_zzz, 0.0);

                std::ranges::fill(buffer_xy_xxx, 0.0);
                std::ranges::fill(buffer_xy_xxy, 0.0);
                std::ranges::fill(buffer_xy_xxz, 0.0);
                std::ranges::fill(buffer_xy_xyy, 0.0);
                std::ranges::fill(buffer_xy_xyz, 0.0);
                std::ranges::fill(buffer_xy_xzz, 0.0);
                std::ranges::fill(buffer_xy_yyy, 0.0);
                std::ranges::fill(buffer_xy_yyz, 0.0);
                std::ranges::fill(buffer_xy_yzz, 0.0);
                std::ranges::fill(buffer_xy_zzz, 0.0);

                std::ranges::fill(buffer_xz_xxx, 0.0);
                std::ranges::fill(buffer_xz_xxy, 0.0);
                std::ranges::fill(buffer_xz_xxz, 0.0);
                std::ranges::fill(buffer_xz_xyy, 0.0);
                std::ranges::fill(buffer_xz_xyz, 0.0);
                std::ranges::fill(buffer_xz_xzz, 0.0);
                std::ranges::fill(buffer_xz_yyy, 0.0);
                std::ranges::fill(buffer_xz_yyz, 0.0);
                std::ranges::fill(buffer_xz_yzz, 0.0);
                std::ranges::fill(buffer_xz_zzz, 0.0);

                std::ranges::fill(buffer_yy_xxx, 0.0);
                std::ranges::fill(buffer_yy_xxy, 0.0);
                std::ranges::fill(buffer_yy_xxz, 0.0);
                std::ranges::fill(buffer_yy_xyy, 0.0);
                std::ranges::fill(buffer_yy_xyz, 0.0);
                std::ranges::fill(buffer_yy_xzz, 0.0);
                std::ranges::fill(buffer_yy_yyy, 0.0);
                std::ranges::fill(buffer_yy_yyz, 0.0);
                std::ranges::fill(buffer_yy_yzz, 0.0);
                std::ranges::fill(buffer_yy_zzz, 0.0);

                std::ranges::fill(buffer_yz_xxx, 0.0);
                std::ranges::fill(buffer_yz_xxy, 0.0);
                std::ranges::fill(buffer_yz_xxz, 0.0);
                std::ranges::fill(buffer_yz_xyy, 0.0);
                std::ranges::fill(buffer_yz_xyz, 0.0);
                std::ranges::fill(buffer_yz_xzz, 0.0);
                std::ranges::fill(buffer_yz_yyy, 0.0);
                std::ranges::fill(buffer_yz_yyz, 0.0);
                std::ranges::fill(buffer_yz_yzz, 0.0);
                std::ranges::fill(buffer_yz_zzz, 0.0);

                std::ranges::fill(buffer_zz_xxx, 0.0);
                std::ranges::fill(buffer_zz_xxy, 0.0);
                std::ranges::fill(buffer_zz_xxz, 0.0);
                std::ranges::fill(buffer_zz_xyy, 0.0);
                std::ranges::fill(buffer_zz_xyz, 0.0);
                std::ranges::fill(buffer_zz_xzz, 0.0);
                std::ranges::fill(buffer_zz_yyy, 0.0);
                std::ranges::fill(buffer_zz_yyz, 0.0);
                std::ranges::fill(buffer_zz_yzz, 0.0);
                std::ranges::fill(buffer_zz_zzz, 0.0);

                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
                        
                        const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const auto fg0 = -2.0 * fexp;

                        const auto fg1 = f00 * fg0;

                        const auto fg2 = fg1 * fg0;

                        ptr_buffer_0_xxx[k] += f00 * gr_x * gr_x * gr_x;

                        ptr_buffer_x_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_x + 3.0 * f00 * gr_x * gr_x;
                        ptr_buffer_y_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_xx_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_x + 7.0 * fg1 * gr_x * gr_x * gr_x + 6.0 * f00 * gr_x;
                        ptr_buffer_xy_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_y + 3.0 * fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_xz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_z + 3.0 * fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_yy_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x * gr_x;
                        ptr_buffer_yz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x * gr_x;

                        ptr_buffer_0_xxy[k] += f00 * gr_x * gr_x * gr_y;

                        ptr_buffer_x_xxy[k] += fg1 * gr_x * gr_x * gr_x * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_y_xxy[k] += fg1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_x * gr_x;
                        ptr_buffer_z_xxy[k] += fg1 * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_xx_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_y + 5.0 * fg1 * gr_x * gr_x * gr_y + 2.0 * f00 * gr_y;
                        ptr_buffer_xy_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x * gr_x + 2.0 * fg1 * gr_x * gr_y * gr_y + 2.0 * f00 * gr_x;
                        ptr_buffer_xz_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_yz_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_zz_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y;

                        ptr_buffer_0_xxz[k] += f00 * gr_x * gr_x * gr_z;

                        ptr_buffer_x_xxz[k] += fg1 * gr_x * gr_x * gr_x * gr_z + 2.0 * f00 * gr_x * gr_z;
                        ptr_buffer_y_xxz[k] += fg1 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xxz[k] += fg1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_x * gr_x;

                        ptr_buffer_xx_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_z + 5.0 * fg1 * gr_x * gr_x * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_xy_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x * gr_x + 2.0 * fg1 * gr_x * gr_z * gr_z + 2.0 * f00 * gr_x;
                        ptr_buffer_yy_xxz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_yz_xxz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_zz_xxz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xyy[k] += f00 * gr_x * gr_y * gr_y;

                        ptr_buffer_x_xyy[k] += fg1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_y * gr_y;
                        ptr_buffer_y_xyy[k] += fg1 * gr_x * gr_y * gr_y * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_z_xyy[k] += fg1 * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_xx_xyy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_xy_xyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_x * gr_x * gr_y + fg1 * gr_y * gr_y * gr_y + 2.0 * f00 * gr_y;
                        ptr_buffer_xz_xyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_yy_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_y + 5.0 * fg1 * gr_x * gr_y * gr_y + 2.0 * f00 * gr_x;
                        ptr_buffer_yz_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y;

                        ptr_buffer_0_xyz[k] += f00 * gr_x * gr_y * gr_z;

                        ptr_buffer_x_xyz[k] += fg1 * gr_x * gr_x * gr_y * gr_z + f00 * gr_y * gr_z;
                        ptr_buffer_y_xyz[k] += fg1 * gr_x * gr_y * gr_y * gr_z + f00 * gr_x * gr_z;
                        ptr_buffer_z_xyz[k] += fg1 * gr_x * gr_y * gr_z * gr_z + f00 * gr_x * gr_y;

                        ptr_buffer_xx_xyz[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xy_xyz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z + fg1 * gr_y * gr_y * gr_z + f00 * gr_z;
                        ptr_buffer_xz_xyz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y + fg1 * gr_y * gr_z * gr_z + f00 * gr_y;
                        ptr_buffer_yy_xyz[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yz_xyz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y + fg1 * gr_x * gr_z * gr_z + f00 * gr_x;
                        ptr_buffer_zz_xyz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xzz[k] += f00 * gr_x * gr_z * gr_z;

                        ptr_buffer_x_xzz[k] += fg1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_y_xzz[k] += fg1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_z_xzz[k] += fg1 * gr_x * gr_z * gr_z * gr_z + 2.0 * f00 * gr_x * gr_z;

                        ptr_buffer_xx_xzz[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_xy_xzz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_xz_xzz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_x * gr_z + fg1 * gr_z * gr_z * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_yy_xzz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_yz_xzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xzz[k] += fg2 * gr_x * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_x * gr_z * gr_z + 2.0 * f00 * gr_x;

                        ptr_buffer_0_yyy[k] += f00 * gr_y * gr_y * gr_y;

                        ptr_buffer_x_yyy[k] += fg1 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_y_yyy[k] += fg1 * gr_y * gr_y * gr_y * gr_y + 3.0 * f00 * gr_y * gr_y;
                        ptr_buffer_z_yyy[k] += fg1 * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_xx_yyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + fg1 * gr_y * gr_y * gr_y;
                        ptr_buffer_xy_yyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_xz_yyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z;
                        ptr_buffer_yy_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_y + 7.0 * fg1 * gr_y * gr_y * gr_y + 6.0 * f00 * gr_y;
                        ptr_buffer_yz_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_zz_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y * gr_y;

                        ptr_buffer_0_yyz[k] += f00 * gr_y * gr_y * gr_z;

                        ptr_buffer_x_yyz[k] += fg1 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_y_yyz[k] += fg1 * gr_y * gr_y * gr_y * gr_z + 2.0 * f00 * gr_y * gr_z;
                        ptr_buffer_z_yyz[k] += fg1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_y * gr_y;

                        ptr_buffer_xx_yyz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_xy_yyz[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_yyz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_yy_yyz[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_z + 5.0 * fg1 * gr_y * gr_y * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_yz_yyz[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_y * gr_z * gr_z + 2.0 * f00 * gr_y;
                        ptr_buffer_zz_yyz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yzz[k] += f00 * gr_y * gr_z * gr_z;

                        ptr_buffer_x_yzz[k] += fg1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_y_yzz[k] += fg1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_z_yzz[k] += fg1 * gr_y * gr_z * gr_z * gr_z + 2.0 * f00 * gr_y * gr_z;

                        ptr_buffer_xx_yzz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_xy_yzz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_xz_yzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_yzz[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_yz_yzz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_y * gr_z + fg1 * gr_z * gr_z * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_zz_yzz[k] += fg2 * gr_y * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_y * gr_z * gr_z + 2.0 * f00 * gr_y;

                        ptr_buffer_0_zzz[k] += f00 * gr_z * gr_z * gr_z;

                        ptr_buffer_x_zzz[k] += fg1 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_y_zzz[k] += fg1 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_z_zzz[k] += fg1 * gr_z * gr_z * gr_z * gr_z + 3.0 * f00 * gr_z * gr_z;

                        ptr_buffer_xx_zzz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + fg1 * gr_z * gr_z * gr_z;
                        ptr_buffer_xy_zzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_xz_zzz[k] += fg2 * gr_x * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_yy_zzz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + fg1 * gr_z * gr_z * gr_z;
                        ptr_buffer_yz_zzz[k] += fg2 * gr_y * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_zz_zzz[k] += fg2 * gr_z * gr_z * gr_z * gr_z * gr_z + 7.0 * fg1 * gr_z * gr_z * gr_z + 6.0 * f00 * gr_z;
                    }
                }
        
                // distribute GTO values into submatrix
        
                for (int isph = 0; isph < nsph; isph++)
                {
                    for (int icomp = 0; icomp < static_cast<int>(sph_cart_coefs[isph].size()); icomp++)
                    {
                        auto icart = sph_cart_coefs[isph][icomp].first;
                        auto fcart = sph_cart_coefs[isph][icomp].second;

                        gtoval::distribute(submat_0, *(buffer_0_refs[icart]), fcart, isph * nrows + irow);

                        gtoval::distribute(submat_x, *(buffer_x_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_y, *(buffer_y_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_z, *(buffer_z_refs[icart]), fcart, isph * nrows + irow);

                        gtoval::distribute(submat_xx, *(buffer_xx_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xy, *(buffer_xy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xz, *(buffer_xz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yy, *(buffer_yy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yz, *(buffer_yz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_zz, *(buffer_zz_refs[icart]), fcart, isph * nrows + irow);
                    }
                }

                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

auto
get_3rd_order_values_rec_f(const CGtoBlock&            gto_block,
                           const std::vector<double>&  grid_coords_x,
                           const std::vector<double>&  grid_coords_y,
                           const std::vector<double>&  grid_coords_z,
                           const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_F * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_F>(isph);
    }

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("3RD_ORDER", 7 * nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x = gto_values.sub_matrix({1, 0});
        auto submat_y = gto_values.sub_matrix({1, 1});
        auto submat_z = gto_values.sub_matrix({1, 2});

        auto submat_xx = gto_values.sub_matrix({2, 0});
        auto submat_xy = gto_values.sub_matrix({2, 1});
        auto submat_xz = gto_values.sub_matrix({2, 2});
        auto submat_yy = gto_values.sub_matrix({2, 3});
        auto submat_yz = gto_values.sub_matrix({2, 4});
        auto submat_zz = gto_values.sub_matrix({2, 5});

        auto submat_xxx = gto_values.sub_matrix({3, 0});
        auto submat_xxy = gto_values.sub_matrix({3, 1});
        auto submat_xxz = gto_values.sub_matrix({3, 2});
        auto submat_xyy = gto_values.sub_matrix({3, 3});
        auto submat_xyz = gto_values.sub_matrix({3, 4});
        auto submat_xzz = gto_values.sub_matrix({3, 5});
        auto submat_yyy = gto_values.sub_matrix({3, 6});
        auto submat_yyz = gto_values.sub_matrix({3, 7});
        auto submat_yzz = gto_values.sub_matrix({3, 8});
        auto submat_zzz = gto_values.sub_matrix({3, 9});
        
        // compute GTO values for F type GTOs on grid
        
        std::vector<double> buffer_0_xxx(ncols);
        std::vector<double> buffer_0_xxy(ncols);
        std::vector<double> buffer_0_xxz(ncols);
        std::vector<double> buffer_0_xyy(ncols);
        std::vector<double> buffer_0_xyz(ncols);
        std::vector<double> buffer_0_xzz(ncols);
        std::vector<double> buffer_0_yyy(ncols);
        std::vector<double> buffer_0_yyz(ncols);
        std::vector<double> buffer_0_yzz(ncols);
        std::vector<double> buffer_0_zzz(ncols);

        std::vector<double> buffer_x_xxx(ncols);
        std::vector<double> buffer_x_xxy(ncols);
        std::vector<double> buffer_x_xxz(ncols);
        std::vector<double> buffer_x_xyy(ncols);
        std::vector<double> buffer_x_xyz(ncols);
        std::vector<double> buffer_x_xzz(ncols);
        std::vector<double> buffer_x_yyy(ncols);
        std::vector<double> buffer_x_yyz(ncols);
        std::vector<double> buffer_x_yzz(ncols);
        std::vector<double> buffer_x_zzz(ncols);

        std::vector<double> buffer_y_xxx(ncols);
        std::vector<double> buffer_y_xxy(ncols);
        std::vector<double> buffer_y_xxz(ncols);
        std::vector<double> buffer_y_xyy(ncols);
        std::vector<double> buffer_y_xyz(ncols);
        std::vector<double> buffer_y_xzz(ncols);
        std::vector<double> buffer_y_yyy(ncols);
        std::vector<double> buffer_y_yyz(ncols);
        std::vector<double> buffer_y_yzz(ncols);
        std::vector<double> buffer_y_zzz(ncols);

        std::vector<double> buffer_z_xxx(ncols);
        std::vector<double> buffer_z_xxy(ncols);
        std::vector<double> buffer_z_xxz(ncols);
        std::vector<double> buffer_z_xyy(ncols);
        std::vector<double> buffer_z_xyz(ncols);
        std::vector<double> buffer_z_xzz(ncols);
        std::vector<double> buffer_z_yyy(ncols);
        std::vector<double> buffer_z_yyz(ncols);
        std::vector<double> buffer_z_yzz(ncols);
        std::vector<double> buffer_z_zzz(ncols);

        std::vector<double> buffer_xx_xxx(ncols);
        std::vector<double> buffer_xx_xxy(ncols);
        std::vector<double> buffer_xx_xxz(ncols);
        std::vector<double> buffer_xx_xyy(ncols);
        std::vector<double> buffer_xx_xyz(ncols);
        std::vector<double> buffer_xx_xzz(ncols);
        std::vector<double> buffer_xx_yyy(ncols);
        std::vector<double> buffer_xx_yyz(ncols);
        std::vector<double> buffer_xx_yzz(ncols);
        std::vector<double> buffer_xx_zzz(ncols);

        std::vector<double> buffer_xy_xxx(ncols);
        std::vector<double> buffer_xy_xxy(ncols);
        std::vector<double> buffer_xy_xxz(ncols);
        std::vector<double> buffer_xy_xyy(ncols);
        std::vector<double> buffer_xy_xyz(ncols);
        std::vector<double> buffer_xy_xzz(ncols);
        std::vector<double> buffer_xy_yyy(ncols);
        std::vector<double> buffer_xy_yyz(ncols);
        std::vector<double> buffer_xy_yzz(ncols);
        std::vector<double> buffer_xy_zzz(ncols);

        std::vector<double> buffer_xz_xxx(ncols);
        std::vector<double> buffer_xz_xxy(ncols);
        std::vector<double> buffer_xz_xxz(ncols);
        std::vector<double> buffer_xz_xyy(ncols);
        std::vector<double> buffer_xz_xyz(ncols);
        std::vector<double> buffer_xz_xzz(ncols);
        std::vector<double> buffer_xz_yyy(ncols);
        std::vector<double> buffer_xz_yyz(ncols);
        std::vector<double> buffer_xz_yzz(ncols);
        std::vector<double> buffer_xz_zzz(ncols);

        std::vector<double> buffer_yy_xxx(ncols);
        std::vector<double> buffer_yy_xxy(ncols);
        std::vector<double> buffer_yy_xxz(ncols);
        std::vector<double> buffer_yy_xyy(ncols);
        std::vector<double> buffer_yy_xyz(ncols);
        std::vector<double> buffer_yy_xzz(ncols);
        std::vector<double> buffer_yy_yyy(ncols);
        std::vector<double> buffer_yy_yyz(ncols);
        std::vector<double> buffer_yy_yzz(ncols);
        std::vector<double> buffer_yy_zzz(ncols);

        std::vector<double> buffer_yz_xxx(ncols);
        std::vector<double> buffer_yz_xxy(ncols);
        std::vector<double> buffer_yz_xxz(ncols);
        std::vector<double> buffer_yz_xyy(ncols);
        std::vector<double> buffer_yz_xyz(ncols);
        std::vector<double> buffer_yz_xzz(ncols);
        std::vector<double> buffer_yz_yyy(ncols);
        std::vector<double> buffer_yz_yyz(ncols);
        std::vector<double> buffer_yz_yzz(ncols);
        std::vector<double> buffer_yz_zzz(ncols);

        std::vector<double> buffer_zz_xxx(ncols);
        std::vector<double> buffer_zz_xxy(ncols);
        std::vector<double> buffer_zz_xxz(ncols);
        std::vector<double> buffer_zz_xyy(ncols);
        std::vector<double> buffer_zz_xyz(ncols);
        std::vector<double> buffer_zz_xzz(ncols);
        std::vector<double> buffer_zz_yyy(ncols);
        std::vector<double> buffer_zz_yyz(ncols);
        std::vector<double> buffer_zz_yzz(ncols);
        std::vector<double> buffer_zz_zzz(ncols);

        std::vector<double> buffer_xxx_xxx(ncols);
        std::vector<double> buffer_xxx_xxy(ncols);
        std::vector<double> buffer_xxx_xxz(ncols);
        std::vector<double> buffer_xxx_xyy(ncols);
        std::vector<double> buffer_xxx_xyz(ncols);
        std::vector<double> buffer_xxx_xzz(ncols);
        std::vector<double> buffer_xxx_yyy(ncols);
        std::vector<double> buffer_xxx_yyz(ncols);
        std::vector<double> buffer_xxx_yzz(ncols);
        std::vector<double> buffer_xxx_zzz(ncols);

        std::vector<double> buffer_xxy_xxx(ncols);
        std::vector<double> buffer_xxy_xxy(ncols);
        std::vector<double> buffer_xxy_xxz(ncols);
        std::vector<double> buffer_xxy_xyy(ncols);
        std::vector<double> buffer_xxy_xyz(ncols);
        std::vector<double> buffer_xxy_xzz(ncols);
        std::vector<double> buffer_xxy_yyy(ncols);
        std::vector<double> buffer_xxy_yyz(ncols);
        std::vector<double> buffer_xxy_yzz(ncols);
        std::vector<double> buffer_xxy_zzz(ncols);

        std::vector<double> buffer_xxz_xxx(ncols);
        std::vector<double> buffer_xxz_xxy(ncols);
        std::vector<double> buffer_xxz_xxz(ncols);
        std::vector<double> buffer_xxz_xyy(ncols);
        std::vector<double> buffer_xxz_xyz(ncols);
        std::vector<double> buffer_xxz_xzz(ncols);
        std::vector<double> buffer_xxz_yyy(ncols);
        std::vector<double> buffer_xxz_yyz(ncols);
        std::vector<double> buffer_xxz_yzz(ncols);
        std::vector<double> buffer_xxz_zzz(ncols);

        std::vector<double> buffer_xyy_xxx(ncols);
        std::vector<double> buffer_xyy_xxy(ncols);
        std::vector<double> buffer_xyy_xxz(ncols);
        std::vector<double> buffer_xyy_xyy(ncols);
        std::vector<double> buffer_xyy_xyz(ncols);
        std::vector<double> buffer_xyy_xzz(ncols);
        std::vector<double> buffer_xyy_yyy(ncols);
        std::vector<double> buffer_xyy_yyz(ncols);
        std::vector<double> buffer_xyy_yzz(ncols);
        std::vector<double> buffer_xyy_zzz(ncols);

        std::vector<double> buffer_xyz_xxx(ncols);
        std::vector<double> buffer_xyz_xxy(ncols);
        std::vector<double> buffer_xyz_xxz(ncols);
        std::vector<double> buffer_xyz_xyy(ncols);
        std::vector<double> buffer_xyz_xyz(ncols);
        std::vector<double> buffer_xyz_xzz(ncols);
        std::vector<double> buffer_xyz_yyy(ncols);
        std::vector<double> buffer_xyz_yyz(ncols);
        std::vector<double> buffer_xyz_yzz(ncols);
        std::vector<double> buffer_xyz_zzz(ncols);

        std::vector<double> buffer_xzz_xxx(ncols);
        std::vector<double> buffer_xzz_xxy(ncols);
        std::vector<double> buffer_xzz_xxz(ncols);
        std::vector<double> buffer_xzz_xyy(ncols);
        std::vector<double> buffer_xzz_xyz(ncols);
        std::vector<double> buffer_xzz_xzz(ncols);
        std::vector<double> buffer_xzz_yyy(ncols);
        std::vector<double> buffer_xzz_yyz(ncols);
        std::vector<double> buffer_xzz_yzz(ncols);
        std::vector<double> buffer_xzz_zzz(ncols);

        std::vector<double> buffer_yyy_xxx(ncols);
        std::vector<double> buffer_yyy_xxy(ncols);
        std::vector<double> buffer_yyy_xxz(ncols);
        std::vector<double> buffer_yyy_xyy(ncols);
        std::vector<double> buffer_yyy_xyz(ncols);
        std::vector<double> buffer_yyy_xzz(ncols);
        std::vector<double> buffer_yyy_yyy(ncols);
        std::vector<double> buffer_yyy_yyz(ncols);
        std::vector<double> buffer_yyy_yzz(ncols);
        std::vector<double> buffer_yyy_zzz(ncols);

        std::vector<double> buffer_yyz_xxx(ncols);
        std::vector<double> buffer_yyz_xxy(ncols);
        std::vector<double> buffer_yyz_xxz(ncols);
        std::vector<double> buffer_yyz_xyy(ncols);
        std::vector<double> buffer_yyz_xyz(ncols);
        std::vector<double> buffer_yyz_xzz(ncols);
        std::vector<double> buffer_yyz_yyy(ncols);
        std::vector<double> buffer_yyz_yyz(ncols);
        std::vector<double> buffer_yyz_yzz(ncols);
        std::vector<double> buffer_yyz_zzz(ncols);

        std::vector<double> buffer_yzz_xxx(ncols);
        std::vector<double> buffer_yzz_xxy(ncols);
        std::vector<double> buffer_yzz_xxz(ncols);
        std::vector<double> buffer_yzz_xyy(ncols);
        std::vector<double> buffer_yzz_xyz(ncols);
        std::vector<double> buffer_yzz_xzz(ncols);
        std::vector<double> buffer_yzz_yyy(ncols);
        std::vector<double> buffer_yzz_yyz(ncols);
        std::vector<double> buffer_yzz_yzz(ncols);
        std::vector<double> buffer_yzz_zzz(ncols);

        std::vector<double> buffer_zzz_xxx(ncols);
        std::vector<double> buffer_zzz_xxy(ncols);
        std::vector<double> buffer_zzz_xxz(ncols);
        std::vector<double> buffer_zzz_xyy(ncols);
        std::vector<double> buffer_zzz_xyz(ncols);
        std::vector<double> buffer_zzz_xzz(ncols);
        std::vector<double> buffer_zzz_yyy(ncols);
        std::vector<double> buffer_zzz_yyz(ncols);
        std::vector<double> buffer_zzz_yzz(ncols);
        std::vector<double> buffer_zzz_zzz(ncols);

        std::vector<std::vector<double>*> buffer_0_refs({
            &buffer_0_xxx, &buffer_0_xxy, &buffer_0_xxz, &buffer_0_xyy, &buffer_0_xyz, &buffer_0_xzz,
            &buffer_0_yyy, &buffer_0_yyz, &buffer_0_yzz, &buffer_0_zzz});

        std::vector<std::vector<double>*> buffer_x_refs({
            &buffer_x_xxx, &buffer_x_xxy, &buffer_x_xxz, &buffer_x_xyy, &buffer_x_xyz, &buffer_x_xzz,
            &buffer_x_yyy, &buffer_x_yyz, &buffer_x_yzz, &buffer_x_zzz});
        std::vector<std::vector<double>*> buffer_y_refs({
            &buffer_y_xxx, &buffer_y_xxy, &buffer_y_xxz, &buffer_y_xyy, &buffer_y_xyz, &buffer_y_xzz,
            &buffer_y_yyy, &buffer_y_yyz, &buffer_y_yzz, &buffer_y_zzz});
        std::vector<std::vector<double>*> buffer_z_refs({
            &buffer_z_xxx, &buffer_z_xxy, &buffer_z_xxz, &buffer_z_xyy, &buffer_z_xyz, &buffer_z_xzz,
            &buffer_z_yyy, &buffer_z_yyz, &buffer_z_yzz, &buffer_z_zzz});

        std::vector<std::vector<double>*> buffer_xx_refs({
            &buffer_xx_xxx, &buffer_xx_xxy, &buffer_xx_xxz, &buffer_xx_xyy, &buffer_xx_xyz, &buffer_xx_xzz,
            &buffer_xx_yyy, &buffer_xx_yyz, &buffer_xx_yzz, &buffer_xx_zzz});
        std::vector<std::vector<double>*> buffer_xy_refs({
            &buffer_xy_xxx, &buffer_xy_xxy, &buffer_xy_xxz, &buffer_xy_xyy, &buffer_xy_xyz, &buffer_xy_xzz,
            &buffer_xy_yyy, &buffer_xy_yyz, &buffer_xy_yzz, &buffer_xy_zzz});
        std::vector<std::vector<double>*> buffer_xz_refs({
            &buffer_xz_xxx, &buffer_xz_xxy, &buffer_xz_xxz, &buffer_xz_xyy, &buffer_xz_xyz, &buffer_xz_xzz,
            &buffer_xz_yyy, &buffer_xz_yyz, &buffer_xz_yzz, &buffer_xz_zzz});
        std::vector<std::vector<double>*> buffer_yy_refs({
            &buffer_yy_xxx, &buffer_yy_xxy, &buffer_yy_xxz, &buffer_yy_xyy, &buffer_yy_xyz, &buffer_yy_xzz,
            &buffer_yy_yyy, &buffer_yy_yyz, &buffer_yy_yzz, &buffer_yy_zzz});
        std::vector<std::vector<double>*> buffer_yz_refs({
            &buffer_yz_xxx, &buffer_yz_xxy, &buffer_yz_xxz, &buffer_yz_xyy, &buffer_yz_xyz, &buffer_yz_xzz,
            &buffer_yz_yyy, &buffer_yz_yyz, &buffer_yz_yzz, &buffer_yz_zzz});
        std::vector<std::vector<double>*> buffer_zz_refs({
            &buffer_zz_xxx, &buffer_zz_xxy, &buffer_zz_xxz, &buffer_zz_xyy, &buffer_zz_xyz, &buffer_zz_xzz,
            &buffer_zz_yyy, &buffer_zz_yyz, &buffer_zz_yzz, &buffer_zz_zzz});

        std::vector<std::vector<double>*> buffer_xxx_refs({
            &buffer_xxx_xxx,
            &buffer_xxx_xxy,
            &buffer_xxx_xxz,
            &buffer_xxx_xyy,
            &buffer_xxx_xyz,
            &buffer_xxx_xzz,
            &buffer_xxx_yyy,
            &buffer_xxx_yyz,
            &buffer_xxx_yzz,
            &buffer_xxx_zzz,
        });
        std::vector<std::vector<double>*> buffer_xxy_refs({
            &buffer_xxy_xxx,
            &buffer_xxy_xxy,
            &buffer_xxy_xxz,
            &buffer_xxy_xyy,
            &buffer_xxy_xyz,
            &buffer_xxy_xzz,
            &buffer_xxy_yyy,
            &buffer_xxy_yyz,
            &buffer_xxy_yzz,
            &buffer_xxy_zzz,
        });
        std::vector<std::vector<double>*> buffer_xxz_refs({
            &buffer_xxz_xxx,
            &buffer_xxz_xxy,
            &buffer_xxz_xxz,
            &buffer_xxz_xyy,
            &buffer_xxz_xyz,
            &buffer_xxz_xzz,
            &buffer_xxz_yyy,
            &buffer_xxz_yyz,
            &buffer_xxz_yzz,
            &buffer_xxz_zzz,
        });
        std::vector<std::vector<double>*> buffer_xyy_refs({
            &buffer_xyy_xxx,
            &buffer_xyy_xxy,
            &buffer_xyy_xxz,
            &buffer_xyy_xyy,
            &buffer_xyy_xyz,
            &buffer_xyy_xzz,
            &buffer_xyy_yyy,
            &buffer_xyy_yyz,
            &buffer_xyy_yzz,
            &buffer_xyy_zzz,
        });
        std::vector<std::vector<double>*> buffer_xyz_refs({
            &buffer_xyz_xxx,
            &buffer_xyz_xxy,
            &buffer_xyz_xxz,
            &buffer_xyz_xyy,
            &buffer_xyz_xyz,
            &buffer_xyz_xzz,
            &buffer_xyz_yyy,
            &buffer_xyz_yyz,
            &buffer_xyz_yzz,
            &buffer_xyz_zzz,
        });
        std::vector<std::vector<double>*> buffer_xzz_refs({
            &buffer_xzz_xxx,
            &buffer_xzz_xxy,
            &buffer_xzz_xxz,
            &buffer_xzz_xyy,
            &buffer_xzz_xyz,
            &buffer_xzz_xzz,
            &buffer_xzz_yyy,
            &buffer_xzz_yyz,
            &buffer_xzz_yzz,
            &buffer_xzz_zzz,
        });
        std::vector<std::vector<double>*> buffer_yyy_refs({
            &buffer_yyy_xxx,
            &buffer_yyy_xxy,
            &buffer_yyy_xxz,
            &buffer_yyy_xyy,
            &buffer_yyy_xyz,
            &buffer_yyy_xzz,
            &buffer_yyy_yyy,
            &buffer_yyy_yyz,
            &buffer_yyy_yzz,
            &buffer_yyy_zzz,
        });
        std::vector<std::vector<double>*> buffer_yyz_refs({
            &buffer_yyz_xxx,
            &buffer_yyz_xxy,
            &buffer_yyz_xxz,
            &buffer_yyz_xyy,
            &buffer_yyz_xyz,
            &buffer_yyz_xzz,
            &buffer_yyz_yyy,
            &buffer_yyz_yyz,
            &buffer_yyz_yzz,
            &buffer_yyz_zzz,
        });
        std::vector<std::vector<double>*> buffer_yzz_refs({
            &buffer_yzz_xxx,
            &buffer_yzz_xxy,
            &buffer_yzz_xxz,
            &buffer_yzz_xyy,
            &buffer_yzz_xyz,
            &buffer_yzz_xzz,
            &buffer_yzz_yyy,
            &buffer_yzz_yyz,
            &buffer_yzz_yzz,
            &buffer_yzz_zzz,
        });
        std::vector<std::vector<double>*> buffer_zzz_refs({
            &buffer_zzz_xxx,
            &buffer_zzz_xxy,
            &buffer_zzz_xxz,
            &buffer_zzz_xyy,
            &buffer_zzz_xyz,
            &buffer_zzz_xzz,
            &buffer_zzz_yyy,
            &buffer_zzz_yyz,
            &buffer_zzz_yzz,
            &buffer_zzz_zzz,
        });

        auto ptr_buffer_0_xxx = buffer_0_xxx.data();
        auto ptr_buffer_0_xxy = buffer_0_xxy.data();
        auto ptr_buffer_0_xxz = buffer_0_xxz.data();
        auto ptr_buffer_0_xyy = buffer_0_xyy.data();
        auto ptr_buffer_0_xyz = buffer_0_xyz.data();
        auto ptr_buffer_0_xzz = buffer_0_xzz.data();
        auto ptr_buffer_0_yyy = buffer_0_yyy.data();
        auto ptr_buffer_0_yyz = buffer_0_yyz.data();
        auto ptr_buffer_0_yzz = buffer_0_yzz.data();
        auto ptr_buffer_0_zzz = buffer_0_zzz.data();

        auto ptr_buffer_x_xxx = buffer_x_xxx.data();
        auto ptr_buffer_x_xxy = buffer_x_xxy.data();
        auto ptr_buffer_x_xxz = buffer_x_xxz.data();
        auto ptr_buffer_x_xyy = buffer_x_xyy.data();
        auto ptr_buffer_x_xyz = buffer_x_xyz.data();
        auto ptr_buffer_x_xzz = buffer_x_xzz.data();
        auto ptr_buffer_x_yyy = buffer_x_yyy.data();
        auto ptr_buffer_x_yyz = buffer_x_yyz.data();
        auto ptr_buffer_x_yzz = buffer_x_yzz.data();
        auto ptr_buffer_x_zzz = buffer_x_zzz.data();

        auto ptr_buffer_y_xxx = buffer_y_xxx.data();
        auto ptr_buffer_y_xxy = buffer_y_xxy.data();
        auto ptr_buffer_y_xxz = buffer_y_xxz.data();
        auto ptr_buffer_y_xyy = buffer_y_xyy.data();
        auto ptr_buffer_y_xyz = buffer_y_xyz.data();
        auto ptr_buffer_y_xzz = buffer_y_xzz.data();
        auto ptr_buffer_y_yyy = buffer_y_yyy.data();
        auto ptr_buffer_y_yyz = buffer_y_yyz.data();
        auto ptr_buffer_y_yzz = buffer_y_yzz.data();
        auto ptr_buffer_y_zzz = buffer_y_zzz.data();

        auto ptr_buffer_z_xxx = buffer_z_xxx.data();
        auto ptr_buffer_z_xxy = buffer_z_xxy.data();
        auto ptr_buffer_z_xxz = buffer_z_xxz.data();
        auto ptr_buffer_z_xyy = buffer_z_xyy.data();
        auto ptr_buffer_z_xyz = buffer_z_xyz.data();
        auto ptr_buffer_z_xzz = buffer_z_xzz.data();
        auto ptr_buffer_z_yyy = buffer_z_yyy.data();
        auto ptr_buffer_z_yyz = buffer_z_yyz.data();
        auto ptr_buffer_z_yzz = buffer_z_yzz.data();
        auto ptr_buffer_z_zzz = buffer_z_zzz.data();

        auto ptr_buffer_xx_xxx = buffer_xx_xxx.data();
        auto ptr_buffer_xx_xxy = buffer_xx_xxy.data();
        auto ptr_buffer_xx_xxz = buffer_xx_xxz.data();
        auto ptr_buffer_xx_xyy = buffer_xx_xyy.data();
        auto ptr_buffer_xx_xyz = buffer_xx_xyz.data();
        auto ptr_buffer_xx_xzz = buffer_xx_xzz.data();
        auto ptr_buffer_xx_yyy = buffer_xx_yyy.data();
        auto ptr_buffer_xx_yyz = buffer_xx_yyz.data();
        auto ptr_buffer_xx_yzz = buffer_xx_yzz.data();
        auto ptr_buffer_xx_zzz = buffer_xx_zzz.data();

        auto ptr_buffer_xy_xxx = buffer_xy_xxx.data();
        auto ptr_buffer_xy_xxy = buffer_xy_xxy.data();
        auto ptr_buffer_xy_xxz = buffer_xy_xxz.data();
        auto ptr_buffer_xy_xyy = buffer_xy_xyy.data();
        auto ptr_buffer_xy_xyz = buffer_xy_xyz.data();
        auto ptr_buffer_xy_xzz = buffer_xy_xzz.data();
        auto ptr_buffer_xy_yyy = buffer_xy_yyy.data();
        auto ptr_buffer_xy_yyz = buffer_xy_yyz.data();
        auto ptr_buffer_xy_yzz = buffer_xy_yzz.data();
        auto ptr_buffer_xy_zzz = buffer_xy_zzz.data();

        auto ptr_buffer_xz_xxx = buffer_xz_xxx.data();
        auto ptr_buffer_xz_xxy = buffer_xz_xxy.data();
        auto ptr_buffer_xz_xxz = buffer_xz_xxz.data();
        auto ptr_buffer_xz_xyy = buffer_xz_xyy.data();
        auto ptr_buffer_xz_xyz = buffer_xz_xyz.data();
        auto ptr_buffer_xz_xzz = buffer_xz_xzz.data();
        auto ptr_buffer_xz_yyy = buffer_xz_yyy.data();
        auto ptr_buffer_xz_yyz = buffer_xz_yyz.data();
        auto ptr_buffer_xz_yzz = buffer_xz_yzz.data();
        auto ptr_buffer_xz_zzz = buffer_xz_zzz.data();

        auto ptr_buffer_yy_xxx = buffer_yy_xxx.data();
        auto ptr_buffer_yy_xxy = buffer_yy_xxy.data();
        auto ptr_buffer_yy_xxz = buffer_yy_xxz.data();
        auto ptr_buffer_yy_xyy = buffer_yy_xyy.data();
        auto ptr_buffer_yy_xyz = buffer_yy_xyz.data();
        auto ptr_buffer_yy_xzz = buffer_yy_xzz.data();
        auto ptr_buffer_yy_yyy = buffer_yy_yyy.data();
        auto ptr_buffer_yy_yyz = buffer_yy_yyz.data();
        auto ptr_buffer_yy_yzz = buffer_yy_yzz.data();
        auto ptr_buffer_yy_zzz = buffer_yy_zzz.data();

        auto ptr_buffer_yz_xxx = buffer_yz_xxx.data();
        auto ptr_buffer_yz_xxy = buffer_yz_xxy.data();
        auto ptr_buffer_yz_xxz = buffer_yz_xxz.data();
        auto ptr_buffer_yz_xyy = buffer_yz_xyy.data();
        auto ptr_buffer_yz_xyz = buffer_yz_xyz.data();
        auto ptr_buffer_yz_xzz = buffer_yz_xzz.data();
        auto ptr_buffer_yz_yyy = buffer_yz_yyy.data();
        auto ptr_buffer_yz_yyz = buffer_yz_yyz.data();
        auto ptr_buffer_yz_yzz = buffer_yz_yzz.data();
        auto ptr_buffer_yz_zzz = buffer_yz_zzz.data();

        auto ptr_buffer_zz_xxx = buffer_zz_xxx.data();
        auto ptr_buffer_zz_xxy = buffer_zz_xxy.data();
        auto ptr_buffer_zz_xxz = buffer_zz_xxz.data();
        auto ptr_buffer_zz_xyy = buffer_zz_xyy.data();
        auto ptr_buffer_zz_xyz = buffer_zz_xyz.data();
        auto ptr_buffer_zz_xzz = buffer_zz_xzz.data();
        auto ptr_buffer_zz_yyy = buffer_zz_yyy.data();
        auto ptr_buffer_zz_yyz = buffer_zz_yyz.data();
        auto ptr_buffer_zz_yzz = buffer_zz_yzz.data();
        auto ptr_buffer_zz_zzz = buffer_zz_zzz.data();

        auto ptr_buffer_xxx_xxx = buffer_xxx_xxx.data();
        auto ptr_buffer_xxx_xxy = buffer_xxx_xxy.data();
        auto ptr_buffer_xxx_xxz = buffer_xxx_xxz.data();
        auto ptr_buffer_xxx_xyy = buffer_xxx_xyy.data();
        auto ptr_buffer_xxx_xyz = buffer_xxx_xyz.data();
        auto ptr_buffer_xxx_xzz = buffer_xxx_xzz.data();
        auto ptr_buffer_xxx_yyy = buffer_xxx_yyy.data();
        auto ptr_buffer_xxx_yyz = buffer_xxx_yyz.data();
        auto ptr_buffer_xxx_yzz = buffer_xxx_yzz.data();
        auto ptr_buffer_xxx_zzz = buffer_xxx_zzz.data();

        auto ptr_buffer_xxy_xxx = buffer_xxy_xxx.data();
        auto ptr_buffer_xxy_xxy = buffer_xxy_xxy.data();
        auto ptr_buffer_xxy_xxz = buffer_xxy_xxz.data();
        auto ptr_buffer_xxy_xyy = buffer_xxy_xyy.data();
        auto ptr_buffer_xxy_xyz = buffer_xxy_xyz.data();
        auto ptr_buffer_xxy_xzz = buffer_xxy_xzz.data();
        auto ptr_buffer_xxy_yyy = buffer_xxy_yyy.data();
        auto ptr_buffer_xxy_yyz = buffer_xxy_yyz.data();
        auto ptr_buffer_xxy_yzz = buffer_xxy_yzz.data();
        auto ptr_buffer_xxy_zzz = buffer_xxy_zzz.data();

        auto ptr_buffer_xxz_xxx = buffer_xxz_xxx.data();
        auto ptr_buffer_xxz_xxy = buffer_xxz_xxy.data();
        auto ptr_buffer_xxz_xxz = buffer_xxz_xxz.data();
        auto ptr_buffer_xxz_xyy = buffer_xxz_xyy.data();
        auto ptr_buffer_xxz_xyz = buffer_xxz_xyz.data();
        auto ptr_buffer_xxz_xzz = buffer_xxz_xzz.data();
        auto ptr_buffer_xxz_yyy = buffer_xxz_yyy.data();
        auto ptr_buffer_xxz_yyz = buffer_xxz_yyz.data();
        auto ptr_buffer_xxz_yzz = buffer_xxz_yzz.data();
        auto ptr_buffer_xxz_zzz = buffer_xxz_zzz.data();

        auto ptr_buffer_xyy_xxx = buffer_xyy_xxx.data();
        auto ptr_buffer_xyy_xxy = buffer_xyy_xxy.data();
        auto ptr_buffer_xyy_xxz = buffer_xyy_xxz.data();
        auto ptr_buffer_xyy_xyy = buffer_xyy_xyy.data();
        auto ptr_buffer_xyy_xyz = buffer_xyy_xyz.data();
        auto ptr_buffer_xyy_xzz = buffer_xyy_xzz.data();
        auto ptr_buffer_xyy_yyy = buffer_xyy_yyy.data();
        auto ptr_buffer_xyy_yyz = buffer_xyy_yyz.data();
        auto ptr_buffer_xyy_yzz = buffer_xyy_yzz.data();
        auto ptr_buffer_xyy_zzz = buffer_xyy_zzz.data();

        auto ptr_buffer_xyz_xxx = buffer_xyz_xxx.data();
        auto ptr_buffer_xyz_xxy = buffer_xyz_xxy.data();
        auto ptr_buffer_xyz_xxz = buffer_xyz_xxz.data();
        auto ptr_buffer_xyz_xyy = buffer_xyz_xyy.data();
        auto ptr_buffer_xyz_xyz = buffer_xyz_xyz.data();
        auto ptr_buffer_xyz_xzz = buffer_xyz_xzz.data();
        auto ptr_buffer_xyz_yyy = buffer_xyz_yyy.data();
        auto ptr_buffer_xyz_yyz = buffer_xyz_yyz.data();
        auto ptr_buffer_xyz_yzz = buffer_xyz_yzz.data();
        auto ptr_buffer_xyz_zzz = buffer_xyz_zzz.data();

        auto ptr_buffer_xzz_xxx = buffer_xzz_xxx.data();
        auto ptr_buffer_xzz_xxy = buffer_xzz_xxy.data();
        auto ptr_buffer_xzz_xxz = buffer_xzz_xxz.data();
        auto ptr_buffer_xzz_xyy = buffer_xzz_xyy.data();
        auto ptr_buffer_xzz_xyz = buffer_xzz_xyz.data();
        auto ptr_buffer_xzz_xzz = buffer_xzz_xzz.data();
        auto ptr_buffer_xzz_yyy = buffer_xzz_yyy.data();
        auto ptr_buffer_xzz_yyz = buffer_xzz_yyz.data();
        auto ptr_buffer_xzz_yzz = buffer_xzz_yzz.data();
        auto ptr_buffer_xzz_zzz = buffer_xzz_zzz.data();

        auto ptr_buffer_yyy_xxx = buffer_yyy_xxx.data();
        auto ptr_buffer_yyy_xxy = buffer_yyy_xxy.data();
        auto ptr_buffer_yyy_xxz = buffer_yyy_xxz.data();
        auto ptr_buffer_yyy_xyy = buffer_yyy_xyy.data();
        auto ptr_buffer_yyy_xyz = buffer_yyy_xyz.data();
        auto ptr_buffer_yyy_xzz = buffer_yyy_xzz.data();
        auto ptr_buffer_yyy_yyy = buffer_yyy_yyy.data();
        auto ptr_buffer_yyy_yyz = buffer_yyy_yyz.data();
        auto ptr_buffer_yyy_yzz = buffer_yyy_yzz.data();
        auto ptr_buffer_yyy_zzz = buffer_yyy_zzz.data();

        auto ptr_buffer_yyz_xxx = buffer_yyz_xxx.data();
        auto ptr_buffer_yyz_xxy = buffer_yyz_xxy.data();
        auto ptr_buffer_yyz_xxz = buffer_yyz_xxz.data();
        auto ptr_buffer_yyz_xyy = buffer_yyz_xyy.data();
        auto ptr_buffer_yyz_xyz = buffer_yyz_xyz.data();
        auto ptr_buffer_yyz_xzz = buffer_yyz_xzz.data();
        auto ptr_buffer_yyz_yyy = buffer_yyz_yyy.data();
        auto ptr_buffer_yyz_yyz = buffer_yyz_yyz.data();
        auto ptr_buffer_yyz_yzz = buffer_yyz_yzz.data();
        auto ptr_buffer_yyz_zzz = buffer_yyz_zzz.data();

        auto ptr_buffer_yzz_xxx = buffer_yzz_xxx.data();
        auto ptr_buffer_yzz_xxy = buffer_yzz_xxy.data();
        auto ptr_buffer_yzz_xxz = buffer_yzz_xxz.data();
        auto ptr_buffer_yzz_xyy = buffer_yzz_xyy.data();
        auto ptr_buffer_yzz_xyz = buffer_yzz_xyz.data();
        auto ptr_buffer_yzz_xzz = buffer_yzz_xzz.data();
        auto ptr_buffer_yzz_yyy = buffer_yzz_yyy.data();
        auto ptr_buffer_yzz_yyz = buffer_yzz_yyz.data();
        auto ptr_buffer_yzz_yzz = buffer_yzz_yzz.data();
        auto ptr_buffer_yzz_zzz = buffer_yzz_zzz.data();

        auto ptr_buffer_zzz_xxx = buffer_zzz_xxx.data();
        auto ptr_buffer_zzz_xxy = buffer_zzz_xxy.data();
        auto ptr_buffer_zzz_xxz = buffer_zzz_xxz.data();
        auto ptr_buffer_zzz_xyy = buffer_zzz_xyy.data();
        auto ptr_buffer_zzz_xyz = buffer_zzz_xyz.data();
        auto ptr_buffer_zzz_xzz = buffer_zzz_xzz.data();
        auto ptr_buffer_zzz_yyy = buffer_zzz_yyy.data();
        auto ptr_buffer_zzz_yyz = buffer_zzz_yyz.data();
        auto ptr_buffer_zzz_yzz = buffer_zzz_yzz.data();
        auto ptr_buffer_zzz_zzz = buffer_zzz_zzz.data();

        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];

                std::ranges::fill(buffer_0_xxx, 0.0);
                std::ranges::fill(buffer_0_xxy, 0.0);
                std::ranges::fill(buffer_0_xxz, 0.0);
                std::ranges::fill(buffer_0_xyy, 0.0);
                std::ranges::fill(buffer_0_xyz, 0.0);
                std::ranges::fill(buffer_0_xzz, 0.0);
                std::ranges::fill(buffer_0_yyy, 0.0);
                std::ranges::fill(buffer_0_yyz, 0.0);
                std::ranges::fill(buffer_0_yzz, 0.0);
                std::ranges::fill(buffer_0_zzz, 0.0);

                std::ranges::fill(buffer_x_xxx, 0.0);
                std::ranges::fill(buffer_x_xxy, 0.0);
                std::ranges::fill(buffer_x_xxz, 0.0);
                std::ranges::fill(buffer_x_xyy, 0.0);
                std::ranges::fill(buffer_x_xyz, 0.0);
                std::ranges::fill(buffer_x_xzz, 0.0);
                std::ranges::fill(buffer_x_yyy, 0.0);
                std::ranges::fill(buffer_x_yyz, 0.0);
                std::ranges::fill(buffer_x_yzz, 0.0);
                std::ranges::fill(buffer_x_zzz, 0.0);

                std::ranges::fill(buffer_y_xxx, 0.0);
                std::ranges::fill(buffer_y_xxy, 0.0);
                std::ranges::fill(buffer_y_xxz, 0.0);
                std::ranges::fill(buffer_y_xyy, 0.0);
                std::ranges::fill(buffer_y_xyz, 0.0);
                std::ranges::fill(buffer_y_xzz, 0.0);
                std::ranges::fill(buffer_y_yyy, 0.0);
                std::ranges::fill(buffer_y_yyz, 0.0);
                std::ranges::fill(buffer_y_yzz, 0.0);
                std::ranges::fill(buffer_y_zzz, 0.0);

                std::ranges::fill(buffer_z_xxx, 0.0);
                std::ranges::fill(buffer_z_xxy, 0.0);
                std::ranges::fill(buffer_z_xxz, 0.0);
                std::ranges::fill(buffer_z_xyy, 0.0);
                std::ranges::fill(buffer_z_xyz, 0.0);
                std::ranges::fill(buffer_z_xzz, 0.0);
                std::ranges::fill(buffer_z_yyy, 0.0);
                std::ranges::fill(buffer_z_yyz, 0.0);
                std::ranges::fill(buffer_z_yzz, 0.0);
                std::ranges::fill(buffer_z_zzz, 0.0);

                std::ranges::fill(buffer_xx_xxx, 0.0);
                std::ranges::fill(buffer_xx_xxy, 0.0);
                std::ranges::fill(buffer_xx_xxz, 0.0);
                std::ranges::fill(buffer_xx_xyy, 0.0);
                std::ranges::fill(buffer_xx_xyz, 0.0);
                std::ranges::fill(buffer_xx_xzz, 0.0);
                std::ranges::fill(buffer_xx_yyy, 0.0);
                std::ranges::fill(buffer_xx_yyz, 0.0);
                std::ranges::fill(buffer_xx_yzz, 0.0);
                std::ranges::fill(buffer_xx_zzz, 0.0);

                std::ranges::fill(buffer_xy_xxx, 0.0);
                std::ranges::fill(buffer_xy_xxy, 0.0);
                std::ranges::fill(buffer_xy_xxz, 0.0);
                std::ranges::fill(buffer_xy_xyy, 0.0);
                std::ranges::fill(buffer_xy_xyz, 0.0);
                std::ranges::fill(buffer_xy_xzz, 0.0);
                std::ranges::fill(buffer_xy_yyy, 0.0);
                std::ranges::fill(buffer_xy_yyz, 0.0);
                std::ranges::fill(buffer_xy_yzz, 0.0);
                std::ranges::fill(buffer_xy_zzz, 0.0);

                std::ranges::fill(buffer_xz_xxx, 0.0);
                std::ranges::fill(buffer_xz_xxy, 0.0);
                std::ranges::fill(buffer_xz_xxz, 0.0);
                std::ranges::fill(buffer_xz_xyy, 0.0);
                std::ranges::fill(buffer_xz_xyz, 0.0);
                std::ranges::fill(buffer_xz_xzz, 0.0);
                std::ranges::fill(buffer_xz_yyy, 0.0);
                std::ranges::fill(buffer_xz_yyz, 0.0);
                std::ranges::fill(buffer_xz_yzz, 0.0);
                std::ranges::fill(buffer_xz_zzz, 0.0);

                std::ranges::fill(buffer_yy_xxx, 0.0);
                std::ranges::fill(buffer_yy_xxy, 0.0);
                std::ranges::fill(buffer_yy_xxz, 0.0);
                std::ranges::fill(buffer_yy_xyy, 0.0);
                std::ranges::fill(buffer_yy_xyz, 0.0);
                std::ranges::fill(buffer_yy_xzz, 0.0);
                std::ranges::fill(buffer_yy_yyy, 0.0);
                std::ranges::fill(buffer_yy_yyz, 0.0);
                std::ranges::fill(buffer_yy_yzz, 0.0);
                std::ranges::fill(buffer_yy_zzz, 0.0);

                std::ranges::fill(buffer_yz_xxx, 0.0);
                std::ranges::fill(buffer_yz_xxy, 0.0);
                std::ranges::fill(buffer_yz_xxz, 0.0);
                std::ranges::fill(buffer_yz_xyy, 0.0);
                std::ranges::fill(buffer_yz_xyz, 0.0);
                std::ranges::fill(buffer_yz_xzz, 0.0);
                std::ranges::fill(buffer_yz_yyy, 0.0);
                std::ranges::fill(buffer_yz_yyz, 0.0);
                std::ranges::fill(buffer_yz_yzz, 0.0);
                std::ranges::fill(buffer_yz_zzz, 0.0);

                std::ranges::fill(buffer_zz_xxx, 0.0);
                std::ranges::fill(buffer_zz_xxy, 0.0);
                std::ranges::fill(buffer_zz_xxz, 0.0);
                std::ranges::fill(buffer_zz_xyy, 0.0);
                std::ranges::fill(buffer_zz_xyz, 0.0);
                std::ranges::fill(buffer_zz_xzz, 0.0);
                std::ranges::fill(buffer_zz_yyy, 0.0);
                std::ranges::fill(buffer_zz_yyz, 0.0);
                std::ranges::fill(buffer_zz_yzz, 0.0);
                std::ranges::fill(buffer_zz_zzz, 0.0);

                std::ranges::fill(buffer_xxx_xxx, 0.0);
                std::ranges::fill(buffer_xxx_xxy, 0.0);
                std::ranges::fill(buffer_xxx_xxz, 0.0);
                std::ranges::fill(buffer_xxx_xyy, 0.0);
                std::ranges::fill(buffer_xxx_xyz, 0.0);
                std::ranges::fill(buffer_xxx_xzz, 0.0);
                std::ranges::fill(buffer_xxx_yyy, 0.0);
                std::ranges::fill(buffer_xxx_yyz, 0.0);
                std::ranges::fill(buffer_xxx_yzz, 0.0);
                std::ranges::fill(buffer_xxx_zzz, 0.0);

                std::ranges::fill(buffer_xxy_xxx, 0.0);
                std::ranges::fill(buffer_xxy_xxy, 0.0);
                std::ranges::fill(buffer_xxy_xxz, 0.0);
                std::ranges::fill(buffer_xxy_xyy, 0.0);
                std::ranges::fill(buffer_xxy_xyz, 0.0);
                std::ranges::fill(buffer_xxy_xzz, 0.0);
                std::ranges::fill(buffer_xxy_yyy, 0.0);
                std::ranges::fill(buffer_xxy_yyz, 0.0);
                std::ranges::fill(buffer_xxy_yzz, 0.0);
                std::ranges::fill(buffer_xxy_zzz, 0.0);

                std::ranges::fill(buffer_xxz_xxx, 0.0);
                std::ranges::fill(buffer_xxz_xxy, 0.0);
                std::ranges::fill(buffer_xxz_xxz, 0.0);
                std::ranges::fill(buffer_xxz_xyy, 0.0);
                std::ranges::fill(buffer_xxz_xyz, 0.0);
                std::ranges::fill(buffer_xxz_xzz, 0.0);
                std::ranges::fill(buffer_xxz_yyy, 0.0);
                std::ranges::fill(buffer_xxz_yyz, 0.0);
                std::ranges::fill(buffer_xxz_yzz, 0.0);
                std::ranges::fill(buffer_xxz_zzz, 0.0);

                std::ranges::fill(buffer_xyy_xxx, 0.0);
                std::ranges::fill(buffer_xyy_xxy, 0.0);
                std::ranges::fill(buffer_xyy_xxz, 0.0);
                std::ranges::fill(buffer_xyy_xyy, 0.0);
                std::ranges::fill(buffer_xyy_xyz, 0.0);
                std::ranges::fill(buffer_xyy_xzz, 0.0);
                std::ranges::fill(buffer_xyy_yyy, 0.0);
                std::ranges::fill(buffer_xyy_yyz, 0.0);
                std::ranges::fill(buffer_xyy_yzz, 0.0);
                std::ranges::fill(buffer_xyy_zzz, 0.0);

                std::ranges::fill(buffer_xyz_xxx, 0.0);
                std::ranges::fill(buffer_xyz_xxy, 0.0);
                std::ranges::fill(buffer_xyz_xxz, 0.0);
                std::ranges::fill(buffer_xyz_xyy, 0.0);
                std::ranges::fill(buffer_xyz_xyz, 0.0);
                std::ranges::fill(buffer_xyz_xzz, 0.0);
                std::ranges::fill(buffer_xyz_yyy, 0.0);
                std::ranges::fill(buffer_xyz_yyz, 0.0);
                std::ranges::fill(buffer_xyz_yzz, 0.0);
                std::ranges::fill(buffer_xyz_zzz, 0.0);

                std::ranges::fill(buffer_xzz_xxx, 0.0);
                std::ranges::fill(buffer_xzz_xxy, 0.0);
                std::ranges::fill(buffer_xzz_xxz, 0.0);
                std::ranges::fill(buffer_xzz_xyy, 0.0);
                std::ranges::fill(buffer_xzz_xyz, 0.0);
                std::ranges::fill(buffer_xzz_xzz, 0.0);
                std::ranges::fill(buffer_xzz_yyy, 0.0);
                std::ranges::fill(buffer_xzz_yyz, 0.0);
                std::ranges::fill(buffer_xzz_yzz, 0.0);
                std::ranges::fill(buffer_xzz_zzz, 0.0);

                std::ranges::fill(buffer_yyy_xxx, 0.0);
                std::ranges::fill(buffer_yyy_xxy, 0.0);
                std::ranges::fill(buffer_yyy_xxz, 0.0);
                std::ranges::fill(buffer_yyy_xyy, 0.0);
                std::ranges::fill(buffer_yyy_xyz, 0.0);
                std::ranges::fill(buffer_yyy_xzz, 0.0);
                std::ranges::fill(buffer_yyy_yyy, 0.0);
                std::ranges::fill(buffer_yyy_yyz, 0.0);
                std::ranges::fill(buffer_yyy_yzz, 0.0);
                std::ranges::fill(buffer_yyy_zzz, 0.0);

                std::ranges::fill(buffer_yyz_xxx, 0.0);
                std::ranges::fill(buffer_yyz_xxy, 0.0);
                std::ranges::fill(buffer_yyz_xxz, 0.0);
                std::ranges::fill(buffer_yyz_xyy, 0.0);
                std::ranges::fill(buffer_yyz_xyz, 0.0);
                std::ranges::fill(buffer_yyz_xzz, 0.0);
                std::ranges::fill(buffer_yyz_yyy, 0.0);
                std::ranges::fill(buffer_yyz_yyz, 0.0);
                std::ranges::fill(buffer_yyz_yzz, 0.0);
                std::ranges::fill(buffer_yyz_zzz, 0.0);

                std::ranges::fill(buffer_yzz_xxx, 0.0);
                std::ranges::fill(buffer_yzz_xxy, 0.0);
                std::ranges::fill(buffer_yzz_xxz, 0.0);
                std::ranges::fill(buffer_yzz_xyy, 0.0);
                std::ranges::fill(buffer_yzz_xyz, 0.0);
                std::ranges::fill(buffer_yzz_xzz, 0.0);
                std::ranges::fill(buffer_yzz_yyy, 0.0);
                std::ranges::fill(buffer_yzz_yyz, 0.0);
                std::ranges::fill(buffer_yzz_yzz, 0.0);
                std::ranges::fill(buffer_yzz_zzz, 0.0);

                std::ranges::fill(buffer_zzz_xxx, 0.0);
                std::ranges::fill(buffer_zzz_xxy, 0.0);
                std::ranges::fill(buffer_zzz_xxz, 0.0);
                std::ranges::fill(buffer_zzz_xyy, 0.0);
                std::ranges::fill(buffer_zzz_xyz, 0.0);
                std::ranges::fill(buffer_zzz_xzz, 0.0);
                std::ranges::fill(buffer_zzz_yyy, 0.0);
                std::ranges::fill(buffer_zzz_yyz, 0.0);
                std::ranges::fill(buffer_zzz_yzz, 0.0);
                std::ranges::fill(buffer_zzz_zzz, 0.0);

                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
                        
                        const auto f00 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                        const auto fg0 = -2.0 * fexp;

                        const auto fg1 = f00 * fg0;

                        const auto fg2 = fg1 * fg0;

                        const auto fg3 = fg2 * fg0;

                        ptr_buffer_0_xxx[k] += f00 * gr_x * gr_x * gr_x;

                        ptr_buffer_x_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_x + 3.0 * f00 * gr_x * gr_x;
                        ptr_buffer_y_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxx[k] += fg1 * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_xx_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_x + 7.0 * fg1 * gr_x * gr_x * gr_x + 6.0 * f00 * gr_x;
                        ptr_buffer_xy_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_y + 3.0 * fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_xz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_z + 3.0 * fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_yy_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x * gr_x;
                        ptr_buffer_yz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xxx[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x * gr_x;

                        ptr_buffer_xxx_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_x * gr_x + 12.0 * fg2 * gr_x * gr_x * gr_x * gr_x + 27.0 * fg1 * gr_x * gr_x + 6.0 * f00;
                        ptr_buffer_xxy_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_x * gr_y + 7.0 * fg2 * gr_x * gr_x * gr_x * gr_y + 6.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xxz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_x * gr_z + 7.0 * fg2 * gr_x * gr_x * gr_x * gr_z + 6.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xyy_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x * gr_x * gr_x + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_x;
                        ptr_buffer_xyz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_xzz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_x + 3.0 * fg2 * gr_x * gr_x * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_x;
                        ptr_buffer_yyy_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_yyz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_x * gr_x * gr_z;
                        ptr_buffer_yzz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_zzz_xxx[k] += fg3 * gr_x * gr_x * gr_x * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xxy[k] += f00 * gr_x * gr_x * gr_y;

                        ptr_buffer_x_xxy[k] += fg1 * gr_x * gr_x * gr_x * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_y_xxy[k] += fg1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_x * gr_x;
                        ptr_buffer_z_xxy[k] += fg1 * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_xx_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_y + 5.0 * fg1 * gr_x * gr_x * gr_y + 2.0 * f00 * gr_y;
                        ptr_buffer_xy_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + fg1 * gr_x * gr_x * gr_x + 2.0 * fg1 * gr_x * gr_y * gr_y + 2.0 * f00 * gr_x;
                        ptr_buffer_xz_xxy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_yz_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_zz_xxy[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y;

                        ptr_buffer_xxx_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_x * gr_y + 9.0 * fg2 * gr_x * gr_x * gr_x * gr_y + 12.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xxy_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x * gr_x * gr_x + 5.0 * fg2 * gr_x * gr_x * gr_y * gr_y + 5.0 * fg1 * gr_x * gr_x + 2.0 * fg1 * gr_y * gr_y + 2.0 * f00;
                        ptr_buffer_xxz_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_z + 5.0 * fg2 * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xyy_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_x * gr_x * gr_y + 2.0 * fg2 * gr_x * gr_y * gr_y * gr_y + 6.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xyz_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_x * gr_x * gr_z + 2.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xzz_xxy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_y + 2.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_yyy_xxy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_y + 6.0 * fg2 * gr_x * gr_x * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_x;
                        ptr_buffer_yyz_xxy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_yzz_xxy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x;
                        ptr_buffer_zzz_xxy[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xxz[k] += f00 * gr_x * gr_x * gr_z;

                        ptr_buffer_x_xxz[k] += fg1 * gr_x * gr_x * gr_x * gr_z + 2.0 * f00 * gr_x * gr_z;
                        ptr_buffer_y_xxz[k] += fg1 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xxz[k] += fg1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_x * gr_x;

                        ptr_buffer_xx_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_x * gr_z + 5.0 * fg1 * gr_x * gr_x * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_xy_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_xxz[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x * gr_x + 2.0 * fg1 * gr_x * gr_z * gr_z + 2.0 * f00 * gr_x;
                        ptr_buffer_yy_xxz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z;
                        ptr_buffer_yz_xxz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y;
                        ptr_buffer_zz_xxz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_x * gr_z;

                        ptr_buffer_xxx_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_x * gr_z + 9.0 * fg2 * gr_x * gr_x * gr_x * gr_z + 12.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xxy_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_z + 5.0 * fg2 * gr_x * gr_x * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xxz_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_x + 5.0 * fg2 * gr_x * gr_x * gr_z * gr_z + 5.0 * fg1 * gr_x * gr_x + 2.0 * fg1 * gr_z * gr_z + 2.0 * f00;
                        ptr_buffer_xyy_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_x * gr_x * gr_z + 2.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xyz_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_y + 2.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xzz_xxz[k] += fg3 * gr_x * gr_x * gr_x * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_x * gr_z + 2.0 * fg2 * gr_x * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yyy_xxz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_yyz_xxz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x;
                        ptr_buffer_yzz_xxz[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_zzz_xxz[k] += fg3 * gr_x * gr_x * gr_z * gr_z * gr_z * gr_z + 6.0 * fg2 * gr_x * gr_x * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_x;

                        ptr_buffer_0_xyy[k] += f00 * gr_x * gr_y * gr_y;

                        ptr_buffer_x_xyy[k] += fg1 * gr_x * gr_x * gr_y * gr_y + f00 * gr_y * gr_y;
                        ptr_buffer_y_xyy[k] += fg1 * gr_x * gr_y * gr_y * gr_y + 2.0 * f00 * gr_x * gr_y;
                        ptr_buffer_z_xyy[k] += fg1 * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_xx_xyy[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_xy_xyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_x * gr_x * gr_y + fg1 * gr_y * gr_y * gr_y + 2.0 * f00 * gr_y;
                        ptr_buffer_xz_xyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_yy_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_y + 5.0 * fg1 * gr_x * gr_y * gr_y + 2.0 * f00 * gr_x;
                        ptr_buffer_yz_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xyy[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y;

                        ptr_buffer_xxx_xyy[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_y + 6.0 * fg2 * gr_x * gr_x * gr_y * gr_y + 3.0 * fg1 * gr_y * gr_y;
                        ptr_buffer_xxy_xyy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_y + 2.0 * fg2 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_y + 6.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xxz_xyy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_xyy_xyy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_y + 5.0 * fg2 * gr_x * gr_x * gr_y * gr_y + 2.0 * fg1 * gr_x * gr_x + fg2 * gr_y * gr_y * gr_y * gr_y + 5.0 * fg1 * gr_y * gr_y + 2.0 * f00;
                        ptr_buffer_xyz_xyy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xzz_xyy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y;
                        ptr_buffer_yyy_xyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_y + 9.0 * fg2 * gr_x * gr_y * gr_y * gr_y + 12.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_yyz_xyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_z + 5.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yzz_xyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y * gr_y * gr_y + 2.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_zzz_xyy[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_0_xyz[k] += f00 * gr_x * gr_y * gr_z;

                        ptr_buffer_x_xyz[k] += fg1 * gr_x * gr_x * gr_y * gr_z + f00 * gr_y * gr_z;
                        ptr_buffer_y_xyz[k] += fg1 * gr_x * gr_y * gr_y * gr_z + f00 * gr_x * gr_z;
                        ptr_buffer_z_xyz[k] += fg1 * gr_x * gr_y * gr_z * gr_z + f00 * gr_x * gr_y;

                        ptr_buffer_xx_xyz[k] += fg2 * gr_x * gr_x * gr_x * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xy_xyz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_x * gr_x * gr_z + fg1 * gr_y * gr_y * gr_z + f00 * gr_z;
                        ptr_buffer_xz_xyz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_x * gr_x * gr_y + fg1 * gr_y * gr_z * gr_z + f00 * gr_y;
                        ptr_buffer_yy_xyz[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yz_xyz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y + fg1 * gr_x * gr_z * gr_z + f00 * gr_x;
                        ptr_buffer_zz_xyz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_y * gr_z;

                        ptr_buffer_xxx_xyz[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_y * gr_z + 6.0 * fg2 * gr_x * gr_x * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xxy_xyz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xxz_xyz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_x * gr_y + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xyy_xyz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xyz_xyz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_x * gr_x * gr_z * gr_z + fg1 * gr_x * gr_x + fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y + fg1 * gr_z * gr_z + f00;
                        ptr_buffer_xzz_xyz[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_yyy_xyz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_z + 6.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yyz_xyz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_yzz_xyz[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_zzz_xyz[k] += fg3 * gr_x * gr_y * gr_z * gr_z * gr_z * gr_z + 6.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_y;

                        ptr_buffer_0_xzz[k] += f00 * gr_x * gr_z * gr_z;

                        ptr_buffer_x_xzz[k] += fg1 * gr_x * gr_x * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_y_xzz[k] += fg1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_z_xzz[k] += fg1 * gr_x * gr_z * gr_z * gr_z + 2.0 * f00 * gr_x * gr_z;

                        ptr_buffer_xx_xzz[k] += fg2 * gr_x * gr_x * gr_x * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_xy_xzz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_xz_xzz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_x * gr_z + fg1 * gr_z * gr_z * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_yy_xzz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_yz_xzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_zz_xzz[k] += fg2 * gr_x * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_x * gr_z * gr_z + 2.0 * f00 * gr_x;

                        ptr_buffer_xxx_xzz[k] += fg3 * gr_x * gr_x * gr_x * gr_x * gr_z * gr_z + 6.0 * fg2 * gr_x * gr_x * gr_z * gr_z + 3.0 * fg1 * gr_z * gr_z;
                        ptr_buffer_xxy_xzz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_xxz_xzz[k] += fg3 * gr_x * gr_x * gr_x * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_x * gr_x * gr_x * gr_z + 3.0 * fg2 * gr_x * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xyy_xzz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_z * gr_z + fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_z * gr_z;
                        ptr_buffer_xyz_xzz[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xzz_xzz[k] += fg3 * gr_x * gr_x * gr_z * gr_z * gr_z * gr_z + 5.0 * fg2 * gr_x * gr_x * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_x + fg2 * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_z * gr_z + 2.0 * f00;
                        ptr_buffer_yyy_xzz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_yyz_xzz[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yzz_xzz[k] += fg3 * gr_x * gr_y * gr_z * gr_z * gr_z * gr_z + 5.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_zzz_xzz[k] += fg3 * gr_x * gr_z * gr_z * gr_z * gr_z * gr_z + 9.0 * fg2 * gr_x * gr_z * gr_z * gr_z + 12.0 * fg1 * gr_x * gr_z;

                        ptr_buffer_0_yyy[k] += f00 * gr_y * gr_y * gr_y;

                        ptr_buffer_x_yyy[k] += fg1 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_y_yyy[k] += fg1 * gr_y * gr_y * gr_y * gr_y + 3.0 * f00 * gr_y * gr_y;
                        ptr_buffer_z_yyy[k] += fg1 * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_xx_yyy[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_y + fg1 * gr_y * gr_y * gr_y;
                        ptr_buffer_xy_yyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_xz_yyy[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z;
                        ptr_buffer_yy_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_y + 7.0 * fg1 * gr_y * gr_y * gr_y + 6.0 * f00 * gr_y;
                        ptr_buffer_yz_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_z + 3.0 * fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_zz_yyy[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y * gr_y;

                        ptr_buffer_xxx_yyy[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_xxy_yyy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_y * gr_y * gr_y * gr_y + 3.0 * fg1 * gr_y * gr_y;
                        ptr_buffer_xxz_yyy[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + fg2 * gr_y * gr_y * gr_y * gr_z;
                        ptr_buffer_xyy_yyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_y + 7.0 * fg2 * gr_x * gr_y * gr_y * gr_y + 6.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xyz_yyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_xzz_yyy[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_yyy_yyy[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_y * gr_y + 12.0 * fg2 * gr_y * gr_y * gr_y * gr_y + 27.0 * fg1 * gr_y * gr_y + 6.0 * f00;
                        ptr_buffer_yyz_yyy[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_y * gr_z + 7.0 * fg2 * gr_y * gr_y * gr_y * gr_z + 6.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_yzz_yyy[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_y * gr_y * gr_y * gr_y + 3.0 * fg2 * gr_y * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_y;
                        ptr_buffer_zzz_yyy[k] += fg3 * gr_y * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yyz[k] += f00 * gr_y * gr_y * gr_z;

                        ptr_buffer_x_yyz[k] += fg1 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_y_yyz[k] += fg1 * gr_y * gr_y * gr_y * gr_z + 2.0 * f00 * gr_y * gr_z;
                        ptr_buffer_z_yyz[k] += fg1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_y * gr_y;

                        ptr_buffer_xx_yyz[k] += fg2 * gr_x * gr_x * gr_y * gr_y * gr_z + fg1 * gr_y * gr_y * gr_z;
                        ptr_buffer_xy_yyz[k] += fg2 * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_xz_yyz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_y * gr_y;
                        ptr_buffer_yy_yyz[k] += fg2 * gr_y * gr_y * gr_y * gr_y * gr_z + 5.0 * fg1 * gr_y * gr_y * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_yz_yyz[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y * gr_y + 2.0 * fg1 * gr_y * gr_z * gr_z + 2.0 * f00 * gr_y;
                        ptr_buffer_zz_yyz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_y * gr_z;

                        ptr_buffer_xxx_yyz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_xxy_yyz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xxz_yyz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_y * gr_y + fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_y * gr_y;
                        ptr_buffer_xyy_yyz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_y * gr_z + 5.0 * fg2 * gr_x * gr_y * gr_y * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xyz_yyz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_y * gr_y * gr_y + 2.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_xzz_yyz[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_yyy_yyz[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_y * gr_z + 9.0 * fg2 * gr_y * gr_y * gr_y * gr_z + 12.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_yyz_yyz[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_z * gr_z + fg2 * gr_y * gr_y * gr_y * gr_y + 5.0 * fg2 * gr_y * gr_y * gr_z * gr_z + 5.0 * fg1 * gr_y * gr_y + 2.0 * fg1 * gr_z * gr_z + 2.0 * f00;
                        ptr_buffer_yzz_yyz[k] += fg3 * gr_y * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_y * gr_y * gr_z + 2.0 * fg2 * gr_y * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_zzz_yyz[k] += fg3 * gr_y * gr_y * gr_z * gr_z * gr_z * gr_z + 6.0 * fg2 * gr_y * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_y;

                        ptr_buffer_0_yzz[k] += f00 * gr_y * gr_z * gr_z;

                        ptr_buffer_x_yzz[k] += fg1 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_y_yzz[k] += fg1 * gr_y * gr_y * gr_z * gr_z + f00 * gr_z * gr_z;
                        ptr_buffer_z_yzz[k] += fg1 * gr_y * gr_z * gr_z * gr_z + 2.0 * f00 * gr_y * gr_z;

                        ptr_buffer_xx_yzz[k] += fg2 * gr_x * gr_x * gr_y * gr_z * gr_z + fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_xy_yzz[k] += fg2 * gr_x * gr_y * gr_y * gr_z * gr_z + fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_xz_yzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y * gr_z;
                        ptr_buffer_yy_yzz[k] += fg2 * gr_y * gr_y * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_yz_yzz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_y * gr_z + fg1 * gr_z * gr_z * gr_z + 2.0 * f00 * gr_z;
                        ptr_buffer_zz_yzz[k] += fg2 * gr_y * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_y * gr_z * gr_z + 2.0 * f00 * gr_y;

                        ptr_buffer_xxx_yzz[k] += fg3 * gr_x * gr_x * gr_x * gr_y * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_xxy_yzz[k] += fg3 * gr_x * gr_x * gr_y * gr_y * gr_z * gr_z + fg2 * gr_x * gr_x * gr_z * gr_z + fg2 * gr_y * gr_y * gr_z * gr_z + fg1 * gr_z * gr_z;
                        ptr_buffer_xxz_yzz[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_x * gr_x * gr_y * gr_z + fg2 * gr_y * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_xyy_yzz[k] += fg3 * gr_x * gr_y * gr_y * gr_y * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_xyz_yzz[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_x * gr_y * gr_y * gr_z + fg2 * gr_x * gr_z * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_xzz_yzz[k] += fg3 * gr_x * gr_y * gr_z * gr_z * gr_z * gr_z + 5.0 * fg2 * gr_x * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_x * gr_y;
                        ptr_buffer_yyy_yzz[k] += fg3 * gr_y * gr_y * gr_y * gr_y * gr_z * gr_z + 6.0 * fg2 * gr_y * gr_y * gr_z * gr_z + 3.0 * fg1 * gr_z * gr_z;
                        ptr_buffer_yyz_yzz[k] += fg3 * gr_y * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * fg2 * gr_y * gr_y * gr_y * gr_z + 3.0 * fg2 * gr_y * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_yzz_yzz[k] += fg3 * gr_y * gr_y * gr_z * gr_z * gr_z * gr_z + 5.0 * fg2 * gr_y * gr_y * gr_z * gr_z + 2.0 * fg1 * gr_y * gr_y + fg2 * gr_z * gr_z * gr_z * gr_z + 5.0 * fg1 * gr_z * gr_z + 2.0 * f00;
                        ptr_buffer_zzz_yzz[k] += fg3 * gr_y * gr_z * gr_z * gr_z * gr_z * gr_z + 9.0 * fg2 * gr_y * gr_z * gr_z * gr_z + 12.0 * fg1 * gr_y * gr_z;

                        ptr_buffer_0_zzz[k] += f00 * gr_z * gr_z * gr_z;

                        ptr_buffer_x_zzz[k] += fg1 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_y_zzz[k] += fg1 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_z_zzz[k] += fg1 * gr_z * gr_z * gr_z * gr_z + 3.0 * f00 * gr_z * gr_z;

                        ptr_buffer_xx_zzz[k] += fg2 * gr_x * gr_x * gr_z * gr_z * gr_z + fg1 * gr_z * gr_z * gr_z;
                        ptr_buffer_xy_zzz[k] += fg2 * gr_x * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_xz_zzz[k] += fg2 * gr_x * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_x * gr_z * gr_z;
                        ptr_buffer_yy_zzz[k] += fg2 * gr_y * gr_y * gr_z * gr_z * gr_z + fg1 * gr_z * gr_z * gr_z;
                        ptr_buffer_yz_zzz[k] += fg2 * gr_y * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_y * gr_z * gr_z;
                        ptr_buffer_zz_zzz[k] += fg2 * gr_z * gr_z * gr_z * gr_z * gr_z + 7.0 * fg1 * gr_z * gr_z * gr_z + 6.0 * f00 * gr_z;

                        ptr_buffer_xxx_zzz[k] += fg3 * gr_x * gr_x * gr_x * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_xxy_zzz[k] += fg3 * gr_x * gr_x * gr_y * gr_z * gr_z * gr_z + fg2 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_xxz_zzz[k] += fg3 * gr_x * gr_x * gr_z * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_x * gr_z * gr_z + fg2 * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_z * gr_z;
                        ptr_buffer_xyy_zzz[k] += fg3 * gr_x * gr_y * gr_y * gr_z * gr_z * gr_z + fg2 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_xyz_zzz[k] += fg3 * gr_x * gr_y * gr_z * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_xzz_zzz[k] += fg3 * gr_x * gr_z * gr_z * gr_z * gr_z * gr_z + 7.0 * fg2 * gr_x * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_x * gr_z;
                        ptr_buffer_yyy_zzz[k] += fg3 * gr_y * gr_y * gr_y * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_yyz_zzz[k] += fg3 * gr_y * gr_y * gr_z * gr_z * gr_z * gr_z + 3.0 * fg2 * gr_y * gr_y * gr_z * gr_z + fg2 * gr_z * gr_z * gr_z * gr_z + 3.0 * fg1 * gr_z * gr_z;
                        ptr_buffer_yzz_zzz[k] += fg3 * gr_y * gr_z * gr_z * gr_z * gr_z * gr_z + 7.0 * fg2 * gr_y * gr_z * gr_z * gr_z + 6.0 * fg1 * gr_y * gr_z;
                        ptr_buffer_zzz_zzz[k] += fg3 * gr_z * gr_z * gr_z * gr_z * gr_z * gr_z + 12.0 * fg2 * gr_z * gr_z * gr_z * gr_z + 27.0 * fg1 * gr_z * gr_z + 6.0 * f00;
                    }
                }
        
                // distribute GTO values into submatrix
        
                for (int isph = 0; isph < nsph; isph++)
                {
                    for (int icomp = 0; icomp < static_cast<int>(sph_cart_coefs[isph].size()); icomp++)
                    {
                        auto icart = sph_cart_coefs[isph][icomp].first;
                        auto fcart = sph_cart_coefs[isph][icomp].second;

                        gtoval::distribute(submat_0, *(buffer_0_refs[icart]), fcart, isph * nrows + irow);

                        gtoval::distribute(submat_x, *(buffer_x_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_y, *(buffer_y_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_z, *(buffer_z_refs[icart]), fcart, isph * nrows + irow);

                        gtoval::distribute(submat_xx, *(buffer_xx_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xy, *(buffer_xy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xz, *(buffer_xz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yy, *(buffer_yy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yz, *(buffer_yz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_zz, *(buffer_zz_refs[icart]), fcart, isph * nrows + irow);

                        gtoval::distribute(submat_xxx, *(buffer_xxx_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xxy, *(buffer_xxy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xxz, *(buffer_xxz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xyy, *(buffer_xyy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xyz, *(buffer_xyz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_xzz, *(buffer_xzz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yyy, *(buffer_yyy_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yyz, *(buffer_yyz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_yzz, *(buffer_yzz_refs[icart]), fcart, isph * nrows + irow);
                        gtoval::distribute(submat_zzz, *(buffer_zzz_refs[icart]), fcart, isph * nrows + irow);
                    }
                }

                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

}  // namespace gtoval

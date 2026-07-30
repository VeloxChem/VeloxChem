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

#include "GtoValuesRecG.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "SphericalMomentum.hpp"

#define ANGULAR_MOMENTUM_G 4

namespace gtoval {  // gtoval namespace

auto
get_lda_values_rec_g(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>& gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_G * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_G>(isph);
    }

    // set up GTO values storage

    const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1));

    const size_t ncols = grid_coords_x.size();

    auto gto_values = matfunc::make_matrix("LDA", 9 * nrows, ncols);

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

    // compute GTO values for G type GTOs on grid

    std::vector<double> buffer_xxxx(ncols);

    std::vector<double> buffer_xxxy(ncols);

    std::vector<double> buffer_xxxz(ncols);

    std::vector<double> buffer_xxyy(ncols);

    std::vector<double> buffer_xxyz(ncols);

    std::vector<double> buffer_xxzz(ncols);

    std::vector<double> buffer_xyyy(ncols);

    std::vector<double> buffer_xyyz(ncols);

    std::vector<double> buffer_xyzz(ncols);

    std::vector<double> buffer_xzzz(ncols);

    std::vector<double> buffer_yyyy(ncols);

    std::vector<double> buffer_yyyz(ncols);

    std::vector<double> buffer_yyzz(ncols);

    std::vector<double> buffer_yzzz(ncols);

    std::vector<double> buffer_zzzz(ncols);

    std::vector<std::vector<double>*> buffer_refs({
        &buffer_xxxx, &buffer_xxxy, &buffer_xxxz, &buffer_xxyy, &buffer_xxyz, &buffer_xxzz,
        &buffer_xyyy, &buffer_xyyz, &buffer_xyzz, &buffer_xzzz, &buffer_yyyy, &buffer_yyyz,
        &buffer_yyzz, &buffer_yzzz, &buffer_zzzz});

    auto ptr_buffer_xxxx = buffer_xxxx.data();

    auto ptr_buffer_xxxy = buffer_xxxy.data();

    auto ptr_buffer_xxxz = buffer_xxxz.data();

    auto ptr_buffer_xxyy = buffer_xxyy.data();

    auto ptr_buffer_xxyz = buffer_xxyz.data();

    auto ptr_buffer_xxzz = buffer_xxzz.data();

    auto ptr_buffer_xyyy = buffer_xyyy.data();

    auto ptr_buffer_xyyz = buffer_xyyz.data();

    auto ptr_buffer_xyzz = buffer_xyzz.data();

    auto ptr_buffer_xzzz = buffer_xzzz.data();

    auto ptr_buffer_yyyy = buffer_yyyy.data();

    auto ptr_buffer_yyyz = buffer_yyyz.data();

    auto ptr_buffer_yyzz = buffer_yyzz.data();

    auto ptr_buffer_yzzz = buffer_yzzz.data();

    auto ptr_buffer_zzzz = buffer_zzzz.data();

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

            std::ranges::fill(buffer_xxxx, 0.0);

            std::ranges::fill(buffer_xxxy, 0.0);

            std::ranges::fill(buffer_xxxz, 0.0);

            std::ranges::fill(buffer_xxyy, 0.0);

            std::ranges::fill(buffer_xxyz, 0.0);

            std::ranges::fill(buffer_xxzz, 0.0);

            std::ranges::fill(buffer_xyyy, 0.0);

            std::ranges::fill(buffer_xyyz, 0.0);

            std::ranges::fill(buffer_xyzz, 0.0);

            std::ranges::fill(buffer_xzzz, 0.0);

            std::ranges::fill(buffer_yyyy, 0.0);

            std::ranges::fill(buffer_yyyz, 0.0);

            std::ranges::fill(buffer_yyzz, 0.0);

            std::ranges::fill(buffer_yzzz, 0.0);

            std::ranges::fill(buffer_zzzz, 0.0);

            for (int j = 0; j < npgtos; j++)
            {
                const auto fexp = gto_exps[j * ncgtos + i];

                const auto fnorm = gto_norms[j * ncgtos + i];

#pragma omp simd
                for (int k = 0; k < static_cast<int>(ncols); k++)
                {
                    const auto gr_x = g_x[k] - r_x;

                    const auto gr_y = g_y[k] - r_y;

                    const auto gr_z = g_z[k] - r_z;

                    const auto fss = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));

                    ptr_buffer_xxxx[k] +=  fss * gr_x * gr_x * gr_x * gr_x;

                    ptr_buffer_xxxy[k] +=  fss * gr_x * gr_x * gr_x * gr_y;

                    ptr_buffer_xxxz[k] +=  fss * gr_x * gr_x * gr_x * gr_z;

                    ptr_buffer_xxyy[k] +=  fss * gr_x * gr_x * gr_y * gr_y;

                    ptr_buffer_xxyz[k] +=  fss * gr_x * gr_x * gr_y * gr_z;

                    ptr_buffer_xxzz[k] +=  fss * gr_x * gr_x * gr_z * gr_z;

                    ptr_buffer_xyyy[k] +=  fss * gr_x * gr_y * gr_y * gr_y;

                    ptr_buffer_xyyz[k] +=  fss * gr_x * gr_y * gr_y * gr_z;

                    ptr_buffer_xyzz[k] +=  fss * gr_x * gr_y * gr_z * gr_z;

                    ptr_buffer_xzzz[k] +=  fss * gr_x * gr_z * gr_z * gr_z;

                    ptr_buffer_yyyy[k] +=  fss * gr_y * gr_y * gr_y * gr_y;

                    ptr_buffer_yyyz[k] +=  fss * gr_y * gr_y * gr_y * gr_z;

                    ptr_buffer_yyzz[k] +=  fss * gr_y * gr_y * gr_z * gr_z;

                    ptr_buffer_yzzz[k] +=  fss * gr_y * gr_z * gr_z * gr_z;

                    ptr_buffer_zzzz[k] +=  fss * gr_z * gr_z * gr_z * gr_z;
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
get_gga_values_rec_g(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // spherical-Cartesian transformation factors

    auto nsph = ANGULAR_MOMENTUM_G * 2 + 1;

    std::vector<std::vector<std::pair<int, double>>> sph_cart_coefs(nsph);

    for (int isph = 0; isph < nsph; isph++)
    {
        sph_cart_coefs[isph] = spher_mom::transformation_factors<ANGULAR_MOMENTUM_G>(isph);
    }

    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", 9 * nrows, ncols);
        
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
        
        // compute GTO values for G type GTOs on grid

        std::vector<double> buffer_0_xxxx(ncols);
        std::vector<double> buffer_0_xxxy(ncols);
        std::vector<double> buffer_0_xxxz(ncols);
        std::vector<double> buffer_0_xxyy(ncols);
        std::vector<double> buffer_0_xxyz(ncols);
        std::vector<double> buffer_0_xxzz(ncols);
        std::vector<double> buffer_0_xyyy(ncols);
        std::vector<double> buffer_0_xyyz(ncols);
        std::vector<double> buffer_0_xyzz(ncols);
        std::vector<double> buffer_0_xzzz(ncols);
        std::vector<double> buffer_0_yyyy(ncols);
        std::vector<double> buffer_0_yyyz(ncols);
        std::vector<double> buffer_0_yyzz(ncols);
        std::vector<double> buffer_0_yzzz(ncols);
        std::vector<double> buffer_0_zzzz(ncols);

        std::vector<double> buffer_x_xxxx(ncols);
        std::vector<double> buffer_x_xxxy(ncols);
        std::vector<double> buffer_x_xxxz(ncols);
        std::vector<double> buffer_x_xxyy(ncols);
        std::vector<double> buffer_x_xxyz(ncols);
        std::vector<double> buffer_x_xxzz(ncols);
        std::vector<double> buffer_x_xyyy(ncols);
        std::vector<double> buffer_x_xyyz(ncols);
        std::vector<double> buffer_x_xyzz(ncols);
        std::vector<double> buffer_x_xzzz(ncols);
        std::vector<double> buffer_x_yyyy(ncols);
        std::vector<double> buffer_x_yyyz(ncols);
        std::vector<double> buffer_x_yyzz(ncols);
        std::vector<double> buffer_x_yzzz(ncols);
        std::vector<double> buffer_x_zzzz(ncols);

        std::vector<double> buffer_y_xxxx(ncols);
        std::vector<double> buffer_y_xxxy(ncols);
        std::vector<double> buffer_y_xxxz(ncols);
        std::vector<double> buffer_y_xxyy(ncols);
        std::vector<double> buffer_y_xxyz(ncols);
        std::vector<double> buffer_y_xxzz(ncols);
        std::vector<double> buffer_y_xyyy(ncols);
        std::vector<double> buffer_y_xyyz(ncols);
        std::vector<double> buffer_y_xyzz(ncols);
        std::vector<double> buffer_y_xzzz(ncols);
        std::vector<double> buffer_y_yyyy(ncols);
        std::vector<double> buffer_y_yyyz(ncols);
        std::vector<double> buffer_y_yyzz(ncols);
        std::vector<double> buffer_y_yzzz(ncols);
        std::vector<double> buffer_y_zzzz(ncols);

        std::vector<double> buffer_z_xxxx(ncols);
        std::vector<double> buffer_z_xxxy(ncols);
        std::vector<double> buffer_z_xxxz(ncols);
        std::vector<double> buffer_z_xxyy(ncols);
        std::vector<double> buffer_z_xxyz(ncols);
        std::vector<double> buffer_z_xxzz(ncols);
        std::vector<double> buffer_z_xyyy(ncols);
        std::vector<double> buffer_z_xyyz(ncols);
        std::vector<double> buffer_z_xyzz(ncols);
        std::vector<double> buffer_z_xzzz(ncols);
        std::vector<double> buffer_z_yyyy(ncols);
        std::vector<double> buffer_z_yyyz(ncols);
        std::vector<double> buffer_z_yyzz(ncols);
        std::vector<double> buffer_z_yzzz(ncols);
        std::vector<double> buffer_z_zzzz(ncols);

        std::vector<std::vector<double>*> buffer_0_refs({
            &buffer_0_xxxx, &buffer_0_xxxy, &buffer_0_xxxz, &buffer_0_xxyy, &buffer_0_xxyz, &buffer_0_xxzz,
            &buffer_0_xyyy, &buffer_0_xyyz, &buffer_0_xyzz, &buffer_0_xzzz, &buffer_0_yyyy, &buffer_0_yyyz,
            &buffer_0_yyzz, &buffer_0_yzzz, &buffer_0_zzzz});

        std::vector<std::vector<double>*> buffer_x_refs({
            &buffer_x_xxxx, &buffer_x_xxxy, &buffer_x_xxxz, &buffer_x_xxyy, &buffer_x_xxyz, &buffer_x_xxzz,
            &buffer_x_xyyy, &buffer_x_xyyz, &buffer_x_xyzz, &buffer_x_xzzz, &buffer_x_yyyy, &buffer_x_yyyz,
            &buffer_x_yyzz, &buffer_x_yzzz, &buffer_x_zzzz});
        std::vector<std::vector<double>*> buffer_y_refs({
            &buffer_y_xxxx, &buffer_y_xxxy, &buffer_y_xxxz, &buffer_y_xxyy, &buffer_y_xxyz, &buffer_y_xxzz,
            &buffer_y_xyyy, &buffer_y_xyyz, &buffer_y_xyzz, &buffer_y_xzzz, &buffer_y_yyyy, &buffer_y_yyyz,
            &buffer_y_yyzz, &buffer_y_yzzz, &buffer_y_zzzz});
        std::vector<std::vector<double>*> buffer_z_refs({
            &buffer_z_xxxx, &buffer_z_xxxy, &buffer_z_xxxz, &buffer_z_xxyy, &buffer_z_xxyz, &buffer_z_xxzz,
            &buffer_z_xyyy, &buffer_z_xyyz, &buffer_z_xyzz, &buffer_z_xzzz, &buffer_z_yyyy, &buffer_z_yyyz,
            &buffer_z_yyzz, &buffer_z_yzzz, &buffer_z_zzzz});

        auto ptr_buffer_0_xxxx = buffer_0_xxxx.data();
        auto ptr_buffer_0_xxxy = buffer_0_xxxy.data();
        auto ptr_buffer_0_xxxz = buffer_0_xxxz.data();
        auto ptr_buffer_0_xxyy = buffer_0_xxyy.data();
        auto ptr_buffer_0_xxyz = buffer_0_xxyz.data();
        auto ptr_buffer_0_xxzz = buffer_0_xxzz.data();
        auto ptr_buffer_0_xyyy = buffer_0_xyyy.data();
        auto ptr_buffer_0_xyyz = buffer_0_xyyz.data();
        auto ptr_buffer_0_xyzz = buffer_0_xyzz.data();
        auto ptr_buffer_0_xzzz = buffer_0_xzzz.data();
        auto ptr_buffer_0_yyyy = buffer_0_yyyy.data();
        auto ptr_buffer_0_yyyz = buffer_0_yyyz.data();
        auto ptr_buffer_0_yyzz = buffer_0_yyzz.data();
        auto ptr_buffer_0_yzzz = buffer_0_yzzz.data();
        auto ptr_buffer_0_zzzz = buffer_0_zzzz.data();

        auto ptr_buffer_x_xxxx = buffer_x_xxxx.data();
        auto ptr_buffer_x_xxxy = buffer_x_xxxy.data();
        auto ptr_buffer_x_xxxz = buffer_x_xxxz.data();
        auto ptr_buffer_x_xxyy = buffer_x_xxyy.data();
        auto ptr_buffer_x_xxyz = buffer_x_xxyz.data();
        auto ptr_buffer_x_xxzz = buffer_x_xxzz.data();
        auto ptr_buffer_x_xyyy = buffer_x_xyyy.data();
        auto ptr_buffer_x_xyyz = buffer_x_xyyz.data();
        auto ptr_buffer_x_xyzz = buffer_x_xyzz.data();
        auto ptr_buffer_x_xzzz = buffer_x_xzzz.data();
        auto ptr_buffer_x_yyyy = buffer_x_yyyy.data();
        auto ptr_buffer_x_yyyz = buffer_x_yyyz.data();
        auto ptr_buffer_x_yyzz = buffer_x_yyzz.data();
        auto ptr_buffer_x_yzzz = buffer_x_yzzz.data();
        auto ptr_buffer_x_zzzz = buffer_x_zzzz.data();

        auto ptr_buffer_y_xxxx = buffer_y_xxxx.data();
        auto ptr_buffer_y_xxxy = buffer_y_xxxy.data();
        auto ptr_buffer_y_xxxz = buffer_y_xxxz.data();
        auto ptr_buffer_y_xxyy = buffer_y_xxyy.data();
        auto ptr_buffer_y_xxyz = buffer_y_xxyz.data();
        auto ptr_buffer_y_xxzz = buffer_y_xxzz.data();
        auto ptr_buffer_y_xyyy = buffer_y_xyyy.data();
        auto ptr_buffer_y_xyyz = buffer_y_xyyz.data();
        auto ptr_buffer_y_xyzz = buffer_y_xyzz.data();
        auto ptr_buffer_y_xzzz = buffer_y_xzzz.data();
        auto ptr_buffer_y_yyyy = buffer_y_yyyy.data();
        auto ptr_buffer_y_yyyz = buffer_y_yyyz.data();
        auto ptr_buffer_y_yyzz = buffer_y_yyzz.data();
        auto ptr_buffer_y_yzzz = buffer_y_yzzz.data();
        auto ptr_buffer_y_zzzz = buffer_y_zzzz.data();

        auto ptr_buffer_z_xxxx = buffer_z_xxxx.data();
        auto ptr_buffer_z_xxxy = buffer_z_xxxy.data();
        auto ptr_buffer_z_xxxz = buffer_z_xxxz.data();
        auto ptr_buffer_z_xxyy = buffer_z_xxyy.data();
        auto ptr_buffer_z_xxyz = buffer_z_xxyz.data();
        auto ptr_buffer_z_xxzz = buffer_z_xxzz.data();
        auto ptr_buffer_z_xyyy = buffer_z_xyyy.data();
        auto ptr_buffer_z_xyyz = buffer_z_xyyz.data();
        auto ptr_buffer_z_xyzz = buffer_z_xyzz.data();
        auto ptr_buffer_z_xzzz = buffer_z_xzzz.data();
        auto ptr_buffer_z_yyyy = buffer_z_yyyy.data();
        auto ptr_buffer_z_yyyz = buffer_z_yyyz.data();
        auto ptr_buffer_z_yyzz = buffer_z_yyzz.data();
        auto ptr_buffer_z_yzzz = buffer_z_yzzz.data();
        auto ptr_buffer_z_zzzz = buffer_z_zzzz.data();

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

                std::ranges::fill(buffer_0_xxxx, 0.0);
                std::ranges::fill(buffer_0_xxxy, 0.0);
                std::ranges::fill(buffer_0_xxxz, 0.0);
                std::ranges::fill(buffer_0_xxyy, 0.0);
                std::ranges::fill(buffer_0_xxyz, 0.0);
                std::ranges::fill(buffer_0_xxzz, 0.0);
                std::ranges::fill(buffer_0_xyyy, 0.0);
                std::ranges::fill(buffer_0_xyyz, 0.0);
                std::ranges::fill(buffer_0_xyzz, 0.0);
                std::ranges::fill(buffer_0_xzzz, 0.0);
                std::ranges::fill(buffer_0_yyyy, 0.0);
                std::ranges::fill(buffer_0_yyyz, 0.0);
                std::ranges::fill(buffer_0_yyzz, 0.0);
                std::ranges::fill(buffer_0_yzzz, 0.0);
                std::ranges::fill(buffer_0_zzzz, 0.0);

                std::ranges::fill(buffer_x_xxxx, 0.0);
                std::ranges::fill(buffer_x_xxxy, 0.0);
                std::ranges::fill(buffer_x_xxxz, 0.0);
                std::ranges::fill(buffer_x_xxyy, 0.0);
                std::ranges::fill(buffer_x_xxyz, 0.0);
                std::ranges::fill(buffer_x_xxzz, 0.0);
                std::ranges::fill(buffer_x_xyyy, 0.0);
                std::ranges::fill(buffer_x_xyyz, 0.0);
                std::ranges::fill(buffer_x_xyzz, 0.0);
                std::ranges::fill(buffer_x_xzzz, 0.0);
                std::ranges::fill(buffer_x_yyyy, 0.0);
                std::ranges::fill(buffer_x_yyyz, 0.0);
                std::ranges::fill(buffer_x_yyzz, 0.0);
                std::ranges::fill(buffer_x_yzzz, 0.0);
                std::ranges::fill(buffer_x_zzzz, 0.0);

                std::ranges::fill(buffer_y_xxxx, 0.0);
                std::ranges::fill(buffer_y_xxxy, 0.0);
                std::ranges::fill(buffer_y_xxxz, 0.0);
                std::ranges::fill(buffer_y_xxyy, 0.0);
                std::ranges::fill(buffer_y_xxyz, 0.0);
                std::ranges::fill(buffer_y_xxzz, 0.0);
                std::ranges::fill(buffer_y_xyyy, 0.0);
                std::ranges::fill(buffer_y_xyyz, 0.0);
                std::ranges::fill(buffer_y_xyzz, 0.0);
                std::ranges::fill(buffer_y_xzzz, 0.0);
                std::ranges::fill(buffer_y_yyyy, 0.0);
                std::ranges::fill(buffer_y_yyyz, 0.0);
                std::ranges::fill(buffer_y_yyzz, 0.0);
                std::ranges::fill(buffer_y_yzzz, 0.0);
                std::ranges::fill(buffer_y_zzzz, 0.0);

                std::ranges::fill(buffer_z_xxxx, 0.0);
                std::ranges::fill(buffer_z_xxxy, 0.0);
                std::ranges::fill(buffer_z_xxxz, 0.0);
                std::ranges::fill(buffer_z_xxyy, 0.0);
                std::ranges::fill(buffer_z_xxyz, 0.0);
                std::ranges::fill(buffer_z_xxzz, 0.0);
                std::ranges::fill(buffer_z_xyyy, 0.0);
                std::ranges::fill(buffer_z_xyyz, 0.0);
                std::ranges::fill(buffer_z_xyzz, 0.0);
                std::ranges::fill(buffer_z_xzzz, 0.0);
                std::ranges::fill(buffer_z_yyyy, 0.0);
                std::ranges::fill(buffer_z_yyyz, 0.0);
                std::ranges::fill(buffer_z_yyzz, 0.0);
                std::ranges::fill(buffer_z_yzzz, 0.0);
                std::ranges::fill(buffer_z_zzzz, 0.0);

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

                        ptr_buffer_0_xxxx[k] += f00 * gr_x * gr_x * gr_x * gr_x;
                        ptr_buffer_x_xxxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_x * gr_x + 4.0 * f00 * gr_x * gr_x * gr_x;
                        ptr_buffer_y_xxxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxxx[k] += fg_1 * gr_x * gr_x * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xxxy[k] += f00 * gr_x * gr_x * gr_x * gr_y;
                        ptr_buffer_x_xxxy[k] += fg_1 * gr_x * gr_x * gr_x * gr_x * gr_y + 3.0 * f00 * gr_x * gr_x * gr_y;
                        ptr_buffer_y_xxxy[k] += fg_1 * gr_x * gr_x * gr_x * gr_y * gr_y + f00 * gr_x * gr_x * gr_x;
                        ptr_buffer_z_xxxy[k] += fg_1 * gr_x * gr_x * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xxxz[k] += f00 * gr_x * gr_x * gr_x * gr_z;
                        ptr_buffer_x_xxxz[k] += fg_1 * gr_x * gr_x * gr_x * gr_x * gr_z + 3.0 * f00 * gr_x * gr_x * gr_z;
                        ptr_buffer_y_xxxz[k] += fg_1 * gr_x * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xxxz[k] += fg_1 * gr_x * gr_x * gr_x * gr_z * gr_z + f00 * gr_x * gr_x * gr_x;

                        ptr_buffer_0_xxyy[k] += f00 * gr_x * gr_x * gr_y * gr_y;
                        ptr_buffer_x_xxyy[k] += fg_1 * gr_x * gr_x * gr_x * gr_y * gr_y + 2.0 * f00 * gr_x * gr_y * gr_y;
                        ptr_buffer_y_xxyy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y * gr_y + 2.0 * f00 * gr_x * gr_x * gr_y;
                        ptr_buffer_z_xxyy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y * gr_z;

                        ptr_buffer_0_xxyz[k] += f00 * gr_x * gr_x * gr_y * gr_z;
                        ptr_buffer_x_xxyz[k] += fg_1 * gr_x * gr_x * gr_x * gr_y * gr_z + 2.0 * f00 * gr_x * gr_y * gr_z;
                        ptr_buffer_y_xxyz[k] += fg_1 * gr_x * gr_x * gr_y * gr_y * gr_z + f00 * gr_x * gr_x * gr_z;
                        ptr_buffer_z_xxyz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z * gr_z + f00 * gr_x * gr_x * gr_y;

                        ptr_buffer_0_xxzz[k] += f00 * gr_x * gr_x * gr_z * gr_z;
                        ptr_buffer_x_xxzz[k] += fg_1 * gr_x * gr_x * gr_x * gr_z * gr_z + 2.0 * f00 * gr_x * gr_z * gr_z;
                        ptr_buffer_y_xxzz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_z_xxzz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z * gr_z + 2.0 * f00 * gr_x * gr_x * gr_z;

                        ptr_buffer_0_xyyy[k] += f00 * gr_x * gr_y * gr_y * gr_y;
                        ptr_buffer_x_xyyy[k] += fg_1 * gr_x * gr_x * gr_y * gr_y * gr_y + f00 * gr_y * gr_y * gr_y;
                        ptr_buffer_y_xyyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y * gr_y + 3.0 * f00 * gr_x * gr_y * gr_y;
                        ptr_buffer_z_xyyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_0_xyyz[k] += f00 * gr_x * gr_y * gr_y * gr_z;
                        ptr_buffer_x_xyyz[k] += fg_1 * gr_x * gr_x * gr_y * gr_y * gr_z + f00 * gr_y * gr_y * gr_z;
                        ptr_buffer_y_xyyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_y * gr_z + 2.0 * f00 * gr_x * gr_y * gr_z;
                        ptr_buffer_z_xyyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z * gr_z + f00 * gr_x * gr_y * gr_y;

                        ptr_buffer_0_xyzz[k] += f00 * gr_x * gr_y * gr_z * gr_z;
                        ptr_buffer_x_xyzz[k] += fg_1 * gr_x * gr_x * gr_y * gr_z * gr_z + f00 * gr_y * gr_z * gr_z;
                        ptr_buffer_y_xyzz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z * gr_z + f00 * gr_x * gr_z * gr_z;
                        ptr_buffer_z_xyzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z * gr_z + 2.0 * f00 * gr_x * gr_y * gr_z;

                        ptr_buffer_0_xzzz[k] += f00 * gr_x * gr_z * gr_z * gr_z;
                        ptr_buffer_x_xzzz[k] += fg_1 * gr_x * gr_x * gr_z * gr_z * gr_z + f00 * gr_z * gr_z * gr_z;
                        ptr_buffer_y_xzzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_z_xzzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z * gr_z + 3.0 * f00 * gr_x * gr_z * gr_z;

                        ptr_buffer_0_yyyy[k] += f00 * gr_y * gr_y * gr_y * gr_y;
                        ptr_buffer_x_yyyy[k] += fg_1 * gr_x * gr_y * gr_y * gr_y * gr_y;
                        ptr_buffer_y_yyyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_y * gr_y + 4.0 * f00 * gr_y * gr_y * gr_y;
                        ptr_buffer_z_yyyy[k] += fg_1 * gr_y * gr_y * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yyyz[k] += f00 * gr_y * gr_y * gr_y * gr_z;
                        ptr_buffer_x_yyyz[k] += fg_1 * gr_x * gr_y * gr_y * gr_y * gr_z;
                        ptr_buffer_y_yyyz[k] += fg_1 * gr_y * gr_y * gr_y * gr_y * gr_z + 3.0 * f00 * gr_y * gr_y * gr_z;
                        ptr_buffer_z_yyyz[k] += fg_1 * gr_y * gr_y * gr_y * gr_z * gr_z + f00 * gr_y * gr_y * gr_y;

                        ptr_buffer_0_yyzz[k] += f00 * gr_y * gr_y * gr_z * gr_z;
                        ptr_buffer_x_yyzz[k] += fg_1 * gr_x * gr_y * gr_y * gr_z * gr_z;
                        ptr_buffer_y_yyzz[k] += fg_1 * gr_y * gr_y * gr_y * gr_z * gr_z + 2.0 * f00 * gr_y * gr_z * gr_z;
                        ptr_buffer_z_yyzz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z * gr_z + 2.0 * f00 * gr_y * gr_y * gr_z;

                        ptr_buffer_0_yzzz[k] += f00 * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_x_yzzz[k] += fg_1 * gr_x * gr_y * gr_z * gr_z * gr_z;
                        ptr_buffer_y_yzzz[k] += fg_1 * gr_y * gr_y * gr_z * gr_z * gr_z + f00 * gr_z * gr_z * gr_z;
                        ptr_buffer_z_yzzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z * gr_z + 3.0 * f00 * gr_y * gr_z * gr_z;

                        ptr_buffer_0_zzzz[k] += f00 * gr_z * gr_z * gr_z * gr_z;
                        ptr_buffer_x_zzzz[k] += fg_1 * gr_x * gr_z * gr_z * gr_z * gr_z;
                        ptr_buffer_y_zzzz[k] += fg_1 * gr_y * gr_z * gr_z * gr_z * gr_z;
                        ptr_buffer_z_zzzz[k] += fg_1 * gr_z * gr_z * gr_z * gr_z * gr_z + 4.0 * f00 * gr_z * gr_z * gr_z;
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

}  // namespace gtoval

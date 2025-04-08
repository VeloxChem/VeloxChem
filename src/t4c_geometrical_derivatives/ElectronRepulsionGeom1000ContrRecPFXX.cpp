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

#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_pfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_pfxx,
                                            const size_t idx_sfxx,
                                            const size_t idx_geom_10_sfxx,
                                            const size_t idx_geom_10_sgxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SFSS

            const auto sf_off = idx_sfxx + i * dcomps + j;

            auto g_0_xxx = cbuffer.data(sf_off + 0 * ccomps * dcomps);

            auto g_0_xxy = cbuffer.data(sf_off + 1 * ccomps * dcomps);

            auto g_0_xxz = cbuffer.data(sf_off + 2 * ccomps * dcomps);

            auto g_0_xyy = cbuffer.data(sf_off + 3 * ccomps * dcomps);

            auto g_0_xyz = cbuffer.data(sf_off + 4 * ccomps * dcomps);

            auto g_0_xzz = cbuffer.data(sf_off + 5 * ccomps * dcomps);

            auto g_0_yyy = cbuffer.data(sf_off + 6 * ccomps * dcomps);

            auto g_0_yyz = cbuffer.data(sf_off + 7 * ccomps * dcomps);

            auto g_0_yzz = cbuffer.data(sf_off + 8 * ccomps * dcomps);

            auto g_0_zzz = cbuffer.data(sf_off + 9 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_10_off = idx_geom_10_sfxx + i * dcomps + j;

            auto g_x_0_0_xxx = cbuffer.data(sf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxy = cbuffer.data(sf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxz = cbuffer.data(sf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xyy = cbuffer.data(sf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xyz = cbuffer.data(sf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xzz = cbuffer.data(sf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_yyy = cbuffer.data(sf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_yyz = cbuffer.data(sf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_yzz = cbuffer.data(sf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_zzz = cbuffer.data(sf_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_0_xxx = cbuffer.data(sf_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_0_xxy = cbuffer.data(sf_geom_10_off + 11 * ccomps * dcomps);

            auto g_y_0_0_xxz = cbuffer.data(sf_geom_10_off + 12 * ccomps * dcomps);

            auto g_y_0_0_xyy = cbuffer.data(sf_geom_10_off + 13 * ccomps * dcomps);

            auto g_y_0_0_xyz = cbuffer.data(sf_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xzz = cbuffer.data(sf_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_0_yyy = cbuffer.data(sf_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_0_yyz = cbuffer.data(sf_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_0_yzz = cbuffer.data(sf_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_0_zzz = cbuffer.data(sf_geom_10_off + 19 * ccomps * dcomps);

            auto g_z_0_0_xxx = cbuffer.data(sf_geom_10_off + 20 * ccomps * dcomps);

            auto g_z_0_0_xxy = cbuffer.data(sf_geom_10_off + 21 * ccomps * dcomps);

            auto g_z_0_0_xxz = cbuffer.data(sf_geom_10_off + 22 * ccomps * dcomps);

            auto g_z_0_0_xyy = cbuffer.data(sf_geom_10_off + 23 * ccomps * dcomps);

            auto g_z_0_0_xyz = cbuffer.data(sf_geom_10_off + 24 * ccomps * dcomps);

            auto g_z_0_0_xzz = cbuffer.data(sf_geom_10_off + 25 * ccomps * dcomps);

            auto g_z_0_0_yyy = cbuffer.data(sf_geom_10_off + 26 * ccomps * dcomps);

            auto g_z_0_0_yyz = cbuffer.data(sf_geom_10_off + 27 * ccomps * dcomps);

            auto g_z_0_0_yzz = cbuffer.data(sf_geom_10_off + 28 * ccomps * dcomps);

            auto g_z_0_0_zzz = cbuffer.data(sf_geom_10_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_10_off = idx_geom_10_sgxx + i * dcomps + j;

            auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 29 * ccomps * dcomps);

            auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * ccomps * dcomps);

            auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 31 * ccomps * dcomps);

            auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 32 * ccomps * dcomps);

            auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 33 * ccomps * dcomps);

            auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 34 * ccomps * dcomps);

            auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 44 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pfxx

            const auto pf_geom_10_off = idx_geom_10_pfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_x_0_0_xxx, g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxx[k] = -g_0_xxx[k] - g_x_0_0_xxx[k] * ab_x + g_x_0_0_xxxx[k];

                g_x_0_x_xxy[k] = -g_0_xxy[k] - g_x_0_0_xxy[k] * ab_x + g_x_0_0_xxxy[k];

                g_x_0_x_xxz[k] = -g_0_xxz[k] - g_x_0_0_xxz[k] * ab_x + g_x_0_0_xxxz[k];

                g_x_0_x_xyy[k] = -g_0_xyy[k] - g_x_0_0_xyy[k] * ab_x + g_x_0_0_xxyy[k];

                g_x_0_x_xyz[k] = -g_0_xyz[k] - g_x_0_0_xyz[k] * ab_x + g_x_0_0_xxyz[k];

                g_x_0_x_xzz[k] = -g_0_xzz[k] - g_x_0_0_xzz[k] * ab_x + g_x_0_0_xxzz[k];

                g_x_0_x_yyy[k] = -g_0_yyy[k] - g_x_0_0_yyy[k] * ab_x + g_x_0_0_xyyy[k];

                g_x_0_x_yyz[k] = -g_0_yyz[k] - g_x_0_0_yyz[k] * ab_x + g_x_0_0_xyyz[k];

                g_x_0_x_yzz[k] = -g_0_yzz[k] - g_x_0_0_yzz[k] * ab_x + g_x_0_0_xyzz[k];

                g_x_0_x_zzz[k] = -g_0_zzz[k] - g_x_0_0_zzz[k] * ab_x + g_x_0_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxxy, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_y_xxx, g_x_0_y_xxy, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yzz, g_x_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxx[k] = -g_x_0_0_xxx[k] * ab_y + g_x_0_0_xxxy[k];

                g_x_0_y_xxy[k] = -g_x_0_0_xxy[k] * ab_y + g_x_0_0_xxyy[k];

                g_x_0_y_xxz[k] = -g_x_0_0_xxz[k] * ab_y + g_x_0_0_xxyz[k];

                g_x_0_y_xyy[k] = -g_x_0_0_xyy[k] * ab_y + g_x_0_0_xyyy[k];

                g_x_0_y_xyz[k] = -g_x_0_0_xyz[k] * ab_y + g_x_0_0_xyyz[k];

                g_x_0_y_xzz[k] = -g_x_0_0_xzz[k] * ab_y + g_x_0_0_xyzz[k];

                g_x_0_y_yyy[k] = -g_x_0_0_yyy[k] * ab_y + g_x_0_0_yyyy[k];

                g_x_0_y_yyz[k] = -g_x_0_0_yyz[k] * ab_y + g_x_0_0_yyyz[k];

                g_x_0_y_yzz[k] = -g_x_0_0_yzz[k] * ab_y + g_x_0_0_yyzz[k];

                g_x_0_y_zzz[k] = -g_x_0_0_zzz[k] * ab_y + g_x_0_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_0_zzzz, g_x_0_z_xxx, g_x_0_z_xxy, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yzz, g_x_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxx[k] = -g_x_0_0_xxx[k] * ab_z + g_x_0_0_xxxz[k];

                g_x_0_z_xxy[k] = -g_x_0_0_xxy[k] * ab_z + g_x_0_0_xxyz[k];

                g_x_0_z_xxz[k] = -g_x_0_0_xxz[k] * ab_z + g_x_0_0_xxzz[k];

                g_x_0_z_xyy[k] = -g_x_0_0_xyy[k] * ab_z + g_x_0_0_xyyz[k];

                g_x_0_z_xyz[k] = -g_x_0_0_xyz[k] * ab_z + g_x_0_0_xyzz[k];

                g_x_0_z_xzz[k] = -g_x_0_0_xzz[k] * ab_z + g_x_0_0_xzzz[k];

                g_x_0_z_yyy[k] = -g_x_0_0_yyy[k] * ab_z + g_x_0_0_yyyz[k];

                g_x_0_z_yyz[k] = -g_x_0_0_yyz[k] * ab_z + g_x_0_0_yyzz[k];

                g_x_0_z_yzz[k] = -g_x_0_0_yzz[k] * ab_z + g_x_0_0_yzzz[k];

                g_x_0_z_zzz[k] = -g_x_0_0_zzz[k] * ab_z + g_x_0_0_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxx[k] = -g_y_0_0_xxx[k] * ab_x + g_y_0_0_xxxx[k];

                g_y_0_x_xxy[k] = -g_y_0_0_xxy[k] * ab_x + g_y_0_0_xxxy[k];

                g_y_0_x_xxz[k] = -g_y_0_0_xxz[k] * ab_x + g_y_0_0_xxxz[k];

                g_y_0_x_xyy[k] = -g_y_0_0_xyy[k] * ab_x + g_y_0_0_xxyy[k];

                g_y_0_x_xyz[k] = -g_y_0_0_xyz[k] * ab_x + g_y_0_0_xxyz[k];

                g_y_0_x_xzz[k] = -g_y_0_0_xzz[k] * ab_x + g_y_0_0_xxzz[k];

                g_y_0_x_yyy[k] = -g_y_0_0_yyy[k] * ab_x + g_y_0_0_xyyy[k];

                g_y_0_x_yyz[k] = -g_y_0_0_yyz[k] * ab_x + g_y_0_0_xyyz[k];

                g_y_0_x_yzz[k] = -g_y_0_0_yzz[k] * ab_x + g_y_0_0_xyzz[k];

                g_y_0_x_zzz[k] = -g_y_0_0_zzz[k] * ab_x + g_y_0_0_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_y_0_0_xxx, g_y_0_0_xxxy, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxx[k] = -g_0_xxx[k] - g_y_0_0_xxx[k] * ab_y + g_y_0_0_xxxy[k];

                g_y_0_y_xxy[k] = -g_0_xxy[k] - g_y_0_0_xxy[k] * ab_y + g_y_0_0_xxyy[k];

                g_y_0_y_xxz[k] = -g_0_xxz[k] - g_y_0_0_xxz[k] * ab_y + g_y_0_0_xxyz[k];

                g_y_0_y_xyy[k] = -g_0_xyy[k] - g_y_0_0_xyy[k] * ab_y + g_y_0_0_xyyy[k];

                g_y_0_y_xyz[k] = -g_0_xyz[k] - g_y_0_0_xyz[k] * ab_y + g_y_0_0_xyyz[k];

                g_y_0_y_xzz[k] = -g_0_xzz[k] - g_y_0_0_xzz[k] * ab_y + g_y_0_0_xyzz[k];

                g_y_0_y_yyy[k] = -g_0_yyy[k] - g_y_0_0_yyy[k] * ab_y + g_y_0_0_yyyy[k];

                g_y_0_y_yyz[k] = -g_0_yyz[k] - g_y_0_0_yyz[k] * ab_y + g_y_0_0_yyyz[k];

                g_y_0_y_yzz[k] = -g_0_yzz[k] - g_y_0_0_yzz[k] * ab_y + g_y_0_0_yyzz[k];

                g_y_0_y_zzz[k] = -g_0_zzz[k] - g_y_0_0_zzz[k] * ab_y + g_y_0_0_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxx, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_0_zzzz, g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxx[k] = -g_y_0_0_xxx[k] * ab_z + g_y_0_0_xxxz[k];

                g_y_0_z_xxy[k] = -g_y_0_0_xxy[k] * ab_z + g_y_0_0_xxyz[k];

                g_y_0_z_xxz[k] = -g_y_0_0_xxz[k] * ab_z + g_y_0_0_xxzz[k];

                g_y_0_z_xyy[k] = -g_y_0_0_xyy[k] * ab_z + g_y_0_0_xyyz[k];

                g_y_0_z_xyz[k] = -g_y_0_0_xyz[k] * ab_z + g_y_0_0_xyzz[k];

                g_y_0_z_xzz[k] = -g_y_0_0_xzz[k] * ab_z + g_y_0_0_xzzz[k];

                g_y_0_z_yyy[k] = -g_y_0_0_yyy[k] * ab_z + g_y_0_0_yyyz[k];

                g_y_0_z_yyz[k] = -g_y_0_0_yyz[k] * ab_z + g_y_0_0_yyzz[k];

                g_y_0_z_yzz[k] = -g_y_0_0_yzz[k] * ab_z + g_y_0_0_yzzz[k];

                g_y_0_z_zzz[k] = -g_y_0_0_zzz[k] * ab_z + g_y_0_0_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 68 * ccomps * dcomps);

            auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxx, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxx[k] = -g_z_0_0_xxx[k] * ab_x + g_z_0_0_xxxx[k];

                g_z_0_x_xxy[k] = -g_z_0_0_xxy[k] * ab_x + g_z_0_0_xxxy[k];

                g_z_0_x_xxz[k] = -g_z_0_0_xxz[k] * ab_x + g_z_0_0_xxxz[k];

                g_z_0_x_xyy[k] = -g_z_0_0_xyy[k] * ab_x + g_z_0_0_xxyy[k];

                g_z_0_x_xyz[k] = -g_z_0_0_xyz[k] * ab_x + g_z_0_0_xxyz[k];

                g_z_0_x_xzz[k] = -g_z_0_0_xzz[k] * ab_x + g_z_0_0_xxzz[k];

                g_z_0_x_yyy[k] = -g_z_0_0_yyy[k] * ab_x + g_z_0_0_xyyy[k];

                g_z_0_x_yyz[k] = -g_z_0_0_yyz[k] * ab_x + g_z_0_0_xyyz[k];

                g_z_0_x_yzz[k] = -g_z_0_0_yzz[k] * ab_x + g_z_0_0_xyzz[k];

                g_z_0_x_zzz[k] = -g_z_0_0_zzz[k] * ab_x + g_z_0_0_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxx, g_z_0_0_xxxy, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxx[k] = -g_z_0_0_xxx[k] * ab_y + g_z_0_0_xxxy[k];

                g_z_0_y_xxy[k] = -g_z_0_0_xxy[k] * ab_y + g_z_0_0_xxyy[k];

                g_z_0_y_xxz[k] = -g_z_0_0_xxz[k] * ab_y + g_z_0_0_xxyz[k];

                g_z_0_y_xyy[k] = -g_z_0_0_xyy[k] * ab_y + g_z_0_0_xyyy[k];

                g_z_0_y_xyz[k] = -g_z_0_0_xyz[k] * ab_y + g_z_0_0_xyyz[k];

                g_z_0_y_xzz[k] = -g_z_0_0_xzz[k] * ab_y + g_z_0_0_xyzz[k];

                g_z_0_y_yyy[k] = -g_z_0_0_yyy[k] * ab_y + g_z_0_0_yyyy[k];

                g_z_0_y_yyz[k] = -g_z_0_0_yyz[k] * ab_y + g_z_0_0_yyyz[k];

                g_z_0_y_yzz[k] = -g_z_0_0_yzz[k] * ab_y + g_z_0_0_yyzz[k];

                g_z_0_y_zzz[k] = -g_z_0_0_zzz[k] * ab_y + g_z_0_0_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 83 * ccomps * dcomps);

            auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_z_0_0_xxx, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_0_zzzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxx[k] = -g_0_xxx[k] - g_z_0_0_xxx[k] * ab_z + g_z_0_0_xxxz[k];

                g_z_0_z_xxy[k] = -g_0_xxy[k] - g_z_0_0_xxy[k] * ab_z + g_z_0_0_xxyz[k];

                g_z_0_z_xxz[k] = -g_0_xxz[k] - g_z_0_0_xxz[k] * ab_z + g_z_0_0_xxzz[k];

                g_z_0_z_xyy[k] = -g_0_xyy[k] - g_z_0_0_xyy[k] * ab_z + g_z_0_0_xyyz[k];

                g_z_0_z_xyz[k] = -g_0_xyz[k] - g_z_0_0_xyz[k] * ab_z + g_z_0_0_xyzz[k];

                g_z_0_z_xzz[k] = -g_0_xzz[k] - g_z_0_0_xzz[k] * ab_z + g_z_0_0_xzzz[k];

                g_z_0_z_yyy[k] = -g_0_yyy[k] - g_z_0_0_yyy[k] * ab_z + g_z_0_0_yyyz[k];

                g_z_0_z_yyz[k] = -g_0_yyz[k] - g_z_0_0_yyz[k] * ab_z + g_z_0_0_yyzz[k];

                g_z_0_z_yzz[k] = -g_0_yzz[k] - g_z_0_0_yzz[k] * ab_z + g_z_0_0_yzzz[k];

                g_z_0_z_zzz[k] = -g_0_zzz[k] - g_z_0_0_zzz[k] * ab_z + g_z_0_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace


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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xpf(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xpf,
                                        const size_t idx_xsf,
                                        const size_t idx_geom_10_xsf,
                                        const size_t idx_geom_10_xsg,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SSF

        const auto sf_off = idx_xsf + i * 10;

        auto g_0_xxx = cbuffer.data(sf_off + 0);

        auto g_0_xxy = cbuffer.data(sf_off + 1);

        auto g_0_xxz = cbuffer.data(sf_off + 2);

        auto g_0_xyy = cbuffer.data(sf_off + 3);

        auto g_0_xyz = cbuffer.data(sf_off + 4);

        auto g_0_xzz = cbuffer.data(sf_off + 5);

        auto g_0_yyy = cbuffer.data(sf_off + 6);

        auto g_0_yyz = cbuffer.data(sf_off + 7);

        auto g_0_yzz = cbuffer.data(sf_off + 8);

        auto g_0_zzz = cbuffer.data(sf_off + 9);

        /// Set up components of auxilary buffer : SSF

        const auto sf_geom_10_off = idx_geom_10_xsf + i * 10;

        auto g_x_0_0_xxx = cbuffer.data(sf_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxy = cbuffer.data(sf_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxz = cbuffer.data(sf_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xyy = cbuffer.data(sf_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xyz = cbuffer.data(sf_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_yyy = cbuffer.data(sf_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_yyz = cbuffer.data(sf_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_yzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_zzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 9);

        auto g_y_0_0_xxx = cbuffer.data(sf_geom_10_off + 10 * acomps + 0);

        auto g_y_0_0_xxy = cbuffer.data(sf_geom_10_off + 10 * acomps + 1);

        auto g_y_0_0_xxz = cbuffer.data(sf_geom_10_off + 10 * acomps + 2);

        auto g_y_0_0_xyy = cbuffer.data(sf_geom_10_off + 10 * acomps + 3);

        auto g_y_0_0_xyz = cbuffer.data(sf_geom_10_off + 10 * acomps + 4);

        auto g_y_0_0_xzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 5);

        auto g_y_0_0_yyy = cbuffer.data(sf_geom_10_off + 10 * acomps + 6);

        auto g_y_0_0_yyz = cbuffer.data(sf_geom_10_off + 10 * acomps + 7);

        auto g_y_0_0_yzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 8);

        auto g_y_0_0_zzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 9);

        auto g_z_0_0_xxx = cbuffer.data(sf_geom_10_off + 20 * acomps + 0);

        auto g_z_0_0_xxy = cbuffer.data(sf_geom_10_off + 20 * acomps + 1);

        auto g_z_0_0_xxz = cbuffer.data(sf_geom_10_off + 20 * acomps + 2);

        auto g_z_0_0_xyy = cbuffer.data(sf_geom_10_off + 20 * acomps + 3);

        auto g_z_0_0_xyz = cbuffer.data(sf_geom_10_off + 20 * acomps + 4);

        auto g_z_0_0_xzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 5);

        auto g_z_0_0_yyy = cbuffer.data(sf_geom_10_off + 20 * acomps + 6);

        auto g_z_0_0_yyz = cbuffer.data(sf_geom_10_off + 20 * acomps + 7);

        auto g_z_0_0_yzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 8);

        auto g_z_0_0_zzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 9);

        /// Set up components of auxilary buffer : SSG

        const auto sg_geom_10_off = idx_geom_10_xsg + i * 15;

        auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 9);

        auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 0 * acomps + 10);

        auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 0 * acomps + 11);

        auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 12);

        auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 13);

        auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 0 * acomps + 14);

        auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * acomps + 0);

        auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 15 * acomps + 1);

        auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 15 * acomps + 2);

        auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 15 * acomps + 3);

        auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 15 * acomps + 4);

        auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 5);

        auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 15 * acomps + 6);

        auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 15 * acomps + 7);

        auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 8);

        auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 9);

        auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 15 * acomps + 10);

        auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 15 * acomps + 11);

        auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 12);

        auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 13);

        auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 15 * acomps + 14);

        auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * acomps + 0);

        auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 30 * acomps + 1);

        auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 30 * acomps + 2);

        auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 30 * acomps + 3);

        auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 30 * acomps + 4);

        auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 5);

        auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 30 * acomps + 6);

        auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 30 * acomps + 7);

        auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 8);

        auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 9);

        auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 30 * acomps + 10);

        auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 30 * acomps + 11);

        auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 12);

        auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 13);

        auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 30 * acomps + 14);

        /// set up bra offset for contr_buffer_xxpf

        const auto pf_geom_10_off = idx_geom_10_xpf + i * 30;

        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_x_0_0_xxx, g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yzz, g_x_0_0_zzz, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_x_xxx[k] = -g_0_xxx[k] - g_x_0_0_xxx[k] * cd_x[k] + g_x_0_0_xxxx[k];

            g_x_0_x_xxy[k] = -g_0_xxy[k] - g_x_0_0_xxy[k] * cd_x[k] + g_x_0_0_xxxy[k];

            g_x_0_x_xxz[k] = -g_0_xxz[k] - g_x_0_0_xxz[k] * cd_x[k] + g_x_0_0_xxxz[k];

            g_x_0_x_xyy[k] = -g_0_xyy[k] - g_x_0_0_xyy[k] * cd_x[k] + g_x_0_0_xxyy[k];

            g_x_0_x_xyz[k] = -g_0_xyz[k] - g_x_0_0_xyz[k] * cd_x[k] + g_x_0_0_xxyz[k];

            g_x_0_x_xzz[k] = -g_0_xzz[k] - g_x_0_0_xzz[k] * cd_x[k] + g_x_0_0_xxzz[k];

            g_x_0_x_yyy[k] = -g_0_yyy[k] - g_x_0_0_yyy[k] * cd_x[k] + g_x_0_0_xyyy[k];

            g_x_0_x_yyz[k] = -g_0_yyz[k] - g_x_0_0_yyz[k] * cd_x[k] + g_x_0_0_xyyz[k];

            g_x_0_x_yzz[k] = -g_0_yzz[k] - g_x_0_0_yzz[k] * cd_x[k] + g_x_0_0_xyzz[k];

            g_x_0_x_zzz[k] = -g_0_zzz[k] - g_x_0_0_zzz[k] * cd_x[k] + g_x_0_0_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 19);

        #pragma omp simd aligned(cd_y, g_x_0_0_xxx, g_x_0_0_xxxy, g_x_0_0_xxy, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xyy, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_yyy, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_y_xxx, g_x_0_y_xxy, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yzz, g_x_0_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_y_xxx[k] = -g_x_0_0_xxx[k] * cd_y[k] + g_x_0_0_xxxy[k];

            g_x_0_y_xxy[k] = -g_x_0_0_xxy[k] * cd_y[k] + g_x_0_0_xxyy[k];

            g_x_0_y_xxz[k] = -g_x_0_0_xxz[k] * cd_y[k] + g_x_0_0_xxyz[k];

            g_x_0_y_xyy[k] = -g_x_0_0_xyy[k] * cd_y[k] + g_x_0_0_xyyy[k];

            g_x_0_y_xyz[k] = -g_x_0_0_xyz[k] * cd_y[k] + g_x_0_0_xyyz[k];

            g_x_0_y_xzz[k] = -g_x_0_0_xzz[k] * cd_y[k] + g_x_0_0_xyzz[k];

            g_x_0_y_yyy[k] = -g_x_0_0_yyy[k] * cd_y[k] + g_x_0_0_yyyy[k];

            g_x_0_y_yyz[k] = -g_x_0_0_yyz[k] * cd_y[k] + g_x_0_0_yyyz[k];

            g_x_0_y_yzz[k] = -g_x_0_0_yzz[k] * cd_y[k] + g_x_0_0_yyzz[k];

            g_x_0_y_zzz[k] = -g_x_0_0_zzz[k] * cd_y[k] + g_x_0_0_yzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_x_0_0_xxx, g_x_0_0_xxxz, g_x_0_0_xxy, g_x_0_0_xxyz, g_x_0_0_xxz, g_x_0_0_xxzz, g_x_0_0_xyy, g_x_0_0_xyyz, g_x_0_0_xyz, g_x_0_0_xyzz, g_x_0_0_xzz, g_x_0_0_xzzz, g_x_0_0_yyy, g_x_0_0_yyyz, g_x_0_0_yyz, g_x_0_0_yyzz, g_x_0_0_yzz, g_x_0_0_yzzz, g_x_0_0_zzz, g_x_0_0_zzzz, g_x_0_z_xxx, g_x_0_z_xxy, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yzz, g_x_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_z_xxx[k] = -g_x_0_0_xxx[k] * cd_z[k] + g_x_0_0_xxxz[k];

            g_x_0_z_xxy[k] = -g_x_0_0_xxy[k] * cd_z[k] + g_x_0_0_xxyz[k];

            g_x_0_z_xxz[k] = -g_x_0_0_xxz[k] * cd_z[k] + g_x_0_0_xxzz[k];

            g_x_0_z_xyy[k] = -g_x_0_0_xyy[k] * cd_z[k] + g_x_0_0_xyyz[k];

            g_x_0_z_xyz[k] = -g_x_0_0_xyz[k] * cd_z[k] + g_x_0_0_xyzz[k];

            g_x_0_z_xzz[k] = -g_x_0_0_xzz[k] * cd_z[k] + g_x_0_0_xzzz[k];

            g_x_0_z_yyy[k] = -g_x_0_0_yyy[k] * cd_z[k] + g_x_0_0_yyyz[k];

            g_x_0_z_yyz[k] = -g_x_0_0_yyz[k] * cd_z[k] + g_x_0_0_yyzz[k];

            g_x_0_z_yzz[k] = -g_x_0_0_yzz[k] * cd_z[k] + g_x_0_0_yzzz[k];

            g_x_0_z_zzz[k] = -g_x_0_0_zzz[k] * cd_z[k] + g_x_0_0_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps  + 0);

        auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 1);

        auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 2);

        auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 3);

        auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 4);

        auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 5);

        auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 6);

        auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 7);

        auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 8);

        auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_y_0_0_xxx, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yzz, g_y_0_0_zzz, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_x_xxx[k] = -g_y_0_0_xxx[k] * cd_x[k] + g_y_0_0_xxxx[k];

            g_y_0_x_xxy[k] = -g_y_0_0_xxy[k] * cd_x[k] + g_y_0_0_xxxy[k];

            g_y_0_x_xxz[k] = -g_y_0_0_xxz[k] * cd_x[k] + g_y_0_0_xxxz[k];

            g_y_0_x_xyy[k] = -g_y_0_0_xyy[k] * cd_x[k] + g_y_0_0_xxyy[k];

            g_y_0_x_xyz[k] = -g_y_0_0_xyz[k] * cd_x[k] + g_y_0_0_xxyz[k];

            g_y_0_x_xzz[k] = -g_y_0_0_xzz[k] * cd_x[k] + g_y_0_0_xxzz[k];

            g_y_0_x_yyy[k] = -g_y_0_0_yyy[k] * cd_x[k] + g_y_0_0_xyyy[k];

            g_y_0_x_yyz[k] = -g_y_0_0_yyz[k] * cd_x[k] + g_y_0_0_xyyz[k];

            g_y_0_x_yzz[k] = -g_y_0_0_yzz[k] * cd_x[k] + g_y_0_0_xyzz[k];

            g_y_0_x_zzz[k] = -g_y_0_0_zzz[k] * cd_x[k] + g_y_0_0_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps  + 10);

        auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 11);

        auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 12);

        auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 13);

        auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 14);

        auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 15);

        auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 16);

        auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 17);

        auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 18);

        auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 19);

        #pragma omp simd aligned(cd_y, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_y_0_0_xxx, g_y_0_0_xxxy, g_y_0_0_xxy, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xyy, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_yyy, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_y_xxx[k] = -g_0_xxx[k] - g_y_0_0_xxx[k] * cd_y[k] + g_y_0_0_xxxy[k];

            g_y_0_y_xxy[k] = -g_0_xxy[k] - g_y_0_0_xxy[k] * cd_y[k] + g_y_0_0_xxyy[k];

            g_y_0_y_xxz[k] = -g_0_xxz[k] - g_y_0_0_xxz[k] * cd_y[k] + g_y_0_0_xxyz[k];

            g_y_0_y_xyy[k] = -g_0_xyy[k] - g_y_0_0_xyy[k] * cd_y[k] + g_y_0_0_xyyy[k];

            g_y_0_y_xyz[k] = -g_0_xyz[k] - g_y_0_0_xyz[k] * cd_y[k] + g_y_0_0_xyyz[k];

            g_y_0_y_xzz[k] = -g_0_xzz[k] - g_y_0_0_xzz[k] * cd_y[k] + g_y_0_0_xyzz[k];

            g_y_0_y_yyy[k] = -g_0_yyy[k] - g_y_0_0_yyy[k] * cd_y[k] + g_y_0_0_yyyy[k];

            g_y_0_y_yyz[k] = -g_0_yyz[k] - g_y_0_0_yyz[k] * cd_y[k] + g_y_0_0_yyyz[k];

            g_y_0_y_yzz[k] = -g_0_yzz[k] - g_y_0_0_yzz[k] * cd_y[k] + g_y_0_0_yyzz[k];

            g_y_0_y_zzz[k] = -g_0_zzz[k] - g_y_0_0_zzz[k] * cd_y[k] + g_y_0_0_yzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps  + 20);

        auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 21);

        auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 22);

        auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 23);

        auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 24);

        auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 25);

        auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps  + 26);

        auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 27);

        auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 28);

        auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_y_0_0_xxx, g_y_0_0_xxxz, g_y_0_0_xxy, g_y_0_0_xxyz, g_y_0_0_xxz, g_y_0_0_xxzz, g_y_0_0_xyy, g_y_0_0_xyyz, g_y_0_0_xyz, g_y_0_0_xyzz, g_y_0_0_xzz, g_y_0_0_xzzz, g_y_0_0_yyy, g_y_0_0_yyyz, g_y_0_0_yyz, g_y_0_0_yyzz, g_y_0_0_yzz, g_y_0_0_yzzz, g_y_0_0_zzz, g_y_0_0_zzzz, g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_z_xxx[k] = -g_y_0_0_xxx[k] * cd_z[k] + g_y_0_0_xxxz[k];

            g_y_0_z_xxy[k] = -g_y_0_0_xxy[k] * cd_z[k] + g_y_0_0_xxyz[k];

            g_y_0_z_xxz[k] = -g_y_0_0_xxz[k] * cd_z[k] + g_y_0_0_xxzz[k];

            g_y_0_z_xyy[k] = -g_y_0_0_xyy[k] * cd_z[k] + g_y_0_0_xyyz[k];

            g_y_0_z_xyz[k] = -g_y_0_0_xyz[k] * cd_z[k] + g_y_0_0_xyzz[k];

            g_y_0_z_xzz[k] = -g_y_0_0_xzz[k] * cd_z[k] + g_y_0_0_xzzz[k];

            g_y_0_z_yyy[k] = -g_y_0_0_yyy[k] * cd_z[k] + g_y_0_0_yyyz[k];

            g_y_0_z_yyz[k] = -g_y_0_0_yyz[k] * cd_z[k] + g_y_0_0_yyzz[k];

            g_y_0_z_yzz[k] = -g_y_0_0_yzz[k] * cd_z[k] + g_y_0_0_yzzz[k];

            g_y_0_z_zzz[k] = -g_y_0_0_zzz[k] * cd_z[k] + g_y_0_0_zzzz[k];
        }
        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps  + 0);

        auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 1);

        auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 2);

        auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 3);

        auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 4);

        auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 5);

        auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 6);

        auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 7);

        auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 8);

        auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 9);

        #pragma omp simd aligned(cd_x, g_z_0_0_xxx, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yzz, g_z_0_0_zzz, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_x_xxx[k] = -g_z_0_0_xxx[k] * cd_x[k] + g_z_0_0_xxxx[k];

            g_z_0_x_xxy[k] = -g_z_0_0_xxy[k] * cd_x[k] + g_z_0_0_xxxy[k];

            g_z_0_x_xxz[k] = -g_z_0_0_xxz[k] * cd_x[k] + g_z_0_0_xxxz[k];

            g_z_0_x_xyy[k] = -g_z_0_0_xyy[k] * cd_x[k] + g_z_0_0_xxyy[k];

            g_z_0_x_xyz[k] = -g_z_0_0_xyz[k] * cd_x[k] + g_z_0_0_xxyz[k];

            g_z_0_x_xzz[k] = -g_z_0_0_xzz[k] * cd_x[k] + g_z_0_0_xxzz[k];

            g_z_0_x_yyy[k] = -g_z_0_0_yyy[k] * cd_x[k] + g_z_0_0_xyyy[k];

            g_z_0_x_yyz[k] = -g_z_0_0_yyz[k] * cd_x[k] + g_z_0_0_xyyz[k];

            g_z_0_x_yzz[k] = -g_z_0_0_yzz[k] * cd_x[k] + g_z_0_0_xyzz[k];

            g_z_0_x_zzz[k] = -g_z_0_0_zzz[k] * cd_x[k] + g_z_0_0_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps  + 10);

        auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 11);

        auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 12);

        auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 13);

        auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 14);

        auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 15);

        auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 16);

        auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 17);

        auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 18);

        auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 19);

        #pragma omp simd aligned(cd_y, g_z_0_0_xxx, g_z_0_0_xxxy, g_z_0_0_xxy, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xyy, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_yyy, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_y_xxx[k] = -g_z_0_0_xxx[k] * cd_y[k] + g_z_0_0_xxxy[k];

            g_z_0_y_xxy[k] = -g_z_0_0_xxy[k] * cd_y[k] + g_z_0_0_xxyy[k];

            g_z_0_y_xxz[k] = -g_z_0_0_xxz[k] * cd_y[k] + g_z_0_0_xxyz[k];

            g_z_0_y_xyy[k] = -g_z_0_0_xyy[k] * cd_y[k] + g_z_0_0_xyyy[k];

            g_z_0_y_xyz[k] = -g_z_0_0_xyz[k] * cd_y[k] + g_z_0_0_xyyz[k];

            g_z_0_y_xzz[k] = -g_z_0_0_xzz[k] * cd_y[k] + g_z_0_0_xyzz[k];

            g_z_0_y_yyy[k] = -g_z_0_0_yyy[k] * cd_y[k] + g_z_0_0_yyyy[k];

            g_z_0_y_yyz[k] = -g_z_0_0_yyz[k] * cd_y[k] + g_z_0_0_yyyz[k];

            g_z_0_y_yzz[k] = -g_z_0_0_yzz[k] * cd_y[k] + g_z_0_0_yyzz[k];

            g_z_0_y_zzz[k] = -g_z_0_0_zzz[k] * cd_y[k] + g_z_0_0_yzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps  + 20);

        auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 21);

        auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 22);

        auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 23);

        auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 24);

        auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 25);

        auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps  + 26);

        auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 27);

        auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 28);

        auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz, g_z_0_0_xxx, g_z_0_0_xxxz, g_z_0_0_xxy, g_z_0_0_xxyz, g_z_0_0_xxz, g_z_0_0_xxzz, g_z_0_0_xyy, g_z_0_0_xyyz, g_z_0_0_xyz, g_z_0_0_xyzz, g_z_0_0_xzz, g_z_0_0_xzzz, g_z_0_0_yyy, g_z_0_0_yyyz, g_z_0_0_yyz, g_z_0_0_yyzz, g_z_0_0_yzz, g_z_0_0_yzzz, g_z_0_0_zzz, g_z_0_0_zzzz, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_z_xxx[k] = -g_0_xxx[k] - g_z_0_0_xxx[k] * cd_z[k] + g_z_0_0_xxxz[k];

            g_z_0_z_xxy[k] = -g_0_xxy[k] - g_z_0_0_xxy[k] * cd_z[k] + g_z_0_0_xxyz[k];

            g_z_0_z_xxz[k] = -g_0_xxz[k] - g_z_0_0_xxz[k] * cd_z[k] + g_z_0_0_xxzz[k];

            g_z_0_z_xyy[k] = -g_0_xyy[k] - g_z_0_0_xyy[k] * cd_z[k] + g_z_0_0_xyyz[k];

            g_z_0_z_xyz[k] = -g_0_xyz[k] - g_z_0_0_xyz[k] * cd_z[k] + g_z_0_0_xyzz[k];

            g_z_0_z_xzz[k] = -g_0_xzz[k] - g_z_0_0_xzz[k] * cd_z[k] + g_z_0_0_xzzz[k];

            g_z_0_z_yyy[k] = -g_0_yyy[k] - g_z_0_0_yyy[k] * cd_z[k] + g_z_0_0_yyyz[k];

            g_z_0_z_yyz[k] = -g_0_yyz[k] - g_z_0_0_yyz[k] * cd_z[k] + g_z_0_0_yyzz[k];

            g_z_0_z_yzz[k] = -g_0_yzz[k] - g_z_0_0_yzz[k] * cd_z[k] + g_z_0_0_yzzz[k];

            g_z_0_z_zzz[k] = -g_0_zzz[k] - g_z_0_0_zzz[k] * cd_z[k] + g_z_0_0_zzzz[k];
        }
    }
}

} // t3ceri namespace


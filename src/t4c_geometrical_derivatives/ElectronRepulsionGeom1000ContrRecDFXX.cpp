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

#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_dfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_dfxx,
                                            const size_t idx_pfxx,
                                            const size_t idx_geom_10_pfxx,
                                            const size_t idx_geom_10_pgxx,
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
            /// Set up components of auxilary buffer : PFSS

            const auto pf_off = idx_pfxx + i * dcomps + j;

            auto g_x_xxx = cbuffer.data(pf_off + 0 * ccomps * dcomps);

            auto g_x_xxy = cbuffer.data(pf_off + 1 * ccomps * dcomps);

            auto g_x_xxz = cbuffer.data(pf_off + 2 * ccomps * dcomps);

            auto g_x_xyy = cbuffer.data(pf_off + 3 * ccomps * dcomps);

            auto g_x_xyz = cbuffer.data(pf_off + 4 * ccomps * dcomps);

            auto g_x_xzz = cbuffer.data(pf_off + 5 * ccomps * dcomps);

            auto g_x_yyy = cbuffer.data(pf_off + 6 * ccomps * dcomps);

            auto g_x_yyz = cbuffer.data(pf_off + 7 * ccomps * dcomps);

            auto g_x_yzz = cbuffer.data(pf_off + 8 * ccomps * dcomps);

            auto g_x_zzz = cbuffer.data(pf_off + 9 * ccomps * dcomps);

            auto g_y_xxx = cbuffer.data(pf_off + 10 * ccomps * dcomps);

            auto g_y_xxy = cbuffer.data(pf_off + 11 * ccomps * dcomps);

            auto g_y_xxz = cbuffer.data(pf_off + 12 * ccomps * dcomps);

            auto g_y_xyy = cbuffer.data(pf_off + 13 * ccomps * dcomps);

            auto g_y_xyz = cbuffer.data(pf_off + 14 * ccomps * dcomps);

            auto g_y_xzz = cbuffer.data(pf_off + 15 * ccomps * dcomps);

            auto g_y_yyy = cbuffer.data(pf_off + 16 * ccomps * dcomps);

            auto g_y_yyz = cbuffer.data(pf_off + 17 * ccomps * dcomps);

            auto g_y_yzz = cbuffer.data(pf_off + 18 * ccomps * dcomps);

            auto g_y_zzz = cbuffer.data(pf_off + 19 * ccomps * dcomps);

            auto g_z_xxx = cbuffer.data(pf_off + 20 * ccomps * dcomps);

            auto g_z_xxy = cbuffer.data(pf_off + 21 * ccomps * dcomps);

            auto g_z_xxz = cbuffer.data(pf_off + 22 * ccomps * dcomps);

            auto g_z_xyy = cbuffer.data(pf_off + 23 * ccomps * dcomps);

            auto g_z_xyz = cbuffer.data(pf_off + 24 * ccomps * dcomps);

            auto g_z_xzz = cbuffer.data(pf_off + 25 * ccomps * dcomps);

            auto g_z_yyy = cbuffer.data(pf_off + 26 * ccomps * dcomps);

            auto g_z_yyz = cbuffer.data(pf_off + 27 * ccomps * dcomps);

            auto g_z_yzz = cbuffer.data(pf_off + 28 * ccomps * dcomps);

            auto g_z_zzz = cbuffer.data(pf_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PFSS

            const auto pf_geom_10_off = idx_geom_10_pfxx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : PGSS

            const auto pg_geom_10_off = idx_geom_10_pgxx + i * dcomps + j;

            auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 89 * ccomps * dcomps);

            auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 95 * ccomps * dcomps);

            auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 107 * ccomps * dcomps);

            auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 108 * ccomps * dcomps);

            auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 109 * ccomps * dcomps);

            auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 110 * ccomps * dcomps);

            auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 111 * ccomps * dcomps);

            auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 112 * ccomps * dcomps);

            auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 113 * ccomps * dcomps);

            auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 114 * ccomps * dcomps);

            auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 115 * ccomps * dcomps);

            auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 116 * ccomps * dcomps);

            auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 117 * ccomps * dcomps);

            auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 118 * ccomps * dcomps);

            auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 134 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dfxx

            const auto df_geom_10_off = idx_geom_10_dfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_xxx = cbuffer.data(df_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxy = cbuffer.data(df_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxz = cbuffer.data(df_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xyy = cbuffer.data(df_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xyz = cbuffer.data(df_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xzz = cbuffer.data(df_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_yyy = cbuffer.data(df_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_yyz = cbuffer.data(df_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_yzz = cbuffer.data(df_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_zzz = cbuffer.data(df_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxx, g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxy, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xxzz, g_x_0_x_xyy, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_xzzz, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yzz, g_x_0_x_zzz, g_x_0_xx_xxx, g_x_0_xx_xxy, g_x_0_xx_xxz, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xzz, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yzz, g_x_0_xx_zzz, g_x_xxx, g_x_xxy, g_x_xxz, g_x_xyy, g_x_xyz, g_x_xzz, g_x_yyy, g_x_yyz, g_x_yzz, g_x_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_xxx[k] = -g_x_xxx[k] - g_x_0_x_xxx[k] * ab_x + g_x_0_x_xxxx[k];

                g_x_0_xx_xxy[k] = -g_x_xxy[k] - g_x_0_x_xxy[k] * ab_x + g_x_0_x_xxxy[k];

                g_x_0_xx_xxz[k] = -g_x_xxz[k] - g_x_0_x_xxz[k] * ab_x + g_x_0_x_xxxz[k];

                g_x_0_xx_xyy[k] = -g_x_xyy[k] - g_x_0_x_xyy[k] * ab_x + g_x_0_x_xxyy[k];

                g_x_0_xx_xyz[k] = -g_x_xyz[k] - g_x_0_x_xyz[k] * ab_x + g_x_0_x_xxyz[k];

                g_x_0_xx_xzz[k] = -g_x_xzz[k] - g_x_0_x_xzz[k] * ab_x + g_x_0_x_xxzz[k];

                g_x_0_xx_yyy[k] = -g_x_yyy[k] - g_x_0_x_yyy[k] * ab_x + g_x_0_x_xyyy[k];

                g_x_0_xx_yyz[k] = -g_x_yyz[k] - g_x_0_x_yyz[k] * ab_x + g_x_0_x_xyyz[k];

                g_x_0_xx_yzz[k] = -g_x_yzz[k] - g_x_0_x_yzz[k] * ab_x + g_x_0_x_xyzz[k];

                g_x_0_xx_zzz[k] = -g_x_zzz[k] - g_x_0_x_zzz[k] * ab_x + g_x_0_x_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_xxx = cbuffer.data(df_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xy_xxy = cbuffer.data(df_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xy_xxz = cbuffer.data(df_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xy_xyy = cbuffer.data(df_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xy_xyz = cbuffer.data(df_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xy_xzz = cbuffer.data(df_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xy_yyy = cbuffer.data(df_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xy_yyz = cbuffer.data(df_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xy_yzz = cbuffer.data(df_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xy_zzz = cbuffer.data(df_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxx, g_x_0_x_xxxy, g_x_0_x_xxy, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xyy, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_yyy, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyz, g_x_0_x_yyzz, g_x_0_x_yzz, g_x_0_x_yzzz, g_x_0_x_zzz, g_x_0_xy_xxx, g_x_0_xy_xxy, g_x_0_xy_xxz, g_x_0_xy_xyy, g_x_0_xy_xyz, g_x_0_xy_xzz, g_x_0_xy_yyy, g_x_0_xy_yyz, g_x_0_xy_yzz, g_x_0_xy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_xxx[k] = -g_x_0_x_xxx[k] * ab_y + g_x_0_x_xxxy[k];

                g_x_0_xy_xxy[k] = -g_x_0_x_xxy[k] * ab_y + g_x_0_x_xxyy[k];

                g_x_0_xy_xxz[k] = -g_x_0_x_xxz[k] * ab_y + g_x_0_x_xxyz[k];

                g_x_0_xy_xyy[k] = -g_x_0_x_xyy[k] * ab_y + g_x_0_x_xyyy[k];

                g_x_0_xy_xyz[k] = -g_x_0_x_xyz[k] * ab_y + g_x_0_x_xyyz[k];

                g_x_0_xy_xzz[k] = -g_x_0_x_xzz[k] * ab_y + g_x_0_x_xyzz[k];

                g_x_0_xy_yyy[k] = -g_x_0_x_yyy[k] * ab_y + g_x_0_x_yyyy[k];

                g_x_0_xy_yyz[k] = -g_x_0_x_yyz[k] * ab_y + g_x_0_x_yyyz[k];

                g_x_0_xy_yzz[k] = -g_x_0_x_yzz[k] * ab_y + g_x_0_x_yyzz[k];

                g_x_0_xy_zzz[k] = -g_x_0_x_zzz[k] * ab_y + g_x_0_x_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_xxx = cbuffer.data(df_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xz_xxy = cbuffer.data(df_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xz_xxz = cbuffer.data(df_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xz_xyy = cbuffer.data(df_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xz_xyz = cbuffer.data(df_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xz_xzz = cbuffer.data(df_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xz_yyy = cbuffer.data(df_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xz_yyz = cbuffer.data(df_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xz_yzz = cbuffer.data(df_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xz_zzz = cbuffer.data(df_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxx, g_x_0_x_xxxz, g_x_0_x_xxy, g_x_0_x_xxyz, g_x_0_x_xxz, g_x_0_x_xxzz, g_x_0_x_xyy, g_x_0_x_xyyz, g_x_0_x_xyz, g_x_0_x_xyzz, g_x_0_x_xzz, g_x_0_x_xzzz, g_x_0_x_yyy, g_x_0_x_yyyz, g_x_0_x_yyz, g_x_0_x_yyzz, g_x_0_x_yzz, g_x_0_x_yzzz, g_x_0_x_zzz, g_x_0_x_zzzz, g_x_0_xz_xxx, g_x_0_xz_xxy, g_x_0_xz_xxz, g_x_0_xz_xyy, g_x_0_xz_xyz, g_x_0_xz_xzz, g_x_0_xz_yyy, g_x_0_xz_yyz, g_x_0_xz_yzz, g_x_0_xz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_xxx[k] = -g_x_0_x_xxx[k] * ab_z + g_x_0_x_xxxz[k];

                g_x_0_xz_xxy[k] = -g_x_0_x_xxy[k] * ab_z + g_x_0_x_xxyz[k];

                g_x_0_xz_xxz[k] = -g_x_0_x_xxz[k] * ab_z + g_x_0_x_xxzz[k];

                g_x_0_xz_xyy[k] = -g_x_0_x_xyy[k] * ab_z + g_x_0_x_xyyz[k];

                g_x_0_xz_xyz[k] = -g_x_0_x_xyz[k] * ab_z + g_x_0_x_xyzz[k];

                g_x_0_xz_xzz[k] = -g_x_0_x_xzz[k] * ab_z + g_x_0_x_xzzz[k];

                g_x_0_xz_yyy[k] = -g_x_0_x_yyy[k] * ab_z + g_x_0_x_yyyz[k];

                g_x_0_xz_yyz[k] = -g_x_0_x_yyz[k] * ab_z + g_x_0_x_yyzz[k];

                g_x_0_xz_yzz[k] = -g_x_0_x_yzz[k] * ab_z + g_x_0_x_yzzz[k];

                g_x_0_xz_zzz[k] = -g_x_0_x_zzz[k] * ab_z + g_x_0_x_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_xxx = cbuffer.data(df_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_yy_xxy = cbuffer.data(df_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_yy_xxz = cbuffer.data(df_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_yy_xyy = cbuffer.data(df_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_yy_xyz = cbuffer.data(df_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_yy_xzz = cbuffer.data(df_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_yy_yyy = cbuffer.data(df_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_yy_yyz = cbuffer.data(df_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_yy_yzz = cbuffer.data(df_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_yy_zzz = cbuffer.data(df_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxx, g_x_0_y_xxxy, g_x_0_y_xxy, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxz, g_x_0_y_xyy, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyz, g_x_0_y_xyzz, g_x_0_y_xzz, g_x_0_y_yyy, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyz, g_x_0_y_yyzz, g_x_0_y_yzz, g_x_0_y_yzzz, g_x_0_y_zzz, g_x_0_yy_xxx, g_x_0_yy_xxy, g_x_0_yy_xxz, g_x_0_yy_xyy, g_x_0_yy_xyz, g_x_0_yy_xzz, g_x_0_yy_yyy, g_x_0_yy_yyz, g_x_0_yy_yzz, g_x_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_xxx[k] = -g_x_0_y_xxx[k] * ab_y + g_x_0_y_xxxy[k];

                g_x_0_yy_xxy[k] = -g_x_0_y_xxy[k] * ab_y + g_x_0_y_xxyy[k];

                g_x_0_yy_xxz[k] = -g_x_0_y_xxz[k] * ab_y + g_x_0_y_xxyz[k];

                g_x_0_yy_xyy[k] = -g_x_0_y_xyy[k] * ab_y + g_x_0_y_xyyy[k];

                g_x_0_yy_xyz[k] = -g_x_0_y_xyz[k] * ab_y + g_x_0_y_xyyz[k];

                g_x_0_yy_xzz[k] = -g_x_0_y_xzz[k] * ab_y + g_x_0_y_xyzz[k];

                g_x_0_yy_yyy[k] = -g_x_0_y_yyy[k] * ab_y + g_x_0_y_yyyy[k];

                g_x_0_yy_yyz[k] = -g_x_0_y_yyz[k] * ab_y + g_x_0_y_yyyz[k];

                g_x_0_yy_yzz[k] = -g_x_0_y_yzz[k] * ab_y + g_x_0_y_yyzz[k];

                g_x_0_yy_zzz[k] = -g_x_0_y_zzz[k] * ab_y + g_x_0_y_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_xxx = cbuffer.data(df_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_yz_xxy = cbuffer.data(df_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_yz_xxz = cbuffer.data(df_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_yz_xyy = cbuffer.data(df_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_yz_xyz = cbuffer.data(df_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yz_xzz = cbuffer.data(df_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yz_yyy = cbuffer.data(df_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yz_yyz = cbuffer.data(df_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yz_yzz = cbuffer.data(df_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yz_zzz = cbuffer.data(df_geom_10_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yz_xxx, g_x_0_yz_xxy, g_x_0_yz_xxz, g_x_0_yz_xyy, g_x_0_yz_xyz, g_x_0_yz_xzz, g_x_0_yz_yyy, g_x_0_yz_yyz, g_x_0_yz_yzz, g_x_0_yz_zzz, g_x_0_z_xxx, g_x_0_z_xxxy, g_x_0_z_xxy, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxz, g_x_0_z_xyy, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyz, g_x_0_z_xyzz, g_x_0_z_xzz, g_x_0_z_yyy, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyz, g_x_0_z_yyzz, g_x_0_z_yzz, g_x_0_z_yzzz, g_x_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_xxx[k] = -g_x_0_z_xxx[k] * ab_y + g_x_0_z_xxxy[k];

                g_x_0_yz_xxy[k] = -g_x_0_z_xxy[k] * ab_y + g_x_0_z_xxyy[k];

                g_x_0_yz_xxz[k] = -g_x_0_z_xxz[k] * ab_y + g_x_0_z_xxyz[k];

                g_x_0_yz_xyy[k] = -g_x_0_z_xyy[k] * ab_y + g_x_0_z_xyyy[k];

                g_x_0_yz_xyz[k] = -g_x_0_z_xyz[k] * ab_y + g_x_0_z_xyyz[k];

                g_x_0_yz_xzz[k] = -g_x_0_z_xzz[k] * ab_y + g_x_0_z_xyzz[k];

                g_x_0_yz_yyy[k] = -g_x_0_z_yyy[k] * ab_y + g_x_0_z_yyyy[k];

                g_x_0_yz_yyz[k] = -g_x_0_z_yyz[k] * ab_y + g_x_0_z_yyyz[k];

                g_x_0_yz_yzz[k] = -g_x_0_z_yzz[k] * ab_y + g_x_0_z_yyzz[k];

                g_x_0_yz_zzz[k] = -g_x_0_z_zzz[k] * ab_y + g_x_0_z_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_xxx = cbuffer.data(df_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_zz_xxy = cbuffer.data(df_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_zz_xxz = cbuffer.data(df_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_zz_xyy = cbuffer.data(df_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_zz_xyz = cbuffer.data(df_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_zz_xzz = cbuffer.data(df_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_zz_yyy = cbuffer.data(df_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_zz_yyz = cbuffer.data(df_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_zz_yzz = cbuffer.data(df_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_zz_zzz = cbuffer.data(df_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxx, g_x_0_z_xxxz, g_x_0_z_xxy, g_x_0_z_xxyz, g_x_0_z_xxz, g_x_0_z_xxzz, g_x_0_z_xyy, g_x_0_z_xyyz, g_x_0_z_xyz, g_x_0_z_xyzz, g_x_0_z_xzz, g_x_0_z_xzzz, g_x_0_z_yyy, g_x_0_z_yyyz, g_x_0_z_yyz, g_x_0_z_yyzz, g_x_0_z_yzz, g_x_0_z_yzzz, g_x_0_z_zzz, g_x_0_z_zzzz, g_x_0_zz_xxx, g_x_0_zz_xxy, g_x_0_zz_xxz, g_x_0_zz_xyy, g_x_0_zz_xyz, g_x_0_zz_xzz, g_x_0_zz_yyy, g_x_0_zz_yyz, g_x_0_zz_yzz, g_x_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_xxx[k] = -g_x_0_z_xxx[k] * ab_z + g_x_0_z_xxxz[k];

                g_x_0_zz_xxy[k] = -g_x_0_z_xxy[k] * ab_z + g_x_0_z_xxyz[k];

                g_x_0_zz_xxz[k] = -g_x_0_z_xxz[k] * ab_z + g_x_0_z_xxzz[k];

                g_x_0_zz_xyy[k] = -g_x_0_z_xyy[k] * ab_z + g_x_0_z_xyyz[k];

                g_x_0_zz_xyz[k] = -g_x_0_z_xyz[k] * ab_z + g_x_0_z_xyzz[k];

                g_x_0_zz_xzz[k] = -g_x_0_z_xzz[k] * ab_z + g_x_0_z_xzzz[k];

                g_x_0_zz_yyy[k] = -g_x_0_z_yyy[k] * ab_z + g_x_0_z_yyyz[k];

                g_x_0_zz_yyz[k] = -g_x_0_z_yyz[k] * ab_z + g_x_0_z_yyzz[k];

                g_x_0_zz_yzz[k] = -g_x_0_z_yzz[k] * ab_z + g_x_0_z_yzzz[k];

                g_x_0_zz_zzz[k] = -g_x_0_z_zzz[k] * ab_z + g_x_0_z_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_xxx = cbuffer.data(df_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_xx_xxy = cbuffer.data(df_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_xx_xxz = cbuffer.data(df_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xx_xyy = cbuffer.data(df_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xx_xyz = cbuffer.data(df_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xx_xzz = cbuffer.data(df_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xx_yyy = cbuffer.data(df_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xx_yyz = cbuffer.data(df_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xx_yzz = cbuffer.data(df_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xx_zzz = cbuffer.data(df_geom_10_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxx, g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxy, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxz, g_y_0_x_xxzz, g_y_0_x_xyy, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyz, g_y_0_x_xyzz, g_y_0_x_xzz, g_y_0_x_xzzz, g_y_0_x_yyy, g_y_0_x_yyz, g_y_0_x_yzz, g_y_0_x_zzz, g_y_0_xx_xxx, g_y_0_xx_xxy, g_y_0_xx_xxz, g_y_0_xx_xyy, g_y_0_xx_xyz, g_y_0_xx_xzz, g_y_0_xx_yyy, g_y_0_xx_yyz, g_y_0_xx_yzz, g_y_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_xxx[k] = -g_y_0_x_xxx[k] * ab_x + g_y_0_x_xxxx[k];

                g_y_0_xx_xxy[k] = -g_y_0_x_xxy[k] * ab_x + g_y_0_x_xxxy[k];

                g_y_0_xx_xxz[k] = -g_y_0_x_xxz[k] * ab_x + g_y_0_x_xxxz[k];

                g_y_0_xx_xyy[k] = -g_y_0_x_xyy[k] * ab_x + g_y_0_x_xxyy[k];

                g_y_0_xx_xyz[k] = -g_y_0_x_xyz[k] * ab_x + g_y_0_x_xxyz[k];

                g_y_0_xx_xzz[k] = -g_y_0_x_xzz[k] * ab_x + g_y_0_x_xxzz[k];

                g_y_0_xx_yyy[k] = -g_y_0_x_yyy[k] * ab_x + g_y_0_x_xyyy[k];

                g_y_0_xx_yyz[k] = -g_y_0_x_yyz[k] * ab_x + g_y_0_x_xyyz[k];

                g_y_0_xx_yzz[k] = -g_y_0_x_yzz[k] * ab_x + g_y_0_x_xyzz[k];

                g_y_0_xx_zzz[k] = -g_y_0_x_zzz[k] * ab_x + g_y_0_x_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_xxx = cbuffer.data(df_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xy_xxy = cbuffer.data(df_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xy_xxz = cbuffer.data(df_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xy_xyy = cbuffer.data(df_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xy_xyz = cbuffer.data(df_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xy_xzz = cbuffer.data(df_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xy_yyy = cbuffer.data(df_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xy_yyz = cbuffer.data(df_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_xy_yzz = cbuffer.data(df_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xy_zzz = cbuffer.data(df_geom_10_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xxx, g_y_0_xy_xxy, g_y_0_xy_xxz, g_y_0_xy_xyy, g_y_0_xy_xyz, g_y_0_xy_xzz, g_y_0_xy_yyy, g_y_0_xy_yyz, g_y_0_xy_yzz, g_y_0_xy_zzz, g_y_0_y_xxx, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxy, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xxzz, g_y_0_y_xyy, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_xzzz, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yzz, g_y_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_xxx[k] = -g_y_0_y_xxx[k] * ab_x + g_y_0_y_xxxx[k];

                g_y_0_xy_xxy[k] = -g_y_0_y_xxy[k] * ab_x + g_y_0_y_xxxy[k];

                g_y_0_xy_xxz[k] = -g_y_0_y_xxz[k] * ab_x + g_y_0_y_xxxz[k];

                g_y_0_xy_xyy[k] = -g_y_0_y_xyy[k] * ab_x + g_y_0_y_xxyy[k];

                g_y_0_xy_xyz[k] = -g_y_0_y_xyz[k] * ab_x + g_y_0_y_xxyz[k];

                g_y_0_xy_xzz[k] = -g_y_0_y_xzz[k] * ab_x + g_y_0_y_xxzz[k];

                g_y_0_xy_yyy[k] = -g_y_0_y_yyy[k] * ab_x + g_y_0_y_xyyy[k];

                g_y_0_xy_yyz[k] = -g_y_0_y_yyz[k] * ab_x + g_y_0_y_xyyz[k];

                g_y_0_xy_yzz[k] = -g_y_0_y_yzz[k] * ab_x + g_y_0_y_xyzz[k];

                g_y_0_xy_zzz[k] = -g_y_0_y_zzz[k] * ab_x + g_y_0_y_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_xxx = cbuffer.data(df_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xz_xxy = cbuffer.data(df_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xz_xxz = cbuffer.data(df_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xz_xyy = cbuffer.data(df_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_xz_xyz = cbuffer.data(df_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xz_xzz = cbuffer.data(df_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xz_yyy = cbuffer.data(df_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xz_yyz = cbuffer.data(df_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xz_yzz = cbuffer.data(df_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xz_zzz = cbuffer.data(df_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xz_xxx, g_y_0_xz_xxy, g_y_0_xz_xxz, g_y_0_xz_xyy, g_y_0_xz_xyz, g_y_0_xz_xzz, g_y_0_xz_yyy, g_y_0_xz_yyz, g_y_0_xz_yzz, g_y_0_xz_zzz, g_y_0_z_xxx, g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxy, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxz, g_y_0_z_xxzz, g_y_0_z_xyy, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyz, g_y_0_z_xyzz, g_y_0_z_xzz, g_y_0_z_xzzz, g_y_0_z_yyy, g_y_0_z_yyz, g_y_0_z_yzz, g_y_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_xxx[k] = -g_y_0_z_xxx[k] * ab_x + g_y_0_z_xxxx[k];

                g_y_0_xz_xxy[k] = -g_y_0_z_xxy[k] * ab_x + g_y_0_z_xxxy[k];

                g_y_0_xz_xxz[k] = -g_y_0_z_xxz[k] * ab_x + g_y_0_z_xxxz[k];

                g_y_0_xz_xyy[k] = -g_y_0_z_xyy[k] * ab_x + g_y_0_z_xxyy[k];

                g_y_0_xz_xyz[k] = -g_y_0_z_xyz[k] * ab_x + g_y_0_z_xxyz[k];

                g_y_0_xz_xzz[k] = -g_y_0_z_xzz[k] * ab_x + g_y_0_z_xxzz[k];

                g_y_0_xz_yyy[k] = -g_y_0_z_yyy[k] * ab_x + g_y_0_z_xyyy[k];

                g_y_0_xz_yyz[k] = -g_y_0_z_yyz[k] * ab_x + g_y_0_z_xyyz[k];

                g_y_0_xz_yzz[k] = -g_y_0_z_yzz[k] * ab_x + g_y_0_z_xyzz[k];

                g_y_0_xz_zzz[k] = -g_y_0_z_zzz[k] * ab_x + g_y_0_z_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_xxx = cbuffer.data(df_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_yy_xxy = cbuffer.data(df_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_yy_xxz = cbuffer.data(df_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_yy_xyy = cbuffer.data(df_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_yy_xyz = cbuffer.data(df_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_yy_xzz = cbuffer.data(df_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_yy_yyy = cbuffer.data(df_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_yy_yyz = cbuffer.data(df_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_yy_yzz = cbuffer.data(df_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_yy_zzz = cbuffer.data(df_geom_10_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxx, g_y_0_y_xxxy, g_y_0_y_xxy, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xyy, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_yyy, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyz, g_y_0_y_yyzz, g_y_0_y_yzz, g_y_0_y_yzzz, g_y_0_y_zzz, g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xzz, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yzz, g_y_0_yy_zzz, g_y_xxx, g_y_xxy, g_y_xxz, g_y_xyy, g_y_xyz, g_y_xzz, g_y_yyy, g_y_yyz, g_y_yzz, g_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_xxx[k] = -g_y_xxx[k] - g_y_0_y_xxx[k] * ab_y + g_y_0_y_xxxy[k];

                g_y_0_yy_xxy[k] = -g_y_xxy[k] - g_y_0_y_xxy[k] * ab_y + g_y_0_y_xxyy[k];

                g_y_0_yy_xxz[k] = -g_y_xxz[k] - g_y_0_y_xxz[k] * ab_y + g_y_0_y_xxyz[k];

                g_y_0_yy_xyy[k] = -g_y_xyy[k] - g_y_0_y_xyy[k] * ab_y + g_y_0_y_xyyy[k];

                g_y_0_yy_xyz[k] = -g_y_xyz[k] - g_y_0_y_xyz[k] * ab_y + g_y_0_y_xyyz[k];

                g_y_0_yy_xzz[k] = -g_y_xzz[k] - g_y_0_y_xzz[k] * ab_y + g_y_0_y_xyzz[k];

                g_y_0_yy_yyy[k] = -g_y_yyy[k] - g_y_0_y_yyy[k] * ab_y + g_y_0_y_yyyy[k];

                g_y_0_yy_yyz[k] = -g_y_yyz[k] - g_y_0_y_yyz[k] * ab_y + g_y_0_y_yyyz[k];

                g_y_0_yy_yzz[k] = -g_y_yzz[k] - g_y_0_y_yzz[k] * ab_y + g_y_0_y_yyzz[k];

                g_y_0_yy_zzz[k] = -g_y_zzz[k] - g_y_0_y_zzz[k] * ab_y + g_y_0_y_yzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_xxx = cbuffer.data(df_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_yz_xxy = cbuffer.data(df_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_yz_xxz = cbuffer.data(df_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_yz_xyy = cbuffer.data(df_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_yz_xyz = cbuffer.data(df_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_yz_xzz = cbuffer.data(df_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_yz_yyy = cbuffer.data(df_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_yz_yyz = cbuffer.data(df_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yz_yzz = cbuffer.data(df_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yz_zzz = cbuffer.data(df_geom_10_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxx, g_y_0_y_xxxz, g_y_0_y_xxy, g_y_0_y_xxyz, g_y_0_y_xxz, g_y_0_y_xxzz, g_y_0_y_xyy, g_y_0_y_xyyz, g_y_0_y_xyz, g_y_0_y_xyzz, g_y_0_y_xzz, g_y_0_y_xzzz, g_y_0_y_yyy, g_y_0_y_yyyz, g_y_0_y_yyz, g_y_0_y_yyzz, g_y_0_y_yzz, g_y_0_y_yzzz, g_y_0_y_zzz, g_y_0_y_zzzz, g_y_0_yz_xxx, g_y_0_yz_xxy, g_y_0_yz_xxz, g_y_0_yz_xyy, g_y_0_yz_xyz, g_y_0_yz_xzz, g_y_0_yz_yyy, g_y_0_yz_yyz, g_y_0_yz_yzz, g_y_0_yz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_xxx[k] = -g_y_0_y_xxx[k] * ab_z + g_y_0_y_xxxz[k];

                g_y_0_yz_xxy[k] = -g_y_0_y_xxy[k] * ab_z + g_y_0_y_xxyz[k];

                g_y_0_yz_xxz[k] = -g_y_0_y_xxz[k] * ab_z + g_y_0_y_xxzz[k];

                g_y_0_yz_xyy[k] = -g_y_0_y_xyy[k] * ab_z + g_y_0_y_xyyz[k];

                g_y_0_yz_xyz[k] = -g_y_0_y_xyz[k] * ab_z + g_y_0_y_xyzz[k];

                g_y_0_yz_xzz[k] = -g_y_0_y_xzz[k] * ab_z + g_y_0_y_xzzz[k];

                g_y_0_yz_yyy[k] = -g_y_0_y_yyy[k] * ab_z + g_y_0_y_yyyz[k];

                g_y_0_yz_yyz[k] = -g_y_0_y_yyz[k] * ab_z + g_y_0_y_yyzz[k];

                g_y_0_yz_yzz[k] = -g_y_0_y_yzz[k] * ab_z + g_y_0_y_yzzz[k];

                g_y_0_yz_zzz[k] = -g_y_0_y_zzz[k] * ab_z + g_y_0_y_zzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_xxx = cbuffer.data(df_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_zz_xxy = cbuffer.data(df_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_zz_xxz = cbuffer.data(df_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_zz_xyy = cbuffer.data(df_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_zz_xyz = cbuffer.data(df_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_zz_xzz = cbuffer.data(df_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_zz_yyy = cbuffer.data(df_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_zz_yyz = cbuffer.data(df_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_zz_yzz = cbuffer.data(df_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_zz_zzz = cbuffer.data(df_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxx, g_y_0_z_xxxz, g_y_0_z_xxy, g_y_0_z_xxyz, g_y_0_z_xxz, g_y_0_z_xxzz, g_y_0_z_xyy, g_y_0_z_xyyz, g_y_0_z_xyz, g_y_0_z_xyzz, g_y_0_z_xzz, g_y_0_z_xzzz, g_y_0_z_yyy, g_y_0_z_yyyz, g_y_0_z_yyz, g_y_0_z_yyzz, g_y_0_z_yzz, g_y_0_z_yzzz, g_y_0_z_zzz, g_y_0_z_zzzz, g_y_0_zz_xxx, g_y_0_zz_xxy, g_y_0_zz_xxz, g_y_0_zz_xyy, g_y_0_zz_xyz, g_y_0_zz_xzz, g_y_0_zz_yyy, g_y_0_zz_yyz, g_y_0_zz_yzz, g_y_0_zz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_xxx[k] = -g_y_0_z_xxx[k] * ab_z + g_y_0_z_xxxz[k];

                g_y_0_zz_xxy[k] = -g_y_0_z_xxy[k] * ab_z + g_y_0_z_xxyz[k];

                g_y_0_zz_xxz[k] = -g_y_0_z_xxz[k] * ab_z + g_y_0_z_xxzz[k];

                g_y_0_zz_xyy[k] = -g_y_0_z_xyy[k] * ab_z + g_y_0_z_xyyz[k];

                g_y_0_zz_xyz[k] = -g_y_0_z_xyz[k] * ab_z + g_y_0_z_xyzz[k];

                g_y_0_zz_xzz[k] = -g_y_0_z_xzz[k] * ab_z + g_y_0_z_xzzz[k];

                g_y_0_zz_yyy[k] = -g_y_0_z_yyy[k] * ab_z + g_y_0_z_yyyz[k];

                g_y_0_zz_yyz[k] = -g_y_0_z_yyz[k] * ab_z + g_y_0_z_yyzz[k];

                g_y_0_zz_yzz[k] = -g_y_0_z_yzz[k] * ab_z + g_y_0_z_yzzz[k];

                g_y_0_zz_zzz[k] = -g_y_0_z_zzz[k] * ab_z + g_y_0_z_zzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_xxx = cbuffer.data(df_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_xx_xxy = cbuffer.data(df_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_xx_xxz = cbuffer.data(df_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_xx_xyy = cbuffer.data(df_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_xx_xyz = cbuffer.data(df_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_xx_xzz = cbuffer.data(df_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_xx_yyy = cbuffer.data(df_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xx_yyz = cbuffer.data(df_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xx_yzz = cbuffer.data(df_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xx_zzz = cbuffer.data(df_geom_10_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxx, g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxy, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxz, g_z_0_x_xxzz, g_z_0_x_xyy, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyz, g_z_0_x_xyzz, g_z_0_x_xzz, g_z_0_x_xzzz, g_z_0_x_yyy, g_z_0_x_yyz, g_z_0_x_yzz, g_z_0_x_zzz, g_z_0_xx_xxx, g_z_0_xx_xxy, g_z_0_xx_xxz, g_z_0_xx_xyy, g_z_0_xx_xyz, g_z_0_xx_xzz, g_z_0_xx_yyy, g_z_0_xx_yyz, g_z_0_xx_yzz, g_z_0_xx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_xxx[k] = -g_z_0_x_xxx[k] * ab_x + g_z_0_x_xxxx[k];

                g_z_0_xx_xxy[k] = -g_z_0_x_xxy[k] * ab_x + g_z_0_x_xxxy[k];

                g_z_0_xx_xxz[k] = -g_z_0_x_xxz[k] * ab_x + g_z_0_x_xxxz[k];

                g_z_0_xx_xyy[k] = -g_z_0_x_xyy[k] * ab_x + g_z_0_x_xxyy[k];

                g_z_0_xx_xyz[k] = -g_z_0_x_xyz[k] * ab_x + g_z_0_x_xxyz[k];

                g_z_0_xx_xzz[k] = -g_z_0_x_xzz[k] * ab_x + g_z_0_x_xxzz[k];

                g_z_0_xx_yyy[k] = -g_z_0_x_yyy[k] * ab_x + g_z_0_x_xyyy[k];

                g_z_0_xx_yyz[k] = -g_z_0_x_yyz[k] * ab_x + g_z_0_x_xyyz[k];

                g_z_0_xx_yzz[k] = -g_z_0_x_yzz[k] * ab_x + g_z_0_x_xyzz[k];

                g_z_0_xx_zzz[k] = -g_z_0_x_zzz[k] * ab_x + g_z_0_x_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_xxx = cbuffer.data(df_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xy_xxy = cbuffer.data(df_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_xy_xxz = cbuffer.data(df_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xy_xyy = cbuffer.data(df_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xy_xyz = cbuffer.data(df_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xy_xzz = cbuffer.data(df_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xy_yyy = cbuffer.data(df_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xy_yyz = cbuffer.data(df_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_xy_yzz = cbuffer.data(df_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xy_zzz = cbuffer.data(df_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xy_xxx, g_z_0_xy_xxy, g_z_0_xy_xxz, g_z_0_xy_xyy, g_z_0_xy_xyz, g_z_0_xy_xzz, g_z_0_xy_yyy, g_z_0_xy_yyz, g_z_0_xy_yzz, g_z_0_xy_zzz, g_z_0_y_xxx, g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxy, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxz, g_z_0_y_xxzz, g_z_0_y_xyy, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyz, g_z_0_y_xyzz, g_z_0_y_xzz, g_z_0_y_xzzz, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yzz, g_z_0_y_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_xxx[k] = -g_z_0_y_xxx[k] * ab_x + g_z_0_y_xxxx[k];

                g_z_0_xy_xxy[k] = -g_z_0_y_xxy[k] * ab_x + g_z_0_y_xxxy[k];

                g_z_0_xy_xxz[k] = -g_z_0_y_xxz[k] * ab_x + g_z_0_y_xxxz[k];

                g_z_0_xy_xyy[k] = -g_z_0_y_xyy[k] * ab_x + g_z_0_y_xxyy[k];

                g_z_0_xy_xyz[k] = -g_z_0_y_xyz[k] * ab_x + g_z_0_y_xxyz[k];

                g_z_0_xy_xzz[k] = -g_z_0_y_xzz[k] * ab_x + g_z_0_y_xxzz[k];

                g_z_0_xy_yyy[k] = -g_z_0_y_yyy[k] * ab_x + g_z_0_y_xyyy[k];

                g_z_0_xy_yyz[k] = -g_z_0_y_yyz[k] * ab_x + g_z_0_y_xyyz[k];

                g_z_0_xy_yzz[k] = -g_z_0_y_yzz[k] * ab_x + g_z_0_y_xyzz[k];

                g_z_0_xy_zzz[k] = -g_z_0_y_zzz[k] * ab_x + g_z_0_y_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_xxx = cbuffer.data(df_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xz_xxy = cbuffer.data(df_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xz_xxz = cbuffer.data(df_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xz_xyy = cbuffer.data(df_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_xz_xyz = cbuffer.data(df_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xz_xzz = cbuffer.data(df_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xz_yyy = cbuffer.data(df_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xz_yyz = cbuffer.data(df_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xz_yzz = cbuffer.data(df_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xz_zzz = cbuffer.data(df_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xxx, g_z_0_xz_xxy, g_z_0_xz_xxz, g_z_0_xz_xyy, g_z_0_xz_xyz, g_z_0_xz_xzz, g_z_0_xz_yyy, g_z_0_xz_yyz, g_z_0_xz_yzz, g_z_0_xz_zzz, g_z_0_z_xxx, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxy, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xxzz, g_z_0_z_xyy, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_xzzz, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_xxx[k] = -g_z_0_z_xxx[k] * ab_x + g_z_0_z_xxxx[k];

                g_z_0_xz_xxy[k] = -g_z_0_z_xxy[k] * ab_x + g_z_0_z_xxxy[k];

                g_z_0_xz_xxz[k] = -g_z_0_z_xxz[k] * ab_x + g_z_0_z_xxxz[k];

                g_z_0_xz_xyy[k] = -g_z_0_z_xyy[k] * ab_x + g_z_0_z_xxyy[k];

                g_z_0_xz_xyz[k] = -g_z_0_z_xyz[k] * ab_x + g_z_0_z_xxyz[k];

                g_z_0_xz_xzz[k] = -g_z_0_z_xzz[k] * ab_x + g_z_0_z_xxzz[k];

                g_z_0_xz_yyy[k] = -g_z_0_z_yyy[k] * ab_x + g_z_0_z_xyyy[k];

                g_z_0_xz_yyz[k] = -g_z_0_z_yyz[k] * ab_x + g_z_0_z_xyyz[k];

                g_z_0_xz_yzz[k] = -g_z_0_z_yzz[k] * ab_x + g_z_0_z_xyzz[k];

                g_z_0_xz_zzz[k] = -g_z_0_z_zzz[k] * ab_x + g_z_0_z_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_xxx = cbuffer.data(df_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_yy_xxy = cbuffer.data(df_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_yy_xxz = cbuffer.data(df_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_yy_xyy = cbuffer.data(df_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_yy_xyz = cbuffer.data(df_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_yy_xzz = cbuffer.data(df_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_yy_yyy = cbuffer.data(df_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_yy_yyz = cbuffer.data(df_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_yy_yzz = cbuffer.data(df_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_yy_zzz = cbuffer.data(df_geom_10_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxx, g_z_0_y_xxxy, g_z_0_y_xxy, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxz, g_z_0_y_xyy, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyz, g_z_0_y_xyzz, g_z_0_y_xzz, g_z_0_y_yyy, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyz, g_z_0_y_yyzz, g_z_0_y_yzz, g_z_0_y_yzzz, g_z_0_y_zzz, g_z_0_yy_xxx, g_z_0_yy_xxy, g_z_0_yy_xxz, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xzz, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yzz, g_z_0_yy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_xxx[k] = -g_z_0_y_xxx[k] * ab_y + g_z_0_y_xxxy[k];

                g_z_0_yy_xxy[k] = -g_z_0_y_xxy[k] * ab_y + g_z_0_y_xxyy[k];

                g_z_0_yy_xxz[k] = -g_z_0_y_xxz[k] * ab_y + g_z_0_y_xxyz[k];

                g_z_0_yy_xyy[k] = -g_z_0_y_xyy[k] * ab_y + g_z_0_y_xyyy[k];

                g_z_0_yy_xyz[k] = -g_z_0_y_xyz[k] * ab_y + g_z_0_y_xyyz[k];

                g_z_0_yy_xzz[k] = -g_z_0_y_xzz[k] * ab_y + g_z_0_y_xyzz[k];

                g_z_0_yy_yyy[k] = -g_z_0_y_yyy[k] * ab_y + g_z_0_y_yyyy[k];

                g_z_0_yy_yyz[k] = -g_z_0_y_yyz[k] * ab_y + g_z_0_y_yyyz[k];

                g_z_0_yy_yzz[k] = -g_z_0_y_yzz[k] * ab_y + g_z_0_y_yyzz[k];

                g_z_0_yy_zzz[k] = -g_z_0_y_zzz[k] * ab_y + g_z_0_y_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_xxx = cbuffer.data(df_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_yz_xxy = cbuffer.data(df_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_yz_xxz = cbuffer.data(df_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_yz_xyy = cbuffer.data(df_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_yz_xyz = cbuffer.data(df_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_yz_xzz = cbuffer.data(df_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_yz_yyy = cbuffer.data(df_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_yz_yyz = cbuffer.data(df_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_yz_yzz = cbuffer.data(df_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_yz_zzz = cbuffer.data(df_geom_10_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xxx, g_z_0_yz_xxy, g_z_0_yz_xxz, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xzz, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yzz, g_z_0_yz_zzz, g_z_0_z_xxx, g_z_0_z_xxxy, g_z_0_z_xxy, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xyy, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_yyy, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyz, g_z_0_z_yyzz, g_z_0_z_yzz, g_z_0_z_yzzz, g_z_0_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_xxx[k] = -g_z_0_z_xxx[k] * ab_y + g_z_0_z_xxxy[k];

                g_z_0_yz_xxy[k] = -g_z_0_z_xxy[k] * ab_y + g_z_0_z_xxyy[k];

                g_z_0_yz_xxz[k] = -g_z_0_z_xxz[k] * ab_y + g_z_0_z_xxyz[k];

                g_z_0_yz_xyy[k] = -g_z_0_z_xyy[k] * ab_y + g_z_0_z_xyyy[k];

                g_z_0_yz_xyz[k] = -g_z_0_z_xyz[k] * ab_y + g_z_0_z_xyyz[k];

                g_z_0_yz_xzz[k] = -g_z_0_z_xzz[k] * ab_y + g_z_0_z_xyzz[k];

                g_z_0_yz_yyy[k] = -g_z_0_z_yyy[k] * ab_y + g_z_0_z_yyyy[k];

                g_z_0_yz_yyz[k] = -g_z_0_z_yyz[k] * ab_y + g_z_0_z_yyyz[k];

                g_z_0_yz_yzz[k] = -g_z_0_z_yzz[k] * ab_y + g_z_0_z_yyzz[k];

                g_z_0_yz_zzz[k] = -g_z_0_z_zzz[k] * ab_y + g_z_0_z_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_xxx = cbuffer.data(df_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_zz_xxy = cbuffer.data(df_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_zz_xxz = cbuffer.data(df_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_zz_xyy = cbuffer.data(df_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_zz_xyz = cbuffer.data(df_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_zz_xzz = cbuffer.data(df_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_zz_yyy = cbuffer.data(df_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_zz_yyz = cbuffer.data(df_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_zz_yzz = cbuffer.data(df_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_zz_zzz = cbuffer.data(df_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxx, g_z_0_z_xxxz, g_z_0_z_xxy, g_z_0_z_xxyz, g_z_0_z_xxz, g_z_0_z_xxzz, g_z_0_z_xyy, g_z_0_z_xyyz, g_z_0_z_xyz, g_z_0_z_xyzz, g_z_0_z_xzz, g_z_0_z_xzzz, g_z_0_z_yyy, g_z_0_z_yyyz, g_z_0_z_yyz, g_z_0_z_yyzz, g_z_0_z_yzz, g_z_0_z_yzzz, g_z_0_z_zzz, g_z_0_z_zzzz, g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xzz, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yzz, g_z_0_zz_zzz, g_z_xxx, g_z_xxy, g_z_xxz, g_z_xyy, g_z_xyz, g_z_xzz, g_z_yyy, g_z_yyz, g_z_yzz, g_z_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_xxx[k] = -g_z_xxx[k] - g_z_0_z_xxx[k] * ab_z + g_z_0_z_xxxz[k];

                g_z_0_zz_xxy[k] = -g_z_xxy[k] - g_z_0_z_xxy[k] * ab_z + g_z_0_z_xxyz[k];

                g_z_0_zz_xxz[k] = -g_z_xxz[k] - g_z_0_z_xxz[k] * ab_z + g_z_0_z_xxzz[k];

                g_z_0_zz_xyy[k] = -g_z_xyy[k] - g_z_0_z_xyy[k] * ab_z + g_z_0_z_xyyz[k];

                g_z_0_zz_xyz[k] = -g_z_xyz[k] - g_z_0_z_xyz[k] * ab_z + g_z_0_z_xyzz[k];

                g_z_0_zz_xzz[k] = -g_z_xzz[k] - g_z_0_z_xzz[k] * ab_z + g_z_0_z_xzzz[k];

                g_z_0_zz_yyy[k] = -g_z_yyy[k] - g_z_0_z_yyy[k] * ab_z + g_z_0_z_yyyz[k];

                g_z_0_zz_yyz[k] = -g_z_yyz[k] - g_z_0_z_yyz[k] * ab_z + g_z_0_z_yyzz[k];

                g_z_0_zz_yzz[k] = -g_z_yzz[k] - g_z_0_z_yzz[k] * ab_z + g_z_0_z_yzzz[k];

                g_z_0_zz_zzz[k] = -g_z_zzz[k] - g_z_0_z_zzz[k] * ab_z + g_z_0_z_zzzz[k];
            }
        }
    }
}

} // erirec namespace


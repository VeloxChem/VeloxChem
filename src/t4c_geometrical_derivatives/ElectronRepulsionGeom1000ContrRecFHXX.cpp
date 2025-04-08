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

#include "ElectronRepulsionGeom1000ContrRecFHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_fhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_fhxx,
                                            const size_t idx_dhxx,
                                            const size_t idx_geom_10_dhxx,
                                            const size_t idx_geom_10_dixx,
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
            /// Set up components of auxilary buffer : DHSS

            const auto dh_off = idx_dhxx + i * dcomps + j;

            auto g_xx_xxxxx = cbuffer.data(dh_off + 0 * ccomps * dcomps);

            auto g_xx_xxxxy = cbuffer.data(dh_off + 1 * ccomps * dcomps);

            auto g_xx_xxxxz = cbuffer.data(dh_off + 2 * ccomps * dcomps);

            auto g_xx_xxxyy = cbuffer.data(dh_off + 3 * ccomps * dcomps);

            auto g_xx_xxxyz = cbuffer.data(dh_off + 4 * ccomps * dcomps);

            auto g_xx_xxxzz = cbuffer.data(dh_off + 5 * ccomps * dcomps);

            auto g_xx_xxyyy = cbuffer.data(dh_off + 6 * ccomps * dcomps);

            auto g_xx_xxyyz = cbuffer.data(dh_off + 7 * ccomps * dcomps);

            auto g_xx_xxyzz = cbuffer.data(dh_off + 8 * ccomps * dcomps);

            auto g_xx_xxzzz = cbuffer.data(dh_off + 9 * ccomps * dcomps);

            auto g_xx_xyyyy = cbuffer.data(dh_off + 10 * ccomps * dcomps);

            auto g_xx_xyyyz = cbuffer.data(dh_off + 11 * ccomps * dcomps);

            auto g_xx_xyyzz = cbuffer.data(dh_off + 12 * ccomps * dcomps);

            auto g_xx_xyzzz = cbuffer.data(dh_off + 13 * ccomps * dcomps);

            auto g_xx_xzzzz = cbuffer.data(dh_off + 14 * ccomps * dcomps);

            auto g_xx_yyyyy = cbuffer.data(dh_off + 15 * ccomps * dcomps);

            auto g_xx_yyyyz = cbuffer.data(dh_off + 16 * ccomps * dcomps);

            auto g_xx_yyyzz = cbuffer.data(dh_off + 17 * ccomps * dcomps);

            auto g_xx_yyzzz = cbuffer.data(dh_off + 18 * ccomps * dcomps);

            auto g_xx_yzzzz = cbuffer.data(dh_off + 19 * ccomps * dcomps);

            auto g_xx_zzzzz = cbuffer.data(dh_off + 20 * ccomps * dcomps);

            auto g_yy_xxxxx = cbuffer.data(dh_off + 63 * ccomps * dcomps);

            auto g_yy_xxxxy = cbuffer.data(dh_off + 64 * ccomps * dcomps);

            auto g_yy_xxxxz = cbuffer.data(dh_off + 65 * ccomps * dcomps);

            auto g_yy_xxxyy = cbuffer.data(dh_off + 66 * ccomps * dcomps);

            auto g_yy_xxxyz = cbuffer.data(dh_off + 67 * ccomps * dcomps);

            auto g_yy_xxxzz = cbuffer.data(dh_off + 68 * ccomps * dcomps);

            auto g_yy_xxyyy = cbuffer.data(dh_off + 69 * ccomps * dcomps);

            auto g_yy_xxyyz = cbuffer.data(dh_off + 70 * ccomps * dcomps);

            auto g_yy_xxyzz = cbuffer.data(dh_off + 71 * ccomps * dcomps);

            auto g_yy_xxzzz = cbuffer.data(dh_off + 72 * ccomps * dcomps);

            auto g_yy_xyyyy = cbuffer.data(dh_off + 73 * ccomps * dcomps);

            auto g_yy_xyyyz = cbuffer.data(dh_off + 74 * ccomps * dcomps);

            auto g_yy_xyyzz = cbuffer.data(dh_off + 75 * ccomps * dcomps);

            auto g_yy_xyzzz = cbuffer.data(dh_off + 76 * ccomps * dcomps);

            auto g_yy_xzzzz = cbuffer.data(dh_off + 77 * ccomps * dcomps);

            auto g_yy_yyyyy = cbuffer.data(dh_off + 78 * ccomps * dcomps);

            auto g_yy_yyyyz = cbuffer.data(dh_off + 79 * ccomps * dcomps);

            auto g_yy_yyyzz = cbuffer.data(dh_off + 80 * ccomps * dcomps);

            auto g_yy_yyzzz = cbuffer.data(dh_off + 81 * ccomps * dcomps);

            auto g_yy_yzzzz = cbuffer.data(dh_off + 82 * ccomps * dcomps);

            auto g_yy_zzzzz = cbuffer.data(dh_off + 83 * ccomps * dcomps);

            auto g_zz_xxxxx = cbuffer.data(dh_off + 105 * ccomps * dcomps);

            auto g_zz_xxxxy = cbuffer.data(dh_off + 106 * ccomps * dcomps);

            auto g_zz_xxxxz = cbuffer.data(dh_off + 107 * ccomps * dcomps);

            auto g_zz_xxxyy = cbuffer.data(dh_off + 108 * ccomps * dcomps);

            auto g_zz_xxxyz = cbuffer.data(dh_off + 109 * ccomps * dcomps);

            auto g_zz_xxxzz = cbuffer.data(dh_off + 110 * ccomps * dcomps);

            auto g_zz_xxyyy = cbuffer.data(dh_off + 111 * ccomps * dcomps);

            auto g_zz_xxyyz = cbuffer.data(dh_off + 112 * ccomps * dcomps);

            auto g_zz_xxyzz = cbuffer.data(dh_off + 113 * ccomps * dcomps);

            auto g_zz_xxzzz = cbuffer.data(dh_off + 114 * ccomps * dcomps);

            auto g_zz_xyyyy = cbuffer.data(dh_off + 115 * ccomps * dcomps);

            auto g_zz_xyyyz = cbuffer.data(dh_off + 116 * ccomps * dcomps);

            auto g_zz_xyyzz = cbuffer.data(dh_off + 117 * ccomps * dcomps);

            auto g_zz_xyzzz = cbuffer.data(dh_off + 118 * ccomps * dcomps);

            auto g_zz_xzzzz = cbuffer.data(dh_off + 119 * ccomps * dcomps);

            auto g_zz_yyyyy = cbuffer.data(dh_off + 120 * ccomps * dcomps);

            auto g_zz_yyyyz = cbuffer.data(dh_off + 121 * ccomps * dcomps);

            auto g_zz_yyyzz = cbuffer.data(dh_off + 122 * ccomps * dcomps);

            auto g_zz_yyzzz = cbuffer.data(dh_off + 123 * ccomps * dcomps);

            auto g_zz_yzzzz = cbuffer.data(dh_off + 124 * ccomps * dcomps);

            auto g_zz_zzzzz = cbuffer.data(dh_off + 125 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DHSS

            const auto dh_geom_10_off = idx_geom_10_dhxx + i * dcomps + j;

            auto g_x_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 251 * ccomps * dcomps);

            auto g_z_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 269 * ccomps * dcomps);

            auto g_z_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 270 * ccomps * dcomps);

            auto g_z_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 271 * ccomps * dcomps);

            auto g_z_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 272 * ccomps * dcomps);

            auto g_z_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 273 * ccomps * dcomps);

            auto g_z_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 274 * ccomps * dcomps);

            auto g_z_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 275 * ccomps * dcomps);

            auto g_z_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 276 * ccomps * dcomps);

            auto g_z_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 277 * ccomps * dcomps);

            auto g_z_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 278 * ccomps * dcomps);

            auto g_z_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 279 * ccomps * dcomps);

            auto g_z_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 280 * ccomps * dcomps);

            auto g_z_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 281 * ccomps * dcomps);

            auto g_z_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 282 * ccomps * dcomps);

            auto g_z_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 283 * ccomps * dcomps);

            auto g_z_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 284 * ccomps * dcomps);

            auto g_z_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 285 * ccomps * dcomps);

            auto g_z_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 286 * ccomps * dcomps);

            auto g_z_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 287 * ccomps * dcomps);

            auto g_z_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 288 * ccomps * dcomps);

            auto g_z_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 289 * ccomps * dcomps);

            auto g_z_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 290 * ccomps * dcomps);

            auto g_z_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 291 * ccomps * dcomps);

            auto g_z_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 292 * ccomps * dcomps);

            auto g_z_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 293 * ccomps * dcomps);

            auto g_z_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 294 * ccomps * dcomps);

            auto g_z_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 295 * ccomps * dcomps);

            auto g_z_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 296 * ccomps * dcomps);

            auto g_z_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 297 * ccomps * dcomps);

            auto g_z_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 298 * ccomps * dcomps);

            auto g_z_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 299 * ccomps * dcomps);

            auto g_z_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 300 * ccomps * dcomps);

            auto g_z_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 301 * ccomps * dcomps);

            auto g_z_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 302 * ccomps * dcomps);

            auto g_z_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 303 * ccomps * dcomps);

            auto g_z_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 304 * ccomps * dcomps);

            auto g_z_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 305 * ccomps * dcomps);

            auto g_z_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 306 * ccomps * dcomps);

            auto g_z_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 307 * ccomps * dcomps);

            auto g_z_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 308 * ccomps * dcomps);

            auto g_z_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 309 * ccomps * dcomps);

            auto g_z_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 310 * ccomps * dcomps);

            auto g_z_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 311 * ccomps * dcomps);

            auto g_z_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 312 * ccomps * dcomps);

            auto g_z_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 313 * ccomps * dcomps);

            auto g_z_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 314 * ccomps * dcomps);

            auto g_z_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 315 * ccomps * dcomps);

            auto g_z_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 316 * ccomps * dcomps);

            auto g_z_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 317 * ccomps * dcomps);

            auto g_z_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 318 * ccomps * dcomps);

            auto g_z_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 319 * ccomps * dcomps);

            auto g_z_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 320 * ccomps * dcomps);

            auto g_z_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 321 * ccomps * dcomps);

            auto g_z_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 322 * ccomps * dcomps);

            auto g_z_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 323 * ccomps * dcomps);

            auto g_z_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 324 * ccomps * dcomps);

            auto g_z_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 325 * ccomps * dcomps);

            auto g_z_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 326 * ccomps * dcomps);

            auto g_z_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 327 * ccomps * dcomps);

            auto g_z_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 328 * ccomps * dcomps);

            auto g_z_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 329 * ccomps * dcomps);

            auto g_z_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 330 * ccomps * dcomps);

            auto g_z_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 331 * ccomps * dcomps);

            auto g_z_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 332 * ccomps * dcomps);

            auto g_z_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 333 * ccomps * dcomps);

            auto g_z_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 334 * ccomps * dcomps);

            auto g_z_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 335 * ccomps * dcomps);

            auto g_z_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 357 * ccomps * dcomps);

            auto g_z_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 358 * ccomps * dcomps);

            auto g_z_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 359 * ccomps * dcomps);

            auto g_z_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 360 * ccomps * dcomps);

            auto g_z_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 361 * ccomps * dcomps);

            auto g_z_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 362 * ccomps * dcomps);

            auto g_z_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 363 * ccomps * dcomps);

            auto g_z_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 377 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DISS

            const auto di_geom_10_off = idx_geom_10_dixx + i * dcomps + j;

            auto g_x_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 335 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 377 * ccomps * dcomps);

            auto g_z_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 378 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 379 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 380 * ccomps * dcomps);

            auto g_z_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 381 * ccomps * dcomps);

            auto g_z_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 382 * ccomps * dcomps);

            auto g_z_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 383 * ccomps * dcomps);

            auto g_z_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 384 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 392 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 393 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 394 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 395 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 396 * ccomps * dcomps);

            auto g_z_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 397 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 398 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 399 * ccomps * dcomps);

            auto g_z_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 400 * ccomps * dcomps);

            auto g_z_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 401 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 402 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 403 * ccomps * dcomps);

            auto g_z_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 404 * ccomps * dcomps);

            auto g_z_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 405 * ccomps * dcomps);

            auto g_z_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 406 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 407 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 408 * ccomps * dcomps);

            auto g_z_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 409 * ccomps * dcomps);

            auto g_z_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 410 * ccomps * dcomps);

            auto g_z_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 411 * ccomps * dcomps);

            auto g_z_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 412 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 503 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fhxx

            const auto fh_geom_10_off = idx_geom_10_fhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxx, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyz, g_x_0_xx_yyyzz, g_x_0_xx_yyzzz, g_x_0_xx_yzzzz, g_x_0_xx_zzzzz, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzzz, g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz, g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxzz, g_xx_xxyyy, g_xx_xxyyz, g_xx_xxyzz, g_xx_xxzzz, g_xx_xyyyy, g_xx_xyyyz, g_xx_xyyzz, g_xx_xyzzz, g_xx_xzzzz, g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz, g_xx_yyzzz, g_xx_yzzzz, g_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_xxxxx[k] = -g_xx_xxxxx[k] - g_x_0_xx_xxxxx[k] * ab_x + g_x_0_xx_xxxxxx[k];

                g_x_0_xxx_xxxxy[k] = -g_xx_xxxxy[k] - g_x_0_xx_xxxxy[k] * ab_x + g_x_0_xx_xxxxxy[k];

                g_x_0_xxx_xxxxz[k] = -g_xx_xxxxz[k] - g_x_0_xx_xxxxz[k] * ab_x + g_x_0_xx_xxxxxz[k];

                g_x_0_xxx_xxxyy[k] = -g_xx_xxxyy[k] - g_x_0_xx_xxxyy[k] * ab_x + g_x_0_xx_xxxxyy[k];

                g_x_0_xxx_xxxyz[k] = -g_xx_xxxyz[k] - g_x_0_xx_xxxyz[k] * ab_x + g_x_0_xx_xxxxyz[k];

                g_x_0_xxx_xxxzz[k] = -g_xx_xxxzz[k] - g_x_0_xx_xxxzz[k] * ab_x + g_x_0_xx_xxxxzz[k];

                g_x_0_xxx_xxyyy[k] = -g_xx_xxyyy[k] - g_x_0_xx_xxyyy[k] * ab_x + g_x_0_xx_xxxyyy[k];

                g_x_0_xxx_xxyyz[k] = -g_xx_xxyyz[k] - g_x_0_xx_xxyyz[k] * ab_x + g_x_0_xx_xxxyyz[k];

                g_x_0_xxx_xxyzz[k] = -g_xx_xxyzz[k] - g_x_0_xx_xxyzz[k] * ab_x + g_x_0_xx_xxxyzz[k];

                g_x_0_xxx_xxzzz[k] = -g_xx_xxzzz[k] - g_x_0_xx_xxzzz[k] * ab_x + g_x_0_xx_xxxzzz[k];

                g_x_0_xxx_xyyyy[k] = -g_xx_xyyyy[k] - g_x_0_xx_xyyyy[k] * ab_x + g_x_0_xx_xxyyyy[k];

                g_x_0_xxx_xyyyz[k] = -g_xx_xyyyz[k] - g_x_0_xx_xyyyz[k] * ab_x + g_x_0_xx_xxyyyz[k];

                g_x_0_xxx_xyyzz[k] = -g_xx_xyyzz[k] - g_x_0_xx_xyyzz[k] * ab_x + g_x_0_xx_xxyyzz[k];

                g_x_0_xxx_xyzzz[k] = -g_xx_xyzzz[k] - g_x_0_xx_xyzzz[k] * ab_x + g_x_0_xx_xxyzzz[k];

                g_x_0_xxx_xzzzz[k] = -g_xx_xzzzz[k] - g_x_0_xx_xzzzz[k] * ab_x + g_x_0_xx_xxzzzz[k];

                g_x_0_xxx_yyyyy[k] = -g_xx_yyyyy[k] - g_x_0_xx_yyyyy[k] * ab_x + g_x_0_xx_xyyyyy[k];

                g_x_0_xxx_yyyyz[k] = -g_xx_yyyyz[k] - g_x_0_xx_yyyyz[k] * ab_x + g_x_0_xx_xyyyyz[k];

                g_x_0_xxx_yyyzz[k] = -g_xx_yyyzz[k] - g_x_0_xx_yyyzz[k] * ab_x + g_x_0_xx_xyyyzz[k];

                g_x_0_xxx_yyzzz[k] = -g_xx_yyzzz[k] - g_x_0_xx_yyzzz[k] * ab_x + g_x_0_xx_xyyzzz[k];

                g_x_0_xxx_yzzzz[k] = -g_xx_yzzzz[k] - g_x_0_xx_yzzzz[k] * ab_x + g_x_0_xx_xyzzzz[k];

                g_x_0_xxx_zzzzz[k] = -g_xx_zzzzz[k] - g_x_0_xx_zzzzz[k] * ab_x + g_x_0_xx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxx, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzz, g_x_0_xxy_xxxxx, g_x_0_xxy_xxxxy, g_x_0_xxy_xxxxz, g_x_0_xxy_xxxyy, g_x_0_xxy_xxxyz, g_x_0_xxy_xxxzz, g_x_0_xxy_xxyyy, g_x_0_xxy_xxyyz, g_x_0_xxy_xxyzz, g_x_0_xxy_xxzzz, g_x_0_xxy_xyyyy, g_x_0_xxy_xyyyz, g_x_0_xxy_xyyzz, g_x_0_xxy_xyzzz, g_x_0_xxy_xzzzz, g_x_0_xxy_yyyyy, g_x_0_xxy_yyyyz, g_x_0_xxy_yyyzz, g_x_0_xxy_yyzzz, g_x_0_xxy_yzzzz, g_x_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_xxxxx[k] = -g_x_0_xx_xxxxx[k] * ab_y + g_x_0_xx_xxxxxy[k];

                g_x_0_xxy_xxxxy[k] = -g_x_0_xx_xxxxy[k] * ab_y + g_x_0_xx_xxxxyy[k];

                g_x_0_xxy_xxxxz[k] = -g_x_0_xx_xxxxz[k] * ab_y + g_x_0_xx_xxxxyz[k];

                g_x_0_xxy_xxxyy[k] = -g_x_0_xx_xxxyy[k] * ab_y + g_x_0_xx_xxxyyy[k];

                g_x_0_xxy_xxxyz[k] = -g_x_0_xx_xxxyz[k] * ab_y + g_x_0_xx_xxxyyz[k];

                g_x_0_xxy_xxxzz[k] = -g_x_0_xx_xxxzz[k] * ab_y + g_x_0_xx_xxxyzz[k];

                g_x_0_xxy_xxyyy[k] = -g_x_0_xx_xxyyy[k] * ab_y + g_x_0_xx_xxyyyy[k];

                g_x_0_xxy_xxyyz[k] = -g_x_0_xx_xxyyz[k] * ab_y + g_x_0_xx_xxyyyz[k];

                g_x_0_xxy_xxyzz[k] = -g_x_0_xx_xxyzz[k] * ab_y + g_x_0_xx_xxyyzz[k];

                g_x_0_xxy_xxzzz[k] = -g_x_0_xx_xxzzz[k] * ab_y + g_x_0_xx_xxyzzz[k];

                g_x_0_xxy_xyyyy[k] = -g_x_0_xx_xyyyy[k] * ab_y + g_x_0_xx_xyyyyy[k];

                g_x_0_xxy_xyyyz[k] = -g_x_0_xx_xyyyz[k] * ab_y + g_x_0_xx_xyyyyz[k];

                g_x_0_xxy_xyyzz[k] = -g_x_0_xx_xyyzz[k] * ab_y + g_x_0_xx_xyyyzz[k];

                g_x_0_xxy_xyzzz[k] = -g_x_0_xx_xyzzz[k] * ab_y + g_x_0_xx_xyyzzz[k];

                g_x_0_xxy_xzzzz[k] = -g_x_0_xx_xzzzz[k] * ab_y + g_x_0_xx_xyzzzz[k];

                g_x_0_xxy_yyyyy[k] = -g_x_0_xx_yyyyy[k] * ab_y + g_x_0_xx_yyyyyy[k];

                g_x_0_xxy_yyyyz[k] = -g_x_0_xx_yyyyz[k] * ab_y + g_x_0_xx_yyyyyz[k];

                g_x_0_xxy_yyyzz[k] = -g_x_0_xx_yyyzz[k] * ab_y + g_x_0_xx_yyyyzz[k];

                g_x_0_xxy_yyzzz[k] = -g_x_0_xx_yyzzz[k] * ab_y + g_x_0_xx_yyyzzz[k];

                g_x_0_xxy_yzzzz[k] = -g_x_0_xx_yzzzz[k] * ab_y + g_x_0_xx_yyzzzz[k];

                g_x_0_xxy_zzzzz[k] = -g_x_0_xx_zzzzz[k] * ab_y + g_x_0_xx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxx, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzz, g_x_0_xx_zzzzzz, g_x_0_xxz_xxxxx, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_xxxxx[k] = -g_x_0_xx_xxxxx[k] * ab_z + g_x_0_xx_xxxxxz[k];

                g_x_0_xxz_xxxxy[k] = -g_x_0_xx_xxxxy[k] * ab_z + g_x_0_xx_xxxxyz[k];

                g_x_0_xxz_xxxxz[k] = -g_x_0_xx_xxxxz[k] * ab_z + g_x_0_xx_xxxxzz[k];

                g_x_0_xxz_xxxyy[k] = -g_x_0_xx_xxxyy[k] * ab_z + g_x_0_xx_xxxyyz[k];

                g_x_0_xxz_xxxyz[k] = -g_x_0_xx_xxxyz[k] * ab_z + g_x_0_xx_xxxyzz[k];

                g_x_0_xxz_xxxzz[k] = -g_x_0_xx_xxxzz[k] * ab_z + g_x_0_xx_xxxzzz[k];

                g_x_0_xxz_xxyyy[k] = -g_x_0_xx_xxyyy[k] * ab_z + g_x_0_xx_xxyyyz[k];

                g_x_0_xxz_xxyyz[k] = -g_x_0_xx_xxyyz[k] * ab_z + g_x_0_xx_xxyyzz[k];

                g_x_0_xxz_xxyzz[k] = -g_x_0_xx_xxyzz[k] * ab_z + g_x_0_xx_xxyzzz[k];

                g_x_0_xxz_xxzzz[k] = -g_x_0_xx_xxzzz[k] * ab_z + g_x_0_xx_xxzzzz[k];

                g_x_0_xxz_xyyyy[k] = -g_x_0_xx_xyyyy[k] * ab_z + g_x_0_xx_xyyyyz[k];

                g_x_0_xxz_xyyyz[k] = -g_x_0_xx_xyyyz[k] * ab_z + g_x_0_xx_xyyyzz[k];

                g_x_0_xxz_xyyzz[k] = -g_x_0_xx_xyyzz[k] * ab_z + g_x_0_xx_xyyzzz[k];

                g_x_0_xxz_xyzzz[k] = -g_x_0_xx_xyzzz[k] * ab_z + g_x_0_xx_xyzzzz[k];

                g_x_0_xxz_xzzzz[k] = -g_x_0_xx_xzzzz[k] * ab_z + g_x_0_xx_xzzzzz[k];

                g_x_0_xxz_yyyyy[k] = -g_x_0_xx_yyyyy[k] * ab_z + g_x_0_xx_yyyyyz[k];

                g_x_0_xxz_yyyyz[k] = -g_x_0_xx_yyyyz[k] * ab_z + g_x_0_xx_yyyyzz[k];

                g_x_0_xxz_yyyzz[k] = -g_x_0_xx_yyyzz[k] * ab_z + g_x_0_xx_yyyzzz[k];

                g_x_0_xxz_yyzzz[k] = -g_x_0_xx_yyzzz[k] * ab_z + g_x_0_xx_yyzzzz[k];

                g_x_0_xxz_yzzzz[k] = -g_x_0_xx_yzzzz[k] * ab_z + g_x_0_xx_yzzzzz[k];

                g_x_0_xxz_zzzzz[k] = -g_x_0_xx_zzzzz[k] * ab_z + g_x_0_xx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xy_xxxxx, g_x_0_xy_xxxxxy, g_x_0_xy_xxxxy, g_x_0_xy_xxxxyy, g_x_0_xy_xxxxyz, g_x_0_xy_xxxxz, g_x_0_xy_xxxyy, g_x_0_xy_xxxyyy, g_x_0_xy_xxxyyz, g_x_0_xy_xxxyz, g_x_0_xy_xxxyzz, g_x_0_xy_xxxzz, g_x_0_xy_xxyyy, g_x_0_xy_xxyyyy, g_x_0_xy_xxyyyz, g_x_0_xy_xxyyz, g_x_0_xy_xxyyzz, g_x_0_xy_xxyzz, g_x_0_xy_xxyzzz, g_x_0_xy_xxzzz, g_x_0_xy_xyyyy, g_x_0_xy_xyyyyy, g_x_0_xy_xyyyyz, g_x_0_xy_xyyyz, g_x_0_xy_xyyyzz, g_x_0_xy_xyyzz, g_x_0_xy_xyyzzz, g_x_0_xy_xyzzz, g_x_0_xy_xyzzzz, g_x_0_xy_xzzzz, g_x_0_xy_yyyyy, g_x_0_xy_yyyyyy, g_x_0_xy_yyyyyz, g_x_0_xy_yyyyz, g_x_0_xy_yyyyzz, g_x_0_xy_yyyzz, g_x_0_xy_yyyzzz, g_x_0_xy_yyzzz, g_x_0_xy_yyzzzz, g_x_0_xy_yzzzz, g_x_0_xy_yzzzzz, g_x_0_xy_zzzzz, g_x_0_xyy_xxxxx, g_x_0_xyy_xxxxy, g_x_0_xyy_xxxxz, g_x_0_xyy_xxxyy, g_x_0_xyy_xxxyz, g_x_0_xyy_xxxzz, g_x_0_xyy_xxyyy, g_x_0_xyy_xxyyz, g_x_0_xyy_xxyzz, g_x_0_xyy_xxzzz, g_x_0_xyy_xyyyy, g_x_0_xyy_xyyyz, g_x_0_xyy_xyyzz, g_x_0_xyy_xyzzz, g_x_0_xyy_xzzzz, g_x_0_xyy_yyyyy, g_x_0_xyy_yyyyz, g_x_0_xyy_yyyzz, g_x_0_xyy_yyzzz, g_x_0_xyy_yzzzz, g_x_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_xxxxx[k] = -g_x_0_xy_xxxxx[k] * ab_y + g_x_0_xy_xxxxxy[k];

                g_x_0_xyy_xxxxy[k] = -g_x_0_xy_xxxxy[k] * ab_y + g_x_0_xy_xxxxyy[k];

                g_x_0_xyy_xxxxz[k] = -g_x_0_xy_xxxxz[k] * ab_y + g_x_0_xy_xxxxyz[k];

                g_x_0_xyy_xxxyy[k] = -g_x_0_xy_xxxyy[k] * ab_y + g_x_0_xy_xxxyyy[k];

                g_x_0_xyy_xxxyz[k] = -g_x_0_xy_xxxyz[k] * ab_y + g_x_0_xy_xxxyyz[k];

                g_x_0_xyy_xxxzz[k] = -g_x_0_xy_xxxzz[k] * ab_y + g_x_0_xy_xxxyzz[k];

                g_x_0_xyy_xxyyy[k] = -g_x_0_xy_xxyyy[k] * ab_y + g_x_0_xy_xxyyyy[k];

                g_x_0_xyy_xxyyz[k] = -g_x_0_xy_xxyyz[k] * ab_y + g_x_0_xy_xxyyyz[k];

                g_x_0_xyy_xxyzz[k] = -g_x_0_xy_xxyzz[k] * ab_y + g_x_0_xy_xxyyzz[k];

                g_x_0_xyy_xxzzz[k] = -g_x_0_xy_xxzzz[k] * ab_y + g_x_0_xy_xxyzzz[k];

                g_x_0_xyy_xyyyy[k] = -g_x_0_xy_xyyyy[k] * ab_y + g_x_0_xy_xyyyyy[k];

                g_x_0_xyy_xyyyz[k] = -g_x_0_xy_xyyyz[k] * ab_y + g_x_0_xy_xyyyyz[k];

                g_x_0_xyy_xyyzz[k] = -g_x_0_xy_xyyzz[k] * ab_y + g_x_0_xy_xyyyzz[k];

                g_x_0_xyy_xyzzz[k] = -g_x_0_xy_xyzzz[k] * ab_y + g_x_0_xy_xyyzzz[k];

                g_x_0_xyy_xzzzz[k] = -g_x_0_xy_xzzzz[k] * ab_y + g_x_0_xy_xyzzzz[k];

                g_x_0_xyy_yyyyy[k] = -g_x_0_xy_yyyyy[k] * ab_y + g_x_0_xy_yyyyyy[k];

                g_x_0_xyy_yyyyz[k] = -g_x_0_xy_yyyyz[k] * ab_y + g_x_0_xy_yyyyyz[k];

                g_x_0_xyy_yyyzz[k] = -g_x_0_xy_yyyzz[k] * ab_y + g_x_0_xy_yyyyzz[k];

                g_x_0_xyy_yyzzz[k] = -g_x_0_xy_yyzzz[k] * ab_y + g_x_0_xy_yyyzzz[k];

                g_x_0_xyy_yzzzz[k] = -g_x_0_xy_yzzzz[k] * ab_y + g_x_0_xy_yyzzzz[k];

                g_x_0_xyy_zzzzz[k] = -g_x_0_xy_zzzzz[k] * ab_y + g_x_0_xy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyz_xxxxx, g_x_0_xyz_xxxxy, g_x_0_xyz_xxxxz, g_x_0_xyz_xxxyy, g_x_0_xyz_xxxyz, g_x_0_xyz_xxxzz, g_x_0_xyz_xxyyy, g_x_0_xyz_xxyyz, g_x_0_xyz_xxyzz, g_x_0_xyz_xxzzz, g_x_0_xyz_xyyyy, g_x_0_xyz_xyyyz, g_x_0_xyz_xyyzz, g_x_0_xyz_xyzzz, g_x_0_xyz_xzzzz, g_x_0_xyz_yyyyy, g_x_0_xyz_yyyyz, g_x_0_xyz_yyyzz, g_x_0_xyz_yyzzz, g_x_0_xyz_yzzzz, g_x_0_xyz_zzzzz, g_x_0_xz_xxxxx, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxy, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxz, g_x_0_xz_xxxyy, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxzz, g_x_0_xz_xxyyy, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxzzz, g_x_0_xz_xyyyy, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xzzzz, g_x_0_xz_yyyyy, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_xxxxx[k] = -g_x_0_xz_xxxxx[k] * ab_y + g_x_0_xz_xxxxxy[k];

                g_x_0_xyz_xxxxy[k] = -g_x_0_xz_xxxxy[k] * ab_y + g_x_0_xz_xxxxyy[k];

                g_x_0_xyz_xxxxz[k] = -g_x_0_xz_xxxxz[k] * ab_y + g_x_0_xz_xxxxyz[k];

                g_x_0_xyz_xxxyy[k] = -g_x_0_xz_xxxyy[k] * ab_y + g_x_0_xz_xxxyyy[k];

                g_x_0_xyz_xxxyz[k] = -g_x_0_xz_xxxyz[k] * ab_y + g_x_0_xz_xxxyyz[k];

                g_x_0_xyz_xxxzz[k] = -g_x_0_xz_xxxzz[k] * ab_y + g_x_0_xz_xxxyzz[k];

                g_x_0_xyz_xxyyy[k] = -g_x_0_xz_xxyyy[k] * ab_y + g_x_0_xz_xxyyyy[k];

                g_x_0_xyz_xxyyz[k] = -g_x_0_xz_xxyyz[k] * ab_y + g_x_0_xz_xxyyyz[k];

                g_x_0_xyz_xxyzz[k] = -g_x_0_xz_xxyzz[k] * ab_y + g_x_0_xz_xxyyzz[k];

                g_x_0_xyz_xxzzz[k] = -g_x_0_xz_xxzzz[k] * ab_y + g_x_0_xz_xxyzzz[k];

                g_x_0_xyz_xyyyy[k] = -g_x_0_xz_xyyyy[k] * ab_y + g_x_0_xz_xyyyyy[k];

                g_x_0_xyz_xyyyz[k] = -g_x_0_xz_xyyyz[k] * ab_y + g_x_0_xz_xyyyyz[k];

                g_x_0_xyz_xyyzz[k] = -g_x_0_xz_xyyzz[k] * ab_y + g_x_0_xz_xyyyzz[k];

                g_x_0_xyz_xyzzz[k] = -g_x_0_xz_xyzzz[k] * ab_y + g_x_0_xz_xyyzzz[k];

                g_x_0_xyz_xzzzz[k] = -g_x_0_xz_xzzzz[k] * ab_y + g_x_0_xz_xyzzzz[k];

                g_x_0_xyz_yyyyy[k] = -g_x_0_xz_yyyyy[k] * ab_y + g_x_0_xz_yyyyyy[k];

                g_x_0_xyz_yyyyz[k] = -g_x_0_xz_yyyyz[k] * ab_y + g_x_0_xz_yyyyyz[k];

                g_x_0_xyz_yyyzz[k] = -g_x_0_xz_yyyzz[k] * ab_y + g_x_0_xz_yyyyzz[k];

                g_x_0_xyz_yyzzz[k] = -g_x_0_xz_yyzzz[k] * ab_y + g_x_0_xz_yyyzzz[k];

                g_x_0_xyz_yzzzz[k] = -g_x_0_xz_yzzzz[k] * ab_y + g_x_0_xz_yyzzzz[k];

                g_x_0_xyz_zzzzz[k] = -g_x_0_xz_zzzzz[k] * ab_y + g_x_0_xz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xz_xxxxx, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxy, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxyy, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxyyy, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xyyyy, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_yyyyy, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_zzzzz, g_x_0_xz_zzzzzz, g_x_0_xzz_xxxxx, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_xxxxx[k] = -g_x_0_xz_xxxxx[k] * ab_z + g_x_0_xz_xxxxxz[k];

                g_x_0_xzz_xxxxy[k] = -g_x_0_xz_xxxxy[k] * ab_z + g_x_0_xz_xxxxyz[k];

                g_x_0_xzz_xxxxz[k] = -g_x_0_xz_xxxxz[k] * ab_z + g_x_0_xz_xxxxzz[k];

                g_x_0_xzz_xxxyy[k] = -g_x_0_xz_xxxyy[k] * ab_z + g_x_0_xz_xxxyyz[k];

                g_x_0_xzz_xxxyz[k] = -g_x_0_xz_xxxyz[k] * ab_z + g_x_0_xz_xxxyzz[k];

                g_x_0_xzz_xxxzz[k] = -g_x_0_xz_xxxzz[k] * ab_z + g_x_0_xz_xxxzzz[k];

                g_x_0_xzz_xxyyy[k] = -g_x_0_xz_xxyyy[k] * ab_z + g_x_0_xz_xxyyyz[k];

                g_x_0_xzz_xxyyz[k] = -g_x_0_xz_xxyyz[k] * ab_z + g_x_0_xz_xxyyzz[k];

                g_x_0_xzz_xxyzz[k] = -g_x_0_xz_xxyzz[k] * ab_z + g_x_0_xz_xxyzzz[k];

                g_x_0_xzz_xxzzz[k] = -g_x_0_xz_xxzzz[k] * ab_z + g_x_0_xz_xxzzzz[k];

                g_x_0_xzz_xyyyy[k] = -g_x_0_xz_xyyyy[k] * ab_z + g_x_0_xz_xyyyyz[k];

                g_x_0_xzz_xyyyz[k] = -g_x_0_xz_xyyyz[k] * ab_z + g_x_0_xz_xyyyzz[k];

                g_x_0_xzz_xyyzz[k] = -g_x_0_xz_xyyzz[k] * ab_z + g_x_0_xz_xyyzzz[k];

                g_x_0_xzz_xyzzz[k] = -g_x_0_xz_xyzzz[k] * ab_z + g_x_0_xz_xyzzzz[k];

                g_x_0_xzz_xzzzz[k] = -g_x_0_xz_xzzzz[k] * ab_z + g_x_0_xz_xzzzzz[k];

                g_x_0_xzz_yyyyy[k] = -g_x_0_xz_yyyyy[k] * ab_z + g_x_0_xz_yyyyyz[k];

                g_x_0_xzz_yyyyz[k] = -g_x_0_xz_yyyyz[k] * ab_z + g_x_0_xz_yyyyzz[k];

                g_x_0_xzz_yyyzz[k] = -g_x_0_xz_yyyzz[k] * ab_z + g_x_0_xz_yyyzzz[k];

                g_x_0_xzz_yyzzz[k] = -g_x_0_xz_yyzzz[k] * ab_z + g_x_0_xz_yyzzzz[k];

                g_x_0_xzz_yzzzz[k] = -g_x_0_xz_yzzzz[k] * ab_z + g_x_0_xz_yzzzzz[k];

                g_x_0_xzz_zzzzz[k] = -g_x_0_xz_zzzzz[k] * ab_z + g_x_0_xz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xxxxx, g_x_0_yy_xxxxxy, g_x_0_yy_xxxxy, g_x_0_yy_xxxxyy, g_x_0_yy_xxxxyz, g_x_0_yy_xxxxz, g_x_0_yy_xxxyy, g_x_0_yy_xxxyyy, g_x_0_yy_xxxyyz, g_x_0_yy_xxxyz, g_x_0_yy_xxxyzz, g_x_0_yy_xxxzz, g_x_0_yy_xxyyy, g_x_0_yy_xxyyyy, g_x_0_yy_xxyyyz, g_x_0_yy_xxyyz, g_x_0_yy_xxyyzz, g_x_0_yy_xxyzz, g_x_0_yy_xxyzzz, g_x_0_yy_xxzzz, g_x_0_yy_xyyyy, g_x_0_yy_xyyyyy, g_x_0_yy_xyyyyz, g_x_0_yy_xyyyz, g_x_0_yy_xyyyzz, g_x_0_yy_xyyzz, g_x_0_yy_xyyzzz, g_x_0_yy_xyzzz, g_x_0_yy_xyzzzz, g_x_0_yy_xzzzz, g_x_0_yy_yyyyy, g_x_0_yy_yyyyyy, g_x_0_yy_yyyyyz, g_x_0_yy_yyyyz, g_x_0_yy_yyyyzz, g_x_0_yy_yyyzz, g_x_0_yy_yyyzzz, g_x_0_yy_yyzzz, g_x_0_yy_yyzzzz, g_x_0_yy_yzzzz, g_x_0_yy_yzzzzz, g_x_0_yy_zzzzz, g_x_0_yyy_xxxxx, g_x_0_yyy_xxxxy, g_x_0_yyy_xxxxz, g_x_0_yyy_xxxyy, g_x_0_yyy_xxxyz, g_x_0_yyy_xxxzz, g_x_0_yyy_xxyyy, g_x_0_yyy_xxyyz, g_x_0_yyy_xxyzz, g_x_0_yyy_xxzzz, g_x_0_yyy_xyyyy, g_x_0_yyy_xyyyz, g_x_0_yyy_xyyzz, g_x_0_yyy_xyzzz, g_x_0_yyy_xzzzz, g_x_0_yyy_yyyyy, g_x_0_yyy_yyyyz, g_x_0_yyy_yyyzz, g_x_0_yyy_yyzzz, g_x_0_yyy_yzzzz, g_x_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_xxxxx[k] = -g_x_0_yy_xxxxx[k] * ab_y + g_x_0_yy_xxxxxy[k];

                g_x_0_yyy_xxxxy[k] = -g_x_0_yy_xxxxy[k] * ab_y + g_x_0_yy_xxxxyy[k];

                g_x_0_yyy_xxxxz[k] = -g_x_0_yy_xxxxz[k] * ab_y + g_x_0_yy_xxxxyz[k];

                g_x_0_yyy_xxxyy[k] = -g_x_0_yy_xxxyy[k] * ab_y + g_x_0_yy_xxxyyy[k];

                g_x_0_yyy_xxxyz[k] = -g_x_0_yy_xxxyz[k] * ab_y + g_x_0_yy_xxxyyz[k];

                g_x_0_yyy_xxxzz[k] = -g_x_0_yy_xxxzz[k] * ab_y + g_x_0_yy_xxxyzz[k];

                g_x_0_yyy_xxyyy[k] = -g_x_0_yy_xxyyy[k] * ab_y + g_x_0_yy_xxyyyy[k];

                g_x_0_yyy_xxyyz[k] = -g_x_0_yy_xxyyz[k] * ab_y + g_x_0_yy_xxyyyz[k];

                g_x_0_yyy_xxyzz[k] = -g_x_0_yy_xxyzz[k] * ab_y + g_x_0_yy_xxyyzz[k];

                g_x_0_yyy_xxzzz[k] = -g_x_0_yy_xxzzz[k] * ab_y + g_x_0_yy_xxyzzz[k];

                g_x_0_yyy_xyyyy[k] = -g_x_0_yy_xyyyy[k] * ab_y + g_x_0_yy_xyyyyy[k];

                g_x_0_yyy_xyyyz[k] = -g_x_0_yy_xyyyz[k] * ab_y + g_x_0_yy_xyyyyz[k];

                g_x_0_yyy_xyyzz[k] = -g_x_0_yy_xyyzz[k] * ab_y + g_x_0_yy_xyyyzz[k];

                g_x_0_yyy_xyzzz[k] = -g_x_0_yy_xyzzz[k] * ab_y + g_x_0_yy_xyyzzz[k];

                g_x_0_yyy_xzzzz[k] = -g_x_0_yy_xzzzz[k] * ab_y + g_x_0_yy_xyzzzz[k];

                g_x_0_yyy_yyyyy[k] = -g_x_0_yy_yyyyy[k] * ab_y + g_x_0_yy_yyyyyy[k];

                g_x_0_yyy_yyyyz[k] = -g_x_0_yy_yyyyz[k] * ab_y + g_x_0_yy_yyyyyz[k];

                g_x_0_yyy_yyyzz[k] = -g_x_0_yy_yyyzz[k] * ab_y + g_x_0_yy_yyyyzz[k];

                g_x_0_yyy_yyzzz[k] = -g_x_0_yy_yyzzz[k] * ab_y + g_x_0_yy_yyyzzz[k];

                g_x_0_yyy_yzzzz[k] = -g_x_0_yy_yzzzz[k] * ab_y + g_x_0_yy_yyzzzz[k];

                g_x_0_yyy_zzzzz[k] = -g_x_0_yy_zzzzz[k] * ab_y + g_x_0_yy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyz_xxxxx, g_x_0_yyz_xxxxy, g_x_0_yyz_xxxxz, g_x_0_yyz_xxxyy, g_x_0_yyz_xxxyz, g_x_0_yyz_xxxzz, g_x_0_yyz_xxyyy, g_x_0_yyz_xxyyz, g_x_0_yyz_xxyzz, g_x_0_yyz_xxzzz, g_x_0_yyz_xyyyy, g_x_0_yyz_xyyyz, g_x_0_yyz_xyyzz, g_x_0_yyz_xyzzz, g_x_0_yyz_xzzzz, g_x_0_yyz_yyyyy, g_x_0_yyz_yyyyz, g_x_0_yyz_yyyzz, g_x_0_yyz_yyzzz, g_x_0_yyz_yzzzz, g_x_0_yyz_zzzzz, g_x_0_yz_xxxxx, g_x_0_yz_xxxxxy, g_x_0_yz_xxxxy, g_x_0_yz_xxxxyy, g_x_0_yz_xxxxyz, g_x_0_yz_xxxxz, g_x_0_yz_xxxyy, g_x_0_yz_xxxyyy, g_x_0_yz_xxxyyz, g_x_0_yz_xxxyz, g_x_0_yz_xxxyzz, g_x_0_yz_xxxzz, g_x_0_yz_xxyyy, g_x_0_yz_xxyyyy, g_x_0_yz_xxyyyz, g_x_0_yz_xxyyz, g_x_0_yz_xxyyzz, g_x_0_yz_xxyzz, g_x_0_yz_xxyzzz, g_x_0_yz_xxzzz, g_x_0_yz_xyyyy, g_x_0_yz_xyyyyy, g_x_0_yz_xyyyyz, g_x_0_yz_xyyyz, g_x_0_yz_xyyyzz, g_x_0_yz_xyyzz, g_x_0_yz_xyyzzz, g_x_0_yz_xyzzz, g_x_0_yz_xyzzzz, g_x_0_yz_xzzzz, g_x_0_yz_yyyyy, g_x_0_yz_yyyyyy, g_x_0_yz_yyyyyz, g_x_0_yz_yyyyz, g_x_0_yz_yyyyzz, g_x_0_yz_yyyzz, g_x_0_yz_yyyzzz, g_x_0_yz_yyzzz, g_x_0_yz_yyzzzz, g_x_0_yz_yzzzz, g_x_0_yz_yzzzzz, g_x_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_xxxxx[k] = -g_x_0_yz_xxxxx[k] * ab_y + g_x_0_yz_xxxxxy[k];

                g_x_0_yyz_xxxxy[k] = -g_x_0_yz_xxxxy[k] * ab_y + g_x_0_yz_xxxxyy[k];

                g_x_0_yyz_xxxxz[k] = -g_x_0_yz_xxxxz[k] * ab_y + g_x_0_yz_xxxxyz[k];

                g_x_0_yyz_xxxyy[k] = -g_x_0_yz_xxxyy[k] * ab_y + g_x_0_yz_xxxyyy[k];

                g_x_0_yyz_xxxyz[k] = -g_x_0_yz_xxxyz[k] * ab_y + g_x_0_yz_xxxyyz[k];

                g_x_0_yyz_xxxzz[k] = -g_x_0_yz_xxxzz[k] * ab_y + g_x_0_yz_xxxyzz[k];

                g_x_0_yyz_xxyyy[k] = -g_x_0_yz_xxyyy[k] * ab_y + g_x_0_yz_xxyyyy[k];

                g_x_0_yyz_xxyyz[k] = -g_x_0_yz_xxyyz[k] * ab_y + g_x_0_yz_xxyyyz[k];

                g_x_0_yyz_xxyzz[k] = -g_x_0_yz_xxyzz[k] * ab_y + g_x_0_yz_xxyyzz[k];

                g_x_0_yyz_xxzzz[k] = -g_x_0_yz_xxzzz[k] * ab_y + g_x_0_yz_xxyzzz[k];

                g_x_0_yyz_xyyyy[k] = -g_x_0_yz_xyyyy[k] * ab_y + g_x_0_yz_xyyyyy[k];

                g_x_0_yyz_xyyyz[k] = -g_x_0_yz_xyyyz[k] * ab_y + g_x_0_yz_xyyyyz[k];

                g_x_0_yyz_xyyzz[k] = -g_x_0_yz_xyyzz[k] * ab_y + g_x_0_yz_xyyyzz[k];

                g_x_0_yyz_xyzzz[k] = -g_x_0_yz_xyzzz[k] * ab_y + g_x_0_yz_xyyzzz[k];

                g_x_0_yyz_xzzzz[k] = -g_x_0_yz_xzzzz[k] * ab_y + g_x_0_yz_xyzzzz[k];

                g_x_0_yyz_yyyyy[k] = -g_x_0_yz_yyyyy[k] * ab_y + g_x_0_yz_yyyyyy[k];

                g_x_0_yyz_yyyyz[k] = -g_x_0_yz_yyyyz[k] * ab_y + g_x_0_yz_yyyyyz[k];

                g_x_0_yyz_yyyzz[k] = -g_x_0_yz_yyyzz[k] * ab_y + g_x_0_yz_yyyyzz[k];

                g_x_0_yyz_yyzzz[k] = -g_x_0_yz_yyzzz[k] * ab_y + g_x_0_yz_yyyzzz[k];

                g_x_0_yyz_yzzzz[k] = -g_x_0_yz_yzzzz[k] * ab_y + g_x_0_yz_yyzzzz[k];

                g_x_0_yyz_zzzzz[k] = -g_x_0_yz_zzzzz[k] * ab_y + g_x_0_yz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzz_xxxxx, g_x_0_yzz_xxxxy, g_x_0_yzz_xxxxz, g_x_0_yzz_xxxyy, g_x_0_yzz_xxxyz, g_x_0_yzz_xxxzz, g_x_0_yzz_xxyyy, g_x_0_yzz_xxyyz, g_x_0_yzz_xxyzz, g_x_0_yzz_xxzzz, g_x_0_yzz_xyyyy, g_x_0_yzz_xyyyz, g_x_0_yzz_xyyzz, g_x_0_yzz_xyzzz, g_x_0_yzz_xzzzz, g_x_0_yzz_yyyyy, g_x_0_yzz_yyyyz, g_x_0_yzz_yyyzz, g_x_0_yzz_yyzzz, g_x_0_yzz_yzzzz, g_x_0_yzz_zzzzz, g_x_0_zz_xxxxx, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxy, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxz, g_x_0_zz_xxxyy, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxzz, g_x_0_zz_xxyyy, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxzzz, g_x_0_zz_xyyyy, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xzzzz, g_x_0_zz_yyyyy, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_xxxxx[k] = -g_x_0_zz_xxxxx[k] * ab_y + g_x_0_zz_xxxxxy[k];

                g_x_0_yzz_xxxxy[k] = -g_x_0_zz_xxxxy[k] * ab_y + g_x_0_zz_xxxxyy[k];

                g_x_0_yzz_xxxxz[k] = -g_x_0_zz_xxxxz[k] * ab_y + g_x_0_zz_xxxxyz[k];

                g_x_0_yzz_xxxyy[k] = -g_x_0_zz_xxxyy[k] * ab_y + g_x_0_zz_xxxyyy[k];

                g_x_0_yzz_xxxyz[k] = -g_x_0_zz_xxxyz[k] * ab_y + g_x_0_zz_xxxyyz[k];

                g_x_0_yzz_xxxzz[k] = -g_x_0_zz_xxxzz[k] * ab_y + g_x_0_zz_xxxyzz[k];

                g_x_0_yzz_xxyyy[k] = -g_x_0_zz_xxyyy[k] * ab_y + g_x_0_zz_xxyyyy[k];

                g_x_0_yzz_xxyyz[k] = -g_x_0_zz_xxyyz[k] * ab_y + g_x_0_zz_xxyyyz[k];

                g_x_0_yzz_xxyzz[k] = -g_x_0_zz_xxyzz[k] * ab_y + g_x_0_zz_xxyyzz[k];

                g_x_0_yzz_xxzzz[k] = -g_x_0_zz_xxzzz[k] * ab_y + g_x_0_zz_xxyzzz[k];

                g_x_0_yzz_xyyyy[k] = -g_x_0_zz_xyyyy[k] * ab_y + g_x_0_zz_xyyyyy[k];

                g_x_0_yzz_xyyyz[k] = -g_x_0_zz_xyyyz[k] * ab_y + g_x_0_zz_xyyyyz[k];

                g_x_0_yzz_xyyzz[k] = -g_x_0_zz_xyyzz[k] * ab_y + g_x_0_zz_xyyyzz[k];

                g_x_0_yzz_xyzzz[k] = -g_x_0_zz_xyzzz[k] * ab_y + g_x_0_zz_xyyzzz[k];

                g_x_0_yzz_xzzzz[k] = -g_x_0_zz_xzzzz[k] * ab_y + g_x_0_zz_xyzzzz[k];

                g_x_0_yzz_yyyyy[k] = -g_x_0_zz_yyyyy[k] * ab_y + g_x_0_zz_yyyyyy[k];

                g_x_0_yzz_yyyyz[k] = -g_x_0_zz_yyyyz[k] * ab_y + g_x_0_zz_yyyyyz[k];

                g_x_0_yzz_yyyzz[k] = -g_x_0_zz_yyyzz[k] * ab_y + g_x_0_zz_yyyyzz[k];

                g_x_0_yzz_yyzzz[k] = -g_x_0_zz_yyzzz[k] * ab_y + g_x_0_zz_yyyzzz[k];

                g_x_0_yzz_yzzzz[k] = -g_x_0_zz_yzzzz[k] * ab_y + g_x_0_zz_yyzzzz[k];

                g_x_0_yzz_zzzzz[k] = -g_x_0_zz_zzzzz[k] * ab_y + g_x_0_zz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xxxxx, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxy, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxyy, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxyyy, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xyyyy, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_yyyyy, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_zzzzz, g_x_0_zz_zzzzzz, g_x_0_zzz_xxxxx, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_xxxxx[k] = -g_x_0_zz_xxxxx[k] * ab_z + g_x_0_zz_xxxxxz[k];

                g_x_0_zzz_xxxxy[k] = -g_x_0_zz_xxxxy[k] * ab_z + g_x_0_zz_xxxxyz[k];

                g_x_0_zzz_xxxxz[k] = -g_x_0_zz_xxxxz[k] * ab_z + g_x_0_zz_xxxxzz[k];

                g_x_0_zzz_xxxyy[k] = -g_x_0_zz_xxxyy[k] * ab_z + g_x_0_zz_xxxyyz[k];

                g_x_0_zzz_xxxyz[k] = -g_x_0_zz_xxxyz[k] * ab_z + g_x_0_zz_xxxyzz[k];

                g_x_0_zzz_xxxzz[k] = -g_x_0_zz_xxxzz[k] * ab_z + g_x_0_zz_xxxzzz[k];

                g_x_0_zzz_xxyyy[k] = -g_x_0_zz_xxyyy[k] * ab_z + g_x_0_zz_xxyyyz[k];

                g_x_0_zzz_xxyyz[k] = -g_x_0_zz_xxyyz[k] * ab_z + g_x_0_zz_xxyyzz[k];

                g_x_0_zzz_xxyzz[k] = -g_x_0_zz_xxyzz[k] * ab_z + g_x_0_zz_xxyzzz[k];

                g_x_0_zzz_xxzzz[k] = -g_x_0_zz_xxzzz[k] * ab_z + g_x_0_zz_xxzzzz[k];

                g_x_0_zzz_xyyyy[k] = -g_x_0_zz_xyyyy[k] * ab_z + g_x_0_zz_xyyyyz[k];

                g_x_0_zzz_xyyyz[k] = -g_x_0_zz_xyyyz[k] * ab_z + g_x_0_zz_xyyyzz[k];

                g_x_0_zzz_xyyzz[k] = -g_x_0_zz_xyyzz[k] * ab_z + g_x_0_zz_xyyzzz[k];

                g_x_0_zzz_xyzzz[k] = -g_x_0_zz_xyzzz[k] * ab_z + g_x_0_zz_xyzzzz[k];

                g_x_0_zzz_xzzzz[k] = -g_x_0_zz_xzzzz[k] * ab_z + g_x_0_zz_xzzzzz[k];

                g_x_0_zzz_yyyyy[k] = -g_x_0_zz_yyyyy[k] * ab_z + g_x_0_zz_yyyyyz[k];

                g_x_0_zzz_yyyyz[k] = -g_x_0_zz_yyyyz[k] * ab_z + g_x_0_zz_yyyyzz[k];

                g_x_0_zzz_yyyzz[k] = -g_x_0_zz_yyyzz[k] * ab_z + g_x_0_zz_yyyzzz[k];

                g_x_0_zzz_yyzzz[k] = -g_x_0_zz_yyzzz[k] * ab_z + g_x_0_zz_yyzzzz[k];

                g_x_0_zzz_yzzzz[k] = -g_x_0_zz_yzzzz[k] * ab_z + g_x_0_zz_yzzzzz[k];

                g_x_0_zzz_zzzzz[k] = -g_x_0_zz_zzzzz[k] * ab_z + g_x_0_zz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_xxxxx, g_y_0_xx_xxxxxx, g_y_0_xx_xxxxxy, g_y_0_xx_xxxxxz, g_y_0_xx_xxxxy, g_y_0_xx_xxxxyy, g_y_0_xx_xxxxyz, g_y_0_xx_xxxxz, g_y_0_xx_xxxxzz, g_y_0_xx_xxxyy, g_y_0_xx_xxxyyy, g_y_0_xx_xxxyyz, g_y_0_xx_xxxyz, g_y_0_xx_xxxyzz, g_y_0_xx_xxxzz, g_y_0_xx_xxxzzz, g_y_0_xx_xxyyy, g_y_0_xx_xxyyyy, g_y_0_xx_xxyyyz, g_y_0_xx_xxyyz, g_y_0_xx_xxyyzz, g_y_0_xx_xxyzz, g_y_0_xx_xxyzzz, g_y_0_xx_xxzzz, g_y_0_xx_xxzzzz, g_y_0_xx_xyyyy, g_y_0_xx_xyyyyy, g_y_0_xx_xyyyyz, g_y_0_xx_xyyyz, g_y_0_xx_xyyyzz, g_y_0_xx_xyyzz, g_y_0_xx_xyyzzz, g_y_0_xx_xyzzz, g_y_0_xx_xyzzzz, g_y_0_xx_xzzzz, g_y_0_xx_xzzzzz, g_y_0_xx_yyyyy, g_y_0_xx_yyyyz, g_y_0_xx_yyyzz, g_y_0_xx_yyzzz, g_y_0_xx_yzzzz, g_y_0_xx_zzzzz, g_y_0_xxx_xxxxx, g_y_0_xxx_xxxxy, g_y_0_xxx_xxxxz, g_y_0_xxx_xxxyy, g_y_0_xxx_xxxyz, g_y_0_xxx_xxxzz, g_y_0_xxx_xxyyy, g_y_0_xxx_xxyyz, g_y_0_xxx_xxyzz, g_y_0_xxx_xxzzz, g_y_0_xxx_xyyyy, g_y_0_xxx_xyyyz, g_y_0_xxx_xyyzz, g_y_0_xxx_xyzzz, g_y_0_xxx_xzzzz, g_y_0_xxx_yyyyy, g_y_0_xxx_yyyyz, g_y_0_xxx_yyyzz, g_y_0_xxx_yyzzz, g_y_0_xxx_yzzzz, g_y_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_xxxxx[k] = -g_y_0_xx_xxxxx[k] * ab_x + g_y_0_xx_xxxxxx[k];

                g_y_0_xxx_xxxxy[k] = -g_y_0_xx_xxxxy[k] * ab_x + g_y_0_xx_xxxxxy[k];

                g_y_0_xxx_xxxxz[k] = -g_y_0_xx_xxxxz[k] * ab_x + g_y_0_xx_xxxxxz[k];

                g_y_0_xxx_xxxyy[k] = -g_y_0_xx_xxxyy[k] * ab_x + g_y_0_xx_xxxxyy[k];

                g_y_0_xxx_xxxyz[k] = -g_y_0_xx_xxxyz[k] * ab_x + g_y_0_xx_xxxxyz[k];

                g_y_0_xxx_xxxzz[k] = -g_y_0_xx_xxxzz[k] * ab_x + g_y_0_xx_xxxxzz[k];

                g_y_0_xxx_xxyyy[k] = -g_y_0_xx_xxyyy[k] * ab_x + g_y_0_xx_xxxyyy[k];

                g_y_0_xxx_xxyyz[k] = -g_y_0_xx_xxyyz[k] * ab_x + g_y_0_xx_xxxyyz[k];

                g_y_0_xxx_xxyzz[k] = -g_y_0_xx_xxyzz[k] * ab_x + g_y_0_xx_xxxyzz[k];

                g_y_0_xxx_xxzzz[k] = -g_y_0_xx_xxzzz[k] * ab_x + g_y_0_xx_xxxzzz[k];

                g_y_0_xxx_xyyyy[k] = -g_y_0_xx_xyyyy[k] * ab_x + g_y_0_xx_xxyyyy[k];

                g_y_0_xxx_xyyyz[k] = -g_y_0_xx_xyyyz[k] * ab_x + g_y_0_xx_xxyyyz[k];

                g_y_0_xxx_xyyzz[k] = -g_y_0_xx_xyyzz[k] * ab_x + g_y_0_xx_xxyyzz[k];

                g_y_0_xxx_xyzzz[k] = -g_y_0_xx_xyzzz[k] * ab_x + g_y_0_xx_xxyzzz[k];

                g_y_0_xxx_xzzzz[k] = -g_y_0_xx_xzzzz[k] * ab_x + g_y_0_xx_xxzzzz[k];

                g_y_0_xxx_yyyyy[k] = -g_y_0_xx_yyyyy[k] * ab_x + g_y_0_xx_xyyyyy[k];

                g_y_0_xxx_yyyyz[k] = -g_y_0_xx_yyyyz[k] * ab_x + g_y_0_xx_xyyyyz[k];

                g_y_0_xxx_yyyzz[k] = -g_y_0_xx_yyyzz[k] * ab_x + g_y_0_xx_xyyyzz[k];

                g_y_0_xxx_yyzzz[k] = -g_y_0_xx_yyzzz[k] * ab_x + g_y_0_xx_xyyzzz[k];

                g_y_0_xxx_yzzzz[k] = -g_y_0_xx_yzzzz[k] * ab_x + g_y_0_xx_xyzzzz[k];

                g_y_0_xxx_zzzzz[k] = -g_y_0_xx_zzzzz[k] * ab_x + g_y_0_xx_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_xxxxx, g_y_0_xxy_xxxxy, g_y_0_xxy_xxxxz, g_y_0_xxy_xxxyy, g_y_0_xxy_xxxyz, g_y_0_xxy_xxxzz, g_y_0_xxy_xxyyy, g_y_0_xxy_xxyyz, g_y_0_xxy_xxyzz, g_y_0_xxy_xxzzz, g_y_0_xxy_xyyyy, g_y_0_xxy_xyyyz, g_y_0_xxy_xyyzz, g_y_0_xxy_xyzzz, g_y_0_xxy_xzzzz, g_y_0_xxy_yyyyy, g_y_0_xxy_yyyyz, g_y_0_xxy_yyyzz, g_y_0_xxy_yyzzz, g_y_0_xxy_yzzzz, g_y_0_xxy_zzzzz, g_y_0_xy_xxxxx, g_y_0_xy_xxxxxx, g_y_0_xy_xxxxxy, g_y_0_xy_xxxxxz, g_y_0_xy_xxxxy, g_y_0_xy_xxxxyy, g_y_0_xy_xxxxyz, g_y_0_xy_xxxxz, g_y_0_xy_xxxxzz, g_y_0_xy_xxxyy, g_y_0_xy_xxxyyy, g_y_0_xy_xxxyyz, g_y_0_xy_xxxyz, g_y_0_xy_xxxyzz, g_y_0_xy_xxxzz, g_y_0_xy_xxxzzz, g_y_0_xy_xxyyy, g_y_0_xy_xxyyyy, g_y_0_xy_xxyyyz, g_y_0_xy_xxyyz, g_y_0_xy_xxyyzz, g_y_0_xy_xxyzz, g_y_0_xy_xxyzzz, g_y_0_xy_xxzzz, g_y_0_xy_xxzzzz, g_y_0_xy_xyyyy, g_y_0_xy_xyyyyy, g_y_0_xy_xyyyyz, g_y_0_xy_xyyyz, g_y_0_xy_xyyyzz, g_y_0_xy_xyyzz, g_y_0_xy_xyyzzz, g_y_0_xy_xyzzz, g_y_0_xy_xyzzzz, g_y_0_xy_xzzzz, g_y_0_xy_xzzzzz, g_y_0_xy_yyyyy, g_y_0_xy_yyyyz, g_y_0_xy_yyyzz, g_y_0_xy_yyzzz, g_y_0_xy_yzzzz, g_y_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_xxxxx[k] = -g_y_0_xy_xxxxx[k] * ab_x + g_y_0_xy_xxxxxx[k];

                g_y_0_xxy_xxxxy[k] = -g_y_0_xy_xxxxy[k] * ab_x + g_y_0_xy_xxxxxy[k];

                g_y_0_xxy_xxxxz[k] = -g_y_0_xy_xxxxz[k] * ab_x + g_y_0_xy_xxxxxz[k];

                g_y_0_xxy_xxxyy[k] = -g_y_0_xy_xxxyy[k] * ab_x + g_y_0_xy_xxxxyy[k];

                g_y_0_xxy_xxxyz[k] = -g_y_0_xy_xxxyz[k] * ab_x + g_y_0_xy_xxxxyz[k];

                g_y_0_xxy_xxxzz[k] = -g_y_0_xy_xxxzz[k] * ab_x + g_y_0_xy_xxxxzz[k];

                g_y_0_xxy_xxyyy[k] = -g_y_0_xy_xxyyy[k] * ab_x + g_y_0_xy_xxxyyy[k];

                g_y_0_xxy_xxyyz[k] = -g_y_0_xy_xxyyz[k] * ab_x + g_y_0_xy_xxxyyz[k];

                g_y_0_xxy_xxyzz[k] = -g_y_0_xy_xxyzz[k] * ab_x + g_y_0_xy_xxxyzz[k];

                g_y_0_xxy_xxzzz[k] = -g_y_0_xy_xxzzz[k] * ab_x + g_y_0_xy_xxxzzz[k];

                g_y_0_xxy_xyyyy[k] = -g_y_0_xy_xyyyy[k] * ab_x + g_y_0_xy_xxyyyy[k];

                g_y_0_xxy_xyyyz[k] = -g_y_0_xy_xyyyz[k] * ab_x + g_y_0_xy_xxyyyz[k];

                g_y_0_xxy_xyyzz[k] = -g_y_0_xy_xyyzz[k] * ab_x + g_y_0_xy_xxyyzz[k];

                g_y_0_xxy_xyzzz[k] = -g_y_0_xy_xyzzz[k] * ab_x + g_y_0_xy_xxyzzz[k];

                g_y_0_xxy_xzzzz[k] = -g_y_0_xy_xzzzz[k] * ab_x + g_y_0_xy_xxzzzz[k];

                g_y_0_xxy_yyyyy[k] = -g_y_0_xy_yyyyy[k] * ab_x + g_y_0_xy_xyyyyy[k];

                g_y_0_xxy_yyyyz[k] = -g_y_0_xy_yyyyz[k] * ab_x + g_y_0_xy_xyyyyz[k];

                g_y_0_xxy_yyyzz[k] = -g_y_0_xy_yyyzz[k] * ab_x + g_y_0_xy_xyyyzz[k];

                g_y_0_xxy_yyzzz[k] = -g_y_0_xy_yyzzz[k] * ab_x + g_y_0_xy_xyyzzz[k];

                g_y_0_xxy_yzzzz[k] = -g_y_0_xy_yzzzz[k] * ab_x + g_y_0_xy_xyzzzz[k];

                g_y_0_xxy_zzzzz[k] = -g_y_0_xy_zzzzz[k] * ab_x + g_y_0_xy_xzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxz_xxxxx, g_y_0_xxz_xxxxy, g_y_0_xxz_xxxxz, g_y_0_xxz_xxxyy, g_y_0_xxz_xxxyz, g_y_0_xxz_xxxzz, g_y_0_xxz_xxyyy, g_y_0_xxz_xxyyz, g_y_0_xxz_xxyzz, g_y_0_xxz_xxzzz, g_y_0_xxz_xyyyy, g_y_0_xxz_xyyyz, g_y_0_xxz_xyyzz, g_y_0_xxz_xyzzz, g_y_0_xxz_xzzzz, g_y_0_xxz_yyyyy, g_y_0_xxz_yyyyz, g_y_0_xxz_yyyzz, g_y_0_xxz_yyzzz, g_y_0_xxz_yzzzz, g_y_0_xxz_zzzzz, g_y_0_xz_xxxxx, g_y_0_xz_xxxxxx, g_y_0_xz_xxxxxy, g_y_0_xz_xxxxxz, g_y_0_xz_xxxxy, g_y_0_xz_xxxxyy, g_y_0_xz_xxxxyz, g_y_0_xz_xxxxz, g_y_0_xz_xxxxzz, g_y_0_xz_xxxyy, g_y_0_xz_xxxyyy, g_y_0_xz_xxxyyz, g_y_0_xz_xxxyz, g_y_0_xz_xxxyzz, g_y_0_xz_xxxzz, g_y_0_xz_xxxzzz, g_y_0_xz_xxyyy, g_y_0_xz_xxyyyy, g_y_0_xz_xxyyyz, g_y_0_xz_xxyyz, g_y_0_xz_xxyyzz, g_y_0_xz_xxyzz, g_y_0_xz_xxyzzz, g_y_0_xz_xxzzz, g_y_0_xz_xxzzzz, g_y_0_xz_xyyyy, g_y_0_xz_xyyyyy, g_y_0_xz_xyyyyz, g_y_0_xz_xyyyz, g_y_0_xz_xyyyzz, g_y_0_xz_xyyzz, g_y_0_xz_xyyzzz, g_y_0_xz_xyzzz, g_y_0_xz_xyzzzz, g_y_0_xz_xzzzz, g_y_0_xz_xzzzzz, g_y_0_xz_yyyyy, g_y_0_xz_yyyyz, g_y_0_xz_yyyzz, g_y_0_xz_yyzzz, g_y_0_xz_yzzzz, g_y_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_xxxxx[k] = -g_y_0_xz_xxxxx[k] * ab_x + g_y_0_xz_xxxxxx[k];

                g_y_0_xxz_xxxxy[k] = -g_y_0_xz_xxxxy[k] * ab_x + g_y_0_xz_xxxxxy[k];

                g_y_0_xxz_xxxxz[k] = -g_y_0_xz_xxxxz[k] * ab_x + g_y_0_xz_xxxxxz[k];

                g_y_0_xxz_xxxyy[k] = -g_y_0_xz_xxxyy[k] * ab_x + g_y_0_xz_xxxxyy[k];

                g_y_0_xxz_xxxyz[k] = -g_y_0_xz_xxxyz[k] * ab_x + g_y_0_xz_xxxxyz[k];

                g_y_0_xxz_xxxzz[k] = -g_y_0_xz_xxxzz[k] * ab_x + g_y_0_xz_xxxxzz[k];

                g_y_0_xxz_xxyyy[k] = -g_y_0_xz_xxyyy[k] * ab_x + g_y_0_xz_xxxyyy[k];

                g_y_0_xxz_xxyyz[k] = -g_y_0_xz_xxyyz[k] * ab_x + g_y_0_xz_xxxyyz[k];

                g_y_0_xxz_xxyzz[k] = -g_y_0_xz_xxyzz[k] * ab_x + g_y_0_xz_xxxyzz[k];

                g_y_0_xxz_xxzzz[k] = -g_y_0_xz_xxzzz[k] * ab_x + g_y_0_xz_xxxzzz[k];

                g_y_0_xxz_xyyyy[k] = -g_y_0_xz_xyyyy[k] * ab_x + g_y_0_xz_xxyyyy[k];

                g_y_0_xxz_xyyyz[k] = -g_y_0_xz_xyyyz[k] * ab_x + g_y_0_xz_xxyyyz[k];

                g_y_0_xxz_xyyzz[k] = -g_y_0_xz_xyyzz[k] * ab_x + g_y_0_xz_xxyyzz[k];

                g_y_0_xxz_xyzzz[k] = -g_y_0_xz_xyzzz[k] * ab_x + g_y_0_xz_xxyzzz[k];

                g_y_0_xxz_xzzzz[k] = -g_y_0_xz_xzzzz[k] * ab_x + g_y_0_xz_xxzzzz[k];

                g_y_0_xxz_yyyyy[k] = -g_y_0_xz_yyyyy[k] * ab_x + g_y_0_xz_xyyyyy[k];

                g_y_0_xxz_yyyyz[k] = -g_y_0_xz_yyyyz[k] * ab_x + g_y_0_xz_xyyyyz[k];

                g_y_0_xxz_yyyzz[k] = -g_y_0_xz_yyyzz[k] * ab_x + g_y_0_xz_xyyyzz[k];

                g_y_0_xxz_yyzzz[k] = -g_y_0_xz_yyzzz[k] * ab_x + g_y_0_xz_xyyzzz[k];

                g_y_0_xxz_yzzzz[k] = -g_y_0_xz_yzzzz[k] * ab_x + g_y_0_xz_xyzzzz[k];

                g_y_0_xxz_zzzzz[k] = -g_y_0_xz_zzzzz[k] * ab_x + g_y_0_xz_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_xxxxx, g_y_0_xyy_xxxxy, g_y_0_xyy_xxxxz, g_y_0_xyy_xxxyy, g_y_0_xyy_xxxyz, g_y_0_xyy_xxxzz, g_y_0_xyy_xxyyy, g_y_0_xyy_xxyyz, g_y_0_xyy_xxyzz, g_y_0_xyy_xxzzz, g_y_0_xyy_xyyyy, g_y_0_xyy_xyyyz, g_y_0_xyy_xyyzz, g_y_0_xyy_xyzzz, g_y_0_xyy_xzzzz, g_y_0_xyy_yyyyy, g_y_0_xyy_yyyyz, g_y_0_xyy_yyyzz, g_y_0_xyy_yyzzz, g_y_0_xyy_yzzzz, g_y_0_xyy_zzzzz, g_y_0_yy_xxxxx, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyz, g_y_0_yy_yyyzz, g_y_0_yy_yyzzz, g_y_0_yy_yzzzz, g_y_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_xxxxx[k] = -g_y_0_yy_xxxxx[k] * ab_x + g_y_0_yy_xxxxxx[k];

                g_y_0_xyy_xxxxy[k] = -g_y_0_yy_xxxxy[k] * ab_x + g_y_0_yy_xxxxxy[k];

                g_y_0_xyy_xxxxz[k] = -g_y_0_yy_xxxxz[k] * ab_x + g_y_0_yy_xxxxxz[k];

                g_y_0_xyy_xxxyy[k] = -g_y_0_yy_xxxyy[k] * ab_x + g_y_0_yy_xxxxyy[k];

                g_y_0_xyy_xxxyz[k] = -g_y_0_yy_xxxyz[k] * ab_x + g_y_0_yy_xxxxyz[k];

                g_y_0_xyy_xxxzz[k] = -g_y_0_yy_xxxzz[k] * ab_x + g_y_0_yy_xxxxzz[k];

                g_y_0_xyy_xxyyy[k] = -g_y_0_yy_xxyyy[k] * ab_x + g_y_0_yy_xxxyyy[k];

                g_y_0_xyy_xxyyz[k] = -g_y_0_yy_xxyyz[k] * ab_x + g_y_0_yy_xxxyyz[k];

                g_y_0_xyy_xxyzz[k] = -g_y_0_yy_xxyzz[k] * ab_x + g_y_0_yy_xxxyzz[k];

                g_y_0_xyy_xxzzz[k] = -g_y_0_yy_xxzzz[k] * ab_x + g_y_0_yy_xxxzzz[k];

                g_y_0_xyy_xyyyy[k] = -g_y_0_yy_xyyyy[k] * ab_x + g_y_0_yy_xxyyyy[k];

                g_y_0_xyy_xyyyz[k] = -g_y_0_yy_xyyyz[k] * ab_x + g_y_0_yy_xxyyyz[k];

                g_y_0_xyy_xyyzz[k] = -g_y_0_yy_xyyzz[k] * ab_x + g_y_0_yy_xxyyzz[k];

                g_y_0_xyy_xyzzz[k] = -g_y_0_yy_xyzzz[k] * ab_x + g_y_0_yy_xxyzzz[k];

                g_y_0_xyy_xzzzz[k] = -g_y_0_yy_xzzzz[k] * ab_x + g_y_0_yy_xxzzzz[k];

                g_y_0_xyy_yyyyy[k] = -g_y_0_yy_yyyyy[k] * ab_x + g_y_0_yy_xyyyyy[k];

                g_y_0_xyy_yyyyz[k] = -g_y_0_yy_yyyyz[k] * ab_x + g_y_0_yy_xyyyyz[k];

                g_y_0_xyy_yyyzz[k] = -g_y_0_yy_yyyzz[k] * ab_x + g_y_0_yy_xyyyzz[k];

                g_y_0_xyy_yyzzz[k] = -g_y_0_yy_yyzzz[k] * ab_x + g_y_0_yy_xyyzzz[k];

                g_y_0_xyy_yzzzz[k] = -g_y_0_yy_yzzzz[k] * ab_x + g_y_0_yy_xyzzzz[k];

                g_y_0_xyy_zzzzz[k] = -g_y_0_yy_zzzzz[k] * ab_x + g_y_0_yy_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyz_xxxxx, g_y_0_xyz_xxxxy, g_y_0_xyz_xxxxz, g_y_0_xyz_xxxyy, g_y_0_xyz_xxxyz, g_y_0_xyz_xxxzz, g_y_0_xyz_xxyyy, g_y_0_xyz_xxyyz, g_y_0_xyz_xxyzz, g_y_0_xyz_xxzzz, g_y_0_xyz_xyyyy, g_y_0_xyz_xyyyz, g_y_0_xyz_xyyzz, g_y_0_xyz_xyzzz, g_y_0_xyz_xzzzz, g_y_0_xyz_yyyyy, g_y_0_xyz_yyyyz, g_y_0_xyz_yyyzz, g_y_0_xyz_yyzzz, g_y_0_xyz_yzzzz, g_y_0_xyz_zzzzz, g_y_0_yz_xxxxx, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxy, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxyy, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxyyy, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xyyyy, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_yyyyy, g_y_0_yz_yyyyz, g_y_0_yz_yyyzz, g_y_0_yz_yyzzz, g_y_0_yz_yzzzz, g_y_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_xxxxx[k] = -g_y_0_yz_xxxxx[k] * ab_x + g_y_0_yz_xxxxxx[k];

                g_y_0_xyz_xxxxy[k] = -g_y_0_yz_xxxxy[k] * ab_x + g_y_0_yz_xxxxxy[k];

                g_y_0_xyz_xxxxz[k] = -g_y_0_yz_xxxxz[k] * ab_x + g_y_0_yz_xxxxxz[k];

                g_y_0_xyz_xxxyy[k] = -g_y_0_yz_xxxyy[k] * ab_x + g_y_0_yz_xxxxyy[k];

                g_y_0_xyz_xxxyz[k] = -g_y_0_yz_xxxyz[k] * ab_x + g_y_0_yz_xxxxyz[k];

                g_y_0_xyz_xxxzz[k] = -g_y_0_yz_xxxzz[k] * ab_x + g_y_0_yz_xxxxzz[k];

                g_y_0_xyz_xxyyy[k] = -g_y_0_yz_xxyyy[k] * ab_x + g_y_0_yz_xxxyyy[k];

                g_y_0_xyz_xxyyz[k] = -g_y_0_yz_xxyyz[k] * ab_x + g_y_0_yz_xxxyyz[k];

                g_y_0_xyz_xxyzz[k] = -g_y_0_yz_xxyzz[k] * ab_x + g_y_0_yz_xxxyzz[k];

                g_y_0_xyz_xxzzz[k] = -g_y_0_yz_xxzzz[k] * ab_x + g_y_0_yz_xxxzzz[k];

                g_y_0_xyz_xyyyy[k] = -g_y_0_yz_xyyyy[k] * ab_x + g_y_0_yz_xxyyyy[k];

                g_y_0_xyz_xyyyz[k] = -g_y_0_yz_xyyyz[k] * ab_x + g_y_0_yz_xxyyyz[k];

                g_y_0_xyz_xyyzz[k] = -g_y_0_yz_xyyzz[k] * ab_x + g_y_0_yz_xxyyzz[k];

                g_y_0_xyz_xyzzz[k] = -g_y_0_yz_xyzzz[k] * ab_x + g_y_0_yz_xxyzzz[k];

                g_y_0_xyz_xzzzz[k] = -g_y_0_yz_xzzzz[k] * ab_x + g_y_0_yz_xxzzzz[k];

                g_y_0_xyz_yyyyy[k] = -g_y_0_yz_yyyyy[k] * ab_x + g_y_0_yz_xyyyyy[k];

                g_y_0_xyz_yyyyz[k] = -g_y_0_yz_yyyyz[k] * ab_x + g_y_0_yz_xyyyyz[k];

                g_y_0_xyz_yyyzz[k] = -g_y_0_yz_yyyzz[k] * ab_x + g_y_0_yz_xyyyzz[k];

                g_y_0_xyz_yyzzz[k] = -g_y_0_yz_yyzzz[k] * ab_x + g_y_0_yz_xyyzzz[k];

                g_y_0_xyz_yzzzz[k] = -g_y_0_yz_yzzzz[k] * ab_x + g_y_0_yz_xyzzzz[k];

                g_y_0_xyz_zzzzz[k] = -g_y_0_yz_zzzzz[k] * ab_x + g_y_0_yz_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzz_xxxxx, g_y_0_xzz_xxxxy, g_y_0_xzz_xxxxz, g_y_0_xzz_xxxyy, g_y_0_xzz_xxxyz, g_y_0_xzz_xxxzz, g_y_0_xzz_xxyyy, g_y_0_xzz_xxyyz, g_y_0_xzz_xxyzz, g_y_0_xzz_xxzzz, g_y_0_xzz_xyyyy, g_y_0_xzz_xyyyz, g_y_0_xzz_xyyzz, g_y_0_xzz_xyzzz, g_y_0_xzz_xzzzz, g_y_0_xzz_yyyyy, g_y_0_xzz_yyyyz, g_y_0_xzz_yyyzz, g_y_0_xzz_yyzzz, g_y_0_xzz_yzzzz, g_y_0_xzz_zzzzz, g_y_0_zz_xxxxx, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxy, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxyy, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxyyy, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xyyyy, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_yyyyy, g_y_0_zz_yyyyz, g_y_0_zz_yyyzz, g_y_0_zz_yyzzz, g_y_0_zz_yzzzz, g_y_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_xxxxx[k] = -g_y_0_zz_xxxxx[k] * ab_x + g_y_0_zz_xxxxxx[k];

                g_y_0_xzz_xxxxy[k] = -g_y_0_zz_xxxxy[k] * ab_x + g_y_0_zz_xxxxxy[k];

                g_y_0_xzz_xxxxz[k] = -g_y_0_zz_xxxxz[k] * ab_x + g_y_0_zz_xxxxxz[k];

                g_y_0_xzz_xxxyy[k] = -g_y_0_zz_xxxyy[k] * ab_x + g_y_0_zz_xxxxyy[k];

                g_y_0_xzz_xxxyz[k] = -g_y_0_zz_xxxyz[k] * ab_x + g_y_0_zz_xxxxyz[k];

                g_y_0_xzz_xxxzz[k] = -g_y_0_zz_xxxzz[k] * ab_x + g_y_0_zz_xxxxzz[k];

                g_y_0_xzz_xxyyy[k] = -g_y_0_zz_xxyyy[k] * ab_x + g_y_0_zz_xxxyyy[k];

                g_y_0_xzz_xxyyz[k] = -g_y_0_zz_xxyyz[k] * ab_x + g_y_0_zz_xxxyyz[k];

                g_y_0_xzz_xxyzz[k] = -g_y_0_zz_xxyzz[k] * ab_x + g_y_0_zz_xxxyzz[k];

                g_y_0_xzz_xxzzz[k] = -g_y_0_zz_xxzzz[k] * ab_x + g_y_0_zz_xxxzzz[k];

                g_y_0_xzz_xyyyy[k] = -g_y_0_zz_xyyyy[k] * ab_x + g_y_0_zz_xxyyyy[k];

                g_y_0_xzz_xyyyz[k] = -g_y_0_zz_xyyyz[k] * ab_x + g_y_0_zz_xxyyyz[k];

                g_y_0_xzz_xyyzz[k] = -g_y_0_zz_xyyzz[k] * ab_x + g_y_0_zz_xxyyzz[k];

                g_y_0_xzz_xyzzz[k] = -g_y_0_zz_xyzzz[k] * ab_x + g_y_0_zz_xxyzzz[k];

                g_y_0_xzz_xzzzz[k] = -g_y_0_zz_xzzzz[k] * ab_x + g_y_0_zz_xxzzzz[k];

                g_y_0_xzz_yyyyy[k] = -g_y_0_zz_yyyyy[k] * ab_x + g_y_0_zz_xyyyyy[k];

                g_y_0_xzz_yyyyz[k] = -g_y_0_zz_yyyyz[k] * ab_x + g_y_0_zz_xyyyyz[k];

                g_y_0_xzz_yyyzz[k] = -g_y_0_zz_yyyzz[k] * ab_x + g_y_0_zz_xyyyzz[k];

                g_y_0_xzz_yyzzz[k] = -g_y_0_zz_yyzzz[k] * ab_x + g_y_0_zz_xyyzzz[k];

                g_y_0_xzz_yzzzz[k] = -g_y_0_zz_yzzzz[k] * ab_x + g_y_0_zz_xyzzzz[k];

                g_y_0_xzz_zzzzz[k] = -g_y_0_zz_zzzzz[k] * ab_x + g_y_0_zz_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxxx, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzz, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzzz, g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz, g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyzz, g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyz, g_yy_xyyzz, g_yy_xyzzz, g_yy_xzzzz, g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz, g_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_xxxxx[k] = -g_yy_xxxxx[k] - g_y_0_yy_xxxxx[k] * ab_y + g_y_0_yy_xxxxxy[k];

                g_y_0_yyy_xxxxy[k] = -g_yy_xxxxy[k] - g_y_0_yy_xxxxy[k] * ab_y + g_y_0_yy_xxxxyy[k];

                g_y_0_yyy_xxxxz[k] = -g_yy_xxxxz[k] - g_y_0_yy_xxxxz[k] * ab_y + g_y_0_yy_xxxxyz[k];

                g_y_0_yyy_xxxyy[k] = -g_yy_xxxyy[k] - g_y_0_yy_xxxyy[k] * ab_y + g_y_0_yy_xxxyyy[k];

                g_y_0_yyy_xxxyz[k] = -g_yy_xxxyz[k] - g_y_0_yy_xxxyz[k] * ab_y + g_y_0_yy_xxxyyz[k];

                g_y_0_yyy_xxxzz[k] = -g_yy_xxxzz[k] - g_y_0_yy_xxxzz[k] * ab_y + g_y_0_yy_xxxyzz[k];

                g_y_0_yyy_xxyyy[k] = -g_yy_xxyyy[k] - g_y_0_yy_xxyyy[k] * ab_y + g_y_0_yy_xxyyyy[k];

                g_y_0_yyy_xxyyz[k] = -g_yy_xxyyz[k] - g_y_0_yy_xxyyz[k] * ab_y + g_y_0_yy_xxyyyz[k];

                g_y_0_yyy_xxyzz[k] = -g_yy_xxyzz[k] - g_y_0_yy_xxyzz[k] * ab_y + g_y_0_yy_xxyyzz[k];

                g_y_0_yyy_xxzzz[k] = -g_yy_xxzzz[k] - g_y_0_yy_xxzzz[k] * ab_y + g_y_0_yy_xxyzzz[k];

                g_y_0_yyy_xyyyy[k] = -g_yy_xyyyy[k] - g_y_0_yy_xyyyy[k] * ab_y + g_y_0_yy_xyyyyy[k];

                g_y_0_yyy_xyyyz[k] = -g_yy_xyyyz[k] - g_y_0_yy_xyyyz[k] * ab_y + g_y_0_yy_xyyyyz[k];

                g_y_0_yyy_xyyzz[k] = -g_yy_xyyzz[k] - g_y_0_yy_xyyzz[k] * ab_y + g_y_0_yy_xyyyzz[k];

                g_y_0_yyy_xyzzz[k] = -g_yy_xyzzz[k] - g_y_0_yy_xyzzz[k] * ab_y + g_y_0_yy_xyyzzz[k];

                g_y_0_yyy_xzzzz[k] = -g_yy_xzzzz[k] - g_y_0_yy_xzzzz[k] * ab_y + g_y_0_yy_xyzzzz[k];

                g_y_0_yyy_yyyyy[k] = -g_yy_yyyyy[k] - g_y_0_yy_yyyyy[k] * ab_y + g_y_0_yy_yyyyyy[k];

                g_y_0_yyy_yyyyz[k] = -g_yy_yyyyz[k] - g_y_0_yy_yyyyz[k] * ab_y + g_y_0_yy_yyyyyz[k];

                g_y_0_yyy_yyyzz[k] = -g_yy_yyyzz[k] - g_y_0_yy_yyyzz[k] * ab_y + g_y_0_yy_yyyyzz[k];

                g_y_0_yyy_yyzzz[k] = -g_yy_yyzzz[k] - g_y_0_yy_yyzzz[k] * ab_y + g_y_0_yy_yyyzzz[k];

                g_y_0_yyy_yzzzz[k] = -g_yy_yzzzz[k] - g_y_0_yy_yzzzz[k] * ab_y + g_y_0_yy_yyzzzz[k];

                g_y_0_yyy_zzzzz[k] = -g_yy_zzzzz[k] - g_y_0_yy_zzzzz[k] * ab_y + g_y_0_yy_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxxx, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzz, g_y_0_yy_zzzzzz, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyyy, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_xxxxx[k] = -g_y_0_yy_xxxxx[k] * ab_z + g_y_0_yy_xxxxxz[k];

                g_y_0_yyz_xxxxy[k] = -g_y_0_yy_xxxxy[k] * ab_z + g_y_0_yy_xxxxyz[k];

                g_y_0_yyz_xxxxz[k] = -g_y_0_yy_xxxxz[k] * ab_z + g_y_0_yy_xxxxzz[k];

                g_y_0_yyz_xxxyy[k] = -g_y_0_yy_xxxyy[k] * ab_z + g_y_0_yy_xxxyyz[k];

                g_y_0_yyz_xxxyz[k] = -g_y_0_yy_xxxyz[k] * ab_z + g_y_0_yy_xxxyzz[k];

                g_y_0_yyz_xxxzz[k] = -g_y_0_yy_xxxzz[k] * ab_z + g_y_0_yy_xxxzzz[k];

                g_y_0_yyz_xxyyy[k] = -g_y_0_yy_xxyyy[k] * ab_z + g_y_0_yy_xxyyyz[k];

                g_y_0_yyz_xxyyz[k] = -g_y_0_yy_xxyyz[k] * ab_z + g_y_0_yy_xxyyzz[k];

                g_y_0_yyz_xxyzz[k] = -g_y_0_yy_xxyzz[k] * ab_z + g_y_0_yy_xxyzzz[k];

                g_y_0_yyz_xxzzz[k] = -g_y_0_yy_xxzzz[k] * ab_z + g_y_0_yy_xxzzzz[k];

                g_y_0_yyz_xyyyy[k] = -g_y_0_yy_xyyyy[k] * ab_z + g_y_0_yy_xyyyyz[k];

                g_y_0_yyz_xyyyz[k] = -g_y_0_yy_xyyyz[k] * ab_z + g_y_0_yy_xyyyzz[k];

                g_y_0_yyz_xyyzz[k] = -g_y_0_yy_xyyzz[k] * ab_z + g_y_0_yy_xyyzzz[k];

                g_y_0_yyz_xyzzz[k] = -g_y_0_yy_xyzzz[k] * ab_z + g_y_0_yy_xyzzzz[k];

                g_y_0_yyz_xzzzz[k] = -g_y_0_yy_xzzzz[k] * ab_z + g_y_0_yy_xzzzzz[k];

                g_y_0_yyz_yyyyy[k] = -g_y_0_yy_yyyyy[k] * ab_z + g_y_0_yy_yyyyyz[k];

                g_y_0_yyz_yyyyz[k] = -g_y_0_yy_yyyyz[k] * ab_z + g_y_0_yy_yyyyzz[k];

                g_y_0_yyz_yyyzz[k] = -g_y_0_yy_yyyzz[k] * ab_z + g_y_0_yy_yyyzzz[k];

                g_y_0_yyz_yyzzz[k] = -g_y_0_yy_yyzzz[k] * ab_z + g_y_0_yy_yyzzzz[k];

                g_y_0_yyz_yzzzz[k] = -g_y_0_yy_yzzzz[k] * ab_z + g_y_0_yy_yzzzzz[k];

                g_y_0_yyz_zzzzz[k] = -g_y_0_yy_zzzzz[k] * ab_z + g_y_0_yy_zzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yz_xxxxx, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxy, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxyy, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxyyy, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xyyyy, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_yyyyy, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_zzzzz, g_y_0_yz_zzzzzz, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyyy, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_xxxxx[k] = -g_y_0_yz_xxxxx[k] * ab_z + g_y_0_yz_xxxxxz[k];

                g_y_0_yzz_xxxxy[k] = -g_y_0_yz_xxxxy[k] * ab_z + g_y_0_yz_xxxxyz[k];

                g_y_0_yzz_xxxxz[k] = -g_y_0_yz_xxxxz[k] * ab_z + g_y_0_yz_xxxxzz[k];

                g_y_0_yzz_xxxyy[k] = -g_y_0_yz_xxxyy[k] * ab_z + g_y_0_yz_xxxyyz[k];

                g_y_0_yzz_xxxyz[k] = -g_y_0_yz_xxxyz[k] * ab_z + g_y_0_yz_xxxyzz[k];

                g_y_0_yzz_xxxzz[k] = -g_y_0_yz_xxxzz[k] * ab_z + g_y_0_yz_xxxzzz[k];

                g_y_0_yzz_xxyyy[k] = -g_y_0_yz_xxyyy[k] * ab_z + g_y_0_yz_xxyyyz[k];

                g_y_0_yzz_xxyyz[k] = -g_y_0_yz_xxyyz[k] * ab_z + g_y_0_yz_xxyyzz[k];

                g_y_0_yzz_xxyzz[k] = -g_y_0_yz_xxyzz[k] * ab_z + g_y_0_yz_xxyzzz[k];

                g_y_0_yzz_xxzzz[k] = -g_y_0_yz_xxzzz[k] * ab_z + g_y_0_yz_xxzzzz[k];

                g_y_0_yzz_xyyyy[k] = -g_y_0_yz_xyyyy[k] * ab_z + g_y_0_yz_xyyyyz[k];

                g_y_0_yzz_xyyyz[k] = -g_y_0_yz_xyyyz[k] * ab_z + g_y_0_yz_xyyyzz[k];

                g_y_0_yzz_xyyzz[k] = -g_y_0_yz_xyyzz[k] * ab_z + g_y_0_yz_xyyzzz[k];

                g_y_0_yzz_xyzzz[k] = -g_y_0_yz_xyzzz[k] * ab_z + g_y_0_yz_xyzzzz[k];

                g_y_0_yzz_xzzzz[k] = -g_y_0_yz_xzzzz[k] * ab_z + g_y_0_yz_xzzzzz[k];

                g_y_0_yzz_yyyyy[k] = -g_y_0_yz_yyyyy[k] * ab_z + g_y_0_yz_yyyyyz[k];

                g_y_0_yzz_yyyyz[k] = -g_y_0_yz_yyyyz[k] * ab_z + g_y_0_yz_yyyyzz[k];

                g_y_0_yzz_yyyzz[k] = -g_y_0_yz_yyyzz[k] * ab_z + g_y_0_yz_yyyzzz[k];

                g_y_0_yzz_yyzzz[k] = -g_y_0_yz_yyzzz[k] * ab_z + g_y_0_yz_yyzzzz[k];

                g_y_0_yzz_yzzzz[k] = -g_y_0_yz_yzzzz[k] * ab_z + g_y_0_yz_yzzzzz[k];

                g_y_0_yzz_zzzzz[k] = -g_y_0_yz_zzzzz[k] * ab_z + g_y_0_yz_zzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xxxxx, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxy, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxyy, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxyyy, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xyyyy, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_yyyyy, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_zzzzz, g_y_0_zz_zzzzzz, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyyy, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_xxxxx[k] = -g_y_0_zz_xxxxx[k] * ab_z + g_y_0_zz_xxxxxz[k];

                g_y_0_zzz_xxxxy[k] = -g_y_0_zz_xxxxy[k] * ab_z + g_y_0_zz_xxxxyz[k];

                g_y_0_zzz_xxxxz[k] = -g_y_0_zz_xxxxz[k] * ab_z + g_y_0_zz_xxxxzz[k];

                g_y_0_zzz_xxxyy[k] = -g_y_0_zz_xxxyy[k] * ab_z + g_y_0_zz_xxxyyz[k];

                g_y_0_zzz_xxxyz[k] = -g_y_0_zz_xxxyz[k] * ab_z + g_y_0_zz_xxxyzz[k];

                g_y_0_zzz_xxxzz[k] = -g_y_0_zz_xxxzz[k] * ab_z + g_y_0_zz_xxxzzz[k];

                g_y_0_zzz_xxyyy[k] = -g_y_0_zz_xxyyy[k] * ab_z + g_y_0_zz_xxyyyz[k];

                g_y_0_zzz_xxyyz[k] = -g_y_0_zz_xxyyz[k] * ab_z + g_y_0_zz_xxyyzz[k];

                g_y_0_zzz_xxyzz[k] = -g_y_0_zz_xxyzz[k] * ab_z + g_y_0_zz_xxyzzz[k];

                g_y_0_zzz_xxzzz[k] = -g_y_0_zz_xxzzz[k] * ab_z + g_y_0_zz_xxzzzz[k];

                g_y_0_zzz_xyyyy[k] = -g_y_0_zz_xyyyy[k] * ab_z + g_y_0_zz_xyyyyz[k];

                g_y_0_zzz_xyyyz[k] = -g_y_0_zz_xyyyz[k] * ab_z + g_y_0_zz_xyyyzz[k];

                g_y_0_zzz_xyyzz[k] = -g_y_0_zz_xyyzz[k] * ab_z + g_y_0_zz_xyyzzz[k];

                g_y_0_zzz_xyzzz[k] = -g_y_0_zz_xyzzz[k] * ab_z + g_y_0_zz_xyzzzz[k];

                g_y_0_zzz_xzzzz[k] = -g_y_0_zz_xzzzz[k] * ab_z + g_y_0_zz_xzzzzz[k];

                g_y_0_zzz_yyyyy[k] = -g_y_0_zz_yyyyy[k] * ab_z + g_y_0_zz_yyyyyz[k];

                g_y_0_zzz_yyyyz[k] = -g_y_0_zz_yyyyz[k] * ab_z + g_y_0_zz_yyyyzz[k];

                g_y_0_zzz_yyyzz[k] = -g_y_0_zz_yyyzz[k] * ab_z + g_y_0_zz_yyyzzz[k];

                g_y_0_zzz_yyzzz[k] = -g_y_0_zz_yyzzz[k] * ab_z + g_y_0_zz_yyzzzz[k];

                g_y_0_zzz_yzzzz[k] = -g_y_0_zz_yzzzz[k] * ab_z + g_y_0_zz_yzzzzz[k];

                g_y_0_zzz_zzzzz[k] = -g_y_0_zz_zzzzz[k] * ab_z + g_y_0_zz_zzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_xxxxx, g_z_0_xx_xxxxxx, g_z_0_xx_xxxxxy, g_z_0_xx_xxxxxz, g_z_0_xx_xxxxy, g_z_0_xx_xxxxyy, g_z_0_xx_xxxxyz, g_z_0_xx_xxxxz, g_z_0_xx_xxxxzz, g_z_0_xx_xxxyy, g_z_0_xx_xxxyyy, g_z_0_xx_xxxyyz, g_z_0_xx_xxxyz, g_z_0_xx_xxxyzz, g_z_0_xx_xxxzz, g_z_0_xx_xxxzzz, g_z_0_xx_xxyyy, g_z_0_xx_xxyyyy, g_z_0_xx_xxyyyz, g_z_0_xx_xxyyz, g_z_0_xx_xxyyzz, g_z_0_xx_xxyzz, g_z_0_xx_xxyzzz, g_z_0_xx_xxzzz, g_z_0_xx_xxzzzz, g_z_0_xx_xyyyy, g_z_0_xx_xyyyyy, g_z_0_xx_xyyyyz, g_z_0_xx_xyyyz, g_z_0_xx_xyyyzz, g_z_0_xx_xyyzz, g_z_0_xx_xyyzzz, g_z_0_xx_xyzzz, g_z_0_xx_xyzzzz, g_z_0_xx_xzzzz, g_z_0_xx_xzzzzz, g_z_0_xx_yyyyy, g_z_0_xx_yyyyz, g_z_0_xx_yyyzz, g_z_0_xx_yyzzz, g_z_0_xx_yzzzz, g_z_0_xx_zzzzz, g_z_0_xxx_xxxxx, g_z_0_xxx_xxxxy, g_z_0_xxx_xxxxz, g_z_0_xxx_xxxyy, g_z_0_xxx_xxxyz, g_z_0_xxx_xxxzz, g_z_0_xxx_xxyyy, g_z_0_xxx_xxyyz, g_z_0_xxx_xxyzz, g_z_0_xxx_xxzzz, g_z_0_xxx_xyyyy, g_z_0_xxx_xyyyz, g_z_0_xxx_xyyzz, g_z_0_xxx_xyzzz, g_z_0_xxx_xzzzz, g_z_0_xxx_yyyyy, g_z_0_xxx_yyyyz, g_z_0_xxx_yyyzz, g_z_0_xxx_yyzzz, g_z_0_xxx_yzzzz, g_z_0_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_xxxxx[k] = -g_z_0_xx_xxxxx[k] * ab_x + g_z_0_xx_xxxxxx[k];

                g_z_0_xxx_xxxxy[k] = -g_z_0_xx_xxxxy[k] * ab_x + g_z_0_xx_xxxxxy[k];

                g_z_0_xxx_xxxxz[k] = -g_z_0_xx_xxxxz[k] * ab_x + g_z_0_xx_xxxxxz[k];

                g_z_0_xxx_xxxyy[k] = -g_z_0_xx_xxxyy[k] * ab_x + g_z_0_xx_xxxxyy[k];

                g_z_0_xxx_xxxyz[k] = -g_z_0_xx_xxxyz[k] * ab_x + g_z_0_xx_xxxxyz[k];

                g_z_0_xxx_xxxzz[k] = -g_z_0_xx_xxxzz[k] * ab_x + g_z_0_xx_xxxxzz[k];

                g_z_0_xxx_xxyyy[k] = -g_z_0_xx_xxyyy[k] * ab_x + g_z_0_xx_xxxyyy[k];

                g_z_0_xxx_xxyyz[k] = -g_z_0_xx_xxyyz[k] * ab_x + g_z_0_xx_xxxyyz[k];

                g_z_0_xxx_xxyzz[k] = -g_z_0_xx_xxyzz[k] * ab_x + g_z_0_xx_xxxyzz[k];

                g_z_0_xxx_xxzzz[k] = -g_z_0_xx_xxzzz[k] * ab_x + g_z_0_xx_xxxzzz[k];

                g_z_0_xxx_xyyyy[k] = -g_z_0_xx_xyyyy[k] * ab_x + g_z_0_xx_xxyyyy[k];

                g_z_0_xxx_xyyyz[k] = -g_z_0_xx_xyyyz[k] * ab_x + g_z_0_xx_xxyyyz[k];

                g_z_0_xxx_xyyzz[k] = -g_z_0_xx_xyyzz[k] * ab_x + g_z_0_xx_xxyyzz[k];

                g_z_0_xxx_xyzzz[k] = -g_z_0_xx_xyzzz[k] * ab_x + g_z_0_xx_xxyzzz[k];

                g_z_0_xxx_xzzzz[k] = -g_z_0_xx_xzzzz[k] * ab_x + g_z_0_xx_xxzzzz[k];

                g_z_0_xxx_yyyyy[k] = -g_z_0_xx_yyyyy[k] * ab_x + g_z_0_xx_xyyyyy[k];

                g_z_0_xxx_yyyyz[k] = -g_z_0_xx_yyyyz[k] * ab_x + g_z_0_xx_xyyyyz[k];

                g_z_0_xxx_yyyzz[k] = -g_z_0_xx_yyyzz[k] * ab_x + g_z_0_xx_xyyyzz[k];

                g_z_0_xxx_yyzzz[k] = -g_z_0_xx_yyzzz[k] * ab_x + g_z_0_xx_xyyzzz[k];

                g_z_0_xxx_yzzzz[k] = -g_z_0_xx_yzzzz[k] * ab_x + g_z_0_xx_xyzzzz[k];

                g_z_0_xxx_zzzzz[k] = -g_z_0_xx_zzzzz[k] * ab_x + g_z_0_xx_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxy_xxxxx, g_z_0_xxy_xxxxy, g_z_0_xxy_xxxxz, g_z_0_xxy_xxxyy, g_z_0_xxy_xxxyz, g_z_0_xxy_xxxzz, g_z_0_xxy_xxyyy, g_z_0_xxy_xxyyz, g_z_0_xxy_xxyzz, g_z_0_xxy_xxzzz, g_z_0_xxy_xyyyy, g_z_0_xxy_xyyyz, g_z_0_xxy_xyyzz, g_z_0_xxy_xyzzz, g_z_0_xxy_xzzzz, g_z_0_xxy_yyyyy, g_z_0_xxy_yyyyz, g_z_0_xxy_yyyzz, g_z_0_xxy_yyzzz, g_z_0_xxy_yzzzz, g_z_0_xxy_zzzzz, g_z_0_xy_xxxxx, g_z_0_xy_xxxxxx, g_z_0_xy_xxxxxy, g_z_0_xy_xxxxxz, g_z_0_xy_xxxxy, g_z_0_xy_xxxxyy, g_z_0_xy_xxxxyz, g_z_0_xy_xxxxz, g_z_0_xy_xxxxzz, g_z_0_xy_xxxyy, g_z_0_xy_xxxyyy, g_z_0_xy_xxxyyz, g_z_0_xy_xxxyz, g_z_0_xy_xxxyzz, g_z_0_xy_xxxzz, g_z_0_xy_xxxzzz, g_z_0_xy_xxyyy, g_z_0_xy_xxyyyy, g_z_0_xy_xxyyyz, g_z_0_xy_xxyyz, g_z_0_xy_xxyyzz, g_z_0_xy_xxyzz, g_z_0_xy_xxyzzz, g_z_0_xy_xxzzz, g_z_0_xy_xxzzzz, g_z_0_xy_xyyyy, g_z_0_xy_xyyyyy, g_z_0_xy_xyyyyz, g_z_0_xy_xyyyz, g_z_0_xy_xyyyzz, g_z_0_xy_xyyzz, g_z_0_xy_xyyzzz, g_z_0_xy_xyzzz, g_z_0_xy_xyzzzz, g_z_0_xy_xzzzz, g_z_0_xy_xzzzzz, g_z_0_xy_yyyyy, g_z_0_xy_yyyyz, g_z_0_xy_yyyzz, g_z_0_xy_yyzzz, g_z_0_xy_yzzzz, g_z_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_xxxxx[k] = -g_z_0_xy_xxxxx[k] * ab_x + g_z_0_xy_xxxxxx[k];

                g_z_0_xxy_xxxxy[k] = -g_z_0_xy_xxxxy[k] * ab_x + g_z_0_xy_xxxxxy[k];

                g_z_0_xxy_xxxxz[k] = -g_z_0_xy_xxxxz[k] * ab_x + g_z_0_xy_xxxxxz[k];

                g_z_0_xxy_xxxyy[k] = -g_z_0_xy_xxxyy[k] * ab_x + g_z_0_xy_xxxxyy[k];

                g_z_0_xxy_xxxyz[k] = -g_z_0_xy_xxxyz[k] * ab_x + g_z_0_xy_xxxxyz[k];

                g_z_0_xxy_xxxzz[k] = -g_z_0_xy_xxxzz[k] * ab_x + g_z_0_xy_xxxxzz[k];

                g_z_0_xxy_xxyyy[k] = -g_z_0_xy_xxyyy[k] * ab_x + g_z_0_xy_xxxyyy[k];

                g_z_0_xxy_xxyyz[k] = -g_z_0_xy_xxyyz[k] * ab_x + g_z_0_xy_xxxyyz[k];

                g_z_0_xxy_xxyzz[k] = -g_z_0_xy_xxyzz[k] * ab_x + g_z_0_xy_xxxyzz[k];

                g_z_0_xxy_xxzzz[k] = -g_z_0_xy_xxzzz[k] * ab_x + g_z_0_xy_xxxzzz[k];

                g_z_0_xxy_xyyyy[k] = -g_z_0_xy_xyyyy[k] * ab_x + g_z_0_xy_xxyyyy[k];

                g_z_0_xxy_xyyyz[k] = -g_z_0_xy_xyyyz[k] * ab_x + g_z_0_xy_xxyyyz[k];

                g_z_0_xxy_xyyzz[k] = -g_z_0_xy_xyyzz[k] * ab_x + g_z_0_xy_xxyyzz[k];

                g_z_0_xxy_xyzzz[k] = -g_z_0_xy_xyzzz[k] * ab_x + g_z_0_xy_xxyzzz[k];

                g_z_0_xxy_xzzzz[k] = -g_z_0_xy_xzzzz[k] * ab_x + g_z_0_xy_xxzzzz[k];

                g_z_0_xxy_yyyyy[k] = -g_z_0_xy_yyyyy[k] * ab_x + g_z_0_xy_xyyyyy[k];

                g_z_0_xxy_yyyyz[k] = -g_z_0_xy_yyyyz[k] * ab_x + g_z_0_xy_xyyyyz[k];

                g_z_0_xxy_yyyzz[k] = -g_z_0_xy_yyyzz[k] * ab_x + g_z_0_xy_xyyyzz[k];

                g_z_0_xxy_yyzzz[k] = -g_z_0_xy_yyzzz[k] * ab_x + g_z_0_xy_xyyzzz[k];

                g_z_0_xxy_yzzzz[k] = -g_z_0_xy_yzzzz[k] * ab_x + g_z_0_xy_xyzzzz[k];

                g_z_0_xxy_zzzzz[k] = -g_z_0_xy_zzzzz[k] * ab_x + g_z_0_xy_xzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_xxxxx, g_z_0_xxz_xxxxy, g_z_0_xxz_xxxxz, g_z_0_xxz_xxxyy, g_z_0_xxz_xxxyz, g_z_0_xxz_xxxzz, g_z_0_xxz_xxyyy, g_z_0_xxz_xxyyz, g_z_0_xxz_xxyzz, g_z_0_xxz_xxzzz, g_z_0_xxz_xyyyy, g_z_0_xxz_xyyyz, g_z_0_xxz_xyyzz, g_z_0_xxz_xyzzz, g_z_0_xxz_xzzzz, g_z_0_xxz_yyyyy, g_z_0_xxz_yyyyz, g_z_0_xxz_yyyzz, g_z_0_xxz_yyzzz, g_z_0_xxz_yzzzz, g_z_0_xxz_zzzzz, g_z_0_xz_xxxxx, g_z_0_xz_xxxxxx, g_z_0_xz_xxxxxy, g_z_0_xz_xxxxxz, g_z_0_xz_xxxxy, g_z_0_xz_xxxxyy, g_z_0_xz_xxxxyz, g_z_0_xz_xxxxz, g_z_0_xz_xxxxzz, g_z_0_xz_xxxyy, g_z_0_xz_xxxyyy, g_z_0_xz_xxxyyz, g_z_0_xz_xxxyz, g_z_0_xz_xxxyzz, g_z_0_xz_xxxzz, g_z_0_xz_xxxzzz, g_z_0_xz_xxyyy, g_z_0_xz_xxyyyy, g_z_0_xz_xxyyyz, g_z_0_xz_xxyyz, g_z_0_xz_xxyyzz, g_z_0_xz_xxyzz, g_z_0_xz_xxyzzz, g_z_0_xz_xxzzz, g_z_0_xz_xxzzzz, g_z_0_xz_xyyyy, g_z_0_xz_xyyyyy, g_z_0_xz_xyyyyz, g_z_0_xz_xyyyz, g_z_0_xz_xyyyzz, g_z_0_xz_xyyzz, g_z_0_xz_xyyzzz, g_z_0_xz_xyzzz, g_z_0_xz_xyzzzz, g_z_0_xz_xzzzz, g_z_0_xz_xzzzzz, g_z_0_xz_yyyyy, g_z_0_xz_yyyyz, g_z_0_xz_yyyzz, g_z_0_xz_yyzzz, g_z_0_xz_yzzzz, g_z_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_xxxxx[k] = -g_z_0_xz_xxxxx[k] * ab_x + g_z_0_xz_xxxxxx[k];

                g_z_0_xxz_xxxxy[k] = -g_z_0_xz_xxxxy[k] * ab_x + g_z_0_xz_xxxxxy[k];

                g_z_0_xxz_xxxxz[k] = -g_z_0_xz_xxxxz[k] * ab_x + g_z_0_xz_xxxxxz[k];

                g_z_0_xxz_xxxyy[k] = -g_z_0_xz_xxxyy[k] * ab_x + g_z_0_xz_xxxxyy[k];

                g_z_0_xxz_xxxyz[k] = -g_z_0_xz_xxxyz[k] * ab_x + g_z_0_xz_xxxxyz[k];

                g_z_0_xxz_xxxzz[k] = -g_z_0_xz_xxxzz[k] * ab_x + g_z_0_xz_xxxxzz[k];

                g_z_0_xxz_xxyyy[k] = -g_z_0_xz_xxyyy[k] * ab_x + g_z_0_xz_xxxyyy[k];

                g_z_0_xxz_xxyyz[k] = -g_z_0_xz_xxyyz[k] * ab_x + g_z_0_xz_xxxyyz[k];

                g_z_0_xxz_xxyzz[k] = -g_z_0_xz_xxyzz[k] * ab_x + g_z_0_xz_xxxyzz[k];

                g_z_0_xxz_xxzzz[k] = -g_z_0_xz_xxzzz[k] * ab_x + g_z_0_xz_xxxzzz[k];

                g_z_0_xxz_xyyyy[k] = -g_z_0_xz_xyyyy[k] * ab_x + g_z_0_xz_xxyyyy[k];

                g_z_0_xxz_xyyyz[k] = -g_z_0_xz_xyyyz[k] * ab_x + g_z_0_xz_xxyyyz[k];

                g_z_0_xxz_xyyzz[k] = -g_z_0_xz_xyyzz[k] * ab_x + g_z_0_xz_xxyyzz[k];

                g_z_0_xxz_xyzzz[k] = -g_z_0_xz_xyzzz[k] * ab_x + g_z_0_xz_xxyzzz[k];

                g_z_0_xxz_xzzzz[k] = -g_z_0_xz_xzzzz[k] * ab_x + g_z_0_xz_xxzzzz[k];

                g_z_0_xxz_yyyyy[k] = -g_z_0_xz_yyyyy[k] * ab_x + g_z_0_xz_xyyyyy[k];

                g_z_0_xxz_yyyyz[k] = -g_z_0_xz_yyyyz[k] * ab_x + g_z_0_xz_xyyyyz[k];

                g_z_0_xxz_yyyzz[k] = -g_z_0_xz_yyyzz[k] * ab_x + g_z_0_xz_xyyyzz[k];

                g_z_0_xxz_yyzzz[k] = -g_z_0_xz_yyzzz[k] * ab_x + g_z_0_xz_xyyzzz[k];

                g_z_0_xxz_yzzzz[k] = -g_z_0_xz_yzzzz[k] * ab_x + g_z_0_xz_xyzzzz[k];

                g_z_0_xxz_zzzzz[k] = -g_z_0_xz_zzzzz[k] * ab_x + g_z_0_xz_xzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyy_xxxxx, g_z_0_xyy_xxxxy, g_z_0_xyy_xxxxz, g_z_0_xyy_xxxyy, g_z_0_xyy_xxxyz, g_z_0_xyy_xxxzz, g_z_0_xyy_xxyyy, g_z_0_xyy_xxyyz, g_z_0_xyy_xxyzz, g_z_0_xyy_xxzzz, g_z_0_xyy_xyyyy, g_z_0_xyy_xyyyz, g_z_0_xyy_xyyzz, g_z_0_xyy_xyzzz, g_z_0_xyy_xzzzz, g_z_0_xyy_yyyyy, g_z_0_xyy_yyyyz, g_z_0_xyy_yyyzz, g_z_0_xyy_yyzzz, g_z_0_xyy_yzzzz, g_z_0_xyy_zzzzz, g_z_0_yy_xxxxx, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxy, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxyy, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxyyy, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xyyyy, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_yyyyy, g_z_0_yy_yyyyz, g_z_0_yy_yyyzz, g_z_0_yy_yyzzz, g_z_0_yy_yzzzz, g_z_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_xxxxx[k] = -g_z_0_yy_xxxxx[k] * ab_x + g_z_0_yy_xxxxxx[k];

                g_z_0_xyy_xxxxy[k] = -g_z_0_yy_xxxxy[k] * ab_x + g_z_0_yy_xxxxxy[k];

                g_z_0_xyy_xxxxz[k] = -g_z_0_yy_xxxxz[k] * ab_x + g_z_0_yy_xxxxxz[k];

                g_z_0_xyy_xxxyy[k] = -g_z_0_yy_xxxyy[k] * ab_x + g_z_0_yy_xxxxyy[k];

                g_z_0_xyy_xxxyz[k] = -g_z_0_yy_xxxyz[k] * ab_x + g_z_0_yy_xxxxyz[k];

                g_z_0_xyy_xxxzz[k] = -g_z_0_yy_xxxzz[k] * ab_x + g_z_0_yy_xxxxzz[k];

                g_z_0_xyy_xxyyy[k] = -g_z_0_yy_xxyyy[k] * ab_x + g_z_0_yy_xxxyyy[k];

                g_z_0_xyy_xxyyz[k] = -g_z_0_yy_xxyyz[k] * ab_x + g_z_0_yy_xxxyyz[k];

                g_z_0_xyy_xxyzz[k] = -g_z_0_yy_xxyzz[k] * ab_x + g_z_0_yy_xxxyzz[k];

                g_z_0_xyy_xxzzz[k] = -g_z_0_yy_xxzzz[k] * ab_x + g_z_0_yy_xxxzzz[k];

                g_z_0_xyy_xyyyy[k] = -g_z_0_yy_xyyyy[k] * ab_x + g_z_0_yy_xxyyyy[k];

                g_z_0_xyy_xyyyz[k] = -g_z_0_yy_xyyyz[k] * ab_x + g_z_0_yy_xxyyyz[k];

                g_z_0_xyy_xyyzz[k] = -g_z_0_yy_xyyzz[k] * ab_x + g_z_0_yy_xxyyzz[k];

                g_z_0_xyy_xyzzz[k] = -g_z_0_yy_xyzzz[k] * ab_x + g_z_0_yy_xxyzzz[k];

                g_z_0_xyy_xzzzz[k] = -g_z_0_yy_xzzzz[k] * ab_x + g_z_0_yy_xxzzzz[k];

                g_z_0_xyy_yyyyy[k] = -g_z_0_yy_yyyyy[k] * ab_x + g_z_0_yy_xyyyyy[k];

                g_z_0_xyy_yyyyz[k] = -g_z_0_yy_yyyyz[k] * ab_x + g_z_0_yy_xyyyyz[k];

                g_z_0_xyy_yyyzz[k] = -g_z_0_yy_yyyzz[k] * ab_x + g_z_0_yy_xyyyzz[k];

                g_z_0_xyy_yyzzz[k] = -g_z_0_yy_yyzzz[k] * ab_x + g_z_0_yy_xyyzzz[k];

                g_z_0_xyy_yzzzz[k] = -g_z_0_yy_yzzzz[k] * ab_x + g_z_0_yy_xyzzzz[k];

                g_z_0_xyy_zzzzz[k] = -g_z_0_yy_zzzzz[k] * ab_x + g_z_0_yy_xzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 509 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyz_xxxxx, g_z_0_xyz_xxxxy, g_z_0_xyz_xxxxz, g_z_0_xyz_xxxyy, g_z_0_xyz_xxxyz, g_z_0_xyz_xxxzz, g_z_0_xyz_xxyyy, g_z_0_xyz_xxyyz, g_z_0_xyz_xxyzz, g_z_0_xyz_xxzzz, g_z_0_xyz_xyyyy, g_z_0_xyz_xyyyz, g_z_0_xyz_xyyzz, g_z_0_xyz_xyzzz, g_z_0_xyz_xzzzz, g_z_0_xyz_yyyyy, g_z_0_xyz_yyyyz, g_z_0_xyz_yyyzz, g_z_0_xyz_yyzzz, g_z_0_xyz_yzzzz, g_z_0_xyz_zzzzz, g_z_0_yz_xxxxx, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxy, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxyy, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxyyy, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xyyyy, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_yyyyy, g_z_0_yz_yyyyz, g_z_0_yz_yyyzz, g_z_0_yz_yyzzz, g_z_0_yz_yzzzz, g_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_xxxxx[k] = -g_z_0_yz_xxxxx[k] * ab_x + g_z_0_yz_xxxxxx[k];

                g_z_0_xyz_xxxxy[k] = -g_z_0_yz_xxxxy[k] * ab_x + g_z_0_yz_xxxxxy[k];

                g_z_0_xyz_xxxxz[k] = -g_z_0_yz_xxxxz[k] * ab_x + g_z_0_yz_xxxxxz[k];

                g_z_0_xyz_xxxyy[k] = -g_z_0_yz_xxxyy[k] * ab_x + g_z_0_yz_xxxxyy[k];

                g_z_0_xyz_xxxyz[k] = -g_z_0_yz_xxxyz[k] * ab_x + g_z_0_yz_xxxxyz[k];

                g_z_0_xyz_xxxzz[k] = -g_z_0_yz_xxxzz[k] * ab_x + g_z_0_yz_xxxxzz[k];

                g_z_0_xyz_xxyyy[k] = -g_z_0_yz_xxyyy[k] * ab_x + g_z_0_yz_xxxyyy[k];

                g_z_0_xyz_xxyyz[k] = -g_z_0_yz_xxyyz[k] * ab_x + g_z_0_yz_xxxyyz[k];

                g_z_0_xyz_xxyzz[k] = -g_z_0_yz_xxyzz[k] * ab_x + g_z_0_yz_xxxyzz[k];

                g_z_0_xyz_xxzzz[k] = -g_z_0_yz_xxzzz[k] * ab_x + g_z_0_yz_xxxzzz[k];

                g_z_0_xyz_xyyyy[k] = -g_z_0_yz_xyyyy[k] * ab_x + g_z_0_yz_xxyyyy[k];

                g_z_0_xyz_xyyyz[k] = -g_z_0_yz_xyyyz[k] * ab_x + g_z_0_yz_xxyyyz[k];

                g_z_0_xyz_xyyzz[k] = -g_z_0_yz_xyyzz[k] * ab_x + g_z_0_yz_xxyyzz[k];

                g_z_0_xyz_xyzzz[k] = -g_z_0_yz_xyzzz[k] * ab_x + g_z_0_yz_xxyzzz[k];

                g_z_0_xyz_xzzzz[k] = -g_z_0_yz_xzzzz[k] * ab_x + g_z_0_yz_xxzzzz[k];

                g_z_0_xyz_yyyyy[k] = -g_z_0_yz_yyyyy[k] * ab_x + g_z_0_yz_xyyyyy[k];

                g_z_0_xyz_yyyyz[k] = -g_z_0_yz_yyyyz[k] * ab_x + g_z_0_yz_xyyyyz[k];

                g_z_0_xyz_yyyzz[k] = -g_z_0_yz_yyyzz[k] * ab_x + g_z_0_yz_xyyyzz[k];

                g_z_0_xyz_yyzzz[k] = -g_z_0_yz_yyzzz[k] * ab_x + g_z_0_yz_xyyzzz[k];

                g_z_0_xyz_yzzzz[k] = -g_z_0_yz_yzzzz[k] * ab_x + g_z_0_yz_xyzzzz[k];

                g_z_0_xyz_zzzzz[k] = -g_z_0_yz_zzzzz[k] * ab_x + g_z_0_yz_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 539 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_xxxxx, g_z_0_xzz_xxxxy, g_z_0_xzz_xxxxz, g_z_0_xzz_xxxyy, g_z_0_xzz_xxxyz, g_z_0_xzz_xxxzz, g_z_0_xzz_xxyyy, g_z_0_xzz_xxyyz, g_z_0_xzz_xxyzz, g_z_0_xzz_xxzzz, g_z_0_xzz_xyyyy, g_z_0_xzz_xyyyz, g_z_0_xzz_xyyzz, g_z_0_xzz_xyzzz, g_z_0_xzz_xzzzz, g_z_0_xzz_yyyyy, g_z_0_xzz_yyyyz, g_z_0_xzz_yyyzz, g_z_0_xzz_yyzzz, g_z_0_xzz_yzzzz, g_z_0_xzz_zzzzz, g_z_0_zz_xxxxx, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyz, g_z_0_zz_yyyzz, g_z_0_zz_yyzzz, g_z_0_zz_yzzzz, g_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_xxxxx[k] = -g_z_0_zz_xxxxx[k] * ab_x + g_z_0_zz_xxxxxx[k];

                g_z_0_xzz_xxxxy[k] = -g_z_0_zz_xxxxy[k] * ab_x + g_z_0_zz_xxxxxy[k];

                g_z_0_xzz_xxxxz[k] = -g_z_0_zz_xxxxz[k] * ab_x + g_z_0_zz_xxxxxz[k];

                g_z_0_xzz_xxxyy[k] = -g_z_0_zz_xxxyy[k] * ab_x + g_z_0_zz_xxxxyy[k];

                g_z_0_xzz_xxxyz[k] = -g_z_0_zz_xxxyz[k] * ab_x + g_z_0_zz_xxxxyz[k];

                g_z_0_xzz_xxxzz[k] = -g_z_0_zz_xxxzz[k] * ab_x + g_z_0_zz_xxxxzz[k];

                g_z_0_xzz_xxyyy[k] = -g_z_0_zz_xxyyy[k] * ab_x + g_z_0_zz_xxxyyy[k];

                g_z_0_xzz_xxyyz[k] = -g_z_0_zz_xxyyz[k] * ab_x + g_z_0_zz_xxxyyz[k];

                g_z_0_xzz_xxyzz[k] = -g_z_0_zz_xxyzz[k] * ab_x + g_z_0_zz_xxxyzz[k];

                g_z_0_xzz_xxzzz[k] = -g_z_0_zz_xxzzz[k] * ab_x + g_z_0_zz_xxxzzz[k];

                g_z_0_xzz_xyyyy[k] = -g_z_0_zz_xyyyy[k] * ab_x + g_z_0_zz_xxyyyy[k];

                g_z_0_xzz_xyyyz[k] = -g_z_0_zz_xyyyz[k] * ab_x + g_z_0_zz_xxyyyz[k];

                g_z_0_xzz_xyyzz[k] = -g_z_0_zz_xyyzz[k] * ab_x + g_z_0_zz_xxyyzz[k];

                g_z_0_xzz_xyzzz[k] = -g_z_0_zz_xyzzz[k] * ab_x + g_z_0_zz_xxyzzz[k];

                g_z_0_xzz_xzzzz[k] = -g_z_0_zz_xzzzz[k] * ab_x + g_z_0_zz_xxzzzz[k];

                g_z_0_xzz_yyyyy[k] = -g_z_0_zz_yyyyy[k] * ab_x + g_z_0_zz_xyyyyy[k];

                g_z_0_xzz_yyyyz[k] = -g_z_0_zz_yyyyz[k] * ab_x + g_z_0_zz_xyyyyz[k];

                g_z_0_xzz_yyyzz[k] = -g_z_0_zz_yyyzz[k] * ab_x + g_z_0_zz_xyyyzz[k];

                g_z_0_xzz_yyzzz[k] = -g_z_0_zz_yyzzz[k] * ab_x + g_z_0_zz_xyyzzz[k];

                g_z_0_xzz_yzzzz[k] = -g_z_0_zz_yzzzz[k] * ab_x + g_z_0_zz_xyzzzz[k];

                g_z_0_xzz_zzzzz[k] = -g_z_0_zz_zzzzz[k] * ab_x + g_z_0_zz_xzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_xxxxx, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxy, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxz, g_z_0_yy_xxxyy, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxzz, g_z_0_yy_xxyyy, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxzzz, g_z_0_yy_xyyyy, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xzzzz, g_z_0_yy_yyyyy, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_zzzzz, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_xxxxx[k] = -g_z_0_yy_xxxxx[k] * ab_y + g_z_0_yy_xxxxxy[k];

                g_z_0_yyy_xxxxy[k] = -g_z_0_yy_xxxxy[k] * ab_y + g_z_0_yy_xxxxyy[k];

                g_z_0_yyy_xxxxz[k] = -g_z_0_yy_xxxxz[k] * ab_y + g_z_0_yy_xxxxyz[k];

                g_z_0_yyy_xxxyy[k] = -g_z_0_yy_xxxyy[k] * ab_y + g_z_0_yy_xxxyyy[k];

                g_z_0_yyy_xxxyz[k] = -g_z_0_yy_xxxyz[k] * ab_y + g_z_0_yy_xxxyyz[k];

                g_z_0_yyy_xxxzz[k] = -g_z_0_yy_xxxzz[k] * ab_y + g_z_0_yy_xxxyzz[k];

                g_z_0_yyy_xxyyy[k] = -g_z_0_yy_xxyyy[k] * ab_y + g_z_0_yy_xxyyyy[k];

                g_z_0_yyy_xxyyz[k] = -g_z_0_yy_xxyyz[k] * ab_y + g_z_0_yy_xxyyyz[k];

                g_z_0_yyy_xxyzz[k] = -g_z_0_yy_xxyzz[k] * ab_y + g_z_0_yy_xxyyzz[k];

                g_z_0_yyy_xxzzz[k] = -g_z_0_yy_xxzzz[k] * ab_y + g_z_0_yy_xxyzzz[k];

                g_z_0_yyy_xyyyy[k] = -g_z_0_yy_xyyyy[k] * ab_y + g_z_0_yy_xyyyyy[k];

                g_z_0_yyy_xyyyz[k] = -g_z_0_yy_xyyyz[k] * ab_y + g_z_0_yy_xyyyyz[k];

                g_z_0_yyy_xyyzz[k] = -g_z_0_yy_xyyzz[k] * ab_y + g_z_0_yy_xyyyzz[k];

                g_z_0_yyy_xyzzz[k] = -g_z_0_yy_xyzzz[k] * ab_y + g_z_0_yy_xyyzzz[k];

                g_z_0_yyy_xzzzz[k] = -g_z_0_yy_xzzzz[k] * ab_y + g_z_0_yy_xyzzzz[k];

                g_z_0_yyy_yyyyy[k] = -g_z_0_yy_yyyyy[k] * ab_y + g_z_0_yy_yyyyyy[k];

                g_z_0_yyy_yyyyz[k] = -g_z_0_yy_yyyyz[k] * ab_y + g_z_0_yy_yyyyyz[k];

                g_z_0_yyy_yyyzz[k] = -g_z_0_yy_yyyzz[k] * ab_y + g_z_0_yy_yyyyzz[k];

                g_z_0_yyy_yyzzz[k] = -g_z_0_yy_yyzzz[k] * ab_y + g_z_0_yy_yyyzzz[k];

                g_z_0_yyy_yzzzz[k] = -g_z_0_yy_yzzzz[k] * ab_y + g_z_0_yy_yyzzzz[k];

                g_z_0_yyy_zzzzz[k] = -g_z_0_yy_zzzzz[k] * ab_y + g_z_0_yy_yzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_zzzzz, g_z_0_yz_xxxxx, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxy, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxz, g_z_0_yz_xxxyy, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxzz, g_z_0_yz_xxyyy, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxzzz, g_z_0_yz_xyyyy, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xzzzz, g_z_0_yz_yyyyy, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_xxxxx[k] = -g_z_0_yz_xxxxx[k] * ab_y + g_z_0_yz_xxxxxy[k];

                g_z_0_yyz_xxxxy[k] = -g_z_0_yz_xxxxy[k] * ab_y + g_z_0_yz_xxxxyy[k];

                g_z_0_yyz_xxxxz[k] = -g_z_0_yz_xxxxz[k] * ab_y + g_z_0_yz_xxxxyz[k];

                g_z_0_yyz_xxxyy[k] = -g_z_0_yz_xxxyy[k] * ab_y + g_z_0_yz_xxxyyy[k];

                g_z_0_yyz_xxxyz[k] = -g_z_0_yz_xxxyz[k] * ab_y + g_z_0_yz_xxxyyz[k];

                g_z_0_yyz_xxxzz[k] = -g_z_0_yz_xxxzz[k] * ab_y + g_z_0_yz_xxxyzz[k];

                g_z_0_yyz_xxyyy[k] = -g_z_0_yz_xxyyy[k] * ab_y + g_z_0_yz_xxyyyy[k];

                g_z_0_yyz_xxyyz[k] = -g_z_0_yz_xxyyz[k] * ab_y + g_z_0_yz_xxyyyz[k];

                g_z_0_yyz_xxyzz[k] = -g_z_0_yz_xxyzz[k] * ab_y + g_z_0_yz_xxyyzz[k];

                g_z_0_yyz_xxzzz[k] = -g_z_0_yz_xxzzz[k] * ab_y + g_z_0_yz_xxyzzz[k];

                g_z_0_yyz_xyyyy[k] = -g_z_0_yz_xyyyy[k] * ab_y + g_z_0_yz_xyyyyy[k];

                g_z_0_yyz_xyyyz[k] = -g_z_0_yz_xyyyz[k] * ab_y + g_z_0_yz_xyyyyz[k];

                g_z_0_yyz_xyyzz[k] = -g_z_0_yz_xyyzz[k] * ab_y + g_z_0_yz_xyyyzz[k];

                g_z_0_yyz_xyzzz[k] = -g_z_0_yz_xyzzz[k] * ab_y + g_z_0_yz_xyyzzz[k];

                g_z_0_yyz_xzzzz[k] = -g_z_0_yz_xzzzz[k] * ab_y + g_z_0_yz_xyzzzz[k];

                g_z_0_yyz_yyyyy[k] = -g_z_0_yz_yyyyy[k] * ab_y + g_z_0_yz_yyyyyy[k];

                g_z_0_yyz_yyyyz[k] = -g_z_0_yz_yyyyz[k] * ab_y + g_z_0_yz_yyyyyz[k];

                g_z_0_yyz_yyyzz[k] = -g_z_0_yz_yyyzz[k] * ab_y + g_z_0_yz_yyyyzz[k];

                g_z_0_yyz_yyzzz[k] = -g_z_0_yz_yyzzz[k] * ab_y + g_z_0_yz_yyyzzz[k];

                g_z_0_yyz_yzzzz[k] = -g_z_0_yz_yzzzz[k] * ab_y + g_z_0_yz_yyzzzz[k];

                g_z_0_yyz_zzzzz[k] = -g_z_0_yz_zzzzz[k] * ab_y + g_z_0_yz_yzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_zzzzz, g_z_0_zz_xxxxx, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_xxxxx[k] = -g_z_0_zz_xxxxx[k] * ab_y + g_z_0_zz_xxxxxy[k];

                g_z_0_yzz_xxxxy[k] = -g_z_0_zz_xxxxy[k] * ab_y + g_z_0_zz_xxxxyy[k];

                g_z_0_yzz_xxxxz[k] = -g_z_0_zz_xxxxz[k] * ab_y + g_z_0_zz_xxxxyz[k];

                g_z_0_yzz_xxxyy[k] = -g_z_0_zz_xxxyy[k] * ab_y + g_z_0_zz_xxxyyy[k];

                g_z_0_yzz_xxxyz[k] = -g_z_0_zz_xxxyz[k] * ab_y + g_z_0_zz_xxxyyz[k];

                g_z_0_yzz_xxxzz[k] = -g_z_0_zz_xxxzz[k] * ab_y + g_z_0_zz_xxxyzz[k];

                g_z_0_yzz_xxyyy[k] = -g_z_0_zz_xxyyy[k] * ab_y + g_z_0_zz_xxyyyy[k];

                g_z_0_yzz_xxyyz[k] = -g_z_0_zz_xxyyz[k] * ab_y + g_z_0_zz_xxyyyz[k];

                g_z_0_yzz_xxyzz[k] = -g_z_0_zz_xxyzz[k] * ab_y + g_z_0_zz_xxyyzz[k];

                g_z_0_yzz_xxzzz[k] = -g_z_0_zz_xxzzz[k] * ab_y + g_z_0_zz_xxyzzz[k];

                g_z_0_yzz_xyyyy[k] = -g_z_0_zz_xyyyy[k] * ab_y + g_z_0_zz_xyyyyy[k];

                g_z_0_yzz_xyyyz[k] = -g_z_0_zz_xyyyz[k] * ab_y + g_z_0_zz_xyyyyz[k];

                g_z_0_yzz_xyyzz[k] = -g_z_0_zz_xyyzz[k] * ab_y + g_z_0_zz_xyyyzz[k];

                g_z_0_yzz_xyzzz[k] = -g_z_0_zz_xyzzz[k] * ab_y + g_z_0_zz_xyyzzz[k];

                g_z_0_yzz_xzzzz[k] = -g_z_0_zz_xzzzz[k] * ab_y + g_z_0_zz_xyzzzz[k];

                g_z_0_yzz_yyyyy[k] = -g_z_0_zz_yyyyy[k] * ab_y + g_z_0_zz_yyyyyy[k];

                g_z_0_yzz_yyyyz[k] = -g_z_0_zz_yyyyz[k] * ab_y + g_z_0_zz_yyyyyz[k];

                g_z_0_yzz_yyyzz[k] = -g_z_0_zz_yyyzz[k] * ab_y + g_z_0_zz_yyyyzz[k];

                g_z_0_yzz_yyzzz[k] = -g_z_0_zz_yyzzz[k] * ab_y + g_z_0_zz_yyyzzz[k];

                g_z_0_yzz_yzzzz[k] = -g_z_0_zz_yzzzz[k] * ab_y + g_z_0_zz_yyzzzz[k];

                g_z_0_yzz_zzzzz[k] = -g_z_0_zz_zzzzz[k] * ab_y + g_z_0_zz_yzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxxxx, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzz, g_z_0_zz_zzzzzz, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzzz, g_zz_xxxxx, g_zz_xxxxy, g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxzz, g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyzz, g_zz_xxzzz, g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyzz, g_zz_xyzzz, g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz, g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_xxxxx[k] = -g_zz_xxxxx[k] - g_z_0_zz_xxxxx[k] * ab_z + g_z_0_zz_xxxxxz[k];

                g_z_0_zzz_xxxxy[k] = -g_zz_xxxxy[k] - g_z_0_zz_xxxxy[k] * ab_z + g_z_0_zz_xxxxyz[k];

                g_z_0_zzz_xxxxz[k] = -g_zz_xxxxz[k] - g_z_0_zz_xxxxz[k] * ab_z + g_z_0_zz_xxxxzz[k];

                g_z_0_zzz_xxxyy[k] = -g_zz_xxxyy[k] - g_z_0_zz_xxxyy[k] * ab_z + g_z_0_zz_xxxyyz[k];

                g_z_0_zzz_xxxyz[k] = -g_zz_xxxyz[k] - g_z_0_zz_xxxyz[k] * ab_z + g_z_0_zz_xxxyzz[k];

                g_z_0_zzz_xxxzz[k] = -g_zz_xxxzz[k] - g_z_0_zz_xxxzz[k] * ab_z + g_z_0_zz_xxxzzz[k];

                g_z_0_zzz_xxyyy[k] = -g_zz_xxyyy[k] - g_z_0_zz_xxyyy[k] * ab_z + g_z_0_zz_xxyyyz[k];

                g_z_0_zzz_xxyyz[k] = -g_zz_xxyyz[k] - g_z_0_zz_xxyyz[k] * ab_z + g_z_0_zz_xxyyzz[k];

                g_z_0_zzz_xxyzz[k] = -g_zz_xxyzz[k] - g_z_0_zz_xxyzz[k] * ab_z + g_z_0_zz_xxyzzz[k];

                g_z_0_zzz_xxzzz[k] = -g_zz_xxzzz[k] - g_z_0_zz_xxzzz[k] * ab_z + g_z_0_zz_xxzzzz[k];

                g_z_0_zzz_xyyyy[k] = -g_zz_xyyyy[k] - g_z_0_zz_xyyyy[k] * ab_z + g_z_0_zz_xyyyyz[k];

                g_z_0_zzz_xyyyz[k] = -g_zz_xyyyz[k] - g_z_0_zz_xyyyz[k] * ab_z + g_z_0_zz_xyyyzz[k];

                g_z_0_zzz_xyyzz[k] = -g_zz_xyyzz[k] - g_z_0_zz_xyyzz[k] * ab_z + g_z_0_zz_xyyzzz[k];

                g_z_0_zzz_xyzzz[k] = -g_zz_xyzzz[k] - g_z_0_zz_xyzzz[k] * ab_z + g_z_0_zz_xyzzzz[k];

                g_z_0_zzz_xzzzz[k] = -g_zz_xzzzz[k] - g_z_0_zz_xzzzz[k] * ab_z + g_z_0_zz_xzzzzz[k];

                g_z_0_zzz_yyyyy[k] = -g_zz_yyyyy[k] - g_z_0_zz_yyyyy[k] * ab_z + g_z_0_zz_yyyyyz[k];

                g_z_0_zzz_yyyyz[k] = -g_zz_yyyyz[k] - g_z_0_zz_yyyyz[k] * ab_z + g_z_0_zz_yyyyzz[k];

                g_z_0_zzz_yyyzz[k] = -g_zz_yyyzz[k] - g_z_0_zz_yyyzz[k] * ab_z + g_z_0_zz_yyyzzz[k];

                g_z_0_zzz_yyzzz[k] = -g_zz_yyzzz[k] - g_z_0_zz_yyzzz[k] * ab_z + g_z_0_zz_yyzzzz[k];

                g_z_0_zzz_yzzzz[k] = -g_zz_yzzzz[k] - g_z_0_zz_yzzzz[k] * ab_z + g_z_0_zz_yzzzzz[k];

                g_z_0_zzz_zzzzz[k] = -g_zz_zzzzz[k] - g_z_0_zz_zzzzz[k] * ab_z + g_z_0_zz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace


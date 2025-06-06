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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXDH.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xdh(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdh,
                                        const size_t idx_xph,
                                        const size_t idx_geom_10_xph,
                                        const size_t idx_geom_10_xpi,
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
        /// Set up components of auxilary buffer : SPH

        const auto ph_off = idx_xph + i * 63;

        auto g_x_xxxxx = cbuffer.data(ph_off + 0);

        auto g_x_xxxxy = cbuffer.data(ph_off + 1);

        auto g_x_xxxxz = cbuffer.data(ph_off + 2);

        auto g_x_xxxyy = cbuffer.data(ph_off + 3);

        auto g_x_xxxyz = cbuffer.data(ph_off + 4);

        auto g_x_xxxzz = cbuffer.data(ph_off + 5);

        auto g_x_xxyyy = cbuffer.data(ph_off + 6);

        auto g_x_xxyyz = cbuffer.data(ph_off + 7);

        auto g_x_xxyzz = cbuffer.data(ph_off + 8);

        auto g_x_xxzzz = cbuffer.data(ph_off + 9);

        auto g_x_xyyyy = cbuffer.data(ph_off + 10);

        auto g_x_xyyyz = cbuffer.data(ph_off + 11);

        auto g_x_xyyzz = cbuffer.data(ph_off + 12);

        auto g_x_xyzzz = cbuffer.data(ph_off + 13);

        auto g_x_xzzzz = cbuffer.data(ph_off + 14);

        auto g_x_yyyyy = cbuffer.data(ph_off + 15);

        auto g_x_yyyyz = cbuffer.data(ph_off + 16);

        auto g_x_yyyzz = cbuffer.data(ph_off + 17);

        auto g_x_yyzzz = cbuffer.data(ph_off + 18);

        auto g_x_yzzzz = cbuffer.data(ph_off + 19);

        auto g_x_zzzzz = cbuffer.data(ph_off + 20);

        auto g_y_xxxxx = cbuffer.data(ph_off + 21);

        auto g_y_xxxxy = cbuffer.data(ph_off + 22);

        auto g_y_xxxxz = cbuffer.data(ph_off + 23);

        auto g_y_xxxyy = cbuffer.data(ph_off + 24);

        auto g_y_xxxyz = cbuffer.data(ph_off + 25);

        auto g_y_xxxzz = cbuffer.data(ph_off + 26);

        auto g_y_xxyyy = cbuffer.data(ph_off + 27);

        auto g_y_xxyyz = cbuffer.data(ph_off + 28);

        auto g_y_xxyzz = cbuffer.data(ph_off + 29);

        auto g_y_xxzzz = cbuffer.data(ph_off + 30);

        auto g_y_xyyyy = cbuffer.data(ph_off + 31);

        auto g_y_xyyyz = cbuffer.data(ph_off + 32);

        auto g_y_xyyzz = cbuffer.data(ph_off + 33);

        auto g_y_xyzzz = cbuffer.data(ph_off + 34);

        auto g_y_xzzzz = cbuffer.data(ph_off + 35);

        auto g_y_yyyyy = cbuffer.data(ph_off + 36);

        auto g_y_yyyyz = cbuffer.data(ph_off + 37);

        auto g_y_yyyzz = cbuffer.data(ph_off + 38);

        auto g_y_yyzzz = cbuffer.data(ph_off + 39);

        auto g_y_yzzzz = cbuffer.data(ph_off + 40);

        auto g_y_zzzzz = cbuffer.data(ph_off + 41);

        auto g_z_xxxxx = cbuffer.data(ph_off + 42);

        auto g_z_xxxxy = cbuffer.data(ph_off + 43);

        auto g_z_xxxxz = cbuffer.data(ph_off + 44);

        auto g_z_xxxyy = cbuffer.data(ph_off + 45);

        auto g_z_xxxyz = cbuffer.data(ph_off + 46);

        auto g_z_xxxzz = cbuffer.data(ph_off + 47);

        auto g_z_xxyyy = cbuffer.data(ph_off + 48);

        auto g_z_xxyyz = cbuffer.data(ph_off + 49);

        auto g_z_xxyzz = cbuffer.data(ph_off + 50);

        auto g_z_xxzzz = cbuffer.data(ph_off + 51);

        auto g_z_xyyyy = cbuffer.data(ph_off + 52);

        auto g_z_xyyyz = cbuffer.data(ph_off + 53);

        auto g_z_xyyzz = cbuffer.data(ph_off + 54);

        auto g_z_xyzzz = cbuffer.data(ph_off + 55);

        auto g_z_xzzzz = cbuffer.data(ph_off + 56);

        auto g_z_yyyyy = cbuffer.data(ph_off + 57);

        auto g_z_yyyyz = cbuffer.data(ph_off + 58);

        auto g_z_yyyzz = cbuffer.data(ph_off + 59);

        auto g_z_yyzzz = cbuffer.data(ph_off + 60);

        auto g_z_yzzzz = cbuffer.data(ph_off + 61);

        auto g_z_zzzzz = cbuffer.data(ph_off + 62);

        /// Set up components of auxilary buffer : SPH

        const auto ph_geom_10_off = idx_geom_10_xph + i * 63;

        auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 9);

        auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 10);

        auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 11);

        auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 12);

        auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 13);

        auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 14);

        auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 15);

        auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 16);

        auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 17);

        auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 18);

        auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 19);

        auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 20);

        auto g_x_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps + 21);

        auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 22);

        auto g_x_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps + 23);

        auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 24);

        auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 25);

        auto g_x_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 26);

        auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 27);

        auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 28);

        auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 29);

        auto g_x_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 30);

        auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 31);

        auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 32);

        auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 33);

        auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 34);

        auto g_x_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 35);

        auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 36);

        auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 37);

        auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 38);

        auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 39);

        auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 40);

        auto g_x_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 41);

        auto g_x_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps + 42);

        auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 43);

        auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps + 44);

        auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 45);

        auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 46);

        auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 47);

        auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 48);

        auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 49);

        auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 50);

        auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 51);

        auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 52);

        auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 53);

        auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 54);

        auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 55);

        auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 56);

        auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 57);

        auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 58);

        auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 59);

        auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 60);

        auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 61);

        auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 62);

        auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 0);

        auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 1);

        auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 2);

        auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 3);

        auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 4);

        auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 5);

        auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 6);

        auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 7);

        auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 8);

        auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 9);

        auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 10);

        auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 11);

        auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 12);

        auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 13);

        auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 14);

        auto g_y_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 15);

        auto g_y_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 16);

        auto g_y_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 17);

        auto g_y_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 18);

        auto g_y_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 19);

        auto g_y_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 20);

        auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 21);

        auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 22);

        auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 23);

        auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 24);

        auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 25);

        auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 26);

        auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 27);

        auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 28);

        auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 29);

        auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 30);

        auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 31);

        auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 32);

        auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 33);

        auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 34);

        auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 35);

        auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 36);

        auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 37);

        auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 38);

        auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 39);

        auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 40);

        auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 41);

        auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 42);

        auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 43);

        auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 44);

        auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 45);

        auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 46);

        auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 47);

        auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 48);

        auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 49);

        auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 50);

        auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 51);

        auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 52);

        auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 53);

        auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 54);

        auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 55);

        auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 56);

        auto g_y_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 57);

        auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 58);

        auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 59);

        auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 60);

        auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 61);

        auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 62);

        auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 0);

        auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 1);

        auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 2);

        auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 3);

        auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 4);

        auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 5);

        auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 6);

        auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 7);

        auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 8);

        auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 9);

        auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 10);

        auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 11);

        auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 12);

        auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 13);

        auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 14);

        auto g_z_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 15);

        auto g_z_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 16);

        auto g_z_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 17);

        auto g_z_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 18);

        auto g_z_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 19);

        auto g_z_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 20);

        auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 21);

        auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 22);

        auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 23);

        auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 24);

        auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 25);

        auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 26);

        auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 27);

        auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 28);

        auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 29);

        auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 30);

        auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 31);

        auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 32);

        auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 33);

        auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 34);

        auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 35);

        auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 36);

        auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 37);

        auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 38);

        auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 39);

        auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 40);

        auto g_z_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 41);

        auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 42);

        auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 43);

        auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 44);

        auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 45);

        auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 46);

        auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 47);

        auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 48);

        auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 49);

        auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 50);

        auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 51);

        auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 52);

        auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 53);

        auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 54);

        auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 55);

        auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 56);

        auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 57);

        auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 58);

        auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 59);

        auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 60);

        auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 61);

        auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 62);

        /// Set up components of auxilary buffer : SPI

        const auto pi_geom_10_off = idx_geom_10_xpi + i * 84;

        auto g_x_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 9);

        auto g_x_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 10);

        auto g_x_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 11);

        auto g_x_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 12);

        auto g_x_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 13);

        auto g_x_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 14);

        auto g_x_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 15);

        auto g_x_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 16);

        auto g_x_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 17);

        auto g_x_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 18);

        auto g_x_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 19);

        auto g_x_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 20);

        auto g_x_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 21);

        auto g_x_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 22);

        auto g_x_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 23);

        auto g_x_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 24);

        auto g_x_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 25);

        auto g_x_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 26);

        auto g_x_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 27);

        auto g_x_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps + 29);

        auto g_x_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 31);

        auto g_x_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 32);

        auto g_x_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 34);

        auto g_x_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 35);

        auto g_x_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 36);

        auto g_x_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 38);

        auto g_x_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 39);

        auto g_x_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 40);

        auto g_x_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 41);

        auto g_x_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 43);

        auto g_x_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 44);

        auto g_x_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 45);

        auto g_x_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 46);

        auto g_x_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 47);

        auto g_x_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 49);

        auto g_x_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 50);

        auto g_x_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 51);

        auto g_x_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 52);

        auto g_x_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 53);

        auto g_x_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 54);

        auto g_x_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps + 57);

        auto g_x_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 0 * acomps + 58);

        auto g_x_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 59);

        auto g_x_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 60);

        auto g_x_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 61);

        auto g_x_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 62);

        auto g_x_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 63);

        auto g_x_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 64);

        auto g_x_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 65);

        auto g_x_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 66);

        auto g_x_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 67);

        auto g_x_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 68);

        auto g_x_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 69);

        auto g_x_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 70);

        auto g_x_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 71);

        auto g_x_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 72);

        auto g_x_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 73);

        auto g_x_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 74);

        auto g_x_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 75);

        auto g_x_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 76);

        auto g_x_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps + 77);

        auto g_x_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps + 78);

        auto g_x_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 79);

        auto g_x_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 80);

        auto g_x_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 81);

        auto g_x_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 82);

        auto g_x_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps + 83);

        auto g_y_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps + 0);

        auto g_y_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps + 1);

        auto g_y_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps + 2);

        auto g_y_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 3);

        auto g_y_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 4);

        auto g_y_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 5);

        auto g_y_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 6);

        auto g_y_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 7);

        auto g_y_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 8);

        auto g_y_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 9);

        auto g_y_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 10);

        auto g_y_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 11);

        auto g_y_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 12);

        auto g_y_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 13);

        auto g_y_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 14);

        auto g_y_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 15);

        auto g_y_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 16);

        auto g_y_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 17);

        auto g_y_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 18);

        auto g_y_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 19);

        auto g_y_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 20);

        auto g_y_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps + 28);

        auto g_y_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps + 29);

        auto g_y_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps + 30);

        auto g_y_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 31);

        auto g_y_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 32);

        auto g_y_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 33);

        auto g_y_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 34);

        auto g_y_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 35);

        auto g_y_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 36);

        auto g_y_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 37);

        auto g_y_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 38);

        auto g_y_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 39);

        auto g_y_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 40);

        auto g_y_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 41);

        auto g_y_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 42);

        auto g_y_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 43);

        auto g_y_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 44);

        auto g_y_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 45);

        auto g_y_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 46);

        auto g_y_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 47);

        auto g_y_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 48);

        auto g_y_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 49);

        auto g_y_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 50);

        auto g_y_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 51);

        auto g_y_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 52);

        auto g_y_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 53);

        auto g_y_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 54);

        auto g_y_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 55);

        auto g_y_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps + 56);

        auto g_y_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps + 57);

        auto g_y_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps + 58);

        auto g_y_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 59);

        auto g_y_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 60);

        auto g_y_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 61);

        auto g_y_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 62);

        auto g_y_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 63);

        auto g_y_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 64);

        auto g_y_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 65);

        auto g_y_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 66);

        auto g_y_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 67);

        auto g_y_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 68);

        auto g_y_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 69);

        auto g_y_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 70);

        auto g_y_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps + 71);

        auto g_y_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 72);

        auto g_y_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 73);

        auto g_y_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 74);

        auto g_y_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 75);

        auto g_y_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 76);

        auto g_y_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps + 78);

        auto g_y_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 79);

        auto g_y_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 80);

        auto g_y_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 81);

        auto g_y_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 82);

        auto g_y_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps + 83);

        auto g_z_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps + 0);

        auto g_z_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps + 1);

        auto g_z_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps + 2);

        auto g_z_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 3);

        auto g_z_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 4);

        auto g_z_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 5);

        auto g_z_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 6);

        auto g_z_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 7);

        auto g_z_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 8);

        auto g_z_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 9);

        auto g_z_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 10);

        auto g_z_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 11);

        auto g_z_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 12);

        auto g_z_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 13);

        auto g_z_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 14);

        auto g_z_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 15);

        auto g_z_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 16);

        auto g_z_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 17);

        auto g_z_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 18);

        auto g_z_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 19);

        auto g_z_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 20);

        auto g_z_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps + 28);

        auto g_z_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps + 29);

        auto g_z_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps + 30);

        auto g_z_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 31);

        auto g_z_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 32);

        auto g_z_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 33);

        auto g_z_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 34);

        auto g_z_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 35);

        auto g_z_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 36);

        auto g_z_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 37);

        auto g_z_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 38);

        auto g_z_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 39);

        auto g_z_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 40);

        auto g_z_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 41);

        auto g_z_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 42);

        auto g_z_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 43);

        auto g_z_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 44);

        auto g_z_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 45);

        auto g_z_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 46);

        auto g_z_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 47);

        auto g_z_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 48);

        auto g_z_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 49);

        auto g_z_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 50);

        auto g_z_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 51);

        auto g_z_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 52);

        auto g_z_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 53);

        auto g_z_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 54);

        auto g_z_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps + 56);

        auto g_z_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps + 57);

        auto g_z_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps + 58);

        auto g_z_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 59);

        auto g_z_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 60);

        auto g_z_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 61);

        auto g_z_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 62);

        auto g_z_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 63);

        auto g_z_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 64);

        auto g_z_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 65);

        auto g_z_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 66);

        auto g_z_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 67);

        auto g_z_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 68);

        auto g_z_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 69);

        auto g_z_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 70);

        auto g_z_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 71);

        auto g_z_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 72);

        auto g_z_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 73);

        auto g_z_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 74);

        auto g_z_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 75);

        auto g_z_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 76);

        auto g_z_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps + 77);

        auto g_z_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps + 78);

        auto g_z_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 79);

        auto g_z_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 80);

        auto g_z_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 81);

        auto g_z_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 82);

        auto g_z_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps + 83);

        /// set up bra offset for contr_buffer_xxdh

        const auto dh_geom_10_off = idx_geom_10_xdh + i * 126;

        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_x_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_x_0_x_xxxxx, g_x_0_x_xxxxxx, g_x_0_x_xxxxxy, g_x_0_x_xxxxxz, g_x_0_x_xxxxy, g_x_0_x_xxxxyy, g_x_0_x_xxxxyz, g_x_0_x_xxxxz, g_x_0_x_xxxxzz, g_x_0_x_xxxyy, g_x_0_x_xxxyyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyy, g_x_0_x_xxyyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyzz, g_x_0_x_yyzzz, g_x_0_x_yzzzz, g_x_0_x_zzzzz, g_x_0_xx_xxxxx, g_x_0_xx_xxxxy, g_x_0_xx_xxxxz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyz, g_x_0_xx_xxxzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyz, g_x_0_xx_xxyzz, g_x_0_xx_xxzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyz, g_x_0_xx_xyyzz, g_x_0_xx_xyzzz, g_x_0_xx_xzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyz, g_x_0_xx_yyyzz, g_x_0_xx_yyzzz, g_x_0_xx_yzzzz, g_x_0_xx_zzzzz, g_x_xxxxx, g_x_xxxxy, g_x_xxxxz, g_x_xxxyy, g_x_xxxyz, g_x_xxxzz, g_x_xxyyy, g_x_xxyyz, g_x_xxyzz, g_x_xxzzz, g_x_xyyyy, g_x_xyyyz, g_x_xyyzz, g_x_xyzzz, g_x_xzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xx_xxxxx[k] = -g_x_xxxxx[k] - g_x_0_x_xxxxx[k] * cd_x[k] + g_x_0_x_xxxxxx[k];

            g_x_0_xx_xxxxy[k] = -g_x_xxxxy[k] - g_x_0_x_xxxxy[k] * cd_x[k] + g_x_0_x_xxxxxy[k];

            g_x_0_xx_xxxxz[k] = -g_x_xxxxz[k] - g_x_0_x_xxxxz[k] * cd_x[k] + g_x_0_x_xxxxxz[k];

            g_x_0_xx_xxxyy[k] = -g_x_xxxyy[k] - g_x_0_x_xxxyy[k] * cd_x[k] + g_x_0_x_xxxxyy[k];

            g_x_0_xx_xxxyz[k] = -g_x_xxxyz[k] - g_x_0_x_xxxyz[k] * cd_x[k] + g_x_0_x_xxxxyz[k];

            g_x_0_xx_xxxzz[k] = -g_x_xxxzz[k] - g_x_0_x_xxxzz[k] * cd_x[k] + g_x_0_x_xxxxzz[k];

            g_x_0_xx_xxyyy[k] = -g_x_xxyyy[k] - g_x_0_x_xxyyy[k] * cd_x[k] + g_x_0_x_xxxyyy[k];

            g_x_0_xx_xxyyz[k] = -g_x_xxyyz[k] - g_x_0_x_xxyyz[k] * cd_x[k] + g_x_0_x_xxxyyz[k];

            g_x_0_xx_xxyzz[k] = -g_x_xxyzz[k] - g_x_0_x_xxyzz[k] * cd_x[k] + g_x_0_x_xxxyzz[k];

            g_x_0_xx_xxzzz[k] = -g_x_xxzzz[k] - g_x_0_x_xxzzz[k] * cd_x[k] + g_x_0_x_xxxzzz[k];

            g_x_0_xx_xyyyy[k] = -g_x_xyyyy[k] - g_x_0_x_xyyyy[k] * cd_x[k] + g_x_0_x_xxyyyy[k];

            g_x_0_xx_xyyyz[k] = -g_x_xyyyz[k] - g_x_0_x_xyyyz[k] * cd_x[k] + g_x_0_x_xxyyyz[k];

            g_x_0_xx_xyyzz[k] = -g_x_xyyzz[k] - g_x_0_x_xyyzz[k] * cd_x[k] + g_x_0_x_xxyyzz[k];

            g_x_0_xx_xyzzz[k] = -g_x_xyzzz[k] - g_x_0_x_xyzzz[k] * cd_x[k] + g_x_0_x_xxyzzz[k];

            g_x_0_xx_xzzzz[k] = -g_x_xzzzz[k] - g_x_0_x_xzzzz[k] * cd_x[k] + g_x_0_x_xxzzzz[k];

            g_x_0_xx_yyyyy[k] = -g_x_yyyyy[k] - g_x_0_x_yyyyy[k] * cd_x[k] + g_x_0_x_xyyyyy[k];

            g_x_0_xx_yyyyz[k] = -g_x_yyyyz[k] - g_x_0_x_yyyyz[k] * cd_x[k] + g_x_0_x_xyyyyz[k];

            g_x_0_xx_yyyzz[k] = -g_x_yyyzz[k] - g_x_0_x_yyyzz[k] * cd_x[k] + g_x_0_x_xyyyzz[k];

            g_x_0_xx_yyzzz[k] = -g_x_yyzzz[k] - g_x_0_x_yyzzz[k] * cd_x[k] + g_x_0_x_xyyzzz[k];

            g_x_0_xx_yzzzz[k] = -g_x_yzzzz[k] - g_x_0_x_yzzzz[k] * cd_x[k] + g_x_0_x_xyzzzz[k];

            g_x_0_xx_zzzzz[k] = -g_x_zzzzz[k] - g_x_0_x_zzzzz[k] * cd_x[k] + g_x_0_x_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_x_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_x_0_x_xxxxx, g_x_0_x_xxxxxy, g_x_0_x_xxxxy, g_x_0_x_xxxxyy, g_x_0_x_xxxxyz, g_x_0_x_xxxxz, g_x_0_x_xxxyy, g_x_0_x_xxxyyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzz, g_x_0_x_xxyyy, g_x_0_x_xxyyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzz, g_x_0_xy_xxxxx, g_x_0_xy_xxxxy, g_x_0_xy_xxxxz, g_x_0_xy_xxxyy, g_x_0_xy_xxxyz, g_x_0_xy_xxxzz, g_x_0_xy_xxyyy, g_x_0_xy_xxyyz, g_x_0_xy_xxyzz, g_x_0_xy_xxzzz, g_x_0_xy_xyyyy, g_x_0_xy_xyyyz, g_x_0_xy_xyyzz, g_x_0_xy_xyzzz, g_x_0_xy_xzzzz, g_x_0_xy_yyyyy, g_x_0_xy_yyyyz, g_x_0_xy_yyyzz, g_x_0_xy_yyzzz, g_x_0_xy_yzzzz, g_x_0_xy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xy_xxxxx[k] = -g_x_0_x_xxxxx[k] * cd_y[k] + g_x_0_x_xxxxxy[k];

            g_x_0_xy_xxxxy[k] = -g_x_0_x_xxxxy[k] * cd_y[k] + g_x_0_x_xxxxyy[k];

            g_x_0_xy_xxxxz[k] = -g_x_0_x_xxxxz[k] * cd_y[k] + g_x_0_x_xxxxyz[k];

            g_x_0_xy_xxxyy[k] = -g_x_0_x_xxxyy[k] * cd_y[k] + g_x_0_x_xxxyyy[k];

            g_x_0_xy_xxxyz[k] = -g_x_0_x_xxxyz[k] * cd_y[k] + g_x_0_x_xxxyyz[k];

            g_x_0_xy_xxxzz[k] = -g_x_0_x_xxxzz[k] * cd_y[k] + g_x_0_x_xxxyzz[k];

            g_x_0_xy_xxyyy[k] = -g_x_0_x_xxyyy[k] * cd_y[k] + g_x_0_x_xxyyyy[k];

            g_x_0_xy_xxyyz[k] = -g_x_0_x_xxyyz[k] * cd_y[k] + g_x_0_x_xxyyyz[k];

            g_x_0_xy_xxyzz[k] = -g_x_0_x_xxyzz[k] * cd_y[k] + g_x_0_x_xxyyzz[k];

            g_x_0_xy_xxzzz[k] = -g_x_0_x_xxzzz[k] * cd_y[k] + g_x_0_x_xxyzzz[k];

            g_x_0_xy_xyyyy[k] = -g_x_0_x_xyyyy[k] * cd_y[k] + g_x_0_x_xyyyyy[k];

            g_x_0_xy_xyyyz[k] = -g_x_0_x_xyyyz[k] * cd_y[k] + g_x_0_x_xyyyyz[k];

            g_x_0_xy_xyyzz[k] = -g_x_0_x_xyyzz[k] * cd_y[k] + g_x_0_x_xyyyzz[k];

            g_x_0_xy_xyzzz[k] = -g_x_0_x_xyzzz[k] * cd_y[k] + g_x_0_x_xyyzzz[k];

            g_x_0_xy_xzzzz[k] = -g_x_0_x_xzzzz[k] * cd_y[k] + g_x_0_x_xyzzzz[k];

            g_x_0_xy_yyyyy[k] = -g_x_0_x_yyyyy[k] * cd_y[k] + g_x_0_x_yyyyyy[k];

            g_x_0_xy_yyyyz[k] = -g_x_0_x_yyyyz[k] * cd_y[k] + g_x_0_x_yyyyyz[k];

            g_x_0_xy_yyyzz[k] = -g_x_0_x_yyyzz[k] * cd_y[k] + g_x_0_x_yyyyzz[k];

            g_x_0_xy_yyzzz[k] = -g_x_0_x_yyzzz[k] * cd_y[k] + g_x_0_x_yyyzzz[k];

            g_x_0_xy_yzzzz[k] = -g_x_0_x_yzzzz[k] * cd_y[k] + g_x_0_x_yyzzzz[k];

            g_x_0_xy_zzzzz[k] = -g_x_0_x_zzzzz[k] * cd_y[k] + g_x_0_x_yzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_x_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 59);

        auto g_x_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 62);

        #pragma omp simd aligned(cd_z, g_x_0_x_xxxxx, g_x_0_x_xxxxxz, g_x_0_x_xxxxy, g_x_0_x_xxxxyz, g_x_0_x_xxxxz, g_x_0_x_xxxxzz, g_x_0_x_xxxyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzz, g_x_0_x_zzzzzz, g_x_0_xz_xxxxx, g_x_0_xz_xxxxy, g_x_0_xz_xxxxz, g_x_0_xz_xxxyy, g_x_0_xz_xxxyz, g_x_0_xz_xxxzz, g_x_0_xz_xxyyy, g_x_0_xz_xxyyz, g_x_0_xz_xxyzz, g_x_0_xz_xxzzz, g_x_0_xz_xyyyy, g_x_0_xz_xyyyz, g_x_0_xz_xyyzz, g_x_0_xz_xyzzz, g_x_0_xz_xzzzz, g_x_0_xz_yyyyy, g_x_0_xz_yyyyz, g_x_0_xz_yyyzz, g_x_0_xz_yyzzz, g_x_0_xz_yzzzz, g_x_0_xz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xz_xxxxx[k] = -g_x_0_x_xxxxx[k] * cd_z[k] + g_x_0_x_xxxxxz[k];

            g_x_0_xz_xxxxy[k] = -g_x_0_x_xxxxy[k] * cd_z[k] + g_x_0_x_xxxxyz[k];

            g_x_0_xz_xxxxz[k] = -g_x_0_x_xxxxz[k] * cd_z[k] + g_x_0_x_xxxxzz[k];

            g_x_0_xz_xxxyy[k] = -g_x_0_x_xxxyy[k] * cd_z[k] + g_x_0_x_xxxyyz[k];

            g_x_0_xz_xxxyz[k] = -g_x_0_x_xxxyz[k] * cd_z[k] + g_x_0_x_xxxyzz[k];

            g_x_0_xz_xxxzz[k] = -g_x_0_x_xxxzz[k] * cd_z[k] + g_x_0_x_xxxzzz[k];

            g_x_0_xz_xxyyy[k] = -g_x_0_x_xxyyy[k] * cd_z[k] + g_x_0_x_xxyyyz[k];

            g_x_0_xz_xxyyz[k] = -g_x_0_x_xxyyz[k] * cd_z[k] + g_x_0_x_xxyyzz[k];

            g_x_0_xz_xxyzz[k] = -g_x_0_x_xxyzz[k] * cd_z[k] + g_x_0_x_xxyzzz[k];

            g_x_0_xz_xxzzz[k] = -g_x_0_x_xxzzz[k] * cd_z[k] + g_x_0_x_xxzzzz[k];

            g_x_0_xz_xyyyy[k] = -g_x_0_x_xyyyy[k] * cd_z[k] + g_x_0_x_xyyyyz[k];

            g_x_0_xz_xyyyz[k] = -g_x_0_x_xyyyz[k] * cd_z[k] + g_x_0_x_xyyyzz[k];

            g_x_0_xz_xyyzz[k] = -g_x_0_x_xyyzz[k] * cd_z[k] + g_x_0_x_xyyzzz[k];

            g_x_0_xz_xyzzz[k] = -g_x_0_x_xyzzz[k] * cd_z[k] + g_x_0_x_xyzzzz[k];

            g_x_0_xz_xzzzz[k] = -g_x_0_x_xzzzz[k] * cd_z[k] + g_x_0_x_xzzzzz[k];

            g_x_0_xz_yyyyy[k] = -g_x_0_x_yyyyy[k] * cd_z[k] + g_x_0_x_yyyyyz[k];

            g_x_0_xz_yyyyz[k] = -g_x_0_x_yyyyz[k] * cd_z[k] + g_x_0_x_yyyyzz[k];

            g_x_0_xz_yyyzz[k] = -g_x_0_x_yyyzz[k] * cd_z[k] + g_x_0_x_yyyzzz[k];

            g_x_0_xz_yyzzz[k] = -g_x_0_x_yyzzz[k] * cd_z[k] + g_x_0_x_yyzzzz[k];

            g_x_0_xz_yzzzz[k] = -g_x_0_x_yzzzz[k] * cd_z[k] + g_x_0_x_yzzzzz[k];

            g_x_0_xz_zzzzz[k] = -g_x_0_x_zzzzz[k] * cd_z[k] + g_x_0_x_zzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_x_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 74);

        auto g_x_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_x_0_y_xxxxx, g_x_0_y_xxxxxy, g_x_0_y_xxxxy, g_x_0_y_xxxxyy, g_x_0_y_xxxxyz, g_x_0_y_xxxxz, g_x_0_y_xxxyy, g_x_0_y_xxxyyy, g_x_0_y_xxxyyz, g_x_0_y_xxxyz, g_x_0_y_xxxyzz, g_x_0_y_xxxzz, g_x_0_y_xxyyy, g_x_0_y_xxyyyy, g_x_0_y_xxyyyz, g_x_0_y_xxyyz, g_x_0_y_xxyyzz, g_x_0_y_xxyzz, g_x_0_y_xxyzzz, g_x_0_y_xxzzz, g_x_0_y_xyyyy, g_x_0_y_xyyyyy, g_x_0_y_xyyyyz, g_x_0_y_xyyyz, g_x_0_y_xyyyzz, g_x_0_y_xyyzz, g_x_0_y_xyyzzz, g_x_0_y_xyzzz, g_x_0_y_xyzzzz, g_x_0_y_xzzzz, g_x_0_y_yyyyy, g_x_0_y_yyyyyy, g_x_0_y_yyyyyz, g_x_0_y_yyyyz, g_x_0_y_yyyyzz, g_x_0_y_yyyzz, g_x_0_y_yyyzzz, g_x_0_y_yyzzz, g_x_0_y_yyzzzz, g_x_0_y_yzzzz, g_x_0_y_yzzzzz, g_x_0_y_zzzzz, g_x_0_yy_xxxxx, g_x_0_yy_xxxxy, g_x_0_yy_xxxxz, g_x_0_yy_xxxyy, g_x_0_yy_xxxyz, g_x_0_yy_xxxzz, g_x_0_yy_xxyyy, g_x_0_yy_xxyyz, g_x_0_yy_xxyzz, g_x_0_yy_xxzzz, g_x_0_yy_xyyyy, g_x_0_yy_xyyyz, g_x_0_yy_xyyzz, g_x_0_yy_xyzzz, g_x_0_yy_xzzzz, g_x_0_yy_yyyyy, g_x_0_yy_yyyyz, g_x_0_yy_yyyzz, g_x_0_yy_yyzzz, g_x_0_yy_yzzzz, g_x_0_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yy_xxxxx[k] = -g_x_0_y_xxxxx[k] * cd_y[k] + g_x_0_y_xxxxxy[k];

            g_x_0_yy_xxxxy[k] = -g_x_0_y_xxxxy[k] * cd_y[k] + g_x_0_y_xxxxyy[k];

            g_x_0_yy_xxxxz[k] = -g_x_0_y_xxxxz[k] * cd_y[k] + g_x_0_y_xxxxyz[k];

            g_x_0_yy_xxxyy[k] = -g_x_0_y_xxxyy[k] * cd_y[k] + g_x_0_y_xxxyyy[k];

            g_x_0_yy_xxxyz[k] = -g_x_0_y_xxxyz[k] * cd_y[k] + g_x_0_y_xxxyyz[k];

            g_x_0_yy_xxxzz[k] = -g_x_0_y_xxxzz[k] * cd_y[k] + g_x_0_y_xxxyzz[k];

            g_x_0_yy_xxyyy[k] = -g_x_0_y_xxyyy[k] * cd_y[k] + g_x_0_y_xxyyyy[k];

            g_x_0_yy_xxyyz[k] = -g_x_0_y_xxyyz[k] * cd_y[k] + g_x_0_y_xxyyyz[k];

            g_x_0_yy_xxyzz[k] = -g_x_0_y_xxyzz[k] * cd_y[k] + g_x_0_y_xxyyzz[k];

            g_x_0_yy_xxzzz[k] = -g_x_0_y_xxzzz[k] * cd_y[k] + g_x_0_y_xxyzzz[k];

            g_x_0_yy_xyyyy[k] = -g_x_0_y_xyyyy[k] * cd_y[k] + g_x_0_y_xyyyyy[k];

            g_x_0_yy_xyyyz[k] = -g_x_0_y_xyyyz[k] * cd_y[k] + g_x_0_y_xyyyyz[k];

            g_x_0_yy_xyyzz[k] = -g_x_0_y_xyyzz[k] * cd_y[k] + g_x_0_y_xyyyzz[k];

            g_x_0_yy_xyzzz[k] = -g_x_0_y_xyzzz[k] * cd_y[k] + g_x_0_y_xyyzzz[k];

            g_x_0_yy_xzzzz[k] = -g_x_0_y_xzzzz[k] * cd_y[k] + g_x_0_y_xyzzzz[k];

            g_x_0_yy_yyyyy[k] = -g_x_0_y_yyyyy[k] * cd_y[k] + g_x_0_y_yyyyyy[k];

            g_x_0_yy_yyyyz[k] = -g_x_0_y_yyyyz[k] * cd_y[k] + g_x_0_y_yyyyyz[k];

            g_x_0_yy_yyyzz[k] = -g_x_0_y_yyyzz[k] * cd_y[k] + g_x_0_y_yyyyzz[k];

            g_x_0_yy_yyzzz[k] = -g_x_0_y_yyzzz[k] * cd_y[k] + g_x_0_y_yyyzzz[k];

            g_x_0_yy_yzzzz[k] = -g_x_0_y_yzzzz[k] * cd_y[k] + g_x_0_y_yyzzzz[k];

            g_x_0_yy_zzzzz[k] = -g_x_0_y_zzzzz[k] * cd_y[k] + g_x_0_y_yzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_x_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 89);

        auto g_x_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 99);

        auto g_x_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 100);

        auto g_x_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 101);

        auto g_x_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 102);

        auto g_x_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 103);

        auto g_x_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_x_0_yz_xxxxx, g_x_0_yz_xxxxy, g_x_0_yz_xxxxz, g_x_0_yz_xxxyy, g_x_0_yz_xxxyz, g_x_0_yz_xxxzz, g_x_0_yz_xxyyy, g_x_0_yz_xxyyz, g_x_0_yz_xxyzz, g_x_0_yz_xxzzz, g_x_0_yz_xyyyy, g_x_0_yz_xyyyz, g_x_0_yz_xyyzz, g_x_0_yz_xyzzz, g_x_0_yz_xzzzz, g_x_0_yz_yyyyy, g_x_0_yz_yyyyz, g_x_0_yz_yyyzz, g_x_0_yz_yyzzz, g_x_0_yz_yzzzz, g_x_0_yz_zzzzz, g_x_0_z_xxxxx, g_x_0_z_xxxxxy, g_x_0_z_xxxxy, g_x_0_z_xxxxyy, g_x_0_z_xxxxyz, g_x_0_z_xxxxz, g_x_0_z_xxxyy, g_x_0_z_xxxyyy, g_x_0_z_xxxyyz, g_x_0_z_xxxyz, g_x_0_z_xxxyzz, g_x_0_z_xxxzz, g_x_0_z_xxyyy, g_x_0_z_xxyyyy, g_x_0_z_xxyyyz, g_x_0_z_xxyyz, g_x_0_z_xxyyzz, g_x_0_z_xxyzz, g_x_0_z_xxyzzz, g_x_0_z_xxzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyyy, g_x_0_z_xyyyyz, g_x_0_z_xyyyz, g_x_0_z_xyyyzz, g_x_0_z_xyyzz, g_x_0_z_xyyzzz, g_x_0_z_xyzzz, g_x_0_z_xyzzzz, g_x_0_z_xzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyyy, g_x_0_z_yyyyyz, g_x_0_z_yyyyz, g_x_0_z_yyyyzz, g_x_0_z_yyyzz, g_x_0_z_yyyzzz, g_x_0_z_yyzzz, g_x_0_z_yyzzzz, g_x_0_z_yzzzz, g_x_0_z_yzzzzz, g_x_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yz_xxxxx[k] = -g_x_0_z_xxxxx[k] * cd_y[k] + g_x_0_z_xxxxxy[k];

            g_x_0_yz_xxxxy[k] = -g_x_0_z_xxxxy[k] * cd_y[k] + g_x_0_z_xxxxyy[k];

            g_x_0_yz_xxxxz[k] = -g_x_0_z_xxxxz[k] * cd_y[k] + g_x_0_z_xxxxyz[k];

            g_x_0_yz_xxxyy[k] = -g_x_0_z_xxxyy[k] * cd_y[k] + g_x_0_z_xxxyyy[k];

            g_x_0_yz_xxxyz[k] = -g_x_0_z_xxxyz[k] * cd_y[k] + g_x_0_z_xxxyyz[k];

            g_x_0_yz_xxxzz[k] = -g_x_0_z_xxxzz[k] * cd_y[k] + g_x_0_z_xxxyzz[k];

            g_x_0_yz_xxyyy[k] = -g_x_0_z_xxyyy[k] * cd_y[k] + g_x_0_z_xxyyyy[k];

            g_x_0_yz_xxyyz[k] = -g_x_0_z_xxyyz[k] * cd_y[k] + g_x_0_z_xxyyyz[k];

            g_x_0_yz_xxyzz[k] = -g_x_0_z_xxyzz[k] * cd_y[k] + g_x_0_z_xxyyzz[k];

            g_x_0_yz_xxzzz[k] = -g_x_0_z_xxzzz[k] * cd_y[k] + g_x_0_z_xxyzzz[k];

            g_x_0_yz_xyyyy[k] = -g_x_0_z_xyyyy[k] * cd_y[k] + g_x_0_z_xyyyyy[k];

            g_x_0_yz_xyyyz[k] = -g_x_0_z_xyyyz[k] * cd_y[k] + g_x_0_z_xyyyyz[k];

            g_x_0_yz_xyyzz[k] = -g_x_0_z_xyyzz[k] * cd_y[k] + g_x_0_z_xyyyzz[k];

            g_x_0_yz_xyzzz[k] = -g_x_0_z_xyzzz[k] * cd_y[k] + g_x_0_z_xyyzzz[k];

            g_x_0_yz_xzzzz[k] = -g_x_0_z_xzzzz[k] * cd_y[k] + g_x_0_z_xyzzzz[k];

            g_x_0_yz_yyyyy[k] = -g_x_0_z_yyyyy[k] * cd_y[k] + g_x_0_z_yyyyyy[k];

            g_x_0_yz_yyyyz[k] = -g_x_0_z_yyyyz[k] * cd_y[k] + g_x_0_z_yyyyyz[k];

            g_x_0_yz_yyyzz[k] = -g_x_0_z_yyyzz[k] * cd_y[k] + g_x_0_z_yyyyzz[k];

            g_x_0_yz_yyzzz[k] = -g_x_0_z_yyzzz[k] * cd_y[k] + g_x_0_z_yyyzzz[k];

            g_x_0_yz_yzzzz[k] = -g_x_0_z_yzzzz[k] * cd_y[k] + g_x_0_z_yyzzzz[k];

            g_x_0_yz_zzzzz[k] = -g_x_0_z_zzzzz[k] * cd_y[k] + g_x_0_z_yzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_x_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps  + 105);

        auto g_x_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 106);

        auto g_x_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 107);

        auto g_x_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 108);

        auto g_x_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 109);

        auto g_x_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 110);

        auto g_x_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 111);

        auto g_x_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 112);

        auto g_x_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 113);

        auto g_x_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 114);

        auto g_x_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 115);

        auto g_x_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 116);

        auto g_x_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 117);

        auto g_x_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 118);

        auto g_x_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 119);

        auto g_x_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps  + 120);

        auto g_x_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 121);

        auto g_x_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 122);

        auto g_x_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 123);

        auto g_x_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 124);

        auto g_x_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps  + 125);

        #pragma omp simd aligned(cd_z, g_x_0_z_xxxxx, g_x_0_z_xxxxxz, g_x_0_z_xxxxy, g_x_0_z_xxxxyz, g_x_0_z_xxxxz, g_x_0_z_xxxxzz, g_x_0_z_xxxyy, g_x_0_z_xxxyyz, g_x_0_z_xxxyz, g_x_0_z_xxxyzz, g_x_0_z_xxxzz, g_x_0_z_xxxzzz, g_x_0_z_xxyyy, g_x_0_z_xxyyyz, g_x_0_z_xxyyz, g_x_0_z_xxyyzz, g_x_0_z_xxyzz, g_x_0_z_xxyzzz, g_x_0_z_xxzzz, g_x_0_z_xxzzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyyz, g_x_0_z_xyyyz, g_x_0_z_xyyyzz, g_x_0_z_xyyzz, g_x_0_z_xyyzzz, g_x_0_z_xyzzz, g_x_0_z_xyzzzz, g_x_0_z_xzzzz, g_x_0_z_xzzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyyz, g_x_0_z_yyyyz, g_x_0_z_yyyyzz, g_x_0_z_yyyzz, g_x_0_z_yyyzzz, g_x_0_z_yyzzz, g_x_0_z_yyzzzz, g_x_0_z_yzzzz, g_x_0_z_yzzzzz, g_x_0_z_zzzzz, g_x_0_z_zzzzzz, g_x_0_zz_xxxxx, g_x_0_zz_xxxxy, g_x_0_zz_xxxxz, g_x_0_zz_xxxyy, g_x_0_zz_xxxyz, g_x_0_zz_xxxzz, g_x_0_zz_xxyyy, g_x_0_zz_xxyyz, g_x_0_zz_xxyzz, g_x_0_zz_xxzzz, g_x_0_zz_xyyyy, g_x_0_zz_xyyyz, g_x_0_zz_xyyzz, g_x_0_zz_xyzzz, g_x_0_zz_xzzzz, g_x_0_zz_yyyyy, g_x_0_zz_yyyyz, g_x_0_zz_yyyzz, g_x_0_zz_yyzzz, g_x_0_zz_yzzzz, g_x_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zz_xxxxx[k] = -g_x_0_z_xxxxx[k] * cd_z[k] + g_x_0_z_xxxxxz[k];

            g_x_0_zz_xxxxy[k] = -g_x_0_z_xxxxy[k] * cd_z[k] + g_x_0_z_xxxxyz[k];

            g_x_0_zz_xxxxz[k] = -g_x_0_z_xxxxz[k] * cd_z[k] + g_x_0_z_xxxxzz[k];

            g_x_0_zz_xxxyy[k] = -g_x_0_z_xxxyy[k] * cd_z[k] + g_x_0_z_xxxyyz[k];

            g_x_0_zz_xxxyz[k] = -g_x_0_z_xxxyz[k] * cd_z[k] + g_x_0_z_xxxyzz[k];

            g_x_0_zz_xxxzz[k] = -g_x_0_z_xxxzz[k] * cd_z[k] + g_x_0_z_xxxzzz[k];

            g_x_0_zz_xxyyy[k] = -g_x_0_z_xxyyy[k] * cd_z[k] + g_x_0_z_xxyyyz[k];

            g_x_0_zz_xxyyz[k] = -g_x_0_z_xxyyz[k] * cd_z[k] + g_x_0_z_xxyyzz[k];

            g_x_0_zz_xxyzz[k] = -g_x_0_z_xxyzz[k] * cd_z[k] + g_x_0_z_xxyzzz[k];

            g_x_0_zz_xxzzz[k] = -g_x_0_z_xxzzz[k] * cd_z[k] + g_x_0_z_xxzzzz[k];

            g_x_0_zz_xyyyy[k] = -g_x_0_z_xyyyy[k] * cd_z[k] + g_x_0_z_xyyyyz[k];

            g_x_0_zz_xyyyz[k] = -g_x_0_z_xyyyz[k] * cd_z[k] + g_x_0_z_xyyyzz[k];

            g_x_0_zz_xyyzz[k] = -g_x_0_z_xyyzz[k] * cd_z[k] + g_x_0_z_xyyzzz[k];

            g_x_0_zz_xyzzz[k] = -g_x_0_z_xyzzz[k] * cd_z[k] + g_x_0_z_xyzzzz[k];

            g_x_0_zz_xzzzz[k] = -g_x_0_z_xzzzz[k] * cd_z[k] + g_x_0_z_xzzzzz[k];

            g_x_0_zz_yyyyy[k] = -g_x_0_z_yyyyy[k] * cd_z[k] + g_x_0_z_yyyyyz[k];

            g_x_0_zz_yyyyz[k] = -g_x_0_z_yyyyz[k] * cd_z[k] + g_x_0_z_yyyyzz[k];

            g_x_0_zz_yyyzz[k] = -g_x_0_z_yyyzz[k] * cd_z[k] + g_x_0_z_yyyzzz[k];

            g_x_0_zz_yyzzz[k] = -g_x_0_z_yyzzz[k] * cd_z[k] + g_x_0_z_yyzzzz[k];

            g_x_0_zz_yzzzz[k] = -g_x_0_z_yzzzz[k] * cd_z[k] + g_x_0_z_yzzzzz[k];

            g_x_0_zz_zzzzz[k] = -g_x_0_z_zzzzz[k] * cd_z[k] + g_x_0_z_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_y_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 0);

        auto g_y_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 1);

        auto g_y_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 2);

        auto g_y_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 3);

        auto g_y_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 4);

        auto g_y_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 5);

        auto g_y_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 6);

        auto g_y_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 7);

        auto g_y_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 8);

        auto g_y_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 9);

        auto g_y_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 10);

        auto g_y_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 11);

        auto g_y_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 12);

        auto g_y_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 13);

        auto g_y_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 14);

        auto g_y_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 15);

        auto g_y_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 16);

        auto g_y_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 17);

        auto g_y_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 18);

        auto g_y_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 19);

        auto g_y_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_y_0_x_xxxxx, g_y_0_x_xxxxxx, g_y_0_x_xxxxxy, g_y_0_x_xxxxxz, g_y_0_x_xxxxy, g_y_0_x_xxxxyy, g_y_0_x_xxxxyz, g_y_0_x_xxxxz, g_y_0_x_xxxxzz, g_y_0_x_xxxyy, g_y_0_x_xxxyyy, g_y_0_x_xxxyyz, g_y_0_x_xxxyz, g_y_0_x_xxxyzz, g_y_0_x_xxxzz, g_y_0_x_xxxzzz, g_y_0_x_xxyyy, g_y_0_x_xxyyyy, g_y_0_x_xxyyyz, g_y_0_x_xxyyz, g_y_0_x_xxyyzz, g_y_0_x_xxyzz, g_y_0_x_xxyzzz, g_y_0_x_xxzzz, g_y_0_x_xxzzzz, g_y_0_x_xyyyy, g_y_0_x_xyyyyy, g_y_0_x_xyyyyz, g_y_0_x_xyyyz, g_y_0_x_xyyyzz, g_y_0_x_xyyzz, g_y_0_x_xyyzzz, g_y_0_x_xyzzz, g_y_0_x_xyzzzz, g_y_0_x_xzzzz, g_y_0_x_xzzzzz, g_y_0_x_yyyyy, g_y_0_x_yyyyz, g_y_0_x_yyyzz, g_y_0_x_yyzzz, g_y_0_x_yzzzz, g_y_0_x_zzzzz, g_y_0_xx_xxxxx, g_y_0_xx_xxxxy, g_y_0_xx_xxxxz, g_y_0_xx_xxxyy, g_y_0_xx_xxxyz, g_y_0_xx_xxxzz, g_y_0_xx_xxyyy, g_y_0_xx_xxyyz, g_y_0_xx_xxyzz, g_y_0_xx_xxzzz, g_y_0_xx_xyyyy, g_y_0_xx_xyyyz, g_y_0_xx_xyyzz, g_y_0_xx_xyzzz, g_y_0_xx_xzzzz, g_y_0_xx_yyyyy, g_y_0_xx_yyyyz, g_y_0_xx_yyyzz, g_y_0_xx_yyzzz, g_y_0_xx_yzzzz, g_y_0_xx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xx_xxxxx[k] = -g_y_0_x_xxxxx[k] * cd_x[k] + g_y_0_x_xxxxxx[k];

            g_y_0_xx_xxxxy[k] = -g_y_0_x_xxxxy[k] * cd_x[k] + g_y_0_x_xxxxxy[k];

            g_y_0_xx_xxxxz[k] = -g_y_0_x_xxxxz[k] * cd_x[k] + g_y_0_x_xxxxxz[k];

            g_y_0_xx_xxxyy[k] = -g_y_0_x_xxxyy[k] * cd_x[k] + g_y_0_x_xxxxyy[k];

            g_y_0_xx_xxxyz[k] = -g_y_0_x_xxxyz[k] * cd_x[k] + g_y_0_x_xxxxyz[k];

            g_y_0_xx_xxxzz[k] = -g_y_0_x_xxxzz[k] * cd_x[k] + g_y_0_x_xxxxzz[k];

            g_y_0_xx_xxyyy[k] = -g_y_0_x_xxyyy[k] * cd_x[k] + g_y_0_x_xxxyyy[k];

            g_y_0_xx_xxyyz[k] = -g_y_0_x_xxyyz[k] * cd_x[k] + g_y_0_x_xxxyyz[k];

            g_y_0_xx_xxyzz[k] = -g_y_0_x_xxyzz[k] * cd_x[k] + g_y_0_x_xxxyzz[k];

            g_y_0_xx_xxzzz[k] = -g_y_0_x_xxzzz[k] * cd_x[k] + g_y_0_x_xxxzzz[k];

            g_y_0_xx_xyyyy[k] = -g_y_0_x_xyyyy[k] * cd_x[k] + g_y_0_x_xxyyyy[k];

            g_y_0_xx_xyyyz[k] = -g_y_0_x_xyyyz[k] * cd_x[k] + g_y_0_x_xxyyyz[k];

            g_y_0_xx_xyyzz[k] = -g_y_0_x_xyyzz[k] * cd_x[k] + g_y_0_x_xxyyzz[k];

            g_y_0_xx_xyzzz[k] = -g_y_0_x_xyzzz[k] * cd_x[k] + g_y_0_x_xxyzzz[k];

            g_y_0_xx_xzzzz[k] = -g_y_0_x_xzzzz[k] * cd_x[k] + g_y_0_x_xxzzzz[k];

            g_y_0_xx_yyyyy[k] = -g_y_0_x_yyyyy[k] * cd_x[k] + g_y_0_x_xyyyyy[k];

            g_y_0_xx_yyyyz[k] = -g_y_0_x_yyyyz[k] * cd_x[k] + g_y_0_x_xyyyyz[k];

            g_y_0_xx_yyyzz[k] = -g_y_0_x_yyyzz[k] * cd_x[k] + g_y_0_x_xyyyzz[k];

            g_y_0_xx_yyzzz[k] = -g_y_0_x_yyzzz[k] * cd_x[k] + g_y_0_x_xyyzzz[k];

            g_y_0_xx_yzzzz[k] = -g_y_0_x_yzzzz[k] * cd_x[k] + g_y_0_x_xyzzzz[k];

            g_y_0_xx_zzzzz[k] = -g_y_0_x_zzzzz[k] * cd_x[k] + g_y_0_x_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_y_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 21);

        auto g_y_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 22);

        auto g_y_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 23);

        auto g_y_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 24);

        auto g_y_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 25);

        auto g_y_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 26);

        auto g_y_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 27);

        auto g_y_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 28);

        auto g_y_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 29);

        auto g_y_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 30);

        auto g_y_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 31);

        auto g_y_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 32);

        auto g_y_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 33);

        auto g_y_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 34);

        auto g_y_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 35);

        auto g_y_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 36);

        auto g_y_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 37);

        auto g_y_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 38);

        auto g_y_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 39);

        auto g_y_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 40);

        auto g_y_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_y_0_xy_xxxxx, g_y_0_xy_xxxxy, g_y_0_xy_xxxxz, g_y_0_xy_xxxyy, g_y_0_xy_xxxyz, g_y_0_xy_xxxzz, g_y_0_xy_xxyyy, g_y_0_xy_xxyyz, g_y_0_xy_xxyzz, g_y_0_xy_xxzzz, g_y_0_xy_xyyyy, g_y_0_xy_xyyyz, g_y_0_xy_xyyzz, g_y_0_xy_xyzzz, g_y_0_xy_xzzzz, g_y_0_xy_yyyyy, g_y_0_xy_yyyyz, g_y_0_xy_yyyzz, g_y_0_xy_yyzzz, g_y_0_xy_yzzzz, g_y_0_xy_zzzzz, g_y_0_y_xxxxx, g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxy, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxz, g_y_0_y_xxxxzz, g_y_0_y_xxxyy, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyy, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xy_xxxxx[k] = -g_y_0_y_xxxxx[k] * cd_x[k] + g_y_0_y_xxxxxx[k];

            g_y_0_xy_xxxxy[k] = -g_y_0_y_xxxxy[k] * cd_x[k] + g_y_0_y_xxxxxy[k];

            g_y_0_xy_xxxxz[k] = -g_y_0_y_xxxxz[k] * cd_x[k] + g_y_0_y_xxxxxz[k];

            g_y_0_xy_xxxyy[k] = -g_y_0_y_xxxyy[k] * cd_x[k] + g_y_0_y_xxxxyy[k];

            g_y_0_xy_xxxyz[k] = -g_y_0_y_xxxyz[k] * cd_x[k] + g_y_0_y_xxxxyz[k];

            g_y_0_xy_xxxzz[k] = -g_y_0_y_xxxzz[k] * cd_x[k] + g_y_0_y_xxxxzz[k];

            g_y_0_xy_xxyyy[k] = -g_y_0_y_xxyyy[k] * cd_x[k] + g_y_0_y_xxxyyy[k];

            g_y_0_xy_xxyyz[k] = -g_y_0_y_xxyyz[k] * cd_x[k] + g_y_0_y_xxxyyz[k];

            g_y_0_xy_xxyzz[k] = -g_y_0_y_xxyzz[k] * cd_x[k] + g_y_0_y_xxxyzz[k];

            g_y_0_xy_xxzzz[k] = -g_y_0_y_xxzzz[k] * cd_x[k] + g_y_0_y_xxxzzz[k];

            g_y_0_xy_xyyyy[k] = -g_y_0_y_xyyyy[k] * cd_x[k] + g_y_0_y_xxyyyy[k];

            g_y_0_xy_xyyyz[k] = -g_y_0_y_xyyyz[k] * cd_x[k] + g_y_0_y_xxyyyz[k];

            g_y_0_xy_xyyzz[k] = -g_y_0_y_xyyzz[k] * cd_x[k] + g_y_0_y_xxyyzz[k];

            g_y_0_xy_xyzzz[k] = -g_y_0_y_xyzzz[k] * cd_x[k] + g_y_0_y_xxyzzz[k];

            g_y_0_xy_xzzzz[k] = -g_y_0_y_xzzzz[k] * cd_x[k] + g_y_0_y_xxzzzz[k];

            g_y_0_xy_yyyyy[k] = -g_y_0_y_yyyyy[k] * cd_x[k] + g_y_0_y_xyyyyy[k];

            g_y_0_xy_yyyyz[k] = -g_y_0_y_yyyyz[k] * cd_x[k] + g_y_0_y_xyyyyz[k];

            g_y_0_xy_yyyzz[k] = -g_y_0_y_yyyzz[k] * cd_x[k] + g_y_0_y_xyyyzz[k];

            g_y_0_xy_yyzzz[k] = -g_y_0_y_yyzzz[k] * cd_x[k] + g_y_0_y_xyyzzz[k];

            g_y_0_xy_yzzzz[k] = -g_y_0_y_yzzzz[k] * cd_x[k] + g_y_0_y_xyzzzz[k];

            g_y_0_xy_zzzzz[k] = -g_y_0_y_zzzzz[k] * cd_x[k] + g_y_0_y_xzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_y_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 42);

        auto g_y_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 43);

        auto g_y_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 44);

        auto g_y_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 45);

        auto g_y_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 46);

        auto g_y_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 47);

        auto g_y_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 48);

        auto g_y_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 49);

        auto g_y_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 50);

        auto g_y_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 51);

        auto g_y_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 52);

        auto g_y_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 53);

        auto g_y_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 54);

        auto g_y_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 55);

        auto g_y_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 56);

        auto g_y_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 57);

        auto g_y_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 58);

        auto g_y_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 59);

        auto g_y_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 60);

        auto g_y_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 61);

        auto g_y_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 62);

        #pragma omp simd aligned(cd_x, g_y_0_xz_xxxxx, g_y_0_xz_xxxxy, g_y_0_xz_xxxxz, g_y_0_xz_xxxyy, g_y_0_xz_xxxyz, g_y_0_xz_xxxzz, g_y_0_xz_xxyyy, g_y_0_xz_xxyyz, g_y_0_xz_xxyzz, g_y_0_xz_xxzzz, g_y_0_xz_xyyyy, g_y_0_xz_xyyyz, g_y_0_xz_xyyzz, g_y_0_xz_xyzzz, g_y_0_xz_xzzzz, g_y_0_xz_yyyyy, g_y_0_xz_yyyyz, g_y_0_xz_yyyzz, g_y_0_xz_yyzzz, g_y_0_xz_yzzzz, g_y_0_xz_zzzzz, g_y_0_z_xxxxx, g_y_0_z_xxxxxx, g_y_0_z_xxxxxy, g_y_0_z_xxxxxz, g_y_0_z_xxxxy, g_y_0_z_xxxxyy, g_y_0_z_xxxxyz, g_y_0_z_xxxxz, g_y_0_z_xxxxzz, g_y_0_z_xxxyy, g_y_0_z_xxxyyy, g_y_0_z_xxxyyz, g_y_0_z_xxxyz, g_y_0_z_xxxyzz, g_y_0_z_xxxzz, g_y_0_z_xxxzzz, g_y_0_z_xxyyy, g_y_0_z_xxyyyy, g_y_0_z_xxyyyz, g_y_0_z_xxyyz, g_y_0_z_xxyyzz, g_y_0_z_xxyzz, g_y_0_z_xxyzzz, g_y_0_z_xxzzz, g_y_0_z_xxzzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyyy, g_y_0_z_xyyyyz, g_y_0_z_xyyyz, g_y_0_z_xyyyzz, g_y_0_z_xyyzz, g_y_0_z_xyyzzz, g_y_0_z_xyzzz, g_y_0_z_xyzzzz, g_y_0_z_xzzzz, g_y_0_z_xzzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyz, g_y_0_z_yyyzz, g_y_0_z_yyzzz, g_y_0_z_yzzzz, g_y_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xz_xxxxx[k] = -g_y_0_z_xxxxx[k] * cd_x[k] + g_y_0_z_xxxxxx[k];

            g_y_0_xz_xxxxy[k] = -g_y_0_z_xxxxy[k] * cd_x[k] + g_y_0_z_xxxxxy[k];

            g_y_0_xz_xxxxz[k] = -g_y_0_z_xxxxz[k] * cd_x[k] + g_y_0_z_xxxxxz[k];

            g_y_0_xz_xxxyy[k] = -g_y_0_z_xxxyy[k] * cd_x[k] + g_y_0_z_xxxxyy[k];

            g_y_0_xz_xxxyz[k] = -g_y_0_z_xxxyz[k] * cd_x[k] + g_y_0_z_xxxxyz[k];

            g_y_0_xz_xxxzz[k] = -g_y_0_z_xxxzz[k] * cd_x[k] + g_y_0_z_xxxxzz[k];

            g_y_0_xz_xxyyy[k] = -g_y_0_z_xxyyy[k] * cd_x[k] + g_y_0_z_xxxyyy[k];

            g_y_0_xz_xxyyz[k] = -g_y_0_z_xxyyz[k] * cd_x[k] + g_y_0_z_xxxyyz[k];

            g_y_0_xz_xxyzz[k] = -g_y_0_z_xxyzz[k] * cd_x[k] + g_y_0_z_xxxyzz[k];

            g_y_0_xz_xxzzz[k] = -g_y_0_z_xxzzz[k] * cd_x[k] + g_y_0_z_xxxzzz[k];

            g_y_0_xz_xyyyy[k] = -g_y_0_z_xyyyy[k] * cd_x[k] + g_y_0_z_xxyyyy[k];

            g_y_0_xz_xyyyz[k] = -g_y_0_z_xyyyz[k] * cd_x[k] + g_y_0_z_xxyyyz[k];

            g_y_0_xz_xyyzz[k] = -g_y_0_z_xyyzz[k] * cd_x[k] + g_y_0_z_xxyyzz[k];

            g_y_0_xz_xyzzz[k] = -g_y_0_z_xyzzz[k] * cd_x[k] + g_y_0_z_xxyzzz[k];

            g_y_0_xz_xzzzz[k] = -g_y_0_z_xzzzz[k] * cd_x[k] + g_y_0_z_xxzzzz[k];

            g_y_0_xz_yyyyy[k] = -g_y_0_z_yyyyy[k] * cd_x[k] + g_y_0_z_xyyyyy[k];

            g_y_0_xz_yyyyz[k] = -g_y_0_z_yyyyz[k] * cd_x[k] + g_y_0_z_xyyyyz[k];

            g_y_0_xz_yyyzz[k] = -g_y_0_z_yyyzz[k] * cd_x[k] + g_y_0_z_xyyyzz[k];

            g_y_0_xz_yyzzz[k] = -g_y_0_z_yyzzz[k] * cd_x[k] + g_y_0_z_xyyzzz[k];

            g_y_0_xz_yzzzz[k] = -g_y_0_z_yzzzz[k] * cd_x[k] + g_y_0_z_xyzzzz[k];

            g_y_0_xz_zzzzz[k] = -g_y_0_z_zzzzz[k] * cd_x[k] + g_y_0_z_xzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_y_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 63);

        auto g_y_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 64);

        auto g_y_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 65);

        auto g_y_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 66);

        auto g_y_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 67);

        auto g_y_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 68);

        auto g_y_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 69);

        auto g_y_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 70);

        auto g_y_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 71);

        auto g_y_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 72);

        auto g_y_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 73);

        auto g_y_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 74);

        auto g_y_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 75);

        auto g_y_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 76);

        auto g_y_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 77);

        auto g_y_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 78);

        auto g_y_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 79);

        auto g_y_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 80);

        auto g_y_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 81);

        auto g_y_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 82);

        auto g_y_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_y_0_y_xxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxy, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzz, g_y_0_yy_xxxxx, g_y_0_yy_xxxxy, g_y_0_yy_xxxxz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyz, g_y_0_yy_xxxzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyz, g_y_0_yy_xxyzz, g_y_0_yy_xxzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyz, g_y_0_yy_xyyzz, g_y_0_yy_xyzzz, g_y_0_yy_xzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyz, g_y_0_yy_yyyzz, g_y_0_yy_yyzzz, g_y_0_yy_yzzzz, g_y_0_yy_zzzzz, g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxyy, g_y_xxxyz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyz, g_y_xxyzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyz, g_y_xyyzz, g_y_xyzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yy_xxxxx[k] = -g_y_xxxxx[k] - g_y_0_y_xxxxx[k] * cd_y[k] + g_y_0_y_xxxxxy[k];

            g_y_0_yy_xxxxy[k] = -g_y_xxxxy[k] - g_y_0_y_xxxxy[k] * cd_y[k] + g_y_0_y_xxxxyy[k];

            g_y_0_yy_xxxxz[k] = -g_y_xxxxz[k] - g_y_0_y_xxxxz[k] * cd_y[k] + g_y_0_y_xxxxyz[k];

            g_y_0_yy_xxxyy[k] = -g_y_xxxyy[k] - g_y_0_y_xxxyy[k] * cd_y[k] + g_y_0_y_xxxyyy[k];

            g_y_0_yy_xxxyz[k] = -g_y_xxxyz[k] - g_y_0_y_xxxyz[k] * cd_y[k] + g_y_0_y_xxxyyz[k];

            g_y_0_yy_xxxzz[k] = -g_y_xxxzz[k] - g_y_0_y_xxxzz[k] * cd_y[k] + g_y_0_y_xxxyzz[k];

            g_y_0_yy_xxyyy[k] = -g_y_xxyyy[k] - g_y_0_y_xxyyy[k] * cd_y[k] + g_y_0_y_xxyyyy[k];

            g_y_0_yy_xxyyz[k] = -g_y_xxyyz[k] - g_y_0_y_xxyyz[k] * cd_y[k] + g_y_0_y_xxyyyz[k];

            g_y_0_yy_xxyzz[k] = -g_y_xxyzz[k] - g_y_0_y_xxyzz[k] * cd_y[k] + g_y_0_y_xxyyzz[k];

            g_y_0_yy_xxzzz[k] = -g_y_xxzzz[k] - g_y_0_y_xxzzz[k] * cd_y[k] + g_y_0_y_xxyzzz[k];

            g_y_0_yy_xyyyy[k] = -g_y_xyyyy[k] - g_y_0_y_xyyyy[k] * cd_y[k] + g_y_0_y_xyyyyy[k];

            g_y_0_yy_xyyyz[k] = -g_y_xyyyz[k] - g_y_0_y_xyyyz[k] * cd_y[k] + g_y_0_y_xyyyyz[k];

            g_y_0_yy_xyyzz[k] = -g_y_xyyzz[k] - g_y_0_y_xyyzz[k] * cd_y[k] + g_y_0_y_xyyyzz[k];

            g_y_0_yy_xyzzz[k] = -g_y_xyzzz[k] - g_y_0_y_xyzzz[k] * cd_y[k] + g_y_0_y_xyyzzz[k];

            g_y_0_yy_xzzzz[k] = -g_y_xzzzz[k] - g_y_0_y_xzzzz[k] * cd_y[k] + g_y_0_y_xyzzzz[k];

            g_y_0_yy_yyyyy[k] = -g_y_yyyyy[k] - g_y_0_y_yyyyy[k] * cd_y[k] + g_y_0_y_yyyyyy[k];

            g_y_0_yy_yyyyz[k] = -g_y_yyyyz[k] - g_y_0_y_yyyyz[k] * cd_y[k] + g_y_0_y_yyyyyz[k];

            g_y_0_yy_yyyzz[k] = -g_y_yyyzz[k] - g_y_0_y_yyyzz[k] * cd_y[k] + g_y_0_y_yyyyzz[k];

            g_y_0_yy_yyzzz[k] = -g_y_yyzzz[k] - g_y_0_y_yyzzz[k] * cd_y[k] + g_y_0_y_yyyzzz[k];

            g_y_0_yy_yzzzz[k] = -g_y_yzzzz[k] - g_y_0_y_yzzzz[k] * cd_y[k] + g_y_0_y_yyzzzz[k];

            g_y_0_yy_zzzzz[k] = -g_y_zzzzz[k] - g_y_0_y_zzzzz[k] * cd_y[k] + g_y_0_y_yzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_y_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 84);

        auto g_y_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 85);

        auto g_y_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 86);

        auto g_y_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 87);

        auto g_y_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 88);

        auto g_y_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 89);

        auto g_y_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 90);

        auto g_y_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 91);

        auto g_y_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 92);

        auto g_y_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 93);

        auto g_y_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 94);

        auto g_y_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 95);

        auto g_y_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 96);

        auto g_y_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 97);

        auto g_y_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 98);

        auto g_y_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 99);

        auto g_y_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 100);

        auto g_y_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 101);

        auto g_y_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 102);

        auto g_y_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 103);

        auto g_y_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 104);

        #pragma omp simd aligned(cd_z, g_y_0_y_xxxxx, g_y_0_y_xxxxxz, g_y_0_y_xxxxy, g_y_0_y_xxxxyz, g_y_0_y_xxxxz, g_y_0_y_xxxxzz, g_y_0_y_xxxyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzz, g_y_0_y_zzzzzz, g_y_0_yz_xxxxx, g_y_0_yz_xxxxy, g_y_0_yz_xxxxz, g_y_0_yz_xxxyy, g_y_0_yz_xxxyz, g_y_0_yz_xxxzz, g_y_0_yz_xxyyy, g_y_0_yz_xxyyz, g_y_0_yz_xxyzz, g_y_0_yz_xxzzz, g_y_0_yz_xyyyy, g_y_0_yz_xyyyz, g_y_0_yz_xyyzz, g_y_0_yz_xyzzz, g_y_0_yz_xzzzz, g_y_0_yz_yyyyy, g_y_0_yz_yyyyz, g_y_0_yz_yyyzz, g_y_0_yz_yyzzz, g_y_0_yz_yzzzz, g_y_0_yz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yz_xxxxx[k] = -g_y_0_y_xxxxx[k] * cd_z[k] + g_y_0_y_xxxxxz[k];

            g_y_0_yz_xxxxy[k] = -g_y_0_y_xxxxy[k] * cd_z[k] + g_y_0_y_xxxxyz[k];

            g_y_0_yz_xxxxz[k] = -g_y_0_y_xxxxz[k] * cd_z[k] + g_y_0_y_xxxxzz[k];

            g_y_0_yz_xxxyy[k] = -g_y_0_y_xxxyy[k] * cd_z[k] + g_y_0_y_xxxyyz[k];

            g_y_0_yz_xxxyz[k] = -g_y_0_y_xxxyz[k] * cd_z[k] + g_y_0_y_xxxyzz[k];

            g_y_0_yz_xxxzz[k] = -g_y_0_y_xxxzz[k] * cd_z[k] + g_y_0_y_xxxzzz[k];

            g_y_0_yz_xxyyy[k] = -g_y_0_y_xxyyy[k] * cd_z[k] + g_y_0_y_xxyyyz[k];

            g_y_0_yz_xxyyz[k] = -g_y_0_y_xxyyz[k] * cd_z[k] + g_y_0_y_xxyyzz[k];

            g_y_0_yz_xxyzz[k] = -g_y_0_y_xxyzz[k] * cd_z[k] + g_y_0_y_xxyzzz[k];

            g_y_0_yz_xxzzz[k] = -g_y_0_y_xxzzz[k] * cd_z[k] + g_y_0_y_xxzzzz[k];

            g_y_0_yz_xyyyy[k] = -g_y_0_y_xyyyy[k] * cd_z[k] + g_y_0_y_xyyyyz[k];

            g_y_0_yz_xyyyz[k] = -g_y_0_y_xyyyz[k] * cd_z[k] + g_y_0_y_xyyyzz[k];

            g_y_0_yz_xyyzz[k] = -g_y_0_y_xyyzz[k] * cd_z[k] + g_y_0_y_xyyzzz[k];

            g_y_0_yz_xyzzz[k] = -g_y_0_y_xyzzz[k] * cd_z[k] + g_y_0_y_xyzzzz[k];

            g_y_0_yz_xzzzz[k] = -g_y_0_y_xzzzz[k] * cd_z[k] + g_y_0_y_xzzzzz[k];

            g_y_0_yz_yyyyy[k] = -g_y_0_y_yyyyy[k] * cd_z[k] + g_y_0_y_yyyyyz[k];

            g_y_0_yz_yyyyz[k] = -g_y_0_y_yyyyz[k] * cd_z[k] + g_y_0_y_yyyyzz[k];

            g_y_0_yz_yyyzz[k] = -g_y_0_y_yyyzz[k] * cd_z[k] + g_y_0_y_yyyzzz[k];

            g_y_0_yz_yyzzz[k] = -g_y_0_y_yyzzz[k] * cd_z[k] + g_y_0_y_yyzzzz[k];

            g_y_0_yz_yzzzz[k] = -g_y_0_y_yzzzz[k] * cd_z[k] + g_y_0_y_yzzzzz[k];

            g_y_0_yz_zzzzz[k] = -g_y_0_y_zzzzz[k] * cd_z[k] + g_y_0_y_zzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_y_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps  + 105);

        auto g_y_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 106);

        auto g_y_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 107);

        auto g_y_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 108);

        auto g_y_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 109);

        auto g_y_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 110);

        auto g_y_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 111);

        auto g_y_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 112);

        auto g_y_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 113);

        auto g_y_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 114);

        auto g_y_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 115);

        auto g_y_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 116);

        auto g_y_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 117);

        auto g_y_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 118);

        auto g_y_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 119);

        auto g_y_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps  + 120);

        auto g_y_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 121);

        auto g_y_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 122);

        auto g_y_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 123);

        auto g_y_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 124);

        auto g_y_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps  + 125);

        #pragma omp simd aligned(cd_z, g_y_0_z_xxxxx, g_y_0_z_xxxxxz, g_y_0_z_xxxxy, g_y_0_z_xxxxyz, g_y_0_z_xxxxz, g_y_0_z_xxxxzz, g_y_0_z_xxxyy, g_y_0_z_xxxyyz, g_y_0_z_xxxyz, g_y_0_z_xxxyzz, g_y_0_z_xxxzz, g_y_0_z_xxxzzz, g_y_0_z_xxyyy, g_y_0_z_xxyyyz, g_y_0_z_xxyyz, g_y_0_z_xxyyzz, g_y_0_z_xxyzz, g_y_0_z_xxyzzz, g_y_0_z_xxzzz, g_y_0_z_xxzzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyyz, g_y_0_z_xyyyz, g_y_0_z_xyyyzz, g_y_0_z_xyyzz, g_y_0_z_xyyzzz, g_y_0_z_xyzzz, g_y_0_z_xyzzzz, g_y_0_z_xzzzz, g_y_0_z_xzzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyyz, g_y_0_z_yyyyz, g_y_0_z_yyyyzz, g_y_0_z_yyyzz, g_y_0_z_yyyzzz, g_y_0_z_yyzzz, g_y_0_z_yyzzzz, g_y_0_z_yzzzz, g_y_0_z_yzzzzz, g_y_0_z_zzzzz, g_y_0_z_zzzzzz, g_y_0_zz_xxxxx, g_y_0_zz_xxxxy, g_y_0_zz_xxxxz, g_y_0_zz_xxxyy, g_y_0_zz_xxxyz, g_y_0_zz_xxxzz, g_y_0_zz_xxyyy, g_y_0_zz_xxyyz, g_y_0_zz_xxyzz, g_y_0_zz_xxzzz, g_y_0_zz_xyyyy, g_y_0_zz_xyyyz, g_y_0_zz_xyyzz, g_y_0_zz_xyzzz, g_y_0_zz_xzzzz, g_y_0_zz_yyyyy, g_y_0_zz_yyyyz, g_y_0_zz_yyyzz, g_y_0_zz_yyzzz, g_y_0_zz_yzzzz, g_y_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zz_xxxxx[k] = -g_y_0_z_xxxxx[k] * cd_z[k] + g_y_0_z_xxxxxz[k];

            g_y_0_zz_xxxxy[k] = -g_y_0_z_xxxxy[k] * cd_z[k] + g_y_0_z_xxxxyz[k];

            g_y_0_zz_xxxxz[k] = -g_y_0_z_xxxxz[k] * cd_z[k] + g_y_0_z_xxxxzz[k];

            g_y_0_zz_xxxyy[k] = -g_y_0_z_xxxyy[k] * cd_z[k] + g_y_0_z_xxxyyz[k];

            g_y_0_zz_xxxyz[k] = -g_y_0_z_xxxyz[k] * cd_z[k] + g_y_0_z_xxxyzz[k];

            g_y_0_zz_xxxzz[k] = -g_y_0_z_xxxzz[k] * cd_z[k] + g_y_0_z_xxxzzz[k];

            g_y_0_zz_xxyyy[k] = -g_y_0_z_xxyyy[k] * cd_z[k] + g_y_0_z_xxyyyz[k];

            g_y_0_zz_xxyyz[k] = -g_y_0_z_xxyyz[k] * cd_z[k] + g_y_0_z_xxyyzz[k];

            g_y_0_zz_xxyzz[k] = -g_y_0_z_xxyzz[k] * cd_z[k] + g_y_0_z_xxyzzz[k];

            g_y_0_zz_xxzzz[k] = -g_y_0_z_xxzzz[k] * cd_z[k] + g_y_0_z_xxzzzz[k];

            g_y_0_zz_xyyyy[k] = -g_y_0_z_xyyyy[k] * cd_z[k] + g_y_0_z_xyyyyz[k];

            g_y_0_zz_xyyyz[k] = -g_y_0_z_xyyyz[k] * cd_z[k] + g_y_0_z_xyyyzz[k];

            g_y_0_zz_xyyzz[k] = -g_y_0_z_xyyzz[k] * cd_z[k] + g_y_0_z_xyyzzz[k];

            g_y_0_zz_xyzzz[k] = -g_y_0_z_xyzzz[k] * cd_z[k] + g_y_0_z_xyzzzz[k];

            g_y_0_zz_xzzzz[k] = -g_y_0_z_xzzzz[k] * cd_z[k] + g_y_0_z_xzzzzz[k];

            g_y_0_zz_yyyyy[k] = -g_y_0_z_yyyyy[k] * cd_z[k] + g_y_0_z_yyyyyz[k];

            g_y_0_zz_yyyyz[k] = -g_y_0_z_yyyyz[k] * cd_z[k] + g_y_0_z_yyyyzz[k];

            g_y_0_zz_yyyzz[k] = -g_y_0_z_yyyzz[k] * cd_z[k] + g_y_0_z_yyyzzz[k];

            g_y_0_zz_yyzzz[k] = -g_y_0_z_yyzzz[k] * cd_z[k] + g_y_0_z_yyzzzz[k];

            g_y_0_zz_yzzzz[k] = -g_y_0_z_yzzzz[k] * cd_z[k] + g_y_0_z_yzzzzz[k];

            g_y_0_zz_zzzzz[k] = -g_y_0_z_zzzzz[k] * cd_z[k] + g_y_0_z_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_z_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 0);

        auto g_z_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 1);

        auto g_z_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 2);

        auto g_z_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 3);

        auto g_z_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 4);

        auto g_z_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 5);

        auto g_z_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 6);

        auto g_z_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 7);

        auto g_z_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 8);

        auto g_z_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 9);

        auto g_z_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 10);

        auto g_z_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 11);

        auto g_z_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 12);

        auto g_z_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 13);

        auto g_z_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 14);

        auto g_z_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 15);

        auto g_z_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 16);

        auto g_z_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 17);

        auto g_z_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 18);

        auto g_z_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 19);

        auto g_z_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_z_0_x_xxxxx, g_z_0_x_xxxxxx, g_z_0_x_xxxxxy, g_z_0_x_xxxxxz, g_z_0_x_xxxxy, g_z_0_x_xxxxyy, g_z_0_x_xxxxyz, g_z_0_x_xxxxz, g_z_0_x_xxxxzz, g_z_0_x_xxxyy, g_z_0_x_xxxyyy, g_z_0_x_xxxyyz, g_z_0_x_xxxyz, g_z_0_x_xxxyzz, g_z_0_x_xxxzz, g_z_0_x_xxxzzz, g_z_0_x_xxyyy, g_z_0_x_xxyyyy, g_z_0_x_xxyyyz, g_z_0_x_xxyyz, g_z_0_x_xxyyzz, g_z_0_x_xxyzz, g_z_0_x_xxyzzz, g_z_0_x_xxzzz, g_z_0_x_xxzzzz, g_z_0_x_xyyyy, g_z_0_x_xyyyyy, g_z_0_x_xyyyyz, g_z_0_x_xyyyz, g_z_0_x_xyyyzz, g_z_0_x_xyyzz, g_z_0_x_xyyzzz, g_z_0_x_xyzzz, g_z_0_x_xyzzzz, g_z_0_x_xzzzz, g_z_0_x_xzzzzz, g_z_0_x_yyyyy, g_z_0_x_yyyyz, g_z_0_x_yyyzz, g_z_0_x_yyzzz, g_z_0_x_yzzzz, g_z_0_x_zzzzz, g_z_0_xx_xxxxx, g_z_0_xx_xxxxy, g_z_0_xx_xxxxz, g_z_0_xx_xxxyy, g_z_0_xx_xxxyz, g_z_0_xx_xxxzz, g_z_0_xx_xxyyy, g_z_0_xx_xxyyz, g_z_0_xx_xxyzz, g_z_0_xx_xxzzz, g_z_0_xx_xyyyy, g_z_0_xx_xyyyz, g_z_0_xx_xyyzz, g_z_0_xx_xyzzz, g_z_0_xx_xzzzz, g_z_0_xx_yyyyy, g_z_0_xx_yyyyz, g_z_0_xx_yyyzz, g_z_0_xx_yyzzz, g_z_0_xx_yzzzz, g_z_0_xx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xx_xxxxx[k] = -g_z_0_x_xxxxx[k] * cd_x[k] + g_z_0_x_xxxxxx[k];

            g_z_0_xx_xxxxy[k] = -g_z_0_x_xxxxy[k] * cd_x[k] + g_z_0_x_xxxxxy[k];

            g_z_0_xx_xxxxz[k] = -g_z_0_x_xxxxz[k] * cd_x[k] + g_z_0_x_xxxxxz[k];

            g_z_0_xx_xxxyy[k] = -g_z_0_x_xxxyy[k] * cd_x[k] + g_z_0_x_xxxxyy[k];

            g_z_0_xx_xxxyz[k] = -g_z_0_x_xxxyz[k] * cd_x[k] + g_z_0_x_xxxxyz[k];

            g_z_0_xx_xxxzz[k] = -g_z_0_x_xxxzz[k] * cd_x[k] + g_z_0_x_xxxxzz[k];

            g_z_0_xx_xxyyy[k] = -g_z_0_x_xxyyy[k] * cd_x[k] + g_z_0_x_xxxyyy[k];

            g_z_0_xx_xxyyz[k] = -g_z_0_x_xxyyz[k] * cd_x[k] + g_z_0_x_xxxyyz[k];

            g_z_0_xx_xxyzz[k] = -g_z_0_x_xxyzz[k] * cd_x[k] + g_z_0_x_xxxyzz[k];

            g_z_0_xx_xxzzz[k] = -g_z_0_x_xxzzz[k] * cd_x[k] + g_z_0_x_xxxzzz[k];

            g_z_0_xx_xyyyy[k] = -g_z_0_x_xyyyy[k] * cd_x[k] + g_z_0_x_xxyyyy[k];

            g_z_0_xx_xyyyz[k] = -g_z_0_x_xyyyz[k] * cd_x[k] + g_z_0_x_xxyyyz[k];

            g_z_0_xx_xyyzz[k] = -g_z_0_x_xyyzz[k] * cd_x[k] + g_z_0_x_xxyyzz[k];

            g_z_0_xx_xyzzz[k] = -g_z_0_x_xyzzz[k] * cd_x[k] + g_z_0_x_xxyzzz[k];

            g_z_0_xx_xzzzz[k] = -g_z_0_x_xzzzz[k] * cd_x[k] + g_z_0_x_xxzzzz[k];

            g_z_0_xx_yyyyy[k] = -g_z_0_x_yyyyy[k] * cd_x[k] + g_z_0_x_xyyyyy[k];

            g_z_0_xx_yyyyz[k] = -g_z_0_x_yyyyz[k] * cd_x[k] + g_z_0_x_xyyyyz[k];

            g_z_0_xx_yyyzz[k] = -g_z_0_x_yyyzz[k] * cd_x[k] + g_z_0_x_xyyyzz[k];

            g_z_0_xx_yyzzz[k] = -g_z_0_x_yyzzz[k] * cd_x[k] + g_z_0_x_xyyzzz[k];

            g_z_0_xx_yzzzz[k] = -g_z_0_x_yzzzz[k] * cd_x[k] + g_z_0_x_xyzzzz[k];

            g_z_0_xx_zzzzz[k] = -g_z_0_x_zzzzz[k] * cd_x[k] + g_z_0_x_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_z_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 21);

        auto g_z_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 22);

        auto g_z_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 23);

        auto g_z_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 24);

        auto g_z_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 25);

        auto g_z_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 26);

        auto g_z_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 27);

        auto g_z_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 28);

        auto g_z_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 29);

        auto g_z_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 30);

        auto g_z_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 31);

        auto g_z_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 32);

        auto g_z_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 33);

        auto g_z_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 34);

        auto g_z_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 35);

        auto g_z_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 36);

        auto g_z_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 37);

        auto g_z_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 38);

        auto g_z_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 39);

        auto g_z_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 40);

        auto g_z_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_z_0_xy_xxxxx, g_z_0_xy_xxxxy, g_z_0_xy_xxxxz, g_z_0_xy_xxxyy, g_z_0_xy_xxxyz, g_z_0_xy_xxxzz, g_z_0_xy_xxyyy, g_z_0_xy_xxyyz, g_z_0_xy_xxyzz, g_z_0_xy_xxzzz, g_z_0_xy_xyyyy, g_z_0_xy_xyyyz, g_z_0_xy_xyyzz, g_z_0_xy_xyzzz, g_z_0_xy_xzzzz, g_z_0_xy_yyyyy, g_z_0_xy_yyyyz, g_z_0_xy_yyyzz, g_z_0_xy_yyzzz, g_z_0_xy_yzzzz, g_z_0_xy_zzzzz, g_z_0_y_xxxxx, g_z_0_y_xxxxxx, g_z_0_y_xxxxxy, g_z_0_y_xxxxxz, g_z_0_y_xxxxy, g_z_0_y_xxxxyy, g_z_0_y_xxxxyz, g_z_0_y_xxxxz, g_z_0_y_xxxxzz, g_z_0_y_xxxyy, g_z_0_y_xxxyyy, g_z_0_y_xxxyyz, g_z_0_y_xxxyz, g_z_0_y_xxxyzz, g_z_0_y_xxxzz, g_z_0_y_xxxzzz, g_z_0_y_xxyyy, g_z_0_y_xxyyyy, g_z_0_y_xxyyyz, g_z_0_y_xxyyz, g_z_0_y_xxyyzz, g_z_0_y_xxyzz, g_z_0_y_xxyzzz, g_z_0_y_xxzzz, g_z_0_y_xxzzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyyy, g_z_0_y_xyyyyz, g_z_0_y_xyyyz, g_z_0_y_xyyyzz, g_z_0_y_xyyzz, g_z_0_y_xyyzzz, g_z_0_y_xyzzz, g_z_0_y_xyzzzz, g_z_0_y_xzzzz, g_z_0_y_xzzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyzz, g_z_0_y_yyzzz, g_z_0_y_yzzzz, g_z_0_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xy_xxxxx[k] = -g_z_0_y_xxxxx[k] * cd_x[k] + g_z_0_y_xxxxxx[k];

            g_z_0_xy_xxxxy[k] = -g_z_0_y_xxxxy[k] * cd_x[k] + g_z_0_y_xxxxxy[k];

            g_z_0_xy_xxxxz[k] = -g_z_0_y_xxxxz[k] * cd_x[k] + g_z_0_y_xxxxxz[k];

            g_z_0_xy_xxxyy[k] = -g_z_0_y_xxxyy[k] * cd_x[k] + g_z_0_y_xxxxyy[k];

            g_z_0_xy_xxxyz[k] = -g_z_0_y_xxxyz[k] * cd_x[k] + g_z_0_y_xxxxyz[k];

            g_z_0_xy_xxxzz[k] = -g_z_0_y_xxxzz[k] * cd_x[k] + g_z_0_y_xxxxzz[k];

            g_z_0_xy_xxyyy[k] = -g_z_0_y_xxyyy[k] * cd_x[k] + g_z_0_y_xxxyyy[k];

            g_z_0_xy_xxyyz[k] = -g_z_0_y_xxyyz[k] * cd_x[k] + g_z_0_y_xxxyyz[k];

            g_z_0_xy_xxyzz[k] = -g_z_0_y_xxyzz[k] * cd_x[k] + g_z_0_y_xxxyzz[k];

            g_z_0_xy_xxzzz[k] = -g_z_0_y_xxzzz[k] * cd_x[k] + g_z_0_y_xxxzzz[k];

            g_z_0_xy_xyyyy[k] = -g_z_0_y_xyyyy[k] * cd_x[k] + g_z_0_y_xxyyyy[k];

            g_z_0_xy_xyyyz[k] = -g_z_0_y_xyyyz[k] * cd_x[k] + g_z_0_y_xxyyyz[k];

            g_z_0_xy_xyyzz[k] = -g_z_0_y_xyyzz[k] * cd_x[k] + g_z_0_y_xxyyzz[k];

            g_z_0_xy_xyzzz[k] = -g_z_0_y_xyzzz[k] * cd_x[k] + g_z_0_y_xxyzzz[k];

            g_z_0_xy_xzzzz[k] = -g_z_0_y_xzzzz[k] * cd_x[k] + g_z_0_y_xxzzzz[k];

            g_z_0_xy_yyyyy[k] = -g_z_0_y_yyyyy[k] * cd_x[k] + g_z_0_y_xyyyyy[k];

            g_z_0_xy_yyyyz[k] = -g_z_0_y_yyyyz[k] * cd_x[k] + g_z_0_y_xyyyyz[k];

            g_z_0_xy_yyyzz[k] = -g_z_0_y_yyyzz[k] * cd_x[k] + g_z_0_y_xyyyzz[k];

            g_z_0_xy_yyzzz[k] = -g_z_0_y_yyzzz[k] * cd_x[k] + g_z_0_y_xyyzzz[k];

            g_z_0_xy_yzzzz[k] = -g_z_0_y_yzzzz[k] * cd_x[k] + g_z_0_y_xyzzzz[k];

            g_z_0_xy_zzzzz[k] = -g_z_0_y_zzzzz[k] * cd_x[k] + g_z_0_y_xzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_z_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 42);

        auto g_z_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 43);

        auto g_z_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 44);

        auto g_z_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 45);

        auto g_z_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 46);

        auto g_z_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 47);

        auto g_z_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 48);

        auto g_z_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 49);

        auto g_z_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 50);

        auto g_z_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 51);

        auto g_z_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 52);

        auto g_z_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 53);

        auto g_z_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 54);

        auto g_z_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 55);

        auto g_z_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 56);

        auto g_z_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 57);

        auto g_z_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 58);

        auto g_z_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 59);

        auto g_z_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 60);

        auto g_z_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 61);

        auto g_z_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 62);

        #pragma omp simd aligned(cd_x, g_z_0_xz_xxxxx, g_z_0_xz_xxxxy, g_z_0_xz_xxxxz, g_z_0_xz_xxxyy, g_z_0_xz_xxxyz, g_z_0_xz_xxxzz, g_z_0_xz_xxyyy, g_z_0_xz_xxyyz, g_z_0_xz_xxyzz, g_z_0_xz_xxzzz, g_z_0_xz_xyyyy, g_z_0_xz_xyyyz, g_z_0_xz_xyyzz, g_z_0_xz_xyzzz, g_z_0_xz_xzzzz, g_z_0_xz_yyyyy, g_z_0_xz_yyyyz, g_z_0_xz_yyyzz, g_z_0_xz_yyzzz, g_z_0_xz_yzzzz, g_z_0_xz_zzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxy, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxz, g_z_0_z_xxxxzz, g_z_0_z_xxxyy, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyy, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xz_xxxxx[k] = -g_z_0_z_xxxxx[k] * cd_x[k] + g_z_0_z_xxxxxx[k];

            g_z_0_xz_xxxxy[k] = -g_z_0_z_xxxxy[k] * cd_x[k] + g_z_0_z_xxxxxy[k];

            g_z_0_xz_xxxxz[k] = -g_z_0_z_xxxxz[k] * cd_x[k] + g_z_0_z_xxxxxz[k];

            g_z_0_xz_xxxyy[k] = -g_z_0_z_xxxyy[k] * cd_x[k] + g_z_0_z_xxxxyy[k];

            g_z_0_xz_xxxyz[k] = -g_z_0_z_xxxyz[k] * cd_x[k] + g_z_0_z_xxxxyz[k];

            g_z_0_xz_xxxzz[k] = -g_z_0_z_xxxzz[k] * cd_x[k] + g_z_0_z_xxxxzz[k];

            g_z_0_xz_xxyyy[k] = -g_z_0_z_xxyyy[k] * cd_x[k] + g_z_0_z_xxxyyy[k];

            g_z_0_xz_xxyyz[k] = -g_z_0_z_xxyyz[k] * cd_x[k] + g_z_0_z_xxxyyz[k];

            g_z_0_xz_xxyzz[k] = -g_z_0_z_xxyzz[k] * cd_x[k] + g_z_0_z_xxxyzz[k];

            g_z_0_xz_xxzzz[k] = -g_z_0_z_xxzzz[k] * cd_x[k] + g_z_0_z_xxxzzz[k];

            g_z_0_xz_xyyyy[k] = -g_z_0_z_xyyyy[k] * cd_x[k] + g_z_0_z_xxyyyy[k];

            g_z_0_xz_xyyyz[k] = -g_z_0_z_xyyyz[k] * cd_x[k] + g_z_0_z_xxyyyz[k];

            g_z_0_xz_xyyzz[k] = -g_z_0_z_xyyzz[k] * cd_x[k] + g_z_0_z_xxyyzz[k];

            g_z_0_xz_xyzzz[k] = -g_z_0_z_xyzzz[k] * cd_x[k] + g_z_0_z_xxyzzz[k];

            g_z_0_xz_xzzzz[k] = -g_z_0_z_xzzzz[k] * cd_x[k] + g_z_0_z_xxzzzz[k];

            g_z_0_xz_yyyyy[k] = -g_z_0_z_yyyyy[k] * cd_x[k] + g_z_0_z_xyyyyy[k];

            g_z_0_xz_yyyyz[k] = -g_z_0_z_yyyyz[k] * cd_x[k] + g_z_0_z_xyyyyz[k];

            g_z_0_xz_yyyzz[k] = -g_z_0_z_yyyzz[k] * cd_x[k] + g_z_0_z_xyyyzz[k];

            g_z_0_xz_yyzzz[k] = -g_z_0_z_yyzzz[k] * cd_x[k] + g_z_0_z_xyyzzz[k];

            g_z_0_xz_yzzzz[k] = -g_z_0_z_yzzzz[k] * cd_x[k] + g_z_0_z_xyzzzz[k];

            g_z_0_xz_zzzzz[k] = -g_z_0_z_zzzzz[k] * cd_x[k] + g_z_0_z_xzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_z_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 63);

        auto g_z_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 64);

        auto g_z_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 65);

        auto g_z_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 66);

        auto g_z_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 67);

        auto g_z_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 68);

        auto g_z_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 69);

        auto g_z_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 70);

        auto g_z_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 71);

        auto g_z_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 72);

        auto g_z_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 73);

        auto g_z_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 74);

        auto g_z_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 75);

        auto g_z_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 76);

        auto g_z_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 77);

        auto g_z_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 78);

        auto g_z_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 79);

        auto g_z_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 80);

        auto g_z_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 81);

        auto g_z_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 82);

        auto g_z_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_z_0_y_xxxxx, g_z_0_y_xxxxxy, g_z_0_y_xxxxy, g_z_0_y_xxxxyy, g_z_0_y_xxxxyz, g_z_0_y_xxxxz, g_z_0_y_xxxyy, g_z_0_y_xxxyyy, g_z_0_y_xxxyyz, g_z_0_y_xxxyz, g_z_0_y_xxxyzz, g_z_0_y_xxxzz, g_z_0_y_xxyyy, g_z_0_y_xxyyyy, g_z_0_y_xxyyyz, g_z_0_y_xxyyz, g_z_0_y_xxyyzz, g_z_0_y_xxyzz, g_z_0_y_xxyzzz, g_z_0_y_xxzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyyy, g_z_0_y_xyyyyz, g_z_0_y_xyyyz, g_z_0_y_xyyyzz, g_z_0_y_xyyzz, g_z_0_y_xyyzzz, g_z_0_y_xyzzz, g_z_0_y_xyzzzz, g_z_0_y_xzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyyy, g_z_0_y_yyyyyz, g_z_0_y_yyyyz, g_z_0_y_yyyyzz, g_z_0_y_yyyzz, g_z_0_y_yyyzzz, g_z_0_y_yyzzz, g_z_0_y_yyzzzz, g_z_0_y_yzzzz, g_z_0_y_yzzzzz, g_z_0_y_zzzzz, g_z_0_yy_xxxxx, g_z_0_yy_xxxxy, g_z_0_yy_xxxxz, g_z_0_yy_xxxyy, g_z_0_yy_xxxyz, g_z_0_yy_xxxzz, g_z_0_yy_xxyyy, g_z_0_yy_xxyyz, g_z_0_yy_xxyzz, g_z_0_yy_xxzzz, g_z_0_yy_xyyyy, g_z_0_yy_xyyyz, g_z_0_yy_xyyzz, g_z_0_yy_xyzzz, g_z_0_yy_xzzzz, g_z_0_yy_yyyyy, g_z_0_yy_yyyyz, g_z_0_yy_yyyzz, g_z_0_yy_yyzzz, g_z_0_yy_yzzzz, g_z_0_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yy_xxxxx[k] = -g_z_0_y_xxxxx[k] * cd_y[k] + g_z_0_y_xxxxxy[k];

            g_z_0_yy_xxxxy[k] = -g_z_0_y_xxxxy[k] * cd_y[k] + g_z_0_y_xxxxyy[k];

            g_z_0_yy_xxxxz[k] = -g_z_0_y_xxxxz[k] * cd_y[k] + g_z_0_y_xxxxyz[k];

            g_z_0_yy_xxxyy[k] = -g_z_0_y_xxxyy[k] * cd_y[k] + g_z_0_y_xxxyyy[k];

            g_z_0_yy_xxxyz[k] = -g_z_0_y_xxxyz[k] * cd_y[k] + g_z_0_y_xxxyyz[k];

            g_z_0_yy_xxxzz[k] = -g_z_0_y_xxxzz[k] * cd_y[k] + g_z_0_y_xxxyzz[k];

            g_z_0_yy_xxyyy[k] = -g_z_0_y_xxyyy[k] * cd_y[k] + g_z_0_y_xxyyyy[k];

            g_z_0_yy_xxyyz[k] = -g_z_0_y_xxyyz[k] * cd_y[k] + g_z_0_y_xxyyyz[k];

            g_z_0_yy_xxyzz[k] = -g_z_0_y_xxyzz[k] * cd_y[k] + g_z_0_y_xxyyzz[k];

            g_z_0_yy_xxzzz[k] = -g_z_0_y_xxzzz[k] * cd_y[k] + g_z_0_y_xxyzzz[k];

            g_z_0_yy_xyyyy[k] = -g_z_0_y_xyyyy[k] * cd_y[k] + g_z_0_y_xyyyyy[k];

            g_z_0_yy_xyyyz[k] = -g_z_0_y_xyyyz[k] * cd_y[k] + g_z_0_y_xyyyyz[k];

            g_z_0_yy_xyyzz[k] = -g_z_0_y_xyyzz[k] * cd_y[k] + g_z_0_y_xyyyzz[k];

            g_z_0_yy_xyzzz[k] = -g_z_0_y_xyzzz[k] * cd_y[k] + g_z_0_y_xyyzzz[k];

            g_z_0_yy_xzzzz[k] = -g_z_0_y_xzzzz[k] * cd_y[k] + g_z_0_y_xyzzzz[k];

            g_z_0_yy_yyyyy[k] = -g_z_0_y_yyyyy[k] * cd_y[k] + g_z_0_y_yyyyyy[k];

            g_z_0_yy_yyyyz[k] = -g_z_0_y_yyyyz[k] * cd_y[k] + g_z_0_y_yyyyyz[k];

            g_z_0_yy_yyyzz[k] = -g_z_0_y_yyyzz[k] * cd_y[k] + g_z_0_y_yyyyzz[k];

            g_z_0_yy_yyzzz[k] = -g_z_0_y_yyzzz[k] * cd_y[k] + g_z_0_y_yyyzzz[k];

            g_z_0_yy_yzzzz[k] = -g_z_0_y_yzzzz[k] * cd_y[k] + g_z_0_y_yyzzzz[k];

            g_z_0_yy_zzzzz[k] = -g_z_0_y_zzzzz[k] * cd_y[k] + g_z_0_y_yzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_z_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 84);

        auto g_z_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 85);

        auto g_z_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 86);

        auto g_z_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 87);

        auto g_z_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 88);

        auto g_z_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 89);

        auto g_z_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 90);

        auto g_z_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 91);

        auto g_z_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 92);

        auto g_z_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 93);

        auto g_z_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 94);

        auto g_z_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 95);

        auto g_z_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 96);

        auto g_z_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 97);

        auto g_z_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 98);

        auto g_z_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 99);

        auto g_z_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 100);

        auto g_z_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 101);

        auto g_z_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 102);

        auto g_z_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 103);

        auto g_z_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_z_0_yz_xxxxx, g_z_0_yz_xxxxy, g_z_0_yz_xxxxz, g_z_0_yz_xxxyy, g_z_0_yz_xxxyz, g_z_0_yz_xxxzz, g_z_0_yz_xxyyy, g_z_0_yz_xxyyz, g_z_0_yz_xxyzz, g_z_0_yz_xxzzz, g_z_0_yz_xyyyy, g_z_0_yz_xyyyz, g_z_0_yz_xyyzz, g_z_0_yz_xyzzz, g_z_0_yz_xzzzz, g_z_0_yz_yyyyy, g_z_0_yz_yyyyz, g_z_0_yz_yyyzz, g_z_0_yz_yyzzz, g_z_0_yz_yzzzz, g_z_0_yz_zzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxy, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yz_xxxxx[k] = -g_z_0_z_xxxxx[k] * cd_y[k] + g_z_0_z_xxxxxy[k];

            g_z_0_yz_xxxxy[k] = -g_z_0_z_xxxxy[k] * cd_y[k] + g_z_0_z_xxxxyy[k];

            g_z_0_yz_xxxxz[k] = -g_z_0_z_xxxxz[k] * cd_y[k] + g_z_0_z_xxxxyz[k];

            g_z_0_yz_xxxyy[k] = -g_z_0_z_xxxyy[k] * cd_y[k] + g_z_0_z_xxxyyy[k];

            g_z_0_yz_xxxyz[k] = -g_z_0_z_xxxyz[k] * cd_y[k] + g_z_0_z_xxxyyz[k];

            g_z_0_yz_xxxzz[k] = -g_z_0_z_xxxzz[k] * cd_y[k] + g_z_0_z_xxxyzz[k];

            g_z_0_yz_xxyyy[k] = -g_z_0_z_xxyyy[k] * cd_y[k] + g_z_0_z_xxyyyy[k];

            g_z_0_yz_xxyyz[k] = -g_z_0_z_xxyyz[k] * cd_y[k] + g_z_0_z_xxyyyz[k];

            g_z_0_yz_xxyzz[k] = -g_z_0_z_xxyzz[k] * cd_y[k] + g_z_0_z_xxyyzz[k];

            g_z_0_yz_xxzzz[k] = -g_z_0_z_xxzzz[k] * cd_y[k] + g_z_0_z_xxyzzz[k];

            g_z_0_yz_xyyyy[k] = -g_z_0_z_xyyyy[k] * cd_y[k] + g_z_0_z_xyyyyy[k];

            g_z_0_yz_xyyyz[k] = -g_z_0_z_xyyyz[k] * cd_y[k] + g_z_0_z_xyyyyz[k];

            g_z_0_yz_xyyzz[k] = -g_z_0_z_xyyzz[k] * cd_y[k] + g_z_0_z_xyyyzz[k];

            g_z_0_yz_xyzzz[k] = -g_z_0_z_xyzzz[k] * cd_y[k] + g_z_0_z_xyyzzz[k];

            g_z_0_yz_xzzzz[k] = -g_z_0_z_xzzzz[k] * cd_y[k] + g_z_0_z_xyzzzz[k];

            g_z_0_yz_yyyyy[k] = -g_z_0_z_yyyyy[k] * cd_y[k] + g_z_0_z_yyyyyy[k];

            g_z_0_yz_yyyyz[k] = -g_z_0_z_yyyyz[k] * cd_y[k] + g_z_0_z_yyyyyz[k];

            g_z_0_yz_yyyzz[k] = -g_z_0_z_yyyzz[k] * cd_y[k] + g_z_0_z_yyyyzz[k];

            g_z_0_yz_yyzzz[k] = -g_z_0_z_yyzzz[k] * cd_y[k] + g_z_0_z_yyyzzz[k];

            g_z_0_yz_yzzzz[k] = -g_z_0_z_yzzzz[k] * cd_y[k] + g_z_0_z_yyzzzz[k];

            g_z_0_yz_zzzzz[k] = -g_z_0_z_zzzzz[k] * cd_y[k] + g_z_0_z_yzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_z_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps  + 105);

        auto g_z_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 106);

        auto g_z_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 107);

        auto g_z_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 108);

        auto g_z_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 109);

        auto g_z_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 110);

        auto g_z_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 111);

        auto g_z_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 112);

        auto g_z_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 113);

        auto g_z_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 114);

        auto g_z_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 115);

        auto g_z_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 116);

        auto g_z_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 117);

        auto g_z_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 118);

        auto g_z_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 119);

        auto g_z_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps  + 120);

        auto g_z_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 121);

        auto g_z_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 122);

        auto g_z_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 123);

        auto g_z_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 124);

        auto g_z_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps  + 125);

        #pragma omp simd aligned(cd_z, g_z_0_z_xxxxx, g_z_0_z_xxxxxz, g_z_0_z_xxxxy, g_z_0_z_xxxxyz, g_z_0_z_xxxxz, g_z_0_z_xxxxzz, g_z_0_z_xxxyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzz, g_z_0_z_zzzzzz, g_z_0_zz_xxxxx, g_z_0_zz_xxxxy, g_z_0_zz_xxxxz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyz, g_z_0_zz_xxxzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyz, g_z_0_zz_xxyzz, g_z_0_zz_xxzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyz, g_z_0_zz_xyyzz, g_z_0_zz_xyzzz, g_z_0_zz_xzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyz, g_z_0_zz_yyyzz, g_z_0_zz_yyzzz, g_z_0_zz_yzzzz, g_z_0_zz_zzzzz, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxyy, g_z_xxxyz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyz, g_z_xxyzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyz, g_z_xyyzz, g_z_xyzzz, g_z_xzzzz, g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz, g_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zz_xxxxx[k] = -g_z_xxxxx[k] - g_z_0_z_xxxxx[k] * cd_z[k] + g_z_0_z_xxxxxz[k];

            g_z_0_zz_xxxxy[k] = -g_z_xxxxy[k] - g_z_0_z_xxxxy[k] * cd_z[k] + g_z_0_z_xxxxyz[k];

            g_z_0_zz_xxxxz[k] = -g_z_xxxxz[k] - g_z_0_z_xxxxz[k] * cd_z[k] + g_z_0_z_xxxxzz[k];

            g_z_0_zz_xxxyy[k] = -g_z_xxxyy[k] - g_z_0_z_xxxyy[k] * cd_z[k] + g_z_0_z_xxxyyz[k];

            g_z_0_zz_xxxyz[k] = -g_z_xxxyz[k] - g_z_0_z_xxxyz[k] * cd_z[k] + g_z_0_z_xxxyzz[k];

            g_z_0_zz_xxxzz[k] = -g_z_xxxzz[k] - g_z_0_z_xxxzz[k] * cd_z[k] + g_z_0_z_xxxzzz[k];

            g_z_0_zz_xxyyy[k] = -g_z_xxyyy[k] - g_z_0_z_xxyyy[k] * cd_z[k] + g_z_0_z_xxyyyz[k];

            g_z_0_zz_xxyyz[k] = -g_z_xxyyz[k] - g_z_0_z_xxyyz[k] * cd_z[k] + g_z_0_z_xxyyzz[k];

            g_z_0_zz_xxyzz[k] = -g_z_xxyzz[k] - g_z_0_z_xxyzz[k] * cd_z[k] + g_z_0_z_xxyzzz[k];

            g_z_0_zz_xxzzz[k] = -g_z_xxzzz[k] - g_z_0_z_xxzzz[k] * cd_z[k] + g_z_0_z_xxzzzz[k];

            g_z_0_zz_xyyyy[k] = -g_z_xyyyy[k] - g_z_0_z_xyyyy[k] * cd_z[k] + g_z_0_z_xyyyyz[k];

            g_z_0_zz_xyyyz[k] = -g_z_xyyyz[k] - g_z_0_z_xyyyz[k] * cd_z[k] + g_z_0_z_xyyyzz[k];

            g_z_0_zz_xyyzz[k] = -g_z_xyyzz[k] - g_z_0_z_xyyzz[k] * cd_z[k] + g_z_0_z_xyyzzz[k];

            g_z_0_zz_xyzzz[k] = -g_z_xyzzz[k] - g_z_0_z_xyzzz[k] * cd_z[k] + g_z_0_z_xyzzzz[k];

            g_z_0_zz_xzzzz[k] = -g_z_xzzzz[k] - g_z_0_z_xzzzz[k] * cd_z[k] + g_z_0_z_xzzzzz[k];

            g_z_0_zz_yyyyy[k] = -g_z_yyyyy[k] - g_z_0_z_yyyyy[k] * cd_z[k] + g_z_0_z_yyyyyz[k];

            g_z_0_zz_yyyyz[k] = -g_z_yyyyz[k] - g_z_0_z_yyyyz[k] * cd_z[k] + g_z_0_z_yyyyzz[k];

            g_z_0_zz_yyyzz[k] = -g_z_yyyzz[k] - g_z_0_z_yyyzz[k] * cd_z[k] + g_z_0_z_yyyzzz[k];

            g_z_0_zz_yyzzz[k] = -g_z_yyzzz[k] - g_z_0_z_yyzzz[k] * cd_z[k] + g_z_0_z_yyzzzz[k];

            g_z_0_zz_yzzzz[k] = -g_z_yzzzz[k] - g_z_0_z_yzzzz[k] * cd_z[k] + g_z_0_z_yzzzzz[k];

            g_z_0_zz_zzzzz[k] = -g_z_zzzzz[k] - g_z_0_z_zzzzz[k] * cd_z[k] + g_z_0_z_zzzzzz[k];
        }
    }
}

} // t3ceri namespace


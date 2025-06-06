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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXPK.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xpk(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xpk,
                                        const size_t idx_xsk,
                                        const size_t idx_geom_10_xsk,
                                        const size_t idx_geom_10_xsl,
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
        /// Set up components of auxilary buffer : SSK

        const auto sk_off = idx_xsk + i * 36;

        auto g_0_xxxxxxx = cbuffer.data(sk_off + 0);

        auto g_0_xxxxxxy = cbuffer.data(sk_off + 1);

        auto g_0_xxxxxxz = cbuffer.data(sk_off + 2);

        auto g_0_xxxxxyy = cbuffer.data(sk_off + 3);

        auto g_0_xxxxxyz = cbuffer.data(sk_off + 4);

        auto g_0_xxxxxzz = cbuffer.data(sk_off + 5);

        auto g_0_xxxxyyy = cbuffer.data(sk_off + 6);

        auto g_0_xxxxyyz = cbuffer.data(sk_off + 7);

        auto g_0_xxxxyzz = cbuffer.data(sk_off + 8);

        auto g_0_xxxxzzz = cbuffer.data(sk_off + 9);

        auto g_0_xxxyyyy = cbuffer.data(sk_off + 10);

        auto g_0_xxxyyyz = cbuffer.data(sk_off + 11);

        auto g_0_xxxyyzz = cbuffer.data(sk_off + 12);

        auto g_0_xxxyzzz = cbuffer.data(sk_off + 13);

        auto g_0_xxxzzzz = cbuffer.data(sk_off + 14);

        auto g_0_xxyyyyy = cbuffer.data(sk_off + 15);

        auto g_0_xxyyyyz = cbuffer.data(sk_off + 16);

        auto g_0_xxyyyzz = cbuffer.data(sk_off + 17);

        auto g_0_xxyyzzz = cbuffer.data(sk_off + 18);

        auto g_0_xxyzzzz = cbuffer.data(sk_off + 19);

        auto g_0_xxzzzzz = cbuffer.data(sk_off + 20);

        auto g_0_xyyyyyy = cbuffer.data(sk_off + 21);

        auto g_0_xyyyyyz = cbuffer.data(sk_off + 22);

        auto g_0_xyyyyzz = cbuffer.data(sk_off + 23);

        auto g_0_xyyyzzz = cbuffer.data(sk_off + 24);

        auto g_0_xyyzzzz = cbuffer.data(sk_off + 25);

        auto g_0_xyzzzzz = cbuffer.data(sk_off + 26);

        auto g_0_xzzzzzz = cbuffer.data(sk_off + 27);

        auto g_0_yyyyyyy = cbuffer.data(sk_off + 28);

        auto g_0_yyyyyyz = cbuffer.data(sk_off + 29);

        auto g_0_yyyyyzz = cbuffer.data(sk_off + 30);

        auto g_0_yyyyzzz = cbuffer.data(sk_off + 31);

        auto g_0_yyyzzzz = cbuffer.data(sk_off + 32);

        auto g_0_yyzzzzz = cbuffer.data(sk_off + 33);

        auto g_0_yzzzzzz = cbuffer.data(sk_off + 34);

        auto g_0_zzzzzzz = cbuffer.data(sk_off + 35);

        /// Set up components of auxilary buffer : SSK

        const auto sk_geom_10_off = idx_geom_10_xsk + i * 36;

        auto g_x_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 9);

        auto g_x_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 10);

        auto g_x_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 11);

        auto g_x_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 12);

        auto g_x_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 13);

        auto g_x_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 14);

        auto g_x_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 15);

        auto g_x_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 16);

        auto g_x_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 17);

        auto g_x_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 18);

        auto g_x_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 19);

        auto g_x_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 20);

        auto g_x_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 21);

        auto g_x_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 22);

        auto g_x_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 23);

        auto g_x_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 24);

        auto g_x_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 25);

        auto g_x_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 26);

        auto g_x_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 27);

        auto g_x_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps + 28);

        auto g_x_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps + 29);

        auto g_x_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 30);

        auto g_x_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 31);

        auto g_x_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 32);

        auto g_x_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 33);

        auto g_x_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 34);

        auto g_x_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps + 35);

        auto g_y_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 36 * acomps + 0);

        auto g_y_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 36 * acomps + 1);

        auto g_y_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 36 * acomps + 2);

        auto g_y_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 3);

        auto g_y_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 4);

        auto g_y_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 5);

        auto g_y_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 6);

        auto g_y_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 7);

        auto g_y_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 8);

        auto g_y_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 9);

        auto g_y_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 10);

        auto g_y_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 11);

        auto g_y_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 12);

        auto g_y_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 13);

        auto g_y_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 14);

        auto g_y_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 15);

        auto g_y_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 16);

        auto g_y_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 17);

        auto g_y_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 18);

        auto g_y_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 19);

        auto g_y_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 20);

        auto g_y_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 21);

        auto g_y_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 22);

        auto g_y_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 23);

        auto g_y_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 24);

        auto g_y_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 25);

        auto g_y_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 26);

        auto g_y_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 27);

        auto g_y_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps + 28);

        auto g_y_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps + 29);

        auto g_y_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 30);

        auto g_y_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 31);

        auto g_y_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 32);

        auto g_y_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 33);

        auto g_y_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 34);

        auto g_y_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps + 35);

        auto g_z_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 72 * acomps + 0);

        auto g_z_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 72 * acomps + 1);

        auto g_z_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 72 * acomps + 2);

        auto g_z_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 3);

        auto g_z_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 4);

        auto g_z_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 5);

        auto g_z_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 6);

        auto g_z_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 7);

        auto g_z_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 8);

        auto g_z_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 9);

        auto g_z_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 10);

        auto g_z_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 11);

        auto g_z_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 12);

        auto g_z_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 13);

        auto g_z_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 14);

        auto g_z_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 15);

        auto g_z_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 16);

        auto g_z_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 17);

        auto g_z_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 18);

        auto g_z_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 19);

        auto g_z_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 20);

        auto g_z_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 21);

        auto g_z_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 22);

        auto g_z_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 23);

        auto g_z_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 24);

        auto g_z_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 25);

        auto g_z_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 26);

        auto g_z_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 27);

        auto g_z_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps + 28);

        auto g_z_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps + 29);

        auto g_z_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 30);

        auto g_z_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 31);

        auto g_z_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 32);

        auto g_z_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 33);

        auto g_z_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 34);

        auto g_z_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps + 35);

        /// Set up components of auxilary buffer : SSL

        const auto sl_geom_10_off = idx_geom_10_xsl + i * 45;

        auto g_x_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 9);

        auto g_x_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 10);

        auto g_x_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 11);

        auto g_x_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 12);

        auto g_x_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 13);

        auto g_x_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 14);

        auto g_x_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 15);

        auto g_x_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 16);

        auto g_x_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 17);

        auto g_x_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 18);

        auto g_x_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 19);

        auto g_x_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 20);

        auto g_x_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 21);

        auto g_x_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 22);

        auto g_x_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 23);

        auto g_x_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 24);

        auto g_x_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 25);

        auto g_x_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 26);

        auto g_x_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 27);

        auto g_x_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 28);

        auto g_x_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 29);

        auto g_x_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 30);

        auto g_x_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 31);

        auto g_x_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 32);

        auto g_x_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 33);

        auto g_x_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 34);

        auto g_x_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 35);

        auto g_x_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps + 36);

        auto g_x_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps + 37);

        auto g_x_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 38);

        auto g_x_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 39);

        auto g_x_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 40);

        auto g_x_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 41);

        auto g_x_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 42);

        auto g_x_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 43);

        auto g_x_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps + 44);

        auto g_y_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 45 * acomps + 0);

        auto g_y_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 45 * acomps + 1);

        auto g_y_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 45 * acomps + 2);

        auto g_y_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 3);

        auto g_y_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 4);

        auto g_y_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 5);

        auto g_y_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 6);

        auto g_y_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 7);

        auto g_y_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 8);

        auto g_y_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 9);

        auto g_y_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 10);

        auto g_y_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 11);

        auto g_y_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 12);

        auto g_y_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 13);

        auto g_y_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 14);

        auto g_y_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 15);

        auto g_y_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 16);

        auto g_y_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 17);

        auto g_y_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 18);

        auto g_y_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 19);

        auto g_y_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 20);

        auto g_y_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 21);

        auto g_y_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 22);

        auto g_y_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 23);

        auto g_y_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 24);

        auto g_y_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 25);

        auto g_y_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 26);

        auto g_y_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 27);

        auto g_y_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 28);

        auto g_y_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 29);

        auto g_y_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 30);

        auto g_y_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 31);

        auto g_y_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 32);

        auto g_y_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 33);

        auto g_y_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 34);

        auto g_y_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 35);

        auto g_y_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps + 36);

        auto g_y_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps + 37);

        auto g_y_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 38);

        auto g_y_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 39);

        auto g_y_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 40);

        auto g_y_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 41);

        auto g_y_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 42);

        auto g_y_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 43);

        auto g_y_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps + 44);

        auto g_z_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 90 * acomps + 0);

        auto g_z_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 90 * acomps + 1);

        auto g_z_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 90 * acomps + 2);

        auto g_z_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 3);

        auto g_z_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 4);

        auto g_z_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 5);

        auto g_z_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 6);

        auto g_z_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 7);

        auto g_z_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 8);

        auto g_z_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 9);

        auto g_z_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 10);

        auto g_z_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 11);

        auto g_z_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 12);

        auto g_z_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 13);

        auto g_z_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 14);

        auto g_z_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 15);

        auto g_z_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 16);

        auto g_z_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 17);

        auto g_z_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 18);

        auto g_z_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 19);

        auto g_z_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 20);

        auto g_z_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 21);

        auto g_z_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 22);

        auto g_z_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 23);

        auto g_z_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 24);

        auto g_z_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 25);

        auto g_z_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 26);

        auto g_z_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 27);

        auto g_z_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 28);

        auto g_z_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 29);

        auto g_z_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 30);

        auto g_z_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 31);

        auto g_z_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 32);

        auto g_z_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 33);

        auto g_z_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 34);

        auto g_z_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 35);

        auto g_z_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps + 36);

        auto g_z_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps + 37);

        auto g_z_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 38);

        auto g_z_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 39);

        auto g_z_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 40);

        auto g_z_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 41);

        auto g_z_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 42);

        auto g_z_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 43);

        auto g_z_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps + 44);

        /// set up bra offset for contr_buffer_xxpk

        const auto pk_geom_10_off = idx_geom_10_xpk + i * 108;

        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_x_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_x_yyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_x_yyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_x_yyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_x_yyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_x_yyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_x_yyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_x_yzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_x_zzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz, g_x_0_0_xxxxxxx, g_x_0_0_xxxxxxxx, g_x_0_0_xxxxxxxy, g_x_0_0_xxxxxxxz, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxxyy, g_x_0_0_xxxxxxyz, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxxzz, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyyy, g_x_0_0_xxxxxyyz, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxyzz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxxzzz, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyyy, g_x_0_0_xxxxyyyz, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyyzz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxyzzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxxzzzz, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyyy, g_x_0_0_xxxyyyyz, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyyzz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyyzzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxyzzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxxzzzzz, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyyy, g_x_0_0_xxyyyyyz, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyyzz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyyzzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyyzzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxyzzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xxzzzzzz, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyyy, g_x_0_0_xyyyyyyz, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyyzz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyyzzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyyzzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyyzzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xyzzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_xzzzzzzz, g_x_0_0_yyyyyyy, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_zzzzzzz, g_x_0_x_xxxxxxx, g_x_0_x_xxxxxxy, g_x_0_x_xxxxxxz, g_x_0_x_xxxxxyy, g_x_0_x_xxxxxyz, g_x_0_x_xxxxxzz, g_x_0_x_xxxxyyy, g_x_0_x_xxxxyyz, g_x_0_x_xxxxyzz, g_x_0_x_xxxxzzz, g_x_0_x_xxxyyyy, g_x_0_x_xxxyyyz, g_x_0_x_xxxyyzz, g_x_0_x_xxxyzzz, g_x_0_x_xxxzzzz, g_x_0_x_xxyyyyy, g_x_0_x_xxyyyyz, g_x_0_x_xxyyyzz, g_x_0_x_xxyyzzz, g_x_0_x_xxyzzzz, g_x_0_x_xxzzzzz, g_x_0_x_xyyyyyy, g_x_0_x_xyyyyyz, g_x_0_x_xyyyyzz, g_x_0_x_xyyyzzz, g_x_0_x_xyyzzzz, g_x_0_x_xyzzzzz, g_x_0_x_xzzzzzz, g_x_0_x_yyyyyyy, g_x_0_x_yyyyyyz, g_x_0_x_yyyyyzz, g_x_0_x_yyyyzzz, g_x_0_x_yyyzzzz, g_x_0_x_yyzzzzz, g_x_0_x_yzzzzzz, g_x_0_x_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_x_xxxxxxx[k] = -g_0_xxxxxxx[k] - g_x_0_0_xxxxxxx[k] * cd_x[k] + g_x_0_0_xxxxxxxx[k];

            g_x_0_x_xxxxxxy[k] = -g_0_xxxxxxy[k] - g_x_0_0_xxxxxxy[k] * cd_x[k] + g_x_0_0_xxxxxxxy[k];

            g_x_0_x_xxxxxxz[k] = -g_0_xxxxxxz[k] - g_x_0_0_xxxxxxz[k] * cd_x[k] + g_x_0_0_xxxxxxxz[k];

            g_x_0_x_xxxxxyy[k] = -g_0_xxxxxyy[k] - g_x_0_0_xxxxxyy[k] * cd_x[k] + g_x_0_0_xxxxxxyy[k];

            g_x_0_x_xxxxxyz[k] = -g_0_xxxxxyz[k] - g_x_0_0_xxxxxyz[k] * cd_x[k] + g_x_0_0_xxxxxxyz[k];

            g_x_0_x_xxxxxzz[k] = -g_0_xxxxxzz[k] - g_x_0_0_xxxxxzz[k] * cd_x[k] + g_x_0_0_xxxxxxzz[k];

            g_x_0_x_xxxxyyy[k] = -g_0_xxxxyyy[k] - g_x_0_0_xxxxyyy[k] * cd_x[k] + g_x_0_0_xxxxxyyy[k];

            g_x_0_x_xxxxyyz[k] = -g_0_xxxxyyz[k] - g_x_0_0_xxxxyyz[k] * cd_x[k] + g_x_0_0_xxxxxyyz[k];

            g_x_0_x_xxxxyzz[k] = -g_0_xxxxyzz[k] - g_x_0_0_xxxxyzz[k] * cd_x[k] + g_x_0_0_xxxxxyzz[k];

            g_x_0_x_xxxxzzz[k] = -g_0_xxxxzzz[k] - g_x_0_0_xxxxzzz[k] * cd_x[k] + g_x_0_0_xxxxxzzz[k];

            g_x_0_x_xxxyyyy[k] = -g_0_xxxyyyy[k] - g_x_0_0_xxxyyyy[k] * cd_x[k] + g_x_0_0_xxxxyyyy[k];

            g_x_0_x_xxxyyyz[k] = -g_0_xxxyyyz[k] - g_x_0_0_xxxyyyz[k] * cd_x[k] + g_x_0_0_xxxxyyyz[k];

            g_x_0_x_xxxyyzz[k] = -g_0_xxxyyzz[k] - g_x_0_0_xxxyyzz[k] * cd_x[k] + g_x_0_0_xxxxyyzz[k];

            g_x_0_x_xxxyzzz[k] = -g_0_xxxyzzz[k] - g_x_0_0_xxxyzzz[k] * cd_x[k] + g_x_0_0_xxxxyzzz[k];

            g_x_0_x_xxxzzzz[k] = -g_0_xxxzzzz[k] - g_x_0_0_xxxzzzz[k] * cd_x[k] + g_x_0_0_xxxxzzzz[k];

            g_x_0_x_xxyyyyy[k] = -g_0_xxyyyyy[k] - g_x_0_0_xxyyyyy[k] * cd_x[k] + g_x_0_0_xxxyyyyy[k];

            g_x_0_x_xxyyyyz[k] = -g_0_xxyyyyz[k] - g_x_0_0_xxyyyyz[k] * cd_x[k] + g_x_0_0_xxxyyyyz[k];

            g_x_0_x_xxyyyzz[k] = -g_0_xxyyyzz[k] - g_x_0_0_xxyyyzz[k] * cd_x[k] + g_x_0_0_xxxyyyzz[k];

            g_x_0_x_xxyyzzz[k] = -g_0_xxyyzzz[k] - g_x_0_0_xxyyzzz[k] * cd_x[k] + g_x_0_0_xxxyyzzz[k];

            g_x_0_x_xxyzzzz[k] = -g_0_xxyzzzz[k] - g_x_0_0_xxyzzzz[k] * cd_x[k] + g_x_0_0_xxxyzzzz[k];

            g_x_0_x_xxzzzzz[k] = -g_0_xxzzzzz[k] - g_x_0_0_xxzzzzz[k] * cd_x[k] + g_x_0_0_xxxzzzzz[k];

            g_x_0_x_xyyyyyy[k] = -g_0_xyyyyyy[k] - g_x_0_0_xyyyyyy[k] * cd_x[k] + g_x_0_0_xxyyyyyy[k];

            g_x_0_x_xyyyyyz[k] = -g_0_xyyyyyz[k] - g_x_0_0_xyyyyyz[k] * cd_x[k] + g_x_0_0_xxyyyyyz[k];

            g_x_0_x_xyyyyzz[k] = -g_0_xyyyyzz[k] - g_x_0_0_xyyyyzz[k] * cd_x[k] + g_x_0_0_xxyyyyzz[k];

            g_x_0_x_xyyyzzz[k] = -g_0_xyyyzzz[k] - g_x_0_0_xyyyzzz[k] * cd_x[k] + g_x_0_0_xxyyyzzz[k];

            g_x_0_x_xyyzzzz[k] = -g_0_xyyzzzz[k] - g_x_0_0_xyyzzzz[k] * cd_x[k] + g_x_0_0_xxyyzzzz[k];

            g_x_0_x_xyzzzzz[k] = -g_0_xyzzzzz[k] - g_x_0_0_xyzzzzz[k] * cd_x[k] + g_x_0_0_xxyzzzzz[k];

            g_x_0_x_xzzzzzz[k] = -g_0_xzzzzzz[k] - g_x_0_0_xzzzzzz[k] * cd_x[k] + g_x_0_0_xxzzzzzz[k];

            g_x_0_x_yyyyyyy[k] = -g_0_yyyyyyy[k] - g_x_0_0_yyyyyyy[k] * cd_x[k] + g_x_0_0_xyyyyyyy[k];

            g_x_0_x_yyyyyyz[k] = -g_0_yyyyyyz[k] - g_x_0_0_yyyyyyz[k] * cd_x[k] + g_x_0_0_xyyyyyyz[k];

            g_x_0_x_yyyyyzz[k] = -g_0_yyyyyzz[k] - g_x_0_0_yyyyyzz[k] * cd_x[k] + g_x_0_0_xyyyyyzz[k];

            g_x_0_x_yyyyzzz[k] = -g_0_yyyyzzz[k] - g_x_0_0_yyyyzzz[k] * cd_x[k] + g_x_0_0_xyyyyzzz[k];

            g_x_0_x_yyyzzzz[k] = -g_0_yyyzzzz[k] - g_x_0_0_yyyzzzz[k] * cd_x[k] + g_x_0_0_xyyyzzzz[k];

            g_x_0_x_yyzzzzz[k] = -g_0_yyzzzzz[k] - g_x_0_0_yyzzzzz[k] * cd_x[k] + g_x_0_0_xyyzzzzz[k];

            g_x_0_x_yzzzzzz[k] = -g_0_yzzzzzz[k] - g_x_0_0_yzzzzzz[k] * cd_x[k] + g_x_0_0_xyzzzzzz[k];

            g_x_0_x_zzzzzzz[k] = -g_0_zzzzzzz[k] - g_x_0_0_zzzzzzz[k] * cd_x[k] + g_x_0_0_xzzzzzzz[k];
        }

        /// Set up 36-72 components of targeted buffer : cbuffer.data(

        auto g_x_0_y_xxxxxxx = cbuffer.data(pk_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_y_xxxxxxz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_y_xxxxxzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_y_xxxxzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_y_xxxzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_y_xxzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 59);

        auto g_x_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_y_xzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_y_zzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 71);

        #pragma omp simd aligned(cd_y, g_x_0_0_xxxxxxx, g_x_0_0_xxxxxxxy, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxxyy, g_x_0_0_xxxxxxyz, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyyy, g_x_0_0_xxxxxyyz, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxyzz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyyy, g_x_0_0_xxxxyyyz, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyyzz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxyzzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyyy, g_x_0_0_xxxyyyyz, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyyzz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyyzzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxyzzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyyy, g_x_0_0_xxyyyyyz, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyyzz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyyzzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyyzzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxyzzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyyy, g_x_0_0_xyyyyyyz, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyyzz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyyzzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyyzzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyyzzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xyzzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_yyyyyyy, g_x_0_0_yyyyyyyy, g_x_0_0_yyyyyyyz, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyyzz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyyzzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyyzzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyyzzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yyzzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_yzzzzzzz, g_x_0_0_zzzzzzz, g_x_0_y_xxxxxxx, g_x_0_y_xxxxxxy, g_x_0_y_xxxxxxz, g_x_0_y_xxxxxyy, g_x_0_y_xxxxxyz, g_x_0_y_xxxxxzz, g_x_0_y_xxxxyyy, g_x_0_y_xxxxyyz, g_x_0_y_xxxxyzz, g_x_0_y_xxxxzzz, g_x_0_y_xxxyyyy, g_x_0_y_xxxyyyz, g_x_0_y_xxxyyzz, g_x_0_y_xxxyzzz, g_x_0_y_xxxzzzz, g_x_0_y_xxyyyyy, g_x_0_y_xxyyyyz, g_x_0_y_xxyyyzz, g_x_0_y_xxyyzzz, g_x_0_y_xxyzzzz, g_x_0_y_xxzzzzz, g_x_0_y_xyyyyyy, g_x_0_y_xyyyyyz, g_x_0_y_xyyyyzz, g_x_0_y_xyyyzzz, g_x_0_y_xyyzzzz, g_x_0_y_xyzzzzz, g_x_0_y_xzzzzzz, g_x_0_y_yyyyyyy, g_x_0_y_yyyyyyz, g_x_0_y_yyyyyzz, g_x_0_y_yyyyzzz, g_x_0_y_yyyzzzz, g_x_0_y_yyzzzzz, g_x_0_y_yzzzzzz, g_x_0_y_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_y_xxxxxxx[k] = -g_x_0_0_xxxxxxx[k] * cd_y[k] + g_x_0_0_xxxxxxxy[k];

            g_x_0_y_xxxxxxy[k] = -g_x_0_0_xxxxxxy[k] * cd_y[k] + g_x_0_0_xxxxxxyy[k];

            g_x_0_y_xxxxxxz[k] = -g_x_0_0_xxxxxxz[k] * cd_y[k] + g_x_0_0_xxxxxxyz[k];

            g_x_0_y_xxxxxyy[k] = -g_x_0_0_xxxxxyy[k] * cd_y[k] + g_x_0_0_xxxxxyyy[k];

            g_x_0_y_xxxxxyz[k] = -g_x_0_0_xxxxxyz[k] * cd_y[k] + g_x_0_0_xxxxxyyz[k];

            g_x_0_y_xxxxxzz[k] = -g_x_0_0_xxxxxzz[k] * cd_y[k] + g_x_0_0_xxxxxyzz[k];

            g_x_0_y_xxxxyyy[k] = -g_x_0_0_xxxxyyy[k] * cd_y[k] + g_x_0_0_xxxxyyyy[k];

            g_x_0_y_xxxxyyz[k] = -g_x_0_0_xxxxyyz[k] * cd_y[k] + g_x_0_0_xxxxyyyz[k];

            g_x_0_y_xxxxyzz[k] = -g_x_0_0_xxxxyzz[k] * cd_y[k] + g_x_0_0_xxxxyyzz[k];

            g_x_0_y_xxxxzzz[k] = -g_x_0_0_xxxxzzz[k] * cd_y[k] + g_x_0_0_xxxxyzzz[k];

            g_x_0_y_xxxyyyy[k] = -g_x_0_0_xxxyyyy[k] * cd_y[k] + g_x_0_0_xxxyyyyy[k];

            g_x_0_y_xxxyyyz[k] = -g_x_0_0_xxxyyyz[k] * cd_y[k] + g_x_0_0_xxxyyyyz[k];

            g_x_0_y_xxxyyzz[k] = -g_x_0_0_xxxyyzz[k] * cd_y[k] + g_x_0_0_xxxyyyzz[k];

            g_x_0_y_xxxyzzz[k] = -g_x_0_0_xxxyzzz[k] * cd_y[k] + g_x_0_0_xxxyyzzz[k];

            g_x_0_y_xxxzzzz[k] = -g_x_0_0_xxxzzzz[k] * cd_y[k] + g_x_0_0_xxxyzzzz[k];

            g_x_0_y_xxyyyyy[k] = -g_x_0_0_xxyyyyy[k] * cd_y[k] + g_x_0_0_xxyyyyyy[k];

            g_x_0_y_xxyyyyz[k] = -g_x_0_0_xxyyyyz[k] * cd_y[k] + g_x_0_0_xxyyyyyz[k];

            g_x_0_y_xxyyyzz[k] = -g_x_0_0_xxyyyzz[k] * cd_y[k] + g_x_0_0_xxyyyyzz[k];

            g_x_0_y_xxyyzzz[k] = -g_x_0_0_xxyyzzz[k] * cd_y[k] + g_x_0_0_xxyyyzzz[k];

            g_x_0_y_xxyzzzz[k] = -g_x_0_0_xxyzzzz[k] * cd_y[k] + g_x_0_0_xxyyzzzz[k];

            g_x_0_y_xxzzzzz[k] = -g_x_0_0_xxzzzzz[k] * cd_y[k] + g_x_0_0_xxyzzzzz[k];

            g_x_0_y_xyyyyyy[k] = -g_x_0_0_xyyyyyy[k] * cd_y[k] + g_x_0_0_xyyyyyyy[k];

            g_x_0_y_xyyyyyz[k] = -g_x_0_0_xyyyyyz[k] * cd_y[k] + g_x_0_0_xyyyyyyz[k];

            g_x_0_y_xyyyyzz[k] = -g_x_0_0_xyyyyzz[k] * cd_y[k] + g_x_0_0_xyyyyyzz[k];

            g_x_0_y_xyyyzzz[k] = -g_x_0_0_xyyyzzz[k] * cd_y[k] + g_x_0_0_xyyyyzzz[k];

            g_x_0_y_xyyzzzz[k] = -g_x_0_0_xyyzzzz[k] * cd_y[k] + g_x_0_0_xyyyzzzz[k];

            g_x_0_y_xyzzzzz[k] = -g_x_0_0_xyzzzzz[k] * cd_y[k] + g_x_0_0_xyyzzzzz[k];

            g_x_0_y_xzzzzzz[k] = -g_x_0_0_xzzzzzz[k] * cd_y[k] + g_x_0_0_xyzzzzzz[k];

            g_x_0_y_yyyyyyy[k] = -g_x_0_0_yyyyyyy[k] * cd_y[k] + g_x_0_0_yyyyyyyy[k];

            g_x_0_y_yyyyyyz[k] = -g_x_0_0_yyyyyyz[k] * cd_y[k] + g_x_0_0_yyyyyyyz[k];

            g_x_0_y_yyyyyzz[k] = -g_x_0_0_yyyyyzz[k] * cd_y[k] + g_x_0_0_yyyyyyzz[k];

            g_x_0_y_yyyyzzz[k] = -g_x_0_0_yyyyzzz[k] * cd_y[k] + g_x_0_0_yyyyyzzz[k];

            g_x_0_y_yyyzzzz[k] = -g_x_0_0_yyyzzzz[k] * cd_y[k] + g_x_0_0_yyyyzzzz[k];

            g_x_0_y_yyzzzzz[k] = -g_x_0_0_yyzzzzz[k] * cd_y[k] + g_x_0_0_yyyzzzzz[k];

            g_x_0_y_yzzzzzz[k] = -g_x_0_0_yzzzzzz[k] * cd_y[k] + g_x_0_0_yyzzzzzz[k];

            g_x_0_y_zzzzzzz[k] = -g_x_0_0_zzzzzzz[k] * cd_y[k] + g_x_0_0_yzzzzzzz[k];
        }

        /// Set up 72-108 components of targeted buffer : cbuffer.data(

        auto g_x_0_z_xxxxxxx = cbuffer.data(pk_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 74);

        auto g_x_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 83);

        auto g_x_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 89);

        auto g_x_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 99);

        auto g_x_0_z_yyyyyyy = cbuffer.data(pk_geom_10_off + 0 * acomps  + 100);

        auto g_x_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 101);

        auto g_x_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 102);

        auto g_x_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 103);

        auto g_x_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 104);

        auto g_x_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 105);

        auto g_x_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 106);

        auto g_x_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 0 * acomps  + 107);

        #pragma omp simd aligned(cd_z, g_x_0_0_xxxxxxx, g_x_0_0_xxxxxxxz, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxxyz, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxxzz, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyyz, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxyzz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxxzzz, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyyz, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyyzz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxyzzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxxzzzz, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyyz, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyyzz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyyzzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxyzzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxxzzzzz, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyyz, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyyzz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyyzzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyyzzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxyzzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xxzzzzzz, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyyz, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyyzz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyyzzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyyzzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyyzzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xyzzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_xzzzzzzz, g_x_0_0_yyyyyyy, g_x_0_0_yyyyyyyz, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyyzz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyyzzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyyzzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyyzzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yyzzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_yzzzzzzz, g_x_0_0_zzzzzzz, g_x_0_0_zzzzzzzz, g_x_0_z_xxxxxxx, g_x_0_z_xxxxxxy, g_x_0_z_xxxxxxz, g_x_0_z_xxxxxyy, g_x_0_z_xxxxxyz, g_x_0_z_xxxxxzz, g_x_0_z_xxxxyyy, g_x_0_z_xxxxyyz, g_x_0_z_xxxxyzz, g_x_0_z_xxxxzzz, g_x_0_z_xxxyyyy, g_x_0_z_xxxyyyz, g_x_0_z_xxxyyzz, g_x_0_z_xxxyzzz, g_x_0_z_xxxzzzz, g_x_0_z_xxyyyyy, g_x_0_z_xxyyyyz, g_x_0_z_xxyyyzz, g_x_0_z_xxyyzzz, g_x_0_z_xxyzzzz, g_x_0_z_xxzzzzz, g_x_0_z_xyyyyyy, g_x_0_z_xyyyyyz, g_x_0_z_xyyyyzz, g_x_0_z_xyyyzzz, g_x_0_z_xyyzzzz, g_x_0_z_xyzzzzz, g_x_0_z_xzzzzzz, g_x_0_z_yyyyyyy, g_x_0_z_yyyyyyz, g_x_0_z_yyyyyzz, g_x_0_z_yyyyzzz, g_x_0_z_yyyzzzz, g_x_0_z_yyzzzzz, g_x_0_z_yzzzzzz, g_x_0_z_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_z_xxxxxxx[k] = -g_x_0_0_xxxxxxx[k] * cd_z[k] + g_x_0_0_xxxxxxxz[k];

            g_x_0_z_xxxxxxy[k] = -g_x_0_0_xxxxxxy[k] * cd_z[k] + g_x_0_0_xxxxxxyz[k];

            g_x_0_z_xxxxxxz[k] = -g_x_0_0_xxxxxxz[k] * cd_z[k] + g_x_0_0_xxxxxxzz[k];

            g_x_0_z_xxxxxyy[k] = -g_x_0_0_xxxxxyy[k] * cd_z[k] + g_x_0_0_xxxxxyyz[k];

            g_x_0_z_xxxxxyz[k] = -g_x_0_0_xxxxxyz[k] * cd_z[k] + g_x_0_0_xxxxxyzz[k];

            g_x_0_z_xxxxxzz[k] = -g_x_0_0_xxxxxzz[k] * cd_z[k] + g_x_0_0_xxxxxzzz[k];

            g_x_0_z_xxxxyyy[k] = -g_x_0_0_xxxxyyy[k] * cd_z[k] + g_x_0_0_xxxxyyyz[k];

            g_x_0_z_xxxxyyz[k] = -g_x_0_0_xxxxyyz[k] * cd_z[k] + g_x_0_0_xxxxyyzz[k];

            g_x_0_z_xxxxyzz[k] = -g_x_0_0_xxxxyzz[k] * cd_z[k] + g_x_0_0_xxxxyzzz[k];

            g_x_0_z_xxxxzzz[k] = -g_x_0_0_xxxxzzz[k] * cd_z[k] + g_x_0_0_xxxxzzzz[k];

            g_x_0_z_xxxyyyy[k] = -g_x_0_0_xxxyyyy[k] * cd_z[k] + g_x_0_0_xxxyyyyz[k];

            g_x_0_z_xxxyyyz[k] = -g_x_0_0_xxxyyyz[k] * cd_z[k] + g_x_0_0_xxxyyyzz[k];

            g_x_0_z_xxxyyzz[k] = -g_x_0_0_xxxyyzz[k] * cd_z[k] + g_x_0_0_xxxyyzzz[k];

            g_x_0_z_xxxyzzz[k] = -g_x_0_0_xxxyzzz[k] * cd_z[k] + g_x_0_0_xxxyzzzz[k];

            g_x_0_z_xxxzzzz[k] = -g_x_0_0_xxxzzzz[k] * cd_z[k] + g_x_0_0_xxxzzzzz[k];

            g_x_0_z_xxyyyyy[k] = -g_x_0_0_xxyyyyy[k] * cd_z[k] + g_x_0_0_xxyyyyyz[k];

            g_x_0_z_xxyyyyz[k] = -g_x_0_0_xxyyyyz[k] * cd_z[k] + g_x_0_0_xxyyyyzz[k];

            g_x_0_z_xxyyyzz[k] = -g_x_0_0_xxyyyzz[k] * cd_z[k] + g_x_0_0_xxyyyzzz[k];

            g_x_0_z_xxyyzzz[k] = -g_x_0_0_xxyyzzz[k] * cd_z[k] + g_x_0_0_xxyyzzzz[k];

            g_x_0_z_xxyzzzz[k] = -g_x_0_0_xxyzzzz[k] * cd_z[k] + g_x_0_0_xxyzzzzz[k];

            g_x_0_z_xxzzzzz[k] = -g_x_0_0_xxzzzzz[k] * cd_z[k] + g_x_0_0_xxzzzzzz[k];

            g_x_0_z_xyyyyyy[k] = -g_x_0_0_xyyyyyy[k] * cd_z[k] + g_x_0_0_xyyyyyyz[k];

            g_x_0_z_xyyyyyz[k] = -g_x_0_0_xyyyyyz[k] * cd_z[k] + g_x_0_0_xyyyyyzz[k];

            g_x_0_z_xyyyyzz[k] = -g_x_0_0_xyyyyzz[k] * cd_z[k] + g_x_0_0_xyyyyzzz[k];

            g_x_0_z_xyyyzzz[k] = -g_x_0_0_xyyyzzz[k] * cd_z[k] + g_x_0_0_xyyyzzzz[k];

            g_x_0_z_xyyzzzz[k] = -g_x_0_0_xyyzzzz[k] * cd_z[k] + g_x_0_0_xyyzzzzz[k];

            g_x_0_z_xyzzzzz[k] = -g_x_0_0_xyzzzzz[k] * cd_z[k] + g_x_0_0_xyzzzzzz[k];

            g_x_0_z_xzzzzzz[k] = -g_x_0_0_xzzzzzz[k] * cd_z[k] + g_x_0_0_xzzzzzzz[k];

            g_x_0_z_yyyyyyy[k] = -g_x_0_0_yyyyyyy[k] * cd_z[k] + g_x_0_0_yyyyyyyz[k];

            g_x_0_z_yyyyyyz[k] = -g_x_0_0_yyyyyyz[k] * cd_z[k] + g_x_0_0_yyyyyyzz[k];

            g_x_0_z_yyyyyzz[k] = -g_x_0_0_yyyyyzz[k] * cd_z[k] + g_x_0_0_yyyyyzzz[k];

            g_x_0_z_yyyyzzz[k] = -g_x_0_0_yyyyzzz[k] * cd_z[k] + g_x_0_0_yyyyzzzz[k];

            g_x_0_z_yyyzzzz[k] = -g_x_0_0_yyyzzzz[k] * cd_z[k] + g_x_0_0_yyyzzzzz[k];

            g_x_0_z_yyzzzzz[k] = -g_x_0_0_yyzzzzz[k] * cd_z[k] + g_x_0_0_yyzzzzzz[k];

            g_x_0_z_yzzzzzz[k] = -g_x_0_0_yzzzzzz[k] * cd_z[k] + g_x_0_0_yzzzzzzz[k];

            g_x_0_z_zzzzzzz[k] = -g_x_0_0_zzzzzzz[k] * cd_z[k] + g_x_0_0_zzzzzzzz[k];
        }
        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_y_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 108 * acomps  + 0);

        auto g_y_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 1);

        auto g_y_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 2);

        auto g_y_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 3);

        auto g_y_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 4);

        auto g_y_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 5);

        auto g_y_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 6);

        auto g_y_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 7);

        auto g_y_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 8);

        auto g_y_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 9);

        auto g_y_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 10);

        auto g_y_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 11);

        auto g_y_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 12);

        auto g_y_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 13);

        auto g_y_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 14);

        auto g_y_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 15);

        auto g_y_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 16);

        auto g_y_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 17);

        auto g_y_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 18);

        auto g_y_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 19);

        auto g_y_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 20);

        auto g_y_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 21);

        auto g_y_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 22);

        auto g_y_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 23);

        auto g_y_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 24);

        auto g_y_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 25);

        auto g_y_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 26);

        auto g_y_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 27);

        auto g_y_0_x_yyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 28);

        auto g_y_0_x_yyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 29);

        auto g_y_0_x_yyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 30);

        auto g_y_0_x_yyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 31);

        auto g_y_0_x_yyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 32);

        auto g_y_0_x_yyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 33);

        auto g_y_0_x_yzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 34);

        auto g_y_0_x_zzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_y_0_0_xxxxxxx, g_y_0_0_xxxxxxxx, g_y_0_0_xxxxxxxy, g_y_0_0_xxxxxxxz, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxxyy, g_y_0_0_xxxxxxyz, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxxzz, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyyy, g_y_0_0_xxxxxyyz, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxyzz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxxzzz, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyyy, g_y_0_0_xxxxyyyz, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyyzz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxyzzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxxzzzz, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyyy, g_y_0_0_xxxyyyyz, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyyzz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyyzzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxyzzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxxzzzzz, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyyy, g_y_0_0_xxyyyyyz, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyyzz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyyzzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyyzzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxyzzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xxzzzzzz, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyyy, g_y_0_0_xyyyyyyz, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyyzz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyyzzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyyzzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyyzzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xyzzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_xzzzzzzz, g_y_0_0_yyyyyyy, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_zzzzzzz, g_y_0_x_xxxxxxx, g_y_0_x_xxxxxxy, g_y_0_x_xxxxxxz, g_y_0_x_xxxxxyy, g_y_0_x_xxxxxyz, g_y_0_x_xxxxxzz, g_y_0_x_xxxxyyy, g_y_0_x_xxxxyyz, g_y_0_x_xxxxyzz, g_y_0_x_xxxxzzz, g_y_0_x_xxxyyyy, g_y_0_x_xxxyyyz, g_y_0_x_xxxyyzz, g_y_0_x_xxxyzzz, g_y_0_x_xxxzzzz, g_y_0_x_xxyyyyy, g_y_0_x_xxyyyyz, g_y_0_x_xxyyyzz, g_y_0_x_xxyyzzz, g_y_0_x_xxyzzzz, g_y_0_x_xxzzzzz, g_y_0_x_xyyyyyy, g_y_0_x_xyyyyyz, g_y_0_x_xyyyyzz, g_y_0_x_xyyyzzz, g_y_0_x_xyyzzzz, g_y_0_x_xyzzzzz, g_y_0_x_xzzzzzz, g_y_0_x_yyyyyyy, g_y_0_x_yyyyyyz, g_y_0_x_yyyyyzz, g_y_0_x_yyyyzzz, g_y_0_x_yyyzzzz, g_y_0_x_yyzzzzz, g_y_0_x_yzzzzzz, g_y_0_x_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_x_xxxxxxx[k] = -g_y_0_0_xxxxxxx[k] * cd_x[k] + g_y_0_0_xxxxxxxx[k];

            g_y_0_x_xxxxxxy[k] = -g_y_0_0_xxxxxxy[k] * cd_x[k] + g_y_0_0_xxxxxxxy[k];

            g_y_0_x_xxxxxxz[k] = -g_y_0_0_xxxxxxz[k] * cd_x[k] + g_y_0_0_xxxxxxxz[k];

            g_y_0_x_xxxxxyy[k] = -g_y_0_0_xxxxxyy[k] * cd_x[k] + g_y_0_0_xxxxxxyy[k];

            g_y_0_x_xxxxxyz[k] = -g_y_0_0_xxxxxyz[k] * cd_x[k] + g_y_0_0_xxxxxxyz[k];

            g_y_0_x_xxxxxzz[k] = -g_y_0_0_xxxxxzz[k] * cd_x[k] + g_y_0_0_xxxxxxzz[k];

            g_y_0_x_xxxxyyy[k] = -g_y_0_0_xxxxyyy[k] * cd_x[k] + g_y_0_0_xxxxxyyy[k];

            g_y_0_x_xxxxyyz[k] = -g_y_0_0_xxxxyyz[k] * cd_x[k] + g_y_0_0_xxxxxyyz[k];

            g_y_0_x_xxxxyzz[k] = -g_y_0_0_xxxxyzz[k] * cd_x[k] + g_y_0_0_xxxxxyzz[k];

            g_y_0_x_xxxxzzz[k] = -g_y_0_0_xxxxzzz[k] * cd_x[k] + g_y_0_0_xxxxxzzz[k];

            g_y_0_x_xxxyyyy[k] = -g_y_0_0_xxxyyyy[k] * cd_x[k] + g_y_0_0_xxxxyyyy[k];

            g_y_0_x_xxxyyyz[k] = -g_y_0_0_xxxyyyz[k] * cd_x[k] + g_y_0_0_xxxxyyyz[k];

            g_y_0_x_xxxyyzz[k] = -g_y_0_0_xxxyyzz[k] * cd_x[k] + g_y_0_0_xxxxyyzz[k];

            g_y_0_x_xxxyzzz[k] = -g_y_0_0_xxxyzzz[k] * cd_x[k] + g_y_0_0_xxxxyzzz[k];

            g_y_0_x_xxxzzzz[k] = -g_y_0_0_xxxzzzz[k] * cd_x[k] + g_y_0_0_xxxxzzzz[k];

            g_y_0_x_xxyyyyy[k] = -g_y_0_0_xxyyyyy[k] * cd_x[k] + g_y_0_0_xxxyyyyy[k];

            g_y_0_x_xxyyyyz[k] = -g_y_0_0_xxyyyyz[k] * cd_x[k] + g_y_0_0_xxxyyyyz[k];

            g_y_0_x_xxyyyzz[k] = -g_y_0_0_xxyyyzz[k] * cd_x[k] + g_y_0_0_xxxyyyzz[k];

            g_y_0_x_xxyyzzz[k] = -g_y_0_0_xxyyzzz[k] * cd_x[k] + g_y_0_0_xxxyyzzz[k];

            g_y_0_x_xxyzzzz[k] = -g_y_0_0_xxyzzzz[k] * cd_x[k] + g_y_0_0_xxxyzzzz[k];

            g_y_0_x_xxzzzzz[k] = -g_y_0_0_xxzzzzz[k] * cd_x[k] + g_y_0_0_xxxzzzzz[k];

            g_y_0_x_xyyyyyy[k] = -g_y_0_0_xyyyyyy[k] * cd_x[k] + g_y_0_0_xxyyyyyy[k];

            g_y_0_x_xyyyyyz[k] = -g_y_0_0_xyyyyyz[k] * cd_x[k] + g_y_0_0_xxyyyyyz[k];

            g_y_0_x_xyyyyzz[k] = -g_y_0_0_xyyyyzz[k] * cd_x[k] + g_y_0_0_xxyyyyzz[k];

            g_y_0_x_xyyyzzz[k] = -g_y_0_0_xyyyzzz[k] * cd_x[k] + g_y_0_0_xxyyyzzz[k];

            g_y_0_x_xyyzzzz[k] = -g_y_0_0_xyyzzzz[k] * cd_x[k] + g_y_0_0_xxyyzzzz[k];

            g_y_0_x_xyzzzzz[k] = -g_y_0_0_xyzzzzz[k] * cd_x[k] + g_y_0_0_xxyzzzzz[k];

            g_y_0_x_xzzzzzz[k] = -g_y_0_0_xzzzzzz[k] * cd_x[k] + g_y_0_0_xxzzzzzz[k];

            g_y_0_x_yyyyyyy[k] = -g_y_0_0_yyyyyyy[k] * cd_x[k] + g_y_0_0_xyyyyyyy[k];

            g_y_0_x_yyyyyyz[k] = -g_y_0_0_yyyyyyz[k] * cd_x[k] + g_y_0_0_xyyyyyyz[k];

            g_y_0_x_yyyyyzz[k] = -g_y_0_0_yyyyyzz[k] * cd_x[k] + g_y_0_0_xyyyyyzz[k];

            g_y_0_x_yyyyzzz[k] = -g_y_0_0_yyyyzzz[k] * cd_x[k] + g_y_0_0_xyyyyzzz[k];

            g_y_0_x_yyyzzzz[k] = -g_y_0_0_yyyzzzz[k] * cd_x[k] + g_y_0_0_xyyyzzzz[k];

            g_y_0_x_yyzzzzz[k] = -g_y_0_0_yyzzzzz[k] * cd_x[k] + g_y_0_0_xyyzzzzz[k];

            g_y_0_x_yzzzzzz[k] = -g_y_0_0_yzzzzzz[k] * cd_x[k] + g_y_0_0_xyzzzzzz[k];

            g_y_0_x_zzzzzzz[k] = -g_y_0_0_zzzzzzz[k] * cd_x[k] + g_y_0_0_xzzzzzzz[k];
        }

        /// Set up 36-72 components of targeted buffer : cbuffer.data(

        auto g_y_0_y_xxxxxxx = cbuffer.data(pk_geom_10_off + 108 * acomps  + 36);

        auto g_y_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 37);

        auto g_y_0_y_xxxxxxz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 38);

        auto g_y_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 39);

        auto g_y_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 40);

        auto g_y_0_y_xxxxxzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 41);

        auto g_y_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 42);

        auto g_y_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 43);

        auto g_y_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 44);

        auto g_y_0_y_xxxxzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 45);

        auto g_y_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 46);

        auto g_y_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 47);

        auto g_y_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 48);

        auto g_y_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 49);

        auto g_y_0_y_xxxzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 50);

        auto g_y_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 51);

        auto g_y_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 52);

        auto g_y_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 53);

        auto g_y_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 54);

        auto g_y_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 55);

        auto g_y_0_y_xxzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 56);

        auto g_y_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 57);

        auto g_y_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 58);

        auto g_y_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 59);

        auto g_y_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 60);

        auto g_y_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 61);

        auto g_y_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 62);

        auto g_y_0_y_xzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 63);

        auto g_y_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 64);

        auto g_y_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 65);

        auto g_y_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 66);

        auto g_y_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 67);

        auto g_y_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 68);

        auto g_y_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 69);

        auto g_y_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 70);

        auto g_y_0_y_zzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 71);

        #pragma omp simd aligned(cd_y, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz, g_y_0_0_xxxxxxx, g_y_0_0_xxxxxxxy, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxxyy, g_y_0_0_xxxxxxyz, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyyy, g_y_0_0_xxxxxyyz, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxyzz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyyy, g_y_0_0_xxxxyyyz, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyyzz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxyzzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyyy, g_y_0_0_xxxyyyyz, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyyzz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyyzzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxyzzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyyy, g_y_0_0_xxyyyyyz, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyyzz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyyzzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyyzzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxyzzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyyy, g_y_0_0_xyyyyyyz, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyyzz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyyzzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyyzzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyyzzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xyzzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_yyyyyyy, g_y_0_0_yyyyyyyy, g_y_0_0_yyyyyyyz, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyyzz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyyzzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyyzzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyyzzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yyzzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_yzzzzzzz, g_y_0_0_zzzzzzz, g_y_0_y_xxxxxxx, g_y_0_y_xxxxxxy, g_y_0_y_xxxxxxz, g_y_0_y_xxxxxyy, g_y_0_y_xxxxxyz, g_y_0_y_xxxxxzz, g_y_0_y_xxxxyyy, g_y_0_y_xxxxyyz, g_y_0_y_xxxxyzz, g_y_0_y_xxxxzzz, g_y_0_y_xxxyyyy, g_y_0_y_xxxyyyz, g_y_0_y_xxxyyzz, g_y_0_y_xxxyzzz, g_y_0_y_xxxzzzz, g_y_0_y_xxyyyyy, g_y_0_y_xxyyyyz, g_y_0_y_xxyyyzz, g_y_0_y_xxyyzzz, g_y_0_y_xxyzzzz, g_y_0_y_xxzzzzz, g_y_0_y_xyyyyyy, g_y_0_y_xyyyyyz, g_y_0_y_xyyyyzz, g_y_0_y_xyyyzzz, g_y_0_y_xyyzzzz, g_y_0_y_xyzzzzz, g_y_0_y_xzzzzzz, g_y_0_y_yyyyyyy, g_y_0_y_yyyyyyz, g_y_0_y_yyyyyzz, g_y_0_y_yyyyzzz, g_y_0_y_yyyzzzz, g_y_0_y_yyzzzzz, g_y_0_y_yzzzzzz, g_y_0_y_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_y_xxxxxxx[k] = -g_0_xxxxxxx[k] - g_y_0_0_xxxxxxx[k] * cd_y[k] + g_y_0_0_xxxxxxxy[k];

            g_y_0_y_xxxxxxy[k] = -g_0_xxxxxxy[k] - g_y_0_0_xxxxxxy[k] * cd_y[k] + g_y_0_0_xxxxxxyy[k];

            g_y_0_y_xxxxxxz[k] = -g_0_xxxxxxz[k] - g_y_0_0_xxxxxxz[k] * cd_y[k] + g_y_0_0_xxxxxxyz[k];

            g_y_0_y_xxxxxyy[k] = -g_0_xxxxxyy[k] - g_y_0_0_xxxxxyy[k] * cd_y[k] + g_y_0_0_xxxxxyyy[k];

            g_y_0_y_xxxxxyz[k] = -g_0_xxxxxyz[k] - g_y_0_0_xxxxxyz[k] * cd_y[k] + g_y_0_0_xxxxxyyz[k];

            g_y_0_y_xxxxxzz[k] = -g_0_xxxxxzz[k] - g_y_0_0_xxxxxzz[k] * cd_y[k] + g_y_0_0_xxxxxyzz[k];

            g_y_0_y_xxxxyyy[k] = -g_0_xxxxyyy[k] - g_y_0_0_xxxxyyy[k] * cd_y[k] + g_y_0_0_xxxxyyyy[k];

            g_y_0_y_xxxxyyz[k] = -g_0_xxxxyyz[k] - g_y_0_0_xxxxyyz[k] * cd_y[k] + g_y_0_0_xxxxyyyz[k];

            g_y_0_y_xxxxyzz[k] = -g_0_xxxxyzz[k] - g_y_0_0_xxxxyzz[k] * cd_y[k] + g_y_0_0_xxxxyyzz[k];

            g_y_0_y_xxxxzzz[k] = -g_0_xxxxzzz[k] - g_y_0_0_xxxxzzz[k] * cd_y[k] + g_y_0_0_xxxxyzzz[k];

            g_y_0_y_xxxyyyy[k] = -g_0_xxxyyyy[k] - g_y_0_0_xxxyyyy[k] * cd_y[k] + g_y_0_0_xxxyyyyy[k];

            g_y_0_y_xxxyyyz[k] = -g_0_xxxyyyz[k] - g_y_0_0_xxxyyyz[k] * cd_y[k] + g_y_0_0_xxxyyyyz[k];

            g_y_0_y_xxxyyzz[k] = -g_0_xxxyyzz[k] - g_y_0_0_xxxyyzz[k] * cd_y[k] + g_y_0_0_xxxyyyzz[k];

            g_y_0_y_xxxyzzz[k] = -g_0_xxxyzzz[k] - g_y_0_0_xxxyzzz[k] * cd_y[k] + g_y_0_0_xxxyyzzz[k];

            g_y_0_y_xxxzzzz[k] = -g_0_xxxzzzz[k] - g_y_0_0_xxxzzzz[k] * cd_y[k] + g_y_0_0_xxxyzzzz[k];

            g_y_0_y_xxyyyyy[k] = -g_0_xxyyyyy[k] - g_y_0_0_xxyyyyy[k] * cd_y[k] + g_y_0_0_xxyyyyyy[k];

            g_y_0_y_xxyyyyz[k] = -g_0_xxyyyyz[k] - g_y_0_0_xxyyyyz[k] * cd_y[k] + g_y_0_0_xxyyyyyz[k];

            g_y_0_y_xxyyyzz[k] = -g_0_xxyyyzz[k] - g_y_0_0_xxyyyzz[k] * cd_y[k] + g_y_0_0_xxyyyyzz[k];

            g_y_0_y_xxyyzzz[k] = -g_0_xxyyzzz[k] - g_y_0_0_xxyyzzz[k] * cd_y[k] + g_y_0_0_xxyyyzzz[k];

            g_y_0_y_xxyzzzz[k] = -g_0_xxyzzzz[k] - g_y_0_0_xxyzzzz[k] * cd_y[k] + g_y_0_0_xxyyzzzz[k];

            g_y_0_y_xxzzzzz[k] = -g_0_xxzzzzz[k] - g_y_0_0_xxzzzzz[k] * cd_y[k] + g_y_0_0_xxyzzzzz[k];

            g_y_0_y_xyyyyyy[k] = -g_0_xyyyyyy[k] - g_y_0_0_xyyyyyy[k] * cd_y[k] + g_y_0_0_xyyyyyyy[k];

            g_y_0_y_xyyyyyz[k] = -g_0_xyyyyyz[k] - g_y_0_0_xyyyyyz[k] * cd_y[k] + g_y_0_0_xyyyyyyz[k];

            g_y_0_y_xyyyyzz[k] = -g_0_xyyyyzz[k] - g_y_0_0_xyyyyzz[k] * cd_y[k] + g_y_0_0_xyyyyyzz[k];

            g_y_0_y_xyyyzzz[k] = -g_0_xyyyzzz[k] - g_y_0_0_xyyyzzz[k] * cd_y[k] + g_y_0_0_xyyyyzzz[k];

            g_y_0_y_xyyzzzz[k] = -g_0_xyyzzzz[k] - g_y_0_0_xyyzzzz[k] * cd_y[k] + g_y_0_0_xyyyzzzz[k];

            g_y_0_y_xyzzzzz[k] = -g_0_xyzzzzz[k] - g_y_0_0_xyzzzzz[k] * cd_y[k] + g_y_0_0_xyyzzzzz[k];

            g_y_0_y_xzzzzzz[k] = -g_0_xzzzzzz[k] - g_y_0_0_xzzzzzz[k] * cd_y[k] + g_y_0_0_xyzzzzzz[k];

            g_y_0_y_yyyyyyy[k] = -g_0_yyyyyyy[k] - g_y_0_0_yyyyyyy[k] * cd_y[k] + g_y_0_0_yyyyyyyy[k];

            g_y_0_y_yyyyyyz[k] = -g_0_yyyyyyz[k] - g_y_0_0_yyyyyyz[k] * cd_y[k] + g_y_0_0_yyyyyyyz[k];

            g_y_0_y_yyyyyzz[k] = -g_0_yyyyyzz[k] - g_y_0_0_yyyyyzz[k] * cd_y[k] + g_y_0_0_yyyyyyzz[k];

            g_y_0_y_yyyyzzz[k] = -g_0_yyyyzzz[k] - g_y_0_0_yyyyzzz[k] * cd_y[k] + g_y_0_0_yyyyyzzz[k];

            g_y_0_y_yyyzzzz[k] = -g_0_yyyzzzz[k] - g_y_0_0_yyyzzzz[k] * cd_y[k] + g_y_0_0_yyyyzzzz[k];

            g_y_0_y_yyzzzzz[k] = -g_0_yyzzzzz[k] - g_y_0_0_yyzzzzz[k] * cd_y[k] + g_y_0_0_yyyzzzzz[k];

            g_y_0_y_yzzzzzz[k] = -g_0_yzzzzzz[k] - g_y_0_0_yzzzzzz[k] * cd_y[k] + g_y_0_0_yyzzzzzz[k];

            g_y_0_y_zzzzzzz[k] = -g_0_zzzzzzz[k] - g_y_0_0_zzzzzzz[k] * cd_y[k] + g_y_0_0_yzzzzzzz[k];
        }

        /// Set up 72-108 components of targeted buffer : cbuffer.data(

        auto g_y_0_z_xxxxxxx = cbuffer.data(pk_geom_10_off + 108 * acomps  + 72);

        auto g_y_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 73);

        auto g_y_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 74);

        auto g_y_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 75);

        auto g_y_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 76);

        auto g_y_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 77);

        auto g_y_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 78);

        auto g_y_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 79);

        auto g_y_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 80);

        auto g_y_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 81);

        auto g_y_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 82);

        auto g_y_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 83);

        auto g_y_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 84);

        auto g_y_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 85);

        auto g_y_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 86);

        auto g_y_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 87);

        auto g_y_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 88);

        auto g_y_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 89);

        auto g_y_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 90);

        auto g_y_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 91);

        auto g_y_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 92);

        auto g_y_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 93);

        auto g_y_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 94);

        auto g_y_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 95);

        auto g_y_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 96);

        auto g_y_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 97);

        auto g_y_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 98);

        auto g_y_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 99);

        auto g_y_0_z_yyyyyyy = cbuffer.data(pk_geom_10_off + 108 * acomps  + 100);

        auto g_y_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 101);

        auto g_y_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 102);

        auto g_y_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 103);

        auto g_y_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 104);

        auto g_y_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 105);

        auto g_y_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 106);

        auto g_y_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 108 * acomps  + 107);

        #pragma omp simd aligned(cd_z, g_y_0_0_xxxxxxx, g_y_0_0_xxxxxxxz, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxxyz, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxxzz, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyyz, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxyzz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxxzzz, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyyz, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyyzz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxyzzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxxzzzz, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyyz, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyyzz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyyzzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxyzzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxxzzzzz, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyyz, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyyzz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyyzzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyyzzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxyzzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xxzzzzzz, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyyz, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyyzz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyyzzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyyzzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyyzzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xyzzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_xzzzzzzz, g_y_0_0_yyyyyyy, g_y_0_0_yyyyyyyz, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyyzz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyyzzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyyzzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyyzzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yyzzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_yzzzzzzz, g_y_0_0_zzzzzzz, g_y_0_0_zzzzzzzz, g_y_0_z_xxxxxxx, g_y_0_z_xxxxxxy, g_y_0_z_xxxxxxz, g_y_0_z_xxxxxyy, g_y_0_z_xxxxxyz, g_y_0_z_xxxxxzz, g_y_0_z_xxxxyyy, g_y_0_z_xxxxyyz, g_y_0_z_xxxxyzz, g_y_0_z_xxxxzzz, g_y_0_z_xxxyyyy, g_y_0_z_xxxyyyz, g_y_0_z_xxxyyzz, g_y_0_z_xxxyzzz, g_y_0_z_xxxzzzz, g_y_0_z_xxyyyyy, g_y_0_z_xxyyyyz, g_y_0_z_xxyyyzz, g_y_0_z_xxyyzzz, g_y_0_z_xxyzzzz, g_y_0_z_xxzzzzz, g_y_0_z_xyyyyyy, g_y_0_z_xyyyyyz, g_y_0_z_xyyyyzz, g_y_0_z_xyyyzzz, g_y_0_z_xyyzzzz, g_y_0_z_xyzzzzz, g_y_0_z_xzzzzzz, g_y_0_z_yyyyyyy, g_y_0_z_yyyyyyz, g_y_0_z_yyyyyzz, g_y_0_z_yyyyzzz, g_y_0_z_yyyzzzz, g_y_0_z_yyzzzzz, g_y_0_z_yzzzzzz, g_y_0_z_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_z_xxxxxxx[k] = -g_y_0_0_xxxxxxx[k] * cd_z[k] + g_y_0_0_xxxxxxxz[k];

            g_y_0_z_xxxxxxy[k] = -g_y_0_0_xxxxxxy[k] * cd_z[k] + g_y_0_0_xxxxxxyz[k];

            g_y_0_z_xxxxxxz[k] = -g_y_0_0_xxxxxxz[k] * cd_z[k] + g_y_0_0_xxxxxxzz[k];

            g_y_0_z_xxxxxyy[k] = -g_y_0_0_xxxxxyy[k] * cd_z[k] + g_y_0_0_xxxxxyyz[k];

            g_y_0_z_xxxxxyz[k] = -g_y_0_0_xxxxxyz[k] * cd_z[k] + g_y_0_0_xxxxxyzz[k];

            g_y_0_z_xxxxxzz[k] = -g_y_0_0_xxxxxzz[k] * cd_z[k] + g_y_0_0_xxxxxzzz[k];

            g_y_0_z_xxxxyyy[k] = -g_y_0_0_xxxxyyy[k] * cd_z[k] + g_y_0_0_xxxxyyyz[k];

            g_y_0_z_xxxxyyz[k] = -g_y_0_0_xxxxyyz[k] * cd_z[k] + g_y_0_0_xxxxyyzz[k];

            g_y_0_z_xxxxyzz[k] = -g_y_0_0_xxxxyzz[k] * cd_z[k] + g_y_0_0_xxxxyzzz[k];

            g_y_0_z_xxxxzzz[k] = -g_y_0_0_xxxxzzz[k] * cd_z[k] + g_y_0_0_xxxxzzzz[k];

            g_y_0_z_xxxyyyy[k] = -g_y_0_0_xxxyyyy[k] * cd_z[k] + g_y_0_0_xxxyyyyz[k];

            g_y_0_z_xxxyyyz[k] = -g_y_0_0_xxxyyyz[k] * cd_z[k] + g_y_0_0_xxxyyyzz[k];

            g_y_0_z_xxxyyzz[k] = -g_y_0_0_xxxyyzz[k] * cd_z[k] + g_y_0_0_xxxyyzzz[k];

            g_y_0_z_xxxyzzz[k] = -g_y_0_0_xxxyzzz[k] * cd_z[k] + g_y_0_0_xxxyzzzz[k];

            g_y_0_z_xxxzzzz[k] = -g_y_0_0_xxxzzzz[k] * cd_z[k] + g_y_0_0_xxxzzzzz[k];

            g_y_0_z_xxyyyyy[k] = -g_y_0_0_xxyyyyy[k] * cd_z[k] + g_y_0_0_xxyyyyyz[k];

            g_y_0_z_xxyyyyz[k] = -g_y_0_0_xxyyyyz[k] * cd_z[k] + g_y_0_0_xxyyyyzz[k];

            g_y_0_z_xxyyyzz[k] = -g_y_0_0_xxyyyzz[k] * cd_z[k] + g_y_0_0_xxyyyzzz[k];

            g_y_0_z_xxyyzzz[k] = -g_y_0_0_xxyyzzz[k] * cd_z[k] + g_y_0_0_xxyyzzzz[k];

            g_y_0_z_xxyzzzz[k] = -g_y_0_0_xxyzzzz[k] * cd_z[k] + g_y_0_0_xxyzzzzz[k];

            g_y_0_z_xxzzzzz[k] = -g_y_0_0_xxzzzzz[k] * cd_z[k] + g_y_0_0_xxzzzzzz[k];

            g_y_0_z_xyyyyyy[k] = -g_y_0_0_xyyyyyy[k] * cd_z[k] + g_y_0_0_xyyyyyyz[k];

            g_y_0_z_xyyyyyz[k] = -g_y_0_0_xyyyyyz[k] * cd_z[k] + g_y_0_0_xyyyyyzz[k];

            g_y_0_z_xyyyyzz[k] = -g_y_0_0_xyyyyzz[k] * cd_z[k] + g_y_0_0_xyyyyzzz[k];

            g_y_0_z_xyyyzzz[k] = -g_y_0_0_xyyyzzz[k] * cd_z[k] + g_y_0_0_xyyyzzzz[k];

            g_y_0_z_xyyzzzz[k] = -g_y_0_0_xyyzzzz[k] * cd_z[k] + g_y_0_0_xyyzzzzz[k];

            g_y_0_z_xyzzzzz[k] = -g_y_0_0_xyzzzzz[k] * cd_z[k] + g_y_0_0_xyzzzzzz[k];

            g_y_0_z_xzzzzzz[k] = -g_y_0_0_xzzzzzz[k] * cd_z[k] + g_y_0_0_xzzzzzzz[k];

            g_y_0_z_yyyyyyy[k] = -g_y_0_0_yyyyyyy[k] * cd_z[k] + g_y_0_0_yyyyyyyz[k];

            g_y_0_z_yyyyyyz[k] = -g_y_0_0_yyyyyyz[k] * cd_z[k] + g_y_0_0_yyyyyyzz[k];

            g_y_0_z_yyyyyzz[k] = -g_y_0_0_yyyyyzz[k] * cd_z[k] + g_y_0_0_yyyyyzzz[k];

            g_y_0_z_yyyyzzz[k] = -g_y_0_0_yyyyzzz[k] * cd_z[k] + g_y_0_0_yyyyzzzz[k];

            g_y_0_z_yyyzzzz[k] = -g_y_0_0_yyyzzzz[k] * cd_z[k] + g_y_0_0_yyyzzzzz[k];

            g_y_0_z_yyzzzzz[k] = -g_y_0_0_yyzzzzz[k] * cd_z[k] + g_y_0_0_yyzzzzzz[k];

            g_y_0_z_yzzzzzz[k] = -g_y_0_0_yzzzzzz[k] * cd_z[k] + g_y_0_0_yzzzzzzz[k];

            g_y_0_z_zzzzzzz[k] = -g_y_0_0_zzzzzzz[k] * cd_z[k] + g_y_0_0_zzzzzzzz[k];
        }
        /// Set up 0-36 components of targeted buffer : cbuffer.data(

        auto g_z_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 216 * acomps  + 0);

        auto g_z_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 1);

        auto g_z_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 2);

        auto g_z_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 3);

        auto g_z_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 4);

        auto g_z_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 5);

        auto g_z_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 6);

        auto g_z_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 7);

        auto g_z_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 8);

        auto g_z_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 9);

        auto g_z_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 10);

        auto g_z_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 11);

        auto g_z_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 12);

        auto g_z_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 13);

        auto g_z_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 14);

        auto g_z_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 15);

        auto g_z_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 16);

        auto g_z_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 17);

        auto g_z_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 18);

        auto g_z_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 19);

        auto g_z_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 20);

        auto g_z_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 21);

        auto g_z_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 22);

        auto g_z_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 23);

        auto g_z_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 24);

        auto g_z_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 25);

        auto g_z_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 26);

        auto g_z_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 27);

        auto g_z_0_x_yyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 28);

        auto g_z_0_x_yyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 29);

        auto g_z_0_x_yyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 30);

        auto g_z_0_x_yyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 31);

        auto g_z_0_x_yyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 32);

        auto g_z_0_x_yyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 33);

        auto g_z_0_x_yzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 34);

        auto g_z_0_x_zzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_z_0_0_xxxxxxx, g_z_0_0_xxxxxxxx, g_z_0_0_xxxxxxxy, g_z_0_0_xxxxxxxz, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxxyy, g_z_0_0_xxxxxxyz, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxxzz, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyyy, g_z_0_0_xxxxxyyz, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxyzz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxxzzz, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyyy, g_z_0_0_xxxxyyyz, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyyzz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxyzzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxxzzzz, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyyy, g_z_0_0_xxxyyyyz, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyyzz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyyzzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxyzzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxxzzzzz, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyyy, g_z_0_0_xxyyyyyz, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyyzz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyyzzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyyzzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxyzzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xxzzzzzz, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyyy, g_z_0_0_xyyyyyyz, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyyzz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyyzzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyyzzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyyzzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xyzzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_xzzzzzzz, g_z_0_0_yyyyyyy, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_zzzzzzz, g_z_0_x_xxxxxxx, g_z_0_x_xxxxxxy, g_z_0_x_xxxxxxz, g_z_0_x_xxxxxyy, g_z_0_x_xxxxxyz, g_z_0_x_xxxxxzz, g_z_0_x_xxxxyyy, g_z_0_x_xxxxyyz, g_z_0_x_xxxxyzz, g_z_0_x_xxxxzzz, g_z_0_x_xxxyyyy, g_z_0_x_xxxyyyz, g_z_0_x_xxxyyzz, g_z_0_x_xxxyzzz, g_z_0_x_xxxzzzz, g_z_0_x_xxyyyyy, g_z_0_x_xxyyyyz, g_z_0_x_xxyyyzz, g_z_0_x_xxyyzzz, g_z_0_x_xxyzzzz, g_z_0_x_xxzzzzz, g_z_0_x_xyyyyyy, g_z_0_x_xyyyyyz, g_z_0_x_xyyyyzz, g_z_0_x_xyyyzzz, g_z_0_x_xyyzzzz, g_z_0_x_xyzzzzz, g_z_0_x_xzzzzzz, g_z_0_x_yyyyyyy, g_z_0_x_yyyyyyz, g_z_0_x_yyyyyzz, g_z_0_x_yyyyzzz, g_z_0_x_yyyzzzz, g_z_0_x_yyzzzzz, g_z_0_x_yzzzzzz, g_z_0_x_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_x_xxxxxxx[k] = -g_z_0_0_xxxxxxx[k] * cd_x[k] + g_z_0_0_xxxxxxxx[k];

            g_z_0_x_xxxxxxy[k] = -g_z_0_0_xxxxxxy[k] * cd_x[k] + g_z_0_0_xxxxxxxy[k];

            g_z_0_x_xxxxxxz[k] = -g_z_0_0_xxxxxxz[k] * cd_x[k] + g_z_0_0_xxxxxxxz[k];

            g_z_0_x_xxxxxyy[k] = -g_z_0_0_xxxxxyy[k] * cd_x[k] + g_z_0_0_xxxxxxyy[k];

            g_z_0_x_xxxxxyz[k] = -g_z_0_0_xxxxxyz[k] * cd_x[k] + g_z_0_0_xxxxxxyz[k];

            g_z_0_x_xxxxxzz[k] = -g_z_0_0_xxxxxzz[k] * cd_x[k] + g_z_0_0_xxxxxxzz[k];

            g_z_0_x_xxxxyyy[k] = -g_z_0_0_xxxxyyy[k] * cd_x[k] + g_z_0_0_xxxxxyyy[k];

            g_z_0_x_xxxxyyz[k] = -g_z_0_0_xxxxyyz[k] * cd_x[k] + g_z_0_0_xxxxxyyz[k];

            g_z_0_x_xxxxyzz[k] = -g_z_0_0_xxxxyzz[k] * cd_x[k] + g_z_0_0_xxxxxyzz[k];

            g_z_0_x_xxxxzzz[k] = -g_z_0_0_xxxxzzz[k] * cd_x[k] + g_z_0_0_xxxxxzzz[k];

            g_z_0_x_xxxyyyy[k] = -g_z_0_0_xxxyyyy[k] * cd_x[k] + g_z_0_0_xxxxyyyy[k];

            g_z_0_x_xxxyyyz[k] = -g_z_0_0_xxxyyyz[k] * cd_x[k] + g_z_0_0_xxxxyyyz[k];

            g_z_0_x_xxxyyzz[k] = -g_z_0_0_xxxyyzz[k] * cd_x[k] + g_z_0_0_xxxxyyzz[k];

            g_z_0_x_xxxyzzz[k] = -g_z_0_0_xxxyzzz[k] * cd_x[k] + g_z_0_0_xxxxyzzz[k];

            g_z_0_x_xxxzzzz[k] = -g_z_0_0_xxxzzzz[k] * cd_x[k] + g_z_0_0_xxxxzzzz[k];

            g_z_0_x_xxyyyyy[k] = -g_z_0_0_xxyyyyy[k] * cd_x[k] + g_z_0_0_xxxyyyyy[k];

            g_z_0_x_xxyyyyz[k] = -g_z_0_0_xxyyyyz[k] * cd_x[k] + g_z_0_0_xxxyyyyz[k];

            g_z_0_x_xxyyyzz[k] = -g_z_0_0_xxyyyzz[k] * cd_x[k] + g_z_0_0_xxxyyyzz[k];

            g_z_0_x_xxyyzzz[k] = -g_z_0_0_xxyyzzz[k] * cd_x[k] + g_z_0_0_xxxyyzzz[k];

            g_z_0_x_xxyzzzz[k] = -g_z_0_0_xxyzzzz[k] * cd_x[k] + g_z_0_0_xxxyzzzz[k];

            g_z_0_x_xxzzzzz[k] = -g_z_0_0_xxzzzzz[k] * cd_x[k] + g_z_0_0_xxxzzzzz[k];

            g_z_0_x_xyyyyyy[k] = -g_z_0_0_xyyyyyy[k] * cd_x[k] + g_z_0_0_xxyyyyyy[k];

            g_z_0_x_xyyyyyz[k] = -g_z_0_0_xyyyyyz[k] * cd_x[k] + g_z_0_0_xxyyyyyz[k];

            g_z_0_x_xyyyyzz[k] = -g_z_0_0_xyyyyzz[k] * cd_x[k] + g_z_0_0_xxyyyyzz[k];

            g_z_0_x_xyyyzzz[k] = -g_z_0_0_xyyyzzz[k] * cd_x[k] + g_z_0_0_xxyyyzzz[k];

            g_z_0_x_xyyzzzz[k] = -g_z_0_0_xyyzzzz[k] * cd_x[k] + g_z_0_0_xxyyzzzz[k];

            g_z_0_x_xyzzzzz[k] = -g_z_0_0_xyzzzzz[k] * cd_x[k] + g_z_0_0_xxyzzzzz[k];

            g_z_0_x_xzzzzzz[k] = -g_z_0_0_xzzzzzz[k] * cd_x[k] + g_z_0_0_xxzzzzzz[k];

            g_z_0_x_yyyyyyy[k] = -g_z_0_0_yyyyyyy[k] * cd_x[k] + g_z_0_0_xyyyyyyy[k];

            g_z_0_x_yyyyyyz[k] = -g_z_0_0_yyyyyyz[k] * cd_x[k] + g_z_0_0_xyyyyyyz[k];

            g_z_0_x_yyyyyzz[k] = -g_z_0_0_yyyyyzz[k] * cd_x[k] + g_z_0_0_xyyyyyzz[k];

            g_z_0_x_yyyyzzz[k] = -g_z_0_0_yyyyzzz[k] * cd_x[k] + g_z_0_0_xyyyyzzz[k];

            g_z_0_x_yyyzzzz[k] = -g_z_0_0_yyyzzzz[k] * cd_x[k] + g_z_0_0_xyyyzzzz[k];

            g_z_0_x_yyzzzzz[k] = -g_z_0_0_yyzzzzz[k] * cd_x[k] + g_z_0_0_xyyzzzzz[k];

            g_z_0_x_yzzzzzz[k] = -g_z_0_0_yzzzzzz[k] * cd_x[k] + g_z_0_0_xyzzzzzz[k];

            g_z_0_x_zzzzzzz[k] = -g_z_0_0_zzzzzzz[k] * cd_x[k] + g_z_0_0_xzzzzzzz[k];
        }

        /// Set up 36-72 components of targeted buffer : cbuffer.data(

        auto g_z_0_y_xxxxxxx = cbuffer.data(pk_geom_10_off + 216 * acomps  + 36);

        auto g_z_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 37);

        auto g_z_0_y_xxxxxxz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 38);

        auto g_z_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 39);

        auto g_z_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 40);

        auto g_z_0_y_xxxxxzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 41);

        auto g_z_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 42);

        auto g_z_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 43);

        auto g_z_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 44);

        auto g_z_0_y_xxxxzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 45);

        auto g_z_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 46);

        auto g_z_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 47);

        auto g_z_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 48);

        auto g_z_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 49);

        auto g_z_0_y_xxxzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 50);

        auto g_z_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 51);

        auto g_z_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 52);

        auto g_z_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 53);

        auto g_z_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 54);

        auto g_z_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 55);

        auto g_z_0_y_xxzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 56);

        auto g_z_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 57);

        auto g_z_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 58);

        auto g_z_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 59);

        auto g_z_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 60);

        auto g_z_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 61);

        auto g_z_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 62);

        auto g_z_0_y_xzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 63);

        auto g_z_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 64);

        auto g_z_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 65);

        auto g_z_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 66);

        auto g_z_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 67);

        auto g_z_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 68);

        auto g_z_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 69);

        auto g_z_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 70);

        auto g_z_0_y_zzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 71);

        #pragma omp simd aligned(cd_y, g_z_0_0_xxxxxxx, g_z_0_0_xxxxxxxy, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxxyy, g_z_0_0_xxxxxxyz, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyyy, g_z_0_0_xxxxxyyz, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxyzz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyyy, g_z_0_0_xxxxyyyz, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyyzz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxyzzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyyy, g_z_0_0_xxxyyyyz, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyyzz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyyzzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxyzzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyyy, g_z_0_0_xxyyyyyz, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyyzz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyyzzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyyzzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxyzzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyyy, g_z_0_0_xyyyyyyz, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyyzz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyyzzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyyzzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyyzzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xyzzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_yyyyyyy, g_z_0_0_yyyyyyyy, g_z_0_0_yyyyyyyz, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyyzz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyyzzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyyzzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyyzzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yyzzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_yzzzzzzz, g_z_0_0_zzzzzzz, g_z_0_y_xxxxxxx, g_z_0_y_xxxxxxy, g_z_0_y_xxxxxxz, g_z_0_y_xxxxxyy, g_z_0_y_xxxxxyz, g_z_0_y_xxxxxzz, g_z_0_y_xxxxyyy, g_z_0_y_xxxxyyz, g_z_0_y_xxxxyzz, g_z_0_y_xxxxzzz, g_z_0_y_xxxyyyy, g_z_0_y_xxxyyyz, g_z_0_y_xxxyyzz, g_z_0_y_xxxyzzz, g_z_0_y_xxxzzzz, g_z_0_y_xxyyyyy, g_z_0_y_xxyyyyz, g_z_0_y_xxyyyzz, g_z_0_y_xxyyzzz, g_z_0_y_xxyzzzz, g_z_0_y_xxzzzzz, g_z_0_y_xyyyyyy, g_z_0_y_xyyyyyz, g_z_0_y_xyyyyzz, g_z_0_y_xyyyzzz, g_z_0_y_xyyzzzz, g_z_0_y_xyzzzzz, g_z_0_y_xzzzzzz, g_z_0_y_yyyyyyy, g_z_0_y_yyyyyyz, g_z_0_y_yyyyyzz, g_z_0_y_yyyyzzz, g_z_0_y_yyyzzzz, g_z_0_y_yyzzzzz, g_z_0_y_yzzzzzz, g_z_0_y_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_y_xxxxxxx[k] = -g_z_0_0_xxxxxxx[k] * cd_y[k] + g_z_0_0_xxxxxxxy[k];

            g_z_0_y_xxxxxxy[k] = -g_z_0_0_xxxxxxy[k] * cd_y[k] + g_z_0_0_xxxxxxyy[k];

            g_z_0_y_xxxxxxz[k] = -g_z_0_0_xxxxxxz[k] * cd_y[k] + g_z_0_0_xxxxxxyz[k];

            g_z_0_y_xxxxxyy[k] = -g_z_0_0_xxxxxyy[k] * cd_y[k] + g_z_0_0_xxxxxyyy[k];

            g_z_0_y_xxxxxyz[k] = -g_z_0_0_xxxxxyz[k] * cd_y[k] + g_z_0_0_xxxxxyyz[k];

            g_z_0_y_xxxxxzz[k] = -g_z_0_0_xxxxxzz[k] * cd_y[k] + g_z_0_0_xxxxxyzz[k];

            g_z_0_y_xxxxyyy[k] = -g_z_0_0_xxxxyyy[k] * cd_y[k] + g_z_0_0_xxxxyyyy[k];

            g_z_0_y_xxxxyyz[k] = -g_z_0_0_xxxxyyz[k] * cd_y[k] + g_z_0_0_xxxxyyyz[k];

            g_z_0_y_xxxxyzz[k] = -g_z_0_0_xxxxyzz[k] * cd_y[k] + g_z_0_0_xxxxyyzz[k];

            g_z_0_y_xxxxzzz[k] = -g_z_0_0_xxxxzzz[k] * cd_y[k] + g_z_0_0_xxxxyzzz[k];

            g_z_0_y_xxxyyyy[k] = -g_z_0_0_xxxyyyy[k] * cd_y[k] + g_z_0_0_xxxyyyyy[k];

            g_z_0_y_xxxyyyz[k] = -g_z_0_0_xxxyyyz[k] * cd_y[k] + g_z_0_0_xxxyyyyz[k];

            g_z_0_y_xxxyyzz[k] = -g_z_0_0_xxxyyzz[k] * cd_y[k] + g_z_0_0_xxxyyyzz[k];

            g_z_0_y_xxxyzzz[k] = -g_z_0_0_xxxyzzz[k] * cd_y[k] + g_z_0_0_xxxyyzzz[k];

            g_z_0_y_xxxzzzz[k] = -g_z_0_0_xxxzzzz[k] * cd_y[k] + g_z_0_0_xxxyzzzz[k];

            g_z_0_y_xxyyyyy[k] = -g_z_0_0_xxyyyyy[k] * cd_y[k] + g_z_0_0_xxyyyyyy[k];

            g_z_0_y_xxyyyyz[k] = -g_z_0_0_xxyyyyz[k] * cd_y[k] + g_z_0_0_xxyyyyyz[k];

            g_z_0_y_xxyyyzz[k] = -g_z_0_0_xxyyyzz[k] * cd_y[k] + g_z_0_0_xxyyyyzz[k];

            g_z_0_y_xxyyzzz[k] = -g_z_0_0_xxyyzzz[k] * cd_y[k] + g_z_0_0_xxyyyzzz[k];

            g_z_0_y_xxyzzzz[k] = -g_z_0_0_xxyzzzz[k] * cd_y[k] + g_z_0_0_xxyyzzzz[k];

            g_z_0_y_xxzzzzz[k] = -g_z_0_0_xxzzzzz[k] * cd_y[k] + g_z_0_0_xxyzzzzz[k];

            g_z_0_y_xyyyyyy[k] = -g_z_0_0_xyyyyyy[k] * cd_y[k] + g_z_0_0_xyyyyyyy[k];

            g_z_0_y_xyyyyyz[k] = -g_z_0_0_xyyyyyz[k] * cd_y[k] + g_z_0_0_xyyyyyyz[k];

            g_z_0_y_xyyyyzz[k] = -g_z_0_0_xyyyyzz[k] * cd_y[k] + g_z_0_0_xyyyyyzz[k];

            g_z_0_y_xyyyzzz[k] = -g_z_0_0_xyyyzzz[k] * cd_y[k] + g_z_0_0_xyyyyzzz[k];

            g_z_0_y_xyyzzzz[k] = -g_z_0_0_xyyzzzz[k] * cd_y[k] + g_z_0_0_xyyyzzzz[k];

            g_z_0_y_xyzzzzz[k] = -g_z_0_0_xyzzzzz[k] * cd_y[k] + g_z_0_0_xyyzzzzz[k];

            g_z_0_y_xzzzzzz[k] = -g_z_0_0_xzzzzzz[k] * cd_y[k] + g_z_0_0_xyzzzzzz[k];

            g_z_0_y_yyyyyyy[k] = -g_z_0_0_yyyyyyy[k] * cd_y[k] + g_z_0_0_yyyyyyyy[k];

            g_z_0_y_yyyyyyz[k] = -g_z_0_0_yyyyyyz[k] * cd_y[k] + g_z_0_0_yyyyyyyz[k];

            g_z_0_y_yyyyyzz[k] = -g_z_0_0_yyyyyzz[k] * cd_y[k] + g_z_0_0_yyyyyyzz[k];

            g_z_0_y_yyyyzzz[k] = -g_z_0_0_yyyyzzz[k] * cd_y[k] + g_z_0_0_yyyyyzzz[k];

            g_z_0_y_yyyzzzz[k] = -g_z_0_0_yyyzzzz[k] * cd_y[k] + g_z_0_0_yyyyzzzz[k];

            g_z_0_y_yyzzzzz[k] = -g_z_0_0_yyzzzzz[k] * cd_y[k] + g_z_0_0_yyyzzzzz[k];

            g_z_0_y_yzzzzzz[k] = -g_z_0_0_yzzzzzz[k] * cd_y[k] + g_z_0_0_yyzzzzzz[k];

            g_z_0_y_zzzzzzz[k] = -g_z_0_0_zzzzzzz[k] * cd_y[k] + g_z_0_0_yzzzzzzz[k];
        }

        /// Set up 72-108 components of targeted buffer : cbuffer.data(

        auto g_z_0_z_xxxxxxx = cbuffer.data(pk_geom_10_off + 216 * acomps  + 72);

        auto g_z_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 73);

        auto g_z_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 74);

        auto g_z_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 75);

        auto g_z_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 76);

        auto g_z_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 77);

        auto g_z_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 78);

        auto g_z_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 79);

        auto g_z_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 80);

        auto g_z_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 81);

        auto g_z_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 82);

        auto g_z_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 83);

        auto g_z_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 84);

        auto g_z_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 85);

        auto g_z_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 86);

        auto g_z_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 87);

        auto g_z_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 88);

        auto g_z_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 89);

        auto g_z_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 90);

        auto g_z_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 91);

        auto g_z_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 92);

        auto g_z_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 93);

        auto g_z_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 94);

        auto g_z_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 95);

        auto g_z_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 96);

        auto g_z_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 97);

        auto g_z_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 98);

        auto g_z_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 99);

        auto g_z_0_z_yyyyyyy = cbuffer.data(pk_geom_10_off + 216 * acomps  + 100);

        auto g_z_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 101);

        auto g_z_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 102);

        auto g_z_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 103);

        auto g_z_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 104);

        auto g_z_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 105);

        auto g_z_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 106);

        auto g_z_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 216 * acomps  + 107);

        #pragma omp simd aligned(cd_z, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz, g_z_0_0_xxxxxxx, g_z_0_0_xxxxxxxz, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxxyz, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxxzz, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyyz, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxyzz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxxzzz, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyyz, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyyzz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxyzzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxxzzzz, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyyz, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyyzz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyyzzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxyzzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxxzzzzz, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyyz, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyyzz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyyzzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyyzzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxyzzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xxzzzzzz, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyyz, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyyzz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyyzzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyyzzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyyzzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xyzzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_xzzzzzzz, g_z_0_0_yyyyyyy, g_z_0_0_yyyyyyyz, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyyzz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyyzzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyyzzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyyzzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yyzzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_yzzzzzzz, g_z_0_0_zzzzzzz, g_z_0_0_zzzzzzzz, g_z_0_z_xxxxxxx, g_z_0_z_xxxxxxy, g_z_0_z_xxxxxxz, g_z_0_z_xxxxxyy, g_z_0_z_xxxxxyz, g_z_0_z_xxxxxzz, g_z_0_z_xxxxyyy, g_z_0_z_xxxxyyz, g_z_0_z_xxxxyzz, g_z_0_z_xxxxzzz, g_z_0_z_xxxyyyy, g_z_0_z_xxxyyyz, g_z_0_z_xxxyyzz, g_z_0_z_xxxyzzz, g_z_0_z_xxxzzzz, g_z_0_z_xxyyyyy, g_z_0_z_xxyyyyz, g_z_0_z_xxyyyzz, g_z_0_z_xxyyzzz, g_z_0_z_xxyzzzz, g_z_0_z_xxzzzzz, g_z_0_z_xyyyyyy, g_z_0_z_xyyyyyz, g_z_0_z_xyyyyzz, g_z_0_z_xyyyzzz, g_z_0_z_xyyzzzz, g_z_0_z_xyzzzzz, g_z_0_z_xzzzzzz, g_z_0_z_yyyyyyy, g_z_0_z_yyyyyyz, g_z_0_z_yyyyyzz, g_z_0_z_yyyyzzz, g_z_0_z_yyyzzzz, g_z_0_z_yyzzzzz, g_z_0_z_yzzzzzz, g_z_0_z_zzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_z_xxxxxxx[k] = -g_0_xxxxxxx[k] - g_z_0_0_xxxxxxx[k] * cd_z[k] + g_z_0_0_xxxxxxxz[k];

            g_z_0_z_xxxxxxy[k] = -g_0_xxxxxxy[k] - g_z_0_0_xxxxxxy[k] * cd_z[k] + g_z_0_0_xxxxxxyz[k];

            g_z_0_z_xxxxxxz[k] = -g_0_xxxxxxz[k] - g_z_0_0_xxxxxxz[k] * cd_z[k] + g_z_0_0_xxxxxxzz[k];

            g_z_0_z_xxxxxyy[k] = -g_0_xxxxxyy[k] - g_z_0_0_xxxxxyy[k] * cd_z[k] + g_z_0_0_xxxxxyyz[k];

            g_z_0_z_xxxxxyz[k] = -g_0_xxxxxyz[k] - g_z_0_0_xxxxxyz[k] * cd_z[k] + g_z_0_0_xxxxxyzz[k];

            g_z_0_z_xxxxxzz[k] = -g_0_xxxxxzz[k] - g_z_0_0_xxxxxzz[k] * cd_z[k] + g_z_0_0_xxxxxzzz[k];

            g_z_0_z_xxxxyyy[k] = -g_0_xxxxyyy[k] - g_z_0_0_xxxxyyy[k] * cd_z[k] + g_z_0_0_xxxxyyyz[k];

            g_z_0_z_xxxxyyz[k] = -g_0_xxxxyyz[k] - g_z_0_0_xxxxyyz[k] * cd_z[k] + g_z_0_0_xxxxyyzz[k];

            g_z_0_z_xxxxyzz[k] = -g_0_xxxxyzz[k] - g_z_0_0_xxxxyzz[k] * cd_z[k] + g_z_0_0_xxxxyzzz[k];

            g_z_0_z_xxxxzzz[k] = -g_0_xxxxzzz[k] - g_z_0_0_xxxxzzz[k] * cd_z[k] + g_z_0_0_xxxxzzzz[k];

            g_z_0_z_xxxyyyy[k] = -g_0_xxxyyyy[k] - g_z_0_0_xxxyyyy[k] * cd_z[k] + g_z_0_0_xxxyyyyz[k];

            g_z_0_z_xxxyyyz[k] = -g_0_xxxyyyz[k] - g_z_0_0_xxxyyyz[k] * cd_z[k] + g_z_0_0_xxxyyyzz[k];

            g_z_0_z_xxxyyzz[k] = -g_0_xxxyyzz[k] - g_z_0_0_xxxyyzz[k] * cd_z[k] + g_z_0_0_xxxyyzzz[k];

            g_z_0_z_xxxyzzz[k] = -g_0_xxxyzzz[k] - g_z_0_0_xxxyzzz[k] * cd_z[k] + g_z_0_0_xxxyzzzz[k];

            g_z_0_z_xxxzzzz[k] = -g_0_xxxzzzz[k] - g_z_0_0_xxxzzzz[k] * cd_z[k] + g_z_0_0_xxxzzzzz[k];

            g_z_0_z_xxyyyyy[k] = -g_0_xxyyyyy[k] - g_z_0_0_xxyyyyy[k] * cd_z[k] + g_z_0_0_xxyyyyyz[k];

            g_z_0_z_xxyyyyz[k] = -g_0_xxyyyyz[k] - g_z_0_0_xxyyyyz[k] * cd_z[k] + g_z_0_0_xxyyyyzz[k];

            g_z_0_z_xxyyyzz[k] = -g_0_xxyyyzz[k] - g_z_0_0_xxyyyzz[k] * cd_z[k] + g_z_0_0_xxyyyzzz[k];

            g_z_0_z_xxyyzzz[k] = -g_0_xxyyzzz[k] - g_z_0_0_xxyyzzz[k] * cd_z[k] + g_z_0_0_xxyyzzzz[k];

            g_z_0_z_xxyzzzz[k] = -g_0_xxyzzzz[k] - g_z_0_0_xxyzzzz[k] * cd_z[k] + g_z_0_0_xxyzzzzz[k];

            g_z_0_z_xxzzzzz[k] = -g_0_xxzzzzz[k] - g_z_0_0_xxzzzzz[k] * cd_z[k] + g_z_0_0_xxzzzzzz[k];

            g_z_0_z_xyyyyyy[k] = -g_0_xyyyyyy[k] - g_z_0_0_xyyyyyy[k] * cd_z[k] + g_z_0_0_xyyyyyyz[k];

            g_z_0_z_xyyyyyz[k] = -g_0_xyyyyyz[k] - g_z_0_0_xyyyyyz[k] * cd_z[k] + g_z_0_0_xyyyyyzz[k];

            g_z_0_z_xyyyyzz[k] = -g_0_xyyyyzz[k] - g_z_0_0_xyyyyzz[k] * cd_z[k] + g_z_0_0_xyyyyzzz[k];

            g_z_0_z_xyyyzzz[k] = -g_0_xyyyzzz[k] - g_z_0_0_xyyyzzz[k] * cd_z[k] + g_z_0_0_xyyyzzzz[k];

            g_z_0_z_xyyzzzz[k] = -g_0_xyyzzzz[k] - g_z_0_0_xyyzzzz[k] * cd_z[k] + g_z_0_0_xyyzzzzz[k];

            g_z_0_z_xyzzzzz[k] = -g_0_xyzzzzz[k] - g_z_0_0_xyzzzzz[k] * cd_z[k] + g_z_0_0_xyzzzzzz[k];

            g_z_0_z_xzzzzzz[k] = -g_0_xzzzzzz[k] - g_z_0_0_xzzzzzz[k] * cd_z[k] + g_z_0_0_xzzzzzzz[k];

            g_z_0_z_yyyyyyy[k] = -g_0_yyyyyyy[k] - g_z_0_0_yyyyyyy[k] * cd_z[k] + g_z_0_0_yyyyyyyz[k];

            g_z_0_z_yyyyyyz[k] = -g_0_yyyyyyz[k] - g_z_0_0_yyyyyyz[k] * cd_z[k] + g_z_0_0_yyyyyyzz[k];

            g_z_0_z_yyyyyzz[k] = -g_0_yyyyyzz[k] - g_z_0_0_yyyyyzz[k] * cd_z[k] + g_z_0_0_yyyyyzzz[k];

            g_z_0_z_yyyyzzz[k] = -g_0_yyyyzzz[k] - g_z_0_0_yyyyzzz[k] * cd_z[k] + g_z_0_0_yyyyzzzz[k];

            g_z_0_z_yyyzzzz[k] = -g_0_yyyzzzz[k] - g_z_0_0_yyyzzzz[k] * cd_z[k] + g_z_0_0_yyyzzzzz[k];

            g_z_0_z_yyzzzzz[k] = -g_0_yyzzzzz[k] - g_z_0_0_yyzzzzz[k] * cd_z[k] + g_z_0_0_yyzzzzzz[k];

            g_z_0_z_yzzzzzz[k] = -g_0_yzzzzzz[k] - g_z_0_0_yzzzzzz[k] * cd_z[k] + g_z_0_0_yzzzzzzz[k];

            g_z_0_z_zzzzzzz[k] = -g_0_zzzzzzz[k] - g_z_0_0_zzzzzzz[k] * cd_z[k] + g_z_0_0_zzzzzzzz[k];
        }
    }
}

} // t3ceri namespace


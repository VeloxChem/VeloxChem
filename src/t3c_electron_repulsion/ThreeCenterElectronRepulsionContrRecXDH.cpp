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

#include "ThreeCenterElectronRepulsionContrRecXDH.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xdh(CSimdArray<double>& cbuffer,
                                const size_t idx_xdh,
                                const size_t idx_xph,
                                const size_t idx_xpi,
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

        /// Set up components of auxilary buffer : SPI

        const auto pi_off = idx_xpi + i * 84;

        auto g_x_xxxxxx = cbuffer.data(pi_off + 0);

        auto g_x_xxxxxy = cbuffer.data(pi_off + 1);

        auto g_x_xxxxxz = cbuffer.data(pi_off + 2);

        auto g_x_xxxxyy = cbuffer.data(pi_off + 3);

        auto g_x_xxxxyz = cbuffer.data(pi_off + 4);

        auto g_x_xxxxzz = cbuffer.data(pi_off + 5);

        auto g_x_xxxyyy = cbuffer.data(pi_off + 6);

        auto g_x_xxxyyz = cbuffer.data(pi_off + 7);

        auto g_x_xxxyzz = cbuffer.data(pi_off + 8);

        auto g_x_xxxzzz = cbuffer.data(pi_off + 9);

        auto g_x_xxyyyy = cbuffer.data(pi_off + 10);

        auto g_x_xxyyyz = cbuffer.data(pi_off + 11);

        auto g_x_xxyyzz = cbuffer.data(pi_off + 12);

        auto g_x_xxyzzz = cbuffer.data(pi_off + 13);

        auto g_x_xxzzzz = cbuffer.data(pi_off + 14);

        auto g_x_xyyyyy = cbuffer.data(pi_off + 15);

        auto g_x_xyyyyz = cbuffer.data(pi_off + 16);

        auto g_x_xyyyzz = cbuffer.data(pi_off + 17);

        auto g_x_xyyzzz = cbuffer.data(pi_off + 18);

        auto g_x_xyzzzz = cbuffer.data(pi_off + 19);

        auto g_x_xzzzzz = cbuffer.data(pi_off + 20);

        auto g_y_xxxxxx = cbuffer.data(pi_off + 28);

        auto g_y_xxxxxy = cbuffer.data(pi_off + 29);

        auto g_y_xxxxxz = cbuffer.data(pi_off + 30);

        auto g_y_xxxxyy = cbuffer.data(pi_off + 31);

        auto g_y_xxxxyz = cbuffer.data(pi_off + 32);

        auto g_y_xxxxzz = cbuffer.data(pi_off + 33);

        auto g_y_xxxyyy = cbuffer.data(pi_off + 34);

        auto g_y_xxxyyz = cbuffer.data(pi_off + 35);

        auto g_y_xxxyzz = cbuffer.data(pi_off + 36);

        auto g_y_xxxzzz = cbuffer.data(pi_off + 37);

        auto g_y_xxyyyy = cbuffer.data(pi_off + 38);

        auto g_y_xxyyyz = cbuffer.data(pi_off + 39);

        auto g_y_xxyyzz = cbuffer.data(pi_off + 40);

        auto g_y_xxyzzz = cbuffer.data(pi_off + 41);

        auto g_y_xxzzzz = cbuffer.data(pi_off + 42);

        auto g_y_xyyyyy = cbuffer.data(pi_off + 43);

        auto g_y_xyyyyz = cbuffer.data(pi_off + 44);

        auto g_y_xyyyzz = cbuffer.data(pi_off + 45);

        auto g_y_xyyzzz = cbuffer.data(pi_off + 46);

        auto g_y_xyzzzz = cbuffer.data(pi_off + 47);

        auto g_y_xzzzzz = cbuffer.data(pi_off + 48);

        auto g_y_yyyyyy = cbuffer.data(pi_off + 49);

        auto g_y_yyyyyz = cbuffer.data(pi_off + 50);

        auto g_y_yyyyzz = cbuffer.data(pi_off + 51);

        auto g_y_yyyzzz = cbuffer.data(pi_off + 52);

        auto g_y_yyzzzz = cbuffer.data(pi_off + 53);

        auto g_y_yzzzzz = cbuffer.data(pi_off + 54);

        auto g_z_xxxxxx = cbuffer.data(pi_off + 56);

        auto g_z_xxxxxy = cbuffer.data(pi_off + 57);

        auto g_z_xxxxxz = cbuffer.data(pi_off + 58);

        auto g_z_xxxxyy = cbuffer.data(pi_off + 59);

        auto g_z_xxxxyz = cbuffer.data(pi_off + 60);

        auto g_z_xxxxzz = cbuffer.data(pi_off + 61);

        auto g_z_xxxyyy = cbuffer.data(pi_off + 62);

        auto g_z_xxxyyz = cbuffer.data(pi_off + 63);

        auto g_z_xxxyzz = cbuffer.data(pi_off + 64);

        auto g_z_xxxzzz = cbuffer.data(pi_off + 65);

        auto g_z_xxyyyy = cbuffer.data(pi_off + 66);

        auto g_z_xxyyyz = cbuffer.data(pi_off + 67);

        auto g_z_xxyyzz = cbuffer.data(pi_off + 68);

        auto g_z_xxyzzz = cbuffer.data(pi_off + 69);

        auto g_z_xxzzzz = cbuffer.data(pi_off + 70);

        auto g_z_xyyyyy = cbuffer.data(pi_off + 71);

        auto g_z_xyyyyz = cbuffer.data(pi_off + 72);

        auto g_z_xyyyzz = cbuffer.data(pi_off + 73);

        auto g_z_xyyzzz = cbuffer.data(pi_off + 74);

        auto g_z_xyzzzz = cbuffer.data(pi_off + 75);

        auto g_z_xzzzzz = cbuffer.data(pi_off + 76);

        auto g_z_yyyyyy = cbuffer.data(pi_off + 77);

        auto g_z_yyyyyz = cbuffer.data(pi_off + 78);

        auto g_z_yyyyzz = cbuffer.data(pi_off + 79);

        auto g_z_yyyzzz = cbuffer.data(pi_off + 80);

        auto g_z_yyzzzz = cbuffer.data(pi_off + 81);

        auto g_z_yzzzzz = cbuffer.data(pi_off + 82);

        auto g_z_zzzzzz = cbuffer.data(pi_off + 83);

        /// set up bra offset for contr_buffer_xdh

        const auto dh_off = idx_xdh + i * 126;

        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_xx_xxxxx = cbuffer.data(dh_off + 0);

        auto g_xx_xxxxy = cbuffer.data(dh_off + 1);

        auto g_xx_xxxxz = cbuffer.data(dh_off + 2);

        auto g_xx_xxxyy = cbuffer.data(dh_off + 3);

        auto g_xx_xxxyz = cbuffer.data(dh_off + 4);

        auto g_xx_xxxzz = cbuffer.data(dh_off + 5);

        auto g_xx_xxyyy = cbuffer.data(dh_off + 6);

        auto g_xx_xxyyz = cbuffer.data(dh_off + 7);

        auto g_xx_xxyzz = cbuffer.data(dh_off + 8);

        auto g_xx_xxzzz = cbuffer.data(dh_off + 9);

        auto g_xx_xyyyy = cbuffer.data(dh_off + 10);

        auto g_xx_xyyyz = cbuffer.data(dh_off + 11);

        auto g_xx_xyyzz = cbuffer.data(dh_off + 12);

        auto g_xx_xyzzz = cbuffer.data(dh_off + 13);

        auto g_xx_xzzzz = cbuffer.data(dh_off + 14);

        auto g_xx_yyyyy = cbuffer.data(dh_off + 15);

        auto g_xx_yyyyz = cbuffer.data(dh_off + 16);

        auto g_xx_yyyzz = cbuffer.data(dh_off + 17);

        auto g_xx_yyzzz = cbuffer.data(dh_off + 18);

        auto g_xx_yzzzz = cbuffer.data(dh_off + 19);

        auto g_xx_zzzzz = cbuffer.data(dh_off + 20);

        #pragma omp simd aligned(cd_x, g_x_xxxxx, g_x_xxxxxx, g_x_xxxxxy, g_x_xxxxxz, g_x_xxxxy, g_x_xxxxyy, g_x_xxxxyz, g_x_xxxxz, g_x_xxxxzz, g_x_xxxyy, g_x_xxxyyy, g_x_xxxyyz, g_x_xxxyz, g_x_xxxyzz, g_x_xxxzz, g_x_xxxzzz, g_x_xxyyy, g_x_xxyyyy, g_x_xxyyyz, g_x_xxyyz, g_x_xxyyzz, g_x_xxyzz, g_x_xxyzzz, g_x_xxzzz, g_x_xxzzzz, g_x_xyyyy, g_x_xyyyyy, g_x_xyyyyz, g_x_xyyyz, g_x_xyyyzz, g_x_xyyzz, g_x_xyyzzz, g_x_xyzzz, g_x_xyzzzz, g_x_xzzzz, g_x_xzzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz, g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz, g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxzz, g_xx_xxyyy, g_xx_xxyyz, g_xx_xxyzz, g_xx_xxzzz, g_xx_xyyyy, g_xx_xyyyz, g_xx_xyyzz, g_xx_xyzzz, g_xx_xzzzz, g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz, g_xx_yyzzz, g_xx_yzzzz, g_xx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xx_xxxxx[k] = -g_x_xxxxx[k] * cd_x[k] + g_x_xxxxxx[k];

            g_xx_xxxxy[k] = -g_x_xxxxy[k] * cd_x[k] + g_x_xxxxxy[k];

            g_xx_xxxxz[k] = -g_x_xxxxz[k] * cd_x[k] + g_x_xxxxxz[k];

            g_xx_xxxyy[k] = -g_x_xxxyy[k] * cd_x[k] + g_x_xxxxyy[k];

            g_xx_xxxyz[k] = -g_x_xxxyz[k] * cd_x[k] + g_x_xxxxyz[k];

            g_xx_xxxzz[k] = -g_x_xxxzz[k] * cd_x[k] + g_x_xxxxzz[k];

            g_xx_xxyyy[k] = -g_x_xxyyy[k] * cd_x[k] + g_x_xxxyyy[k];

            g_xx_xxyyz[k] = -g_x_xxyyz[k] * cd_x[k] + g_x_xxxyyz[k];

            g_xx_xxyzz[k] = -g_x_xxyzz[k] * cd_x[k] + g_x_xxxyzz[k];

            g_xx_xxzzz[k] = -g_x_xxzzz[k] * cd_x[k] + g_x_xxxzzz[k];

            g_xx_xyyyy[k] = -g_x_xyyyy[k] * cd_x[k] + g_x_xxyyyy[k];

            g_xx_xyyyz[k] = -g_x_xyyyz[k] * cd_x[k] + g_x_xxyyyz[k];

            g_xx_xyyzz[k] = -g_x_xyyzz[k] * cd_x[k] + g_x_xxyyzz[k];

            g_xx_xyzzz[k] = -g_x_xyzzz[k] * cd_x[k] + g_x_xxyzzz[k];

            g_xx_xzzzz[k] = -g_x_xzzzz[k] * cd_x[k] + g_x_xxzzzz[k];

            g_xx_yyyyy[k] = -g_x_yyyyy[k] * cd_x[k] + g_x_xyyyyy[k];

            g_xx_yyyyz[k] = -g_x_yyyyz[k] * cd_x[k] + g_x_xyyyyz[k];

            g_xx_yyyzz[k] = -g_x_yyyzz[k] * cd_x[k] + g_x_xyyyzz[k];

            g_xx_yyzzz[k] = -g_x_yyzzz[k] * cd_x[k] + g_x_xyyzzz[k];

            g_xx_yzzzz[k] = -g_x_yzzzz[k] * cd_x[k] + g_x_xyzzzz[k];

            g_xx_zzzzz[k] = -g_x_zzzzz[k] * cd_x[k] + g_x_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_xy_xxxxx = cbuffer.data(dh_off + 21);

        auto g_xy_xxxxy = cbuffer.data(dh_off + 22);

        auto g_xy_xxxxz = cbuffer.data(dh_off + 23);

        auto g_xy_xxxyy = cbuffer.data(dh_off + 24);

        auto g_xy_xxxyz = cbuffer.data(dh_off + 25);

        auto g_xy_xxxzz = cbuffer.data(dh_off + 26);

        auto g_xy_xxyyy = cbuffer.data(dh_off + 27);

        auto g_xy_xxyyz = cbuffer.data(dh_off + 28);

        auto g_xy_xxyzz = cbuffer.data(dh_off + 29);

        auto g_xy_xxzzz = cbuffer.data(dh_off + 30);

        auto g_xy_xyyyy = cbuffer.data(dh_off + 31);

        auto g_xy_xyyyz = cbuffer.data(dh_off + 32);

        auto g_xy_xyyzz = cbuffer.data(dh_off + 33);

        auto g_xy_xyzzz = cbuffer.data(dh_off + 34);

        auto g_xy_xzzzz = cbuffer.data(dh_off + 35);

        auto g_xy_yyyyy = cbuffer.data(dh_off + 36);

        auto g_xy_yyyyz = cbuffer.data(dh_off + 37);

        auto g_xy_yyyzz = cbuffer.data(dh_off + 38);

        auto g_xy_yyzzz = cbuffer.data(dh_off + 39);

        auto g_xy_yzzzz = cbuffer.data(dh_off + 40);

        auto g_xy_zzzzz = cbuffer.data(dh_off + 41);

        #pragma omp simd aligned(cd_x, g_xy_xxxxx, g_xy_xxxxy, g_xy_xxxxz, g_xy_xxxyy, g_xy_xxxyz, g_xy_xxxzz, g_xy_xxyyy, g_xy_xxyyz, g_xy_xxyzz, g_xy_xxzzz, g_xy_xyyyy, g_xy_xyyyz, g_xy_xyyzz, g_xy_xyzzz, g_xy_xzzzz, g_xy_yyyyy, g_xy_yyyyz, g_xy_yyyzz, g_xy_yyzzz, g_xy_yzzzz, g_xy_zzzzz, g_y_xxxxx, g_y_xxxxxx, g_y_xxxxxy, g_y_xxxxxz, g_y_xxxxy, g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxz, g_y_xxxxzz, g_y_xxxyy, g_y_xxxyyy, g_y_xxxyyz, g_y_xxxyz, g_y_xxxyzz, g_y_xxxzz, g_y_xxxzzz, g_y_xxyyy, g_y_xxyyyy, g_y_xxyyyz, g_y_xxyyz, g_y_xxyyzz, g_y_xxyzz, g_y_xxyzzz, g_y_xxzzz, g_y_xxzzzz, g_y_xyyyy, g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyz, g_y_xyyyzz, g_y_xyyzz, g_y_xyyzzz, g_y_xyzzz, g_y_xyzzzz, g_y_xzzzz, g_y_xzzzzz, g_y_yyyyy, g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xy_xxxxx[k] = -g_y_xxxxx[k] * cd_x[k] + g_y_xxxxxx[k];

            g_xy_xxxxy[k] = -g_y_xxxxy[k] * cd_x[k] + g_y_xxxxxy[k];

            g_xy_xxxxz[k] = -g_y_xxxxz[k] * cd_x[k] + g_y_xxxxxz[k];

            g_xy_xxxyy[k] = -g_y_xxxyy[k] * cd_x[k] + g_y_xxxxyy[k];

            g_xy_xxxyz[k] = -g_y_xxxyz[k] * cd_x[k] + g_y_xxxxyz[k];

            g_xy_xxxzz[k] = -g_y_xxxzz[k] * cd_x[k] + g_y_xxxxzz[k];

            g_xy_xxyyy[k] = -g_y_xxyyy[k] * cd_x[k] + g_y_xxxyyy[k];

            g_xy_xxyyz[k] = -g_y_xxyyz[k] * cd_x[k] + g_y_xxxyyz[k];

            g_xy_xxyzz[k] = -g_y_xxyzz[k] * cd_x[k] + g_y_xxxyzz[k];

            g_xy_xxzzz[k] = -g_y_xxzzz[k] * cd_x[k] + g_y_xxxzzz[k];

            g_xy_xyyyy[k] = -g_y_xyyyy[k] * cd_x[k] + g_y_xxyyyy[k];

            g_xy_xyyyz[k] = -g_y_xyyyz[k] * cd_x[k] + g_y_xxyyyz[k];

            g_xy_xyyzz[k] = -g_y_xyyzz[k] * cd_x[k] + g_y_xxyyzz[k];

            g_xy_xyzzz[k] = -g_y_xyzzz[k] * cd_x[k] + g_y_xxyzzz[k];

            g_xy_xzzzz[k] = -g_y_xzzzz[k] * cd_x[k] + g_y_xxzzzz[k];

            g_xy_yyyyy[k] = -g_y_yyyyy[k] * cd_x[k] + g_y_xyyyyy[k];

            g_xy_yyyyz[k] = -g_y_yyyyz[k] * cd_x[k] + g_y_xyyyyz[k];

            g_xy_yyyzz[k] = -g_y_yyyzz[k] * cd_x[k] + g_y_xyyyzz[k];

            g_xy_yyzzz[k] = -g_y_yyzzz[k] * cd_x[k] + g_y_xyyzzz[k];

            g_xy_yzzzz[k] = -g_y_yzzzz[k] * cd_x[k] + g_y_xyzzzz[k];

            g_xy_zzzzz[k] = -g_y_zzzzz[k] * cd_x[k] + g_y_xzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_xz_xxxxx = cbuffer.data(dh_off + 42);

        auto g_xz_xxxxy = cbuffer.data(dh_off + 43);

        auto g_xz_xxxxz = cbuffer.data(dh_off + 44);

        auto g_xz_xxxyy = cbuffer.data(dh_off + 45);

        auto g_xz_xxxyz = cbuffer.data(dh_off + 46);

        auto g_xz_xxxzz = cbuffer.data(dh_off + 47);

        auto g_xz_xxyyy = cbuffer.data(dh_off + 48);

        auto g_xz_xxyyz = cbuffer.data(dh_off + 49);

        auto g_xz_xxyzz = cbuffer.data(dh_off + 50);

        auto g_xz_xxzzz = cbuffer.data(dh_off + 51);

        auto g_xz_xyyyy = cbuffer.data(dh_off + 52);

        auto g_xz_xyyyz = cbuffer.data(dh_off + 53);

        auto g_xz_xyyzz = cbuffer.data(dh_off + 54);

        auto g_xz_xyzzz = cbuffer.data(dh_off + 55);

        auto g_xz_xzzzz = cbuffer.data(dh_off + 56);

        auto g_xz_yyyyy = cbuffer.data(dh_off + 57);

        auto g_xz_yyyyz = cbuffer.data(dh_off + 58);

        auto g_xz_yyyzz = cbuffer.data(dh_off + 59);

        auto g_xz_yyzzz = cbuffer.data(dh_off + 60);

        auto g_xz_yzzzz = cbuffer.data(dh_off + 61);

        auto g_xz_zzzzz = cbuffer.data(dh_off + 62);

        #pragma omp simd aligned(cd_x, g_xz_xxxxx, g_xz_xxxxy, g_xz_xxxxz, g_xz_xxxyy, g_xz_xxxyz, g_xz_xxxzz, g_xz_xxyyy, g_xz_xxyyz, g_xz_xxyzz, g_xz_xxzzz, g_xz_xyyyy, g_xz_xyyyz, g_xz_xyyzz, g_xz_xyzzz, g_xz_xzzzz, g_xz_yyyyy, g_xz_yyyyz, g_xz_yyyzz, g_xz_yyzzz, g_xz_yzzzz, g_xz_zzzzz, g_z_xxxxx, g_z_xxxxxx, g_z_xxxxxy, g_z_xxxxxz, g_z_xxxxy, g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxz, g_z_xxxxzz, g_z_xxxyy, g_z_xxxyyy, g_z_xxxyyz, g_z_xxxyz, g_z_xxxyzz, g_z_xxxzz, g_z_xxxzzz, g_z_xxyyy, g_z_xxyyyy, g_z_xxyyyz, g_z_xxyyz, g_z_xxyyzz, g_z_xxyzz, g_z_xxyzzz, g_z_xxzzz, g_z_xxzzzz, g_z_xyyyy, g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyz, g_z_xyyyzz, g_z_xyyzz, g_z_xyyzzz, g_z_xyzzz, g_z_xyzzzz, g_z_xzzzz, g_z_xzzzzz, g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz, g_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xz_xxxxx[k] = -g_z_xxxxx[k] * cd_x[k] + g_z_xxxxxx[k];

            g_xz_xxxxy[k] = -g_z_xxxxy[k] * cd_x[k] + g_z_xxxxxy[k];

            g_xz_xxxxz[k] = -g_z_xxxxz[k] * cd_x[k] + g_z_xxxxxz[k];

            g_xz_xxxyy[k] = -g_z_xxxyy[k] * cd_x[k] + g_z_xxxxyy[k];

            g_xz_xxxyz[k] = -g_z_xxxyz[k] * cd_x[k] + g_z_xxxxyz[k];

            g_xz_xxxzz[k] = -g_z_xxxzz[k] * cd_x[k] + g_z_xxxxzz[k];

            g_xz_xxyyy[k] = -g_z_xxyyy[k] * cd_x[k] + g_z_xxxyyy[k];

            g_xz_xxyyz[k] = -g_z_xxyyz[k] * cd_x[k] + g_z_xxxyyz[k];

            g_xz_xxyzz[k] = -g_z_xxyzz[k] * cd_x[k] + g_z_xxxyzz[k];

            g_xz_xxzzz[k] = -g_z_xxzzz[k] * cd_x[k] + g_z_xxxzzz[k];

            g_xz_xyyyy[k] = -g_z_xyyyy[k] * cd_x[k] + g_z_xxyyyy[k];

            g_xz_xyyyz[k] = -g_z_xyyyz[k] * cd_x[k] + g_z_xxyyyz[k];

            g_xz_xyyzz[k] = -g_z_xyyzz[k] * cd_x[k] + g_z_xxyyzz[k];

            g_xz_xyzzz[k] = -g_z_xyzzz[k] * cd_x[k] + g_z_xxyzzz[k];

            g_xz_xzzzz[k] = -g_z_xzzzz[k] * cd_x[k] + g_z_xxzzzz[k];

            g_xz_yyyyy[k] = -g_z_yyyyy[k] * cd_x[k] + g_z_xyyyyy[k];

            g_xz_yyyyz[k] = -g_z_yyyyz[k] * cd_x[k] + g_z_xyyyyz[k];

            g_xz_yyyzz[k] = -g_z_yyyzz[k] * cd_x[k] + g_z_xyyyzz[k];

            g_xz_yyzzz[k] = -g_z_yyzzz[k] * cd_x[k] + g_z_xyyzzz[k];

            g_xz_yzzzz[k] = -g_z_yzzzz[k] * cd_x[k] + g_z_xyzzzz[k];

            g_xz_zzzzz[k] = -g_z_zzzzz[k] * cd_x[k] + g_z_xzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_yy_xxxxx = cbuffer.data(dh_off + 63);

        auto g_yy_xxxxy = cbuffer.data(dh_off + 64);

        auto g_yy_xxxxz = cbuffer.data(dh_off + 65);

        auto g_yy_xxxyy = cbuffer.data(dh_off + 66);

        auto g_yy_xxxyz = cbuffer.data(dh_off + 67);

        auto g_yy_xxxzz = cbuffer.data(dh_off + 68);

        auto g_yy_xxyyy = cbuffer.data(dh_off + 69);

        auto g_yy_xxyyz = cbuffer.data(dh_off + 70);

        auto g_yy_xxyzz = cbuffer.data(dh_off + 71);

        auto g_yy_xxzzz = cbuffer.data(dh_off + 72);

        auto g_yy_xyyyy = cbuffer.data(dh_off + 73);

        auto g_yy_xyyyz = cbuffer.data(dh_off + 74);

        auto g_yy_xyyzz = cbuffer.data(dh_off + 75);

        auto g_yy_xyzzz = cbuffer.data(dh_off + 76);

        auto g_yy_xzzzz = cbuffer.data(dh_off + 77);

        auto g_yy_yyyyy = cbuffer.data(dh_off + 78);

        auto g_yy_yyyyz = cbuffer.data(dh_off + 79);

        auto g_yy_yyyzz = cbuffer.data(dh_off + 80);

        auto g_yy_yyzzz = cbuffer.data(dh_off + 81);

        auto g_yy_yzzzz = cbuffer.data(dh_off + 82);

        auto g_yy_zzzzz = cbuffer.data(dh_off + 83);

        #pragma omp simd aligned(cd_y, g_y_xxxxx, g_y_xxxxxy, g_y_xxxxy, g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxz, g_y_xxxyy, g_y_xxxyyy, g_y_xxxyyz, g_y_xxxyz, g_y_xxxyzz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyyy, g_y_xxyyyz, g_y_xxyyz, g_y_xxyyzz, g_y_xxyzz, g_y_xxyzzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyz, g_y_xyyyzz, g_y_xyyzz, g_y_xyyzzz, g_y_xyzzz, g_y_xyzzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyz, g_y_yyyyzz, g_y_yyyzz, g_y_yyyzzz, g_y_yyzzz, g_y_yyzzzz, g_y_yzzzz, g_y_yzzzzz, g_y_zzzzz, g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz, g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyzz, g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyz, g_yy_xyyzz, g_yy_xyzzz, g_yy_xzzzz, g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz, g_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yy_xxxxx[k] = -g_y_xxxxx[k] * cd_y[k] + g_y_xxxxxy[k];

            g_yy_xxxxy[k] = -g_y_xxxxy[k] * cd_y[k] + g_y_xxxxyy[k];

            g_yy_xxxxz[k] = -g_y_xxxxz[k] * cd_y[k] + g_y_xxxxyz[k];

            g_yy_xxxyy[k] = -g_y_xxxyy[k] * cd_y[k] + g_y_xxxyyy[k];

            g_yy_xxxyz[k] = -g_y_xxxyz[k] * cd_y[k] + g_y_xxxyyz[k];

            g_yy_xxxzz[k] = -g_y_xxxzz[k] * cd_y[k] + g_y_xxxyzz[k];

            g_yy_xxyyy[k] = -g_y_xxyyy[k] * cd_y[k] + g_y_xxyyyy[k];

            g_yy_xxyyz[k] = -g_y_xxyyz[k] * cd_y[k] + g_y_xxyyyz[k];

            g_yy_xxyzz[k] = -g_y_xxyzz[k] * cd_y[k] + g_y_xxyyzz[k];

            g_yy_xxzzz[k] = -g_y_xxzzz[k] * cd_y[k] + g_y_xxyzzz[k];

            g_yy_xyyyy[k] = -g_y_xyyyy[k] * cd_y[k] + g_y_xyyyyy[k];

            g_yy_xyyyz[k] = -g_y_xyyyz[k] * cd_y[k] + g_y_xyyyyz[k];

            g_yy_xyyzz[k] = -g_y_xyyzz[k] * cd_y[k] + g_y_xyyyzz[k];

            g_yy_xyzzz[k] = -g_y_xyzzz[k] * cd_y[k] + g_y_xyyzzz[k];

            g_yy_xzzzz[k] = -g_y_xzzzz[k] * cd_y[k] + g_y_xyzzzz[k];

            g_yy_yyyyy[k] = -g_y_yyyyy[k] * cd_y[k] + g_y_yyyyyy[k];

            g_yy_yyyyz[k] = -g_y_yyyyz[k] * cd_y[k] + g_y_yyyyyz[k];

            g_yy_yyyzz[k] = -g_y_yyyzz[k] * cd_y[k] + g_y_yyyyzz[k];

            g_yy_yyzzz[k] = -g_y_yyzzz[k] * cd_y[k] + g_y_yyyzzz[k];

            g_yy_yzzzz[k] = -g_y_yzzzz[k] * cd_y[k] + g_y_yyzzzz[k];

            g_yy_zzzzz[k] = -g_y_zzzzz[k] * cd_y[k] + g_y_yzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_yz_xxxxx = cbuffer.data(dh_off + 84);

        auto g_yz_xxxxy = cbuffer.data(dh_off + 85);

        auto g_yz_xxxxz = cbuffer.data(dh_off + 86);

        auto g_yz_xxxyy = cbuffer.data(dh_off + 87);

        auto g_yz_xxxyz = cbuffer.data(dh_off + 88);

        auto g_yz_xxxzz = cbuffer.data(dh_off + 89);

        auto g_yz_xxyyy = cbuffer.data(dh_off + 90);

        auto g_yz_xxyyz = cbuffer.data(dh_off + 91);

        auto g_yz_xxyzz = cbuffer.data(dh_off + 92);

        auto g_yz_xxzzz = cbuffer.data(dh_off + 93);

        auto g_yz_xyyyy = cbuffer.data(dh_off + 94);

        auto g_yz_xyyyz = cbuffer.data(dh_off + 95);

        auto g_yz_xyyzz = cbuffer.data(dh_off + 96);

        auto g_yz_xyzzz = cbuffer.data(dh_off + 97);

        auto g_yz_xzzzz = cbuffer.data(dh_off + 98);

        auto g_yz_yyyyy = cbuffer.data(dh_off + 99);

        auto g_yz_yyyyz = cbuffer.data(dh_off + 100);

        auto g_yz_yyyzz = cbuffer.data(dh_off + 101);

        auto g_yz_yyzzz = cbuffer.data(dh_off + 102);

        auto g_yz_yzzzz = cbuffer.data(dh_off + 103);

        auto g_yz_zzzzz = cbuffer.data(dh_off + 104);

        #pragma omp simd aligned(cd_y, g_yz_xxxxx, g_yz_xxxxy, g_yz_xxxxz, g_yz_xxxyy, g_yz_xxxyz, g_yz_xxxzz, g_yz_xxyyy, g_yz_xxyyz, g_yz_xxyzz, g_yz_xxzzz, g_yz_xyyyy, g_yz_xyyyz, g_yz_xyyzz, g_yz_xyzzz, g_yz_xzzzz, g_yz_yyyyy, g_yz_yyyyz, g_yz_yyyzz, g_yz_yyzzz, g_yz_yzzzz, g_yz_zzzzz, g_z_xxxxx, g_z_xxxxxy, g_z_xxxxy, g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxz, g_z_xxxyy, g_z_xxxyyy, g_z_xxxyyz, g_z_xxxyz, g_z_xxxyzz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyyy, g_z_xxyyyz, g_z_xxyyz, g_z_xxyyzz, g_z_xxyzz, g_z_xxyzzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyz, g_z_xyyyzz, g_z_xyyzz, g_z_xyyzzz, g_z_xyzzz, g_z_xyzzzz, g_z_xzzzz, g_z_yyyyy, g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyz, g_z_yyyyzz, g_z_yyyzz, g_z_yyyzzz, g_z_yyzzz, g_z_yyzzzz, g_z_yzzzz, g_z_yzzzzz, g_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yz_xxxxx[k] = -g_z_xxxxx[k] * cd_y[k] + g_z_xxxxxy[k];

            g_yz_xxxxy[k] = -g_z_xxxxy[k] * cd_y[k] + g_z_xxxxyy[k];

            g_yz_xxxxz[k] = -g_z_xxxxz[k] * cd_y[k] + g_z_xxxxyz[k];

            g_yz_xxxyy[k] = -g_z_xxxyy[k] * cd_y[k] + g_z_xxxyyy[k];

            g_yz_xxxyz[k] = -g_z_xxxyz[k] * cd_y[k] + g_z_xxxyyz[k];

            g_yz_xxxzz[k] = -g_z_xxxzz[k] * cd_y[k] + g_z_xxxyzz[k];

            g_yz_xxyyy[k] = -g_z_xxyyy[k] * cd_y[k] + g_z_xxyyyy[k];

            g_yz_xxyyz[k] = -g_z_xxyyz[k] * cd_y[k] + g_z_xxyyyz[k];

            g_yz_xxyzz[k] = -g_z_xxyzz[k] * cd_y[k] + g_z_xxyyzz[k];

            g_yz_xxzzz[k] = -g_z_xxzzz[k] * cd_y[k] + g_z_xxyzzz[k];

            g_yz_xyyyy[k] = -g_z_xyyyy[k] * cd_y[k] + g_z_xyyyyy[k];

            g_yz_xyyyz[k] = -g_z_xyyyz[k] * cd_y[k] + g_z_xyyyyz[k];

            g_yz_xyyzz[k] = -g_z_xyyzz[k] * cd_y[k] + g_z_xyyyzz[k];

            g_yz_xyzzz[k] = -g_z_xyzzz[k] * cd_y[k] + g_z_xyyzzz[k];

            g_yz_xzzzz[k] = -g_z_xzzzz[k] * cd_y[k] + g_z_xyzzzz[k];

            g_yz_yyyyy[k] = -g_z_yyyyy[k] * cd_y[k] + g_z_yyyyyy[k];

            g_yz_yyyyz[k] = -g_z_yyyyz[k] * cd_y[k] + g_z_yyyyyz[k];

            g_yz_yyyzz[k] = -g_z_yyyzz[k] * cd_y[k] + g_z_yyyyzz[k];

            g_yz_yyzzz[k] = -g_z_yyzzz[k] * cd_y[k] + g_z_yyyzzz[k];

            g_yz_yzzzz[k] = -g_z_yzzzz[k] * cd_y[k] + g_z_yyzzzz[k];

            g_yz_zzzzz[k] = -g_z_zzzzz[k] * cd_y[k] + g_z_yzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_zz_xxxxx = cbuffer.data(dh_off + 105);

        auto g_zz_xxxxy = cbuffer.data(dh_off + 106);

        auto g_zz_xxxxz = cbuffer.data(dh_off + 107);

        auto g_zz_xxxyy = cbuffer.data(dh_off + 108);

        auto g_zz_xxxyz = cbuffer.data(dh_off + 109);

        auto g_zz_xxxzz = cbuffer.data(dh_off + 110);

        auto g_zz_xxyyy = cbuffer.data(dh_off + 111);

        auto g_zz_xxyyz = cbuffer.data(dh_off + 112);

        auto g_zz_xxyzz = cbuffer.data(dh_off + 113);

        auto g_zz_xxzzz = cbuffer.data(dh_off + 114);

        auto g_zz_xyyyy = cbuffer.data(dh_off + 115);

        auto g_zz_xyyyz = cbuffer.data(dh_off + 116);

        auto g_zz_xyyzz = cbuffer.data(dh_off + 117);

        auto g_zz_xyzzz = cbuffer.data(dh_off + 118);

        auto g_zz_xzzzz = cbuffer.data(dh_off + 119);

        auto g_zz_yyyyy = cbuffer.data(dh_off + 120);

        auto g_zz_yyyyz = cbuffer.data(dh_off + 121);

        auto g_zz_yyyzz = cbuffer.data(dh_off + 122);

        auto g_zz_yyzzz = cbuffer.data(dh_off + 123);

        auto g_zz_yzzzz = cbuffer.data(dh_off + 124);

        auto g_zz_zzzzz = cbuffer.data(dh_off + 125);

        #pragma omp simd aligned(cd_z, g_z_xxxxx, g_z_xxxxxz, g_z_xxxxy, g_z_xxxxyz, g_z_xxxxz, g_z_xxxxzz, g_z_xxxyy, g_z_xxxyyz, g_z_xxxyz, g_z_xxxyzz, g_z_xxxzz, g_z_xxxzzz, g_z_xxyyy, g_z_xxyyyz, g_z_xxyyz, g_z_xxyyzz, g_z_xxyzz, g_z_xxyzzz, g_z_xxzzz, g_z_xxzzzz, g_z_xyyyy, g_z_xyyyyz, g_z_xyyyz, g_z_xyyyzz, g_z_xyyzz, g_z_xyyzzz, g_z_xyzzz, g_z_xyzzzz, g_z_xzzzz, g_z_xzzzzz, g_z_yyyyy, g_z_yyyyyz, g_z_yyyyz, g_z_yyyyzz, g_z_yyyzz, g_z_yyyzzz, g_z_yyzzz, g_z_yyzzzz, g_z_yzzzz, g_z_yzzzzz, g_z_zzzzz, g_z_zzzzzz, g_zz_xxxxx, g_zz_xxxxy, g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxzz, g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyzz, g_zz_xxzzz, g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyzz, g_zz_xyzzz, g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz, g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zz_xxxxx[k] = -g_z_xxxxx[k] * cd_z[k] + g_z_xxxxxz[k];

            g_zz_xxxxy[k] = -g_z_xxxxy[k] * cd_z[k] + g_z_xxxxyz[k];

            g_zz_xxxxz[k] = -g_z_xxxxz[k] * cd_z[k] + g_z_xxxxzz[k];

            g_zz_xxxyy[k] = -g_z_xxxyy[k] * cd_z[k] + g_z_xxxyyz[k];

            g_zz_xxxyz[k] = -g_z_xxxyz[k] * cd_z[k] + g_z_xxxyzz[k];

            g_zz_xxxzz[k] = -g_z_xxxzz[k] * cd_z[k] + g_z_xxxzzz[k];

            g_zz_xxyyy[k] = -g_z_xxyyy[k] * cd_z[k] + g_z_xxyyyz[k];

            g_zz_xxyyz[k] = -g_z_xxyyz[k] * cd_z[k] + g_z_xxyyzz[k];

            g_zz_xxyzz[k] = -g_z_xxyzz[k] * cd_z[k] + g_z_xxyzzz[k];

            g_zz_xxzzz[k] = -g_z_xxzzz[k] * cd_z[k] + g_z_xxzzzz[k];

            g_zz_xyyyy[k] = -g_z_xyyyy[k] * cd_z[k] + g_z_xyyyyz[k];

            g_zz_xyyyz[k] = -g_z_xyyyz[k] * cd_z[k] + g_z_xyyyzz[k];

            g_zz_xyyzz[k] = -g_z_xyyzz[k] * cd_z[k] + g_z_xyyzzz[k];

            g_zz_xyzzz[k] = -g_z_xyzzz[k] * cd_z[k] + g_z_xyzzzz[k];

            g_zz_xzzzz[k] = -g_z_xzzzz[k] * cd_z[k] + g_z_xzzzzz[k];

            g_zz_yyyyy[k] = -g_z_yyyyy[k] * cd_z[k] + g_z_yyyyyz[k];

            g_zz_yyyyz[k] = -g_z_yyyyz[k] * cd_z[k] + g_z_yyyyzz[k];

            g_zz_yyyzz[k] = -g_z_yyyzz[k] * cd_z[k] + g_z_yyyzzz[k];

            g_zz_yyzzz[k] = -g_z_yyzzz[k] * cd_z[k] + g_z_yyzzzz[k];

            g_zz_yzzzz[k] = -g_z_yzzzz[k] * cd_z[k] + g_z_yzzzzz[k];

            g_zz_zzzzz[k] = -g_z_zzzzz[k] * cd_z[k] + g_z_zzzzzz[k];
        }
    }
}

} // t3ceri namespace


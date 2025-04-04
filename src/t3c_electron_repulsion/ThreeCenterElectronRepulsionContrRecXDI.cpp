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

#include "ThreeCenterElectronRepulsionContrRecXDI.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xdi(CSimdArray<double>& cbuffer,
                                const size_t idx_xdi,
                                const size_t idx_xpi,
                                const size_t idx_xpk,
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

        auto g_x_yyyyyy = cbuffer.data(pi_off + 21);

        auto g_x_yyyyyz = cbuffer.data(pi_off + 22);

        auto g_x_yyyyzz = cbuffer.data(pi_off + 23);

        auto g_x_yyyzzz = cbuffer.data(pi_off + 24);

        auto g_x_yyzzzz = cbuffer.data(pi_off + 25);

        auto g_x_yzzzzz = cbuffer.data(pi_off + 26);

        auto g_x_zzzzzz = cbuffer.data(pi_off + 27);

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

        auto g_y_zzzzzz = cbuffer.data(pi_off + 55);

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

        /// Set up components of auxilary buffer : SPK

        const auto pk_off = idx_xpk + i * 108;

        auto g_x_xxxxxxx = cbuffer.data(pk_off + 0);

        auto g_x_xxxxxxy = cbuffer.data(pk_off + 1);

        auto g_x_xxxxxxz = cbuffer.data(pk_off + 2);

        auto g_x_xxxxxyy = cbuffer.data(pk_off + 3);

        auto g_x_xxxxxyz = cbuffer.data(pk_off + 4);

        auto g_x_xxxxxzz = cbuffer.data(pk_off + 5);

        auto g_x_xxxxyyy = cbuffer.data(pk_off + 6);

        auto g_x_xxxxyyz = cbuffer.data(pk_off + 7);

        auto g_x_xxxxyzz = cbuffer.data(pk_off + 8);

        auto g_x_xxxxzzz = cbuffer.data(pk_off + 9);

        auto g_x_xxxyyyy = cbuffer.data(pk_off + 10);

        auto g_x_xxxyyyz = cbuffer.data(pk_off + 11);

        auto g_x_xxxyyzz = cbuffer.data(pk_off + 12);

        auto g_x_xxxyzzz = cbuffer.data(pk_off + 13);

        auto g_x_xxxzzzz = cbuffer.data(pk_off + 14);

        auto g_x_xxyyyyy = cbuffer.data(pk_off + 15);

        auto g_x_xxyyyyz = cbuffer.data(pk_off + 16);

        auto g_x_xxyyyzz = cbuffer.data(pk_off + 17);

        auto g_x_xxyyzzz = cbuffer.data(pk_off + 18);

        auto g_x_xxyzzzz = cbuffer.data(pk_off + 19);

        auto g_x_xxzzzzz = cbuffer.data(pk_off + 20);

        auto g_x_xyyyyyy = cbuffer.data(pk_off + 21);

        auto g_x_xyyyyyz = cbuffer.data(pk_off + 22);

        auto g_x_xyyyyzz = cbuffer.data(pk_off + 23);

        auto g_x_xyyyzzz = cbuffer.data(pk_off + 24);

        auto g_x_xyyzzzz = cbuffer.data(pk_off + 25);

        auto g_x_xyzzzzz = cbuffer.data(pk_off + 26);

        auto g_x_xzzzzzz = cbuffer.data(pk_off + 27);

        auto g_y_xxxxxxx = cbuffer.data(pk_off + 36);

        auto g_y_xxxxxxy = cbuffer.data(pk_off + 37);

        auto g_y_xxxxxxz = cbuffer.data(pk_off + 38);

        auto g_y_xxxxxyy = cbuffer.data(pk_off + 39);

        auto g_y_xxxxxyz = cbuffer.data(pk_off + 40);

        auto g_y_xxxxxzz = cbuffer.data(pk_off + 41);

        auto g_y_xxxxyyy = cbuffer.data(pk_off + 42);

        auto g_y_xxxxyyz = cbuffer.data(pk_off + 43);

        auto g_y_xxxxyzz = cbuffer.data(pk_off + 44);

        auto g_y_xxxxzzz = cbuffer.data(pk_off + 45);

        auto g_y_xxxyyyy = cbuffer.data(pk_off + 46);

        auto g_y_xxxyyyz = cbuffer.data(pk_off + 47);

        auto g_y_xxxyyzz = cbuffer.data(pk_off + 48);

        auto g_y_xxxyzzz = cbuffer.data(pk_off + 49);

        auto g_y_xxxzzzz = cbuffer.data(pk_off + 50);

        auto g_y_xxyyyyy = cbuffer.data(pk_off + 51);

        auto g_y_xxyyyyz = cbuffer.data(pk_off + 52);

        auto g_y_xxyyyzz = cbuffer.data(pk_off + 53);

        auto g_y_xxyyzzz = cbuffer.data(pk_off + 54);

        auto g_y_xxyzzzz = cbuffer.data(pk_off + 55);

        auto g_y_xxzzzzz = cbuffer.data(pk_off + 56);

        auto g_y_xyyyyyy = cbuffer.data(pk_off + 57);

        auto g_y_xyyyyyz = cbuffer.data(pk_off + 58);

        auto g_y_xyyyyzz = cbuffer.data(pk_off + 59);

        auto g_y_xyyyzzz = cbuffer.data(pk_off + 60);

        auto g_y_xyyzzzz = cbuffer.data(pk_off + 61);

        auto g_y_xyzzzzz = cbuffer.data(pk_off + 62);

        auto g_y_xzzzzzz = cbuffer.data(pk_off + 63);

        auto g_y_yyyyyyy = cbuffer.data(pk_off + 64);

        auto g_y_yyyyyyz = cbuffer.data(pk_off + 65);

        auto g_y_yyyyyzz = cbuffer.data(pk_off + 66);

        auto g_y_yyyyzzz = cbuffer.data(pk_off + 67);

        auto g_y_yyyzzzz = cbuffer.data(pk_off + 68);

        auto g_y_yyzzzzz = cbuffer.data(pk_off + 69);

        auto g_y_yzzzzzz = cbuffer.data(pk_off + 70);

        auto g_z_xxxxxxx = cbuffer.data(pk_off + 72);

        auto g_z_xxxxxxy = cbuffer.data(pk_off + 73);

        auto g_z_xxxxxxz = cbuffer.data(pk_off + 74);

        auto g_z_xxxxxyy = cbuffer.data(pk_off + 75);

        auto g_z_xxxxxyz = cbuffer.data(pk_off + 76);

        auto g_z_xxxxxzz = cbuffer.data(pk_off + 77);

        auto g_z_xxxxyyy = cbuffer.data(pk_off + 78);

        auto g_z_xxxxyyz = cbuffer.data(pk_off + 79);

        auto g_z_xxxxyzz = cbuffer.data(pk_off + 80);

        auto g_z_xxxxzzz = cbuffer.data(pk_off + 81);

        auto g_z_xxxyyyy = cbuffer.data(pk_off + 82);

        auto g_z_xxxyyyz = cbuffer.data(pk_off + 83);

        auto g_z_xxxyyzz = cbuffer.data(pk_off + 84);

        auto g_z_xxxyzzz = cbuffer.data(pk_off + 85);

        auto g_z_xxxzzzz = cbuffer.data(pk_off + 86);

        auto g_z_xxyyyyy = cbuffer.data(pk_off + 87);

        auto g_z_xxyyyyz = cbuffer.data(pk_off + 88);

        auto g_z_xxyyyzz = cbuffer.data(pk_off + 89);

        auto g_z_xxyyzzz = cbuffer.data(pk_off + 90);

        auto g_z_xxyzzzz = cbuffer.data(pk_off + 91);

        auto g_z_xxzzzzz = cbuffer.data(pk_off + 92);

        auto g_z_xyyyyyy = cbuffer.data(pk_off + 93);

        auto g_z_xyyyyyz = cbuffer.data(pk_off + 94);

        auto g_z_xyyyyzz = cbuffer.data(pk_off + 95);

        auto g_z_xyyyzzz = cbuffer.data(pk_off + 96);

        auto g_z_xyyzzzz = cbuffer.data(pk_off + 97);

        auto g_z_xyzzzzz = cbuffer.data(pk_off + 98);

        auto g_z_xzzzzzz = cbuffer.data(pk_off + 99);

        auto g_z_yyyyyyy = cbuffer.data(pk_off + 100);

        auto g_z_yyyyyyz = cbuffer.data(pk_off + 101);

        auto g_z_yyyyyzz = cbuffer.data(pk_off + 102);

        auto g_z_yyyyzzz = cbuffer.data(pk_off + 103);

        auto g_z_yyyzzzz = cbuffer.data(pk_off + 104);

        auto g_z_yyzzzzz = cbuffer.data(pk_off + 105);

        auto g_z_yzzzzzz = cbuffer.data(pk_off + 106);

        auto g_z_zzzzzzz = cbuffer.data(pk_off + 107);

        /// set up bra offset for contr_buffer_xdi

        const auto di_off = idx_xdi + i * 168;

        /// Set up 0-28 components of targeted buffer : cbuffer.data(

        auto g_xx_xxxxxx = cbuffer.data(di_off + 0);

        auto g_xx_xxxxxy = cbuffer.data(di_off + 1);

        auto g_xx_xxxxxz = cbuffer.data(di_off + 2);

        auto g_xx_xxxxyy = cbuffer.data(di_off + 3);

        auto g_xx_xxxxyz = cbuffer.data(di_off + 4);

        auto g_xx_xxxxzz = cbuffer.data(di_off + 5);

        auto g_xx_xxxyyy = cbuffer.data(di_off + 6);

        auto g_xx_xxxyyz = cbuffer.data(di_off + 7);

        auto g_xx_xxxyzz = cbuffer.data(di_off + 8);

        auto g_xx_xxxzzz = cbuffer.data(di_off + 9);

        auto g_xx_xxyyyy = cbuffer.data(di_off + 10);

        auto g_xx_xxyyyz = cbuffer.data(di_off + 11);

        auto g_xx_xxyyzz = cbuffer.data(di_off + 12);

        auto g_xx_xxyzzz = cbuffer.data(di_off + 13);

        auto g_xx_xxzzzz = cbuffer.data(di_off + 14);

        auto g_xx_xyyyyy = cbuffer.data(di_off + 15);

        auto g_xx_xyyyyz = cbuffer.data(di_off + 16);

        auto g_xx_xyyyzz = cbuffer.data(di_off + 17);

        auto g_xx_xyyzzz = cbuffer.data(di_off + 18);

        auto g_xx_xyzzzz = cbuffer.data(di_off + 19);

        auto g_xx_xzzzzz = cbuffer.data(di_off + 20);

        auto g_xx_yyyyyy = cbuffer.data(di_off + 21);

        auto g_xx_yyyyyz = cbuffer.data(di_off + 22);

        auto g_xx_yyyyzz = cbuffer.data(di_off + 23);

        auto g_xx_yyyzzz = cbuffer.data(di_off + 24);

        auto g_xx_yyzzzz = cbuffer.data(di_off + 25);

        auto g_xx_yzzzzz = cbuffer.data(di_off + 26);

        auto g_xx_zzzzzz = cbuffer.data(di_off + 27);

        #pragma omp simd aligned(cd_x, g_x_xxxxxx, g_x_xxxxxxx, g_x_xxxxxxy, g_x_xxxxxxz, g_x_xxxxxy, g_x_xxxxxyy, g_x_xxxxxyz, g_x_xxxxxz, g_x_xxxxxzz, g_x_xxxxyy, g_x_xxxxyyy, g_x_xxxxyyz, g_x_xxxxyz, g_x_xxxxyzz, g_x_xxxxzz, g_x_xxxxzzz, g_x_xxxyyy, g_x_xxxyyyy, g_x_xxxyyyz, g_x_xxxyyz, g_x_xxxyyzz, g_x_xxxyzz, g_x_xxxyzzz, g_x_xxxzzz, g_x_xxxzzzz, g_x_xxyyyy, g_x_xxyyyyy, g_x_xxyyyyz, g_x_xxyyyz, g_x_xxyyyzz, g_x_xxyyzz, g_x_xxyyzzz, g_x_xxyzzz, g_x_xxyzzzz, g_x_xxzzzz, g_x_xxzzzzz, g_x_xyyyyy, g_x_xyyyyyy, g_x_xyyyyyz, g_x_xyyyyz, g_x_xyyyyzz, g_x_xyyyzz, g_x_xyyyzzz, g_x_xyyzzz, g_x_xyyzzzz, g_x_xyzzzz, g_x_xyzzzzz, g_x_xzzzzz, g_x_xzzzzzz, g_x_yyyyyy, g_x_yyyyyz, g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz, g_x_yzzzzz, g_x_zzzzzz, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz, g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz, g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy, g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz, g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz, g_xx_xyyzzz, g_xx_xyzzzz, g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz, g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz, g_xx_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xx_xxxxxx[k] = -g_x_xxxxxx[k] * cd_x[k] + g_x_xxxxxxx[k];

            g_xx_xxxxxy[k] = -g_x_xxxxxy[k] * cd_x[k] + g_x_xxxxxxy[k];

            g_xx_xxxxxz[k] = -g_x_xxxxxz[k] * cd_x[k] + g_x_xxxxxxz[k];

            g_xx_xxxxyy[k] = -g_x_xxxxyy[k] * cd_x[k] + g_x_xxxxxyy[k];

            g_xx_xxxxyz[k] = -g_x_xxxxyz[k] * cd_x[k] + g_x_xxxxxyz[k];

            g_xx_xxxxzz[k] = -g_x_xxxxzz[k] * cd_x[k] + g_x_xxxxxzz[k];

            g_xx_xxxyyy[k] = -g_x_xxxyyy[k] * cd_x[k] + g_x_xxxxyyy[k];

            g_xx_xxxyyz[k] = -g_x_xxxyyz[k] * cd_x[k] + g_x_xxxxyyz[k];

            g_xx_xxxyzz[k] = -g_x_xxxyzz[k] * cd_x[k] + g_x_xxxxyzz[k];

            g_xx_xxxzzz[k] = -g_x_xxxzzz[k] * cd_x[k] + g_x_xxxxzzz[k];

            g_xx_xxyyyy[k] = -g_x_xxyyyy[k] * cd_x[k] + g_x_xxxyyyy[k];

            g_xx_xxyyyz[k] = -g_x_xxyyyz[k] * cd_x[k] + g_x_xxxyyyz[k];

            g_xx_xxyyzz[k] = -g_x_xxyyzz[k] * cd_x[k] + g_x_xxxyyzz[k];

            g_xx_xxyzzz[k] = -g_x_xxyzzz[k] * cd_x[k] + g_x_xxxyzzz[k];

            g_xx_xxzzzz[k] = -g_x_xxzzzz[k] * cd_x[k] + g_x_xxxzzzz[k];

            g_xx_xyyyyy[k] = -g_x_xyyyyy[k] * cd_x[k] + g_x_xxyyyyy[k];

            g_xx_xyyyyz[k] = -g_x_xyyyyz[k] * cd_x[k] + g_x_xxyyyyz[k];

            g_xx_xyyyzz[k] = -g_x_xyyyzz[k] * cd_x[k] + g_x_xxyyyzz[k];

            g_xx_xyyzzz[k] = -g_x_xyyzzz[k] * cd_x[k] + g_x_xxyyzzz[k];

            g_xx_xyzzzz[k] = -g_x_xyzzzz[k] * cd_x[k] + g_x_xxyzzzz[k];

            g_xx_xzzzzz[k] = -g_x_xzzzzz[k] * cd_x[k] + g_x_xxzzzzz[k];

            g_xx_yyyyyy[k] = -g_x_yyyyyy[k] * cd_x[k] + g_x_xyyyyyy[k];

            g_xx_yyyyyz[k] = -g_x_yyyyyz[k] * cd_x[k] + g_x_xyyyyyz[k];

            g_xx_yyyyzz[k] = -g_x_yyyyzz[k] * cd_x[k] + g_x_xyyyyzz[k];

            g_xx_yyyzzz[k] = -g_x_yyyzzz[k] * cd_x[k] + g_x_xyyyzzz[k];

            g_xx_yyzzzz[k] = -g_x_yyzzzz[k] * cd_x[k] + g_x_xyyzzzz[k];

            g_xx_yzzzzz[k] = -g_x_yzzzzz[k] * cd_x[k] + g_x_xyzzzzz[k];

            g_xx_zzzzzz[k] = -g_x_zzzzzz[k] * cd_x[k] + g_x_xzzzzzz[k];
        }

        /// Set up 28-56 components of targeted buffer : cbuffer.data(

        auto g_xy_xxxxxx = cbuffer.data(di_off + 28);

        auto g_xy_xxxxxy = cbuffer.data(di_off + 29);

        auto g_xy_xxxxxz = cbuffer.data(di_off + 30);

        auto g_xy_xxxxyy = cbuffer.data(di_off + 31);

        auto g_xy_xxxxyz = cbuffer.data(di_off + 32);

        auto g_xy_xxxxzz = cbuffer.data(di_off + 33);

        auto g_xy_xxxyyy = cbuffer.data(di_off + 34);

        auto g_xy_xxxyyz = cbuffer.data(di_off + 35);

        auto g_xy_xxxyzz = cbuffer.data(di_off + 36);

        auto g_xy_xxxzzz = cbuffer.data(di_off + 37);

        auto g_xy_xxyyyy = cbuffer.data(di_off + 38);

        auto g_xy_xxyyyz = cbuffer.data(di_off + 39);

        auto g_xy_xxyyzz = cbuffer.data(di_off + 40);

        auto g_xy_xxyzzz = cbuffer.data(di_off + 41);

        auto g_xy_xxzzzz = cbuffer.data(di_off + 42);

        auto g_xy_xyyyyy = cbuffer.data(di_off + 43);

        auto g_xy_xyyyyz = cbuffer.data(di_off + 44);

        auto g_xy_xyyyzz = cbuffer.data(di_off + 45);

        auto g_xy_xyyzzz = cbuffer.data(di_off + 46);

        auto g_xy_xyzzzz = cbuffer.data(di_off + 47);

        auto g_xy_xzzzzz = cbuffer.data(di_off + 48);

        auto g_xy_yyyyyy = cbuffer.data(di_off + 49);

        auto g_xy_yyyyyz = cbuffer.data(di_off + 50);

        auto g_xy_yyyyzz = cbuffer.data(di_off + 51);

        auto g_xy_yyyzzz = cbuffer.data(di_off + 52);

        auto g_xy_yyzzzz = cbuffer.data(di_off + 53);

        auto g_xy_yzzzzz = cbuffer.data(di_off + 54);

        auto g_xy_zzzzzz = cbuffer.data(di_off + 55);

        #pragma omp simd aligned(cd_x, g_xy_xxxxxx, g_xy_xxxxxy, g_xy_xxxxxz, g_xy_xxxxyy, g_xy_xxxxyz, g_xy_xxxxzz, g_xy_xxxyyy, g_xy_xxxyyz, g_xy_xxxyzz, g_xy_xxxzzz, g_xy_xxyyyy, g_xy_xxyyyz, g_xy_xxyyzz, g_xy_xxyzzz, g_xy_xxzzzz, g_xy_xyyyyy, g_xy_xyyyyz, g_xy_xyyyzz, g_xy_xyyzzz, g_xy_xyzzzz, g_xy_xzzzzz, g_xy_yyyyyy, g_xy_yyyyyz, g_xy_yyyyzz, g_xy_yyyzzz, g_xy_yyzzzz, g_xy_yzzzzz, g_xy_zzzzzz, g_y_xxxxxx, g_y_xxxxxxx, g_y_xxxxxxy, g_y_xxxxxxz, g_y_xxxxxy, g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxz, g_y_xxxxxzz, g_y_xxxxyy, g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyz, g_y_xxxxyzz, g_y_xxxxzz, g_y_xxxxzzz, g_y_xxxyyy, g_y_xxxyyyy, g_y_xxxyyyz, g_y_xxxyyz, g_y_xxxyyzz, g_y_xxxyzz, g_y_xxxyzzz, g_y_xxxzzz, g_y_xxxzzzz, g_y_xxyyyy, g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyz, g_y_xxyyyzz, g_y_xxyyzz, g_y_xxyyzzz, g_y_xxyzzz, g_y_xxyzzzz, g_y_xxzzzz, g_y_xxzzzzz, g_y_xyyyyy, g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyz, g_y_xyyyyzz, g_y_xyyyzz, g_y_xyyyzzz, g_y_xyyzzz, g_y_xyyzzzz, g_y_xyzzzz, g_y_xyzzzzz, g_y_xzzzzz, g_y_xzzzzzz, g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyzz, g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz, g_y_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xy_xxxxxx[k] = -g_y_xxxxxx[k] * cd_x[k] + g_y_xxxxxxx[k];

            g_xy_xxxxxy[k] = -g_y_xxxxxy[k] * cd_x[k] + g_y_xxxxxxy[k];

            g_xy_xxxxxz[k] = -g_y_xxxxxz[k] * cd_x[k] + g_y_xxxxxxz[k];

            g_xy_xxxxyy[k] = -g_y_xxxxyy[k] * cd_x[k] + g_y_xxxxxyy[k];

            g_xy_xxxxyz[k] = -g_y_xxxxyz[k] * cd_x[k] + g_y_xxxxxyz[k];

            g_xy_xxxxzz[k] = -g_y_xxxxzz[k] * cd_x[k] + g_y_xxxxxzz[k];

            g_xy_xxxyyy[k] = -g_y_xxxyyy[k] * cd_x[k] + g_y_xxxxyyy[k];

            g_xy_xxxyyz[k] = -g_y_xxxyyz[k] * cd_x[k] + g_y_xxxxyyz[k];

            g_xy_xxxyzz[k] = -g_y_xxxyzz[k] * cd_x[k] + g_y_xxxxyzz[k];

            g_xy_xxxzzz[k] = -g_y_xxxzzz[k] * cd_x[k] + g_y_xxxxzzz[k];

            g_xy_xxyyyy[k] = -g_y_xxyyyy[k] * cd_x[k] + g_y_xxxyyyy[k];

            g_xy_xxyyyz[k] = -g_y_xxyyyz[k] * cd_x[k] + g_y_xxxyyyz[k];

            g_xy_xxyyzz[k] = -g_y_xxyyzz[k] * cd_x[k] + g_y_xxxyyzz[k];

            g_xy_xxyzzz[k] = -g_y_xxyzzz[k] * cd_x[k] + g_y_xxxyzzz[k];

            g_xy_xxzzzz[k] = -g_y_xxzzzz[k] * cd_x[k] + g_y_xxxzzzz[k];

            g_xy_xyyyyy[k] = -g_y_xyyyyy[k] * cd_x[k] + g_y_xxyyyyy[k];

            g_xy_xyyyyz[k] = -g_y_xyyyyz[k] * cd_x[k] + g_y_xxyyyyz[k];

            g_xy_xyyyzz[k] = -g_y_xyyyzz[k] * cd_x[k] + g_y_xxyyyzz[k];

            g_xy_xyyzzz[k] = -g_y_xyyzzz[k] * cd_x[k] + g_y_xxyyzzz[k];

            g_xy_xyzzzz[k] = -g_y_xyzzzz[k] * cd_x[k] + g_y_xxyzzzz[k];

            g_xy_xzzzzz[k] = -g_y_xzzzzz[k] * cd_x[k] + g_y_xxzzzzz[k];

            g_xy_yyyyyy[k] = -g_y_yyyyyy[k] * cd_x[k] + g_y_xyyyyyy[k];

            g_xy_yyyyyz[k] = -g_y_yyyyyz[k] * cd_x[k] + g_y_xyyyyyz[k];

            g_xy_yyyyzz[k] = -g_y_yyyyzz[k] * cd_x[k] + g_y_xyyyyzz[k];

            g_xy_yyyzzz[k] = -g_y_yyyzzz[k] * cd_x[k] + g_y_xyyyzzz[k];

            g_xy_yyzzzz[k] = -g_y_yyzzzz[k] * cd_x[k] + g_y_xyyzzzz[k];

            g_xy_yzzzzz[k] = -g_y_yzzzzz[k] * cd_x[k] + g_y_xyzzzzz[k];

            g_xy_zzzzzz[k] = -g_y_zzzzzz[k] * cd_x[k] + g_y_xzzzzzz[k];
        }

        /// Set up 56-84 components of targeted buffer : cbuffer.data(

        auto g_xz_xxxxxx = cbuffer.data(di_off + 56);

        auto g_xz_xxxxxy = cbuffer.data(di_off + 57);

        auto g_xz_xxxxxz = cbuffer.data(di_off + 58);

        auto g_xz_xxxxyy = cbuffer.data(di_off + 59);

        auto g_xz_xxxxyz = cbuffer.data(di_off + 60);

        auto g_xz_xxxxzz = cbuffer.data(di_off + 61);

        auto g_xz_xxxyyy = cbuffer.data(di_off + 62);

        auto g_xz_xxxyyz = cbuffer.data(di_off + 63);

        auto g_xz_xxxyzz = cbuffer.data(di_off + 64);

        auto g_xz_xxxzzz = cbuffer.data(di_off + 65);

        auto g_xz_xxyyyy = cbuffer.data(di_off + 66);

        auto g_xz_xxyyyz = cbuffer.data(di_off + 67);

        auto g_xz_xxyyzz = cbuffer.data(di_off + 68);

        auto g_xz_xxyzzz = cbuffer.data(di_off + 69);

        auto g_xz_xxzzzz = cbuffer.data(di_off + 70);

        auto g_xz_xyyyyy = cbuffer.data(di_off + 71);

        auto g_xz_xyyyyz = cbuffer.data(di_off + 72);

        auto g_xz_xyyyzz = cbuffer.data(di_off + 73);

        auto g_xz_xyyzzz = cbuffer.data(di_off + 74);

        auto g_xz_xyzzzz = cbuffer.data(di_off + 75);

        auto g_xz_xzzzzz = cbuffer.data(di_off + 76);

        auto g_xz_yyyyyy = cbuffer.data(di_off + 77);

        auto g_xz_yyyyyz = cbuffer.data(di_off + 78);

        auto g_xz_yyyyzz = cbuffer.data(di_off + 79);

        auto g_xz_yyyzzz = cbuffer.data(di_off + 80);

        auto g_xz_yyzzzz = cbuffer.data(di_off + 81);

        auto g_xz_yzzzzz = cbuffer.data(di_off + 82);

        auto g_xz_zzzzzz = cbuffer.data(di_off + 83);

        #pragma omp simd aligned(cd_x, g_xz_xxxxxx, g_xz_xxxxxy, g_xz_xxxxxz, g_xz_xxxxyy, g_xz_xxxxyz, g_xz_xxxxzz, g_xz_xxxyyy, g_xz_xxxyyz, g_xz_xxxyzz, g_xz_xxxzzz, g_xz_xxyyyy, g_xz_xxyyyz, g_xz_xxyyzz, g_xz_xxyzzz, g_xz_xxzzzz, g_xz_xyyyyy, g_xz_xyyyyz, g_xz_xyyyzz, g_xz_xyyzzz, g_xz_xyzzzz, g_xz_xzzzzz, g_xz_yyyyyy, g_xz_yyyyyz, g_xz_yyyyzz, g_xz_yyyzzz, g_xz_yyzzzz, g_xz_yzzzzz, g_xz_zzzzzz, g_z_xxxxxx, g_z_xxxxxxx, g_z_xxxxxxy, g_z_xxxxxxz, g_z_xxxxxy, g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxxzz, g_z_xxxxyy, g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxxzzz, g_z_xxxyyy, g_z_xxxyyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxxzzzz, g_z_xxyyyy, g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xxzzzzz, g_z_xyyyyy, g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_xzzzzzz, g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz, g_z_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xz_xxxxxx[k] = -g_z_xxxxxx[k] * cd_x[k] + g_z_xxxxxxx[k];

            g_xz_xxxxxy[k] = -g_z_xxxxxy[k] * cd_x[k] + g_z_xxxxxxy[k];

            g_xz_xxxxxz[k] = -g_z_xxxxxz[k] * cd_x[k] + g_z_xxxxxxz[k];

            g_xz_xxxxyy[k] = -g_z_xxxxyy[k] * cd_x[k] + g_z_xxxxxyy[k];

            g_xz_xxxxyz[k] = -g_z_xxxxyz[k] * cd_x[k] + g_z_xxxxxyz[k];

            g_xz_xxxxzz[k] = -g_z_xxxxzz[k] * cd_x[k] + g_z_xxxxxzz[k];

            g_xz_xxxyyy[k] = -g_z_xxxyyy[k] * cd_x[k] + g_z_xxxxyyy[k];

            g_xz_xxxyyz[k] = -g_z_xxxyyz[k] * cd_x[k] + g_z_xxxxyyz[k];

            g_xz_xxxyzz[k] = -g_z_xxxyzz[k] * cd_x[k] + g_z_xxxxyzz[k];

            g_xz_xxxzzz[k] = -g_z_xxxzzz[k] * cd_x[k] + g_z_xxxxzzz[k];

            g_xz_xxyyyy[k] = -g_z_xxyyyy[k] * cd_x[k] + g_z_xxxyyyy[k];

            g_xz_xxyyyz[k] = -g_z_xxyyyz[k] * cd_x[k] + g_z_xxxyyyz[k];

            g_xz_xxyyzz[k] = -g_z_xxyyzz[k] * cd_x[k] + g_z_xxxyyzz[k];

            g_xz_xxyzzz[k] = -g_z_xxyzzz[k] * cd_x[k] + g_z_xxxyzzz[k];

            g_xz_xxzzzz[k] = -g_z_xxzzzz[k] * cd_x[k] + g_z_xxxzzzz[k];

            g_xz_xyyyyy[k] = -g_z_xyyyyy[k] * cd_x[k] + g_z_xxyyyyy[k];

            g_xz_xyyyyz[k] = -g_z_xyyyyz[k] * cd_x[k] + g_z_xxyyyyz[k];

            g_xz_xyyyzz[k] = -g_z_xyyyzz[k] * cd_x[k] + g_z_xxyyyzz[k];

            g_xz_xyyzzz[k] = -g_z_xyyzzz[k] * cd_x[k] + g_z_xxyyzzz[k];

            g_xz_xyzzzz[k] = -g_z_xyzzzz[k] * cd_x[k] + g_z_xxyzzzz[k];

            g_xz_xzzzzz[k] = -g_z_xzzzzz[k] * cd_x[k] + g_z_xxzzzzz[k];

            g_xz_yyyyyy[k] = -g_z_yyyyyy[k] * cd_x[k] + g_z_xyyyyyy[k];

            g_xz_yyyyyz[k] = -g_z_yyyyyz[k] * cd_x[k] + g_z_xyyyyyz[k];

            g_xz_yyyyzz[k] = -g_z_yyyyzz[k] * cd_x[k] + g_z_xyyyyzz[k];

            g_xz_yyyzzz[k] = -g_z_yyyzzz[k] * cd_x[k] + g_z_xyyyzzz[k];

            g_xz_yyzzzz[k] = -g_z_yyzzzz[k] * cd_x[k] + g_z_xyyzzzz[k];

            g_xz_yzzzzz[k] = -g_z_yzzzzz[k] * cd_x[k] + g_z_xyzzzzz[k];

            g_xz_zzzzzz[k] = -g_z_zzzzzz[k] * cd_x[k] + g_z_xzzzzzz[k];
        }

        /// Set up 84-112 components of targeted buffer : cbuffer.data(

        auto g_yy_xxxxxx = cbuffer.data(di_off + 84);

        auto g_yy_xxxxxy = cbuffer.data(di_off + 85);

        auto g_yy_xxxxxz = cbuffer.data(di_off + 86);

        auto g_yy_xxxxyy = cbuffer.data(di_off + 87);

        auto g_yy_xxxxyz = cbuffer.data(di_off + 88);

        auto g_yy_xxxxzz = cbuffer.data(di_off + 89);

        auto g_yy_xxxyyy = cbuffer.data(di_off + 90);

        auto g_yy_xxxyyz = cbuffer.data(di_off + 91);

        auto g_yy_xxxyzz = cbuffer.data(di_off + 92);

        auto g_yy_xxxzzz = cbuffer.data(di_off + 93);

        auto g_yy_xxyyyy = cbuffer.data(di_off + 94);

        auto g_yy_xxyyyz = cbuffer.data(di_off + 95);

        auto g_yy_xxyyzz = cbuffer.data(di_off + 96);

        auto g_yy_xxyzzz = cbuffer.data(di_off + 97);

        auto g_yy_xxzzzz = cbuffer.data(di_off + 98);

        auto g_yy_xyyyyy = cbuffer.data(di_off + 99);

        auto g_yy_xyyyyz = cbuffer.data(di_off + 100);

        auto g_yy_xyyyzz = cbuffer.data(di_off + 101);

        auto g_yy_xyyzzz = cbuffer.data(di_off + 102);

        auto g_yy_xyzzzz = cbuffer.data(di_off + 103);

        auto g_yy_xzzzzz = cbuffer.data(di_off + 104);

        auto g_yy_yyyyyy = cbuffer.data(di_off + 105);

        auto g_yy_yyyyyz = cbuffer.data(di_off + 106);

        auto g_yy_yyyyzz = cbuffer.data(di_off + 107);

        auto g_yy_yyyzzz = cbuffer.data(di_off + 108);

        auto g_yy_yyzzzz = cbuffer.data(di_off + 109);

        auto g_yy_yzzzzz = cbuffer.data(di_off + 110);

        auto g_yy_zzzzzz = cbuffer.data(di_off + 111);

        #pragma omp simd aligned(cd_y, g_y_xxxxxx, g_y_xxxxxxy, g_y_xxxxxy, g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxz, g_y_xxxxyy, g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyz, g_y_xxxxyzz, g_y_xxxxzz, g_y_xxxyyy, g_y_xxxyyyy, g_y_xxxyyyz, g_y_xxxyyz, g_y_xxxyyzz, g_y_xxxyzz, g_y_xxxyzzz, g_y_xxxzzz, g_y_xxyyyy, g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyz, g_y_xxyyyzz, g_y_xxyyzz, g_y_xxyyzzz, g_y_xxyzzz, g_y_xxyzzzz, g_y_xxzzzz, g_y_xyyyyy, g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyz, g_y_xyyyyzz, g_y_xyyyzz, g_y_xyyyzzz, g_y_xyyzzz, g_y_xyyzzzz, g_y_xyzzzz, g_y_xyzzzzz, g_y_xzzzzz, g_y_yyyyyy, g_y_yyyyyyy, g_y_yyyyyyz, g_y_yyyyyz, g_y_yyyyyzz, g_y_yyyyzz, g_y_yyyyzzz, g_y_yyyzzz, g_y_yyyzzzz, g_y_yyzzzz, g_y_yyzzzzz, g_y_yzzzzz, g_y_yzzzzzz, g_y_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz, g_yy_xyyzzz, g_yy_xyzzzz, g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz, g_yy_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yy_xxxxxx[k] = -g_y_xxxxxx[k] * cd_y[k] + g_y_xxxxxxy[k];

            g_yy_xxxxxy[k] = -g_y_xxxxxy[k] * cd_y[k] + g_y_xxxxxyy[k];

            g_yy_xxxxxz[k] = -g_y_xxxxxz[k] * cd_y[k] + g_y_xxxxxyz[k];

            g_yy_xxxxyy[k] = -g_y_xxxxyy[k] * cd_y[k] + g_y_xxxxyyy[k];

            g_yy_xxxxyz[k] = -g_y_xxxxyz[k] * cd_y[k] + g_y_xxxxyyz[k];

            g_yy_xxxxzz[k] = -g_y_xxxxzz[k] * cd_y[k] + g_y_xxxxyzz[k];

            g_yy_xxxyyy[k] = -g_y_xxxyyy[k] * cd_y[k] + g_y_xxxyyyy[k];

            g_yy_xxxyyz[k] = -g_y_xxxyyz[k] * cd_y[k] + g_y_xxxyyyz[k];

            g_yy_xxxyzz[k] = -g_y_xxxyzz[k] * cd_y[k] + g_y_xxxyyzz[k];

            g_yy_xxxzzz[k] = -g_y_xxxzzz[k] * cd_y[k] + g_y_xxxyzzz[k];

            g_yy_xxyyyy[k] = -g_y_xxyyyy[k] * cd_y[k] + g_y_xxyyyyy[k];

            g_yy_xxyyyz[k] = -g_y_xxyyyz[k] * cd_y[k] + g_y_xxyyyyz[k];

            g_yy_xxyyzz[k] = -g_y_xxyyzz[k] * cd_y[k] + g_y_xxyyyzz[k];

            g_yy_xxyzzz[k] = -g_y_xxyzzz[k] * cd_y[k] + g_y_xxyyzzz[k];

            g_yy_xxzzzz[k] = -g_y_xxzzzz[k] * cd_y[k] + g_y_xxyzzzz[k];

            g_yy_xyyyyy[k] = -g_y_xyyyyy[k] * cd_y[k] + g_y_xyyyyyy[k];

            g_yy_xyyyyz[k] = -g_y_xyyyyz[k] * cd_y[k] + g_y_xyyyyyz[k];

            g_yy_xyyyzz[k] = -g_y_xyyyzz[k] * cd_y[k] + g_y_xyyyyzz[k];

            g_yy_xyyzzz[k] = -g_y_xyyzzz[k] * cd_y[k] + g_y_xyyyzzz[k];

            g_yy_xyzzzz[k] = -g_y_xyzzzz[k] * cd_y[k] + g_y_xyyzzzz[k];

            g_yy_xzzzzz[k] = -g_y_xzzzzz[k] * cd_y[k] + g_y_xyzzzzz[k];

            g_yy_yyyyyy[k] = -g_y_yyyyyy[k] * cd_y[k] + g_y_yyyyyyy[k];

            g_yy_yyyyyz[k] = -g_y_yyyyyz[k] * cd_y[k] + g_y_yyyyyyz[k];

            g_yy_yyyyzz[k] = -g_y_yyyyzz[k] * cd_y[k] + g_y_yyyyyzz[k];

            g_yy_yyyzzz[k] = -g_y_yyyzzz[k] * cd_y[k] + g_y_yyyyzzz[k];

            g_yy_yyzzzz[k] = -g_y_yyzzzz[k] * cd_y[k] + g_y_yyyzzzz[k];

            g_yy_yzzzzz[k] = -g_y_yzzzzz[k] * cd_y[k] + g_y_yyzzzzz[k];

            g_yy_zzzzzz[k] = -g_y_zzzzzz[k] * cd_y[k] + g_y_yzzzzzz[k];
        }

        /// Set up 112-140 components of targeted buffer : cbuffer.data(

        auto g_yz_xxxxxx = cbuffer.data(di_off + 112);

        auto g_yz_xxxxxy = cbuffer.data(di_off + 113);

        auto g_yz_xxxxxz = cbuffer.data(di_off + 114);

        auto g_yz_xxxxyy = cbuffer.data(di_off + 115);

        auto g_yz_xxxxyz = cbuffer.data(di_off + 116);

        auto g_yz_xxxxzz = cbuffer.data(di_off + 117);

        auto g_yz_xxxyyy = cbuffer.data(di_off + 118);

        auto g_yz_xxxyyz = cbuffer.data(di_off + 119);

        auto g_yz_xxxyzz = cbuffer.data(di_off + 120);

        auto g_yz_xxxzzz = cbuffer.data(di_off + 121);

        auto g_yz_xxyyyy = cbuffer.data(di_off + 122);

        auto g_yz_xxyyyz = cbuffer.data(di_off + 123);

        auto g_yz_xxyyzz = cbuffer.data(di_off + 124);

        auto g_yz_xxyzzz = cbuffer.data(di_off + 125);

        auto g_yz_xxzzzz = cbuffer.data(di_off + 126);

        auto g_yz_xyyyyy = cbuffer.data(di_off + 127);

        auto g_yz_xyyyyz = cbuffer.data(di_off + 128);

        auto g_yz_xyyyzz = cbuffer.data(di_off + 129);

        auto g_yz_xyyzzz = cbuffer.data(di_off + 130);

        auto g_yz_xyzzzz = cbuffer.data(di_off + 131);

        auto g_yz_xzzzzz = cbuffer.data(di_off + 132);

        auto g_yz_yyyyyy = cbuffer.data(di_off + 133);

        auto g_yz_yyyyyz = cbuffer.data(di_off + 134);

        auto g_yz_yyyyzz = cbuffer.data(di_off + 135);

        auto g_yz_yyyzzz = cbuffer.data(di_off + 136);

        auto g_yz_yyzzzz = cbuffer.data(di_off + 137);

        auto g_yz_yzzzzz = cbuffer.data(di_off + 138);

        auto g_yz_zzzzzz = cbuffer.data(di_off + 139);

        #pragma omp simd aligned(cd_y, g_yz_xxxxxx, g_yz_xxxxxy, g_yz_xxxxxz, g_yz_xxxxyy, g_yz_xxxxyz, g_yz_xxxxzz, g_yz_xxxyyy, g_yz_xxxyyz, g_yz_xxxyzz, g_yz_xxxzzz, g_yz_xxyyyy, g_yz_xxyyyz, g_yz_xxyyzz, g_yz_xxyzzz, g_yz_xxzzzz, g_yz_xyyyyy, g_yz_xyyyyz, g_yz_xyyyzz, g_yz_xyyzzz, g_yz_xyzzzz, g_yz_xzzzzz, g_yz_yyyyyy, g_yz_yyyyyz, g_yz_yyyyzz, g_yz_yyyzzz, g_yz_yyzzzz, g_yz_yzzzzz, g_yz_zzzzzz, g_z_xxxxxx, g_z_xxxxxxy, g_z_xxxxxy, g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxyy, g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxyyy, g_z_xxxyyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxyyyy, g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xyyyyy, g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyyy, g_z_yyyyyyz, g_z_yyyyyz, g_z_yyyyyzz, g_z_yyyyzz, g_z_yyyyzzz, g_z_yyyzzz, g_z_yyyzzzz, g_z_yyzzzz, g_z_yyzzzzz, g_z_yzzzzz, g_z_yzzzzzz, g_z_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yz_xxxxxx[k] = -g_z_xxxxxx[k] * cd_y[k] + g_z_xxxxxxy[k];

            g_yz_xxxxxy[k] = -g_z_xxxxxy[k] * cd_y[k] + g_z_xxxxxyy[k];

            g_yz_xxxxxz[k] = -g_z_xxxxxz[k] * cd_y[k] + g_z_xxxxxyz[k];

            g_yz_xxxxyy[k] = -g_z_xxxxyy[k] * cd_y[k] + g_z_xxxxyyy[k];

            g_yz_xxxxyz[k] = -g_z_xxxxyz[k] * cd_y[k] + g_z_xxxxyyz[k];

            g_yz_xxxxzz[k] = -g_z_xxxxzz[k] * cd_y[k] + g_z_xxxxyzz[k];

            g_yz_xxxyyy[k] = -g_z_xxxyyy[k] * cd_y[k] + g_z_xxxyyyy[k];

            g_yz_xxxyyz[k] = -g_z_xxxyyz[k] * cd_y[k] + g_z_xxxyyyz[k];

            g_yz_xxxyzz[k] = -g_z_xxxyzz[k] * cd_y[k] + g_z_xxxyyzz[k];

            g_yz_xxxzzz[k] = -g_z_xxxzzz[k] * cd_y[k] + g_z_xxxyzzz[k];

            g_yz_xxyyyy[k] = -g_z_xxyyyy[k] * cd_y[k] + g_z_xxyyyyy[k];

            g_yz_xxyyyz[k] = -g_z_xxyyyz[k] * cd_y[k] + g_z_xxyyyyz[k];

            g_yz_xxyyzz[k] = -g_z_xxyyzz[k] * cd_y[k] + g_z_xxyyyzz[k];

            g_yz_xxyzzz[k] = -g_z_xxyzzz[k] * cd_y[k] + g_z_xxyyzzz[k];

            g_yz_xxzzzz[k] = -g_z_xxzzzz[k] * cd_y[k] + g_z_xxyzzzz[k];

            g_yz_xyyyyy[k] = -g_z_xyyyyy[k] * cd_y[k] + g_z_xyyyyyy[k];

            g_yz_xyyyyz[k] = -g_z_xyyyyz[k] * cd_y[k] + g_z_xyyyyyz[k];

            g_yz_xyyyzz[k] = -g_z_xyyyzz[k] * cd_y[k] + g_z_xyyyyzz[k];

            g_yz_xyyzzz[k] = -g_z_xyyzzz[k] * cd_y[k] + g_z_xyyyzzz[k];

            g_yz_xyzzzz[k] = -g_z_xyzzzz[k] * cd_y[k] + g_z_xyyzzzz[k];

            g_yz_xzzzzz[k] = -g_z_xzzzzz[k] * cd_y[k] + g_z_xyzzzzz[k];

            g_yz_yyyyyy[k] = -g_z_yyyyyy[k] * cd_y[k] + g_z_yyyyyyy[k];

            g_yz_yyyyyz[k] = -g_z_yyyyyz[k] * cd_y[k] + g_z_yyyyyyz[k];

            g_yz_yyyyzz[k] = -g_z_yyyyzz[k] * cd_y[k] + g_z_yyyyyzz[k];

            g_yz_yyyzzz[k] = -g_z_yyyzzz[k] * cd_y[k] + g_z_yyyyzzz[k];

            g_yz_yyzzzz[k] = -g_z_yyzzzz[k] * cd_y[k] + g_z_yyyzzzz[k];

            g_yz_yzzzzz[k] = -g_z_yzzzzz[k] * cd_y[k] + g_z_yyzzzzz[k];

            g_yz_zzzzzz[k] = -g_z_zzzzzz[k] * cd_y[k] + g_z_yzzzzzz[k];
        }

        /// Set up 140-168 components of targeted buffer : cbuffer.data(

        auto g_zz_xxxxxx = cbuffer.data(di_off + 140);

        auto g_zz_xxxxxy = cbuffer.data(di_off + 141);

        auto g_zz_xxxxxz = cbuffer.data(di_off + 142);

        auto g_zz_xxxxyy = cbuffer.data(di_off + 143);

        auto g_zz_xxxxyz = cbuffer.data(di_off + 144);

        auto g_zz_xxxxzz = cbuffer.data(di_off + 145);

        auto g_zz_xxxyyy = cbuffer.data(di_off + 146);

        auto g_zz_xxxyyz = cbuffer.data(di_off + 147);

        auto g_zz_xxxyzz = cbuffer.data(di_off + 148);

        auto g_zz_xxxzzz = cbuffer.data(di_off + 149);

        auto g_zz_xxyyyy = cbuffer.data(di_off + 150);

        auto g_zz_xxyyyz = cbuffer.data(di_off + 151);

        auto g_zz_xxyyzz = cbuffer.data(di_off + 152);

        auto g_zz_xxyzzz = cbuffer.data(di_off + 153);

        auto g_zz_xxzzzz = cbuffer.data(di_off + 154);

        auto g_zz_xyyyyy = cbuffer.data(di_off + 155);

        auto g_zz_xyyyyz = cbuffer.data(di_off + 156);

        auto g_zz_xyyyzz = cbuffer.data(di_off + 157);

        auto g_zz_xyyzzz = cbuffer.data(di_off + 158);

        auto g_zz_xyzzzz = cbuffer.data(di_off + 159);

        auto g_zz_xzzzzz = cbuffer.data(di_off + 160);

        auto g_zz_yyyyyy = cbuffer.data(di_off + 161);

        auto g_zz_yyyyyz = cbuffer.data(di_off + 162);

        auto g_zz_yyyyzz = cbuffer.data(di_off + 163);

        auto g_zz_yyyzzz = cbuffer.data(di_off + 164);

        auto g_zz_yyzzzz = cbuffer.data(di_off + 165);

        auto g_zz_yzzzzz = cbuffer.data(di_off + 166);

        auto g_zz_zzzzzz = cbuffer.data(di_off + 167);

        #pragma omp simd aligned(cd_z, g_z_xxxxxx, g_z_xxxxxxz, g_z_xxxxxy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxxzz, g_z_xxxxyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxxzzz, g_z_xxxyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxxzzzz, g_z_xxyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xxzzzzz, g_z_xyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_xzzzzzz, g_z_yyyyyy, g_z_yyyyyyz, g_z_yyyyyz, g_z_yyyyyzz, g_z_yyyyzz, g_z_yyyyzzz, g_z_yyyzzz, g_z_yyyzzzz, g_z_yyzzzz, g_z_yyzzzzz, g_z_yzzzzz, g_z_yzzzzzz, g_z_zzzzzz, g_z_zzzzzzz, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz, g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyzz, g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz, g_zz_zzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zz_xxxxxx[k] = -g_z_xxxxxx[k] * cd_z[k] + g_z_xxxxxxz[k];

            g_zz_xxxxxy[k] = -g_z_xxxxxy[k] * cd_z[k] + g_z_xxxxxyz[k];

            g_zz_xxxxxz[k] = -g_z_xxxxxz[k] * cd_z[k] + g_z_xxxxxzz[k];

            g_zz_xxxxyy[k] = -g_z_xxxxyy[k] * cd_z[k] + g_z_xxxxyyz[k];

            g_zz_xxxxyz[k] = -g_z_xxxxyz[k] * cd_z[k] + g_z_xxxxyzz[k];

            g_zz_xxxxzz[k] = -g_z_xxxxzz[k] * cd_z[k] + g_z_xxxxzzz[k];

            g_zz_xxxyyy[k] = -g_z_xxxyyy[k] * cd_z[k] + g_z_xxxyyyz[k];

            g_zz_xxxyyz[k] = -g_z_xxxyyz[k] * cd_z[k] + g_z_xxxyyzz[k];

            g_zz_xxxyzz[k] = -g_z_xxxyzz[k] * cd_z[k] + g_z_xxxyzzz[k];

            g_zz_xxxzzz[k] = -g_z_xxxzzz[k] * cd_z[k] + g_z_xxxzzzz[k];

            g_zz_xxyyyy[k] = -g_z_xxyyyy[k] * cd_z[k] + g_z_xxyyyyz[k];

            g_zz_xxyyyz[k] = -g_z_xxyyyz[k] * cd_z[k] + g_z_xxyyyzz[k];

            g_zz_xxyyzz[k] = -g_z_xxyyzz[k] * cd_z[k] + g_z_xxyyzzz[k];

            g_zz_xxyzzz[k] = -g_z_xxyzzz[k] * cd_z[k] + g_z_xxyzzzz[k];

            g_zz_xxzzzz[k] = -g_z_xxzzzz[k] * cd_z[k] + g_z_xxzzzzz[k];

            g_zz_xyyyyy[k] = -g_z_xyyyyy[k] * cd_z[k] + g_z_xyyyyyz[k];

            g_zz_xyyyyz[k] = -g_z_xyyyyz[k] * cd_z[k] + g_z_xyyyyzz[k];

            g_zz_xyyyzz[k] = -g_z_xyyyzz[k] * cd_z[k] + g_z_xyyyzzz[k];

            g_zz_xyyzzz[k] = -g_z_xyyzzz[k] * cd_z[k] + g_z_xyyzzzz[k];

            g_zz_xyzzzz[k] = -g_z_xyzzzz[k] * cd_z[k] + g_z_xyzzzzz[k];

            g_zz_xzzzzz[k] = -g_z_xzzzzz[k] * cd_z[k] + g_z_xzzzzzz[k];

            g_zz_yyyyyy[k] = -g_z_yyyyyy[k] * cd_z[k] + g_z_yyyyyyz[k];

            g_zz_yyyyyz[k] = -g_z_yyyyyz[k] * cd_z[k] + g_z_yyyyyzz[k];

            g_zz_yyyyzz[k] = -g_z_yyyyzz[k] * cd_z[k] + g_z_yyyyzzz[k];

            g_zz_yyyzzz[k] = -g_z_yyyzzz[k] * cd_z[k] + g_z_yyyzzzz[k];

            g_zz_yyzzzz[k] = -g_z_yyzzzz[k] * cd_z[k] + g_z_yyzzzzz[k];

            g_zz_yzzzzz[k] = -g_z_yzzzzz[k] * cd_z[k] + g_z_yzzzzzz[k];

            g_zz_zzzzzz[k] = -g_z_zzzzzz[k] * cd_z[k] + g_z_zzzzzzz[k];
        }
    }
}

} // t3ceri namespace


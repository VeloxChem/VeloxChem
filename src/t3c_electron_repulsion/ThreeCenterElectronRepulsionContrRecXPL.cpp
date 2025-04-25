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

#include "ThreeCenterElectronRepulsionContrRecXPL.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xpl(CSimdArray<double>& cbuffer,
                                const size_t idx_xpl,
                                const size_t idx_xsl,
                                const size_t idx_xsm,
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
        /// Set up components of auxilary buffer : SSL

        const auto sl_off = idx_xsl + i * 45;

        auto g_0_xxxxxxxx = cbuffer.data(sl_off + 0);

        auto g_0_xxxxxxxy = cbuffer.data(sl_off + 1);

        auto g_0_xxxxxxxz = cbuffer.data(sl_off + 2);

        auto g_0_xxxxxxyy = cbuffer.data(sl_off + 3);

        auto g_0_xxxxxxyz = cbuffer.data(sl_off + 4);

        auto g_0_xxxxxxzz = cbuffer.data(sl_off + 5);

        auto g_0_xxxxxyyy = cbuffer.data(sl_off + 6);

        auto g_0_xxxxxyyz = cbuffer.data(sl_off + 7);

        auto g_0_xxxxxyzz = cbuffer.data(sl_off + 8);

        auto g_0_xxxxxzzz = cbuffer.data(sl_off + 9);

        auto g_0_xxxxyyyy = cbuffer.data(sl_off + 10);

        auto g_0_xxxxyyyz = cbuffer.data(sl_off + 11);

        auto g_0_xxxxyyzz = cbuffer.data(sl_off + 12);

        auto g_0_xxxxyzzz = cbuffer.data(sl_off + 13);

        auto g_0_xxxxzzzz = cbuffer.data(sl_off + 14);

        auto g_0_xxxyyyyy = cbuffer.data(sl_off + 15);

        auto g_0_xxxyyyyz = cbuffer.data(sl_off + 16);

        auto g_0_xxxyyyzz = cbuffer.data(sl_off + 17);

        auto g_0_xxxyyzzz = cbuffer.data(sl_off + 18);

        auto g_0_xxxyzzzz = cbuffer.data(sl_off + 19);

        auto g_0_xxxzzzzz = cbuffer.data(sl_off + 20);

        auto g_0_xxyyyyyy = cbuffer.data(sl_off + 21);

        auto g_0_xxyyyyyz = cbuffer.data(sl_off + 22);

        auto g_0_xxyyyyzz = cbuffer.data(sl_off + 23);

        auto g_0_xxyyyzzz = cbuffer.data(sl_off + 24);

        auto g_0_xxyyzzzz = cbuffer.data(sl_off + 25);

        auto g_0_xxyzzzzz = cbuffer.data(sl_off + 26);

        auto g_0_xxzzzzzz = cbuffer.data(sl_off + 27);

        auto g_0_xyyyyyyy = cbuffer.data(sl_off + 28);

        auto g_0_xyyyyyyz = cbuffer.data(sl_off + 29);

        auto g_0_xyyyyyzz = cbuffer.data(sl_off + 30);

        auto g_0_xyyyyzzz = cbuffer.data(sl_off + 31);

        auto g_0_xyyyzzzz = cbuffer.data(sl_off + 32);

        auto g_0_xyyzzzzz = cbuffer.data(sl_off + 33);

        auto g_0_xyzzzzzz = cbuffer.data(sl_off + 34);

        auto g_0_xzzzzzzz = cbuffer.data(sl_off + 35);

        auto g_0_yyyyyyyy = cbuffer.data(sl_off + 36);

        auto g_0_yyyyyyyz = cbuffer.data(sl_off + 37);

        auto g_0_yyyyyyzz = cbuffer.data(sl_off + 38);

        auto g_0_yyyyyzzz = cbuffer.data(sl_off + 39);

        auto g_0_yyyyzzzz = cbuffer.data(sl_off + 40);

        auto g_0_yyyzzzzz = cbuffer.data(sl_off + 41);

        auto g_0_yyzzzzzz = cbuffer.data(sl_off + 42);

        auto g_0_yzzzzzzz = cbuffer.data(sl_off + 43);

        auto g_0_zzzzzzzz = cbuffer.data(sl_off + 44);

        /// Set up components of auxilary buffer : SSM

        const auto sm_off = idx_xsm + i * 55;

        auto g_0_xxxxxxxxx = cbuffer.data(sm_off + 0);

        auto g_0_xxxxxxxxy = cbuffer.data(sm_off + 1);

        auto g_0_xxxxxxxxz = cbuffer.data(sm_off + 2);

        auto g_0_xxxxxxxyy = cbuffer.data(sm_off + 3);

        auto g_0_xxxxxxxyz = cbuffer.data(sm_off + 4);

        auto g_0_xxxxxxxzz = cbuffer.data(sm_off + 5);

        auto g_0_xxxxxxyyy = cbuffer.data(sm_off + 6);

        auto g_0_xxxxxxyyz = cbuffer.data(sm_off + 7);

        auto g_0_xxxxxxyzz = cbuffer.data(sm_off + 8);

        auto g_0_xxxxxxzzz = cbuffer.data(sm_off + 9);

        auto g_0_xxxxxyyyy = cbuffer.data(sm_off + 10);

        auto g_0_xxxxxyyyz = cbuffer.data(sm_off + 11);

        auto g_0_xxxxxyyzz = cbuffer.data(sm_off + 12);

        auto g_0_xxxxxyzzz = cbuffer.data(sm_off + 13);

        auto g_0_xxxxxzzzz = cbuffer.data(sm_off + 14);

        auto g_0_xxxxyyyyy = cbuffer.data(sm_off + 15);

        auto g_0_xxxxyyyyz = cbuffer.data(sm_off + 16);

        auto g_0_xxxxyyyzz = cbuffer.data(sm_off + 17);

        auto g_0_xxxxyyzzz = cbuffer.data(sm_off + 18);

        auto g_0_xxxxyzzzz = cbuffer.data(sm_off + 19);

        auto g_0_xxxxzzzzz = cbuffer.data(sm_off + 20);

        auto g_0_xxxyyyyyy = cbuffer.data(sm_off + 21);

        auto g_0_xxxyyyyyz = cbuffer.data(sm_off + 22);

        auto g_0_xxxyyyyzz = cbuffer.data(sm_off + 23);

        auto g_0_xxxyyyzzz = cbuffer.data(sm_off + 24);

        auto g_0_xxxyyzzzz = cbuffer.data(sm_off + 25);

        auto g_0_xxxyzzzzz = cbuffer.data(sm_off + 26);

        auto g_0_xxxzzzzzz = cbuffer.data(sm_off + 27);

        auto g_0_xxyyyyyyy = cbuffer.data(sm_off + 28);

        auto g_0_xxyyyyyyz = cbuffer.data(sm_off + 29);

        auto g_0_xxyyyyyzz = cbuffer.data(sm_off + 30);

        auto g_0_xxyyyyzzz = cbuffer.data(sm_off + 31);

        auto g_0_xxyyyzzzz = cbuffer.data(sm_off + 32);

        auto g_0_xxyyzzzzz = cbuffer.data(sm_off + 33);

        auto g_0_xxyzzzzzz = cbuffer.data(sm_off + 34);

        auto g_0_xxzzzzzzz = cbuffer.data(sm_off + 35);

        auto g_0_xyyyyyyyy = cbuffer.data(sm_off + 36);

        auto g_0_xyyyyyyyz = cbuffer.data(sm_off + 37);

        auto g_0_xyyyyyyzz = cbuffer.data(sm_off + 38);

        auto g_0_xyyyyyzzz = cbuffer.data(sm_off + 39);

        auto g_0_xyyyyzzzz = cbuffer.data(sm_off + 40);

        auto g_0_xyyyzzzzz = cbuffer.data(sm_off + 41);

        auto g_0_xyyzzzzzz = cbuffer.data(sm_off + 42);

        auto g_0_xyzzzzzzz = cbuffer.data(sm_off + 43);

        auto g_0_xzzzzzzzz = cbuffer.data(sm_off + 44);

        auto g_0_yyyyyyyyy = cbuffer.data(sm_off + 45);

        auto g_0_yyyyyyyyz = cbuffer.data(sm_off + 46);

        auto g_0_yyyyyyyzz = cbuffer.data(sm_off + 47);

        auto g_0_yyyyyyzzz = cbuffer.data(sm_off + 48);

        auto g_0_yyyyyzzzz = cbuffer.data(sm_off + 49);

        auto g_0_yyyyzzzzz = cbuffer.data(sm_off + 50);

        auto g_0_yyyzzzzzz = cbuffer.data(sm_off + 51);

        auto g_0_yyzzzzzzz = cbuffer.data(sm_off + 52);

        auto g_0_yzzzzzzzz = cbuffer.data(sm_off + 53);

        auto g_0_zzzzzzzzz = cbuffer.data(sm_off + 54);

        /// set up bra offset for contr_buffer_xpl

        const auto pl_off = idx_xpl + i * 135;

        /// Set up 0-45 components of targeted buffer : cbuffer.data(

        auto g_x_xxxxxxxx = cbuffer.data(pl_off + 0);

        auto g_x_xxxxxxxy = cbuffer.data(pl_off + 1);

        auto g_x_xxxxxxxz = cbuffer.data(pl_off + 2);

        auto g_x_xxxxxxyy = cbuffer.data(pl_off + 3);

        auto g_x_xxxxxxyz = cbuffer.data(pl_off + 4);

        auto g_x_xxxxxxzz = cbuffer.data(pl_off + 5);

        auto g_x_xxxxxyyy = cbuffer.data(pl_off + 6);

        auto g_x_xxxxxyyz = cbuffer.data(pl_off + 7);

        auto g_x_xxxxxyzz = cbuffer.data(pl_off + 8);

        auto g_x_xxxxxzzz = cbuffer.data(pl_off + 9);

        auto g_x_xxxxyyyy = cbuffer.data(pl_off + 10);

        auto g_x_xxxxyyyz = cbuffer.data(pl_off + 11);

        auto g_x_xxxxyyzz = cbuffer.data(pl_off + 12);

        auto g_x_xxxxyzzz = cbuffer.data(pl_off + 13);

        auto g_x_xxxxzzzz = cbuffer.data(pl_off + 14);

        auto g_x_xxxyyyyy = cbuffer.data(pl_off + 15);

        auto g_x_xxxyyyyz = cbuffer.data(pl_off + 16);

        auto g_x_xxxyyyzz = cbuffer.data(pl_off + 17);

        auto g_x_xxxyyzzz = cbuffer.data(pl_off + 18);

        auto g_x_xxxyzzzz = cbuffer.data(pl_off + 19);

        auto g_x_xxxzzzzz = cbuffer.data(pl_off + 20);

        auto g_x_xxyyyyyy = cbuffer.data(pl_off + 21);

        auto g_x_xxyyyyyz = cbuffer.data(pl_off + 22);

        auto g_x_xxyyyyzz = cbuffer.data(pl_off + 23);

        auto g_x_xxyyyzzz = cbuffer.data(pl_off + 24);

        auto g_x_xxyyzzzz = cbuffer.data(pl_off + 25);

        auto g_x_xxyzzzzz = cbuffer.data(pl_off + 26);

        auto g_x_xxzzzzzz = cbuffer.data(pl_off + 27);

        auto g_x_xyyyyyyy = cbuffer.data(pl_off + 28);

        auto g_x_xyyyyyyz = cbuffer.data(pl_off + 29);

        auto g_x_xyyyyyzz = cbuffer.data(pl_off + 30);

        auto g_x_xyyyyzzz = cbuffer.data(pl_off + 31);

        auto g_x_xyyyzzzz = cbuffer.data(pl_off + 32);

        auto g_x_xyyzzzzz = cbuffer.data(pl_off + 33);

        auto g_x_xyzzzzzz = cbuffer.data(pl_off + 34);

        auto g_x_xzzzzzzz = cbuffer.data(pl_off + 35);

        auto g_x_yyyyyyyy = cbuffer.data(pl_off + 36);

        auto g_x_yyyyyyyz = cbuffer.data(pl_off + 37);

        auto g_x_yyyyyyzz = cbuffer.data(pl_off + 38);

        auto g_x_yyyyyzzz = cbuffer.data(pl_off + 39);

        auto g_x_yyyyzzzz = cbuffer.data(pl_off + 40);

        auto g_x_yyyzzzzz = cbuffer.data(pl_off + 41);

        auto g_x_yyzzzzzz = cbuffer.data(pl_off + 42);

        auto g_x_yzzzzzzz = cbuffer.data(pl_off + 43);

        auto g_x_zzzzzzzz = cbuffer.data(pl_off + 44);

        #pragma omp simd aligned(cd_x, g_0_xxxxxxxx, g_0_xxxxxxxxx, g_0_xxxxxxxxy, g_0_xxxxxxxxz, g_0_xxxxxxxy, g_0_xxxxxxxyy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxxzz, g_0_xxxxxxyy, g_0_xxxxxxyyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxxzzz, g_0_xxxxxyyy, g_0_xxxxxyyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxxzzzz, g_0_xxxxyyyy, g_0_xxxxyyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxxzzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxxzzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xxzzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_xzzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyz, g_0_yyyyyyzz, g_0_yyyyyzzz, g_0_yyyyzzzz, g_0_yyyzzzzz, g_0_yyzzzzzz, g_0_yzzzzzzz, g_0_zzzzzzzz, g_x_xxxxxxxx, g_x_xxxxxxxy, g_x_xxxxxxxz, g_x_xxxxxxyy, g_x_xxxxxxyz, g_x_xxxxxxzz, g_x_xxxxxyyy, g_x_xxxxxyyz, g_x_xxxxxyzz, g_x_xxxxxzzz, g_x_xxxxyyyy, g_x_xxxxyyyz, g_x_xxxxyyzz, g_x_xxxxyzzz, g_x_xxxxzzzz, g_x_xxxyyyyy, g_x_xxxyyyyz, g_x_xxxyyyzz, g_x_xxxyyzzz, g_x_xxxyzzzz, g_x_xxxzzzzz, g_x_xxyyyyyy, g_x_xxyyyyyz, g_x_xxyyyyzz, g_x_xxyyyzzz, g_x_xxyyzzzz, g_x_xxyzzzzz, g_x_xxzzzzzz, g_x_xyyyyyyy, g_x_xyyyyyyz, g_x_xyyyyyzz, g_x_xyyyyzzz, g_x_xyyyzzzz, g_x_xyyzzzzz, g_x_xyzzzzzz, g_x_xzzzzzzz, g_x_yyyyyyyy, g_x_yyyyyyyz, g_x_yyyyyyzz, g_x_yyyyyzzz, g_x_yyyyzzzz, g_x_yyyzzzzz, g_x_yyzzzzzz, g_x_yzzzzzzz, g_x_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_x[k] + g_0_xxxxxxxxx[k];

            g_x_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_x[k] + g_0_xxxxxxxxy[k];

            g_x_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_x[k] + g_0_xxxxxxxxz[k];

            g_x_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_x[k] + g_0_xxxxxxxyy[k];

            g_x_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_x[k] + g_0_xxxxxxxyz[k];

            g_x_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_x[k] + g_0_xxxxxxxzz[k];

            g_x_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_x[k] + g_0_xxxxxxyyy[k];

            g_x_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_x[k] + g_0_xxxxxxyyz[k];

            g_x_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_x[k] + g_0_xxxxxxyzz[k];

            g_x_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_x[k] + g_0_xxxxxxzzz[k];

            g_x_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_x[k] + g_0_xxxxxyyyy[k];

            g_x_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_x[k] + g_0_xxxxxyyyz[k];

            g_x_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_x[k] + g_0_xxxxxyyzz[k];

            g_x_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_x[k] + g_0_xxxxxyzzz[k];

            g_x_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_x[k] + g_0_xxxxxzzzz[k];

            g_x_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_x[k] + g_0_xxxxyyyyy[k];

            g_x_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_x[k] + g_0_xxxxyyyyz[k];

            g_x_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_x[k] + g_0_xxxxyyyzz[k];

            g_x_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_x[k] + g_0_xxxxyyzzz[k];

            g_x_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_x[k] + g_0_xxxxyzzzz[k];

            g_x_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_x[k] + g_0_xxxxzzzzz[k];

            g_x_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_x[k] + g_0_xxxyyyyyy[k];

            g_x_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_x[k] + g_0_xxxyyyyyz[k];

            g_x_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_x[k] + g_0_xxxyyyyzz[k];

            g_x_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_x[k] + g_0_xxxyyyzzz[k];

            g_x_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_x[k] + g_0_xxxyyzzzz[k];

            g_x_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_x[k] + g_0_xxxyzzzzz[k];

            g_x_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_x[k] + g_0_xxxzzzzzz[k];

            g_x_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_x[k] + g_0_xxyyyyyyy[k];

            g_x_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_x[k] + g_0_xxyyyyyyz[k];

            g_x_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_x[k] + g_0_xxyyyyyzz[k];

            g_x_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_x[k] + g_0_xxyyyyzzz[k];

            g_x_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_x[k] + g_0_xxyyyzzzz[k];

            g_x_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_x[k] + g_0_xxyyzzzzz[k];

            g_x_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_x[k] + g_0_xxyzzzzzz[k];

            g_x_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_x[k] + g_0_xxzzzzzzz[k];

            g_x_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_x[k] + g_0_xyyyyyyyy[k];

            g_x_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_x[k] + g_0_xyyyyyyyz[k];

            g_x_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_x[k] + g_0_xyyyyyyzz[k];

            g_x_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_x[k] + g_0_xyyyyyzzz[k];

            g_x_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_x[k] + g_0_xyyyyzzzz[k];

            g_x_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_x[k] + g_0_xyyyzzzzz[k];

            g_x_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_x[k] + g_0_xyyzzzzzz[k];

            g_x_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_x[k] + g_0_xyzzzzzzz[k];

            g_x_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_x[k] + g_0_xzzzzzzzz[k];
        }

        /// Set up 45-90 components of targeted buffer : cbuffer.data(

        auto g_y_xxxxxxxx = cbuffer.data(pl_off + 45);

        auto g_y_xxxxxxxy = cbuffer.data(pl_off + 46);

        auto g_y_xxxxxxxz = cbuffer.data(pl_off + 47);

        auto g_y_xxxxxxyy = cbuffer.data(pl_off + 48);

        auto g_y_xxxxxxyz = cbuffer.data(pl_off + 49);

        auto g_y_xxxxxxzz = cbuffer.data(pl_off + 50);

        auto g_y_xxxxxyyy = cbuffer.data(pl_off + 51);

        auto g_y_xxxxxyyz = cbuffer.data(pl_off + 52);

        auto g_y_xxxxxyzz = cbuffer.data(pl_off + 53);

        auto g_y_xxxxxzzz = cbuffer.data(pl_off + 54);

        auto g_y_xxxxyyyy = cbuffer.data(pl_off + 55);

        auto g_y_xxxxyyyz = cbuffer.data(pl_off + 56);

        auto g_y_xxxxyyzz = cbuffer.data(pl_off + 57);

        auto g_y_xxxxyzzz = cbuffer.data(pl_off + 58);

        auto g_y_xxxxzzzz = cbuffer.data(pl_off + 59);

        auto g_y_xxxyyyyy = cbuffer.data(pl_off + 60);

        auto g_y_xxxyyyyz = cbuffer.data(pl_off + 61);

        auto g_y_xxxyyyzz = cbuffer.data(pl_off + 62);

        auto g_y_xxxyyzzz = cbuffer.data(pl_off + 63);

        auto g_y_xxxyzzzz = cbuffer.data(pl_off + 64);

        auto g_y_xxxzzzzz = cbuffer.data(pl_off + 65);

        auto g_y_xxyyyyyy = cbuffer.data(pl_off + 66);

        auto g_y_xxyyyyyz = cbuffer.data(pl_off + 67);

        auto g_y_xxyyyyzz = cbuffer.data(pl_off + 68);

        auto g_y_xxyyyzzz = cbuffer.data(pl_off + 69);

        auto g_y_xxyyzzzz = cbuffer.data(pl_off + 70);

        auto g_y_xxyzzzzz = cbuffer.data(pl_off + 71);

        auto g_y_xxzzzzzz = cbuffer.data(pl_off + 72);

        auto g_y_xyyyyyyy = cbuffer.data(pl_off + 73);

        auto g_y_xyyyyyyz = cbuffer.data(pl_off + 74);

        auto g_y_xyyyyyzz = cbuffer.data(pl_off + 75);

        auto g_y_xyyyyzzz = cbuffer.data(pl_off + 76);

        auto g_y_xyyyzzzz = cbuffer.data(pl_off + 77);

        auto g_y_xyyzzzzz = cbuffer.data(pl_off + 78);

        auto g_y_xyzzzzzz = cbuffer.data(pl_off + 79);

        auto g_y_xzzzzzzz = cbuffer.data(pl_off + 80);

        auto g_y_yyyyyyyy = cbuffer.data(pl_off + 81);

        auto g_y_yyyyyyyz = cbuffer.data(pl_off + 82);

        auto g_y_yyyyyyzz = cbuffer.data(pl_off + 83);

        auto g_y_yyyyyzzz = cbuffer.data(pl_off + 84);

        auto g_y_yyyyzzzz = cbuffer.data(pl_off + 85);

        auto g_y_yyyzzzzz = cbuffer.data(pl_off + 86);

        auto g_y_yyzzzzzz = cbuffer.data(pl_off + 87);

        auto g_y_yzzzzzzz = cbuffer.data(pl_off + 88);

        auto g_y_zzzzzzzz = cbuffer.data(pl_off + 89);

        #pragma omp simd aligned(cd_y, g_0_xxxxxxxx, g_0_xxxxxxxxy, g_0_xxxxxxxy, g_0_xxxxxxxyy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxyy, g_0_xxxxxxyyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxyyy, g_0_xxxxxyyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxyyyy, g_0_xxxxyyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyyy, g_0_yyyyyyyyz, g_0_yyyyyyyz, g_0_yyyyyyyzz, g_0_yyyyyyzz, g_0_yyyyyyzzz, g_0_yyyyyzzz, g_0_yyyyyzzzz, g_0_yyyyzzzz, g_0_yyyyzzzzz, g_0_yyyzzzzz, g_0_yyyzzzzzz, g_0_yyzzzzzz, g_0_yyzzzzzzz, g_0_yzzzzzzz, g_0_yzzzzzzzz, g_0_zzzzzzzz, g_y_xxxxxxxx, g_y_xxxxxxxy, g_y_xxxxxxxz, g_y_xxxxxxyy, g_y_xxxxxxyz, g_y_xxxxxxzz, g_y_xxxxxyyy, g_y_xxxxxyyz, g_y_xxxxxyzz, g_y_xxxxxzzz, g_y_xxxxyyyy, g_y_xxxxyyyz, g_y_xxxxyyzz, g_y_xxxxyzzz, g_y_xxxxzzzz, g_y_xxxyyyyy, g_y_xxxyyyyz, g_y_xxxyyyzz, g_y_xxxyyzzz, g_y_xxxyzzzz, g_y_xxxzzzzz, g_y_xxyyyyyy, g_y_xxyyyyyz, g_y_xxyyyyzz, g_y_xxyyyzzz, g_y_xxyyzzzz, g_y_xxyzzzzz, g_y_xxzzzzzz, g_y_xyyyyyyy, g_y_xyyyyyyz, g_y_xyyyyyzz, g_y_xyyyyzzz, g_y_xyyyzzzz, g_y_xyyzzzzz, g_y_xyzzzzzz, g_y_xzzzzzzz, g_y_yyyyyyyy, g_y_yyyyyyyz, g_y_yyyyyyzz, g_y_yyyyyzzz, g_y_yyyyzzzz, g_y_yyyzzzzz, g_y_yyzzzzzz, g_y_yzzzzzzz, g_y_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_y[k] + g_0_xxxxxxxxy[k];

            g_y_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_y[k] + g_0_xxxxxxxyy[k];

            g_y_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_y[k] + g_0_xxxxxxxyz[k];

            g_y_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_y[k] + g_0_xxxxxxyyy[k];

            g_y_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_y[k] + g_0_xxxxxxyyz[k];

            g_y_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_y[k] + g_0_xxxxxxyzz[k];

            g_y_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_y[k] + g_0_xxxxxyyyy[k];

            g_y_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_y[k] + g_0_xxxxxyyyz[k];

            g_y_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_y[k] + g_0_xxxxxyyzz[k];

            g_y_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_y[k] + g_0_xxxxxyzzz[k];

            g_y_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_y[k] + g_0_xxxxyyyyy[k];

            g_y_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_y[k] + g_0_xxxxyyyyz[k];

            g_y_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_y[k] + g_0_xxxxyyyzz[k];

            g_y_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_y[k] + g_0_xxxxyyzzz[k];

            g_y_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_y[k] + g_0_xxxxyzzzz[k];

            g_y_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_y[k] + g_0_xxxyyyyyy[k];

            g_y_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_y[k] + g_0_xxxyyyyyz[k];

            g_y_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_y[k] + g_0_xxxyyyyzz[k];

            g_y_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_y[k] + g_0_xxxyyyzzz[k];

            g_y_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_y[k] + g_0_xxxyyzzzz[k];

            g_y_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_y[k] + g_0_xxxyzzzzz[k];

            g_y_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_y[k] + g_0_xxyyyyyyy[k];

            g_y_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_y[k] + g_0_xxyyyyyyz[k];

            g_y_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_y[k] + g_0_xxyyyyyzz[k];

            g_y_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_y[k] + g_0_xxyyyyzzz[k];

            g_y_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_y[k] + g_0_xxyyyzzzz[k];

            g_y_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_y[k] + g_0_xxyyzzzzz[k];

            g_y_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_y[k] + g_0_xxyzzzzzz[k];

            g_y_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_y[k] + g_0_xyyyyyyyy[k];

            g_y_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_y[k] + g_0_xyyyyyyyz[k];

            g_y_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_y[k] + g_0_xyyyyyyzz[k];

            g_y_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_y[k] + g_0_xyyyyyzzz[k];

            g_y_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_y[k] + g_0_xyyyyzzzz[k];

            g_y_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_y[k] + g_0_xyyyzzzzz[k];

            g_y_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_y[k] + g_0_xyyzzzzzz[k];

            g_y_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_y[k] + g_0_xyzzzzzzz[k];

            g_y_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_y[k] + g_0_yyyyyyyyy[k];

            g_y_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_y[k] + g_0_yyyyyyyyz[k];

            g_y_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_y[k] + g_0_yyyyyyyzz[k];

            g_y_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_y[k] + g_0_yyyyyyzzz[k];

            g_y_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_y[k] + g_0_yyyyyzzzz[k];

            g_y_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_y[k] + g_0_yyyyzzzzz[k];

            g_y_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_y[k] + g_0_yyyzzzzzz[k];

            g_y_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_y[k] + g_0_yyzzzzzzz[k];

            g_y_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_y[k] + g_0_yzzzzzzzz[k];
        }

        /// Set up 90-135 components of targeted buffer : cbuffer.data(

        auto g_z_xxxxxxxx = cbuffer.data(pl_off + 90);

        auto g_z_xxxxxxxy = cbuffer.data(pl_off + 91);

        auto g_z_xxxxxxxz = cbuffer.data(pl_off + 92);

        auto g_z_xxxxxxyy = cbuffer.data(pl_off + 93);

        auto g_z_xxxxxxyz = cbuffer.data(pl_off + 94);

        auto g_z_xxxxxxzz = cbuffer.data(pl_off + 95);

        auto g_z_xxxxxyyy = cbuffer.data(pl_off + 96);

        auto g_z_xxxxxyyz = cbuffer.data(pl_off + 97);

        auto g_z_xxxxxyzz = cbuffer.data(pl_off + 98);

        auto g_z_xxxxxzzz = cbuffer.data(pl_off + 99);

        auto g_z_xxxxyyyy = cbuffer.data(pl_off + 100);

        auto g_z_xxxxyyyz = cbuffer.data(pl_off + 101);

        auto g_z_xxxxyyzz = cbuffer.data(pl_off + 102);

        auto g_z_xxxxyzzz = cbuffer.data(pl_off + 103);

        auto g_z_xxxxzzzz = cbuffer.data(pl_off + 104);

        auto g_z_xxxyyyyy = cbuffer.data(pl_off + 105);

        auto g_z_xxxyyyyz = cbuffer.data(pl_off + 106);

        auto g_z_xxxyyyzz = cbuffer.data(pl_off + 107);

        auto g_z_xxxyyzzz = cbuffer.data(pl_off + 108);

        auto g_z_xxxyzzzz = cbuffer.data(pl_off + 109);

        auto g_z_xxxzzzzz = cbuffer.data(pl_off + 110);

        auto g_z_xxyyyyyy = cbuffer.data(pl_off + 111);

        auto g_z_xxyyyyyz = cbuffer.data(pl_off + 112);

        auto g_z_xxyyyyzz = cbuffer.data(pl_off + 113);

        auto g_z_xxyyyzzz = cbuffer.data(pl_off + 114);

        auto g_z_xxyyzzzz = cbuffer.data(pl_off + 115);

        auto g_z_xxyzzzzz = cbuffer.data(pl_off + 116);

        auto g_z_xxzzzzzz = cbuffer.data(pl_off + 117);

        auto g_z_xyyyyyyy = cbuffer.data(pl_off + 118);

        auto g_z_xyyyyyyz = cbuffer.data(pl_off + 119);

        auto g_z_xyyyyyzz = cbuffer.data(pl_off + 120);

        auto g_z_xyyyyzzz = cbuffer.data(pl_off + 121);

        auto g_z_xyyyzzzz = cbuffer.data(pl_off + 122);

        auto g_z_xyyzzzzz = cbuffer.data(pl_off + 123);

        auto g_z_xyzzzzzz = cbuffer.data(pl_off + 124);

        auto g_z_xzzzzzzz = cbuffer.data(pl_off + 125);

        auto g_z_yyyyyyyy = cbuffer.data(pl_off + 126);

        auto g_z_yyyyyyyz = cbuffer.data(pl_off + 127);

        auto g_z_yyyyyyzz = cbuffer.data(pl_off + 128);

        auto g_z_yyyyyzzz = cbuffer.data(pl_off + 129);

        auto g_z_yyyyzzzz = cbuffer.data(pl_off + 130);

        auto g_z_yyyzzzzz = cbuffer.data(pl_off + 131);

        auto g_z_yyzzzzzz = cbuffer.data(pl_off + 132);

        auto g_z_yzzzzzzz = cbuffer.data(pl_off + 133);

        auto g_z_zzzzzzzz = cbuffer.data(pl_off + 134);

        #pragma omp simd aligned(cd_z, g_0_xxxxxxxx, g_0_xxxxxxxxz, g_0_xxxxxxxy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxxzz, g_0_xxxxxxyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxxzzz, g_0_xxxxxyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxxzzzz, g_0_xxxxyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxxzzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxxzzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xxzzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_xzzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyyz, g_0_yyyyyyyz, g_0_yyyyyyyzz, g_0_yyyyyyzz, g_0_yyyyyyzzz, g_0_yyyyyzzz, g_0_yyyyyzzzz, g_0_yyyyzzzz, g_0_yyyyzzzzz, g_0_yyyzzzzz, g_0_yyyzzzzzz, g_0_yyzzzzzz, g_0_yyzzzzzzz, g_0_yzzzzzzz, g_0_yzzzzzzzz, g_0_zzzzzzzz, g_0_zzzzzzzzz, g_z_xxxxxxxx, g_z_xxxxxxxy, g_z_xxxxxxxz, g_z_xxxxxxyy, g_z_xxxxxxyz, g_z_xxxxxxzz, g_z_xxxxxyyy, g_z_xxxxxyyz, g_z_xxxxxyzz, g_z_xxxxxzzz, g_z_xxxxyyyy, g_z_xxxxyyyz, g_z_xxxxyyzz, g_z_xxxxyzzz, g_z_xxxxzzzz, g_z_xxxyyyyy, g_z_xxxyyyyz, g_z_xxxyyyzz, g_z_xxxyyzzz, g_z_xxxyzzzz, g_z_xxxzzzzz, g_z_xxyyyyyy, g_z_xxyyyyyz, g_z_xxyyyyzz, g_z_xxyyyzzz, g_z_xxyyzzzz, g_z_xxyzzzzz, g_z_xxzzzzzz, g_z_xyyyyyyy, g_z_xyyyyyyz, g_z_xyyyyyzz, g_z_xyyyyzzz, g_z_xyyyzzzz, g_z_xyyzzzzz, g_z_xyzzzzzz, g_z_xzzzzzzz, g_z_yyyyyyyy, g_z_yyyyyyyz, g_z_yyyyyyzz, g_z_yyyyyzzz, g_z_yyyyzzzz, g_z_yyyzzzzz, g_z_yyzzzzzz, g_z_yzzzzzzz, g_z_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_z[k] + g_0_xxxxxxxxz[k];

            g_z_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_z[k] + g_0_xxxxxxxyz[k];

            g_z_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_z[k] + g_0_xxxxxxxzz[k];

            g_z_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_z[k] + g_0_xxxxxxyyz[k];

            g_z_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_z[k] + g_0_xxxxxxyzz[k];

            g_z_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_z[k] + g_0_xxxxxxzzz[k];

            g_z_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_z[k] + g_0_xxxxxyyyz[k];

            g_z_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_z[k] + g_0_xxxxxyyzz[k];

            g_z_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_z[k] + g_0_xxxxxyzzz[k];

            g_z_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_z[k] + g_0_xxxxxzzzz[k];

            g_z_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_z[k] + g_0_xxxxyyyyz[k];

            g_z_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_z[k] + g_0_xxxxyyyzz[k];

            g_z_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_z[k] + g_0_xxxxyyzzz[k];

            g_z_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_z[k] + g_0_xxxxyzzzz[k];

            g_z_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_z[k] + g_0_xxxxzzzzz[k];

            g_z_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_z[k] + g_0_xxxyyyyyz[k];

            g_z_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_z[k] + g_0_xxxyyyyzz[k];

            g_z_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_z[k] + g_0_xxxyyyzzz[k];

            g_z_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_z[k] + g_0_xxxyyzzzz[k];

            g_z_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_z[k] + g_0_xxxyzzzzz[k];

            g_z_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_z[k] + g_0_xxxzzzzzz[k];

            g_z_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_z[k] + g_0_xxyyyyyyz[k];

            g_z_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_z[k] + g_0_xxyyyyyzz[k];

            g_z_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_z[k] + g_0_xxyyyyzzz[k];

            g_z_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_z[k] + g_0_xxyyyzzzz[k];

            g_z_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_z[k] + g_0_xxyyzzzzz[k];

            g_z_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_z[k] + g_0_xxyzzzzzz[k];

            g_z_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_z[k] + g_0_xxzzzzzzz[k];

            g_z_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_z[k] + g_0_xyyyyyyyz[k];

            g_z_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_z[k] + g_0_xyyyyyyzz[k];

            g_z_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_z[k] + g_0_xyyyyyzzz[k];

            g_z_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_z[k] + g_0_xyyyyzzzz[k];

            g_z_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_z[k] + g_0_xyyyzzzzz[k];

            g_z_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_z[k] + g_0_xyyzzzzzz[k];

            g_z_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_z[k] + g_0_xyzzzzzzz[k];

            g_z_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_z[k] + g_0_xzzzzzzzz[k];

            g_z_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_z[k] + g_0_yyyyyyyyz[k];

            g_z_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_z[k] + g_0_yyyyyyyzz[k];

            g_z_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_z[k] + g_0_yyyyyyzzz[k];

            g_z_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_z[k] + g_0_yyyyyzzzz[k];

            g_z_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_z[k] + g_0_yyyyzzzzz[k];

            g_z_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_z[k] + g_0_yyyzzzzzz[k];

            g_z_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_z[k] + g_0_yyzzzzzzz[k];

            g_z_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_z[k] + g_0_yzzzzzzzz[k];

            g_z_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_z[k] + g_0_zzzzzzzzz[k];
        }
    }
}

} // t3ceri namespace


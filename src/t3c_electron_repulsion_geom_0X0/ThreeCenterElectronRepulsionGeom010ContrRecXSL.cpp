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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXSL.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xsl(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xsl,
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

        /// set up bra offset for contr_buffer_xxsl

        const auto sl_geom_10_off = idx_geom_10_xsl + i * 45;

        /// Set up 0-45 components of targeted buffer : cbuffer.data(

        auto g_x_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 0 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_0_xxxxxxxx, g_0_xxxxxxxxx, g_0_xxxxxxxxy, g_0_xxxxxxxxz, g_0_xxxxxxxy, g_0_xxxxxxxyy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxxzz, g_0_xxxxxxyy, g_0_xxxxxxyyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxxzzz, g_0_xxxxxyyy, g_0_xxxxxyyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxxzzzz, g_0_xxxxyyyy, g_0_xxxxyyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxxzzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxxzzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xxzzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_xzzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyz, g_0_yyyyyyzz, g_0_yyyyyzzz, g_0_yyyyzzzz, g_0_yyyzzzzz, g_0_yyzzzzzz, g_0_yzzzzzzz, g_0_zzzzzzzz, g_x_0_0_xxxxxxxx, g_x_0_0_xxxxxxxy, g_x_0_0_xxxxxxxz, g_x_0_0_xxxxxxyy, g_x_0_0_xxxxxxyz, g_x_0_0_xxxxxxzz, g_x_0_0_xxxxxyyy, g_x_0_0_xxxxxyyz, g_x_0_0_xxxxxyzz, g_x_0_0_xxxxxzzz, g_x_0_0_xxxxyyyy, g_x_0_0_xxxxyyyz, g_x_0_0_xxxxyyzz, g_x_0_0_xxxxyzzz, g_x_0_0_xxxxzzzz, g_x_0_0_xxxyyyyy, g_x_0_0_xxxyyyyz, g_x_0_0_xxxyyyzz, g_x_0_0_xxxyyzzz, g_x_0_0_xxxyzzzz, g_x_0_0_xxxzzzzz, g_x_0_0_xxyyyyyy, g_x_0_0_xxyyyyyz, g_x_0_0_xxyyyyzz, g_x_0_0_xxyyyzzz, g_x_0_0_xxyyzzzz, g_x_0_0_xxyzzzzz, g_x_0_0_xxzzzzzz, g_x_0_0_xyyyyyyy, g_x_0_0_xyyyyyyz, g_x_0_0_xyyyyyzz, g_x_0_0_xyyyyzzz, g_x_0_0_xyyyzzzz, g_x_0_0_xyyzzzzz, g_x_0_0_xyzzzzzz, g_x_0_0_xzzzzzzz, g_x_0_0_yyyyyyyy, g_x_0_0_yyyyyyyz, g_x_0_0_yyyyyyzz, g_x_0_0_yyyyyzzz, g_x_0_0_yyyyzzzz, g_x_0_0_yyyzzzzz, g_x_0_0_yyzzzzzz, g_x_0_0_yzzzzzzz, g_x_0_0_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_0_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_x[k] + g_0_xxxxxxxxx[k];

            g_x_0_0_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_x[k] + g_0_xxxxxxxxy[k];

            g_x_0_0_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_x[k] + g_0_xxxxxxxxz[k];

            g_x_0_0_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_x[k] + g_0_xxxxxxxyy[k];

            g_x_0_0_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_x[k] + g_0_xxxxxxxyz[k];

            g_x_0_0_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_x[k] + g_0_xxxxxxxzz[k];

            g_x_0_0_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_x[k] + g_0_xxxxxxyyy[k];

            g_x_0_0_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_x[k] + g_0_xxxxxxyyz[k];

            g_x_0_0_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_x[k] + g_0_xxxxxxyzz[k];

            g_x_0_0_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_x[k] + g_0_xxxxxxzzz[k];

            g_x_0_0_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_x[k] + g_0_xxxxxyyyy[k];

            g_x_0_0_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_x[k] + g_0_xxxxxyyyz[k];

            g_x_0_0_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_x[k] + g_0_xxxxxyyzz[k];

            g_x_0_0_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_x[k] + g_0_xxxxxyzzz[k];

            g_x_0_0_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_x[k] + g_0_xxxxxzzzz[k];

            g_x_0_0_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_x[k] + g_0_xxxxyyyyy[k];

            g_x_0_0_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_x[k] + g_0_xxxxyyyyz[k];

            g_x_0_0_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_x[k] + g_0_xxxxyyyzz[k];

            g_x_0_0_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_x[k] + g_0_xxxxyyzzz[k];

            g_x_0_0_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_x[k] + g_0_xxxxyzzzz[k];

            g_x_0_0_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_x[k] + g_0_xxxxzzzzz[k];

            g_x_0_0_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_x[k] + g_0_xxxyyyyyy[k];

            g_x_0_0_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_x[k] + g_0_xxxyyyyyz[k];

            g_x_0_0_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_x[k] + g_0_xxxyyyyzz[k];

            g_x_0_0_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_x[k] + g_0_xxxyyyzzz[k];

            g_x_0_0_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_x[k] + g_0_xxxyyzzzz[k];

            g_x_0_0_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_x[k] + g_0_xxxyzzzzz[k];

            g_x_0_0_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_x[k] + g_0_xxxzzzzzz[k];

            g_x_0_0_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_x[k] + g_0_xxyyyyyyy[k];

            g_x_0_0_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_x[k] + g_0_xxyyyyyyz[k];

            g_x_0_0_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_x[k] + g_0_xxyyyyyzz[k];

            g_x_0_0_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_x[k] + g_0_xxyyyyzzz[k];

            g_x_0_0_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_x[k] + g_0_xxyyyzzzz[k];

            g_x_0_0_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_x[k] + g_0_xxyyzzzzz[k];

            g_x_0_0_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_x[k] + g_0_xxyzzzzzz[k];

            g_x_0_0_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_x[k] + g_0_xxzzzzzzz[k];

            g_x_0_0_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_x[k] + g_0_xyyyyyyyy[k];

            g_x_0_0_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_x[k] + g_0_xyyyyyyyz[k];

            g_x_0_0_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_x[k] + g_0_xyyyyyyzz[k];

            g_x_0_0_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_x[k] + g_0_xyyyyyzzz[k];

            g_x_0_0_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_x[k] + g_0_xyyyyzzzz[k];

            g_x_0_0_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_x[k] + g_0_xyyyzzzzz[k];

            g_x_0_0_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_x[k] + g_0_xyyzzzzzz[k];

            g_x_0_0_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_x[k] + g_0_xyzzzzzzz[k];

            g_x_0_0_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_x[k] + g_0_xzzzzzzzz[k];
        }
        /// Set up 0-45 components of targeted buffer : cbuffer.data(

        auto g_y_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 45 * acomps  + 0);

        auto g_y_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 1);

        auto g_y_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 2);

        auto g_y_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 3);

        auto g_y_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 4);

        auto g_y_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 5);

        auto g_y_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 6);

        auto g_y_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 7);

        auto g_y_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 8);

        auto g_y_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 9);

        auto g_y_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 10);

        auto g_y_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 11);

        auto g_y_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 12);

        auto g_y_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 13);

        auto g_y_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 14);

        auto g_y_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 15);

        auto g_y_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 16);

        auto g_y_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 17);

        auto g_y_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 18);

        auto g_y_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 19);

        auto g_y_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 20);

        auto g_y_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 21);

        auto g_y_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 22);

        auto g_y_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 23);

        auto g_y_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 24);

        auto g_y_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 25);

        auto g_y_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 26);

        auto g_y_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 27);

        auto g_y_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 28);

        auto g_y_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 29);

        auto g_y_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 30);

        auto g_y_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 31);

        auto g_y_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 32);

        auto g_y_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 33);

        auto g_y_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 34);

        auto g_y_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 35);

        auto g_y_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 45 * acomps  + 36);

        auto g_y_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 37);

        auto g_y_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 38);

        auto g_y_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 39);

        auto g_y_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 40);

        auto g_y_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 41);

        auto g_y_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 42);

        auto g_y_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 43);

        auto g_y_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 45 * acomps  + 44);

        #pragma omp simd aligned(cd_y, g_0_xxxxxxxx, g_0_xxxxxxxxy, g_0_xxxxxxxy, g_0_xxxxxxxyy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxyy, g_0_xxxxxxyyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxyyy, g_0_xxxxxyyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxyyyy, g_0_xxxxyyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyyy, g_0_yyyyyyyyz, g_0_yyyyyyyz, g_0_yyyyyyyzz, g_0_yyyyyyzz, g_0_yyyyyyzzz, g_0_yyyyyzzz, g_0_yyyyyzzzz, g_0_yyyyzzzz, g_0_yyyyzzzzz, g_0_yyyzzzzz, g_0_yyyzzzzzz, g_0_yyzzzzzz, g_0_yyzzzzzzz, g_0_yzzzzzzz, g_0_yzzzzzzzz, g_0_zzzzzzzz, g_y_0_0_xxxxxxxx, g_y_0_0_xxxxxxxy, g_y_0_0_xxxxxxxz, g_y_0_0_xxxxxxyy, g_y_0_0_xxxxxxyz, g_y_0_0_xxxxxxzz, g_y_0_0_xxxxxyyy, g_y_0_0_xxxxxyyz, g_y_0_0_xxxxxyzz, g_y_0_0_xxxxxzzz, g_y_0_0_xxxxyyyy, g_y_0_0_xxxxyyyz, g_y_0_0_xxxxyyzz, g_y_0_0_xxxxyzzz, g_y_0_0_xxxxzzzz, g_y_0_0_xxxyyyyy, g_y_0_0_xxxyyyyz, g_y_0_0_xxxyyyzz, g_y_0_0_xxxyyzzz, g_y_0_0_xxxyzzzz, g_y_0_0_xxxzzzzz, g_y_0_0_xxyyyyyy, g_y_0_0_xxyyyyyz, g_y_0_0_xxyyyyzz, g_y_0_0_xxyyyzzz, g_y_0_0_xxyyzzzz, g_y_0_0_xxyzzzzz, g_y_0_0_xxzzzzzz, g_y_0_0_xyyyyyyy, g_y_0_0_xyyyyyyz, g_y_0_0_xyyyyyzz, g_y_0_0_xyyyyzzz, g_y_0_0_xyyyzzzz, g_y_0_0_xyyzzzzz, g_y_0_0_xyzzzzzz, g_y_0_0_xzzzzzzz, g_y_0_0_yyyyyyyy, g_y_0_0_yyyyyyyz, g_y_0_0_yyyyyyzz, g_y_0_0_yyyyyzzz, g_y_0_0_yyyyzzzz, g_y_0_0_yyyzzzzz, g_y_0_0_yyzzzzzz, g_y_0_0_yzzzzzzz, g_y_0_0_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_0_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_y[k] + g_0_xxxxxxxxy[k];

            g_y_0_0_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_y[k] + g_0_xxxxxxxyy[k];

            g_y_0_0_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_y[k] + g_0_xxxxxxxyz[k];

            g_y_0_0_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_y[k] + g_0_xxxxxxyyy[k];

            g_y_0_0_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_y[k] + g_0_xxxxxxyyz[k];

            g_y_0_0_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_y[k] + g_0_xxxxxxyzz[k];

            g_y_0_0_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_y[k] + g_0_xxxxxyyyy[k];

            g_y_0_0_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_y[k] + g_0_xxxxxyyyz[k];

            g_y_0_0_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_y[k] + g_0_xxxxxyyzz[k];

            g_y_0_0_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_y[k] + g_0_xxxxxyzzz[k];

            g_y_0_0_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_y[k] + g_0_xxxxyyyyy[k];

            g_y_0_0_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_y[k] + g_0_xxxxyyyyz[k];

            g_y_0_0_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_y[k] + g_0_xxxxyyyzz[k];

            g_y_0_0_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_y[k] + g_0_xxxxyyzzz[k];

            g_y_0_0_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_y[k] + g_0_xxxxyzzzz[k];

            g_y_0_0_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_y[k] + g_0_xxxyyyyyy[k];

            g_y_0_0_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_y[k] + g_0_xxxyyyyyz[k];

            g_y_0_0_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_y[k] + g_0_xxxyyyyzz[k];

            g_y_0_0_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_y[k] + g_0_xxxyyyzzz[k];

            g_y_0_0_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_y[k] + g_0_xxxyyzzzz[k];

            g_y_0_0_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_y[k] + g_0_xxxyzzzzz[k];

            g_y_0_0_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_y[k] + g_0_xxyyyyyyy[k];

            g_y_0_0_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_y[k] + g_0_xxyyyyyyz[k];

            g_y_0_0_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_y[k] + g_0_xxyyyyyzz[k];

            g_y_0_0_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_y[k] + g_0_xxyyyyzzz[k];

            g_y_0_0_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_y[k] + g_0_xxyyyzzzz[k];

            g_y_0_0_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_y[k] + g_0_xxyyzzzzz[k];

            g_y_0_0_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_y[k] + g_0_xxyzzzzzz[k];

            g_y_0_0_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_y[k] + g_0_xyyyyyyyy[k];

            g_y_0_0_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_y[k] + g_0_xyyyyyyyz[k];

            g_y_0_0_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_y[k] + g_0_xyyyyyyzz[k];

            g_y_0_0_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_y[k] + g_0_xyyyyyzzz[k];

            g_y_0_0_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_y[k] + g_0_xyyyyzzzz[k];

            g_y_0_0_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_y[k] + g_0_xyyyzzzzz[k];

            g_y_0_0_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_y[k] + g_0_xyyzzzzzz[k];

            g_y_0_0_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_y[k] + g_0_xyzzzzzzz[k];

            g_y_0_0_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_y[k] + g_0_yyyyyyyyy[k];

            g_y_0_0_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_y[k] + g_0_yyyyyyyyz[k];

            g_y_0_0_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_y[k] + g_0_yyyyyyyzz[k];

            g_y_0_0_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_y[k] + g_0_yyyyyyzzz[k];

            g_y_0_0_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_y[k] + g_0_yyyyyzzzz[k];

            g_y_0_0_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_y[k] + g_0_yyyyzzzzz[k];

            g_y_0_0_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_y[k] + g_0_yyyzzzzzz[k];

            g_y_0_0_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_y[k] + g_0_yyzzzzzzz[k];

            g_y_0_0_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_y[k] + g_0_yzzzzzzzz[k];
        }
        /// Set up 0-45 components of targeted buffer : cbuffer.data(

        auto g_z_0_0_xxxxxxxx = cbuffer.data(sl_geom_10_off + 90 * acomps  + 0);

        auto g_z_0_0_xxxxxxxy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 1);

        auto g_z_0_0_xxxxxxxz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 2);

        auto g_z_0_0_xxxxxxyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 3);

        auto g_z_0_0_xxxxxxyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 4);

        auto g_z_0_0_xxxxxxzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 5);

        auto g_z_0_0_xxxxxyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 6);

        auto g_z_0_0_xxxxxyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 7);

        auto g_z_0_0_xxxxxyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 8);

        auto g_z_0_0_xxxxxzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 9);

        auto g_z_0_0_xxxxyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 10);

        auto g_z_0_0_xxxxyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 11);

        auto g_z_0_0_xxxxyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 12);

        auto g_z_0_0_xxxxyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 13);

        auto g_z_0_0_xxxxzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 14);

        auto g_z_0_0_xxxyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 15);

        auto g_z_0_0_xxxyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 16);

        auto g_z_0_0_xxxyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 17);

        auto g_z_0_0_xxxyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 18);

        auto g_z_0_0_xxxyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 19);

        auto g_z_0_0_xxxzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 20);

        auto g_z_0_0_xxyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 21);

        auto g_z_0_0_xxyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 22);

        auto g_z_0_0_xxyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 23);

        auto g_z_0_0_xxyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 24);

        auto g_z_0_0_xxyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 25);

        auto g_z_0_0_xxyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 26);

        auto g_z_0_0_xxzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 27);

        auto g_z_0_0_xyyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 28);

        auto g_z_0_0_xyyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 29);

        auto g_z_0_0_xyyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 30);

        auto g_z_0_0_xyyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 31);

        auto g_z_0_0_xyyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 32);

        auto g_z_0_0_xyyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 33);

        auto g_z_0_0_xyzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 34);

        auto g_z_0_0_xzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 35);

        auto g_z_0_0_yyyyyyyy = cbuffer.data(sl_geom_10_off + 90 * acomps  + 36);

        auto g_z_0_0_yyyyyyyz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 37);

        auto g_z_0_0_yyyyyyzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 38);

        auto g_z_0_0_yyyyyzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 39);

        auto g_z_0_0_yyyyzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 40);

        auto g_z_0_0_yyyzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 41);

        auto g_z_0_0_yyzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 42);

        auto g_z_0_0_yzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 43);

        auto g_z_0_0_zzzzzzzz = cbuffer.data(sl_geom_10_off + 90 * acomps  + 44);

        #pragma omp simd aligned(cd_z, g_0_xxxxxxxx, g_0_xxxxxxxxz, g_0_xxxxxxxy, g_0_xxxxxxxyz, g_0_xxxxxxxz, g_0_xxxxxxxzz, g_0_xxxxxxyy, g_0_xxxxxxyyz, g_0_xxxxxxyz, g_0_xxxxxxyzz, g_0_xxxxxxzz, g_0_xxxxxxzzz, g_0_xxxxxyyy, g_0_xxxxxyyyz, g_0_xxxxxyyz, g_0_xxxxxyyzz, g_0_xxxxxyzz, g_0_xxxxxyzzz, g_0_xxxxxzzz, g_0_xxxxxzzzz, g_0_xxxxyyyy, g_0_xxxxyyyyz, g_0_xxxxyyyz, g_0_xxxxyyyzz, g_0_xxxxyyzz, g_0_xxxxyyzzz, g_0_xxxxyzzz, g_0_xxxxyzzzz, g_0_xxxxzzzz, g_0_xxxxzzzzz, g_0_xxxyyyyy, g_0_xxxyyyyyz, g_0_xxxyyyyz, g_0_xxxyyyyzz, g_0_xxxyyyzz, g_0_xxxyyyzzz, g_0_xxxyyzzz, g_0_xxxyyzzzz, g_0_xxxyzzzz, g_0_xxxyzzzzz, g_0_xxxzzzzz, g_0_xxxzzzzzz, g_0_xxyyyyyy, g_0_xxyyyyyyz, g_0_xxyyyyyz, g_0_xxyyyyyzz, g_0_xxyyyyzz, g_0_xxyyyyzzz, g_0_xxyyyzzz, g_0_xxyyyzzzz, g_0_xxyyzzzz, g_0_xxyyzzzzz, g_0_xxyzzzzz, g_0_xxyzzzzzz, g_0_xxzzzzzz, g_0_xxzzzzzzz, g_0_xyyyyyyy, g_0_xyyyyyyyz, g_0_xyyyyyyz, g_0_xyyyyyyzz, g_0_xyyyyyzz, g_0_xyyyyyzzz, g_0_xyyyyzzz, g_0_xyyyyzzzz, g_0_xyyyzzzz, g_0_xyyyzzzzz, g_0_xyyzzzzz, g_0_xyyzzzzzz, g_0_xyzzzzzz, g_0_xyzzzzzzz, g_0_xzzzzzzz, g_0_xzzzzzzzz, g_0_yyyyyyyy, g_0_yyyyyyyyz, g_0_yyyyyyyz, g_0_yyyyyyyzz, g_0_yyyyyyzz, g_0_yyyyyyzzz, g_0_yyyyyzzz, g_0_yyyyyzzzz, g_0_yyyyzzzz, g_0_yyyyzzzzz, g_0_yyyzzzzz, g_0_yyyzzzzzz, g_0_yyzzzzzz, g_0_yyzzzzzzz, g_0_yzzzzzzz, g_0_yzzzzzzzz, g_0_zzzzzzzz, g_0_zzzzzzzzz, g_z_0_0_xxxxxxxx, g_z_0_0_xxxxxxxy, g_z_0_0_xxxxxxxz, g_z_0_0_xxxxxxyy, g_z_0_0_xxxxxxyz, g_z_0_0_xxxxxxzz, g_z_0_0_xxxxxyyy, g_z_0_0_xxxxxyyz, g_z_0_0_xxxxxyzz, g_z_0_0_xxxxxzzz, g_z_0_0_xxxxyyyy, g_z_0_0_xxxxyyyz, g_z_0_0_xxxxyyzz, g_z_0_0_xxxxyzzz, g_z_0_0_xxxxzzzz, g_z_0_0_xxxyyyyy, g_z_0_0_xxxyyyyz, g_z_0_0_xxxyyyzz, g_z_0_0_xxxyyzzz, g_z_0_0_xxxyzzzz, g_z_0_0_xxxzzzzz, g_z_0_0_xxyyyyyy, g_z_0_0_xxyyyyyz, g_z_0_0_xxyyyyzz, g_z_0_0_xxyyyzzz, g_z_0_0_xxyyzzzz, g_z_0_0_xxyzzzzz, g_z_0_0_xxzzzzzz, g_z_0_0_xyyyyyyy, g_z_0_0_xyyyyyyz, g_z_0_0_xyyyyyzz, g_z_0_0_xyyyyzzz, g_z_0_0_xyyyzzzz, g_z_0_0_xyyzzzzz, g_z_0_0_xyzzzzzz, g_z_0_0_xzzzzzzz, g_z_0_0_yyyyyyyy, g_z_0_0_yyyyyyyz, g_z_0_0_yyyyyyzz, g_z_0_0_yyyyyzzz, g_z_0_0_yyyyzzzz, g_z_0_0_yyyzzzzz, g_z_0_0_yyzzzzzz, g_z_0_0_yzzzzzzz, g_z_0_0_zzzzzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_0_xxxxxxxx[k] = -g_0_xxxxxxxx[k] * cd_z[k] + g_0_xxxxxxxxz[k];

            g_z_0_0_xxxxxxxy[k] = -g_0_xxxxxxxy[k] * cd_z[k] + g_0_xxxxxxxyz[k];

            g_z_0_0_xxxxxxxz[k] = -g_0_xxxxxxxz[k] * cd_z[k] + g_0_xxxxxxxzz[k];

            g_z_0_0_xxxxxxyy[k] = -g_0_xxxxxxyy[k] * cd_z[k] + g_0_xxxxxxyyz[k];

            g_z_0_0_xxxxxxyz[k] = -g_0_xxxxxxyz[k] * cd_z[k] + g_0_xxxxxxyzz[k];

            g_z_0_0_xxxxxxzz[k] = -g_0_xxxxxxzz[k] * cd_z[k] + g_0_xxxxxxzzz[k];

            g_z_0_0_xxxxxyyy[k] = -g_0_xxxxxyyy[k] * cd_z[k] + g_0_xxxxxyyyz[k];

            g_z_0_0_xxxxxyyz[k] = -g_0_xxxxxyyz[k] * cd_z[k] + g_0_xxxxxyyzz[k];

            g_z_0_0_xxxxxyzz[k] = -g_0_xxxxxyzz[k] * cd_z[k] + g_0_xxxxxyzzz[k];

            g_z_0_0_xxxxxzzz[k] = -g_0_xxxxxzzz[k] * cd_z[k] + g_0_xxxxxzzzz[k];

            g_z_0_0_xxxxyyyy[k] = -g_0_xxxxyyyy[k] * cd_z[k] + g_0_xxxxyyyyz[k];

            g_z_0_0_xxxxyyyz[k] = -g_0_xxxxyyyz[k] * cd_z[k] + g_0_xxxxyyyzz[k];

            g_z_0_0_xxxxyyzz[k] = -g_0_xxxxyyzz[k] * cd_z[k] + g_0_xxxxyyzzz[k];

            g_z_0_0_xxxxyzzz[k] = -g_0_xxxxyzzz[k] * cd_z[k] + g_0_xxxxyzzzz[k];

            g_z_0_0_xxxxzzzz[k] = -g_0_xxxxzzzz[k] * cd_z[k] + g_0_xxxxzzzzz[k];

            g_z_0_0_xxxyyyyy[k] = -g_0_xxxyyyyy[k] * cd_z[k] + g_0_xxxyyyyyz[k];

            g_z_0_0_xxxyyyyz[k] = -g_0_xxxyyyyz[k] * cd_z[k] + g_0_xxxyyyyzz[k];

            g_z_0_0_xxxyyyzz[k] = -g_0_xxxyyyzz[k] * cd_z[k] + g_0_xxxyyyzzz[k];

            g_z_0_0_xxxyyzzz[k] = -g_0_xxxyyzzz[k] * cd_z[k] + g_0_xxxyyzzzz[k];

            g_z_0_0_xxxyzzzz[k] = -g_0_xxxyzzzz[k] * cd_z[k] + g_0_xxxyzzzzz[k];

            g_z_0_0_xxxzzzzz[k] = -g_0_xxxzzzzz[k] * cd_z[k] + g_0_xxxzzzzzz[k];

            g_z_0_0_xxyyyyyy[k] = -g_0_xxyyyyyy[k] * cd_z[k] + g_0_xxyyyyyyz[k];

            g_z_0_0_xxyyyyyz[k] = -g_0_xxyyyyyz[k] * cd_z[k] + g_0_xxyyyyyzz[k];

            g_z_0_0_xxyyyyzz[k] = -g_0_xxyyyyzz[k] * cd_z[k] + g_0_xxyyyyzzz[k];

            g_z_0_0_xxyyyzzz[k] = -g_0_xxyyyzzz[k] * cd_z[k] + g_0_xxyyyzzzz[k];

            g_z_0_0_xxyyzzzz[k] = -g_0_xxyyzzzz[k] * cd_z[k] + g_0_xxyyzzzzz[k];

            g_z_0_0_xxyzzzzz[k] = -g_0_xxyzzzzz[k] * cd_z[k] + g_0_xxyzzzzzz[k];

            g_z_0_0_xxzzzzzz[k] = -g_0_xxzzzzzz[k] * cd_z[k] + g_0_xxzzzzzzz[k];

            g_z_0_0_xyyyyyyy[k] = -g_0_xyyyyyyy[k] * cd_z[k] + g_0_xyyyyyyyz[k];

            g_z_0_0_xyyyyyyz[k] = -g_0_xyyyyyyz[k] * cd_z[k] + g_0_xyyyyyyzz[k];

            g_z_0_0_xyyyyyzz[k] = -g_0_xyyyyyzz[k] * cd_z[k] + g_0_xyyyyyzzz[k];

            g_z_0_0_xyyyyzzz[k] = -g_0_xyyyyzzz[k] * cd_z[k] + g_0_xyyyyzzzz[k];

            g_z_0_0_xyyyzzzz[k] = -g_0_xyyyzzzz[k] * cd_z[k] + g_0_xyyyzzzzz[k];

            g_z_0_0_xyyzzzzz[k] = -g_0_xyyzzzzz[k] * cd_z[k] + g_0_xyyzzzzzz[k];

            g_z_0_0_xyzzzzzz[k] = -g_0_xyzzzzzz[k] * cd_z[k] + g_0_xyzzzzzzz[k];

            g_z_0_0_xzzzzzzz[k] = -g_0_xzzzzzzz[k] * cd_z[k] + g_0_xzzzzzzzz[k];

            g_z_0_0_yyyyyyyy[k] = -g_0_yyyyyyyy[k] * cd_z[k] + g_0_yyyyyyyyz[k];

            g_z_0_0_yyyyyyyz[k] = -g_0_yyyyyyyz[k] * cd_z[k] + g_0_yyyyyyyzz[k];

            g_z_0_0_yyyyyyzz[k] = -g_0_yyyyyyzz[k] * cd_z[k] + g_0_yyyyyyzzz[k];

            g_z_0_0_yyyyyzzz[k] = -g_0_yyyyyzzz[k] * cd_z[k] + g_0_yyyyyzzzz[k];

            g_z_0_0_yyyyzzzz[k] = -g_0_yyyyzzzz[k] * cd_z[k] + g_0_yyyyzzzzz[k];

            g_z_0_0_yyyzzzzz[k] = -g_0_yyyzzzzz[k] * cd_z[k] + g_0_yyyzzzzzz[k];

            g_z_0_0_yyzzzzzz[k] = -g_0_yyzzzzzz[k] * cd_z[k] + g_0_yyzzzzzzz[k];

            g_z_0_0_yzzzzzzz[k] = -g_0_yzzzzzzz[k] * cd_z[k] + g_0_yzzzzzzzz[k];

            g_z_0_0_zzzzzzzz[k] = -g_0_zzzzzzzz[k] * cd_z[k] + g_0_zzzzzzzzz[k];
        }
    }
}

} // t3ceri namespace


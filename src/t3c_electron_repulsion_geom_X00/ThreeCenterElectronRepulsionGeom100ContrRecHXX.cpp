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

#include "ThreeCenterElectronRepulsionGeom100ContrRecHXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_hxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_hxx,
                                      const size_t idx_gxx,
                                      const size_t idx_ixx,
                                      const int c_angmom,
                                      const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : GSS

            const auto g_off = idx_gxx + i * dcomps + j;

            auto g_xxxx_0 = cbuffer.data(g_off + 0 * ccomps * dcomps);

            auto g_xxxy_0 = cbuffer.data(g_off + 1 * ccomps * dcomps);

            auto g_xxxz_0 = cbuffer.data(g_off + 2 * ccomps * dcomps);

            auto g_xxyy_0 = cbuffer.data(g_off + 3 * ccomps * dcomps);

            auto g_xxyz_0 = cbuffer.data(g_off + 4 * ccomps * dcomps);

            auto g_xxzz_0 = cbuffer.data(g_off + 5 * ccomps * dcomps);

            auto g_xyyy_0 = cbuffer.data(g_off + 6 * ccomps * dcomps);

            auto g_xyyz_0 = cbuffer.data(g_off + 7 * ccomps * dcomps);

            auto g_xyzz_0 = cbuffer.data(g_off + 8 * ccomps * dcomps);

            auto g_xzzz_0 = cbuffer.data(g_off + 9 * ccomps * dcomps);

            auto g_yyyy_0 = cbuffer.data(g_off + 10 * ccomps * dcomps);

            auto g_yyyz_0 = cbuffer.data(g_off + 11 * ccomps * dcomps);

            auto g_yyzz_0 = cbuffer.data(g_off + 12 * ccomps * dcomps);

            auto g_yzzz_0 = cbuffer.data(g_off + 13 * ccomps * dcomps);

            auto g_zzzz_0 = cbuffer.data(g_off + 14 * ccomps * dcomps);

            /// Set up components of auxilary buffer : ISS

            const auto i_off = idx_ixx + i * dcomps + j;

            auto g_xxxxxx_0 = cbuffer.data(i_off + 0 * ccomps * dcomps);

            auto g_xxxxxy_0 = cbuffer.data(i_off + 1 * ccomps * dcomps);

            auto g_xxxxxz_0 = cbuffer.data(i_off + 2 * ccomps * dcomps);

            auto g_xxxxyy_0 = cbuffer.data(i_off + 3 * ccomps * dcomps);

            auto g_xxxxyz_0 = cbuffer.data(i_off + 4 * ccomps * dcomps);

            auto g_xxxxzz_0 = cbuffer.data(i_off + 5 * ccomps * dcomps);

            auto g_xxxyyy_0 = cbuffer.data(i_off + 6 * ccomps * dcomps);

            auto g_xxxyyz_0 = cbuffer.data(i_off + 7 * ccomps * dcomps);

            auto g_xxxyzz_0 = cbuffer.data(i_off + 8 * ccomps * dcomps);

            auto g_xxxzzz_0 = cbuffer.data(i_off + 9 * ccomps * dcomps);

            auto g_xxyyyy_0 = cbuffer.data(i_off + 10 * ccomps * dcomps);

            auto g_xxyyyz_0 = cbuffer.data(i_off + 11 * ccomps * dcomps);

            auto g_xxyyzz_0 = cbuffer.data(i_off + 12 * ccomps * dcomps);

            auto g_xxyzzz_0 = cbuffer.data(i_off + 13 * ccomps * dcomps);

            auto g_xxzzzz_0 = cbuffer.data(i_off + 14 * ccomps * dcomps);

            auto g_xyyyyy_0 = cbuffer.data(i_off + 15 * ccomps * dcomps);

            auto g_xyyyyz_0 = cbuffer.data(i_off + 16 * ccomps * dcomps);

            auto g_xyyyzz_0 = cbuffer.data(i_off + 17 * ccomps * dcomps);

            auto g_xyyzzz_0 = cbuffer.data(i_off + 18 * ccomps * dcomps);

            auto g_xyzzzz_0 = cbuffer.data(i_off + 19 * ccomps * dcomps);

            auto g_xzzzzz_0 = cbuffer.data(i_off + 20 * ccomps * dcomps);

            auto g_yyyyyy_0 = cbuffer.data(i_off + 21 * ccomps * dcomps);

            auto g_yyyyyz_0 = cbuffer.data(i_off + 22 * ccomps * dcomps);

            auto g_yyyyzz_0 = cbuffer.data(i_off + 23 * ccomps * dcomps);

            auto g_yyyzzz_0 = cbuffer.data(i_off + 24 * ccomps * dcomps);

            auto g_yyzzzz_0 = cbuffer.data(i_off + 25 * ccomps * dcomps);

            auto g_yzzzzz_0 = cbuffer.data(i_off + 26 * ccomps * dcomps);

            auto g_zzzzzz_0 = cbuffer.data(i_off + 27 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hxx

            const auto h_geom_100_off = idx_geom_100_hxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxxxx_0 = cbuffer.data(h_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxy_0 = cbuffer.data(h_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxz_0 = cbuffer.data(h_geom_100_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxyy_0 = cbuffer.data(h_geom_100_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxyz_0 = cbuffer.data(h_geom_100_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxzz_0 = cbuffer.data(h_geom_100_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxyyy_0 = cbuffer.data(h_geom_100_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxyyz_0 = cbuffer.data(h_geom_100_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxyzz_0 = cbuffer.data(h_geom_100_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxzzz_0 = cbuffer.data(h_geom_100_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xyyyy_0 = cbuffer.data(h_geom_100_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xyyyz_0 = cbuffer.data(h_geom_100_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xyyzz_0 = cbuffer.data(h_geom_100_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xyzzz_0 = cbuffer.data(h_geom_100_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xzzzz_0 = cbuffer.data(h_geom_100_off + 14 * ccomps * dcomps);

            auto g_x_0_0_yyyyy_0 = cbuffer.data(h_geom_100_off + 15 * ccomps * dcomps);

            auto g_x_0_0_yyyyz_0 = cbuffer.data(h_geom_100_off + 16 * ccomps * dcomps);

            auto g_x_0_0_yyyzz_0 = cbuffer.data(h_geom_100_off + 17 * ccomps * dcomps);

            auto g_x_0_0_yyzzz_0 = cbuffer.data(h_geom_100_off + 18 * ccomps * dcomps);

            auto g_x_0_0_yzzzz_0 = cbuffer.data(h_geom_100_off + 19 * ccomps * dcomps);

            auto g_x_0_0_zzzzz_0 = cbuffer.data(h_geom_100_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxx_0, g_x_0_0_xxxxy_0, g_x_0_0_xxxxz_0, g_x_0_0_xxxyy_0, g_x_0_0_xxxyz_0, g_x_0_0_xxxzz_0, g_x_0_0_xxyyy_0, g_x_0_0_xxyyz_0, g_x_0_0_xxyzz_0, g_x_0_0_xxzzz_0, g_x_0_0_xyyyy_0, g_x_0_0_xyyyz_0, g_x_0_0_xyyzz_0, g_x_0_0_xyzzz_0, g_x_0_0_xzzzz_0, g_x_0_0_yyyyy_0, g_x_0_0_yyyyz_0, g_x_0_0_yyyzz_0, g_x_0_0_yyzzz_0, g_x_0_0_yzzzz_0, g_x_0_0_zzzzz_0, g_xxxx_0, g_xxxxxx_0, g_xxxxxy_0, g_xxxxxz_0, g_xxxxyy_0, g_xxxxyz_0, g_xxxxzz_0, g_xxxy_0, g_xxxyyy_0, g_xxxyyz_0, g_xxxyzz_0, g_xxxz_0, g_xxxzzz_0, g_xxyy_0, g_xxyyyy_0, g_xxyyyz_0, g_xxyyzz_0, g_xxyz_0, g_xxyzzz_0, g_xxzz_0, g_xxzzzz_0, g_xyyy_0, g_xyyyyy_0, g_xyyyyz_0, g_xyyyzz_0, g_xyyz_0, g_xyyzzz_0, g_xyzz_0, g_xyzzzz_0, g_xzzz_0, g_xzzzzz_0, g_yyyy_0, g_yyyz_0, g_yyzz_0, g_yzzz_0, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxxx_0[k] = -5.0 * g_xxxx_0[k] + g_xxxxxx_0[k];

                g_x_0_0_xxxxy_0[k] = -4.0 * g_xxxy_0[k] + g_xxxxxy_0[k];

                g_x_0_0_xxxxz_0[k] = -4.0 * g_xxxz_0[k] + g_xxxxxz_0[k];

                g_x_0_0_xxxyy_0[k] = -3.0 * g_xxyy_0[k] + g_xxxxyy_0[k];

                g_x_0_0_xxxyz_0[k] = -3.0 * g_xxyz_0[k] + g_xxxxyz_0[k];

                g_x_0_0_xxxzz_0[k] = -3.0 * g_xxzz_0[k] + g_xxxxzz_0[k];

                g_x_0_0_xxyyy_0[k] = -2.0 * g_xyyy_0[k] + g_xxxyyy_0[k];

                g_x_0_0_xxyyz_0[k] = -2.0 * g_xyyz_0[k] + g_xxxyyz_0[k];

                g_x_0_0_xxyzz_0[k] = -2.0 * g_xyzz_0[k] + g_xxxyzz_0[k];

                g_x_0_0_xxzzz_0[k] = -2.0 * g_xzzz_0[k] + g_xxxzzz_0[k];

                g_x_0_0_xyyyy_0[k] = -g_yyyy_0[k] + g_xxyyyy_0[k];

                g_x_0_0_xyyyz_0[k] = -g_yyyz_0[k] + g_xxyyyz_0[k];

                g_x_0_0_xyyzz_0[k] = -g_yyzz_0[k] + g_xxyyzz_0[k];

                g_x_0_0_xyzzz_0[k] = -g_yzzz_0[k] + g_xxyzzz_0[k];

                g_x_0_0_xzzzz_0[k] = -g_zzzz_0[k] + g_xxzzzz_0[k];

                g_x_0_0_yyyyy_0[k] = g_xyyyyy_0[k];

                g_x_0_0_yyyyz_0[k] = g_xyyyyz_0[k];

                g_x_0_0_yyyzz_0[k] = g_xyyyzz_0[k];

                g_x_0_0_yyzzz_0[k] = g_xyyzzz_0[k];

                g_x_0_0_yzzzz_0[k] = g_xyzzzz_0[k];

                g_x_0_0_zzzzz_0[k] = g_xzzzzz_0[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxxxx_0 = cbuffer.data(h_geom_100_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xxxxy_0 = cbuffer.data(h_geom_100_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xxxxz_0 = cbuffer.data(h_geom_100_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xxxyy_0 = cbuffer.data(h_geom_100_off + 24 * ccomps * dcomps);

            auto g_y_0_0_xxxyz_0 = cbuffer.data(h_geom_100_off + 25 * ccomps * dcomps);

            auto g_y_0_0_xxxzz_0 = cbuffer.data(h_geom_100_off + 26 * ccomps * dcomps);

            auto g_y_0_0_xxyyy_0 = cbuffer.data(h_geom_100_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxyyz_0 = cbuffer.data(h_geom_100_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxyzz_0 = cbuffer.data(h_geom_100_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxzzz_0 = cbuffer.data(h_geom_100_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xyyyy_0 = cbuffer.data(h_geom_100_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xyyyz_0 = cbuffer.data(h_geom_100_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xyyzz_0 = cbuffer.data(h_geom_100_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xyzzz_0 = cbuffer.data(h_geom_100_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xzzzz_0 = cbuffer.data(h_geom_100_off + 35 * ccomps * dcomps);

            auto g_y_0_0_yyyyy_0 = cbuffer.data(h_geom_100_off + 36 * ccomps * dcomps);

            auto g_y_0_0_yyyyz_0 = cbuffer.data(h_geom_100_off + 37 * ccomps * dcomps);

            auto g_y_0_0_yyyzz_0 = cbuffer.data(h_geom_100_off + 38 * ccomps * dcomps);

            auto g_y_0_0_yyzzz_0 = cbuffer.data(h_geom_100_off + 39 * ccomps * dcomps);

            auto g_y_0_0_yzzzz_0 = cbuffer.data(h_geom_100_off + 40 * ccomps * dcomps);

            auto g_y_0_0_zzzzz_0 = cbuffer.data(h_geom_100_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxx_0, g_xxxxxy_0, g_xxxxyy_0, g_xxxxyz_0, g_xxxy_0, g_xxxyyy_0, g_xxxyyz_0, g_xxxyzz_0, g_xxxz_0, g_xxyy_0, g_xxyyyy_0, g_xxyyyz_0, g_xxyyzz_0, g_xxyz_0, g_xxyzzz_0, g_xxzz_0, g_xyyy_0, g_xyyyyy_0, g_xyyyyz_0, g_xyyyzz_0, g_xyyz_0, g_xyyzzz_0, g_xyzz_0, g_xyzzzz_0, g_xzzz_0, g_y_0_0_xxxxx_0, g_y_0_0_xxxxy_0, g_y_0_0_xxxxz_0, g_y_0_0_xxxyy_0, g_y_0_0_xxxyz_0, g_y_0_0_xxxzz_0, g_y_0_0_xxyyy_0, g_y_0_0_xxyyz_0, g_y_0_0_xxyzz_0, g_y_0_0_xxzzz_0, g_y_0_0_xyyyy_0, g_y_0_0_xyyyz_0, g_y_0_0_xyyzz_0, g_y_0_0_xyzzz_0, g_y_0_0_xzzzz_0, g_y_0_0_yyyyy_0, g_y_0_0_yyyyz_0, g_y_0_0_yyyzz_0, g_y_0_0_yyzzz_0, g_y_0_0_yzzzz_0, g_y_0_0_zzzzz_0, g_yyyy_0, g_yyyyyy_0, g_yyyyyz_0, g_yyyyzz_0, g_yyyz_0, g_yyyzzz_0, g_yyzz_0, g_yyzzzz_0, g_yzzz_0, g_yzzzzz_0, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxxx_0[k] = g_xxxxxy_0[k];

                g_y_0_0_xxxxy_0[k] = -g_xxxx_0[k] + g_xxxxyy_0[k];

                g_y_0_0_xxxxz_0[k] = g_xxxxyz_0[k];

                g_y_0_0_xxxyy_0[k] = -2.0 * g_xxxy_0[k] + g_xxxyyy_0[k];

                g_y_0_0_xxxyz_0[k] = -g_xxxz_0[k] + g_xxxyyz_0[k];

                g_y_0_0_xxxzz_0[k] = g_xxxyzz_0[k];

                g_y_0_0_xxyyy_0[k] = -3.0 * g_xxyy_0[k] + g_xxyyyy_0[k];

                g_y_0_0_xxyyz_0[k] = -2.0 * g_xxyz_0[k] + g_xxyyyz_0[k];

                g_y_0_0_xxyzz_0[k] = -g_xxzz_0[k] + g_xxyyzz_0[k];

                g_y_0_0_xxzzz_0[k] = g_xxyzzz_0[k];

                g_y_0_0_xyyyy_0[k] = -4.0 * g_xyyy_0[k] + g_xyyyyy_0[k];

                g_y_0_0_xyyyz_0[k] = -3.0 * g_xyyz_0[k] + g_xyyyyz_0[k];

                g_y_0_0_xyyzz_0[k] = -2.0 * g_xyzz_0[k] + g_xyyyzz_0[k];

                g_y_0_0_xyzzz_0[k] = -g_xzzz_0[k] + g_xyyzzz_0[k];

                g_y_0_0_xzzzz_0[k] = g_xyzzzz_0[k];

                g_y_0_0_yyyyy_0[k] = -5.0 * g_yyyy_0[k] + g_yyyyyy_0[k];

                g_y_0_0_yyyyz_0[k] = -4.0 * g_yyyz_0[k] + g_yyyyyz_0[k];

                g_y_0_0_yyyzz_0[k] = -3.0 * g_yyzz_0[k] + g_yyyyzz_0[k];

                g_y_0_0_yyzzz_0[k] = -2.0 * g_yzzz_0[k] + g_yyyzzz_0[k];

                g_y_0_0_yzzzz_0[k] = -g_zzzz_0[k] + g_yyzzzz_0[k];

                g_y_0_0_zzzzz_0[k] = g_yzzzzz_0[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxxxx_0 = cbuffer.data(h_geom_100_off + 42 * ccomps * dcomps);

            auto g_z_0_0_xxxxy_0 = cbuffer.data(h_geom_100_off + 43 * ccomps * dcomps);

            auto g_z_0_0_xxxxz_0 = cbuffer.data(h_geom_100_off + 44 * ccomps * dcomps);

            auto g_z_0_0_xxxyy_0 = cbuffer.data(h_geom_100_off + 45 * ccomps * dcomps);

            auto g_z_0_0_xxxyz_0 = cbuffer.data(h_geom_100_off + 46 * ccomps * dcomps);

            auto g_z_0_0_xxxzz_0 = cbuffer.data(h_geom_100_off + 47 * ccomps * dcomps);

            auto g_z_0_0_xxyyy_0 = cbuffer.data(h_geom_100_off + 48 * ccomps * dcomps);

            auto g_z_0_0_xxyyz_0 = cbuffer.data(h_geom_100_off + 49 * ccomps * dcomps);

            auto g_z_0_0_xxyzz_0 = cbuffer.data(h_geom_100_off + 50 * ccomps * dcomps);

            auto g_z_0_0_xxzzz_0 = cbuffer.data(h_geom_100_off + 51 * ccomps * dcomps);

            auto g_z_0_0_xyyyy_0 = cbuffer.data(h_geom_100_off + 52 * ccomps * dcomps);

            auto g_z_0_0_xyyyz_0 = cbuffer.data(h_geom_100_off + 53 * ccomps * dcomps);

            auto g_z_0_0_xyyzz_0 = cbuffer.data(h_geom_100_off + 54 * ccomps * dcomps);

            auto g_z_0_0_xyzzz_0 = cbuffer.data(h_geom_100_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xzzzz_0 = cbuffer.data(h_geom_100_off + 56 * ccomps * dcomps);

            auto g_z_0_0_yyyyy_0 = cbuffer.data(h_geom_100_off + 57 * ccomps * dcomps);

            auto g_z_0_0_yyyyz_0 = cbuffer.data(h_geom_100_off + 58 * ccomps * dcomps);

            auto g_z_0_0_yyyzz_0 = cbuffer.data(h_geom_100_off + 59 * ccomps * dcomps);

            auto g_z_0_0_yyzzz_0 = cbuffer.data(h_geom_100_off + 60 * ccomps * dcomps);

            auto g_z_0_0_yzzzz_0 = cbuffer.data(h_geom_100_off + 61 * ccomps * dcomps);

            auto g_z_0_0_zzzzz_0 = cbuffer.data(h_geom_100_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxx_0, g_xxxxxz_0, g_xxxxyz_0, g_xxxxzz_0, g_xxxy_0, g_xxxyyz_0, g_xxxyzz_0, g_xxxz_0, g_xxxzzz_0, g_xxyy_0, g_xxyyyz_0, g_xxyyzz_0, g_xxyz_0, g_xxyzzz_0, g_xxzz_0, g_xxzzzz_0, g_xyyy_0, g_xyyyyz_0, g_xyyyzz_0, g_xyyz_0, g_xyyzzz_0, g_xyzz_0, g_xyzzzz_0, g_xzzz_0, g_xzzzzz_0, g_yyyy_0, g_yyyyyz_0, g_yyyyzz_0, g_yyyz_0, g_yyyzzz_0, g_yyzz_0, g_yyzzzz_0, g_yzzz_0, g_yzzzzz_0, g_z_0_0_xxxxx_0, g_z_0_0_xxxxy_0, g_z_0_0_xxxxz_0, g_z_0_0_xxxyy_0, g_z_0_0_xxxyz_0, g_z_0_0_xxxzz_0, g_z_0_0_xxyyy_0, g_z_0_0_xxyyz_0, g_z_0_0_xxyzz_0, g_z_0_0_xxzzz_0, g_z_0_0_xyyyy_0, g_z_0_0_xyyyz_0, g_z_0_0_xyyzz_0, g_z_0_0_xyzzz_0, g_z_0_0_xzzzz_0, g_z_0_0_yyyyy_0, g_z_0_0_yyyyz_0, g_z_0_0_yyyzz_0, g_z_0_0_yyzzz_0, g_z_0_0_yzzzz_0, g_z_0_0_zzzzz_0, g_zzzz_0, g_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxxx_0[k] = g_xxxxxz_0[k];

                g_z_0_0_xxxxy_0[k] = g_xxxxyz_0[k];

                g_z_0_0_xxxxz_0[k] = -g_xxxx_0[k] + g_xxxxzz_0[k];

                g_z_0_0_xxxyy_0[k] = g_xxxyyz_0[k];

                g_z_0_0_xxxyz_0[k] = -g_xxxy_0[k] + g_xxxyzz_0[k];

                g_z_0_0_xxxzz_0[k] = -2.0 * g_xxxz_0[k] + g_xxxzzz_0[k];

                g_z_0_0_xxyyy_0[k] = g_xxyyyz_0[k];

                g_z_0_0_xxyyz_0[k] = -g_xxyy_0[k] + g_xxyyzz_0[k];

                g_z_0_0_xxyzz_0[k] = -2.0 * g_xxyz_0[k] + g_xxyzzz_0[k];

                g_z_0_0_xxzzz_0[k] = -3.0 * g_xxzz_0[k] + g_xxzzzz_0[k];

                g_z_0_0_xyyyy_0[k] = g_xyyyyz_0[k];

                g_z_0_0_xyyyz_0[k] = -g_xyyy_0[k] + g_xyyyzz_0[k];

                g_z_0_0_xyyzz_0[k] = -2.0 * g_xyyz_0[k] + g_xyyzzz_0[k];

                g_z_0_0_xyzzz_0[k] = -3.0 * g_xyzz_0[k] + g_xyzzzz_0[k];

                g_z_0_0_xzzzz_0[k] = -4.0 * g_xzzz_0[k] + g_xzzzzz_0[k];

                g_z_0_0_yyyyy_0[k] = g_yyyyyz_0[k];

                g_z_0_0_yyyyz_0[k] = -g_yyyy_0[k] + g_yyyyzz_0[k];

                g_z_0_0_yyyzz_0[k] = -2.0 * g_yyyz_0[k] + g_yyyzzz_0[k];

                g_z_0_0_yyzzz_0[k] = -3.0 * g_yyzz_0[k] + g_yyzzzz_0[k];

                g_z_0_0_yzzzz_0[k] = -4.0 * g_yzzz_0[k] + g_yzzzzz_0[k];

                g_z_0_0_zzzzz_0[k] = -5.0 * g_zzzz_0[k] + g_zzzzzz_0[k];
            }
        }
    }
}

} // t3ceri namespace


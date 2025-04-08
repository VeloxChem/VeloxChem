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

#include "ElectronRepulsionGeom1000ContrRecIHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_ihxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ihxx,
                                            const size_t idx_hhxx,
                                            const size_t idx_geom_10_hhxx,
                                            const size_t idx_geom_10_hixx,
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
            /// Set up components of auxilary buffer : HHSS

            const auto hh_off = idx_hhxx + i * dcomps + j;

            auto g_xxxxx_xxxxx = cbuffer.data(hh_off + 0 * ccomps * dcomps);

            auto g_xxxxx_xxxxy = cbuffer.data(hh_off + 1 * ccomps * dcomps);

            auto g_xxxxx_xxxxz = cbuffer.data(hh_off + 2 * ccomps * dcomps);

            auto g_xxxxx_xxxyy = cbuffer.data(hh_off + 3 * ccomps * dcomps);

            auto g_xxxxx_xxxyz = cbuffer.data(hh_off + 4 * ccomps * dcomps);

            auto g_xxxxx_xxxzz = cbuffer.data(hh_off + 5 * ccomps * dcomps);

            auto g_xxxxx_xxyyy = cbuffer.data(hh_off + 6 * ccomps * dcomps);

            auto g_xxxxx_xxyyz = cbuffer.data(hh_off + 7 * ccomps * dcomps);

            auto g_xxxxx_xxyzz = cbuffer.data(hh_off + 8 * ccomps * dcomps);

            auto g_xxxxx_xxzzz = cbuffer.data(hh_off + 9 * ccomps * dcomps);

            auto g_xxxxx_xyyyy = cbuffer.data(hh_off + 10 * ccomps * dcomps);

            auto g_xxxxx_xyyyz = cbuffer.data(hh_off + 11 * ccomps * dcomps);

            auto g_xxxxx_xyyzz = cbuffer.data(hh_off + 12 * ccomps * dcomps);

            auto g_xxxxx_xyzzz = cbuffer.data(hh_off + 13 * ccomps * dcomps);

            auto g_xxxxx_xzzzz = cbuffer.data(hh_off + 14 * ccomps * dcomps);

            auto g_xxxxx_yyyyy = cbuffer.data(hh_off + 15 * ccomps * dcomps);

            auto g_xxxxx_yyyyz = cbuffer.data(hh_off + 16 * ccomps * dcomps);

            auto g_xxxxx_yyyzz = cbuffer.data(hh_off + 17 * ccomps * dcomps);

            auto g_xxxxx_yyzzz = cbuffer.data(hh_off + 18 * ccomps * dcomps);

            auto g_xxxxx_yzzzz = cbuffer.data(hh_off + 19 * ccomps * dcomps);

            auto g_xxxxx_zzzzz = cbuffer.data(hh_off + 20 * ccomps * dcomps);

            auto g_yyyyy_xxxxx = cbuffer.data(hh_off + 315 * ccomps * dcomps);

            auto g_yyyyy_xxxxy = cbuffer.data(hh_off + 316 * ccomps * dcomps);

            auto g_yyyyy_xxxxz = cbuffer.data(hh_off + 317 * ccomps * dcomps);

            auto g_yyyyy_xxxyy = cbuffer.data(hh_off + 318 * ccomps * dcomps);

            auto g_yyyyy_xxxyz = cbuffer.data(hh_off + 319 * ccomps * dcomps);

            auto g_yyyyy_xxxzz = cbuffer.data(hh_off + 320 * ccomps * dcomps);

            auto g_yyyyy_xxyyy = cbuffer.data(hh_off + 321 * ccomps * dcomps);

            auto g_yyyyy_xxyyz = cbuffer.data(hh_off + 322 * ccomps * dcomps);

            auto g_yyyyy_xxyzz = cbuffer.data(hh_off + 323 * ccomps * dcomps);

            auto g_yyyyy_xxzzz = cbuffer.data(hh_off + 324 * ccomps * dcomps);

            auto g_yyyyy_xyyyy = cbuffer.data(hh_off + 325 * ccomps * dcomps);

            auto g_yyyyy_xyyyz = cbuffer.data(hh_off + 326 * ccomps * dcomps);

            auto g_yyyyy_xyyzz = cbuffer.data(hh_off + 327 * ccomps * dcomps);

            auto g_yyyyy_xyzzz = cbuffer.data(hh_off + 328 * ccomps * dcomps);

            auto g_yyyyy_xzzzz = cbuffer.data(hh_off + 329 * ccomps * dcomps);

            auto g_yyyyy_yyyyy = cbuffer.data(hh_off + 330 * ccomps * dcomps);

            auto g_yyyyy_yyyyz = cbuffer.data(hh_off + 331 * ccomps * dcomps);

            auto g_yyyyy_yyyzz = cbuffer.data(hh_off + 332 * ccomps * dcomps);

            auto g_yyyyy_yyzzz = cbuffer.data(hh_off + 333 * ccomps * dcomps);

            auto g_yyyyy_yzzzz = cbuffer.data(hh_off + 334 * ccomps * dcomps);

            auto g_yyyyy_zzzzz = cbuffer.data(hh_off + 335 * ccomps * dcomps);

            auto g_zzzzz_xxxxx = cbuffer.data(hh_off + 420 * ccomps * dcomps);

            auto g_zzzzz_xxxxy = cbuffer.data(hh_off + 421 * ccomps * dcomps);

            auto g_zzzzz_xxxxz = cbuffer.data(hh_off + 422 * ccomps * dcomps);

            auto g_zzzzz_xxxyy = cbuffer.data(hh_off + 423 * ccomps * dcomps);

            auto g_zzzzz_xxxyz = cbuffer.data(hh_off + 424 * ccomps * dcomps);

            auto g_zzzzz_xxxzz = cbuffer.data(hh_off + 425 * ccomps * dcomps);

            auto g_zzzzz_xxyyy = cbuffer.data(hh_off + 426 * ccomps * dcomps);

            auto g_zzzzz_xxyyz = cbuffer.data(hh_off + 427 * ccomps * dcomps);

            auto g_zzzzz_xxyzz = cbuffer.data(hh_off + 428 * ccomps * dcomps);

            auto g_zzzzz_xxzzz = cbuffer.data(hh_off + 429 * ccomps * dcomps);

            auto g_zzzzz_xyyyy = cbuffer.data(hh_off + 430 * ccomps * dcomps);

            auto g_zzzzz_xyyyz = cbuffer.data(hh_off + 431 * ccomps * dcomps);

            auto g_zzzzz_xyyzz = cbuffer.data(hh_off + 432 * ccomps * dcomps);

            auto g_zzzzz_xyzzz = cbuffer.data(hh_off + 433 * ccomps * dcomps);

            auto g_zzzzz_xzzzz = cbuffer.data(hh_off + 434 * ccomps * dcomps);

            auto g_zzzzz_yyyyy = cbuffer.data(hh_off + 435 * ccomps * dcomps);

            auto g_zzzzz_yyyyz = cbuffer.data(hh_off + 436 * ccomps * dcomps);

            auto g_zzzzz_yyyzz = cbuffer.data(hh_off + 437 * ccomps * dcomps);

            auto g_zzzzz_yyzzz = cbuffer.data(hh_off + 438 * ccomps * dcomps);

            auto g_zzzzz_yzzzz = cbuffer.data(hh_off + 439 * ccomps * dcomps);

            auto g_zzzzz_zzzzz = cbuffer.data(hh_off + 440 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HHSS

            const auto hh_geom_10_off = idx_geom_10_hhxx + i * dcomps + j;

            auto g_x_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 419 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 434 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 755 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 839 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 860 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 861 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 862 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 863 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 867 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 944 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 951 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 959 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 974 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 979 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 989 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 1004 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 1007 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1019 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1034 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1035 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1049 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1063 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1064 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1079 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1091 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1094 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1108 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1109 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1110 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1111 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1112 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1113 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1114 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1115 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1119 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1124 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1139 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1144 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1145 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1146 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1147 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1148 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1149 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1150 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1151 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1154 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1169 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1259 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1280 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1281 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1282 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1283 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1284 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1285 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1286 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1287 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1288 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1289 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1290 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1291 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1292 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1293 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1294 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1295 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1301 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1309 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1310 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1311 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1312 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1313 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1314 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1315 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1322 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HISS

            const auto hi_geom_10_off = idx_geom_10_hixx + i * dcomps + j;

            auto g_x_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 419 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 441 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 446 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 452 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 459 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 461 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 467 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 470 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 473 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 477 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 482 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 488 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 495 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 539 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 540 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 542 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 543 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 544 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 545 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 547 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 548 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 549 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 550 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 551 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 553 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 554 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 555 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 556 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 557 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 558 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 561 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 562 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 563 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 564 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 565 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 566 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 567 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 568 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 569 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 570 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 571 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 572 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 573 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 574 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 575 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 576 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 577 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 578 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 579 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 580 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 581 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 582 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 583 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 584 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 585 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 586 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 860 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 881 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 896 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 897 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 898 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 899 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 902 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 928 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 929 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 930 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 931 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 932 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 933 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 934 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 935 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 944 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 964 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 965 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 966 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 967 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 968 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 969 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 970 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 971 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 986 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1000 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1028 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1035 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1036 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1037 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1038 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1039 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1040 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1041 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1042 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1043 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1049 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1063 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1070 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1072 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1073 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1074 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1075 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1076 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1077 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1078 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1079 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1080 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1081 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1082 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1083 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1084 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1086 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1087 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1088 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1089 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1090 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1091 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1092 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1093 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1094 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1095 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1096 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1097 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1098 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1099 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1100 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1101 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1102 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1103 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1104 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1105 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1106 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1107 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1108 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1109 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1110 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1111 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1112 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1114 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1115 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1116 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1117 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1118 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1119 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1120 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1121 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1122 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1123 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1124 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1125 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1126 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1127 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1128 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1129 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1130 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1131 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1132 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1133 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1134 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1135 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1136 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1137 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1138 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1139 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1140 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1142 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1143 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1144 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1145 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1146 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1147 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1148 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1149 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1150 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1151 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1152 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1153 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1154 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1155 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1156 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1157 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1158 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1159 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1160 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1161 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1162 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1163 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1164 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1165 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1166 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1167 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1168 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1170 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1171 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1172 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1173 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1174 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1280 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1288 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1289 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1290 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1291 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1292 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1293 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1294 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1295 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1301 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1322 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1323 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1324 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1325 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1326 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1327 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1328 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1329 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1330 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1331 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1332 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1333 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1334 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1335 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1336 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1344 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1345 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1346 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1347 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1348 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1349 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1350 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1351 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1352 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1353 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1354 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1355 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1356 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1357 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1358 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1359 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1360 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1361 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1362 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1363 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1364 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1372 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1373 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1374 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1375 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1376 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1377 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1378 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1379 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1380 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1381 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1382 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1383 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1384 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1385 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1386 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1387 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1388 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1389 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1390 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1391 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1392 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1400 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1401 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1402 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1403 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1404 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1405 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1406 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1407 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1408 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1409 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1410 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1411 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1412 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1413 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1414 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1415 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1416 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1417 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1418 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1419 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1420 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1428 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1429 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1430 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1431 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1432 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1433 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1434 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1435 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1436 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1437 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1438 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1439 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1440 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1441 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1442 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1443 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1444 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1445 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1446 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1447 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1448 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1456 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1457 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1458 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1459 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1460 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1461 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1462 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1463 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1464 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1465 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1466 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1467 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1468 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1469 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1470 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1471 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1472 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1473 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1474 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1475 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1476 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1484 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1485 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1486 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1487 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1488 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1489 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1490 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1491 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1492 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1493 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1494 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1495 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1496 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1497 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1498 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1499 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1500 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1501 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1502 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1503 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1504 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1512 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1513 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1514 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1515 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1516 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1517 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1518 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1519 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1520 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1521 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1522 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1523 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1524 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1525 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1526 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1527 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1528 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1529 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1530 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1531 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1532 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1540 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1541 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1542 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1543 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1544 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1545 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1546 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1547 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1548 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1549 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1550 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1551 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1552 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1553 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1554 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1555 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1556 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1557 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1558 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1559 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1560 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1574 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1576 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1577 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1578 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1579 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1580 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1581 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1582 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1583 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1612 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1613 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1614 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1615 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1616 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1617 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1618 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1619 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1620 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1621 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1622 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1624 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1625 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1626 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1627 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1628 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1629 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1630 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1631 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1632 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1633 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1634 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1635 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1636 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1637 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1638 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1639 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1640 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1641 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1642 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1643 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1644 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1645 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1646 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1647 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1648 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1649 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1650 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1652 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1653 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1654 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1655 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1656 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1657 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1658 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1659 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1660 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1661 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1662 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1663 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1664 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1665 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1666 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1667 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1668 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1669 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1670 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1671 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1672 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1673 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1674 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1675 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1676 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1677 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1678 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1680 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1681 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1682 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1683 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1684 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1685 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1686 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1687 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1688 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1689 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1690 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1691 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1692 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1693 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1694 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1695 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1696 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1697 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1698 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1699 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1700 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1701 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1702 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1703 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1704 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1705 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1706 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1708 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1709 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1710 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1711 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1712 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1713 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1714 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1715 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1716 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1717 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1718 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1719 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1720 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1721 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1722 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1723 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1724 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1725 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1726 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1727 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1728 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1729 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1730 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1731 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1732 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1733 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1734 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1736 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1737 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1738 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1739 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1740 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1741 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1742 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1743 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1744 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1745 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1746 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1747 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1748 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1749 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1750 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1751 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1752 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1753 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1754 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1755 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1756 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1757 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1758 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1759 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1760 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1761 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1762 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1763 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ihxx

            const auto ih_geom_10_off = idx_geom_10_ihxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxx_xxxxx, g_x_0_xxxxxx_xxxxy, g_x_0_xxxxxx_xxxxz, g_x_0_xxxxxx_xxxyy, g_x_0_xxxxxx_xxxyz, g_x_0_xxxxxx_xxxzz, g_x_0_xxxxxx_xxyyy, g_x_0_xxxxxx_xxyyz, g_x_0_xxxxxx_xxyzz, g_x_0_xxxxxx_xxzzz, g_x_0_xxxxxx_xyyyy, g_x_0_xxxxxx_xyyyz, g_x_0_xxxxxx_xyyzz, g_x_0_xxxxxx_xyzzz, g_x_0_xxxxxx_xzzzz, g_x_0_xxxxxx_yyyyy, g_x_0_xxxxxx_yyyyz, g_x_0_xxxxxx_yyyzz, g_x_0_xxxxxx_yyzzz, g_x_0_xxxxxx_yzzzz, g_x_0_xxxxxx_zzzzz, g_xxxxx_xxxxx, g_xxxxx_xxxxy, g_xxxxx_xxxxz, g_xxxxx_xxxyy, g_xxxxx_xxxyz, g_xxxxx_xxxzz, g_xxxxx_xxyyy, g_xxxxx_xxyyz, g_xxxxx_xxyzz, g_xxxxx_xxzzz, g_xxxxx_xyyyy, g_xxxxx_xyyyz, g_xxxxx_xyyzz, g_xxxxx_xyzzz, g_xxxxx_xzzzz, g_xxxxx_yyyyy, g_xxxxx_yyyyz, g_xxxxx_yyyzz, g_xxxxx_yyzzz, g_xxxxx_yzzzz, g_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxxxx[k] = -g_xxxxx_xxxxx[k] - g_x_0_xxxxx_xxxxx[k] * ab_x + g_x_0_xxxxx_xxxxxx[k];

                g_x_0_xxxxxx_xxxxy[k] = -g_xxxxx_xxxxy[k] - g_x_0_xxxxx_xxxxy[k] * ab_x + g_x_0_xxxxx_xxxxxy[k];

                g_x_0_xxxxxx_xxxxz[k] = -g_xxxxx_xxxxz[k] - g_x_0_xxxxx_xxxxz[k] * ab_x + g_x_0_xxxxx_xxxxxz[k];

                g_x_0_xxxxxx_xxxyy[k] = -g_xxxxx_xxxyy[k] - g_x_0_xxxxx_xxxyy[k] * ab_x + g_x_0_xxxxx_xxxxyy[k];

                g_x_0_xxxxxx_xxxyz[k] = -g_xxxxx_xxxyz[k] - g_x_0_xxxxx_xxxyz[k] * ab_x + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxx_xxxzz[k] = -g_xxxxx_xxxzz[k] - g_x_0_xxxxx_xxxzz[k] * ab_x + g_x_0_xxxxx_xxxxzz[k];

                g_x_0_xxxxxx_xxyyy[k] = -g_xxxxx_xxyyy[k] - g_x_0_xxxxx_xxyyy[k] * ab_x + g_x_0_xxxxx_xxxyyy[k];

                g_x_0_xxxxxx_xxyyz[k] = -g_xxxxx_xxyyz[k] - g_x_0_xxxxx_xxyyz[k] * ab_x + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxx_xxyzz[k] = -g_xxxxx_xxyzz[k] - g_x_0_xxxxx_xxyzz[k] * ab_x + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxx_xxzzz[k] = -g_xxxxx_xxzzz[k] - g_x_0_xxxxx_xxzzz[k] * ab_x + g_x_0_xxxxx_xxxzzz[k];

                g_x_0_xxxxxx_xyyyy[k] = -g_xxxxx_xyyyy[k] - g_x_0_xxxxx_xyyyy[k] * ab_x + g_x_0_xxxxx_xxyyyy[k];

                g_x_0_xxxxxx_xyyyz[k] = -g_xxxxx_xyyyz[k] - g_x_0_xxxxx_xyyyz[k] * ab_x + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxx_xyyzz[k] = -g_xxxxx_xyyzz[k] - g_x_0_xxxxx_xyyzz[k] * ab_x + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxx_xyzzz[k] = -g_xxxxx_xyzzz[k] - g_x_0_xxxxx_xyzzz[k] * ab_x + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxx_xzzzz[k] = -g_xxxxx_xzzzz[k] - g_x_0_xxxxx_xzzzz[k] * ab_x + g_x_0_xxxxx_xxzzzz[k];

                g_x_0_xxxxxx_yyyyy[k] = -g_xxxxx_yyyyy[k] - g_x_0_xxxxx_yyyyy[k] * ab_x + g_x_0_xxxxx_xyyyyy[k];

                g_x_0_xxxxxx_yyyyz[k] = -g_xxxxx_yyyyz[k] - g_x_0_xxxxx_yyyyz[k] * ab_x + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxx_yyyzz[k] = -g_xxxxx_yyyzz[k] - g_x_0_xxxxx_yyyzz[k] * ab_x + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxx_yyzzz[k] = -g_xxxxx_yyzzz[k] - g_x_0_xxxxx_yyzzz[k] * ab_x + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxx_yzzzz[k] = -g_xxxxx_yzzzz[k] - g_x_0_xxxxx_yzzzz[k] * ab_x + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxx_zzzzz[k] = -g_xxxxx_zzzzz[k] - g_x_0_xxxxx_zzzzz[k] * ab_x + g_x_0_xxxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxy_xxxxx, g_x_0_xxxxxy_xxxxy, g_x_0_xxxxxy_xxxxz, g_x_0_xxxxxy_xxxyy, g_x_0_xxxxxy_xxxyz, g_x_0_xxxxxy_xxxzz, g_x_0_xxxxxy_xxyyy, g_x_0_xxxxxy_xxyyz, g_x_0_xxxxxy_xxyzz, g_x_0_xxxxxy_xxzzz, g_x_0_xxxxxy_xyyyy, g_x_0_xxxxxy_xyyyz, g_x_0_xxxxxy_xyyzz, g_x_0_xxxxxy_xyzzz, g_x_0_xxxxxy_xzzzz, g_x_0_xxxxxy_yyyyy, g_x_0_xxxxxy_yyyyz, g_x_0_xxxxxy_yyyzz, g_x_0_xxxxxy_yyzzz, g_x_0_xxxxxy_yzzzz, g_x_0_xxxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxxxx[k] = -g_x_0_xxxxx_xxxxx[k] * ab_y + g_x_0_xxxxx_xxxxxy[k];

                g_x_0_xxxxxy_xxxxy[k] = -g_x_0_xxxxx_xxxxy[k] * ab_y + g_x_0_xxxxx_xxxxyy[k];

                g_x_0_xxxxxy_xxxxz[k] = -g_x_0_xxxxx_xxxxz[k] * ab_y + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxy_xxxyy[k] = -g_x_0_xxxxx_xxxyy[k] * ab_y + g_x_0_xxxxx_xxxyyy[k];

                g_x_0_xxxxxy_xxxyz[k] = -g_x_0_xxxxx_xxxyz[k] * ab_y + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxy_xxxzz[k] = -g_x_0_xxxxx_xxxzz[k] * ab_y + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxy_xxyyy[k] = -g_x_0_xxxxx_xxyyy[k] * ab_y + g_x_0_xxxxx_xxyyyy[k];

                g_x_0_xxxxxy_xxyyz[k] = -g_x_0_xxxxx_xxyyz[k] * ab_y + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxy_xxyzz[k] = -g_x_0_xxxxx_xxyzz[k] * ab_y + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxy_xxzzz[k] = -g_x_0_xxxxx_xxzzz[k] * ab_y + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxy_xyyyy[k] = -g_x_0_xxxxx_xyyyy[k] * ab_y + g_x_0_xxxxx_xyyyyy[k];

                g_x_0_xxxxxy_xyyyz[k] = -g_x_0_xxxxx_xyyyz[k] * ab_y + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxy_xyyzz[k] = -g_x_0_xxxxx_xyyzz[k] * ab_y + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxy_xyzzz[k] = -g_x_0_xxxxx_xyzzz[k] * ab_y + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxy_xzzzz[k] = -g_x_0_xxxxx_xzzzz[k] * ab_y + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxy_yyyyy[k] = -g_x_0_xxxxx_yyyyy[k] * ab_y + g_x_0_xxxxx_yyyyyy[k];

                g_x_0_xxxxxy_yyyyz[k] = -g_x_0_xxxxx_yyyyz[k] * ab_y + g_x_0_xxxxx_yyyyyz[k];

                g_x_0_xxxxxy_yyyzz[k] = -g_x_0_xxxxx_yyyzz[k] * ab_y + g_x_0_xxxxx_yyyyzz[k];

                g_x_0_xxxxxy_yyzzz[k] = -g_x_0_xxxxx_yyzzz[k] * ab_y + g_x_0_xxxxx_yyyzzz[k];

                g_x_0_xxxxxy_yzzzz[k] = -g_x_0_xxxxx_yzzzz[k] * ab_y + g_x_0_xxxxx_yyzzzz[k];

                g_x_0_xxxxxy_zzzzz[k] = -g_x_0_xxxxx_zzzzz[k] * ab_y + g_x_0_xxxxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxx_zzzzzz, g_x_0_xxxxxz_xxxxx, g_x_0_xxxxxz_xxxxy, g_x_0_xxxxxz_xxxxz, g_x_0_xxxxxz_xxxyy, g_x_0_xxxxxz_xxxyz, g_x_0_xxxxxz_xxxzz, g_x_0_xxxxxz_xxyyy, g_x_0_xxxxxz_xxyyz, g_x_0_xxxxxz_xxyzz, g_x_0_xxxxxz_xxzzz, g_x_0_xxxxxz_xyyyy, g_x_0_xxxxxz_xyyyz, g_x_0_xxxxxz_xyyzz, g_x_0_xxxxxz_xyzzz, g_x_0_xxxxxz_xzzzz, g_x_0_xxxxxz_yyyyy, g_x_0_xxxxxz_yyyyz, g_x_0_xxxxxz_yyyzz, g_x_0_xxxxxz_yyzzz, g_x_0_xxxxxz_yzzzz, g_x_0_xxxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxxxx[k] = -g_x_0_xxxxx_xxxxx[k] * ab_z + g_x_0_xxxxx_xxxxxz[k];

                g_x_0_xxxxxz_xxxxy[k] = -g_x_0_xxxxx_xxxxy[k] * ab_z + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxz_xxxxz[k] = -g_x_0_xxxxx_xxxxz[k] * ab_z + g_x_0_xxxxx_xxxxzz[k];

                g_x_0_xxxxxz_xxxyy[k] = -g_x_0_xxxxx_xxxyy[k] * ab_z + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxz_xxxyz[k] = -g_x_0_xxxxx_xxxyz[k] * ab_z + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxz_xxxzz[k] = -g_x_0_xxxxx_xxxzz[k] * ab_z + g_x_0_xxxxx_xxxzzz[k];

                g_x_0_xxxxxz_xxyyy[k] = -g_x_0_xxxxx_xxyyy[k] * ab_z + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxz_xxyyz[k] = -g_x_0_xxxxx_xxyyz[k] * ab_z + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxz_xxyzz[k] = -g_x_0_xxxxx_xxyzz[k] * ab_z + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxz_xxzzz[k] = -g_x_0_xxxxx_xxzzz[k] * ab_z + g_x_0_xxxxx_xxzzzz[k];

                g_x_0_xxxxxz_xyyyy[k] = -g_x_0_xxxxx_xyyyy[k] * ab_z + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxz_xyyyz[k] = -g_x_0_xxxxx_xyyyz[k] * ab_z + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxz_xyyzz[k] = -g_x_0_xxxxx_xyyzz[k] * ab_z + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxz_xyzzz[k] = -g_x_0_xxxxx_xyzzz[k] * ab_z + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxz_xzzzz[k] = -g_x_0_xxxxx_xzzzz[k] * ab_z + g_x_0_xxxxx_xzzzzz[k];

                g_x_0_xxxxxz_yyyyy[k] = -g_x_0_xxxxx_yyyyy[k] * ab_z + g_x_0_xxxxx_yyyyyz[k];

                g_x_0_xxxxxz_yyyyz[k] = -g_x_0_xxxxx_yyyyz[k] * ab_z + g_x_0_xxxxx_yyyyzz[k];

                g_x_0_xxxxxz_yyyzz[k] = -g_x_0_xxxxx_yyyzz[k] * ab_z + g_x_0_xxxxx_yyyzzz[k];

                g_x_0_xxxxxz_yyzzz[k] = -g_x_0_xxxxx_yyzzz[k] * ab_z + g_x_0_xxxxx_yyzzzz[k];

                g_x_0_xxxxxz_yzzzz[k] = -g_x_0_xxxxx_yzzzz[k] * ab_z + g_x_0_xxxxx_yzzzzz[k];

                g_x_0_xxxxxz_zzzzz[k] = -g_x_0_xxxxx_zzzzz[k] * ab_z + g_x_0_xxxxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxy_xxxxx, g_x_0_xxxxy_xxxxxy, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxxyy, g_x_0_xxxxy_xxxxyz, g_x_0_xxxxy_xxxxz, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyyy, g_x_0_xxxxy_xxxyyz, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxyzz, g_x_0_xxxxy_xxxzz, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyyy, g_x_0_xxxxy_xxyyyz, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyyzz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxyzzz, g_x_0_xxxxy_xxzzz, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyyy, g_x_0_xxxxy_xyyyyz, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyyzz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyyzzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xyzzzz, g_x_0_xxxxy_xzzzz, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyyy, g_x_0_xxxxy_yyyyyz, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyyzz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyyzzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yyzzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_yzzzzz, g_x_0_xxxxy_zzzzz, g_x_0_xxxxyy_xxxxx, g_x_0_xxxxyy_xxxxy, g_x_0_xxxxyy_xxxxz, g_x_0_xxxxyy_xxxyy, g_x_0_xxxxyy_xxxyz, g_x_0_xxxxyy_xxxzz, g_x_0_xxxxyy_xxyyy, g_x_0_xxxxyy_xxyyz, g_x_0_xxxxyy_xxyzz, g_x_0_xxxxyy_xxzzz, g_x_0_xxxxyy_xyyyy, g_x_0_xxxxyy_xyyyz, g_x_0_xxxxyy_xyyzz, g_x_0_xxxxyy_xyzzz, g_x_0_xxxxyy_xzzzz, g_x_0_xxxxyy_yyyyy, g_x_0_xxxxyy_yyyyz, g_x_0_xxxxyy_yyyzz, g_x_0_xxxxyy_yyzzz, g_x_0_xxxxyy_yzzzz, g_x_0_xxxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxxxx[k] = -g_x_0_xxxxy_xxxxx[k] * ab_y + g_x_0_xxxxy_xxxxxy[k];

                g_x_0_xxxxyy_xxxxy[k] = -g_x_0_xxxxy_xxxxy[k] * ab_y + g_x_0_xxxxy_xxxxyy[k];

                g_x_0_xxxxyy_xxxxz[k] = -g_x_0_xxxxy_xxxxz[k] * ab_y + g_x_0_xxxxy_xxxxyz[k];

                g_x_0_xxxxyy_xxxyy[k] = -g_x_0_xxxxy_xxxyy[k] * ab_y + g_x_0_xxxxy_xxxyyy[k];

                g_x_0_xxxxyy_xxxyz[k] = -g_x_0_xxxxy_xxxyz[k] * ab_y + g_x_0_xxxxy_xxxyyz[k];

                g_x_0_xxxxyy_xxxzz[k] = -g_x_0_xxxxy_xxxzz[k] * ab_y + g_x_0_xxxxy_xxxyzz[k];

                g_x_0_xxxxyy_xxyyy[k] = -g_x_0_xxxxy_xxyyy[k] * ab_y + g_x_0_xxxxy_xxyyyy[k];

                g_x_0_xxxxyy_xxyyz[k] = -g_x_0_xxxxy_xxyyz[k] * ab_y + g_x_0_xxxxy_xxyyyz[k];

                g_x_0_xxxxyy_xxyzz[k] = -g_x_0_xxxxy_xxyzz[k] * ab_y + g_x_0_xxxxy_xxyyzz[k];

                g_x_0_xxxxyy_xxzzz[k] = -g_x_0_xxxxy_xxzzz[k] * ab_y + g_x_0_xxxxy_xxyzzz[k];

                g_x_0_xxxxyy_xyyyy[k] = -g_x_0_xxxxy_xyyyy[k] * ab_y + g_x_0_xxxxy_xyyyyy[k];

                g_x_0_xxxxyy_xyyyz[k] = -g_x_0_xxxxy_xyyyz[k] * ab_y + g_x_0_xxxxy_xyyyyz[k];

                g_x_0_xxxxyy_xyyzz[k] = -g_x_0_xxxxy_xyyzz[k] * ab_y + g_x_0_xxxxy_xyyyzz[k];

                g_x_0_xxxxyy_xyzzz[k] = -g_x_0_xxxxy_xyzzz[k] * ab_y + g_x_0_xxxxy_xyyzzz[k];

                g_x_0_xxxxyy_xzzzz[k] = -g_x_0_xxxxy_xzzzz[k] * ab_y + g_x_0_xxxxy_xyzzzz[k];

                g_x_0_xxxxyy_yyyyy[k] = -g_x_0_xxxxy_yyyyy[k] * ab_y + g_x_0_xxxxy_yyyyyy[k];

                g_x_0_xxxxyy_yyyyz[k] = -g_x_0_xxxxy_yyyyz[k] * ab_y + g_x_0_xxxxy_yyyyyz[k];

                g_x_0_xxxxyy_yyyzz[k] = -g_x_0_xxxxy_yyyzz[k] * ab_y + g_x_0_xxxxy_yyyyzz[k];

                g_x_0_xxxxyy_yyzzz[k] = -g_x_0_xxxxy_yyzzz[k] * ab_y + g_x_0_xxxxy_yyyzzz[k];

                g_x_0_xxxxyy_yzzzz[k] = -g_x_0_xxxxy_yzzzz[k] * ab_y + g_x_0_xxxxy_yyzzzz[k];

                g_x_0_xxxxyy_zzzzz[k] = -g_x_0_xxxxy_zzzzz[k] * ab_y + g_x_0_xxxxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxyz_xxxxx, g_x_0_xxxxyz_xxxxy, g_x_0_xxxxyz_xxxxz, g_x_0_xxxxyz_xxxyy, g_x_0_xxxxyz_xxxyz, g_x_0_xxxxyz_xxxzz, g_x_0_xxxxyz_xxyyy, g_x_0_xxxxyz_xxyyz, g_x_0_xxxxyz_xxyzz, g_x_0_xxxxyz_xxzzz, g_x_0_xxxxyz_xyyyy, g_x_0_xxxxyz_xyyyz, g_x_0_xxxxyz_xyyzz, g_x_0_xxxxyz_xyzzz, g_x_0_xxxxyz_xzzzz, g_x_0_xxxxyz_yyyyy, g_x_0_xxxxyz_yyyyz, g_x_0_xxxxyz_yyyzz, g_x_0_xxxxyz_yyzzz, g_x_0_xxxxyz_yzzzz, g_x_0_xxxxyz_zzzzz, g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxxxx[k] = -g_x_0_xxxxz_xxxxx[k] * ab_y + g_x_0_xxxxz_xxxxxy[k];

                g_x_0_xxxxyz_xxxxy[k] = -g_x_0_xxxxz_xxxxy[k] * ab_y + g_x_0_xxxxz_xxxxyy[k];

                g_x_0_xxxxyz_xxxxz[k] = -g_x_0_xxxxz_xxxxz[k] * ab_y + g_x_0_xxxxz_xxxxyz[k];

                g_x_0_xxxxyz_xxxyy[k] = -g_x_0_xxxxz_xxxyy[k] * ab_y + g_x_0_xxxxz_xxxyyy[k];

                g_x_0_xxxxyz_xxxyz[k] = -g_x_0_xxxxz_xxxyz[k] * ab_y + g_x_0_xxxxz_xxxyyz[k];

                g_x_0_xxxxyz_xxxzz[k] = -g_x_0_xxxxz_xxxzz[k] * ab_y + g_x_0_xxxxz_xxxyzz[k];

                g_x_0_xxxxyz_xxyyy[k] = -g_x_0_xxxxz_xxyyy[k] * ab_y + g_x_0_xxxxz_xxyyyy[k];

                g_x_0_xxxxyz_xxyyz[k] = -g_x_0_xxxxz_xxyyz[k] * ab_y + g_x_0_xxxxz_xxyyyz[k];

                g_x_0_xxxxyz_xxyzz[k] = -g_x_0_xxxxz_xxyzz[k] * ab_y + g_x_0_xxxxz_xxyyzz[k];

                g_x_0_xxxxyz_xxzzz[k] = -g_x_0_xxxxz_xxzzz[k] * ab_y + g_x_0_xxxxz_xxyzzz[k];

                g_x_0_xxxxyz_xyyyy[k] = -g_x_0_xxxxz_xyyyy[k] * ab_y + g_x_0_xxxxz_xyyyyy[k];

                g_x_0_xxxxyz_xyyyz[k] = -g_x_0_xxxxz_xyyyz[k] * ab_y + g_x_0_xxxxz_xyyyyz[k];

                g_x_0_xxxxyz_xyyzz[k] = -g_x_0_xxxxz_xyyzz[k] * ab_y + g_x_0_xxxxz_xyyyzz[k];

                g_x_0_xxxxyz_xyzzz[k] = -g_x_0_xxxxz_xyzzz[k] * ab_y + g_x_0_xxxxz_xyyzzz[k];

                g_x_0_xxxxyz_xzzzz[k] = -g_x_0_xxxxz_xzzzz[k] * ab_y + g_x_0_xxxxz_xyzzzz[k];

                g_x_0_xxxxyz_yyyyy[k] = -g_x_0_xxxxz_yyyyy[k] * ab_y + g_x_0_xxxxz_yyyyyy[k];

                g_x_0_xxxxyz_yyyyz[k] = -g_x_0_xxxxz_yyyyz[k] * ab_y + g_x_0_xxxxz_yyyyyz[k];

                g_x_0_xxxxyz_yyyzz[k] = -g_x_0_xxxxz_yyyzz[k] * ab_y + g_x_0_xxxxz_yyyyzz[k];

                g_x_0_xxxxyz_yyzzz[k] = -g_x_0_xxxxz_yyzzz[k] * ab_y + g_x_0_xxxxz_yyyzzz[k];

                g_x_0_xxxxyz_yzzzz[k] = -g_x_0_xxxxz_yzzzz[k] * ab_y + g_x_0_xxxxz_yyzzzz[k];

                g_x_0_xxxxyz_zzzzz[k] = -g_x_0_xxxxz_zzzzz[k] * ab_y + g_x_0_xxxxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzz, g_x_0_xxxxz_zzzzzz, g_x_0_xxxxzz_xxxxx, g_x_0_xxxxzz_xxxxy, g_x_0_xxxxzz_xxxxz, g_x_0_xxxxzz_xxxyy, g_x_0_xxxxzz_xxxyz, g_x_0_xxxxzz_xxxzz, g_x_0_xxxxzz_xxyyy, g_x_0_xxxxzz_xxyyz, g_x_0_xxxxzz_xxyzz, g_x_0_xxxxzz_xxzzz, g_x_0_xxxxzz_xyyyy, g_x_0_xxxxzz_xyyyz, g_x_0_xxxxzz_xyyzz, g_x_0_xxxxzz_xyzzz, g_x_0_xxxxzz_xzzzz, g_x_0_xxxxzz_yyyyy, g_x_0_xxxxzz_yyyyz, g_x_0_xxxxzz_yyyzz, g_x_0_xxxxzz_yyzzz, g_x_0_xxxxzz_yzzzz, g_x_0_xxxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxxxx[k] = -g_x_0_xxxxz_xxxxx[k] * ab_z + g_x_0_xxxxz_xxxxxz[k];

                g_x_0_xxxxzz_xxxxy[k] = -g_x_0_xxxxz_xxxxy[k] * ab_z + g_x_0_xxxxz_xxxxyz[k];

                g_x_0_xxxxzz_xxxxz[k] = -g_x_0_xxxxz_xxxxz[k] * ab_z + g_x_0_xxxxz_xxxxzz[k];

                g_x_0_xxxxzz_xxxyy[k] = -g_x_0_xxxxz_xxxyy[k] * ab_z + g_x_0_xxxxz_xxxyyz[k];

                g_x_0_xxxxzz_xxxyz[k] = -g_x_0_xxxxz_xxxyz[k] * ab_z + g_x_0_xxxxz_xxxyzz[k];

                g_x_0_xxxxzz_xxxzz[k] = -g_x_0_xxxxz_xxxzz[k] * ab_z + g_x_0_xxxxz_xxxzzz[k];

                g_x_0_xxxxzz_xxyyy[k] = -g_x_0_xxxxz_xxyyy[k] * ab_z + g_x_0_xxxxz_xxyyyz[k];

                g_x_0_xxxxzz_xxyyz[k] = -g_x_0_xxxxz_xxyyz[k] * ab_z + g_x_0_xxxxz_xxyyzz[k];

                g_x_0_xxxxzz_xxyzz[k] = -g_x_0_xxxxz_xxyzz[k] * ab_z + g_x_0_xxxxz_xxyzzz[k];

                g_x_0_xxxxzz_xxzzz[k] = -g_x_0_xxxxz_xxzzz[k] * ab_z + g_x_0_xxxxz_xxzzzz[k];

                g_x_0_xxxxzz_xyyyy[k] = -g_x_0_xxxxz_xyyyy[k] * ab_z + g_x_0_xxxxz_xyyyyz[k];

                g_x_0_xxxxzz_xyyyz[k] = -g_x_0_xxxxz_xyyyz[k] * ab_z + g_x_0_xxxxz_xyyyzz[k];

                g_x_0_xxxxzz_xyyzz[k] = -g_x_0_xxxxz_xyyzz[k] * ab_z + g_x_0_xxxxz_xyyzzz[k];

                g_x_0_xxxxzz_xyzzz[k] = -g_x_0_xxxxz_xyzzz[k] * ab_z + g_x_0_xxxxz_xyzzzz[k];

                g_x_0_xxxxzz_xzzzz[k] = -g_x_0_xxxxz_xzzzz[k] * ab_z + g_x_0_xxxxz_xzzzzz[k];

                g_x_0_xxxxzz_yyyyy[k] = -g_x_0_xxxxz_yyyyy[k] * ab_z + g_x_0_xxxxz_yyyyyz[k];

                g_x_0_xxxxzz_yyyyz[k] = -g_x_0_xxxxz_yyyyz[k] * ab_z + g_x_0_xxxxz_yyyyzz[k];

                g_x_0_xxxxzz_yyyzz[k] = -g_x_0_xxxxz_yyyzz[k] * ab_z + g_x_0_xxxxz_yyyzzz[k];

                g_x_0_xxxxzz_yyzzz[k] = -g_x_0_xxxxz_yyzzz[k] * ab_z + g_x_0_xxxxz_yyzzzz[k];

                g_x_0_xxxxzz_yzzzz[k] = -g_x_0_xxxxz_yzzzz[k] * ab_z + g_x_0_xxxxz_yzzzzz[k];

                g_x_0_xxxxzz_zzzzz[k] = -g_x_0_xxxxz_zzzzz[k] * ab_z + g_x_0_xxxxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyy_xxxxx, g_x_0_xxxyy_xxxxxy, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxxyy, g_x_0_xxxyy_xxxxyz, g_x_0_xxxyy_xxxxz, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyyy, g_x_0_xxxyy_xxxyyz, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxyzz, g_x_0_xxxyy_xxxzz, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyyy, g_x_0_xxxyy_xxyyyz, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyyzz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxyzzz, g_x_0_xxxyy_xxzzz, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyyy, g_x_0_xxxyy_xyyyyz, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyyzz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyyzzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xyzzzz, g_x_0_xxxyy_xzzzz, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyyy, g_x_0_xxxyy_yyyyyz, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyyzz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyyzzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yyzzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_yzzzzz, g_x_0_xxxyy_zzzzz, g_x_0_xxxyyy_xxxxx, g_x_0_xxxyyy_xxxxy, g_x_0_xxxyyy_xxxxz, g_x_0_xxxyyy_xxxyy, g_x_0_xxxyyy_xxxyz, g_x_0_xxxyyy_xxxzz, g_x_0_xxxyyy_xxyyy, g_x_0_xxxyyy_xxyyz, g_x_0_xxxyyy_xxyzz, g_x_0_xxxyyy_xxzzz, g_x_0_xxxyyy_xyyyy, g_x_0_xxxyyy_xyyyz, g_x_0_xxxyyy_xyyzz, g_x_0_xxxyyy_xyzzz, g_x_0_xxxyyy_xzzzz, g_x_0_xxxyyy_yyyyy, g_x_0_xxxyyy_yyyyz, g_x_0_xxxyyy_yyyzz, g_x_0_xxxyyy_yyzzz, g_x_0_xxxyyy_yzzzz, g_x_0_xxxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxxxx[k] = -g_x_0_xxxyy_xxxxx[k] * ab_y + g_x_0_xxxyy_xxxxxy[k];

                g_x_0_xxxyyy_xxxxy[k] = -g_x_0_xxxyy_xxxxy[k] * ab_y + g_x_0_xxxyy_xxxxyy[k];

                g_x_0_xxxyyy_xxxxz[k] = -g_x_0_xxxyy_xxxxz[k] * ab_y + g_x_0_xxxyy_xxxxyz[k];

                g_x_0_xxxyyy_xxxyy[k] = -g_x_0_xxxyy_xxxyy[k] * ab_y + g_x_0_xxxyy_xxxyyy[k];

                g_x_0_xxxyyy_xxxyz[k] = -g_x_0_xxxyy_xxxyz[k] * ab_y + g_x_0_xxxyy_xxxyyz[k];

                g_x_0_xxxyyy_xxxzz[k] = -g_x_0_xxxyy_xxxzz[k] * ab_y + g_x_0_xxxyy_xxxyzz[k];

                g_x_0_xxxyyy_xxyyy[k] = -g_x_0_xxxyy_xxyyy[k] * ab_y + g_x_0_xxxyy_xxyyyy[k];

                g_x_0_xxxyyy_xxyyz[k] = -g_x_0_xxxyy_xxyyz[k] * ab_y + g_x_0_xxxyy_xxyyyz[k];

                g_x_0_xxxyyy_xxyzz[k] = -g_x_0_xxxyy_xxyzz[k] * ab_y + g_x_0_xxxyy_xxyyzz[k];

                g_x_0_xxxyyy_xxzzz[k] = -g_x_0_xxxyy_xxzzz[k] * ab_y + g_x_0_xxxyy_xxyzzz[k];

                g_x_0_xxxyyy_xyyyy[k] = -g_x_0_xxxyy_xyyyy[k] * ab_y + g_x_0_xxxyy_xyyyyy[k];

                g_x_0_xxxyyy_xyyyz[k] = -g_x_0_xxxyy_xyyyz[k] * ab_y + g_x_0_xxxyy_xyyyyz[k];

                g_x_0_xxxyyy_xyyzz[k] = -g_x_0_xxxyy_xyyzz[k] * ab_y + g_x_0_xxxyy_xyyyzz[k];

                g_x_0_xxxyyy_xyzzz[k] = -g_x_0_xxxyy_xyzzz[k] * ab_y + g_x_0_xxxyy_xyyzzz[k];

                g_x_0_xxxyyy_xzzzz[k] = -g_x_0_xxxyy_xzzzz[k] * ab_y + g_x_0_xxxyy_xyzzzz[k];

                g_x_0_xxxyyy_yyyyy[k] = -g_x_0_xxxyy_yyyyy[k] * ab_y + g_x_0_xxxyy_yyyyyy[k];

                g_x_0_xxxyyy_yyyyz[k] = -g_x_0_xxxyy_yyyyz[k] * ab_y + g_x_0_xxxyy_yyyyyz[k];

                g_x_0_xxxyyy_yyyzz[k] = -g_x_0_xxxyy_yyyzz[k] * ab_y + g_x_0_xxxyy_yyyyzz[k];

                g_x_0_xxxyyy_yyzzz[k] = -g_x_0_xxxyy_yyzzz[k] * ab_y + g_x_0_xxxyy_yyyzzz[k];

                g_x_0_xxxyyy_yzzzz[k] = -g_x_0_xxxyy_yzzzz[k] * ab_y + g_x_0_xxxyy_yyzzzz[k];

                g_x_0_xxxyyy_zzzzz[k] = -g_x_0_xxxyy_zzzzz[k] * ab_y + g_x_0_xxxyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyyz_xxxxx, g_x_0_xxxyyz_xxxxy, g_x_0_xxxyyz_xxxxz, g_x_0_xxxyyz_xxxyy, g_x_0_xxxyyz_xxxyz, g_x_0_xxxyyz_xxxzz, g_x_0_xxxyyz_xxyyy, g_x_0_xxxyyz_xxyyz, g_x_0_xxxyyz_xxyzz, g_x_0_xxxyyz_xxzzz, g_x_0_xxxyyz_xyyyy, g_x_0_xxxyyz_xyyyz, g_x_0_xxxyyz_xyyzz, g_x_0_xxxyyz_xyzzz, g_x_0_xxxyyz_xzzzz, g_x_0_xxxyyz_yyyyy, g_x_0_xxxyyz_yyyyz, g_x_0_xxxyyz_yyyzz, g_x_0_xxxyyz_yyzzz, g_x_0_xxxyyz_yzzzz, g_x_0_xxxyyz_zzzzz, g_x_0_xxxyz_xxxxx, g_x_0_xxxyz_xxxxxy, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxxyy, g_x_0_xxxyz_xxxxyz, g_x_0_xxxyz_xxxxz, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyyy, g_x_0_xxxyz_xxxyyz, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxyzz, g_x_0_xxxyz_xxxzz, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyyy, g_x_0_xxxyz_xxyyyz, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyyzz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxyzzz, g_x_0_xxxyz_xxzzz, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyyy, g_x_0_xxxyz_xyyyyz, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyyzz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyyzzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xyzzzz, g_x_0_xxxyz_xzzzz, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyyy, g_x_0_xxxyz_yyyyyz, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyyzz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyyzzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yyzzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_yzzzzz, g_x_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxxxx[k] = -g_x_0_xxxyz_xxxxx[k] * ab_y + g_x_0_xxxyz_xxxxxy[k];

                g_x_0_xxxyyz_xxxxy[k] = -g_x_0_xxxyz_xxxxy[k] * ab_y + g_x_0_xxxyz_xxxxyy[k];

                g_x_0_xxxyyz_xxxxz[k] = -g_x_0_xxxyz_xxxxz[k] * ab_y + g_x_0_xxxyz_xxxxyz[k];

                g_x_0_xxxyyz_xxxyy[k] = -g_x_0_xxxyz_xxxyy[k] * ab_y + g_x_0_xxxyz_xxxyyy[k];

                g_x_0_xxxyyz_xxxyz[k] = -g_x_0_xxxyz_xxxyz[k] * ab_y + g_x_0_xxxyz_xxxyyz[k];

                g_x_0_xxxyyz_xxxzz[k] = -g_x_0_xxxyz_xxxzz[k] * ab_y + g_x_0_xxxyz_xxxyzz[k];

                g_x_0_xxxyyz_xxyyy[k] = -g_x_0_xxxyz_xxyyy[k] * ab_y + g_x_0_xxxyz_xxyyyy[k];

                g_x_0_xxxyyz_xxyyz[k] = -g_x_0_xxxyz_xxyyz[k] * ab_y + g_x_0_xxxyz_xxyyyz[k];

                g_x_0_xxxyyz_xxyzz[k] = -g_x_0_xxxyz_xxyzz[k] * ab_y + g_x_0_xxxyz_xxyyzz[k];

                g_x_0_xxxyyz_xxzzz[k] = -g_x_0_xxxyz_xxzzz[k] * ab_y + g_x_0_xxxyz_xxyzzz[k];

                g_x_0_xxxyyz_xyyyy[k] = -g_x_0_xxxyz_xyyyy[k] * ab_y + g_x_0_xxxyz_xyyyyy[k];

                g_x_0_xxxyyz_xyyyz[k] = -g_x_0_xxxyz_xyyyz[k] * ab_y + g_x_0_xxxyz_xyyyyz[k];

                g_x_0_xxxyyz_xyyzz[k] = -g_x_0_xxxyz_xyyzz[k] * ab_y + g_x_0_xxxyz_xyyyzz[k];

                g_x_0_xxxyyz_xyzzz[k] = -g_x_0_xxxyz_xyzzz[k] * ab_y + g_x_0_xxxyz_xyyzzz[k];

                g_x_0_xxxyyz_xzzzz[k] = -g_x_0_xxxyz_xzzzz[k] * ab_y + g_x_0_xxxyz_xyzzzz[k];

                g_x_0_xxxyyz_yyyyy[k] = -g_x_0_xxxyz_yyyyy[k] * ab_y + g_x_0_xxxyz_yyyyyy[k];

                g_x_0_xxxyyz_yyyyz[k] = -g_x_0_xxxyz_yyyyz[k] * ab_y + g_x_0_xxxyz_yyyyyz[k];

                g_x_0_xxxyyz_yyyzz[k] = -g_x_0_xxxyz_yyyzz[k] * ab_y + g_x_0_xxxyz_yyyyzz[k];

                g_x_0_xxxyyz_yyzzz[k] = -g_x_0_xxxyz_yyzzz[k] * ab_y + g_x_0_xxxyz_yyyzzz[k];

                g_x_0_xxxyyz_yzzzz[k] = -g_x_0_xxxyz_yzzzz[k] * ab_y + g_x_0_xxxyz_yyzzzz[k];

                g_x_0_xxxyyz_zzzzz[k] = -g_x_0_xxxyz_zzzzz[k] * ab_y + g_x_0_xxxyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyzz_xxxxx, g_x_0_xxxyzz_xxxxy, g_x_0_xxxyzz_xxxxz, g_x_0_xxxyzz_xxxyy, g_x_0_xxxyzz_xxxyz, g_x_0_xxxyzz_xxxzz, g_x_0_xxxyzz_xxyyy, g_x_0_xxxyzz_xxyyz, g_x_0_xxxyzz_xxyzz, g_x_0_xxxyzz_xxzzz, g_x_0_xxxyzz_xyyyy, g_x_0_xxxyzz_xyyyz, g_x_0_xxxyzz_xyyzz, g_x_0_xxxyzz_xyzzz, g_x_0_xxxyzz_xzzzz, g_x_0_xxxyzz_yyyyy, g_x_0_xxxyzz_yyyyz, g_x_0_xxxyzz_yyyzz, g_x_0_xxxyzz_yyzzz, g_x_0_xxxyzz_yzzzz, g_x_0_xxxyzz_zzzzz, g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxxxx[k] = -g_x_0_xxxzz_xxxxx[k] * ab_y + g_x_0_xxxzz_xxxxxy[k];

                g_x_0_xxxyzz_xxxxy[k] = -g_x_0_xxxzz_xxxxy[k] * ab_y + g_x_0_xxxzz_xxxxyy[k];

                g_x_0_xxxyzz_xxxxz[k] = -g_x_0_xxxzz_xxxxz[k] * ab_y + g_x_0_xxxzz_xxxxyz[k];

                g_x_0_xxxyzz_xxxyy[k] = -g_x_0_xxxzz_xxxyy[k] * ab_y + g_x_0_xxxzz_xxxyyy[k];

                g_x_0_xxxyzz_xxxyz[k] = -g_x_0_xxxzz_xxxyz[k] * ab_y + g_x_0_xxxzz_xxxyyz[k];

                g_x_0_xxxyzz_xxxzz[k] = -g_x_0_xxxzz_xxxzz[k] * ab_y + g_x_0_xxxzz_xxxyzz[k];

                g_x_0_xxxyzz_xxyyy[k] = -g_x_0_xxxzz_xxyyy[k] * ab_y + g_x_0_xxxzz_xxyyyy[k];

                g_x_0_xxxyzz_xxyyz[k] = -g_x_0_xxxzz_xxyyz[k] * ab_y + g_x_0_xxxzz_xxyyyz[k];

                g_x_0_xxxyzz_xxyzz[k] = -g_x_0_xxxzz_xxyzz[k] * ab_y + g_x_0_xxxzz_xxyyzz[k];

                g_x_0_xxxyzz_xxzzz[k] = -g_x_0_xxxzz_xxzzz[k] * ab_y + g_x_0_xxxzz_xxyzzz[k];

                g_x_0_xxxyzz_xyyyy[k] = -g_x_0_xxxzz_xyyyy[k] * ab_y + g_x_0_xxxzz_xyyyyy[k];

                g_x_0_xxxyzz_xyyyz[k] = -g_x_0_xxxzz_xyyyz[k] * ab_y + g_x_0_xxxzz_xyyyyz[k];

                g_x_0_xxxyzz_xyyzz[k] = -g_x_0_xxxzz_xyyzz[k] * ab_y + g_x_0_xxxzz_xyyyzz[k];

                g_x_0_xxxyzz_xyzzz[k] = -g_x_0_xxxzz_xyzzz[k] * ab_y + g_x_0_xxxzz_xyyzzz[k];

                g_x_0_xxxyzz_xzzzz[k] = -g_x_0_xxxzz_xzzzz[k] * ab_y + g_x_0_xxxzz_xyzzzz[k];

                g_x_0_xxxyzz_yyyyy[k] = -g_x_0_xxxzz_yyyyy[k] * ab_y + g_x_0_xxxzz_yyyyyy[k];

                g_x_0_xxxyzz_yyyyz[k] = -g_x_0_xxxzz_yyyyz[k] * ab_y + g_x_0_xxxzz_yyyyyz[k];

                g_x_0_xxxyzz_yyyzz[k] = -g_x_0_xxxzz_yyyzz[k] * ab_y + g_x_0_xxxzz_yyyyzz[k];

                g_x_0_xxxyzz_yyzzz[k] = -g_x_0_xxxzz_yyzzz[k] * ab_y + g_x_0_xxxzz_yyyzzz[k];

                g_x_0_xxxyzz_yzzzz[k] = -g_x_0_xxxzz_yzzzz[k] * ab_y + g_x_0_xxxzz_yyzzzz[k];

                g_x_0_xxxyzz_zzzzz[k] = -g_x_0_xxxzz_zzzzz[k] * ab_y + g_x_0_xxxzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzz, g_x_0_xxxzz_zzzzzz, g_x_0_xxxzzz_xxxxx, g_x_0_xxxzzz_xxxxy, g_x_0_xxxzzz_xxxxz, g_x_0_xxxzzz_xxxyy, g_x_0_xxxzzz_xxxyz, g_x_0_xxxzzz_xxxzz, g_x_0_xxxzzz_xxyyy, g_x_0_xxxzzz_xxyyz, g_x_0_xxxzzz_xxyzz, g_x_0_xxxzzz_xxzzz, g_x_0_xxxzzz_xyyyy, g_x_0_xxxzzz_xyyyz, g_x_0_xxxzzz_xyyzz, g_x_0_xxxzzz_xyzzz, g_x_0_xxxzzz_xzzzz, g_x_0_xxxzzz_yyyyy, g_x_0_xxxzzz_yyyyz, g_x_0_xxxzzz_yyyzz, g_x_0_xxxzzz_yyzzz, g_x_0_xxxzzz_yzzzz, g_x_0_xxxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxxxx[k] = -g_x_0_xxxzz_xxxxx[k] * ab_z + g_x_0_xxxzz_xxxxxz[k];

                g_x_0_xxxzzz_xxxxy[k] = -g_x_0_xxxzz_xxxxy[k] * ab_z + g_x_0_xxxzz_xxxxyz[k];

                g_x_0_xxxzzz_xxxxz[k] = -g_x_0_xxxzz_xxxxz[k] * ab_z + g_x_0_xxxzz_xxxxzz[k];

                g_x_0_xxxzzz_xxxyy[k] = -g_x_0_xxxzz_xxxyy[k] * ab_z + g_x_0_xxxzz_xxxyyz[k];

                g_x_0_xxxzzz_xxxyz[k] = -g_x_0_xxxzz_xxxyz[k] * ab_z + g_x_0_xxxzz_xxxyzz[k];

                g_x_0_xxxzzz_xxxzz[k] = -g_x_0_xxxzz_xxxzz[k] * ab_z + g_x_0_xxxzz_xxxzzz[k];

                g_x_0_xxxzzz_xxyyy[k] = -g_x_0_xxxzz_xxyyy[k] * ab_z + g_x_0_xxxzz_xxyyyz[k];

                g_x_0_xxxzzz_xxyyz[k] = -g_x_0_xxxzz_xxyyz[k] * ab_z + g_x_0_xxxzz_xxyyzz[k];

                g_x_0_xxxzzz_xxyzz[k] = -g_x_0_xxxzz_xxyzz[k] * ab_z + g_x_0_xxxzz_xxyzzz[k];

                g_x_0_xxxzzz_xxzzz[k] = -g_x_0_xxxzz_xxzzz[k] * ab_z + g_x_0_xxxzz_xxzzzz[k];

                g_x_0_xxxzzz_xyyyy[k] = -g_x_0_xxxzz_xyyyy[k] * ab_z + g_x_0_xxxzz_xyyyyz[k];

                g_x_0_xxxzzz_xyyyz[k] = -g_x_0_xxxzz_xyyyz[k] * ab_z + g_x_0_xxxzz_xyyyzz[k];

                g_x_0_xxxzzz_xyyzz[k] = -g_x_0_xxxzz_xyyzz[k] * ab_z + g_x_0_xxxzz_xyyzzz[k];

                g_x_0_xxxzzz_xyzzz[k] = -g_x_0_xxxzz_xyzzz[k] * ab_z + g_x_0_xxxzz_xyzzzz[k];

                g_x_0_xxxzzz_xzzzz[k] = -g_x_0_xxxzz_xzzzz[k] * ab_z + g_x_0_xxxzz_xzzzzz[k];

                g_x_0_xxxzzz_yyyyy[k] = -g_x_0_xxxzz_yyyyy[k] * ab_z + g_x_0_xxxzz_yyyyyz[k];

                g_x_0_xxxzzz_yyyyz[k] = -g_x_0_xxxzz_yyyyz[k] * ab_z + g_x_0_xxxzz_yyyyzz[k];

                g_x_0_xxxzzz_yyyzz[k] = -g_x_0_xxxzz_yyyzz[k] * ab_z + g_x_0_xxxzz_yyyzzz[k];

                g_x_0_xxxzzz_yyzzz[k] = -g_x_0_xxxzz_yyzzz[k] * ab_z + g_x_0_xxxzz_yyzzzz[k];

                g_x_0_xxxzzz_yzzzz[k] = -g_x_0_xxxzz_yzzzz[k] * ab_z + g_x_0_xxxzz_yzzzzz[k];

                g_x_0_xxxzzz_zzzzz[k] = -g_x_0_xxxzz_zzzzz[k] * ab_z + g_x_0_xxxzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyy_xxxxx, g_x_0_xxyyy_xxxxxy, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxxyy, g_x_0_xxyyy_xxxxyz, g_x_0_xxyyy_xxxxz, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyyy, g_x_0_xxyyy_xxxyyz, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxyzz, g_x_0_xxyyy_xxxzz, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyyy, g_x_0_xxyyy_xxyyyz, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyyzz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxyzzz, g_x_0_xxyyy_xxzzz, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyyy, g_x_0_xxyyy_xyyyyz, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyyzz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyyzzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xyzzzz, g_x_0_xxyyy_xzzzz, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyyy, g_x_0_xxyyy_yyyyyz, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyyzz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyyzzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yyzzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_yzzzzz, g_x_0_xxyyy_zzzzz, g_x_0_xxyyyy_xxxxx, g_x_0_xxyyyy_xxxxy, g_x_0_xxyyyy_xxxxz, g_x_0_xxyyyy_xxxyy, g_x_0_xxyyyy_xxxyz, g_x_0_xxyyyy_xxxzz, g_x_0_xxyyyy_xxyyy, g_x_0_xxyyyy_xxyyz, g_x_0_xxyyyy_xxyzz, g_x_0_xxyyyy_xxzzz, g_x_0_xxyyyy_xyyyy, g_x_0_xxyyyy_xyyyz, g_x_0_xxyyyy_xyyzz, g_x_0_xxyyyy_xyzzz, g_x_0_xxyyyy_xzzzz, g_x_0_xxyyyy_yyyyy, g_x_0_xxyyyy_yyyyz, g_x_0_xxyyyy_yyyzz, g_x_0_xxyyyy_yyzzz, g_x_0_xxyyyy_yzzzz, g_x_0_xxyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxxxx[k] = -g_x_0_xxyyy_xxxxx[k] * ab_y + g_x_0_xxyyy_xxxxxy[k];

                g_x_0_xxyyyy_xxxxy[k] = -g_x_0_xxyyy_xxxxy[k] * ab_y + g_x_0_xxyyy_xxxxyy[k];

                g_x_0_xxyyyy_xxxxz[k] = -g_x_0_xxyyy_xxxxz[k] * ab_y + g_x_0_xxyyy_xxxxyz[k];

                g_x_0_xxyyyy_xxxyy[k] = -g_x_0_xxyyy_xxxyy[k] * ab_y + g_x_0_xxyyy_xxxyyy[k];

                g_x_0_xxyyyy_xxxyz[k] = -g_x_0_xxyyy_xxxyz[k] * ab_y + g_x_0_xxyyy_xxxyyz[k];

                g_x_0_xxyyyy_xxxzz[k] = -g_x_0_xxyyy_xxxzz[k] * ab_y + g_x_0_xxyyy_xxxyzz[k];

                g_x_0_xxyyyy_xxyyy[k] = -g_x_0_xxyyy_xxyyy[k] * ab_y + g_x_0_xxyyy_xxyyyy[k];

                g_x_0_xxyyyy_xxyyz[k] = -g_x_0_xxyyy_xxyyz[k] * ab_y + g_x_0_xxyyy_xxyyyz[k];

                g_x_0_xxyyyy_xxyzz[k] = -g_x_0_xxyyy_xxyzz[k] * ab_y + g_x_0_xxyyy_xxyyzz[k];

                g_x_0_xxyyyy_xxzzz[k] = -g_x_0_xxyyy_xxzzz[k] * ab_y + g_x_0_xxyyy_xxyzzz[k];

                g_x_0_xxyyyy_xyyyy[k] = -g_x_0_xxyyy_xyyyy[k] * ab_y + g_x_0_xxyyy_xyyyyy[k];

                g_x_0_xxyyyy_xyyyz[k] = -g_x_0_xxyyy_xyyyz[k] * ab_y + g_x_0_xxyyy_xyyyyz[k];

                g_x_0_xxyyyy_xyyzz[k] = -g_x_0_xxyyy_xyyzz[k] * ab_y + g_x_0_xxyyy_xyyyzz[k];

                g_x_0_xxyyyy_xyzzz[k] = -g_x_0_xxyyy_xyzzz[k] * ab_y + g_x_0_xxyyy_xyyzzz[k];

                g_x_0_xxyyyy_xzzzz[k] = -g_x_0_xxyyy_xzzzz[k] * ab_y + g_x_0_xxyyy_xyzzzz[k];

                g_x_0_xxyyyy_yyyyy[k] = -g_x_0_xxyyy_yyyyy[k] * ab_y + g_x_0_xxyyy_yyyyyy[k];

                g_x_0_xxyyyy_yyyyz[k] = -g_x_0_xxyyy_yyyyz[k] * ab_y + g_x_0_xxyyy_yyyyyz[k];

                g_x_0_xxyyyy_yyyzz[k] = -g_x_0_xxyyy_yyyzz[k] * ab_y + g_x_0_xxyyy_yyyyzz[k];

                g_x_0_xxyyyy_yyzzz[k] = -g_x_0_xxyyy_yyzzz[k] * ab_y + g_x_0_xxyyy_yyyzzz[k];

                g_x_0_xxyyyy_yzzzz[k] = -g_x_0_xxyyy_yzzzz[k] * ab_y + g_x_0_xxyyy_yyzzzz[k];

                g_x_0_xxyyyy_zzzzz[k] = -g_x_0_xxyyy_zzzzz[k] * ab_y + g_x_0_xxyyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyyz_xxxxx, g_x_0_xxyyyz_xxxxy, g_x_0_xxyyyz_xxxxz, g_x_0_xxyyyz_xxxyy, g_x_0_xxyyyz_xxxyz, g_x_0_xxyyyz_xxxzz, g_x_0_xxyyyz_xxyyy, g_x_0_xxyyyz_xxyyz, g_x_0_xxyyyz_xxyzz, g_x_0_xxyyyz_xxzzz, g_x_0_xxyyyz_xyyyy, g_x_0_xxyyyz_xyyyz, g_x_0_xxyyyz_xyyzz, g_x_0_xxyyyz_xyzzz, g_x_0_xxyyyz_xzzzz, g_x_0_xxyyyz_yyyyy, g_x_0_xxyyyz_yyyyz, g_x_0_xxyyyz_yyyzz, g_x_0_xxyyyz_yyzzz, g_x_0_xxyyyz_yzzzz, g_x_0_xxyyyz_zzzzz, g_x_0_xxyyz_xxxxx, g_x_0_xxyyz_xxxxxy, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxxyy, g_x_0_xxyyz_xxxxyz, g_x_0_xxyyz_xxxxz, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyyy, g_x_0_xxyyz_xxxyyz, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxyzz, g_x_0_xxyyz_xxxzz, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyyy, g_x_0_xxyyz_xxyyyz, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyyzz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxyzzz, g_x_0_xxyyz_xxzzz, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyyy, g_x_0_xxyyz_xyyyyz, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyyzz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyyzzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xyzzzz, g_x_0_xxyyz_xzzzz, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyyy, g_x_0_xxyyz_yyyyyz, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyyzz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyyzzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yyzzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_yzzzzz, g_x_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxxxx[k] = -g_x_0_xxyyz_xxxxx[k] * ab_y + g_x_0_xxyyz_xxxxxy[k];

                g_x_0_xxyyyz_xxxxy[k] = -g_x_0_xxyyz_xxxxy[k] * ab_y + g_x_0_xxyyz_xxxxyy[k];

                g_x_0_xxyyyz_xxxxz[k] = -g_x_0_xxyyz_xxxxz[k] * ab_y + g_x_0_xxyyz_xxxxyz[k];

                g_x_0_xxyyyz_xxxyy[k] = -g_x_0_xxyyz_xxxyy[k] * ab_y + g_x_0_xxyyz_xxxyyy[k];

                g_x_0_xxyyyz_xxxyz[k] = -g_x_0_xxyyz_xxxyz[k] * ab_y + g_x_0_xxyyz_xxxyyz[k];

                g_x_0_xxyyyz_xxxzz[k] = -g_x_0_xxyyz_xxxzz[k] * ab_y + g_x_0_xxyyz_xxxyzz[k];

                g_x_0_xxyyyz_xxyyy[k] = -g_x_0_xxyyz_xxyyy[k] * ab_y + g_x_0_xxyyz_xxyyyy[k];

                g_x_0_xxyyyz_xxyyz[k] = -g_x_0_xxyyz_xxyyz[k] * ab_y + g_x_0_xxyyz_xxyyyz[k];

                g_x_0_xxyyyz_xxyzz[k] = -g_x_0_xxyyz_xxyzz[k] * ab_y + g_x_0_xxyyz_xxyyzz[k];

                g_x_0_xxyyyz_xxzzz[k] = -g_x_0_xxyyz_xxzzz[k] * ab_y + g_x_0_xxyyz_xxyzzz[k];

                g_x_0_xxyyyz_xyyyy[k] = -g_x_0_xxyyz_xyyyy[k] * ab_y + g_x_0_xxyyz_xyyyyy[k];

                g_x_0_xxyyyz_xyyyz[k] = -g_x_0_xxyyz_xyyyz[k] * ab_y + g_x_0_xxyyz_xyyyyz[k];

                g_x_0_xxyyyz_xyyzz[k] = -g_x_0_xxyyz_xyyzz[k] * ab_y + g_x_0_xxyyz_xyyyzz[k];

                g_x_0_xxyyyz_xyzzz[k] = -g_x_0_xxyyz_xyzzz[k] * ab_y + g_x_0_xxyyz_xyyzzz[k];

                g_x_0_xxyyyz_xzzzz[k] = -g_x_0_xxyyz_xzzzz[k] * ab_y + g_x_0_xxyyz_xyzzzz[k];

                g_x_0_xxyyyz_yyyyy[k] = -g_x_0_xxyyz_yyyyy[k] * ab_y + g_x_0_xxyyz_yyyyyy[k];

                g_x_0_xxyyyz_yyyyz[k] = -g_x_0_xxyyz_yyyyz[k] * ab_y + g_x_0_xxyyz_yyyyyz[k];

                g_x_0_xxyyyz_yyyzz[k] = -g_x_0_xxyyz_yyyzz[k] * ab_y + g_x_0_xxyyz_yyyyzz[k];

                g_x_0_xxyyyz_yyzzz[k] = -g_x_0_xxyyz_yyzzz[k] * ab_y + g_x_0_xxyyz_yyyzzz[k];

                g_x_0_xxyyyz_yzzzz[k] = -g_x_0_xxyyz_yzzzz[k] * ab_y + g_x_0_xxyyz_yyzzzz[k];

                g_x_0_xxyyyz_zzzzz[k] = -g_x_0_xxyyz_zzzzz[k] * ab_y + g_x_0_xxyyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyzz_xxxxx, g_x_0_xxyyzz_xxxxy, g_x_0_xxyyzz_xxxxz, g_x_0_xxyyzz_xxxyy, g_x_0_xxyyzz_xxxyz, g_x_0_xxyyzz_xxxzz, g_x_0_xxyyzz_xxyyy, g_x_0_xxyyzz_xxyyz, g_x_0_xxyyzz_xxyzz, g_x_0_xxyyzz_xxzzz, g_x_0_xxyyzz_xyyyy, g_x_0_xxyyzz_xyyyz, g_x_0_xxyyzz_xyyzz, g_x_0_xxyyzz_xyzzz, g_x_0_xxyyzz_xzzzz, g_x_0_xxyyzz_yyyyy, g_x_0_xxyyzz_yyyyz, g_x_0_xxyyzz_yyyzz, g_x_0_xxyyzz_yyzzz, g_x_0_xxyyzz_yzzzz, g_x_0_xxyyzz_zzzzz, g_x_0_xxyzz_xxxxx, g_x_0_xxyzz_xxxxxy, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxxyy, g_x_0_xxyzz_xxxxyz, g_x_0_xxyzz_xxxxz, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyyy, g_x_0_xxyzz_xxxyyz, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxyzz, g_x_0_xxyzz_xxxzz, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyyy, g_x_0_xxyzz_xxyyyz, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyyzz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxyzzz, g_x_0_xxyzz_xxzzz, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyyy, g_x_0_xxyzz_xyyyyz, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyyzz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyyzzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xyzzzz, g_x_0_xxyzz_xzzzz, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyyy, g_x_0_xxyzz_yyyyyz, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyyzz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyyzzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yyzzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_yzzzzz, g_x_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxxxx[k] = -g_x_0_xxyzz_xxxxx[k] * ab_y + g_x_0_xxyzz_xxxxxy[k];

                g_x_0_xxyyzz_xxxxy[k] = -g_x_0_xxyzz_xxxxy[k] * ab_y + g_x_0_xxyzz_xxxxyy[k];

                g_x_0_xxyyzz_xxxxz[k] = -g_x_0_xxyzz_xxxxz[k] * ab_y + g_x_0_xxyzz_xxxxyz[k];

                g_x_0_xxyyzz_xxxyy[k] = -g_x_0_xxyzz_xxxyy[k] * ab_y + g_x_0_xxyzz_xxxyyy[k];

                g_x_0_xxyyzz_xxxyz[k] = -g_x_0_xxyzz_xxxyz[k] * ab_y + g_x_0_xxyzz_xxxyyz[k];

                g_x_0_xxyyzz_xxxzz[k] = -g_x_0_xxyzz_xxxzz[k] * ab_y + g_x_0_xxyzz_xxxyzz[k];

                g_x_0_xxyyzz_xxyyy[k] = -g_x_0_xxyzz_xxyyy[k] * ab_y + g_x_0_xxyzz_xxyyyy[k];

                g_x_0_xxyyzz_xxyyz[k] = -g_x_0_xxyzz_xxyyz[k] * ab_y + g_x_0_xxyzz_xxyyyz[k];

                g_x_0_xxyyzz_xxyzz[k] = -g_x_0_xxyzz_xxyzz[k] * ab_y + g_x_0_xxyzz_xxyyzz[k];

                g_x_0_xxyyzz_xxzzz[k] = -g_x_0_xxyzz_xxzzz[k] * ab_y + g_x_0_xxyzz_xxyzzz[k];

                g_x_0_xxyyzz_xyyyy[k] = -g_x_0_xxyzz_xyyyy[k] * ab_y + g_x_0_xxyzz_xyyyyy[k];

                g_x_0_xxyyzz_xyyyz[k] = -g_x_0_xxyzz_xyyyz[k] * ab_y + g_x_0_xxyzz_xyyyyz[k];

                g_x_0_xxyyzz_xyyzz[k] = -g_x_0_xxyzz_xyyzz[k] * ab_y + g_x_0_xxyzz_xyyyzz[k];

                g_x_0_xxyyzz_xyzzz[k] = -g_x_0_xxyzz_xyzzz[k] * ab_y + g_x_0_xxyzz_xyyzzz[k];

                g_x_0_xxyyzz_xzzzz[k] = -g_x_0_xxyzz_xzzzz[k] * ab_y + g_x_0_xxyzz_xyzzzz[k];

                g_x_0_xxyyzz_yyyyy[k] = -g_x_0_xxyzz_yyyyy[k] * ab_y + g_x_0_xxyzz_yyyyyy[k];

                g_x_0_xxyyzz_yyyyz[k] = -g_x_0_xxyzz_yyyyz[k] * ab_y + g_x_0_xxyzz_yyyyyz[k];

                g_x_0_xxyyzz_yyyzz[k] = -g_x_0_xxyzz_yyyzz[k] * ab_y + g_x_0_xxyzz_yyyyzz[k];

                g_x_0_xxyyzz_yyzzz[k] = -g_x_0_xxyzz_yyzzz[k] * ab_y + g_x_0_xxyzz_yyyzzz[k];

                g_x_0_xxyyzz_yzzzz[k] = -g_x_0_xxyzz_yzzzz[k] * ab_y + g_x_0_xxyzz_yyzzzz[k];

                g_x_0_xxyyzz_zzzzz[k] = -g_x_0_xxyzz_zzzzz[k] * ab_y + g_x_0_xxyzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzzz_xxxxx, g_x_0_xxyzzz_xxxxy, g_x_0_xxyzzz_xxxxz, g_x_0_xxyzzz_xxxyy, g_x_0_xxyzzz_xxxyz, g_x_0_xxyzzz_xxxzz, g_x_0_xxyzzz_xxyyy, g_x_0_xxyzzz_xxyyz, g_x_0_xxyzzz_xxyzz, g_x_0_xxyzzz_xxzzz, g_x_0_xxyzzz_xyyyy, g_x_0_xxyzzz_xyyyz, g_x_0_xxyzzz_xyyzz, g_x_0_xxyzzz_xyzzz, g_x_0_xxyzzz_xzzzz, g_x_0_xxyzzz_yyyyy, g_x_0_xxyzzz_yyyyz, g_x_0_xxyzzz_yyyzz, g_x_0_xxyzzz_yyzzz, g_x_0_xxyzzz_yzzzz, g_x_0_xxyzzz_zzzzz, g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxxxx[k] = -g_x_0_xxzzz_xxxxx[k] * ab_y + g_x_0_xxzzz_xxxxxy[k];

                g_x_0_xxyzzz_xxxxy[k] = -g_x_0_xxzzz_xxxxy[k] * ab_y + g_x_0_xxzzz_xxxxyy[k];

                g_x_0_xxyzzz_xxxxz[k] = -g_x_0_xxzzz_xxxxz[k] * ab_y + g_x_0_xxzzz_xxxxyz[k];

                g_x_0_xxyzzz_xxxyy[k] = -g_x_0_xxzzz_xxxyy[k] * ab_y + g_x_0_xxzzz_xxxyyy[k];

                g_x_0_xxyzzz_xxxyz[k] = -g_x_0_xxzzz_xxxyz[k] * ab_y + g_x_0_xxzzz_xxxyyz[k];

                g_x_0_xxyzzz_xxxzz[k] = -g_x_0_xxzzz_xxxzz[k] * ab_y + g_x_0_xxzzz_xxxyzz[k];

                g_x_0_xxyzzz_xxyyy[k] = -g_x_0_xxzzz_xxyyy[k] * ab_y + g_x_0_xxzzz_xxyyyy[k];

                g_x_0_xxyzzz_xxyyz[k] = -g_x_0_xxzzz_xxyyz[k] * ab_y + g_x_0_xxzzz_xxyyyz[k];

                g_x_0_xxyzzz_xxyzz[k] = -g_x_0_xxzzz_xxyzz[k] * ab_y + g_x_0_xxzzz_xxyyzz[k];

                g_x_0_xxyzzz_xxzzz[k] = -g_x_0_xxzzz_xxzzz[k] * ab_y + g_x_0_xxzzz_xxyzzz[k];

                g_x_0_xxyzzz_xyyyy[k] = -g_x_0_xxzzz_xyyyy[k] * ab_y + g_x_0_xxzzz_xyyyyy[k];

                g_x_0_xxyzzz_xyyyz[k] = -g_x_0_xxzzz_xyyyz[k] * ab_y + g_x_0_xxzzz_xyyyyz[k];

                g_x_0_xxyzzz_xyyzz[k] = -g_x_0_xxzzz_xyyzz[k] * ab_y + g_x_0_xxzzz_xyyyzz[k];

                g_x_0_xxyzzz_xyzzz[k] = -g_x_0_xxzzz_xyzzz[k] * ab_y + g_x_0_xxzzz_xyyzzz[k];

                g_x_0_xxyzzz_xzzzz[k] = -g_x_0_xxzzz_xzzzz[k] * ab_y + g_x_0_xxzzz_xyzzzz[k];

                g_x_0_xxyzzz_yyyyy[k] = -g_x_0_xxzzz_yyyyy[k] * ab_y + g_x_0_xxzzz_yyyyyy[k];

                g_x_0_xxyzzz_yyyyz[k] = -g_x_0_xxzzz_yyyyz[k] * ab_y + g_x_0_xxzzz_yyyyyz[k];

                g_x_0_xxyzzz_yyyzz[k] = -g_x_0_xxzzz_yyyzz[k] * ab_y + g_x_0_xxzzz_yyyyzz[k];

                g_x_0_xxyzzz_yyzzz[k] = -g_x_0_xxzzz_yyzzz[k] * ab_y + g_x_0_xxzzz_yyyzzz[k];

                g_x_0_xxyzzz_yzzzz[k] = -g_x_0_xxzzz_yzzzz[k] * ab_y + g_x_0_xxzzz_yyzzzz[k];

                g_x_0_xxyzzz_zzzzz[k] = -g_x_0_xxzzz_zzzzz[k] * ab_y + g_x_0_xxzzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzz, g_x_0_xxzzz_zzzzzz, g_x_0_xxzzzz_xxxxx, g_x_0_xxzzzz_xxxxy, g_x_0_xxzzzz_xxxxz, g_x_0_xxzzzz_xxxyy, g_x_0_xxzzzz_xxxyz, g_x_0_xxzzzz_xxxzz, g_x_0_xxzzzz_xxyyy, g_x_0_xxzzzz_xxyyz, g_x_0_xxzzzz_xxyzz, g_x_0_xxzzzz_xxzzz, g_x_0_xxzzzz_xyyyy, g_x_0_xxzzzz_xyyyz, g_x_0_xxzzzz_xyyzz, g_x_0_xxzzzz_xyzzz, g_x_0_xxzzzz_xzzzz, g_x_0_xxzzzz_yyyyy, g_x_0_xxzzzz_yyyyz, g_x_0_xxzzzz_yyyzz, g_x_0_xxzzzz_yyzzz, g_x_0_xxzzzz_yzzzz, g_x_0_xxzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxxxx[k] = -g_x_0_xxzzz_xxxxx[k] * ab_z + g_x_0_xxzzz_xxxxxz[k];

                g_x_0_xxzzzz_xxxxy[k] = -g_x_0_xxzzz_xxxxy[k] * ab_z + g_x_0_xxzzz_xxxxyz[k];

                g_x_0_xxzzzz_xxxxz[k] = -g_x_0_xxzzz_xxxxz[k] * ab_z + g_x_0_xxzzz_xxxxzz[k];

                g_x_0_xxzzzz_xxxyy[k] = -g_x_0_xxzzz_xxxyy[k] * ab_z + g_x_0_xxzzz_xxxyyz[k];

                g_x_0_xxzzzz_xxxyz[k] = -g_x_0_xxzzz_xxxyz[k] * ab_z + g_x_0_xxzzz_xxxyzz[k];

                g_x_0_xxzzzz_xxxzz[k] = -g_x_0_xxzzz_xxxzz[k] * ab_z + g_x_0_xxzzz_xxxzzz[k];

                g_x_0_xxzzzz_xxyyy[k] = -g_x_0_xxzzz_xxyyy[k] * ab_z + g_x_0_xxzzz_xxyyyz[k];

                g_x_0_xxzzzz_xxyyz[k] = -g_x_0_xxzzz_xxyyz[k] * ab_z + g_x_0_xxzzz_xxyyzz[k];

                g_x_0_xxzzzz_xxyzz[k] = -g_x_0_xxzzz_xxyzz[k] * ab_z + g_x_0_xxzzz_xxyzzz[k];

                g_x_0_xxzzzz_xxzzz[k] = -g_x_0_xxzzz_xxzzz[k] * ab_z + g_x_0_xxzzz_xxzzzz[k];

                g_x_0_xxzzzz_xyyyy[k] = -g_x_0_xxzzz_xyyyy[k] * ab_z + g_x_0_xxzzz_xyyyyz[k];

                g_x_0_xxzzzz_xyyyz[k] = -g_x_0_xxzzz_xyyyz[k] * ab_z + g_x_0_xxzzz_xyyyzz[k];

                g_x_0_xxzzzz_xyyzz[k] = -g_x_0_xxzzz_xyyzz[k] * ab_z + g_x_0_xxzzz_xyyzzz[k];

                g_x_0_xxzzzz_xyzzz[k] = -g_x_0_xxzzz_xyzzz[k] * ab_z + g_x_0_xxzzz_xyzzzz[k];

                g_x_0_xxzzzz_xzzzz[k] = -g_x_0_xxzzz_xzzzz[k] * ab_z + g_x_0_xxzzz_xzzzzz[k];

                g_x_0_xxzzzz_yyyyy[k] = -g_x_0_xxzzz_yyyyy[k] * ab_z + g_x_0_xxzzz_yyyyyz[k];

                g_x_0_xxzzzz_yyyyz[k] = -g_x_0_xxzzz_yyyyz[k] * ab_z + g_x_0_xxzzz_yyyyzz[k];

                g_x_0_xxzzzz_yyyzz[k] = -g_x_0_xxzzz_yyyzz[k] * ab_z + g_x_0_xxzzz_yyyzzz[k];

                g_x_0_xxzzzz_yyzzz[k] = -g_x_0_xxzzz_yyzzz[k] * ab_z + g_x_0_xxzzz_yyzzzz[k];

                g_x_0_xxzzzz_yzzzz[k] = -g_x_0_xxzzz_yzzzz[k] * ab_z + g_x_0_xxzzz_yzzzzz[k];

                g_x_0_xxzzzz_zzzzz[k] = -g_x_0_xxzzz_zzzzz[k] * ab_z + g_x_0_xxzzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyy_xxxxx, g_x_0_xyyyy_xxxxxy, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxxyy, g_x_0_xyyyy_xxxxyz, g_x_0_xyyyy_xxxxz, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyyy, g_x_0_xyyyy_xxxyyz, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxyzz, g_x_0_xyyyy_xxxzz, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyyy, g_x_0_xyyyy_xxyyyz, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyyzz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxyzzz, g_x_0_xyyyy_xxzzz, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyyy, g_x_0_xyyyy_xyyyyz, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyyzz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyyzzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xyzzzz, g_x_0_xyyyy_xzzzz, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyyy, g_x_0_xyyyy_yyyyyz, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyyzz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyyzzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yyzzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_yzzzzz, g_x_0_xyyyy_zzzzz, g_x_0_xyyyyy_xxxxx, g_x_0_xyyyyy_xxxxy, g_x_0_xyyyyy_xxxxz, g_x_0_xyyyyy_xxxyy, g_x_0_xyyyyy_xxxyz, g_x_0_xyyyyy_xxxzz, g_x_0_xyyyyy_xxyyy, g_x_0_xyyyyy_xxyyz, g_x_0_xyyyyy_xxyzz, g_x_0_xyyyyy_xxzzz, g_x_0_xyyyyy_xyyyy, g_x_0_xyyyyy_xyyyz, g_x_0_xyyyyy_xyyzz, g_x_0_xyyyyy_xyzzz, g_x_0_xyyyyy_xzzzz, g_x_0_xyyyyy_yyyyy, g_x_0_xyyyyy_yyyyz, g_x_0_xyyyyy_yyyzz, g_x_0_xyyyyy_yyzzz, g_x_0_xyyyyy_yzzzz, g_x_0_xyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxxxx[k] = -g_x_0_xyyyy_xxxxx[k] * ab_y + g_x_0_xyyyy_xxxxxy[k];

                g_x_0_xyyyyy_xxxxy[k] = -g_x_0_xyyyy_xxxxy[k] * ab_y + g_x_0_xyyyy_xxxxyy[k];

                g_x_0_xyyyyy_xxxxz[k] = -g_x_0_xyyyy_xxxxz[k] * ab_y + g_x_0_xyyyy_xxxxyz[k];

                g_x_0_xyyyyy_xxxyy[k] = -g_x_0_xyyyy_xxxyy[k] * ab_y + g_x_0_xyyyy_xxxyyy[k];

                g_x_0_xyyyyy_xxxyz[k] = -g_x_0_xyyyy_xxxyz[k] * ab_y + g_x_0_xyyyy_xxxyyz[k];

                g_x_0_xyyyyy_xxxzz[k] = -g_x_0_xyyyy_xxxzz[k] * ab_y + g_x_0_xyyyy_xxxyzz[k];

                g_x_0_xyyyyy_xxyyy[k] = -g_x_0_xyyyy_xxyyy[k] * ab_y + g_x_0_xyyyy_xxyyyy[k];

                g_x_0_xyyyyy_xxyyz[k] = -g_x_0_xyyyy_xxyyz[k] * ab_y + g_x_0_xyyyy_xxyyyz[k];

                g_x_0_xyyyyy_xxyzz[k] = -g_x_0_xyyyy_xxyzz[k] * ab_y + g_x_0_xyyyy_xxyyzz[k];

                g_x_0_xyyyyy_xxzzz[k] = -g_x_0_xyyyy_xxzzz[k] * ab_y + g_x_0_xyyyy_xxyzzz[k];

                g_x_0_xyyyyy_xyyyy[k] = -g_x_0_xyyyy_xyyyy[k] * ab_y + g_x_0_xyyyy_xyyyyy[k];

                g_x_0_xyyyyy_xyyyz[k] = -g_x_0_xyyyy_xyyyz[k] * ab_y + g_x_0_xyyyy_xyyyyz[k];

                g_x_0_xyyyyy_xyyzz[k] = -g_x_0_xyyyy_xyyzz[k] * ab_y + g_x_0_xyyyy_xyyyzz[k];

                g_x_0_xyyyyy_xyzzz[k] = -g_x_0_xyyyy_xyzzz[k] * ab_y + g_x_0_xyyyy_xyyzzz[k];

                g_x_0_xyyyyy_xzzzz[k] = -g_x_0_xyyyy_xzzzz[k] * ab_y + g_x_0_xyyyy_xyzzzz[k];

                g_x_0_xyyyyy_yyyyy[k] = -g_x_0_xyyyy_yyyyy[k] * ab_y + g_x_0_xyyyy_yyyyyy[k];

                g_x_0_xyyyyy_yyyyz[k] = -g_x_0_xyyyy_yyyyz[k] * ab_y + g_x_0_xyyyy_yyyyyz[k];

                g_x_0_xyyyyy_yyyzz[k] = -g_x_0_xyyyy_yyyzz[k] * ab_y + g_x_0_xyyyy_yyyyzz[k];

                g_x_0_xyyyyy_yyzzz[k] = -g_x_0_xyyyy_yyzzz[k] * ab_y + g_x_0_xyyyy_yyyzzz[k];

                g_x_0_xyyyyy_yzzzz[k] = -g_x_0_xyyyy_yzzzz[k] * ab_y + g_x_0_xyyyy_yyzzzz[k];

                g_x_0_xyyyyy_zzzzz[k] = -g_x_0_xyyyy_zzzzz[k] * ab_y + g_x_0_xyyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyyz_xxxxx, g_x_0_xyyyyz_xxxxy, g_x_0_xyyyyz_xxxxz, g_x_0_xyyyyz_xxxyy, g_x_0_xyyyyz_xxxyz, g_x_0_xyyyyz_xxxzz, g_x_0_xyyyyz_xxyyy, g_x_0_xyyyyz_xxyyz, g_x_0_xyyyyz_xxyzz, g_x_0_xyyyyz_xxzzz, g_x_0_xyyyyz_xyyyy, g_x_0_xyyyyz_xyyyz, g_x_0_xyyyyz_xyyzz, g_x_0_xyyyyz_xyzzz, g_x_0_xyyyyz_xzzzz, g_x_0_xyyyyz_yyyyy, g_x_0_xyyyyz_yyyyz, g_x_0_xyyyyz_yyyzz, g_x_0_xyyyyz_yyzzz, g_x_0_xyyyyz_yzzzz, g_x_0_xyyyyz_zzzzz, g_x_0_xyyyz_xxxxx, g_x_0_xyyyz_xxxxxy, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxxyy, g_x_0_xyyyz_xxxxyz, g_x_0_xyyyz_xxxxz, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyyy, g_x_0_xyyyz_xxxyyz, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxyzz, g_x_0_xyyyz_xxxzz, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyyy, g_x_0_xyyyz_xxyyyz, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyyzz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxyzzz, g_x_0_xyyyz_xxzzz, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyyy, g_x_0_xyyyz_xyyyyz, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyyzz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyyzzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xyzzzz, g_x_0_xyyyz_xzzzz, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyyy, g_x_0_xyyyz_yyyyyz, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyyzz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyyzzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yyzzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_yzzzzz, g_x_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxxxx[k] = -g_x_0_xyyyz_xxxxx[k] * ab_y + g_x_0_xyyyz_xxxxxy[k];

                g_x_0_xyyyyz_xxxxy[k] = -g_x_0_xyyyz_xxxxy[k] * ab_y + g_x_0_xyyyz_xxxxyy[k];

                g_x_0_xyyyyz_xxxxz[k] = -g_x_0_xyyyz_xxxxz[k] * ab_y + g_x_0_xyyyz_xxxxyz[k];

                g_x_0_xyyyyz_xxxyy[k] = -g_x_0_xyyyz_xxxyy[k] * ab_y + g_x_0_xyyyz_xxxyyy[k];

                g_x_0_xyyyyz_xxxyz[k] = -g_x_0_xyyyz_xxxyz[k] * ab_y + g_x_0_xyyyz_xxxyyz[k];

                g_x_0_xyyyyz_xxxzz[k] = -g_x_0_xyyyz_xxxzz[k] * ab_y + g_x_0_xyyyz_xxxyzz[k];

                g_x_0_xyyyyz_xxyyy[k] = -g_x_0_xyyyz_xxyyy[k] * ab_y + g_x_0_xyyyz_xxyyyy[k];

                g_x_0_xyyyyz_xxyyz[k] = -g_x_0_xyyyz_xxyyz[k] * ab_y + g_x_0_xyyyz_xxyyyz[k];

                g_x_0_xyyyyz_xxyzz[k] = -g_x_0_xyyyz_xxyzz[k] * ab_y + g_x_0_xyyyz_xxyyzz[k];

                g_x_0_xyyyyz_xxzzz[k] = -g_x_0_xyyyz_xxzzz[k] * ab_y + g_x_0_xyyyz_xxyzzz[k];

                g_x_0_xyyyyz_xyyyy[k] = -g_x_0_xyyyz_xyyyy[k] * ab_y + g_x_0_xyyyz_xyyyyy[k];

                g_x_0_xyyyyz_xyyyz[k] = -g_x_0_xyyyz_xyyyz[k] * ab_y + g_x_0_xyyyz_xyyyyz[k];

                g_x_0_xyyyyz_xyyzz[k] = -g_x_0_xyyyz_xyyzz[k] * ab_y + g_x_0_xyyyz_xyyyzz[k];

                g_x_0_xyyyyz_xyzzz[k] = -g_x_0_xyyyz_xyzzz[k] * ab_y + g_x_0_xyyyz_xyyzzz[k];

                g_x_0_xyyyyz_xzzzz[k] = -g_x_0_xyyyz_xzzzz[k] * ab_y + g_x_0_xyyyz_xyzzzz[k];

                g_x_0_xyyyyz_yyyyy[k] = -g_x_0_xyyyz_yyyyy[k] * ab_y + g_x_0_xyyyz_yyyyyy[k];

                g_x_0_xyyyyz_yyyyz[k] = -g_x_0_xyyyz_yyyyz[k] * ab_y + g_x_0_xyyyz_yyyyyz[k];

                g_x_0_xyyyyz_yyyzz[k] = -g_x_0_xyyyz_yyyzz[k] * ab_y + g_x_0_xyyyz_yyyyzz[k];

                g_x_0_xyyyyz_yyzzz[k] = -g_x_0_xyyyz_yyzzz[k] * ab_y + g_x_0_xyyyz_yyyzzz[k];

                g_x_0_xyyyyz_yzzzz[k] = -g_x_0_xyyyz_yzzzz[k] * ab_y + g_x_0_xyyyz_yyzzzz[k];

                g_x_0_xyyyyz_zzzzz[k] = -g_x_0_xyyyz_zzzzz[k] * ab_y + g_x_0_xyyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyzz_xxxxx, g_x_0_xyyyzz_xxxxy, g_x_0_xyyyzz_xxxxz, g_x_0_xyyyzz_xxxyy, g_x_0_xyyyzz_xxxyz, g_x_0_xyyyzz_xxxzz, g_x_0_xyyyzz_xxyyy, g_x_0_xyyyzz_xxyyz, g_x_0_xyyyzz_xxyzz, g_x_0_xyyyzz_xxzzz, g_x_0_xyyyzz_xyyyy, g_x_0_xyyyzz_xyyyz, g_x_0_xyyyzz_xyyzz, g_x_0_xyyyzz_xyzzz, g_x_0_xyyyzz_xzzzz, g_x_0_xyyyzz_yyyyy, g_x_0_xyyyzz_yyyyz, g_x_0_xyyyzz_yyyzz, g_x_0_xyyyzz_yyzzz, g_x_0_xyyyzz_yzzzz, g_x_0_xyyyzz_zzzzz, g_x_0_xyyzz_xxxxx, g_x_0_xyyzz_xxxxxy, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxxyy, g_x_0_xyyzz_xxxxyz, g_x_0_xyyzz_xxxxz, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyyy, g_x_0_xyyzz_xxxyyz, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxyzz, g_x_0_xyyzz_xxxzz, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyyy, g_x_0_xyyzz_xxyyyz, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyyzz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxyzzz, g_x_0_xyyzz_xxzzz, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyyy, g_x_0_xyyzz_xyyyyz, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyyzz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyyzzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xyzzzz, g_x_0_xyyzz_xzzzz, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyyy, g_x_0_xyyzz_yyyyyz, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyyzz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyyzzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yyzzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_yzzzzz, g_x_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxxxx[k] = -g_x_0_xyyzz_xxxxx[k] * ab_y + g_x_0_xyyzz_xxxxxy[k];

                g_x_0_xyyyzz_xxxxy[k] = -g_x_0_xyyzz_xxxxy[k] * ab_y + g_x_0_xyyzz_xxxxyy[k];

                g_x_0_xyyyzz_xxxxz[k] = -g_x_0_xyyzz_xxxxz[k] * ab_y + g_x_0_xyyzz_xxxxyz[k];

                g_x_0_xyyyzz_xxxyy[k] = -g_x_0_xyyzz_xxxyy[k] * ab_y + g_x_0_xyyzz_xxxyyy[k];

                g_x_0_xyyyzz_xxxyz[k] = -g_x_0_xyyzz_xxxyz[k] * ab_y + g_x_0_xyyzz_xxxyyz[k];

                g_x_0_xyyyzz_xxxzz[k] = -g_x_0_xyyzz_xxxzz[k] * ab_y + g_x_0_xyyzz_xxxyzz[k];

                g_x_0_xyyyzz_xxyyy[k] = -g_x_0_xyyzz_xxyyy[k] * ab_y + g_x_0_xyyzz_xxyyyy[k];

                g_x_0_xyyyzz_xxyyz[k] = -g_x_0_xyyzz_xxyyz[k] * ab_y + g_x_0_xyyzz_xxyyyz[k];

                g_x_0_xyyyzz_xxyzz[k] = -g_x_0_xyyzz_xxyzz[k] * ab_y + g_x_0_xyyzz_xxyyzz[k];

                g_x_0_xyyyzz_xxzzz[k] = -g_x_0_xyyzz_xxzzz[k] * ab_y + g_x_0_xyyzz_xxyzzz[k];

                g_x_0_xyyyzz_xyyyy[k] = -g_x_0_xyyzz_xyyyy[k] * ab_y + g_x_0_xyyzz_xyyyyy[k];

                g_x_0_xyyyzz_xyyyz[k] = -g_x_0_xyyzz_xyyyz[k] * ab_y + g_x_0_xyyzz_xyyyyz[k];

                g_x_0_xyyyzz_xyyzz[k] = -g_x_0_xyyzz_xyyzz[k] * ab_y + g_x_0_xyyzz_xyyyzz[k];

                g_x_0_xyyyzz_xyzzz[k] = -g_x_0_xyyzz_xyzzz[k] * ab_y + g_x_0_xyyzz_xyyzzz[k];

                g_x_0_xyyyzz_xzzzz[k] = -g_x_0_xyyzz_xzzzz[k] * ab_y + g_x_0_xyyzz_xyzzzz[k];

                g_x_0_xyyyzz_yyyyy[k] = -g_x_0_xyyzz_yyyyy[k] * ab_y + g_x_0_xyyzz_yyyyyy[k];

                g_x_0_xyyyzz_yyyyz[k] = -g_x_0_xyyzz_yyyyz[k] * ab_y + g_x_0_xyyzz_yyyyyz[k];

                g_x_0_xyyyzz_yyyzz[k] = -g_x_0_xyyzz_yyyzz[k] * ab_y + g_x_0_xyyzz_yyyyzz[k];

                g_x_0_xyyyzz_yyzzz[k] = -g_x_0_xyyzz_yyzzz[k] * ab_y + g_x_0_xyyzz_yyyzzz[k];

                g_x_0_xyyyzz_yzzzz[k] = -g_x_0_xyyzz_yzzzz[k] * ab_y + g_x_0_xyyzz_yyzzzz[k];

                g_x_0_xyyyzz_zzzzz[k] = -g_x_0_xyyzz_zzzzz[k] * ab_y + g_x_0_xyyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzzz_xxxxx, g_x_0_xyyzzz_xxxxy, g_x_0_xyyzzz_xxxxz, g_x_0_xyyzzz_xxxyy, g_x_0_xyyzzz_xxxyz, g_x_0_xyyzzz_xxxzz, g_x_0_xyyzzz_xxyyy, g_x_0_xyyzzz_xxyyz, g_x_0_xyyzzz_xxyzz, g_x_0_xyyzzz_xxzzz, g_x_0_xyyzzz_xyyyy, g_x_0_xyyzzz_xyyyz, g_x_0_xyyzzz_xyyzz, g_x_0_xyyzzz_xyzzz, g_x_0_xyyzzz_xzzzz, g_x_0_xyyzzz_yyyyy, g_x_0_xyyzzz_yyyyz, g_x_0_xyyzzz_yyyzz, g_x_0_xyyzzz_yyzzz, g_x_0_xyyzzz_yzzzz, g_x_0_xyyzzz_zzzzz, g_x_0_xyzzz_xxxxx, g_x_0_xyzzz_xxxxxy, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxxyy, g_x_0_xyzzz_xxxxyz, g_x_0_xyzzz_xxxxz, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyyy, g_x_0_xyzzz_xxxyyz, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxyzz, g_x_0_xyzzz_xxxzz, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyyy, g_x_0_xyzzz_xxyyyz, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyyzz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxyzzz, g_x_0_xyzzz_xxzzz, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyyy, g_x_0_xyzzz_xyyyyz, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyyzz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyyzzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xyzzzz, g_x_0_xyzzz_xzzzz, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyyy, g_x_0_xyzzz_yyyyyz, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyyzz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyyzzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yyzzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_yzzzzz, g_x_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxxxx[k] = -g_x_0_xyzzz_xxxxx[k] * ab_y + g_x_0_xyzzz_xxxxxy[k];

                g_x_0_xyyzzz_xxxxy[k] = -g_x_0_xyzzz_xxxxy[k] * ab_y + g_x_0_xyzzz_xxxxyy[k];

                g_x_0_xyyzzz_xxxxz[k] = -g_x_0_xyzzz_xxxxz[k] * ab_y + g_x_0_xyzzz_xxxxyz[k];

                g_x_0_xyyzzz_xxxyy[k] = -g_x_0_xyzzz_xxxyy[k] * ab_y + g_x_0_xyzzz_xxxyyy[k];

                g_x_0_xyyzzz_xxxyz[k] = -g_x_0_xyzzz_xxxyz[k] * ab_y + g_x_0_xyzzz_xxxyyz[k];

                g_x_0_xyyzzz_xxxzz[k] = -g_x_0_xyzzz_xxxzz[k] * ab_y + g_x_0_xyzzz_xxxyzz[k];

                g_x_0_xyyzzz_xxyyy[k] = -g_x_0_xyzzz_xxyyy[k] * ab_y + g_x_0_xyzzz_xxyyyy[k];

                g_x_0_xyyzzz_xxyyz[k] = -g_x_0_xyzzz_xxyyz[k] * ab_y + g_x_0_xyzzz_xxyyyz[k];

                g_x_0_xyyzzz_xxyzz[k] = -g_x_0_xyzzz_xxyzz[k] * ab_y + g_x_0_xyzzz_xxyyzz[k];

                g_x_0_xyyzzz_xxzzz[k] = -g_x_0_xyzzz_xxzzz[k] * ab_y + g_x_0_xyzzz_xxyzzz[k];

                g_x_0_xyyzzz_xyyyy[k] = -g_x_0_xyzzz_xyyyy[k] * ab_y + g_x_0_xyzzz_xyyyyy[k];

                g_x_0_xyyzzz_xyyyz[k] = -g_x_0_xyzzz_xyyyz[k] * ab_y + g_x_0_xyzzz_xyyyyz[k];

                g_x_0_xyyzzz_xyyzz[k] = -g_x_0_xyzzz_xyyzz[k] * ab_y + g_x_0_xyzzz_xyyyzz[k];

                g_x_0_xyyzzz_xyzzz[k] = -g_x_0_xyzzz_xyzzz[k] * ab_y + g_x_0_xyzzz_xyyzzz[k];

                g_x_0_xyyzzz_xzzzz[k] = -g_x_0_xyzzz_xzzzz[k] * ab_y + g_x_0_xyzzz_xyzzzz[k];

                g_x_0_xyyzzz_yyyyy[k] = -g_x_0_xyzzz_yyyyy[k] * ab_y + g_x_0_xyzzz_yyyyyy[k];

                g_x_0_xyyzzz_yyyyz[k] = -g_x_0_xyzzz_yyyyz[k] * ab_y + g_x_0_xyzzz_yyyyyz[k];

                g_x_0_xyyzzz_yyyzz[k] = -g_x_0_xyzzz_yyyzz[k] * ab_y + g_x_0_xyzzz_yyyyzz[k];

                g_x_0_xyyzzz_yyzzz[k] = -g_x_0_xyzzz_yyzzz[k] * ab_y + g_x_0_xyzzz_yyyzzz[k];

                g_x_0_xyyzzz_yzzzz[k] = -g_x_0_xyzzz_yzzzz[k] * ab_y + g_x_0_xyzzz_yyzzzz[k];

                g_x_0_xyyzzz_zzzzz[k] = -g_x_0_xyzzz_zzzzz[k] * ab_y + g_x_0_xyzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzzz_xxxxx, g_x_0_xyzzzz_xxxxy, g_x_0_xyzzzz_xxxxz, g_x_0_xyzzzz_xxxyy, g_x_0_xyzzzz_xxxyz, g_x_0_xyzzzz_xxxzz, g_x_0_xyzzzz_xxyyy, g_x_0_xyzzzz_xxyyz, g_x_0_xyzzzz_xxyzz, g_x_0_xyzzzz_xxzzz, g_x_0_xyzzzz_xyyyy, g_x_0_xyzzzz_xyyyz, g_x_0_xyzzzz_xyyzz, g_x_0_xyzzzz_xyzzz, g_x_0_xyzzzz_xzzzz, g_x_0_xyzzzz_yyyyy, g_x_0_xyzzzz_yyyyz, g_x_0_xyzzzz_yyyzz, g_x_0_xyzzzz_yyzzz, g_x_0_xyzzzz_yzzzz, g_x_0_xyzzzz_zzzzz, g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxxxx[k] = -g_x_0_xzzzz_xxxxx[k] * ab_y + g_x_0_xzzzz_xxxxxy[k];

                g_x_0_xyzzzz_xxxxy[k] = -g_x_0_xzzzz_xxxxy[k] * ab_y + g_x_0_xzzzz_xxxxyy[k];

                g_x_0_xyzzzz_xxxxz[k] = -g_x_0_xzzzz_xxxxz[k] * ab_y + g_x_0_xzzzz_xxxxyz[k];

                g_x_0_xyzzzz_xxxyy[k] = -g_x_0_xzzzz_xxxyy[k] * ab_y + g_x_0_xzzzz_xxxyyy[k];

                g_x_0_xyzzzz_xxxyz[k] = -g_x_0_xzzzz_xxxyz[k] * ab_y + g_x_0_xzzzz_xxxyyz[k];

                g_x_0_xyzzzz_xxxzz[k] = -g_x_0_xzzzz_xxxzz[k] * ab_y + g_x_0_xzzzz_xxxyzz[k];

                g_x_0_xyzzzz_xxyyy[k] = -g_x_0_xzzzz_xxyyy[k] * ab_y + g_x_0_xzzzz_xxyyyy[k];

                g_x_0_xyzzzz_xxyyz[k] = -g_x_0_xzzzz_xxyyz[k] * ab_y + g_x_0_xzzzz_xxyyyz[k];

                g_x_0_xyzzzz_xxyzz[k] = -g_x_0_xzzzz_xxyzz[k] * ab_y + g_x_0_xzzzz_xxyyzz[k];

                g_x_0_xyzzzz_xxzzz[k] = -g_x_0_xzzzz_xxzzz[k] * ab_y + g_x_0_xzzzz_xxyzzz[k];

                g_x_0_xyzzzz_xyyyy[k] = -g_x_0_xzzzz_xyyyy[k] * ab_y + g_x_0_xzzzz_xyyyyy[k];

                g_x_0_xyzzzz_xyyyz[k] = -g_x_0_xzzzz_xyyyz[k] * ab_y + g_x_0_xzzzz_xyyyyz[k];

                g_x_0_xyzzzz_xyyzz[k] = -g_x_0_xzzzz_xyyzz[k] * ab_y + g_x_0_xzzzz_xyyyzz[k];

                g_x_0_xyzzzz_xyzzz[k] = -g_x_0_xzzzz_xyzzz[k] * ab_y + g_x_0_xzzzz_xyyzzz[k];

                g_x_0_xyzzzz_xzzzz[k] = -g_x_0_xzzzz_xzzzz[k] * ab_y + g_x_0_xzzzz_xyzzzz[k];

                g_x_0_xyzzzz_yyyyy[k] = -g_x_0_xzzzz_yyyyy[k] * ab_y + g_x_0_xzzzz_yyyyyy[k];

                g_x_0_xyzzzz_yyyyz[k] = -g_x_0_xzzzz_yyyyz[k] * ab_y + g_x_0_xzzzz_yyyyyz[k];

                g_x_0_xyzzzz_yyyzz[k] = -g_x_0_xzzzz_yyyzz[k] * ab_y + g_x_0_xzzzz_yyyyzz[k];

                g_x_0_xyzzzz_yyzzz[k] = -g_x_0_xzzzz_yyzzz[k] * ab_y + g_x_0_xzzzz_yyyzzz[k];

                g_x_0_xyzzzz_yzzzz[k] = -g_x_0_xzzzz_yzzzz[k] * ab_y + g_x_0_xzzzz_yyzzzz[k];

                g_x_0_xyzzzz_zzzzz[k] = -g_x_0_xzzzz_zzzzz[k] * ab_y + g_x_0_xzzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 434 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzz, g_x_0_xzzzz_zzzzzz, g_x_0_xzzzzz_xxxxx, g_x_0_xzzzzz_xxxxy, g_x_0_xzzzzz_xxxxz, g_x_0_xzzzzz_xxxyy, g_x_0_xzzzzz_xxxyz, g_x_0_xzzzzz_xxxzz, g_x_0_xzzzzz_xxyyy, g_x_0_xzzzzz_xxyyz, g_x_0_xzzzzz_xxyzz, g_x_0_xzzzzz_xxzzz, g_x_0_xzzzzz_xyyyy, g_x_0_xzzzzz_xyyyz, g_x_0_xzzzzz_xyyzz, g_x_0_xzzzzz_xyzzz, g_x_0_xzzzzz_xzzzz, g_x_0_xzzzzz_yyyyy, g_x_0_xzzzzz_yyyyz, g_x_0_xzzzzz_yyyzz, g_x_0_xzzzzz_yyzzz, g_x_0_xzzzzz_yzzzz, g_x_0_xzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxxxx[k] = -g_x_0_xzzzz_xxxxx[k] * ab_z + g_x_0_xzzzz_xxxxxz[k];

                g_x_0_xzzzzz_xxxxy[k] = -g_x_0_xzzzz_xxxxy[k] * ab_z + g_x_0_xzzzz_xxxxyz[k];

                g_x_0_xzzzzz_xxxxz[k] = -g_x_0_xzzzz_xxxxz[k] * ab_z + g_x_0_xzzzz_xxxxzz[k];

                g_x_0_xzzzzz_xxxyy[k] = -g_x_0_xzzzz_xxxyy[k] * ab_z + g_x_0_xzzzz_xxxyyz[k];

                g_x_0_xzzzzz_xxxyz[k] = -g_x_0_xzzzz_xxxyz[k] * ab_z + g_x_0_xzzzz_xxxyzz[k];

                g_x_0_xzzzzz_xxxzz[k] = -g_x_0_xzzzz_xxxzz[k] * ab_z + g_x_0_xzzzz_xxxzzz[k];

                g_x_0_xzzzzz_xxyyy[k] = -g_x_0_xzzzz_xxyyy[k] * ab_z + g_x_0_xzzzz_xxyyyz[k];

                g_x_0_xzzzzz_xxyyz[k] = -g_x_0_xzzzz_xxyyz[k] * ab_z + g_x_0_xzzzz_xxyyzz[k];

                g_x_0_xzzzzz_xxyzz[k] = -g_x_0_xzzzz_xxyzz[k] * ab_z + g_x_0_xzzzz_xxyzzz[k];

                g_x_0_xzzzzz_xxzzz[k] = -g_x_0_xzzzz_xxzzz[k] * ab_z + g_x_0_xzzzz_xxzzzz[k];

                g_x_0_xzzzzz_xyyyy[k] = -g_x_0_xzzzz_xyyyy[k] * ab_z + g_x_0_xzzzz_xyyyyz[k];

                g_x_0_xzzzzz_xyyyz[k] = -g_x_0_xzzzz_xyyyz[k] * ab_z + g_x_0_xzzzz_xyyyzz[k];

                g_x_0_xzzzzz_xyyzz[k] = -g_x_0_xzzzz_xyyzz[k] * ab_z + g_x_0_xzzzz_xyyzzz[k];

                g_x_0_xzzzzz_xyzzz[k] = -g_x_0_xzzzz_xyzzz[k] * ab_z + g_x_0_xzzzz_xyzzzz[k];

                g_x_0_xzzzzz_xzzzz[k] = -g_x_0_xzzzz_xzzzz[k] * ab_z + g_x_0_xzzzz_xzzzzz[k];

                g_x_0_xzzzzz_yyyyy[k] = -g_x_0_xzzzz_yyyyy[k] * ab_z + g_x_0_xzzzz_yyyyyz[k];

                g_x_0_xzzzzz_yyyyz[k] = -g_x_0_xzzzz_yyyyz[k] * ab_z + g_x_0_xzzzz_yyyyzz[k];

                g_x_0_xzzzzz_yyyzz[k] = -g_x_0_xzzzz_yyyzz[k] * ab_z + g_x_0_xzzzz_yyyzzz[k];

                g_x_0_xzzzzz_yyzzz[k] = -g_x_0_xzzzz_yyzzz[k] * ab_z + g_x_0_xzzzz_yyzzzz[k];

                g_x_0_xzzzzz_yzzzz[k] = -g_x_0_xzzzz_yzzzz[k] * ab_z + g_x_0_xzzzz_yzzzzz[k];

                g_x_0_xzzzzz_zzzzz[k] = -g_x_0_xzzzz_zzzzz[k] * ab_z + g_x_0_xzzzz_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 441 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 446 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 447 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 448 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 450 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 452 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 453 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 457 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 459 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_xxxxx, g_x_0_yyyyy_xxxxxy, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxxyy, g_x_0_yyyyy_xxxxyz, g_x_0_yyyyy_xxxxz, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyyy, g_x_0_yyyyy_xxxyyz, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxyzz, g_x_0_yyyyy_xxxzz, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyyy, g_x_0_yyyyy_xxyyyz, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyyzz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxyzzz, g_x_0_yyyyy_xxzzz, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyyy, g_x_0_yyyyy_xyyyyz, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyyzz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyyzzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xyzzzz, g_x_0_yyyyy_xzzzz, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyyy, g_x_0_yyyyy_yyyyyz, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyyzz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyyzzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yyzzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_yzzzzz, g_x_0_yyyyy_zzzzz, g_x_0_yyyyyy_xxxxx, g_x_0_yyyyyy_xxxxy, g_x_0_yyyyyy_xxxxz, g_x_0_yyyyyy_xxxyy, g_x_0_yyyyyy_xxxyz, g_x_0_yyyyyy_xxxzz, g_x_0_yyyyyy_xxyyy, g_x_0_yyyyyy_xxyyz, g_x_0_yyyyyy_xxyzz, g_x_0_yyyyyy_xxzzz, g_x_0_yyyyyy_xyyyy, g_x_0_yyyyyy_xyyyz, g_x_0_yyyyyy_xyyzz, g_x_0_yyyyyy_xyzzz, g_x_0_yyyyyy_xzzzz, g_x_0_yyyyyy_yyyyy, g_x_0_yyyyyy_yyyyz, g_x_0_yyyyyy_yyyzz, g_x_0_yyyyyy_yyzzz, g_x_0_yyyyyy_yzzzz, g_x_0_yyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxxxx[k] = -g_x_0_yyyyy_xxxxx[k] * ab_y + g_x_0_yyyyy_xxxxxy[k];

                g_x_0_yyyyyy_xxxxy[k] = -g_x_0_yyyyy_xxxxy[k] * ab_y + g_x_0_yyyyy_xxxxyy[k];

                g_x_0_yyyyyy_xxxxz[k] = -g_x_0_yyyyy_xxxxz[k] * ab_y + g_x_0_yyyyy_xxxxyz[k];

                g_x_0_yyyyyy_xxxyy[k] = -g_x_0_yyyyy_xxxyy[k] * ab_y + g_x_0_yyyyy_xxxyyy[k];

                g_x_0_yyyyyy_xxxyz[k] = -g_x_0_yyyyy_xxxyz[k] * ab_y + g_x_0_yyyyy_xxxyyz[k];

                g_x_0_yyyyyy_xxxzz[k] = -g_x_0_yyyyy_xxxzz[k] * ab_y + g_x_0_yyyyy_xxxyzz[k];

                g_x_0_yyyyyy_xxyyy[k] = -g_x_0_yyyyy_xxyyy[k] * ab_y + g_x_0_yyyyy_xxyyyy[k];

                g_x_0_yyyyyy_xxyyz[k] = -g_x_0_yyyyy_xxyyz[k] * ab_y + g_x_0_yyyyy_xxyyyz[k];

                g_x_0_yyyyyy_xxyzz[k] = -g_x_0_yyyyy_xxyzz[k] * ab_y + g_x_0_yyyyy_xxyyzz[k];

                g_x_0_yyyyyy_xxzzz[k] = -g_x_0_yyyyy_xxzzz[k] * ab_y + g_x_0_yyyyy_xxyzzz[k];

                g_x_0_yyyyyy_xyyyy[k] = -g_x_0_yyyyy_xyyyy[k] * ab_y + g_x_0_yyyyy_xyyyyy[k];

                g_x_0_yyyyyy_xyyyz[k] = -g_x_0_yyyyy_xyyyz[k] * ab_y + g_x_0_yyyyy_xyyyyz[k];

                g_x_0_yyyyyy_xyyzz[k] = -g_x_0_yyyyy_xyyzz[k] * ab_y + g_x_0_yyyyy_xyyyzz[k];

                g_x_0_yyyyyy_xyzzz[k] = -g_x_0_yyyyy_xyzzz[k] * ab_y + g_x_0_yyyyy_xyyzzz[k];

                g_x_0_yyyyyy_xzzzz[k] = -g_x_0_yyyyy_xzzzz[k] * ab_y + g_x_0_yyyyy_xyzzzz[k];

                g_x_0_yyyyyy_yyyyy[k] = -g_x_0_yyyyy_yyyyy[k] * ab_y + g_x_0_yyyyy_yyyyyy[k];

                g_x_0_yyyyyy_yyyyz[k] = -g_x_0_yyyyy_yyyyz[k] * ab_y + g_x_0_yyyyy_yyyyyz[k];

                g_x_0_yyyyyy_yyyzz[k] = -g_x_0_yyyyy_yyyzz[k] * ab_y + g_x_0_yyyyy_yyyyzz[k];

                g_x_0_yyyyyy_yyzzz[k] = -g_x_0_yyyyy_yyzzz[k] * ab_y + g_x_0_yyyyy_yyyzzz[k];

                g_x_0_yyyyyy_yzzzz[k] = -g_x_0_yyyyy_yzzzz[k] * ab_y + g_x_0_yyyyy_yyzzzz[k];

                g_x_0_yyyyyy_zzzzz[k] = -g_x_0_yyyyy_zzzzz[k] * ab_y + g_x_0_yyyyy_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 462 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 467 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 468 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 470 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 473 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 475 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 476 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 477 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 478 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 481 * ccomps * dcomps);

            auto g_x_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyyz_xxxxx, g_x_0_yyyyyz_xxxxy, g_x_0_yyyyyz_xxxxz, g_x_0_yyyyyz_xxxyy, g_x_0_yyyyyz_xxxyz, g_x_0_yyyyyz_xxxzz, g_x_0_yyyyyz_xxyyy, g_x_0_yyyyyz_xxyyz, g_x_0_yyyyyz_xxyzz, g_x_0_yyyyyz_xxzzz, g_x_0_yyyyyz_xyyyy, g_x_0_yyyyyz_xyyyz, g_x_0_yyyyyz_xyyzz, g_x_0_yyyyyz_xyzzz, g_x_0_yyyyyz_xzzzz, g_x_0_yyyyyz_yyyyy, g_x_0_yyyyyz_yyyyz, g_x_0_yyyyyz_yyyzz, g_x_0_yyyyyz_yyzzz, g_x_0_yyyyyz_yzzzz, g_x_0_yyyyyz_zzzzz, g_x_0_yyyyz_xxxxx, g_x_0_yyyyz_xxxxxy, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxxyy, g_x_0_yyyyz_xxxxyz, g_x_0_yyyyz_xxxxz, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyyy, g_x_0_yyyyz_xxxyyz, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxyzz, g_x_0_yyyyz_xxxzz, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyyy, g_x_0_yyyyz_xxyyyz, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyyzz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxyzzz, g_x_0_yyyyz_xxzzz, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyyy, g_x_0_yyyyz_xyyyyz, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyyzz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyyzzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xyzzzz, g_x_0_yyyyz_xzzzz, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyyy, g_x_0_yyyyz_yyyyyz, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyyzz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyyzzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yyzzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_yzzzzz, g_x_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxxxx[k] = -g_x_0_yyyyz_xxxxx[k] * ab_y + g_x_0_yyyyz_xxxxxy[k];

                g_x_0_yyyyyz_xxxxy[k] = -g_x_0_yyyyz_xxxxy[k] * ab_y + g_x_0_yyyyz_xxxxyy[k];

                g_x_0_yyyyyz_xxxxz[k] = -g_x_0_yyyyz_xxxxz[k] * ab_y + g_x_0_yyyyz_xxxxyz[k];

                g_x_0_yyyyyz_xxxyy[k] = -g_x_0_yyyyz_xxxyy[k] * ab_y + g_x_0_yyyyz_xxxyyy[k];

                g_x_0_yyyyyz_xxxyz[k] = -g_x_0_yyyyz_xxxyz[k] * ab_y + g_x_0_yyyyz_xxxyyz[k];

                g_x_0_yyyyyz_xxxzz[k] = -g_x_0_yyyyz_xxxzz[k] * ab_y + g_x_0_yyyyz_xxxyzz[k];

                g_x_0_yyyyyz_xxyyy[k] = -g_x_0_yyyyz_xxyyy[k] * ab_y + g_x_0_yyyyz_xxyyyy[k];

                g_x_0_yyyyyz_xxyyz[k] = -g_x_0_yyyyz_xxyyz[k] * ab_y + g_x_0_yyyyz_xxyyyz[k];

                g_x_0_yyyyyz_xxyzz[k] = -g_x_0_yyyyz_xxyzz[k] * ab_y + g_x_0_yyyyz_xxyyzz[k];

                g_x_0_yyyyyz_xxzzz[k] = -g_x_0_yyyyz_xxzzz[k] * ab_y + g_x_0_yyyyz_xxyzzz[k];

                g_x_0_yyyyyz_xyyyy[k] = -g_x_0_yyyyz_xyyyy[k] * ab_y + g_x_0_yyyyz_xyyyyy[k];

                g_x_0_yyyyyz_xyyyz[k] = -g_x_0_yyyyz_xyyyz[k] * ab_y + g_x_0_yyyyz_xyyyyz[k];

                g_x_0_yyyyyz_xyyzz[k] = -g_x_0_yyyyz_xyyzz[k] * ab_y + g_x_0_yyyyz_xyyyzz[k];

                g_x_0_yyyyyz_xyzzz[k] = -g_x_0_yyyyz_xyzzz[k] * ab_y + g_x_0_yyyyz_xyyzzz[k];

                g_x_0_yyyyyz_xzzzz[k] = -g_x_0_yyyyz_xzzzz[k] * ab_y + g_x_0_yyyyz_xyzzzz[k];

                g_x_0_yyyyyz_yyyyy[k] = -g_x_0_yyyyz_yyyyy[k] * ab_y + g_x_0_yyyyz_yyyyyy[k];

                g_x_0_yyyyyz_yyyyz[k] = -g_x_0_yyyyz_yyyyz[k] * ab_y + g_x_0_yyyyz_yyyyyz[k];

                g_x_0_yyyyyz_yyyzz[k] = -g_x_0_yyyyz_yyyzz[k] * ab_y + g_x_0_yyyyz_yyyyzz[k];

                g_x_0_yyyyyz_yyzzz[k] = -g_x_0_yyyyz_yyzzz[k] * ab_y + g_x_0_yyyyz_yyyzzz[k];

                g_x_0_yyyyyz_yzzzz[k] = -g_x_0_yyyyz_yzzzz[k] * ab_y + g_x_0_yyyyz_yyzzzz[k];

                g_x_0_yyyyyz_zzzzz[k] = -g_x_0_yyyyz_zzzzz[k] * ab_y + g_x_0_yyyyz_yzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 485 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 488 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 490 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 495 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 496 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyzz_xxxxx, g_x_0_yyyyzz_xxxxy, g_x_0_yyyyzz_xxxxz, g_x_0_yyyyzz_xxxyy, g_x_0_yyyyzz_xxxyz, g_x_0_yyyyzz_xxxzz, g_x_0_yyyyzz_xxyyy, g_x_0_yyyyzz_xxyyz, g_x_0_yyyyzz_xxyzz, g_x_0_yyyyzz_xxzzz, g_x_0_yyyyzz_xyyyy, g_x_0_yyyyzz_xyyyz, g_x_0_yyyyzz_xyyzz, g_x_0_yyyyzz_xyzzz, g_x_0_yyyyzz_xzzzz, g_x_0_yyyyzz_yyyyy, g_x_0_yyyyzz_yyyyz, g_x_0_yyyyzz_yyyzz, g_x_0_yyyyzz_yyzzz, g_x_0_yyyyzz_yzzzz, g_x_0_yyyyzz_zzzzz, g_x_0_yyyzz_xxxxx, g_x_0_yyyzz_xxxxxy, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxxyy, g_x_0_yyyzz_xxxxyz, g_x_0_yyyzz_xxxxz, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyyy, g_x_0_yyyzz_xxxyyz, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxyzz, g_x_0_yyyzz_xxxzz, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyyy, g_x_0_yyyzz_xxyyyz, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyyzz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxyzzz, g_x_0_yyyzz_xxzzz, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyyy, g_x_0_yyyzz_xyyyyz, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyyzz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyyzzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xyzzzz, g_x_0_yyyzz_xzzzz, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyyy, g_x_0_yyyzz_yyyyyz, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyyzz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyyzzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yyzzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_yzzzzz, g_x_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxxxx[k] = -g_x_0_yyyzz_xxxxx[k] * ab_y + g_x_0_yyyzz_xxxxxy[k];

                g_x_0_yyyyzz_xxxxy[k] = -g_x_0_yyyzz_xxxxy[k] * ab_y + g_x_0_yyyzz_xxxxyy[k];

                g_x_0_yyyyzz_xxxxz[k] = -g_x_0_yyyzz_xxxxz[k] * ab_y + g_x_0_yyyzz_xxxxyz[k];

                g_x_0_yyyyzz_xxxyy[k] = -g_x_0_yyyzz_xxxyy[k] * ab_y + g_x_0_yyyzz_xxxyyy[k];

                g_x_0_yyyyzz_xxxyz[k] = -g_x_0_yyyzz_xxxyz[k] * ab_y + g_x_0_yyyzz_xxxyyz[k];

                g_x_0_yyyyzz_xxxzz[k] = -g_x_0_yyyzz_xxxzz[k] * ab_y + g_x_0_yyyzz_xxxyzz[k];

                g_x_0_yyyyzz_xxyyy[k] = -g_x_0_yyyzz_xxyyy[k] * ab_y + g_x_0_yyyzz_xxyyyy[k];

                g_x_0_yyyyzz_xxyyz[k] = -g_x_0_yyyzz_xxyyz[k] * ab_y + g_x_0_yyyzz_xxyyyz[k];

                g_x_0_yyyyzz_xxyzz[k] = -g_x_0_yyyzz_xxyzz[k] * ab_y + g_x_0_yyyzz_xxyyzz[k];

                g_x_0_yyyyzz_xxzzz[k] = -g_x_0_yyyzz_xxzzz[k] * ab_y + g_x_0_yyyzz_xxyzzz[k];

                g_x_0_yyyyzz_xyyyy[k] = -g_x_0_yyyzz_xyyyy[k] * ab_y + g_x_0_yyyzz_xyyyyy[k];

                g_x_0_yyyyzz_xyyyz[k] = -g_x_0_yyyzz_xyyyz[k] * ab_y + g_x_0_yyyzz_xyyyyz[k];

                g_x_0_yyyyzz_xyyzz[k] = -g_x_0_yyyzz_xyyzz[k] * ab_y + g_x_0_yyyzz_xyyyzz[k];

                g_x_0_yyyyzz_xyzzz[k] = -g_x_0_yyyzz_xyzzz[k] * ab_y + g_x_0_yyyzz_xyyzzz[k];

                g_x_0_yyyyzz_xzzzz[k] = -g_x_0_yyyzz_xzzzz[k] * ab_y + g_x_0_yyyzz_xyzzzz[k];

                g_x_0_yyyyzz_yyyyy[k] = -g_x_0_yyyzz_yyyyy[k] * ab_y + g_x_0_yyyzz_yyyyyy[k];

                g_x_0_yyyyzz_yyyyz[k] = -g_x_0_yyyzz_yyyyz[k] * ab_y + g_x_0_yyyzz_yyyyyz[k];

                g_x_0_yyyyzz_yyyzz[k] = -g_x_0_yyyzz_yyyzz[k] * ab_y + g_x_0_yyyzz_yyyyzz[k];

                g_x_0_yyyyzz_yyzzz[k] = -g_x_0_yyyzz_yyzzz[k] * ab_y + g_x_0_yyyzz_yyyzzz[k];

                g_x_0_yyyyzz_yzzzz[k] = -g_x_0_yyyzz_yzzzz[k] * ab_y + g_x_0_yyyzz_yyzzzz[k];

                g_x_0_yyyyzz_zzzzz[k] = -g_x_0_yyyzz_zzzzz[k] * ab_y + g_x_0_yyyzz_yzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 504 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 506 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 509 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 513 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 518 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzzz_xxxxx, g_x_0_yyyzzz_xxxxy, g_x_0_yyyzzz_xxxxz, g_x_0_yyyzzz_xxxyy, g_x_0_yyyzzz_xxxyz, g_x_0_yyyzzz_xxxzz, g_x_0_yyyzzz_xxyyy, g_x_0_yyyzzz_xxyyz, g_x_0_yyyzzz_xxyzz, g_x_0_yyyzzz_xxzzz, g_x_0_yyyzzz_xyyyy, g_x_0_yyyzzz_xyyyz, g_x_0_yyyzzz_xyyzz, g_x_0_yyyzzz_xyzzz, g_x_0_yyyzzz_xzzzz, g_x_0_yyyzzz_yyyyy, g_x_0_yyyzzz_yyyyz, g_x_0_yyyzzz_yyyzz, g_x_0_yyyzzz_yyzzz, g_x_0_yyyzzz_yzzzz, g_x_0_yyyzzz_zzzzz, g_x_0_yyzzz_xxxxx, g_x_0_yyzzz_xxxxxy, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxxyy, g_x_0_yyzzz_xxxxyz, g_x_0_yyzzz_xxxxz, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyyy, g_x_0_yyzzz_xxxyyz, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxyzz, g_x_0_yyzzz_xxxzz, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyyy, g_x_0_yyzzz_xxyyyz, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyyzz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxyzzz, g_x_0_yyzzz_xxzzz, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyyy, g_x_0_yyzzz_xyyyyz, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyyzz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyyzzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xyzzzz, g_x_0_yyzzz_xzzzz, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyyy, g_x_0_yyzzz_yyyyyz, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyyzz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyyzzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yyzzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_yzzzzz, g_x_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxxxx[k] = -g_x_0_yyzzz_xxxxx[k] * ab_y + g_x_0_yyzzz_xxxxxy[k];

                g_x_0_yyyzzz_xxxxy[k] = -g_x_0_yyzzz_xxxxy[k] * ab_y + g_x_0_yyzzz_xxxxyy[k];

                g_x_0_yyyzzz_xxxxz[k] = -g_x_0_yyzzz_xxxxz[k] * ab_y + g_x_0_yyzzz_xxxxyz[k];

                g_x_0_yyyzzz_xxxyy[k] = -g_x_0_yyzzz_xxxyy[k] * ab_y + g_x_0_yyzzz_xxxyyy[k];

                g_x_0_yyyzzz_xxxyz[k] = -g_x_0_yyzzz_xxxyz[k] * ab_y + g_x_0_yyzzz_xxxyyz[k];

                g_x_0_yyyzzz_xxxzz[k] = -g_x_0_yyzzz_xxxzz[k] * ab_y + g_x_0_yyzzz_xxxyzz[k];

                g_x_0_yyyzzz_xxyyy[k] = -g_x_0_yyzzz_xxyyy[k] * ab_y + g_x_0_yyzzz_xxyyyy[k];

                g_x_0_yyyzzz_xxyyz[k] = -g_x_0_yyzzz_xxyyz[k] * ab_y + g_x_0_yyzzz_xxyyyz[k];

                g_x_0_yyyzzz_xxyzz[k] = -g_x_0_yyzzz_xxyzz[k] * ab_y + g_x_0_yyzzz_xxyyzz[k];

                g_x_0_yyyzzz_xxzzz[k] = -g_x_0_yyzzz_xxzzz[k] * ab_y + g_x_0_yyzzz_xxyzzz[k];

                g_x_0_yyyzzz_xyyyy[k] = -g_x_0_yyzzz_xyyyy[k] * ab_y + g_x_0_yyzzz_xyyyyy[k];

                g_x_0_yyyzzz_xyyyz[k] = -g_x_0_yyzzz_xyyyz[k] * ab_y + g_x_0_yyzzz_xyyyyz[k];

                g_x_0_yyyzzz_xyyzz[k] = -g_x_0_yyzzz_xyyzz[k] * ab_y + g_x_0_yyzzz_xyyyzz[k];

                g_x_0_yyyzzz_xyzzz[k] = -g_x_0_yyzzz_xyzzz[k] * ab_y + g_x_0_yyzzz_xyyzzz[k];

                g_x_0_yyyzzz_xzzzz[k] = -g_x_0_yyzzz_xzzzz[k] * ab_y + g_x_0_yyzzz_xyzzzz[k];

                g_x_0_yyyzzz_yyyyy[k] = -g_x_0_yyzzz_yyyyy[k] * ab_y + g_x_0_yyzzz_yyyyyy[k];

                g_x_0_yyyzzz_yyyyz[k] = -g_x_0_yyzzz_yyyyz[k] * ab_y + g_x_0_yyzzz_yyyyyz[k];

                g_x_0_yyyzzz_yyyzz[k] = -g_x_0_yyzzz_yyyzz[k] * ab_y + g_x_0_yyzzz_yyyyzz[k];

                g_x_0_yyyzzz_yyzzz[k] = -g_x_0_yyzzz_yyzzz[k] * ab_y + g_x_0_yyzzz_yyyzzz[k];

                g_x_0_yyyzzz_yzzzz[k] = -g_x_0_yyzzz_yzzzz[k] * ab_y + g_x_0_yyzzz_yyzzzz[k];

                g_x_0_yyyzzz_zzzzz[k] = -g_x_0_yyzzz_zzzzz[k] * ab_y + g_x_0_yyzzz_yzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 531 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 532 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 534 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 537 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 539 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 540 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 541 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 542 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 543 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 544 * ccomps * dcomps);

            auto g_x_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzzz_xxxxx, g_x_0_yyzzzz_xxxxy, g_x_0_yyzzzz_xxxxz, g_x_0_yyzzzz_xxxyy, g_x_0_yyzzzz_xxxyz, g_x_0_yyzzzz_xxxzz, g_x_0_yyzzzz_xxyyy, g_x_0_yyzzzz_xxyyz, g_x_0_yyzzzz_xxyzz, g_x_0_yyzzzz_xxzzz, g_x_0_yyzzzz_xyyyy, g_x_0_yyzzzz_xyyyz, g_x_0_yyzzzz_xyyzz, g_x_0_yyzzzz_xyzzz, g_x_0_yyzzzz_xzzzz, g_x_0_yyzzzz_yyyyy, g_x_0_yyzzzz_yyyyz, g_x_0_yyzzzz_yyyzz, g_x_0_yyzzzz_yyzzz, g_x_0_yyzzzz_yzzzz, g_x_0_yyzzzz_zzzzz, g_x_0_yzzzz_xxxxx, g_x_0_yzzzz_xxxxxy, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxxyy, g_x_0_yzzzz_xxxxyz, g_x_0_yzzzz_xxxxz, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyyy, g_x_0_yzzzz_xxxyyz, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxyzz, g_x_0_yzzzz_xxxzz, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyyy, g_x_0_yzzzz_xxyyyz, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyyzz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxyzzz, g_x_0_yzzzz_xxzzz, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyyy, g_x_0_yzzzz_xyyyyz, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyyzz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyyzzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xyzzzz, g_x_0_yzzzz_xzzzz, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyyy, g_x_0_yzzzz_yyyyyz, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyyzz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyyzzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yyzzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_yzzzzz, g_x_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxxxx[k] = -g_x_0_yzzzz_xxxxx[k] * ab_y + g_x_0_yzzzz_xxxxxy[k];

                g_x_0_yyzzzz_xxxxy[k] = -g_x_0_yzzzz_xxxxy[k] * ab_y + g_x_0_yzzzz_xxxxyy[k];

                g_x_0_yyzzzz_xxxxz[k] = -g_x_0_yzzzz_xxxxz[k] * ab_y + g_x_0_yzzzz_xxxxyz[k];

                g_x_0_yyzzzz_xxxyy[k] = -g_x_0_yzzzz_xxxyy[k] * ab_y + g_x_0_yzzzz_xxxyyy[k];

                g_x_0_yyzzzz_xxxyz[k] = -g_x_0_yzzzz_xxxyz[k] * ab_y + g_x_0_yzzzz_xxxyyz[k];

                g_x_0_yyzzzz_xxxzz[k] = -g_x_0_yzzzz_xxxzz[k] * ab_y + g_x_0_yzzzz_xxxyzz[k];

                g_x_0_yyzzzz_xxyyy[k] = -g_x_0_yzzzz_xxyyy[k] * ab_y + g_x_0_yzzzz_xxyyyy[k];

                g_x_0_yyzzzz_xxyyz[k] = -g_x_0_yzzzz_xxyyz[k] * ab_y + g_x_0_yzzzz_xxyyyz[k];

                g_x_0_yyzzzz_xxyzz[k] = -g_x_0_yzzzz_xxyzz[k] * ab_y + g_x_0_yzzzz_xxyyzz[k];

                g_x_0_yyzzzz_xxzzz[k] = -g_x_0_yzzzz_xxzzz[k] * ab_y + g_x_0_yzzzz_xxyzzz[k];

                g_x_0_yyzzzz_xyyyy[k] = -g_x_0_yzzzz_xyyyy[k] * ab_y + g_x_0_yzzzz_xyyyyy[k];

                g_x_0_yyzzzz_xyyyz[k] = -g_x_0_yzzzz_xyyyz[k] * ab_y + g_x_0_yzzzz_xyyyyz[k];

                g_x_0_yyzzzz_xyyzz[k] = -g_x_0_yzzzz_xyyzz[k] * ab_y + g_x_0_yzzzz_xyyyzz[k];

                g_x_0_yyzzzz_xyzzz[k] = -g_x_0_yzzzz_xyzzz[k] * ab_y + g_x_0_yzzzz_xyyzzz[k];

                g_x_0_yyzzzz_xzzzz[k] = -g_x_0_yzzzz_xzzzz[k] * ab_y + g_x_0_yzzzz_xyzzzz[k];

                g_x_0_yyzzzz_yyyyy[k] = -g_x_0_yzzzz_yyyyy[k] * ab_y + g_x_0_yzzzz_yyyyyy[k];

                g_x_0_yyzzzz_yyyyz[k] = -g_x_0_yzzzz_yyyyz[k] * ab_y + g_x_0_yzzzz_yyyyyz[k];

                g_x_0_yyzzzz_yyyzz[k] = -g_x_0_yzzzz_yyyzz[k] * ab_y + g_x_0_yzzzz_yyyyzz[k];

                g_x_0_yyzzzz_yyzzz[k] = -g_x_0_yzzzz_yyzzz[k] * ab_y + g_x_0_yzzzz_yyyzzz[k];

                g_x_0_yyzzzz_yzzzz[k] = -g_x_0_yzzzz_yzzzz[k] * ab_y + g_x_0_yzzzz_yyzzzz[k];

                g_x_0_yyzzzz_zzzzz[k] = -g_x_0_yzzzz_zzzzz[k] * ab_y + g_x_0_yzzzz_yzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 546 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 547 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 548 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 549 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 550 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 551 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 552 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 553 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 554 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 555 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 556 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 557 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 558 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 559 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 560 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 561 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 562 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 563 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 564 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 565 * ccomps * dcomps);

            auto g_x_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzzz_xxxxx, g_x_0_yzzzzz_xxxxy, g_x_0_yzzzzz_xxxxz, g_x_0_yzzzzz_xxxyy, g_x_0_yzzzzz_xxxyz, g_x_0_yzzzzz_xxxzz, g_x_0_yzzzzz_xxyyy, g_x_0_yzzzzz_xxyyz, g_x_0_yzzzzz_xxyzz, g_x_0_yzzzzz_xxzzz, g_x_0_yzzzzz_xyyyy, g_x_0_yzzzzz_xyyyz, g_x_0_yzzzzz_xyyzz, g_x_0_yzzzzz_xyzzz, g_x_0_yzzzzz_xzzzz, g_x_0_yzzzzz_yyyyy, g_x_0_yzzzzz_yyyyz, g_x_0_yzzzzz_yyyzz, g_x_0_yzzzzz_yyzzz, g_x_0_yzzzzz_yzzzz, g_x_0_yzzzzz_zzzzz, g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxxxx[k] = -g_x_0_zzzzz_xxxxx[k] * ab_y + g_x_0_zzzzz_xxxxxy[k];

                g_x_0_yzzzzz_xxxxy[k] = -g_x_0_zzzzz_xxxxy[k] * ab_y + g_x_0_zzzzz_xxxxyy[k];

                g_x_0_yzzzzz_xxxxz[k] = -g_x_0_zzzzz_xxxxz[k] * ab_y + g_x_0_zzzzz_xxxxyz[k];

                g_x_0_yzzzzz_xxxyy[k] = -g_x_0_zzzzz_xxxyy[k] * ab_y + g_x_0_zzzzz_xxxyyy[k];

                g_x_0_yzzzzz_xxxyz[k] = -g_x_0_zzzzz_xxxyz[k] * ab_y + g_x_0_zzzzz_xxxyyz[k];

                g_x_0_yzzzzz_xxxzz[k] = -g_x_0_zzzzz_xxxzz[k] * ab_y + g_x_0_zzzzz_xxxyzz[k];

                g_x_0_yzzzzz_xxyyy[k] = -g_x_0_zzzzz_xxyyy[k] * ab_y + g_x_0_zzzzz_xxyyyy[k];

                g_x_0_yzzzzz_xxyyz[k] = -g_x_0_zzzzz_xxyyz[k] * ab_y + g_x_0_zzzzz_xxyyyz[k];

                g_x_0_yzzzzz_xxyzz[k] = -g_x_0_zzzzz_xxyzz[k] * ab_y + g_x_0_zzzzz_xxyyzz[k];

                g_x_0_yzzzzz_xxzzz[k] = -g_x_0_zzzzz_xxzzz[k] * ab_y + g_x_0_zzzzz_xxyzzz[k];

                g_x_0_yzzzzz_xyyyy[k] = -g_x_0_zzzzz_xyyyy[k] * ab_y + g_x_0_zzzzz_xyyyyy[k];

                g_x_0_yzzzzz_xyyyz[k] = -g_x_0_zzzzz_xyyyz[k] * ab_y + g_x_0_zzzzz_xyyyyz[k];

                g_x_0_yzzzzz_xyyzz[k] = -g_x_0_zzzzz_xyyzz[k] * ab_y + g_x_0_zzzzz_xyyyzz[k];

                g_x_0_yzzzzz_xyzzz[k] = -g_x_0_zzzzz_xyzzz[k] * ab_y + g_x_0_zzzzz_xyyzzz[k];

                g_x_0_yzzzzz_xzzzz[k] = -g_x_0_zzzzz_xzzzz[k] * ab_y + g_x_0_zzzzz_xyzzzz[k];

                g_x_0_yzzzzz_yyyyy[k] = -g_x_0_zzzzz_yyyyy[k] * ab_y + g_x_0_zzzzz_yyyyyy[k];

                g_x_0_yzzzzz_yyyyz[k] = -g_x_0_zzzzz_yyyyz[k] * ab_y + g_x_0_zzzzz_yyyyyz[k];

                g_x_0_yzzzzz_yyyzz[k] = -g_x_0_zzzzz_yyyzz[k] * ab_y + g_x_0_zzzzz_yyyyzz[k];

                g_x_0_yzzzzz_yyzzz[k] = -g_x_0_zzzzz_yyzzz[k] * ab_y + g_x_0_zzzzz_yyyzzz[k];

                g_x_0_yzzzzz_yzzzz[k] = -g_x_0_zzzzz_yzzzz[k] * ab_y + g_x_0_zzzzz_yyzzzz[k];

                g_x_0_yzzzzz_zzzzz[k] = -g_x_0_zzzzz_zzzzz[k] * ab_y + g_x_0_zzzzz_yzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 567 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 568 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 569 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 570 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 571 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 572 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 573 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 574 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 575 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 576 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 577 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 578 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 579 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 580 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 581 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 582 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 583 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 584 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 585 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 586 * ccomps * dcomps);

            auto g_x_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzz, g_x_0_zzzzz_zzzzzz, g_x_0_zzzzzz_xxxxx, g_x_0_zzzzzz_xxxxy, g_x_0_zzzzzz_xxxxz, g_x_0_zzzzzz_xxxyy, g_x_0_zzzzzz_xxxyz, g_x_0_zzzzzz_xxxzz, g_x_0_zzzzzz_xxyyy, g_x_0_zzzzzz_xxyyz, g_x_0_zzzzzz_xxyzz, g_x_0_zzzzzz_xxzzz, g_x_0_zzzzzz_xyyyy, g_x_0_zzzzzz_xyyyz, g_x_0_zzzzzz_xyyzz, g_x_0_zzzzzz_xyzzz, g_x_0_zzzzzz_xzzzz, g_x_0_zzzzzz_yyyyy, g_x_0_zzzzzz_yyyyz, g_x_0_zzzzzz_yyyzz, g_x_0_zzzzzz_yyzzz, g_x_0_zzzzzz_yzzzz, g_x_0_zzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxxxx[k] = -g_x_0_zzzzz_xxxxx[k] * ab_z + g_x_0_zzzzz_xxxxxz[k];

                g_x_0_zzzzzz_xxxxy[k] = -g_x_0_zzzzz_xxxxy[k] * ab_z + g_x_0_zzzzz_xxxxyz[k];

                g_x_0_zzzzzz_xxxxz[k] = -g_x_0_zzzzz_xxxxz[k] * ab_z + g_x_0_zzzzz_xxxxzz[k];

                g_x_0_zzzzzz_xxxyy[k] = -g_x_0_zzzzz_xxxyy[k] * ab_z + g_x_0_zzzzz_xxxyyz[k];

                g_x_0_zzzzzz_xxxyz[k] = -g_x_0_zzzzz_xxxyz[k] * ab_z + g_x_0_zzzzz_xxxyzz[k];

                g_x_0_zzzzzz_xxxzz[k] = -g_x_0_zzzzz_xxxzz[k] * ab_z + g_x_0_zzzzz_xxxzzz[k];

                g_x_0_zzzzzz_xxyyy[k] = -g_x_0_zzzzz_xxyyy[k] * ab_z + g_x_0_zzzzz_xxyyyz[k];

                g_x_0_zzzzzz_xxyyz[k] = -g_x_0_zzzzz_xxyyz[k] * ab_z + g_x_0_zzzzz_xxyyzz[k];

                g_x_0_zzzzzz_xxyzz[k] = -g_x_0_zzzzz_xxyzz[k] * ab_z + g_x_0_zzzzz_xxyzzz[k];

                g_x_0_zzzzzz_xxzzz[k] = -g_x_0_zzzzz_xxzzz[k] * ab_z + g_x_0_zzzzz_xxzzzz[k];

                g_x_0_zzzzzz_xyyyy[k] = -g_x_0_zzzzz_xyyyy[k] * ab_z + g_x_0_zzzzz_xyyyyz[k];

                g_x_0_zzzzzz_xyyyz[k] = -g_x_0_zzzzz_xyyyz[k] * ab_z + g_x_0_zzzzz_xyyyzz[k];

                g_x_0_zzzzzz_xyyzz[k] = -g_x_0_zzzzz_xyyzz[k] * ab_z + g_x_0_zzzzz_xyyzzz[k];

                g_x_0_zzzzzz_xyzzz[k] = -g_x_0_zzzzz_xyzzz[k] * ab_z + g_x_0_zzzzz_xyzzzz[k];

                g_x_0_zzzzzz_xzzzz[k] = -g_x_0_zzzzz_xzzzz[k] * ab_z + g_x_0_zzzzz_xzzzzz[k];

                g_x_0_zzzzzz_yyyyy[k] = -g_x_0_zzzzz_yyyyy[k] * ab_z + g_x_0_zzzzz_yyyyyz[k];

                g_x_0_zzzzzz_yyyyz[k] = -g_x_0_zzzzz_yyyyz[k] * ab_z + g_x_0_zzzzz_yyyyzz[k];

                g_x_0_zzzzzz_yyyzz[k] = -g_x_0_zzzzz_yyyzz[k] * ab_z + g_x_0_zzzzz_yyyzzz[k];

                g_x_0_zzzzzz_yyzzz[k] = -g_x_0_zzzzz_yyzzz[k] * ab_z + g_x_0_zzzzz_yyzzzz[k];

                g_x_0_zzzzzz_yzzzz[k] = -g_x_0_zzzzz_yzzzz[k] * ab_z + g_x_0_zzzzz_yzzzzz[k];

                g_x_0_zzzzzz_zzzzz[k] = -g_x_0_zzzzz_zzzzz[k] * ab_z + g_x_0_zzzzz_zzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxxx, g_y_0_xxxxx_xxxxxy, g_y_0_xxxxx_xxxxxz, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxyy, g_y_0_xxxxx_xxxxyz, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxxzz, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyyy, g_y_0_xxxxx_xxxyyz, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxyzz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxxzzz, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyyy, g_y_0_xxxxx_xxyyyz, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyyzz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxyzzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xxzzzz, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyyy, g_y_0_xxxxx_xyyyyz, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyyzz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyyzzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xyzzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_xzzzzz, g_y_0_xxxxx_yyyyy, g_y_0_xxxxx_yyyyz, g_y_0_xxxxx_yyyzz, g_y_0_xxxxx_yyzzz, g_y_0_xxxxx_yzzzz, g_y_0_xxxxx_zzzzz, g_y_0_xxxxxx_xxxxx, g_y_0_xxxxxx_xxxxy, g_y_0_xxxxxx_xxxxz, g_y_0_xxxxxx_xxxyy, g_y_0_xxxxxx_xxxyz, g_y_0_xxxxxx_xxxzz, g_y_0_xxxxxx_xxyyy, g_y_0_xxxxxx_xxyyz, g_y_0_xxxxxx_xxyzz, g_y_0_xxxxxx_xxzzz, g_y_0_xxxxxx_xyyyy, g_y_0_xxxxxx_xyyyz, g_y_0_xxxxxx_xyyzz, g_y_0_xxxxxx_xyzzz, g_y_0_xxxxxx_xzzzz, g_y_0_xxxxxx_yyyyy, g_y_0_xxxxxx_yyyyz, g_y_0_xxxxxx_yyyzz, g_y_0_xxxxxx_yyzzz, g_y_0_xxxxxx_yzzzz, g_y_0_xxxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxxxx[k] = -g_y_0_xxxxx_xxxxx[k] * ab_x + g_y_0_xxxxx_xxxxxx[k];

                g_y_0_xxxxxx_xxxxy[k] = -g_y_0_xxxxx_xxxxy[k] * ab_x + g_y_0_xxxxx_xxxxxy[k];

                g_y_0_xxxxxx_xxxxz[k] = -g_y_0_xxxxx_xxxxz[k] * ab_x + g_y_0_xxxxx_xxxxxz[k];

                g_y_0_xxxxxx_xxxyy[k] = -g_y_0_xxxxx_xxxyy[k] * ab_x + g_y_0_xxxxx_xxxxyy[k];

                g_y_0_xxxxxx_xxxyz[k] = -g_y_0_xxxxx_xxxyz[k] * ab_x + g_y_0_xxxxx_xxxxyz[k];

                g_y_0_xxxxxx_xxxzz[k] = -g_y_0_xxxxx_xxxzz[k] * ab_x + g_y_0_xxxxx_xxxxzz[k];

                g_y_0_xxxxxx_xxyyy[k] = -g_y_0_xxxxx_xxyyy[k] * ab_x + g_y_0_xxxxx_xxxyyy[k];

                g_y_0_xxxxxx_xxyyz[k] = -g_y_0_xxxxx_xxyyz[k] * ab_x + g_y_0_xxxxx_xxxyyz[k];

                g_y_0_xxxxxx_xxyzz[k] = -g_y_0_xxxxx_xxyzz[k] * ab_x + g_y_0_xxxxx_xxxyzz[k];

                g_y_0_xxxxxx_xxzzz[k] = -g_y_0_xxxxx_xxzzz[k] * ab_x + g_y_0_xxxxx_xxxzzz[k];

                g_y_0_xxxxxx_xyyyy[k] = -g_y_0_xxxxx_xyyyy[k] * ab_x + g_y_0_xxxxx_xxyyyy[k];

                g_y_0_xxxxxx_xyyyz[k] = -g_y_0_xxxxx_xyyyz[k] * ab_x + g_y_0_xxxxx_xxyyyz[k];

                g_y_0_xxxxxx_xyyzz[k] = -g_y_0_xxxxx_xyyzz[k] * ab_x + g_y_0_xxxxx_xxyyzz[k];

                g_y_0_xxxxxx_xyzzz[k] = -g_y_0_xxxxx_xyzzz[k] * ab_x + g_y_0_xxxxx_xxyzzz[k];

                g_y_0_xxxxxx_xzzzz[k] = -g_y_0_xxxxx_xzzzz[k] * ab_x + g_y_0_xxxxx_xxzzzz[k];

                g_y_0_xxxxxx_yyyyy[k] = -g_y_0_xxxxx_yyyyy[k] * ab_x + g_y_0_xxxxx_xyyyyy[k];

                g_y_0_xxxxxx_yyyyz[k] = -g_y_0_xxxxx_yyyyz[k] * ab_x + g_y_0_xxxxx_xyyyyz[k];

                g_y_0_xxxxxx_yyyzz[k] = -g_y_0_xxxxx_yyyzz[k] * ab_x + g_y_0_xxxxx_xyyyzz[k];

                g_y_0_xxxxxx_yyzzz[k] = -g_y_0_xxxxx_yyzzz[k] * ab_x + g_y_0_xxxxx_xyyzzz[k];

                g_y_0_xxxxxx_yzzzz[k] = -g_y_0_xxxxx_yzzzz[k] * ab_x + g_y_0_xxxxx_xyzzzz[k];

                g_y_0_xxxxxx_zzzzz[k] = -g_y_0_xxxxx_zzzzz[k] * ab_x + g_y_0_xxxxx_xzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxy_xxxxx, g_y_0_xxxxxy_xxxxy, g_y_0_xxxxxy_xxxxz, g_y_0_xxxxxy_xxxyy, g_y_0_xxxxxy_xxxyz, g_y_0_xxxxxy_xxxzz, g_y_0_xxxxxy_xxyyy, g_y_0_xxxxxy_xxyyz, g_y_0_xxxxxy_xxyzz, g_y_0_xxxxxy_xxzzz, g_y_0_xxxxxy_xyyyy, g_y_0_xxxxxy_xyyyz, g_y_0_xxxxxy_xyyzz, g_y_0_xxxxxy_xyzzz, g_y_0_xxxxxy_xzzzz, g_y_0_xxxxxy_yyyyy, g_y_0_xxxxxy_yyyyz, g_y_0_xxxxxy_yyyzz, g_y_0_xxxxxy_yyzzz, g_y_0_xxxxxy_yzzzz, g_y_0_xxxxxy_zzzzz, g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxxx, g_y_0_xxxxy_xxxxxy, g_y_0_xxxxy_xxxxxz, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxyy, g_y_0_xxxxy_xxxxyz, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxxzz, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyyy, g_y_0_xxxxy_xxxyyz, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxyzz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxxzzz, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyyy, g_y_0_xxxxy_xxyyyz, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyyzz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxyzzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xxzzzz, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyyy, g_y_0_xxxxy_xyyyyz, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyyzz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyyzzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xyzzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_xzzzzz, g_y_0_xxxxy_yyyyy, g_y_0_xxxxy_yyyyz, g_y_0_xxxxy_yyyzz, g_y_0_xxxxy_yyzzz, g_y_0_xxxxy_yzzzz, g_y_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxxxx[k] = -g_y_0_xxxxy_xxxxx[k] * ab_x + g_y_0_xxxxy_xxxxxx[k];

                g_y_0_xxxxxy_xxxxy[k] = -g_y_0_xxxxy_xxxxy[k] * ab_x + g_y_0_xxxxy_xxxxxy[k];

                g_y_0_xxxxxy_xxxxz[k] = -g_y_0_xxxxy_xxxxz[k] * ab_x + g_y_0_xxxxy_xxxxxz[k];

                g_y_0_xxxxxy_xxxyy[k] = -g_y_0_xxxxy_xxxyy[k] * ab_x + g_y_0_xxxxy_xxxxyy[k];

                g_y_0_xxxxxy_xxxyz[k] = -g_y_0_xxxxy_xxxyz[k] * ab_x + g_y_0_xxxxy_xxxxyz[k];

                g_y_0_xxxxxy_xxxzz[k] = -g_y_0_xxxxy_xxxzz[k] * ab_x + g_y_0_xxxxy_xxxxzz[k];

                g_y_0_xxxxxy_xxyyy[k] = -g_y_0_xxxxy_xxyyy[k] * ab_x + g_y_0_xxxxy_xxxyyy[k];

                g_y_0_xxxxxy_xxyyz[k] = -g_y_0_xxxxy_xxyyz[k] * ab_x + g_y_0_xxxxy_xxxyyz[k];

                g_y_0_xxxxxy_xxyzz[k] = -g_y_0_xxxxy_xxyzz[k] * ab_x + g_y_0_xxxxy_xxxyzz[k];

                g_y_0_xxxxxy_xxzzz[k] = -g_y_0_xxxxy_xxzzz[k] * ab_x + g_y_0_xxxxy_xxxzzz[k];

                g_y_0_xxxxxy_xyyyy[k] = -g_y_0_xxxxy_xyyyy[k] * ab_x + g_y_0_xxxxy_xxyyyy[k];

                g_y_0_xxxxxy_xyyyz[k] = -g_y_0_xxxxy_xyyyz[k] * ab_x + g_y_0_xxxxy_xxyyyz[k];

                g_y_0_xxxxxy_xyyzz[k] = -g_y_0_xxxxy_xyyzz[k] * ab_x + g_y_0_xxxxy_xxyyzz[k];

                g_y_0_xxxxxy_xyzzz[k] = -g_y_0_xxxxy_xyzzz[k] * ab_x + g_y_0_xxxxy_xxyzzz[k];

                g_y_0_xxxxxy_xzzzz[k] = -g_y_0_xxxxy_xzzzz[k] * ab_x + g_y_0_xxxxy_xxzzzz[k];

                g_y_0_xxxxxy_yyyyy[k] = -g_y_0_xxxxy_yyyyy[k] * ab_x + g_y_0_xxxxy_xyyyyy[k];

                g_y_0_xxxxxy_yyyyz[k] = -g_y_0_xxxxy_yyyyz[k] * ab_x + g_y_0_xxxxy_xyyyyz[k];

                g_y_0_xxxxxy_yyyzz[k] = -g_y_0_xxxxy_yyyzz[k] * ab_x + g_y_0_xxxxy_xyyyzz[k];

                g_y_0_xxxxxy_yyzzz[k] = -g_y_0_xxxxy_yyzzz[k] * ab_x + g_y_0_xxxxy_xyyzzz[k];

                g_y_0_xxxxxy_yzzzz[k] = -g_y_0_xxxxy_yzzzz[k] * ab_x + g_y_0_xxxxy_xyzzzz[k];

                g_y_0_xxxxxy_zzzzz[k] = -g_y_0_xxxxy_zzzzz[k] * ab_x + g_y_0_xxxxy_xzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxz_xxxxx, g_y_0_xxxxxz_xxxxy, g_y_0_xxxxxz_xxxxz, g_y_0_xxxxxz_xxxyy, g_y_0_xxxxxz_xxxyz, g_y_0_xxxxxz_xxxzz, g_y_0_xxxxxz_xxyyy, g_y_0_xxxxxz_xxyyz, g_y_0_xxxxxz_xxyzz, g_y_0_xxxxxz_xxzzz, g_y_0_xxxxxz_xyyyy, g_y_0_xxxxxz_xyyyz, g_y_0_xxxxxz_xyyzz, g_y_0_xxxxxz_xyzzz, g_y_0_xxxxxz_xzzzz, g_y_0_xxxxxz_yyyyy, g_y_0_xxxxxz_yyyyz, g_y_0_xxxxxz_yyyzz, g_y_0_xxxxxz_yyzzz, g_y_0_xxxxxz_yzzzz, g_y_0_xxxxxz_zzzzz, g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxxx, g_y_0_xxxxz_xxxxxy, g_y_0_xxxxz_xxxxxz, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxyy, g_y_0_xxxxz_xxxxyz, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxxzz, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyyy, g_y_0_xxxxz_xxxyyz, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxyzz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxxzzz, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyyy, g_y_0_xxxxz_xxyyyz, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyyzz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxyzzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xxzzzz, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyyy, g_y_0_xxxxz_xyyyyz, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyyzz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyyzzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xyzzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_xzzzzz, g_y_0_xxxxz_yyyyy, g_y_0_xxxxz_yyyyz, g_y_0_xxxxz_yyyzz, g_y_0_xxxxz_yyzzz, g_y_0_xxxxz_yzzzz, g_y_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxxxx[k] = -g_y_0_xxxxz_xxxxx[k] * ab_x + g_y_0_xxxxz_xxxxxx[k];

                g_y_0_xxxxxz_xxxxy[k] = -g_y_0_xxxxz_xxxxy[k] * ab_x + g_y_0_xxxxz_xxxxxy[k];

                g_y_0_xxxxxz_xxxxz[k] = -g_y_0_xxxxz_xxxxz[k] * ab_x + g_y_0_xxxxz_xxxxxz[k];

                g_y_0_xxxxxz_xxxyy[k] = -g_y_0_xxxxz_xxxyy[k] * ab_x + g_y_0_xxxxz_xxxxyy[k];

                g_y_0_xxxxxz_xxxyz[k] = -g_y_0_xxxxz_xxxyz[k] * ab_x + g_y_0_xxxxz_xxxxyz[k];

                g_y_0_xxxxxz_xxxzz[k] = -g_y_0_xxxxz_xxxzz[k] * ab_x + g_y_0_xxxxz_xxxxzz[k];

                g_y_0_xxxxxz_xxyyy[k] = -g_y_0_xxxxz_xxyyy[k] * ab_x + g_y_0_xxxxz_xxxyyy[k];

                g_y_0_xxxxxz_xxyyz[k] = -g_y_0_xxxxz_xxyyz[k] * ab_x + g_y_0_xxxxz_xxxyyz[k];

                g_y_0_xxxxxz_xxyzz[k] = -g_y_0_xxxxz_xxyzz[k] * ab_x + g_y_0_xxxxz_xxxyzz[k];

                g_y_0_xxxxxz_xxzzz[k] = -g_y_0_xxxxz_xxzzz[k] * ab_x + g_y_0_xxxxz_xxxzzz[k];

                g_y_0_xxxxxz_xyyyy[k] = -g_y_0_xxxxz_xyyyy[k] * ab_x + g_y_0_xxxxz_xxyyyy[k];

                g_y_0_xxxxxz_xyyyz[k] = -g_y_0_xxxxz_xyyyz[k] * ab_x + g_y_0_xxxxz_xxyyyz[k];

                g_y_0_xxxxxz_xyyzz[k] = -g_y_0_xxxxz_xyyzz[k] * ab_x + g_y_0_xxxxz_xxyyzz[k];

                g_y_0_xxxxxz_xyzzz[k] = -g_y_0_xxxxz_xyzzz[k] * ab_x + g_y_0_xxxxz_xxyzzz[k];

                g_y_0_xxxxxz_xzzzz[k] = -g_y_0_xxxxz_xzzzz[k] * ab_x + g_y_0_xxxxz_xxzzzz[k];

                g_y_0_xxxxxz_yyyyy[k] = -g_y_0_xxxxz_yyyyy[k] * ab_x + g_y_0_xxxxz_xyyyyy[k];

                g_y_0_xxxxxz_yyyyz[k] = -g_y_0_xxxxz_yyyyz[k] * ab_x + g_y_0_xxxxz_xyyyyz[k];

                g_y_0_xxxxxz_yyyzz[k] = -g_y_0_xxxxz_yyyzz[k] * ab_x + g_y_0_xxxxz_xyyyzz[k];

                g_y_0_xxxxxz_yyzzz[k] = -g_y_0_xxxxz_yyzzz[k] * ab_x + g_y_0_xxxxz_xyyzzz[k];

                g_y_0_xxxxxz_yzzzz[k] = -g_y_0_xxxxz_yzzzz[k] * ab_x + g_y_0_xxxxz_xyzzzz[k];

                g_y_0_xxxxxz_zzzzz[k] = -g_y_0_xxxxz_zzzzz[k] * ab_x + g_y_0_xxxxz_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyy_xxxxx, g_y_0_xxxxyy_xxxxy, g_y_0_xxxxyy_xxxxz, g_y_0_xxxxyy_xxxyy, g_y_0_xxxxyy_xxxyz, g_y_0_xxxxyy_xxxzz, g_y_0_xxxxyy_xxyyy, g_y_0_xxxxyy_xxyyz, g_y_0_xxxxyy_xxyzz, g_y_0_xxxxyy_xxzzz, g_y_0_xxxxyy_xyyyy, g_y_0_xxxxyy_xyyyz, g_y_0_xxxxyy_xyyzz, g_y_0_xxxxyy_xyzzz, g_y_0_xxxxyy_xzzzz, g_y_0_xxxxyy_yyyyy, g_y_0_xxxxyy_yyyyz, g_y_0_xxxxyy_yyyzz, g_y_0_xxxxyy_yyzzz, g_y_0_xxxxyy_yzzzz, g_y_0_xxxxyy_zzzzz, g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxxx, g_y_0_xxxyy_xxxxxy, g_y_0_xxxyy_xxxxxz, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxyy, g_y_0_xxxyy_xxxxyz, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxxzz, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyyy, g_y_0_xxxyy_xxxyyz, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxyzz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxxzzz, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyyy, g_y_0_xxxyy_xxyyyz, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyyzz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxyzzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xxzzzz, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyyy, g_y_0_xxxyy_xyyyyz, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyyzz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyyzzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xyzzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_xzzzzz, g_y_0_xxxyy_yyyyy, g_y_0_xxxyy_yyyyz, g_y_0_xxxyy_yyyzz, g_y_0_xxxyy_yyzzz, g_y_0_xxxyy_yzzzz, g_y_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxxxx[k] = -g_y_0_xxxyy_xxxxx[k] * ab_x + g_y_0_xxxyy_xxxxxx[k];

                g_y_0_xxxxyy_xxxxy[k] = -g_y_0_xxxyy_xxxxy[k] * ab_x + g_y_0_xxxyy_xxxxxy[k];

                g_y_0_xxxxyy_xxxxz[k] = -g_y_0_xxxyy_xxxxz[k] * ab_x + g_y_0_xxxyy_xxxxxz[k];

                g_y_0_xxxxyy_xxxyy[k] = -g_y_0_xxxyy_xxxyy[k] * ab_x + g_y_0_xxxyy_xxxxyy[k];

                g_y_0_xxxxyy_xxxyz[k] = -g_y_0_xxxyy_xxxyz[k] * ab_x + g_y_0_xxxyy_xxxxyz[k];

                g_y_0_xxxxyy_xxxzz[k] = -g_y_0_xxxyy_xxxzz[k] * ab_x + g_y_0_xxxyy_xxxxzz[k];

                g_y_0_xxxxyy_xxyyy[k] = -g_y_0_xxxyy_xxyyy[k] * ab_x + g_y_0_xxxyy_xxxyyy[k];

                g_y_0_xxxxyy_xxyyz[k] = -g_y_0_xxxyy_xxyyz[k] * ab_x + g_y_0_xxxyy_xxxyyz[k];

                g_y_0_xxxxyy_xxyzz[k] = -g_y_0_xxxyy_xxyzz[k] * ab_x + g_y_0_xxxyy_xxxyzz[k];

                g_y_0_xxxxyy_xxzzz[k] = -g_y_0_xxxyy_xxzzz[k] * ab_x + g_y_0_xxxyy_xxxzzz[k];

                g_y_0_xxxxyy_xyyyy[k] = -g_y_0_xxxyy_xyyyy[k] * ab_x + g_y_0_xxxyy_xxyyyy[k];

                g_y_0_xxxxyy_xyyyz[k] = -g_y_0_xxxyy_xyyyz[k] * ab_x + g_y_0_xxxyy_xxyyyz[k];

                g_y_0_xxxxyy_xyyzz[k] = -g_y_0_xxxyy_xyyzz[k] * ab_x + g_y_0_xxxyy_xxyyzz[k];

                g_y_0_xxxxyy_xyzzz[k] = -g_y_0_xxxyy_xyzzz[k] * ab_x + g_y_0_xxxyy_xxyzzz[k];

                g_y_0_xxxxyy_xzzzz[k] = -g_y_0_xxxyy_xzzzz[k] * ab_x + g_y_0_xxxyy_xxzzzz[k];

                g_y_0_xxxxyy_yyyyy[k] = -g_y_0_xxxyy_yyyyy[k] * ab_x + g_y_0_xxxyy_xyyyyy[k];

                g_y_0_xxxxyy_yyyyz[k] = -g_y_0_xxxyy_yyyyz[k] * ab_x + g_y_0_xxxyy_xyyyyz[k];

                g_y_0_xxxxyy_yyyzz[k] = -g_y_0_xxxyy_yyyzz[k] * ab_x + g_y_0_xxxyy_xyyyzz[k];

                g_y_0_xxxxyy_yyzzz[k] = -g_y_0_xxxyy_yyzzz[k] * ab_x + g_y_0_xxxyy_xyyzzz[k];

                g_y_0_xxxxyy_yzzzz[k] = -g_y_0_xxxyy_yzzzz[k] * ab_x + g_y_0_xxxyy_xyzzzz[k];

                g_y_0_xxxxyy_zzzzz[k] = -g_y_0_xxxyy_zzzzz[k] * ab_x + g_y_0_xxxyy_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyz_xxxxx, g_y_0_xxxxyz_xxxxy, g_y_0_xxxxyz_xxxxz, g_y_0_xxxxyz_xxxyy, g_y_0_xxxxyz_xxxyz, g_y_0_xxxxyz_xxxzz, g_y_0_xxxxyz_xxyyy, g_y_0_xxxxyz_xxyyz, g_y_0_xxxxyz_xxyzz, g_y_0_xxxxyz_xxzzz, g_y_0_xxxxyz_xyyyy, g_y_0_xxxxyz_xyyyz, g_y_0_xxxxyz_xyyzz, g_y_0_xxxxyz_xyzzz, g_y_0_xxxxyz_xzzzz, g_y_0_xxxxyz_yyyyy, g_y_0_xxxxyz_yyyyz, g_y_0_xxxxyz_yyyzz, g_y_0_xxxxyz_yyzzz, g_y_0_xxxxyz_yzzzz, g_y_0_xxxxyz_zzzzz, g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxxx, g_y_0_xxxyz_xxxxxy, g_y_0_xxxyz_xxxxxz, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxyy, g_y_0_xxxyz_xxxxyz, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxxzz, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyyy, g_y_0_xxxyz_xxxyyz, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxyzz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxxzzz, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyyy, g_y_0_xxxyz_xxyyyz, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyyzz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxyzzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xxzzzz, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyyy, g_y_0_xxxyz_xyyyyz, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyyzz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyyzzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xyzzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_xzzzzz, g_y_0_xxxyz_yyyyy, g_y_0_xxxyz_yyyyz, g_y_0_xxxyz_yyyzz, g_y_0_xxxyz_yyzzz, g_y_0_xxxyz_yzzzz, g_y_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxxxx[k] = -g_y_0_xxxyz_xxxxx[k] * ab_x + g_y_0_xxxyz_xxxxxx[k];

                g_y_0_xxxxyz_xxxxy[k] = -g_y_0_xxxyz_xxxxy[k] * ab_x + g_y_0_xxxyz_xxxxxy[k];

                g_y_0_xxxxyz_xxxxz[k] = -g_y_0_xxxyz_xxxxz[k] * ab_x + g_y_0_xxxyz_xxxxxz[k];

                g_y_0_xxxxyz_xxxyy[k] = -g_y_0_xxxyz_xxxyy[k] * ab_x + g_y_0_xxxyz_xxxxyy[k];

                g_y_0_xxxxyz_xxxyz[k] = -g_y_0_xxxyz_xxxyz[k] * ab_x + g_y_0_xxxyz_xxxxyz[k];

                g_y_0_xxxxyz_xxxzz[k] = -g_y_0_xxxyz_xxxzz[k] * ab_x + g_y_0_xxxyz_xxxxzz[k];

                g_y_0_xxxxyz_xxyyy[k] = -g_y_0_xxxyz_xxyyy[k] * ab_x + g_y_0_xxxyz_xxxyyy[k];

                g_y_0_xxxxyz_xxyyz[k] = -g_y_0_xxxyz_xxyyz[k] * ab_x + g_y_0_xxxyz_xxxyyz[k];

                g_y_0_xxxxyz_xxyzz[k] = -g_y_0_xxxyz_xxyzz[k] * ab_x + g_y_0_xxxyz_xxxyzz[k];

                g_y_0_xxxxyz_xxzzz[k] = -g_y_0_xxxyz_xxzzz[k] * ab_x + g_y_0_xxxyz_xxxzzz[k];

                g_y_0_xxxxyz_xyyyy[k] = -g_y_0_xxxyz_xyyyy[k] * ab_x + g_y_0_xxxyz_xxyyyy[k];

                g_y_0_xxxxyz_xyyyz[k] = -g_y_0_xxxyz_xyyyz[k] * ab_x + g_y_0_xxxyz_xxyyyz[k];

                g_y_0_xxxxyz_xyyzz[k] = -g_y_0_xxxyz_xyyzz[k] * ab_x + g_y_0_xxxyz_xxyyzz[k];

                g_y_0_xxxxyz_xyzzz[k] = -g_y_0_xxxyz_xyzzz[k] * ab_x + g_y_0_xxxyz_xxyzzz[k];

                g_y_0_xxxxyz_xzzzz[k] = -g_y_0_xxxyz_xzzzz[k] * ab_x + g_y_0_xxxyz_xxzzzz[k];

                g_y_0_xxxxyz_yyyyy[k] = -g_y_0_xxxyz_yyyyy[k] * ab_x + g_y_0_xxxyz_xyyyyy[k];

                g_y_0_xxxxyz_yyyyz[k] = -g_y_0_xxxyz_yyyyz[k] * ab_x + g_y_0_xxxyz_xyyyyz[k];

                g_y_0_xxxxyz_yyyzz[k] = -g_y_0_xxxyz_yyyzz[k] * ab_x + g_y_0_xxxyz_xyyyzz[k];

                g_y_0_xxxxyz_yyzzz[k] = -g_y_0_xxxyz_yyzzz[k] * ab_x + g_y_0_xxxyz_xyyzzz[k];

                g_y_0_xxxxyz_yzzzz[k] = -g_y_0_xxxyz_yzzzz[k] * ab_x + g_y_0_xxxyz_xyzzzz[k];

                g_y_0_xxxxyz_zzzzz[k] = -g_y_0_xxxyz_zzzzz[k] * ab_x + g_y_0_xxxyz_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxzz_xxxxx, g_y_0_xxxxzz_xxxxy, g_y_0_xxxxzz_xxxxz, g_y_0_xxxxzz_xxxyy, g_y_0_xxxxzz_xxxyz, g_y_0_xxxxzz_xxxzz, g_y_0_xxxxzz_xxyyy, g_y_0_xxxxzz_xxyyz, g_y_0_xxxxzz_xxyzz, g_y_0_xxxxzz_xxzzz, g_y_0_xxxxzz_xyyyy, g_y_0_xxxxzz_xyyyz, g_y_0_xxxxzz_xyyzz, g_y_0_xxxxzz_xyzzz, g_y_0_xxxxzz_xzzzz, g_y_0_xxxxzz_yyyyy, g_y_0_xxxxzz_yyyyz, g_y_0_xxxxzz_yyyzz, g_y_0_xxxxzz_yyzzz, g_y_0_xxxxzz_yzzzz, g_y_0_xxxxzz_zzzzz, g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxxx, g_y_0_xxxzz_xxxxxy, g_y_0_xxxzz_xxxxxz, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxyy, g_y_0_xxxzz_xxxxyz, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxxzz, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyyy, g_y_0_xxxzz_xxxyyz, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxyzz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxxzzz, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyyy, g_y_0_xxxzz_xxyyyz, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyyzz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxyzzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xxzzzz, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyyy, g_y_0_xxxzz_xyyyyz, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyyzz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyyzzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xyzzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_xzzzzz, g_y_0_xxxzz_yyyyy, g_y_0_xxxzz_yyyyz, g_y_0_xxxzz_yyyzz, g_y_0_xxxzz_yyzzz, g_y_0_xxxzz_yzzzz, g_y_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxxxx[k] = -g_y_0_xxxzz_xxxxx[k] * ab_x + g_y_0_xxxzz_xxxxxx[k];

                g_y_0_xxxxzz_xxxxy[k] = -g_y_0_xxxzz_xxxxy[k] * ab_x + g_y_0_xxxzz_xxxxxy[k];

                g_y_0_xxxxzz_xxxxz[k] = -g_y_0_xxxzz_xxxxz[k] * ab_x + g_y_0_xxxzz_xxxxxz[k];

                g_y_0_xxxxzz_xxxyy[k] = -g_y_0_xxxzz_xxxyy[k] * ab_x + g_y_0_xxxzz_xxxxyy[k];

                g_y_0_xxxxzz_xxxyz[k] = -g_y_0_xxxzz_xxxyz[k] * ab_x + g_y_0_xxxzz_xxxxyz[k];

                g_y_0_xxxxzz_xxxzz[k] = -g_y_0_xxxzz_xxxzz[k] * ab_x + g_y_0_xxxzz_xxxxzz[k];

                g_y_0_xxxxzz_xxyyy[k] = -g_y_0_xxxzz_xxyyy[k] * ab_x + g_y_0_xxxzz_xxxyyy[k];

                g_y_0_xxxxzz_xxyyz[k] = -g_y_0_xxxzz_xxyyz[k] * ab_x + g_y_0_xxxzz_xxxyyz[k];

                g_y_0_xxxxzz_xxyzz[k] = -g_y_0_xxxzz_xxyzz[k] * ab_x + g_y_0_xxxzz_xxxyzz[k];

                g_y_0_xxxxzz_xxzzz[k] = -g_y_0_xxxzz_xxzzz[k] * ab_x + g_y_0_xxxzz_xxxzzz[k];

                g_y_0_xxxxzz_xyyyy[k] = -g_y_0_xxxzz_xyyyy[k] * ab_x + g_y_0_xxxzz_xxyyyy[k];

                g_y_0_xxxxzz_xyyyz[k] = -g_y_0_xxxzz_xyyyz[k] * ab_x + g_y_0_xxxzz_xxyyyz[k];

                g_y_0_xxxxzz_xyyzz[k] = -g_y_0_xxxzz_xyyzz[k] * ab_x + g_y_0_xxxzz_xxyyzz[k];

                g_y_0_xxxxzz_xyzzz[k] = -g_y_0_xxxzz_xyzzz[k] * ab_x + g_y_0_xxxzz_xxyzzz[k];

                g_y_0_xxxxzz_xzzzz[k] = -g_y_0_xxxzz_xzzzz[k] * ab_x + g_y_0_xxxzz_xxzzzz[k];

                g_y_0_xxxxzz_yyyyy[k] = -g_y_0_xxxzz_yyyyy[k] * ab_x + g_y_0_xxxzz_xyyyyy[k];

                g_y_0_xxxxzz_yyyyz[k] = -g_y_0_xxxzz_yyyyz[k] * ab_x + g_y_0_xxxzz_xyyyyz[k];

                g_y_0_xxxxzz_yyyzz[k] = -g_y_0_xxxzz_yyyzz[k] * ab_x + g_y_0_xxxzz_xyyyzz[k];

                g_y_0_xxxxzz_yyzzz[k] = -g_y_0_xxxzz_yyzzz[k] * ab_x + g_y_0_xxxzz_xyyzzz[k];

                g_y_0_xxxxzz_yzzzz[k] = -g_y_0_xxxzz_yzzzz[k] * ab_x + g_y_0_xxxzz_xyzzzz[k];

                g_y_0_xxxxzz_zzzzz[k] = -g_y_0_xxxzz_zzzzz[k] * ab_x + g_y_0_xxxzz_xzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyy_xxxxx, g_y_0_xxxyyy_xxxxy, g_y_0_xxxyyy_xxxxz, g_y_0_xxxyyy_xxxyy, g_y_0_xxxyyy_xxxyz, g_y_0_xxxyyy_xxxzz, g_y_0_xxxyyy_xxyyy, g_y_0_xxxyyy_xxyyz, g_y_0_xxxyyy_xxyzz, g_y_0_xxxyyy_xxzzz, g_y_0_xxxyyy_xyyyy, g_y_0_xxxyyy_xyyyz, g_y_0_xxxyyy_xyyzz, g_y_0_xxxyyy_xyzzz, g_y_0_xxxyyy_xzzzz, g_y_0_xxxyyy_yyyyy, g_y_0_xxxyyy_yyyyz, g_y_0_xxxyyy_yyyzz, g_y_0_xxxyyy_yyzzz, g_y_0_xxxyyy_yzzzz, g_y_0_xxxyyy_zzzzz, g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxxx, g_y_0_xxyyy_xxxxxy, g_y_0_xxyyy_xxxxxz, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxyy, g_y_0_xxyyy_xxxxyz, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxxzz, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyyy, g_y_0_xxyyy_xxxyyz, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxyzz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxxzzz, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyyy, g_y_0_xxyyy_xxyyyz, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyyzz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxyzzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xxzzzz, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyyy, g_y_0_xxyyy_xyyyyz, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyyzz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyyzzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xyzzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_xzzzzz, g_y_0_xxyyy_yyyyy, g_y_0_xxyyy_yyyyz, g_y_0_xxyyy_yyyzz, g_y_0_xxyyy_yyzzz, g_y_0_xxyyy_yzzzz, g_y_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxxxx[k] = -g_y_0_xxyyy_xxxxx[k] * ab_x + g_y_0_xxyyy_xxxxxx[k];

                g_y_0_xxxyyy_xxxxy[k] = -g_y_0_xxyyy_xxxxy[k] * ab_x + g_y_0_xxyyy_xxxxxy[k];

                g_y_0_xxxyyy_xxxxz[k] = -g_y_0_xxyyy_xxxxz[k] * ab_x + g_y_0_xxyyy_xxxxxz[k];

                g_y_0_xxxyyy_xxxyy[k] = -g_y_0_xxyyy_xxxyy[k] * ab_x + g_y_0_xxyyy_xxxxyy[k];

                g_y_0_xxxyyy_xxxyz[k] = -g_y_0_xxyyy_xxxyz[k] * ab_x + g_y_0_xxyyy_xxxxyz[k];

                g_y_0_xxxyyy_xxxzz[k] = -g_y_0_xxyyy_xxxzz[k] * ab_x + g_y_0_xxyyy_xxxxzz[k];

                g_y_0_xxxyyy_xxyyy[k] = -g_y_0_xxyyy_xxyyy[k] * ab_x + g_y_0_xxyyy_xxxyyy[k];

                g_y_0_xxxyyy_xxyyz[k] = -g_y_0_xxyyy_xxyyz[k] * ab_x + g_y_0_xxyyy_xxxyyz[k];

                g_y_0_xxxyyy_xxyzz[k] = -g_y_0_xxyyy_xxyzz[k] * ab_x + g_y_0_xxyyy_xxxyzz[k];

                g_y_0_xxxyyy_xxzzz[k] = -g_y_0_xxyyy_xxzzz[k] * ab_x + g_y_0_xxyyy_xxxzzz[k];

                g_y_0_xxxyyy_xyyyy[k] = -g_y_0_xxyyy_xyyyy[k] * ab_x + g_y_0_xxyyy_xxyyyy[k];

                g_y_0_xxxyyy_xyyyz[k] = -g_y_0_xxyyy_xyyyz[k] * ab_x + g_y_0_xxyyy_xxyyyz[k];

                g_y_0_xxxyyy_xyyzz[k] = -g_y_0_xxyyy_xyyzz[k] * ab_x + g_y_0_xxyyy_xxyyzz[k];

                g_y_0_xxxyyy_xyzzz[k] = -g_y_0_xxyyy_xyzzz[k] * ab_x + g_y_0_xxyyy_xxyzzz[k];

                g_y_0_xxxyyy_xzzzz[k] = -g_y_0_xxyyy_xzzzz[k] * ab_x + g_y_0_xxyyy_xxzzzz[k];

                g_y_0_xxxyyy_yyyyy[k] = -g_y_0_xxyyy_yyyyy[k] * ab_x + g_y_0_xxyyy_xyyyyy[k];

                g_y_0_xxxyyy_yyyyz[k] = -g_y_0_xxyyy_yyyyz[k] * ab_x + g_y_0_xxyyy_xyyyyz[k];

                g_y_0_xxxyyy_yyyzz[k] = -g_y_0_xxyyy_yyyzz[k] * ab_x + g_y_0_xxyyy_xyyyzz[k];

                g_y_0_xxxyyy_yyzzz[k] = -g_y_0_xxyyy_yyzzz[k] * ab_x + g_y_0_xxyyy_xyyzzz[k];

                g_y_0_xxxyyy_yzzzz[k] = -g_y_0_xxyyy_yzzzz[k] * ab_x + g_y_0_xxyyy_xyzzzz[k];

                g_y_0_xxxyyy_zzzzz[k] = -g_y_0_xxyyy_zzzzz[k] * ab_x + g_y_0_xxyyy_xzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyz_xxxxx, g_y_0_xxxyyz_xxxxy, g_y_0_xxxyyz_xxxxz, g_y_0_xxxyyz_xxxyy, g_y_0_xxxyyz_xxxyz, g_y_0_xxxyyz_xxxzz, g_y_0_xxxyyz_xxyyy, g_y_0_xxxyyz_xxyyz, g_y_0_xxxyyz_xxyzz, g_y_0_xxxyyz_xxzzz, g_y_0_xxxyyz_xyyyy, g_y_0_xxxyyz_xyyyz, g_y_0_xxxyyz_xyyzz, g_y_0_xxxyyz_xyzzz, g_y_0_xxxyyz_xzzzz, g_y_0_xxxyyz_yyyyy, g_y_0_xxxyyz_yyyyz, g_y_0_xxxyyz_yyyzz, g_y_0_xxxyyz_yyzzz, g_y_0_xxxyyz_yzzzz, g_y_0_xxxyyz_zzzzz, g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxxx, g_y_0_xxyyz_xxxxxy, g_y_0_xxyyz_xxxxxz, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxyy, g_y_0_xxyyz_xxxxyz, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxxzz, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyyy, g_y_0_xxyyz_xxxyyz, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxyzz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxxzzz, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyyy, g_y_0_xxyyz_xxyyyz, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyyzz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxyzzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xxzzzz, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyyy, g_y_0_xxyyz_xyyyyz, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyyzz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyyzzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xyzzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_xzzzzz, g_y_0_xxyyz_yyyyy, g_y_0_xxyyz_yyyyz, g_y_0_xxyyz_yyyzz, g_y_0_xxyyz_yyzzz, g_y_0_xxyyz_yzzzz, g_y_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxxxx[k] = -g_y_0_xxyyz_xxxxx[k] * ab_x + g_y_0_xxyyz_xxxxxx[k];

                g_y_0_xxxyyz_xxxxy[k] = -g_y_0_xxyyz_xxxxy[k] * ab_x + g_y_0_xxyyz_xxxxxy[k];

                g_y_0_xxxyyz_xxxxz[k] = -g_y_0_xxyyz_xxxxz[k] * ab_x + g_y_0_xxyyz_xxxxxz[k];

                g_y_0_xxxyyz_xxxyy[k] = -g_y_0_xxyyz_xxxyy[k] * ab_x + g_y_0_xxyyz_xxxxyy[k];

                g_y_0_xxxyyz_xxxyz[k] = -g_y_0_xxyyz_xxxyz[k] * ab_x + g_y_0_xxyyz_xxxxyz[k];

                g_y_0_xxxyyz_xxxzz[k] = -g_y_0_xxyyz_xxxzz[k] * ab_x + g_y_0_xxyyz_xxxxzz[k];

                g_y_0_xxxyyz_xxyyy[k] = -g_y_0_xxyyz_xxyyy[k] * ab_x + g_y_0_xxyyz_xxxyyy[k];

                g_y_0_xxxyyz_xxyyz[k] = -g_y_0_xxyyz_xxyyz[k] * ab_x + g_y_0_xxyyz_xxxyyz[k];

                g_y_0_xxxyyz_xxyzz[k] = -g_y_0_xxyyz_xxyzz[k] * ab_x + g_y_0_xxyyz_xxxyzz[k];

                g_y_0_xxxyyz_xxzzz[k] = -g_y_0_xxyyz_xxzzz[k] * ab_x + g_y_0_xxyyz_xxxzzz[k];

                g_y_0_xxxyyz_xyyyy[k] = -g_y_0_xxyyz_xyyyy[k] * ab_x + g_y_0_xxyyz_xxyyyy[k];

                g_y_0_xxxyyz_xyyyz[k] = -g_y_0_xxyyz_xyyyz[k] * ab_x + g_y_0_xxyyz_xxyyyz[k];

                g_y_0_xxxyyz_xyyzz[k] = -g_y_0_xxyyz_xyyzz[k] * ab_x + g_y_0_xxyyz_xxyyzz[k];

                g_y_0_xxxyyz_xyzzz[k] = -g_y_0_xxyyz_xyzzz[k] * ab_x + g_y_0_xxyyz_xxyzzz[k];

                g_y_0_xxxyyz_xzzzz[k] = -g_y_0_xxyyz_xzzzz[k] * ab_x + g_y_0_xxyyz_xxzzzz[k];

                g_y_0_xxxyyz_yyyyy[k] = -g_y_0_xxyyz_yyyyy[k] * ab_x + g_y_0_xxyyz_xyyyyy[k];

                g_y_0_xxxyyz_yyyyz[k] = -g_y_0_xxyyz_yyyyz[k] * ab_x + g_y_0_xxyyz_xyyyyz[k];

                g_y_0_xxxyyz_yyyzz[k] = -g_y_0_xxyyz_yyyzz[k] * ab_x + g_y_0_xxyyz_xyyyzz[k];

                g_y_0_xxxyyz_yyzzz[k] = -g_y_0_xxyyz_yyzzz[k] * ab_x + g_y_0_xxyyz_xyyzzz[k];

                g_y_0_xxxyyz_yzzzz[k] = -g_y_0_xxyyz_yzzzz[k] * ab_x + g_y_0_xxyyz_xyzzzz[k];

                g_y_0_xxxyyz_zzzzz[k] = -g_y_0_xxyyz_zzzzz[k] * ab_x + g_y_0_xxyyz_xzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 776 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyzz_xxxxx, g_y_0_xxxyzz_xxxxy, g_y_0_xxxyzz_xxxxz, g_y_0_xxxyzz_xxxyy, g_y_0_xxxyzz_xxxyz, g_y_0_xxxyzz_xxxzz, g_y_0_xxxyzz_xxyyy, g_y_0_xxxyzz_xxyyz, g_y_0_xxxyzz_xxyzz, g_y_0_xxxyzz_xxzzz, g_y_0_xxxyzz_xyyyy, g_y_0_xxxyzz_xyyyz, g_y_0_xxxyzz_xyyzz, g_y_0_xxxyzz_xyzzz, g_y_0_xxxyzz_xzzzz, g_y_0_xxxyzz_yyyyy, g_y_0_xxxyzz_yyyyz, g_y_0_xxxyzz_yyyzz, g_y_0_xxxyzz_yyzzz, g_y_0_xxxyzz_yzzzz, g_y_0_xxxyzz_zzzzz, g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxxx, g_y_0_xxyzz_xxxxxy, g_y_0_xxyzz_xxxxxz, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxyy, g_y_0_xxyzz_xxxxyz, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxxzz, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyyy, g_y_0_xxyzz_xxxyyz, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxyzz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxxzzz, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyyy, g_y_0_xxyzz_xxyyyz, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyyzz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxyzzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xxzzzz, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyyy, g_y_0_xxyzz_xyyyyz, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyyzz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyyzzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xyzzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_xzzzzz, g_y_0_xxyzz_yyyyy, g_y_0_xxyzz_yyyyz, g_y_0_xxyzz_yyyzz, g_y_0_xxyzz_yyzzz, g_y_0_xxyzz_yzzzz, g_y_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxxxx[k] = -g_y_0_xxyzz_xxxxx[k] * ab_x + g_y_0_xxyzz_xxxxxx[k];

                g_y_0_xxxyzz_xxxxy[k] = -g_y_0_xxyzz_xxxxy[k] * ab_x + g_y_0_xxyzz_xxxxxy[k];

                g_y_0_xxxyzz_xxxxz[k] = -g_y_0_xxyzz_xxxxz[k] * ab_x + g_y_0_xxyzz_xxxxxz[k];

                g_y_0_xxxyzz_xxxyy[k] = -g_y_0_xxyzz_xxxyy[k] * ab_x + g_y_0_xxyzz_xxxxyy[k];

                g_y_0_xxxyzz_xxxyz[k] = -g_y_0_xxyzz_xxxyz[k] * ab_x + g_y_0_xxyzz_xxxxyz[k];

                g_y_0_xxxyzz_xxxzz[k] = -g_y_0_xxyzz_xxxzz[k] * ab_x + g_y_0_xxyzz_xxxxzz[k];

                g_y_0_xxxyzz_xxyyy[k] = -g_y_0_xxyzz_xxyyy[k] * ab_x + g_y_0_xxyzz_xxxyyy[k];

                g_y_0_xxxyzz_xxyyz[k] = -g_y_0_xxyzz_xxyyz[k] * ab_x + g_y_0_xxyzz_xxxyyz[k];

                g_y_0_xxxyzz_xxyzz[k] = -g_y_0_xxyzz_xxyzz[k] * ab_x + g_y_0_xxyzz_xxxyzz[k];

                g_y_0_xxxyzz_xxzzz[k] = -g_y_0_xxyzz_xxzzz[k] * ab_x + g_y_0_xxyzz_xxxzzz[k];

                g_y_0_xxxyzz_xyyyy[k] = -g_y_0_xxyzz_xyyyy[k] * ab_x + g_y_0_xxyzz_xxyyyy[k];

                g_y_0_xxxyzz_xyyyz[k] = -g_y_0_xxyzz_xyyyz[k] * ab_x + g_y_0_xxyzz_xxyyyz[k];

                g_y_0_xxxyzz_xyyzz[k] = -g_y_0_xxyzz_xyyzz[k] * ab_x + g_y_0_xxyzz_xxyyzz[k];

                g_y_0_xxxyzz_xyzzz[k] = -g_y_0_xxyzz_xyzzz[k] * ab_x + g_y_0_xxyzz_xxyzzz[k];

                g_y_0_xxxyzz_xzzzz[k] = -g_y_0_xxyzz_xzzzz[k] * ab_x + g_y_0_xxyzz_xxzzzz[k];

                g_y_0_xxxyzz_yyyyy[k] = -g_y_0_xxyzz_yyyyy[k] * ab_x + g_y_0_xxyzz_xyyyyy[k];

                g_y_0_xxxyzz_yyyyz[k] = -g_y_0_xxyzz_yyyyz[k] * ab_x + g_y_0_xxyzz_xyyyyz[k];

                g_y_0_xxxyzz_yyyzz[k] = -g_y_0_xxyzz_yyyzz[k] * ab_x + g_y_0_xxyzz_xyyyzz[k];

                g_y_0_xxxyzz_yyzzz[k] = -g_y_0_xxyzz_yyzzz[k] * ab_x + g_y_0_xxyzz_xyyzzz[k];

                g_y_0_xxxyzz_yzzzz[k] = -g_y_0_xxyzz_yzzzz[k] * ab_x + g_y_0_xxyzz_xyzzzz[k];

                g_y_0_xxxyzz_zzzzz[k] = -g_y_0_xxyzz_zzzzz[k] * ab_x + g_y_0_xxyzz_xzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzzz_xxxxx, g_y_0_xxxzzz_xxxxy, g_y_0_xxxzzz_xxxxz, g_y_0_xxxzzz_xxxyy, g_y_0_xxxzzz_xxxyz, g_y_0_xxxzzz_xxxzz, g_y_0_xxxzzz_xxyyy, g_y_0_xxxzzz_xxyyz, g_y_0_xxxzzz_xxyzz, g_y_0_xxxzzz_xxzzz, g_y_0_xxxzzz_xyyyy, g_y_0_xxxzzz_xyyyz, g_y_0_xxxzzz_xyyzz, g_y_0_xxxzzz_xyzzz, g_y_0_xxxzzz_xzzzz, g_y_0_xxxzzz_yyyyy, g_y_0_xxxzzz_yyyyz, g_y_0_xxxzzz_yyyzz, g_y_0_xxxzzz_yyzzz, g_y_0_xxxzzz_yzzzz, g_y_0_xxxzzz_zzzzz, g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxxx, g_y_0_xxzzz_xxxxxy, g_y_0_xxzzz_xxxxxz, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxyy, g_y_0_xxzzz_xxxxyz, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxxzz, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyyy, g_y_0_xxzzz_xxxyyz, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxyzz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxxzzz, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyyy, g_y_0_xxzzz_xxyyyz, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyyzz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxyzzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xxzzzz, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyyy, g_y_0_xxzzz_xyyyyz, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyyzz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyyzzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xyzzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_xzzzzz, g_y_0_xxzzz_yyyyy, g_y_0_xxzzz_yyyyz, g_y_0_xxzzz_yyyzz, g_y_0_xxzzz_yyzzz, g_y_0_xxzzz_yzzzz, g_y_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxxxx[k] = -g_y_0_xxzzz_xxxxx[k] * ab_x + g_y_0_xxzzz_xxxxxx[k];

                g_y_0_xxxzzz_xxxxy[k] = -g_y_0_xxzzz_xxxxy[k] * ab_x + g_y_0_xxzzz_xxxxxy[k];

                g_y_0_xxxzzz_xxxxz[k] = -g_y_0_xxzzz_xxxxz[k] * ab_x + g_y_0_xxzzz_xxxxxz[k];

                g_y_0_xxxzzz_xxxyy[k] = -g_y_0_xxzzz_xxxyy[k] * ab_x + g_y_0_xxzzz_xxxxyy[k];

                g_y_0_xxxzzz_xxxyz[k] = -g_y_0_xxzzz_xxxyz[k] * ab_x + g_y_0_xxzzz_xxxxyz[k];

                g_y_0_xxxzzz_xxxzz[k] = -g_y_0_xxzzz_xxxzz[k] * ab_x + g_y_0_xxzzz_xxxxzz[k];

                g_y_0_xxxzzz_xxyyy[k] = -g_y_0_xxzzz_xxyyy[k] * ab_x + g_y_0_xxzzz_xxxyyy[k];

                g_y_0_xxxzzz_xxyyz[k] = -g_y_0_xxzzz_xxyyz[k] * ab_x + g_y_0_xxzzz_xxxyyz[k];

                g_y_0_xxxzzz_xxyzz[k] = -g_y_0_xxzzz_xxyzz[k] * ab_x + g_y_0_xxzzz_xxxyzz[k];

                g_y_0_xxxzzz_xxzzz[k] = -g_y_0_xxzzz_xxzzz[k] * ab_x + g_y_0_xxzzz_xxxzzz[k];

                g_y_0_xxxzzz_xyyyy[k] = -g_y_0_xxzzz_xyyyy[k] * ab_x + g_y_0_xxzzz_xxyyyy[k];

                g_y_0_xxxzzz_xyyyz[k] = -g_y_0_xxzzz_xyyyz[k] * ab_x + g_y_0_xxzzz_xxyyyz[k];

                g_y_0_xxxzzz_xyyzz[k] = -g_y_0_xxzzz_xyyzz[k] * ab_x + g_y_0_xxzzz_xxyyzz[k];

                g_y_0_xxxzzz_xyzzz[k] = -g_y_0_xxzzz_xyzzz[k] * ab_x + g_y_0_xxzzz_xxyzzz[k];

                g_y_0_xxxzzz_xzzzz[k] = -g_y_0_xxzzz_xzzzz[k] * ab_x + g_y_0_xxzzz_xxzzzz[k];

                g_y_0_xxxzzz_yyyyy[k] = -g_y_0_xxzzz_yyyyy[k] * ab_x + g_y_0_xxzzz_xyyyyy[k];

                g_y_0_xxxzzz_yyyyz[k] = -g_y_0_xxzzz_yyyyz[k] * ab_x + g_y_0_xxzzz_xyyyyz[k];

                g_y_0_xxxzzz_yyyzz[k] = -g_y_0_xxzzz_yyyzz[k] * ab_x + g_y_0_xxzzz_xyyyzz[k];

                g_y_0_xxxzzz_yyzzz[k] = -g_y_0_xxzzz_yyzzz[k] * ab_x + g_y_0_xxzzz_xyyzzz[k];

                g_y_0_xxxzzz_yzzzz[k] = -g_y_0_xxzzz_yzzzz[k] * ab_x + g_y_0_xxzzz_xyzzzz[k];

                g_y_0_xxxzzz_zzzzz[k] = -g_y_0_xxzzz_zzzzz[k] * ab_x + g_y_0_xxzzz_xzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 818 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyy_xxxxx, g_y_0_xxyyyy_xxxxy, g_y_0_xxyyyy_xxxxz, g_y_0_xxyyyy_xxxyy, g_y_0_xxyyyy_xxxyz, g_y_0_xxyyyy_xxxzz, g_y_0_xxyyyy_xxyyy, g_y_0_xxyyyy_xxyyz, g_y_0_xxyyyy_xxyzz, g_y_0_xxyyyy_xxzzz, g_y_0_xxyyyy_xyyyy, g_y_0_xxyyyy_xyyyz, g_y_0_xxyyyy_xyyzz, g_y_0_xxyyyy_xyzzz, g_y_0_xxyyyy_xzzzz, g_y_0_xxyyyy_yyyyy, g_y_0_xxyyyy_yyyyz, g_y_0_xxyyyy_yyyzz, g_y_0_xxyyyy_yyzzz, g_y_0_xxyyyy_yzzzz, g_y_0_xxyyyy_zzzzz, g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxxx, g_y_0_xyyyy_xxxxxy, g_y_0_xyyyy_xxxxxz, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxyy, g_y_0_xyyyy_xxxxyz, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxxzz, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyyy, g_y_0_xyyyy_xxxyyz, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxyzz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxxzzz, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyyy, g_y_0_xyyyy_xxyyyz, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyyzz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxyzzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xxzzzz, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyyy, g_y_0_xyyyy_xyyyyz, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyyzz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyyzzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xyzzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_xzzzzz, g_y_0_xyyyy_yyyyy, g_y_0_xyyyy_yyyyz, g_y_0_xyyyy_yyyzz, g_y_0_xyyyy_yyzzz, g_y_0_xyyyy_yzzzz, g_y_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxxxx[k] = -g_y_0_xyyyy_xxxxx[k] * ab_x + g_y_0_xyyyy_xxxxxx[k];

                g_y_0_xxyyyy_xxxxy[k] = -g_y_0_xyyyy_xxxxy[k] * ab_x + g_y_0_xyyyy_xxxxxy[k];

                g_y_0_xxyyyy_xxxxz[k] = -g_y_0_xyyyy_xxxxz[k] * ab_x + g_y_0_xyyyy_xxxxxz[k];

                g_y_0_xxyyyy_xxxyy[k] = -g_y_0_xyyyy_xxxyy[k] * ab_x + g_y_0_xyyyy_xxxxyy[k];

                g_y_0_xxyyyy_xxxyz[k] = -g_y_0_xyyyy_xxxyz[k] * ab_x + g_y_0_xyyyy_xxxxyz[k];

                g_y_0_xxyyyy_xxxzz[k] = -g_y_0_xyyyy_xxxzz[k] * ab_x + g_y_0_xyyyy_xxxxzz[k];

                g_y_0_xxyyyy_xxyyy[k] = -g_y_0_xyyyy_xxyyy[k] * ab_x + g_y_0_xyyyy_xxxyyy[k];

                g_y_0_xxyyyy_xxyyz[k] = -g_y_0_xyyyy_xxyyz[k] * ab_x + g_y_0_xyyyy_xxxyyz[k];

                g_y_0_xxyyyy_xxyzz[k] = -g_y_0_xyyyy_xxyzz[k] * ab_x + g_y_0_xyyyy_xxxyzz[k];

                g_y_0_xxyyyy_xxzzz[k] = -g_y_0_xyyyy_xxzzz[k] * ab_x + g_y_0_xyyyy_xxxzzz[k];

                g_y_0_xxyyyy_xyyyy[k] = -g_y_0_xyyyy_xyyyy[k] * ab_x + g_y_0_xyyyy_xxyyyy[k];

                g_y_0_xxyyyy_xyyyz[k] = -g_y_0_xyyyy_xyyyz[k] * ab_x + g_y_0_xyyyy_xxyyyz[k];

                g_y_0_xxyyyy_xyyzz[k] = -g_y_0_xyyyy_xyyzz[k] * ab_x + g_y_0_xyyyy_xxyyzz[k];

                g_y_0_xxyyyy_xyzzz[k] = -g_y_0_xyyyy_xyzzz[k] * ab_x + g_y_0_xyyyy_xxyzzz[k];

                g_y_0_xxyyyy_xzzzz[k] = -g_y_0_xyyyy_xzzzz[k] * ab_x + g_y_0_xyyyy_xxzzzz[k];

                g_y_0_xxyyyy_yyyyy[k] = -g_y_0_xyyyy_yyyyy[k] * ab_x + g_y_0_xyyyy_xyyyyy[k];

                g_y_0_xxyyyy_yyyyz[k] = -g_y_0_xyyyy_yyyyz[k] * ab_x + g_y_0_xyyyy_xyyyyz[k];

                g_y_0_xxyyyy_yyyzz[k] = -g_y_0_xyyyy_yyyzz[k] * ab_x + g_y_0_xyyyy_xyyyzz[k];

                g_y_0_xxyyyy_yyzzz[k] = -g_y_0_xyyyy_yyzzz[k] * ab_x + g_y_0_xyyyy_xyyzzz[k];

                g_y_0_xxyyyy_yzzzz[k] = -g_y_0_xyyyy_yzzzz[k] * ab_x + g_y_0_xyyyy_xyzzzz[k];

                g_y_0_xxyyyy_zzzzz[k] = -g_y_0_xyyyy_zzzzz[k] * ab_x + g_y_0_xyyyy_xzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyz_xxxxx, g_y_0_xxyyyz_xxxxy, g_y_0_xxyyyz_xxxxz, g_y_0_xxyyyz_xxxyy, g_y_0_xxyyyz_xxxyz, g_y_0_xxyyyz_xxxzz, g_y_0_xxyyyz_xxyyy, g_y_0_xxyyyz_xxyyz, g_y_0_xxyyyz_xxyzz, g_y_0_xxyyyz_xxzzz, g_y_0_xxyyyz_xyyyy, g_y_0_xxyyyz_xyyyz, g_y_0_xxyyyz_xyyzz, g_y_0_xxyyyz_xyzzz, g_y_0_xxyyyz_xzzzz, g_y_0_xxyyyz_yyyyy, g_y_0_xxyyyz_yyyyz, g_y_0_xxyyyz_yyyzz, g_y_0_xxyyyz_yyzzz, g_y_0_xxyyyz_yzzzz, g_y_0_xxyyyz_zzzzz, g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxxx, g_y_0_xyyyz_xxxxxy, g_y_0_xyyyz_xxxxxz, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxyy, g_y_0_xyyyz_xxxxyz, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxxzz, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyyy, g_y_0_xyyyz_xxxyyz, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxyzz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxxzzz, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyyy, g_y_0_xyyyz_xxyyyz, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyyzz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxyzzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xxzzzz, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyyy, g_y_0_xyyyz_xyyyyz, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyyzz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyyzzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xyzzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_xzzzzz, g_y_0_xyyyz_yyyyy, g_y_0_xyyyz_yyyyz, g_y_0_xyyyz_yyyzz, g_y_0_xyyyz_yyzzz, g_y_0_xyyyz_yzzzz, g_y_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxxxx[k] = -g_y_0_xyyyz_xxxxx[k] * ab_x + g_y_0_xyyyz_xxxxxx[k];

                g_y_0_xxyyyz_xxxxy[k] = -g_y_0_xyyyz_xxxxy[k] * ab_x + g_y_0_xyyyz_xxxxxy[k];

                g_y_0_xxyyyz_xxxxz[k] = -g_y_0_xyyyz_xxxxz[k] * ab_x + g_y_0_xyyyz_xxxxxz[k];

                g_y_0_xxyyyz_xxxyy[k] = -g_y_0_xyyyz_xxxyy[k] * ab_x + g_y_0_xyyyz_xxxxyy[k];

                g_y_0_xxyyyz_xxxyz[k] = -g_y_0_xyyyz_xxxyz[k] * ab_x + g_y_0_xyyyz_xxxxyz[k];

                g_y_0_xxyyyz_xxxzz[k] = -g_y_0_xyyyz_xxxzz[k] * ab_x + g_y_0_xyyyz_xxxxzz[k];

                g_y_0_xxyyyz_xxyyy[k] = -g_y_0_xyyyz_xxyyy[k] * ab_x + g_y_0_xyyyz_xxxyyy[k];

                g_y_0_xxyyyz_xxyyz[k] = -g_y_0_xyyyz_xxyyz[k] * ab_x + g_y_0_xyyyz_xxxyyz[k];

                g_y_0_xxyyyz_xxyzz[k] = -g_y_0_xyyyz_xxyzz[k] * ab_x + g_y_0_xyyyz_xxxyzz[k];

                g_y_0_xxyyyz_xxzzz[k] = -g_y_0_xyyyz_xxzzz[k] * ab_x + g_y_0_xyyyz_xxxzzz[k];

                g_y_0_xxyyyz_xyyyy[k] = -g_y_0_xyyyz_xyyyy[k] * ab_x + g_y_0_xyyyz_xxyyyy[k];

                g_y_0_xxyyyz_xyyyz[k] = -g_y_0_xyyyz_xyyyz[k] * ab_x + g_y_0_xyyyz_xxyyyz[k];

                g_y_0_xxyyyz_xyyzz[k] = -g_y_0_xyyyz_xyyzz[k] * ab_x + g_y_0_xyyyz_xxyyzz[k];

                g_y_0_xxyyyz_xyzzz[k] = -g_y_0_xyyyz_xyzzz[k] * ab_x + g_y_0_xyyyz_xxyzzz[k];

                g_y_0_xxyyyz_xzzzz[k] = -g_y_0_xyyyz_xzzzz[k] * ab_x + g_y_0_xyyyz_xxzzzz[k];

                g_y_0_xxyyyz_yyyyy[k] = -g_y_0_xyyyz_yyyyy[k] * ab_x + g_y_0_xyyyz_xyyyyy[k];

                g_y_0_xxyyyz_yyyyz[k] = -g_y_0_xyyyz_yyyyz[k] * ab_x + g_y_0_xyyyz_xyyyyz[k];

                g_y_0_xxyyyz_yyyzz[k] = -g_y_0_xyyyz_yyyzz[k] * ab_x + g_y_0_xyyyz_xyyyzz[k];

                g_y_0_xxyyyz_yyzzz[k] = -g_y_0_xyyyz_yyzzz[k] * ab_x + g_y_0_xyyyz_xyyzzz[k];

                g_y_0_xxyyyz_yzzzz[k] = -g_y_0_xyyyz_yzzzz[k] * ab_x + g_y_0_xyyyz_xyzzzz[k];

                g_y_0_xxyyyz_zzzzz[k] = -g_y_0_xyyyz_zzzzz[k] * ab_x + g_y_0_xyyyz_xzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 860 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyzz_xxxxx, g_y_0_xxyyzz_xxxxy, g_y_0_xxyyzz_xxxxz, g_y_0_xxyyzz_xxxyy, g_y_0_xxyyzz_xxxyz, g_y_0_xxyyzz_xxxzz, g_y_0_xxyyzz_xxyyy, g_y_0_xxyyzz_xxyyz, g_y_0_xxyyzz_xxyzz, g_y_0_xxyyzz_xxzzz, g_y_0_xxyyzz_xyyyy, g_y_0_xxyyzz_xyyyz, g_y_0_xxyyzz_xyyzz, g_y_0_xxyyzz_xyzzz, g_y_0_xxyyzz_xzzzz, g_y_0_xxyyzz_yyyyy, g_y_0_xxyyzz_yyyyz, g_y_0_xxyyzz_yyyzz, g_y_0_xxyyzz_yyzzz, g_y_0_xxyyzz_yzzzz, g_y_0_xxyyzz_zzzzz, g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxxx, g_y_0_xyyzz_xxxxxy, g_y_0_xyyzz_xxxxxz, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxyy, g_y_0_xyyzz_xxxxyz, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxxzz, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyyy, g_y_0_xyyzz_xxxyyz, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxyzz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxxzzz, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyyy, g_y_0_xyyzz_xxyyyz, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyyzz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxyzzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xxzzzz, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyyy, g_y_0_xyyzz_xyyyyz, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyyzz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyyzzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xyzzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_xzzzzz, g_y_0_xyyzz_yyyyy, g_y_0_xyyzz_yyyyz, g_y_0_xyyzz_yyyzz, g_y_0_xyyzz_yyzzz, g_y_0_xyyzz_yzzzz, g_y_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxxxx[k] = -g_y_0_xyyzz_xxxxx[k] * ab_x + g_y_0_xyyzz_xxxxxx[k];

                g_y_0_xxyyzz_xxxxy[k] = -g_y_0_xyyzz_xxxxy[k] * ab_x + g_y_0_xyyzz_xxxxxy[k];

                g_y_0_xxyyzz_xxxxz[k] = -g_y_0_xyyzz_xxxxz[k] * ab_x + g_y_0_xyyzz_xxxxxz[k];

                g_y_0_xxyyzz_xxxyy[k] = -g_y_0_xyyzz_xxxyy[k] * ab_x + g_y_0_xyyzz_xxxxyy[k];

                g_y_0_xxyyzz_xxxyz[k] = -g_y_0_xyyzz_xxxyz[k] * ab_x + g_y_0_xyyzz_xxxxyz[k];

                g_y_0_xxyyzz_xxxzz[k] = -g_y_0_xyyzz_xxxzz[k] * ab_x + g_y_0_xyyzz_xxxxzz[k];

                g_y_0_xxyyzz_xxyyy[k] = -g_y_0_xyyzz_xxyyy[k] * ab_x + g_y_0_xyyzz_xxxyyy[k];

                g_y_0_xxyyzz_xxyyz[k] = -g_y_0_xyyzz_xxyyz[k] * ab_x + g_y_0_xyyzz_xxxyyz[k];

                g_y_0_xxyyzz_xxyzz[k] = -g_y_0_xyyzz_xxyzz[k] * ab_x + g_y_0_xyyzz_xxxyzz[k];

                g_y_0_xxyyzz_xxzzz[k] = -g_y_0_xyyzz_xxzzz[k] * ab_x + g_y_0_xyyzz_xxxzzz[k];

                g_y_0_xxyyzz_xyyyy[k] = -g_y_0_xyyzz_xyyyy[k] * ab_x + g_y_0_xyyzz_xxyyyy[k];

                g_y_0_xxyyzz_xyyyz[k] = -g_y_0_xyyzz_xyyyz[k] * ab_x + g_y_0_xyyzz_xxyyyz[k];

                g_y_0_xxyyzz_xyyzz[k] = -g_y_0_xyyzz_xyyzz[k] * ab_x + g_y_0_xyyzz_xxyyzz[k];

                g_y_0_xxyyzz_xyzzz[k] = -g_y_0_xyyzz_xyzzz[k] * ab_x + g_y_0_xyyzz_xxyzzz[k];

                g_y_0_xxyyzz_xzzzz[k] = -g_y_0_xyyzz_xzzzz[k] * ab_x + g_y_0_xyyzz_xxzzzz[k];

                g_y_0_xxyyzz_yyyyy[k] = -g_y_0_xyyzz_yyyyy[k] * ab_x + g_y_0_xyyzz_xyyyyy[k];

                g_y_0_xxyyzz_yyyyz[k] = -g_y_0_xyyzz_yyyyz[k] * ab_x + g_y_0_xyyzz_xyyyyz[k];

                g_y_0_xxyyzz_yyyzz[k] = -g_y_0_xyyzz_yyyzz[k] * ab_x + g_y_0_xyyzz_xyyyzz[k];

                g_y_0_xxyyzz_yyzzz[k] = -g_y_0_xyyzz_yyzzz[k] * ab_x + g_y_0_xyyzz_xyyzzz[k];

                g_y_0_xxyyzz_yzzzz[k] = -g_y_0_xyyzz_yzzzz[k] * ab_x + g_y_0_xyyzz_xyzzzz[k];

                g_y_0_xxyyzz_zzzzz[k] = -g_y_0_xyyzz_zzzzz[k] * ab_x + g_y_0_xyyzz_xzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 861 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 862 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 863 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 867 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzzz_xxxxx, g_y_0_xxyzzz_xxxxy, g_y_0_xxyzzz_xxxxz, g_y_0_xxyzzz_xxxyy, g_y_0_xxyzzz_xxxyz, g_y_0_xxyzzz_xxxzz, g_y_0_xxyzzz_xxyyy, g_y_0_xxyzzz_xxyyz, g_y_0_xxyzzz_xxyzz, g_y_0_xxyzzz_xxzzz, g_y_0_xxyzzz_xyyyy, g_y_0_xxyzzz_xyyyz, g_y_0_xxyzzz_xyyzz, g_y_0_xxyzzz_xyzzz, g_y_0_xxyzzz_xzzzz, g_y_0_xxyzzz_yyyyy, g_y_0_xxyzzz_yyyyz, g_y_0_xxyzzz_yyyzz, g_y_0_xxyzzz_yyzzz, g_y_0_xxyzzz_yzzzz, g_y_0_xxyzzz_zzzzz, g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxxx, g_y_0_xyzzz_xxxxxy, g_y_0_xyzzz_xxxxxz, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxyy, g_y_0_xyzzz_xxxxyz, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxxzz, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyyy, g_y_0_xyzzz_xxxyyz, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxyzz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxxzzz, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyyy, g_y_0_xyzzz_xxyyyz, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyyzz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxyzzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xxzzzz, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyyy, g_y_0_xyzzz_xyyyyz, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyyzz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyyzzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xyzzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_xzzzzz, g_y_0_xyzzz_yyyyy, g_y_0_xyzzz_yyyyz, g_y_0_xyzzz_yyyzz, g_y_0_xyzzz_yyzzz, g_y_0_xyzzz_yzzzz, g_y_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxxxx[k] = -g_y_0_xyzzz_xxxxx[k] * ab_x + g_y_0_xyzzz_xxxxxx[k];

                g_y_0_xxyzzz_xxxxy[k] = -g_y_0_xyzzz_xxxxy[k] * ab_x + g_y_0_xyzzz_xxxxxy[k];

                g_y_0_xxyzzz_xxxxz[k] = -g_y_0_xyzzz_xxxxz[k] * ab_x + g_y_0_xyzzz_xxxxxz[k];

                g_y_0_xxyzzz_xxxyy[k] = -g_y_0_xyzzz_xxxyy[k] * ab_x + g_y_0_xyzzz_xxxxyy[k];

                g_y_0_xxyzzz_xxxyz[k] = -g_y_0_xyzzz_xxxyz[k] * ab_x + g_y_0_xyzzz_xxxxyz[k];

                g_y_0_xxyzzz_xxxzz[k] = -g_y_0_xyzzz_xxxzz[k] * ab_x + g_y_0_xyzzz_xxxxzz[k];

                g_y_0_xxyzzz_xxyyy[k] = -g_y_0_xyzzz_xxyyy[k] * ab_x + g_y_0_xyzzz_xxxyyy[k];

                g_y_0_xxyzzz_xxyyz[k] = -g_y_0_xyzzz_xxyyz[k] * ab_x + g_y_0_xyzzz_xxxyyz[k];

                g_y_0_xxyzzz_xxyzz[k] = -g_y_0_xyzzz_xxyzz[k] * ab_x + g_y_0_xyzzz_xxxyzz[k];

                g_y_0_xxyzzz_xxzzz[k] = -g_y_0_xyzzz_xxzzz[k] * ab_x + g_y_0_xyzzz_xxxzzz[k];

                g_y_0_xxyzzz_xyyyy[k] = -g_y_0_xyzzz_xyyyy[k] * ab_x + g_y_0_xyzzz_xxyyyy[k];

                g_y_0_xxyzzz_xyyyz[k] = -g_y_0_xyzzz_xyyyz[k] * ab_x + g_y_0_xyzzz_xxyyyz[k];

                g_y_0_xxyzzz_xyyzz[k] = -g_y_0_xyzzz_xyyzz[k] * ab_x + g_y_0_xyzzz_xxyyzz[k];

                g_y_0_xxyzzz_xyzzz[k] = -g_y_0_xyzzz_xyzzz[k] * ab_x + g_y_0_xyzzz_xxyzzz[k];

                g_y_0_xxyzzz_xzzzz[k] = -g_y_0_xyzzz_xzzzz[k] * ab_x + g_y_0_xyzzz_xxzzzz[k];

                g_y_0_xxyzzz_yyyyy[k] = -g_y_0_xyzzz_yyyyy[k] * ab_x + g_y_0_xyzzz_xyyyyy[k];

                g_y_0_xxyzzz_yyyyz[k] = -g_y_0_xyzzz_yyyyz[k] * ab_x + g_y_0_xyzzz_xyyyyz[k];

                g_y_0_xxyzzz_yyyzz[k] = -g_y_0_xyzzz_yyyzz[k] * ab_x + g_y_0_xyzzz_xyyyzz[k];

                g_y_0_xxyzzz_yyzzz[k] = -g_y_0_xyzzz_yyzzz[k] * ab_x + g_y_0_xyzzz_xyyzzz[k];

                g_y_0_xxyzzz_yzzzz[k] = -g_y_0_xyzzz_yzzzz[k] * ab_x + g_y_0_xyzzz_xyzzzz[k];

                g_y_0_xxyzzz_zzzzz[k] = -g_y_0_xyzzz_zzzzz[k] * ab_x + g_y_0_xyzzz_xzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 889 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 890 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 891 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 892 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 893 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 894 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 895 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 896 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 897 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 898 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 899 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 902 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzzz_xxxxx, g_y_0_xxzzzz_xxxxy, g_y_0_xxzzzz_xxxxz, g_y_0_xxzzzz_xxxyy, g_y_0_xxzzzz_xxxyz, g_y_0_xxzzzz_xxxzz, g_y_0_xxzzzz_xxyyy, g_y_0_xxzzzz_xxyyz, g_y_0_xxzzzz_xxyzz, g_y_0_xxzzzz_xxzzz, g_y_0_xxzzzz_xyyyy, g_y_0_xxzzzz_xyyyz, g_y_0_xxzzzz_xyyzz, g_y_0_xxzzzz_xyzzz, g_y_0_xxzzzz_xzzzz, g_y_0_xxzzzz_yyyyy, g_y_0_xxzzzz_yyyyz, g_y_0_xxzzzz_yyyzz, g_y_0_xxzzzz_yyzzz, g_y_0_xxzzzz_yzzzz, g_y_0_xxzzzz_zzzzz, g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxxx, g_y_0_xzzzz_xxxxxy, g_y_0_xzzzz_xxxxxz, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxyy, g_y_0_xzzzz_xxxxyz, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxxzz, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyyy, g_y_0_xzzzz_xxxyyz, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxyzz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxxzzz, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyyy, g_y_0_xzzzz_xxyyyz, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyyzz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxyzzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xxzzzz, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyyy, g_y_0_xzzzz_xyyyyz, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyyzz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyyzzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xyzzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_xzzzzz, g_y_0_xzzzz_yyyyy, g_y_0_xzzzz_yyyyz, g_y_0_xzzzz_yyyzz, g_y_0_xzzzz_yyzzz, g_y_0_xzzzz_yzzzz, g_y_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxxxx[k] = -g_y_0_xzzzz_xxxxx[k] * ab_x + g_y_0_xzzzz_xxxxxx[k];

                g_y_0_xxzzzz_xxxxy[k] = -g_y_0_xzzzz_xxxxy[k] * ab_x + g_y_0_xzzzz_xxxxxy[k];

                g_y_0_xxzzzz_xxxxz[k] = -g_y_0_xzzzz_xxxxz[k] * ab_x + g_y_0_xzzzz_xxxxxz[k];

                g_y_0_xxzzzz_xxxyy[k] = -g_y_0_xzzzz_xxxyy[k] * ab_x + g_y_0_xzzzz_xxxxyy[k];

                g_y_0_xxzzzz_xxxyz[k] = -g_y_0_xzzzz_xxxyz[k] * ab_x + g_y_0_xzzzz_xxxxyz[k];

                g_y_0_xxzzzz_xxxzz[k] = -g_y_0_xzzzz_xxxzz[k] * ab_x + g_y_0_xzzzz_xxxxzz[k];

                g_y_0_xxzzzz_xxyyy[k] = -g_y_0_xzzzz_xxyyy[k] * ab_x + g_y_0_xzzzz_xxxyyy[k];

                g_y_0_xxzzzz_xxyyz[k] = -g_y_0_xzzzz_xxyyz[k] * ab_x + g_y_0_xzzzz_xxxyyz[k];

                g_y_0_xxzzzz_xxyzz[k] = -g_y_0_xzzzz_xxyzz[k] * ab_x + g_y_0_xzzzz_xxxyzz[k];

                g_y_0_xxzzzz_xxzzz[k] = -g_y_0_xzzzz_xxzzz[k] * ab_x + g_y_0_xzzzz_xxxzzz[k];

                g_y_0_xxzzzz_xyyyy[k] = -g_y_0_xzzzz_xyyyy[k] * ab_x + g_y_0_xzzzz_xxyyyy[k];

                g_y_0_xxzzzz_xyyyz[k] = -g_y_0_xzzzz_xyyyz[k] * ab_x + g_y_0_xzzzz_xxyyyz[k];

                g_y_0_xxzzzz_xyyzz[k] = -g_y_0_xzzzz_xyyzz[k] * ab_x + g_y_0_xzzzz_xxyyzz[k];

                g_y_0_xxzzzz_xyzzz[k] = -g_y_0_xzzzz_xyzzz[k] * ab_x + g_y_0_xzzzz_xxyzzz[k];

                g_y_0_xxzzzz_xzzzz[k] = -g_y_0_xzzzz_xzzzz[k] * ab_x + g_y_0_xzzzz_xxzzzz[k];

                g_y_0_xxzzzz_yyyyy[k] = -g_y_0_xzzzz_yyyyy[k] * ab_x + g_y_0_xzzzz_xyyyyy[k];

                g_y_0_xxzzzz_yyyyz[k] = -g_y_0_xzzzz_yyyyz[k] * ab_x + g_y_0_xzzzz_xyyyyz[k];

                g_y_0_xxzzzz_yyyzz[k] = -g_y_0_xzzzz_yyyzz[k] * ab_x + g_y_0_xzzzz_xyyyzz[k];

                g_y_0_xxzzzz_yyzzz[k] = -g_y_0_xzzzz_yyzzz[k] * ab_x + g_y_0_xzzzz_xyyzzz[k];

                g_y_0_xxzzzz_yzzzz[k] = -g_y_0_xzzzz_yzzzz[k] * ab_x + g_y_0_xzzzz_xyzzzz[k];

                g_y_0_xxzzzz_zzzzz[k] = -g_y_0_xzzzz_zzzzz[k] * ab_x + g_y_0_xzzzz_xzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 917 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 918 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 919 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 920 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 921 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 922 * ccomps * dcomps);

            auto g_y_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyy_xxxxx, g_y_0_xyyyyy_xxxxy, g_y_0_xyyyyy_xxxxz, g_y_0_xyyyyy_xxxyy, g_y_0_xyyyyy_xxxyz, g_y_0_xyyyyy_xxxzz, g_y_0_xyyyyy_xxyyy, g_y_0_xyyyyy_xxyyz, g_y_0_xyyyyy_xxyzz, g_y_0_xyyyyy_xxzzz, g_y_0_xyyyyy_xyyyy, g_y_0_xyyyyy_xyyyz, g_y_0_xyyyyy_xyyzz, g_y_0_xyyyyy_xyzzz, g_y_0_xyyyyy_xzzzz, g_y_0_xyyyyy_yyyyy, g_y_0_xyyyyy_yyyyz, g_y_0_xyyyyy_yyyzz, g_y_0_xyyyyy_yyzzz, g_y_0_xyyyyy_yzzzz, g_y_0_xyyyyy_zzzzz, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxxxx[k] = -g_y_0_yyyyy_xxxxx[k] * ab_x + g_y_0_yyyyy_xxxxxx[k];

                g_y_0_xyyyyy_xxxxy[k] = -g_y_0_yyyyy_xxxxy[k] * ab_x + g_y_0_yyyyy_xxxxxy[k];

                g_y_0_xyyyyy_xxxxz[k] = -g_y_0_yyyyy_xxxxz[k] * ab_x + g_y_0_yyyyy_xxxxxz[k];

                g_y_0_xyyyyy_xxxyy[k] = -g_y_0_yyyyy_xxxyy[k] * ab_x + g_y_0_yyyyy_xxxxyy[k];

                g_y_0_xyyyyy_xxxyz[k] = -g_y_0_yyyyy_xxxyz[k] * ab_x + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_xyyyyy_xxxzz[k] = -g_y_0_yyyyy_xxxzz[k] * ab_x + g_y_0_yyyyy_xxxxzz[k];

                g_y_0_xyyyyy_xxyyy[k] = -g_y_0_yyyyy_xxyyy[k] * ab_x + g_y_0_yyyyy_xxxyyy[k];

                g_y_0_xyyyyy_xxyyz[k] = -g_y_0_yyyyy_xxyyz[k] * ab_x + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_xyyyyy_xxyzz[k] = -g_y_0_yyyyy_xxyzz[k] * ab_x + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_xyyyyy_xxzzz[k] = -g_y_0_yyyyy_xxzzz[k] * ab_x + g_y_0_yyyyy_xxxzzz[k];

                g_y_0_xyyyyy_xyyyy[k] = -g_y_0_yyyyy_xyyyy[k] * ab_x + g_y_0_yyyyy_xxyyyy[k];

                g_y_0_xyyyyy_xyyyz[k] = -g_y_0_yyyyy_xyyyz[k] * ab_x + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_xyyyyy_xyyzz[k] = -g_y_0_yyyyy_xyyzz[k] * ab_x + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_xyyyyy_xyzzz[k] = -g_y_0_yyyyy_xyzzz[k] * ab_x + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_xyyyyy_xzzzz[k] = -g_y_0_yyyyy_xzzzz[k] * ab_x + g_y_0_yyyyy_xxzzzz[k];

                g_y_0_xyyyyy_yyyyy[k] = -g_y_0_yyyyy_yyyyy[k] * ab_x + g_y_0_yyyyy_xyyyyy[k];

                g_y_0_xyyyyy_yyyyz[k] = -g_y_0_yyyyy_yyyyz[k] * ab_x + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_xyyyyy_yyyzz[k] = -g_y_0_yyyyy_yyyzz[k] * ab_x + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_xyyyyy_yyzzz[k] = -g_y_0_yyyyy_yyzzz[k] * ab_x + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_xyyyyy_yzzzz[k] = -g_y_0_yyyyy_yzzzz[k] * ab_x + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_xyyyyy_zzzzz[k] = -g_y_0_yyyyy_zzzzz[k] * ab_x + g_y_0_yyyyy_xzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 928 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 929 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 930 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 931 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 932 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 933 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 934 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 935 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyz_xxxxx, g_y_0_xyyyyz_xxxxy, g_y_0_xyyyyz_xxxxz, g_y_0_xyyyyz_xxxyy, g_y_0_xyyyyz_xxxyz, g_y_0_xyyyyz_xxxzz, g_y_0_xyyyyz_xxyyy, g_y_0_xyyyyz_xxyyz, g_y_0_xyyyyz_xxyzz, g_y_0_xyyyyz_xxzzz, g_y_0_xyyyyz_xyyyy, g_y_0_xyyyyz_xyyyz, g_y_0_xyyyyz_xyyzz, g_y_0_xyyyyz_xyzzz, g_y_0_xyyyyz_xzzzz, g_y_0_xyyyyz_yyyyy, g_y_0_xyyyyz_yyyyz, g_y_0_xyyyyz_yyyzz, g_y_0_xyyyyz_yyzzz, g_y_0_xyyyyz_yzzzz, g_y_0_xyyyyz_zzzzz, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxxxx[k] = -g_y_0_yyyyz_xxxxx[k] * ab_x + g_y_0_yyyyz_xxxxxx[k];

                g_y_0_xyyyyz_xxxxy[k] = -g_y_0_yyyyz_xxxxy[k] * ab_x + g_y_0_yyyyz_xxxxxy[k];

                g_y_0_xyyyyz_xxxxz[k] = -g_y_0_yyyyz_xxxxz[k] * ab_x + g_y_0_yyyyz_xxxxxz[k];

                g_y_0_xyyyyz_xxxyy[k] = -g_y_0_yyyyz_xxxyy[k] * ab_x + g_y_0_yyyyz_xxxxyy[k];

                g_y_0_xyyyyz_xxxyz[k] = -g_y_0_yyyyz_xxxyz[k] * ab_x + g_y_0_yyyyz_xxxxyz[k];

                g_y_0_xyyyyz_xxxzz[k] = -g_y_0_yyyyz_xxxzz[k] * ab_x + g_y_0_yyyyz_xxxxzz[k];

                g_y_0_xyyyyz_xxyyy[k] = -g_y_0_yyyyz_xxyyy[k] * ab_x + g_y_0_yyyyz_xxxyyy[k];

                g_y_0_xyyyyz_xxyyz[k] = -g_y_0_yyyyz_xxyyz[k] * ab_x + g_y_0_yyyyz_xxxyyz[k];

                g_y_0_xyyyyz_xxyzz[k] = -g_y_0_yyyyz_xxyzz[k] * ab_x + g_y_0_yyyyz_xxxyzz[k];

                g_y_0_xyyyyz_xxzzz[k] = -g_y_0_yyyyz_xxzzz[k] * ab_x + g_y_0_yyyyz_xxxzzz[k];

                g_y_0_xyyyyz_xyyyy[k] = -g_y_0_yyyyz_xyyyy[k] * ab_x + g_y_0_yyyyz_xxyyyy[k];

                g_y_0_xyyyyz_xyyyz[k] = -g_y_0_yyyyz_xyyyz[k] * ab_x + g_y_0_yyyyz_xxyyyz[k];

                g_y_0_xyyyyz_xyyzz[k] = -g_y_0_yyyyz_xyyzz[k] * ab_x + g_y_0_yyyyz_xxyyzz[k];

                g_y_0_xyyyyz_xyzzz[k] = -g_y_0_yyyyz_xyzzz[k] * ab_x + g_y_0_yyyyz_xxyzzz[k];

                g_y_0_xyyyyz_xzzzz[k] = -g_y_0_yyyyz_xzzzz[k] * ab_x + g_y_0_yyyyz_xxzzzz[k];

                g_y_0_xyyyyz_yyyyy[k] = -g_y_0_yyyyz_yyyyy[k] * ab_x + g_y_0_yyyyz_xyyyyy[k];

                g_y_0_xyyyyz_yyyyz[k] = -g_y_0_yyyyz_yyyyz[k] * ab_x + g_y_0_yyyyz_xyyyyz[k];

                g_y_0_xyyyyz_yyyzz[k] = -g_y_0_yyyyz_yyyzz[k] * ab_x + g_y_0_yyyyz_xyyyzz[k];

                g_y_0_xyyyyz_yyzzz[k] = -g_y_0_yyyyz_yyzzz[k] * ab_x + g_y_0_yyyyz_xyyzzz[k];

                g_y_0_xyyyyz_yzzzz[k] = -g_y_0_yyyyz_yzzzz[k] * ab_x + g_y_0_yyyyz_xyzzzz[k];

                g_y_0_xyyyyz_zzzzz[k] = -g_y_0_yyyyz_zzzzz[k] * ab_x + g_y_0_yyyyz_xzzzzz[k];
            }

            /// Set up 945-966 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 945 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 946 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 947 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 948 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 949 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 950 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 951 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 964 * ccomps * dcomps);

            auto g_y_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyzz_xxxxx, g_y_0_xyyyzz_xxxxy, g_y_0_xyyyzz_xxxxz, g_y_0_xyyyzz_xxxyy, g_y_0_xyyyzz_xxxyz, g_y_0_xyyyzz_xxxzz, g_y_0_xyyyzz_xxyyy, g_y_0_xyyyzz_xxyyz, g_y_0_xyyyzz_xxyzz, g_y_0_xyyyzz_xxzzz, g_y_0_xyyyzz_xyyyy, g_y_0_xyyyzz_xyyyz, g_y_0_xyyyzz_xyyzz, g_y_0_xyyyzz_xyzzz, g_y_0_xyyyzz_xzzzz, g_y_0_xyyyzz_yyyyy, g_y_0_xyyyzz_yyyyz, g_y_0_xyyyzz_yyyzz, g_y_0_xyyyzz_yyzzz, g_y_0_xyyyzz_yzzzz, g_y_0_xyyyzz_zzzzz, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxxxx[k] = -g_y_0_yyyzz_xxxxx[k] * ab_x + g_y_0_yyyzz_xxxxxx[k];

                g_y_0_xyyyzz_xxxxy[k] = -g_y_0_yyyzz_xxxxy[k] * ab_x + g_y_0_yyyzz_xxxxxy[k];

                g_y_0_xyyyzz_xxxxz[k] = -g_y_0_yyyzz_xxxxz[k] * ab_x + g_y_0_yyyzz_xxxxxz[k];

                g_y_0_xyyyzz_xxxyy[k] = -g_y_0_yyyzz_xxxyy[k] * ab_x + g_y_0_yyyzz_xxxxyy[k];

                g_y_0_xyyyzz_xxxyz[k] = -g_y_0_yyyzz_xxxyz[k] * ab_x + g_y_0_yyyzz_xxxxyz[k];

                g_y_0_xyyyzz_xxxzz[k] = -g_y_0_yyyzz_xxxzz[k] * ab_x + g_y_0_yyyzz_xxxxzz[k];

                g_y_0_xyyyzz_xxyyy[k] = -g_y_0_yyyzz_xxyyy[k] * ab_x + g_y_0_yyyzz_xxxyyy[k];

                g_y_0_xyyyzz_xxyyz[k] = -g_y_0_yyyzz_xxyyz[k] * ab_x + g_y_0_yyyzz_xxxyyz[k];

                g_y_0_xyyyzz_xxyzz[k] = -g_y_0_yyyzz_xxyzz[k] * ab_x + g_y_0_yyyzz_xxxyzz[k];

                g_y_0_xyyyzz_xxzzz[k] = -g_y_0_yyyzz_xxzzz[k] * ab_x + g_y_0_yyyzz_xxxzzz[k];

                g_y_0_xyyyzz_xyyyy[k] = -g_y_0_yyyzz_xyyyy[k] * ab_x + g_y_0_yyyzz_xxyyyy[k];

                g_y_0_xyyyzz_xyyyz[k] = -g_y_0_yyyzz_xyyyz[k] * ab_x + g_y_0_yyyzz_xxyyyz[k];

                g_y_0_xyyyzz_xyyzz[k] = -g_y_0_yyyzz_xyyzz[k] * ab_x + g_y_0_yyyzz_xxyyzz[k];

                g_y_0_xyyyzz_xyzzz[k] = -g_y_0_yyyzz_xyzzz[k] * ab_x + g_y_0_yyyzz_xxyzzz[k];

                g_y_0_xyyyzz_xzzzz[k] = -g_y_0_yyyzz_xzzzz[k] * ab_x + g_y_0_yyyzz_xxzzzz[k];

                g_y_0_xyyyzz_yyyyy[k] = -g_y_0_yyyzz_yyyyy[k] * ab_x + g_y_0_yyyzz_xyyyyy[k];

                g_y_0_xyyyzz_yyyyz[k] = -g_y_0_yyyzz_yyyyz[k] * ab_x + g_y_0_yyyzz_xyyyyz[k];

                g_y_0_xyyyzz_yyyzz[k] = -g_y_0_yyyzz_yyyzz[k] * ab_x + g_y_0_yyyzz_xyyyzz[k];

                g_y_0_xyyyzz_yyzzz[k] = -g_y_0_yyyzz_yyzzz[k] * ab_x + g_y_0_yyyzz_xyyzzz[k];

                g_y_0_xyyyzz_yzzzz[k] = -g_y_0_yyyzz_yzzzz[k] * ab_x + g_y_0_yyyzz_xyzzzz[k];

                g_y_0_xyyyzz_zzzzz[k] = -g_y_0_yyyzz_zzzzz[k] * ab_x + g_y_0_yyyzz_xzzzzz[k];
            }

            /// Set up 966-987 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 966 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 967 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 968 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 969 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 970 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 971 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 973 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 974 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 975 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 976 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 977 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 978 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 979 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 986 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzzz_xxxxx, g_y_0_xyyzzz_xxxxy, g_y_0_xyyzzz_xxxxz, g_y_0_xyyzzz_xxxyy, g_y_0_xyyzzz_xxxyz, g_y_0_xyyzzz_xxxzz, g_y_0_xyyzzz_xxyyy, g_y_0_xyyzzz_xxyyz, g_y_0_xyyzzz_xxyzz, g_y_0_xyyzzz_xxzzz, g_y_0_xyyzzz_xyyyy, g_y_0_xyyzzz_xyyyz, g_y_0_xyyzzz_xyyzz, g_y_0_xyyzzz_xyzzz, g_y_0_xyyzzz_xzzzz, g_y_0_xyyzzz_yyyyy, g_y_0_xyyzzz_yyyyz, g_y_0_xyyzzz_yyyzz, g_y_0_xyyzzz_yyzzz, g_y_0_xyyzzz_yzzzz, g_y_0_xyyzzz_zzzzz, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxxxx[k] = -g_y_0_yyzzz_xxxxx[k] * ab_x + g_y_0_yyzzz_xxxxxx[k];

                g_y_0_xyyzzz_xxxxy[k] = -g_y_0_yyzzz_xxxxy[k] * ab_x + g_y_0_yyzzz_xxxxxy[k];

                g_y_0_xyyzzz_xxxxz[k] = -g_y_0_yyzzz_xxxxz[k] * ab_x + g_y_0_yyzzz_xxxxxz[k];

                g_y_0_xyyzzz_xxxyy[k] = -g_y_0_yyzzz_xxxyy[k] * ab_x + g_y_0_yyzzz_xxxxyy[k];

                g_y_0_xyyzzz_xxxyz[k] = -g_y_0_yyzzz_xxxyz[k] * ab_x + g_y_0_yyzzz_xxxxyz[k];

                g_y_0_xyyzzz_xxxzz[k] = -g_y_0_yyzzz_xxxzz[k] * ab_x + g_y_0_yyzzz_xxxxzz[k];

                g_y_0_xyyzzz_xxyyy[k] = -g_y_0_yyzzz_xxyyy[k] * ab_x + g_y_0_yyzzz_xxxyyy[k];

                g_y_0_xyyzzz_xxyyz[k] = -g_y_0_yyzzz_xxyyz[k] * ab_x + g_y_0_yyzzz_xxxyyz[k];

                g_y_0_xyyzzz_xxyzz[k] = -g_y_0_yyzzz_xxyzz[k] * ab_x + g_y_0_yyzzz_xxxyzz[k];

                g_y_0_xyyzzz_xxzzz[k] = -g_y_0_yyzzz_xxzzz[k] * ab_x + g_y_0_yyzzz_xxxzzz[k];

                g_y_0_xyyzzz_xyyyy[k] = -g_y_0_yyzzz_xyyyy[k] * ab_x + g_y_0_yyzzz_xxyyyy[k];

                g_y_0_xyyzzz_xyyyz[k] = -g_y_0_yyzzz_xyyyz[k] * ab_x + g_y_0_yyzzz_xxyyyz[k];

                g_y_0_xyyzzz_xyyzz[k] = -g_y_0_yyzzz_xyyzz[k] * ab_x + g_y_0_yyzzz_xxyyzz[k];

                g_y_0_xyyzzz_xyzzz[k] = -g_y_0_yyzzz_xyzzz[k] * ab_x + g_y_0_yyzzz_xxyzzz[k];

                g_y_0_xyyzzz_xzzzz[k] = -g_y_0_yyzzz_xzzzz[k] * ab_x + g_y_0_yyzzz_xxzzzz[k];

                g_y_0_xyyzzz_yyyyy[k] = -g_y_0_yyzzz_yyyyy[k] * ab_x + g_y_0_yyzzz_xyyyyy[k];

                g_y_0_xyyzzz_yyyyz[k] = -g_y_0_yyzzz_yyyyz[k] * ab_x + g_y_0_yyzzz_xyyyyz[k];

                g_y_0_xyyzzz_yyyzz[k] = -g_y_0_yyzzz_yyyzz[k] * ab_x + g_y_0_yyzzz_xyyyzz[k];

                g_y_0_xyyzzz_yyzzz[k] = -g_y_0_yyzzz_yyzzz[k] * ab_x + g_y_0_yyzzz_xyyzzz[k];

                g_y_0_xyyzzz_yzzzz[k] = -g_y_0_yyzzz_yzzzz[k] * ab_x + g_y_0_yyzzz_xyzzzz[k];

                g_y_0_xyyzzz_zzzzz[k] = -g_y_0_yyzzz_zzzzz[k] * ab_x + g_y_0_yyzzz_xzzzzz[k];
            }

            /// Set up 987-1008 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1000 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1001 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1002 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1003 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1004 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1005 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1006 * ccomps * dcomps);

            auto g_y_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzzz_xxxxx, g_y_0_xyzzzz_xxxxy, g_y_0_xyzzzz_xxxxz, g_y_0_xyzzzz_xxxyy, g_y_0_xyzzzz_xxxyz, g_y_0_xyzzzz_xxxzz, g_y_0_xyzzzz_xxyyy, g_y_0_xyzzzz_xxyyz, g_y_0_xyzzzz_xxyzz, g_y_0_xyzzzz_xxzzz, g_y_0_xyzzzz_xyyyy, g_y_0_xyzzzz_xyyyz, g_y_0_xyzzzz_xyyzz, g_y_0_xyzzzz_xyzzz, g_y_0_xyzzzz_xzzzz, g_y_0_xyzzzz_yyyyy, g_y_0_xyzzzz_yyyyz, g_y_0_xyzzzz_yyyzz, g_y_0_xyzzzz_yyzzz, g_y_0_xyzzzz_yzzzz, g_y_0_xyzzzz_zzzzz, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxxxx[k] = -g_y_0_yzzzz_xxxxx[k] * ab_x + g_y_0_yzzzz_xxxxxx[k];

                g_y_0_xyzzzz_xxxxy[k] = -g_y_0_yzzzz_xxxxy[k] * ab_x + g_y_0_yzzzz_xxxxxy[k];

                g_y_0_xyzzzz_xxxxz[k] = -g_y_0_yzzzz_xxxxz[k] * ab_x + g_y_0_yzzzz_xxxxxz[k];

                g_y_0_xyzzzz_xxxyy[k] = -g_y_0_yzzzz_xxxyy[k] * ab_x + g_y_0_yzzzz_xxxxyy[k];

                g_y_0_xyzzzz_xxxyz[k] = -g_y_0_yzzzz_xxxyz[k] * ab_x + g_y_0_yzzzz_xxxxyz[k];

                g_y_0_xyzzzz_xxxzz[k] = -g_y_0_yzzzz_xxxzz[k] * ab_x + g_y_0_yzzzz_xxxxzz[k];

                g_y_0_xyzzzz_xxyyy[k] = -g_y_0_yzzzz_xxyyy[k] * ab_x + g_y_0_yzzzz_xxxyyy[k];

                g_y_0_xyzzzz_xxyyz[k] = -g_y_0_yzzzz_xxyyz[k] * ab_x + g_y_0_yzzzz_xxxyyz[k];

                g_y_0_xyzzzz_xxyzz[k] = -g_y_0_yzzzz_xxyzz[k] * ab_x + g_y_0_yzzzz_xxxyzz[k];

                g_y_0_xyzzzz_xxzzz[k] = -g_y_0_yzzzz_xxzzz[k] * ab_x + g_y_0_yzzzz_xxxzzz[k];

                g_y_0_xyzzzz_xyyyy[k] = -g_y_0_yzzzz_xyyyy[k] * ab_x + g_y_0_yzzzz_xxyyyy[k];

                g_y_0_xyzzzz_xyyyz[k] = -g_y_0_yzzzz_xyyyz[k] * ab_x + g_y_0_yzzzz_xxyyyz[k];

                g_y_0_xyzzzz_xyyzz[k] = -g_y_0_yzzzz_xyyzz[k] * ab_x + g_y_0_yzzzz_xxyyzz[k];

                g_y_0_xyzzzz_xyzzz[k] = -g_y_0_yzzzz_xyzzz[k] * ab_x + g_y_0_yzzzz_xxyzzz[k];

                g_y_0_xyzzzz_xzzzz[k] = -g_y_0_yzzzz_xzzzz[k] * ab_x + g_y_0_yzzzz_xxzzzz[k];

                g_y_0_xyzzzz_yyyyy[k] = -g_y_0_yzzzz_yyyyy[k] * ab_x + g_y_0_yzzzz_xyyyyy[k];

                g_y_0_xyzzzz_yyyyz[k] = -g_y_0_yzzzz_yyyyz[k] * ab_x + g_y_0_yzzzz_xyyyyz[k];

                g_y_0_xyzzzz_yyyzz[k] = -g_y_0_yzzzz_yyyzz[k] * ab_x + g_y_0_yzzzz_xyyyzz[k];

                g_y_0_xyzzzz_yyzzz[k] = -g_y_0_yzzzz_yyzzz[k] * ab_x + g_y_0_yzzzz_xyyzzz[k];

                g_y_0_xyzzzz_yzzzz[k] = -g_y_0_yzzzz_yzzzz[k] * ab_x + g_y_0_yzzzz_xyzzzz[k];

                g_y_0_xyzzzz_zzzzz[k] = -g_y_0_yzzzz_zzzzz[k] * ab_x + g_y_0_yzzzz_xzzzzz[k];
            }

            /// Set up 1008-1029 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1028 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzzz_xxxxx, g_y_0_xzzzzz_xxxxy, g_y_0_xzzzzz_xxxxz, g_y_0_xzzzzz_xxxyy, g_y_0_xzzzzz_xxxyz, g_y_0_xzzzzz_xxxzz, g_y_0_xzzzzz_xxyyy, g_y_0_xzzzzz_xxyyz, g_y_0_xzzzzz_xxyzz, g_y_0_xzzzzz_xxzzz, g_y_0_xzzzzz_xyyyy, g_y_0_xzzzzz_xyyyz, g_y_0_xzzzzz_xyyzz, g_y_0_xzzzzz_xyzzz, g_y_0_xzzzzz_xzzzz, g_y_0_xzzzzz_yyyyy, g_y_0_xzzzzz_yyyyz, g_y_0_xzzzzz_yyyzz, g_y_0_xzzzzz_yyzzz, g_y_0_xzzzzz_yzzzz, g_y_0_xzzzzz_zzzzz, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxxxx[k] = -g_y_0_zzzzz_xxxxx[k] * ab_x + g_y_0_zzzzz_xxxxxx[k];

                g_y_0_xzzzzz_xxxxy[k] = -g_y_0_zzzzz_xxxxy[k] * ab_x + g_y_0_zzzzz_xxxxxy[k];

                g_y_0_xzzzzz_xxxxz[k] = -g_y_0_zzzzz_xxxxz[k] * ab_x + g_y_0_zzzzz_xxxxxz[k];

                g_y_0_xzzzzz_xxxyy[k] = -g_y_0_zzzzz_xxxyy[k] * ab_x + g_y_0_zzzzz_xxxxyy[k];

                g_y_0_xzzzzz_xxxyz[k] = -g_y_0_zzzzz_xxxyz[k] * ab_x + g_y_0_zzzzz_xxxxyz[k];

                g_y_0_xzzzzz_xxxzz[k] = -g_y_0_zzzzz_xxxzz[k] * ab_x + g_y_0_zzzzz_xxxxzz[k];

                g_y_0_xzzzzz_xxyyy[k] = -g_y_0_zzzzz_xxyyy[k] * ab_x + g_y_0_zzzzz_xxxyyy[k];

                g_y_0_xzzzzz_xxyyz[k] = -g_y_0_zzzzz_xxyyz[k] * ab_x + g_y_0_zzzzz_xxxyyz[k];

                g_y_0_xzzzzz_xxyzz[k] = -g_y_0_zzzzz_xxyzz[k] * ab_x + g_y_0_zzzzz_xxxyzz[k];

                g_y_0_xzzzzz_xxzzz[k] = -g_y_0_zzzzz_xxzzz[k] * ab_x + g_y_0_zzzzz_xxxzzz[k];

                g_y_0_xzzzzz_xyyyy[k] = -g_y_0_zzzzz_xyyyy[k] * ab_x + g_y_0_zzzzz_xxyyyy[k];

                g_y_0_xzzzzz_xyyyz[k] = -g_y_0_zzzzz_xyyyz[k] * ab_x + g_y_0_zzzzz_xxyyyz[k];

                g_y_0_xzzzzz_xyyzz[k] = -g_y_0_zzzzz_xyyzz[k] * ab_x + g_y_0_zzzzz_xxyyzz[k];

                g_y_0_xzzzzz_xyzzz[k] = -g_y_0_zzzzz_xyzzz[k] * ab_x + g_y_0_zzzzz_xxyzzz[k];

                g_y_0_xzzzzz_xzzzz[k] = -g_y_0_zzzzz_xzzzz[k] * ab_x + g_y_0_zzzzz_xxzzzz[k];

                g_y_0_xzzzzz_yyyyy[k] = -g_y_0_zzzzz_yyyyy[k] * ab_x + g_y_0_zzzzz_xyyyyy[k];

                g_y_0_xzzzzz_yyyyz[k] = -g_y_0_zzzzz_yyyyz[k] * ab_x + g_y_0_zzzzz_xyyyyz[k];

                g_y_0_xzzzzz_yyyzz[k] = -g_y_0_zzzzz_yyyzz[k] * ab_x + g_y_0_zzzzz_xyyyzz[k];

                g_y_0_xzzzzz_yyzzz[k] = -g_y_0_zzzzz_yyzzz[k] * ab_x + g_y_0_zzzzz_xyyzzz[k];

                g_y_0_xzzzzz_yzzzz[k] = -g_y_0_zzzzz_yzzzz[k] * ab_x + g_y_0_zzzzz_xyzzzz[k];

                g_y_0_xzzzzz_zzzzz[k] = -g_y_0_zzzzz_zzzzz[k] * ab_x + g_y_0_zzzzz_xzzzzz[k];
            }

            /// Set up 1029-1050 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1035 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1036 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1037 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1038 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1039 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1040 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1041 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1042 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1043 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyyy_xxxxx, g_y_0_yyyyyy_xxxxy, g_y_0_yyyyyy_xxxxz, g_y_0_yyyyyy_xxxyy, g_y_0_yyyyyy_xxxyz, g_y_0_yyyyyy_xxxzz, g_y_0_yyyyyy_xxyyy, g_y_0_yyyyyy_xxyyz, g_y_0_yyyyyy_xxyzz, g_y_0_yyyyyy_xxzzz, g_y_0_yyyyyy_xyyyy, g_y_0_yyyyyy_xyyyz, g_y_0_yyyyyy_xyyzz, g_y_0_yyyyyy_xyzzz, g_y_0_yyyyyy_xzzzz, g_y_0_yyyyyy_yyyyy, g_y_0_yyyyyy_yyyyz, g_y_0_yyyyyy_yyyzz, g_y_0_yyyyyy_yyzzz, g_y_0_yyyyyy_yzzzz, g_y_0_yyyyyy_zzzzz, g_yyyyy_xxxxx, g_yyyyy_xxxxy, g_yyyyy_xxxxz, g_yyyyy_xxxyy, g_yyyyy_xxxyz, g_yyyyy_xxxzz, g_yyyyy_xxyyy, g_yyyyy_xxyyz, g_yyyyy_xxyzz, g_yyyyy_xxzzz, g_yyyyy_xyyyy, g_yyyyy_xyyyz, g_yyyyy_xyyzz, g_yyyyy_xyzzz, g_yyyyy_xzzzz, g_yyyyy_yyyyy, g_yyyyy_yyyyz, g_yyyyy_yyyzz, g_yyyyy_yyzzz, g_yyyyy_yzzzz, g_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxxxx[k] = -g_yyyyy_xxxxx[k] - g_y_0_yyyyy_xxxxx[k] * ab_y + g_y_0_yyyyy_xxxxxy[k];

                g_y_0_yyyyyy_xxxxy[k] = -g_yyyyy_xxxxy[k] - g_y_0_yyyyy_xxxxy[k] * ab_y + g_y_0_yyyyy_xxxxyy[k];

                g_y_0_yyyyyy_xxxxz[k] = -g_yyyyy_xxxxz[k] - g_y_0_yyyyy_xxxxz[k] * ab_y + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_yyyyyy_xxxyy[k] = -g_yyyyy_xxxyy[k] - g_y_0_yyyyy_xxxyy[k] * ab_y + g_y_0_yyyyy_xxxyyy[k];

                g_y_0_yyyyyy_xxxyz[k] = -g_yyyyy_xxxyz[k] - g_y_0_yyyyy_xxxyz[k] * ab_y + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_yyyyyy_xxxzz[k] = -g_yyyyy_xxxzz[k] - g_y_0_yyyyy_xxxzz[k] * ab_y + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_yyyyyy_xxyyy[k] = -g_yyyyy_xxyyy[k] - g_y_0_yyyyy_xxyyy[k] * ab_y + g_y_0_yyyyy_xxyyyy[k];

                g_y_0_yyyyyy_xxyyz[k] = -g_yyyyy_xxyyz[k] - g_y_0_yyyyy_xxyyz[k] * ab_y + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_yyyyyy_xxyzz[k] = -g_yyyyy_xxyzz[k] - g_y_0_yyyyy_xxyzz[k] * ab_y + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_yyyyyy_xxzzz[k] = -g_yyyyy_xxzzz[k] - g_y_0_yyyyy_xxzzz[k] * ab_y + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_yyyyyy_xyyyy[k] = -g_yyyyy_xyyyy[k] - g_y_0_yyyyy_xyyyy[k] * ab_y + g_y_0_yyyyy_xyyyyy[k];

                g_y_0_yyyyyy_xyyyz[k] = -g_yyyyy_xyyyz[k] - g_y_0_yyyyy_xyyyz[k] * ab_y + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_yyyyyy_xyyzz[k] = -g_yyyyy_xyyzz[k] - g_y_0_yyyyy_xyyzz[k] * ab_y + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_yyyyyy_xyzzz[k] = -g_yyyyy_xyzzz[k] - g_y_0_yyyyy_xyzzz[k] * ab_y + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_yyyyyy_xzzzz[k] = -g_yyyyy_xzzzz[k] - g_y_0_yyyyy_xzzzz[k] * ab_y + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_yyyyyy_yyyyy[k] = -g_yyyyy_yyyyy[k] - g_y_0_yyyyy_yyyyy[k] * ab_y + g_y_0_yyyyy_yyyyyy[k];

                g_y_0_yyyyyy_yyyyz[k] = -g_yyyyy_yyyyz[k] - g_y_0_yyyyy_yyyyz[k] * ab_y + g_y_0_yyyyy_yyyyyz[k];

                g_y_0_yyyyyy_yyyzz[k] = -g_yyyyy_yyyzz[k] - g_y_0_yyyyy_yyyzz[k] * ab_y + g_y_0_yyyyy_yyyyzz[k];

                g_y_0_yyyyyy_yyzzz[k] = -g_yyyyy_yyzzz[k] - g_y_0_yyyyy_yyzzz[k] * ab_y + g_y_0_yyyyy_yyyzzz[k];

                g_y_0_yyyyyy_yzzzz[k] = -g_yyyyy_yzzzz[k] - g_y_0_yyyyy_yzzzz[k] * ab_y + g_y_0_yyyyy_yyzzzz[k];

                g_y_0_yyyyyy_zzzzz[k] = -g_yyyyy_zzzzz[k] - g_y_0_yyyyy_zzzzz[k] * ab_y + g_y_0_yyyyy_yzzzzz[k];
            }

            /// Set up 1050-1071 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1057 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1063 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1070 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyy_zzzzzz, g_y_0_yyyyyz_xxxxx, g_y_0_yyyyyz_xxxxy, g_y_0_yyyyyz_xxxxz, g_y_0_yyyyyz_xxxyy, g_y_0_yyyyyz_xxxyz, g_y_0_yyyyyz_xxxzz, g_y_0_yyyyyz_xxyyy, g_y_0_yyyyyz_xxyyz, g_y_0_yyyyyz_xxyzz, g_y_0_yyyyyz_xxzzz, g_y_0_yyyyyz_xyyyy, g_y_0_yyyyyz_xyyyz, g_y_0_yyyyyz_xyyzz, g_y_0_yyyyyz_xyzzz, g_y_0_yyyyyz_xzzzz, g_y_0_yyyyyz_yyyyy, g_y_0_yyyyyz_yyyyz, g_y_0_yyyyyz_yyyzz, g_y_0_yyyyyz_yyzzz, g_y_0_yyyyyz_yzzzz, g_y_0_yyyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxxxx[k] = -g_y_0_yyyyy_xxxxx[k] * ab_z + g_y_0_yyyyy_xxxxxz[k];

                g_y_0_yyyyyz_xxxxy[k] = -g_y_0_yyyyy_xxxxy[k] * ab_z + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_yyyyyz_xxxxz[k] = -g_y_0_yyyyy_xxxxz[k] * ab_z + g_y_0_yyyyy_xxxxzz[k];

                g_y_0_yyyyyz_xxxyy[k] = -g_y_0_yyyyy_xxxyy[k] * ab_z + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_yyyyyz_xxxyz[k] = -g_y_0_yyyyy_xxxyz[k] * ab_z + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_yyyyyz_xxxzz[k] = -g_y_0_yyyyy_xxxzz[k] * ab_z + g_y_0_yyyyy_xxxzzz[k];

                g_y_0_yyyyyz_xxyyy[k] = -g_y_0_yyyyy_xxyyy[k] * ab_z + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_yyyyyz_xxyyz[k] = -g_y_0_yyyyy_xxyyz[k] * ab_z + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_yyyyyz_xxyzz[k] = -g_y_0_yyyyy_xxyzz[k] * ab_z + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_yyyyyz_xxzzz[k] = -g_y_0_yyyyy_xxzzz[k] * ab_z + g_y_0_yyyyy_xxzzzz[k];

                g_y_0_yyyyyz_xyyyy[k] = -g_y_0_yyyyy_xyyyy[k] * ab_z + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_yyyyyz_xyyyz[k] = -g_y_0_yyyyy_xyyyz[k] * ab_z + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_yyyyyz_xyyzz[k] = -g_y_0_yyyyy_xyyzz[k] * ab_z + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_yyyyyz_xyzzz[k] = -g_y_0_yyyyy_xyzzz[k] * ab_z + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_yyyyyz_xzzzz[k] = -g_y_0_yyyyy_xzzzz[k] * ab_z + g_y_0_yyyyy_xzzzzz[k];

                g_y_0_yyyyyz_yyyyy[k] = -g_y_0_yyyyy_yyyyy[k] * ab_z + g_y_0_yyyyy_yyyyyz[k];

                g_y_0_yyyyyz_yyyyz[k] = -g_y_0_yyyyy_yyyyz[k] * ab_z + g_y_0_yyyyy_yyyyzz[k];

                g_y_0_yyyyyz_yyyzz[k] = -g_y_0_yyyyy_yyyzz[k] * ab_z + g_y_0_yyyyy_yyyzzz[k];

                g_y_0_yyyyyz_yyzzz[k] = -g_y_0_yyyyy_yyzzz[k] * ab_z + g_y_0_yyyyy_yyzzzz[k];

                g_y_0_yyyyyz_yzzzz[k] = -g_y_0_yyyyy_yzzzz[k] * ab_z + g_y_0_yyyyy_yzzzzz[k];

                g_y_0_yyyyyz_zzzzz[k] = -g_y_0_yyyyy_zzzzz[k] * ab_z + g_y_0_yyyyy_zzzzzz[k];
            }

            /// Set up 1071-1092 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1072 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1073 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1074 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1075 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1076 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1077 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1078 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1079 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1080 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1081 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1082 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1083 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1084 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1085 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1086 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1087 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1088 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1089 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1090 * ccomps * dcomps);

            auto g_y_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_zzzzz, g_y_0_yyyyz_zzzzzz, g_y_0_yyyyzz_xxxxx, g_y_0_yyyyzz_xxxxy, g_y_0_yyyyzz_xxxxz, g_y_0_yyyyzz_xxxyy, g_y_0_yyyyzz_xxxyz, g_y_0_yyyyzz_xxxzz, g_y_0_yyyyzz_xxyyy, g_y_0_yyyyzz_xxyyz, g_y_0_yyyyzz_xxyzz, g_y_0_yyyyzz_xxzzz, g_y_0_yyyyzz_xyyyy, g_y_0_yyyyzz_xyyyz, g_y_0_yyyyzz_xyyzz, g_y_0_yyyyzz_xyzzz, g_y_0_yyyyzz_xzzzz, g_y_0_yyyyzz_yyyyy, g_y_0_yyyyzz_yyyyz, g_y_0_yyyyzz_yyyzz, g_y_0_yyyyzz_yyzzz, g_y_0_yyyyzz_yzzzz, g_y_0_yyyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxxxx[k] = -g_y_0_yyyyz_xxxxx[k] * ab_z + g_y_0_yyyyz_xxxxxz[k];

                g_y_0_yyyyzz_xxxxy[k] = -g_y_0_yyyyz_xxxxy[k] * ab_z + g_y_0_yyyyz_xxxxyz[k];

                g_y_0_yyyyzz_xxxxz[k] = -g_y_0_yyyyz_xxxxz[k] * ab_z + g_y_0_yyyyz_xxxxzz[k];

                g_y_0_yyyyzz_xxxyy[k] = -g_y_0_yyyyz_xxxyy[k] * ab_z + g_y_0_yyyyz_xxxyyz[k];

                g_y_0_yyyyzz_xxxyz[k] = -g_y_0_yyyyz_xxxyz[k] * ab_z + g_y_0_yyyyz_xxxyzz[k];

                g_y_0_yyyyzz_xxxzz[k] = -g_y_0_yyyyz_xxxzz[k] * ab_z + g_y_0_yyyyz_xxxzzz[k];

                g_y_0_yyyyzz_xxyyy[k] = -g_y_0_yyyyz_xxyyy[k] * ab_z + g_y_0_yyyyz_xxyyyz[k];

                g_y_0_yyyyzz_xxyyz[k] = -g_y_0_yyyyz_xxyyz[k] * ab_z + g_y_0_yyyyz_xxyyzz[k];

                g_y_0_yyyyzz_xxyzz[k] = -g_y_0_yyyyz_xxyzz[k] * ab_z + g_y_0_yyyyz_xxyzzz[k];

                g_y_0_yyyyzz_xxzzz[k] = -g_y_0_yyyyz_xxzzz[k] * ab_z + g_y_0_yyyyz_xxzzzz[k];

                g_y_0_yyyyzz_xyyyy[k] = -g_y_0_yyyyz_xyyyy[k] * ab_z + g_y_0_yyyyz_xyyyyz[k];

                g_y_0_yyyyzz_xyyyz[k] = -g_y_0_yyyyz_xyyyz[k] * ab_z + g_y_0_yyyyz_xyyyzz[k];

                g_y_0_yyyyzz_xyyzz[k] = -g_y_0_yyyyz_xyyzz[k] * ab_z + g_y_0_yyyyz_xyyzzz[k];

                g_y_0_yyyyzz_xyzzz[k] = -g_y_0_yyyyz_xyzzz[k] * ab_z + g_y_0_yyyyz_xyzzzz[k];

                g_y_0_yyyyzz_xzzzz[k] = -g_y_0_yyyyz_xzzzz[k] * ab_z + g_y_0_yyyyz_xzzzzz[k];

                g_y_0_yyyyzz_yyyyy[k] = -g_y_0_yyyyz_yyyyy[k] * ab_z + g_y_0_yyyyz_yyyyyz[k];

                g_y_0_yyyyzz_yyyyz[k] = -g_y_0_yyyyz_yyyyz[k] * ab_z + g_y_0_yyyyz_yyyyzz[k];

                g_y_0_yyyyzz_yyyzz[k] = -g_y_0_yyyyz_yyyzz[k] * ab_z + g_y_0_yyyyz_yyyzzz[k];

                g_y_0_yyyyzz_yyzzz[k] = -g_y_0_yyyyz_yyzzz[k] * ab_z + g_y_0_yyyyz_yyzzzz[k];

                g_y_0_yyyyzz_yzzzz[k] = -g_y_0_yyyyz_yzzzz[k] * ab_z + g_y_0_yyyyz_yzzzzz[k];

                g_y_0_yyyyzz_zzzzz[k] = -g_y_0_yyyyz_zzzzz[k] * ab_z + g_y_0_yyyyz_zzzzzz[k];
            }

            /// Set up 1092-1113 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1092 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1093 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1094 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1095 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1096 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1097 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1098 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1099 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1100 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1101 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1102 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1103 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1104 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1105 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1106 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1107 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1108 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1109 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1110 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1111 * ccomps * dcomps);

            auto g_y_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1112 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_zzzzz, g_y_0_yyyzz_zzzzzz, g_y_0_yyyzzz_xxxxx, g_y_0_yyyzzz_xxxxy, g_y_0_yyyzzz_xxxxz, g_y_0_yyyzzz_xxxyy, g_y_0_yyyzzz_xxxyz, g_y_0_yyyzzz_xxxzz, g_y_0_yyyzzz_xxyyy, g_y_0_yyyzzz_xxyyz, g_y_0_yyyzzz_xxyzz, g_y_0_yyyzzz_xxzzz, g_y_0_yyyzzz_xyyyy, g_y_0_yyyzzz_xyyyz, g_y_0_yyyzzz_xyyzz, g_y_0_yyyzzz_xyzzz, g_y_0_yyyzzz_xzzzz, g_y_0_yyyzzz_yyyyy, g_y_0_yyyzzz_yyyyz, g_y_0_yyyzzz_yyyzz, g_y_0_yyyzzz_yyzzz, g_y_0_yyyzzz_yzzzz, g_y_0_yyyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxxxx[k] = -g_y_0_yyyzz_xxxxx[k] * ab_z + g_y_0_yyyzz_xxxxxz[k];

                g_y_0_yyyzzz_xxxxy[k] = -g_y_0_yyyzz_xxxxy[k] * ab_z + g_y_0_yyyzz_xxxxyz[k];

                g_y_0_yyyzzz_xxxxz[k] = -g_y_0_yyyzz_xxxxz[k] * ab_z + g_y_0_yyyzz_xxxxzz[k];

                g_y_0_yyyzzz_xxxyy[k] = -g_y_0_yyyzz_xxxyy[k] * ab_z + g_y_0_yyyzz_xxxyyz[k];

                g_y_0_yyyzzz_xxxyz[k] = -g_y_0_yyyzz_xxxyz[k] * ab_z + g_y_0_yyyzz_xxxyzz[k];

                g_y_0_yyyzzz_xxxzz[k] = -g_y_0_yyyzz_xxxzz[k] * ab_z + g_y_0_yyyzz_xxxzzz[k];

                g_y_0_yyyzzz_xxyyy[k] = -g_y_0_yyyzz_xxyyy[k] * ab_z + g_y_0_yyyzz_xxyyyz[k];

                g_y_0_yyyzzz_xxyyz[k] = -g_y_0_yyyzz_xxyyz[k] * ab_z + g_y_0_yyyzz_xxyyzz[k];

                g_y_0_yyyzzz_xxyzz[k] = -g_y_0_yyyzz_xxyzz[k] * ab_z + g_y_0_yyyzz_xxyzzz[k];

                g_y_0_yyyzzz_xxzzz[k] = -g_y_0_yyyzz_xxzzz[k] * ab_z + g_y_0_yyyzz_xxzzzz[k];

                g_y_0_yyyzzz_xyyyy[k] = -g_y_0_yyyzz_xyyyy[k] * ab_z + g_y_0_yyyzz_xyyyyz[k];

                g_y_0_yyyzzz_xyyyz[k] = -g_y_0_yyyzz_xyyyz[k] * ab_z + g_y_0_yyyzz_xyyyzz[k];

                g_y_0_yyyzzz_xyyzz[k] = -g_y_0_yyyzz_xyyzz[k] * ab_z + g_y_0_yyyzz_xyyzzz[k];

                g_y_0_yyyzzz_xyzzz[k] = -g_y_0_yyyzz_xyzzz[k] * ab_z + g_y_0_yyyzz_xyzzzz[k];

                g_y_0_yyyzzz_xzzzz[k] = -g_y_0_yyyzz_xzzzz[k] * ab_z + g_y_0_yyyzz_xzzzzz[k];

                g_y_0_yyyzzz_yyyyy[k] = -g_y_0_yyyzz_yyyyy[k] * ab_z + g_y_0_yyyzz_yyyyyz[k];

                g_y_0_yyyzzz_yyyyz[k] = -g_y_0_yyyzz_yyyyz[k] * ab_z + g_y_0_yyyzz_yyyyzz[k];

                g_y_0_yyyzzz_yyyzz[k] = -g_y_0_yyyzz_yyyzz[k] * ab_z + g_y_0_yyyzz_yyyzzz[k];

                g_y_0_yyyzzz_yyzzz[k] = -g_y_0_yyyzz_yyzzz[k] * ab_z + g_y_0_yyyzz_yyzzzz[k];

                g_y_0_yyyzzz_yzzzz[k] = -g_y_0_yyyzz_yzzzz[k] * ab_z + g_y_0_yyyzz_yzzzzz[k];

                g_y_0_yyyzzz_zzzzz[k] = -g_y_0_yyyzz_zzzzz[k] * ab_z + g_y_0_yyyzz_zzzzzz[k];
            }

            /// Set up 1113-1134 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1113 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1114 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1115 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1116 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1117 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1118 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1119 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1120 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1121 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1122 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1123 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1124 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1125 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1126 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1127 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1128 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1129 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1130 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1131 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1132 * ccomps * dcomps);

            auto g_y_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_zzzzz, g_y_0_yyzzz_zzzzzz, g_y_0_yyzzzz_xxxxx, g_y_0_yyzzzz_xxxxy, g_y_0_yyzzzz_xxxxz, g_y_0_yyzzzz_xxxyy, g_y_0_yyzzzz_xxxyz, g_y_0_yyzzzz_xxxzz, g_y_0_yyzzzz_xxyyy, g_y_0_yyzzzz_xxyyz, g_y_0_yyzzzz_xxyzz, g_y_0_yyzzzz_xxzzz, g_y_0_yyzzzz_xyyyy, g_y_0_yyzzzz_xyyyz, g_y_0_yyzzzz_xyyzz, g_y_0_yyzzzz_xyzzz, g_y_0_yyzzzz_xzzzz, g_y_0_yyzzzz_yyyyy, g_y_0_yyzzzz_yyyyz, g_y_0_yyzzzz_yyyzz, g_y_0_yyzzzz_yyzzz, g_y_0_yyzzzz_yzzzz, g_y_0_yyzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxxxx[k] = -g_y_0_yyzzz_xxxxx[k] * ab_z + g_y_0_yyzzz_xxxxxz[k];

                g_y_0_yyzzzz_xxxxy[k] = -g_y_0_yyzzz_xxxxy[k] * ab_z + g_y_0_yyzzz_xxxxyz[k];

                g_y_0_yyzzzz_xxxxz[k] = -g_y_0_yyzzz_xxxxz[k] * ab_z + g_y_0_yyzzz_xxxxzz[k];

                g_y_0_yyzzzz_xxxyy[k] = -g_y_0_yyzzz_xxxyy[k] * ab_z + g_y_0_yyzzz_xxxyyz[k];

                g_y_0_yyzzzz_xxxyz[k] = -g_y_0_yyzzz_xxxyz[k] * ab_z + g_y_0_yyzzz_xxxyzz[k];

                g_y_0_yyzzzz_xxxzz[k] = -g_y_0_yyzzz_xxxzz[k] * ab_z + g_y_0_yyzzz_xxxzzz[k];

                g_y_0_yyzzzz_xxyyy[k] = -g_y_0_yyzzz_xxyyy[k] * ab_z + g_y_0_yyzzz_xxyyyz[k];

                g_y_0_yyzzzz_xxyyz[k] = -g_y_0_yyzzz_xxyyz[k] * ab_z + g_y_0_yyzzz_xxyyzz[k];

                g_y_0_yyzzzz_xxyzz[k] = -g_y_0_yyzzz_xxyzz[k] * ab_z + g_y_0_yyzzz_xxyzzz[k];

                g_y_0_yyzzzz_xxzzz[k] = -g_y_0_yyzzz_xxzzz[k] * ab_z + g_y_0_yyzzz_xxzzzz[k];

                g_y_0_yyzzzz_xyyyy[k] = -g_y_0_yyzzz_xyyyy[k] * ab_z + g_y_0_yyzzz_xyyyyz[k];

                g_y_0_yyzzzz_xyyyz[k] = -g_y_0_yyzzz_xyyyz[k] * ab_z + g_y_0_yyzzz_xyyyzz[k];

                g_y_0_yyzzzz_xyyzz[k] = -g_y_0_yyzzz_xyyzz[k] * ab_z + g_y_0_yyzzz_xyyzzz[k];

                g_y_0_yyzzzz_xyzzz[k] = -g_y_0_yyzzz_xyzzz[k] * ab_z + g_y_0_yyzzz_xyzzzz[k];

                g_y_0_yyzzzz_xzzzz[k] = -g_y_0_yyzzz_xzzzz[k] * ab_z + g_y_0_yyzzz_xzzzzz[k];

                g_y_0_yyzzzz_yyyyy[k] = -g_y_0_yyzzz_yyyyy[k] * ab_z + g_y_0_yyzzz_yyyyyz[k];

                g_y_0_yyzzzz_yyyyz[k] = -g_y_0_yyzzz_yyyyz[k] * ab_z + g_y_0_yyzzz_yyyyzz[k];

                g_y_0_yyzzzz_yyyzz[k] = -g_y_0_yyzzz_yyyzz[k] * ab_z + g_y_0_yyzzz_yyyzzz[k];

                g_y_0_yyzzzz_yyzzz[k] = -g_y_0_yyzzz_yyzzz[k] * ab_z + g_y_0_yyzzz_yyzzzz[k];

                g_y_0_yyzzzz_yzzzz[k] = -g_y_0_yyzzz_yzzzz[k] * ab_z + g_y_0_yyzzz_yzzzzz[k];

                g_y_0_yyzzzz_zzzzz[k] = -g_y_0_yyzzz_zzzzz[k] * ab_z + g_y_0_yyzzz_zzzzzz[k];
            }

            /// Set up 1134-1155 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1134 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1135 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1136 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1137 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1138 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1139 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1140 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1141 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1142 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1143 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1144 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1145 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1146 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1147 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1148 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1149 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1150 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1151 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1152 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1153 * ccomps * dcomps);

            auto g_y_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_zzzzz, g_y_0_yzzzz_zzzzzz, g_y_0_yzzzzz_xxxxx, g_y_0_yzzzzz_xxxxy, g_y_0_yzzzzz_xxxxz, g_y_0_yzzzzz_xxxyy, g_y_0_yzzzzz_xxxyz, g_y_0_yzzzzz_xxxzz, g_y_0_yzzzzz_xxyyy, g_y_0_yzzzzz_xxyyz, g_y_0_yzzzzz_xxyzz, g_y_0_yzzzzz_xxzzz, g_y_0_yzzzzz_xyyyy, g_y_0_yzzzzz_xyyyz, g_y_0_yzzzzz_xyyzz, g_y_0_yzzzzz_xyzzz, g_y_0_yzzzzz_xzzzz, g_y_0_yzzzzz_yyyyy, g_y_0_yzzzzz_yyyyz, g_y_0_yzzzzz_yyyzz, g_y_0_yzzzzz_yyzzz, g_y_0_yzzzzz_yzzzz, g_y_0_yzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxxxx[k] = -g_y_0_yzzzz_xxxxx[k] * ab_z + g_y_0_yzzzz_xxxxxz[k];

                g_y_0_yzzzzz_xxxxy[k] = -g_y_0_yzzzz_xxxxy[k] * ab_z + g_y_0_yzzzz_xxxxyz[k];

                g_y_0_yzzzzz_xxxxz[k] = -g_y_0_yzzzz_xxxxz[k] * ab_z + g_y_0_yzzzz_xxxxzz[k];

                g_y_0_yzzzzz_xxxyy[k] = -g_y_0_yzzzz_xxxyy[k] * ab_z + g_y_0_yzzzz_xxxyyz[k];

                g_y_0_yzzzzz_xxxyz[k] = -g_y_0_yzzzz_xxxyz[k] * ab_z + g_y_0_yzzzz_xxxyzz[k];

                g_y_0_yzzzzz_xxxzz[k] = -g_y_0_yzzzz_xxxzz[k] * ab_z + g_y_0_yzzzz_xxxzzz[k];

                g_y_0_yzzzzz_xxyyy[k] = -g_y_0_yzzzz_xxyyy[k] * ab_z + g_y_0_yzzzz_xxyyyz[k];

                g_y_0_yzzzzz_xxyyz[k] = -g_y_0_yzzzz_xxyyz[k] * ab_z + g_y_0_yzzzz_xxyyzz[k];

                g_y_0_yzzzzz_xxyzz[k] = -g_y_0_yzzzz_xxyzz[k] * ab_z + g_y_0_yzzzz_xxyzzz[k];

                g_y_0_yzzzzz_xxzzz[k] = -g_y_0_yzzzz_xxzzz[k] * ab_z + g_y_0_yzzzz_xxzzzz[k];

                g_y_0_yzzzzz_xyyyy[k] = -g_y_0_yzzzz_xyyyy[k] * ab_z + g_y_0_yzzzz_xyyyyz[k];

                g_y_0_yzzzzz_xyyyz[k] = -g_y_0_yzzzz_xyyyz[k] * ab_z + g_y_0_yzzzz_xyyyzz[k];

                g_y_0_yzzzzz_xyyzz[k] = -g_y_0_yzzzz_xyyzz[k] * ab_z + g_y_0_yzzzz_xyyzzz[k];

                g_y_0_yzzzzz_xyzzz[k] = -g_y_0_yzzzz_xyzzz[k] * ab_z + g_y_0_yzzzz_xyzzzz[k];

                g_y_0_yzzzzz_xzzzz[k] = -g_y_0_yzzzz_xzzzz[k] * ab_z + g_y_0_yzzzz_xzzzzz[k];

                g_y_0_yzzzzz_yyyyy[k] = -g_y_0_yzzzz_yyyyy[k] * ab_z + g_y_0_yzzzz_yyyyyz[k];

                g_y_0_yzzzzz_yyyyz[k] = -g_y_0_yzzzz_yyyyz[k] * ab_z + g_y_0_yzzzz_yyyyzz[k];

                g_y_0_yzzzzz_yyyzz[k] = -g_y_0_yzzzz_yyyzz[k] * ab_z + g_y_0_yzzzz_yyyzzz[k];

                g_y_0_yzzzzz_yyzzz[k] = -g_y_0_yzzzz_yyzzz[k] * ab_z + g_y_0_yzzzz_yyzzzz[k];

                g_y_0_yzzzzz_yzzzz[k] = -g_y_0_yzzzz_yzzzz[k] * ab_z + g_y_0_yzzzz_yzzzzz[k];

                g_y_0_yzzzzz_zzzzz[k] = -g_y_0_yzzzz_zzzzz[k] * ab_z + g_y_0_yzzzz_zzzzzz[k];
            }

            /// Set up 1155-1176 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1155 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1156 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1157 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1158 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1159 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1160 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1161 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1162 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1163 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1164 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1165 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1166 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1167 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1168 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1169 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1170 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1171 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1172 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1173 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1174 * ccomps * dcomps);

            auto g_y_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_zzzzz, g_y_0_zzzzz_zzzzzz, g_y_0_zzzzzz_xxxxx, g_y_0_zzzzzz_xxxxy, g_y_0_zzzzzz_xxxxz, g_y_0_zzzzzz_xxxyy, g_y_0_zzzzzz_xxxyz, g_y_0_zzzzzz_xxxzz, g_y_0_zzzzzz_xxyyy, g_y_0_zzzzzz_xxyyz, g_y_0_zzzzzz_xxyzz, g_y_0_zzzzzz_xxzzz, g_y_0_zzzzzz_xyyyy, g_y_0_zzzzzz_xyyyz, g_y_0_zzzzzz_xyyzz, g_y_0_zzzzzz_xyzzz, g_y_0_zzzzzz_xzzzz, g_y_0_zzzzzz_yyyyy, g_y_0_zzzzzz_yyyyz, g_y_0_zzzzzz_yyyzz, g_y_0_zzzzzz_yyzzz, g_y_0_zzzzzz_yzzzz, g_y_0_zzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxxxx[k] = -g_y_0_zzzzz_xxxxx[k] * ab_z + g_y_0_zzzzz_xxxxxz[k];

                g_y_0_zzzzzz_xxxxy[k] = -g_y_0_zzzzz_xxxxy[k] * ab_z + g_y_0_zzzzz_xxxxyz[k];

                g_y_0_zzzzzz_xxxxz[k] = -g_y_0_zzzzz_xxxxz[k] * ab_z + g_y_0_zzzzz_xxxxzz[k];

                g_y_0_zzzzzz_xxxyy[k] = -g_y_0_zzzzz_xxxyy[k] * ab_z + g_y_0_zzzzz_xxxyyz[k];

                g_y_0_zzzzzz_xxxyz[k] = -g_y_0_zzzzz_xxxyz[k] * ab_z + g_y_0_zzzzz_xxxyzz[k];

                g_y_0_zzzzzz_xxxzz[k] = -g_y_0_zzzzz_xxxzz[k] * ab_z + g_y_0_zzzzz_xxxzzz[k];

                g_y_0_zzzzzz_xxyyy[k] = -g_y_0_zzzzz_xxyyy[k] * ab_z + g_y_0_zzzzz_xxyyyz[k];

                g_y_0_zzzzzz_xxyyz[k] = -g_y_0_zzzzz_xxyyz[k] * ab_z + g_y_0_zzzzz_xxyyzz[k];

                g_y_0_zzzzzz_xxyzz[k] = -g_y_0_zzzzz_xxyzz[k] * ab_z + g_y_0_zzzzz_xxyzzz[k];

                g_y_0_zzzzzz_xxzzz[k] = -g_y_0_zzzzz_xxzzz[k] * ab_z + g_y_0_zzzzz_xxzzzz[k];

                g_y_0_zzzzzz_xyyyy[k] = -g_y_0_zzzzz_xyyyy[k] * ab_z + g_y_0_zzzzz_xyyyyz[k];

                g_y_0_zzzzzz_xyyyz[k] = -g_y_0_zzzzz_xyyyz[k] * ab_z + g_y_0_zzzzz_xyyyzz[k];

                g_y_0_zzzzzz_xyyzz[k] = -g_y_0_zzzzz_xyyzz[k] * ab_z + g_y_0_zzzzz_xyyzzz[k];

                g_y_0_zzzzzz_xyzzz[k] = -g_y_0_zzzzz_xyzzz[k] * ab_z + g_y_0_zzzzz_xyzzzz[k];

                g_y_0_zzzzzz_xzzzz[k] = -g_y_0_zzzzz_xzzzz[k] * ab_z + g_y_0_zzzzz_xzzzzz[k];

                g_y_0_zzzzzz_yyyyy[k] = -g_y_0_zzzzz_yyyyy[k] * ab_z + g_y_0_zzzzz_yyyyyz[k];

                g_y_0_zzzzzz_yyyyz[k] = -g_y_0_zzzzz_yyyyz[k] * ab_z + g_y_0_zzzzz_yyyyzz[k];

                g_y_0_zzzzzz_yyyzz[k] = -g_y_0_zzzzz_yyyzz[k] * ab_z + g_y_0_zzzzz_yyyzzz[k];

                g_y_0_zzzzzz_yyzzz[k] = -g_y_0_zzzzz_yyzzz[k] * ab_z + g_y_0_zzzzz_yyzzzz[k];

                g_y_0_zzzzzz_yzzzz[k] = -g_y_0_zzzzz_yzzzz[k] * ab_z + g_y_0_zzzzz_yzzzzz[k];

                g_y_0_zzzzzz_zzzzz[k] = -g_y_0_zzzzz_zzzzz[k] * ab_z + g_y_0_zzzzz_zzzzzz[k];
            }

            /// Set up 1176-1197 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 1196 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxxx, g_z_0_xxxxx_xxxxxy, g_z_0_xxxxx_xxxxxz, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxyy, g_z_0_xxxxx_xxxxyz, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxxzz, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyyy, g_z_0_xxxxx_xxxyyz, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxyzz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxxzzz, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyyy, g_z_0_xxxxx_xxyyyz, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyyzz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxyzzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xxzzzz, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyyy, g_z_0_xxxxx_xyyyyz, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyyzz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyyzzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xyzzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_xzzzzz, g_z_0_xxxxx_yyyyy, g_z_0_xxxxx_yyyyz, g_z_0_xxxxx_yyyzz, g_z_0_xxxxx_yyzzz, g_z_0_xxxxx_yzzzz, g_z_0_xxxxx_zzzzz, g_z_0_xxxxxx_xxxxx, g_z_0_xxxxxx_xxxxy, g_z_0_xxxxxx_xxxxz, g_z_0_xxxxxx_xxxyy, g_z_0_xxxxxx_xxxyz, g_z_0_xxxxxx_xxxzz, g_z_0_xxxxxx_xxyyy, g_z_0_xxxxxx_xxyyz, g_z_0_xxxxxx_xxyzz, g_z_0_xxxxxx_xxzzz, g_z_0_xxxxxx_xyyyy, g_z_0_xxxxxx_xyyyz, g_z_0_xxxxxx_xyyzz, g_z_0_xxxxxx_xyzzz, g_z_0_xxxxxx_xzzzz, g_z_0_xxxxxx_yyyyy, g_z_0_xxxxxx_yyyyz, g_z_0_xxxxxx_yyyzz, g_z_0_xxxxxx_yyzzz, g_z_0_xxxxxx_yzzzz, g_z_0_xxxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxxxx[k] = -g_z_0_xxxxx_xxxxx[k] * ab_x + g_z_0_xxxxx_xxxxxx[k];

                g_z_0_xxxxxx_xxxxy[k] = -g_z_0_xxxxx_xxxxy[k] * ab_x + g_z_0_xxxxx_xxxxxy[k];

                g_z_0_xxxxxx_xxxxz[k] = -g_z_0_xxxxx_xxxxz[k] * ab_x + g_z_0_xxxxx_xxxxxz[k];

                g_z_0_xxxxxx_xxxyy[k] = -g_z_0_xxxxx_xxxyy[k] * ab_x + g_z_0_xxxxx_xxxxyy[k];

                g_z_0_xxxxxx_xxxyz[k] = -g_z_0_xxxxx_xxxyz[k] * ab_x + g_z_0_xxxxx_xxxxyz[k];

                g_z_0_xxxxxx_xxxzz[k] = -g_z_0_xxxxx_xxxzz[k] * ab_x + g_z_0_xxxxx_xxxxzz[k];

                g_z_0_xxxxxx_xxyyy[k] = -g_z_0_xxxxx_xxyyy[k] * ab_x + g_z_0_xxxxx_xxxyyy[k];

                g_z_0_xxxxxx_xxyyz[k] = -g_z_0_xxxxx_xxyyz[k] * ab_x + g_z_0_xxxxx_xxxyyz[k];

                g_z_0_xxxxxx_xxyzz[k] = -g_z_0_xxxxx_xxyzz[k] * ab_x + g_z_0_xxxxx_xxxyzz[k];

                g_z_0_xxxxxx_xxzzz[k] = -g_z_0_xxxxx_xxzzz[k] * ab_x + g_z_0_xxxxx_xxxzzz[k];

                g_z_0_xxxxxx_xyyyy[k] = -g_z_0_xxxxx_xyyyy[k] * ab_x + g_z_0_xxxxx_xxyyyy[k];

                g_z_0_xxxxxx_xyyyz[k] = -g_z_0_xxxxx_xyyyz[k] * ab_x + g_z_0_xxxxx_xxyyyz[k];

                g_z_0_xxxxxx_xyyzz[k] = -g_z_0_xxxxx_xyyzz[k] * ab_x + g_z_0_xxxxx_xxyyzz[k];

                g_z_0_xxxxxx_xyzzz[k] = -g_z_0_xxxxx_xyzzz[k] * ab_x + g_z_0_xxxxx_xxyzzz[k];

                g_z_0_xxxxxx_xzzzz[k] = -g_z_0_xxxxx_xzzzz[k] * ab_x + g_z_0_xxxxx_xxzzzz[k];

                g_z_0_xxxxxx_yyyyy[k] = -g_z_0_xxxxx_yyyyy[k] * ab_x + g_z_0_xxxxx_xyyyyy[k];

                g_z_0_xxxxxx_yyyyz[k] = -g_z_0_xxxxx_yyyyz[k] * ab_x + g_z_0_xxxxx_xyyyyz[k];

                g_z_0_xxxxxx_yyyzz[k] = -g_z_0_xxxxx_yyyzz[k] * ab_x + g_z_0_xxxxx_xyyyzz[k];

                g_z_0_xxxxxx_yyzzz[k] = -g_z_0_xxxxx_yyzzz[k] * ab_x + g_z_0_xxxxx_xyyzzz[k];

                g_z_0_xxxxxx_yzzzz[k] = -g_z_0_xxxxx_yzzzz[k] * ab_x + g_z_0_xxxxx_xyzzzz[k];

                g_z_0_xxxxxx_zzzzz[k] = -g_z_0_xxxxx_zzzzz[k] * ab_x + g_z_0_xxxxx_xzzzzz[k];
            }

            /// Set up 1197-1218 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 1217 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxy_xxxxx, g_z_0_xxxxxy_xxxxy, g_z_0_xxxxxy_xxxxz, g_z_0_xxxxxy_xxxyy, g_z_0_xxxxxy_xxxyz, g_z_0_xxxxxy_xxxzz, g_z_0_xxxxxy_xxyyy, g_z_0_xxxxxy_xxyyz, g_z_0_xxxxxy_xxyzz, g_z_0_xxxxxy_xxzzz, g_z_0_xxxxxy_xyyyy, g_z_0_xxxxxy_xyyyz, g_z_0_xxxxxy_xyyzz, g_z_0_xxxxxy_xyzzz, g_z_0_xxxxxy_xzzzz, g_z_0_xxxxxy_yyyyy, g_z_0_xxxxxy_yyyyz, g_z_0_xxxxxy_yyyzz, g_z_0_xxxxxy_yyzzz, g_z_0_xxxxxy_yzzzz, g_z_0_xxxxxy_zzzzz, g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxxx, g_z_0_xxxxy_xxxxxy, g_z_0_xxxxy_xxxxxz, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxyy, g_z_0_xxxxy_xxxxyz, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxxzz, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyyy, g_z_0_xxxxy_xxxyyz, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxyzz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxxzzz, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyyy, g_z_0_xxxxy_xxyyyz, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyyzz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxyzzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xxzzzz, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyyy, g_z_0_xxxxy_xyyyyz, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyyzz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyyzzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xyzzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_xzzzzz, g_z_0_xxxxy_yyyyy, g_z_0_xxxxy_yyyyz, g_z_0_xxxxy_yyyzz, g_z_0_xxxxy_yyzzz, g_z_0_xxxxy_yzzzz, g_z_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxxxx[k] = -g_z_0_xxxxy_xxxxx[k] * ab_x + g_z_0_xxxxy_xxxxxx[k];

                g_z_0_xxxxxy_xxxxy[k] = -g_z_0_xxxxy_xxxxy[k] * ab_x + g_z_0_xxxxy_xxxxxy[k];

                g_z_0_xxxxxy_xxxxz[k] = -g_z_0_xxxxy_xxxxz[k] * ab_x + g_z_0_xxxxy_xxxxxz[k];

                g_z_0_xxxxxy_xxxyy[k] = -g_z_0_xxxxy_xxxyy[k] * ab_x + g_z_0_xxxxy_xxxxyy[k];

                g_z_0_xxxxxy_xxxyz[k] = -g_z_0_xxxxy_xxxyz[k] * ab_x + g_z_0_xxxxy_xxxxyz[k];

                g_z_0_xxxxxy_xxxzz[k] = -g_z_0_xxxxy_xxxzz[k] * ab_x + g_z_0_xxxxy_xxxxzz[k];

                g_z_0_xxxxxy_xxyyy[k] = -g_z_0_xxxxy_xxyyy[k] * ab_x + g_z_0_xxxxy_xxxyyy[k];

                g_z_0_xxxxxy_xxyyz[k] = -g_z_0_xxxxy_xxyyz[k] * ab_x + g_z_0_xxxxy_xxxyyz[k];

                g_z_0_xxxxxy_xxyzz[k] = -g_z_0_xxxxy_xxyzz[k] * ab_x + g_z_0_xxxxy_xxxyzz[k];

                g_z_0_xxxxxy_xxzzz[k] = -g_z_0_xxxxy_xxzzz[k] * ab_x + g_z_0_xxxxy_xxxzzz[k];

                g_z_0_xxxxxy_xyyyy[k] = -g_z_0_xxxxy_xyyyy[k] * ab_x + g_z_0_xxxxy_xxyyyy[k];

                g_z_0_xxxxxy_xyyyz[k] = -g_z_0_xxxxy_xyyyz[k] * ab_x + g_z_0_xxxxy_xxyyyz[k];

                g_z_0_xxxxxy_xyyzz[k] = -g_z_0_xxxxy_xyyzz[k] * ab_x + g_z_0_xxxxy_xxyyzz[k];

                g_z_0_xxxxxy_xyzzz[k] = -g_z_0_xxxxy_xyzzz[k] * ab_x + g_z_0_xxxxy_xxyzzz[k];

                g_z_0_xxxxxy_xzzzz[k] = -g_z_0_xxxxy_xzzzz[k] * ab_x + g_z_0_xxxxy_xxzzzz[k];

                g_z_0_xxxxxy_yyyyy[k] = -g_z_0_xxxxy_yyyyy[k] * ab_x + g_z_0_xxxxy_xyyyyy[k];

                g_z_0_xxxxxy_yyyyz[k] = -g_z_0_xxxxy_yyyyz[k] * ab_x + g_z_0_xxxxy_xyyyyz[k];

                g_z_0_xxxxxy_yyyzz[k] = -g_z_0_xxxxy_yyyzz[k] * ab_x + g_z_0_xxxxy_xyyyzz[k];

                g_z_0_xxxxxy_yyzzz[k] = -g_z_0_xxxxy_yyzzz[k] * ab_x + g_z_0_xxxxy_xyyzzz[k];

                g_z_0_xxxxxy_yzzzz[k] = -g_z_0_xxxxy_yzzzz[k] * ab_x + g_z_0_xxxxy_xyzzzz[k];

                g_z_0_xxxxxy_zzzzz[k] = -g_z_0_xxxxy_zzzzz[k] * ab_x + g_z_0_xxxxy_xzzzzz[k];
            }

            /// Set up 1218-1239 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 1238 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxz_xxxxx, g_z_0_xxxxxz_xxxxy, g_z_0_xxxxxz_xxxxz, g_z_0_xxxxxz_xxxyy, g_z_0_xxxxxz_xxxyz, g_z_0_xxxxxz_xxxzz, g_z_0_xxxxxz_xxyyy, g_z_0_xxxxxz_xxyyz, g_z_0_xxxxxz_xxyzz, g_z_0_xxxxxz_xxzzz, g_z_0_xxxxxz_xyyyy, g_z_0_xxxxxz_xyyyz, g_z_0_xxxxxz_xyyzz, g_z_0_xxxxxz_xyzzz, g_z_0_xxxxxz_xzzzz, g_z_0_xxxxxz_yyyyy, g_z_0_xxxxxz_yyyyz, g_z_0_xxxxxz_yyyzz, g_z_0_xxxxxz_yyzzz, g_z_0_xxxxxz_yzzzz, g_z_0_xxxxxz_zzzzz, g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxxx, g_z_0_xxxxz_xxxxxy, g_z_0_xxxxz_xxxxxz, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxyy, g_z_0_xxxxz_xxxxyz, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxxzz, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyyy, g_z_0_xxxxz_xxxyyz, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxyzz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxxzzz, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyyy, g_z_0_xxxxz_xxyyyz, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyyzz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxyzzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xxzzzz, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyyy, g_z_0_xxxxz_xyyyyz, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyyzz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyyzzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xyzzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_xzzzzz, g_z_0_xxxxz_yyyyy, g_z_0_xxxxz_yyyyz, g_z_0_xxxxz_yyyzz, g_z_0_xxxxz_yyzzz, g_z_0_xxxxz_yzzzz, g_z_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxxxx[k] = -g_z_0_xxxxz_xxxxx[k] * ab_x + g_z_0_xxxxz_xxxxxx[k];

                g_z_0_xxxxxz_xxxxy[k] = -g_z_0_xxxxz_xxxxy[k] * ab_x + g_z_0_xxxxz_xxxxxy[k];

                g_z_0_xxxxxz_xxxxz[k] = -g_z_0_xxxxz_xxxxz[k] * ab_x + g_z_0_xxxxz_xxxxxz[k];

                g_z_0_xxxxxz_xxxyy[k] = -g_z_0_xxxxz_xxxyy[k] * ab_x + g_z_0_xxxxz_xxxxyy[k];

                g_z_0_xxxxxz_xxxyz[k] = -g_z_0_xxxxz_xxxyz[k] * ab_x + g_z_0_xxxxz_xxxxyz[k];

                g_z_0_xxxxxz_xxxzz[k] = -g_z_0_xxxxz_xxxzz[k] * ab_x + g_z_0_xxxxz_xxxxzz[k];

                g_z_0_xxxxxz_xxyyy[k] = -g_z_0_xxxxz_xxyyy[k] * ab_x + g_z_0_xxxxz_xxxyyy[k];

                g_z_0_xxxxxz_xxyyz[k] = -g_z_0_xxxxz_xxyyz[k] * ab_x + g_z_0_xxxxz_xxxyyz[k];

                g_z_0_xxxxxz_xxyzz[k] = -g_z_0_xxxxz_xxyzz[k] * ab_x + g_z_0_xxxxz_xxxyzz[k];

                g_z_0_xxxxxz_xxzzz[k] = -g_z_0_xxxxz_xxzzz[k] * ab_x + g_z_0_xxxxz_xxxzzz[k];

                g_z_0_xxxxxz_xyyyy[k] = -g_z_0_xxxxz_xyyyy[k] * ab_x + g_z_0_xxxxz_xxyyyy[k];

                g_z_0_xxxxxz_xyyyz[k] = -g_z_0_xxxxz_xyyyz[k] * ab_x + g_z_0_xxxxz_xxyyyz[k];

                g_z_0_xxxxxz_xyyzz[k] = -g_z_0_xxxxz_xyyzz[k] * ab_x + g_z_0_xxxxz_xxyyzz[k];

                g_z_0_xxxxxz_xyzzz[k] = -g_z_0_xxxxz_xyzzz[k] * ab_x + g_z_0_xxxxz_xxyzzz[k];

                g_z_0_xxxxxz_xzzzz[k] = -g_z_0_xxxxz_xzzzz[k] * ab_x + g_z_0_xxxxz_xxzzzz[k];

                g_z_0_xxxxxz_yyyyy[k] = -g_z_0_xxxxz_yyyyy[k] * ab_x + g_z_0_xxxxz_xyyyyy[k];

                g_z_0_xxxxxz_yyyyz[k] = -g_z_0_xxxxz_yyyyz[k] * ab_x + g_z_0_xxxxz_xyyyyz[k];

                g_z_0_xxxxxz_yyyzz[k] = -g_z_0_xxxxz_yyyzz[k] * ab_x + g_z_0_xxxxz_xyyyzz[k];

                g_z_0_xxxxxz_yyzzz[k] = -g_z_0_xxxxz_yyzzz[k] * ab_x + g_z_0_xxxxz_xyyzzz[k];

                g_z_0_xxxxxz_yzzzz[k] = -g_z_0_xxxxz_yzzzz[k] * ab_x + g_z_0_xxxxz_xyzzzz[k];

                g_z_0_xxxxxz_zzzzz[k] = -g_z_0_xxxxz_zzzzz[k] * ab_x + g_z_0_xxxxz_xzzzzz[k];
            }

            /// Set up 1239-1260 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyy_xxxxx, g_z_0_xxxxyy_xxxxy, g_z_0_xxxxyy_xxxxz, g_z_0_xxxxyy_xxxyy, g_z_0_xxxxyy_xxxyz, g_z_0_xxxxyy_xxxzz, g_z_0_xxxxyy_xxyyy, g_z_0_xxxxyy_xxyyz, g_z_0_xxxxyy_xxyzz, g_z_0_xxxxyy_xxzzz, g_z_0_xxxxyy_xyyyy, g_z_0_xxxxyy_xyyyz, g_z_0_xxxxyy_xyyzz, g_z_0_xxxxyy_xyzzz, g_z_0_xxxxyy_xzzzz, g_z_0_xxxxyy_yyyyy, g_z_0_xxxxyy_yyyyz, g_z_0_xxxxyy_yyyzz, g_z_0_xxxxyy_yyzzz, g_z_0_xxxxyy_yzzzz, g_z_0_xxxxyy_zzzzz, g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxxx, g_z_0_xxxyy_xxxxxy, g_z_0_xxxyy_xxxxxz, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxyy, g_z_0_xxxyy_xxxxyz, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxxzz, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyyy, g_z_0_xxxyy_xxxyyz, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxyzz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxxzzz, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyyy, g_z_0_xxxyy_xxyyyz, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyyzz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxyzzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xxzzzz, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyyy, g_z_0_xxxyy_xyyyyz, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyyzz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyyzzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xyzzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_xzzzzz, g_z_0_xxxyy_yyyyy, g_z_0_xxxyy_yyyyz, g_z_0_xxxyy_yyyzz, g_z_0_xxxyy_yyzzz, g_z_0_xxxyy_yzzzz, g_z_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxxxx[k] = -g_z_0_xxxyy_xxxxx[k] * ab_x + g_z_0_xxxyy_xxxxxx[k];

                g_z_0_xxxxyy_xxxxy[k] = -g_z_0_xxxyy_xxxxy[k] * ab_x + g_z_0_xxxyy_xxxxxy[k];

                g_z_0_xxxxyy_xxxxz[k] = -g_z_0_xxxyy_xxxxz[k] * ab_x + g_z_0_xxxyy_xxxxxz[k];

                g_z_0_xxxxyy_xxxyy[k] = -g_z_0_xxxyy_xxxyy[k] * ab_x + g_z_0_xxxyy_xxxxyy[k];

                g_z_0_xxxxyy_xxxyz[k] = -g_z_0_xxxyy_xxxyz[k] * ab_x + g_z_0_xxxyy_xxxxyz[k];

                g_z_0_xxxxyy_xxxzz[k] = -g_z_0_xxxyy_xxxzz[k] * ab_x + g_z_0_xxxyy_xxxxzz[k];

                g_z_0_xxxxyy_xxyyy[k] = -g_z_0_xxxyy_xxyyy[k] * ab_x + g_z_0_xxxyy_xxxyyy[k];

                g_z_0_xxxxyy_xxyyz[k] = -g_z_0_xxxyy_xxyyz[k] * ab_x + g_z_0_xxxyy_xxxyyz[k];

                g_z_0_xxxxyy_xxyzz[k] = -g_z_0_xxxyy_xxyzz[k] * ab_x + g_z_0_xxxyy_xxxyzz[k];

                g_z_0_xxxxyy_xxzzz[k] = -g_z_0_xxxyy_xxzzz[k] * ab_x + g_z_0_xxxyy_xxxzzz[k];

                g_z_0_xxxxyy_xyyyy[k] = -g_z_0_xxxyy_xyyyy[k] * ab_x + g_z_0_xxxyy_xxyyyy[k];

                g_z_0_xxxxyy_xyyyz[k] = -g_z_0_xxxyy_xyyyz[k] * ab_x + g_z_0_xxxyy_xxyyyz[k];

                g_z_0_xxxxyy_xyyzz[k] = -g_z_0_xxxyy_xyyzz[k] * ab_x + g_z_0_xxxyy_xxyyzz[k];

                g_z_0_xxxxyy_xyzzz[k] = -g_z_0_xxxyy_xyzzz[k] * ab_x + g_z_0_xxxyy_xxyzzz[k];

                g_z_0_xxxxyy_xzzzz[k] = -g_z_0_xxxyy_xzzzz[k] * ab_x + g_z_0_xxxyy_xxzzzz[k];

                g_z_0_xxxxyy_yyyyy[k] = -g_z_0_xxxyy_yyyyy[k] * ab_x + g_z_0_xxxyy_xyyyyy[k];

                g_z_0_xxxxyy_yyyyz[k] = -g_z_0_xxxyy_yyyyz[k] * ab_x + g_z_0_xxxyy_xyyyyz[k];

                g_z_0_xxxxyy_yyyzz[k] = -g_z_0_xxxyy_yyyzz[k] * ab_x + g_z_0_xxxyy_xyyyzz[k];

                g_z_0_xxxxyy_yyzzz[k] = -g_z_0_xxxyy_yyzzz[k] * ab_x + g_z_0_xxxyy_xyyzzz[k];

                g_z_0_xxxxyy_yzzzz[k] = -g_z_0_xxxyy_yzzzz[k] * ab_x + g_z_0_xxxyy_xyzzzz[k];

                g_z_0_xxxxyy_zzzzz[k] = -g_z_0_xxxyy_zzzzz[k] * ab_x + g_z_0_xxxyy_xzzzzz[k];
            }

            /// Set up 1260-1281 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 1280 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyz_xxxxx, g_z_0_xxxxyz_xxxxy, g_z_0_xxxxyz_xxxxz, g_z_0_xxxxyz_xxxyy, g_z_0_xxxxyz_xxxyz, g_z_0_xxxxyz_xxxzz, g_z_0_xxxxyz_xxyyy, g_z_0_xxxxyz_xxyyz, g_z_0_xxxxyz_xxyzz, g_z_0_xxxxyz_xxzzz, g_z_0_xxxxyz_xyyyy, g_z_0_xxxxyz_xyyyz, g_z_0_xxxxyz_xyyzz, g_z_0_xxxxyz_xyzzz, g_z_0_xxxxyz_xzzzz, g_z_0_xxxxyz_yyyyy, g_z_0_xxxxyz_yyyyz, g_z_0_xxxxyz_yyyzz, g_z_0_xxxxyz_yyzzz, g_z_0_xxxxyz_yzzzz, g_z_0_xxxxyz_zzzzz, g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxxx, g_z_0_xxxyz_xxxxxy, g_z_0_xxxyz_xxxxxz, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxyy, g_z_0_xxxyz_xxxxyz, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxxzz, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyyy, g_z_0_xxxyz_xxxyyz, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxyzz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxxzzz, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyyy, g_z_0_xxxyz_xxyyyz, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyyzz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxyzzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xxzzzz, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyyy, g_z_0_xxxyz_xyyyyz, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyyzz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyyzzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xyzzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_xzzzzz, g_z_0_xxxyz_yyyyy, g_z_0_xxxyz_yyyyz, g_z_0_xxxyz_yyyzz, g_z_0_xxxyz_yyzzz, g_z_0_xxxyz_yzzzz, g_z_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxxxx[k] = -g_z_0_xxxyz_xxxxx[k] * ab_x + g_z_0_xxxyz_xxxxxx[k];

                g_z_0_xxxxyz_xxxxy[k] = -g_z_0_xxxyz_xxxxy[k] * ab_x + g_z_0_xxxyz_xxxxxy[k];

                g_z_0_xxxxyz_xxxxz[k] = -g_z_0_xxxyz_xxxxz[k] * ab_x + g_z_0_xxxyz_xxxxxz[k];

                g_z_0_xxxxyz_xxxyy[k] = -g_z_0_xxxyz_xxxyy[k] * ab_x + g_z_0_xxxyz_xxxxyy[k];

                g_z_0_xxxxyz_xxxyz[k] = -g_z_0_xxxyz_xxxyz[k] * ab_x + g_z_0_xxxyz_xxxxyz[k];

                g_z_0_xxxxyz_xxxzz[k] = -g_z_0_xxxyz_xxxzz[k] * ab_x + g_z_0_xxxyz_xxxxzz[k];

                g_z_0_xxxxyz_xxyyy[k] = -g_z_0_xxxyz_xxyyy[k] * ab_x + g_z_0_xxxyz_xxxyyy[k];

                g_z_0_xxxxyz_xxyyz[k] = -g_z_0_xxxyz_xxyyz[k] * ab_x + g_z_0_xxxyz_xxxyyz[k];

                g_z_0_xxxxyz_xxyzz[k] = -g_z_0_xxxyz_xxyzz[k] * ab_x + g_z_0_xxxyz_xxxyzz[k];

                g_z_0_xxxxyz_xxzzz[k] = -g_z_0_xxxyz_xxzzz[k] * ab_x + g_z_0_xxxyz_xxxzzz[k];

                g_z_0_xxxxyz_xyyyy[k] = -g_z_0_xxxyz_xyyyy[k] * ab_x + g_z_0_xxxyz_xxyyyy[k];

                g_z_0_xxxxyz_xyyyz[k] = -g_z_0_xxxyz_xyyyz[k] * ab_x + g_z_0_xxxyz_xxyyyz[k];

                g_z_0_xxxxyz_xyyzz[k] = -g_z_0_xxxyz_xyyzz[k] * ab_x + g_z_0_xxxyz_xxyyzz[k];

                g_z_0_xxxxyz_xyzzz[k] = -g_z_0_xxxyz_xyzzz[k] * ab_x + g_z_0_xxxyz_xxyzzz[k];

                g_z_0_xxxxyz_xzzzz[k] = -g_z_0_xxxyz_xzzzz[k] * ab_x + g_z_0_xxxyz_xxzzzz[k];

                g_z_0_xxxxyz_yyyyy[k] = -g_z_0_xxxyz_yyyyy[k] * ab_x + g_z_0_xxxyz_xyyyyy[k];

                g_z_0_xxxxyz_yyyyz[k] = -g_z_0_xxxyz_yyyyz[k] * ab_x + g_z_0_xxxyz_xyyyyz[k];

                g_z_0_xxxxyz_yyyzz[k] = -g_z_0_xxxyz_yyyzz[k] * ab_x + g_z_0_xxxyz_xyyyzz[k];

                g_z_0_xxxxyz_yyzzz[k] = -g_z_0_xxxyz_yyzzz[k] * ab_x + g_z_0_xxxyz_xyyzzz[k];

                g_z_0_xxxxyz_yzzzz[k] = -g_z_0_xxxyz_yzzzz[k] * ab_x + g_z_0_xxxyz_xyzzzz[k];

                g_z_0_xxxxyz_zzzzz[k] = -g_z_0_xxxyz_zzzzz[k] * ab_x + g_z_0_xxxyz_xzzzzz[k];
            }

            /// Set up 1281-1302 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 1281 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 1282 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 1283 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 1284 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 1285 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 1286 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 1287 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 1288 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 1289 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 1290 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 1291 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 1292 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 1293 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 1294 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 1295 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 1301 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxzz_xxxxx, g_z_0_xxxxzz_xxxxy, g_z_0_xxxxzz_xxxxz, g_z_0_xxxxzz_xxxyy, g_z_0_xxxxzz_xxxyz, g_z_0_xxxxzz_xxxzz, g_z_0_xxxxzz_xxyyy, g_z_0_xxxxzz_xxyyz, g_z_0_xxxxzz_xxyzz, g_z_0_xxxxzz_xxzzz, g_z_0_xxxxzz_xyyyy, g_z_0_xxxxzz_xyyyz, g_z_0_xxxxzz_xyyzz, g_z_0_xxxxzz_xyzzz, g_z_0_xxxxzz_xzzzz, g_z_0_xxxxzz_yyyyy, g_z_0_xxxxzz_yyyyz, g_z_0_xxxxzz_yyyzz, g_z_0_xxxxzz_yyzzz, g_z_0_xxxxzz_yzzzz, g_z_0_xxxxzz_zzzzz, g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxxx, g_z_0_xxxzz_xxxxxy, g_z_0_xxxzz_xxxxxz, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxyy, g_z_0_xxxzz_xxxxyz, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxxzz, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyyy, g_z_0_xxxzz_xxxyyz, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxyzz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxxzzz, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyyy, g_z_0_xxxzz_xxyyyz, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyyzz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxyzzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xxzzzz, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyyy, g_z_0_xxxzz_xyyyyz, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyyzz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyyzzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xyzzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_xzzzzz, g_z_0_xxxzz_yyyyy, g_z_0_xxxzz_yyyyz, g_z_0_xxxzz_yyyzz, g_z_0_xxxzz_yyzzz, g_z_0_xxxzz_yzzzz, g_z_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxxxx[k] = -g_z_0_xxxzz_xxxxx[k] * ab_x + g_z_0_xxxzz_xxxxxx[k];

                g_z_0_xxxxzz_xxxxy[k] = -g_z_0_xxxzz_xxxxy[k] * ab_x + g_z_0_xxxzz_xxxxxy[k];

                g_z_0_xxxxzz_xxxxz[k] = -g_z_0_xxxzz_xxxxz[k] * ab_x + g_z_0_xxxzz_xxxxxz[k];

                g_z_0_xxxxzz_xxxyy[k] = -g_z_0_xxxzz_xxxyy[k] * ab_x + g_z_0_xxxzz_xxxxyy[k];

                g_z_0_xxxxzz_xxxyz[k] = -g_z_0_xxxzz_xxxyz[k] * ab_x + g_z_0_xxxzz_xxxxyz[k];

                g_z_0_xxxxzz_xxxzz[k] = -g_z_0_xxxzz_xxxzz[k] * ab_x + g_z_0_xxxzz_xxxxzz[k];

                g_z_0_xxxxzz_xxyyy[k] = -g_z_0_xxxzz_xxyyy[k] * ab_x + g_z_0_xxxzz_xxxyyy[k];

                g_z_0_xxxxzz_xxyyz[k] = -g_z_0_xxxzz_xxyyz[k] * ab_x + g_z_0_xxxzz_xxxyyz[k];

                g_z_0_xxxxzz_xxyzz[k] = -g_z_0_xxxzz_xxyzz[k] * ab_x + g_z_0_xxxzz_xxxyzz[k];

                g_z_0_xxxxzz_xxzzz[k] = -g_z_0_xxxzz_xxzzz[k] * ab_x + g_z_0_xxxzz_xxxzzz[k];

                g_z_0_xxxxzz_xyyyy[k] = -g_z_0_xxxzz_xyyyy[k] * ab_x + g_z_0_xxxzz_xxyyyy[k];

                g_z_0_xxxxzz_xyyyz[k] = -g_z_0_xxxzz_xyyyz[k] * ab_x + g_z_0_xxxzz_xxyyyz[k];

                g_z_0_xxxxzz_xyyzz[k] = -g_z_0_xxxzz_xyyzz[k] * ab_x + g_z_0_xxxzz_xxyyzz[k];

                g_z_0_xxxxzz_xyzzz[k] = -g_z_0_xxxzz_xyzzz[k] * ab_x + g_z_0_xxxzz_xxyzzz[k];

                g_z_0_xxxxzz_xzzzz[k] = -g_z_0_xxxzz_xzzzz[k] * ab_x + g_z_0_xxxzz_xxzzzz[k];

                g_z_0_xxxxzz_yyyyy[k] = -g_z_0_xxxzz_yyyyy[k] * ab_x + g_z_0_xxxzz_xyyyyy[k];

                g_z_0_xxxxzz_yyyyz[k] = -g_z_0_xxxzz_yyyyz[k] * ab_x + g_z_0_xxxzz_xyyyyz[k];

                g_z_0_xxxxzz_yyyzz[k] = -g_z_0_xxxzz_yyyzz[k] * ab_x + g_z_0_xxxzz_xyyyzz[k];

                g_z_0_xxxxzz_yyzzz[k] = -g_z_0_xxxzz_yyzzz[k] * ab_x + g_z_0_xxxzz_xyyzzz[k];

                g_z_0_xxxxzz_yzzzz[k] = -g_z_0_xxxzz_yzzzz[k] * ab_x + g_z_0_xxxzz_xyzzzz[k];

                g_z_0_xxxxzz_zzzzz[k] = -g_z_0_xxxzz_zzzzz[k] * ab_x + g_z_0_xxxzz_xzzzzz[k];
            }

            /// Set up 1302-1323 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1309 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1310 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1311 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1312 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1313 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1314 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1315 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1322 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyy_xxxxx, g_z_0_xxxyyy_xxxxy, g_z_0_xxxyyy_xxxxz, g_z_0_xxxyyy_xxxyy, g_z_0_xxxyyy_xxxyz, g_z_0_xxxyyy_xxxzz, g_z_0_xxxyyy_xxyyy, g_z_0_xxxyyy_xxyyz, g_z_0_xxxyyy_xxyzz, g_z_0_xxxyyy_xxzzz, g_z_0_xxxyyy_xyyyy, g_z_0_xxxyyy_xyyyz, g_z_0_xxxyyy_xyyzz, g_z_0_xxxyyy_xyzzz, g_z_0_xxxyyy_xzzzz, g_z_0_xxxyyy_yyyyy, g_z_0_xxxyyy_yyyyz, g_z_0_xxxyyy_yyyzz, g_z_0_xxxyyy_yyzzz, g_z_0_xxxyyy_yzzzz, g_z_0_xxxyyy_zzzzz, g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxxx, g_z_0_xxyyy_xxxxxy, g_z_0_xxyyy_xxxxxz, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxyy, g_z_0_xxyyy_xxxxyz, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxxzz, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyyy, g_z_0_xxyyy_xxxyyz, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxyzz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxxzzz, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyyy, g_z_0_xxyyy_xxyyyz, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyyzz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxyzzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xxzzzz, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyyy, g_z_0_xxyyy_xyyyyz, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyyzz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyyzzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xyzzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_xzzzzz, g_z_0_xxyyy_yyyyy, g_z_0_xxyyy_yyyyz, g_z_0_xxyyy_yyyzz, g_z_0_xxyyy_yyzzz, g_z_0_xxyyy_yzzzz, g_z_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxxxx[k] = -g_z_0_xxyyy_xxxxx[k] * ab_x + g_z_0_xxyyy_xxxxxx[k];

                g_z_0_xxxyyy_xxxxy[k] = -g_z_0_xxyyy_xxxxy[k] * ab_x + g_z_0_xxyyy_xxxxxy[k];

                g_z_0_xxxyyy_xxxxz[k] = -g_z_0_xxyyy_xxxxz[k] * ab_x + g_z_0_xxyyy_xxxxxz[k];

                g_z_0_xxxyyy_xxxyy[k] = -g_z_0_xxyyy_xxxyy[k] * ab_x + g_z_0_xxyyy_xxxxyy[k];

                g_z_0_xxxyyy_xxxyz[k] = -g_z_0_xxyyy_xxxyz[k] * ab_x + g_z_0_xxyyy_xxxxyz[k];

                g_z_0_xxxyyy_xxxzz[k] = -g_z_0_xxyyy_xxxzz[k] * ab_x + g_z_0_xxyyy_xxxxzz[k];

                g_z_0_xxxyyy_xxyyy[k] = -g_z_0_xxyyy_xxyyy[k] * ab_x + g_z_0_xxyyy_xxxyyy[k];

                g_z_0_xxxyyy_xxyyz[k] = -g_z_0_xxyyy_xxyyz[k] * ab_x + g_z_0_xxyyy_xxxyyz[k];

                g_z_0_xxxyyy_xxyzz[k] = -g_z_0_xxyyy_xxyzz[k] * ab_x + g_z_0_xxyyy_xxxyzz[k];

                g_z_0_xxxyyy_xxzzz[k] = -g_z_0_xxyyy_xxzzz[k] * ab_x + g_z_0_xxyyy_xxxzzz[k];

                g_z_0_xxxyyy_xyyyy[k] = -g_z_0_xxyyy_xyyyy[k] * ab_x + g_z_0_xxyyy_xxyyyy[k];

                g_z_0_xxxyyy_xyyyz[k] = -g_z_0_xxyyy_xyyyz[k] * ab_x + g_z_0_xxyyy_xxyyyz[k];

                g_z_0_xxxyyy_xyyzz[k] = -g_z_0_xxyyy_xyyzz[k] * ab_x + g_z_0_xxyyy_xxyyzz[k];

                g_z_0_xxxyyy_xyzzz[k] = -g_z_0_xxyyy_xyzzz[k] * ab_x + g_z_0_xxyyy_xxyzzz[k];

                g_z_0_xxxyyy_xzzzz[k] = -g_z_0_xxyyy_xzzzz[k] * ab_x + g_z_0_xxyyy_xxzzzz[k];

                g_z_0_xxxyyy_yyyyy[k] = -g_z_0_xxyyy_yyyyy[k] * ab_x + g_z_0_xxyyy_xyyyyy[k];

                g_z_0_xxxyyy_yyyyz[k] = -g_z_0_xxyyy_yyyyz[k] * ab_x + g_z_0_xxyyy_xyyyyz[k];

                g_z_0_xxxyyy_yyyzz[k] = -g_z_0_xxyyy_yyyzz[k] * ab_x + g_z_0_xxyyy_xyyyzz[k];

                g_z_0_xxxyyy_yyzzz[k] = -g_z_0_xxyyy_yyzzz[k] * ab_x + g_z_0_xxyyy_xyyzzz[k];

                g_z_0_xxxyyy_yzzzz[k] = -g_z_0_xxyyy_yzzzz[k] * ab_x + g_z_0_xxyyy_xyzzzz[k];

                g_z_0_xxxyyy_zzzzz[k] = -g_z_0_xxyyy_zzzzz[k] * ab_x + g_z_0_xxyyy_xzzzzz[k];
            }

            /// Set up 1323-1344 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1323 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1324 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1325 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1326 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1327 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1328 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1329 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1330 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1331 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1332 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1333 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1334 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1335 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1336 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1337 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1338 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1339 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1340 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1341 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1342 * ccomps * dcomps);

            auto g_z_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyz_xxxxx, g_z_0_xxxyyz_xxxxy, g_z_0_xxxyyz_xxxxz, g_z_0_xxxyyz_xxxyy, g_z_0_xxxyyz_xxxyz, g_z_0_xxxyyz_xxxzz, g_z_0_xxxyyz_xxyyy, g_z_0_xxxyyz_xxyyz, g_z_0_xxxyyz_xxyzz, g_z_0_xxxyyz_xxzzz, g_z_0_xxxyyz_xyyyy, g_z_0_xxxyyz_xyyyz, g_z_0_xxxyyz_xyyzz, g_z_0_xxxyyz_xyzzz, g_z_0_xxxyyz_xzzzz, g_z_0_xxxyyz_yyyyy, g_z_0_xxxyyz_yyyyz, g_z_0_xxxyyz_yyyzz, g_z_0_xxxyyz_yyzzz, g_z_0_xxxyyz_yzzzz, g_z_0_xxxyyz_zzzzz, g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxxx, g_z_0_xxyyz_xxxxxy, g_z_0_xxyyz_xxxxxz, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxyy, g_z_0_xxyyz_xxxxyz, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxxzz, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyyy, g_z_0_xxyyz_xxxyyz, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxyzz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxxzzz, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyyy, g_z_0_xxyyz_xxyyyz, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyyzz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxyzzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xxzzzz, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyyy, g_z_0_xxyyz_xyyyyz, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyyzz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyyzzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xyzzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_xzzzzz, g_z_0_xxyyz_yyyyy, g_z_0_xxyyz_yyyyz, g_z_0_xxyyz_yyyzz, g_z_0_xxyyz_yyzzz, g_z_0_xxyyz_yzzzz, g_z_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxxxx[k] = -g_z_0_xxyyz_xxxxx[k] * ab_x + g_z_0_xxyyz_xxxxxx[k];

                g_z_0_xxxyyz_xxxxy[k] = -g_z_0_xxyyz_xxxxy[k] * ab_x + g_z_0_xxyyz_xxxxxy[k];

                g_z_0_xxxyyz_xxxxz[k] = -g_z_0_xxyyz_xxxxz[k] * ab_x + g_z_0_xxyyz_xxxxxz[k];

                g_z_0_xxxyyz_xxxyy[k] = -g_z_0_xxyyz_xxxyy[k] * ab_x + g_z_0_xxyyz_xxxxyy[k];

                g_z_0_xxxyyz_xxxyz[k] = -g_z_0_xxyyz_xxxyz[k] * ab_x + g_z_0_xxyyz_xxxxyz[k];

                g_z_0_xxxyyz_xxxzz[k] = -g_z_0_xxyyz_xxxzz[k] * ab_x + g_z_0_xxyyz_xxxxzz[k];

                g_z_0_xxxyyz_xxyyy[k] = -g_z_0_xxyyz_xxyyy[k] * ab_x + g_z_0_xxyyz_xxxyyy[k];

                g_z_0_xxxyyz_xxyyz[k] = -g_z_0_xxyyz_xxyyz[k] * ab_x + g_z_0_xxyyz_xxxyyz[k];

                g_z_0_xxxyyz_xxyzz[k] = -g_z_0_xxyyz_xxyzz[k] * ab_x + g_z_0_xxyyz_xxxyzz[k];

                g_z_0_xxxyyz_xxzzz[k] = -g_z_0_xxyyz_xxzzz[k] * ab_x + g_z_0_xxyyz_xxxzzz[k];

                g_z_0_xxxyyz_xyyyy[k] = -g_z_0_xxyyz_xyyyy[k] * ab_x + g_z_0_xxyyz_xxyyyy[k];

                g_z_0_xxxyyz_xyyyz[k] = -g_z_0_xxyyz_xyyyz[k] * ab_x + g_z_0_xxyyz_xxyyyz[k];

                g_z_0_xxxyyz_xyyzz[k] = -g_z_0_xxyyz_xyyzz[k] * ab_x + g_z_0_xxyyz_xxyyzz[k];

                g_z_0_xxxyyz_xyzzz[k] = -g_z_0_xxyyz_xyzzz[k] * ab_x + g_z_0_xxyyz_xxyzzz[k];

                g_z_0_xxxyyz_xzzzz[k] = -g_z_0_xxyyz_xzzzz[k] * ab_x + g_z_0_xxyyz_xxzzzz[k];

                g_z_0_xxxyyz_yyyyy[k] = -g_z_0_xxyyz_yyyyy[k] * ab_x + g_z_0_xxyyz_xyyyyy[k];

                g_z_0_xxxyyz_yyyyz[k] = -g_z_0_xxyyz_yyyyz[k] * ab_x + g_z_0_xxyyz_xyyyyz[k];

                g_z_0_xxxyyz_yyyzz[k] = -g_z_0_xxyyz_yyyzz[k] * ab_x + g_z_0_xxyyz_xyyyzz[k];

                g_z_0_xxxyyz_yyzzz[k] = -g_z_0_xxyyz_yyzzz[k] * ab_x + g_z_0_xxyyz_xyyzzz[k];

                g_z_0_xxxyyz_yzzzz[k] = -g_z_0_xxyyz_yzzzz[k] * ab_x + g_z_0_xxyyz_xyzzzz[k];

                g_z_0_xxxyyz_zzzzz[k] = -g_z_0_xxyyz_zzzzz[k] * ab_x + g_z_0_xxyyz_xzzzzz[k];
            }

            /// Set up 1344-1365 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1344 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1345 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1346 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1347 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1348 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1349 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1350 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1351 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1352 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1353 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1354 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1355 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1356 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1357 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1358 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1359 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1360 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1361 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1362 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1363 * ccomps * dcomps);

            auto g_z_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1364 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyzz_xxxxx, g_z_0_xxxyzz_xxxxy, g_z_0_xxxyzz_xxxxz, g_z_0_xxxyzz_xxxyy, g_z_0_xxxyzz_xxxyz, g_z_0_xxxyzz_xxxzz, g_z_0_xxxyzz_xxyyy, g_z_0_xxxyzz_xxyyz, g_z_0_xxxyzz_xxyzz, g_z_0_xxxyzz_xxzzz, g_z_0_xxxyzz_xyyyy, g_z_0_xxxyzz_xyyyz, g_z_0_xxxyzz_xyyzz, g_z_0_xxxyzz_xyzzz, g_z_0_xxxyzz_xzzzz, g_z_0_xxxyzz_yyyyy, g_z_0_xxxyzz_yyyyz, g_z_0_xxxyzz_yyyzz, g_z_0_xxxyzz_yyzzz, g_z_0_xxxyzz_yzzzz, g_z_0_xxxyzz_zzzzz, g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxxx, g_z_0_xxyzz_xxxxxy, g_z_0_xxyzz_xxxxxz, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxyy, g_z_0_xxyzz_xxxxyz, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxxzz, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyyy, g_z_0_xxyzz_xxxyyz, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxyzz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxxzzz, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyyy, g_z_0_xxyzz_xxyyyz, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyyzz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxyzzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xxzzzz, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyyy, g_z_0_xxyzz_xyyyyz, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyyzz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyyzzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xyzzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_xzzzzz, g_z_0_xxyzz_yyyyy, g_z_0_xxyzz_yyyyz, g_z_0_xxyzz_yyyzz, g_z_0_xxyzz_yyzzz, g_z_0_xxyzz_yzzzz, g_z_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxxxx[k] = -g_z_0_xxyzz_xxxxx[k] * ab_x + g_z_0_xxyzz_xxxxxx[k];

                g_z_0_xxxyzz_xxxxy[k] = -g_z_0_xxyzz_xxxxy[k] * ab_x + g_z_0_xxyzz_xxxxxy[k];

                g_z_0_xxxyzz_xxxxz[k] = -g_z_0_xxyzz_xxxxz[k] * ab_x + g_z_0_xxyzz_xxxxxz[k];

                g_z_0_xxxyzz_xxxyy[k] = -g_z_0_xxyzz_xxxyy[k] * ab_x + g_z_0_xxyzz_xxxxyy[k];

                g_z_0_xxxyzz_xxxyz[k] = -g_z_0_xxyzz_xxxyz[k] * ab_x + g_z_0_xxyzz_xxxxyz[k];

                g_z_0_xxxyzz_xxxzz[k] = -g_z_0_xxyzz_xxxzz[k] * ab_x + g_z_0_xxyzz_xxxxzz[k];

                g_z_0_xxxyzz_xxyyy[k] = -g_z_0_xxyzz_xxyyy[k] * ab_x + g_z_0_xxyzz_xxxyyy[k];

                g_z_0_xxxyzz_xxyyz[k] = -g_z_0_xxyzz_xxyyz[k] * ab_x + g_z_0_xxyzz_xxxyyz[k];

                g_z_0_xxxyzz_xxyzz[k] = -g_z_0_xxyzz_xxyzz[k] * ab_x + g_z_0_xxyzz_xxxyzz[k];

                g_z_0_xxxyzz_xxzzz[k] = -g_z_0_xxyzz_xxzzz[k] * ab_x + g_z_0_xxyzz_xxxzzz[k];

                g_z_0_xxxyzz_xyyyy[k] = -g_z_0_xxyzz_xyyyy[k] * ab_x + g_z_0_xxyzz_xxyyyy[k];

                g_z_0_xxxyzz_xyyyz[k] = -g_z_0_xxyzz_xyyyz[k] * ab_x + g_z_0_xxyzz_xxyyyz[k];

                g_z_0_xxxyzz_xyyzz[k] = -g_z_0_xxyzz_xyyzz[k] * ab_x + g_z_0_xxyzz_xxyyzz[k];

                g_z_0_xxxyzz_xyzzz[k] = -g_z_0_xxyzz_xyzzz[k] * ab_x + g_z_0_xxyzz_xxyzzz[k];

                g_z_0_xxxyzz_xzzzz[k] = -g_z_0_xxyzz_xzzzz[k] * ab_x + g_z_0_xxyzz_xxzzzz[k];

                g_z_0_xxxyzz_yyyyy[k] = -g_z_0_xxyzz_yyyyy[k] * ab_x + g_z_0_xxyzz_xyyyyy[k];

                g_z_0_xxxyzz_yyyyz[k] = -g_z_0_xxyzz_yyyyz[k] * ab_x + g_z_0_xxyzz_xyyyyz[k];

                g_z_0_xxxyzz_yyyzz[k] = -g_z_0_xxyzz_yyyzz[k] * ab_x + g_z_0_xxyzz_xyyyzz[k];

                g_z_0_xxxyzz_yyzzz[k] = -g_z_0_xxyzz_yyzzz[k] * ab_x + g_z_0_xxyzz_xyyzzz[k];

                g_z_0_xxxyzz_yzzzz[k] = -g_z_0_xxyzz_yzzzz[k] * ab_x + g_z_0_xxyzz_xyzzzz[k];

                g_z_0_xxxyzz_zzzzz[k] = -g_z_0_xxyzz_zzzzz[k] * ab_x + g_z_0_xxyzz_xzzzzz[k];
            }

            /// Set up 1365-1386 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1365 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1366 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1367 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1368 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1369 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1370 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1371 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1372 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1373 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1374 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1375 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1376 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1377 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1378 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1379 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1380 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1381 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1382 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1383 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1384 * ccomps * dcomps);

            auto g_z_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1385 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzzz_xxxxx, g_z_0_xxxzzz_xxxxy, g_z_0_xxxzzz_xxxxz, g_z_0_xxxzzz_xxxyy, g_z_0_xxxzzz_xxxyz, g_z_0_xxxzzz_xxxzz, g_z_0_xxxzzz_xxyyy, g_z_0_xxxzzz_xxyyz, g_z_0_xxxzzz_xxyzz, g_z_0_xxxzzz_xxzzz, g_z_0_xxxzzz_xyyyy, g_z_0_xxxzzz_xyyyz, g_z_0_xxxzzz_xyyzz, g_z_0_xxxzzz_xyzzz, g_z_0_xxxzzz_xzzzz, g_z_0_xxxzzz_yyyyy, g_z_0_xxxzzz_yyyyz, g_z_0_xxxzzz_yyyzz, g_z_0_xxxzzz_yyzzz, g_z_0_xxxzzz_yzzzz, g_z_0_xxxzzz_zzzzz, g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxxx, g_z_0_xxzzz_xxxxxy, g_z_0_xxzzz_xxxxxz, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxyy, g_z_0_xxzzz_xxxxyz, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxxzz, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyyy, g_z_0_xxzzz_xxxyyz, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxyzz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxxzzz, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyyy, g_z_0_xxzzz_xxyyyz, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyyzz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxyzzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xxzzzz, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyyy, g_z_0_xxzzz_xyyyyz, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyyzz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyyzzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xyzzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_xzzzzz, g_z_0_xxzzz_yyyyy, g_z_0_xxzzz_yyyyz, g_z_0_xxzzz_yyyzz, g_z_0_xxzzz_yyzzz, g_z_0_xxzzz_yzzzz, g_z_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxxxx[k] = -g_z_0_xxzzz_xxxxx[k] * ab_x + g_z_0_xxzzz_xxxxxx[k];

                g_z_0_xxxzzz_xxxxy[k] = -g_z_0_xxzzz_xxxxy[k] * ab_x + g_z_0_xxzzz_xxxxxy[k];

                g_z_0_xxxzzz_xxxxz[k] = -g_z_0_xxzzz_xxxxz[k] * ab_x + g_z_0_xxzzz_xxxxxz[k];

                g_z_0_xxxzzz_xxxyy[k] = -g_z_0_xxzzz_xxxyy[k] * ab_x + g_z_0_xxzzz_xxxxyy[k];

                g_z_0_xxxzzz_xxxyz[k] = -g_z_0_xxzzz_xxxyz[k] * ab_x + g_z_0_xxzzz_xxxxyz[k];

                g_z_0_xxxzzz_xxxzz[k] = -g_z_0_xxzzz_xxxzz[k] * ab_x + g_z_0_xxzzz_xxxxzz[k];

                g_z_0_xxxzzz_xxyyy[k] = -g_z_0_xxzzz_xxyyy[k] * ab_x + g_z_0_xxzzz_xxxyyy[k];

                g_z_0_xxxzzz_xxyyz[k] = -g_z_0_xxzzz_xxyyz[k] * ab_x + g_z_0_xxzzz_xxxyyz[k];

                g_z_0_xxxzzz_xxyzz[k] = -g_z_0_xxzzz_xxyzz[k] * ab_x + g_z_0_xxzzz_xxxyzz[k];

                g_z_0_xxxzzz_xxzzz[k] = -g_z_0_xxzzz_xxzzz[k] * ab_x + g_z_0_xxzzz_xxxzzz[k];

                g_z_0_xxxzzz_xyyyy[k] = -g_z_0_xxzzz_xyyyy[k] * ab_x + g_z_0_xxzzz_xxyyyy[k];

                g_z_0_xxxzzz_xyyyz[k] = -g_z_0_xxzzz_xyyyz[k] * ab_x + g_z_0_xxzzz_xxyyyz[k];

                g_z_0_xxxzzz_xyyzz[k] = -g_z_0_xxzzz_xyyzz[k] * ab_x + g_z_0_xxzzz_xxyyzz[k];

                g_z_0_xxxzzz_xyzzz[k] = -g_z_0_xxzzz_xyzzz[k] * ab_x + g_z_0_xxzzz_xxyzzz[k];

                g_z_0_xxxzzz_xzzzz[k] = -g_z_0_xxzzz_xzzzz[k] * ab_x + g_z_0_xxzzz_xxzzzz[k];

                g_z_0_xxxzzz_yyyyy[k] = -g_z_0_xxzzz_yyyyy[k] * ab_x + g_z_0_xxzzz_xyyyyy[k];

                g_z_0_xxxzzz_yyyyz[k] = -g_z_0_xxzzz_yyyyz[k] * ab_x + g_z_0_xxzzz_xyyyyz[k];

                g_z_0_xxxzzz_yyyzz[k] = -g_z_0_xxzzz_yyyzz[k] * ab_x + g_z_0_xxzzz_xyyyzz[k];

                g_z_0_xxxzzz_yyzzz[k] = -g_z_0_xxzzz_yyzzz[k] * ab_x + g_z_0_xxzzz_xyyzzz[k];

                g_z_0_xxxzzz_yzzzz[k] = -g_z_0_xxzzz_yzzzz[k] * ab_x + g_z_0_xxzzz_xyzzzz[k];

                g_z_0_xxxzzz_zzzzz[k] = -g_z_0_xxzzz_zzzzz[k] * ab_x + g_z_0_xxzzz_xzzzzz[k];
            }

            /// Set up 1386-1407 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1386 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1387 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1388 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1389 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1390 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1391 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1392 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1393 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1394 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1395 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1396 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1397 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1398 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1399 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1400 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1401 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1402 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1403 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1404 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1405 * ccomps * dcomps);

            auto g_z_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1406 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyy_xxxxx, g_z_0_xxyyyy_xxxxy, g_z_0_xxyyyy_xxxxz, g_z_0_xxyyyy_xxxyy, g_z_0_xxyyyy_xxxyz, g_z_0_xxyyyy_xxxzz, g_z_0_xxyyyy_xxyyy, g_z_0_xxyyyy_xxyyz, g_z_0_xxyyyy_xxyzz, g_z_0_xxyyyy_xxzzz, g_z_0_xxyyyy_xyyyy, g_z_0_xxyyyy_xyyyz, g_z_0_xxyyyy_xyyzz, g_z_0_xxyyyy_xyzzz, g_z_0_xxyyyy_xzzzz, g_z_0_xxyyyy_yyyyy, g_z_0_xxyyyy_yyyyz, g_z_0_xxyyyy_yyyzz, g_z_0_xxyyyy_yyzzz, g_z_0_xxyyyy_yzzzz, g_z_0_xxyyyy_zzzzz, g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxxx, g_z_0_xyyyy_xxxxxy, g_z_0_xyyyy_xxxxxz, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxyy, g_z_0_xyyyy_xxxxyz, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxxzz, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyyy, g_z_0_xyyyy_xxxyyz, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxyzz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxxzzz, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyyy, g_z_0_xyyyy_xxyyyz, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyyzz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxyzzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xxzzzz, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyyy, g_z_0_xyyyy_xyyyyz, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyyzz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyyzzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xyzzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_xzzzzz, g_z_0_xyyyy_yyyyy, g_z_0_xyyyy_yyyyz, g_z_0_xyyyy_yyyzz, g_z_0_xyyyy_yyzzz, g_z_0_xyyyy_yzzzz, g_z_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxxxx[k] = -g_z_0_xyyyy_xxxxx[k] * ab_x + g_z_0_xyyyy_xxxxxx[k];

                g_z_0_xxyyyy_xxxxy[k] = -g_z_0_xyyyy_xxxxy[k] * ab_x + g_z_0_xyyyy_xxxxxy[k];

                g_z_0_xxyyyy_xxxxz[k] = -g_z_0_xyyyy_xxxxz[k] * ab_x + g_z_0_xyyyy_xxxxxz[k];

                g_z_0_xxyyyy_xxxyy[k] = -g_z_0_xyyyy_xxxyy[k] * ab_x + g_z_0_xyyyy_xxxxyy[k];

                g_z_0_xxyyyy_xxxyz[k] = -g_z_0_xyyyy_xxxyz[k] * ab_x + g_z_0_xyyyy_xxxxyz[k];

                g_z_0_xxyyyy_xxxzz[k] = -g_z_0_xyyyy_xxxzz[k] * ab_x + g_z_0_xyyyy_xxxxzz[k];

                g_z_0_xxyyyy_xxyyy[k] = -g_z_0_xyyyy_xxyyy[k] * ab_x + g_z_0_xyyyy_xxxyyy[k];

                g_z_0_xxyyyy_xxyyz[k] = -g_z_0_xyyyy_xxyyz[k] * ab_x + g_z_0_xyyyy_xxxyyz[k];

                g_z_0_xxyyyy_xxyzz[k] = -g_z_0_xyyyy_xxyzz[k] * ab_x + g_z_0_xyyyy_xxxyzz[k];

                g_z_0_xxyyyy_xxzzz[k] = -g_z_0_xyyyy_xxzzz[k] * ab_x + g_z_0_xyyyy_xxxzzz[k];

                g_z_0_xxyyyy_xyyyy[k] = -g_z_0_xyyyy_xyyyy[k] * ab_x + g_z_0_xyyyy_xxyyyy[k];

                g_z_0_xxyyyy_xyyyz[k] = -g_z_0_xyyyy_xyyyz[k] * ab_x + g_z_0_xyyyy_xxyyyz[k];

                g_z_0_xxyyyy_xyyzz[k] = -g_z_0_xyyyy_xyyzz[k] * ab_x + g_z_0_xyyyy_xxyyzz[k];

                g_z_0_xxyyyy_xyzzz[k] = -g_z_0_xyyyy_xyzzz[k] * ab_x + g_z_0_xyyyy_xxyzzz[k];

                g_z_0_xxyyyy_xzzzz[k] = -g_z_0_xyyyy_xzzzz[k] * ab_x + g_z_0_xyyyy_xxzzzz[k];

                g_z_0_xxyyyy_yyyyy[k] = -g_z_0_xyyyy_yyyyy[k] * ab_x + g_z_0_xyyyy_xyyyyy[k];

                g_z_0_xxyyyy_yyyyz[k] = -g_z_0_xyyyy_yyyyz[k] * ab_x + g_z_0_xyyyy_xyyyyz[k];

                g_z_0_xxyyyy_yyyzz[k] = -g_z_0_xyyyy_yyyzz[k] * ab_x + g_z_0_xyyyy_xyyyzz[k];

                g_z_0_xxyyyy_yyzzz[k] = -g_z_0_xyyyy_yyzzz[k] * ab_x + g_z_0_xyyyy_xyyzzz[k];

                g_z_0_xxyyyy_yzzzz[k] = -g_z_0_xyyyy_yzzzz[k] * ab_x + g_z_0_xyyyy_xyzzzz[k];

                g_z_0_xxyyyy_zzzzz[k] = -g_z_0_xyyyy_zzzzz[k] * ab_x + g_z_0_xyyyy_xzzzzz[k];
            }

            /// Set up 1407-1428 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1407 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1408 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1409 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1410 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1411 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1412 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1413 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1414 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1415 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1416 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1417 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1418 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1419 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1420 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1421 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1422 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1423 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1424 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1425 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1426 * ccomps * dcomps);

            auto g_z_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyz_xxxxx, g_z_0_xxyyyz_xxxxy, g_z_0_xxyyyz_xxxxz, g_z_0_xxyyyz_xxxyy, g_z_0_xxyyyz_xxxyz, g_z_0_xxyyyz_xxxzz, g_z_0_xxyyyz_xxyyy, g_z_0_xxyyyz_xxyyz, g_z_0_xxyyyz_xxyzz, g_z_0_xxyyyz_xxzzz, g_z_0_xxyyyz_xyyyy, g_z_0_xxyyyz_xyyyz, g_z_0_xxyyyz_xyyzz, g_z_0_xxyyyz_xyzzz, g_z_0_xxyyyz_xzzzz, g_z_0_xxyyyz_yyyyy, g_z_0_xxyyyz_yyyyz, g_z_0_xxyyyz_yyyzz, g_z_0_xxyyyz_yyzzz, g_z_0_xxyyyz_yzzzz, g_z_0_xxyyyz_zzzzz, g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxxx, g_z_0_xyyyz_xxxxxy, g_z_0_xyyyz_xxxxxz, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxyy, g_z_0_xyyyz_xxxxyz, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxxzz, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyyy, g_z_0_xyyyz_xxxyyz, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxyzz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxxzzz, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyyy, g_z_0_xyyyz_xxyyyz, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyyzz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxyzzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xxzzzz, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyyy, g_z_0_xyyyz_xyyyyz, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyyzz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyyzzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xyzzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_xzzzzz, g_z_0_xyyyz_yyyyy, g_z_0_xyyyz_yyyyz, g_z_0_xyyyz_yyyzz, g_z_0_xyyyz_yyzzz, g_z_0_xyyyz_yzzzz, g_z_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxxxx[k] = -g_z_0_xyyyz_xxxxx[k] * ab_x + g_z_0_xyyyz_xxxxxx[k];

                g_z_0_xxyyyz_xxxxy[k] = -g_z_0_xyyyz_xxxxy[k] * ab_x + g_z_0_xyyyz_xxxxxy[k];

                g_z_0_xxyyyz_xxxxz[k] = -g_z_0_xyyyz_xxxxz[k] * ab_x + g_z_0_xyyyz_xxxxxz[k];

                g_z_0_xxyyyz_xxxyy[k] = -g_z_0_xyyyz_xxxyy[k] * ab_x + g_z_0_xyyyz_xxxxyy[k];

                g_z_0_xxyyyz_xxxyz[k] = -g_z_0_xyyyz_xxxyz[k] * ab_x + g_z_0_xyyyz_xxxxyz[k];

                g_z_0_xxyyyz_xxxzz[k] = -g_z_0_xyyyz_xxxzz[k] * ab_x + g_z_0_xyyyz_xxxxzz[k];

                g_z_0_xxyyyz_xxyyy[k] = -g_z_0_xyyyz_xxyyy[k] * ab_x + g_z_0_xyyyz_xxxyyy[k];

                g_z_0_xxyyyz_xxyyz[k] = -g_z_0_xyyyz_xxyyz[k] * ab_x + g_z_0_xyyyz_xxxyyz[k];

                g_z_0_xxyyyz_xxyzz[k] = -g_z_0_xyyyz_xxyzz[k] * ab_x + g_z_0_xyyyz_xxxyzz[k];

                g_z_0_xxyyyz_xxzzz[k] = -g_z_0_xyyyz_xxzzz[k] * ab_x + g_z_0_xyyyz_xxxzzz[k];

                g_z_0_xxyyyz_xyyyy[k] = -g_z_0_xyyyz_xyyyy[k] * ab_x + g_z_0_xyyyz_xxyyyy[k];

                g_z_0_xxyyyz_xyyyz[k] = -g_z_0_xyyyz_xyyyz[k] * ab_x + g_z_0_xyyyz_xxyyyz[k];

                g_z_0_xxyyyz_xyyzz[k] = -g_z_0_xyyyz_xyyzz[k] * ab_x + g_z_0_xyyyz_xxyyzz[k];

                g_z_0_xxyyyz_xyzzz[k] = -g_z_0_xyyyz_xyzzz[k] * ab_x + g_z_0_xyyyz_xxyzzz[k];

                g_z_0_xxyyyz_xzzzz[k] = -g_z_0_xyyyz_xzzzz[k] * ab_x + g_z_0_xyyyz_xxzzzz[k];

                g_z_0_xxyyyz_yyyyy[k] = -g_z_0_xyyyz_yyyyy[k] * ab_x + g_z_0_xyyyz_xyyyyy[k];

                g_z_0_xxyyyz_yyyyz[k] = -g_z_0_xyyyz_yyyyz[k] * ab_x + g_z_0_xyyyz_xyyyyz[k];

                g_z_0_xxyyyz_yyyzz[k] = -g_z_0_xyyyz_yyyzz[k] * ab_x + g_z_0_xyyyz_xyyyzz[k];

                g_z_0_xxyyyz_yyzzz[k] = -g_z_0_xyyyz_yyzzz[k] * ab_x + g_z_0_xyyyz_xyyzzz[k];

                g_z_0_xxyyyz_yzzzz[k] = -g_z_0_xyyyz_yzzzz[k] * ab_x + g_z_0_xyyyz_xyzzzz[k];

                g_z_0_xxyyyz_zzzzz[k] = -g_z_0_xyyyz_zzzzz[k] * ab_x + g_z_0_xyyyz_xzzzzz[k];
            }

            /// Set up 1428-1449 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1428 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1429 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1430 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1431 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1432 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1433 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1434 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1435 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1436 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1437 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1438 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1439 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1440 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1441 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1442 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1443 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1444 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1445 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1446 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1447 * ccomps * dcomps);

            auto g_z_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1448 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyzz_xxxxx, g_z_0_xxyyzz_xxxxy, g_z_0_xxyyzz_xxxxz, g_z_0_xxyyzz_xxxyy, g_z_0_xxyyzz_xxxyz, g_z_0_xxyyzz_xxxzz, g_z_0_xxyyzz_xxyyy, g_z_0_xxyyzz_xxyyz, g_z_0_xxyyzz_xxyzz, g_z_0_xxyyzz_xxzzz, g_z_0_xxyyzz_xyyyy, g_z_0_xxyyzz_xyyyz, g_z_0_xxyyzz_xyyzz, g_z_0_xxyyzz_xyzzz, g_z_0_xxyyzz_xzzzz, g_z_0_xxyyzz_yyyyy, g_z_0_xxyyzz_yyyyz, g_z_0_xxyyzz_yyyzz, g_z_0_xxyyzz_yyzzz, g_z_0_xxyyzz_yzzzz, g_z_0_xxyyzz_zzzzz, g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxxx, g_z_0_xyyzz_xxxxxy, g_z_0_xyyzz_xxxxxz, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxyy, g_z_0_xyyzz_xxxxyz, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxxzz, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyyy, g_z_0_xyyzz_xxxyyz, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxyzz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxxzzz, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyyy, g_z_0_xyyzz_xxyyyz, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyyzz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxyzzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xxzzzz, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyyy, g_z_0_xyyzz_xyyyyz, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyyzz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyyzzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xyzzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_xzzzzz, g_z_0_xyyzz_yyyyy, g_z_0_xyyzz_yyyyz, g_z_0_xyyzz_yyyzz, g_z_0_xyyzz_yyzzz, g_z_0_xyyzz_yzzzz, g_z_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxxxx[k] = -g_z_0_xyyzz_xxxxx[k] * ab_x + g_z_0_xyyzz_xxxxxx[k];

                g_z_0_xxyyzz_xxxxy[k] = -g_z_0_xyyzz_xxxxy[k] * ab_x + g_z_0_xyyzz_xxxxxy[k];

                g_z_0_xxyyzz_xxxxz[k] = -g_z_0_xyyzz_xxxxz[k] * ab_x + g_z_0_xyyzz_xxxxxz[k];

                g_z_0_xxyyzz_xxxyy[k] = -g_z_0_xyyzz_xxxyy[k] * ab_x + g_z_0_xyyzz_xxxxyy[k];

                g_z_0_xxyyzz_xxxyz[k] = -g_z_0_xyyzz_xxxyz[k] * ab_x + g_z_0_xyyzz_xxxxyz[k];

                g_z_0_xxyyzz_xxxzz[k] = -g_z_0_xyyzz_xxxzz[k] * ab_x + g_z_0_xyyzz_xxxxzz[k];

                g_z_0_xxyyzz_xxyyy[k] = -g_z_0_xyyzz_xxyyy[k] * ab_x + g_z_0_xyyzz_xxxyyy[k];

                g_z_0_xxyyzz_xxyyz[k] = -g_z_0_xyyzz_xxyyz[k] * ab_x + g_z_0_xyyzz_xxxyyz[k];

                g_z_0_xxyyzz_xxyzz[k] = -g_z_0_xyyzz_xxyzz[k] * ab_x + g_z_0_xyyzz_xxxyzz[k];

                g_z_0_xxyyzz_xxzzz[k] = -g_z_0_xyyzz_xxzzz[k] * ab_x + g_z_0_xyyzz_xxxzzz[k];

                g_z_0_xxyyzz_xyyyy[k] = -g_z_0_xyyzz_xyyyy[k] * ab_x + g_z_0_xyyzz_xxyyyy[k];

                g_z_0_xxyyzz_xyyyz[k] = -g_z_0_xyyzz_xyyyz[k] * ab_x + g_z_0_xyyzz_xxyyyz[k];

                g_z_0_xxyyzz_xyyzz[k] = -g_z_0_xyyzz_xyyzz[k] * ab_x + g_z_0_xyyzz_xxyyzz[k];

                g_z_0_xxyyzz_xyzzz[k] = -g_z_0_xyyzz_xyzzz[k] * ab_x + g_z_0_xyyzz_xxyzzz[k];

                g_z_0_xxyyzz_xzzzz[k] = -g_z_0_xyyzz_xzzzz[k] * ab_x + g_z_0_xyyzz_xxzzzz[k];

                g_z_0_xxyyzz_yyyyy[k] = -g_z_0_xyyzz_yyyyy[k] * ab_x + g_z_0_xyyzz_xyyyyy[k];

                g_z_0_xxyyzz_yyyyz[k] = -g_z_0_xyyzz_yyyyz[k] * ab_x + g_z_0_xyyzz_xyyyyz[k];

                g_z_0_xxyyzz_yyyzz[k] = -g_z_0_xyyzz_yyyzz[k] * ab_x + g_z_0_xyyzz_xyyyzz[k];

                g_z_0_xxyyzz_yyzzz[k] = -g_z_0_xyyzz_yyzzz[k] * ab_x + g_z_0_xyyzz_xyyzzz[k];

                g_z_0_xxyyzz_yzzzz[k] = -g_z_0_xyyzz_yzzzz[k] * ab_x + g_z_0_xyyzz_xyzzzz[k];

                g_z_0_xxyyzz_zzzzz[k] = -g_z_0_xyyzz_zzzzz[k] * ab_x + g_z_0_xyyzz_xzzzzz[k];
            }

            /// Set up 1449-1470 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1449 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1450 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1451 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1452 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1453 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1454 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1455 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1456 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1457 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1458 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1459 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1460 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1461 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1462 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1463 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1464 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1465 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1466 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1467 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1468 * ccomps * dcomps);

            auto g_z_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzzz_xxxxx, g_z_0_xxyzzz_xxxxy, g_z_0_xxyzzz_xxxxz, g_z_0_xxyzzz_xxxyy, g_z_0_xxyzzz_xxxyz, g_z_0_xxyzzz_xxxzz, g_z_0_xxyzzz_xxyyy, g_z_0_xxyzzz_xxyyz, g_z_0_xxyzzz_xxyzz, g_z_0_xxyzzz_xxzzz, g_z_0_xxyzzz_xyyyy, g_z_0_xxyzzz_xyyyz, g_z_0_xxyzzz_xyyzz, g_z_0_xxyzzz_xyzzz, g_z_0_xxyzzz_xzzzz, g_z_0_xxyzzz_yyyyy, g_z_0_xxyzzz_yyyyz, g_z_0_xxyzzz_yyyzz, g_z_0_xxyzzz_yyzzz, g_z_0_xxyzzz_yzzzz, g_z_0_xxyzzz_zzzzz, g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxxx, g_z_0_xyzzz_xxxxxy, g_z_0_xyzzz_xxxxxz, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxyy, g_z_0_xyzzz_xxxxyz, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxxzz, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyyy, g_z_0_xyzzz_xxxyyz, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxyzz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxxzzz, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyyy, g_z_0_xyzzz_xxyyyz, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyyzz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxyzzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xxzzzz, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyyy, g_z_0_xyzzz_xyyyyz, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyyzz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyyzzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xyzzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_xzzzzz, g_z_0_xyzzz_yyyyy, g_z_0_xyzzz_yyyyz, g_z_0_xyzzz_yyyzz, g_z_0_xyzzz_yyzzz, g_z_0_xyzzz_yzzzz, g_z_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxxxx[k] = -g_z_0_xyzzz_xxxxx[k] * ab_x + g_z_0_xyzzz_xxxxxx[k];

                g_z_0_xxyzzz_xxxxy[k] = -g_z_0_xyzzz_xxxxy[k] * ab_x + g_z_0_xyzzz_xxxxxy[k];

                g_z_0_xxyzzz_xxxxz[k] = -g_z_0_xyzzz_xxxxz[k] * ab_x + g_z_0_xyzzz_xxxxxz[k];

                g_z_0_xxyzzz_xxxyy[k] = -g_z_0_xyzzz_xxxyy[k] * ab_x + g_z_0_xyzzz_xxxxyy[k];

                g_z_0_xxyzzz_xxxyz[k] = -g_z_0_xyzzz_xxxyz[k] * ab_x + g_z_0_xyzzz_xxxxyz[k];

                g_z_0_xxyzzz_xxxzz[k] = -g_z_0_xyzzz_xxxzz[k] * ab_x + g_z_0_xyzzz_xxxxzz[k];

                g_z_0_xxyzzz_xxyyy[k] = -g_z_0_xyzzz_xxyyy[k] * ab_x + g_z_0_xyzzz_xxxyyy[k];

                g_z_0_xxyzzz_xxyyz[k] = -g_z_0_xyzzz_xxyyz[k] * ab_x + g_z_0_xyzzz_xxxyyz[k];

                g_z_0_xxyzzz_xxyzz[k] = -g_z_0_xyzzz_xxyzz[k] * ab_x + g_z_0_xyzzz_xxxyzz[k];

                g_z_0_xxyzzz_xxzzz[k] = -g_z_0_xyzzz_xxzzz[k] * ab_x + g_z_0_xyzzz_xxxzzz[k];

                g_z_0_xxyzzz_xyyyy[k] = -g_z_0_xyzzz_xyyyy[k] * ab_x + g_z_0_xyzzz_xxyyyy[k];

                g_z_0_xxyzzz_xyyyz[k] = -g_z_0_xyzzz_xyyyz[k] * ab_x + g_z_0_xyzzz_xxyyyz[k];

                g_z_0_xxyzzz_xyyzz[k] = -g_z_0_xyzzz_xyyzz[k] * ab_x + g_z_0_xyzzz_xxyyzz[k];

                g_z_0_xxyzzz_xyzzz[k] = -g_z_0_xyzzz_xyzzz[k] * ab_x + g_z_0_xyzzz_xxyzzz[k];

                g_z_0_xxyzzz_xzzzz[k] = -g_z_0_xyzzz_xzzzz[k] * ab_x + g_z_0_xyzzz_xxzzzz[k];

                g_z_0_xxyzzz_yyyyy[k] = -g_z_0_xyzzz_yyyyy[k] * ab_x + g_z_0_xyzzz_xyyyyy[k];

                g_z_0_xxyzzz_yyyyz[k] = -g_z_0_xyzzz_yyyyz[k] * ab_x + g_z_0_xyzzz_xyyyyz[k];

                g_z_0_xxyzzz_yyyzz[k] = -g_z_0_xyzzz_yyyzz[k] * ab_x + g_z_0_xyzzz_xyyyzz[k];

                g_z_0_xxyzzz_yyzzz[k] = -g_z_0_xyzzz_yyzzz[k] * ab_x + g_z_0_xyzzz_xyyzzz[k];

                g_z_0_xxyzzz_yzzzz[k] = -g_z_0_xyzzz_yzzzz[k] * ab_x + g_z_0_xyzzz_xyzzzz[k];

                g_z_0_xxyzzz_zzzzz[k] = -g_z_0_xyzzz_zzzzz[k] * ab_x + g_z_0_xyzzz_xzzzzz[k];
            }

            /// Set up 1470-1491 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1470 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1471 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1472 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1473 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1474 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1475 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1476 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1477 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1478 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1479 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1480 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1481 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1482 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1483 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1484 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1485 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1486 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1487 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1488 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1489 * ccomps * dcomps);

            auto g_z_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1490 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzzz_xxxxx, g_z_0_xxzzzz_xxxxy, g_z_0_xxzzzz_xxxxz, g_z_0_xxzzzz_xxxyy, g_z_0_xxzzzz_xxxyz, g_z_0_xxzzzz_xxxzz, g_z_0_xxzzzz_xxyyy, g_z_0_xxzzzz_xxyyz, g_z_0_xxzzzz_xxyzz, g_z_0_xxzzzz_xxzzz, g_z_0_xxzzzz_xyyyy, g_z_0_xxzzzz_xyyyz, g_z_0_xxzzzz_xyyzz, g_z_0_xxzzzz_xyzzz, g_z_0_xxzzzz_xzzzz, g_z_0_xxzzzz_yyyyy, g_z_0_xxzzzz_yyyyz, g_z_0_xxzzzz_yyyzz, g_z_0_xxzzzz_yyzzz, g_z_0_xxzzzz_yzzzz, g_z_0_xxzzzz_zzzzz, g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxxx, g_z_0_xzzzz_xxxxxy, g_z_0_xzzzz_xxxxxz, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxyy, g_z_0_xzzzz_xxxxyz, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxxzz, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyyy, g_z_0_xzzzz_xxxyyz, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxyzz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxxzzz, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyyy, g_z_0_xzzzz_xxyyyz, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyyzz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxyzzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xxzzzz, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyyy, g_z_0_xzzzz_xyyyyz, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyyzz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyyzzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xyzzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_xzzzzz, g_z_0_xzzzz_yyyyy, g_z_0_xzzzz_yyyyz, g_z_0_xzzzz_yyyzz, g_z_0_xzzzz_yyzzz, g_z_0_xzzzz_yzzzz, g_z_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxxxx[k] = -g_z_0_xzzzz_xxxxx[k] * ab_x + g_z_0_xzzzz_xxxxxx[k];

                g_z_0_xxzzzz_xxxxy[k] = -g_z_0_xzzzz_xxxxy[k] * ab_x + g_z_0_xzzzz_xxxxxy[k];

                g_z_0_xxzzzz_xxxxz[k] = -g_z_0_xzzzz_xxxxz[k] * ab_x + g_z_0_xzzzz_xxxxxz[k];

                g_z_0_xxzzzz_xxxyy[k] = -g_z_0_xzzzz_xxxyy[k] * ab_x + g_z_0_xzzzz_xxxxyy[k];

                g_z_0_xxzzzz_xxxyz[k] = -g_z_0_xzzzz_xxxyz[k] * ab_x + g_z_0_xzzzz_xxxxyz[k];

                g_z_0_xxzzzz_xxxzz[k] = -g_z_0_xzzzz_xxxzz[k] * ab_x + g_z_0_xzzzz_xxxxzz[k];

                g_z_0_xxzzzz_xxyyy[k] = -g_z_0_xzzzz_xxyyy[k] * ab_x + g_z_0_xzzzz_xxxyyy[k];

                g_z_0_xxzzzz_xxyyz[k] = -g_z_0_xzzzz_xxyyz[k] * ab_x + g_z_0_xzzzz_xxxyyz[k];

                g_z_0_xxzzzz_xxyzz[k] = -g_z_0_xzzzz_xxyzz[k] * ab_x + g_z_0_xzzzz_xxxyzz[k];

                g_z_0_xxzzzz_xxzzz[k] = -g_z_0_xzzzz_xxzzz[k] * ab_x + g_z_0_xzzzz_xxxzzz[k];

                g_z_0_xxzzzz_xyyyy[k] = -g_z_0_xzzzz_xyyyy[k] * ab_x + g_z_0_xzzzz_xxyyyy[k];

                g_z_0_xxzzzz_xyyyz[k] = -g_z_0_xzzzz_xyyyz[k] * ab_x + g_z_0_xzzzz_xxyyyz[k];

                g_z_0_xxzzzz_xyyzz[k] = -g_z_0_xzzzz_xyyzz[k] * ab_x + g_z_0_xzzzz_xxyyzz[k];

                g_z_0_xxzzzz_xyzzz[k] = -g_z_0_xzzzz_xyzzz[k] * ab_x + g_z_0_xzzzz_xxyzzz[k];

                g_z_0_xxzzzz_xzzzz[k] = -g_z_0_xzzzz_xzzzz[k] * ab_x + g_z_0_xzzzz_xxzzzz[k];

                g_z_0_xxzzzz_yyyyy[k] = -g_z_0_xzzzz_yyyyy[k] * ab_x + g_z_0_xzzzz_xyyyyy[k];

                g_z_0_xxzzzz_yyyyz[k] = -g_z_0_xzzzz_yyyyz[k] * ab_x + g_z_0_xzzzz_xyyyyz[k];

                g_z_0_xxzzzz_yyyzz[k] = -g_z_0_xzzzz_yyyzz[k] * ab_x + g_z_0_xzzzz_xyyyzz[k];

                g_z_0_xxzzzz_yyzzz[k] = -g_z_0_xzzzz_yyzzz[k] * ab_x + g_z_0_xzzzz_xyyzzz[k];

                g_z_0_xxzzzz_yzzzz[k] = -g_z_0_xzzzz_yzzzz[k] * ab_x + g_z_0_xzzzz_xyzzzz[k];

                g_z_0_xxzzzz_zzzzz[k] = -g_z_0_xzzzz_zzzzz[k] * ab_x + g_z_0_xzzzz_xzzzzz[k];
            }

            /// Set up 1491-1512 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1491 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1492 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1493 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1494 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1495 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1496 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1497 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1498 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1499 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1500 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1501 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1502 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1503 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1504 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1505 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1506 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1507 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1508 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1509 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1510 * ccomps * dcomps);

            auto g_z_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyy_xxxxx, g_z_0_xyyyyy_xxxxy, g_z_0_xyyyyy_xxxxz, g_z_0_xyyyyy_xxxyy, g_z_0_xyyyyy_xxxyz, g_z_0_xyyyyy_xxxzz, g_z_0_xyyyyy_xxyyy, g_z_0_xyyyyy_xxyyz, g_z_0_xyyyyy_xxyzz, g_z_0_xyyyyy_xxzzz, g_z_0_xyyyyy_xyyyy, g_z_0_xyyyyy_xyyyz, g_z_0_xyyyyy_xyyzz, g_z_0_xyyyyy_xyzzz, g_z_0_xyyyyy_xzzzz, g_z_0_xyyyyy_yyyyy, g_z_0_xyyyyy_yyyyz, g_z_0_xyyyyy_yyyzz, g_z_0_xyyyyy_yyzzz, g_z_0_xyyyyy_yzzzz, g_z_0_xyyyyy_zzzzz, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxxxx[k] = -g_z_0_yyyyy_xxxxx[k] * ab_x + g_z_0_yyyyy_xxxxxx[k];

                g_z_0_xyyyyy_xxxxy[k] = -g_z_0_yyyyy_xxxxy[k] * ab_x + g_z_0_yyyyy_xxxxxy[k];

                g_z_0_xyyyyy_xxxxz[k] = -g_z_0_yyyyy_xxxxz[k] * ab_x + g_z_0_yyyyy_xxxxxz[k];

                g_z_0_xyyyyy_xxxyy[k] = -g_z_0_yyyyy_xxxyy[k] * ab_x + g_z_0_yyyyy_xxxxyy[k];

                g_z_0_xyyyyy_xxxyz[k] = -g_z_0_yyyyy_xxxyz[k] * ab_x + g_z_0_yyyyy_xxxxyz[k];

                g_z_0_xyyyyy_xxxzz[k] = -g_z_0_yyyyy_xxxzz[k] * ab_x + g_z_0_yyyyy_xxxxzz[k];

                g_z_0_xyyyyy_xxyyy[k] = -g_z_0_yyyyy_xxyyy[k] * ab_x + g_z_0_yyyyy_xxxyyy[k];

                g_z_0_xyyyyy_xxyyz[k] = -g_z_0_yyyyy_xxyyz[k] * ab_x + g_z_0_yyyyy_xxxyyz[k];

                g_z_0_xyyyyy_xxyzz[k] = -g_z_0_yyyyy_xxyzz[k] * ab_x + g_z_0_yyyyy_xxxyzz[k];

                g_z_0_xyyyyy_xxzzz[k] = -g_z_0_yyyyy_xxzzz[k] * ab_x + g_z_0_yyyyy_xxxzzz[k];

                g_z_0_xyyyyy_xyyyy[k] = -g_z_0_yyyyy_xyyyy[k] * ab_x + g_z_0_yyyyy_xxyyyy[k];

                g_z_0_xyyyyy_xyyyz[k] = -g_z_0_yyyyy_xyyyz[k] * ab_x + g_z_0_yyyyy_xxyyyz[k];

                g_z_0_xyyyyy_xyyzz[k] = -g_z_0_yyyyy_xyyzz[k] * ab_x + g_z_0_yyyyy_xxyyzz[k];

                g_z_0_xyyyyy_xyzzz[k] = -g_z_0_yyyyy_xyzzz[k] * ab_x + g_z_0_yyyyy_xxyzzz[k];

                g_z_0_xyyyyy_xzzzz[k] = -g_z_0_yyyyy_xzzzz[k] * ab_x + g_z_0_yyyyy_xxzzzz[k];

                g_z_0_xyyyyy_yyyyy[k] = -g_z_0_yyyyy_yyyyy[k] * ab_x + g_z_0_yyyyy_xyyyyy[k];

                g_z_0_xyyyyy_yyyyz[k] = -g_z_0_yyyyy_yyyyz[k] * ab_x + g_z_0_yyyyy_xyyyyz[k];

                g_z_0_xyyyyy_yyyzz[k] = -g_z_0_yyyyy_yyyzz[k] * ab_x + g_z_0_yyyyy_xyyyzz[k];

                g_z_0_xyyyyy_yyzzz[k] = -g_z_0_yyyyy_yyzzz[k] * ab_x + g_z_0_yyyyy_xyyzzz[k];

                g_z_0_xyyyyy_yzzzz[k] = -g_z_0_yyyyy_yzzzz[k] * ab_x + g_z_0_yyyyy_xyzzzz[k];

                g_z_0_xyyyyy_zzzzz[k] = -g_z_0_yyyyy_zzzzz[k] * ab_x + g_z_0_yyyyy_xzzzzz[k];
            }

            /// Set up 1512-1533 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1512 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1513 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1514 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1515 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1516 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1517 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1518 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1519 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1520 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1521 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1522 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1523 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1524 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1525 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1526 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1527 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1528 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1529 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1530 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1531 * ccomps * dcomps);

            auto g_z_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1532 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyz_xxxxx, g_z_0_xyyyyz_xxxxy, g_z_0_xyyyyz_xxxxz, g_z_0_xyyyyz_xxxyy, g_z_0_xyyyyz_xxxyz, g_z_0_xyyyyz_xxxzz, g_z_0_xyyyyz_xxyyy, g_z_0_xyyyyz_xxyyz, g_z_0_xyyyyz_xxyzz, g_z_0_xyyyyz_xxzzz, g_z_0_xyyyyz_xyyyy, g_z_0_xyyyyz_xyyyz, g_z_0_xyyyyz_xyyzz, g_z_0_xyyyyz_xyzzz, g_z_0_xyyyyz_xzzzz, g_z_0_xyyyyz_yyyyy, g_z_0_xyyyyz_yyyyz, g_z_0_xyyyyz_yyyzz, g_z_0_xyyyyz_yyzzz, g_z_0_xyyyyz_yzzzz, g_z_0_xyyyyz_zzzzz, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxxxx[k] = -g_z_0_yyyyz_xxxxx[k] * ab_x + g_z_0_yyyyz_xxxxxx[k];

                g_z_0_xyyyyz_xxxxy[k] = -g_z_0_yyyyz_xxxxy[k] * ab_x + g_z_0_yyyyz_xxxxxy[k];

                g_z_0_xyyyyz_xxxxz[k] = -g_z_0_yyyyz_xxxxz[k] * ab_x + g_z_0_yyyyz_xxxxxz[k];

                g_z_0_xyyyyz_xxxyy[k] = -g_z_0_yyyyz_xxxyy[k] * ab_x + g_z_0_yyyyz_xxxxyy[k];

                g_z_0_xyyyyz_xxxyz[k] = -g_z_0_yyyyz_xxxyz[k] * ab_x + g_z_0_yyyyz_xxxxyz[k];

                g_z_0_xyyyyz_xxxzz[k] = -g_z_0_yyyyz_xxxzz[k] * ab_x + g_z_0_yyyyz_xxxxzz[k];

                g_z_0_xyyyyz_xxyyy[k] = -g_z_0_yyyyz_xxyyy[k] * ab_x + g_z_0_yyyyz_xxxyyy[k];

                g_z_0_xyyyyz_xxyyz[k] = -g_z_0_yyyyz_xxyyz[k] * ab_x + g_z_0_yyyyz_xxxyyz[k];

                g_z_0_xyyyyz_xxyzz[k] = -g_z_0_yyyyz_xxyzz[k] * ab_x + g_z_0_yyyyz_xxxyzz[k];

                g_z_0_xyyyyz_xxzzz[k] = -g_z_0_yyyyz_xxzzz[k] * ab_x + g_z_0_yyyyz_xxxzzz[k];

                g_z_0_xyyyyz_xyyyy[k] = -g_z_0_yyyyz_xyyyy[k] * ab_x + g_z_0_yyyyz_xxyyyy[k];

                g_z_0_xyyyyz_xyyyz[k] = -g_z_0_yyyyz_xyyyz[k] * ab_x + g_z_0_yyyyz_xxyyyz[k];

                g_z_0_xyyyyz_xyyzz[k] = -g_z_0_yyyyz_xyyzz[k] * ab_x + g_z_0_yyyyz_xxyyzz[k];

                g_z_0_xyyyyz_xyzzz[k] = -g_z_0_yyyyz_xyzzz[k] * ab_x + g_z_0_yyyyz_xxyzzz[k];

                g_z_0_xyyyyz_xzzzz[k] = -g_z_0_yyyyz_xzzzz[k] * ab_x + g_z_0_yyyyz_xxzzzz[k];

                g_z_0_xyyyyz_yyyyy[k] = -g_z_0_yyyyz_yyyyy[k] * ab_x + g_z_0_yyyyz_xyyyyy[k];

                g_z_0_xyyyyz_yyyyz[k] = -g_z_0_yyyyz_yyyyz[k] * ab_x + g_z_0_yyyyz_xyyyyz[k];

                g_z_0_xyyyyz_yyyzz[k] = -g_z_0_yyyyz_yyyzz[k] * ab_x + g_z_0_yyyyz_xyyyzz[k];

                g_z_0_xyyyyz_yyzzz[k] = -g_z_0_yyyyz_yyzzz[k] * ab_x + g_z_0_yyyyz_xyyzzz[k];

                g_z_0_xyyyyz_yzzzz[k] = -g_z_0_yyyyz_yzzzz[k] * ab_x + g_z_0_yyyyz_xyzzzz[k];

                g_z_0_xyyyyz_zzzzz[k] = -g_z_0_yyyyz_zzzzz[k] * ab_x + g_z_0_yyyyz_xzzzzz[k];
            }

            /// Set up 1533-1554 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1533 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1534 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1535 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1536 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1537 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1538 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1539 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1540 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1541 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1542 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1543 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1544 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1545 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1546 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1547 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1548 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1549 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1550 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1551 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1552 * ccomps * dcomps);

            auto g_z_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1553 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyzz_xxxxx, g_z_0_xyyyzz_xxxxy, g_z_0_xyyyzz_xxxxz, g_z_0_xyyyzz_xxxyy, g_z_0_xyyyzz_xxxyz, g_z_0_xyyyzz_xxxzz, g_z_0_xyyyzz_xxyyy, g_z_0_xyyyzz_xxyyz, g_z_0_xyyyzz_xxyzz, g_z_0_xyyyzz_xxzzz, g_z_0_xyyyzz_xyyyy, g_z_0_xyyyzz_xyyyz, g_z_0_xyyyzz_xyyzz, g_z_0_xyyyzz_xyzzz, g_z_0_xyyyzz_xzzzz, g_z_0_xyyyzz_yyyyy, g_z_0_xyyyzz_yyyyz, g_z_0_xyyyzz_yyyzz, g_z_0_xyyyzz_yyzzz, g_z_0_xyyyzz_yzzzz, g_z_0_xyyyzz_zzzzz, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxxxx[k] = -g_z_0_yyyzz_xxxxx[k] * ab_x + g_z_0_yyyzz_xxxxxx[k];

                g_z_0_xyyyzz_xxxxy[k] = -g_z_0_yyyzz_xxxxy[k] * ab_x + g_z_0_yyyzz_xxxxxy[k];

                g_z_0_xyyyzz_xxxxz[k] = -g_z_0_yyyzz_xxxxz[k] * ab_x + g_z_0_yyyzz_xxxxxz[k];

                g_z_0_xyyyzz_xxxyy[k] = -g_z_0_yyyzz_xxxyy[k] * ab_x + g_z_0_yyyzz_xxxxyy[k];

                g_z_0_xyyyzz_xxxyz[k] = -g_z_0_yyyzz_xxxyz[k] * ab_x + g_z_0_yyyzz_xxxxyz[k];

                g_z_0_xyyyzz_xxxzz[k] = -g_z_0_yyyzz_xxxzz[k] * ab_x + g_z_0_yyyzz_xxxxzz[k];

                g_z_0_xyyyzz_xxyyy[k] = -g_z_0_yyyzz_xxyyy[k] * ab_x + g_z_0_yyyzz_xxxyyy[k];

                g_z_0_xyyyzz_xxyyz[k] = -g_z_0_yyyzz_xxyyz[k] * ab_x + g_z_0_yyyzz_xxxyyz[k];

                g_z_0_xyyyzz_xxyzz[k] = -g_z_0_yyyzz_xxyzz[k] * ab_x + g_z_0_yyyzz_xxxyzz[k];

                g_z_0_xyyyzz_xxzzz[k] = -g_z_0_yyyzz_xxzzz[k] * ab_x + g_z_0_yyyzz_xxxzzz[k];

                g_z_0_xyyyzz_xyyyy[k] = -g_z_0_yyyzz_xyyyy[k] * ab_x + g_z_0_yyyzz_xxyyyy[k];

                g_z_0_xyyyzz_xyyyz[k] = -g_z_0_yyyzz_xyyyz[k] * ab_x + g_z_0_yyyzz_xxyyyz[k];

                g_z_0_xyyyzz_xyyzz[k] = -g_z_0_yyyzz_xyyzz[k] * ab_x + g_z_0_yyyzz_xxyyzz[k];

                g_z_0_xyyyzz_xyzzz[k] = -g_z_0_yyyzz_xyzzz[k] * ab_x + g_z_0_yyyzz_xxyzzz[k];

                g_z_0_xyyyzz_xzzzz[k] = -g_z_0_yyyzz_xzzzz[k] * ab_x + g_z_0_yyyzz_xxzzzz[k];

                g_z_0_xyyyzz_yyyyy[k] = -g_z_0_yyyzz_yyyyy[k] * ab_x + g_z_0_yyyzz_xyyyyy[k];

                g_z_0_xyyyzz_yyyyz[k] = -g_z_0_yyyzz_yyyyz[k] * ab_x + g_z_0_yyyzz_xyyyyz[k];

                g_z_0_xyyyzz_yyyzz[k] = -g_z_0_yyyzz_yyyzz[k] * ab_x + g_z_0_yyyzz_xyyyzz[k];

                g_z_0_xyyyzz_yyzzz[k] = -g_z_0_yyyzz_yyzzz[k] * ab_x + g_z_0_yyyzz_xyyzzz[k];

                g_z_0_xyyyzz_yzzzz[k] = -g_z_0_yyyzz_yzzzz[k] * ab_x + g_z_0_yyyzz_xyzzzz[k];

                g_z_0_xyyyzz_zzzzz[k] = -g_z_0_yyyzz_zzzzz[k] * ab_x + g_z_0_yyyzz_xzzzzz[k];
            }

            /// Set up 1554-1575 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1554 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1555 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1556 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1557 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1558 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1559 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1560 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1561 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1562 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1563 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1564 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1565 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1566 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1567 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1574 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzzz_xxxxx, g_z_0_xyyzzz_xxxxy, g_z_0_xyyzzz_xxxxz, g_z_0_xyyzzz_xxxyy, g_z_0_xyyzzz_xxxyz, g_z_0_xyyzzz_xxxzz, g_z_0_xyyzzz_xxyyy, g_z_0_xyyzzz_xxyyz, g_z_0_xyyzzz_xxyzz, g_z_0_xyyzzz_xxzzz, g_z_0_xyyzzz_xyyyy, g_z_0_xyyzzz_xyyyz, g_z_0_xyyzzz_xyyzz, g_z_0_xyyzzz_xyzzz, g_z_0_xyyzzz_xzzzz, g_z_0_xyyzzz_yyyyy, g_z_0_xyyzzz_yyyyz, g_z_0_xyyzzz_yyyzz, g_z_0_xyyzzz_yyzzz, g_z_0_xyyzzz_yzzzz, g_z_0_xyyzzz_zzzzz, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxxxx[k] = -g_z_0_yyzzz_xxxxx[k] * ab_x + g_z_0_yyzzz_xxxxxx[k];

                g_z_0_xyyzzz_xxxxy[k] = -g_z_0_yyzzz_xxxxy[k] * ab_x + g_z_0_yyzzz_xxxxxy[k];

                g_z_0_xyyzzz_xxxxz[k] = -g_z_0_yyzzz_xxxxz[k] * ab_x + g_z_0_yyzzz_xxxxxz[k];

                g_z_0_xyyzzz_xxxyy[k] = -g_z_0_yyzzz_xxxyy[k] * ab_x + g_z_0_yyzzz_xxxxyy[k];

                g_z_0_xyyzzz_xxxyz[k] = -g_z_0_yyzzz_xxxyz[k] * ab_x + g_z_0_yyzzz_xxxxyz[k];

                g_z_0_xyyzzz_xxxzz[k] = -g_z_0_yyzzz_xxxzz[k] * ab_x + g_z_0_yyzzz_xxxxzz[k];

                g_z_0_xyyzzz_xxyyy[k] = -g_z_0_yyzzz_xxyyy[k] * ab_x + g_z_0_yyzzz_xxxyyy[k];

                g_z_0_xyyzzz_xxyyz[k] = -g_z_0_yyzzz_xxyyz[k] * ab_x + g_z_0_yyzzz_xxxyyz[k];

                g_z_0_xyyzzz_xxyzz[k] = -g_z_0_yyzzz_xxyzz[k] * ab_x + g_z_0_yyzzz_xxxyzz[k];

                g_z_0_xyyzzz_xxzzz[k] = -g_z_0_yyzzz_xxzzz[k] * ab_x + g_z_0_yyzzz_xxxzzz[k];

                g_z_0_xyyzzz_xyyyy[k] = -g_z_0_yyzzz_xyyyy[k] * ab_x + g_z_0_yyzzz_xxyyyy[k];

                g_z_0_xyyzzz_xyyyz[k] = -g_z_0_yyzzz_xyyyz[k] * ab_x + g_z_0_yyzzz_xxyyyz[k];

                g_z_0_xyyzzz_xyyzz[k] = -g_z_0_yyzzz_xyyzz[k] * ab_x + g_z_0_yyzzz_xxyyzz[k];

                g_z_0_xyyzzz_xyzzz[k] = -g_z_0_yyzzz_xyzzz[k] * ab_x + g_z_0_yyzzz_xxyzzz[k];

                g_z_0_xyyzzz_xzzzz[k] = -g_z_0_yyzzz_xzzzz[k] * ab_x + g_z_0_yyzzz_xxzzzz[k];

                g_z_0_xyyzzz_yyyyy[k] = -g_z_0_yyzzz_yyyyy[k] * ab_x + g_z_0_yyzzz_xyyyyy[k];

                g_z_0_xyyzzz_yyyyz[k] = -g_z_0_yyzzz_yyyyz[k] * ab_x + g_z_0_yyzzz_xyyyyz[k];

                g_z_0_xyyzzz_yyyzz[k] = -g_z_0_yyzzz_yyyzz[k] * ab_x + g_z_0_yyzzz_xyyyzz[k];

                g_z_0_xyyzzz_yyzzz[k] = -g_z_0_yyzzz_yyzzz[k] * ab_x + g_z_0_yyzzz_xyyzzz[k];

                g_z_0_xyyzzz_yzzzz[k] = -g_z_0_yyzzz_yzzzz[k] * ab_x + g_z_0_yyzzz_xyzzzz[k];

                g_z_0_xyyzzz_zzzzz[k] = -g_z_0_yyzzz_zzzzz[k] * ab_x + g_z_0_yyzzz_xzzzzz[k];
            }

            /// Set up 1575-1596 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1576 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1577 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1578 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1579 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1580 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1581 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1582 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1583 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1589 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1590 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1591 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1592 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1593 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1594 * ccomps * dcomps);

            auto g_z_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1595 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzzz_xxxxx, g_z_0_xyzzzz_xxxxy, g_z_0_xyzzzz_xxxxz, g_z_0_xyzzzz_xxxyy, g_z_0_xyzzzz_xxxyz, g_z_0_xyzzzz_xxxzz, g_z_0_xyzzzz_xxyyy, g_z_0_xyzzzz_xxyyz, g_z_0_xyzzzz_xxyzz, g_z_0_xyzzzz_xxzzz, g_z_0_xyzzzz_xyyyy, g_z_0_xyzzzz_xyyyz, g_z_0_xyzzzz_xyyzz, g_z_0_xyzzzz_xyzzz, g_z_0_xyzzzz_xzzzz, g_z_0_xyzzzz_yyyyy, g_z_0_xyzzzz_yyyyz, g_z_0_xyzzzz_yyyzz, g_z_0_xyzzzz_yyzzz, g_z_0_xyzzzz_yzzzz, g_z_0_xyzzzz_zzzzz, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxxxx[k] = -g_z_0_yzzzz_xxxxx[k] * ab_x + g_z_0_yzzzz_xxxxxx[k];

                g_z_0_xyzzzz_xxxxy[k] = -g_z_0_yzzzz_xxxxy[k] * ab_x + g_z_0_yzzzz_xxxxxy[k];

                g_z_0_xyzzzz_xxxxz[k] = -g_z_0_yzzzz_xxxxz[k] * ab_x + g_z_0_yzzzz_xxxxxz[k];

                g_z_0_xyzzzz_xxxyy[k] = -g_z_0_yzzzz_xxxyy[k] * ab_x + g_z_0_yzzzz_xxxxyy[k];

                g_z_0_xyzzzz_xxxyz[k] = -g_z_0_yzzzz_xxxyz[k] * ab_x + g_z_0_yzzzz_xxxxyz[k];

                g_z_0_xyzzzz_xxxzz[k] = -g_z_0_yzzzz_xxxzz[k] * ab_x + g_z_0_yzzzz_xxxxzz[k];

                g_z_0_xyzzzz_xxyyy[k] = -g_z_0_yzzzz_xxyyy[k] * ab_x + g_z_0_yzzzz_xxxyyy[k];

                g_z_0_xyzzzz_xxyyz[k] = -g_z_0_yzzzz_xxyyz[k] * ab_x + g_z_0_yzzzz_xxxyyz[k];

                g_z_0_xyzzzz_xxyzz[k] = -g_z_0_yzzzz_xxyzz[k] * ab_x + g_z_0_yzzzz_xxxyzz[k];

                g_z_0_xyzzzz_xxzzz[k] = -g_z_0_yzzzz_xxzzz[k] * ab_x + g_z_0_yzzzz_xxxzzz[k];

                g_z_0_xyzzzz_xyyyy[k] = -g_z_0_yzzzz_xyyyy[k] * ab_x + g_z_0_yzzzz_xxyyyy[k];

                g_z_0_xyzzzz_xyyyz[k] = -g_z_0_yzzzz_xyyyz[k] * ab_x + g_z_0_yzzzz_xxyyyz[k];

                g_z_0_xyzzzz_xyyzz[k] = -g_z_0_yzzzz_xyyzz[k] * ab_x + g_z_0_yzzzz_xxyyzz[k];

                g_z_0_xyzzzz_xyzzz[k] = -g_z_0_yzzzz_xyzzz[k] * ab_x + g_z_0_yzzzz_xxyzzz[k];

                g_z_0_xyzzzz_xzzzz[k] = -g_z_0_yzzzz_xzzzz[k] * ab_x + g_z_0_yzzzz_xxzzzz[k];

                g_z_0_xyzzzz_yyyyy[k] = -g_z_0_yzzzz_yyyyy[k] * ab_x + g_z_0_yzzzz_xyyyyy[k];

                g_z_0_xyzzzz_yyyyz[k] = -g_z_0_yzzzz_yyyyz[k] * ab_x + g_z_0_yzzzz_xyyyyz[k];

                g_z_0_xyzzzz_yyyzz[k] = -g_z_0_yzzzz_yyyzz[k] * ab_x + g_z_0_yzzzz_xyyyzz[k];

                g_z_0_xyzzzz_yyzzz[k] = -g_z_0_yzzzz_yyzzz[k] * ab_x + g_z_0_yzzzz_xyyzzz[k];

                g_z_0_xyzzzz_yzzzz[k] = -g_z_0_yzzzz_yzzzz[k] * ab_x + g_z_0_yzzzz_xyzzzz[k];

                g_z_0_xyzzzz_zzzzz[k] = -g_z_0_yzzzz_zzzzz[k] * ab_x + g_z_0_yzzzz_xzzzzz[k];
            }

            /// Set up 1596-1617 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1612 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1613 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1614 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1615 * ccomps * dcomps);

            auto g_z_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1616 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzzz_xxxxx, g_z_0_xzzzzz_xxxxy, g_z_0_xzzzzz_xxxxz, g_z_0_xzzzzz_xxxyy, g_z_0_xzzzzz_xxxyz, g_z_0_xzzzzz_xxxzz, g_z_0_xzzzzz_xxyyy, g_z_0_xzzzzz_xxyyz, g_z_0_xzzzzz_xxyzz, g_z_0_xzzzzz_xxzzz, g_z_0_xzzzzz_xyyyy, g_z_0_xzzzzz_xyyyz, g_z_0_xzzzzz_xyyzz, g_z_0_xzzzzz_xyzzz, g_z_0_xzzzzz_xzzzz, g_z_0_xzzzzz_yyyyy, g_z_0_xzzzzz_yyyyz, g_z_0_xzzzzz_yyyzz, g_z_0_xzzzzz_yyzzz, g_z_0_xzzzzz_yzzzz, g_z_0_xzzzzz_zzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxxxx[k] = -g_z_0_zzzzz_xxxxx[k] * ab_x + g_z_0_zzzzz_xxxxxx[k];

                g_z_0_xzzzzz_xxxxy[k] = -g_z_0_zzzzz_xxxxy[k] * ab_x + g_z_0_zzzzz_xxxxxy[k];

                g_z_0_xzzzzz_xxxxz[k] = -g_z_0_zzzzz_xxxxz[k] * ab_x + g_z_0_zzzzz_xxxxxz[k];

                g_z_0_xzzzzz_xxxyy[k] = -g_z_0_zzzzz_xxxyy[k] * ab_x + g_z_0_zzzzz_xxxxyy[k];

                g_z_0_xzzzzz_xxxyz[k] = -g_z_0_zzzzz_xxxyz[k] * ab_x + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_xzzzzz_xxxzz[k] = -g_z_0_zzzzz_xxxzz[k] * ab_x + g_z_0_zzzzz_xxxxzz[k];

                g_z_0_xzzzzz_xxyyy[k] = -g_z_0_zzzzz_xxyyy[k] * ab_x + g_z_0_zzzzz_xxxyyy[k];

                g_z_0_xzzzzz_xxyyz[k] = -g_z_0_zzzzz_xxyyz[k] * ab_x + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_xzzzzz_xxyzz[k] = -g_z_0_zzzzz_xxyzz[k] * ab_x + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_xzzzzz_xxzzz[k] = -g_z_0_zzzzz_xxzzz[k] * ab_x + g_z_0_zzzzz_xxxzzz[k];

                g_z_0_xzzzzz_xyyyy[k] = -g_z_0_zzzzz_xyyyy[k] * ab_x + g_z_0_zzzzz_xxyyyy[k];

                g_z_0_xzzzzz_xyyyz[k] = -g_z_0_zzzzz_xyyyz[k] * ab_x + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_xzzzzz_xyyzz[k] = -g_z_0_zzzzz_xyyzz[k] * ab_x + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_xzzzzz_xyzzz[k] = -g_z_0_zzzzz_xyzzz[k] * ab_x + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_xzzzzz_xzzzz[k] = -g_z_0_zzzzz_xzzzz[k] * ab_x + g_z_0_zzzzz_xxzzzz[k];

                g_z_0_xzzzzz_yyyyy[k] = -g_z_0_zzzzz_yyyyy[k] * ab_x + g_z_0_zzzzz_xyyyyy[k];

                g_z_0_xzzzzz_yyyyz[k] = -g_z_0_zzzzz_yyyyz[k] * ab_x + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_xzzzzz_yyyzz[k] = -g_z_0_zzzzz_yyyzz[k] * ab_x + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_xzzzzz_yyzzz[k] = -g_z_0_zzzzz_yyzzz[k] * ab_x + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_xzzzzz_yzzzz[k] = -g_z_0_zzzzz_yzzzz[k] * ab_x + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_xzzzzz_zzzzz[k] = -g_z_0_zzzzz_zzzzz[k] * ab_x + g_z_0_zzzzz_xzzzzz[k];
            }

            /// Set up 1617-1638 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1617 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1618 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1619 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1620 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1621 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1622 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1623 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1624 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1625 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1626 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1627 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1628 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1629 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1630 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1631 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1632 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1633 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1634 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1635 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1636 * ccomps * dcomps);

            auto g_z_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1637 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_zzzzz, g_z_0_yyyyyy_xxxxx, g_z_0_yyyyyy_xxxxy, g_z_0_yyyyyy_xxxxz, g_z_0_yyyyyy_xxxyy, g_z_0_yyyyyy_xxxyz, g_z_0_yyyyyy_xxxzz, g_z_0_yyyyyy_xxyyy, g_z_0_yyyyyy_xxyyz, g_z_0_yyyyyy_xxyzz, g_z_0_yyyyyy_xxzzz, g_z_0_yyyyyy_xyyyy, g_z_0_yyyyyy_xyyyz, g_z_0_yyyyyy_xyyzz, g_z_0_yyyyyy_xyzzz, g_z_0_yyyyyy_xzzzz, g_z_0_yyyyyy_yyyyy, g_z_0_yyyyyy_yyyyz, g_z_0_yyyyyy_yyyzz, g_z_0_yyyyyy_yyzzz, g_z_0_yyyyyy_yzzzz, g_z_0_yyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxxxx[k] = -g_z_0_yyyyy_xxxxx[k] * ab_y + g_z_0_yyyyy_xxxxxy[k];

                g_z_0_yyyyyy_xxxxy[k] = -g_z_0_yyyyy_xxxxy[k] * ab_y + g_z_0_yyyyy_xxxxyy[k];

                g_z_0_yyyyyy_xxxxz[k] = -g_z_0_yyyyy_xxxxz[k] * ab_y + g_z_0_yyyyy_xxxxyz[k];

                g_z_0_yyyyyy_xxxyy[k] = -g_z_0_yyyyy_xxxyy[k] * ab_y + g_z_0_yyyyy_xxxyyy[k];

                g_z_0_yyyyyy_xxxyz[k] = -g_z_0_yyyyy_xxxyz[k] * ab_y + g_z_0_yyyyy_xxxyyz[k];

                g_z_0_yyyyyy_xxxzz[k] = -g_z_0_yyyyy_xxxzz[k] * ab_y + g_z_0_yyyyy_xxxyzz[k];

                g_z_0_yyyyyy_xxyyy[k] = -g_z_0_yyyyy_xxyyy[k] * ab_y + g_z_0_yyyyy_xxyyyy[k];

                g_z_0_yyyyyy_xxyyz[k] = -g_z_0_yyyyy_xxyyz[k] * ab_y + g_z_0_yyyyy_xxyyyz[k];

                g_z_0_yyyyyy_xxyzz[k] = -g_z_0_yyyyy_xxyzz[k] * ab_y + g_z_0_yyyyy_xxyyzz[k];

                g_z_0_yyyyyy_xxzzz[k] = -g_z_0_yyyyy_xxzzz[k] * ab_y + g_z_0_yyyyy_xxyzzz[k];

                g_z_0_yyyyyy_xyyyy[k] = -g_z_0_yyyyy_xyyyy[k] * ab_y + g_z_0_yyyyy_xyyyyy[k];

                g_z_0_yyyyyy_xyyyz[k] = -g_z_0_yyyyy_xyyyz[k] * ab_y + g_z_0_yyyyy_xyyyyz[k];

                g_z_0_yyyyyy_xyyzz[k] = -g_z_0_yyyyy_xyyzz[k] * ab_y + g_z_0_yyyyy_xyyyzz[k];

                g_z_0_yyyyyy_xyzzz[k] = -g_z_0_yyyyy_xyzzz[k] * ab_y + g_z_0_yyyyy_xyyzzz[k];

                g_z_0_yyyyyy_xzzzz[k] = -g_z_0_yyyyy_xzzzz[k] * ab_y + g_z_0_yyyyy_xyzzzz[k];

                g_z_0_yyyyyy_yyyyy[k] = -g_z_0_yyyyy_yyyyy[k] * ab_y + g_z_0_yyyyy_yyyyyy[k];

                g_z_0_yyyyyy_yyyyz[k] = -g_z_0_yyyyy_yyyyz[k] * ab_y + g_z_0_yyyyy_yyyyyz[k];

                g_z_0_yyyyyy_yyyzz[k] = -g_z_0_yyyyy_yyyzz[k] * ab_y + g_z_0_yyyyy_yyyyzz[k];

                g_z_0_yyyyyy_yyzzz[k] = -g_z_0_yyyyy_yyzzz[k] * ab_y + g_z_0_yyyyy_yyyzzz[k];

                g_z_0_yyyyyy_yzzzz[k] = -g_z_0_yyyyy_yzzzz[k] * ab_y + g_z_0_yyyyy_yyzzzz[k];

                g_z_0_yyyyyy_zzzzz[k] = -g_z_0_yyyyy_zzzzz[k] * ab_y + g_z_0_yyyyy_yzzzzz[k];
            }

            /// Set up 1638-1659 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1638 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1639 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1640 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1641 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1642 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1643 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1644 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1645 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1646 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1647 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1648 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1649 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1650 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1651 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1652 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1653 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1654 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1655 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1656 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1657 * ccomps * dcomps);

            auto g_z_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1658 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyyz_xxxxx, g_z_0_yyyyyz_xxxxy, g_z_0_yyyyyz_xxxxz, g_z_0_yyyyyz_xxxyy, g_z_0_yyyyyz_xxxyz, g_z_0_yyyyyz_xxxzz, g_z_0_yyyyyz_xxyyy, g_z_0_yyyyyz_xxyyz, g_z_0_yyyyyz_xxyzz, g_z_0_yyyyyz_xxzzz, g_z_0_yyyyyz_xyyyy, g_z_0_yyyyyz_xyyyz, g_z_0_yyyyyz_xyyzz, g_z_0_yyyyyz_xyzzz, g_z_0_yyyyyz_xzzzz, g_z_0_yyyyyz_yyyyy, g_z_0_yyyyyz_yyyyz, g_z_0_yyyyyz_yyyzz, g_z_0_yyyyyz_yyzzz, g_z_0_yyyyyz_yzzzz, g_z_0_yyyyyz_zzzzz, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxxxx[k] = -g_z_0_yyyyz_xxxxx[k] * ab_y + g_z_0_yyyyz_xxxxxy[k];

                g_z_0_yyyyyz_xxxxy[k] = -g_z_0_yyyyz_xxxxy[k] * ab_y + g_z_0_yyyyz_xxxxyy[k];

                g_z_0_yyyyyz_xxxxz[k] = -g_z_0_yyyyz_xxxxz[k] * ab_y + g_z_0_yyyyz_xxxxyz[k];

                g_z_0_yyyyyz_xxxyy[k] = -g_z_0_yyyyz_xxxyy[k] * ab_y + g_z_0_yyyyz_xxxyyy[k];

                g_z_0_yyyyyz_xxxyz[k] = -g_z_0_yyyyz_xxxyz[k] * ab_y + g_z_0_yyyyz_xxxyyz[k];

                g_z_0_yyyyyz_xxxzz[k] = -g_z_0_yyyyz_xxxzz[k] * ab_y + g_z_0_yyyyz_xxxyzz[k];

                g_z_0_yyyyyz_xxyyy[k] = -g_z_0_yyyyz_xxyyy[k] * ab_y + g_z_0_yyyyz_xxyyyy[k];

                g_z_0_yyyyyz_xxyyz[k] = -g_z_0_yyyyz_xxyyz[k] * ab_y + g_z_0_yyyyz_xxyyyz[k];

                g_z_0_yyyyyz_xxyzz[k] = -g_z_0_yyyyz_xxyzz[k] * ab_y + g_z_0_yyyyz_xxyyzz[k];

                g_z_0_yyyyyz_xxzzz[k] = -g_z_0_yyyyz_xxzzz[k] * ab_y + g_z_0_yyyyz_xxyzzz[k];

                g_z_0_yyyyyz_xyyyy[k] = -g_z_0_yyyyz_xyyyy[k] * ab_y + g_z_0_yyyyz_xyyyyy[k];

                g_z_0_yyyyyz_xyyyz[k] = -g_z_0_yyyyz_xyyyz[k] * ab_y + g_z_0_yyyyz_xyyyyz[k];

                g_z_0_yyyyyz_xyyzz[k] = -g_z_0_yyyyz_xyyzz[k] * ab_y + g_z_0_yyyyz_xyyyzz[k];

                g_z_0_yyyyyz_xyzzz[k] = -g_z_0_yyyyz_xyzzz[k] * ab_y + g_z_0_yyyyz_xyyzzz[k];

                g_z_0_yyyyyz_xzzzz[k] = -g_z_0_yyyyz_xzzzz[k] * ab_y + g_z_0_yyyyz_xyzzzz[k];

                g_z_0_yyyyyz_yyyyy[k] = -g_z_0_yyyyz_yyyyy[k] * ab_y + g_z_0_yyyyz_yyyyyy[k];

                g_z_0_yyyyyz_yyyyz[k] = -g_z_0_yyyyz_yyyyz[k] * ab_y + g_z_0_yyyyz_yyyyyz[k];

                g_z_0_yyyyyz_yyyzz[k] = -g_z_0_yyyyz_yyyzz[k] * ab_y + g_z_0_yyyyz_yyyyzz[k];

                g_z_0_yyyyyz_yyzzz[k] = -g_z_0_yyyyz_yyzzz[k] * ab_y + g_z_0_yyyyz_yyyzzz[k];

                g_z_0_yyyyyz_yzzzz[k] = -g_z_0_yyyyz_yzzzz[k] * ab_y + g_z_0_yyyyz_yyzzzz[k];

                g_z_0_yyyyyz_zzzzz[k] = -g_z_0_yyyyz_zzzzz[k] * ab_y + g_z_0_yyyyz_yzzzzz[k];
            }

            /// Set up 1659-1680 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1659 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1660 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1661 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1662 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1663 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1664 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1665 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1666 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1667 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1668 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1669 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1670 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1671 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1672 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1673 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1674 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1675 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1676 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1677 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1678 * ccomps * dcomps);

            auto g_z_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyzz_xxxxx, g_z_0_yyyyzz_xxxxy, g_z_0_yyyyzz_xxxxz, g_z_0_yyyyzz_xxxyy, g_z_0_yyyyzz_xxxyz, g_z_0_yyyyzz_xxxzz, g_z_0_yyyyzz_xxyyy, g_z_0_yyyyzz_xxyyz, g_z_0_yyyyzz_xxyzz, g_z_0_yyyyzz_xxzzz, g_z_0_yyyyzz_xyyyy, g_z_0_yyyyzz_xyyyz, g_z_0_yyyyzz_xyyzz, g_z_0_yyyyzz_xyzzz, g_z_0_yyyyzz_xzzzz, g_z_0_yyyyzz_yyyyy, g_z_0_yyyyzz_yyyyz, g_z_0_yyyyzz_yyyzz, g_z_0_yyyyzz_yyzzz, g_z_0_yyyyzz_yzzzz, g_z_0_yyyyzz_zzzzz, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxxxx[k] = -g_z_0_yyyzz_xxxxx[k] * ab_y + g_z_0_yyyzz_xxxxxy[k];

                g_z_0_yyyyzz_xxxxy[k] = -g_z_0_yyyzz_xxxxy[k] * ab_y + g_z_0_yyyzz_xxxxyy[k];

                g_z_0_yyyyzz_xxxxz[k] = -g_z_0_yyyzz_xxxxz[k] * ab_y + g_z_0_yyyzz_xxxxyz[k];

                g_z_0_yyyyzz_xxxyy[k] = -g_z_0_yyyzz_xxxyy[k] * ab_y + g_z_0_yyyzz_xxxyyy[k];

                g_z_0_yyyyzz_xxxyz[k] = -g_z_0_yyyzz_xxxyz[k] * ab_y + g_z_0_yyyzz_xxxyyz[k];

                g_z_0_yyyyzz_xxxzz[k] = -g_z_0_yyyzz_xxxzz[k] * ab_y + g_z_0_yyyzz_xxxyzz[k];

                g_z_0_yyyyzz_xxyyy[k] = -g_z_0_yyyzz_xxyyy[k] * ab_y + g_z_0_yyyzz_xxyyyy[k];

                g_z_0_yyyyzz_xxyyz[k] = -g_z_0_yyyzz_xxyyz[k] * ab_y + g_z_0_yyyzz_xxyyyz[k];

                g_z_0_yyyyzz_xxyzz[k] = -g_z_0_yyyzz_xxyzz[k] * ab_y + g_z_0_yyyzz_xxyyzz[k];

                g_z_0_yyyyzz_xxzzz[k] = -g_z_0_yyyzz_xxzzz[k] * ab_y + g_z_0_yyyzz_xxyzzz[k];

                g_z_0_yyyyzz_xyyyy[k] = -g_z_0_yyyzz_xyyyy[k] * ab_y + g_z_0_yyyzz_xyyyyy[k];

                g_z_0_yyyyzz_xyyyz[k] = -g_z_0_yyyzz_xyyyz[k] * ab_y + g_z_0_yyyzz_xyyyyz[k];

                g_z_0_yyyyzz_xyyzz[k] = -g_z_0_yyyzz_xyyzz[k] * ab_y + g_z_0_yyyzz_xyyyzz[k];

                g_z_0_yyyyzz_xyzzz[k] = -g_z_0_yyyzz_xyzzz[k] * ab_y + g_z_0_yyyzz_xyyzzz[k];

                g_z_0_yyyyzz_xzzzz[k] = -g_z_0_yyyzz_xzzzz[k] * ab_y + g_z_0_yyyzz_xyzzzz[k];

                g_z_0_yyyyzz_yyyyy[k] = -g_z_0_yyyzz_yyyyy[k] * ab_y + g_z_0_yyyzz_yyyyyy[k];

                g_z_0_yyyyzz_yyyyz[k] = -g_z_0_yyyzz_yyyyz[k] * ab_y + g_z_0_yyyzz_yyyyyz[k];

                g_z_0_yyyyzz_yyyzz[k] = -g_z_0_yyyzz_yyyzz[k] * ab_y + g_z_0_yyyzz_yyyyzz[k];

                g_z_0_yyyyzz_yyzzz[k] = -g_z_0_yyyzz_yyzzz[k] * ab_y + g_z_0_yyyzz_yyyzzz[k];

                g_z_0_yyyyzz_yzzzz[k] = -g_z_0_yyyzz_yzzzz[k] * ab_y + g_z_0_yyyzz_yyzzzz[k];

                g_z_0_yyyyzz_zzzzz[k] = -g_z_0_yyyzz_zzzzz[k] * ab_y + g_z_0_yyyzz_yzzzzz[k];
            }

            /// Set up 1680-1701 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1680 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1681 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1682 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1683 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1684 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1685 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1686 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1687 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1688 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1689 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1690 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1691 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1692 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1693 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1694 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1695 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1696 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1697 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1698 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1699 * ccomps * dcomps);

            auto g_z_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1700 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzzz_xxxxx, g_z_0_yyyzzz_xxxxy, g_z_0_yyyzzz_xxxxz, g_z_0_yyyzzz_xxxyy, g_z_0_yyyzzz_xxxyz, g_z_0_yyyzzz_xxxzz, g_z_0_yyyzzz_xxyyy, g_z_0_yyyzzz_xxyyz, g_z_0_yyyzzz_xxyzz, g_z_0_yyyzzz_xxzzz, g_z_0_yyyzzz_xyyyy, g_z_0_yyyzzz_xyyyz, g_z_0_yyyzzz_xyyzz, g_z_0_yyyzzz_xyzzz, g_z_0_yyyzzz_xzzzz, g_z_0_yyyzzz_yyyyy, g_z_0_yyyzzz_yyyyz, g_z_0_yyyzzz_yyyzz, g_z_0_yyyzzz_yyzzz, g_z_0_yyyzzz_yzzzz, g_z_0_yyyzzz_zzzzz, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxxxx[k] = -g_z_0_yyzzz_xxxxx[k] * ab_y + g_z_0_yyzzz_xxxxxy[k];

                g_z_0_yyyzzz_xxxxy[k] = -g_z_0_yyzzz_xxxxy[k] * ab_y + g_z_0_yyzzz_xxxxyy[k];

                g_z_0_yyyzzz_xxxxz[k] = -g_z_0_yyzzz_xxxxz[k] * ab_y + g_z_0_yyzzz_xxxxyz[k];

                g_z_0_yyyzzz_xxxyy[k] = -g_z_0_yyzzz_xxxyy[k] * ab_y + g_z_0_yyzzz_xxxyyy[k];

                g_z_0_yyyzzz_xxxyz[k] = -g_z_0_yyzzz_xxxyz[k] * ab_y + g_z_0_yyzzz_xxxyyz[k];

                g_z_0_yyyzzz_xxxzz[k] = -g_z_0_yyzzz_xxxzz[k] * ab_y + g_z_0_yyzzz_xxxyzz[k];

                g_z_0_yyyzzz_xxyyy[k] = -g_z_0_yyzzz_xxyyy[k] * ab_y + g_z_0_yyzzz_xxyyyy[k];

                g_z_0_yyyzzz_xxyyz[k] = -g_z_0_yyzzz_xxyyz[k] * ab_y + g_z_0_yyzzz_xxyyyz[k];

                g_z_0_yyyzzz_xxyzz[k] = -g_z_0_yyzzz_xxyzz[k] * ab_y + g_z_0_yyzzz_xxyyzz[k];

                g_z_0_yyyzzz_xxzzz[k] = -g_z_0_yyzzz_xxzzz[k] * ab_y + g_z_0_yyzzz_xxyzzz[k];

                g_z_0_yyyzzz_xyyyy[k] = -g_z_0_yyzzz_xyyyy[k] * ab_y + g_z_0_yyzzz_xyyyyy[k];

                g_z_0_yyyzzz_xyyyz[k] = -g_z_0_yyzzz_xyyyz[k] * ab_y + g_z_0_yyzzz_xyyyyz[k];

                g_z_0_yyyzzz_xyyzz[k] = -g_z_0_yyzzz_xyyzz[k] * ab_y + g_z_0_yyzzz_xyyyzz[k];

                g_z_0_yyyzzz_xyzzz[k] = -g_z_0_yyzzz_xyzzz[k] * ab_y + g_z_0_yyzzz_xyyzzz[k];

                g_z_0_yyyzzz_xzzzz[k] = -g_z_0_yyzzz_xzzzz[k] * ab_y + g_z_0_yyzzz_xyzzzz[k];

                g_z_0_yyyzzz_yyyyy[k] = -g_z_0_yyzzz_yyyyy[k] * ab_y + g_z_0_yyzzz_yyyyyy[k];

                g_z_0_yyyzzz_yyyyz[k] = -g_z_0_yyzzz_yyyyz[k] * ab_y + g_z_0_yyzzz_yyyyyz[k];

                g_z_0_yyyzzz_yyyzz[k] = -g_z_0_yyzzz_yyyzz[k] * ab_y + g_z_0_yyzzz_yyyyzz[k];

                g_z_0_yyyzzz_yyzzz[k] = -g_z_0_yyzzz_yyzzz[k] * ab_y + g_z_0_yyzzz_yyyzzz[k];

                g_z_0_yyyzzz_yzzzz[k] = -g_z_0_yyzzz_yzzzz[k] * ab_y + g_z_0_yyzzz_yyzzzz[k];

                g_z_0_yyyzzz_zzzzz[k] = -g_z_0_yyzzz_zzzzz[k] * ab_y + g_z_0_yyzzz_yzzzzz[k];
            }

            /// Set up 1701-1722 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1701 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1702 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1703 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1704 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1705 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1706 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1707 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1708 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1709 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1710 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1711 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1712 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1713 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1714 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1715 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1716 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1717 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1718 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1719 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1720 * ccomps * dcomps);

            auto g_z_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1721 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzzz_xxxxx, g_z_0_yyzzzz_xxxxy, g_z_0_yyzzzz_xxxxz, g_z_0_yyzzzz_xxxyy, g_z_0_yyzzzz_xxxyz, g_z_0_yyzzzz_xxxzz, g_z_0_yyzzzz_xxyyy, g_z_0_yyzzzz_xxyyz, g_z_0_yyzzzz_xxyzz, g_z_0_yyzzzz_xxzzz, g_z_0_yyzzzz_xyyyy, g_z_0_yyzzzz_xyyyz, g_z_0_yyzzzz_xyyzz, g_z_0_yyzzzz_xyzzz, g_z_0_yyzzzz_xzzzz, g_z_0_yyzzzz_yyyyy, g_z_0_yyzzzz_yyyyz, g_z_0_yyzzzz_yyyzz, g_z_0_yyzzzz_yyzzz, g_z_0_yyzzzz_yzzzz, g_z_0_yyzzzz_zzzzz, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxxxx[k] = -g_z_0_yzzzz_xxxxx[k] * ab_y + g_z_0_yzzzz_xxxxxy[k];

                g_z_0_yyzzzz_xxxxy[k] = -g_z_0_yzzzz_xxxxy[k] * ab_y + g_z_0_yzzzz_xxxxyy[k];

                g_z_0_yyzzzz_xxxxz[k] = -g_z_0_yzzzz_xxxxz[k] * ab_y + g_z_0_yzzzz_xxxxyz[k];

                g_z_0_yyzzzz_xxxyy[k] = -g_z_0_yzzzz_xxxyy[k] * ab_y + g_z_0_yzzzz_xxxyyy[k];

                g_z_0_yyzzzz_xxxyz[k] = -g_z_0_yzzzz_xxxyz[k] * ab_y + g_z_0_yzzzz_xxxyyz[k];

                g_z_0_yyzzzz_xxxzz[k] = -g_z_0_yzzzz_xxxzz[k] * ab_y + g_z_0_yzzzz_xxxyzz[k];

                g_z_0_yyzzzz_xxyyy[k] = -g_z_0_yzzzz_xxyyy[k] * ab_y + g_z_0_yzzzz_xxyyyy[k];

                g_z_0_yyzzzz_xxyyz[k] = -g_z_0_yzzzz_xxyyz[k] * ab_y + g_z_0_yzzzz_xxyyyz[k];

                g_z_0_yyzzzz_xxyzz[k] = -g_z_0_yzzzz_xxyzz[k] * ab_y + g_z_0_yzzzz_xxyyzz[k];

                g_z_0_yyzzzz_xxzzz[k] = -g_z_0_yzzzz_xxzzz[k] * ab_y + g_z_0_yzzzz_xxyzzz[k];

                g_z_0_yyzzzz_xyyyy[k] = -g_z_0_yzzzz_xyyyy[k] * ab_y + g_z_0_yzzzz_xyyyyy[k];

                g_z_0_yyzzzz_xyyyz[k] = -g_z_0_yzzzz_xyyyz[k] * ab_y + g_z_0_yzzzz_xyyyyz[k];

                g_z_0_yyzzzz_xyyzz[k] = -g_z_0_yzzzz_xyyzz[k] * ab_y + g_z_0_yzzzz_xyyyzz[k];

                g_z_0_yyzzzz_xyzzz[k] = -g_z_0_yzzzz_xyzzz[k] * ab_y + g_z_0_yzzzz_xyyzzz[k];

                g_z_0_yyzzzz_xzzzz[k] = -g_z_0_yzzzz_xzzzz[k] * ab_y + g_z_0_yzzzz_xyzzzz[k];

                g_z_0_yyzzzz_yyyyy[k] = -g_z_0_yzzzz_yyyyy[k] * ab_y + g_z_0_yzzzz_yyyyyy[k];

                g_z_0_yyzzzz_yyyyz[k] = -g_z_0_yzzzz_yyyyz[k] * ab_y + g_z_0_yzzzz_yyyyyz[k];

                g_z_0_yyzzzz_yyyzz[k] = -g_z_0_yzzzz_yyyzz[k] * ab_y + g_z_0_yzzzz_yyyyzz[k];

                g_z_0_yyzzzz_yyzzz[k] = -g_z_0_yzzzz_yyzzz[k] * ab_y + g_z_0_yzzzz_yyyzzz[k];

                g_z_0_yyzzzz_yzzzz[k] = -g_z_0_yzzzz_yzzzz[k] * ab_y + g_z_0_yzzzz_yyzzzz[k];

                g_z_0_yyzzzz_zzzzz[k] = -g_z_0_yzzzz_zzzzz[k] * ab_y + g_z_0_yzzzz_yzzzzz[k];
            }

            /// Set up 1722-1743 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1722 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1723 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1724 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1725 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1726 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1727 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1728 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1729 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1730 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1731 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1732 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1733 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1734 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1735 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1736 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1737 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1738 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1739 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1740 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1741 * ccomps * dcomps);

            auto g_z_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1742 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzzz_xxxxx, g_z_0_yzzzzz_xxxxy, g_z_0_yzzzzz_xxxxz, g_z_0_yzzzzz_xxxyy, g_z_0_yzzzzz_xxxyz, g_z_0_yzzzzz_xxxzz, g_z_0_yzzzzz_xxyyy, g_z_0_yzzzzz_xxyyz, g_z_0_yzzzzz_xxyzz, g_z_0_yzzzzz_xxzzz, g_z_0_yzzzzz_xyyyy, g_z_0_yzzzzz_xyyyz, g_z_0_yzzzzz_xyyzz, g_z_0_yzzzzz_xyzzz, g_z_0_yzzzzz_xzzzz, g_z_0_yzzzzz_yyyyy, g_z_0_yzzzzz_yyyyz, g_z_0_yzzzzz_yyyzz, g_z_0_yzzzzz_yyzzz, g_z_0_yzzzzz_yzzzz, g_z_0_yzzzzz_zzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxxxx[k] = -g_z_0_zzzzz_xxxxx[k] * ab_y + g_z_0_zzzzz_xxxxxy[k];

                g_z_0_yzzzzz_xxxxy[k] = -g_z_0_zzzzz_xxxxy[k] * ab_y + g_z_0_zzzzz_xxxxyy[k];

                g_z_0_yzzzzz_xxxxz[k] = -g_z_0_zzzzz_xxxxz[k] * ab_y + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_yzzzzz_xxxyy[k] = -g_z_0_zzzzz_xxxyy[k] * ab_y + g_z_0_zzzzz_xxxyyy[k];

                g_z_0_yzzzzz_xxxyz[k] = -g_z_0_zzzzz_xxxyz[k] * ab_y + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_yzzzzz_xxxzz[k] = -g_z_0_zzzzz_xxxzz[k] * ab_y + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_yzzzzz_xxyyy[k] = -g_z_0_zzzzz_xxyyy[k] * ab_y + g_z_0_zzzzz_xxyyyy[k];

                g_z_0_yzzzzz_xxyyz[k] = -g_z_0_zzzzz_xxyyz[k] * ab_y + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_yzzzzz_xxyzz[k] = -g_z_0_zzzzz_xxyzz[k] * ab_y + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_yzzzzz_xxzzz[k] = -g_z_0_zzzzz_xxzzz[k] * ab_y + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_yzzzzz_xyyyy[k] = -g_z_0_zzzzz_xyyyy[k] * ab_y + g_z_0_zzzzz_xyyyyy[k];

                g_z_0_yzzzzz_xyyyz[k] = -g_z_0_zzzzz_xyyyz[k] * ab_y + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_yzzzzz_xyyzz[k] = -g_z_0_zzzzz_xyyzz[k] * ab_y + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_yzzzzz_xyzzz[k] = -g_z_0_zzzzz_xyzzz[k] * ab_y + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_yzzzzz_xzzzz[k] = -g_z_0_zzzzz_xzzzz[k] * ab_y + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_yzzzzz_yyyyy[k] = -g_z_0_zzzzz_yyyyy[k] * ab_y + g_z_0_zzzzz_yyyyyy[k];

                g_z_0_yzzzzz_yyyyz[k] = -g_z_0_zzzzz_yyyyz[k] * ab_y + g_z_0_zzzzz_yyyyyz[k];

                g_z_0_yzzzzz_yyyzz[k] = -g_z_0_zzzzz_yyyzz[k] * ab_y + g_z_0_zzzzz_yyyyzz[k];

                g_z_0_yzzzzz_yyzzz[k] = -g_z_0_zzzzz_yyzzz[k] * ab_y + g_z_0_zzzzz_yyyzzz[k];

                g_z_0_yzzzzz_yzzzz[k] = -g_z_0_zzzzz_yzzzz[k] * ab_y + g_z_0_zzzzz_yyzzzz[k];

                g_z_0_yzzzzz_zzzzz[k] = -g_z_0_zzzzz_zzzzz[k] * ab_y + g_z_0_zzzzz_yzzzzz[k];
            }

            /// Set up 1743-1764 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1743 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1744 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1745 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1746 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1747 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1748 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1749 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1750 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1751 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1752 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1753 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1754 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1755 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1756 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1757 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1758 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1759 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1760 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1761 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1762 * ccomps * dcomps);

            auto g_z_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1763 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzz, g_z_0_zzzzz_zzzzzz, g_z_0_zzzzzz_xxxxx, g_z_0_zzzzzz_xxxxy, g_z_0_zzzzzz_xxxxz, g_z_0_zzzzzz_xxxyy, g_z_0_zzzzzz_xxxyz, g_z_0_zzzzzz_xxxzz, g_z_0_zzzzzz_xxyyy, g_z_0_zzzzzz_xxyyz, g_z_0_zzzzzz_xxyzz, g_z_0_zzzzzz_xxzzz, g_z_0_zzzzzz_xyyyy, g_z_0_zzzzzz_xyyyz, g_z_0_zzzzzz_xyyzz, g_z_0_zzzzzz_xyzzz, g_z_0_zzzzzz_xzzzz, g_z_0_zzzzzz_yyyyy, g_z_0_zzzzzz_yyyyz, g_z_0_zzzzzz_yyyzz, g_z_0_zzzzzz_yyzzz, g_z_0_zzzzzz_yzzzz, g_z_0_zzzzzz_zzzzz, g_zzzzz_xxxxx, g_zzzzz_xxxxy, g_zzzzz_xxxxz, g_zzzzz_xxxyy, g_zzzzz_xxxyz, g_zzzzz_xxxzz, g_zzzzz_xxyyy, g_zzzzz_xxyyz, g_zzzzz_xxyzz, g_zzzzz_xxzzz, g_zzzzz_xyyyy, g_zzzzz_xyyyz, g_zzzzz_xyyzz, g_zzzzz_xyzzz, g_zzzzz_xzzzz, g_zzzzz_yyyyy, g_zzzzz_yyyyz, g_zzzzz_yyyzz, g_zzzzz_yyzzz, g_zzzzz_yzzzz, g_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxxxx[k] = -g_zzzzz_xxxxx[k] - g_z_0_zzzzz_xxxxx[k] * ab_z + g_z_0_zzzzz_xxxxxz[k];

                g_z_0_zzzzzz_xxxxy[k] = -g_zzzzz_xxxxy[k] - g_z_0_zzzzz_xxxxy[k] * ab_z + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_zzzzzz_xxxxz[k] = -g_zzzzz_xxxxz[k] - g_z_0_zzzzz_xxxxz[k] * ab_z + g_z_0_zzzzz_xxxxzz[k];

                g_z_0_zzzzzz_xxxyy[k] = -g_zzzzz_xxxyy[k] - g_z_0_zzzzz_xxxyy[k] * ab_z + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_zzzzzz_xxxyz[k] = -g_zzzzz_xxxyz[k] - g_z_0_zzzzz_xxxyz[k] * ab_z + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_zzzzzz_xxxzz[k] = -g_zzzzz_xxxzz[k] - g_z_0_zzzzz_xxxzz[k] * ab_z + g_z_0_zzzzz_xxxzzz[k];

                g_z_0_zzzzzz_xxyyy[k] = -g_zzzzz_xxyyy[k] - g_z_0_zzzzz_xxyyy[k] * ab_z + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_zzzzzz_xxyyz[k] = -g_zzzzz_xxyyz[k] - g_z_0_zzzzz_xxyyz[k] * ab_z + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_zzzzzz_xxyzz[k] = -g_zzzzz_xxyzz[k] - g_z_0_zzzzz_xxyzz[k] * ab_z + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_zzzzzz_xxzzz[k] = -g_zzzzz_xxzzz[k] - g_z_0_zzzzz_xxzzz[k] * ab_z + g_z_0_zzzzz_xxzzzz[k];

                g_z_0_zzzzzz_xyyyy[k] = -g_zzzzz_xyyyy[k] - g_z_0_zzzzz_xyyyy[k] * ab_z + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_zzzzzz_xyyyz[k] = -g_zzzzz_xyyyz[k] - g_z_0_zzzzz_xyyyz[k] * ab_z + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_zzzzzz_xyyzz[k] = -g_zzzzz_xyyzz[k] - g_z_0_zzzzz_xyyzz[k] * ab_z + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_zzzzzz_xyzzz[k] = -g_zzzzz_xyzzz[k] - g_z_0_zzzzz_xyzzz[k] * ab_z + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_zzzzzz_xzzzz[k] = -g_zzzzz_xzzzz[k] - g_z_0_zzzzz_xzzzz[k] * ab_z + g_z_0_zzzzz_xzzzzz[k];

                g_z_0_zzzzzz_yyyyy[k] = -g_zzzzz_yyyyy[k] - g_z_0_zzzzz_yyyyy[k] * ab_z + g_z_0_zzzzz_yyyyyz[k];

                g_z_0_zzzzzz_yyyyz[k] = -g_zzzzz_yyyyz[k] - g_z_0_zzzzz_yyyyz[k] * ab_z + g_z_0_zzzzz_yyyyzz[k];

                g_z_0_zzzzzz_yyyzz[k] = -g_zzzzz_yyyzz[k] - g_z_0_zzzzz_yyyzz[k] * ab_z + g_z_0_zzzzz_yyyzzz[k];

                g_z_0_zzzzzz_yyzzz[k] = -g_zzzzz_yyzzz[k] - g_z_0_zzzzz_yyzzz[k] * ab_z + g_z_0_zzzzz_yyzzzz[k];

                g_z_0_zzzzzz_yzzzz[k] = -g_zzzzz_yzzzz[k] - g_z_0_zzzzz_yzzzz[k] * ab_z + g_z_0_zzzzz_yzzzzz[k];

                g_z_0_zzzzzz_zzzzz[k] = -g_zzzzz_zzzzz[k] - g_z_0_zzzzz_zzzzz[k] * ab_z + g_z_0_zzzzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace


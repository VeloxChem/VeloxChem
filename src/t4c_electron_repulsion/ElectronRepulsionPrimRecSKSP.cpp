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

#include "ElectronRepulsionPrimRecSKSP.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sksp(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sksp,
                                  size_t                idx_eri_0_shsp,
                                  size_t                idx_eri_1_shsp,
                                  size_t                idx_eri_1_siss,
                                  size_t                idx_eri_0_sisp,
                                  size_t                idx_eri_1_sisp,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SHSP

    auto g_0_xxxxx_0_x_0 = pbuffer.data(idx_eri_0_shsp);

    auto g_0_xxxxx_0_y_0 = pbuffer.data(idx_eri_0_shsp + 1);

    auto g_0_xxxxx_0_z_0 = pbuffer.data(idx_eri_0_shsp + 2);

    auto g_0_xxxxy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 3);

    auto g_0_xxxxz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 6);

    auto g_0_xxxyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 9);

    auto g_0_xxxyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 10);

    auto g_0_xxxyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 11);

    auto g_0_xxxzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 15);

    auto g_0_xxxzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 16);

    auto g_0_xxxzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 17);

    auto g_0_xxyyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 18);

    auto g_0_xxyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 19);

    auto g_0_xxyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 20);

    auto g_0_xxyzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 24);

    auto g_0_xxzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 27);

    auto g_0_xxzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 28);

    auto g_0_xxzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 29);

    auto g_0_xyyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 31);

    auto g_0_xyyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 32);

    auto g_0_xyyzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 37);

    auto g_0_xyyzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 38);

    auto g_0_xzzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 43);

    auto g_0_xzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 44);

    auto g_0_yyyyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 45);

    auto g_0_yyyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 46);

    auto g_0_yyyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 47);

    auto g_0_yyyyz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 49);

    auto g_0_yyyzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 51);

    auto g_0_yyyzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 52);

    auto g_0_yyyzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 53);

    auto g_0_yyzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 54);

    auto g_0_yyzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 55);

    auto g_0_yyzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 56);

    auto g_0_yzzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 57);

    auto g_0_yzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 59);

    auto g_0_zzzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 60);

    auto g_0_zzzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 61);

    auto g_0_zzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 62);

    /// Set up components of auxilary buffer : SHSP

    auto g_0_xxxxx_0_x_1 = pbuffer.data(idx_eri_1_shsp);

    auto g_0_xxxxx_0_y_1 = pbuffer.data(idx_eri_1_shsp + 1);

    auto g_0_xxxxx_0_z_1 = pbuffer.data(idx_eri_1_shsp + 2);

    auto g_0_xxxxy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 3);

    auto g_0_xxxxz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 6);

    auto g_0_xxxyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 9);

    auto g_0_xxxyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 10);

    auto g_0_xxxyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 11);

    auto g_0_xxxzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 15);

    auto g_0_xxxzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 16);

    auto g_0_xxxzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 17);

    auto g_0_xxyyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 18);

    auto g_0_xxyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 19);

    auto g_0_xxyyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 20);

    auto g_0_xxyzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 24);

    auto g_0_xxzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 27);

    auto g_0_xxzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 28);

    auto g_0_xxzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 29);

    auto g_0_xyyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 31);

    auto g_0_xyyyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 32);

    auto g_0_xyyzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 37);

    auto g_0_xyyzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 38);

    auto g_0_xzzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 43);

    auto g_0_xzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 44);

    auto g_0_yyyyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 45);

    auto g_0_yyyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 46);

    auto g_0_yyyyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 47);

    auto g_0_yyyyz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 49);

    auto g_0_yyyzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 51);

    auto g_0_yyyzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 52);

    auto g_0_yyyzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 53);

    auto g_0_yyzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 54);

    auto g_0_yyzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 55);

    auto g_0_yyzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 56);

    auto g_0_yzzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 57);

    auto g_0_yzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 59);

    auto g_0_zzzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 60);

    auto g_0_zzzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 61);

    auto g_0_zzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 62);

    /// Set up components of auxilary buffer : SISS

    auto g_0_xxxxxx_0_0_1 = pbuffer.data(idx_eri_1_siss);

    auto g_0_xxxxyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 3);

    auto g_0_xxxxzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 5);

    auto g_0_xxxyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 6);

    auto g_0_xxxzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 9);

    auto g_0_xxyyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 10);

    auto g_0_xxzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 14);

    auto g_0_yyyyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 21);

    auto g_0_yyyyzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 23);

    auto g_0_yyyzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 24);

    auto g_0_yyzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 25);

    auto g_0_zzzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 27);

    /// Set up components of auxilary buffer : SISP

    auto g_0_xxxxxx_0_x_0 = pbuffer.data(idx_eri_0_sisp);

    auto g_0_xxxxxx_0_y_0 = pbuffer.data(idx_eri_0_sisp + 1);

    auto g_0_xxxxxx_0_z_0 = pbuffer.data(idx_eri_0_sisp + 2);

    auto g_0_xxxxxy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 3);

    auto g_0_xxxxxy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 4);

    auto g_0_xxxxxz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 6);

    auto g_0_xxxxxz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 8);

    auto g_0_xxxxyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 9);

    auto g_0_xxxxyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 10);

    auto g_0_xxxxyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 11);

    auto g_0_xxxxzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 15);

    auto g_0_xxxxzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 16);

    auto g_0_xxxxzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 17);

    auto g_0_xxxyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 18);

    auto g_0_xxxyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 19);

    auto g_0_xxxyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 20);

    auto g_0_xxxyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 24);

    auto g_0_xxxzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 27);

    auto g_0_xxxzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 28);

    auto g_0_xxxzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 29);

    auto g_0_xxyyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 30);

    auto g_0_xxyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 31);

    auto g_0_xxyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 32);

    auto g_0_xxyyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 36);

    auto g_0_xxyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 37);

    auto g_0_xxyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 38);

    auto g_0_xxyzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 39);

    auto g_0_xxzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 42);

    auto g_0_xxzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 43);

    auto g_0_xxzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 44);

    auto g_0_xyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 45);

    auto g_0_xyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 46);

    auto g_0_xyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 47);

    auto g_0_xyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 52);

    auto g_0_xyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 53);

    auto g_0_xyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 55);

    auto g_0_xyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 56);

    auto g_0_xzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 60);

    auto g_0_xzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 61);

    auto g_0_xzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 62);

    auto g_0_yyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 63);

    auto g_0_yyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 64);

    auto g_0_yyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 65);

    auto g_0_yyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 67);

    auto g_0_yyyyyz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 68);

    auto g_0_yyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 69);

    auto g_0_yyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 70);

    auto g_0_yyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 71);

    auto g_0_yyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 72);

    auto g_0_yyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 73);

    auto g_0_yyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 74);

    auto g_0_yyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 75);

    auto g_0_yyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 76);

    auto g_0_yyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 77);

    auto g_0_yzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 78);

    auto g_0_yzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 79);

    auto g_0_yzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 80);

    auto g_0_zzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 81);

    auto g_0_zzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 82);

    auto g_0_zzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 83);

    /// Set up components of auxilary buffer : SISP

    auto g_0_xxxxxx_0_x_1 = pbuffer.data(idx_eri_1_sisp);

    auto g_0_xxxxxx_0_y_1 = pbuffer.data(idx_eri_1_sisp + 1);

    auto g_0_xxxxxx_0_z_1 = pbuffer.data(idx_eri_1_sisp + 2);

    auto g_0_xxxxxy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 3);

    auto g_0_xxxxxy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 4);

    auto g_0_xxxxxz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 6);

    auto g_0_xxxxxz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 8);

    auto g_0_xxxxyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 9);

    auto g_0_xxxxyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 10);

    auto g_0_xxxxyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 11);

    auto g_0_xxxxzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 15);

    auto g_0_xxxxzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 16);

    auto g_0_xxxxzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 17);

    auto g_0_xxxyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 18);

    auto g_0_xxxyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 19);

    auto g_0_xxxyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 20);

    auto g_0_xxxyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 24);

    auto g_0_xxxzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 27);

    auto g_0_xxxzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 28);

    auto g_0_xxxzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 29);

    auto g_0_xxyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 30);

    auto g_0_xxyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 31);

    auto g_0_xxyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 32);

    auto g_0_xxyyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 36);

    auto g_0_xxyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 37);

    auto g_0_xxyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 38);

    auto g_0_xxyzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 39);

    auto g_0_xxzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 42);

    auto g_0_xxzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 43);

    auto g_0_xxzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 44);

    auto g_0_xyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 45);

    auto g_0_xyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 46);

    auto g_0_xyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 47);

    auto g_0_xyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 52);

    auto g_0_xyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 53);

    auto g_0_xyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 55);

    auto g_0_xyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 56);

    auto g_0_xzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 60);

    auto g_0_xzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 61);

    auto g_0_xzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 62);

    auto g_0_yyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 63);

    auto g_0_yyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 64);

    auto g_0_yyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 65);

    auto g_0_yyyyyz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 67);

    auto g_0_yyyyyz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 68);

    auto g_0_yyyyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 69);

    auto g_0_yyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 70);

    auto g_0_yyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 71);

    auto g_0_yyyzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 72);

    auto g_0_yyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 73);

    auto g_0_yyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 74);

    auto g_0_yyzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 75);

    auto g_0_yyzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 76);

    auto g_0_yyzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 77);

    auto g_0_yzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 78);

    auto g_0_yzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 79);

    auto g_0_yzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 80);

    auto g_0_zzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 81);

    auto g_0_zzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 82);

    auto g_0_zzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 83);

    /// Set up 0-3 components of targeted buffer : SKSP

    auto g_0_xxxxxxx_0_x_0 = pbuffer.data(idx_eri_0_sksp);

    auto g_0_xxxxxxx_0_y_0 = pbuffer.data(idx_eri_0_sksp + 1);

    auto g_0_xxxxxxx_0_z_0 = pbuffer.data(idx_eri_0_sksp + 2);

#pragma omp simd aligned(g_0_xxxxx_0_x_0,       \
                             g_0_xxxxx_0_x_1,   \
                             g_0_xxxxx_0_y_0,   \
                             g_0_xxxxx_0_y_1,   \
                             g_0_xxxxx_0_z_0,   \
                             g_0_xxxxx_0_z_1,   \
                             g_0_xxxxxx_0_0_1,  \
                             g_0_xxxxxx_0_x_0,  \
                             g_0_xxxxxx_0_x_1,  \
                             g_0_xxxxxx_0_y_0,  \
                             g_0_xxxxxx_0_y_1,  \
                             g_0_xxxxxx_0_z_0,  \
                             g_0_xxxxxx_0_z_1,  \
                             g_0_xxxxxxx_0_x_0, \
                             g_0_xxxxxxx_0_y_0, \
                             g_0_xxxxxxx_0_z_0, \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_x_0[i] = 6.0 * g_0_xxxxx_0_x_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxx_0_0_1[i] * fi_abcd_0 +
                               g_0_xxxxxx_0_x_0[i] * pb_x + g_0_xxxxxx_0_x_1[i] * wp_x[i];

        g_0_xxxxxxx_0_y_0[i] =
            6.0 * g_0_xxxxx_0_y_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_y_1[i] * fti_ab_0 + g_0_xxxxxx_0_y_0[i] * pb_x + g_0_xxxxxx_0_y_1[i] * wp_x[i];

        g_0_xxxxxxx_0_z_0[i] =
            6.0 * g_0_xxxxx_0_z_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_z_1[i] * fti_ab_0 + g_0_xxxxxx_0_z_0[i] * pb_x + g_0_xxxxxx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : SKSP

    auto g_0_xxxxxxy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 3);

    auto g_0_xxxxxxy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 4);

    auto g_0_xxxxxxy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 5);

#pragma omp simd aligned(g_0_xxxxxx_0_0_1,      \
                             g_0_xxxxxx_0_x_0,  \
                             g_0_xxxxxx_0_x_1,  \
                             g_0_xxxxxx_0_y_0,  \
                             g_0_xxxxxx_0_y_1,  \
                             g_0_xxxxxx_0_z_0,  \
                             g_0_xxxxxx_0_z_1,  \
                             g_0_xxxxxxy_0_x_0, \
                             g_0_xxxxxxy_0_y_0, \
                             g_0_xxxxxxy_0_z_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_x_0[i] = g_0_xxxxxx_0_x_0[i] * pb_y + g_0_xxxxxx_0_x_1[i] * wp_y[i];

        g_0_xxxxxxy_0_y_0[i] = g_0_xxxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxxx_0_y_0[i] * pb_y + g_0_xxxxxx_0_y_1[i] * wp_y[i];

        g_0_xxxxxxy_0_z_0[i] = g_0_xxxxxx_0_z_0[i] * pb_y + g_0_xxxxxx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : SKSP

    auto g_0_xxxxxxz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 6);

    auto g_0_xxxxxxz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 7);

    auto g_0_xxxxxxz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 8);

#pragma omp simd aligned(g_0_xxxxxx_0_0_1,      \
                             g_0_xxxxxx_0_x_0,  \
                             g_0_xxxxxx_0_x_1,  \
                             g_0_xxxxxx_0_y_0,  \
                             g_0_xxxxxx_0_y_1,  \
                             g_0_xxxxxx_0_z_0,  \
                             g_0_xxxxxx_0_z_1,  \
                             g_0_xxxxxxz_0_x_0, \
                             g_0_xxxxxxz_0_y_0, \
                             g_0_xxxxxxz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_x_0[i] = g_0_xxxxxx_0_x_0[i] * pb_z + g_0_xxxxxx_0_x_1[i] * wp_z[i];

        g_0_xxxxxxz_0_y_0[i] = g_0_xxxxxx_0_y_0[i] * pb_z + g_0_xxxxxx_0_y_1[i] * wp_z[i];

        g_0_xxxxxxz_0_z_0[i] = g_0_xxxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxxx_0_z_0[i] * pb_z + g_0_xxxxxx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : SKSP

    auto g_0_xxxxxyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 9);

    auto g_0_xxxxxyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 10);

    auto g_0_xxxxxyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 11);

#pragma omp simd aligned(g_0_xxxxx_0_x_0,       \
                             g_0_xxxxx_0_x_1,   \
                             g_0_xxxxxy_0_x_0,  \
                             g_0_xxxxxy_0_x_1,  \
                             g_0_xxxxxyy_0_x_0, \
                             g_0_xxxxxyy_0_y_0, \
                             g_0_xxxxxyy_0_z_0, \
                             g_0_xxxxyy_0_y_0,  \
                             g_0_xxxxyy_0_y_1,  \
                             g_0_xxxxyy_0_z_0,  \
                             g_0_xxxxyy_0_z_1,  \
                             g_0_xxxyy_0_y_0,   \
                             g_0_xxxyy_0_y_1,   \
                             g_0_xxxyy_0_z_0,   \
                             g_0_xxxyy_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_x_0[i] =
            g_0_xxxxx_0_x_0[i] * fi_ab_0 - g_0_xxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxy_0_x_0[i] * pb_y + g_0_xxxxxy_0_x_1[i] * wp_y[i];

        g_0_xxxxxyy_0_y_0[i] =
            4.0 * g_0_xxxyy_0_y_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_y_1[i] * fti_ab_0 + g_0_xxxxyy_0_y_0[i] * pb_x + g_0_xxxxyy_0_y_1[i] * wp_x[i];

        g_0_xxxxxyy_0_z_0[i] =
            4.0 * g_0_xxxyy_0_z_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_z_1[i] * fti_ab_0 + g_0_xxxxyy_0_z_0[i] * pb_x + g_0_xxxxyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : SKSP

    auto g_0_xxxxxyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 12);

    auto g_0_xxxxxyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 13);

    auto g_0_xxxxxyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 14);

#pragma omp simd aligned(g_0_xxxxxy_0_y_0,      \
                             g_0_xxxxxy_0_y_1,  \
                             g_0_xxxxxyz_0_x_0, \
                             g_0_xxxxxyz_0_y_0, \
                             g_0_xxxxxyz_0_z_0, \
                             g_0_xxxxxz_0_x_0,  \
                             g_0_xxxxxz_0_x_1,  \
                             g_0_xxxxxz_0_z_0,  \
                             g_0_xxxxxz_0_z_1,  \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xxxxxyz_0_x_0[i] = g_0_xxxxxz_0_x_0[i] * pb_y + g_0_xxxxxz_0_x_1[i] * wp_y[i];

        g_0_xxxxxyz_0_y_0[i] = g_0_xxxxxy_0_y_0[i] * pb_z + g_0_xxxxxy_0_y_1[i] * wp_z[i];

        g_0_xxxxxyz_0_z_0[i] = g_0_xxxxxz_0_z_0[i] * pb_y + g_0_xxxxxz_0_z_1[i] * wp_y[i];
    }

    /// Set up 15-18 components of targeted buffer : SKSP

    auto g_0_xxxxxzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 15);

    auto g_0_xxxxxzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 16);

    auto g_0_xxxxxzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 17);

#pragma omp simd aligned(g_0_xxxxx_0_x_0,       \
                             g_0_xxxxx_0_x_1,   \
                             g_0_xxxxxz_0_x_0,  \
                             g_0_xxxxxz_0_x_1,  \
                             g_0_xxxxxzz_0_x_0, \
                             g_0_xxxxxzz_0_y_0, \
                             g_0_xxxxxzz_0_z_0, \
                             g_0_xxxxzz_0_y_0,  \
                             g_0_xxxxzz_0_y_1,  \
                             g_0_xxxxzz_0_z_0,  \
                             g_0_xxxxzz_0_z_1,  \
                             g_0_xxxzz_0_y_0,   \
                             g_0_xxxzz_0_y_1,   \
                             g_0_xxxzz_0_z_0,   \
                             g_0_xxxzz_0_z_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_x_0[i] =
            g_0_xxxxx_0_x_0[i] * fi_ab_0 - g_0_xxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxz_0_x_0[i] * pb_z + g_0_xxxxxz_0_x_1[i] * wp_z[i];

        g_0_xxxxxzz_0_y_0[i] =
            4.0 * g_0_xxxzz_0_y_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_y_1[i] * fti_ab_0 + g_0_xxxxzz_0_y_0[i] * pb_x + g_0_xxxxzz_0_y_1[i] * wp_x[i];

        g_0_xxxxxzz_0_z_0[i] =
            4.0 * g_0_xxxzz_0_z_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_z_1[i] * fti_ab_0 + g_0_xxxxzz_0_z_0[i] * pb_x + g_0_xxxxzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : SKSP

    auto g_0_xxxxyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 18);

    auto g_0_xxxxyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 19);

    auto g_0_xxxxyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 20);

#pragma omp simd aligned(g_0_xxxxy_0_x_0,       \
                             g_0_xxxxy_0_x_1,   \
                             g_0_xxxxyy_0_x_0,  \
                             g_0_xxxxyy_0_x_1,  \
                             g_0_xxxxyyy_0_x_0, \
                             g_0_xxxxyyy_0_y_0, \
                             g_0_xxxxyyy_0_z_0, \
                             g_0_xxxyyy_0_y_0,  \
                             g_0_xxxyyy_0_y_1,  \
                             g_0_xxxyyy_0_z_0,  \
                             g_0_xxxyyy_0_z_1,  \
                             g_0_xxyyy_0_y_0,   \
                             g_0_xxyyy_0_y_1,   \
                             g_0_xxyyy_0_z_0,   \
                             g_0_xxyyy_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_x_0[i] =
            2.0 * g_0_xxxxy_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_x_1[i] * fti_ab_0 + g_0_xxxxyy_0_x_0[i] * pb_y + g_0_xxxxyy_0_x_1[i] * wp_y[i];

        g_0_xxxxyyy_0_y_0[i] =
            3.0 * g_0_xxyyy_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_y_1[i] * fti_ab_0 + g_0_xxxyyy_0_y_0[i] * pb_x + g_0_xxxyyy_0_y_1[i] * wp_x[i];

        g_0_xxxxyyy_0_z_0[i] =
            3.0 * g_0_xxyyy_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_z_1[i] * fti_ab_0 + g_0_xxxyyy_0_z_0[i] * pb_x + g_0_xxxyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 21-24 components of targeted buffer : SKSP

    auto g_0_xxxxyyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 21);

    auto g_0_xxxxyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 22);

    auto g_0_xxxxyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 23);

#pragma omp simd aligned(g_0_xxxxyy_0_0_1,      \
                             g_0_xxxxyy_0_x_0,  \
                             g_0_xxxxyy_0_x_1,  \
                             g_0_xxxxyy_0_y_0,  \
                             g_0_xxxxyy_0_y_1,  \
                             g_0_xxxxyy_0_z_0,  \
                             g_0_xxxxyy_0_z_1,  \
                             g_0_xxxxyyz_0_x_0, \
                             g_0_xxxxyyz_0_y_0, \
                             g_0_xxxxyyz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_x_0[i] = g_0_xxxxyy_0_x_0[i] * pb_z + g_0_xxxxyy_0_x_1[i] * wp_z[i];

        g_0_xxxxyyz_0_y_0[i] = g_0_xxxxyy_0_y_0[i] * pb_z + g_0_xxxxyy_0_y_1[i] * wp_z[i];

        g_0_xxxxyyz_0_z_0[i] = g_0_xxxxyy_0_0_1[i] * fi_abcd_0 + g_0_xxxxyy_0_z_0[i] * pb_z + g_0_xxxxyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 24-27 components of targeted buffer : SKSP

    auto g_0_xxxxyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 24);

    auto g_0_xxxxyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 25);

    auto g_0_xxxxyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 26);

#pragma omp simd aligned(g_0_xxxxyzz_0_x_0,     \
                             g_0_xxxxyzz_0_y_0, \
                             g_0_xxxxyzz_0_z_0, \
                             g_0_xxxxzz_0_0_1,  \
                             g_0_xxxxzz_0_x_0,  \
                             g_0_xxxxzz_0_x_1,  \
                             g_0_xxxxzz_0_y_0,  \
                             g_0_xxxxzz_0_y_1,  \
                             g_0_xxxxzz_0_z_0,  \
                             g_0_xxxxzz_0_z_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_x_0[i] = g_0_xxxxzz_0_x_0[i] * pb_y + g_0_xxxxzz_0_x_1[i] * wp_y[i];

        g_0_xxxxyzz_0_y_0[i] = g_0_xxxxzz_0_0_1[i] * fi_abcd_0 + g_0_xxxxzz_0_y_0[i] * pb_y + g_0_xxxxzz_0_y_1[i] * wp_y[i];

        g_0_xxxxyzz_0_z_0[i] = g_0_xxxxzz_0_z_0[i] * pb_y + g_0_xxxxzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 27-30 components of targeted buffer : SKSP

    auto g_0_xxxxzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 27);

    auto g_0_xxxxzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 28);

    auto g_0_xxxxzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 29);

#pragma omp simd aligned(g_0_xxxxz_0_x_0,       \
                             g_0_xxxxz_0_x_1,   \
                             g_0_xxxxzz_0_x_0,  \
                             g_0_xxxxzz_0_x_1,  \
                             g_0_xxxxzzz_0_x_0, \
                             g_0_xxxxzzz_0_y_0, \
                             g_0_xxxxzzz_0_z_0, \
                             g_0_xxxzzz_0_y_0,  \
                             g_0_xxxzzz_0_y_1,  \
                             g_0_xxxzzz_0_z_0,  \
                             g_0_xxxzzz_0_z_1,  \
                             g_0_xxzzz_0_y_0,   \
                             g_0_xxzzz_0_y_1,   \
                             g_0_xxzzz_0_z_0,   \
                             g_0_xxzzz_0_z_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_x_0[i] =
            2.0 * g_0_xxxxz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_x_1[i] * fti_ab_0 + g_0_xxxxzz_0_x_0[i] * pb_z + g_0_xxxxzz_0_x_1[i] * wp_z[i];

        g_0_xxxxzzz_0_y_0[i] =
            3.0 * g_0_xxzzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_y_1[i] * fti_ab_0 + g_0_xxxzzz_0_y_0[i] * pb_x + g_0_xxxzzz_0_y_1[i] * wp_x[i];

        g_0_xxxxzzz_0_z_0[i] =
            3.0 * g_0_xxzzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_z_1[i] * fti_ab_0 + g_0_xxxzzz_0_z_0[i] * pb_x + g_0_xxxzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 30-33 components of targeted buffer : SKSP

    auto g_0_xxxyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 30);

    auto g_0_xxxyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 31);

    auto g_0_xxxyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 32);

#pragma omp simd aligned(g_0_xxxyy_0_x_0,       \
                             g_0_xxxyy_0_x_1,   \
                             g_0_xxxyyy_0_x_0,  \
                             g_0_xxxyyy_0_x_1,  \
                             g_0_xxxyyyy_0_x_0, \
                             g_0_xxxyyyy_0_y_0, \
                             g_0_xxxyyyy_0_z_0, \
                             g_0_xxyyyy_0_y_0,  \
                             g_0_xxyyyy_0_y_1,  \
                             g_0_xxyyyy_0_z_0,  \
                             g_0_xxyyyy_0_z_1,  \
                             g_0_xyyyy_0_y_0,   \
                             g_0_xyyyy_0_y_1,   \
                             g_0_xyyyy_0_z_0,   \
                             g_0_xyyyy_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_x_0[i] =
            3.0 * g_0_xxxyy_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_x_1[i] * fti_ab_0 + g_0_xxxyyy_0_x_0[i] * pb_y + g_0_xxxyyy_0_x_1[i] * wp_y[i];

        g_0_xxxyyyy_0_y_0[i] =
            2.0 * g_0_xyyyy_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_y_1[i] * fti_ab_0 + g_0_xxyyyy_0_y_0[i] * pb_x + g_0_xxyyyy_0_y_1[i] * wp_x[i];

        g_0_xxxyyyy_0_z_0[i] =
            2.0 * g_0_xyyyy_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_z_1[i] * fti_ab_0 + g_0_xxyyyy_0_z_0[i] * pb_x + g_0_xxyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 33-36 components of targeted buffer : SKSP

    auto g_0_xxxyyyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 33);

    auto g_0_xxxyyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 34);

    auto g_0_xxxyyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 35);

#pragma omp simd aligned(g_0_xxxyyy_0_0_1,      \
                             g_0_xxxyyy_0_x_0,  \
                             g_0_xxxyyy_0_x_1,  \
                             g_0_xxxyyy_0_y_0,  \
                             g_0_xxxyyy_0_y_1,  \
                             g_0_xxxyyy_0_z_0,  \
                             g_0_xxxyyy_0_z_1,  \
                             g_0_xxxyyyz_0_x_0, \
                             g_0_xxxyyyz_0_y_0, \
                             g_0_xxxyyyz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_x_0[i] = g_0_xxxyyy_0_x_0[i] * pb_z + g_0_xxxyyy_0_x_1[i] * wp_z[i];

        g_0_xxxyyyz_0_y_0[i] = g_0_xxxyyy_0_y_0[i] * pb_z + g_0_xxxyyy_0_y_1[i] * wp_z[i];

        g_0_xxxyyyz_0_z_0[i] = g_0_xxxyyy_0_0_1[i] * fi_abcd_0 + g_0_xxxyyy_0_z_0[i] * pb_z + g_0_xxxyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 36-39 components of targeted buffer : SKSP

    auto g_0_xxxyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 36);

    auto g_0_xxxyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 37);

    auto g_0_xxxyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 38);

#pragma omp simd aligned(g_0_xxxyyzz_0_x_0,     \
                             g_0_xxxyyzz_0_y_0, \
                             g_0_xxxyyzz_0_z_0, \
                             g_0_xxxyzz_0_x_0,  \
                             g_0_xxxyzz_0_x_1,  \
                             g_0_xxxzz_0_x_0,   \
                             g_0_xxxzz_0_x_1,   \
                             g_0_xxyyzz_0_y_0,  \
                             g_0_xxyyzz_0_y_1,  \
                             g_0_xxyyzz_0_z_0,  \
                             g_0_xxyyzz_0_z_1,  \
                             g_0_xyyzz_0_y_0,   \
                             g_0_xyyzz_0_y_1,   \
                             g_0_xyyzz_0_z_0,   \
                             g_0_xyyzz_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_x_0[i] =
            g_0_xxxzz_0_x_0[i] * fi_ab_0 - g_0_xxxzz_0_x_1[i] * fti_ab_0 + g_0_xxxyzz_0_x_0[i] * pb_y + g_0_xxxyzz_0_x_1[i] * wp_y[i];

        g_0_xxxyyzz_0_y_0[i] =
            2.0 * g_0_xyyzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_y_1[i] * fti_ab_0 + g_0_xxyyzz_0_y_0[i] * pb_x + g_0_xxyyzz_0_y_1[i] * wp_x[i];

        g_0_xxxyyzz_0_z_0[i] =
            2.0 * g_0_xyyzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_z_1[i] * fti_ab_0 + g_0_xxyyzz_0_z_0[i] * pb_x + g_0_xxyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 39-42 components of targeted buffer : SKSP

    auto g_0_xxxyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 39);

    auto g_0_xxxyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 40);

    auto g_0_xxxyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 41);

#pragma omp simd aligned(g_0_xxxyzzz_0_x_0,     \
                             g_0_xxxyzzz_0_y_0, \
                             g_0_xxxyzzz_0_z_0, \
                             g_0_xxxzzz_0_0_1,  \
                             g_0_xxxzzz_0_x_0,  \
                             g_0_xxxzzz_0_x_1,  \
                             g_0_xxxzzz_0_y_0,  \
                             g_0_xxxzzz_0_y_1,  \
                             g_0_xxxzzz_0_z_0,  \
                             g_0_xxxzzz_0_z_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_x_0[i] = g_0_xxxzzz_0_x_0[i] * pb_y + g_0_xxxzzz_0_x_1[i] * wp_y[i];

        g_0_xxxyzzz_0_y_0[i] = g_0_xxxzzz_0_0_1[i] * fi_abcd_0 + g_0_xxxzzz_0_y_0[i] * pb_y + g_0_xxxzzz_0_y_1[i] * wp_y[i];

        g_0_xxxyzzz_0_z_0[i] = g_0_xxxzzz_0_z_0[i] * pb_y + g_0_xxxzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 42-45 components of targeted buffer : SKSP

    auto g_0_xxxzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 42);

    auto g_0_xxxzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 43);

    auto g_0_xxxzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 44);

#pragma omp simd aligned(g_0_xxxzz_0_x_0,       \
                             g_0_xxxzz_0_x_1,   \
                             g_0_xxxzzz_0_x_0,  \
                             g_0_xxxzzz_0_x_1,  \
                             g_0_xxxzzzz_0_x_0, \
                             g_0_xxxzzzz_0_y_0, \
                             g_0_xxxzzzz_0_z_0, \
                             g_0_xxzzzz_0_y_0,  \
                             g_0_xxzzzz_0_y_1,  \
                             g_0_xxzzzz_0_z_0,  \
                             g_0_xxzzzz_0_z_1,  \
                             g_0_xzzzz_0_y_0,   \
                             g_0_xzzzz_0_y_1,   \
                             g_0_xzzzz_0_z_0,   \
                             g_0_xzzzz_0_z_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_x_0[i] =
            3.0 * g_0_xxxzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_x_1[i] * fti_ab_0 + g_0_xxxzzz_0_x_0[i] * pb_z + g_0_xxxzzz_0_x_1[i] * wp_z[i];

        g_0_xxxzzzz_0_y_0[i] =
            2.0 * g_0_xzzzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_y_1[i] * fti_ab_0 + g_0_xxzzzz_0_y_0[i] * pb_x + g_0_xxzzzz_0_y_1[i] * wp_x[i];

        g_0_xxxzzzz_0_z_0[i] =
            2.0 * g_0_xzzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_z_1[i] * fti_ab_0 + g_0_xxzzzz_0_z_0[i] * pb_x + g_0_xxzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 45-48 components of targeted buffer : SKSP

    auto g_0_xxyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 45);

    auto g_0_xxyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 46);

    auto g_0_xxyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 47);

#pragma omp simd aligned(g_0_xxyyy_0_x_0,       \
                             g_0_xxyyy_0_x_1,   \
                             g_0_xxyyyy_0_x_0,  \
                             g_0_xxyyyy_0_x_1,  \
                             g_0_xxyyyyy_0_x_0, \
                             g_0_xxyyyyy_0_y_0, \
                             g_0_xxyyyyy_0_z_0, \
                             g_0_xyyyyy_0_y_0,  \
                             g_0_xyyyyy_0_y_1,  \
                             g_0_xyyyyy_0_z_0,  \
                             g_0_xyyyyy_0_z_1,  \
                             g_0_yyyyy_0_y_0,   \
                             g_0_yyyyy_0_y_1,   \
                             g_0_yyyyy_0_z_0,   \
                             g_0_yyyyy_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_x_0[i] =
            4.0 * g_0_xxyyy_0_x_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_x_1[i] * fti_ab_0 + g_0_xxyyyy_0_x_0[i] * pb_y + g_0_xxyyyy_0_x_1[i] * wp_y[i];

        g_0_xxyyyyy_0_y_0[i] =
            g_0_yyyyy_0_y_0[i] * fi_ab_0 - g_0_yyyyy_0_y_1[i] * fti_ab_0 + g_0_xyyyyy_0_y_0[i] * pb_x + g_0_xyyyyy_0_y_1[i] * wp_x[i];

        g_0_xxyyyyy_0_z_0[i] =
            g_0_yyyyy_0_z_0[i] * fi_ab_0 - g_0_yyyyy_0_z_1[i] * fti_ab_0 + g_0_xyyyyy_0_z_0[i] * pb_x + g_0_xyyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 48-51 components of targeted buffer : SKSP

    auto g_0_xxyyyyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 48);

    auto g_0_xxyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 49);

    auto g_0_xxyyyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 50);

#pragma omp simd aligned(g_0_xxyyyy_0_0_1,      \
                             g_0_xxyyyy_0_x_0,  \
                             g_0_xxyyyy_0_x_1,  \
                             g_0_xxyyyy_0_y_0,  \
                             g_0_xxyyyy_0_y_1,  \
                             g_0_xxyyyy_0_z_0,  \
                             g_0_xxyyyy_0_z_1,  \
                             g_0_xxyyyyz_0_x_0, \
                             g_0_xxyyyyz_0_y_0, \
                             g_0_xxyyyyz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_x_0[i] = g_0_xxyyyy_0_x_0[i] * pb_z + g_0_xxyyyy_0_x_1[i] * wp_z[i];

        g_0_xxyyyyz_0_y_0[i] = g_0_xxyyyy_0_y_0[i] * pb_z + g_0_xxyyyy_0_y_1[i] * wp_z[i];

        g_0_xxyyyyz_0_z_0[i] = g_0_xxyyyy_0_0_1[i] * fi_abcd_0 + g_0_xxyyyy_0_z_0[i] * pb_z + g_0_xxyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 51-54 components of targeted buffer : SKSP

    auto g_0_xxyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 51);

    auto g_0_xxyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 52);

    auto g_0_xxyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 53);

#pragma omp simd aligned(g_0_xxyyyzz_0_x_0,     \
                             g_0_xxyyyzz_0_y_0, \
                             g_0_xxyyyzz_0_z_0, \
                             g_0_xxyyzz_0_x_0,  \
                             g_0_xxyyzz_0_x_1,  \
                             g_0_xxyzz_0_x_0,   \
                             g_0_xxyzz_0_x_1,   \
                             g_0_xyyyzz_0_y_0,  \
                             g_0_xyyyzz_0_y_1,  \
                             g_0_xyyyzz_0_z_0,  \
                             g_0_xyyyzz_0_z_1,  \
                             g_0_yyyzz_0_y_0,   \
                             g_0_yyyzz_0_y_1,   \
                             g_0_yyyzz_0_z_0,   \
                             g_0_yyyzz_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_x_0[i] =
            2.0 * g_0_xxyzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_x_1[i] * fti_ab_0 + g_0_xxyyzz_0_x_0[i] * pb_y + g_0_xxyyzz_0_x_1[i] * wp_y[i];

        g_0_xxyyyzz_0_y_0[i] =
            g_0_yyyzz_0_y_0[i] * fi_ab_0 - g_0_yyyzz_0_y_1[i] * fti_ab_0 + g_0_xyyyzz_0_y_0[i] * pb_x + g_0_xyyyzz_0_y_1[i] * wp_x[i];

        g_0_xxyyyzz_0_z_0[i] =
            g_0_yyyzz_0_z_0[i] * fi_ab_0 - g_0_yyyzz_0_z_1[i] * fti_ab_0 + g_0_xyyyzz_0_z_0[i] * pb_x + g_0_xyyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 54-57 components of targeted buffer : SKSP

    auto g_0_xxyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 54);

    auto g_0_xxyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 55);

    auto g_0_xxyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 56);

#pragma omp simd aligned(g_0_xxyyzzz_0_x_0,     \
                             g_0_xxyyzzz_0_y_0, \
                             g_0_xxyyzzz_0_z_0, \
                             g_0_xxyzzz_0_x_0,  \
                             g_0_xxyzzz_0_x_1,  \
                             g_0_xxzzz_0_x_0,   \
                             g_0_xxzzz_0_x_1,   \
                             g_0_xyyzzz_0_y_0,  \
                             g_0_xyyzzz_0_y_1,  \
                             g_0_xyyzzz_0_z_0,  \
                             g_0_xyyzzz_0_z_1,  \
                             g_0_yyzzz_0_y_0,   \
                             g_0_yyzzz_0_y_1,   \
                             g_0_yyzzz_0_z_0,   \
                             g_0_yyzzz_0_z_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_x_0[i] =
            g_0_xxzzz_0_x_0[i] * fi_ab_0 - g_0_xxzzz_0_x_1[i] * fti_ab_0 + g_0_xxyzzz_0_x_0[i] * pb_y + g_0_xxyzzz_0_x_1[i] * wp_y[i];

        g_0_xxyyzzz_0_y_0[i] =
            g_0_yyzzz_0_y_0[i] * fi_ab_0 - g_0_yyzzz_0_y_1[i] * fti_ab_0 + g_0_xyyzzz_0_y_0[i] * pb_x + g_0_xyyzzz_0_y_1[i] * wp_x[i];

        g_0_xxyyzzz_0_z_0[i] =
            g_0_yyzzz_0_z_0[i] * fi_ab_0 - g_0_yyzzz_0_z_1[i] * fti_ab_0 + g_0_xyyzzz_0_z_0[i] * pb_x + g_0_xyyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 57-60 components of targeted buffer : SKSP

    auto g_0_xxyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 57);

    auto g_0_xxyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 58);

    auto g_0_xxyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 59);

#pragma omp simd aligned(g_0_xxyzzzz_0_x_0,     \
                             g_0_xxyzzzz_0_y_0, \
                             g_0_xxyzzzz_0_z_0, \
                             g_0_xxzzzz_0_0_1,  \
                             g_0_xxzzzz_0_x_0,  \
                             g_0_xxzzzz_0_x_1,  \
                             g_0_xxzzzz_0_y_0,  \
                             g_0_xxzzzz_0_y_1,  \
                             g_0_xxzzzz_0_z_0,  \
                             g_0_xxzzzz_0_z_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_x_0[i] = g_0_xxzzzz_0_x_0[i] * pb_y + g_0_xxzzzz_0_x_1[i] * wp_y[i];

        g_0_xxyzzzz_0_y_0[i] = g_0_xxzzzz_0_0_1[i] * fi_abcd_0 + g_0_xxzzzz_0_y_0[i] * pb_y + g_0_xxzzzz_0_y_1[i] * wp_y[i];

        g_0_xxyzzzz_0_z_0[i] = g_0_xxzzzz_0_z_0[i] * pb_y + g_0_xxzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 60-63 components of targeted buffer : SKSP

    auto g_0_xxzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 60);

    auto g_0_xxzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 61);

    auto g_0_xxzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 62);

#pragma omp simd aligned(g_0_xxzzz_0_x_0,       \
                             g_0_xxzzz_0_x_1,   \
                             g_0_xxzzzz_0_x_0,  \
                             g_0_xxzzzz_0_x_1,  \
                             g_0_xxzzzzz_0_x_0, \
                             g_0_xxzzzzz_0_y_0, \
                             g_0_xxzzzzz_0_z_0, \
                             g_0_xzzzzz_0_y_0,  \
                             g_0_xzzzzz_0_y_1,  \
                             g_0_xzzzzz_0_z_0,  \
                             g_0_xzzzzz_0_z_1,  \
                             g_0_zzzzz_0_y_0,   \
                             g_0_zzzzz_0_y_1,   \
                             g_0_zzzzz_0_z_0,   \
                             g_0_zzzzz_0_z_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_x_0[i] =
            4.0 * g_0_xxzzz_0_x_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_x_1[i] * fti_ab_0 + g_0_xxzzzz_0_x_0[i] * pb_z + g_0_xxzzzz_0_x_1[i] * wp_z[i];

        g_0_xxzzzzz_0_y_0[i] =
            g_0_zzzzz_0_y_0[i] * fi_ab_0 - g_0_zzzzz_0_y_1[i] * fti_ab_0 + g_0_xzzzzz_0_y_0[i] * pb_x + g_0_xzzzzz_0_y_1[i] * wp_x[i];

        g_0_xxzzzzz_0_z_0[i] =
            g_0_zzzzz_0_z_0[i] * fi_ab_0 - g_0_zzzzz_0_z_1[i] * fti_ab_0 + g_0_xzzzzz_0_z_0[i] * pb_x + g_0_xzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 63-66 components of targeted buffer : SKSP

    auto g_0_xyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 63);

    auto g_0_xyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 64);

    auto g_0_xyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 65);

#pragma omp simd aligned(g_0_xyyyyyy_0_x_0,     \
                             g_0_xyyyyyy_0_y_0, \
                             g_0_xyyyyyy_0_z_0, \
                             g_0_yyyyyy_0_0_1,  \
                             g_0_yyyyyy_0_x_0,  \
                             g_0_yyyyyy_0_x_1,  \
                             g_0_yyyyyy_0_y_0,  \
                             g_0_yyyyyy_0_y_1,  \
                             g_0_yyyyyy_0_z_0,  \
                             g_0_yyyyyy_0_z_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_x_0[i] = g_0_yyyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyyy_0_x_0[i] * pb_x + g_0_yyyyyy_0_x_1[i] * wp_x[i];

        g_0_xyyyyyy_0_y_0[i] = g_0_yyyyyy_0_y_0[i] * pb_x + g_0_yyyyyy_0_y_1[i] * wp_x[i];

        g_0_xyyyyyy_0_z_0[i] = g_0_yyyyyy_0_z_0[i] * pb_x + g_0_yyyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 66-69 components of targeted buffer : SKSP

    auto g_0_xyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 66);

    auto g_0_xyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 67);

    auto g_0_xyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 68);

#pragma omp simd aligned(g_0_xyyyyy_0_x_0,      \
                             g_0_xyyyyy_0_x_1,  \
                             g_0_xyyyyyz_0_x_0, \
                             g_0_xyyyyyz_0_y_0, \
                             g_0_xyyyyyz_0_z_0, \
                             g_0_yyyyyz_0_y_0,  \
                             g_0_yyyyyz_0_y_1,  \
                             g_0_yyyyyz_0_z_0,  \
                             g_0_yyyyyz_0_z_1,  \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyyyyyz_0_x_0[i] = g_0_xyyyyy_0_x_0[i] * pb_z + g_0_xyyyyy_0_x_1[i] * wp_z[i];

        g_0_xyyyyyz_0_y_0[i] = g_0_yyyyyz_0_y_0[i] * pb_x + g_0_yyyyyz_0_y_1[i] * wp_x[i];

        g_0_xyyyyyz_0_z_0[i] = g_0_yyyyyz_0_z_0[i] * pb_x + g_0_yyyyyz_0_z_1[i] * wp_x[i];
    }

    /// Set up 69-72 components of targeted buffer : SKSP

    auto g_0_xyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 69);

    auto g_0_xyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 70);

    auto g_0_xyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 71);

#pragma omp simd aligned(g_0_xyyyyzz_0_x_0,     \
                             g_0_xyyyyzz_0_y_0, \
                             g_0_xyyyyzz_0_z_0, \
                             g_0_yyyyzz_0_0_1,  \
                             g_0_yyyyzz_0_x_0,  \
                             g_0_yyyyzz_0_x_1,  \
                             g_0_yyyyzz_0_y_0,  \
                             g_0_yyyyzz_0_y_1,  \
                             g_0_yyyyzz_0_z_0,  \
                             g_0_yyyyzz_0_z_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_x_0[i] = g_0_yyyyzz_0_0_1[i] * fi_abcd_0 + g_0_yyyyzz_0_x_0[i] * pb_x + g_0_yyyyzz_0_x_1[i] * wp_x[i];

        g_0_xyyyyzz_0_y_0[i] = g_0_yyyyzz_0_y_0[i] * pb_x + g_0_yyyyzz_0_y_1[i] * wp_x[i];

        g_0_xyyyyzz_0_z_0[i] = g_0_yyyyzz_0_z_0[i] * pb_x + g_0_yyyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 72-75 components of targeted buffer : SKSP

    auto g_0_xyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 72);

    auto g_0_xyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 73);

    auto g_0_xyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 74);

#pragma omp simd aligned(g_0_xyyyzzz_0_x_0,     \
                             g_0_xyyyzzz_0_y_0, \
                             g_0_xyyyzzz_0_z_0, \
                             g_0_yyyzzz_0_0_1,  \
                             g_0_yyyzzz_0_x_0,  \
                             g_0_yyyzzz_0_x_1,  \
                             g_0_yyyzzz_0_y_0,  \
                             g_0_yyyzzz_0_y_1,  \
                             g_0_yyyzzz_0_z_0,  \
                             g_0_yyyzzz_0_z_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_x_0[i] = g_0_yyyzzz_0_0_1[i] * fi_abcd_0 + g_0_yyyzzz_0_x_0[i] * pb_x + g_0_yyyzzz_0_x_1[i] * wp_x[i];

        g_0_xyyyzzz_0_y_0[i] = g_0_yyyzzz_0_y_0[i] * pb_x + g_0_yyyzzz_0_y_1[i] * wp_x[i];

        g_0_xyyyzzz_0_z_0[i] = g_0_yyyzzz_0_z_0[i] * pb_x + g_0_yyyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 75-78 components of targeted buffer : SKSP

    auto g_0_xyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 75);

    auto g_0_xyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 76);

    auto g_0_xyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 77);

#pragma omp simd aligned(g_0_xyyzzzz_0_x_0,     \
                             g_0_xyyzzzz_0_y_0, \
                             g_0_xyyzzzz_0_z_0, \
                             g_0_yyzzzz_0_0_1,  \
                             g_0_yyzzzz_0_x_0,  \
                             g_0_yyzzzz_0_x_1,  \
                             g_0_yyzzzz_0_y_0,  \
                             g_0_yyzzzz_0_y_1,  \
                             g_0_yyzzzz_0_z_0,  \
                             g_0_yyzzzz_0_z_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_x_0[i] = g_0_yyzzzz_0_0_1[i] * fi_abcd_0 + g_0_yyzzzz_0_x_0[i] * pb_x + g_0_yyzzzz_0_x_1[i] * wp_x[i];

        g_0_xyyzzzz_0_y_0[i] = g_0_yyzzzz_0_y_0[i] * pb_x + g_0_yyzzzz_0_y_1[i] * wp_x[i];

        g_0_xyyzzzz_0_z_0[i] = g_0_yyzzzz_0_z_0[i] * pb_x + g_0_yyzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 78-81 components of targeted buffer : SKSP

    auto g_0_xyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 78);

    auto g_0_xyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 79);

    auto g_0_xyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 80);

#pragma omp simd aligned(g_0_xyzzzzz_0_x_0,     \
                             g_0_xyzzzzz_0_y_0, \
                             g_0_xyzzzzz_0_z_0, \
                             g_0_xzzzzz_0_x_0,  \
                             g_0_xzzzzz_0_x_1,  \
                             g_0_yzzzzz_0_y_0,  \
                             g_0_yzzzzz_0_y_1,  \
                             g_0_yzzzzz_0_z_0,  \
                             g_0_yzzzzz_0_z_1,  \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyzzzzz_0_x_0[i] = g_0_xzzzzz_0_x_0[i] * pb_y + g_0_xzzzzz_0_x_1[i] * wp_y[i];

        g_0_xyzzzzz_0_y_0[i] = g_0_yzzzzz_0_y_0[i] * pb_x + g_0_yzzzzz_0_y_1[i] * wp_x[i];

        g_0_xyzzzzz_0_z_0[i] = g_0_yzzzzz_0_z_0[i] * pb_x + g_0_yzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 81-84 components of targeted buffer : SKSP

    auto g_0_xzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 81);

    auto g_0_xzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 82);

    auto g_0_xzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 83);

#pragma omp simd aligned(g_0_xzzzzzz_0_x_0,     \
                             g_0_xzzzzzz_0_y_0, \
                             g_0_xzzzzzz_0_z_0, \
                             g_0_zzzzzz_0_0_1,  \
                             g_0_zzzzzz_0_x_0,  \
                             g_0_zzzzzz_0_x_1,  \
                             g_0_zzzzzz_0_y_0,  \
                             g_0_zzzzzz_0_y_1,  \
                             g_0_zzzzzz_0_z_0,  \
                             g_0_zzzzzz_0_z_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_x_0[i] = g_0_zzzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzzz_0_x_0[i] * pb_x + g_0_zzzzzz_0_x_1[i] * wp_x[i];

        g_0_xzzzzzz_0_y_0[i] = g_0_zzzzzz_0_y_0[i] * pb_x + g_0_zzzzzz_0_y_1[i] * wp_x[i];

        g_0_xzzzzzz_0_z_0[i] = g_0_zzzzzz_0_z_0[i] * pb_x + g_0_zzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 84-87 components of targeted buffer : SKSP

    auto g_0_yyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 84);

    auto g_0_yyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 85);

    auto g_0_yyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 86);

#pragma omp simd aligned(g_0_yyyyy_0_x_0,       \
                             g_0_yyyyy_0_x_1,   \
                             g_0_yyyyy_0_y_0,   \
                             g_0_yyyyy_0_y_1,   \
                             g_0_yyyyy_0_z_0,   \
                             g_0_yyyyy_0_z_1,   \
                             g_0_yyyyyy_0_0_1,  \
                             g_0_yyyyyy_0_x_0,  \
                             g_0_yyyyyy_0_x_1,  \
                             g_0_yyyyyy_0_y_0,  \
                             g_0_yyyyyy_0_y_1,  \
                             g_0_yyyyyy_0_z_0,  \
                             g_0_yyyyyy_0_z_1,  \
                             g_0_yyyyyyy_0_x_0, \
                             g_0_yyyyyyy_0_y_0, \
                             g_0_yyyyyyy_0_z_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_x_0[i] =
            6.0 * g_0_yyyyy_0_x_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_x_1[i] * fti_ab_0 + g_0_yyyyyy_0_x_0[i] * pb_y + g_0_yyyyyy_0_x_1[i] * wp_y[i];

        g_0_yyyyyyy_0_y_0[i] = 6.0 * g_0_yyyyy_0_y_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyyy_0_0_1[i] * fi_abcd_0 +
                               g_0_yyyyyy_0_y_0[i] * pb_y + g_0_yyyyyy_0_y_1[i] * wp_y[i];

        g_0_yyyyyyy_0_z_0[i] =
            6.0 * g_0_yyyyy_0_z_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_z_1[i] * fti_ab_0 + g_0_yyyyyy_0_z_0[i] * pb_y + g_0_yyyyyy_0_z_1[i] * wp_y[i];
    }

    /// Set up 87-90 components of targeted buffer : SKSP

    auto g_0_yyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 87);

    auto g_0_yyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 88);

    auto g_0_yyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 89);

#pragma omp simd aligned(g_0_yyyyyy_0_0_1,      \
                             g_0_yyyyyy_0_x_0,  \
                             g_0_yyyyyy_0_x_1,  \
                             g_0_yyyyyy_0_y_0,  \
                             g_0_yyyyyy_0_y_1,  \
                             g_0_yyyyyy_0_z_0,  \
                             g_0_yyyyyy_0_z_1,  \
                             g_0_yyyyyyz_0_x_0, \
                             g_0_yyyyyyz_0_y_0, \
                             g_0_yyyyyyz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_x_0[i] = g_0_yyyyyy_0_x_0[i] * pb_z + g_0_yyyyyy_0_x_1[i] * wp_z[i];

        g_0_yyyyyyz_0_y_0[i] = g_0_yyyyyy_0_y_0[i] * pb_z + g_0_yyyyyy_0_y_1[i] * wp_z[i];

        g_0_yyyyyyz_0_z_0[i] = g_0_yyyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyyy_0_z_0[i] * pb_z + g_0_yyyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 90-93 components of targeted buffer : SKSP

    auto g_0_yyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 90);

    auto g_0_yyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 91);

    auto g_0_yyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 92);

#pragma omp simd aligned(g_0_yyyyy_0_y_0,       \
                             g_0_yyyyy_0_y_1,   \
                             g_0_yyyyyz_0_y_0,  \
                             g_0_yyyyyz_0_y_1,  \
                             g_0_yyyyyzz_0_x_0, \
                             g_0_yyyyyzz_0_y_0, \
                             g_0_yyyyyzz_0_z_0, \
                             g_0_yyyyzz_0_x_0,  \
                             g_0_yyyyzz_0_x_1,  \
                             g_0_yyyyzz_0_z_0,  \
                             g_0_yyyyzz_0_z_1,  \
                             g_0_yyyzz_0_x_0,   \
                             g_0_yyyzz_0_x_1,   \
                             g_0_yyyzz_0_z_0,   \
                             g_0_yyyzz_0_z_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_x_0[i] =
            4.0 * g_0_yyyzz_0_x_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_x_1[i] * fti_ab_0 + g_0_yyyyzz_0_x_0[i] * pb_y + g_0_yyyyzz_0_x_1[i] * wp_y[i];

        g_0_yyyyyzz_0_y_0[i] =
            g_0_yyyyy_0_y_0[i] * fi_ab_0 - g_0_yyyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyyz_0_y_0[i] * pb_z + g_0_yyyyyz_0_y_1[i] * wp_z[i];

        g_0_yyyyyzz_0_z_0[i] =
            4.0 * g_0_yyyzz_0_z_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_z_1[i] * fti_ab_0 + g_0_yyyyzz_0_z_0[i] * pb_y + g_0_yyyyzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 93-96 components of targeted buffer : SKSP

    auto g_0_yyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 93);

    auto g_0_yyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 94);

    auto g_0_yyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 95);

#pragma omp simd aligned(g_0_yyyyz_0_y_0,       \
                             g_0_yyyyz_0_y_1,   \
                             g_0_yyyyzz_0_y_0,  \
                             g_0_yyyyzz_0_y_1,  \
                             g_0_yyyyzzz_0_x_0, \
                             g_0_yyyyzzz_0_y_0, \
                             g_0_yyyyzzz_0_z_0, \
                             g_0_yyyzzz_0_x_0,  \
                             g_0_yyyzzz_0_x_1,  \
                             g_0_yyyzzz_0_z_0,  \
                             g_0_yyyzzz_0_z_1,  \
                             g_0_yyzzz_0_x_0,   \
                             g_0_yyzzz_0_x_1,   \
                             g_0_yyzzz_0_z_0,   \
                             g_0_yyzzz_0_z_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_x_0[i] =
            3.0 * g_0_yyzzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_x_1[i] * fti_ab_0 + g_0_yyyzzz_0_x_0[i] * pb_y + g_0_yyyzzz_0_x_1[i] * wp_y[i];

        g_0_yyyyzzz_0_y_0[i] =
            2.0 * g_0_yyyyz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_y_1[i] * fti_ab_0 + g_0_yyyyzz_0_y_0[i] * pb_z + g_0_yyyyzz_0_y_1[i] * wp_z[i];

        g_0_yyyyzzz_0_z_0[i] =
            3.0 * g_0_yyzzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_z_1[i] * fti_ab_0 + g_0_yyyzzz_0_z_0[i] * pb_y + g_0_yyyzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 96-99 components of targeted buffer : SKSP

    auto g_0_yyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 96);

    auto g_0_yyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 97);

    auto g_0_yyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 98);

#pragma omp simd aligned(g_0_yyyzz_0_y_0,       \
                             g_0_yyyzz_0_y_1,   \
                             g_0_yyyzzz_0_y_0,  \
                             g_0_yyyzzz_0_y_1,  \
                             g_0_yyyzzzz_0_x_0, \
                             g_0_yyyzzzz_0_y_0, \
                             g_0_yyyzzzz_0_z_0, \
                             g_0_yyzzzz_0_x_0,  \
                             g_0_yyzzzz_0_x_1,  \
                             g_0_yyzzzz_0_z_0,  \
                             g_0_yyzzzz_0_z_1,  \
                             g_0_yzzzz_0_x_0,   \
                             g_0_yzzzz_0_x_1,   \
                             g_0_yzzzz_0_z_0,   \
                             g_0_yzzzz_0_z_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_x_0[i] =
            2.0 * g_0_yzzzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_x_1[i] * fti_ab_0 + g_0_yyzzzz_0_x_0[i] * pb_y + g_0_yyzzzz_0_x_1[i] * wp_y[i];

        g_0_yyyzzzz_0_y_0[i] =
            3.0 * g_0_yyyzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_y_1[i] * fti_ab_0 + g_0_yyyzzz_0_y_0[i] * pb_z + g_0_yyyzzz_0_y_1[i] * wp_z[i];

        g_0_yyyzzzz_0_z_0[i] =
            2.0 * g_0_yzzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_z_1[i] * fti_ab_0 + g_0_yyzzzz_0_z_0[i] * pb_y + g_0_yyzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 99-102 components of targeted buffer : SKSP

    auto g_0_yyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 99);

    auto g_0_yyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 100);

    auto g_0_yyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 101);

#pragma omp simd aligned(g_0_yyzzz_0_y_0,       \
                             g_0_yyzzz_0_y_1,   \
                             g_0_yyzzzz_0_y_0,  \
                             g_0_yyzzzz_0_y_1,  \
                             g_0_yyzzzzz_0_x_0, \
                             g_0_yyzzzzz_0_y_0, \
                             g_0_yyzzzzz_0_z_0, \
                             g_0_yzzzzz_0_x_0,  \
                             g_0_yzzzzz_0_x_1,  \
                             g_0_yzzzzz_0_z_0,  \
                             g_0_yzzzzz_0_z_1,  \
                             g_0_zzzzz_0_x_0,   \
                             g_0_zzzzz_0_x_1,   \
                             g_0_zzzzz_0_z_0,   \
                             g_0_zzzzz_0_z_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_x_0[i] =
            g_0_zzzzz_0_x_0[i] * fi_ab_0 - g_0_zzzzz_0_x_1[i] * fti_ab_0 + g_0_yzzzzz_0_x_0[i] * pb_y + g_0_yzzzzz_0_x_1[i] * wp_y[i];

        g_0_yyzzzzz_0_y_0[i] =
            4.0 * g_0_yyzzz_0_y_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_y_1[i] * fti_ab_0 + g_0_yyzzzz_0_y_0[i] * pb_z + g_0_yyzzzz_0_y_1[i] * wp_z[i];

        g_0_yyzzzzz_0_z_0[i] =
            g_0_zzzzz_0_z_0[i] * fi_ab_0 - g_0_zzzzz_0_z_1[i] * fti_ab_0 + g_0_yzzzzz_0_z_0[i] * pb_y + g_0_yzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 102-105 components of targeted buffer : SKSP

    auto g_0_yzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 102);

    auto g_0_yzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 103);

    auto g_0_yzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 104);

#pragma omp simd aligned(g_0_yzzzzzz_0_x_0,     \
                             g_0_yzzzzzz_0_y_0, \
                             g_0_yzzzzzz_0_z_0, \
                             g_0_zzzzzz_0_0_1,  \
                             g_0_zzzzzz_0_x_0,  \
                             g_0_zzzzzz_0_x_1,  \
                             g_0_zzzzzz_0_y_0,  \
                             g_0_zzzzzz_0_y_1,  \
                             g_0_zzzzzz_0_z_0,  \
                             g_0_zzzzzz_0_z_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_x_0[i] = g_0_zzzzzz_0_x_0[i] * pb_y + g_0_zzzzzz_0_x_1[i] * wp_y[i];

        g_0_yzzzzzz_0_y_0[i] = g_0_zzzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzzz_0_y_0[i] * pb_y + g_0_zzzzzz_0_y_1[i] * wp_y[i];

        g_0_yzzzzzz_0_z_0[i] = g_0_zzzzzz_0_z_0[i] * pb_y + g_0_zzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 105-108 components of targeted buffer : SKSP

    auto g_0_zzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 105);

    auto g_0_zzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 106);

    auto g_0_zzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 107);

#pragma omp simd aligned(g_0_zzzzz_0_x_0,       \
                             g_0_zzzzz_0_x_1,   \
                             g_0_zzzzz_0_y_0,   \
                             g_0_zzzzz_0_y_1,   \
                             g_0_zzzzz_0_z_0,   \
                             g_0_zzzzz_0_z_1,   \
                             g_0_zzzzzz_0_0_1,  \
                             g_0_zzzzzz_0_x_0,  \
                             g_0_zzzzzz_0_x_1,  \
                             g_0_zzzzzz_0_y_0,  \
                             g_0_zzzzzz_0_y_1,  \
                             g_0_zzzzzz_0_z_0,  \
                             g_0_zzzzzz_0_z_1,  \
                             g_0_zzzzzzz_0_x_0, \
                             g_0_zzzzzzz_0_y_0, \
                             g_0_zzzzzzz_0_z_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_x_0[i] =
            6.0 * g_0_zzzzz_0_x_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_x_1[i] * fti_ab_0 + g_0_zzzzzz_0_x_0[i] * pb_z + g_0_zzzzzz_0_x_1[i] * wp_z[i];

        g_0_zzzzzzz_0_y_0[i] =
            6.0 * g_0_zzzzz_0_y_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_y_1[i] * fti_ab_0 + g_0_zzzzzz_0_y_0[i] * pb_z + g_0_zzzzzz_0_y_1[i] * wp_z[i];

        g_0_zzzzzzz_0_z_0[i] = 6.0 * g_0_zzzzz_0_z_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_z_1[i] * fti_ab_0 + g_0_zzzzzz_0_0_1[i] * fi_abcd_0 +
                               g_0_zzzzzz_0_z_0[i] * pb_z + g_0_zzzzzz_0_z_1[i] * wp_z[i];
    }
}

}  // namespace erirec

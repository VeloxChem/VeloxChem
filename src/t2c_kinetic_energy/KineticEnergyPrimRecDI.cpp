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

#include "KineticEnergyPrimRecDI.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_di(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_di,
                            const size_t              idx_ovl_si,
                            const size_t              idx_kin_si,
                            const size_t              idx_kin_ph,
                            const size_t              idx_kin_pi,
                            const size_t              idx_ovl_di,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SI

    auto ts_0_xxxxxx = pbuffer.data(idx_ovl_si);

    auto ts_0_xxxxxy = pbuffer.data(idx_ovl_si + 1);

    auto ts_0_xxxxxz = pbuffer.data(idx_ovl_si + 2);

    auto ts_0_xxxxyy = pbuffer.data(idx_ovl_si + 3);

    auto ts_0_xxxxyz = pbuffer.data(idx_ovl_si + 4);

    auto ts_0_xxxxzz = pbuffer.data(idx_ovl_si + 5);

    auto ts_0_xxxyyy = pbuffer.data(idx_ovl_si + 6);

    auto ts_0_xxxyyz = pbuffer.data(idx_ovl_si + 7);

    auto ts_0_xxxyzz = pbuffer.data(idx_ovl_si + 8);

    auto ts_0_xxxzzz = pbuffer.data(idx_ovl_si + 9);

    auto ts_0_xxyyyy = pbuffer.data(idx_ovl_si + 10);

    auto ts_0_xxyyyz = pbuffer.data(idx_ovl_si + 11);

    auto ts_0_xxyyzz = pbuffer.data(idx_ovl_si + 12);

    auto ts_0_xxyzzz = pbuffer.data(idx_ovl_si + 13);

    auto ts_0_xxzzzz = pbuffer.data(idx_ovl_si + 14);

    auto ts_0_xyyyyy = pbuffer.data(idx_ovl_si + 15);

    auto ts_0_xyyyyz = pbuffer.data(idx_ovl_si + 16);

    auto ts_0_xyyyzz = pbuffer.data(idx_ovl_si + 17);

    auto ts_0_xyyzzz = pbuffer.data(idx_ovl_si + 18);

    auto ts_0_xyzzzz = pbuffer.data(idx_ovl_si + 19);

    auto ts_0_xzzzzz = pbuffer.data(idx_ovl_si + 20);

    auto ts_0_yyyyyy = pbuffer.data(idx_ovl_si + 21);

    auto ts_0_yyyyyz = pbuffer.data(idx_ovl_si + 22);

    auto ts_0_yyyyzz = pbuffer.data(idx_ovl_si + 23);

    auto ts_0_yyyzzz = pbuffer.data(idx_ovl_si + 24);

    auto ts_0_yyzzzz = pbuffer.data(idx_ovl_si + 25);

    auto ts_0_yzzzzz = pbuffer.data(idx_ovl_si + 26);

    auto ts_0_zzzzzz = pbuffer.data(idx_ovl_si + 27);

    // Set up components of auxiliary buffer : SI

    auto tk_0_xxxxxx = pbuffer.data(idx_kin_si);

    auto tk_0_xxxxxy = pbuffer.data(idx_kin_si + 1);

    auto tk_0_xxxxxz = pbuffer.data(idx_kin_si + 2);

    auto tk_0_xxxxyy = pbuffer.data(idx_kin_si + 3);

    auto tk_0_xxxxyz = pbuffer.data(idx_kin_si + 4);

    auto tk_0_xxxxzz = pbuffer.data(idx_kin_si + 5);

    auto tk_0_xxxyyy = pbuffer.data(idx_kin_si + 6);

    auto tk_0_xxxyyz = pbuffer.data(idx_kin_si + 7);

    auto tk_0_xxxyzz = pbuffer.data(idx_kin_si + 8);

    auto tk_0_xxxzzz = pbuffer.data(idx_kin_si + 9);

    auto tk_0_xxyyyy = pbuffer.data(idx_kin_si + 10);

    auto tk_0_xxyyyz = pbuffer.data(idx_kin_si + 11);

    auto tk_0_xxyyzz = pbuffer.data(idx_kin_si + 12);

    auto tk_0_xxyzzz = pbuffer.data(idx_kin_si + 13);

    auto tk_0_xxzzzz = pbuffer.data(idx_kin_si + 14);

    auto tk_0_xyyyyy = pbuffer.data(idx_kin_si + 15);

    auto tk_0_xyyyyz = pbuffer.data(idx_kin_si + 16);

    auto tk_0_xyyyzz = pbuffer.data(idx_kin_si + 17);

    auto tk_0_xyyzzz = pbuffer.data(idx_kin_si + 18);

    auto tk_0_xyzzzz = pbuffer.data(idx_kin_si + 19);

    auto tk_0_xzzzzz = pbuffer.data(idx_kin_si + 20);

    auto tk_0_yyyyyy = pbuffer.data(idx_kin_si + 21);

    auto tk_0_yyyyyz = pbuffer.data(idx_kin_si + 22);

    auto tk_0_yyyyzz = pbuffer.data(idx_kin_si + 23);

    auto tk_0_yyyzzz = pbuffer.data(idx_kin_si + 24);

    auto tk_0_yyzzzz = pbuffer.data(idx_kin_si + 25);

    auto tk_0_yzzzzz = pbuffer.data(idx_kin_si + 26);

    auto tk_0_zzzzzz = pbuffer.data(idx_kin_si + 27);

    // Set up components of auxiliary buffer : PH

    auto tk_x_xxxxx = pbuffer.data(idx_kin_ph);

    auto tk_x_xxxxy = pbuffer.data(idx_kin_ph + 1);

    auto tk_x_xxxxz = pbuffer.data(idx_kin_ph + 2);

    auto tk_x_xxxyy = pbuffer.data(idx_kin_ph + 3);

    auto tk_x_xxxyz = pbuffer.data(idx_kin_ph + 4);

    auto tk_x_xxxzz = pbuffer.data(idx_kin_ph + 5);

    auto tk_x_xxyyy = pbuffer.data(idx_kin_ph + 6);

    auto tk_x_xxyyz = pbuffer.data(idx_kin_ph + 7);

    auto tk_x_xxyzz = pbuffer.data(idx_kin_ph + 8);

    auto tk_x_xxzzz = pbuffer.data(idx_kin_ph + 9);

    auto tk_x_xyyyy = pbuffer.data(idx_kin_ph + 10);

    auto tk_x_xyyyz = pbuffer.data(idx_kin_ph + 11);

    auto tk_x_xyyzz = pbuffer.data(idx_kin_ph + 12);

    auto tk_x_xyzzz = pbuffer.data(idx_kin_ph + 13);

    auto tk_x_xzzzz = pbuffer.data(idx_kin_ph + 14);

    auto tk_x_yyyyy = pbuffer.data(idx_kin_ph + 15);

    auto tk_x_yyyyz = pbuffer.data(idx_kin_ph + 16);

    auto tk_x_yyyzz = pbuffer.data(idx_kin_ph + 17);

    auto tk_x_yyzzz = pbuffer.data(idx_kin_ph + 18);

    auto tk_x_yzzzz = pbuffer.data(idx_kin_ph + 19);

    auto tk_x_zzzzz = pbuffer.data(idx_kin_ph + 20);

    auto tk_y_xxxxx = pbuffer.data(idx_kin_ph + 21);

    auto tk_y_xxxxy = pbuffer.data(idx_kin_ph + 22);

    auto tk_y_xxxxz = pbuffer.data(idx_kin_ph + 23);

    auto tk_y_xxxyy = pbuffer.data(idx_kin_ph + 24);

    auto tk_y_xxxyz = pbuffer.data(idx_kin_ph + 25);

    auto tk_y_xxxzz = pbuffer.data(idx_kin_ph + 26);

    auto tk_y_xxyyy = pbuffer.data(idx_kin_ph + 27);

    auto tk_y_xxyyz = pbuffer.data(idx_kin_ph + 28);

    auto tk_y_xxyzz = pbuffer.data(idx_kin_ph + 29);

    auto tk_y_xxzzz = pbuffer.data(idx_kin_ph + 30);

    auto tk_y_xyyyy = pbuffer.data(idx_kin_ph + 31);

    auto tk_y_xyyyz = pbuffer.data(idx_kin_ph + 32);

    auto tk_y_xyyzz = pbuffer.data(idx_kin_ph + 33);

    auto tk_y_xyzzz = pbuffer.data(idx_kin_ph + 34);

    auto tk_y_xzzzz = pbuffer.data(idx_kin_ph + 35);

    auto tk_y_yyyyy = pbuffer.data(idx_kin_ph + 36);

    auto tk_y_yyyyz = pbuffer.data(idx_kin_ph + 37);

    auto tk_y_yyyzz = pbuffer.data(idx_kin_ph + 38);

    auto tk_y_yyzzz = pbuffer.data(idx_kin_ph + 39);

    auto tk_y_yzzzz = pbuffer.data(idx_kin_ph + 40);

    auto tk_y_zzzzz = pbuffer.data(idx_kin_ph + 41);

    auto tk_z_xxxxx = pbuffer.data(idx_kin_ph + 42);

    auto tk_z_xxxxy = pbuffer.data(idx_kin_ph + 43);

    auto tk_z_xxxxz = pbuffer.data(idx_kin_ph + 44);

    auto tk_z_xxxyy = pbuffer.data(idx_kin_ph + 45);

    auto tk_z_xxxyz = pbuffer.data(idx_kin_ph + 46);

    auto tk_z_xxxzz = pbuffer.data(idx_kin_ph + 47);

    auto tk_z_xxyyy = pbuffer.data(idx_kin_ph + 48);

    auto tk_z_xxyyz = pbuffer.data(idx_kin_ph + 49);

    auto tk_z_xxyzz = pbuffer.data(idx_kin_ph + 50);

    auto tk_z_xxzzz = pbuffer.data(idx_kin_ph + 51);

    auto tk_z_xyyyy = pbuffer.data(idx_kin_ph + 52);

    auto tk_z_xyyyz = pbuffer.data(idx_kin_ph + 53);

    auto tk_z_xyyzz = pbuffer.data(idx_kin_ph + 54);

    auto tk_z_xyzzz = pbuffer.data(idx_kin_ph + 55);

    auto tk_z_xzzzz = pbuffer.data(idx_kin_ph + 56);

    auto tk_z_yyyyy = pbuffer.data(idx_kin_ph + 57);

    auto tk_z_yyyyz = pbuffer.data(idx_kin_ph + 58);

    auto tk_z_yyyzz = pbuffer.data(idx_kin_ph + 59);

    auto tk_z_yyzzz = pbuffer.data(idx_kin_ph + 60);

    auto tk_z_yzzzz = pbuffer.data(idx_kin_ph + 61);

    auto tk_z_zzzzz = pbuffer.data(idx_kin_ph + 62);

    // Set up components of auxiliary buffer : PI

    auto tk_x_xxxxxx = pbuffer.data(idx_kin_pi);

    auto tk_x_xxxxxy = pbuffer.data(idx_kin_pi + 1);

    auto tk_x_xxxxxz = pbuffer.data(idx_kin_pi + 2);

    auto tk_x_xxxxyy = pbuffer.data(idx_kin_pi + 3);

    auto tk_x_xxxxyz = pbuffer.data(idx_kin_pi + 4);

    auto tk_x_xxxxzz = pbuffer.data(idx_kin_pi + 5);

    auto tk_x_xxxyyy = pbuffer.data(idx_kin_pi + 6);

    auto tk_x_xxxyyz = pbuffer.data(idx_kin_pi + 7);

    auto tk_x_xxxyzz = pbuffer.data(idx_kin_pi + 8);

    auto tk_x_xxxzzz = pbuffer.data(idx_kin_pi + 9);

    auto tk_x_xxyyyy = pbuffer.data(idx_kin_pi + 10);

    auto tk_x_xxyyyz = pbuffer.data(idx_kin_pi + 11);

    auto tk_x_xxyyzz = pbuffer.data(idx_kin_pi + 12);

    auto tk_x_xxyzzz = pbuffer.data(idx_kin_pi + 13);

    auto tk_x_xxzzzz = pbuffer.data(idx_kin_pi + 14);

    auto tk_x_xyyyyy = pbuffer.data(idx_kin_pi + 15);

    auto tk_x_xyyyyz = pbuffer.data(idx_kin_pi + 16);

    auto tk_x_xyyyzz = pbuffer.data(idx_kin_pi + 17);

    auto tk_x_xyyzzz = pbuffer.data(idx_kin_pi + 18);

    auto tk_x_xyzzzz = pbuffer.data(idx_kin_pi + 19);

    auto tk_x_xzzzzz = pbuffer.data(idx_kin_pi + 20);

    auto tk_x_yyyyyy = pbuffer.data(idx_kin_pi + 21);

    auto tk_x_yyyyyz = pbuffer.data(idx_kin_pi + 22);

    auto tk_x_yyyyzz = pbuffer.data(idx_kin_pi + 23);

    auto tk_x_yyyzzz = pbuffer.data(idx_kin_pi + 24);

    auto tk_x_yyzzzz = pbuffer.data(idx_kin_pi + 25);

    auto tk_x_yzzzzz = pbuffer.data(idx_kin_pi + 26);

    auto tk_x_zzzzzz = pbuffer.data(idx_kin_pi + 27);

    auto tk_y_xxxxxx = pbuffer.data(idx_kin_pi + 28);

    auto tk_y_xxxxxy = pbuffer.data(idx_kin_pi + 29);

    auto tk_y_xxxxxz = pbuffer.data(idx_kin_pi + 30);

    auto tk_y_xxxxyy = pbuffer.data(idx_kin_pi + 31);

    auto tk_y_xxxxyz = pbuffer.data(idx_kin_pi + 32);

    auto tk_y_xxxxzz = pbuffer.data(idx_kin_pi + 33);

    auto tk_y_xxxyyy = pbuffer.data(idx_kin_pi + 34);

    auto tk_y_xxxyyz = pbuffer.data(idx_kin_pi + 35);

    auto tk_y_xxxyzz = pbuffer.data(idx_kin_pi + 36);

    auto tk_y_xxxzzz = pbuffer.data(idx_kin_pi + 37);

    auto tk_y_xxyyyy = pbuffer.data(idx_kin_pi + 38);

    auto tk_y_xxyyyz = pbuffer.data(idx_kin_pi + 39);

    auto tk_y_xxyyzz = pbuffer.data(idx_kin_pi + 40);

    auto tk_y_xxyzzz = pbuffer.data(idx_kin_pi + 41);

    auto tk_y_xxzzzz = pbuffer.data(idx_kin_pi + 42);

    auto tk_y_xyyyyy = pbuffer.data(idx_kin_pi + 43);

    auto tk_y_xyyyyz = pbuffer.data(idx_kin_pi + 44);

    auto tk_y_xyyyzz = pbuffer.data(idx_kin_pi + 45);

    auto tk_y_xyyzzz = pbuffer.data(idx_kin_pi + 46);

    auto tk_y_xyzzzz = pbuffer.data(idx_kin_pi + 47);

    auto tk_y_xzzzzz = pbuffer.data(idx_kin_pi + 48);

    auto tk_y_yyyyyy = pbuffer.data(idx_kin_pi + 49);

    auto tk_y_yyyyyz = pbuffer.data(idx_kin_pi + 50);

    auto tk_y_yyyyzz = pbuffer.data(idx_kin_pi + 51);

    auto tk_y_yyyzzz = pbuffer.data(idx_kin_pi + 52);

    auto tk_y_yyzzzz = pbuffer.data(idx_kin_pi + 53);

    auto tk_y_yzzzzz = pbuffer.data(idx_kin_pi + 54);

    auto tk_y_zzzzzz = pbuffer.data(idx_kin_pi + 55);

    auto tk_z_xxxxxx = pbuffer.data(idx_kin_pi + 56);

    auto tk_z_xxxxxy = pbuffer.data(idx_kin_pi + 57);

    auto tk_z_xxxxxz = pbuffer.data(idx_kin_pi + 58);

    auto tk_z_xxxxyy = pbuffer.data(idx_kin_pi + 59);

    auto tk_z_xxxxyz = pbuffer.data(idx_kin_pi + 60);

    auto tk_z_xxxxzz = pbuffer.data(idx_kin_pi + 61);

    auto tk_z_xxxyyy = pbuffer.data(idx_kin_pi + 62);

    auto tk_z_xxxyyz = pbuffer.data(idx_kin_pi + 63);

    auto tk_z_xxxyzz = pbuffer.data(idx_kin_pi + 64);

    auto tk_z_xxxzzz = pbuffer.data(idx_kin_pi + 65);

    auto tk_z_xxyyyy = pbuffer.data(idx_kin_pi + 66);

    auto tk_z_xxyyyz = pbuffer.data(idx_kin_pi + 67);

    auto tk_z_xxyyzz = pbuffer.data(idx_kin_pi + 68);

    auto tk_z_xxyzzz = pbuffer.data(idx_kin_pi + 69);

    auto tk_z_xxzzzz = pbuffer.data(idx_kin_pi + 70);

    auto tk_z_xyyyyy = pbuffer.data(idx_kin_pi + 71);

    auto tk_z_xyyyyz = pbuffer.data(idx_kin_pi + 72);

    auto tk_z_xyyyzz = pbuffer.data(idx_kin_pi + 73);

    auto tk_z_xyyzzz = pbuffer.data(idx_kin_pi + 74);

    auto tk_z_xyzzzz = pbuffer.data(idx_kin_pi + 75);

    auto tk_z_xzzzzz = pbuffer.data(idx_kin_pi + 76);

    auto tk_z_yyyyyy = pbuffer.data(idx_kin_pi + 77);

    auto tk_z_yyyyyz = pbuffer.data(idx_kin_pi + 78);

    auto tk_z_yyyyzz = pbuffer.data(idx_kin_pi + 79);

    auto tk_z_yyyzzz = pbuffer.data(idx_kin_pi + 80);

    auto tk_z_yyzzzz = pbuffer.data(idx_kin_pi + 81);

    auto tk_z_yzzzzz = pbuffer.data(idx_kin_pi + 82);

    auto tk_z_zzzzzz = pbuffer.data(idx_kin_pi + 83);

    // Set up components of auxiliary buffer : DI

    auto ts_xx_xxxxxx = pbuffer.data(idx_ovl_di);

    auto ts_xx_xxxxxy = pbuffer.data(idx_ovl_di + 1);

    auto ts_xx_xxxxxz = pbuffer.data(idx_ovl_di + 2);

    auto ts_xx_xxxxyy = pbuffer.data(idx_ovl_di + 3);

    auto ts_xx_xxxxyz = pbuffer.data(idx_ovl_di + 4);

    auto ts_xx_xxxxzz = pbuffer.data(idx_ovl_di + 5);

    auto ts_xx_xxxyyy = pbuffer.data(idx_ovl_di + 6);

    auto ts_xx_xxxyyz = pbuffer.data(idx_ovl_di + 7);

    auto ts_xx_xxxyzz = pbuffer.data(idx_ovl_di + 8);

    auto ts_xx_xxxzzz = pbuffer.data(idx_ovl_di + 9);

    auto ts_xx_xxyyyy = pbuffer.data(idx_ovl_di + 10);

    auto ts_xx_xxyyyz = pbuffer.data(idx_ovl_di + 11);

    auto ts_xx_xxyyzz = pbuffer.data(idx_ovl_di + 12);

    auto ts_xx_xxyzzz = pbuffer.data(idx_ovl_di + 13);

    auto ts_xx_xxzzzz = pbuffer.data(idx_ovl_di + 14);

    auto ts_xx_xyyyyy = pbuffer.data(idx_ovl_di + 15);

    auto ts_xx_xyyyyz = pbuffer.data(idx_ovl_di + 16);

    auto ts_xx_xyyyzz = pbuffer.data(idx_ovl_di + 17);

    auto ts_xx_xyyzzz = pbuffer.data(idx_ovl_di + 18);

    auto ts_xx_xyzzzz = pbuffer.data(idx_ovl_di + 19);

    auto ts_xx_xzzzzz = pbuffer.data(idx_ovl_di + 20);

    auto ts_xx_yyyyyy = pbuffer.data(idx_ovl_di + 21);

    auto ts_xx_yyyyyz = pbuffer.data(idx_ovl_di + 22);

    auto ts_xx_yyyyzz = pbuffer.data(idx_ovl_di + 23);

    auto ts_xx_yyyzzz = pbuffer.data(idx_ovl_di + 24);

    auto ts_xx_yyzzzz = pbuffer.data(idx_ovl_di + 25);

    auto ts_xx_yzzzzz = pbuffer.data(idx_ovl_di + 26);

    auto ts_xx_zzzzzz = pbuffer.data(idx_ovl_di + 27);

    auto ts_xy_xxxxxx = pbuffer.data(idx_ovl_di + 28);

    auto ts_xy_xxxxxy = pbuffer.data(idx_ovl_di + 29);

    auto ts_xy_xxxxxz = pbuffer.data(idx_ovl_di + 30);

    auto ts_xy_xxxxyy = pbuffer.data(idx_ovl_di + 31);

    auto ts_xy_xxxxyz = pbuffer.data(idx_ovl_di + 32);

    auto ts_xy_xxxxzz = pbuffer.data(idx_ovl_di + 33);

    auto ts_xy_xxxyyy = pbuffer.data(idx_ovl_di + 34);

    auto ts_xy_xxxyyz = pbuffer.data(idx_ovl_di + 35);

    auto ts_xy_xxxyzz = pbuffer.data(idx_ovl_di + 36);

    auto ts_xy_xxxzzz = pbuffer.data(idx_ovl_di + 37);

    auto ts_xy_xxyyyy = pbuffer.data(idx_ovl_di + 38);

    auto ts_xy_xxyyyz = pbuffer.data(idx_ovl_di + 39);

    auto ts_xy_xxyyzz = pbuffer.data(idx_ovl_di + 40);

    auto ts_xy_xxyzzz = pbuffer.data(idx_ovl_di + 41);

    auto ts_xy_xxzzzz = pbuffer.data(idx_ovl_di + 42);

    auto ts_xy_xyyyyy = pbuffer.data(idx_ovl_di + 43);

    auto ts_xy_xyyyyz = pbuffer.data(idx_ovl_di + 44);

    auto ts_xy_xyyyzz = pbuffer.data(idx_ovl_di + 45);

    auto ts_xy_xyyzzz = pbuffer.data(idx_ovl_di + 46);

    auto ts_xy_xyzzzz = pbuffer.data(idx_ovl_di + 47);

    auto ts_xy_xzzzzz = pbuffer.data(idx_ovl_di + 48);

    auto ts_xy_yyyyyy = pbuffer.data(idx_ovl_di + 49);

    auto ts_xy_yyyyyz = pbuffer.data(idx_ovl_di + 50);

    auto ts_xy_yyyyzz = pbuffer.data(idx_ovl_di + 51);

    auto ts_xy_yyyzzz = pbuffer.data(idx_ovl_di + 52);

    auto ts_xy_yyzzzz = pbuffer.data(idx_ovl_di + 53);

    auto ts_xy_yzzzzz = pbuffer.data(idx_ovl_di + 54);

    auto ts_xy_zzzzzz = pbuffer.data(idx_ovl_di + 55);

    auto ts_xz_xxxxxx = pbuffer.data(idx_ovl_di + 56);

    auto ts_xz_xxxxxy = pbuffer.data(idx_ovl_di + 57);

    auto ts_xz_xxxxxz = pbuffer.data(idx_ovl_di + 58);

    auto ts_xz_xxxxyy = pbuffer.data(idx_ovl_di + 59);

    auto ts_xz_xxxxyz = pbuffer.data(idx_ovl_di + 60);

    auto ts_xz_xxxxzz = pbuffer.data(idx_ovl_di + 61);

    auto ts_xz_xxxyyy = pbuffer.data(idx_ovl_di + 62);

    auto ts_xz_xxxyyz = pbuffer.data(idx_ovl_di + 63);

    auto ts_xz_xxxyzz = pbuffer.data(idx_ovl_di + 64);

    auto ts_xz_xxxzzz = pbuffer.data(idx_ovl_di + 65);

    auto ts_xz_xxyyyy = pbuffer.data(idx_ovl_di + 66);

    auto ts_xz_xxyyyz = pbuffer.data(idx_ovl_di + 67);

    auto ts_xz_xxyyzz = pbuffer.data(idx_ovl_di + 68);

    auto ts_xz_xxyzzz = pbuffer.data(idx_ovl_di + 69);

    auto ts_xz_xxzzzz = pbuffer.data(idx_ovl_di + 70);

    auto ts_xz_xyyyyy = pbuffer.data(idx_ovl_di + 71);

    auto ts_xz_xyyyyz = pbuffer.data(idx_ovl_di + 72);

    auto ts_xz_xyyyzz = pbuffer.data(idx_ovl_di + 73);

    auto ts_xz_xyyzzz = pbuffer.data(idx_ovl_di + 74);

    auto ts_xz_xyzzzz = pbuffer.data(idx_ovl_di + 75);

    auto ts_xz_xzzzzz = pbuffer.data(idx_ovl_di + 76);

    auto ts_xz_yyyyyy = pbuffer.data(idx_ovl_di + 77);

    auto ts_xz_yyyyyz = pbuffer.data(idx_ovl_di + 78);

    auto ts_xz_yyyyzz = pbuffer.data(idx_ovl_di + 79);

    auto ts_xz_yyyzzz = pbuffer.data(idx_ovl_di + 80);

    auto ts_xz_yyzzzz = pbuffer.data(idx_ovl_di + 81);

    auto ts_xz_yzzzzz = pbuffer.data(idx_ovl_di + 82);

    auto ts_xz_zzzzzz = pbuffer.data(idx_ovl_di + 83);

    auto ts_yy_xxxxxx = pbuffer.data(idx_ovl_di + 84);

    auto ts_yy_xxxxxy = pbuffer.data(idx_ovl_di + 85);

    auto ts_yy_xxxxxz = pbuffer.data(idx_ovl_di + 86);

    auto ts_yy_xxxxyy = pbuffer.data(idx_ovl_di + 87);

    auto ts_yy_xxxxyz = pbuffer.data(idx_ovl_di + 88);

    auto ts_yy_xxxxzz = pbuffer.data(idx_ovl_di + 89);

    auto ts_yy_xxxyyy = pbuffer.data(idx_ovl_di + 90);

    auto ts_yy_xxxyyz = pbuffer.data(idx_ovl_di + 91);

    auto ts_yy_xxxyzz = pbuffer.data(idx_ovl_di + 92);

    auto ts_yy_xxxzzz = pbuffer.data(idx_ovl_di + 93);

    auto ts_yy_xxyyyy = pbuffer.data(idx_ovl_di + 94);

    auto ts_yy_xxyyyz = pbuffer.data(idx_ovl_di + 95);

    auto ts_yy_xxyyzz = pbuffer.data(idx_ovl_di + 96);

    auto ts_yy_xxyzzz = pbuffer.data(idx_ovl_di + 97);

    auto ts_yy_xxzzzz = pbuffer.data(idx_ovl_di + 98);

    auto ts_yy_xyyyyy = pbuffer.data(idx_ovl_di + 99);

    auto ts_yy_xyyyyz = pbuffer.data(idx_ovl_di + 100);

    auto ts_yy_xyyyzz = pbuffer.data(idx_ovl_di + 101);

    auto ts_yy_xyyzzz = pbuffer.data(idx_ovl_di + 102);

    auto ts_yy_xyzzzz = pbuffer.data(idx_ovl_di + 103);

    auto ts_yy_xzzzzz = pbuffer.data(idx_ovl_di + 104);

    auto ts_yy_yyyyyy = pbuffer.data(idx_ovl_di + 105);

    auto ts_yy_yyyyyz = pbuffer.data(idx_ovl_di + 106);

    auto ts_yy_yyyyzz = pbuffer.data(idx_ovl_di + 107);

    auto ts_yy_yyyzzz = pbuffer.data(idx_ovl_di + 108);

    auto ts_yy_yyzzzz = pbuffer.data(idx_ovl_di + 109);

    auto ts_yy_yzzzzz = pbuffer.data(idx_ovl_di + 110);

    auto ts_yy_zzzzzz = pbuffer.data(idx_ovl_di + 111);

    auto ts_yz_xxxxxx = pbuffer.data(idx_ovl_di + 112);

    auto ts_yz_xxxxxy = pbuffer.data(idx_ovl_di + 113);

    auto ts_yz_xxxxxz = pbuffer.data(idx_ovl_di + 114);

    auto ts_yz_xxxxyy = pbuffer.data(idx_ovl_di + 115);

    auto ts_yz_xxxxyz = pbuffer.data(idx_ovl_di + 116);

    auto ts_yz_xxxxzz = pbuffer.data(idx_ovl_di + 117);

    auto ts_yz_xxxyyy = pbuffer.data(idx_ovl_di + 118);

    auto ts_yz_xxxyyz = pbuffer.data(idx_ovl_di + 119);

    auto ts_yz_xxxyzz = pbuffer.data(idx_ovl_di + 120);

    auto ts_yz_xxxzzz = pbuffer.data(idx_ovl_di + 121);

    auto ts_yz_xxyyyy = pbuffer.data(idx_ovl_di + 122);

    auto ts_yz_xxyyyz = pbuffer.data(idx_ovl_di + 123);

    auto ts_yz_xxyyzz = pbuffer.data(idx_ovl_di + 124);

    auto ts_yz_xxyzzz = pbuffer.data(idx_ovl_di + 125);

    auto ts_yz_xxzzzz = pbuffer.data(idx_ovl_di + 126);

    auto ts_yz_xyyyyy = pbuffer.data(idx_ovl_di + 127);

    auto ts_yz_xyyyyz = pbuffer.data(idx_ovl_di + 128);

    auto ts_yz_xyyyzz = pbuffer.data(idx_ovl_di + 129);

    auto ts_yz_xyyzzz = pbuffer.data(idx_ovl_di + 130);

    auto ts_yz_xyzzzz = pbuffer.data(idx_ovl_di + 131);

    auto ts_yz_xzzzzz = pbuffer.data(idx_ovl_di + 132);

    auto ts_yz_yyyyyy = pbuffer.data(idx_ovl_di + 133);

    auto ts_yz_yyyyyz = pbuffer.data(idx_ovl_di + 134);

    auto ts_yz_yyyyzz = pbuffer.data(idx_ovl_di + 135);

    auto ts_yz_yyyzzz = pbuffer.data(idx_ovl_di + 136);

    auto ts_yz_yyzzzz = pbuffer.data(idx_ovl_di + 137);

    auto ts_yz_yzzzzz = pbuffer.data(idx_ovl_di + 138);

    auto ts_yz_zzzzzz = pbuffer.data(idx_ovl_di + 139);

    auto ts_zz_xxxxxx = pbuffer.data(idx_ovl_di + 140);

    auto ts_zz_xxxxxy = pbuffer.data(idx_ovl_di + 141);

    auto ts_zz_xxxxxz = pbuffer.data(idx_ovl_di + 142);

    auto ts_zz_xxxxyy = pbuffer.data(idx_ovl_di + 143);

    auto ts_zz_xxxxyz = pbuffer.data(idx_ovl_di + 144);

    auto ts_zz_xxxxzz = pbuffer.data(idx_ovl_di + 145);

    auto ts_zz_xxxyyy = pbuffer.data(idx_ovl_di + 146);

    auto ts_zz_xxxyyz = pbuffer.data(idx_ovl_di + 147);

    auto ts_zz_xxxyzz = pbuffer.data(idx_ovl_di + 148);

    auto ts_zz_xxxzzz = pbuffer.data(idx_ovl_di + 149);

    auto ts_zz_xxyyyy = pbuffer.data(idx_ovl_di + 150);

    auto ts_zz_xxyyyz = pbuffer.data(idx_ovl_di + 151);

    auto ts_zz_xxyyzz = pbuffer.data(idx_ovl_di + 152);

    auto ts_zz_xxyzzz = pbuffer.data(idx_ovl_di + 153);

    auto ts_zz_xxzzzz = pbuffer.data(idx_ovl_di + 154);

    auto ts_zz_xyyyyy = pbuffer.data(idx_ovl_di + 155);

    auto ts_zz_xyyyyz = pbuffer.data(idx_ovl_di + 156);

    auto ts_zz_xyyyzz = pbuffer.data(idx_ovl_di + 157);

    auto ts_zz_xyyzzz = pbuffer.data(idx_ovl_di + 158);

    auto ts_zz_xyzzzz = pbuffer.data(idx_ovl_di + 159);

    auto ts_zz_xzzzzz = pbuffer.data(idx_ovl_di + 160);

    auto ts_zz_yyyyyy = pbuffer.data(idx_ovl_di + 161);

    auto ts_zz_yyyyyz = pbuffer.data(idx_ovl_di + 162);

    auto ts_zz_yyyyzz = pbuffer.data(idx_ovl_di + 163);

    auto ts_zz_yyyzzz = pbuffer.data(idx_ovl_di + 164);

    auto ts_zz_yyzzzz = pbuffer.data(idx_ovl_di + 165);

    auto ts_zz_yzzzzz = pbuffer.data(idx_ovl_di + 166);

    auto ts_zz_zzzzzz = pbuffer.data(idx_ovl_di + 167);

    // Set up 0-28 components of targeted buffer : DI

    auto tk_xx_xxxxxx = pbuffer.data(idx_kin_di);

    auto tk_xx_xxxxxy = pbuffer.data(idx_kin_di + 1);

    auto tk_xx_xxxxxz = pbuffer.data(idx_kin_di + 2);

    auto tk_xx_xxxxyy = pbuffer.data(idx_kin_di + 3);

    auto tk_xx_xxxxyz = pbuffer.data(idx_kin_di + 4);

    auto tk_xx_xxxxzz = pbuffer.data(idx_kin_di + 5);

    auto tk_xx_xxxyyy = pbuffer.data(idx_kin_di + 6);

    auto tk_xx_xxxyyz = pbuffer.data(idx_kin_di + 7);

    auto tk_xx_xxxyzz = pbuffer.data(idx_kin_di + 8);

    auto tk_xx_xxxzzz = pbuffer.data(idx_kin_di + 9);

    auto tk_xx_xxyyyy = pbuffer.data(idx_kin_di + 10);

    auto tk_xx_xxyyyz = pbuffer.data(idx_kin_di + 11);

    auto tk_xx_xxyyzz = pbuffer.data(idx_kin_di + 12);

    auto tk_xx_xxyzzz = pbuffer.data(idx_kin_di + 13);

    auto tk_xx_xxzzzz = pbuffer.data(idx_kin_di + 14);

    auto tk_xx_xyyyyy = pbuffer.data(idx_kin_di + 15);

    auto tk_xx_xyyyyz = pbuffer.data(idx_kin_di + 16);

    auto tk_xx_xyyyzz = pbuffer.data(idx_kin_di + 17);

    auto tk_xx_xyyzzz = pbuffer.data(idx_kin_di + 18);

    auto tk_xx_xyzzzz = pbuffer.data(idx_kin_di + 19);

    auto tk_xx_xzzzzz = pbuffer.data(idx_kin_di + 20);

    auto tk_xx_yyyyyy = pbuffer.data(idx_kin_di + 21);

    auto tk_xx_yyyyyz = pbuffer.data(idx_kin_di + 22);

    auto tk_xx_yyyyzz = pbuffer.data(idx_kin_di + 23);

    auto tk_xx_yyyzzz = pbuffer.data(idx_kin_di + 24);

    auto tk_xx_yyzzzz = pbuffer.data(idx_kin_di + 25);

    auto tk_xx_yzzzzz = pbuffer.data(idx_kin_di + 26);

    auto tk_xx_zzzzzz = pbuffer.data(idx_kin_di + 27);

#pragma omp simd aligned(pa_x,             \
                             tk_0_xxxxxx,  \
                             tk_0_xxxxxy,  \
                             tk_0_xxxxxz,  \
                             tk_0_xxxxyy,  \
                             tk_0_xxxxyz,  \
                             tk_0_xxxxzz,  \
                             tk_0_xxxyyy,  \
                             tk_0_xxxyyz,  \
                             tk_0_xxxyzz,  \
                             tk_0_xxxzzz,  \
                             tk_0_xxyyyy,  \
                             tk_0_xxyyyz,  \
                             tk_0_xxyyzz,  \
                             tk_0_xxyzzz,  \
                             tk_0_xxzzzz,  \
                             tk_0_xyyyyy,  \
                             tk_0_xyyyyz,  \
                             tk_0_xyyyzz,  \
                             tk_0_xyyzzz,  \
                             tk_0_xyzzzz,  \
                             tk_0_xzzzzz,  \
                             tk_0_yyyyyy,  \
                             tk_0_yyyyyz,  \
                             tk_0_yyyyzz,  \
                             tk_0_yyyzzz,  \
                             tk_0_yyzzzz,  \
                             tk_0_yzzzzz,  \
                             tk_0_zzzzzz,  \
                             tk_x_xxxxx,   \
                             tk_x_xxxxxx,  \
                             tk_x_xxxxxy,  \
                             tk_x_xxxxxz,  \
                             tk_x_xxxxy,   \
                             tk_x_xxxxyy,  \
                             tk_x_xxxxyz,  \
                             tk_x_xxxxz,   \
                             tk_x_xxxxzz,  \
                             tk_x_xxxyy,   \
                             tk_x_xxxyyy,  \
                             tk_x_xxxyyz,  \
                             tk_x_xxxyz,   \
                             tk_x_xxxyzz,  \
                             tk_x_xxxzz,   \
                             tk_x_xxxzzz,  \
                             tk_x_xxyyy,   \
                             tk_x_xxyyyy,  \
                             tk_x_xxyyyz,  \
                             tk_x_xxyyz,   \
                             tk_x_xxyyzz,  \
                             tk_x_xxyzz,   \
                             tk_x_xxyzzz,  \
                             tk_x_xxzzz,   \
                             tk_x_xxzzzz,  \
                             tk_x_xyyyy,   \
                             tk_x_xyyyyy,  \
                             tk_x_xyyyyz,  \
                             tk_x_xyyyz,   \
                             tk_x_xyyyzz,  \
                             tk_x_xyyzz,   \
                             tk_x_xyyzzz,  \
                             tk_x_xyzzz,   \
                             tk_x_xyzzzz,  \
                             tk_x_xzzzz,   \
                             tk_x_xzzzzz,  \
                             tk_x_yyyyy,   \
                             tk_x_yyyyyy,  \
                             tk_x_yyyyyz,  \
                             tk_x_yyyyz,   \
                             tk_x_yyyyzz,  \
                             tk_x_yyyzz,   \
                             tk_x_yyyzzz,  \
                             tk_x_yyzzz,   \
                             tk_x_yyzzzz,  \
                             tk_x_yzzzz,   \
                             tk_x_yzzzzz,  \
                             tk_x_zzzzz,   \
                             tk_x_zzzzzz,  \
                             tk_xx_xxxxxx, \
                             tk_xx_xxxxxy, \
                             tk_xx_xxxxxz, \
                             tk_xx_xxxxyy, \
                             tk_xx_xxxxyz, \
                             tk_xx_xxxxzz, \
                             tk_xx_xxxyyy, \
                             tk_xx_xxxyyz, \
                             tk_xx_xxxyzz, \
                             tk_xx_xxxzzz, \
                             tk_xx_xxyyyy, \
                             tk_xx_xxyyyz, \
                             tk_xx_xxyyzz, \
                             tk_xx_xxyzzz, \
                             tk_xx_xxzzzz, \
                             tk_xx_xyyyyy, \
                             tk_xx_xyyyyz, \
                             tk_xx_xyyyzz, \
                             tk_xx_xyyzzz, \
                             tk_xx_xyzzzz, \
                             tk_xx_xzzzzz, \
                             tk_xx_yyyyyy, \
                             tk_xx_yyyyyz, \
                             tk_xx_yyyyzz, \
                             tk_xx_yyyzzz, \
                             tk_xx_yyzzzz, \
                             tk_xx_yzzzzz, \
                             tk_xx_zzzzzz, \
                             ts_0_xxxxxx,  \
                             ts_0_xxxxxy,  \
                             ts_0_xxxxxz,  \
                             ts_0_xxxxyy,  \
                             ts_0_xxxxyz,  \
                             ts_0_xxxxzz,  \
                             ts_0_xxxyyy,  \
                             ts_0_xxxyyz,  \
                             ts_0_xxxyzz,  \
                             ts_0_xxxzzz,  \
                             ts_0_xxyyyy,  \
                             ts_0_xxyyyz,  \
                             ts_0_xxyyzz,  \
                             ts_0_xxyzzz,  \
                             ts_0_xxzzzz,  \
                             ts_0_xyyyyy,  \
                             ts_0_xyyyyz,  \
                             ts_0_xyyyzz,  \
                             ts_0_xyyzzz,  \
                             ts_0_xyzzzz,  \
                             ts_0_xzzzzz,  \
                             ts_0_yyyyyy,  \
                             ts_0_yyyyyz,  \
                             ts_0_yyyyzz,  \
                             ts_0_yyyzzz,  \
                             ts_0_yyzzzz,  \
                             ts_0_yzzzzz,  \
                             ts_0_zzzzzz,  \
                             ts_xx_xxxxxx, \
                             ts_xx_xxxxxy, \
                             ts_xx_xxxxxz, \
                             ts_xx_xxxxyy, \
                             ts_xx_xxxxyz, \
                             ts_xx_xxxxzz, \
                             ts_xx_xxxyyy, \
                             ts_xx_xxxyyz, \
                             ts_xx_xxxyzz, \
                             ts_xx_xxxzzz, \
                             ts_xx_xxyyyy, \
                             ts_xx_xxyyyz, \
                             ts_xx_xxyyzz, \
                             ts_xx_xxyzzz, \
                             ts_xx_xxzzzz, \
                             ts_xx_xyyyyy, \
                             ts_xx_xyyyyz, \
                             ts_xx_xyyyzz, \
                             ts_xx_xyyzzz, \
                             ts_xx_xyzzzz, \
                             ts_xx_xzzzzz, \
                             ts_xx_yyyyyy, \
                             ts_xx_yyyyyz, \
                             ts_xx_yyyyzz, \
                             ts_xx_yyyzzz, \
                             ts_xx_yyzzzz, \
                             ts_xx_yzzzzz, \
                             ts_xx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_xxxxxx[i] = -2.0 * ts_0_xxxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxxx[i] * fe_0 + 6.0 * tk_x_xxxxx[i] * fe_0 + tk_x_xxxxxx[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxxx[i] * fz_0;

        tk_xx_xxxxxy[i] = -2.0 * ts_0_xxxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxxy[i] * fe_0 + 5.0 * tk_x_xxxxy[i] * fe_0 + tk_x_xxxxxy[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxxy[i] * fz_0;

        tk_xx_xxxxxz[i] = -2.0 * ts_0_xxxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxxz[i] * fe_0 + 5.0 * tk_x_xxxxz[i] * fe_0 + tk_x_xxxxxz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxxz[i] * fz_0;

        tk_xx_xxxxyy[i] = -2.0 * ts_0_xxxxyy[i] * fbe_0 * fz_0 + tk_0_xxxxyy[i] * fe_0 + 4.0 * tk_x_xxxyy[i] * fe_0 + tk_x_xxxxyy[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxyy[i] * fz_0;

        tk_xx_xxxxyz[i] = -2.0 * ts_0_xxxxyz[i] * fbe_0 * fz_0 + tk_0_xxxxyz[i] * fe_0 + 4.0 * tk_x_xxxyz[i] * fe_0 + tk_x_xxxxyz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxyz[i] * fz_0;

        tk_xx_xxxxzz[i] = -2.0 * ts_0_xxxxzz[i] * fbe_0 * fz_0 + tk_0_xxxxzz[i] * fe_0 + 4.0 * tk_x_xxxzz[i] * fe_0 + tk_x_xxxxzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxxzz[i] * fz_0;

        tk_xx_xxxyyy[i] = -2.0 * ts_0_xxxyyy[i] * fbe_0 * fz_0 + tk_0_xxxyyy[i] * fe_0 + 3.0 * tk_x_xxyyy[i] * fe_0 + tk_x_xxxyyy[i] * pa_x[i] +
                          2.0 * ts_xx_xxxyyy[i] * fz_0;

        tk_xx_xxxyyz[i] = -2.0 * ts_0_xxxyyz[i] * fbe_0 * fz_0 + tk_0_xxxyyz[i] * fe_0 + 3.0 * tk_x_xxyyz[i] * fe_0 + tk_x_xxxyyz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxyyz[i] * fz_0;

        tk_xx_xxxyzz[i] = -2.0 * ts_0_xxxyzz[i] * fbe_0 * fz_0 + tk_0_xxxyzz[i] * fe_0 + 3.0 * tk_x_xxyzz[i] * fe_0 + tk_x_xxxyzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxyzz[i] * fz_0;

        tk_xx_xxxzzz[i] = -2.0 * ts_0_xxxzzz[i] * fbe_0 * fz_0 + tk_0_xxxzzz[i] * fe_0 + 3.0 * tk_x_xxzzz[i] * fe_0 + tk_x_xxxzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxxzzz[i] * fz_0;

        tk_xx_xxyyyy[i] = -2.0 * ts_0_xxyyyy[i] * fbe_0 * fz_0 + tk_0_xxyyyy[i] * fe_0 + 2.0 * tk_x_xyyyy[i] * fe_0 + tk_x_xxyyyy[i] * pa_x[i] +
                          2.0 * ts_xx_xxyyyy[i] * fz_0;

        tk_xx_xxyyyz[i] = -2.0 * ts_0_xxyyyz[i] * fbe_0 * fz_0 + tk_0_xxyyyz[i] * fe_0 + 2.0 * tk_x_xyyyz[i] * fe_0 + tk_x_xxyyyz[i] * pa_x[i] +
                          2.0 * ts_xx_xxyyyz[i] * fz_0;

        tk_xx_xxyyzz[i] = -2.0 * ts_0_xxyyzz[i] * fbe_0 * fz_0 + tk_0_xxyyzz[i] * fe_0 + 2.0 * tk_x_xyyzz[i] * fe_0 + tk_x_xxyyzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxyyzz[i] * fz_0;

        tk_xx_xxyzzz[i] = -2.0 * ts_0_xxyzzz[i] * fbe_0 * fz_0 + tk_0_xxyzzz[i] * fe_0 + 2.0 * tk_x_xyzzz[i] * fe_0 + tk_x_xxyzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxyzzz[i] * fz_0;

        tk_xx_xxzzzz[i] = -2.0 * ts_0_xxzzzz[i] * fbe_0 * fz_0 + tk_0_xxzzzz[i] * fe_0 + 2.0 * tk_x_xzzzz[i] * fe_0 + tk_x_xxzzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xxzzzz[i] * fz_0;

        tk_xx_xyyyyy[i] = -2.0 * ts_0_xyyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyyy[i] * fe_0 + tk_x_yyyyy[i] * fe_0 + tk_x_xyyyyy[i] * pa_x[i] +
                          2.0 * ts_xx_xyyyyy[i] * fz_0;

        tk_xx_xyyyyz[i] = -2.0 * ts_0_xyyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyyz[i] * fe_0 + tk_x_yyyyz[i] * fe_0 + tk_x_xyyyyz[i] * pa_x[i] +
                          2.0 * ts_xx_xyyyyz[i] * fz_0;

        tk_xx_xyyyzz[i] = -2.0 * ts_0_xyyyzz[i] * fbe_0 * fz_0 + tk_0_xyyyzz[i] * fe_0 + tk_x_yyyzz[i] * fe_0 + tk_x_xyyyzz[i] * pa_x[i] +
                          2.0 * ts_xx_xyyyzz[i] * fz_0;

        tk_xx_xyyzzz[i] = -2.0 * ts_0_xyyzzz[i] * fbe_0 * fz_0 + tk_0_xyyzzz[i] * fe_0 + tk_x_yyzzz[i] * fe_0 + tk_x_xyyzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xyyzzz[i] * fz_0;

        tk_xx_xyzzzz[i] = -2.0 * ts_0_xyzzzz[i] * fbe_0 * fz_0 + tk_0_xyzzzz[i] * fe_0 + tk_x_yzzzz[i] * fe_0 + tk_x_xyzzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xyzzzz[i] * fz_0;

        tk_xx_xzzzzz[i] = -2.0 * ts_0_xzzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzzz[i] * fe_0 + tk_x_zzzzz[i] * fe_0 + tk_x_xzzzzz[i] * pa_x[i] +
                          2.0 * ts_xx_xzzzzz[i] * fz_0;

        tk_xx_yyyyyy[i] = -2.0 * ts_0_yyyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyyy[i] * fe_0 + tk_x_yyyyyy[i] * pa_x[i] + 2.0 * ts_xx_yyyyyy[i] * fz_0;

        tk_xx_yyyyyz[i] = -2.0 * ts_0_yyyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyyz[i] * fe_0 + tk_x_yyyyyz[i] * pa_x[i] + 2.0 * ts_xx_yyyyyz[i] * fz_0;

        tk_xx_yyyyzz[i] = -2.0 * ts_0_yyyyzz[i] * fbe_0 * fz_0 + tk_0_yyyyzz[i] * fe_0 + tk_x_yyyyzz[i] * pa_x[i] + 2.0 * ts_xx_yyyyzz[i] * fz_0;

        tk_xx_yyyzzz[i] = -2.0 * ts_0_yyyzzz[i] * fbe_0 * fz_0 + tk_0_yyyzzz[i] * fe_0 + tk_x_yyyzzz[i] * pa_x[i] + 2.0 * ts_xx_yyyzzz[i] * fz_0;

        tk_xx_yyzzzz[i] = -2.0 * ts_0_yyzzzz[i] * fbe_0 * fz_0 + tk_0_yyzzzz[i] * fe_0 + tk_x_yyzzzz[i] * pa_x[i] + 2.0 * ts_xx_yyzzzz[i] * fz_0;

        tk_xx_yzzzzz[i] = -2.0 * ts_0_yzzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzzz[i] * fe_0 + tk_x_yzzzzz[i] * pa_x[i] + 2.0 * ts_xx_yzzzzz[i] * fz_0;

        tk_xx_zzzzzz[i] = -2.0 * ts_0_zzzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzzz[i] * fe_0 + tk_x_zzzzzz[i] * pa_x[i] + 2.0 * ts_xx_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : DI

    auto tk_xy_xxxxxx = pbuffer.data(idx_kin_di + 28);

    auto tk_xy_xxxxxy = pbuffer.data(idx_kin_di + 29);

    auto tk_xy_xxxxxz = pbuffer.data(idx_kin_di + 30);

    auto tk_xy_xxxxyy = pbuffer.data(idx_kin_di + 31);

    auto tk_xy_xxxxyz = pbuffer.data(idx_kin_di + 32);

    auto tk_xy_xxxxzz = pbuffer.data(idx_kin_di + 33);

    auto tk_xy_xxxyyy = pbuffer.data(idx_kin_di + 34);

    auto tk_xy_xxxyyz = pbuffer.data(idx_kin_di + 35);

    auto tk_xy_xxxyzz = pbuffer.data(idx_kin_di + 36);

    auto tk_xy_xxxzzz = pbuffer.data(idx_kin_di + 37);

    auto tk_xy_xxyyyy = pbuffer.data(idx_kin_di + 38);

    auto tk_xy_xxyyyz = pbuffer.data(idx_kin_di + 39);

    auto tk_xy_xxyyzz = pbuffer.data(idx_kin_di + 40);

    auto tk_xy_xxyzzz = pbuffer.data(idx_kin_di + 41);

    auto tk_xy_xxzzzz = pbuffer.data(idx_kin_di + 42);

    auto tk_xy_xyyyyy = pbuffer.data(idx_kin_di + 43);

    auto tk_xy_xyyyyz = pbuffer.data(idx_kin_di + 44);

    auto tk_xy_xyyyzz = pbuffer.data(idx_kin_di + 45);

    auto tk_xy_xyyzzz = pbuffer.data(idx_kin_di + 46);

    auto tk_xy_xyzzzz = pbuffer.data(idx_kin_di + 47);

    auto tk_xy_xzzzzz = pbuffer.data(idx_kin_di + 48);

    auto tk_xy_yyyyyy = pbuffer.data(idx_kin_di + 49);

    auto tk_xy_yyyyyz = pbuffer.data(idx_kin_di + 50);

    auto tk_xy_yyyyzz = pbuffer.data(idx_kin_di + 51);

    auto tk_xy_yyyzzz = pbuffer.data(idx_kin_di + 52);

    auto tk_xy_yyzzzz = pbuffer.data(idx_kin_di + 53);

    auto tk_xy_yzzzzz = pbuffer.data(idx_kin_di + 54);

    auto tk_xy_zzzzzz = pbuffer.data(idx_kin_di + 55);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_x_xxxxxx,  \
                             tk_x_xxxxxz,  \
                             tk_x_xxxxzz,  \
                             tk_x_xxxzzz,  \
                             tk_x_xxzzzz,  \
                             tk_x_xzzzzz,  \
                             tk_xy_xxxxxx, \
                             tk_xy_xxxxxy, \
                             tk_xy_xxxxxz, \
                             tk_xy_xxxxyy, \
                             tk_xy_xxxxyz, \
                             tk_xy_xxxxzz, \
                             tk_xy_xxxyyy, \
                             tk_xy_xxxyyz, \
                             tk_xy_xxxyzz, \
                             tk_xy_xxxzzz, \
                             tk_xy_xxyyyy, \
                             tk_xy_xxyyyz, \
                             tk_xy_xxyyzz, \
                             tk_xy_xxyzzz, \
                             tk_xy_xxzzzz, \
                             tk_xy_xyyyyy, \
                             tk_xy_xyyyyz, \
                             tk_xy_xyyyzz, \
                             tk_xy_xyyzzz, \
                             tk_xy_xyzzzz, \
                             tk_xy_xzzzzz, \
                             tk_xy_yyyyyy, \
                             tk_xy_yyyyyz, \
                             tk_xy_yyyyzz, \
                             tk_xy_yyyzzz, \
                             tk_xy_yyzzzz, \
                             tk_xy_yzzzzz, \
                             tk_xy_zzzzzz, \
                             tk_y_xxxxxy,  \
                             tk_y_xxxxy,   \
                             tk_y_xxxxyy,  \
                             tk_y_xxxxyz,  \
                             tk_y_xxxyy,   \
                             tk_y_xxxyyy,  \
                             tk_y_xxxyyz,  \
                             tk_y_xxxyz,   \
                             tk_y_xxxyzz,  \
                             tk_y_xxyyy,   \
                             tk_y_xxyyyy,  \
                             tk_y_xxyyyz,  \
                             tk_y_xxyyz,   \
                             tk_y_xxyyzz,  \
                             tk_y_xxyzz,   \
                             tk_y_xxyzzz,  \
                             tk_y_xyyyy,   \
                             tk_y_xyyyyy,  \
                             tk_y_xyyyyz,  \
                             tk_y_xyyyz,   \
                             tk_y_xyyyzz,  \
                             tk_y_xyyzz,   \
                             tk_y_xyyzzz,  \
                             tk_y_xyzzz,   \
                             tk_y_xyzzzz,  \
                             tk_y_yyyyy,   \
                             tk_y_yyyyyy,  \
                             tk_y_yyyyyz,  \
                             tk_y_yyyyz,   \
                             tk_y_yyyyzz,  \
                             tk_y_yyyzz,   \
                             tk_y_yyyzzz,  \
                             tk_y_yyzzz,   \
                             tk_y_yyzzzz,  \
                             tk_y_yzzzz,   \
                             tk_y_yzzzzz,  \
                             tk_y_zzzzzz,  \
                             ts_xy_xxxxxx, \
                             ts_xy_xxxxxy, \
                             ts_xy_xxxxxz, \
                             ts_xy_xxxxyy, \
                             ts_xy_xxxxyz, \
                             ts_xy_xxxxzz, \
                             ts_xy_xxxyyy, \
                             ts_xy_xxxyyz, \
                             ts_xy_xxxyzz, \
                             ts_xy_xxxzzz, \
                             ts_xy_xxyyyy, \
                             ts_xy_xxyyyz, \
                             ts_xy_xxyyzz, \
                             ts_xy_xxyzzz, \
                             ts_xy_xxzzzz, \
                             ts_xy_xyyyyy, \
                             ts_xy_xyyyyz, \
                             ts_xy_xyyyzz, \
                             ts_xy_xyyzzz, \
                             ts_xy_xyzzzz, \
                             ts_xy_xzzzzz, \
                             ts_xy_yyyyyy, \
                             ts_xy_yyyyyz, \
                             ts_xy_yyyyzz, \
                             ts_xy_yyyzzz, \
                             ts_xy_yyzzzz, \
                             ts_xy_yzzzzz, \
                             ts_xy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xy_xxxxxx[i] = tk_x_xxxxxx[i] * pa_y[i] + 2.0 * ts_xy_xxxxxx[i] * fz_0;

        tk_xy_xxxxxy[i] = 5.0 * tk_y_xxxxy[i] * fe_0 + tk_y_xxxxxy[i] * pa_x[i] + 2.0 * ts_xy_xxxxxy[i] * fz_0;

        tk_xy_xxxxxz[i] = tk_x_xxxxxz[i] * pa_y[i] + 2.0 * ts_xy_xxxxxz[i] * fz_0;

        tk_xy_xxxxyy[i] = 4.0 * tk_y_xxxyy[i] * fe_0 + tk_y_xxxxyy[i] * pa_x[i] + 2.0 * ts_xy_xxxxyy[i] * fz_0;

        tk_xy_xxxxyz[i] = 4.0 * tk_y_xxxyz[i] * fe_0 + tk_y_xxxxyz[i] * pa_x[i] + 2.0 * ts_xy_xxxxyz[i] * fz_0;

        tk_xy_xxxxzz[i] = tk_x_xxxxzz[i] * pa_y[i] + 2.0 * ts_xy_xxxxzz[i] * fz_0;

        tk_xy_xxxyyy[i] = 3.0 * tk_y_xxyyy[i] * fe_0 + tk_y_xxxyyy[i] * pa_x[i] + 2.0 * ts_xy_xxxyyy[i] * fz_0;

        tk_xy_xxxyyz[i] = 3.0 * tk_y_xxyyz[i] * fe_0 + tk_y_xxxyyz[i] * pa_x[i] + 2.0 * ts_xy_xxxyyz[i] * fz_0;

        tk_xy_xxxyzz[i] = 3.0 * tk_y_xxyzz[i] * fe_0 + tk_y_xxxyzz[i] * pa_x[i] + 2.0 * ts_xy_xxxyzz[i] * fz_0;

        tk_xy_xxxzzz[i] = tk_x_xxxzzz[i] * pa_y[i] + 2.0 * ts_xy_xxxzzz[i] * fz_0;

        tk_xy_xxyyyy[i] = 2.0 * tk_y_xyyyy[i] * fe_0 + tk_y_xxyyyy[i] * pa_x[i] + 2.0 * ts_xy_xxyyyy[i] * fz_0;

        tk_xy_xxyyyz[i] = 2.0 * tk_y_xyyyz[i] * fe_0 + tk_y_xxyyyz[i] * pa_x[i] + 2.0 * ts_xy_xxyyyz[i] * fz_0;

        tk_xy_xxyyzz[i] = 2.0 * tk_y_xyyzz[i] * fe_0 + tk_y_xxyyzz[i] * pa_x[i] + 2.0 * ts_xy_xxyyzz[i] * fz_0;

        tk_xy_xxyzzz[i] = 2.0 * tk_y_xyzzz[i] * fe_0 + tk_y_xxyzzz[i] * pa_x[i] + 2.0 * ts_xy_xxyzzz[i] * fz_0;

        tk_xy_xxzzzz[i] = tk_x_xxzzzz[i] * pa_y[i] + 2.0 * ts_xy_xxzzzz[i] * fz_0;

        tk_xy_xyyyyy[i] = tk_y_yyyyy[i] * fe_0 + tk_y_xyyyyy[i] * pa_x[i] + 2.0 * ts_xy_xyyyyy[i] * fz_0;

        tk_xy_xyyyyz[i] = tk_y_yyyyz[i] * fe_0 + tk_y_xyyyyz[i] * pa_x[i] + 2.0 * ts_xy_xyyyyz[i] * fz_0;

        tk_xy_xyyyzz[i] = tk_y_yyyzz[i] * fe_0 + tk_y_xyyyzz[i] * pa_x[i] + 2.0 * ts_xy_xyyyzz[i] * fz_0;

        tk_xy_xyyzzz[i] = tk_y_yyzzz[i] * fe_0 + tk_y_xyyzzz[i] * pa_x[i] + 2.0 * ts_xy_xyyzzz[i] * fz_0;

        tk_xy_xyzzzz[i] = tk_y_yzzzz[i] * fe_0 + tk_y_xyzzzz[i] * pa_x[i] + 2.0 * ts_xy_xyzzzz[i] * fz_0;

        tk_xy_xzzzzz[i] = tk_x_xzzzzz[i] * pa_y[i] + 2.0 * ts_xy_xzzzzz[i] * fz_0;

        tk_xy_yyyyyy[i] = tk_y_yyyyyy[i] * pa_x[i] + 2.0 * ts_xy_yyyyyy[i] * fz_0;

        tk_xy_yyyyyz[i] = tk_y_yyyyyz[i] * pa_x[i] + 2.0 * ts_xy_yyyyyz[i] * fz_0;

        tk_xy_yyyyzz[i] = tk_y_yyyyzz[i] * pa_x[i] + 2.0 * ts_xy_yyyyzz[i] * fz_0;

        tk_xy_yyyzzz[i] = tk_y_yyyzzz[i] * pa_x[i] + 2.0 * ts_xy_yyyzzz[i] * fz_0;

        tk_xy_yyzzzz[i] = tk_y_yyzzzz[i] * pa_x[i] + 2.0 * ts_xy_yyzzzz[i] * fz_0;

        tk_xy_yzzzzz[i] = tk_y_yzzzzz[i] * pa_x[i] + 2.0 * ts_xy_yzzzzz[i] * fz_0;

        tk_xy_zzzzzz[i] = tk_y_zzzzzz[i] * pa_x[i] + 2.0 * ts_xy_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : DI

    auto tk_xz_xxxxxx = pbuffer.data(idx_kin_di + 56);

    auto tk_xz_xxxxxy = pbuffer.data(idx_kin_di + 57);

    auto tk_xz_xxxxxz = pbuffer.data(idx_kin_di + 58);

    auto tk_xz_xxxxyy = pbuffer.data(idx_kin_di + 59);

    auto tk_xz_xxxxyz = pbuffer.data(idx_kin_di + 60);

    auto tk_xz_xxxxzz = pbuffer.data(idx_kin_di + 61);

    auto tk_xz_xxxyyy = pbuffer.data(idx_kin_di + 62);

    auto tk_xz_xxxyyz = pbuffer.data(idx_kin_di + 63);

    auto tk_xz_xxxyzz = pbuffer.data(idx_kin_di + 64);

    auto tk_xz_xxxzzz = pbuffer.data(idx_kin_di + 65);

    auto tk_xz_xxyyyy = pbuffer.data(idx_kin_di + 66);

    auto tk_xz_xxyyyz = pbuffer.data(idx_kin_di + 67);

    auto tk_xz_xxyyzz = pbuffer.data(idx_kin_di + 68);

    auto tk_xz_xxyzzz = pbuffer.data(idx_kin_di + 69);

    auto tk_xz_xxzzzz = pbuffer.data(idx_kin_di + 70);

    auto tk_xz_xyyyyy = pbuffer.data(idx_kin_di + 71);

    auto tk_xz_xyyyyz = pbuffer.data(idx_kin_di + 72);

    auto tk_xz_xyyyzz = pbuffer.data(idx_kin_di + 73);

    auto tk_xz_xyyzzz = pbuffer.data(idx_kin_di + 74);

    auto tk_xz_xyzzzz = pbuffer.data(idx_kin_di + 75);

    auto tk_xz_xzzzzz = pbuffer.data(idx_kin_di + 76);

    auto tk_xz_yyyyyy = pbuffer.data(idx_kin_di + 77);

    auto tk_xz_yyyyyz = pbuffer.data(idx_kin_di + 78);

    auto tk_xz_yyyyzz = pbuffer.data(idx_kin_di + 79);

    auto tk_xz_yyyzzz = pbuffer.data(idx_kin_di + 80);

    auto tk_xz_yyzzzz = pbuffer.data(idx_kin_di + 81);

    auto tk_xz_yzzzzz = pbuffer.data(idx_kin_di + 82);

    auto tk_xz_zzzzzz = pbuffer.data(idx_kin_di + 83);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_x_xxxxxx,  \
                             tk_x_xxxxxy,  \
                             tk_x_xxxxyy,  \
                             tk_x_xxxyyy,  \
                             tk_x_xxyyyy,  \
                             tk_x_xyyyyy,  \
                             tk_xz_xxxxxx, \
                             tk_xz_xxxxxy, \
                             tk_xz_xxxxxz, \
                             tk_xz_xxxxyy, \
                             tk_xz_xxxxyz, \
                             tk_xz_xxxxzz, \
                             tk_xz_xxxyyy, \
                             tk_xz_xxxyyz, \
                             tk_xz_xxxyzz, \
                             tk_xz_xxxzzz, \
                             tk_xz_xxyyyy, \
                             tk_xz_xxyyyz, \
                             tk_xz_xxyyzz, \
                             tk_xz_xxyzzz, \
                             tk_xz_xxzzzz, \
                             tk_xz_xyyyyy, \
                             tk_xz_xyyyyz, \
                             tk_xz_xyyyzz, \
                             tk_xz_xyyzzz, \
                             tk_xz_xyzzzz, \
                             tk_xz_xzzzzz, \
                             tk_xz_yyyyyy, \
                             tk_xz_yyyyyz, \
                             tk_xz_yyyyzz, \
                             tk_xz_yyyzzz, \
                             tk_xz_yyzzzz, \
                             tk_xz_yzzzzz, \
                             tk_xz_zzzzzz, \
                             tk_z_xxxxxz,  \
                             tk_z_xxxxyz,  \
                             tk_z_xxxxz,   \
                             tk_z_xxxxzz,  \
                             tk_z_xxxyyz,  \
                             tk_z_xxxyz,   \
                             tk_z_xxxyzz,  \
                             tk_z_xxxzz,   \
                             tk_z_xxxzzz,  \
                             tk_z_xxyyyz,  \
                             tk_z_xxyyz,   \
                             tk_z_xxyyzz,  \
                             tk_z_xxyzz,   \
                             tk_z_xxyzzz,  \
                             tk_z_xxzzz,   \
                             tk_z_xxzzzz,  \
                             tk_z_xyyyyz,  \
                             tk_z_xyyyz,   \
                             tk_z_xyyyzz,  \
                             tk_z_xyyzz,   \
                             tk_z_xyyzzz,  \
                             tk_z_xyzzz,   \
                             tk_z_xyzzzz,  \
                             tk_z_xzzzz,   \
                             tk_z_xzzzzz,  \
                             tk_z_yyyyyy,  \
                             tk_z_yyyyyz,  \
                             tk_z_yyyyz,   \
                             tk_z_yyyyzz,  \
                             tk_z_yyyzz,   \
                             tk_z_yyyzzz,  \
                             tk_z_yyzzz,   \
                             tk_z_yyzzzz,  \
                             tk_z_yzzzz,   \
                             tk_z_yzzzzz,  \
                             tk_z_zzzzz,   \
                             tk_z_zzzzzz,  \
                             ts_xz_xxxxxx, \
                             ts_xz_xxxxxy, \
                             ts_xz_xxxxxz, \
                             ts_xz_xxxxyy, \
                             ts_xz_xxxxyz, \
                             ts_xz_xxxxzz, \
                             ts_xz_xxxyyy, \
                             ts_xz_xxxyyz, \
                             ts_xz_xxxyzz, \
                             ts_xz_xxxzzz, \
                             ts_xz_xxyyyy, \
                             ts_xz_xxyyyz, \
                             ts_xz_xxyyzz, \
                             ts_xz_xxyzzz, \
                             ts_xz_xxzzzz, \
                             ts_xz_xyyyyy, \
                             ts_xz_xyyyyz, \
                             ts_xz_xyyyzz, \
                             ts_xz_xyyzzz, \
                             ts_xz_xyzzzz, \
                             ts_xz_xzzzzz, \
                             ts_xz_yyyyyy, \
                             ts_xz_yyyyyz, \
                             ts_xz_yyyyzz, \
                             ts_xz_yyyzzz, \
                             ts_xz_yyzzzz, \
                             ts_xz_yzzzzz, \
                             ts_xz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xz_xxxxxx[i] = tk_x_xxxxxx[i] * pa_z[i] + 2.0 * ts_xz_xxxxxx[i] * fz_0;

        tk_xz_xxxxxy[i] = tk_x_xxxxxy[i] * pa_z[i] + 2.0 * ts_xz_xxxxxy[i] * fz_0;

        tk_xz_xxxxxz[i] = 5.0 * tk_z_xxxxz[i] * fe_0 + tk_z_xxxxxz[i] * pa_x[i] + 2.0 * ts_xz_xxxxxz[i] * fz_0;

        tk_xz_xxxxyy[i] = tk_x_xxxxyy[i] * pa_z[i] + 2.0 * ts_xz_xxxxyy[i] * fz_0;

        tk_xz_xxxxyz[i] = 4.0 * tk_z_xxxyz[i] * fe_0 + tk_z_xxxxyz[i] * pa_x[i] + 2.0 * ts_xz_xxxxyz[i] * fz_0;

        tk_xz_xxxxzz[i] = 4.0 * tk_z_xxxzz[i] * fe_0 + tk_z_xxxxzz[i] * pa_x[i] + 2.0 * ts_xz_xxxxzz[i] * fz_0;

        tk_xz_xxxyyy[i] = tk_x_xxxyyy[i] * pa_z[i] + 2.0 * ts_xz_xxxyyy[i] * fz_0;

        tk_xz_xxxyyz[i] = 3.0 * tk_z_xxyyz[i] * fe_0 + tk_z_xxxyyz[i] * pa_x[i] + 2.0 * ts_xz_xxxyyz[i] * fz_0;

        tk_xz_xxxyzz[i] = 3.0 * tk_z_xxyzz[i] * fe_0 + tk_z_xxxyzz[i] * pa_x[i] + 2.0 * ts_xz_xxxyzz[i] * fz_0;

        tk_xz_xxxzzz[i] = 3.0 * tk_z_xxzzz[i] * fe_0 + tk_z_xxxzzz[i] * pa_x[i] + 2.0 * ts_xz_xxxzzz[i] * fz_0;

        tk_xz_xxyyyy[i] = tk_x_xxyyyy[i] * pa_z[i] + 2.0 * ts_xz_xxyyyy[i] * fz_0;

        tk_xz_xxyyyz[i] = 2.0 * tk_z_xyyyz[i] * fe_0 + tk_z_xxyyyz[i] * pa_x[i] + 2.0 * ts_xz_xxyyyz[i] * fz_0;

        tk_xz_xxyyzz[i] = 2.0 * tk_z_xyyzz[i] * fe_0 + tk_z_xxyyzz[i] * pa_x[i] + 2.0 * ts_xz_xxyyzz[i] * fz_0;

        tk_xz_xxyzzz[i] = 2.0 * tk_z_xyzzz[i] * fe_0 + tk_z_xxyzzz[i] * pa_x[i] + 2.0 * ts_xz_xxyzzz[i] * fz_0;

        tk_xz_xxzzzz[i] = 2.0 * tk_z_xzzzz[i] * fe_0 + tk_z_xxzzzz[i] * pa_x[i] + 2.0 * ts_xz_xxzzzz[i] * fz_0;

        tk_xz_xyyyyy[i] = tk_x_xyyyyy[i] * pa_z[i] + 2.0 * ts_xz_xyyyyy[i] * fz_0;

        tk_xz_xyyyyz[i] = tk_z_yyyyz[i] * fe_0 + tk_z_xyyyyz[i] * pa_x[i] + 2.0 * ts_xz_xyyyyz[i] * fz_0;

        tk_xz_xyyyzz[i] = tk_z_yyyzz[i] * fe_0 + tk_z_xyyyzz[i] * pa_x[i] + 2.0 * ts_xz_xyyyzz[i] * fz_0;

        tk_xz_xyyzzz[i] = tk_z_yyzzz[i] * fe_0 + tk_z_xyyzzz[i] * pa_x[i] + 2.0 * ts_xz_xyyzzz[i] * fz_0;

        tk_xz_xyzzzz[i] = tk_z_yzzzz[i] * fe_0 + tk_z_xyzzzz[i] * pa_x[i] + 2.0 * ts_xz_xyzzzz[i] * fz_0;

        tk_xz_xzzzzz[i] = tk_z_zzzzz[i] * fe_0 + tk_z_xzzzzz[i] * pa_x[i] + 2.0 * ts_xz_xzzzzz[i] * fz_0;

        tk_xz_yyyyyy[i] = tk_z_yyyyyy[i] * pa_x[i] + 2.0 * ts_xz_yyyyyy[i] * fz_0;

        tk_xz_yyyyyz[i] = tk_z_yyyyyz[i] * pa_x[i] + 2.0 * ts_xz_yyyyyz[i] * fz_0;

        tk_xz_yyyyzz[i] = tk_z_yyyyzz[i] * pa_x[i] + 2.0 * ts_xz_yyyyzz[i] * fz_0;

        tk_xz_yyyzzz[i] = tk_z_yyyzzz[i] * pa_x[i] + 2.0 * ts_xz_yyyzzz[i] * fz_0;

        tk_xz_yyzzzz[i] = tk_z_yyzzzz[i] * pa_x[i] + 2.0 * ts_xz_yyzzzz[i] * fz_0;

        tk_xz_yzzzzz[i] = tk_z_yzzzzz[i] * pa_x[i] + 2.0 * ts_xz_yzzzzz[i] * fz_0;

        tk_xz_zzzzzz[i] = tk_z_zzzzzz[i] * pa_x[i] + 2.0 * ts_xz_zzzzzz[i] * fz_0;
    }

    // Set up 84-112 components of targeted buffer : DI

    auto tk_yy_xxxxxx = pbuffer.data(idx_kin_di + 84);

    auto tk_yy_xxxxxy = pbuffer.data(idx_kin_di + 85);

    auto tk_yy_xxxxxz = pbuffer.data(idx_kin_di + 86);

    auto tk_yy_xxxxyy = pbuffer.data(idx_kin_di + 87);

    auto tk_yy_xxxxyz = pbuffer.data(idx_kin_di + 88);

    auto tk_yy_xxxxzz = pbuffer.data(idx_kin_di + 89);

    auto tk_yy_xxxyyy = pbuffer.data(idx_kin_di + 90);

    auto tk_yy_xxxyyz = pbuffer.data(idx_kin_di + 91);

    auto tk_yy_xxxyzz = pbuffer.data(idx_kin_di + 92);

    auto tk_yy_xxxzzz = pbuffer.data(idx_kin_di + 93);

    auto tk_yy_xxyyyy = pbuffer.data(idx_kin_di + 94);

    auto tk_yy_xxyyyz = pbuffer.data(idx_kin_di + 95);

    auto tk_yy_xxyyzz = pbuffer.data(idx_kin_di + 96);

    auto tk_yy_xxyzzz = pbuffer.data(idx_kin_di + 97);

    auto tk_yy_xxzzzz = pbuffer.data(idx_kin_di + 98);

    auto tk_yy_xyyyyy = pbuffer.data(idx_kin_di + 99);

    auto tk_yy_xyyyyz = pbuffer.data(idx_kin_di + 100);

    auto tk_yy_xyyyzz = pbuffer.data(idx_kin_di + 101);

    auto tk_yy_xyyzzz = pbuffer.data(idx_kin_di + 102);

    auto tk_yy_xyzzzz = pbuffer.data(idx_kin_di + 103);

    auto tk_yy_xzzzzz = pbuffer.data(idx_kin_di + 104);

    auto tk_yy_yyyyyy = pbuffer.data(idx_kin_di + 105);

    auto tk_yy_yyyyyz = pbuffer.data(idx_kin_di + 106);

    auto tk_yy_yyyyzz = pbuffer.data(idx_kin_di + 107);

    auto tk_yy_yyyzzz = pbuffer.data(idx_kin_di + 108);

    auto tk_yy_yyzzzz = pbuffer.data(idx_kin_di + 109);

    auto tk_yy_yzzzzz = pbuffer.data(idx_kin_di + 110);

    auto tk_yy_zzzzzz = pbuffer.data(idx_kin_di + 111);

#pragma omp simd aligned(pa_y,             \
                             tk_0_xxxxxx,  \
                             tk_0_xxxxxy,  \
                             tk_0_xxxxxz,  \
                             tk_0_xxxxyy,  \
                             tk_0_xxxxyz,  \
                             tk_0_xxxxzz,  \
                             tk_0_xxxyyy,  \
                             tk_0_xxxyyz,  \
                             tk_0_xxxyzz,  \
                             tk_0_xxxzzz,  \
                             tk_0_xxyyyy,  \
                             tk_0_xxyyyz,  \
                             tk_0_xxyyzz,  \
                             tk_0_xxyzzz,  \
                             tk_0_xxzzzz,  \
                             tk_0_xyyyyy,  \
                             tk_0_xyyyyz,  \
                             tk_0_xyyyzz,  \
                             tk_0_xyyzzz,  \
                             tk_0_xyzzzz,  \
                             tk_0_xzzzzz,  \
                             tk_0_yyyyyy,  \
                             tk_0_yyyyyz,  \
                             tk_0_yyyyzz,  \
                             tk_0_yyyzzz,  \
                             tk_0_yyzzzz,  \
                             tk_0_yzzzzz,  \
                             tk_0_zzzzzz,  \
                             tk_y_xxxxx,   \
                             tk_y_xxxxxx,  \
                             tk_y_xxxxxy,  \
                             tk_y_xxxxxz,  \
                             tk_y_xxxxy,   \
                             tk_y_xxxxyy,  \
                             tk_y_xxxxyz,  \
                             tk_y_xxxxz,   \
                             tk_y_xxxxzz,  \
                             tk_y_xxxyy,   \
                             tk_y_xxxyyy,  \
                             tk_y_xxxyyz,  \
                             tk_y_xxxyz,   \
                             tk_y_xxxyzz,  \
                             tk_y_xxxzz,   \
                             tk_y_xxxzzz,  \
                             tk_y_xxyyy,   \
                             tk_y_xxyyyy,  \
                             tk_y_xxyyyz,  \
                             tk_y_xxyyz,   \
                             tk_y_xxyyzz,  \
                             tk_y_xxyzz,   \
                             tk_y_xxyzzz,  \
                             tk_y_xxzzz,   \
                             tk_y_xxzzzz,  \
                             tk_y_xyyyy,   \
                             tk_y_xyyyyy,  \
                             tk_y_xyyyyz,  \
                             tk_y_xyyyz,   \
                             tk_y_xyyyzz,  \
                             tk_y_xyyzz,   \
                             tk_y_xyyzzz,  \
                             tk_y_xyzzz,   \
                             tk_y_xyzzzz,  \
                             tk_y_xzzzz,   \
                             tk_y_xzzzzz,  \
                             tk_y_yyyyy,   \
                             tk_y_yyyyyy,  \
                             tk_y_yyyyyz,  \
                             tk_y_yyyyz,   \
                             tk_y_yyyyzz,  \
                             tk_y_yyyzz,   \
                             tk_y_yyyzzz,  \
                             tk_y_yyzzz,   \
                             tk_y_yyzzzz,  \
                             tk_y_yzzzz,   \
                             tk_y_yzzzzz,  \
                             tk_y_zzzzz,   \
                             tk_y_zzzzzz,  \
                             tk_yy_xxxxxx, \
                             tk_yy_xxxxxy, \
                             tk_yy_xxxxxz, \
                             tk_yy_xxxxyy, \
                             tk_yy_xxxxyz, \
                             tk_yy_xxxxzz, \
                             tk_yy_xxxyyy, \
                             tk_yy_xxxyyz, \
                             tk_yy_xxxyzz, \
                             tk_yy_xxxzzz, \
                             tk_yy_xxyyyy, \
                             tk_yy_xxyyyz, \
                             tk_yy_xxyyzz, \
                             tk_yy_xxyzzz, \
                             tk_yy_xxzzzz, \
                             tk_yy_xyyyyy, \
                             tk_yy_xyyyyz, \
                             tk_yy_xyyyzz, \
                             tk_yy_xyyzzz, \
                             tk_yy_xyzzzz, \
                             tk_yy_xzzzzz, \
                             tk_yy_yyyyyy, \
                             tk_yy_yyyyyz, \
                             tk_yy_yyyyzz, \
                             tk_yy_yyyzzz, \
                             tk_yy_yyzzzz, \
                             tk_yy_yzzzzz, \
                             tk_yy_zzzzzz, \
                             ts_0_xxxxxx,  \
                             ts_0_xxxxxy,  \
                             ts_0_xxxxxz,  \
                             ts_0_xxxxyy,  \
                             ts_0_xxxxyz,  \
                             ts_0_xxxxzz,  \
                             ts_0_xxxyyy,  \
                             ts_0_xxxyyz,  \
                             ts_0_xxxyzz,  \
                             ts_0_xxxzzz,  \
                             ts_0_xxyyyy,  \
                             ts_0_xxyyyz,  \
                             ts_0_xxyyzz,  \
                             ts_0_xxyzzz,  \
                             ts_0_xxzzzz,  \
                             ts_0_xyyyyy,  \
                             ts_0_xyyyyz,  \
                             ts_0_xyyyzz,  \
                             ts_0_xyyzzz,  \
                             ts_0_xyzzzz,  \
                             ts_0_xzzzzz,  \
                             ts_0_yyyyyy,  \
                             ts_0_yyyyyz,  \
                             ts_0_yyyyzz,  \
                             ts_0_yyyzzz,  \
                             ts_0_yyzzzz,  \
                             ts_0_yzzzzz,  \
                             ts_0_zzzzzz,  \
                             ts_yy_xxxxxx, \
                             ts_yy_xxxxxy, \
                             ts_yy_xxxxxz, \
                             ts_yy_xxxxyy, \
                             ts_yy_xxxxyz, \
                             ts_yy_xxxxzz, \
                             ts_yy_xxxyyy, \
                             ts_yy_xxxyyz, \
                             ts_yy_xxxyzz, \
                             ts_yy_xxxzzz, \
                             ts_yy_xxyyyy, \
                             ts_yy_xxyyyz, \
                             ts_yy_xxyyzz, \
                             ts_yy_xxyzzz, \
                             ts_yy_xxzzzz, \
                             ts_yy_xyyyyy, \
                             ts_yy_xyyyyz, \
                             ts_yy_xyyyzz, \
                             ts_yy_xyyzzz, \
                             ts_yy_xyzzzz, \
                             ts_yy_xzzzzz, \
                             ts_yy_yyyyyy, \
                             ts_yy_yyyyyz, \
                             ts_yy_yyyyzz, \
                             ts_yy_yyyzzz, \
                             ts_yy_yyzzzz, \
                             ts_yy_yzzzzz, \
                             ts_yy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_xxxxxx[i] = -2.0 * ts_0_xxxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxxx[i] * fe_0 + tk_y_xxxxxx[i] * pa_y[i] + 2.0 * ts_yy_xxxxxx[i] * fz_0;

        tk_yy_xxxxxy[i] = -2.0 * ts_0_xxxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxxy[i] * fe_0 + tk_y_xxxxx[i] * fe_0 + tk_y_xxxxxy[i] * pa_y[i] +
                          2.0 * ts_yy_xxxxxy[i] * fz_0;

        tk_yy_xxxxxz[i] = -2.0 * ts_0_xxxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxxz[i] * fe_0 + tk_y_xxxxxz[i] * pa_y[i] + 2.0 * ts_yy_xxxxxz[i] * fz_0;

        tk_yy_xxxxyy[i] = -2.0 * ts_0_xxxxyy[i] * fbe_0 * fz_0 + tk_0_xxxxyy[i] * fe_0 + 2.0 * tk_y_xxxxy[i] * fe_0 + tk_y_xxxxyy[i] * pa_y[i] +
                          2.0 * ts_yy_xxxxyy[i] * fz_0;

        tk_yy_xxxxyz[i] = -2.0 * ts_0_xxxxyz[i] * fbe_0 * fz_0 + tk_0_xxxxyz[i] * fe_0 + tk_y_xxxxz[i] * fe_0 + tk_y_xxxxyz[i] * pa_y[i] +
                          2.0 * ts_yy_xxxxyz[i] * fz_0;

        tk_yy_xxxxzz[i] = -2.0 * ts_0_xxxxzz[i] * fbe_0 * fz_0 + tk_0_xxxxzz[i] * fe_0 + tk_y_xxxxzz[i] * pa_y[i] + 2.0 * ts_yy_xxxxzz[i] * fz_0;

        tk_yy_xxxyyy[i] = -2.0 * ts_0_xxxyyy[i] * fbe_0 * fz_0 + tk_0_xxxyyy[i] * fe_0 + 3.0 * tk_y_xxxyy[i] * fe_0 + tk_y_xxxyyy[i] * pa_y[i] +
                          2.0 * ts_yy_xxxyyy[i] * fz_0;

        tk_yy_xxxyyz[i] = -2.0 * ts_0_xxxyyz[i] * fbe_0 * fz_0 + tk_0_xxxyyz[i] * fe_0 + 2.0 * tk_y_xxxyz[i] * fe_0 + tk_y_xxxyyz[i] * pa_y[i] +
                          2.0 * ts_yy_xxxyyz[i] * fz_0;

        tk_yy_xxxyzz[i] = -2.0 * ts_0_xxxyzz[i] * fbe_0 * fz_0 + tk_0_xxxyzz[i] * fe_0 + tk_y_xxxzz[i] * fe_0 + tk_y_xxxyzz[i] * pa_y[i] +
                          2.0 * ts_yy_xxxyzz[i] * fz_0;

        tk_yy_xxxzzz[i] = -2.0 * ts_0_xxxzzz[i] * fbe_0 * fz_0 + tk_0_xxxzzz[i] * fe_0 + tk_y_xxxzzz[i] * pa_y[i] + 2.0 * ts_yy_xxxzzz[i] * fz_0;

        tk_yy_xxyyyy[i] = -2.0 * ts_0_xxyyyy[i] * fbe_0 * fz_0 + tk_0_xxyyyy[i] * fe_0 + 4.0 * tk_y_xxyyy[i] * fe_0 + tk_y_xxyyyy[i] * pa_y[i] +
                          2.0 * ts_yy_xxyyyy[i] * fz_0;

        tk_yy_xxyyyz[i] = -2.0 * ts_0_xxyyyz[i] * fbe_0 * fz_0 + tk_0_xxyyyz[i] * fe_0 + 3.0 * tk_y_xxyyz[i] * fe_0 + tk_y_xxyyyz[i] * pa_y[i] +
                          2.0 * ts_yy_xxyyyz[i] * fz_0;

        tk_yy_xxyyzz[i] = -2.0 * ts_0_xxyyzz[i] * fbe_0 * fz_0 + tk_0_xxyyzz[i] * fe_0 + 2.0 * tk_y_xxyzz[i] * fe_0 + tk_y_xxyyzz[i] * pa_y[i] +
                          2.0 * ts_yy_xxyyzz[i] * fz_0;

        tk_yy_xxyzzz[i] = -2.0 * ts_0_xxyzzz[i] * fbe_0 * fz_0 + tk_0_xxyzzz[i] * fe_0 + tk_y_xxzzz[i] * fe_0 + tk_y_xxyzzz[i] * pa_y[i] +
                          2.0 * ts_yy_xxyzzz[i] * fz_0;

        tk_yy_xxzzzz[i] = -2.0 * ts_0_xxzzzz[i] * fbe_0 * fz_0 + tk_0_xxzzzz[i] * fe_0 + tk_y_xxzzzz[i] * pa_y[i] + 2.0 * ts_yy_xxzzzz[i] * fz_0;

        tk_yy_xyyyyy[i] = -2.0 * ts_0_xyyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyyy[i] * fe_0 + 5.0 * tk_y_xyyyy[i] * fe_0 + tk_y_xyyyyy[i] * pa_y[i] +
                          2.0 * ts_yy_xyyyyy[i] * fz_0;

        tk_yy_xyyyyz[i] = -2.0 * ts_0_xyyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyyz[i] * fe_0 + 4.0 * tk_y_xyyyz[i] * fe_0 + tk_y_xyyyyz[i] * pa_y[i] +
                          2.0 * ts_yy_xyyyyz[i] * fz_0;

        tk_yy_xyyyzz[i] = -2.0 * ts_0_xyyyzz[i] * fbe_0 * fz_0 + tk_0_xyyyzz[i] * fe_0 + 3.0 * tk_y_xyyzz[i] * fe_0 + tk_y_xyyyzz[i] * pa_y[i] +
                          2.0 * ts_yy_xyyyzz[i] * fz_0;

        tk_yy_xyyzzz[i] = -2.0 * ts_0_xyyzzz[i] * fbe_0 * fz_0 + tk_0_xyyzzz[i] * fe_0 + 2.0 * tk_y_xyzzz[i] * fe_0 + tk_y_xyyzzz[i] * pa_y[i] +
                          2.0 * ts_yy_xyyzzz[i] * fz_0;

        tk_yy_xyzzzz[i] = -2.0 * ts_0_xyzzzz[i] * fbe_0 * fz_0 + tk_0_xyzzzz[i] * fe_0 + tk_y_xzzzz[i] * fe_0 + tk_y_xyzzzz[i] * pa_y[i] +
                          2.0 * ts_yy_xyzzzz[i] * fz_0;

        tk_yy_xzzzzz[i] = -2.0 * ts_0_xzzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzzz[i] * fe_0 + tk_y_xzzzzz[i] * pa_y[i] + 2.0 * ts_yy_xzzzzz[i] * fz_0;

        tk_yy_yyyyyy[i] = -2.0 * ts_0_yyyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyyy[i] * fe_0 + 6.0 * tk_y_yyyyy[i] * fe_0 + tk_y_yyyyyy[i] * pa_y[i] +
                          2.0 * ts_yy_yyyyyy[i] * fz_0;

        tk_yy_yyyyyz[i] = -2.0 * ts_0_yyyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyyz[i] * fe_0 + 5.0 * tk_y_yyyyz[i] * fe_0 + tk_y_yyyyyz[i] * pa_y[i] +
                          2.0 * ts_yy_yyyyyz[i] * fz_0;

        tk_yy_yyyyzz[i] = -2.0 * ts_0_yyyyzz[i] * fbe_0 * fz_0 + tk_0_yyyyzz[i] * fe_0 + 4.0 * tk_y_yyyzz[i] * fe_0 + tk_y_yyyyzz[i] * pa_y[i] +
                          2.0 * ts_yy_yyyyzz[i] * fz_0;

        tk_yy_yyyzzz[i] = -2.0 * ts_0_yyyzzz[i] * fbe_0 * fz_0 + tk_0_yyyzzz[i] * fe_0 + 3.0 * tk_y_yyzzz[i] * fe_0 + tk_y_yyyzzz[i] * pa_y[i] +
                          2.0 * ts_yy_yyyzzz[i] * fz_0;

        tk_yy_yyzzzz[i] = -2.0 * ts_0_yyzzzz[i] * fbe_0 * fz_0 + tk_0_yyzzzz[i] * fe_0 + 2.0 * tk_y_yzzzz[i] * fe_0 + tk_y_yyzzzz[i] * pa_y[i] +
                          2.0 * ts_yy_yyzzzz[i] * fz_0;

        tk_yy_yzzzzz[i] = -2.0 * ts_0_yzzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzzz[i] * fe_0 + tk_y_zzzzz[i] * fe_0 + tk_y_yzzzzz[i] * pa_y[i] +
                          2.0 * ts_yy_yzzzzz[i] * fz_0;

        tk_yy_zzzzzz[i] = -2.0 * ts_0_zzzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzzz[i] * fe_0 + tk_y_zzzzzz[i] * pa_y[i] + 2.0 * ts_yy_zzzzzz[i] * fz_0;
    }

    // Set up 112-140 components of targeted buffer : DI

    auto tk_yz_xxxxxx = pbuffer.data(idx_kin_di + 112);

    auto tk_yz_xxxxxy = pbuffer.data(idx_kin_di + 113);

    auto tk_yz_xxxxxz = pbuffer.data(idx_kin_di + 114);

    auto tk_yz_xxxxyy = pbuffer.data(idx_kin_di + 115);

    auto tk_yz_xxxxyz = pbuffer.data(idx_kin_di + 116);

    auto tk_yz_xxxxzz = pbuffer.data(idx_kin_di + 117);

    auto tk_yz_xxxyyy = pbuffer.data(idx_kin_di + 118);

    auto tk_yz_xxxyyz = pbuffer.data(idx_kin_di + 119);

    auto tk_yz_xxxyzz = pbuffer.data(idx_kin_di + 120);

    auto tk_yz_xxxzzz = pbuffer.data(idx_kin_di + 121);

    auto tk_yz_xxyyyy = pbuffer.data(idx_kin_di + 122);

    auto tk_yz_xxyyyz = pbuffer.data(idx_kin_di + 123);

    auto tk_yz_xxyyzz = pbuffer.data(idx_kin_di + 124);

    auto tk_yz_xxyzzz = pbuffer.data(idx_kin_di + 125);

    auto tk_yz_xxzzzz = pbuffer.data(idx_kin_di + 126);

    auto tk_yz_xyyyyy = pbuffer.data(idx_kin_di + 127);

    auto tk_yz_xyyyyz = pbuffer.data(idx_kin_di + 128);

    auto tk_yz_xyyyzz = pbuffer.data(idx_kin_di + 129);

    auto tk_yz_xyyzzz = pbuffer.data(idx_kin_di + 130);

    auto tk_yz_xyzzzz = pbuffer.data(idx_kin_di + 131);

    auto tk_yz_xzzzzz = pbuffer.data(idx_kin_di + 132);

    auto tk_yz_yyyyyy = pbuffer.data(idx_kin_di + 133);

    auto tk_yz_yyyyyz = pbuffer.data(idx_kin_di + 134);

    auto tk_yz_yyyyzz = pbuffer.data(idx_kin_di + 135);

    auto tk_yz_yyyzzz = pbuffer.data(idx_kin_di + 136);

    auto tk_yz_yyzzzz = pbuffer.data(idx_kin_di + 137);

    auto tk_yz_yzzzzz = pbuffer.data(idx_kin_di + 138);

    auto tk_yz_zzzzzz = pbuffer.data(idx_kin_di + 139);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_y_xxxxxy,  \
                             tk_y_xxxxyy,  \
                             tk_y_xxxyyy,  \
                             tk_y_xxyyyy,  \
                             tk_y_xyyyyy,  \
                             tk_y_yyyyyy,  \
                             tk_yz_xxxxxx, \
                             tk_yz_xxxxxy, \
                             tk_yz_xxxxxz, \
                             tk_yz_xxxxyy, \
                             tk_yz_xxxxyz, \
                             tk_yz_xxxxzz, \
                             tk_yz_xxxyyy, \
                             tk_yz_xxxyyz, \
                             tk_yz_xxxyzz, \
                             tk_yz_xxxzzz, \
                             tk_yz_xxyyyy, \
                             tk_yz_xxyyyz, \
                             tk_yz_xxyyzz, \
                             tk_yz_xxyzzz, \
                             tk_yz_xxzzzz, \
                             tk_yz_xyyyyy, \
                             tk_yz_xyyyyz, \
                             tk_yz_xyyyzz, \
                             tk_yz_xyyzzz, \
                             tk_yz_xyzzzz, \
                             tk_yz_xzzzzz, \
                             tk_yz_yyyyyy, \
                             tk_yz_yyyyyz, \
                             tk_yz_yyyyzz, \
                             tk_yz_yyyzzz, \
                             tk_yz_yyzzzz, \
                             tk_yz_yzzzzz, \
                             tk_yz_zzzzzz, \
                             tk_z_xxxxxx,  \
                             tk_z_xxxxxz,  \
                             tk_z_xxxxyz,  \
                             tk_z_xxxxz,   \
                             tk_z_xxxxzz,  \
                             tk_z_xxxyyz,  \
                             tk_z_xxxyz,   \
                             tk_z_xxxyzz,  \
                             tk_z_xxxzz,   \
                             tk_z_xxxzzz,  \
                             tk_z_xxyyyz,  \
                             tk_z_xxyyz,   \
                             tk_z_xxyyzz,  \
                             tk_z_xxyzz,   \
                             tk_z_xxyzzz,  \
                             tk_z_xxzzz,   \
                             tk_z_xxzzzz,  \
                             tk_z_xyyyyz,  \
                             tk_z_xyyyz,   \
                             tk_z_xyyyzz,  \
                             tk_z_xyyzz,   \
                             tk_z_xyyzzz,  \
                             tk_z_xyzzz,   \
                             tk_z_xyzzzz,  \
                             tk_z_xzzzz,   \
                             tk_z_xzzzzz,  \
                             tk_z_yyyyyz,  \
                             tk_z_yyyyz,   \
                             tk_z_yyyyzz,  \
                             tk_z_yyyzz,   \
                             tk_z_yyyzzz,  \
                             tk_z_yyzzz,   \
                             tk_z_yyzzzz,  \
                             tk_z_yzzzz,   \
                             tk_z_yzzzzz,  \
                             tk_z_zzzzz,   \
                             tk_z_zzzzzz,  \
                             ts_yz_xxxxxx, \
                             ts_yz_xxxxxy, \
                             ts_yz_xxxxxz, \
                             ts_yz_xxxxyy, \
                             ts_yz_xxxxyz, \
                             ts_yz_xxxxzz, \
                             ts_yz_xxxyyy, \
                             ts_yz_xxxyyz, \
                             ts_yz_xxxyzz, \
                             ts_yz_xxxzzz, \
                             ts_yz_xxyyyy, \
                             ts_yz_xxyyyz, \
                             ts_yz_xxyyzz, \
                             ts_yz_xxyzzz, \
                             ts_yz_xxzzzz, \
                             ts_yz_xyyyyy, \
                             ts_yz_xyyyyz, \
                             ts_yz_xyyyzz, \
                             ts_yz_xyyzzz, \
                             ts_yz_xyzzzz, \
                             ts_yz_xzzzzz, \
                             ts_yz_yyyyyy, \
                             ts_yz_yyyyyz, \
                             ts_yz_yyyyzz, \
                             ts_yz_yyyzzz, \
                             ts_yz_yyzzzz, \
                             ts_yz_yzzzzz, \
                             ts_yz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yz_xxxxxx[i] = tk_z_xxxxxx[i] * pa_y[i] + 2.0 * ts_yz_xxxxxx[i] * fz_0;

        tk_yz_xxxxxy[i] = tk_y_xxxxxy[i] * pa_z[i] + 2.0 * ts_yz_xxxxxy[i] * fz_0;

        tk_yz_xxxxxz[i] = tk_z_xxxxxz[i] * pa_y[i] + 2.0 * ts_yz_xxxxxz[i] * fz_0;

        tk_yz_xxxxyy[i] = tk_y_xxxxyy[i] * pa_z[i] + 2.0 * ts_yz_xxxxyy[i] * fz_0;

        tk_yz_xxxxyz[i] = tk_z_xxxxz[i] * fe_0 + tk_z_xxxxyz[i] * pa_y[i] + 2.0 * ts_yz_xxxxyz[i] * fz_0;

        tk_yz_xxxxzz[i] = tk_z_xxxxzz[i] * pa_y[i] + 2.0 * ts_yz_xxxxzz[i] * fz_0;

        tk_yz_xxxyyy[i] = tk_y_xxxyyy[i] * pa_z[i] + 2.0 * ts_yz_xxxyyy[i] * fz_0;

        tk_yz_xxxyyz[i] = 2.0 * tk_z_xxxyz[i] * fe_0 + tk_z_xxxyyz[i] * pa_y[i] + 2.0 * ts_yz_xxxyyz[i] * fz_0;

        tk_yz_xxxyzz[i] = tk_z_xxxzz[i] * fe_0 + tk_z_xxxyzz[i] * pa_y[i] + 2.0 * ts_yz_xxxyzz[i] * fz_0;

        tk_yz_xxxzzz[i] = tk_z_xxxzzz[i] * pa_y[i] + 2.0 * ts_yz_xxxzzz[i] * fz_0;

        tk_yz_xxyyyy[i] = tk_y_xxyyyy[i] * pa_z[i] + 2.0 * ts_yz_xxyyyy[i] * fz_0;

        tk_yz_xxyyyz[i] = 3.0 * tk_z_xxyyz[i] * fe_0 + tk_z_xxyyyz[i] * pa_y[i] + 2.0 * ts_yz_xxyyyz[i] * fz_0;

        tk_yz_xxyyzz[i] = 2.0 * tk_z_xxyzz[i] * fe_0 + tk_z_xxyyzz[i] * pa_y[i] + 2.0 * ts_yz_xxyyzz[i] * fz_0;

        tk_yz_xxyzzz[i] = tk_z_xxzzz[i] * fe_0 + tk_z_xxyzzz[i] * pa_y[i] + 2.0 * ts_yz_xxyzzz[i] * fz_0;

        tk_yz_xxzzzz[i] = tk_z_xxzzzz[i] * pa_y[i] + 2.0 * ts_yz_xxzzzz[i] * fz_0;

        tk_yz_xyyyyy[i] = tk_y_xyyyyy[i] * pa_z[i] + 2.0 * ts_yz_xyyyyy[i] * fz_0;

        tk_yz_xyyyyz[i] = 4.0 * tk_z_xyyyz[i] * fe_0 + tk_z_xyyyyz[i] * pa_y[i] + 2.0 * ts_yz_xyyyyz[i] * fz_0;

        tk_yz_xyyyzz[i] = 3.0 * tk_z_xyyzz[i] * fe_0 + tk_z_xyyyzz[i] * pa_y[i] + 2.0 * ts_yz_xyyyzz[i] * fz_0;

        tk_yz_xyyzzz[i] = 2.0 * tk_z_xyzzz[i] * fe_0 + tk_z_xyyzzz[i] * pa_y[i] + 2.0 * ts_yz_xyyzzz[i] * fz_0;

        tk_yz_xyzzzz[i] = tk_z_xzzzz[i] * fe_0 + tk_z_xyzzzz[i] * pa_y[i] + 2.0 * ts_yz_xyzzzz[i] * fz_0;

        tk_yz_xzzzzz[i] = tk_z_xzzzzz[i] * pa_y[i] + 2.0 * ts_yz_xzzzzz[i] * fz_0;

        tk_yz_yyyyyy[i] = tk_y_yyyyyy[i] * pa_z[i] + 2.0 * ts_yz_yyyyyy[i] * fz_0;

        tk_yz_yyyyyz[i] = 5.0 * tk_z_yyyyz[i] * fe_0 + tk_z_yyyyyz[i] * pa_y[i] + 2.0 * ts_yz_yyyyyz[i] * fz_0;

        tk_yz_yyyyzz[i] = 4.0 * tk_z_yyyzz[i] * fe_0 + tk_z_yyyyzz[i] * pa_y[i] + 2.0 * ts_yz_yyyyzz[i] * fz_0;

        tk_yz_yyyzzz[i] = 3.0 * tk_z_yyzzz[i] * fe_0 + tk_z_yyyzzz[i] * pa_y[i] + 2.0 * ts_yz_yyyzzz[i] * fz_0;

        tk_yz_yyzzzz[i] = 2.0 * tk_z_yzzzz[i] * fe_0 + tk_z_yyzzzz[i] * pa_y[i] + 2.0 * ts_yz_yyzzzz[i] * fz_0;

        tk_yz_yzzzzz[i] = tk_z_zzzzz[i] * fe_0 + tk_z_yzzzzz[i] * pa_y[i] + 2.0 * ts_yz_yzzzzz[i] * fz_0;

        tk_yz_zzzzzz[i] = tk_z_zzzzzz[i] * pa_y[i] + 2.0 * ts_yz_zzzzzz[i] * fz_0;
    }

    // Set up 140-168 components of targeted buffer : DI

    auto tk_zz_xxxxxx = pbuffer.data(idx_kin_di + 140);

    auto tk_zz_xxxxxy = pbuffer.data(idx_kin_di + 141);

    auto tk_zz_xxxxxz = pbuffer.data(idx_kin_di + 142);

    auto tk_zz_xxxxyy = pbuffer.data(idx_kin_di + 143);

    auto tk_zz_xxxxyz = pbuffer.data(idx_kin_di + 144);

    auto tk_zz_xxxxzz = pbuffer.data(idx_kin_di + 145);

    auto tk_zz_xxxyyy = pbuffer.data(idx_kin_di + 146);

    auto tk_zz_xxxyyz = pbuffer.data(idx_kin_di + 147);

    auto tk_zz_xxxyzz = pbuffer.data(idx_kin_di + 148);

    auto tk_zz_xxxzzz = pbuffer.data(idx_kin_di + 149);

    auto tk_zz_xxyyyy = pbuffer.data(idx_kin_di + 150);

    auto tk_zz_xxyyyz = pbuffer.data(idx_kin_di + 151);

    auto tk_zz_xxyyzz = pbuffer.data(idx_kin_di + 152);

    auto tk_zz_xxyzzz = pbuffer.data(idx_kin_di + 153);

    auto tk_zz_xxzzzz = pbuffer.data(idx_kin_di + 154);

    auto tk_zz_xyyyyy = pbuffer.data(idx_kin_di + 155);

    auto tk_zz_xyyyyz = pbuffer.data(idx_kin_di + 156);

    auto tk_zz_xyyyzz = pbuffer.data(idx_kin_di + 157);

    auto tk_zz_xyyzzz = pbuffer.data(idx_kin_di + 158);

    auto tk_zz_xyzzzz = pbuffer.data(idx_kin_di + 159);

    auto tk_zz_xzzzzz = pbuffer.data(idx_kin_di + 160);

    auto tk_zz_yyyyyy = pbuffer.data(idx_kin_di + 161);

    auto tk_zz_yyyyyz = pbuffer.data(idx_kin_di + 162);

    auto tk_zz_yyyyzz = pbuffer.data(idx_kin_di + 163);

    auto tk_zz_yyyzzz = pbuffer.data(idx_kin_di + 164);

    auto tk_zz_yyzzzz = pbuffer.data(idx_kin_di + 165);

    auto tk_zz_yzzzzz = pbuffer.data(idx_kin_di + 166);

    auto tk_zz_zzzzzz = pbuffer.data(idx_kin_di + 167);

#pragma omp simd aligned(pa_z,             \
                             tk_0_xxxxxx,  \
                             tk_0_xxxxxy,  \
                             tk_0_xxxxxz,  \
                             tk_0_xxxxyy,  \
                             tk_0_xxxxyz,  \
                             tk_0_xxxxzz,  \
                             tk_0_xxxyyy,  \
                             tk_0_xxxyyz,  \
                             tk_0_xxxyzz,  \
                             tk_0_xxxzzz,  \
                             tk_0_xxyyyy,  \
                             tk_0_xxyyyz,  \
                             tk_0_xxyyzz,  \
                             tk_0_xxyzzz,  \
                             tk_0_xxzzzz,  \
                             tk_0_xyyyyy,  \
                             tk_0_xyyyyz,  \
                             tk_0_xyyyzz,  \
                             tk_0_xyyzzz,  \
                             tk_0_xyzzzz,  \
                             tk_0_xzzzzz,  \
                             tk_0_yyyyyy,  \
                             tk_0_yyyyyz,  \
                             tk_0_yyyyzz,  \
                             tk_0_yyyzzz,  \
                             tk_0_yyzzzz,  \
                             tk_0_yzzzzz,  \
                             tk_0_zzzzzz,  \
                             tk_z_xxxxx,   \
                             tk_z_xxxxxx,  \
                             tk_z_xxxxxy,  \
                             tk_z_xxxxxz,  \
                             tk_z_xxxxy,   \
                             tk_z_xxxxyy,  \
                             tk_z_xxxxyz,  \
                             tk_z_xxxxz,   \
                             tk_z_xxxxzz,  \
                             tk_z_xxxyy,   \
                             tk_z_xxxyyy,  \
                             tk_z_xxxyyz,  \
                             tk_z_xxxyz,   \
                             tk_z_xxxyzz,  \
                             tk_z_xxxzz,   \
                             tk_z_xxxzzz,  \
                             tk_z_xxyyy,   \
                             tk_z_xxyyyy,  \
                             tk_z_xxyyyz,  \
                             tk_z_xxyyz,   \
                             tk_z_xxyyzz,  \
                             tk_z_xxyzz,   \
                             tk_z_xxyzzz,  \
                             tk_z_xxzzz,   \
                             tk_z_xxzzzz,  \
                             tk_z_xyyyy,   \
                             tk_z_xyyyyy,  \
                             tk_z_xyyyyz,  \
                             tk_z_xyyyz,   \
                             tk_z_xyyyzz,  \
                             tk_z_xyyzz,   \
                             tk_z_xyyzzz,  \
                             tk_z_xyzzz,   \
                             tk_z_xyzzzz,  \
                             tk_z_xzzzz,   \
                             tk_z_xzzzzz,  \
                             tk_z_yyyyy,   \
                             tk_z_yyyyyy,  \
                             tk_z_yyyyyz,  \
                             tk_z_yyyyz,   \
                             tk_z_yyyyzz,  \
                             tk_z_yyyzz,   \
                             tk_z_yyyzzz,  \
                             tk_z_yyzzz,   \
                             tk_z_yyzzzz,  \
                             tk_z_yzzzz,   \
                             tk_z_yzzzzz,  \
                             tk_z_zzzzz,   \
                             tk_z_zzzzzz,  \
                             tk_zz_xxxxxx, \
                             tk_zz_xxxxxy, \
                             tk_zz_xxxxxz, \
                             tk_zz_xxxxyy, \
                             tk_zz_xxxxyz, \
                             tk_zz_xxxxzz, \
                             tk_zz_xxxyyy, \
                             tk_zz_xxxyyz, \
                             tk_zz_xxxyzz, \
                             tk_zz_xxxzzz, \
                             tk_zz_xxyyyy, \
                             tk_zz_xxyyyz, \
                             tk_zz_xxyyzz, \
                             tk_zz_xxyzzz, \
                             tk_zz_xxzzzz, \
                             tk_zz_xyyyyy, \
                             tk_zz_xyyyyz, \
                             tk_zz_xyyyzz, \
                             tk_zz_xyyzzz, \
                             tk_zz_xyzzzz, \
                             tk_zz_xzzzzz, \
                             tk_zz_yyyyyy, \
                             tk_zz_yyyyyz, \
                             tk_zz_yyyyzz, \
                             tk_zz_yyyzzz, \
                             tk_zz_yyzzzz, \
                             tk_zz_yzzzzz, \
                             tk_zz_zzzzzz, \
                             ts_0_xxxxxx,  \
                             ts_0_xxxxxy,  \
                             ts_0_xxxxxz,  \
                             ts_0_xxxxyy,  \
                             ts_0_xxxxyz,  \
                             ts_0_xxxxzz,  \
                             ts_0_xxxyyy,  \
                             ts_0_xxxyyz,  \
                             ts_0_xxxyzz,  \
                             ts_0_xxxzzz,  \
                             ts_0_xxyyyy,  \
                             ts_0_xxyyyz,  \
                             ts_0_xxyyzz,  \
                             ts_0_xxyzzz,  \
                             ts_0_xxzzzz,  \
                             ts_0_xyyyyy,  \
                             ts_0_xyyyyz,  \
                             ts_0_xyyyzz,  \
                             ts_0_xyyzzz,  \
                             ts_0_xyzzzz,  \
                             ts_0_xzzzzz,  \
                             ts_0_yyyyyy,  \
                             ts_0_yyyyyz,  \
                             ts_0_yyyyzz,  \
                             ts_0_yyyzzz,  \
                             ts_0_yyzzzz,  \
                             ts_0_yzzzzz,  \
                             ts_0_zzzzzz,  \
                             ts_zz_xxxxxx, \
                             ts_zz_xxxxxy, \
                             ts_zz_xxxxxz, \
                             ts_zz_xxxxyy, \
                             ts_zz_xxxxyz, \
                             ts_zz_xxxxzz, \
                             ts_zz_xxxyyy, \
                             ts_zz_xxxyyz, \
                             ts_zz_xxxyzz, \
                             ts_zz_xxxzzz, \
                             ts_zz_xxyyyy, \
                             ts_zz_xxyyyz, \
                             ts_zz_xxyyzz, \
                             ts_zz_xxyzzz, \
                             ts_zz_xxzzzz, \
                             ts_zz_xyyyyy, \
                             ts_zz_xyyyyz, \
                             ts_zz_xyyyzz, \
                             ts_zz_xyyzzz, \
                             ts_zz_xyzzzz, \
                             ts_zz_xzzzzz, \
                             ts_zz_yyyyyy, \
                             ts_zz_yyyyyz, \
                             ts_zz_yyyyzz, \
                             ts_zz_yyyzzz, \
                             ts_zz_yyzzzz, \
                             ts_zz_yzzzzz, \
                             ts_zz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_xxxxxx[i] = -2.0 * ts_0_xxxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxxx[i] * fe_0 + tk_z_xxxxxx[i] * pa_z[i] + 2.0 * ts_zz_xxxxxx[i] * fz_0;

        tk_zz_xxxxxy[i] = -2.0 * ts_0_xxxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxxy[i] * fe_0 + tk_z_xxxxxy[i] * pa_z[i] + 2.0 * ts_zz_xxxxxy[i] * fz_0;

        tk_zz_xxxxxz[i] = -2.0 * ts_0_xxxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxxz[i] * fe_0 + tk_z_xxxxx[i] * fe_0 + tk_z_xxxxxz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxxxz[i] * fz_0;

        tk_zz_xxxxyy[i] = -2.0 * ts_0_xxxxyy[i] * fbe_0 * fz_0 + tk_0_xxxxyy[i] * fe_0 + tk_z_xxxxyy[i] * pa_z[i] + 2.0 * ts_zz_xxxxyy[i] * fz_0;

        tk_zz_xxxxyz[i] = -2.0 * ts_0_xxxxyz[i] * fbe_0 * fz_0 + tk_0_xxxxyz[i] * fe_0 + tk_z_xxxxy[i] * fe_0 + tk_z_xxxxyz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxxyz[i] * fz_0;

        tk_zz_xxxxzz[i] = -2.0 * ts_0_xxxxzz[i] * fbe_0 * fz_0 + tk_0_xxxxzz[i] * fe_0 + 2.0 * tk_z_xxxxz[i] * fe_0 + tk_z_xxxxzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxxzz[i] * fz_0;

        tk_zz_xxxyyy[i] = -2.0 * ts_0_xxxyyy[i] * fbe_0 * fz_0 + tk_0_xxxyyy[i] * fe_0 + tk_z_xxxyyy[i] * pa_z[i] + 2.0 * ts_zz_xxxyyy[i] * fz_0;

        tk_zz_xxxyyz[i] = -2.0 * ts_0_xxxyyz[i] * fbe_0 * fz_0 + tk_0_xxxyyz[i] * fe_0 + tk_z_xxxyy[i] * fe_0 + tk_z_xxxyyz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxyyz[i] * fz_0;

        tk_zz_xxxyzz[i] = -2.0 * ts_0_xxxyzz[i] * fbe_0 * fz_0 + tk_0_xxxyzz[i] * fe_0 + 2.0 * tk_z_xxxyz[i] * fe_0 + tk_z_xxxyzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxyzz[i] * fz_0;

        tk_zz_xxxzzz[i] = -2.0 * ts_0_xxxzzz[i] * fbe_0 * fz_0 + tk_0_xxxzzz[i] * fe_0 + 3.0 * tk_z_xxxzz[i] * fe_0 + tk_z_xxxzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxxzzz[i] * fz_0;

        tk_zz_xxyyyy[i] = -2.0 * ts_0_xxyyyy[i] * fbe_0 * fz_0 + tk_0_xxyyyy[i] * fe_0 + tk_z_xxyyyy[i] * pa_z[i] + 2.0 * ts_zz_xxyyyy[i] * fz_0;

        tk_zz_xxyyyz[i] = -2.0 * ts_0_xxyyyz[i] * fbe_0 * fz_0 + tk_0_xxyyyz[i] * fe_0 + tk_z_xxyyy[i] * fe_0 + tk_z_xxyyyz[i] * pa_z[i] +
                          2.0 * ts_zz_xxyyyz[i] * fz_0;

        tk_zz_xxyyzz[i] = -2.0 * ts_0_xxyyzz[i] * fbe_0 * fz_0 + tk_0_xxyyzz[i] * fe_0 + 2.0 * tk_z_xxyyz[i] * fe_0 + tk_z_xxyyzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxyyzz[i] * fz_0;

        tk_zz_xxyzzz[i] = -2.0 * ts_0_xxyzzz[i] * fbe_0 * fz_0 + tk_0_xxyzzz[i] * fe_0 + 3.0 * tk_z_xxyzz[i] * fe_0 + tk_z_xxyzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxyzzz[i] * fz_0;

        tk_zz_xxzzzz[i] = -2.0 * ts_0_xxzzzz[i] * fbe_0 * fz_0 + tk_0_xxzzzz[i] * fe_0 + 4.0 * tk_z_xxzzz[i] * fe_0 + tk_z_xxzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xxzzzz[i] * fz_0;

        tk_zz_xyyyyy[i] = -2.0 * ts_0_xyyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyyy[i] * fe_0 + tk_z_xyyyyy[i] * pa_z[i] + 2.0 * ts_zz_xyyyyy[i] * fz_0;

        tk_zz_xyyyyz[i] = -2.0 * ts_0_xyyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyyz[i] * fe_0 + tk_z_xyyyy[i] * fe_0 + tk_z_xyyyyz[i] * pa_z[i] +
                          2.0 * ts_zz_xyyyyz[i] * fz_0;

        tk_zz_xyyyzz[i] = -2.0 * ts_0_xyyyzz[i] * fbe_0 * fz_0 + tk_0_xyyyzz[i] * fe_0 + 2.0 * tk_z_xyyyz[i] * fe_0 + tk_z_xyyyzz[i] * pa_z[i] +
                          2.0 * ts_zz_xyyyzz[i] * fz_0;

        tk_zz_xyyzzz[i] = -2.0 * ts_0_xyyzzz[i] * fbe_0 * fz_0 + tk_0_xyyzzz[i] * fe_0 + 3.0 * tk_z_xyyzz[i] * fe_0 + tk_z_xyyzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xyyzzz[i] * fz_0;

        tk_zz_xyzzzz[i] = -2.0 * ts_0_xyzzzz[i] * fbe_0 * fz_0 + tk_0_xyzzzz[i] * fe_0 + 4.0 * tk_z_xyzzz[i] * fe_0 + tk_z_xyzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xyzzzz[i] * fz_0;

        tk_zz_xzzzzz[i] = -2.0 * ts_0_xzzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzzz[i] * fe_0 + 5.0 * tk_z_xzzzz[i] * fe_0 + tk_z_xzzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_xzzzzz[i] * fz_0;

        tk_zz_yyyyyy[i] = -2.0 * ts_0_yyyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyyy[i] * fe_0 + tk_z_yyyyyy[i] * pa_z[i] + 2.0 * ts_zz_yyyyyy[i] * fz_0;

        tk_zz_yyyyyz[i] = -2.0 * ts_0_yyyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyyz[i] * fe_0 + tk_z_yyyyy[i] * fe_0 + tk_z_yyyyyz[i] * pa_z[i] +
                          2.0 * ts_zz_yyyyyz[i] * fz_0;

        tk_zz_yyyyzz[i] = -2.0 * ts_0_yyyyzz[i] * fbe_0 * fz_0 + tk_0_yyyyzz[i] * fe_0 + 2.0 * tk_z_yyyyz[i] * fe_0 + tk_z_yyyyzz[i] * pa_z[i] +
                          2.0 * ts_zz_yyyyzz[i] * fz_0;

        tk_zz_yyyzzz[i] = -2.0 * ts_0_yyyzzz[i] * fbe_0 * fz_0 + tk_0_yyyzzz[i] * fe_0 + 3.0 * tk_z_yyyzz[i] * fe_0 + tk_z_yyyzzz[i] * pa_z[i] +
                          2.0 * ts_zz_yyyzzz[i] * fz_0;

        tk_zz_yyzzzz[i] = -2.0 * ts_0_yyzzzz[i] * fbe_0 * fz_0 + tk_0_yyzzzz[i] * fe_0 + 4.0 * tk_z_yyzzz[i] * fe_0 + tk_z_yyzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_yyzzzz[i] * fz_0;

        tk_zz_yzzzzz[i] = -2.0 * ts_0_yzzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzzz[i] * fe_0 + 5.0 * tk_z_yzzzz[i] * fe_0 + tk_z_yzzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_yzzzzz[i] * fz_0;

        tk_zz_zzzzzz[i] = -2.0 * ts_0_zzzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzzz[i] * fe_0 + 6.0 * tk_z_zzzzz[i] * fe_0 + tk_z_zzzzzz[i] * pa_z[i] +
                          2.0 * ts_zz_zzzzzz[i] * fz_0;
    }
}

}  // namespace kinrec

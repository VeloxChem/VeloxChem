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

#include "ThreeCenterElectronRepulsionPrimRecPSM.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psm,
                                 size_t idx_eri_1_ssl,
                                 size_t idx_eri_1_ssm,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_ssl);

    auto g_0_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_ssl + 1);

    auto g_0_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_ssl + 2);

    auto g_0_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_ssl + 3);

    auto g_0_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_ssl + 4);

    auto g_0_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_ssl + 5);

    auto g_0_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_ssl + 6);

    auto g_0_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_ssl + 7);

    auto g_0_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_ssl + 8);

    auto g_0_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_ssl + 9);

    auto g_0_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_ssl + 10);

    auto g_0_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_ssl + 11);

    auto g_0_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_ssl + 12);

    auto g_0_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_ssl + 13);

    auto g_0_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_ssl + 14);

    auto g_0_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 15);

    auto g_0_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 16);

    auto g_0_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 17);

    auto g_0_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 18);

    auto g_0_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 19);

    auto g_0_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 20);

    auto g_0_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 21);

    auto g_0_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 22);

    auto g_0_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 23);

    auto g_0_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 24);

    auto g_0_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 25);

    auto g_0_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 26);

    auto g_0_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 27);

    auto g_0_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 28);

    auto g_0_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 29);

    auto g_0_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 30);

    auto g_0_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 31);

    auto g_0_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 32);

    auto g_0_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 33);

    auto g_0_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 34);

    auto g_0_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 35);

    auto g_0_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 36);

    auto g_0_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 37);

    auto g_0_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 38);

    auto g_0_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 39);

    auto g_0_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 40);

    auto g_0_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 41);

    auto g_0_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 42);

    auto g_0_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 43);

    auto g_0_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 44);

    /// Set up components of auxilary buffer : SSM

    auto g_0_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_ssm);

    auto g_0_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_ssm + 1);

    auto g_0_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_ssm + 2);

    auto g_0_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_ssm + 3);

    auto g_0_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_ssm + 4);

    auto g_0_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_ssm + 5);

    auto g_0_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_ssm + 6);

    auto g_0_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_ssm + 7);

    auto g_0_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_ssm + 8);

    auto g_0_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_ssm + 9);

    auto g_0_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_ssm + 10);

    auto g_0_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_ssm + 11);

    auto g_0_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_ssm + 12);

    auto g_0_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_ssm + 13);

    auto g_0_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_ssm + 14);

    auto g_0_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 15);

    auto g_0_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 16);

    auto g_0_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 17);

    auto g_0_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 18);

    auto g_0_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 19);

    auto g_0_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 20);

    auto g_0_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 21);

    auto g_0_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 22);

    auto g_0_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 23);

    auto g_0_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 24);

    auto g_0_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 25);

    auto g_0_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 26);

    auto g_0_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 27);

    auto g_0_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 28);

    auto g_0_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 29);

    auto g_0_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 30);

    auto g_0_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 31);

    auto g_0_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 32);

    auto g_0_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 33);

    auto g_0_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 34);

    auto g_0_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 35);

    auto g_0_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 36);

    auto g_0_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 37);

    auto g_0_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 38);

    auto g_0_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 39);

    auto g_0_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 40);

    auto g_0_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 41);

    auto g_0_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 42);

    auto g_0_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 43);

    auto g_0_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 44);

    auto g_0_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_ssm + 45);

    auto g_0_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_ssm + 46);

    auto g_0_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_ssm + 47);

    auto g_0_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_ssm + 48);

    auto g_0_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_ssm + 49);

    auto g_0_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 50);

    auto g_0_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 51);

    auto g_0_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 52);

    auto g_0_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 53);

    auto g_0_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_ssm + 54);

    /// Set up 0-55 components of targeted buffer : PSM

    auto g_x_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm);

    auto g_x_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 1);

    auto g_x_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 2);

    auto g_x_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 3);

    auto g_x_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 4);

    auto g_x_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 5);

    auto g_x_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 6);

    auto g_x_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 7);

    auto g_x_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 8);

    auto g_x_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 9);

    auto g_x_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 10);

    auto g_x_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 11);

    auto g_x_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 12);

    auto g_x_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 13);

    auto g_x_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 14);

    auto g_x_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 15);

    auto g_x_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 16);

    auto g_x_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 17);

    auto g_x_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 18);

    auto g_x_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 19);

    auto g_x_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 20);

    auto g_x_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 21);

    auto g_x_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 22);

    auto g_x_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 23);

    auto g_x_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 24);

    auto g_x_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 25);

    auto g_x_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 26);

    auto g_x_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 27);

    auto g_x_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 28);

    auto g_x_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 29);

    auto g_x_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 30);

    auto g_x_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 31);

    auto g_x_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 32);

    auto g_x_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 33);

    auto g_x_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 34);

    auto g_x_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 35);

    auto g_x_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 36);

    auto g_x_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 37);

    auto g_x_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 38);

    auto g_x_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 39);

    auto g_x_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 40);

    auto g_x_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 41);

    auto g_x_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 42);

    auto g_x_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 43);

    auto g_x_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 44);

    auto g_x_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 45);

    auto g_x_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 46);

    auto g_x_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 47);

    auto g_x_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 48);

    auto g_x_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 49);

    auto g_x_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 50);

    auto g_x_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 51);

    auto g_x_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 52);

    auto g_x_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 53);

    auto g_x_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 54);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzz_1, g_0_0_zzzzzzzzz_1, g_x_0_xxxxxxxxx_0, g_x_0_xxxxxxxxy_0, g_x_0_xxxxxxxxz_0, g_x_0_xxxxxxxyy_0, g_x_0_xxxxxxxyz_0, g_x_0_xxxxxxxzz_0, g_x_0_xxxxxxyyy_0, g_x_0_xxxxxxyyz_0, g_x_0_xxxxxxyzz_0, g_x_0_xxxxxxzzz_0, g_x_0_xxxxxyyyy_0, g_x_0_xxxxxyyyz_0, g_x_0_xxxxxyyzz_0, g_x_0_xxxxxyzzz_0, g_x_0_xxxxxzzzz_0, g_x_0_xxxxyyyyy_0, g_x_0_xxxxyyyyz_0, g_x_0_xxxxyyyzz_0, g_x_0_xxxxyyzzz_0, g_x_0_xxxxyzzzz_0, g_x_0_xxxxzzzzz_0, g_x_0_xxxyyyyyy_0, g_x_0_xxxyyyyyz_0, g_x_0_xxxyyyyzz_0, g_x_0_xxxyyyzzz_0, g_x_0_xxxyyzzzz_0, g_x_0_xxxyzzzzz_0, g_x_0_xxxzzzzzz_0, g_x_0_xxyyyyyyy_0, g_x_0_xxyyyyyyz_0, g_x_0_xxyyyyyzz_0, g_x_0_xxyyyyzzz_0, g_x_0_xxyyyzzzz_0, g_x_0_xxyyzzzzz_0, g_x_0_xxyzzzzzz_0, g_x_0_xxzzzzzzz_0, g_x_0_xyyyyyyyy_0, g_x_0_xyyyyyyyz_0, g_x_0_xyyyyyyzz_0, g_x_0_xyyyyyzzz_0, g_x_0_xyyyyzzzz_0, g_x_0_xyyyzzzzz_0, g_x_0_xyyzzzzzz_0, g_x_0_xyzzzzzzz_0, g_x_0_xzzzzzzzz_0, g_x_0_yyyyyyyyy_0, g_x_0_yyyyyyyyz_0, g_x_0_yyyyyyyzz_0, g_x_0_yyyyyyzzz_0, g_x_0_yyyyyzzzz_0, g_x_0_yyyyzzzzz_0, g_x_0_yyyzzzzzz_0, g_x_0_yyzzzzzzz_0, g_x_0_yzzzzzzzz_0, g_x_0_zzzzzzzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxxxxxxx_0[i] = 9.0 * g_0_0_xxxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxxx_1[i] * wa_x[i];

        g_x_0_xxxxxxxxy_0[i] = 8.0 * g_0_0_xxxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxxxy_1[i] * wa_x[i];

        g_x_0_xxxxxxxxz_0[i] = 8.0 * g_0_0_xxxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxxxz_1[i] * wa_x[i];

        g_x_0_xxxxxxxyy_0[i] = 7.0 * g_0_0_xxxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxxxyy_1[i] * wa_x[i];

        g_x_0_xxxxxxxyz_0[i] = 7.0 * g_0_0_xxxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxxxyz_1[i] * wa_x[i];

        g_x_0_xxxxxxxzz_0[i] = 7.0 * g_0_0_xxxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxxxzz_1[i] * wa_x[i];

        g_x_0_xxxxxxyyy_0[i] = 6.0 * g_0_0_xxxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxxxyyy_1[i] * wa_x[i];

        g_x_0_xxxxxxyyz_0[i] = 6.0 * g_0_0_xxxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxxxyyz_1[i] * wa_x[i];

        g_x_0_xxxxxxyzz_0[i] = 6.0 * g_0_0_xxxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxxxyzz_1[i] * wa_x[i];

        g_x_0_xxxxxxzzz_0[i] = 6.0 * g_0_0_xxxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxxxzzz_1[i] * wa_x[i];

        g_x_0_xxxxxyyyy_0[i] = 5.0 * g_0_0_xxxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyyy_1[i] * wa_x[i];

        g_x_0_xxxxxyyyz_0[i] = 5.0 * g_0_0_xxxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxxxyyyz_1[i] * wa_x[i];

        g_x_0_xxxxxyyzz_0[i] = 5.0 * g_0_0_xxxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxxxyyzz_1[i] * wa_x[i];

        g_x_0_xxxxxyzzz_0[i] = 5.0 * g_0_0_xxxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxxxyzzz_1[i] * wa_x[i];

        g_x_0_xxxxxzzzz_0[i] = 5.0 * g_0_0_xxxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxxxzzzz_1[i] * wa_x[i];

        g_x_0_xxxxyyyyy_0[i] = 4.0 * g_0_0_xxxyyyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyyy_1[i] * wa_x[i];

        g_x_0_xxxxyyyyz_0[i] = 4.0 * g_0_0_xxxyyyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyyyz_1[i] * wa_x[i];

        g_x_0_xxxxyyyzz_0[i] = 4.0 * g_0_0_xxxyyyzz_1[i] * fi_acd_0 + g_0_0_xxxxyyyzz_1[i] * wa_x[i];

        g_x_0_xxxxyyzzz_0[i] = 4.0 * g_0_0_xxxyyzzz_1[i] * fi_acd_0 + g_0_0_xxxxyyzzz_1[i] * wa_x[i];

        g_x_0_xxxxyzzzz_0[i] = 4.0 * g_0_0_xxxyzzzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzzz_1[i] * wa_x[i];

        g_x_0_xxxxzzzzz_0[i] = 4.0 * g_0_0_xxxzzzzz_1[i] * fi_acd_0 + g_0_0_xxxxzzzzz_1[i] * wa_x[i];

        g_x_0_xxxyyyyyy_0[i] = 3.0 * g_0_0_xxyyyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyyy_1[i] * wa_x[i];

        g_x_0_xxxyyyyyz_0[i] = 3.0 * g_0_0_xxyyyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyyyz_1[i] * wa_x[i];

        g_x_0_xxxyyyyzz_0[i] = 3.0 * g_0_0_xxyyyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyyyzz_1[i] * wa_x[i];

        g_x_0_xxxyyyzzz_0[i] = 3.0 * g_0_0_xxyyyzzz_1[i] * fi_acd_0 + g_0_0_xxxyyyzzz_1[i] * wa_x[i];

        g_x_0_xxxyyzzzz_0[i] = 3.0 * g_0_0_xxyyzzzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzzz_1[i] * wa_x[i];

        g_x_0_xxxyzzzzz_0[i] = 3.0 * g_0_0_xxyzzzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzzz_1[i] * wa_x[i];

        g_x_0_xxxzzzzzz_0[i] = 3.0 * g_0_0_xxzzzzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzzzz_1[i] * wa_x[i];

        g_x_0_xxyyyyyyy_0[i] = 2.0 * g_0_0_xyyyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyyy_1[i] * wa_x[i];

        g_x_0_xxyyyyyyz_0[i] = 2.0 * g_0_0_xyyyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyyyz_1[i] * wa_x[i];

        g_x_0_xxyyyyyzz_0[i] = 2.0 * g_0_0_xyyyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyyyzz_1[i] * wa_x[i];

        g_x_0_xxyyyyzzz_0[i] = 2.0 * g_0_0_xyyyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyyyzzz_1[i] * wa_x[i];

        g_x_0_xxyyyzzzz_0[i] = 2.0 * g_0_0_xyyyzzzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzzz_1[i] * wa_x[i];

        g_x_0_xxyyzzzzz_0[i] = 2.0 * g_0_0_xyyzzzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzzz_1[i] * wa_x[i];

        g_x_0_xxyzzzzzz_0[i] = 2.0 * g_0_0_xyzzzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzzz_1[i] * wa_x[i];

        g_x_0_xxzzzzzzz_0[i] = 2.0 * g_0_0_xzzzzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzzzz_1[i] * wa_x[i];

        g_x_0_xyyyyyyyy_0[i] = g_0_0_yyyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyyy_1[i] * wa_x[i];

        g_x_0_xyyyyyyyz_0[i] = g_0_0_yyyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyyyz_1[i] * wa_x[i];

        g_x_0_xyyyyyyzz_0[i] = g_0_0_yyyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyyyzz_1[i] * wa_x[i];

        g_x_0_xyyyyyzzz_0[i] = g_0_0_yyyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyyyzzz_1[i] * wa_x[i];

        g_x_0_xyyyyzzzz_0[i] = g_0_0_yyyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzzz_1[i] * wa_x[i];

        g_x_0_xyyyzzzzz_0[i] = g_0_0_yyyzzzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzzz_1[i] * wa_x[i];

        g_x_0_xyyzzzzzz_0[i] = g_0_0_yyzzzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzzz_1[i] * wa_x[i];

        g_x_0_xyzzzzzzz_0[i] = g_0_0_yzzzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzzz_1[i] * wa_x[i];

        g_x_0_xzzzzzzzz_0[i] = g_0_0_zzzzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzzzz_1[i] * wa_x[i];

        g_x_0_yyyyyyyyy_0[i] = g_0_0_yyyyyyyyy_1[i] * wa_x[i];

        g_x_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyyz_1[i] * wa_x[i];

        g_x_0_yyyyyyyzz_0[i] = g_0_0_yyyyyyyzz_1[i] * wa_x[i];

        g_x_0_yyyyyyzzz_0[i] = g_0_0_yyyyyyzzz_1[i] * wa_x[i];

        g_x_0_yyyyyzzzz_0[i] = g_0_0_yyyyyzzzz_1[i] * wa_x[i];

        g_x_0_yyyyzzzzz_0[i] = g_0_0_yyyyzzzzz_1[i] * wa_x[i];

        g_x_0_yyyzzzzzz_0[i] = g_0_0_yyyzzzzzz_1[i] * wa_x[i];

        g_x_0_yyzzzzzzz_0[i] = g_0_0_yyzzzzzzz_1[i] * wa_x[i];

        g_x_0_yzzzzzzzz_0[i] = g_0_0_yzzzzzzzz_1[i] * wa_x[i];

        g_x_0_zzzzzzzzz_0[i] = g_0_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 55-110 components of targeted buffer : PSM

    auto g_y_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm + 55);

    auto g_y_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 56);

    auto g_y_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 57);

    auto g_y_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 58);

    auto g_y_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 59);

    auto g_y_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 60);

    auto g_y_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 61);

    auto g_y_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 62);

    auto g_y_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 63);

    auto g_y_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 64);

    auto g_y_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 65);

    auto g_y_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 66);

    auto g_y_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 67);

    auto g_y_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 68);

    auto g_y_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 69);

    auto g_y_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 70);

    auto g_y_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 71);

    auto g_y_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 72);

    auto g_y_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 73);

    auto g_y_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 74);

    auto g_y_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 75);

    auto g_y_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 76);

    auto g_y_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 77);

    auto g_y_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 78);

    auto g_y_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 79);

    auto g_y_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 80);

    auto g_y_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 81);

    auto g_y_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 82);

    auto g_y_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 83);

    auto g_y_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 84);

    auto g_y_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 85);

    auto g_y_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 86);

    auto g_y_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 87);

    auto g_y_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 88);

    auto g_y_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 89);

    auto g_y_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 90);

    auto g_y_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 91);

    auto g_y_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 92);

    auto g_y_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 93);

    auto g_y_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 94);

    auto g_y_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 95);

    auto g_y_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 96);

    auto g_y_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 97);

    auto g_y_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 98);

    auto g_y_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 99);

    auto g_y_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 100);

    auto g_y_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 101);

    auto g_y_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 102);

    auto g_y_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 103);

    auto g_y_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 104);

    auto g_y_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 105);

    auto g_y_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 106);

    auto g_y_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 107);

    auto g_y_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 108);

    auto g_y_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 109);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzz_1, g_0_0_zzzzzzzzz_1, g_y_0_xxxxxxxxx_0, g_y_0_xxxxxxxxy_0, g_y_0_xxxxxxxxz_0, g_y_0_xxxxxxxyy_0, g_y_0_xxxxxxxyz_0, g_y_0_xxxxxxxzz_0, g_y_0_xxxxxxyyy_0, g_y_0_xxxxxxyyz_0, g_y_0_xxxxxxyzz_0, g_y_0_xxxxxxzzz_0, g_y_0_xxxxxyyyy_0, g_y_0_xxxxxyyyz_0, g_y_0_xxxxxyyzz_0, g_y_0_xxxxxyzzz_0, g_y_0_xxxxxzzzz_0, g_y_0_xxxxyyyyy_0, g_y_0_xxxxyyyyz_0, g_y_0_xxxxyyyzz_0, g_y_0_xxxxyyzzz_0, g_y_0_xxxxyzzzz_0, g_y_0_xxxxzzzzz_0, g_y_0_xxxyyyyyy_0, g_y_0_xxxyyyyyz_0, g_y_0_xxxyyyyzz_0, g_y_0_xxxyyyzzz_0, g_y_0_xxxyyzzzz_0, g_y_0_xxxyzzzzz_0, g_y_0_xxxzzzzzz_0, g_y_0_xxyyyyyyy_0, g_y_0_xxyyyyyyz_0, g_y_0_xxyyyyyzz_0, g_y_0_xxyyyyzzz_0, g_y_0_xxyyyzzzz_0, g_y_0_xxyyzzzzz_0, g_y_0_xxyzzzzzz_0, g_y_0_xxzzzzzzz_0, g_y_0_xyyyyyyyy_0, g_y_0_xyyyyyyyz_0, g_y_0_xyyyyyyzz_0, g_y_0_xyyyyyzzz_0, g_y_0_xyyyyzzzz_0, g_y_0_xyyyzzzzz_0, g_y_0_xyyzzzzzz_0, g_y_0_xyzzzzzzz_0, g_y_0_xzzzzzzzz_0, g_y_0_yyyyyyyyy_0, g_y_0_yyyyyyyyz_0, g_y_0_yyyyyyyzz_0, g_y_0_yyyyyyzzz_0, g_y_0_yyyyyzzzz_0, g_y_0_yyyyzzzzz_0, g_y_0_yyyzzzzzz_0, g_y_0_yyzzzzzzz_0, g_y_0_yzzzzzzzz_0, g_y_0_zzzzzzzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxxxxxxx_0[i] = g_0_0_xxxxxxxxx_1[i] * wa_y[i];

        g_y_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxxy_1[i] * wa_y[i];

        g_y_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxxz_1[i] * wa_y[i];

        g_y_0_xxxxxxxyy_0[i] = 2.0 * g_0_0_xxxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxxyy_1[i] * wa_y[i];

        g_y_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxxyz_1[i] * wa_y[i];

        g_y_0_xxxxxxxzz_0[i] = g_0_0_xxxxxxxzz_1[i] * wa_y[i];

        g_y_0_xxxxxxyyy_0[i] = 3.0 * g_0_0_xxxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxxyyy_1[i] * wa_y[i];

        g_y_0_xxxxxxyyz_0[i] = 2.0 * g_0_0_xxxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxxyyz_1[i] * wa_y[i];

        g_y_0_xxxxxxyzz_0[i] = g_0_0_xxxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxxyzz_1[i] * wa_y[i];

        g_y_0_xxxxxxzzz_0[i] = g_0_0_xxxxxxzzz_1[i] * wa_y[i];

        g_y_0_xxxxxyyyy_0[i] = 4.0 * g_0_0_xxxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyyy_1[i] * wa_y[i];

        g_y_0_xxxxxyyyz_0[i] = 3.0 * g_0_0_xxxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxxyyyz_1[i] * wa_y[i];

        g_y_0_xxxxxyyzz_0[i] = 2.0 * g_0_0_xxxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxxyyzz_1[i] * wa_y[i];

        g_y_0_xxxxxyzzz_0[i] = g_0_0_xxxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxxyzzz_1[i] * wa_y[i];

        g_y_0_xxxxxzzzz_0[i] = g_0_0_xxxxxzzzz_1[i] * wa_y[i];

        g_y_0_xxxxyyyyy_0[i] = 5.0 * g_0_0_xxxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyyy_1[i] * wa_y[i];

        g_y_0_xxxxyyyyz_0[i] = 4.0 * g_0_0_xxxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyyyz_1[i] * wa_y[i];

        g_y_0_xxxxyyyzz_0[i] = 3.0 * g_0_0_xxxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxxyyyzz_1[i] * wa_y[i];

        g_y_0_xxxxyyzzz_0[i] = 2.0 * g_0_0_xxxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxxyyzzz_1[i] * wa_y[i];

        g_y_0_xxxxyzzzz_0[i] = g_0_0_xxxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzzz_1[i] * wa_y[i];

        g_y_0_xxxxzzzzz_0[i] = g_0_0_xxxxzzzzz_1[i] * wa_y[i];

        g_y_0_xxxyyyyyy_0[i] = 6.0 * g_0_0_xxxyyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyyy_1[i] * wa_y[i];

        g_y_0_xxxyyyyyz_0[i] = 5.0 * g_0_0_xxxyyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyyyz_1[i] * wa_y[i];

        g_y_0_xxxyyyyzz_0[i] = 4.0 * g_0_0_xxxyyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyyyzz_1[i] * wa_y[i];

        g_y_0_xxxyyyzzz_0[i] = 3.0 * g_0_0_xxxyyzzz_1[i] * fi_acd_0 + g_0_0_xxxyyyzzz_1[i] * wa_y[i];

        g_y_0_xxxyyzzzz_0[i] = 2.0 * g_0_0_xxxyzzzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzzz_1[i] * wa_y[i];

        g_y_0_xxxyzzzzz_0[i] = g_0_0_xxxzzzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzzz_1[i] * wa_y[i];

        g_y_0_xxxzzzzzz_0[i] = g_0_0_xxxzzzzzz_1[i] * wa_y[i];

        g_y_0_xxyyyyyyy_0[i] = 7.0 * g_0_0_xxyyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyyy_1[i] * wa_y[i];

        g_y_0_xxyyyyyyz_0[i] = 6.0 * g_0_0_xxyyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyyyz_1[i] * wa_y[i];

        g_y_0_xxyyyyyzz_0[i] = 5.0 * g_0_0_xxyyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyyyzz_1[i] * wa_y[i];

        g_y_0_xxyyyyzzz_0[i] = 4.0 * g_0_0_xxyyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyyyzzz_1[i] * wa_y[i];

        g_y_0_xxyyyzzzz_0[i] = 3.0 * g_0_0_xxyyzzzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzzz_1[i] * wa_y[i];

        g_y_0_xxyyzzzzz_0[i] = 2.0 * g_0_0_xxyzzzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzzz_1[i] * wa_y[i];

        g_y_0_xxyzzzzzz_0[i] = g_0_0_xxzzzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzzz_1[i] * wa_y[i];

        g_y_0_xxzzzzzzz_0[i] = g_0_0_xxzzzzzzz_1[i] * wa_y[i];

        g_y_0_xyyyyyyyy_0[i] = 8.0 * g_0_0_xyyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyyy_1[i] * wa_y[i];

        g_y_0_xyyyyyyyz_0[i] = 7.0 * g_0_0_xyyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyyyz_1[i] * wa_y[i];

        g_y_0_xyyyyyyzz_0[i] = 6.0 * g_0_0_xyyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyyyzz_1[i] * wa_y[i];

        g_y_0_xyyyyyzzz_0[i] = 5.0 * g_0_0_xyyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyyyzzz_1[i] * wa_y[i];

        g_y_0_xyyyyzzzz_0[i] = 4.0 * g_0_0_xyyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzzz_1[i] * wa_y[i];

        g_y_0_xyyyzzzzz_0[i] = 3.0 * g_0_0_xyyzzzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzzz_1[i] * wa_y[i];

        g_y_0_xyyzzzzzz_0[i] = 2.0 * g_0_0_xyzzzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzzz_1[i] * wa_y[i];

        g_y_0_xyzzzzzzz_0[i] = g_0_0_xzzzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzzz_1[i] * wa_y[i];

        g_y_0_xzzzzzzzz_0[i] = g_0_0_xzzzzzzzz_1[i] * wa_y[i];

        g_y_0_yyyyyyyyy_0[i] = 9.0 * g_0_0_yyyyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyyyy_1[i] * wa_y[i];

        g_y_0_yyyyyyyyz_0[i] = 8.0 * g_0_0_yyyyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyyyyz_1[i] * wa_y[i];

        g_y_0_yyyyyyyzz_0[i] = 7.0 * g_0_0_yyyyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyyyyzz_1[i] * wa_y[i];

        g_y_0_yyyyyyzzz_0[i] = 6.0 * g_0_0_yyyyyzzz_1[i] * fi_acd_0 + g_0_0_yyyyyyzzz_1[i] * wa_y[i];

        g_y_0_yyyyyzzzz_0[i] = 5.0 * g_0_0_yyyyzzzz_1[i] * fi_acd_0 + g_0_0_yyyyyzzzz_1[i] * wa_y[i];

        g_y_0_yyyyzzzzz_0[i] = 4.0 * g_0_0_yyyzzzzz_1[i] * fi_acd_0 + g_0_0_yyyyzzzzz_1[i] * wa_y[i];

        g_y_0_yyyzzzzzz_0[i] = 3.0 * g_0_0_yyzzzzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzzzz_1[i] * wa_y[i];

        g_y_0_yyzzzzzzz_0[i] = 2.0 * g_0_0_yzzzzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzzzz_1[i] * wa_y[i];

        g_y_0_yzzzzzzzz_0[i] = g_0_0_zzzzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzzzz_1[i] * wa_y[i];

        g_y_0_zzzzzzzzz_0[i] = g_0_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 110-165 components of targeted buffer : PSM

    auto g_z_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_psm + 110);

    auto g_z_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_psm + 111);

    auto g_z_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_psm + 112);

    auto g_z_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_psm + 113);

    auto g_z_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_psm + 114);

    auto g_z_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_psm + 115);

    auto g_z_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_psm + 116);

    auto g_z_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_psm + 117);

    auto g_z_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_psm + 118);

    auto g_z_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_psm + 119);

    auto g_z_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_psm + 120);

    auto g_z_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_psm + 121);

    auto g_z_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_psm + 122);

    auto g_z_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_psm + 123);

    auto g_z_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_psm + 124);

    auto g_z_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_psm + 125);

    auto g_z_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_psm + 126);

    auto g_z_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_psm + 127);

    auto g_z_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_psm + 128);

    auto g_z_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_psm + 129);

    auto g_z_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_psm + 130);

    auto g_z_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 131);

    auto g_z_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 132);

    auto g_z_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 133);

    auto g_z_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 134);

    auto g_z_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 135);

    auto g_z_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 136);

    auto g_z_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 137);

    auto g_z_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 138);

    auto g_z_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 139);

    auto g_z_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 140);

    auto g_z_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 141);

    auto g_z_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 142);

    auto g_z_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 143);

    auto g_z_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 144);

    auto g_z_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 145);

    auto g_z_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 146);

    auto g_z_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 147);

    auto g_z_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 148);

    auto g_z_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 149);

    auto g_z_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 150);

    auto g_z_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 151);

    auto g_z_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 152);

    auto g_z_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 153);

    auto g_z_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 154);

    auto g_z_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_psm + 155);

    auto g_z_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_psm + 156);

    auto g_z_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_psm + 157);

    auto g_z_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_psm + 158);

    auto g_z_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_psm + 159);

    auto g_z_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_psm + 160);

    auto g_z_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 161);

    auto g_z_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 162);

    auto g_z_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 163);

    auto g_z_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_psm + 164);

    #pragma omp simd aligned(g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxxx_1, g_0_0_xxxxxxxxy_1, g_0_0_xxxxxxxxz_1, g_0_0_xxxxxxxy_1, g_0_0_xxxxxxxyy_1, g_0_0_xxxxxxxyz_1, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxxzz_1, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyyy_1, g_0_0_xxxxxxyyz_1, g_0_0_xxxxxxyz_1, g_0_0_xxxxxxyzz_1, g_0_0_xxxxxxzz_1, g_0_0_xxxxxxzzz_1, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyyy_1, g_0_0_xxxxxyyyz_1, g_0_0_xxxxxyyz_1, g_0_0_xxxxxyyzz_1, g_0_0_xxxxxyzz_1, g_0_0_xxxxxyzzz_1, g_0_0_xxxxxzzz_1, g_0_0_xxxxxzzzz_1, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyyy_1, g_0_0_xxxxyyyyz_1, g_0_0_xxxxyyyz_1, g_0_0_xxxxyyyzz_1, g_0_0_xxxxyyzz_1, g_0_0_xxxxyyzzz_1, g_0_0_xxxxyzzz_1, g_0_0_xxxxyzzzz_1, g_0_0_xxxxzzzz_1, g_0_0_xxxxzzzzz_1, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyyy_1, g_0_0_xxxyyyyyz_1, g_0_0_xxxyyyyz_1, g_0_0_xxxyyyyzz_1, g_0_0_xxxyyyzz_1, g_0_0_xxxyyyzzz_1, g_0_0_xxxyyzzz_1, g_0_0_xxxyyzzzz_1, g_0_0_xxxyzzzz_1, g_0_0_xxxyzzzzz_1, g_0_0_xxxzzzzz_1, g_0_0_xxxzzzzzz_1, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyyy_1, g_0_0_xxyyyyyyz_1, g_0_0_xxyyyyyz_1, g_0_0_xxyyyyyzz_1, g_0_0_xxyyyyzz_1, g_0_0_xxyyyyzzz_1, g_0_0_xxyyyzzz_1, g_0_0_xxyyyzzzz_1, g_0_0_xxyyzzzz_1, g_0_0_xxyyzzzzz_1, g_0_0_xxyzzzzz_1, g_0_0_xxyzzzzzz_1, g_0_0_xxzzzzzz_1, g_0_0_xxzzzzzzz_1, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyyy_1, g_0_0_xyyyyyyyz_1, g_0_0_xyyyyyyz_1, g_0_0_xyyyyyyzz_1, g_0_0_xyyyyyzz_1, g_0_0_xyyyyyzzz_1, g_0_0_xyyyyzzz_1, g_0_0_xyyyyzzzz_1, g_0_0_xyyyzzzz_1, g_0_0_xyyyzzzzz_1, g_0_0_xyyzzzzz_1, g_0_0_xyyzzzzzz_1, g_0_0_xyzzzzzz_1, g_0_0_xyzzzzzzz_1, g_0_0_xzzzzzzz_1, g_0_0_xzzzzzzzz_1, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyyy_1, g_0_0_yyyyyyyyz_1, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyyzz_1, g_0_0_yyyyyyzz_1, g_0_0_yyyyyyzzz_1, g_0_0_yyyyyzzz_1, g_0_0_yyyyyzzzz_1, g_0_0_yyyyzzzz_1, g_0_0_yyyyzzzzz_1, g_0_0_yyyzzzzz_1, g_0_0_yyyzzzzzz_1, g_0_0_yyzzzzzz_1, g_0_0_yyzzzzzzz_1, g_0_0_yzzzzzzz_1, g_0_0_yzzzzzzzz_1, g_0_0_zzzzzzzz_1, g_0_0_zzzzzzzzz_1, g_z_0_xxxxxxxxx_0, g_z_0_xxxxxxxxy_0, g_z_0_xxxxxxxxz_0, g_z_0_xxxxxxxyy_0, g_z_0_xxxxxxxyz_0, g_z_0_xxxxxxxzz_0, g_z_0_xxxxxxyyy_0, g_z_0_xxxxxxyyz_0, g_z_0_xxxxxxyzz_0, g_z_0_xxxxxxzzz_0, g_z_0_xxxxxyyyy_0, g_z_0_xxxxxyyyz_0, g_z_0_xxxxxyyzz_0, g_z_0_xxxxxyzzz_0, g_z_0_xxxxxzzzz_0, g_z_0_xxxxyyyyy_0, g_z_0_xxxxyyyyz_0, g_z_0_xxxxyyyzz_0, g_z_0_xxxxyyzzz_0, g_z_0_xxxxyzzzz_0, g_z_0_xxxxzzzzz_0, g_z_0_xxxyyyyyy_0, g_z_0_xxxyyyyyz_0, g_z_0_xxxyyyyzz_0, g_z_0_xxxyyyzzz_0, g_z_0_xxxyyzzzz_0, g_z_0_xxxyzzzzz_0, g_z_0_xxxzzzzzz_0, g_z_0_xxyyyyyyy_0, g_z_0_xxyyyyyyz_0, g_z_0_xxyyyyyzz_0, g_z_0_xxyyyyzzz_0, g_z_0_xxyyyzzzz_0, g_z_0_xxyyzzzzz_0, g_z_0_xxyzzzzzz_0, g_z_0_xxzzzzzzz_0, g_z_0_xyyyyyyyy_0, g_z_0_xyyyyyyyz_0, g_z_0_xyyyyyyzz_0, g_z_0_xyyyyyzzz_0, g_z_0_xyyyyzzzz_0, g_z_0_xyyyzzzzz_0, g_z_0_xyyzzzzzz_0, g_z_0_xyzzzzzzz_0, g_z_0_xzzzzzzzz_0, g_z_0_yyyyyyyyy_0, g_z_0_yyyyyyyyz_0, g_z_0_yyyyyyyzz_0, g_z_0_yyyyyyzzz_0, g_z_0_yyyyyzzzz_0, g_z_0_yyyyzzzzz_0, g_z_0_yyyzzzzzz_0, g_z_0_yyzzzzzzz_0, g_z_0_yzzzzzzzz_0, g_z_0_zzzzzzzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxxxxxxx_0[i] = g_0_0_xxxxxxxxx_1[i] * wa_z[i];

        g_z_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxxy_1[i] * wa_z[i];

        g_z_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxx_1[i] * fi_acd_0 + g_0_0_xxxxxxxxz_1[i] * wa_z[i];

        g_z_0_xxxxxxxyy_0[i] = g_0_0_xxxxxxxyy_1[i] * wa_z[i];

        g_z_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxy_1[i] * fi_acd_0 + g_0_0_xxxxxxxyz_1[i] * wa_z[i];

        g_z_0_xxxxxxxzz_0[i] = 2.0 * g_0_0_xxxxxxxz_1[i] * fi_acd_0 + g_0_0_xxxxxxxzz_1[i] * wa_z[i];

        g_z_0_xxxxxxyyy_0[i] = g_0_0_xxxxxxyyy_1[i] * wa_z[i];

        g_z_0_xxxxxxyyz_0[i] = g_0_0_xxxxxxyy_1[i] * fi_acd_0 + g_0_0_xxxxxxyyz_1[i] * wa_z[i];

        g_z_0_xxxxxxyzz_0[i] = 2.0 * g_0_0_xxxxxxyz_1[i] * fi_acd_0 + g_0_0_xxxxxxyzz_1[i] * wa_z[i];

        g_z_0_xxxxxxzzz_0[i] = 3.0 * g_0_0_xxxxxxzz_1[i] * fi_acd_0 + g_0_0_xxxxxxzzz_1[i] * wa_z[i];

        g_z_0_xxxxxyyyy_0[i] = g_0_0_xxxxxyyyy_1[i] * wa_z[i];

        g_z_0_xxxxxyyyz_0[i] = g_0_0_xxxxxyyy_1[i] * fi_acd_0 + g_0_0_xxxxxyyyz_1[i] * wa_z[i];

        g_z_0_xxxxxyyzz_0[i] = 2.0 * g_0_0_xxxxxyyz_1[i] * fi_acd_0 + g_0_0_xxxxxyyzz_1[i] * wa_z[i];

        g_z_0_xxxxxyzzz_0[i] = 3.0 * g_0_0_xxxxxyzz_1[i] * fi_acd_0 + g_0_0_xxxxxyzzz_1[i] * wa_z[i];

        g_z_0_xxxxxzzzz_0[i] = 4.0 * g_0_0_xxxxxzzz_1[i] * fi_acd_0 + g_0_0_xxxxxzzzz_1[i] * wa_z[i];

        g_z_0_xxxxyyyyy_0[i] = g_0_0_xxxxyyyyy_1[i] * wa_z[i];

        g_z_0_xxxxyyyyz_0[i] = g_0_0_xxxxyyyy_1[i] * fi_acd_0 + g_0_0_xxxxyyyyz_1[i] * wa_z[i];

        g_z_0_xxxxyyyzz_0[i] = 2.0 * g_0_0_xxxxyyyz_1[i] * fi_acd_0 + g_0_0_xxxxyyyzz_1[i] * wa_z[i];

        g_z_0_xxxxyyzzz_0[i] = 3.0 * g_0_0_xxxxyyzz_1[i] * fi_acd_0 + g_0_0_xxxxyyzzz_1[i] * wa_z[i];

        g_z_0_xxxxyzzzz_0[i] = 4.0 * g_0_0_xxxxyzzz_1[i] * fi_acd_0 + g_0_0_xxxxyzzzz_1[i] * wa_z[i];

        g_z_0_xxxxzzzzz_0[i] = 5.0 * g_0_0_xxxxzzzz_1[i] * fi_acd_0 + g_0_0_xxxxzzzzz_1[i] * wa_z[i];

        g_z_0_xxxyyyyyy_0[i] = g_0_0_xxxyyyyyy_1[i] * wa_z[i];

        g_z_0_xxxyyyyyz_0[i] = g_0_0_xxxyyyyy_1[i] * fi_acd_0 + g_0_0_xxxyyyyyz_1[i] * wa_z[i];

        g_z_0_xxxyyyyzz_0[i] = 2.0 * g_0_0_xxxyyyyz_1[i] * fi_acd_0 + g_0_0_xxxyyyyzz_1[i] * wa_z[i];

        g_z_0_xxxyyyzzz_0[i] = 3.0 * g_0_0_xxxyyyzz_1[i] * fi_acd_0 + g_0_0_xxxyyyzzz_1[i] * wa_z[i];

        g_z_0_xxxyyzzzz_0[i] = 4.0 * g_0_0_xxxyyzzz_1[i] * fi_acd_0 + g_0_0_xxxyyzzzz_1[i] * wa_z[i];

        g_z_0_xxxyzzzzz_0[i] = 5.0 * g_0_0_xxxyzzzz_1[i] * fi_acd_0 + g_0_0_xxxyzzzzz_1[i] * wa_z[i];

        g_z_0_xxxzzzzzz_0[i] = 6.0 * g_0_0_xxxzzzzz_1[i] * fi_acd_0 + g_0_0_xxxzzzzzz_1[i] * wa_z[i];

        g_z_0_xxyyyyyyy_0[i] = g_0_0_xxyyyyyyy_1[i] * wa_z[i];

        g_z_0_xxyyyyyyz_0[i] = g_0_0_xxyyyyyy_1[i] * fi_acd_0 + g_0_0_xxyyyyyyz_1[i] * wa_z[i];

        g_z_0_xxyyyyyzz_0[i] = 2.0 * g_0_0_xxyyyyyz_1[i] * fi_acd_0 + g_0_0_xxyyyyyzz_1[i] * wa_z[i];

        g_z_0_xxyyyyzzz_0[i] = 3.0 * g_0_0_xxyyyyzz_1[i] * fi_acd_0 + g_0_0_xxyyyyzzz_1[i] * wa_z[i];

        g_z_0_xxyyyzzzz_0[i] = 4.0 * g_0_0_xxyyyzzz_1[i] * fi_acd_0 + g_0_0_xxyyyzzzz_1[i] * wa_z[i];

        g_z_0_xxyyzzzzz_0[i] = 5.0 * g_0_0_xxyyzzzz_1[i] * fi_acd_0 + g_0_0_xxyyzzzzz_1[i] * wa_z[i];

        g_z_0_xxyzzzzzz_0[i] = 6.0 * g_0_0_xxyzzzzz_1[i] * fi_acd_0 + g_0_0_xxyzzzzzz_1[i] * wa_z[i];

        g_z_0_xxzzzzzzz_0[i] = 7.0 * g_0_0_xxzzzzzz_1[i] * fi_acd_0 + g_0_0_xxzzzzzzz_1[i] * wa_z[i];

        g_z_0_xyyyyyyyy_0[i] = g_0_0_xyyyyyyyy_1[i] * wa_z[i];

        g_z_0_xyyyyyyyz_0[i] = g_0_0_xyyyyyyy_1[i] * fi_acd_0 + g_0_0_xyyyyyyyz_1[i] * wa_z[i];

        g_z_0_xyyyyyyzz_0[i] = 2.0 * g_0_0_xyyyyyyz_1[i] * fi_acd_0 + g_0_0_xyyyyyyzz_1[i] * wa_z[i];

        g_z_0_xyyyyyzzz_0[i] = 3.0 * g_0_0_xyyyyyzz_1[i] * fi_acd_0 + g_0_0_xyyyyyzzz_1[i] * wa_z[i];

        g_z_0_xyyyyzzzz_0[i] = 4.0 * g_0_0_xyyyyzzz_1[i] * fi_acd_0 + g_0_0_xyyyyzzzz_1[i] * wa_z[i];

        g_z_0_xyyyzzzzz_0[i] = 5.0 * g_0_0_xyyyzzzz_1[i] * fi_acd_0 + g_0_0_xyyyzzzzz_1[i] * wa_z[i];

        g_z_0_xyyzzzzzz_0[i] = 6.0 * g_0_0_xyyzzzzz_1[i] * fi_acd_0 + g_0_0_xyyzzzzzz_1[i] * wa_z[i];

        g_z_0_xyzzzzzzz_0[i] = 7.0 * g_0_0_xyzzzzzz_1[i] * fi_acd_0 + g_0_0_xyzzzzzzz_1[i] * wa_z[i];

        g_z_0_xzzzzzzzz_0[i] = 8.0 * g_0_0_xzzzzzzz_1[i] * fi_acd_0 + g_0_0_xzzzzzzzz_1[i] * wa_z[i];

        g_z_0_yyyyyyyyy_0[i] = g_0_0_yyyyyyyyy_1[i] * wa_z[i];

        g_z_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyy_1[i] * fi_acd_0 + g_0_0_yyyyyyyyz_1[i] * wa_z[i];

        g_z_0_yyyyyyyzz_0[i] = 2.0 * g_0_0_yyyyyyyz_1[i] * fi_acd_0 + g_0_0_yyyyyyyzz_1[i] * wa_z[i];

        g_z_0_yyyyyyzzz_0[i] = 3.0 * g_0_0_yyyyyyzz_1[i] * fi_acd_0 + g_0_0_yyyyyyzzz_1[i] * wa_z[i];

        g_z_0_yyyyyzzzz_0[i] = 4.0 * g_0_0_yyyyyzzz_1[i] * fi_acd_0 + g_0_0_yyyyyzzzz_1[i] * wa_z[i];

        g_z_0_yyyyzzzzz_0[i] = 5.0 * g_0_0_yyyyzzzz_1[i] * fi_acd_0 + g_0_0_yyyyzzzzz_1[i] * wa_z[i];

        g_z_0_yyyzzzzzz_0[i] = 6.0 * g_0_0_yyyzzzzz_1[i] * fi_acd_0 + g_0_0_yyyzzzzzz_1[i] * wa_z[i];

        g_z_0_yyzzzzzzz_0[i] = 7.0 * g_0_0_yyzzzzzz_1[i] * fi_acd_0 + g_0_0_yyzzzzzzz_1[i] * wa_z[i];

        g_z_0_yzzzzzzzz_0[i] = 8.0 * g_0_0_yzzzzzzz_1[i] * fi_acd_0 + g_0_0_yzzzzzzzz_1[i] * wa_z[i];

        g_z_0_zzzzzzzzz_0[i] = 9.0 * g_0_0_zzzzzzzz_1[i] * fi_acd_0 + g_0_0_zzzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace


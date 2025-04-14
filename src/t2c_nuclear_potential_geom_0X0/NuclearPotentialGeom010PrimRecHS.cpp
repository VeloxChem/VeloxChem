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

#include "NuclearPotentialGeom010PrimRecHS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_hs(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_hs,
                                        const size_t              idx_npot_geom_010_0_fs,
                                        const size_t              idx_npot_geom_010_1_fs,
                                        const size_t              idx_npot_1_gs,
                                        const size_t              idx_npot_geom_010_0_gs,
                                        const size_t              idx_npot_geom_010_1_gs,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs);

    auto ta1_x_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 1);

    auto ta1_x_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 2);

    auto ta1_x_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 6);

    auto ta1_x_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 8);

    auto ta1_x_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 9);

    auto ta1_y_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 10);

    auto ta1_y_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 13);

    auto ta1_y_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 15);

    auto ta1_y_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 16);

    auto ta1_y_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 17);

    auto ta1_y_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 19);

    auto ta1_z_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 20);

    auto ta1_z_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 23);

    auto ta1_z_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 25);

    auto ta1_z_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 26);

    auto ta1_z_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 28);

    auto ta1_z_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 29);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs);

    auto ta1_x_xxy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 1);

    auto ta1_x_xxz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 2);

    auto ta1_x_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 6);

    auto ta1_x_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 8);

    auto ta1_x_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 9);

    auto ta1_y_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 10);

    auto ta1_y_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 13);

    auto ta1_y_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 15);

    auto ta1_y_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 16);

    auto ta1_y_yyz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 17);

    auto ta1_y_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 19);

    auto ta1_z_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 20);

    auto ta1_z_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 23);

    auto ta1_z_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 25);

    auto ta1_z_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 26);

    auto ta1_z_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 28);

    auto ta1_z_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 29);

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_1 = pbuffer.data(idx_npot_1_gs);

    auto ta_xxyy_0_1 = pbuffer.data(idx_npot_1_gs + 3);

    auto ta_xxzz_0_1 = pbuffer.data(idx_npot_1_gs + 5);

    auto ta_yyyy_0_1 = pbuffer.data(idx_npot_1_gs + 10);

    auto ta_yyzz_0_1 = pbuffer.data(idx_npot_1_gs + 12);

    auto ta_zzzz_0_1 = pbuffer.data(idx_npot_1_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs);

    auto ta1_x_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 1);

    auto ta1_x_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 2);

    auto ta1_x_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 3);

    auto ta1_x_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 5);

    auto ta1_x_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 6);

    auto ta1_x_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 9);

    auto ta1_x_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 10);

    auto ta1_x_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 12);

    auto ta1_x_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 13);

    auto ta1_x_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 14);

    auto ta1_y_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 15);

    auto ta1_y_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 16);

    auto ta1_y_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 18);

    auto ta1_y_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 20);

    auto ta1_y_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 21);

    auto ta1_y_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 24);

    auto ta1_y_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 25);

    auto ta1_y_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 26);

    auto ta1_y_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 27);

    auto ta1_y_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 28);

    auto ta1_y_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 29);

    auto ta1_z_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 30);

    auto ta1_z_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 32);

    auto ta1_z_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 33);

    auto ta1_z_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 35);

    auto ta1_z_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 36);

    auto ta1_z_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 39);

    auto ta1_z_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 40);

    auto ta1_z_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 41);

    auto ta1_z_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 42);

    auto ta1_z_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 43);

    auto ta1_z_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 44);

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs);

    auto ta1_x_xxxy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 1);

    auto ta1_x_xxxz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 2);

    auto ta1_x_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 3);

    auto ta1_x_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 5);

    auto ta1_x_xyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 6);

    auto ta1_x_xzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 9);

    auto ta1_x_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 10);

    auto ta1_x_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 12);

    auto ta1_x_yzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 13);

    auto ta1_x_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 14);

    auto ta1_y_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 15);

    auto ta1_y_xxxy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 16);

    auto ta1_y_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 18);

    auto ta1_y_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 20);

    auto ta1_y_xyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 21);

    auto ta1_y_xzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 24);

    auto ta1_y_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 25);

    auto ta1_y_yyyz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 26);

    auto ta1_y_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 27);

    auto ta1_y_yzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 28);

    auto ta1_y_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 29);

    auto ta1_z_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 30);

    auto ta1_z_xxxz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 32);

    auto ta1_z_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 33);

    auto ta1_z_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 35);

    auto ta1_z_xyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 36);

    auto ta1_z_xzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 39);

    auto ta1_z_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 40);

    auto ta1_z_yyyz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 41);

    auto ta1_z_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 42);

    auto ta1_z_yzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 43);

    auto ta1_z_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 44);

    // Set up components of targeted buffer : HS

    auto ta1_x_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs);

    auto ta1_x_xxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 1);

    auto ta1_x_xxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 2);

    auto ta1_x_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 3);

    auto ta1_x_xxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 4);

    auto ta1_x_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 5);

    auto ta1_x_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 6);

    auto ta1_x_xxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 7);

    auto ta1_x_xxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 8);

    auto ta1_x_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 9);

    auto ta1_x_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 10);

    auto ta1_x_xyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 11);

    auto ta1_x_xyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 12);

    auto ta1_x_xyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 13);

    auto ta1_x_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 14);

    auto ta1_x_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 15);

    auto ta1_x_yyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 16);

    auto ta1_x_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 17);

    auto ta1_x_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 18);

    auto ta1_x_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 19);

    auto ta1_x_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 20);

    auto ta1_y_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 21);

    auto ta1_y_xxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 22);

    auto ta1_y_xxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 23);

    auto ta1_y_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 24);

    auto ta1_y_xxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 25);

    auto ta1_y_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 26);

    auto ta1_y_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 27);

    auto ta1_y_xxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 28);

    auto ta1_y_xxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 29);

    auto ta1_y_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 30);

    auto ta1_y_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 31);

    auto ta1_y_xyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 32);

    auto ta1_y_xyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 33);

    auto ta1_y_xyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 34);

    auto ta1_y_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 35);

    auto ta1_y_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 36);

    auto ta1_y_yyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 37);

    auto ta1_y_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 38);

    auto ta1_y_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 39);

    auto ta1_y_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 40);

    auto ta1_y_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 41);

    auto ta1_z_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 42);

    auto ta1_z_xxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 43);

    auto ta1_z_xxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 44);

    auto ta1_z_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 45);

    auto ta1_z_xxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 46);

    auto ta1_z_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 47);

    auto ta1_z_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 48);

    auto ta1_z_xxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 49);

    auto ta1_z_xxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 50);

    auto ta1_z_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 51);

    auto ta1_z_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 52);

    auto ta1_z_xyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 53);

    auto ta1_z_xyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 54);

    auto ta1_z_xyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 55);

    auto ta1_z_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 56);

    auto ta1_z_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 57);

    auto ta1_z_yyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 58);

    auto ta1_z_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 59);

    auto ta1_z_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 60);

    auto ta1_z_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 61);

    auto ta1_z_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 62);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxx_0_0,   \
                             ta1_x_xxx_0_1,   \
                             ta1_x_xxxx_0_0,  \
                             ta1_x_xxxx_0_1,  \
                             ta1_x_xxxxx_0_0, \
                             ta1_x_xxxxy_0_0, \
                             ta1_x_xxxxz_0_0, \
                             ta1_x_xxxy_0_0,  \
                             ta1_x_xxxy_0_1,  \
                             ta1_x_xxxyy_0_0, \
                             ta1_x_xxxyz_0_0, \
                             ta1_x_xxxz_0_0,  \
                             ta1_x_xxxz_0_1,  \
                             ta1_x_xxxzz_0_0, \
                             ta1_x_xxy_0_0,   \
                             ta1_x_xxy_0_1,   \
                             ta1_x_xxyy_0_0,  \
                             ta1_x_xxyy_0_1,  \
                             ta1_x_xxyyy_0_0, \
                             ta1_x_xxyyz_0_0, \
                             ta1_x_xxyzz_0_0, \
                             ta1_x_xxz_0_0,   \
                             ta1_x_xxz_0_1,   \
                             ta1_x_xxzz_0_0,  \
                             ta1_x_xxzz_0_1,  \
                             ta1_x_xxzzz_0_0, \
                             ta1_x_xyyy_0_0,  \
                             ta1_x_xyyy_0_1,  \
                             ta1_x_xyyyy_0_0, \
                             ta1_x_xyyyz_0_0, \
                             ta1_x_xyyzz_0_0, \
                             ta1_x_xyzzz_0_0, \
                             ta1_x_xzzz_0_0,  \
                             ta1_x_xzzz_0_1,  \
                             ta1_x_xzzzz_0_0, \
                             ta1_x_yyy_0_0,   \
                             ta1_x_yyy_0_1,   \
                             ta1_x_yyyy_0_0,  \
                             ta1_x_yyyy_0_1,  \
                             ta1_x_yyyyy_0_0, \
                             ta1_x_yyyyz_0_0, \
                             ta1_x_yyyzz_0_0, \
                             ta1_x_yyzz_0_0,  \
                             ta1_x_yyzz_0_1,  \
                             ta1_x_yyzzz_0_0, \
                             ta1_x_yzz_0_0,   \
                             ta1_x_yzz_0_1,   \
                             ta1_x_yzzz_0_0,  \
                             ta1_x_yzzz_0_1,  \
                             ta1_x_yzzzz_0_0, \
                             ta1_x_zzz_0_0,   \
                             ta1_x_zzz_0_1,   \
                             ta1_x_zzzz_0_0,  \
                             ta1_x_zzzz_0_1,  \
                             ta1_x_zzzzz_0_0, \
                             ta1_y_xxx_0_0,   \
                             ta1_y_xxx_0_1,   \
                             ta1_y_xxxx_0_0,  \
                             ta1_y_xxxx_0_1,  \
                             ta1_y_xxxxx_0_0, \
                             ta1_y_xxxxy_0_0, \
                             ta1_y_xxxxz_0_0, \
                             ta1_y_xxxy_0_0,  \
                             ta1_y_xxxy_0_1,  \
                             ta1_y_xxxyy_0_0, \
                             ta1_y_xxxyz_0_0, \
                             ta1_y_xxxzz_0_0, \
                             ta1_y_xxyy_0_0,  \
                             ta1_y_xxyy_0_1,  \
                             ta1_y_xxyyy_0_0, \
                             ta1_y_xxyyz_0_0, \
                             ta1_y_xxyzz_0_0, \
                             ta1_y_xxzz_0_0,  \
                             ta1_y_xxzz_0_1,  \
                             ta1_y_xxzzz_0_0, \
                             ta1_y_xyy_0_0,   \
                             ta1_y_xyy_0_1,   \
                             ta1_y_xyyy_0_0,  \
                             ta1_y_xyyy_0_1,  \
                             ta1_y_xyyyy_0_0, \
                             ta1_y_xyyyz_0_0, \
                             ta1_y_xyyzz_0_0, \
                             ta1_y_xyzzz_0_0, \
                             ta1_y_xzz_0_0,   \
                             ta1_y_xzz_0_1,   \
                             ta1_y_xzzz_0_0,  \
                             ta1_y_xzzz_0_1,  \
                             ta1_y_xzzzz_0_0, \
                             ta1_y_yyy_0_0,   \
                             ta1_y_yyy_0_1,   \
                             ta1_y_yyyy_0_0,  \
                             ta1_y_yyyy_0_1,  \
                             ta1_y_yyyyy_0_0, \
                             ta1_y_yyyyz_0_0, \
                             ta1_y_yyyz_0_0,  \
                             ta1_y_yyyz_0_1,  \
                             ta1_y_yyyzz_0_0, \
                             ta1_y_yyz_0_0,   \
                             ta1_y_yyz_0_1,   \
                             ta1_y_yyzz_0_0,  \
                             ta1_y_yyzz_0_1,  \
                             ta1_y_yyzzz_0_0, \
                             ta1_y_yzzz_0_0,  \
                             ta1_y_yzzz_0_1,  \
                             ta1_y_yzzzz_0_0, \
                             ta1_y_zzz_0_0,   \
                             ta1_y_zzz_0_1,   \
                             ta1_y_zzzz_0_0,  \
                             ta1_y_zzzz_0_1,  \
                             ta1_y_zzzzz_0_0, \
                             ta1_z_xxx_0_0,   \
                             ta1_z_xxx_0_1,   \
                             ta1_z_xxxx_0_0,  \
                             ta1_z_xxxx_0_1,  \
                             ta1_z_xxxxx_0_0, \
                             ta1_z_xxxxy_0_0, \
                             ta1_z_xxxxz_0_0, \
                             ta1_z_xxxyy_0_0, \
                             ta1_z_xxxyz_0_0, \
                             ta1_z_xxxz_0_0,  \
                             ta1_z_xxxz_0_1,  \
                             ta1_z_xxxzz_0_0, \
                             ta1_z_xxyy_0_0,  \
                             ta1_z_xxyy_0_1,  \
                             ta1_z_xxyyy_0_0, \
                             ta1_z_xxyyz_0_0, \
                             ta1_z_xxyzz_0_0, \
                             ta1_z_xxzz_0_0,  \
                             ta1_z_xxzz_0_1,  \
                             ta1_z_xxzzz_0_0, \
                             ta1_z_xyy_0_0,   \
                             ta1_z_xyy_0_1,   \
                             ta1_z_xyyy_0_0,  \
                             ta1_z_xyyy_0_1,  \
                             ta1_z_xyyyy_0_0, \
                             ta1_z_xyyyz_0_0, \
                             ta1_z_xyyzz_0_0, \
                             ta1_z_xyzzz_0_0, \
                             ta1_z_xzz_0_0,   \
                             ta1_z_xzz_0_1,   \
                             ta1_z_xzzz_0_0,  \
                             ta1_z_xzzz_0_1,  \
                             ta1_z_xzzzz_0_0, \
                             ta1_z_yyy_0_0,   \
                             ta1_z_yyy_0_1,   \
                             ta1_z_yyyy_0_0,  \
                             ta1_z_yyyy_0_1,  \
                             ta1_z_yyyyy_0_0, \
                             ta1_z_yyyyz_0_0, \
                             ta1_z_yyyz_0_0,  \
                             ta1_z_yyyz_0_1,  \
                             ta1_z_yyyzz_0_0, \
                             ta1_z_yyzz_0_0,  \
                             ta1_z_yyzz_0_1,  \
                             ta1_z_yyzzz_0_0, \
                             ta1_z_yzz_0_0,   \
                             ta1_z_yzz_0_1,   \
                             ta1_z_yzzz_0_0,  \
                             ta1_z_yzzz_0_1,  \
                             ta1_z_yzzzz_0_0, \
                             ta1_z_zzz_0_0,   \
                             ta1_z_zzz_0_1,   \
                             ta1_z_zzzz_0_0,  \
                             ta1_z_zzzz_0_1,  \
                             ta1_z_zzzzz_0_0, \
                             ta_xxxx_0_1,     \
                             ta_xxyy_0_1,     \
                             ta_xxzz_0_1,     \
                             ta_yyyy_0_1,     \
                             ta_yyzz_0_1,     \
                             ta_zzzz_0_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxx_0_0[i] = 4.0 * ta1_x_xxx_0_0[i] * fe_0 - 4.0 * ta1_x_xxx_0_1[i] * fe_0 + ta_xxxx_0_1[i] + ta1_x_xxxx_0_0[i] * pa_x[i] -
                             ta1_x_xxxx_0_1[i] * pc_x[i];

        ta1_x_xxxxy_0_0[i] = ta1_x_xxxx_0_0[i] * pa_y[i] - ta1_x_xxxx_0_1[i] * pc_y[i];

        ta1_x_xxxxz_0_0[i] = ta1_x_xxxx_0_0[i] * pa_z[i] - ta1_x_xxxx_0_1[i] * pc_z[i];

        ta1_x_xxxyy_0_0[i] = ta1_x_xxx_0_0[i] * fe_0 - ta1_x_xxx_0_1[i] * fe_0 + ta1_x_xxxy_0_0[i] * pa_y[i] - ta1_x_xxxy_0_1[i] * pc_y[i];

        ta1_x_xxxyz_0_0[i] = ta1_x_xxxz_0_0[i] * pa_y[i] - ta1_x_xxxz_0_1[i] * pc_y[i];

        ta1_x_xxxzz_0_0[i] = ta1_x_xxx_0_0[i] * fe_0 - ta1_x_xxx_0_1[i] * fe_0 + ta1_x_xxxz_0_0[i] * pa_z[i] - ta1_x_xxxz_0_1[i] * pc_z[i];

        ta1_x_xxyyy_0_0[i] =
            2.0 * ta1_x_xxy_0_0[i] * fe_0 - 2.0 * ta1_x_xxy_0_1[i] * fe_0 + ta1_x_xxyy_0_0[i] * pa_y[i] - ta1_x_xxyy_0_1[i] * pc_y[i];

        ta1_x_xxyyz_0_0[i] = ta1_x_xxyy_0_0[i] * pa_z[i] - ta1_x_xxyy_0_1[i] * pc_z[i];

        ta1_x_xxyzz_0_0[i] = ta1_x_xxzz_0_0[i] * pa_y[i] - ta1_x_xxzz_0_1[i] * pc_y[i];

        ta1_x_xxzzz_0_0[i] =
            2.0 * ta1_x_xxz_0_0[i] * fe_0 - 2.0 * ta1_x_xxz_0_1[i] * fe_0 + ta1_x_xxzz_0_0[i] * pa_z[i] - ta1_x_xxzz_0_1[i] * pc_z[i];

        ta1_x_xyyyy_0_0[i] = ta_yyyy_0_1[i] + ta1_x_yyyy_0_0[i] * pa_x[i] - ta1_x_yyyy_0_1[i] * pc_x[i];

        ta1_x_xyyyz_0_0[i] = ta1_x_xyyy_0_0[i] * pa_z[i] - ta1_x_xyyy_0_1[i] * pc_z[i];

        ta1_x_xyyzz_0_0[i] = ta_yyzz_0_1[i] + ta1_x_yyzz_0_0[i] * pa_x[i] - ta1_x_yyzz_0_1[i] * pc_x[i];

        ta1_x_xyzzz_0_0[i] = ta1_x_xzzz_0_0[i] * pa_y[i] - ta1_x_xzzz_0_1[i] * pc_y[i];

        ta1_x_xzzzz_0_0[i] = ta_zzzz_0_1[i] + ta1_x_zzzz_0_0[i] * pa_x[i] - ta1_x_zzzz_0_1[i] * pc_x[i];

        ta1_x_yyyyy_0_0[i] =
            4.0 * ta1_x_yyy_0_0[i] * fe_0 - 4.0 * ta1_x_yyy_0_1[i] * fe_0 + ta1_x_yyyy_0_0[i] * pa_y[i] - ta1_x_yyyy_0_1[i] * pc_y[i];

        ta1_x_yyyyz_0_0[i] = ta1_x_yyyy_0_0[i] * pa_z[i] - ta1_x_yyyy_0_1[i] * pc_z[i];

        ta1_x_yyyzz_0_0[i] =
            2.0 * ta1_x_yzz_0_0[i] * fe_0 - 2.0 * ta1_x_yzz_0_1[i] * fe_0 + ta1_x_yyzz_0_0[i] * pa_y[i] - ta1_x_yyzz_0_1[i] * pc_y[i];

        ta1_x_yyzzz_0_0[i] = ta1_x_zzz_0_0[i] * fe_0 - ta1_x_zzz_0_1[i] * fe_0 + ta1_x_yzzz_0_0[i] * pa_y[i] - ta1_x_yzzz_0_1[i] * pc_y[i];

        ta1_x_yzzzz_0_0[i] = ta1_x_zzzz_0_0[i] * pa_y[i] - ta1_x_zzzz_0_1[i] * pc_y[i];

        ta1_x_zzzzz_0_0[i] =
            4.0 * ta1_x_zzz_0_0[i] * fe_0 - 4.0 * ta1_x_zzz_0_1[i] * fe_0 + ta1_x_zzzz_0_0[i] * pa_z[i] - ta1_x_zzzz_0_1[i] * pc_z[i];

        ta1_y_xxxxx_0_0[i] =
            4.0 * ta1_y_xxx_0_0[i] * fe_0 - 4.0 * ta1_y_xxx_0_1[i] * fe_0 + ta1_y_xxxx_0_0[i] * pa_x[i] - ta1_y_xxxx_0_1[i] * pc_x[i];

        ta1_y_xxxxy_0_0[i] = ta_xxxx_0_1[i] + ta1_y_xxxx_0_0[i] * pa_y[i] - ta1_y_xxxx_0_1[i] * pc_y[i];

        ta1_y_xxxxz_0_0[i] = ta1_y_xxxx_0_0[i] * pa_z[i] - ta1_y_xxxx_0_1[i] * pc_z[i];

        ta1_y_xxxyy_0_0[i] =
            2.0 * ta1_y_xyy_0_0[i] * fe_0 - 2.0 * ta1_y_xyy_0_1[i] * fe_0 + ta1_y_xxyy_0_0[i] * pa_x[i] - ta1_y_xxyy_0_1[i] * pc_x[i];

        ta1_y_xxxyz_0_0[i] = ta1_y_xxxy_0_0[i] * pa_z[i] - ta1_y_xxxy_0_1[i] * pc_z[i];

        ta1_y_xxxzz_0_0[i] =
            2.0 * ta1_y_xzz_0_0[i] * fe_0 - 2.0 * ta1_y_xzz_0_1[i] * fe_0 + ta1_y_xxzz_0_0[i] * pa_x[i] - ta1_y_xxzz_0_1[i] * pc_x[i];

        ta1_y_xxyyy_0_0[i] = ta1_y_yyy_0_0[i] * fe_0 - ta1_y_yyy_0_1[i] * fe_0 + ta1_y_xyyy_0_0[i] * pa_x[i] - ta1_y_xyyy_0_1[i] * pc_x[i];

        ta1_y_xxyyz_0_0[i] = ta1_y_xxyy_0_0[i] * pa_z[i] - ta1_y_xxyy_0_1[i] * pc_z[i];

        ta1_y_xxyzz_0_0[i] = ta_xxzz_0_1[i] + ta1_y_xxzz_0_0[i] * pa_y[i] - ta1_y_xxzz_0_1[i] * pc_y[i];

        ta1_y_xxzzz_0_0[i] = ta1_y_zzz_0_0[i] * fe_0 - ta1_y_zzz_0_1[i] * fe_0 + ta1_y_xzzz_0_0[i] * pa_x[i] - ta1_y_xzzz_0_1[i] * pc_x[i];

        ta1_y_xyyyy_0_0[i] = ta1_y_yyyy_0_0[i] * pa_x[i] - ta1_y_yyyy_0_1[i] * pc_x[i];

        ta1_y_xyyyz_0_0[i] = ta1_y_yyyz_0_0[i] * pa_x[i] - ta1_y_yyyz_0_1[i] * pc_x[i];

        ta1_y_xyyzz_0_0[i] = ta1_y_yyzz_0_0[i] * pa_x[i] - ta1_y_yyzz_0_1[i] * pc_x[i];

        ta1_y_xyzzz_0_0[i] = ta1_y_yzzz_0_0[i] * pa_x[i] - ta1_y_yzzz_0_1[i] * pc_x[i];

        ta1_y_xzzzz_0_0[i] = ta1_y_zzzz_0_0[i] * pa_x[i] - ta1_y_zzzz_0_1[i] * pc_x[i];

        ta1_y_yyyyy_0_0[i] = 4.0 * ta1_y_yyy_0_0[i] * fe_0 - 4.0 * ta1_y_yyy_0_1[i] * fe_0 + ta_yyyy_0_1[i] + ta1_y_yyyy_0_0[i] * pa_y[i] -
                             ta1_y_yyyy_0_1[i] * pc_y[i];

        ta1_y_yyyyz_0_0[i] = ta1_y_yyyy_0_0[i] * pa_z[i] - ta1_y_yyyy_0_1[i] * pc_z[i];

        ta1_y_yyyzz_0_0[i] = ta1_y_yyy_0_0[i] * fe_0 - ta1_y_yyy_0_1[i] * fe_0 + ta1_y_yyyz_0_0[i] * pa_z[i] - ta1_y_yyyz_0_1[i] * pc_z[i];

        ta1_y_yyzzz_0_0[i] =
            2.0 * ta1_y_yyz_0_0[i] * fe_0 - 2.0 * ta1_y_yyz_0_1[i] * fe_0 + ta1_y_yyzz_0_0[i] * pa_z[i] - ta1_y_yyzz_0_1[i] * pc_z[i];

        ta1_y_yzzzz_0_0[i] = ta_zzzz_0_1[i] + ta1_y_zzzz_0_0[i] * pa_y[i] - ta1_y_zzzz_0_1[i] * pc_y[i];

        ta1_y_zzzzz_0_0[i] =
            4.0 * ta1_y_zzz_0_0[i] * fe_0 - 4.0 * ta1_y_zzz_0_1[i] * fe_0 + ta1_y_zzzz_0_0[i] * pa_z[i] - ta1_y_zzzz_0_1[i] * pc_z[i];

        ta1_z_xxxxx_0_0[i] =
            4.0 * ta1_z_xxx_0_0[i] * fe_0 - 4.0 * ta1_z_xxx_0_1[i] * fe_0 + ta1_z_xxxx_0_0[i] * pa_x[i] - ta1_z_xxxx_0_1[i] * pc_x[i];

        ta1_z_xxxxy_0_0[i] = ta1_z_xxxx_0_0[i] * pa_y[i] - ta1_z_xxxx_0_1[i] * pc_y[i];

        ta1_z_xxxxz_0_0[i] = ta_xxxx_0_1[i] + ta1_z_xxxx_0_0[i] * pa_z[i] - ta1_z_xxxx_0_1[i] * pc_z[i];

        ta1_z_xxxyy_0_0[i] =
            2.0 * ta1_z_xyy_0_0[i] * fe_0 - 2.0 * ta1_z_xyy_0_1[i] * fe_0 + ta1_z_xxyy_0_0[i] * pa_x[i] - ta1_z_xxyy_0_1[i] * pc_x[i];

        ta1_z_xxxyz_0_0[i] = ta1_z_xxxz_0_0[i] * pa_y[i] - ta1_z_xxxz_0_1[i] * pc_y[i];

        ta1_z_xxxzz_0_0[i] =
            2.0 * ta1_z_xzz_0_0[i] * fe_0 - 2.0 * ta1_z_xzz_0_1[i] * fe_0 + ta1_z_xxzz_0_0[i] * pa_x[i] - ta1_z_xxzz_0_1[i] * pc_x[i];

        ta1_z_xxyyy_0_0[i] = ta1_z_yyy_0_0[i] * fe_0 - ta1_z_yyy_0_1[i] * fe_0 + ta1_z_xyyy_0_0[i] * pa_x[i] - ta1_z_xyyy_0_1[i] * pc_x[i];

        ta1_z_xxyyz_0_0[i] = ta_xxyy_0_1[i] + ta1_z_xxyy_0_0[i] * pa_z[i] - ta1_z_xxyy_0_1[i] * pc_z[i];

        ta1_z_xxyzz_0_0[i] = ta1_z_xxzz_0_0[i] * pa_y[i] - ta1_z_xxzz_0_1[i] * pc_y[i];

        ta1_z_xxzzz_0_0[i] = ta1_z_zzz_0_0[i] * fe_0 - ta1_z_zzz_0_1[i] * fe_0 + ta1_z_xzzz_0_0[i] * pa_x[i] - ta1_z_xzzz_0_1[i] * pc_x[i];

        ta1_z_xyyyy_0_0[i] = ta1_z_yyyy_0_0[i] * pa_x[i] - ta1_z_yyyy_0_1[i] * pc_x[i];

        ta1_z_xyyyz_0_0[i] = ta1_z_yyyz_0_0[i] * pa_x[i] - ta1_z_yyyz_0_1[i] * pc_x[i];

        ta1_z_xyyzz_0_0[i] = ta1_z_yyzz_0_0[i] * pa_x[i] - ta1_z_yyzz_0_1[i] * pc_x[i];

        ta1_z_xyzzz_0_0[i] = ta1_z_yzzz_0_0[i] * pa_x[i] - ta1_z_yzzz_0_1[i] * pc_x[i];

        ta1_z_xzzzz_0_0[i] = ta1_z_zzzz_0_0[i] * pa_x[i] - ta1_z_zzzz_0_1[i] * pc_x[i];

        ta1_z_yyyyy_0_0[i] =
            4.0 * ta1_z_yyy_0_0[i] * fe_0 - 4.0 * ta1_z_yyy_0_1[i] * fe_0 + ta1_z_yyyy_0_0[i] * pa_y[i] - ta1_z_yyyy_0_1[i] * pc_y[i];

        ta1_z_yyyyz_0_0[i] = ta_yyyy_0_1[i] + ta1_z_yyyy_0_0[i] * pa_z[i] - ta1_z_yyyy_0_1[i] * pc_z[i];

        ta1_z_yyyzz_0_0[i] =
            2.0 * ta1_z_yzz_0_0[i] * fe_0 - 2.0 * ta1_z_yzz_0_1[i] * fe_0 + ta1_z_yyzz_0_0[i] * pa_y[i] - ta1_z_yyzz_0_1[i] * pc_y[i];

        ta1_z_yyzzz_0_0[i] = ta1_z_zzz_0_0[i] * fe_0 - ta1_z_zzz_0_1[i] * fe_0 + ta1_z_yzzz_0_0[i] * pa_y[i] - ta1_z_yzzz_0_1[i] * pc_y[i];

        ta1_z_yzzzz_0_0[i] = ta1_z_zzzz_0_0[i] * pa_y[i] - ta1_z_zzzz_0_1[i] * pc_y[i];

        ta1_z_zzzzz_0_0[i] = 4.0 * ta1_z_zzz_0_0[i] * fe_0 - 4.0 * ta1_z_zzz_0_1[i] * fe_0 + ta_zzzz_0_1[i] + ta1_z_zzzz_0_0[i] * pa_z[i] -
                             ta1_z_zzzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec

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

#include "ThreeCenterElectronRepulsionContrRecXGF.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xgf(CSimdArray<double>& cbuffer,
                                const size_t idx_xgf,
                                const size_t idx_xff,
                                const size_t idx_xfg,
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
        /// Set up components of auxilary buffer : SFF

        const auto ff_off = idx_xff + i * 100;

        auto g_xxx_xxx = cbuffer.data(ff_off + 0);

        auto g_xxx_xxy = cbuffer.data(ff_off + 1);

        auto g_xxx_xxz = cbuffer.data(ff_off + 2);

        auto g_xxx_xyy = cbuffer.data(ff_off + 3);

        auto g_xxx_xyz = cbuffer.data(ff_off + 4);

        auto g_xxx_xzz = cbuffer.data(ff_off + 5);

        auto g_xxx_yyy = cbuffer.data(ff_off + 6);

        auto g_xxx_yyz = cbuffer.data(ff_off + 7);

        auto g_xxx_yzz = cbuffer.data(ff_off + 8);

        auto g_xxx_zzz = cbuffer.data(ff_off + 9);

        auto g_xxy_xxx = cbuffer.data(ff_off + 10);

        auto g_xxy_xxy = cbuffer.data(ff_off + 11);

        auto g_xxy_xxz = cbuffer.data(ff_off + 12);

        auto g_xxy_xyy = cbuffer.data(ff_off + 13);

        auto g_xxy_xyz = cbuffer.data(ff_off + 14);

        auto g_xxy_xzz = cbuffer.data(ff_off + 15);

        auto g_xxy_yyy = cbuffer.data(ff_off + 16);

        auto g_xxy_yyz = cbuffer.data(ff_off + 17);

        auto g_xxy_yzz = cbuffer.data(ff_off + 18);

        auto g_xxy_zzz = cbuffer.data(ff_off + 19);

        auto g_xxz_xxx = cbuffer.data(ff_off + 20);

        auto g_xxz_xxy = cbuffer.data(ff_off + 21);

        auto g_xxz_xxz = cbuffer.data(ff_off + 22);

        auto g_xxz_xyy = cbuffer.data(ff_off + 23);

        auto g_xxz_xyz = cbuffer.data(ff_off + 24);

        auto g_xxz_xzz = cbuffer.data(ff_off + 25);

        auto g_xxz_yyy = cbuffer.data(ff_off + 26);

        auto g_xxz_yyz = cbuffer.data(ff_off + 27);

        auto g_xxz_yzz = cbuffer.data(ff_off + 28);

        auto g_xxz_zzz = cbuffer.data(ff_off + 29);

        auto g_xyy_xxx = cbuffer.data(ff_off + 30);

        auto g_xyy_xxy = cbuffer.data(ff_off + 31);

        auto g_xyy_xxz = cbuffer.data(ff_off + 32);

        auto g_xyy_xyy = cbuffer.data(ff_off + 33);

        auto g_xyy_xyz = cbuffer.data(ff_off + 34);

        auto g_xyy_xzz = cbuffer.data(ff_off + 35);

        auto g_xyy_yyy = cbuffer.data(ff_off + 36);

        auto g_xyy_yyz = cbuffer.data(ff_off + 37);

        auto g_xyy_yzz = cbuffer.data(ff_off + 38);

        auto g_xyy_zzz = cbuffer.data(ff_off + 39);

        auto g_xyz_xxx = cbuffer.data(ff_off + 40);

        auto g_xyz_xxy = cbuffer.data(ff_off + 41);

        auto g_xyz_xxz = cbuffer.data(ff_off + 42);

        auto g_xyz_xyy = cbuffer.data(ff_off + 43);

        auto g_xyz_xyz = cbuffer.data(ff_off + 44);

        auto g_xyz_xzz = cbuffer.data(ff_off + 45);

        auto g_xyz_yyy = cbuffer.data(ff_off + 46);

        auto g_xyz_yyz = cbuffer.data(ff_off + 47);

        auto g_xyz_yzz = cbuffer.data(ff_off + 48);

        auto g_xyz_zzz = cbuffer.data(ff_off + 49);

        auto g_xzz_xxx = cbuffer.data(ff_off + 50);

        auto g_xzz_xxy = cbuffer.data(ff_off + 51);

        auto g_xzz_xxz = cbuffer.data(ff_off + 52);

        auto g_xzz_xyy = cbuffer.data(ff_off + 53);

        auto g_xzz_xyz = cbuffer.data(ff_off + 54);

        auto g_xzz_xzz = cbuffer.data(ff_off + 55);

        auto g_xzz_yyy = cbuffer.data(ff_off + 56);

        auto g_xzz_yyz = cbuffer.data(ff_off + 57);

        auto g_xzz_yzz = cbuffer.data(ff_off + 58);

        auto g_xzz_zzz = cbuffer.data(ff_off + 59);

        auto g_yyy_xxx = cbuffer.data(ff_off + 60);

        auto g_yyy_xxy = cbuffer.data(ff_off + 61);

        auto g_yyy_xxz = cbuffer.data(ff_off + 62);

        auto g_yyy_xyy = cbuffer.data(ff_off + 63);

        auto g_yyy_xyz = cbuffer.data(ff_off + 64);

        auto g_yyy_xzz = cbuffer.data(ff_off + 65);

        auto g_yyy_yyy = cbuffer.data(ff_off + 66);

        auto g_yyy_yyz = cbuffer.data(ff_off + 67);

        auto g_yyy_yzz = cbuffer.data(ff_off + 68);

        auto g_yyy_zzz = cbuffer.data(ff_off + 69);

        auto g_yyz_xxx = cbuffer.data(ff_off + 70);

        auto g_yyz_xxy = cbuffer.data(ff_off + 71);

        auto g_yyz_xxz = cbuffer.data(ff_off + 72);

        auto g_yyz_xyy = cbuffer.data(ff_off + 73);

        auto g_yyz_xyz = cbuffer.data(ff_off + 74);

        auto g_yyz_xzz = cbuffer.data(ff_off + 75);

        auto g_yyz_yyy = cbuffer.data(ff_off + 76);

        auto g_yyz_yyz = cbuffer.data(ff_off + 77);

        auto g_yyz_yzz = cbuffer.data(ff_off + 78);

        auto g_yyz_zzz = cbuffer.data(ff_off + 79);

        auto g_yzz_xxx = cbuffer.data(ff_off + 80);

        auto g_yzz_xxy = cbuffer.data(ff_off + 81);

        auto g_yzz_xxz = cbuffer.data(ff_off + 82);

        auto g_yzz_xyy = cbuffer.data(ff_off + 83);

        auto g_yzz_xyz = cbuffer.data(ff_off + 84);

        auto g_yzz_xzz = cbuffer.data(ff_off + 85);

        auto g_yzz_yyy = cbuffer.data(ff_off + 86);

        auto g_yzz_yyz = cbuffer.data(ff_off + 87);

        auto g_yzz_yzz = cbuffer.data(ff_off + 88);

        auto g_yzz_zzz = cbuffer.data(ff_off + 89);

        auto g_zzz_xxx = cbuffer.data(ff_off + 90);

        auto g_zzz_xxy = cbuffer.data(ff_off + 91);

        auto g_zzz_xxz = cbuffer.data(ff_off + 92);

        auto g_zzz_xyy = cbuffer.data(ff_off + 93);

        auto g_zzz_xyz = cbuffer.data(ff_off + 94);

        auto g_zzz_xzz = cbuffer.data(ff_off + 95);

        auto g_zzz_yyy = cbuffer.data(ff_off + 96);

        auto g_zzz_yyz = cbuffer.data(ff_off + 97);

        auto g_zzz_yzz = cbuffer.data(ff_off + 98);

        auto g_zzz_zzz = cbuffer.data(ff_off + 99);

        /// Set up components of auxilary buffer : SFG

        const auto fg_off = idx_xfg + i * 150;

        auto g_xxx_xxxx = cbuffer.data(fg_off + 0);

        auto g_xxx_xxxy = cbuffer.data(fg_off + 1);

        auto g_xxx_xxxz = cbuffer.data(fg_off + 2);

        auto g_xxx_xxyy = cbuffer.data(fg_off + 3);

        auto g_xxx_xxyz = cbuffer.data(fg_off + 4);

        auto g_xxx_xxzz = cbuffer.data(fg_off + 5);

        auto g_xxx_xyyy = cbuffer.data(fg_off + 6);

        auto g_xxx_xyyz = cbuffer.data(fg_off + 7);

        auto g_xxx_xyzz = cbuffer.data(fg_off + 8);

        auto g_xxx_xzzz = cbuffer.data(fg_off + 9);

        auto g_xxy_xxxx = cbuffer.data(fg_off + 15);

        auto g_xxy_xxxy = cbuffer.data(fg_off + 16);

        auto g_xxy_xxxz = cbuffer.data(fg_off + 17);

        auto g_xxy_xxyy = cbuffer.data(fg_off + 18);

        auto g_xxy_xxyz = cbuffer.data(fg_off + 19);

        auto g_xxy_xxzz = cbuffer.data(fg_off + 20);

        auto g_xxy_xyyy = cbuffer.data(fg_off + 21);

        auto g_xxy_xyyz = cbuffer.data(fg_off + 22);

        auto g_xxy_xyzz = cbuffer.data(fg_off + 23);

        auto g_xxy_xzzz = cbuffer.data(fg_off + 24);

        auto g_xxz_xxxx = cbuffer.data(fg_off + 30);

        auto g_xxz_xxxy = cbuffer.data(fg_off + 31);

        auto g_xxz_xxxz = cbuffer.data(fg_off + 32);

        auto g_xxz_xxyy = cbuffer.data(fg_off + 33);

        auto g_xxz_xxyz = cbuffer.data(fg_off + 34);

        auto g_xxz_xxzz = cbuffer.data(fg_off + 35);

        auto g_xxz_xyyy = cbuffer.data(fg_off + 36);

        auto g_xxz_xyyz = cbuffer.data(fg_off + 37);

        auto g_xxz_xyzz = cbuffer.data(fg_off + 38);

        auto g_xxz_xzzz = cbuffer.data(fg_off + 39);

        auto g_xyy_xxxx = cbuffer.data(fg_off + 45);

        auto g_xyy_xxxy = cbuffer.data(fg_off + 46);

        auto g_xyy_xxxz = cbuffer.data(fg_off + 47);

        auto g_xyy_xxyy = cbuffer.data(fg_off + 48);

        auto g_xyy_xxyz = cbuffer.data(fg_off + 49);

        auto g_xyy_xxzz = cbuffer.data(fg_off + 50);

        auto g_xyy_xyyy = cbuffer.data(fg_off + 51);

        auto g_xyy_xyyz = cbuffer.data(fg_off + 52);

        auto g_xyy_xyzz = cbuffer.data(fg_off + 53);

        auto g_xyy_xzzz = cbuffer.data(fg_off + 54);

        auto g_xyz_xxxx = cbuffer.data(fg_off + 60);

        auto g_xyz_xxxy = cbuffer.data(fg_off + 61);

        auto g_xyz_xxxz = cbuffer.data(fg_off + 62);

        auto g_xyz_xxyy = cbuffer.data(fg_off + 63);

        auto g_xyz_xxyz = cbuffer.data(fg_off + 64);

        auto g_xyz_xxzz = cbuffer.data(fg_off + 65);

        auto g_xyz_xyyy = cbuffer.data(fg_off + 66);

        auto g_xyz_xyyz = cbuffer.data(fg_off + 67);

        auto g_xyz_xyzz = cbuffer.data(fg_off + 68);

        auto g_xyz_xzzz = cbuffer.data(fg_off + 69);

        auto g_xzz_xxxx = cbuffer.data(fg_off + 75);

        auto g_xzz_xxxy = cbuffer.data(fg_off + 76);

        auto g_xzz_xxxz = cbuffer.data(fg_off + 77);

        auto g_xzz_xxyy = cbuffer.data(fg_off + 78);

        auto g_xzz_xxyz = cbuffer.data(fg_off + 79);

        auto g_xzz_xxzz = cbuffer.data(fg_off + 80);

        auto g_xzz_xyyy = cbuffer.data(fg_off + 81);

        auto g_xzz_xyyz = cbuffer.data(fg_off + 82);

        auto g_xzz_xyzz = cbuffer.data(fg_off + 83);

        auto g_xzz_xzzz = cbuffer.data(fg_off + 84);

        auto g_yyy_xxxx = cbuffer.data(fg_off + 90);

        auto g_yyy_xxxy = cbuffer.data(fg_off + 91);

        auto g_yyy_xxxz = cbuffer.data(fg_off + 92);

        auto g_yyy_xxyy = cbuffer.data(fg_off + 93);

        auto g_yyy_xxyz = cbuffer.data(fg_off + 94);

        auto g_yyy_xxzz = cbuffer.data(fg_off + 95);

        auto g_yyy_xyyy = cbuffer.data(fg_off + 96);

        auto g_yyy_xyyz = cbuffer.data(fg_off + 97);

        auto g_yyy_xyzz = cbuffer.data(fg_off + 98);

        auto g_yyy_xzzz = cbuffer.data(fg_off + 99);

        auto g_yyy_yyyy = cbuffer.data(fg_off + 100);

        auto g_yyy_yyyz = cbuffer.data(fg_off + 101);

        auto g_yyy_yyzz = cbuffer.data(fg_off + 102);

        auto g_yyy_yzzz = cbuffer.data(fg_off + 103);

        auto g_yyz_xxxx = cbuffer.data(fg_off + 105);

        auto g_yyz_xxxy = cbuffer.data(fg_off + 106);

        auto g_yyz_xxxz = cbuffer.data(fg_off + 107);

        auto g_yyz_xxyy = cbuffer.data(fg_off + 108);

        auto g_yyz_xxyz = cbuffer.data(fg_off + 109);

        auto g_yyz_xxzz = cbuffer.data(fg_off + 110);

        auto g_yyz_xyyy = cbuffer.data(fg_off + 111);

        auto g_yyz_xyyz = cbuffer.data(fg_off + 112);

        auto g_yyz_xyzz = cbuffer.data(fg_off + 113);

        auto g_yyz_xzzz = cbuffer.data(fg_off + 114);

        auto g_yyz_yyyy = cbuffer.data(fg_off + 115);

        auto g_yyz_yyyz = cbuffer.data(fg_off + 116);

        auto g_yyz_yyzz = cbuffer.data(fg_off + 117);

        auto g_yyz_yzzz = cbuffer.data(fg_off + 118);

        auto g_yzz_xxxx = cbuffer.data(fg_off + 120);

        auto g_yzz_xxxy = cbuffer.data(fg_off + 121);

        auto g_yzz_xxxz = cbuffer.data(fg_off + 122);

        auto g_yzz_xxyy = cbuffer.data(fg_off + 123);

        auto g_yzz_xxyz = cbuffer.data(fg_off + 124);

        auto g_yzz_xxzz = cbuffer.data(fg_off + 125);

        auto g_yzz_xyyy = cbuffer.data(fg_off + 126);

        auto g_yzz_xyyz = cbuffer.data(fg_off + 127);

        auto g_yzz_xyzz = cbuffer.data(fg_off + 128);

        auto g_yzz_xzzz = cbuffer.data(fg_off + 129);

        auto g_yzz_yyyy = cbuffer.data(fg_off + 130);

        auto g_yzz_yyyz = cbuffer.data(fg_off + 131);

        auto g_yzz_yyzz = cbuffer.data(fg_off + 132);

        auto g_yzz_yzzz = cbuffer.data(fg_off + 133);

        auto g_zzz_xxxx = cbuffer.data(fg_off + 135);

        auto g_zzz_xxxy = cbuffer.data(fg_off + 136);

        auto g_zzz_xxxz = cbuffer.data(fg_off + 137);

        auto g_zzz_xxyy = cbuffer.data(fg_off + 138);

        auto g_zzz_xxyz = cbuffer.data(fg_off + 139);

        auto g_zzz_xxzz = cbuffer.data(fg_off + 140);

        auto g_zzz_xyyy = cbuffer.data(fg_off + 141);

        auto g_zzz_xyyz = cbuffer.data(fg_off + 142);

        auto g_zzz_xyzz = cbuffer.data(fg_off + 143);

        auto g_zzz_xzzz = cbuffer.data(fg_off + 144);

        auto g_zzz_yyyy = cbuffer.data(fg_off + 145);

        auto g_zzz_yyyz = cbuffer.data(fg_off + 146);

        auto g_zzz_yyzz = cbuffer.data(fg_off + 147);

        auto g_zzz_yzzz = cbuffer.data(fg_off + 148);

        auto g_zzz_zzzz = cbuffer.data(fg_off + 149);

        /// set up bra offset for contr_buffer_xgf

        const auto gf_off = idx_xgf + i * 150;

        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_xxxx_xxx = cbuffer.data(gf_off + 0);

        auto g_xxxx_xxy = cbuffer.data(gf_off + 1);

        auto g_xxxx_xxz = cbuffer.data(gf_off + 2);

        auto g_xxxx_xyy = cbuffer.data(gf_off + 3);

        auto g_xxxx_xyz = cbuffer.data(gf_off + 4);

        auto g_xxxx_xzz = cbuffer.data(gf_off + 5);

        auto g_xxxx_yyy = cbuffer.data(gf_off + 6);

        auto g_xxxx_yyz = cbuffer.data(gf_off + 7);

        auto g_xxxx_yzz = cbuffer.data(gf_off + 8);

        auto g_xxxx_zzz = cbuffer.data(gf_off + 9);

        #pragma omp simd aligned(cd_x, g_xxx_xxx, g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxy, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxz, g_xxx_xxzz, g_xxx_xyy, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyz, g_xxx_xyzz, g_xxx_xzz, g_xxx_xzzz, g_xxx_yyy, g_xxx_yyz, g_xxx_yzz, g_xxx_zzz, g_xxxx_xxx, g_xxxx_xxy, g_xxxx_xxz, g_xxxx_xyy, g_xxxx_xyz, g_xxxx_xzz, g_xxxx_yyy, g_xxxx_yyz, g_xxxx_yzz, g_xxxx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxx_xxx[k] = -g_xxx_xxx[k] * cd_x[k] + g_xxx_xxxx[k];

            g_xxxx_xxy[k] = -g_xxx_xxy[k] * cd_x[k] + g_xxx_xxxy[k];

            g_xxxx_xxz[k] = -g_xxx_xxz[k] * cd_x[k] + g_xxx_xxxz[k];

            g_xxxx_xyy[k] = -g_xxx_xyy[k] * cd_x[k] + g_xxx_xxyy[k];

            g_xxxx_xyz[k] = -g_xxx_xyz[k] * cd_x[k] + g_xxx_xxyz[k];

            g_xxxx_xzz[k] = -g_xxx_xzz[k] * cd_x[k] + g_xxx_xxzz[k];

            g_xxxx_yyy[k] = -g_xxx_yyy[k] * cd_x[k] + g_xxx_xyyy[k];

            g_xxxx_yyz[k] = -g_xxx_yyz[k] * cd_x[k] + g_xxx_xyyz[k];

            g_xxxx_yzz[k] = -g_xxx_yzz[k] * cd_x[k] + g_xxx_xyzz[k];

            g_xxxx_zzz[k] = -g_xxx_zzz[k] * cd_x[k] + g_xxx_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_xxxy_xxx = cbuffer.data(gf_off + 10);

        auto g_xxxy_xxy = cbuffer.data(gf_off + 11);

        auto g_xxxy_xxz = cbuffer.data(gf_off + 12);

        auto g_xxxy_xyy = cbuffer.data(gf_off + 13);

        auto g_xxxy_xyz = cbuffer.data(gf_off + 14);

        auto g_xxxy_xzz = cbuffer.data(gf_off + 15);

        auto g_xxxy_yyy = cbuffer.data(gf_off + 16);

        auto g_xxxy_yyz = cbuffer.data(gf_off + 17);

        auto g_xxxy_yzz = cbuffer.data(gf_off + 18);

        auto g_xxxy_zzz = cbuffer.data(gf_off + 19);

        #pragma omp simd aligned(cd_x, g_xxxy_xxx, g_xxxy_xxy, g_xxxy_xxz, g_xxxy_xyy, g_xxxy_xyz, g_xxxy_xzz, g_xxxy_yyy, g_xxxy_yyz, g_xxxy_yzz, g_xxxy_zzz, g_xxy_xxx, g_xxy_xxxx, g_xxy_xxxy, g_xxy_xxxz, g_xxy_xxy, g_xxy_xxyy, g_xxy_xxyz, g_xxy_xxz, g_xxy_xxzz, g_xxy_xyy, g_xxy_xyyy, g_xxy_xyyz, g_xxy_xyz, g_xxy_xyzz, g_xxy_xzz, g_xxy_xzzz, g_xxy_yyy, g_xxy_yyz, g_xxy_yzz, g_xxy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxy_xxx[k] = -g_xxy_xxx[k] * cd_x[k] + g_xxy_xxxx[k];

            g_xxxy_xxy[k] = -g_xxy_xxy[k] * cd_x[k] + g_xxy_xxxy[k];

            g_xxxy_xxz[k] = -g_xxy_xxz[k] * cd_x[k] + g_xxy_xxxz[k];

            g_xxxy_xyy[k] = -g_xxy_xyy[k] * cd_x[k] + g_xxy_xxyy[k];

            g_xxxy_xyz[k] = -g_xxy_xyz[k] * cd_x[k] + g_xxy_xxyz[k];

            g_xxxy_xzz[k] = -g_xxy_xzz[k] * cd_x[k] + g_xxy_xxzz[k];

            g_xxxy_yyy[k] = -g_xxy_yyy[k] * cd_x[k] + g_xxy_xyyy[k];

            g_xxxy_yyz[k] = -g_xxy_yyz[k] * cd_x[k] + g_xxy_xyyz[k];

            g_xxxy_yzz[k] = -g_xxy_yzz[k] * cd_x[k] + g_xxy_xyzz[k];

            g_xxxy_zzz[k] = -g_xxy_zzz[k] * cd_x[k] + g_xxy_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_xxxz_xxx = cbuffer.data(gf_off + 20);

        auto g_xxxz_xxy = cbuffer.data(gf_off + 21);

        auto g_xxxz_xxz = cbuffer.data(gf_off + 22);

        auto g_xxxz_xyy = cbuffer.data(gf_off + 23);

        auto g_xxxz_xyz = cbuffer.data(gf_off + 24);

        auto g_xxxz_xzz = cbuffer.data(gf_off + 25);

        auto g_xxxz_yyy = cbuffer.data(gf_off + 26);

        auto g_xxxz_yyz = cbuffer.data(gf_off + 27);

        auto g_xxxz_yzz = cbuffer.data(gf_off + 28);

        auto g_xxxz_zzz = cbuffer.data(gf_off + 29);

        #pragma omp simd aligned(cd_x, g_xxxz_xxx, g_xxxz_xxy, g_xxxz_xxz, g_xxxz_xyy, g_xxxz_xyz, g_xxxz_xzz, g_xxxz_yyy, g_xxxz_yyz, g_xxxz_yzz, g_xxxz_zzz, g_xxz_xxx, g_xxz_xxxx, g_xxz_xxxy, g_xxz_xxxz, g_xxz_xxy, g_xxz_xxyy, g_xxz_xxyz, g_xxz_xxz, g_xxz_xxzz, g_xxz_xyy, g_xxz_xyyy, g_xxz_xyyz, g_xxz_xyz, g_xxz_xyzz, g_xxz_xzz, g_xxz_xzzz, g_xxz_yyy, g_xxz_yyz, g_xxz_yzz, g_xxz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxz_xxx[k] = -g_xxz_xxx[k] * cd_x[k] + g_xxz_xxxx[k];

            g_xxxz_xxy[k] = -g_xxz_xxy[k] * cd_x[k] + g_xxz_xxxy[k];

            g_xxxz_xxz[k] = -g_xxz_xxz[k] * cd_x[k] + g_xxz_xxxz[k];

            g_xxxz_xyy[k] = -g_xxz_xyy[k] * cd_x[k] + g_xxz_xxyy[k];

            g_xxxz_xyz[k] = -g_xxz_xyz[k] * cd_x[k] + g_xxz_xxyz[k];

            g_xxxz_xzz[k] = -g_xxz_xzz[k] * cd_x[k] + g_xxz_xxzz[k];

            g_xxxz_yyy[k] = -g_xxz_yyy[k] * cd_x[k] + g_xxz_xyyy[k];

            g_xxxz_yyz[k] = -g_xxz_yyz[k] * cd_x[k] + g_xxz_xyyz[k];

            g_xxxz_yzz[k] = -g_xxz_yzz[k] * cd_x[k] + g_xxz_xyzz[k];

            g_xxxz_zzz[k] = -g_xxz_zzz[k] * cd_x[k] + g_xxz_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_xxyy_xxx = cbuffer.data(gf_off + 30);

        auto g_xxyy_xxy = cbuffer.data(gf_off + 31);

        auto g_xxyy_xxz = cbuffer.data(gf_off + 32);

        auto g_xxyy_xyy = cbuffer.data(gf_off + 33);

        auto g_xxyy_xyz = cbuffer.data(gf_off + 34);

        auto g_xxyy_xzz = cbuffer.data(gf_off + 35);

        auto g_xxyy_yyy = cbuffer.data(gf_off + 36);

        auto g_xxyy_yyz = cbuffer.data(gf_off + 37);

        auto g_xxyy_yzz = cbuffer.data(gf_off + 38);

        auto g_xxyy_zzz = cbuffer.data(gf_off + 39);

        #pragma omp simd aligned(cd_x, g_xxyy_xxx, g_xxyy_xxy, g_xxyy_xxz, g_xxyy_xyy, g_xxyy_xyz, g_xxyy_xzz, g_xxyy_yyy, g_xxyy_yyz, g_xxyy_yzz, g_xxyy_zzz, g_xyy_xxx, g_xyy_xxxx, g_xyy_xxxy, g_xyy_xxxz, g_xyy_xxy, g_xyy_xxyy, g_xyy_xxyz, g_xyy_xxz, g_xyy_xxzz, g_xyy_xyy, g_xyy_xyyy, g_xyy_xyyz, g_xyy_xyz, g_xyy_xyzz, g_xyy_xzz, g_xyy_xzzz, g_xyy_yyy, g_xyy_yyz, g_xyy_yzz, g_xyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyy_xxx[k] = -g_xyy_xxx[k] * cd_x[k] + g_xyy_xxxx[k];

            g_xxyy_xxy[k] = -g_xyy_xxy[k] * cd_x[k] + g_xyy_xxxy[k];

            g_xxyy_xxz[k] = -g_xyy_xxz[k] * cd_x[k] + g_xyy_xxxz[k];

            g_xxyy_xyy[k] = -g_xyy_xyy[k] * cd_x[k] + g_xyy_xxyy[k];

            g_xxyy_xyz[k] = -g_xyy_xyz[k] * cd_x[k] + g_xyy_xxyz[k];

            g_xxyy_xzz[k] = -g_xyy_xzz[k] * cd_x[k] + g_xyy_xxzz[k];

            g_xxyy_yyy[k] = -g_xyy_yyy[k] * cd_x[k] + g_xyy_xyyy[k];

            g_xxyy_yyz[k] = -g_xyy_yyz[k] * cd_x[k] + g_xyy_xyyz[k];

            g_xxyy_yzz[k] = -g_xyy_yzz[k] * cd_x[k] + g_xyy_xyzz[k];

            g_xxyy_zzz[k] = -g_xyy_zzz[k] * cd_x[k] + g_xyy_xzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_xxyz_xxx = cbuffer.data(gf_off + 40);

        auto g_xxyz_xxy = cbuffer.data(gf_off + 41);

        auto g_xxyz_xxz = cbuffer.data(gf_off + 42);

        auto g_xxyz_xyy = cbuffer.data(gf_off + 43);

        auto g_xxyz_xyz = cbuffer.data(gf_off + 44);

        auto g_xxyz_xzz = cbuffer.data(gf_off + 45);

        auto g_xxyz_yyy = cbuffer.data(gf_off + 46);

        auto g_xxyz_yyz = cbuffer.data(gf_off + 47);

        auto g_xxyz_yzz = cbuffer.data(gf_off + 48);

        auto g_xxyz_zzz = cbuffer.data(gf_off + 49);

        #pragma omp simd aligned(cd_x, g_xxyz_xxx, g_xxyz_xxy, g_xxyz_xxz, g_xxyz_xyy, g_xxyz_xyz, g_xxyz_xzz, g_xxyz_yyy, g_xxyz_yyz, g_xxyz_yzz, g_xxyz_zzz, g_xyz_xxx, g_xyz_xxxx, g_xyz_xxxy, g_xyz_xxxz, g_xyz_xxy, g_xyz_xxyy, g_xyz_xxyz, g_xyz_xxz, g_xyz_xxzz, g_xyz_xyy, g_xyz_xyyy, g_xyz_xyyz, g_xyz_xyz, g_xyz_xyzz, g_xyz_xzz, g_xyz_xzzz, g_xyz_yyy, g_xyz_yyz, g_xyz_yzz, g_xyz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyz_xxx[k] = -g_xyz_xxx[k] * cd_x[k] + g_xyz_xxxx[k];

            g_xxyz_xxy[k] = -g_xyz_xxy[k] * cd_x[k] + g_xyz_xxxy[k];

            g_xxyz_xxz[k] = -g_xyz_xxz[k] * cd_x[k] + g_xyz_xxxz[k];

            g_xxyz_xyy[k] = -g_xyz_xyy[k] * cd_x[k] + g_xyz_xxyy[k];

            g_xxyz_xyz[k] = -g_xyz_xyz[k] * cd_x[k] + g_xyz_xxyz[k];

            g_xxyz_xzz[k] = -g_xyz_xzz[k] * cd_x[k] + g_xyz_xxzz[k];

            g_xxyz_yyy[k] = -g_xyz_yyy[k] * cd_x[k] + g_xyz_xyyy[k];

            g_xxyz_yyz[k] = -g_xyz_yyz[k] * cd_x[k] + g_xyz_xyyz[k];

            g_xxyz_yzz[k] = -g_xyz_yzz[k] * cd_x[k] + g_xyz_xyzz[k];

            g_xxyz_zzz[k] = -g_xyz_zzz[k] * cd_x[k] + g_xyz_xzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_xxzz_xxx = cbuffer.data(gf_off + 50);

        auto g_xxzz_xxy = cbuffer.data(gf_off + 51);

        auto g_xxzz_xxz = cbuffer.data(gf_off + 52);

        auto g_xxzz_xyy = cbuffer.data(gf_off + 53);

        auto g_xxzz_xyz = cbuffer.data(gf_off + 54);

        auto g_xxzz_xzz = cbuffer.data(gf_off + 55);

        auto g_xxzz_yyy = cbuffer.data(gf_off + 56);

        auto g_xxzz_yyz = cbuffer.data(gf_off + 57);

        auto g_xxzz_yzz = cbuffer.data(gf_off + 58);

        auto g_xxzz_zzz = cbuffer.data(gf_off + 59);

        #pragma omp simd aligned(cd_x, g_xxzz_xxx, g_xxzz_xxy, g_xxzz_xxz, g_xxzz_xyy, g_xxzz_xyz, g_xxzz_xzz, g_xxzz_yyy, g_xxzz_yyz, g_xxzz_yzz, g_xxzz_zzz, g_xzz_xxx, g_xzz_xxxx, g_xzz_xxxy, g_xzz_xxxz, g_xzz_xxy, g_xzz_xxyy, g_xzz_xxyz, g_xzz_xxz, g_xzz_xxzz, g_xzz_xyy, g_xzz_xyyy, g_xzz_xyyz, g_xzz_xyz, g_xzz_xyzz, g_xzz_xzz, g_xzz_xzzz, g_xzz_yyy, g_xzz_yyz, g_xzz_yzz, g_xzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxzz_xxx[k] = -g_xzz_xxx[k] * cd_x[k] + g_xzz_xxxx[k];

            g_xxzz_xxy[k] = -g_xzz_xxy[k] * cd_x[k] + g_xzz_xxxy[k];

            g_xxzz_xxz[k] = -g_xzz_xxz[k] * cd_x[k] + g_xzz_xxxz[k];

            g_xxzz_xyy[k] = -g_xzz_xyy[k] * cd_x[k] + g_xzz_xxyy[k];

            g_xxzz_xyz[k] = -g_xzz_xyz[k] * cd_x[k] + g_xzz_xxyz[k];

            g_xxzz_xzz[k] = -g_xzz_xzz[k] * cd_x[k] + g_xzz_xxzz[k];

            g_xxzz_yyy[k] = -g_xzz_yyy[k] * cd_x[k] + g_xzz_xyyy[k];

            g_xxzz_yyz[k] = -g_xzz_yyz[k] * cd_x[k] + g_xzz_xyyz[k];

            g_xxzz_yzz[k] = -g_xzz_yzz[k] * cd_x[k] + g_xzz_xyzz[k];

            g_xxzz_zzz[k] = -g_xzz_zzz[k] * cd_x[k] + g_xzz_xzzz[k];
        }

        /// Set up 60-70 components of targeted buffer : cbuffer.data(

        auto g_xyyy_xxx = cbuffer.data(gf_off + 60);

        auto g_xyyy_xxy = cbuffer.data(gf_off + 61);

        auto g_xyyy_xxz = cbuffer.data(gf_off + 62);

        auto g_xyyy_xyy = cbuffer.data(gf_off + 63);

        auto g_xyyy_xyz = cbuffer.data(gf_off + 64);

        auto g_xyyy_xzz = cbuffer.data(gf_off + 65);

        auto g_xyyy_yyy = cbuffer.data(gf_off + 66);

        auto g_xyyy_yyz = cbuffer.data(gf_off + 67);

        auto g_xyyy_yzz = cbuffer.data(gf_off + 68);

        auto g_xyyy_zzz = cbuffer.data(gf_off + 69);

        #pragma omp simd aligned(cd_x, g_xyyy_xxx, g_xyyy_xxy, g_xyyy_xxz, g_xyyy_xyy, g_xyyy_xyz, g_xyyy_xzz, g_xyyy_yyy, g_xyyy_yyz, g_xyyy_yzz, g_xyyy_zzz, g_yyy_xxx, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxy, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxz, g_yyy_xxzz, g_yyy_xyy, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyz, g_yyy_xyzz, g_yyy_xzz, g_yyy_xzzz, g_yyy_yyy, g_yyy_yyz, g_yyy_yzz, g_yyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyy_xxx[k] = -g_yyy_xxx[k] * cd_x[k] + g_yyy_xxxx[k];

            g_xyyy_xxy[k] = -g_yyy_xxy[k] * cd_x[k] + g_yyy_xxxy[k];

            g_xyyy_xxz[k] = -g_yyy_xxz[k] * cd_x[k] + g_yyy_xxxz[k];

            g_xyyy_xyy[k] = -g_yyy_xyy[k] * cd_x[k] + g_yyy_xxyy[k];

            g_xyyy_xyz[k] = -g_yyy_xyz[k] * cd_x[k] + g_yyy_xxyz[k];

            g_xyyy_xzz[k] = -g_yyy_xzz[k] * cd_x[k] + g_yyy_xxzz[k];

            g_xyyy_yyy[k] = -g_yyy_yyy[k] * cd_x[k] + g_yyy_xyyy[k];

            g_xyyy_yyz[k] = -g_yyy_yyz[k] * cd_x[k] + g_yyy_xyyz[k];

            g_xyyy_yzz[k] = -g_yyy_yzz[k] * cd_x[k] + g_yyy_xyzz[k];

            g_xyyy_zzz[k] = -g_yyy_zzz[k] * cd_x[k] + g_yyy_xzzz[k];
        }

        /// Set up 70-80 components of targeted buffer : cbuffer.data(

        auto g_xyyz_xxx = cbuffer.data(gf_off + 70);

        auto g_xyyz_xxy = cbuffer.data(gf_off + 71);

        auto g_xyyz_xxz = cbuffer.data(gf_off + 72);

        auto g_xyyz_xyy = cbuffer.data(gf_off + 73);

        auto g_xyyz_xyz = cbuffer.data(gf_off + 74);

        auto g_xyyz_xzz = cbuffer.data(gf_off + 75);

        auto g_xyyz_yyy = cbuffer.data(gf_off + 76);

        auto g_xyyz_yyz = cbuffer.data(gf_off + 77);

        auto g_xyyz_yzz = cbuffer.data(gf_off + 78);

        auto g_xyyz_zzz = cbuffer.data(gf_off + 79);

        #pragma omp simd aligned(cd_x, g_xyyz_xxx, g_xyyz_xxy, g_xyyz_xxz, g_xyyz_xyy, g_xyyz_xyz, g_xyyz_xzz, g_xyyz_yyy, g_xyyz_yyz, g_xyyz_yzz, g_xyyz_zzz, g_yyz_xxx, g_yyz_xxxx, g_yyz_xxxy, g_yyz_xxxz, g_yyz_xxy, g_yyz_xxyy, g_yyz_xxyz, g_yyz_xxz, g_yyz_xxzz, g_yyz_xyy, g_yyz_xyyy, g_yyz_xyyz, g_yyz_xyz, g_yyz_xyzz, g_yyz_xzz, g_yyz_xzzz, g_yyz_yyy, g_yyz_yyz, g_yyz_yzz, g_yyz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyz_xxx[k] = -g_yyz_xxx[k] * cd_x[k] + g_yyz_xxxx[k];

            g_xyyz_xxy[k] = -g_yyz_xxy[k] * cd_x[k] + g_yyz_xxxy[k];

            g_xyyz_xxz[k] = -g_yyz_xxz[k] * cd_x[k] + g_yyz_xxxz[k];

            g_xyyz_xyy[k] = -g_yyz_xyy[k] * cd_x[k] + g_yyz_xxyy[k];

            g_xyyz_xyz[k] = -g_yyz_xyz[k] * cd_x[k] + g_yyz_xxyz[k];

            g_xyyz_xzz[k] = -g_yyz_xzz[k] * cd_x[k] + g_yyz_xxzz[k];

            g_xyyz_yyy[k] = -g_yyz_yyy[k] * cd_x[k] + g_yyz_xyyy[k];

            g_xyyz_yyz[k] = -g_yyz_yyz[k] * cd_x[k] + g_yyz_xyyz[k];

            g_xyyz_yzz[k] = -g_yyz_yzz[k] * cd_x[k] + g_yyz_xyzz[k];

            g_xyyz_zzz[k] = -g_yyz_zzz[k] * cd_x[k] + g_yyz_xzzz[k];
        }

        /// Set up 80-90 components of targeted buffer : cbuffer.data(

        auto g_xyzz_xxx = cbuffer.data(gf_off + 80);

        auto g_xyzz_xxy = cbuffer.data(gf_off + 81);

        auto g_xyzz_xxz = cbuffer.data(gf_off + 82);

        auto g_xyzz_xyy = cbuffer.data(gf_off + 83);

        auto g_xyzz_xyz = cbuffer.data(gf_off + 84);

        auto g_xyzz_xzz = cbuffer.data(gf_off + 85);

        auto g_xyzz_yyy = cbuffer.data(gf_off + 86);

        auto g_xyzz_yyz = cbuffer.data(gf_off + 87);

        auto g_xyzz_yzz = cbuffer.data(gf_off + 88);

        auto g_xyzz_zzz = cbuffer.data(gf_off + 89);

        #pragma omp simd aligned(cd_x, g_xyzz_xxx, g_xyzz_xxy, g_xyzz_xxz, g_xyzz_xyy, g_xyzz_xyz, g_xyzz_xzz, g_xyzz_yyy, g_xyzz_yyz, g_xyzz_yzz, g_xyzz_zzz, g_yzz_xxx, g_yzz_xxxx, g_yzz_xxxy, g_yzz_xxxz, g_yzz_xxy, g_yzz_xxyy, g_yzz_xxyz, g_yzz_xxz, g_yzz_xxzz, g_yzz_xyy, g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyz, g_yzz_xyzz, g_yzz_xzz, g_yzz_xzzz, g_yzz_yyy, g_yzz_yyz, g_yzz_yzz, g_yzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyzz_xxx[k] = -g_yzz_xxx[k] * cd_x[k] + g_yzz_xxxx[k];

            g_xyzz_xxy[k] = -g_yzz_xxy[k] * cd_x[k] + g_yzz_xxxy[k];

            g_xyzz_xxz[k] = -g_yzz_xxz[k] * cd_x[k] + g_yzz_xxxz[k];

            g_xyzz_xyy[k] = -g_yzz_xyy[k] * cd_x[k] + g_yzz_xxyy[k];

            g_xyzz_xyz[k] = -g_yzz_xyz[k] * cd_x[k] + g_yzz_xxyz[k];

            g_xyzz_xzz[k] = -g_yzz_xzz[k] * cd_x[k] + g_yzz_xxzz[k];

            g_xyzz_yyy[k] = -g_yzz_yyy[k] * cd_x[k] + g_yzz_xyyy[k];

            g_xyzz_yyz[k] = -g_yzz_yyz[k] * cd_x[k] + g_yzz_xyyz[k];

            g_xyzz_yzz[k] = -g_yzz_yzz[k] * cd_x[k] + g_yzz_xyzz[k];

            g_xyzz_zzz[k] = -g_yzz_zzz[k] * cd_x[k] + g_yzz_xzzz[k];
        }

        /// Set up 90-100 components of targeted buffer : cbuffer.data(

        auto g_xzzz_xxx = cbuffer.data(gf_off + 90);

        auto g_xzzz_xxy = cbuffer.data(gf_off + 91);

        auto g_xzzz_xxz = cbuffer.data(gf_off + 92);

        auto g_xzzz_xyy = cbuffer.data(gf_off + 93);

        auto g_xzzz_xyz = cbuffer.data(gf_off + 94);

        auto g_xzzz_xzz = cbuffer.data(gf_off + 95);

        auto g_xzzz_yyy = cbuffer.data(gf_off + 96);

        auto g_xzzz_yyz = cbuffer.data(gf_off + 97);

        auto g_xzzz_yzz = cbuffer.data(gf_off + 98);

        auto g_xzzz_zzz = cbuffer.data(gf_off + 99);

        #pragma omp simd aligned(cd_x, g_xzzz_xxx, g_xzzz_xxy, g_xzzz_xxz, g_xzzz_xyy, g_xzzz_xyz, g_xzzz_xzz, g_xzzz_yyy, g_xzzz_yyz, g_xzzz_yzz, g_xzzz_zzz, g_zzz_xxx, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz, g_zzz_xxy, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxz, g_zzz_xxzz, g_zzz_xyy, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyz, g_zzz_xyzz, g_zzz_xzz, g_zzz_xzzz, g_zzz_yyy, g_zzz_yyz, g_zzz_yzz, g_zzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzzz_xxx[k] = -g_zzz_xxx[k] * cd_x[k] + g_zzz_xxxx[k];

            g_xzzz_xxy[k] = -g_zzz_xxy[k] * cd_x[k] + g_zzz_xxxy[k];

            g_xzzz_xxz[k] = -g_zzz_xxz[k] * cd_x[k] + g_zzz_xxxz[k];

            g_xzzz_xyy[k] = -g_zzz_xyy[k] * cd_x[k] + g_zzz_xxyy[k];

            g_xzzz_xyz[k] = -g_zzz_xyz[k] * cd_x[k] + g_zzz_xxyz[k];

            g_xzzz_xzz[k] = -g_zzz_xzz[k] * cd_x[k] + g_zzz_xxzz[k];

            g_xzzz_yyy[k] = -g_zzz_yyy[k] * cd_x[k] + g_zzz_xyyy[k];

            g_xzzz_yyz[k] = -g_zzz_yyz[k] * cd_x[k] + g_zzz_xyyz[k];

            g_xzzz_yzz[k] = -g_zzz_yzz[k] * cd_x[k] + g_zzz_xyzz[k];

            g_xzzz_zzz[k] = -g_zzz_zzz[k] * cd_x[k] + g_zzz_xzzz[k];
        }

        /// Set up 100-110 components of targeted buffer : cbuffer.data(

        auto g_yyyy_xxx = cbuffer.data(gf_off + 100);

        auto g_yyyy_xxy = cbuffer.data(gf_off + 101);

        auto g_yyyy_xxz = cbuffer.data(gf_off + 102);

        auto g_yyyy_xyy = cbuffer.data(gf_off + 103);

        auto g_yyyy_xyz = cbuffer.data(gf_off + 104);

        auto g_yyyy_xzz = cbuffer.data(gf_off + 105);

        auto g_yyyy_yyy = cbuffer.data(gf_off + 106);

        auto g_yyyy_yyz = cbuffer.data(gf_off + 107);

        auto g_yyyy_yzz = cbuffer.data(gf_off + 108);

        auto g_yyyy_zzz = cbuffer.data(gf_off + 109);

        #pragma omp simd aligned(cd_y, g_yyy_xxx, g_yyy_xxxy, g_yyy_xxy, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxz, g_yyy_xyy, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyz, g_yyy_xyzz, g_yyy_xzz, g_yyy_yyy, g_yyy_yyyy, g_yyy_yyyz, g_yyy_yyz, g_yyy_yyzz, g_yyy_yzz, g_yyy_yzzz, g_yyy_zzz, g_yyyy_xxx, g_yyyy_xxy, g_yyyy_xxz, g_yyyy_xyy, g_yyyy_xyz, g_yyyy_xzz, g_yyyy_yyy, g_yyyy_yyz, g_yyyy_yzz, g_yyyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyy_xxx[k] = -g_yyy_xxx[k] * cd_y[k] + g_yyy_xxxy[k];

            g_yyyy_xxy[k] = -g_yyy_xxy[k] * cd_y[k] + g_yyy_xxyy[k];

            g_yyyy_xxz[k] = -g_yyy_xxz[k] * cd_y[k] + g_yyy_xxyz[k];

            g_yyyy_xyy[k] = -g_yyy_xyy[k] * cd_y[k] + g_yyy_xyyy[k];

            g_yyyy_xyz[k] = -g_yyy_xyz[k] * cd_y[k] + g_yyy_xyyz[k];

            g_yyyy_xzz[k] = -g_yyy_xzz[k] * cd_y[k] + g_yyy_xyzz[k];

            g_yyyy_yyy[k] = -g_yyy_yyy[k] * cd_y[k] + g_yyy_yyyy[k];

            g_yyyy_yyz[k] = -g_yyy_yyz[k] * cd_y[k] + g_yyy_yyyz[k];

            g_yyyy_yzz[k] = -g_yyy_yzz[k] * cd_y[k] + g_yyy_yyzz[k];

            g_yyyy_zzz[k] = -g_yyy_zzz[k] * cd_y[k] + g_yyy_yzzz[k];
        }

        /// Set up 110-120 components of targeted buffer : cbuffer.data(

        auto g_yyyz_xxx = cbuffer.data(gf_off + 110);

        auto g_yyyz_xxy = cbuffer.data(gf_off + 111);

        auto g_yyyz_xxz = cbuffer.data(gf_off + 112);

        auto g_yyyz_xyy = cbuffer.data(gf_off + 113);

        auto g_yyyz_xyz = cbuffer.data(gf_off + 114);

        auto g_yyyz_xzz = cbuffer.data(gf_off + 115);

        auto g_yyyz_yyy = cbuffer.data(gf_off + 116);

        auto g_yyyz_yyz = cbuffer.data(gf_off + 117);

        auto g_yyyz_yzz = cbuffer.data(gf_off + 118);

        auto g_yyyz_zzz = cbuffer.data(gf_off + 119);

        #pragma omp simd aligned(cd_y, g_yyyz_xxx, g_yyyz_xxy, g_yyyz_xxz, g_yyyz_xyy, g_yyyz_xyz, g_yyyz_xzz, g_yyyz_yyy, g_yyyz_yyz, g_yyyz_yzz, g_yyyz_zzz, g_yyz_xxx, g_yyz_xxxy, g_yyz_xxy, g_yyz_xxyy, g_yyz_xxyz, g_yyz_xxz, g_yyz_xyy, g_yyz_xyyy, g_yyz_xyyz, g_yyz_xyz, g_yyz_xyzz, g_yyz_xzz, g_yyz_yyy, g_yyz_yyyy, g_yyz_yyyz, g_yyz_yyz, g_yyz_yyzz, g_yyz_yzz, g_yyz_yzzz, g_yyz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyz_xxx[k] = -g_yyz_xxx[k] * cd_y[k] + g_yyz_xxxy[k];

            g_yyyz_xxy[k] = -g_yyz_xxy[k] * cd_y[k] + g_yyz_xxyy[k];

            g_yyyz_xxz[k] = -g_yyz_xxz[k] * cd_y[k] + g_yyz_xxyz[k];

            g_yyyz_xyy[k] = -g_yyz_xyy[k] * cd_y[k] + g_yyz_xyyy[k];

            g_yyyz_xyz[k] = -g_yyz_xyz[k] * cd_y[k] + g_yyz_xyyz[k];

            g_yyyz_xzz[k] = -g_yyz_xzz[k] * cd_y[k] + g_yyz_xyzz[k];

            g_yyyz_yyy[k] = -g_yyz_yyy[k] * cd_y[k] + g_yyz_yyyy[k];

            g_yyyz_yyz[k] = -g_yyz_yyz[k] * cd_y[k] + g_yyz_yyyz[k];

            g_yyyz_yzz[k] = -g_yyz_yzz[k] * cd_y[k] + g_yyz_yyzz[k];

            g_yyyz_zzz[k] = -g_yyz_zzz[k] * cd_y[k] + g_yyz_yzzz[k];
        }

        /// Set up 120-130 components of targeted buffer : cbuffer.data(

        auto g_yyzz_xxx = cbuffer.data(gf_off + 120);

        auto g_yyzz_xxy = cbuffer.data(gf_off + 121);

        auto g_yyzz_xxz = cbuffer.data(gf_off + 122);

        auto g_yyzz_xyy = cbuffer.data(gf_off + 123);

        auto g_yyzz_xyz = cbuffer.data(gf_off + 124);

        auto g_yyzz_xzz = cbuffer.data(gf_off + 125);

        auto g_yyzz_yyy = cbuffer.data(gf_off + 126);

        auto g_yyzz_yyz = cbuffer.data(gf_off + 127);

        auto g_yyzz_yzz = cbuffer.data(gf_off + 128);

        auto g_yyzz_zzz = cbuffer.data(gf_off + 129);

        #pragma omp simd aligned(cd_y, g_yyzz_xxx, g_yyzz_xxy, g_yyzz_xxz, g_yyzz_xyy, g_yyzz_xyz, g_yyzz_xzz, g_yyzz_yyy, g_yyzz_yyz, g_yyzz_yzz, g_yyzz_zzz, g_yzz_xxx, g_yzz_xxxy, g_yzz_xxy, g_yzz_xxyy, g_yzz_xxyz, g_yzz_xxz, g_yzz_xyy, g_yzz_xyyy, g_yzz_xyyz, g_yzz_xyz, g_yzz_xyzz, g_yzz_xzz, g_yzz_yyy, g_yzz_yyyy, g_yzz_yyyz, g_yzz_yyz, g_yzz_yyzz, g_yzz_yzz, g_yzz_yzzz, g_yzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyzz_xxx[k] = -g_yzz_xxx[k] * cd_y[k] + g_yzz_xxxy[k];

            g_yyzz_xxy[k] = -g_yzz_xxy[k] * cd_y[k] + g_yzz_xxyy[k];

            g_yyzz_xxz[k] = -g_yzz_xxz[k] * cd_y[k] + g_yzz_xxyz[k];

            g_yyzz_xyy[k] = -g_yzz_xyy[k] * cd_y[k] + g_yzz_xyyy[k];

            g_yyzz_xyz[k] = -g_yzz_xyz[k] * cd_y[k] + g_yzz_xyyz[k];

            g_yyzz_xzz[k] = -g_yzz_xzz[k] * cd_y[k] + g_yzz_xyzz[k];

            g_yyzz_yyy[k] = -g_yzz_yyy[k] * cd_y[k] + g_yzz_yyyy[k];

            g_yyzz_yyz[k] = -g_yzz_yyz[k] * cd_y[k] + g_yzz_yyyz[k];

            g_yyzz_yzz[k] = -g_yzz_yzz[k] * cd_y[k] + g_yzz_yyzz[k];

            g_yyzz_zzz[k] = -g_yzz_zzz[k] * cd_y[k] + g_yzz_yzzz[k];
        }

        /// Set up 130-140 components of targeted buffer : cbuffer.data(

        auto g_yzzz_xxx = cbuffer.data(gf_off + 130);

        auto g_yzzz_xxy = cbuffer.data(gf_off + 131);

        auto g_yzzz_xxz = cbuffer.data(gf_off + 132);

        auto g_yzzz_xyy = cbuffer.data(gf_off + 133);

        auto g_yzzz_xyz = cbuffer.data(gf_off + 134);

        auto g_yzzz_xzz = cbuffer.data(gf_off + 135);

        auto g_yzzz_yyy = cbuffer.data(gf_off + 136);

        auto g_yzzz_yyz = cbuffer.data(gf_off + 137);

        auto g_yzzz_yzz = cbuffer.data(gf_off + 138);

        auto g_yzzz_zzz = cbuffer.data(gf_off + 139);

        #pragma omp simd aligned(cd_y, g_yzzz_xxx, g_yzzz_xxy, g_yzzz_xxz, g_yzzz_xyy, g_yzzz_xyz, g_yzzz_xzz, g_yzzz_yyy, g_yzzz_yyz, g_yzzz_yzz, g_yzzz_zzz, g_zzz_xxx, g_zzz_xxxy, g_zzz_xxy, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxz, g_zzz_xyy, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyz, g_zzz_xyzz, g_zzz_xzz, g_zzz_yyy, g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyz, g_zzz_yyzz, g_zzz_yzz, g_zzz_yzzz, g_zzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzzz_xxx[k] = -g_zzz_xxx[k] * cd_y[k] + g_zzz_xxxy[k];

            g_yzzz_xxy[k] = -g_zzz_xxy[k] * cd_y[k] + g_zzz_xxyy[k];

            g_yzzz_xxz[k] = -g_zzz_xxz[k] * cd_y[k] + g_zzz_xxyz[k];

            g_yzzz_xyy[k] = -g_zzz_xyy[k] * cd_y[k] + g_zzz_xyyy[k];

            g_yzzz_xyz[k] = -g_zzz_xyz[k] * cd_y[k] + g_zzz_xyyz[k];

            g_yzzz_xzz[k] = -g_zzz_xzz[k] * cd_y[k] + g_zzz_xyzz[k];

            g_yzzz_yyy[k] = -g_zzz_yyy[k] * cd_y[k] + g_zzz_yyyy[k];

            g_yzzz_yyz[k] = -g_zzz_yyz[k] * cd_y[k] + g_zzz_yyyz[k];

            g_yzzz_yzz[k] = -g_zzz_yzz[k] * cd_y[k] + g_zzz_yyzz[k];

            g_yzzz_zzz[k] = -g_zzz_zzz[k] * cd_y[k] + g_zzz_yzzz[k];
        }

        /// Set up 140-150 components of targeted buffer : cbuffer.data(

        auto g_zzzz_xxx = cbuffer.data(gf_off + 140);

        auto g_zzzz_xxy = cbuffer.data(gf_off + 141);

        auto g_zzzz_xxz = cbuffer.data(gf_off + 142);

        auto g_zzzz_xyy = cbuffer.data(gf_off + 143);

        auto g_zzzz_xyz = cbuffer.data(gf_off + 144);

        auto g_zzzz_xzz = cbuffer.data(gf_off + 145);

        auto g_zzzz_yyy = cbuffer.data(gf_off + 146);

        auto g_zzzz_yyz = cbuffer.data(gf_off + 147);

        auto g_zzzz_yzz = cbuffer.data(gf_off + 148);

        auto g_zzzz_zzz = cbuffer.data(gf_off + 149);

        #pragma omp simd aligned(cd_z, g_zzz_xxx, g_zzz_xxxz, g_zzz_xxy, g_zzz_xxyz, g_zzz_xxz, g_zzz_xxzz, g_zzz_xyy, g_zzz_xyyz, g_zzz_xyz, g_zzz_xyzz, g_zzz_xzz, g_zzz_xzzz, g_zzz_yyy, g_zzz_yyyz, g_zzz_yyz, g_zzz_yyzz, g_zzz_yzz, g_zzz_yzzz, g_zzz_zzz, g_zzz_zzzz, g_zzzz_xxx, g_zzzz_xxy, g_zzzz_xxz, g_zzzz_xyy, g_zzzz_xyz, g_zzzz_xzz, g_zzzz_yyy, g_zzzz_yyz, g_zzzz_yzz, g_zzzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzzz_xxx[k] = -g_zzz_xxx[k] * cd_z[k] + g_zzz_xxxz[k];

            g_zzzz_xxy[k] = -g_zzz_xxy[k] * cd_z[k] + g_zzz_xxyz[k];

            g_zzzz_xxz[k] = -g_zzz_xxz[k] * cd_z[k] + g_zzz_xxzz[k];

            g_zzzz_xyy[k] = -g_zzz_xyy[k] * cd_z[k] + g_zzz_xyyz[k];

            g_zzzz_xyz[k] = -g_zzz_xyz[k] * cd_z[k] + g_zzz_xyzz[k];

            g_zzzz_xzz[k] = -g_zzz_xzz[k] * cd_z[k] + g_zzz_xzzz[k];

            g_zzzz_yyy[k] = -g_zzz_yyy[k] * cd_z[k] + g_zzz_yyyz[k];

            g_zzzz_yyz[k] = -g_zzz_yyz[k] * cd_z[k] + g_zzz_yyzz[k];

            g_zzzz_yzz[k] = -g_zzz_yzz[k] * cd_z[k] + g_zzz_yzzz[k];

            g_zzzz_zzz[k] = -g_zzz_zzz[k] * cd_z[k] + g_zzz_zzzz[k];
        }
    }
}

} // t3ceri namespace


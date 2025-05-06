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

#include "ElectronRepulsionContrRecXXPI.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxpi(CSimdArray<double>&       cbuffer,
                                     const size_t              idx_xxpi,
                                     CSimdArray<double>&       pbuffer,
                                     const size_t              idx_xxsi,
                                     const size_t              idx_xxsk,
                                     const CSimdArray<double>& factors,
                                     const size_t              idx_cd,
                                     const int                 a_angmom,
                                     const int                 b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{
        a_angmom,
    });

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{
        b_angmom,
    });

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSI

            const auto si_off = idx_xxsi + (i * bcomps + j) * 28;

            auto g_0_xxxxxx = pbuffer.data(si_off + 0);

            auto g_0_xxxxxy = pbuffer.data(si_off + 1);

            auto g_0_xxxxxz = pbuffer.data(si_off + 2);

            auto g_0_xxxxyy = pbuffer.data(si_off + 3);

            auto g_0_xxxxyz = pbuffer.data(si_off + 4);

            auto g_0_xxxxzz = pbuffer.data(si_off + 5);

            auto g_0_xxxyyy = pbuffer.data(si_off + 6);

            auto g_0_xxxyyz = pbuffer.data(si_off + 7);

            auto g_0_xxxyzz = pbuffer.data(si_off + 8);

            auto g_0_xxxzzz = pbuffer.data(si_off + 9);

            auto g_0_xxyyyy = pbuffer.data(si_off + 10);

            auto g_0_xxyyyz = pbuffer.data(si_off + 11);

            auto g_0_xxyyzz = pbuffer.data(si_off + 12);

            auto g_0_xxyzzz = pbuffer.data(si_off + 13);

            auto g_0_xxzzzz = pbuffer.data(si_off + 14);

            auto g_0_xyyyyy = pbuffer.data(si_off + 15);

            auto g_0_xyyyyz = pbuffer.data(si_off + 16);

            auto g_0_xyyyzz = pbuffer.data(si_off + 17);

            auto g_0_xyyzzz = pbuffer.data(si_off + 18);

            auto g_0_xyzzzz = pbuffer.data(si_off + 19);

            auto g_0_xzzzzz = pbuffer.data(si_off + 20);

            auto g_0_yyyyyy = pbuffer.data(si_off + 21);

            auto g_0_yyyyyz = pbuffer.data(si_off + 22);

            auto g_0_yyyyzz = pbuffer.data(si_off + 23);

            auto g_0_yyyzzz = pbuffer.data(si_off + 24);

            auto g_0_yyzzzz = pbuffer.data(si_off + 25);

            auto g_0_yzzzzz = pbuffer.data(si_off + 26);

            auto g_0_zzzzzz = pbuffer.data(si_off + 27);

            /// Set up components of auxilary buffer : SSSK

            const auto sk_off = idx_xxsk + (i * bcomps + j) * 36;

            auto g_0_xxxxxxx = pbuffer.data(sk_off + 0);

            auto g_0_xxxxxxy = pbuffer.data(sk_off + 1);

            auto g_0_xxxxxxz = pbuffer.data(sk_off + 2);

            auto g_0_xxxxxyy = pbuffer.data(sk_off + 3);

            auto g_0_xxxxxyz = pbuffer.data(sk_off + 4);

            auto g_0_xxxxxzz = pbuffer.data(sk_off + 5);

            auto g_0_xxxxyyy = pbuffer.data(sk_off + 6);

            auto g_0_xxxxyyz = pbuffer.data(sk_off + 7);

            auto g_0_xxxxyzz = pbuffer.data(sk_off + 8);

            auto g_0_xxxxzzz = pbuffer.data(sk_off + 9);

            auto g_0_xxxyyyy = pbuffer.data(sk_off + 10);

            auto g_0_xxxyyyz = pbuffer.data(sk_off + 11);

            auto g_0_xxxyyzz = pbuffer.data(sk_off + 12);

            auto g_0_xxxyzzz = pbuffer.data(sk_off + 13);

            auto g_0_xxxzzzz = pbuffer.data(sk_off + 14);

            auto g_0_xxyyyyy = pbuffer.data(sk_off + 15);

            auto g_0_xxyyyyz = pbuffer.data(sk_off + 16);

            auto g_0_xxyyyzz = pbuffer.data(sk_off + 17);

            auto g_0_xxyyzzz = pbuffer.data(sk_off + 18);

            auto g_0_xxyzzzz = pbuffer.data(sk_off + 19);

            auto g_0_xxzzzzz = pbuffer.data(sk_off + 20);

            auto g_0_xyyyyyy = pbuffer.data(sk_off + 21);

            auto g_0_xyyyyyz = pbuffer.data(sk_off + 22);

            auto g_0_xyyyyzz = pbuffer.data(sk_off + 23);

            auto g_0_xyyyzzz = pbuffer.data(sk_off + 24);

            auto g_0_xyyzzzz = pbuffer.data(sk_off + 25);

            auto g_0_xyzzzzz = pbuffer.data(sk_off + 26);

            auto g_0_xzzzzzz = pbuffer.data(sk_off + 27);

            auto g_0_yyyyyyy = pbuffer.data(sk_off + 28);

            auto g_0_yyyyyyz = pbuffer.data(sk_off + 29);

            auto g_0_yyyyyzz = pbuffer.data(sk_off + 30);

            auto g_0_yyyyzzz = pbuffer.data(sk_off + 31);

            auto g_0_yyyzzzz = pbuffer.data(sk_off + 32);

            auto g_0_yyzzzzz = pbuffer.data(sk_off + 33);

            auto g_0_yzzzzzz = pbuffer.data(sk_off + 34);

            auto g_0_zzzzzzz = pbuffer.data(sk_off + 35);

            /// set up bra offset for contr_buffer_xxpi

            const auto pi_off = idx_xxpi + (i * bcomps + j) * 84;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

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

#pragma omp simd aligned(cd_x,            \
                             g_0_xxxxxx,  \
                             g_0_xxxxxxx, \
                             g_0_xxxxxxy, \
                             g_0_xxxxxxz, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyy, \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxxzz, \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyy, \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxxzzz, \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyy, \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxxzzzz, \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyy, \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xxzzzzz, \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyy, \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_xzzzzzz, \
                             g_0_yyyyyy,  \
                             g_0_yyyyyz,  \
                             g_0_yyyyzz,  \
                             g_0_yyyzzz,  \
                             g_0_yyzzzz,  \
                             g_0_yzzzzz,  \
                             g_0_zzzzzz,  \
                             g_x_xxxxxx,  \
                             g_x_xxxxxy,  \
                             g_x_xxxxxz,  \
                             g_x_xxxxyy,  \
                             g_x_xxxxyz,  \
                             g_x_xxxxzz,  \
                             g_x_xxxyyy,  \
                             g_x_xxxyyz,  \
                             g_x_xxxyzz,  \
                             g_x_xxxzzz,  \
                             g_x_xxyyyy,  \
                             g_x_xxyyyz,  \
                             g_x_xxyyzz,  \
                             g_x_xxyzzz,  \
                             g_x_xxzzzz,  \
                             g_x_xyyyyy,  \
                             g_x_xyyyyz,  \
                             g_x_xyyyzz,  \
                             g_x_xyyzzz,  \
                             g_x_xyzzzz,  \
                             g_x_xzzzzz,  \
                             g_x_yyyyyy,  \
                             g_x_yyyyyz,  \
                             g_x_yyyyzz,  \
                             g_x_yyyzzz,  \
                             g_x_yyzzzz,  \
                             g_x_yzzzzz,  \
                             g_x_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xxxxxx[k] = -g_0_xxxxxx[k] * cd_x[k] + g_0_xxxxxxx[k];

                g_x_xxxxxy[k] = -g_0_xxxxxy[k] * cd_x[k] + g_0_xxxxxxy[k];

                g_x_xxxxxz[k] = -g_0_xxxxxz[k] * cd_x[k] + g_0_xxxxxxz[k];

                g_x_xxxxyy[k] = -g_0_xxxxyy[k] * cd_x[k] + g_0_xxxxxyy[k];

                g_x_xxxxyz[k] = -g_0_xxxxyz[k] * cd_x[k] + g_0_xxxxxyz[k];

                g_x_xxxxzz[k] = -g_0_xxxxzz[k] * cd_x[k] + g_0_xxxxxzz[k];

                g_x_xxxyyy[k] = -g_0_xxxyyy[k] * cd_x[k] + g_0_xxxxyyy[k];

                g_x_xxxyyz[k] = -g_0_xxxyyz[k] * cd_x[k] + g_0_xxxxyyz[k];

                g_x_xxxyzz[k] = -g_0_xxxyzz[k] * cd_x[k] + g_0_xxxxyzz[k];

                g_x_xxxzzz[k] = -g_0_xxxzzz[k] * cd_x[k] + g_0_xxxxzzz[k];

                g_x_xxyyyy[k] = -g_0_xxyyyy[k] * cd_x[k] + g_0_xxxyyyy[k];

                g_x_xxyyyz[k] = -g_0_xxyyyz[k] * cd_x[k] + g_0_xxxyyyz[k];

                g_x_xxyyzz[k] = -g_0_xxyyzz[k] * cd_x[k] + g_0_xxxyyzz[k];

                g_x_xxyzzz[k] = -g_0_xxyzzz[k] * cd_x[k] + g_0_xxxyzzz[k];

                g_x_xxzzzz[k] = -g_0_xxzzzz[k] * cd_x[k] + g_0_xxxzzzz[k];

                g_x_xyyyyy[k] = -g_0_xyyyyy[k] * cd_x[k] + g_0_xxyyyyy[k];

                g_x_xyyyyz[k] = -g_0_xyyyyz[k] * cd_x[k] + g_0_xxyyyyz[k];

                g_x_xyyyzz[k] = -g_0_xyyyzz[k] * cd_x[k] + g_0_xxyyyzz[k];

                g_x_xyyzzz[k] = -g_0_xyyzzz[k] * cd_x[k] + g_0_xxyyzzz[k];

                g_x_xyzzzz[k] = -g_0_xyzzzz[k] * cd_x[k] + g_0_xxyzzzz[k];

                g_x_xzzzzz[k] = -g_0_xzzzzz[k] * cd_x[k] + g_0_xxzzzzz[k];

                g_x_yyyyyy[k] = -g_0_yyyyyy[k] * cd_x[k] + g_0_xyyyyyy[k];

                g_x_yyyyyz[k] = -g_0_yyyyyz[k] * cd_x[k] + g_0_xyyyyyz[k];

                g_x_yyyyzz[k] = -g_0_yyyyzz[k] * cd_x[k] + g_0_xyyyyzz[k];

                g_x_yyyzzz[k] = -g_0_yyyzzz[k] * cd_x[k] + g_0_xyyyzzz[k];

                g_x_yyzzzz[k] = -g_0_yyzzzz[k] * cd_x[k] + g_0_xyyzzzz[k];

                g_x_yzzzzz[k] = -g_0_yzzzzz[k] * cd_x[k] + g_0_xyzzzzz[k];

                g_x_zzzzzz[k] = -g_0_zzzzzz[k] * cd_x[k] + g_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

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

#pragma omp simd aligned(cd_y,            \
                             g_0_xxxxxx,  \
                             g_0_xxxxxxy, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyy, \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyy, \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyy, \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyy, \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyy, \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_yyyyyy,  \
                             g_0_yyyyyyy, \
                             g_0_yyyyyyz, \
                             g_0_yyyyyz,  \
                             g_0_yyyyyzz, \
                             g_0_yyyyzz,  \
                             g_0_yyyyzzz, \
                             g_0_yyyzzz,  \
                             g_0_yyyzzzz, \
                             g_0_yyzzzz,  \
                             g_0_yyzzzzz, \
                             g_0_yzzzzz,  \
                             g_0_yzzzzzz, \
                             g_0_zzzzzz,  \
                             g_y_xxxxxx,  \
                             g_y_xxxxxy,  \
                             g_y_xxxxxz,  \
                             g_y_xxxxyy,  \
                             g_y_xxxxyz,  \
                             g_y_xxxxzz,  \
                             g_y_xxxyyy,  \
                             g_y_xxxyyz,  \
                             g_y_xxxyzz,  \
                             g_y_xxxzzz,  \
                             g_y_xxyyyy,  \
                             g_y_xxyyyz,  \
                             g_y_xxyyzz,  \
                             g_y_xxyzzz,  \
                             g_y_xxzzzz,  \
                             g_y_xyyyyy,  \
                             g_y_xyyyyz,  \
                             g_y_xyyyzz,  \
                             g_y_xyyzzz,  \
                             g_y_xyzzzz,  \
                             g_y_xzzzzz,  \
                             g_y_yyyyyy,  \
                             g_y_yyyyyz,  \
                             g_y_yyyyzz,  \
                             g_y_yyyzzz,  \
                             g_y_yyzzzz,  \
                             g_y_yzzzzz,  \
                             g_y_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xxxxxx[k] = -g_0_xxxxxx[k] * cd_y[k] + g_0_xxxxxxy[k];

                g_y_xxxxxy[k] = -g_0_xxxxxy[k] * cd_y[k] + g_0_xxxxxyy[k];

                g_y_xxxxxz[k] = -g_0_xxxxxz[k] * cd_y[k] + g_0_xxxxxyz[k];

                g_y_xxxxyy[k] = -g_0_xxxxyy[k] * cd_y[k] + g_0_xxxxyyy[k];

                g_y_xxxxyz[k] = -g_0_xxxxyz[k] * cd_y[k] + g_0_xxxxyyz[k];

                g_y_xxxxzz[k] = -g_0_xxxxzz[k] * cd_y[k] + g_0_xxxxyzz[k];

                g_y_xxxyyy[k] = -g_0_xxxyyy[k] * cd_y[k] + g_0_xxxyyyy[k];

                g_y_xxxyyz[k] = -g_0_xxxyyz[k] * cd_y[k] + g_0_xxxyyyz[k];

                g_y_xxxyzz[k] = -g_0_xxxyzz[k] * cd_y[k] + g_0_xxxyyzz[k];

                g_y_xxxzzz[k] = -g_0_xxxzzz[k] * cd_y[k] + g_0_xxxyzzz[k];

                g_y_xxyyyy[k] = -g_0_xxyyyy[k] * cd_y[k] + g_0_xxyyyyy[k];

                g_y_xxyyyz[k] = -g_0_xxyyyz[k] * cd_y[k] + g_0_xxyyyyz[k];

                g_y_xxyyzz[k] = -g_0_xxyyzz[k] * cd_y[k] + g_0_xxyyyzz[k];

                g_y_xxyzzz[k] = -g_0_xxyzzz[k] * cd_y[k] + g_0_xxyyzzz[k];

                g_y_xxzzzz[k] = -g_0_xxzzzz[k] * cd_y[k] + g_0_xxyzzzz[k];

                g_y_xyyyyy[k] = -g_0_xyyyyy[k] * cd_y[k] + g_0_xyyyyyy[k];

                g_y_xyyyyz[k] = -g_0_xyyyyz[k] * cd_y[k] + g_0_xyyyyyz[k];

                g_y_xyyyzz[k] = -g_0_xyyyzz[k] * cd_y[k] + g_0_xyyyyzz[k];

                g_y_xyyzzz[k] = -g_0_xyyzzz[k] * cd_y[k] + g_0_xyyyzzz[k];

                g_y_xyzzzz[k] = -g_0_xyzzzz[k] * cd_y[k] + g_0_xyyzzzz[k];

                g_y_xzzzzz[k] = -g_0_xzzzzz[k] * cd_y[k] + g_0_xyzzzzz[k];

                g_y_yyyyyy[k] = -g_0_yyyyyy[k] * cd_y[k] + g_0_yyyyyyy[k];

                g_y_yyyyyz[k] = -g_0_yyyyyz[k] * cd_y[k] + g_0_yyyyyyz[k];

                g_y_yyyyzz[k] = -g_0_yyyyzz[k] * cd_y[k] + g_0_yyyyyzz[k];

                g_y_yyyzzz[k] = -g_0_yyyzzz[k] * cd_y[k] + g_0_yyyyzzz[k];

                g_y_yyzzzz[k] = -g_0_yyzzzz[k] * cd_y[k] + g_0_yyyzzzz[k];

                g_y_yzzzzz[k] = -g_0_yzzzzz[k] * cd_y[k] + g_0_yyzzzzz[k];

                g_y_zzzzzz[k] = -g_0_zzzzzz[k] * cd_y[k] + g_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

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

#pragma omp simd aligned(cd_z,            \
                             g_0_xxxxxx,  \
                             g_0_xxxxxxz, \
                             g_0_xxxxxy,  \
                             g_0_xxxxxyz, \
                             g_0_xxxxxz,  \
                             g_0_xxxxxzz, \
                             g_0_xxxxyy,  \
                             g_0_xxxxyyz, \
                             g_0_xxxxyz,  \
                             g_0_xxxxyzz, \
                             g_0_xxxxzz,  \
                             g_0_xxxxzzz, \
                             g_0_xxxyyy,  \
                             g_0_xxxyyyz, \
                             g_0_xxxyyz,  \
                             g_0_xxxyyzz, \
                             g_0_xxxyzz,  \
                             g_0_xxxyzzz, \
                             g_0_xxxzzz,  \
                             g_0_xxxzzzz, \
                             g_0_xxyyyy,  \
                             g_0_xxyyyyz, \
                             g_0_xxyyyz,  \
                             g_0_xxyyyzz, \
                             g_0_xxyyzz,  \
                             g_0_xxyyzzz, \
                             g_0_xxyzzz,  \
                             g_0_xxyzzzz, \
                             g_0_xxzzzz,  \
                             g_0_xxzzzzz, \
                             g_0_xyyyyy,  \
                             g_0_xyyyyyz, \
                             g_0_xyyyyz,  \
                             g_0_xyyyyzz, \
                             g_0_xyyyzz,  \
                             g_0_xyyyzzz, \
                             g_0_xyyzzz,  \
                             g_0_xyyzzzz, \
                             g_0_xyzzzz,  \
                             g_0_xyzzzzz, \
                             g_0_xzzzzz,  \
                             g_0_xzzzzzz, \
                             g_0_yyyyyy,  \
                             g_0_yyyyyyz, \
                             g_0_yyyyyz,  \
                             g_0_yyyyyzz, \
                             g_0_yyyyzz,  \
                             g_0_yyyyzzz, \
                             g_0_yyyzzz,  \
                             g_0_yyyzzzz, \
                             g_0_yyzzzz,  \
                             g_0_yyzzzzz, \
                             g_0_yzzzzz,  \
                             g_0_yzzzzzz, \
                             g_0_zzzzzz,  \
                             g_0_zzzzzzz, \
                             g_z_xxxxxx,  \
                             g_z_xxxxxy,  \
                             g_z_xxxxxz,  \
                             g_z_xxxxyy,  \
                             g_z_xxxxyz,  \
                             g_z_xxxxzz,  \
                             g_z_xxxyyy,  \
                             g_z_xxxyyz,  \
                             g_z_xxxyzz,  \
                             g_z_xxxzzz,  \
                             g_z_xxyyyy,  \
                             g_z_xxyyyz,  \
                             g_z_xxyyzz,  \
                             g_z_xxyzzz,  \
                             g_z_xxzzzz,  \
                             g_z_xyyyyy,  \
                             g_z_xyyyyz,  \
                             g_z_xyyyzz,  \
                             g_z_xyyzzz,  \
                             g_z_xyzzzz,  \
                             g_z_xzzzzz,  \
                             g_z_yyyyyy,  \
                             g_z_yyyyyz,  \
                             g_z_yyyyzz,  \
                             g_z_yyyzzz,  \
                             g_z_yyzzzz,  \
                             g_z_yzzzzz,  \
                             g_z_zzzzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xxxxxx[k] = -g_0_xxxxxx[k] * cd_z[k] + g_0_xxxxxxz[k];

                g_z_xxxxxy[k] = -g_0_xxxxxy[k] * cd_z[k] + g_0_xxxxxyz[k];

                g_z_xxxxxz[k] = -g_0_xxxxxz[k] * cd_z[k] + g_0_xxxxxzz[k];

                g_z_xxxxyy[k] = -g_0_xxxxyy[k] * cd_z[k] + g_0_xxxxyyz[k];

                g_z_xxxxyz[k] = -g_0_xxxxyz[k] * cd_z[k] + g_0_xxxxyzz[k];

                g_z_xxxxzz[k] = -g_0_xxxxzz[k] * cd_z[k] + g_0_xxxxzzz[k];

                g_z_xxxyyy[k] = -g_0_xxxyyy[k] * cd_z[k] + g_0_xxxyyyz[k];

                g_z_xxxyyz[k] = -g_0_xxxyyz[k] * cd_z[k] + g_0_xxxyyzz[k];

                g_z_xxxyzz[k] = -g_0_xxxyzz[k] * cd_z[k] + g_0_xxxyzzz[k];

                g_z_xxxzzz[k] = -g_0_xxxzzz[k] * cd_z[k] + g_0_xxxzzzz[k];

                g_z_xxyyyy[k] = -g_0_xxyyyy[k] * cd_z[k] + g_0_xxyyyyz[k];

                g_z_xxyyyz[k] = -g_0_xxyyyz[k] * cd_z[k] + g_0_xxyyyzz[k];

                g_z_xxyyzz[k] = -g_0_xxyyzz[k] * cd_z[k] + g_0_xxyyzzz[k];

                g_z_xxyzzz[k] = -g_0_xxyzzz[k] * cd_z[k] + g_0_xxyzzzz[k];

                g_z_xxzzzz[k] = -g_0_xxzzzz[k] * cd_z[k] + g_0_xxzzzzz[k];

                g_z_xyyyyy[k] = -g_0_xyyyyy[k] * cd_z[k] + g_0_xyyyyyz[k];

                g_z_xyyyyz[k] = -g_0_xyyyyz[k] * cd_z[k] + g_0_xyyyyzz[k];

                g_z_xyyyzz[k] = -g_0_xyyyzz[k] * cd_z[k] + g_0_xyyyzzz[k];

                g_z_xyyzzz[k] = -g_0_xyyzzz[k] * cd_z[k] + g_0_xyyzzzz[k];

                g_z_xyzzzz[k] = -g_0_xyzzzz[k] * cd_z[k] + g_0_xyzzzzz[k];

                g_z_xzzzzz[k] = -g_0_xzzzzz[k] * cd_z[k] + g_0_xzzzzzz[k];

                g_z_yyyyyy[k] = -g_0_yyyyyy[k] * cd_z[k] + g_0_yyyyyyz[k];

                g_z_yyyyyz[k] = -g_0_yyyyyz[k] * cd_z[k] + g_0_yyyyyzz[k];

                g_z_yyyyzz[k] = -g_0_yyyyzz[k] * cd_z[k] + g_0_yyyyzzz[k];

                g_z_yyyzzz[k] = -g_0_yyyzzz[k] * cd_z[k] + g_0_yyyzzzz[k];

                g_z_yyzzzz[k] = -g_0_yyzzzz[k] * cd_z[k] + g_0_yyzzzzz[k];

                g_z_yzzzzz[k] = -g_0_yzzzzz[k] * cd_z[k] + g_0_yzzzzzz[k];

                g_z_zzzzzz[k] = -g_0_zzzzzz[k] * cd_z[k] + g_0_zzzzzzz[k];
            }
        }
    }
}

}  // namespace erirec

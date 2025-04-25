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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXDG.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xdg(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdg,
                                        const size_t idx_xpg,
                                        const size_t idx_geom_10_xpg,
                                        const size_t idx_geom_10_xph,
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
        /// Set up components of auxilary buffer : SPG

        const auto pg_off = idx_xpg + i * 45;

        auto g_x_xxxx = cbuffer.data(pg_off + 0);

        auto g_x_xxxy = cbuffer.data(pg_off + 1);

        auto g_x_xxxz = cbuffer.data(pg_off + 2);

        auto g_x_xxyy = cbuffer.data(pg_off + 3);

        auto g_x_xxyz = cbuffer.data(pg_off + 4);

        auto g_x_xxzz = cbuffer.data(pg_off + 5);

        auto g_x_xyyy = cbuffer.data(pg_off + 6);

        auto g_x_xyyz = cbuffer.data(pg_off + 7);

        auto g_x_xyzz = cbuffer.data(pg_off + 8);

        auto g_x_xzzz = cbuffer.data(pg_off + 9);

        auto g_x_yyyy = cbuffer.data(pg_off + 10);

        auto g_x_yyyz = cbuffer.data(pg_off + 11);

        auto g_x_yyzz = cbuffer.data(pg_off + 12);

        auto g_x_yzzz = cbuffer.data(pg_off + 13);

        auto g_x_zzzz = cbuffer.data(pg_off + 14);

        auto g_y_xxxx = cbuffer.data(pg_off + 15);

        auto g_y_xxxy = cbuffer.data(pg_off + 16);

        auto g_y_xxxz = cbuffer.data(pg_off + 17);

        auto g_y_xxyy = cbuffer.data(pg_off + 18);

        auto g_y_xxyz = cbuffer.data(pg_off + 19);

        auto g_y_xxzz = cbuffer.data(pg_off + 20);

        auto g_y_xyyy = cbuffer.data(pg_off + 21);

        auto g_y_xyyz = cbuffer.data(pg_off + 22);

        auto g_y_xyzz = cbuffer.data(pg_off + 23);

        auto g_y_xzzz = cbuffer.data(pg_off + 24);

        auto g_y_yyyy = cbuffer.data(pg_off + 25);

        auto g_y_yyyz = cbuffer.data(pg_off + 26);

        auto g_y_yyzz = cbuffer.data(pg_off + 27);

        auto g_y_yzzz = cbuffer.data(pg_off + 28);

        auto g_y_zzzz = cbuffer.data(pg_off + 29);

        auto g_z_xxxx = cbuffer.data(pg_off + 30);

        auto g_z_xxxy = cbuffer.data(pg_off + 31);

        auto g_z_xxxz = cbuffer.data(pg_off + 32);

        auto g_z_xxyy = cbuffer.data(pg_off + 33);

        auto g_z_xxyz = cbuffer.data(pg_off + 34);

        auto g_z_xxzz = cbuffer.data(pg_off + 35);

        auto g_z_xyyy = cbuffer.data(pg_off + 36);

        auto g_z_xyyz = cbuffer.data(pg_off + 37);

        auto g_z_xyzz = cbuffer.data(pg_off + 38);

        auto g_z_xzzz = cbuffer.data(pg_off + 39);

        auto g_z_yyyy = cbuffer.data(pg_off + 40);

        auto g_z_yyyz = cbuffer.data(pg_off + 41);

        auto g_z_yyzz = cbuffer.data(pg_off + 42);

        auto g_z_yzzz = cbuffer.data(pg_off + 43);

        auto g_z_zzzz = cbuffer.data(pg_off + 44);

        /// Set up components of auxilary buffer : SPG

        const auto pg_geom_10_off = idx_geom_10_xpg + i * 45;

        auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 9);

        auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 10);

        auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 11);

        auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 12);

        auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 13);

        auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 14);

        auto g_x_0_y_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 15);

        auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 16);

        auto g_x_0_y_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 17);

        auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 18);

        auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 19);

        auto g_x_0_y_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 20);

        auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 21);

        auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 22);

        auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 23);

        auto g_x_0_y_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 24);

        auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 25);

        auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 26);

        auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 27);

        auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 28);

        auto g_x_0_y_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 29);

        auto g_x_0_z_xxxx = cbuffer.data(pg_geom_10_off + 0 * acomps + 30);

        auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 0 * acomps + 31);

        auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 0 * acomps + 32);

        auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 33);

        auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 34);

        auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 35);

        auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 36);

        auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 37);

        auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 38);

        auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 39);

        auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 0 * acomps + 40);

        auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 0 * acomps + 41);

        auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 42);

        auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 43);

        auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 0 * acomps + 44);

        auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 0);

        auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 1);

        auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 2);

        auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 3);

        auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 4);

        auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 5);

        auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 6);

        auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 7);

        auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 8);

        auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 9);

        auto g_y_0_x_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 10);

        auto g_y_0_x_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 11);

        auto g_y_0_x_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 12);

        auto g_y_0_x_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 13);

        auto g_y_0_x_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 14);

        auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 15);

        auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 16);

        auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 17);

        auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 18);

        auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 19);

        auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 20);

        auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 21);

        auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 22);

        auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 23);

        auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 24);

        auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 25);

        auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 26);

        auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 27);

        auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 28);

        auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 29);

        auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 45 * acomps + 30);

        auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 45 * acomps + 31);

        auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 45 * acomps + 32);

        auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 33);

        auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 34);

        auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 35);

        auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 36);

        auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 37);

        auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 38);

        auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 39);

        auto g_y_0_z_yyyy = cbuffer.data(pg_geom_10_off + 45 * acomps + 40);

        auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 45 * acomps + 41);

        auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 42);

        auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 43);

        auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 45 * acomps + 44);

        auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 0);

        auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 1);

        auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 2);

        auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 3);

        auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 4);

        auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 5);

        auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 6);

        auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 7);

        auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 8);

        auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 9);

        auto g_z_0_x_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 10);

        auto g_z_0_x_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 11);

        auto g_z_0_x_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 12);

        auto g_z_0_x_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 13);

        auto g_z_0_x_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 14);

        auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 15);

        auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 16);

        auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 17);

        auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 18);

        auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 19);

        auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 20);

        auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 21);

        auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 22);

        auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 23);

        auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 24);

        auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 25);

        auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 26);

        auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 27);

        auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 28);

        auto g_z_0_y_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 29);

        auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 90 * acomps + 30);

        auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 90 * acomps + 31);

        auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 90 * acomps + 32);

        auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 33);

        auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 34);

        auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 35);

        auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 36);

        auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 37);

        auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 38);

        auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 39);

        auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 90 * acomps + 40);

        auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 90 * acomps + 41);

        auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 42);

        auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 43);

        auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 90 * acomps + 44);

        /// Set up components of auxilary buffer : SPH

        const auto ph_geom_10_off = idx_geom_10_xph + i * 63;

        auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 9);

        auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 10);

        auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 11);

        auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 12);

        auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 13);

        auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 14);

        auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 15);

        auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 16);

        auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 17);

        auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 18);

        auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 19);

        auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 20);

        auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 22);

        auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 24);

        auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 25);

        auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 27);

        auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 28);

        auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 29);

        auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 31);

        auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 32);

        auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 33);

        auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 34);

        auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 36);

        auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 37);

        auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 38);

        auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 39);

        auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 40);

        auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps + 43);

        auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps + 44);

        auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 45);

        auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 46);

        auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 47);

        auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 48);

        auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 49);

        auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 50);

        auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 51);

        auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 52);

        auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 53);

        auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 54);

        auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 55);

        auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 56);

        auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps + 57);

        auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps + 58);

        auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 59);

        auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 60);

        auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 61);

        auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps + 62);

        auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 0);

        auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 1);

        auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 2);

        auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 3);

        auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 4);

        auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 5);

        auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 6);

        auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 7);

        auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 8);

        auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 9);

        auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 10);

        auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 11);

        auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 12);

        auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 13);

        auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 14);

        auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 21);

        auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 22);

        auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 23);

        auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 24);

        auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 25);

        auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 26);

        auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 27);

        auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 28);

        auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 29);

        auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 30);

        auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 31);

        auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 32);

        auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 33);

        auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 34);

        auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 35);

        auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 36);

        auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 37);

        auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 38);

        auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 39);

        auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 40);

        auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 41);

        auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps + 42);

        auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps + 43);

        auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps + 44);

        auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 45);

        auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 46);

        auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 47);

        auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 48);

        auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 49);

        auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 50);

        auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 51);

        auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps + 52);

        auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 53);

        auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 54);

        auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 55);

        auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 56);

        auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps + 58);

        auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 59);

        auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 60);

        auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 61);

        auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps + 62);

        auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 0);

        auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 1);

        auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 2);

        auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 3);

        auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 4);

        auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 5);

        auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 6);

        auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 7);

        auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 8);

        auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 9);

        auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 10);

        auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 11);

        auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 12);

        auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 13);

        auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 14);

        auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 21);

        auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 22);

        auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 23);

        auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 24);

        auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 25);

        auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 26);

        auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 27);

        auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 28);

        auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 29);

        auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 30);

        auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 31);

        auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 32);

        auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 33);

        auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 34);

        auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 35);

        auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 36);

        auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 37);

        auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 38);

        auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 39);

        auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 40);

        auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps + 42);

        auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps + 43);

        auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps + 44);

        auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 45);

        auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 46);

        auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 47);

        auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 48);

        auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 49);

        auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 50);

        auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 51);

        auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 52);

        auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 53);

        auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 54);

        auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 55);

        auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 56);

        auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps + 57);

        auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps + 58);

        auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 59);

        auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 60);

        auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 61);

        auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps + 62);

        /// set up bra offset for contr_buffer_xxdg

        const auto dg_geom_10_off = idx_geom_10_xdg + i * 90;

        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_x_0_x_xxxx, g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxy, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxxzz, g_x_0_x_xxyy, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xxzzz, g_x_0_x_xyyy, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_xzzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz, g_x_0_xx_xxxx, g_x_0_xx_xxxy, g_x_0_xx_xxxz, g_x_0_xx_xxyy, g_x_0_xx_xxyz, g_x_0_xx_xxzz, g_x_0_xx_xyyy, g_x_0_xx_xyyz, g_x_0_xx_xyzz, g_x_0_xx_xzzz, g_x_0_xx_yyyy, g_x_0_xx_yyyz, g_x_0_xx_yyzz, g_x_0_xx_yzzz, g_x_0_xx_zzzz, g_x_xxxx, g_x_xxxy, g_x_xxxz, g_x_xxyy, g_x_xxyz, g_x_xxzz, g_x_xyyy, g_x_xyyz, g_x_xyzz, g_x_xzzz, g_x_yyyy, g_x_yyyz, g_x_yyzz, g_x_yzzz, g_x_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xx_xxxx[k] = -g_x_xxxx[k] - g_x_0_x_xxxx[k] * cd_x[k] + g_x_0_x_xxxxx[k];

            g_x_0_xx_xxxy[k] = -g_x_xxxy[k] - g_x_0_x_xxxy[k] * cd_x[k] + g_x_0_x_xxxxy[k];

            g_x_0_xx_xxxz[k] = -g_x_xxxz[k] - g_x_0_x_xxxz[k] * cd_x[k] + g_x_0_x_xxxxz[k];

            g_x_0_xx_xxyy[k] = -g_x_xxyy[k] - g_x_0_x_xxyy[k] * cd_x[k] + g_x_0_x_xxxyy[k];

            g_x_0_xx_xxyz[k] = -g_x_xxyz[k] - g_x_0_x_xxyz[k] * cd_x[k] + g_x_0_x_xxxyz[k];

            g_x_0_xx_xxzz[k] = -g_x_xxzz[k] - g_x_0_x_xxzz[k] * cd_x[k] + g_x_0_x_xxxzz[k];

            g_x_0_xx_xyyy[k] = -g_x_xyyy[k] - g_x_0_x_xyyy[k] * cd_x[k] + g_x_0_x_xxyyy[k];

            g_x_0_xx_xyyz[k] = -g_x_xyyz[k] - g_x_0_x_xyyz[k] * cd_x[k] + g_x_0_x_xxyyz[k];

            g_x_0_xx_xyzz[k] = -g_x_xyzz[k] - g_x_0_x_xyzz[k] * cd_x[k] + g_x_0_x_xxyzz[k];

            g_x_0_xx_xzzz[k] = -g_x_xzzz[k] - g_x_0_x_xzzz[k] * cd_x[k] + g_x_0_x_xxzzz[k];

            g_x_0_xx_yyyy[k] = -g_x_yyyy[k] - g_x_0_x_yyyy[k] * cd_x[k] + g_x_0_x_xyyyy[k];

            g_x_0_xx_yyyz[k] = -g_x_yyyz[k] - g_x_0_x_yyyz[k] * cd_x[k] + g_x_0_x_xyyyz[k];

            g_x_0_xx_yyzz[k] = -g_x_yyzz[k] - g_x_0_x_yyzz[k] * cd_x[k] + g_x_0_x_xyyzz[k];

            g_x_0_xx_yzzz[k] = -g_x_yzzz[k] - g_x_0_x_yzzz[k] * cd_x[k] + g_x_0_x_xyzzz[k];

            g_x_0_xx_zzzz[k] = -g_x_zzzz[k] - g_x_0_x_zzzz[k] * cd_x[k] + g_x_0_x_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_x_xxxx, g_x_0_x_xxxxy, g_x_0_x_xxxy, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyz, g_x_0_x_yyyzz, g_x_0_x_yyzz, g_x_0_x_yyzzz, g_x_0_x_yzzz, g_x_0_x_yzzzz, g_x_0_x_zzzz, g_x_0_xy_xxxx, g_x_0_xy_xxxy, g_x_0_xy_xxxz, g_x_0_xy_xxyy, g_x_0_xy_xxyz, g_x_0_xy_xxzz, g_x_0_xy_xyyy, g_x_0_xy_xyyz, g_x_0_xy_xyzz, g_x_0_xy_xzzz, g_x_0_xy_yyyy, g_x_0_xy_yyyz, g_x_0_xy_yyzz, g_x_0_xy_yzzz, g_x_0_xy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xy_xxxx[k] = -g_x_0_x_xxxx[k] * cd_y[k] + g_x_0_x_xxxxy[k];

            g_x_0_xy_xxxy[k] = -g_x_0_x_xxxy[k] * cd_y[k] + g_x_0_x_xxxyy[k];

            g_x_0_xy_xxxz[k] = -g_x_0_x_xxxz[k] * cd_y[k] + g_x_0_x_xxxyz[k];

            g_x_0_xy_xxyy[k] = -g_x_0_x_xxyy[k] * cd_y[k] + g_x_0_x_xxyyy[k];

            g_x_0_xy_xxyz[k] = -g_x_0_x_xxyz[k] * cd_y[k] + g_x_0_x_xxyyz[k];

            g_x_0_xy_xxzz[k] = -g_x_0_x_xxzz[k] * cd_y[k] + g_x_0_x_xxyzz[k];

            g_x_0_xy_xyyy[k] = -g_x_0_x_xyyy[k] * cd_y[k] + g_x_0_x_xyyyy[k];

            g_x_0_xy_xyyz[k] = -g_x_0_x_xyyz[k] * cd_y[k] + g_x_0_x_xyyyz[k];

            g_x_0_xy_xyzz[k] = -g_x_0_x_xyzz[k] * cd_y[k] + g_x_0_x_xyyzz[k];

            g_x_0_xy_xzzz[k] = -g_x_0_x_xzzz[k] * cd_y[k] + g_x_0_x_xyzzz[k];

            g_x_0_xy_yyyy[k] = -g_x_0_x_yyyy[k] * cd_y[k] + g_x_0_x_yyyyy[k];

            g_x_0_xy_yyyz[k] = -g_x_0_x_yyyz[k] * cd_y[k] + g_x_0_x_yyyyz[k];

            g_x_0_xy_yyzz[k] = -g_x_0_x_yyzz[k] * cd_y[k] + g_x_0_x_yyyzz[k];

            g_x_0_xy_yzzz[k] = -g_x_0_x_yzzz[k] * cd_y[k] + g_x_0_x_yyzzz[k];

            g_x_0_xy_zzzz[k] = -g_x_0_x_zzzz[k] * cd_y[k] + g_x_0_x_yzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_x_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 44);

        #pragma omp simd aligned(cd_z, g_x_0_x_xxxx, g_x_0_x_xxxxz, g_x_0_x_xxxy, g_x_0_x_xxxyz, g_x_0_x_xxxz, g_x_0_x_xxxzz, g_x_0_x_xxyy, g_x_0_x_xxyyz, g_x_0_x_xxyz, g_x_0_x_xxyzz, g_x_0_x_xxzz, g_x_0_x_xxzzz, g_x_0_x_xyyy, g_x_0_x_xyyyz, g_x_0_x_xyyz, g_x_0_x_xyyzz, g_x_0_x_xyzz, g_x_0_x_xyzzz, g_x_0_x_xzzz, g_x_0_x_xzzzz, g_x_0_x_yyyy, g_x_0_x_yyyyz, g_x_0_x_yyyz, g_x_0_x_yyyzz, g_x_0_x_yyzz, g_x_0_x_yyzzz, g_x_0_x_yzzz, g_x_0_x_yzzzz, g_x_0_x_zzzz, g_x_0_x_zzzzz, g_x_0_xz_xxxx, g_x_0_xz_xxxy, g_x_0_xz_xxxz, g_x_0_xz_xxyy, g_x_0_xz_xxyz, g_x_0_xz_xxzz, g_x_0_xz_xyyy, g_x_0_xz_xyyz, g_x_0_xz_xyzz, g_x_0_xz_xzzz, g_x_0_xz_yyyy, g_x_0_xz_yyyz, g_x_0_xz_yyzz, g_x_0_xz_yzzz, g_x_0_xz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xz_xxxx[k] = -g_x_0_x_xxxx[k] * cd_z[k] + g_x_0_x_xxxxz[k];

            g_x_0_xz_xxxy[k] = -g_x_0_x_xxxy[k] * cd_z[k] + g_x_0_x_xxxyz[k];

            g_x_0_xz_xxxz[k] = -g_x_0_x_xxxz[k] * cd_z[k] + g_x_0_x_xxxzz[k];

            g_x_0_xz_xxyy[k] = -g_x_0_x_xxyy[k] * cd_z[k] + g_x_0_x_xxyyz[k];

            g_x_0_xz_xxyz[k] = -g_x_0_x_xxyz[k] * cd_z[k] + g_x_0_x_xxyzz[k];

            g_x_0_xz_xxzz[k] = -g_x_0_x_xxzz[k] * cd_z[k] + g_x_0_x_xxzzz[k];

            g_x_0_xz_xyyy[k] = -g_x_0_x_xyyy[k] * cd_z[k] + g_x_0_x_xyyyz[k];

            g_x_0_xz_xyyz[k] = -g_x_0_x_xyyz[k] * cd_z[k] + g_x_0_x_xyyzz[k];

            g_x_0_xz_xyzz[k] = -g_x_0_x_xyzz[k] * cd_z[k] + g_x_0_x_xyzzz[k];

            g_x_0_xz_xzzz[k] = -g_x_0_x_xzzz[k] * cd_z[k] + g_x_0_x_xzzzz[k];

            g_x_0_xz_yyyy[k] = -g_x_0_x_yyyy[k] * cd_z[k] + g_x_0_x_yyyyz[k];

            g_x_0_xz_yyyz[k] = -g_x_0_x_yyyz[k] * cd_z[k] + g_x_0_x_yyyzz[k];

            g_x_0_xz_yyzz[k] = -g_x_0_x_yyzz[k] * cd_z[k] + g_x_0_x_yyzzz[k];

            g_x_0_xz_yzzz[k] = -g_x_0_x_yzzz[k] * cd_z[k] + g_x_0_x_yzzzz[k];

            g_x_0_xz_zzzz[k] = -g_x_0_x_zzzz[k] * cd_z[k] + g_x_0_x_zzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_y, g_x_0_y_xxxx, g_x_0_y_xxxxy, g_x_0_y_xxxy, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyz, g_x_0_y_xxyzz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyz, g_x_0_y_xyyzz, g_x_0_y_xyzz, g_x_0_y_xyzzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyz, g_x_0_y_yyyzz, g_x_0_y_yyzz, g_x_0_y_yyzzz, g_x_0_y_yzzz, g_x_0_y_yzzzz, g_x_0_y_zzzz, g_x_0_yy_xxxx, g_x_0_yy_xxxy, g_x_0_yy_xxxz, g_x_0_yy_xxyy, g_x_0_yy_xxyz, g_x_0_yy_xxzz, g_x_0_yy_xyyy, g_x_0_yy_xyyz, g_x_0_yy_xyzz, g_x_0_yy_xzzz, g_x_0_yy_yyyy, g_x_0_yy_yyyz, g_x_0_yy_yyzz, g_x_0_yy_yzzz, g_x_0_yy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yy_xxxx[k] = -g_x_0_y_xxxx[k] * cd_y[k] + g_x_0_y_xxxxy[k];

            g_x_0_yy_xxxy[k] = -g_x_0_y_xxxy[k] * cd_y[k] + g_x_0_y_xxxyy[k];

            g_x_0_yy_xxxz[k] = -g_x_0_y_xxxz[k] * cd_y[k] + g_x_0_y_xxxyz[k];

            g_x_0_yy_xxyy[k] = -g_x_0_y_xxyy[k] * cd_y[k] + g_x_0_y_xxyyy[k];

            g_x_0_yy_xxyz[k] = -g_x_0_y_xxyz[k] * cd_y[k] + g_x_0_y_xxyyz[k];

            g_x_0_yy_xxzz[k] = -g_x_0_y_xxzz[k] * cd_y[k] + g_x_0_y_xxyzz[k];

            g_x_0_yy_xyyy[k] = -g_x_0_y_xyyy[k] * cd_y[k] + g_x_0_y_xyyyy[k];

            g_x_0_yy_xyyz[k] = -g_x_0_y_xyyz[k] * cd_y[k] + g_x_0_y_xyyyz[k];

            g_x_0_yy_xyzz[k] = -g_x_0_y_xyzz[k] * cd_y[k] + g_x_0_y_xyyzz[k];

            g_x_0_yy_xzzz[k] = -g_x_0_y_xzzz[k] * cd_y[k] + g_x_0_y_xyzzz[k];

            g_x_0_yy_yyyy[k] = -g_x_0_y_yyyy[k] * cd_y[k] + g_x_0_y_yyyyy[k];

            g_x_0_yy_yyyz[k] = -g_x_0_y_yyyz[k] * cd_y[k] + g_x_0_y_yyyyz[k];

            g_x_0_yy_yyzz[k] = -g_x_0_y_yyzz[k] * cd_y[k] + g_x_0_y_yyyzz[k];

            g_x_0_yy_yzzz[k] = -g_x_0_y_yzzz[k] * cd_y[k] + g_x_0_y_yyzzz[k];

            g_x_0_yy_zzzz[k] = -g_x_0_y_zzzz[k] * cd_y[k] + g_x_0_y_yzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_x_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 74);

        #pragma omp simd aligned(cd_y, g_x_0_yz_xxxx, g_x_0_yz_xxxy, g_x_0_yz_xxxz, g_x_0_yz_xxyy, g_x_0_yz_xxyz, g_x_0_yz_xxzz, g_x_0_yz_xyyy, g_x_0_yz_xyyz, g_x_0_yz_xyzz, g_x_0_yz_xzzz, g_x_0_yz_yyyy, g_x_0_yz_yyyz, g_x_0_yz_yyzz, g_x_0_yz_yzzz, g_x_0_yz_zzzz, g_x_0_z_xxxx, g_x_0_z_xxxxy, g_x_0_z_xxxy, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyz, g_x_0_z_xxyzz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyz, g_x_0_z_xyyzz, g_x_0_z_xyzz, g_x_0_z_xyzzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyz, g_x_0_z_yyyzz, g_x_0_z_yyzz, g_x_0_z_yyzzz, g_x_0_z_yzzz, g_x_0_z_yzzzz, g_x_0_z_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yz_xxxx[k] = -g_x_0_z_xxxx[k] * cd_y[k] + g_x_0_z_xxxxy[k];

            g_x_0_yz_xxxy[k] = -g_x_0_z_xxxy[k] * cd_y[k] + g_x_0_z_xxxyy[k];

            g_x_0_yz_xxxz[k] = -g_x_0_z_xxxz[k] * cd_y[k] + g_x_0_z_xxxyz[k];

            g_x_0_yz_xxyy[k] = -g_x_0_z_xxyy[k] * cd_y[k] + g_x_0_z_xxyyy[k];

            g_x_0_yz_xxyz[k] = -g_x_0_z_xxyz[k] * cd_y[k] + g_x_0_z_xxyyz[k];

            g_x_0_yz_xxzz[k] = -g_x_0_z_xxzz[k] * cd_y[k] + g_x_0_z_xxyzz[k];

            g_x_0_yz_xyyy[k] = -g_x_0_z_xyyy[k] * cd_y[k] + g_x_0_z_xyyyy[k];

            g_x_0_yz_xyyz[k] = -g_x_0_z_xyyz[k] * cd_y[k] + g_x_0_z_xyyyz[k];

            g_x_0_yz_xyzz[k] = -g_x_0_z_xyzz[k] * cd_y[k] + g_x_0_z_xyyzz[k];

            g_x_0_yz_xzzz[k] = -g_x_0_z_xzzz[k] * cd_y[k] + g_x_0_z_xyzzz[k];

            g_x_0_yz_yyyy[k] = -g_x_0_z_yyyy[k] * cd_y[k] + g_x_0_z_yyyyy[k];

            g_x_0_yz_yyyz[k] = -g_x_0_z_yyyz[k] * cd_y[k] + g_x_0_z_yyyyz[k];

            g_x_0_yz_yyzz[k] = -g_x_0_z_yyzz[k] * cd_y[k] + g_x_0_z_yyyzz[k];

            g_x_0_yz_yzzz[k] = -g_x_0_z_yzzz[k] * cd_y[k] + g_x_0_z_yyzzz[k];

            g_x_0_yz_zzzz[k] = -g_x_0_z_zzzz[k] * cd_y[k] + g_x_0_z_yzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_x_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 83);

        auto g_x_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 0 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_x_0_z_xxxx, g_x_0_z_xxxxz, g_x_0_z_xxxy, g_x_0_z_xxxyz, g_x_0_z_xxxz, g_x_0_z_xxxzz, g_x_0_z_xxyy, g_x_0_z_xxyyz, g_x_0_z_xxyz, g_x_0_z_xxyzz, g_x_0_z_xxzz, g_x_0_z_xxzzz, g_x_0_z_xyyy, g_x_0_z_xyyyz, g_x_0_z_xyyz, g_x_0_z_xyyzz, g_x_0_z_xyzz, g_x_0_z_xyzzz, g_x_0_z_xzzz, g_x_0_z_xzzzz, g_x_0_z_yyyy, g_x_0_z_yyyyz, g_x_0_z_yyyz, g_x_0_z_yyyzz, g_x_0_z_yyzz, g_x_0_z_yyzzz, g_x_0_z_yzzz, g_x_0_z_yzzzz, g_x_0_z_zzzz, g_x_0_z_zzzzz, g_x_0_zz_xxxx, g_x_0_zz_xxxy, g_x_0_zz_xxxz, g_x_0_zz_xxyy, g_x_0_zz_xxyz, g_x_0_zz_xxzz, g_x_0_zz_xyyy, g_x_0_zz_xyyz, g_x_0_zz_xyzz, g_x_0_zz_xzzz, g_x_0_zz_yyyy, g_x_0_zz_yyyz, g_x_0_zz_yyzz, g_x_0_zz_yzzz, g_x_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zz_xxxx[k] = -g_x_0_z_xxxx[k] * cd_z[k] + g_x_0_z_xxxxz[k];

            g_x_0_zz_xxxy[k] = -g_x_0_z_xxxy[k] * cd_z[k] + g_x_0_z_xxxyz[k];

            g_x_0_zz_xxxz[k] = -g_x_0_z_xxxz[k] * cd_z[k] + g_x_0_z_xxxzz[k];

            g_x_0_zz_xxyy[k] = -g_x_0_z_xxyy[k] * cd_z[k] + g_x_0_z_xxyyz[k];

            g_x_0_zz_xxyz[k] = -g_x_0_z_xxyz[k] * cd_z[k] + g_x_0_z_xxyzz[k];

            g_x_0_zz_xxzz[k] = -g_x_0_z_xxzz[k] * cd_z[k] + g_x_0_z_xxzzz[k];

            g_x_0_zz_xyyy[k] = -g_x_0_z_xyyy[k] * cd_z[k] + g_x_0_z_xyyyz[k];

            g_x_0_zz_xyyz[k] = -g_x_0_z_xyyz[k] * cd_z[k] + g_x_0_z_xyyzz[k];

            g_x_0_zz_xyzz[k] = -g_x_0_z_xyzz[k] * cd_z[k] + g_x_0_z_xyzzz[k];

            g_x_0_zz_xzzz[k] = -g_x_0_z_xzzz[k] * cd_z[k] + g_x_0_z_xzzzz[k];

            g_x_0_zz_yyyy[k] = -g_x_0_z_yyyy[k] * cd_z[k] + g_x_0_z_yyyyz[k];

            g_x_0_zz_yyyz[k] = -g_x_0_z_yyyz[k] * cd_z[k] + g_x_0_z_yyyzz[k];

            g_x_0_zz_yyzz[k] = -g_x_0_z_yyzz[k] * cd_z[k] + g_x_0_z_yyzzz[k];

            g_x_0_zz_yzzz[k] = -g_x_0_z_yzzz[k] * cd_z[k] + g_x_0_z_yzzzz[k];

            g_x_0_zz_zzzz[k] = -g_x_0_z_zzzz[k] * cd_z[k] + g_x_0_z_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 0);

        auto g_y_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 1);

        auto g_y_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 2);

        auto g_y_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 3);

        auto g_y_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 4);

        auto g_y_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 5);

        auto g_y_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 6);

        auto g_y_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 7);

        auto g_y_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 8);

        auto g_y_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 9);

        auto g_y_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 10);

        auto g_y_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 11);

        auto g_y_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 12);

        auto g_y_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 13);

        auto g_y_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_y_0_x_xxxx, g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxy, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxz, g_y_0_x_xxxzz, g_y_0_x_xxyy, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyz, g_y_0_x_xxyzz, g_y_0_x_xxzz, g_y_0_x_xxzzz, g_y_0_x_xyyy, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyz, g_y_0_x_xyyzz, g_y_0_x_xyzz, g_y_0_x_xyzzz, g_y_0_x_xzzz, g_y_0_x_xzzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz, g_y_0_xx_xxxx, g_y_0_xx_xxxy, g_y_0_xx_xxxz, g_y_0_xx_xxyy, g_y_0_xx_xxyz, g_y_0_xx_xxzz, g_y_0_xx_xyyy, g_y_0_xx_xyyz, g_y_0_xx_xyzz, g_y_0_xx_xzzz, g_y_0_xx_yyyy, g_y_0_xx_yyyz, g_y_0_xx_yyzz, g_y_0_xx_yzzz, g_y_0_xx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xx_xxxx[k] = -g_y_0_x_xxxx[k] * cd_x[k] + g_y_0_x_xxxxx[k];

            g_y_0_xx_xxxy[k] = -g_y_0_x_xxxy[k] * cd_x[k] + g_y_0_x_xxxxy[k];

            g_y_0_xx_xxxz[k] = -g_y_0_x_xxxz[k] * cd_x[k] + g_y_0_x_xxxxz[k];

            g_y_0_xx_xxyy[k] = -g_y_0_x_xxyy[k] * cd_x[k] + g_y_0_x_xxxyy[k];

            g_y_0_xx_xxyz[k] = -g_y_0_x_xxyz[k] * cd_x[k] + g_y_0_x_xxxyz[k];

            g_y_0_xx_xxzz[k] = -g_y_0_x_xxzz[k] * cd_x[k] + g_y_0_x_xxxzz[k];

            g_y_0_xx_xyyy[k] = -g_y_0_x_xyyy[k] * cd_x[k] + g_y_0_x_xxyyy[k];

            g_y_0_xx_xyyz[k] = -g_y_0_x_xyyz[k] * cd_x[k] + g_y_0_x_xxyyz[k];

            g_y_0_xx_xyzz[k] = -g_y_0_x_xyzz[k] * cd_x[k] + g_y_0_x_xxyzz[k];

            g_y_0_xx_xzzz[k] = -g_y_0_x_xzzz[k] * cd_x[k] + g_y_0_x_xxzzz[k];

            g_y_0_xx_yyyy[k] = -g_y_0_x_yyyy[k] * cd_x[k] + g_y_0_x_xyyyy[k];

            g_y_0_xx_yyyz[k] = -g_y_0_x_yyyz[k] * cd_x[k] + g_y_0_x_xyyyz[k];

            g_y_0_xx_yyzz[k] = -g_y_0_x_yyzz[k] * cd_x[k] + g_y_0_x_xyyzz[k];

            g_y_0_xx_yzzz[k] = -g_y_0_x_yzzz[k] * cd_x[k] + g_y_0_x_xyzzz[k];

            g_y_0_xx_zzzz[k] = -g_y_0_x_zzzz[k] * cd_x[k] + g_y_0_x_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 15);

        auto g_y_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 16);

        auto g_y_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 17);

        auto g_y_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 18);

        auto g_y_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 19);

        auto g_y_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 20);

        auto g_y_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 21);

        auto g_y_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 22);

        auto g_y_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 23);

        auto g_y_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 24);

        auto g_y_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 25);

        auto g_y_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 26);

        auto g_y_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 27);

        auto g_y_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 28);

        auto g_y_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xy_xxxx, g_y_0_xy_xxxy, g_y_0_xy_xxxz, g_y_0_xy_xxyy, g_y_0_xy_xxyz, g_y_0_xy_xxzz, g_y_0_xy_xyyy, g_y_0_xy_xyyz, g_y_0_xy_xyzz, g_y_0_xy_xzzz, g_y_0_xy_yyyy, g_y_0_xy_yyyz, g_y_0_xy_yyzz, g_y_0_xy_yzzz, g_y_0_xy_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxy, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxxzz, g_y_0_y_xxyy, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xxzzz, g_y_0_y_xyyy, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_xzzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xy_xxxx[k] = -g_y_0_y_xxxx[k] * cd_x[k] + g_y_0_y_xxxxx[k];

            g_y_0_xy_xxxy[k] = -g_y_0_y_xxxy[k] * cd_x[k] + g_y_0_y_xxxxy[k];

            g_y_0_xy_xxxz[k] = -g_y_0_y_xxxz[k] * cd_x[k] + g_y_0_y_xxxxz[k];

            g_y_0_xy_xxyy[k] = -g_y_0_y_xxyy[k] * cd_x[k] + g_y_0_y_xxxyy[k];

            g_y_0_xy_xxyz[k] = -g_y_0_y_xxyz[k] * cd_x[k] + g_y_0_y_xxxyz[k];

            g_y_0_xy_xxzz[k] = -g_y_0_y_xxzz[k] * cd_x[k] + g_y_0_y_xxxzz[k];

            g_y_0_xy_xyyy[k] = -g_y_0_y_xyyy[k] * cd_x[k] + g_y_0_y_xxyyy[k];

            g_y_0_xy_xyyz[k] = -g_y_0_y_xyyz[k] * cd_x[k] + g_y_0_y_xxyyz[k];

            g_y_0_xy_xyzz[k] = -g_y_0_y_xyzz[k] * cd_x[k] + g_y_0_y_xxyzz[k];

            g_y_0_xy_xzzz[k] = -g_y_0_y_xzzz[k] * cd_x[k] + g_y_0_y_xxzzz[k];

            g_y_0_xy_yyyy[k] = -g_y_0_y_yyyy[k] * cd_x[k] + g_y_0_y_xyyyy[k];

            g_y_0_xy_yyyz[k] = -g_y_0_y_yyyz[k] * cd_x[k] + g_y_0_y_xyyyz[k];

            g_y_0_xy_yyzz[k] = -g_y_0_y_yyzz[k] * cd_x[k] + g_y_0_y_xyyzz[k];

            g_y_0_xy_yzzz[k] = -g_y_0_y_yzzz[k] * cd_x[k] + g_y_0_y_xyzzz[k];

            g_y_0_xy_zzzz[k] = -g_y_0_y_zzzz[k] * cd_x[k] + g_y_0_y_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_y_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 30);

        auto g_y_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 31);

        auto g_y_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 32);

        auto g_y_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 33);

        auto g_y_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 34);

        auto g_y_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 35);

        auto g_y_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 36);

        auto g_y_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 37);

        auto g_y_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 38);

        auto g_y_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 39);

        auto g_y_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 40);

        auto g_y_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 41);

        auto g_y_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 42);

        auto g_y_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 43);

        auto g_y_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_y_0_xz_xxxx, g_y_0_xz_xxxy, g_y_0_xz_xxxz, g_y_0_xz_xxyy, g_y_0_xz_xxyz, g_y_0_xz_xxzz, g_y_0_xz_xyyy, g_y_0_xz_xyyz, g_y_0_xz_xyzz, g_y_0_xz_xzzz, g_y_0_xz_yyyy, g_y_0_xz_yyyz, g_y_0_xz_yyzz, g_y_0_xz_yzzz, g_y_0_xz_zzzz, g_y_0_z_xxxx, g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxy, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxz, g_y_0_z_xxxzz, g_y_0_z_xxyy, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyz, g_y_0_z_xxyzz, g_y_0_z_xxzz, g_y_0_z_xxzzz, g_y_0_z_xyyy, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyz, g_y_0_z_xyyzz, g_y_0_z_xyzz, g_y_0_z_xyzzz, g_y_0_z_xzzz, g_y_0_z_xzzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xz_xxxx[k] = -g_y_0_z_xxxx[k] * cd_x[k] + g_y_0_z_xxxxx[k];

            g_y_0_xz_xxxy[k] = -g_y_0_z_xxxy[k] * cd_x[k] + g_y_0_z_xxxxy[k];

            g_y_0_xz_xxxz[k] = -g_y_0_z_xxxz[k] * cd_x[k] + g_y_0_z_xxxxz[k];

            g_y_0_xz_xxyy[k] = -g_y_0_z_xxyy[k] * cd_x[k] + g_y_0_z_xxxyy[k];

            g_y_0_xz_xxyz[k] = -g_y_0_z_xxyz[k] * cd_x[k] + g_y_0_z_xxxyz[k];

            g_y_0_xz_xxzz[k] = -g_y_0_z_xxzz[k] * cd_x[k] + g_y_0_z_xxxzz[k];

            g_y_0_xz_xyyy[k] = -g_y_0_z_xyyy[k] * cd_x[k] + g_y_0_z_xxyyy[k];

            g_y_0_xz_xyyz[k] = -g_y_0_z_xyyz[k] * cd_x[k] + g_y_0_z_xxyyz[k];

            g_y_0_xz_xyzz[k] = -g_y_0_z_xyzz[k] * cd_x[k] + g_y_0_z_xxyzz[k];

            g_y_0_xz_xzzz[k] = -g_y_0_z_xzzz[k] * cd_x[k] + g_y_0_z_xxzzz[k];

            g_y_0_xz_yyyy[k] = -g_y_0_z_yyyy[k] * cd_x[k] + g_y_0_z_xyyyy[k];

            g_y_0_xz_yyyz[k] = -g_y_0_z_yyyz[k] * cd_x[k] + g_y_0_z_xyyyz[k];

            g_y_0_xz_yyzz[k] = -g_y_0_z_yyzz[k] * cd_x[k] + g_y_0_z_xyyzz[k];

            g_y_0_xz_yzzz[k] = -g_y_0_z_yzzz[k] * cd_x[k] + g_y_0_z_xyzzz[k];

            g_y_0_xz_zzzz[k] = -g_y_0_z_zzzz[k] * cd_x[k] + g_y_0_z_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 45);

        auto g_y_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 46);

        auto g_y_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 47);

        auto g_y_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 48);

        auto g_y_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 49);

        auto g_y_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 50);

        auto g_y_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 51);

        auto g_y_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 52);

        auto g_y_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 53);

        auto g_y_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 54);

        auto g_y_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 55);

        auto g_y_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 56);

        auto g_y_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 57);

        auto g_y_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 58);

        auto g_y_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 59);

        #pragma omp simd aligned(cd_y, g_y_0_y_xxxx, g_y_0_y_xxxxy, g_y_0_y_xxxy, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyz, g_y_0_y_yyyzz, g_y_0_y_yyzz, g_y_0_y_yyzzz, g_y_0_y_yzzz, g_y_0_y_yzzzz, g_y_0_y_zzzz, g_y_0_yy_xxxx, g_y_0_yy_xxxy, g_y_0_yy_xxxz, g_y_0_yy_xxyy, g_y_0_yy_xxyz, g_y_0_yy_xxzz, g_y_0_yy_xyyy, g_y_0_yy_xyyz, g_y_0_yy_xyzz, g_y_0_yy_xzzz, g_y_0_yy_yyyy, g_y_0_yy_yyyz, g_y_0_yy_yyzz, g_y_0_yy_yzzz, g_y_0_yy_zzzz, g_y_xxxx, g_y_xxxy, g_y_xxxz, g_y_xxyy, g_y_xxyz, g_y_xxzz, g_y_xyyy, g_y_xyyz, g_y_xyzz, g_y_xzzz, g_y_yyyy, g_y_yyyz, g_y_yyzz, g_y_yzzz, g_y_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yy_xxxx[k] = -g_y_xxxx[k] - g_y_0_y_xxxx[k] * cd_y[k] + g_y_0_y_xxxxy[k];

            g_y_0_yy_xxxy[k] = -g_y_xxxy[k] - g_y_0_y_xxxy[k] * cd_y[k] + g_y_0_y_xxxyy[k];

            g_y_0_yy_xxxz[k] = -g_y_xxxz[k] - g_y_0_y_xxxz[k] * cd_y[k] + g_y_0_y_xxxyz[k];

            g_y_0_yy_xxyy[k] = -g_y_xxyy[k] - g_y_0_y_xxyy[k] * cd_y[k] + g_y_0_y_xxyyy[k];

            g_y_0_yy_xxyz[k] = -g_y_xxyz[k] - g_y_0_y_xxyz[k] * cd_y[k] + g_y_0_y_xxyyz[k];

            g_y_0_yy_xxzz[k] = -g_y_xxzz[k] - g_y_0_y_xxzz[k] * cd_y[k] + g_y_0_y_xxyzz[k];

            g_y_0_yy_xyyy[k] = -g_y_xyyy[k] - g_y_0_y_xyyy[k] * cd_y[k] + g_y_0_y_xyyyy[k];

            g_y_0_yy_xyyz[k] = -g_y_xyyz[k] - g_y_0_y_xyyz[k] * cd_y[k] + g_y_0_y_xyyyz[k];

            g_y_0_yy_xyzz[k] = -g_y_xyzz[k] - g_y_0_y_xyzz[k] * cd_y[k] + g_y_0_y_xyyzz[k];

            g_y_0_yy_xzzz[k] = -g_y_xzzz[k] - g_y_0_y_xzzz[k] * cd_y[k] + g_y_0_y_xyzzz[k];

            g_y_0_yy_yyyy[k] = -g_y_yyyy[k] - g_y_0_y_yyyy[k] * cd_y[k] + g_y_0_y_yyyyy[k];

            g_y_0_yy_yyyz[k] = -g_y_yyyz[k] - g_y_0_y_yyyz[k] * cd_y[k] + g_y_0_y_yyyyz[k];

            g_y_0_yy_yyzz[k] = -g_y_yyzz[k] - g_y_0_y_yyzz[k] * cd_y[k] + g_y_0_y_yyyzz[k];

            g_y_0_yy_yzzz[k] = -g_y_yzzz[k] - g_y_0_y_yzzz[k] * cd_y[k] + g_y_0_y_yyzzz[k];

            g_y_0_yy_zzzz[k] = -g_y_zzzz[k] - g_y_0_y_zzzz[k] * cd_y[k] + g_y_0_y_yzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_y_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 60);

        auto g_y_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 61);

        auto g_y_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 62);

        auto g_y_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 63);

        auto g_y_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 64);

        auto g_y_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 65);

        auto g_y_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 66);

        auto g_y_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 67);

        auto g_y_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 68);

        auto g_y_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 69);

        auto g_y_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 70);

        auto g_y_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 71);

        auto g_y_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 72);

        auto g_y_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 73);

        auto g_y_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 74);

        #pragma omp simd aligned(cd_z, g_y_0_y_xxxx, g_y_0_y_xxxxz, g_y_0_y_xxxy, g_y_0_y_xxxyz, g_y_0_y_xxxz, g_y_0_y_xxxzz, g_y_0_y_xxyy, g_y_0_y_xxyyz, g_y_0_y_xxyz, g_y_0_y_xxyzz, g_y_0_y_xxzz, g_y_0_y_xxzzz, g_y_0_y_xyyy, g_y_0_y_xyyyz, g_y_0_y_xyyz, g_y_0_y_xyyzz, g_y_0_y_xyzz, g_y_0_y_xyzzz, g_y_0_y_xzzz, g_y_0_y_xzzzz, g_y_0_y_yyyy, g_y_0_y_yyyyz, g_y_0_y_yyyz, g_y_0_y_yyyzz, g_y_0_y_yyzz, g_y_0_y_yyzzz, g_y_0_y_yzzz, g_y_0_y_yzzzz, g_y_0_y_zzzz, g_y_0_y_zzzzz, g_y_0_yz_xxxx, g_y_0_yz_xxxy, g_y_0_yz_xxxz, g_y_0_yz_xxyy, g_y_0_yz_xxyz, g_y_0_yz_xxzz, g_y_0_yz_xyyy, g_y_0_yz_xyyz, g_y_0_yz_xyzz, g_y_0_yz_xzzz, g_y_0_yz_yyyy, g_y_0_yz_yyyz, g_y_0_yz_yyzz, g_y_0_yz_yzzz, g_y_0_yz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yz_xxxx[k] = -g_y_0_y_xxxx[k] * cd_z[k] + g_y_0_y_xxxxz[k];

            g_y_0_yz_xxxy[k] = -g_y_0_y_xxxy[k] * cd_z[k] + g_y_0_y_xxxyz[k];

            g_y_0_yz_xxxz[k] = -g_y_0_y_xxxz[k] * cd_z[k] + g_y_0_y_xxxzz[k];

            g_y_0_yz_xxyy[k] = -g_y_0_y_xxyy[k] * cd_z[k] + g_y_0_y_xxyyz[k];

            g_y_0_yz_xxyz[k] = -g_y_0_y_xxyz[k] * cd_z[k] + g_y_0_y_xxyzz[k];

            g_y_0_yz_xxzz[k] = -g_y_0_y_xxzz[k] * cd_z[k] + g_y_0_y_xxzzz[k];

            g_y_0_yz_xyyy[k] = -g_y_0_y_xyyy[k] * cd_z[k] + g_y_0_y_xyyyz[k];

            g_y_0_yz_xyyz[k] = -g_y_0_y_xyyz[k] * cd_z[k] + g_y_0_y_xyyzz[k];

            g_y_0_yz_xyzz[k] = -g_y_0_y_xyzz[k] * cd_z[k] + g_y_0_y_xyzzz[k];

            g_y_0_yz_xzzz[k] = -g_y_0_y_xzzz[k] * cd_z[k] + g_y_0_y_xzzzz[k];

            g_y_0_yz_yyyy[k] = -g_y_0_y_yyyy[k] * cd_z[k] + g_y_0_y_yyyyz[k];

            g_y_0_yz_yyyz[k] = -g_y_0_y_yyyz[k] * cd_z[k] + g_y_0_y_yyyzz[k];

            g_y_0_yz_yyzz[k] = -g_y_0_y_yyzz[k] * cd_z[k] + g_y_0_y_yyzzz[k];

            g_y_0_yz_yzzz[k] = -g_y_0_y_yzzz[k] * cd_z[k] + g_y_0_y_yzzzz[k];

            g_y_0_yz_zzzz[k] = -g_y_0_y_zzzz[k] * cd_z[k] + g_y_0_y_zzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_y_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 90 * acomps  + 75);

        auto g_y_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 76);

        auto g_y_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 77);

        auto g_y_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 78);

        auto g_y_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 79);

        auto g_y_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 80);

        auto g_y_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 81);

        auto g_y_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 82);

        auto g_y_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 83);

        auto g_y_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 84);

        auto g_y_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 90 * acomps  + 85);

        auto g_y_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 86);

        auto g_y_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 87);

        auto g_y_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 88);

        auto g_y_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 90 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_y_0_z_xxxx, g_y_0_z_xxxxz, g_y_0_z_xxxy, g_y_0_z_xxxyz, g_y_0_z_xxxz, g_y_0_z_xxxzz, g_y_0_z_xxyy, g_y_0_z_xxyyz, g_y_0_z_xxyz, g_y_0_z_xxyzz, g_y_0_z_xxzz, g_y_0_z_xxzzz, g_y_0_z_xyyy, g_y_0_z_xyyyz, g_y_0_z_xyyz, g_y_0_z_xyyzz, g_y_0_z_xyzz, g_y_0_z_xyzzz, g_y_0_z_xzzz, g_y_0_z_xzzzz, g_y_0_z_yyyy, g_y_0_z_yyyyz, g_y_0_z_yyyz, g_y_0_z_yyyzz, g_y_0_z_yyzz, g_y_0_z_yyzzz, g_y_0_z_yzzz, g_y_0_z_yzzzz, g_y_0_z_zzzz, g_y_0_z_zzzzz, g_y_0_zz_xxxx, g_y_0_zz_xxxy, g_y_0_zz_xxxz, g_y_0_zz_xxyy, g_y_0_zz_xxyz, g_y_0_zz_xxzz, g_y_0_zz_xyyy, g_y_0_zz_xyyz, g_y_0_zz_xyzz, g_y_0_zz_xzzz, g_y_0_zz_yyyy, g_y_0_zz_yyyz, g_y_0_zz_yyzz, g_y_0_zz_yzzz, g_y_0_zz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zz_xxxx[k] = -g_y_0_z_xxxx[k] * cd_z[k] + g_y_0_z_xxxxz[k];

            g_y_0_zz_xxxy[k] = -g_y_0_z_xxxy[k] * cd_z[k] + g_y_0_z_xxxyz[k];

            g_y_0_zz_xxxz[k] = -g_y_0_z_xxxz[k] * cd_z[k] + g_y_0_z_xxxzz[k];

            g_y_0_zz_xxyy[k] = -g_y_0_z_xxyy[k] * cd_z[k] + g_y_0_z_xxyyz[k];

            g_y_0_zz_xxyz[k] = -g_y_0_z_xxyz[k] * cd_z[k] + g_y_0_z_xxyzz[k];

            g_y_0_zz_xxzz[k] = -g_y_0_z_xxzz[k] * cd_z[k] + g_y_0_z_xxzzz[k];

            g_y_0_zz_xyyy[k] = -g_y_0_z_xyyy[k] * cd_z[k] + g_y_0_z_xyyyz[k];

            g_y_0_zz_xyyz[k] = -g_y_0_z_xyyz[k] * cd_z[k] + g_y_0_z_xyyzz[k];

            g_y_0_zz_xyzz[k] = -g_y_0_z_xyzz[k] * cd_z[k] + g_y_0_z_xyzzz[k];

            g_y_0_zz_xzzz[k] = -g_y_0_z_xzzz[k] * cd_z[k] + g_y_0_z_xzzzz[k];

            g_y_0_zz_yyyy[k] = -g_y_0_z_yyyy[k] * cd_z[k] + g_y_0_z_yyyyz[k];

            g_y_0_zz_yyyz[k] = -g_y_0_z_yyyz[k] * cd_z[k] + g_y_0_z_yyyzz[k];

            g_y_0_zz_yyzz[k] = -g_y_0_z_yyzz[k] * cd_z[k] + g_y_0_z_yyzzz[k];

            g_y_0_zz_yzzz[k] = -g_y_0_z_yzzz[k] * cd_z[k] + g_y_0_z_yzzzz[k];

            g_y_0_zz_zzzz[k] = -g_y_0_z_zzzz[k] * cd_z[k] + g_y_0_z_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_xx_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 0);

        auto g_z_0_xx_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 1);

        auto g_z_0_xx_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 2);

        auto g_z_0_xx_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 3);

        auto g_z_0_xx_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 4);

        auto g_z_0_xx_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 5);

        auto g_z_0_xx_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 6);

        auto g_z_0_xx_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 7);

        auto g_z_0_xx_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 8);

        auto g_z_0_xx_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 9);

        auto g_z_0_xx_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 10);

        auto g_z_0_xx_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 11);

        auto g_z_0_xx_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 12);

        auto g_z_0_xx_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 13);

        auto g_z_0_xx_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_z_0_x_xxxx, g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxy, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxz, g_z_0_x_xxxzz, g_z_0_x_xxyy, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyz, g_z_0_x_xxyzz, g_z_0_x_xxzz, g_z_0_x_xxzzz, g_z_0_x_xyyy, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyz, g_z_0_x_xyyzz, g_z_0_x_xyzz, g_z_0_x_xyzzz, g_z_0_x_xzzz, g_z_0_x_xzzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz, g_z_0_xx_xxxx, g_z_0_xx_xxxy, g_z_0_xx_xxxz, g_z_0_xx_xxyy, g_z_0_xx_xxyz, g_z_0_xx_xxzz, g_z_0_xx_xyyy, g_z_0_xx_xyyz, g_z_0_xx_xyzz, g_z_0_xx_xzzz, g_z_0_xx_yyyy, g_z_0_xx_yyyz, g_z_0_xx_yyzz, g_z_0_xx_yzzz, g_z_0_xx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xx_xxxx[k] = -g_z_0_x_xxxx[k] * cd_x[k] + g_z_0_x_xxxxx[k];

            g_z_0_xx_xxxy[k] = -g_z_0_x_xxxy[k] * cd_x[k] + g_z_0_x_xxxxy[k];

            g_z_0_xx_xxxz[k] = -g_z_0_x_xxxz[k] * cd_x[k] + g_z_0_x_xxxxz[k];

            g_z_0_xx_xxyy[k] = -g_z_0_x_xxyy[k] * cd_x[k] + g_z_0_x_xxxyy[k];

            g_z_0_xx_xxyz[k] = -g_z_0_x_xxyz[k] * cd_x[k] + g_z_0_x_xxxyz[k];

            g_z_0_xx_xxzz[k] = -g_z_0_x_xxzz[k] * cd_x[k] + g_z_0_x_xxxzz[k];

            g_z_0_xx_xyyy[k] = -g_z_0_x_xyyy[k] * cd_x[k] + g_z_0_x_xxyyy[k];

            g_z_0_xx_xyyz[k] = -g_z_0_x_xyyz[k] * cd_x[k] + g_z_0_x_xxyyz[k];

            g_z_0_xx_xyzz[k] = -g_z_0_x_xyzz[k] * cd_x[k] + g_z_0_x_xxyzz[k];

            g_z_0_xx_xzzz[k] = -g_z_0_x_xzzz[k] * cd_x[k] + g_z_0_x_xxzzz[k];

            g_z_0_xx_yyyy[k] = -g_z_0_x_yyyy[k] * cd_x[k] + g_z_0_x_xyyyy[k];

            g_z_0_xx_yyyz[k] = -g_z_0_x_yyyz[k] * cd_x[k] + g_z_0_x_xyyyz[k];

            g_z_0_xx_yyzz[k] = -g_z_0_x_yyzz[k] * cd_x[k] + g_z_0_x_xyyzz[k];

            g_z_0_xx_yzzz[k] = -g_z_0_x_yzzz[k] * cd_x[k] + g_z_0_x_xyzzz[k];

            g_z_0_xx_zzzz[k] = -g_z_0_x_zzzz[k] * cd_x[k] + g_z_0_x_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xy_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 15);

        auto g_z_0_xy_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 16);

        auto g_z_0_xy_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 17);

        auto g_z_0_xy_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 18);

        auto g_z_0_xy_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 19);

        auto g_z_0_xy_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 20);

        auto g_z_0_xy_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 21);

        auto g_z_0_xy_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 22);

        auto g_z_0_xy_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 23);

        auto g_z_0_xy_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 24);

        auto g_z_0_xy_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 25);

        auto g_z_0_xy_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 26);

        auto g_z_0_xy_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 27);

        auto g_z_0_xy_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 28);

        auto g_z_0_xy_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xy_xxxx, g_z_0_xy_xxxy, g_z_0_xy_xxxz, g_z_0_xy_xxyy, g_z_0_xy_xxyz, g_z_0_xy_xxzz, g_z_0_xy_xyyy, g_z_0_xy_xyyz, g_z_0_xy_xyzz, g_z_0_xy_xzzz, g_z_0_xy_yyyy, g_z_0_xy_yyyz, g_z_0_xy_yyzz, g_z_0_xy_yzzz, g_z_0_xy_zzzz, g_z_0_y_xxxx, g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxy, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxz, g_z_0_y_xxxzz, g_z_0_y_xxyy, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyz, g_z_0_y_xxyzz, g_z_0_y_xxzz, g_z_0_y_xxzzz, g_z_0_y_xyyy, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyz, g_z_0_y_xyyzz, g_z_0_y_xyzz, g_z_0_y_xyzzz, g_z_0_y_xzzz, g_z_0_y_xzzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xy_xxxx[k] = -g_z_0_y_xxxx[k] * cd_x[k] + g_z_0_y_xxxxx[k];

            g_z_0_xy_xxxy[k] = -g_z_0_y_xxxy[k] * cd_x[k] + g_z_0_y_xxxxy[k];

            g_z_0_xy_xxxz[k] = -g_z_0_y_xxxz[k] * cd_x[k] + g_z_0_y_xxxxz[k];

            g_z_0_xy_xxyy[k] = -g_z_0_y_xxyy[k] * cd_x[k] + g_z_0_y_xxxyy[k];

            g_z_0_xy_xxyz[k] = -g_z_0_y_xxyz[k] * cd_x[k] + g_z_0_y_xxxyz[k];

            g_z_0_xy_xxzz[k] = -g_z_0_y_xxzz[k] * cd_x[k] + g_z_0_y_xxxzz[k];

            g_z_0_xy_xyyy[k] = -g_z_0_y_xyyy[k] * cd_x[k] + g_z_0_y_xxyyy[k];

            g_z_0_xy_xyyz[k] = -g_z_0_y_xyyz[k] * cd_x[k] + g_z_0_y_xxyyz[k];

            g_z_0_xy_xyzz[k] = -g_z_0_y_xyzz[k] * cd_x[k] + g_z_0_y_xxyzz[k];

            g_z_0_xy_xzzz[k] = -g_z_0_y_xzzz[k] * cd_x[k] + g_z_0_y_xxzzz[k];

            g_z_0_xy_yyyy[k] = -g_z_0_y_yyyy[k] * cd_x[k] + g_z_0_y_xyyyy[k];

            g_z_0_xy_yyyz[k] = -g_z_0_y_yyyz[k] * cd_x[k] + g_z_0_y_xyyyz[k];

            g_z_0_xy_yyzz[k] = -g_z_0_y_yyzz[k] * cd_x[k] + g_z_0_y_xyyzz[k];

            g_z_0_xy_yzzz[k] = -g_z_0_y_yzzz[k] * cd_x[k] + g_z_0_y_xyzzz[k];

            g_z_0_xy_zzzz[k] = -g_z_0_y_zzzz[k] * cd_x[k] + g_z_0_y_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_z_0_xz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 30);

        auto g_z_0_xz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 31);

        auto g_z_0_xz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 32);

        auto g_z_0_xz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 33);

        auto g_z_0_xz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 34);

        auto g_z_0_xz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 35);

        auto g_z_0_xz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 36);

        auto g_z_0_xz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 37);

        auto g_z_0_xz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 38);

        auto g_z_0_xz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 39);

        auto g_z_0_xz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 40);

        auto g_z_0_xz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 41);

        auto g_z_0_xz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 42);

        auto g_z_0_xz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 43);

        auto g_z_0_xz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_z_0_xz_xxxx, g_z_0_xz_xxxy, g_z_0_xz_xxxz, g_z_0_xz_xxyy, g_z_0_xz_xxyz, g_z_0_xz_xxzz, g_z_0_xz_xyyy, g_z_0_xz_xyyz, g_z_0_xz_xyzz, g_z_0_xz_xzzz, g_z_0_xz_yyyy, g_z_0_xz_yyyz, g_z_0_xz_yyzz, g_z_0_xz_yzzz, g_z_0_xz_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxy, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxxzz, g_z_0_z_xxyy, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xxzzz, g_z_0_z_xyyy, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_xzzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xz_xxxx[k] = -g_z_0_z_xxxx[k] * cd_x[k] + g_z_0_z_xxxxx[k];

            g_z_0_xz_xxxy[k] = -g_z_0_z_xxxy[k] * cd_x[k] + g_z_0_z_xxxxy[k];

            g_z_0_xz_xxxz[k] = -g_z_0_z_xxxz[k] * cd_x[k] + g_z_0_z_xxxxz[k];

            g_z_0_xz_xxyy[k] = -g_z_0_z_xxyy[k] * cd_x[k] + g_z_0_z_xxxyy[k];

            g_z_0_xz_xxyz[k] = -g_z_0_z_xxyz[k] * cd_x[k] + g_z_0_z_xxxyz[k];

            g_z_0_xz_xxzz[k] = -g_z_0_z_xxzz[k] * cd_x[k] + g_z_0_z_xxxzz[k];

            g_z_0_xz_xyyy[k] = -g_z_0_z_xyyy[k] * cd_x[k] + g_z_0_z_xxyyy[k];

            g_z_0_xz_xyyz[k] = -g_z_0_z_xyyz[k] * cd_x[k] + g_z_0_z_xxyyz[k];

            g_z_0_xz_xyzz[k] = -g_z_0_z_xyzz[k] * cd_x[k] + g_z_0_z_xxyzz[k];

            g_z_0_xz_xzzz[k] = -g_z_0_z_xzzz[k] * cd_x[k] + g_z_0_z_xxzzz[k];

            g_z_0_xz_yyyy[k] = -g_z_0_z_yyyy[k] * cd_x[k] + g_z_0_z_xyyyy[k];

            g_z_0_xz_yyyz[k] = -g_z_0_z_yyyz[k] * cd_x[k] + g_z_0_z_xyyyz[k];

            g_z_0_xz_yyzz[k] = -g_z_0_z_yyzz[k] * cd_x[k] + g_z_0_z_xyyzz[k];

            g_z_0_xz_yzzz[k] = -g_z_0_z_yzzz[k] * cd_x[k] + g_z_0_z_xyzzz[k];

            g_z_0_xz_zzzz[k] = -g_z_0_z_zzzz[k] * cd_x[k] + g_z_0_z_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_yy_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 45);

        auto g_z_0_yy_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 46);

        auto g_z_0_yy_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 47);

        auto g_z_0_yy_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 48);

        auto g_z_0_yy_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 49);

        auto g_z_0_yy_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 50);

        auto g_z_0_yy_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 51);

        auto g_z_0_yy_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 52);

        auto g_z_0_yy_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 53);

        auto g_z_0_yy_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 54);

        auto g_z_0_yy_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 55);

        auto g_z_0_yy_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 56);

        auto g_z_0_yy_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 57);

        auto g_z_0_yy_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 58);

        auto g_z_0_yy_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 59);

        #pragma omp simd aligned(cd_y, g_z_0_y_xxxx, g_z_0_y_xxxxy, g_z_0_y_xxxy, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyz, g_z_0_y_xxyzz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyz, g_z_0_y_xyyzz, g_z_0_y_xyzz, g_z_0_y_xyzzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyz, g_z_0_y_yyyzz, g_z_0_y_yyzz, g_z_0_y_yyzzz, g_z_0_y_yzzz, g_z_0_y_yzzzz, g_z_0_y_zzzz, g_z_0_yy_xxxx, g_z_0_yy_xxxy, g_z_0_yy_xxxz, g_z_0_yy_xxyy, g_z_0_yy_xxyz, g_z_0_yy_xxzz, g_z_0_yy_xyyy, g_z_0_yy_xyyz, g_z_0_yy_xyzz, g_z_0_yy_xzzz, g_z_0_yy_yyyy, g_z_0_yy_yyyz, g_z_0_yy_yyzz, g_z_0_yy_yzzz, g_z_0_yy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yy_xxxx[k] = -g_z_0_y_xxxx[k] * cd_y[k] + g_z_0_y_xxxxy[k];

            g_z_0_yy_xxxy[k] = -g_z_0_y_xxxy[k] * cd_y[k] + g_z_0_y_xxxyy[k];

            g_z_0_yy_xxxz[k] = -g_z_0_y_xxxz[k] * cd_y[k] + g_z_0_y_xxxyz[k];

            g_z_0_yy_xxyy[k] = -g_z_0_y_xxyy[k] * cd_y[k] + g_z_0_y_xxyyy[k];

            g_z_0_yy_xxyz[k] = -g_z_0_y_xxyz[k] * cd_y[k] + g_z_0_y_xxyyz[k];

            g_z_0_yy_xxzz[k] = -g_z_0_y_xxzz[k] * cd_y[k] + g_z_0_y_xxyzz[k];

            g_z_0_yy_xyyy[k] = -g_z_0_y_xyyy[k] * cd_y[k] + g_z_0_y_xyyyy[k];

            g_z_0_yy_xyyz[k] = -g_z_0_y_xyyz[k] * cd_y[k] + g_z_0_y_xyyyz[k];

            g_z_0_yy_xyzz[k] = -g_z_0_y_xyzz[k] * cd_y[k] + g_z_0_y_xyyzz[k];

            g_z_0_yy_xzzz[k] = -g_z_0_y_xzzz[k] * cd_y[k] + g_z_0_y_xyzzz[k];

            g_z_0_yy_yyyy[k] = -g_z_0_y_yyyy[k] * cd_y[k] + g_z_0_y_yyyyy[k];

            g_z_0_yy_yyyz[k] = -g_z_0_y_yyyz[k] * cd_y[k] + g_z_0_y_yyyyz[k];

            g_z_0_yy_yyzz[k] = -g_z_0_y_yyzz[k] * cd_y[k] + g_z_0_y_yyyzz[k];

            g_z_0_yy_yzzz[k] = -g_z_0_y_yzzz[k] * cd_y[k] + g_z_0_y_yyzzz[k];

            g_z_0_yy_zzzz[k] = -g_z_0_y_zzzz[k] * cd_y[k] + g_z_0_y_yzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_z_0_yz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 60);

        auto g_z_0_yz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 61);

        auto g_z_0_yz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 62);

        auto g_z_0_yz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 63);

        auto g_z_0_yz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 64);

        auto g_z_0_yz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 65);

        auto g_z_0_yz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 66);

        auto g_z_0_yz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 67);

        auto g_z_0_yz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 68);

        auto g_z_0_yz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 69);

        auto g_z_0_yz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 70);

        auto g_z_0_yz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 71);

        auto g_z_0_yz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 72);

        auto g_z_0_yz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 73);

        auto g_z_0_yz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 74);

        #pragma omp simd aligned(cd_y, g_z_0_yz_xxxx, g_z_0_yz_xxxy, g_z_0_yz_xxxz, g_z_0_yz_xxyy, g_z_0_yz_xxyz, g_z_0_yz_xxzz, g_z_0_yz_xyyy, g_z_0_yz_xyyz, g_z_0_yz_xyzz, g_z_0_yz_xzzz, g_z_0_yz_yyyy, g_z_0_yz_yyyz, g_z_0_yz_yyzz, g_z_0_yz_yzzz, g_z_0_yz_zzzz, g_z_0_z_xxxx, g_z_0_z_xxxxy, g_z_0_z_xxxy, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyz, g_z_0_z_yyyzz, g_z_0_z_yyzz, g_z_0_z_yyzzz, g_z_0_z_yzzz, g_z_0_z_yzzzz, g_z_0_z_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yz_xxxx[k] = -g_z_0_z_xxxx[k] * cd_y[k] + g_z_0_z_xxxxy[k];

            g_z_0_yz_xxxy[k] = -g_z_0_z_xxxy[k] * cd_y[k] + g_z_0_z_xxxyy[k];

            g_z_0_yz_xxxz[k] = -g_z_0_z_xxxz[k] * cd_y[k] + g_z_0_z_xxxyz[k];

            g_z_0_yz_xxyy[k] = -g_z_0_z_xxyy[k] * cd_y[k] + g_z_0_z_xxyyy[k];

            g_z_0_yz_xxyz[k] = -g_z_0_z_xxyz[k] * cd_y[k] + g_z_0_z_xxyyz[k];

            g_z_0_yz_xxzz[k] = -g_z_0_z_xxzz[k] * cd_y[k] + g_z_0_z_xxyzz[k];

            g_z_0_yz_xyyy[k] = -g_z_0_z_xyyy[k] * cd_y[k] + g_z_0_z_xyyyy[k];

            g_z_0_yz_xyyz[k] = -g_z_0_z_xyyz[k] * cd_y[k] + g_z_0_z_xyyyz[k];

            g_z_0_yz_xyzz[k] = -g_z_0_z_xyzz[k] * cd_y[k] + g_z_0_z_xyyzz[k];

            g_z_0_yz_xzzz[k] = -g_z_0_z_xzzz[k] * cd_y[k] + g_z_0_z_xyzzz[k];

            g_z_0_yz_yyyy[k] = -g_z_0_z_yyyy[k] * cd_y[k] + g_z_0_z_yyyyy[k];

            g_z_0_yz_yyyz[k] = -g_z_0_z_yyyz[k] * cd_y[k] + g_z_0_z_yyyyz[k];

            g_z_0_yz_yyzz[k] = -g_z_0_z_yyzz[k] * cd_y[k] + g_z_0_z_yyyzz[k];

            g_z_0_yz_yzzz[k] = -g_z_0_z_yzzz[k] * cd_y[k] + g_z_0_z_yyzzz[k];

            g_z_0_yz_zzzz[k] = -g_z_0_z_zzzz[k] * cd_y[k] + g_z_0_z_yzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_z_0_zz_xxxx = cbuffer.data(dg_geom_10_off + 180 * acomps  + 75);

        auto g_z_0_zz_xxxy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 76);

        auto g_z_0_zz_xxxz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 77);

        auto g_z_0_zz_xxyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 78);

        auto g_z_0_zz_xxyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 79);

        auto g_z_0_zz_xxzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 80);

        auto g_z_0_zz_xyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 81);

        auto g_z_0_zz_xyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 82);

        auto g_z_0_zz_xyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 83);

        auto g_z_0_zz_xzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 84);

        auto g_z_0_zz_yyyy = cbuffer.data(dg_geom_10_off + 180 * acomps  + 85);

        auto g_z_0_zz_yyyz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 86);

        auto g_z_0_zz_yyzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 87);

        auto g_z_0_zz_yzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 88);

        auto g_z_0_zz_zzzz = cbuffer.data(dg_geom_10_off + 180 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_z_0_z_xxxx, g_z_0_z_xxxxz, g_z_0_z_xxxy, g_z_0_z_xxxyz, g_z_0_z_xxxz, g_z_0_z_xxxzz, g_z_0_z_xxyy, g_z_0_z_xxyyz, g_z_0_z_xxyz, g_z_0_z_xxyzz, g_z_0_z_xxzz, g_z_0_z_xxzzz, g_z_0_z_xyyy, g_z_0_z_xyyyz, g_z_0_z_xyyz, g_z_0_z_xyyzz, g_z_0_z_xyzz, g_z_0_z_xyzzz, g_z_0_z_xzzz, g_z_0_z_xzzzz, g_z_0_z_yyyy, g_z_0_z_yyyyz, g_z_0_z_yyyz, g_z_0_z_yyyzz, g_z_0_z_yyzz, g_z_0_z_yyzzz, g_z_0_z_yzzz, g_z_0_z_yzzzz, g_z_0_z_zzzz, g_z_0_z_zzzzz, g_z_0_zz_xxxx, g_z_0_zz_xxxy, g_z_0_zz_xxxz, g_z_0_zz_xxyy, g_z_0_zz_xxyz, g_z_0_zz_xxzz, g_z_0_zz_xyyy, g_z_0_zz_xyyz, g_z_0_zz_xyzz, g_z_0_zz_xzzz, g_z_0_zz_yyyy, g_z_0_zz_yyyz, g_z_0_zz_yyzz, g_z_0_zz_yzzz, g_z_0_zz_zzzz, g_z_xxxx, g_z_xxxy, g_z_xxxz, g_z_xxyy, g_z_xxyz, g_z_xxzz, g_z_xyyy, g_z_xyyz, g_z_xyzz, g_z_xzzz, g_z_yyyy, g_z_yyyz, g_z_yyzz, g_z_yzzz, g_z_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zz_xxxx[k] = -g_z_xxxx[k] - g_z_0_z_xxxx[k] * cd_z[k] + g_z_0_z_xxxxz[k];

            g_z_0_zz_xxxy[k] = -g_z_xxxy[k] - g_z_0_z_xxxy[k] * cd_z[k] + g_z_0_z_xxxyz[k];

            g_z_0_zz_xxxz[k] = -g_z_xxxz[k] - g_z_0_z_xxxz[k] * cd_z[k] + g_z_0_z_xxxzz[k];

            g_z_0_zz_xxyy[k] = -g_z_xxyy[k] - g_z_0_z_xxyy[k] * cd_z[k] + g_z_0_z_xxyyz[k];

            g_z_0_zz_xxyz[k] = -g_z_xxyz[k] - g_z_0_z_xxyz[k] * cd_z[k] + g_z_0_z_xxyzz[k];

            g_z_0_zz_xxzz[k] = -g_z_xxzz[k] - g_z_0_z_xxzz[k] * cd_z[k] + g_z_0_z_xxzzz[k];

            g_z_0_zz_xyyy[k] = -g_z_xyyy[k] - g_z_0_z_xyyy[k] * cd_z[k] + g_z_0_z_xyyyz[k];

            g_z_0_zz_xyyz[k] = -g_z_xyyz[k] - g_z_0_z_xyyz[k] * cd_z[k] + g_z_0_z_xyyzz[k];

            g_z_0_zz_xyzz[k] = -g_z_xyzz[k] - g_z_0_z_xyzz[k] * cd_z[k] + g_z_0_z_xyzzz[k];

            g_z_0_zz_xzzz[k] = -g_z_xzzz[k] - g_z_0_z_xzzz[k] * cd_z[k] + g_z_0_z_xzzzz[k];

            g_z_0_zz_yyyy[k] = -g_z_yyyy[k] - g_z_0_z_yyyy[k] * cd_z[k] + g_z_0_z_yyyyz[k];

            g_z_0_zz_yyyz[k] = -g_z_yyyz[k] - g_z_0_z_yyyz[k] * cd_z[k] + g_z_0_z_yyyzz[k];

            g_z_0_zz_yyzz[k] = -g_z_yyzz[k] - g_z_0_z_yyzz[k] * cd_z[k] + g_z_0_z_yyzzz[k];

            g_z_0_zz_yzzz[k] = -g_z_yzzz[k] - g_z_0_z_yzzz[k] * cd_z[k] + g_z_0_z_yzzzz[k];

            g_z_0_zz_zzzz[k] = -g_z_zzzz[k] - g_z_0_z_zzzz[k] * cd_z[k] + g_z_0_z_zzzzz[k];
        }
    }
}

} // t3ceri namespace


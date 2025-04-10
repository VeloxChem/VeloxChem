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

#include "ElectronRepulsionContrRecPGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_hrr_electron_repulsion_pgxx(CSimdArray<double>&   cbuffer,
                                     const size_t          idx_pgxx,
                                     const size_t          idx_sgxx,
                                     const size_t          idx_shxx,
                                     const TPoint<double>& r_ab,
                                     const int             c_angmom,
                                     const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SGSS

            const auto sg_off = idx_sgxx + i * dcomps + j;

            auto g_0_xxxx = cbuffer.data(sg_off + 0 * ccomps * dcomps);

            auto g_0_xxxy = cbuffer.data(sg_off + 1 * ccomps * dcomps);

            auto g_0_xxxz = cbuffer.data(sg_off + 2 * ccomps * dcomps);

            auto g_0_xxyy = cbuffer.data(sg_off + 3 * ccomps * dcomps);

            auto g_0_xxyz = cbuffer.data(sg_off + 4 * ccomps * dcomps);

            auto g_0_xxzz = cbuffer.data(sg_off + 5 * ccomps * dcomps);

            auto g_0_xyyy = cbuffer.data(sg_off + 6 * ccomps * dcomps);

            auto g_0_xyyz = cbuffer.data(sg_off + 7 * ccomps * dcomps);

            auto g_0_xyzz = cbuffer.data(sg_off + 8 * ccomps * dcomps);

            auto g_0_xzzz = cbuffer.data(sg_off + 9 * ccomps * dcomps);

            auto g_0_yyyy = cbuffer.data(sg_off + 10 * ccomps * dcomps);

            auto g_0_yyyz = cbuffer.data(sg_off + 11 * ccomps * dcomps);

            auto g_0_yyzz = cbuffer.data(sg_off + 12 * ccomps * dcomps);

            auto g_0_yzzz = cbuffer.data(sg_off + 13 * ccomps * dcomps);

            auto g_0_zzzz = cbuffer.data(sg_off + 14 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SHSS

            const auto sh_off = idx_shxx + i * dcomps + j;

            auto g_0_xxxxx = cbuffer.data(sh_off + 0 * ccomps * dcomps);

            auto g_0_xxxxy = cbuffer.data(sh_off + 1 * ccomps * dcomps);

            auto g_0_xxxxz = cbuffer.data(sh_off + 2 * ccomps * dcomps);

            auto g_0_xxxyy = cbuffer.data(sh_off + 3 * ccomps * dcomps);

            auto g_0_xxxyz = cbuffer.data(sh_off + 4 * ccomps * dcomps);

            auto g_0_xxxzz = cbuffer.data(sh_off + 5 * ccomps * dcomps);

            auto g_0_xxyyy = cbuffer.data(sh_off + 6 * ccomps * dcomps);

            auto g_0_xxyyz = cbuffer.data(sh_off + 7 * ccomps * dcomps);

            auto g_0_xxyzz = cbuffer.data(sh_off + 8 * ccomps * dcomps);

            auto g_0_xxzzz = cbuffer.data(sh_off + 9 * ccomps * dcomps);

            auto g_0_xyyyy = cbuffer.data(sh_off + 10 * ccomps * dcomps);

            auto g_0_xyyyz = cbuffer.data(sh_off + 11 * ccomps * dcomps);

            auto g_0_xyyzz = cbuffer.data(sh_off + 12 * ccomps * dcomps);

            auto g_0_xyzzz = cbuffer.data(sh_off + 13 * ccomps * dcomps);

            auto g_0_xzzzz = cbuffer.data(sh_off + 14 * ccomps * dcomps);

            auto g_0_yyyyy = cbuffer.data(sh_off + 15 * ccomps * dcomps);

            auto g_0_yyyyz = cbuffer.data(sh_off + 16 * ccomps * dcomps);

            auto g_0_yyyzz = cbuffer.data(sh_off + 17 * ccomps * dcomps);

            auto g_0_yyzzz = cbuffer.data(sh_off + 18 * ccomps * dcomps);

            auto g_0_yzzzz = cbuffer.data(sh_off + 19 * ccomps * dcomps);

            auto g_0_zzzzz = cbuffer.data(sh_off + 20 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pgxx

            const auto pg_off = idx_pgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_xxxx = cbuffer.data(pg_off + 0 * ccomps * dcomps);

            auto g_x_xxxy = cbuffer.data(pg_off + 1 * ccomps * dcomps);

            auto g_x_xxxz = cbuffer.data(pg_off + 2 * ccomps * dcomps);

            auto g_x_xxyy = cbuffer.data(pg_off + 3 * ccomps * dcomps);

            auto g_x_xxyz = cbuffer.data(pg_off + 4 * ccomps * dcomps);

            auto g_x_xxzz = cbuffer.data(pg_off + 5 * ccomps * dcomps);

            auto g_x_xyyy = cbuffer.data(pg_off + 6 * ccomps * dcomps);

            auto g_x_xyyz = cbuffer.data(pg_off + 7 * ccomps * dcomps);

            auto g_x_xyzz = cbuffer.data(pg_off + 8 * ccomps * dcomps);

            auto g_x_xzzz = cbuffer.data(pg_off + 9 * ccomps * dcomps);

            auto g_x_yyyy = cbuffer.data(pg_off + 10 * ccomps * dcomps);

            auto g_x_yyyz = cbuffer.data(pg_off + 11 * ccomps * dcomps);

            auto g_x_yyzz = cbuffer.data(pg_off + 12 * ccomps * dcomps);

            auto g_x_yzzz = cbuffer.data(pg_off + 13 * ccomps * dcomps);

            auto g_x_zzzz = cbuffer.data(pg_off + 14 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxx,      \
                             g_0_xxxxx, \
                             g_0_xxxxy, \
                             g_0_xxxxz, \
                             g_0_xxxy,  \
                             g_0_xxxyy, \
                             g_0_xxxyz, \
                             g_0_xxxz,  \
                             g_0_xxxzz, \
                             g_0_xxyy,  \
                             g_0_xxyyy, \
                             g_0_xxyyz, \
                             g_0_xxyz,  \
                             g_0_xxyzz, \
                             g_0_xxzz,  \
                             g_0_xxzzz, \
                             g_0_xyyy,  \
                             g_0_xyyyy, \
                             g_0_xyyyz, \
                             g_0_xyyz,  \
                             g_0_xyyzz, \
                             g_0_xyzz,  \
                             g_0_xyzzz, \
                             g_0_xzzz,  \
                             g_0_xzzzz, \
                             g_0_yyyy,  \
                             g_0_yyyz,  \
                             g_0_yyzz,  \
                             g_0_yzzz,  \
                             g_0_zzzz,  \
                             g_x_xxxx,  \
                             g_x_xxxy,  \
                             g_x_xxxz,  \
                             g_x_xxyy,  \
                             g_x_xxyz,  \
                             g_x_xxzz,  \
                             g_x_xyyy,  \
                             g_x_xyyz,  \
                             g_x_xyzz,  \
                             g_x_xzzz,  \
                             g_x_yyyy,  \
                             g_x_yyyz,  \
                             g_x_yyzz,  \
                             g_x_yzzz,  \
                             g_x_zzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xxxx[k] = -g_0_xxxx[k] * ab_x + g_0_xxxxx[k];

                g_x_xxxy[k] = -g_0_xxxy[k] * ab_x + g_0_xxxxy[k];

                g_x_xxxz[k] = -g_0_xxxz[k] * ab_x + g_0_xxxxz[k];

                g_x_xxyy[k] = -g_0_xxyy[k] * ab_x + g_0_xxxyy[k];

                g_x_xxyz[k] = -g_0_xxyz[k] * ab_x + g_0_xxxyz[k];

                g_x_xxzz[k] = -g_0_xxzz[k] * ab_x + g_0_xxxzz[k];

                g_x_xyyy[k] = -g_0_xyyy[k] * ab_x + g_0_xxyyy[k];

                g_x_xyyz[k] = -g_0_xyyz[k] * ab_x + g_0_xxyyz[k];

                g_x_xyzz[k] = -g_0_xyzz[k] * ab_x + g_0_xxyzz[k];

                g_x_xzzz[k] = -g_0_xzzz[k] * ab_x + g_0_xxzzz[k];

                g_x_yyyy[k] = -g_0_yyyy[k] * ab_x + g_0_xyyyy[k];

                g_x_yyyz[k] = -g_0_yyyz[k] * ab_x + g_0_xyyyz[k];

                g_x_yyzz[k] = -g_0_yyzz[k] * ab_x + g_0_xyyzz[k];

                g_x_yzzz[k] = -g_0_yzzz[k] * ab_x + g_0_xyzzz[k];

                g_x_zzzz[k] = -g_0_zzzz[k] * ab_x + g_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_y_xxxx = cbuffer.data(pg_off + 15 * ccomps * dcomps);

            auto g_y_xxxy = cbuffer.data(pg_off + 16 * ccomps * dcomps);

            auto g_y_xxxz = cbuffer.data(pg_off + 17 * ccomps * dcomps);

            auto g_y_xxyy = cbuffer.data(pg_off + 18 * ccomps * dcomps);

            auto g_y_xxyz = cbuffer.data(pg_off + 19 * ccomps * dcomps);

            auto g_y_xxzz = cbuffer.data(pg_off + 20 * ccomps * dcomps);

            auto g_y_xyyy = cbuffer.data(pg_off + 21 * ccomps * dcomps);

            auto g_y_xyyz = cbuffer.data(pg_off + 22 * ccomps * dcomps);

            auto g_y_xyzz = cbuffer.data(pg_off + 23 * ccomps * dcomps);

            auto g_y_xzzz = cbuffer.data(pg_off + 24 * ccomps * dcomps);

            auto g_y_yyyy = cbuffer.data(pg_off + 25 * ccomps * dcomps);

            auto g_y_yyyz = cbuffer.data(pg_off + 26 * ccomps * dcomps);

            auto g_y_yyzz = cbuffer.data(pg_off + 27 * ccomps * dcomps);

            auto g_y_yzzz = cbuffer.data(pg_off + 28 * ccomps * dcomps);

            auto g_y_zzzz = cbuffer.data(pg_off + 29 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxx,      \
                             g_0_xxxxy, \
                             g_0_xxxy,  \
                             g_0_xxxyy, \
                             g_0_xxxyz, \
                             g_0_xxxz,  \
                             g_0_xxyy,  \
                             g_0_xxyyy, \
                             g_0_xxyyz, \
                             g_0_xxyz,  \
                             g_0_xxyzz, \
                             g_0_xxzz,  \
                             g_0_xyyy,  \
                             g_0_xyyyy, \
                             g_0_xyyyz, \
                             g_0_xyyz,  \
                             g_0_xyyzz, \
                             g_0_xyzz,  \
                             g_0_xyzzz, \
                             g_0_xzzz,  \
                             g_0_yyyy,  \
                             g_0_yyyyy, \
                             g_0_yyyyz, \
                             g_0_yyyz,  \
                             g_0_yyyzz, \
                             g_0_yyzz,  \
                             g_0_yyzzz, \
                             g_0_yzzz,  \
                             g_0_yzzzz, \
                             g_0_zzzz,  \
                             g_y_xxxx,  \
                             g_y_xxxy,  \
                             g_y_xxxz,  \
                             g_y_xxyy,  \
                             g_y_xxyz,  \
                             g_y_xxzz,  \
                             g_y_xyyy,  \
                             g_y_xyyz,  \
                             g_y_xyzz,  \
                             g_y_xzzz,  \
                             g_y_yyyy,  \
                             g_y_yyyz,  \
                             g_y_yyzz,  \
                             g_y_yzzz,  \
                             g_y_zzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xxxx[k] = -g_0_xxxx[k] * ab_y + g_0_xxxxy[k];

                g_y_xxxy[k] = -g_0_xxxy[k] * ab_y + g_0_xxxyy[k];

                g_y_xxxz[k] = -g_0_xxxz[k] * ab_y + g_0_xxxyz[k];

                g_y_xxyy[k] = -g_0_xxyy[k] * ab_y + g_0_xxyyy[k];

                g_y_xxyz[k] = -g_0_xxyz[k] * ab_y + g_0_xxyyz[k];

                g_y_xxzz[k] = -g_0_xxzz[k] * ab_y + g_0_xxyzz[k];

                g_y_xyyy[k] = -g_0_xyyy[k] * ab_y + g_0_xyyyy[k];

                g_y_xyyz[k] = -g_0_xyyz[k] * ab_y + g_0_xyyyz[k];

                g_y_xyzz[k] = -g_0_xyzz[k] * ab_y + g_0_xyyzz[k];

                g_y_xzzz[k] = -g_0_xzzz[k] * ab_y + g_0_xyzzz[k];

                g_y_yyyy[k] = -g_0_yyyy[k] * ab_y + g_0_yyyyy[k];

                g_y_yyyz[k] = -g_0_yyyz[k] * ab_y + g_0_yyyyz[k];

                g_y_yyzz[k] = -g_0_yyzz[k] * ab_y + g_0_yyyzz[k];

                g_y_yzzz[k] = -g_0_yzzz[k] * ab_y + g_0_yyzzz[k];

                g_y_zzzz[k] = -g_0_zzzz[k] * ab_y + g_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_z_xxxx = cbuffer.data(pg_off + 30 * ccomps * dcomps);

            auto g_z_xxxy = cbuffer.data(pg_off + 31 * ccomps * dcomps);

            auto g_z_xxxz = cbuffer.data(pg_off + 32 * ccomps * dcomps);

            auto g_z_xxyy = cbuffer.data(pg_off + 33 * ccomps * dcomps);

            auto g_z_xxyz = cbuffer.data(pg_off + 34 * ccomps * dcomps);

            auto g_z_xxzz = cbuffer.data(pg_off + 35 * ccomps * dcomps);

            auto g_z_xyyy = cbuffer.data(pg_off + 36 * ccomps * dcomps);

            auto g_z_xyyz = cbuffer.data(pg_off + 37 * ccomps * dcomps);

            auto g_z_xyzz = cbuffer.data(pg_off + 38 * ccomps * dcomps);

            auto g_z_xzzz = cbuffer.data(pg_off + 39 * ccomps * dcomps);

            auto g_z_yyyy = cbuffer.data(pg_off + 40 * ccomps * dcomps);

            auto g_z_yyyz = cbuffer.data(pg_off + 41 * ccomps * dcomps);

            auto g_z_yyzz = cbuffer.data(pg_off + 42 * ccomps * dcomps);

            auto g_z_yzzz = cbuffer.data(pg_off + 43 * ccomps * dcomps);

            auto g_z_zzzz = cbuffer.data(pg_off + 44 * ccomps * dcomps);

#pragma omp simd aligned(g_0_xxxx,      \
                             g_0_xxxxz, \
                             g_0_xxxy,  \
                             g_0_xxxyz, \
                             g_0_xxxz,  \
                             g_0_xxxzz, \
                             g_0_xxyy,  \
                             g_0_xxyyz, \
                             g_0_xxyz,  \
                             g_0_xxyzz, \
                             g_0_xxzz,  \
                             g_0_xxzzz, \
                             g_0_xyyy,  \
                             g_0_xyyyz, \
                             g_0_xyyz,  \
                             g_0_xyyzz, \
                             g_0_xyzz,  \
                             g_0_xyzzz, \
                             g_0_xzzz,  \
                             g_0_xzzzz, \
                             g_0_yyyy,  \
                             g_0_yyyyz, \
                             g_0_yyyz,  \
                             g_0_yyyzz, \
                             g_0_yyzz,  \
                             g_0_yyzzz, \
                             g_0_yzzz,  \
                             g_0_yzzzz, \
                             g_0_zzzz,  \
                             g_0_zzzzz, \
                             g_z_xxxx,  \
                             g_z_xxxy,  \
                             g_z_xxxz,  \
                             g_z_xxyy,  \
                             g_z_xxyz,  \
                             g_z_xxzz,  \
                             g_z_xyyy,  \
                             g_z_xyyz,  \
                             g_z_xyzz,  \
                             g_z_xzzz,  \
                             g_z_yyyy,  \
                             g_z_yyyz,  \
                             g_z_yyzz,  \
                             g_z_yzzz,  \
                             g_z_zzzz : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xxxx[k] = -g_0_xxxx[k] * ab_z + g_0_xxxxz[k];

                g_z_xxxy[k] = -g_0_xxxy[k] * ab_z + g_0_xxxyz[k];

                g_z_xxxz[k] = -g_0_xxxz[k] * ab_z + g_0_xxxzz[k];

                g_z_xxyy[k] = -g_0_xxyy[k] * ab_z + g_0_xxyyz[k];

                g_z_xxyz[k] = -g_0_xxyz[k] * ab_z + g_0_xxyzz[k];

                g_z_xxzz[k] = -g_0_xxzz[k] * ab_z + g_0_xxzzz[k];

                g_z_xyyy[k] = -g_0_xyyy[k] * ab_z + g_0_xyyyz[k];

                g_z_xyyz[k] = -g_0_xyyz[k] * ab_z + g_0_xyyzz[k];

                g_z_xyzz[k] = -g_0_xyzz[k] * ab_z + g_0_xyzzz[k];

                g_z_xzzz[k] = -g_0_xzzz[k] * ab_z + g_0_xzzzz[k];

                g_z_yyyy[k] = -g_0_yyyy[k] * ab_z + g_0_yyyyz[k];

                g_z_yyyz[k] = -g_0_yyyz[k] * ab_z + g_0_yyyzz[k];

                g_z_yyzz[k] = -g_0_yyzz[k] * ab_z + g_0_yyzzz[k];

                g_z_yzzz[k] = -g_0_yzzz[k] * ab_z + g_0_yzzzz[k];

                g_z_zzzz[k] = -g_0_zzzz[k] * ab_z + g_0_zzzzz[k];
            }
        }
    }
}

}  // namespace erirec

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

#include "ElectronRepulsionPrimRecSSSG.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sssg(CSimdArray<double>& pbuffer,
                                  const size_t        idx_eri_0_sssg,
                                  size_t              idx_eri_0_sssd,
                                  size_t              idx_eri_1_sssd,
                                  size_t              idx_eri_0_sssf,
                                  size_t              idx_eri_1_sssf,
                                  CSimdArray<double>& factors,
                                  const size_t        idx_qd,
                                  const size_t        idx_wq,
                                  const double        a_exp,
                                  const double        b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSSD

    auto g_0_0_0_xx_0 = pbuffer.data(idx_eri_0_sssd);

    auto g_0_0_0_yy_0 = pbuffer.data(idx_eri_0_sssd + 3);

    auto g_0_0_0_zz_0 = pbuffer.data(idx_eri_0_sssd + 5);

    /// Set up components of auxilary buffer : SSSD

    auto g_0_0_0_xx_1 = pbuffer.data(idx_eri_1_sssd);

    auto g_0_0_0_yy_1 = pbuffer.data(idx_eri_1_sssd + 3);

    auto g_0_0_0_zz_1 = pbuffer.data(idx_eri_1_sssd + 5);

    /// Set up components of auxilary buffer : SSSF

    auto g_0_0_0_xxx_0 = pbuffer.data(idx_eri_0_sssf);

    auto g_0_0_0_xxz_0 = pbuffer.data(idx_eri_0_sssf + 2);

    auto g_0_0_0_xyy_0 = pbuffer.data(idx_eri_0_sssf + 3);

    auto g_0_0_0_xzz_0 = pbuffer.data(idx_eri_0_sssf + 5);

    auto g_0_0_0_yyy_0 = pbuffer.data(idx_eri_0_sssf + 6);

    auto g_0_0_0_yyz_0 = pbuffer.data(idx_eri_0_sssf + 7);

    auto g_0_0_0_yzz_0 = pbuffer.data(idx_eri_0_sssf + 8);

    auto g_0_0_0_zzz_0 = pbuffer.data(idx_eri_0_sssf + 9);

    /// Set up components of auxilary buffer : SSSF

    auto g_0_0_0_xxx_1 = pbuffer.data(idx_eri_1_sssf);

    auto g_0_0_0_xxz_1 = pbuffer.data(idx_eri_1_sssf + 2);

    auto g_0_0_0_xyy_1 = pbuffer.data(idx_eri_1_sssf + 3);

    auto g_0_0_0_xzz_1 = pbuffer.data(idx_eri_1_sssf + 5);

    auto g_0_0_0_yyy_1 = pbuffer.data(idx_eri_1_sssf + 6);

    auto g_0_0_0_yyz_1 = pbuffer.data(idx_eri_1_sssf + 7);

    auto g_0_0_0_yzz_1 = pbuffer.data(idx_eri_1_sssf + 8);

    auto g_0_0_0_zzz_1 = pbuffer.data(idx_eri_1_sssf + 9);

    /// Set up components of targeted buffer : SSSG

    auto g_0_0_0_xxxx_0 = pbuffer.data(idx_eri_0_sssg);

    auto g_0_0_0_xxxy_0 = pbuffer.data(idx_eri_0_sssg + 1);

    auto g_0_0_0_xxxz_0 = pbuffer.data(idx_eri_0_sssg + 2);

    auto g_0_0_0_xxyy_0 = pbuffer.data(idx_eri_0_sssg + 3);

    auto g_0_0_0_xxyz_0 = pbuffer.data(idx_eri_0_sssg + 4);

    auto g_0_0_0_xxzz_0 = pbuffer.data(idx_eri_0_sssg + 5);

    auto g_0_0_0_xyyy_0 = pbuffer.data(idx_eri_0_sssg + 6);

    auto g_0_0_0_xyyz_0 = pbuffer.data(idx_eri_0_sssg + 7);

    auto g_0_0_0_xyzz_0 = pbuffer.data(idx_eri_0_sssg + 8);

    auto g_0_0_0_xzzz_0 = pbuffer.data(idx_eri_0_sssg + 9);

    auto g_0_0_0_yyyy_0 = pbuffer.data(idx_eri_0_sssg + 10);

    auto g_0_0_0_yyyz_0 = pbuffer.data(idx_eri_0_sssg + 11);

    auto g_0_0_0_yyzz_0 = pbuffer.data(idx_eri_0_sssg + 12);

    auto g_0_0_0_yzzz_0 = pbuffer.data(idx_eri_0_sssg + 13);

    auto g_0_0_0_zzzz_0 = pbuffer.data(idx_eri_0_sssg + 14);

#pragma omp simd aligned(g_0_0_0_xx_0,       \
                             g_0_0_0_xx_1,   \
                             g_0_0_0_xxx_0,  \
                             g_0_0_0_xxx_1,  \
                             g_0_0_0_xxxx_0, \
                             g_0_0_0_xxxy_0, \
                             g_0_0_0_xxxz_0, \
                             g_0_0_0_xxyy_0, \
                             g_0_0_0_xxyz_0, \
                             g_0_0_0_xxz_0,  \
                             g_0_0_0_xxz_1,  \
                             g_0_0_0_xxzz_0, \
                             g_0_0_0_xyy_0,  \
                             g_0_0_0_xyy_1,  \
                             g_0_0_0_xyyy_0, \
                             g_0_0_0_xyyz_0, \
                             g_0_0_0_xyzz_0, \
                             g_0_0_0_xzz_0,  \
                             g_0_0_0_xzz_1,  \
                             g_0_0_0_xzzz_0, \
                             g_0_0_0_yy_0,   \
                             g_0_0_0_yy_1,   \
                             g_0_0_0_yyy_0,  \
                             g_0_0_0_yyy_1,  \
                             g_0_0_0_yyyy_0, \
                             g_0_0_0_yyyz_0, \
                             g_0_0_0_yyz_0,  \
                             g_0_0_0_yyz_1,  \
                             g_0_0_0_yyzz_0, \
                             g_0_0_0_yzz_0,  \
                             g_0_0_0_yzz_1,  \
                             g_0_0_0_yzzz_0, \
                             g_0_0_0_zz_0,   \
                             g_0_0_0_zz_1,   \
                             g_0_0_0_zzz_0,  \
                             g_0_0_0_zzz_1,  \
                             g_0_0_0_zzzz_0, \
                             qd_x,           \
                             qd_y,           \
                             qd_z,           \
                             wq_x,           \
                             wq_y,           \
                             wq_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 = fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxx_0[i] =
            3.0 * g_0_0_0_xx_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xx_1[i] * fti_cd_0 + g_0_0_0_xxx_0[i] * qd_x[i] + g_0_0_0_xxx_1[i] * wq_x[i];

        g_0_0_0_xxxy_0[i] = g_0_0_0_xxx_0[i] * qd_y[i] + g_0_0_0_xxx_1[i] * wq_y[i];

        g_0_0_0_xxxz_0[i] = g_0_0_0_xxx_0[i] * qd_z[i] + g_0_0_0_xxx_1[i] * wq_z[i];

        g_0_0_0_xxyy_0[i] = g_0_0_0_yy_0[i] * fi_cd_0 - g_0_0_0_yy_1[i] * fti_cd_0 + g_0_0_0_xyy_0[i] * qd_x[i] + g_0_0_0_xyy_1[i] * wq_x[i];

        g_0_0_0_xxyz_0[i] = g_0_0_0_xxz_0[i] * qd_y[i] + g_0_0_0_xxz_1[i] * wq_y[i];

        g_0_0_0_xxzz_0[i] = g_0_0_0_zz_0[i] * fi_cd_0 - g_0_0_0_zz_1[i] * fti_cd_0 + g_0_0_0_xzz_0[i] * qd_x[i] + g_0_0_0_xzz_1[i] * wq_x[i];

        g_0_0_0_xyyy_0[i] = g_0_0_0_yyy_0[i] * qd_x[i] + g_0_0_0_yyy_1[i] * wq_x[i];

        g_0_0_0_xyyz_0[i] = g_0_0_0_yyz_0[i] * qd_x[i] + g_0_0_0_yyz_1[i] * wq_x[i];

        g_0_0_0_xyzz_0[i] = g_0_0_0_yzz_0[i] * qd_x[i] + g_0_0_0_yzz_1[i] * wq_x[i];

        g_0_0_0_xzzz_0[i] = g_0_0_0_zzz_0[i] * qd_x[i] + g_0_0_0_zzz_1[i] * wq_x[i];

        g_0_0_0_yyyy_0[i] =
            3.0 * g_0_0_0_yy_0[i] * fi_cd_0 - 3.0 * g_0_0_0_yy_1[i] * fti_cd_0 + g_0_0_0_yyy_0[i] * qd_y[i] + g_0_0_0_yyy_1[i] * wq_y[i];

        g_0_0_0_yyyz_0[i] = g_0_0_0_yyy_0[i] * qd_z[i] + g_0_0_0_yyy_1[i] * wq_z[i];

        g_0_0_0_yyzz_0[i] = g_0_0_0_zz_0[i] * fi_cd_0 - g_0_0_0_zz_1[i] * fti_cd_0 + g_0_0_0_yzz_0[i] * qd_y[i] + g_0_0_0_yzz_1[i] * wq_y[i];

        g_0_0_0_yzzz_0[i] = g_0_0_0_zzz_0[i] * qd_y[i] + g_0_0_0_zzz_1[i] * wq_y[i];

        g_0_0_0_zzzz_0[i] =
            3.0 * g_0_0_0_zz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_zz_1[i] * fti_cd_0 + g_0_0_0_zzz_0[i] * qd_z[i] + g_0_0_0_zzz_1[i] * wq_z[i];
    }
}

}  // namespace erirec

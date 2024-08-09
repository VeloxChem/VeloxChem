#include "ElectronRepulsionPrimRecSSSH.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssh(CSimdArray<double>& prim_buffer_0_sssh,
                                  const CSimdArray<double>& prim_buffer_0_sssf,
                                  const CSimdArray<double>& prim_buffer_1_sssf,
                                  const CSimdArray<double>& prim_buffer_0_sssg,
                                  const CSimdArray<double>& prim_buffer_1_sssg,
                                  const double* qd_x,
                                  const double* qd_y,
                                  const double* qd_z,
                                  const double* wq_x,
                                  const double* wq_y,
                                  const double* wq_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sssh.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssf

    auto g_0_0_0_xxx_0 = prim_buffer_0_sssf[0];

    auto g_0_0_0_xyy_0 = prim_buffer_0_sssf[3];

    auto g_0_0_0_xzz_0 = prim_buffer_0_sssf[5];

    auto g_0_0_0_yyy_0 = prim_buffer_0_sssf[6];

    auto g_0_0_0_yzz_0 = prim_buffer_0_sssf[8];

    auto g_0_0_0_zzz_0 = prim_buffer_0_sssf[9];

    /// Set up components of auxilary buffer : prim_buffer_1_sssf

    auto g_0_0_0_xxx_1 = prim_buffer_1_sssf[0];

    auto g_0_0_0_xyy_1 = prim_buffer_1_sssf[3];

    auto g_0_0_0_xzz_1 = prim_buffer_1_sssf[5];

    auto g_0_0_0_yyy_1 = prim_buffer_1_sssf[6];

    auto g_0_0_0_yzz_1 = prim_buffer_1_sssf[8];

    auto g_0_0_0_zzz_1 = prim_buffer_1_sssf[9];

    /// Set up components of auxilary buffer : prim_buffer_0_sssg

    auto g_0_0_0_xxxx_0 = prim_buffer_0_sssg[0];

    auto g_0_0_0_xxxz_0 = prim_buffer_0_sssg[2];

    auto g_0_0_0_xxyy_0 = prim_buffer_0_sssg[3];

    auto g_0_0_0_xxzz_0 = prim_buffer_0_sssg[5];

    auto g_0_0_0_xyyy_0 = prim_buffer_0_sssg[6];

    auto g_0_0_0_xzzz_0 = prim_buffer_0_sssg[9];

    auto g_0_0_0_yyyy_0 = prim_buffer_0_sssg[10];

    auto g_0_0_0_yyyz_0 = prim_buffer_0_sssg[11];

    auto g_0_0_0_yyzz_0 = prim_buffer_0_sssg[12];

    auto g_0_0_0_yzzz_0 = prim_buffer_0_sssg[13];

    auto g_0_0_0_zzzz_0 = prim_buffer_0_sssg[14];

    /// Set up components of auxilary buffer : prim_buffer_1_sssg

    auto g_0_0_0_xxxx_1 = prim_buffer_1_sssg[0];

    auto g_0_0_0_xxxz_1 = prim_buffer_1_sssg[2];

    auto g_0_0_0_xxyy_1 = prim_buffer_1_sssg[3];

    auto g_0_0_0_xxzz_1 = prim_buffer_1_sssg[5];

    auto g_0_0_0_xyyy_1 = prim_buffer_1_sssg[6];

    auto g_0_0_0_xzzz_1 = prim_buffer_1_sssg[9];

    auto g_0_0_0_yyyy_1 = prim_buffer_1_sssg[10];

    auto g_0_0_0_yyyz_1 = prim_buffer_1_sssg[11];

    auto g_0_0_0_yyzz_1 = prim_buffer_1_sssg[12];

    auto g_0_0_0_yzzz_1 = prim_buffer_1_sssg[13];

    auto g_0_0_0_zzzz_1 = prim_buffer_1_sssg[14];

    /// Set up components of targeted buffer : prim_buffer_0_sssh

    auto g_0_0_0_xxxxx_0 = prim_buffer_0_sssh[0];

    auto g_0_0_0_xxxxy_0 = prim_buffer_0_sssh[1];

    auto g_0_0_0_xxxxz_0 = prim_buffer_0_sssh[2];

    auto g_0_0_0_xxxyy_0 = prim_buffer_0_sssh[3];

    auto g_0_0_0_xxxyz_0 = prim_buffer_0_sssh[4];

    auto g_0_0_0_xxxzz_0 = prim_buffer_0_sssh[5];

    auto g_0_0_0_xxyyy_0 = prim_buffer_0_sssh[6];

    auto g_0_0_0_xxyyz_0 = prim_buffer_0_sssh[7];

    auto g_0_0_0_xxyzz_0 = prim_buffer_0_sssh[8];

    auto g_0_0_0_xxzzz_0 = prim_buffer_0_sssh[9];

    auto g_0_0_0_xyyyy_0 = prim_buffer_0_sssh[10];

    auto g_0_0_0_xyyyz_0 = prim_buffer_0_sssh[11];

    auto g_0_0_0_xyyzz_0 = prim_buffer_0_sssh[12];

    auto g_0_0_0_xyzzz_0 = prim_buffer_0_sssh[13];

    auto g_0_0_0_xzzzz_0 = prim_buffer_0_sssh[14];

    auto g_0_0_0_yyyyy_0 = prim_buffer_0_sssh[15];

    auto g_0_0_0_yyyyz_0 = prim_buffer_0_sssh[16];

    auto g_0_0_0_yyyzz_0 = prim_buffer_0_sssh[17];

    auto g_0_0_0_yyzzz_0 = prim_buffer_0_sssh[18];

    auto g_0_0_0_yzzzz_0 = prim_buffer_0_sssh[19];

    auto g_0_0_0_zzzzz_0 = prim_buffer_0_sssh[20];

    #pragma omp simd aligned(g_0_0_0_xxx_0, g_0_0_0_xxx_1, g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxy_0, g_0_0_0_xxxxz_0, g_0_0_0_xxxyy_0, g_0_0_0_xxxyz_0, g_0_0_0_xxxz_0, g_0_0_0_xxxz_1, g_0_0_0_xxxzz_0, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyz_0, g_0_0_0_xxyzz_0, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xyy_0, g_0_0_0_xyy_1, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyz_0, g_0_0_0_xyyzz_0, g_0_0_0_xyzzz_0, g_0_0_0_xzz_0, g_0_0_0_xzz_1, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_yyy_0, g_0_0_0_yyy_1, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyz_0, g_0_0_0_yyyz_0, g_0_0_0_yyyz_1, g_0_0_0_yyyzz_0, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yzz_0, g_0_0_0_yzz_1, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_zzz_0, g_0_0_0_zzz_1, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxxx_0[i] = 4.0 * g_0_0_0_xxx_0[i] * fi_cd_0 - 4.0 * g_0_0_0_xxx_1[i] * fti_cd_0 + g_0_0_0_xxxx_0[i] * qd_x[i] + g_0_0_0_xxxx_1[i] * wq_x[i];

        g_0_0_0_xxxxy_0[i] = g_0_0_0_xxxx_0[i] * qd_y[i] + g_0_0_0_xxxx_1[i] * wq_y[i];

        g_0_0_0_xxxxz_0[i] = g_0_0_0_xxxx_0[i] * qd_z[i] + g_0_0_0_xxxx_1[i] * wq_z[i];

        g_0_0_0_xxxyy_0[i] = 2.0 * g_0_0_0_xyy_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyy_1[i] * fti_cd_0 + g_0_0_0_xxyy_0[i] * qd_x[i] + g_0_0_0_xxyy_1[i] * wq_x[i];

        g_0_0_0_xxxyz_0[i] = g_0_0_0_xxxz_0[i] * qd_y[i] + g_0_0_0_xxxz_1[i] * wq_y[i];

        g_0_0_0_xxxzz_0[i] = 2.0 * g_0_0_0_xzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xzz_1[i] * fti_cd_0 + g_0_0_0_xxzz_0[i] * qd_x[i] + g_0_0_0_xxzz_1[i] * wq_x[i];

        g_0_0_0_xxyyy_0[i] = g_0_0_0_yyy_0[i] * fi_cd_0 - g_0_0_0_yyy_1[i] * fti_cd_0 + g_0_0_0_xyyy_0[i] * qd_x[i] + g_0_0_0_xyyy_1[i] * wq_x[i];

        g_0_0_0_xxyyz_0[i] = g_0_0_0_xxyy_0[i] * qd_z[i] + g_0_0_0_xxyy_1[i] * wq_z[i];

        g_0_0_0_xxyzz_0[i] = g_0_0_0_xxzz_0[i] * qd_y[i] + g_0_0_0_xxzz_1[i] * wq_y[i];

        g_0_0_0_xxzzz_0[i] = g_0_0_0_zzz_0[i] * fi_cd_0 - g_0_0_0_zzz_1[i] * fti_cd_0 + g_0_0_0_xzzz_0[i] * qd_x[i] + g_0_0_0_xzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyy_0[i] = g_0_0_0_yyyy_0[i] * qd_x[i] + g_0_0_0_yyyy_1[i] * wq_x[i];

        g_0_0_0_xyyyz_0[i] = g_0_0_0_yyyz_0[i] * qd_x[i] + g_0_0_0_yyyz_1[i] * wq_x[i];

        g_0_0_0_xyyzz_0[i] = g_0_0_0_yyzz_0[i] * qd_x[i] + g_0_0_0_yyzz_1[i] * wq_x[i];

        g_0_0_0_xyzzz_0[i] = g_0_0_0_yzzz_0[i] * qd_x[i] + g_0_0_0_yzzz_1[i] * wq_x[i];

        g_0_0_0_xzzzz_0[i] = g_0_0_0_zzzz_0[i] * qd_x[i] + g_0_0_0_zzzz_1[i] * wq_x[i];

        g_0_0_0_yyyyy_0[i] = 4.0 * g_0_0_0_yyy_0[i] * fi_cd_0 - 4.0 * g_0_0_0_yyy_1[i] * fti_cd_0 + g_0_0_0_yyyy_0[i] * qd_y[i] + g_0_0_0_yyyy_1[i] * wq_y[i];

        g_0_0_0_yyyyz_0[i] = g_0_0_0_yyyy_0[i] * qd_z[i] + g_0_0_0_yyyy_1[i] * wq_z[i];

        g_0_0_0_yyyzz_0[i] = 2.0 * g_0_0_0_yzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_yzz_1[i] * fti_cd_0 + g_0_0_0_yyzz_0[i] * qd_y[i] + g_0_0_0_yyzz_1[i] * wq_y[i];

        g_0_0_0_yyzzz_0[i] = g_0_0_0_zzz_0[i] * fi_cd_0 - g_0_0_0_zzz_1[i] * fti_cd_0 + g_0_0_0_yzzz_0[i] * qd_y[i] + g_0_0_0_yzzz_1[i] * wq_y[i];

        g_0_0_0_yzzzz_0[i] = g_0_0_0_zzzz_0[i] * qd_y[i] + g_0_0_0_zzzz_1[i] * wq_y[i];

        g_0_0_0_zzzzz_0[i] = 4.0 * g_0_0_0_zzz_0[i] * fi_cd_0 - 4.0 * g_0_0_0_zzz_1[i] * fti_cd_0 + g_0_0_0_zzzz_0[i] * qd_z[i] + g_0_0_0_zzzz_1[i] * wq_z[i];
    }
}

} // erirec namespace


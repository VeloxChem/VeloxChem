#include "ElectronRepulsionPrimRecSSSI.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssi(CSimdArray<double>& prim_buffer_0_sssi,
                                  const CSimdArray<double>& prim_buffer_0_sssg,
                                  const CSimdArray<double>& prim_buffer_1_sssg,
                                  const CSimdArray<double>& prim_buffer_0_sssh,
                                  const CSimdArray<double>& prim_buffer_1_sssh,
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
    const auto ndims = prim_buffer_0_sssi.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssg

    auto g_0_0_0_xxxx_0 = prim_buffer_0_sssg[0];

    auto g_0_0_0_xxyy_0 = prim_buffer_0_sssg[3];

    auto g_0_0_0_xxzz_0 = prim_buffer_0_sssg[5];

    auto g_0_0_0_xyyy_0 = prim_buffer_0_sssg[6];

    auto g_0_0_0_xzzz_0 = prim_buffer_0_sssg[9];

    auto g_0_0_0_yyyy_0 = prim_buffer_0_sssg[10];

    auto g_0_0_0_yyzz_0 = prim_buffer_0_sssg[12];

    auto g_0_0_0_yzzz_0 = prim_buffer_0_sssg[13];

    auto g_0_0_0_zzzz_0 = prim_buffer_0_sssg[14];

    /// Set up components of auxilary buffer : prim_buffer_1_sssg

    auto g_0_0_0_xxxx_1 = prim_buffer_1_sssg[0];

    auto g_0_0_0_xxyy_1 = prim_buffer_1_sssg[3];

    auto g_0_0_0_xxzz_1 = prim_buffer_1_sssg[5];

    auto g_0_0_0_xyyy_1 = prim_buffer_1_sssg[6];

    auto g_0_0_0_xzzz_1 = prim_buffer_1_sssg[9];

    auto g_0_0_0_yyyy_1 = prim_buffer_1_sssg[10];

    auto g_0_0_0_yyzz_1 = prim_buffer_1_sssg[12];

    auto g_0_0_0_yzzz_1 = prim_buffer_1_sssg[13];

    auto g_0_0_0_zzzz_1 = prim_buffer_1_sssg[14];

    /// Set up components of auxilary buffer : prim_buffer_0_sssh

    auto g_0_0_0_xxxxx_0 = prim_buffer_0_sssh[0];

    auto g_0_0_0_xxxxz_0 = prim_buffer_0_sssh[2];

    auto g_0_0_0_xxxyy_0 = prim_buffer_0_sssh[3];

    auto g_0_0_0_xxxzz_0 = prim_buffer_0_sssh[5];

    auto g_0_0_0_xxyyy_0 = prim_buffer_0_sssh[6];

    auto g_0_0_0_xxzzz_0 = prim_buffer_0_sssh[9];

    auto g_0_0_0_xyyyy_0 = prim_buffer_0_sssh[10];

    auto g_0_0_0_xyyzz_0 = prim_buffer_0_sssh[12];

    auto g_0_0_0_xzzzz_0 = prim_buffer_0_sssh[14];

    auto g_0_0_0_yyyyy_0 = prim_buffer_0_sssh[15];

    auto g_0_0_0_yyyyz_0 = prim_buffer_0_sssh[16];

    auto g_0_0_0_yyyzz_0 = prim_buffer_0_sssh[17];

    auto g_0_0_0_yyzzz_0 = prim_buffer_0_sssh[18];

    auto g_0_0_0_yzzzz_0 = prim_buffer_0_sssh[19];

    auto g_0_0_0_zzzzz_0 = prim_buffer_0_sssh[20];

    /// Set up components of auxilary buffer : prim_buffer_1_sssh

    auto g_0_0_0_xxxxx_1 = prim_buffer_1_sssh[0];

    auto g_0_0_0_xxxxz_1 = prim_buffer_1_sssh[2];

    auto g_0_0_0_xxxyy_1 = prim_buffer_1_sssh[3];

    auto g_0_0_0_xxxzz_1 = prim_buffer_1_sssh[5];

    auto g_0_0_0_xxyyy_1 = prim_buffer_1_sssh[6];

    auto g_0_0_0_xxzzz_1 = prim_buffer_1_sssh[9];

    auto g_0_0_0_xyyyy_1 = prim_buffer_1_sssh[10];

    auto g_0_0_0_xyyzz_1 = prim_buffer_1_sssh[12];

    auto g_0_0_0_xzzzz_1 = prim_buffer_1_sssh[14];

    auto g_0_0_0_yyyyy_1 = prim_buffer_1_sssh[15];

    auto g_0_0_0_yyyyz_1 = prim_buffer_1_sssh[16];

    auto g_0_0_0_yyyzz_1 = prim_buffer_1_sssh[17];

    auto g_0_0_0_yyzzz_1 = prim_buffer_1_sssh[18];

    auto g_0_0_0_yzzzz_1 = prim_buffer_1_sssh[19];

    auto g_0_0_0_zzzzz_1 = prim_buffer_1_sssh[20];

    /// Set up components of targeted buffer : prim_buffer_0_sssi

    auto g_0_0_0_xxxxxx_0 = prim_buffer_0_sssi[0];

    auto g_0_0_0_xxxxxy_0 = prim_buffer_0_sssi[1];

    auto g_0_0_0_xxxxxz_0 = prim_buffer_0_sssi[2];

    auto g_0_0_0_xxxxyy_0 = prim_buffer_0_sssi[3];

    auto g_0_0_0_xxxxyz_0 = prim_buffer_0_sssi[4];

    auto g_0_0_0_xxxxzz_0 = prim_buffer_0_sssi[5];

    auto g_0_0_0_xxxyyy_0 = prim_buffer_0_sssi[6];

    auto g_0_0_0_xxxyyz_0 = prim_buffer_0_sssi[7];

    auto g_0_0_0_xxxyzz_0 = prim_buffer_0_sssi[8];

    auto g_0_0_0_xxxzzz_0 = prim_buffer_0_sssi[9];

    auto g_0_0_0_xxyyyy_0 = prim_buffer_0_sssi[10];

    auto g_0_0_0_xxyyyz_0 = prim_buffer_0_sssi[11];

    auto g_0_0_0_xxyyzz_0 = prim_buffer_0_sssi[12];

    auto g_0_0_0_xxyzzz_0 = prim_buffer_0_sssi[13];

    auto g_0_0_0_xxzzzz_0 = prim_buffer_0_sssi[14];

    auto g_0_0_0_xyyyyy_0 = prim_buffer_0_sssi[15];

    auto g_0_0_0_xyyyyz_0 = prim_buffer_0_sssi[16];

    auto g_0_0_0_xyyyzz_0 = prim_buffer_0_sssi[17];

    auto g_0_0_0_xyyzzz_0 = prim_buffer_0_sssi[18];

    auto g_0_0_0_xyzzzz_0 = prim_buffer_0_sssi[19];

    auto g_0_0_0_xzzzzz_0 = prim_buffer_0_sssi[20];

    auto g_0_0_0_yyyyyy_0 = prim_buffer_0_sssi[21];

    auto g_0_0_0_yyyyyz_0 = prim_buffer_0_sssi[22];

    auto g_0_0_0_yyyyzz_0 = prim_buffer_0_sssi[23];

    auto g_0_0_0_yyyzzz_0 = prim_buffer_0_sssi[24];

    auto g_0_0_0_yyzzzz_0 = prim_buffer_0_sssi[25];

    auto g_0_0_0_yzzzzz_0 = prim_buffer_0_sssi[26];

    auto g_0_0_0_zzzzzz_0 = prim_buffer_0_sssi[27];

    #pragma omp simd aligned(g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxx_1, g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxy_0, g_0_0_0_xxxxxz_0, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyz_0, g_0_0_0_xxxxz_0, g_0_0_0_xxxxz_1, g_0_0_0_xxxxzz_0, g_0_0_0_xxxyy_0, g_0_0_0_xxxyy_1, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyz_0, g_0_0_0_xxxyzz_0, g_0_0_0_xxxzz_0, g_0_0_0_xxxzz_1, g_0_0_0_xxxzzz_0, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyy_1, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyz_0, g_0_0_0_xxyyzz_0, g_0_0_0_xxyzzz_0, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xxzzz_1, g_0_0_0_xxzzzz_0, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyy_1, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyz_0, g_0_0_0_xyyyzz_0, g_0_0_0_xyyzz_0, g_0_0_0_xyyzz_1, g_0_0_0_xyyzzz_0, g_0_0_0_xyzzzz_0, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_xzzzz_1, g_0_0_0_xzzzzz_0, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyy_1, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyz_0, g_0_0_0_yyyyz_0, g_0_0_0_yyyyz_1, g_0_0_0_yyyyzz_0, g_0_0_0_yyyzz_0, g_0_0_0_yyyzz_1, g_0_0_0_yyyzzz_0, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yyzzz_1, g_0_0_0_yyzzzz_0, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_yzzzz_1, g_0_0_0_yzzzzz_0, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, g_0_0_0_zzzzz_1, g_0_0_0_zzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxxxx_0[i] = 5.0 * g_0_0_0_xxxx_0[i] * fi_cd_0 - 5.0 * g_0_0_0_xxxx_1[i] * fti_cd_0 + g_0_0_0_xxxxx_0[i] * qd_x[i] + g_0_0_0_xxxxx_1[i] * wq_x[i];

        g_0_0_0_xxxxxy_0[i] = g_0_0_0_xxxxx_0[i] * qd_y[i] + g_0_0_0_xxxxx_1[i] * wq_y[i];

        g_0_0_0_xxxxxz_0[i] = g_0_0_0_xxxxx_0[i] * qd_z[i] + g_0_0_0_xxxxx_1[i] * wq_z[i];

        g_0_0_0_xxxxyy_0[i] = 3.0 * g_0_0_0_xxyy_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxyy_1[i] * fti_cd_0 + g_0_0_0_xxxyy_0[i] * qd_x[i] + g_0_0_0_xxxyy_1[i] * wq_x[i];

        g_0_0_0_xxxxyz_0[i] = g_0_0_0_xxxxz_0[i] * qd_y[i] + g_0_0_0_xxxxz_1[i] * wq_y[i];

        g_0_0_0_xxxxzz_0[i] = 3.0 * g_0_0_0_xxzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxzz_1[i] * fti_cd_0 + g_0_0_0_xxxzz_0[i] * qd_x[i] + g_0_0_0_xxxzz_1[i] * wq_x[i];

        g_0_0_0_xxxyyy_0[i] = 2.0 * g_0_0_0_xyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyy_1[i] * fti_cd_0 + g_0_0_0_xxyyy_0[i] * qd_x[i] + g_0_0_0_xxyyy_1[i] * wq_x[i];

        g_0_0_0_xxxyyz_0[i] = g_0_0_0_xxxyy_0[i] * qd_z[i] + g_0_0_0_xxxyy_1[i] * wq_z[i];

        g_0_0_0_xxxyzz_0[i] = g_0_0_0_xxxzz_0[i] * qd_y[i] + g_0_0_0_xxxzz_1[i] * wq_y[i];

        g_0_0_0_xxxzzz_0[i] = 2.0 * g_0_0_0_xzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xzzz_1[i] * fti_cd_0 + g_0_0_0_xxzzz_0[i] * qd_x[i] + g_0_0_0_xxzzz_1[i] * wq_x[i];

        g_0_0_0_xxyyyy_0[i] = g_0_0_0_yyyy_0[i] * fi_cd_0 - g_0_0_0_yyyy_1[i] * fti_cd_0 + g_0_0_0_xyyyy_0[i] * qd_x[i] + g_0_0_0_xyyyy_1[i] * wq_x[i];

        g_0_0_0_xxyyyz_0[i] = g_0_0_0_xxyyy_0[i] * qd_z[i] + g_0_0_0_xxyyy_1[i] * wq_z[i];

        g_0_0_0_xxyyzz_0[i] = g_0_0_0_yyzz_0[i] * fi_cd_0 - g_0_0_0_yyzz_1[i] * fti_cd_0 + g_0_0_0_xyyzz_0[i] * qd_x[i] + g_0_0_0_xyyzz_1[i] * wq_x[i];

        g_0_0_0_xxyzzz_0[i] = g_0_0_0_xxzzz_0[i] * qd_y[i] + g_0_0_0_xxzzz_1[i] * wq_y[i];

        g_0_0_0_xxzzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_xzzzz_0[i] * qd_x[i] + g_0_0_0_xzzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyyy_0[i] = g_0_0_0_yyyyy_0[i] * qd_x[i] + g_0_0_0_yyyyy_1[i] * wq_x[i];

        g_0_0_0_xyyyyz_0[i] = g_0_0_0_yyyyz_0[i] * qd_x[i] + g_0_0_0_yyyyz_1[i] * wq_x[i];

        g_0_0_0_xyyyzz_0[i] = g_0_0_0_yyyzz_0[i] * qd_x[i] + g_0_0_0_yyyzz_1[i] * wq_x[i];

        g_0_0_0_xyyzzz_0[i] = g_0_0_0_yyzzz_0[i] * qd_x[i] + g_0_0_0_yyzzz_1[i] * wq_x[i];

        g_0_0_0_xyzzzz_0[i] = g_0_0_0_yzzzz_0[i] * qd_x[i] + g_0_0_0_yzzzz_1[i] * wq_x[i];

        g_0_0_0_xzzzzz_0[i] = g_0_0_0_zzzzz_0[i] * qd_x[i] + g_0_0_0_zzzzz_1[i] * wq_x[i];

        g_0_0_0_yyyyyy_0[i] = 5.0 * g_0_0_0_yyyy_0[i] * fi_cd_0 - 5.0 * g_0_0_0_yyyy_1[i] * fti_cd_0 + g_0_0_0_yyyyy_0[i] * qd_y[i] + g_0_0_0_yyyyy_1[i] * wq_y[i];

        g_0_0_0_yyyyyz_0[i] = g_0_0_0_yyyyy_0[i] * qd_z[i] + g_0_0_0_yyyyy_1[i] * wq_z[i];

        g_0_0_0_yyyyzz_0[i] = 3.0 * g_0_0_0_yyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_yyzz_1[i] * fti_cd_0 + g_0_0_0_yyyzz_0[i] * qd_y[i] + g_0_0_0_yyyzz_1[i] * wq_y[i];

        g_0_0_0_yyyzzz_0[i] = 2.0 * g_0_0_0_yzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_yzzz_1[i] * fti_cd_0 + g_0_0_0_yyzzz_0[i] * qd_y[i] + g_0_0_0_yyzzz_1[i] * wq_y[i];

        g_0_0_0_yyzzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_yzzzz_0[i] * qd_y[i] + g_0_0_0_yzzzz_1[i] * wq_y[i];

        g_0_0_0_yzzzzz_0[i] = g_0_0_0_zzzzz_0[i] * qd_y[i] + g_0_0_0_zzzzz_1[i] * wq_y[i];

        g_0_0_0_zzzzzz_0[i] = 5.0 * g_0_0_0_zzzz_0[i] * fi_cd_0 - 5.0 * g_0_0_0_zzzz_1[i] * fti_cd_0 + g_0_0_0_zzzzz_0[i] * qd_z[i] + g_0_0_0_zzzzz_1[i] * wq_z[i];
    }
}

} // erirec namespace


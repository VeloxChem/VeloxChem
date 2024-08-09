#include "ElectronRepulsionPrimRecSISS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_siss(CSimdArray<double>& prim_buffer_0_siss,
                                  const CSimdArray<double>& prim_buffer_0_sgss,
                                  const CSimdArray<double>& prim_buffer_1_sgss,
                                  const CSimdArray<double>& prim_buffer_0_shss,
                                  const CSimdArray<double>& prim_buffer_1_shss,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_siss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sgss

    auto g_0_xxxx_0_0_0 = prim_buffer_0_sgss[0];

    auto g_0_xxyy_0_0_0 = prim_buffer_0_sgss[3];

    auto g_0_xxzz_0_0_0 = prim_buffer_0_sgss[5];

    auto g_0_xyyy_0_0_0 = prim_buffer_0_sgss[6];

    auto g_0_xzzz_0_0_0 = prim_buffer_0_sgss[9];

    auto g_0_yyyy_0_0_0 = prim_buffer_0_sgss[10];

    auto g_0_yyzz_0_0_0 = prim_buffer_0_sgss[12];

    auto g_0_yzzz_0_0_0 = prim_buffer_0_sgss[13];

    auto g_0_zzzz_0_0_0 = prim_buffer_0_sgss[14];

    /// Set up components of auxilary buffer : prim_buffer_1_sgss

    auto g_0_xxxx_0_0_1 = prim_buffer_1_sgss[0];

    auto g_0_xxyy_0_0_1 = prim_buffer_1_sgss[3];

    auto g_0_xxzz_0_0_1 = prim_buffer_1_sgss[5];

    auto g_0_xyyy_0_0_1 = prim_buffer_1_sgss[6];

    auto g_0_xzzz_0_0_1 = prim_buffer_1_sgss[9];

    auto g_0_yyyy_0_0_1 = prim_buffer_1_sgss[10];

    auto g_0_yyzz_0_0_1 = prim_buffer_1_sgss[12];

    auto g_0_yzzz_0_0_1 = prim_buffer_1_sgss[13];

    auto g_0_zzzz_0_0_1 = prim_buffer_1_sgss[14];

    /// Set up components of auxilary buffer : prim_buffer_0_shss

    auto g_0_xxxxx_0_0_0 = prim_buffer_0_shss[0];

    auto g_0_xxxxz_0_0_0 = prim_buffer_0_shss[2];

    auto g_0_xxxyy_0_0_0 = prim_buffer_0_shss[3];

    auto g_0_xxxzz_0_0_0 = prim_buffer_0_shss[5];

    auto g_0_xxyyy_0_0_0 = prim_buffer_0_shss[6];

    auto g_0_xxzzz_0_0_0 = prim_buffer_0_shss[9];

    auto g_0_xyyyy_0_0_0 = prim_buffer_0_shss[10];

    auto g_0_xyyzz_0_0_0 = prim_buffer_0_shss[12];

    auto g_0_xzzzz_0_0_0 = prim_buffer_0_shss[14];

    auto g_0_yyyyy_0_0_0 = prim_buffer_0_shss[15];

    auto g_0_yyyyz_0_0_0 = prim_buffer_0_shss[16];

    auto g_0_yyyzz_0_0_0 = prim_buffer_0_shss[17];

    auto g_0_yyzzz_0_0_0 = prim_buffer_0_shss[18];

    auto g_0_yzzzz_0_0_0 = prim_buffer_0_shss[19];

    auto g_0_zzzzz_0_0_0 = prim_buffer_0_shss[20];

    /// Set up components of auxilary buffer : prim_buffer_1_shss

    auto g_0_xxxxx_0_0_1 = prim_buffer_1_shss[0];

    auto g_0_xxxxz_0_0_1 = prim_buffer_1_shss[2];

    auto g_0_xxxyy_0_0_1 = prim_buffer_1_shss[3];

    auto g_0_xxxzz_0_0_1 = prim_buffer_1_shss[5];

    auto g_0_xxyyy_0_0_1 = prim_buffer_1_shss[6];

    auto g_0_xxzzz_0_0_1 = prim_buffer_1_shss[9];

    auto g_0_xyyyy_0_0_1 = prim_buffer_1_shss[10];

    auto g_0_xyyzz_0_0_1 = prim_buffer_1_shss[12];

    auto g_0_xzzzz_0_0_1 = prim_buffer_1_shss[14];

    auto g_0_yyyyy_0_0_1 = prim_buffer_1_shss[15];

    auto g_0_yyyyz_0_0_1 = prim_buffer_1_shss[16];

    auto g_0_yyyzz_0_0_1 = prim_buffer_1_shss[17];

    auto g_0_yyzzz_0_0_1 = prim_buffer_1_shss[18];

    auto g_0_yzzzz_0_0_1 = prim_buffer_1_shss[19];

    auto g_0_zzzzz_0_0_1 = prim_buffer_1_shss[20];

    /// Set up components of targeted buffer : prim_buffer_0_siss

    auto g_0_xxxxxx_0_0_0 = prim_buffer_0_siss[0];

    auto g_0_xxxxxy_0_0_0 = prim_buffer_0_siss[1];

    auto g_0_xxxxxz_0_0_0 = prim_buffer_0_siss[2];

    auto g_0_xxxxyy_0_0_0 = prim_buffer_0_siss[3];

    auto g_0_xxxxyz_0_0_0 = prim_buffer_0_siss[4];

    auto g_0_xxxxzz_0_0_0 = prim_buffer_0_siss[5];

    auto g_0_xxxyyy_0_0_0 = prim_buffer_0_siss[6];

    auto g_0_xxxyyz_0_0_0 = prim_buffer_0_siss[7];

    auto g_0_xxxyzz_0_0_0 = prim_buffer_0_siss[8];

    auto g_0_xxxzzz_0_0_0 = prim_buffer_0_siss[9];

    auto g_0_xxyyyy_0_0_0 = prim_buffer_0_siss[10];

    auto g_0_xxyyyz_0_0_0 = prim_buffer_0_siss[11];

    auto g_0_xxyyzz_0_0_0 = prim_buffer_0_siss[12];

    auto g_0_xxyzzz_0_0_0 = prim_buffer_0_siss[13];

    auto g_0_xxzzzz_0_0_0 = prim_buffer_0_siss[14];

    auto g_0_xyyyyy_0_0_0 = prim_buffer_0_siss[15];

    auto g_0_xyyyyz_0_0_0 = prim_buffer_0_siss[16];

    auto g_0_xyyyzz_0_0_0 = prim_buffer_0_siss[17];

    auto g_0_xyyzzz_0_0_0 = prim_buffer_0_siss[18];

    auto g_0_xyzzzz_0_0_0 = prim_buffer_0_siss[19];

    auto g_0_xzzzzz_0_0_0 = prim_buffer_0_siss[20];

    auto g_0_yyyyyy_0_0_0 = prim_buffer_0_siss[21];

    auto g_0_yyyyyz_0_0_0 = prim_buffer_0_siss[22];

    auto g_0_yyyyzz_0_0_0 = prim_buffer_0_siss[23];

    auto g_0_yyyzzz_0_0_0 = prim_buffer_0_siss[24];

    auto g_0_yyzzzz_0_0_0 = prim_buffer_0_siss[25];

    auto g_0_yzzzzz_0_0_0 = prim_buffer_0_siss[26];

    auto g_0_zzzzzz_0_0_0 = prim_buffer_0_siss[27];

    #pragma omp simd aligned(g_0_xxxx_0_0_0, g_0_xxxx_0_0_1, g_0_xxxxx_0_0_0, g_0_xxxxx_0_0_1, g_0_xxxxxx_0_0_0, g_0_xxxxxy_0_0_0, g_0_xxxxxz_0_0_0, g_0_xxxxyy_0_0_0, g_0_xxxxyz_0_0_0, g_0_xxxxz_0_0_0, g_0_xxxxz_0_0_1, g_0_xxxxzz_0_0_0, g_0_xxxyy_0_0_0, g_0_xxxyy_0_0_1, g_0_xxxyyy_0_0_0, g_0_xxxyyz_0_0_0, g_0_xxxyzz_0_0_0, g_0_xxxzz_0_0_0, g_0_xxxzz_0_0_1, g_0_xxxzzz_0_0_0, g_0_xxyy_0_0_0, g_0_xxyy_0_0_1, g_0_xxyyy_0_0_0, g_0_xxyyy_0_0_1, g_0_xxyyyy_0_0_0, g_0_xxyyyz_0_0_0, g_0_xxyyzz_0_0_0, g_0_xxyzzz_0_0_0, g_0_xxzz_0_0_0, g_0_xxzz_0_0_1, g_0_xxzzz_0_0_0, g_0_xxzzz_0_0_1, g_0_xxzzzz_0_0_0, g_0_xyyy_0_0_0, g_0_xyyy_0_0_1, g_0_xyyyy_0_0_0, g_0_xyyyy_0_0_1, g_0_xyyyyy_0_0_0, g_0_xyyyyz_0_0_0, g_0_xyyyzz_0_0_0, g_0_xyyzz_0_0_0, g_0_xyyzz_0_0_1, g_0_xyyzzz_0_0_0, g_0_xyzzzz_0_0_0, g_0_xzzz_0_0_0, g_0_xzzz_0_0_1, g_0_xzzzz_0_0_0, g_0_xzzzz_0_0_1, g_0_xzzzzz_0_0_0, g_0_yyyy_0_0_0, g_0_yyyy_0_0_1, g_0_yyyyy_0_0_0, g_0_yyyyy_0_0_1, g_0_yyyyyy_0_0_0, g_0_yyyyyz_0_0_0, g_0_yyyyz_0_0_0, g_0_yyyyz_0_0_1, g_0_yyyyzz_0_0_0, g_0_yyyzz_0_0_0, g_0_yyyzz_0_0_1, g_0_yyyzzz_0_0_0, g_0_yyzz_0_0_0, g_0_yyzz_0_0_1, g_0_yyzzz_0_0_0, g_0_yyzzz_0_0_1, g_0_yyzzzz_0_0_0, g_0_yzzz_0_0_0, g_0_yzzz_0_0_1, g_0_yzzzz_0_0_0, g_0_yzzzz_0_0_1, g_0_yzzzzz_0_0_0, g_0_zzzz_0_0_0, g_0_zzzz_0_0_1, g_0_zzzzz_0_0_0, g_0_zzzzz_0_0_1, g_0_zzzzzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_0_0[i] = 5.0 * g_0_xxxx_0_0_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_0_1[i] * fti_ab_0 + g_0_xxxxx_0_0_0[i] * pb_x + g_0_xxxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxxy_0_0_0[i] = g_0_xxxxx_0_0_0[i] * pb_y + g_0_xxxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxxz_0_0_0[i] = g_0_xxxxx_0_0_0[i] * pb_z + g_0_xxxxx_0_0_1[i] * wp_z[i];

        g_0_xxxxyy_0_0_0[i] = 3.0 * g_0_xxyy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_0_1[i] * fti_ab_0 + g_0_xxxyy_0_0_0[i] * pb_x + g_0_xxxyy_0_0_1[i] * wp_x[i];

        g_0_xxxxyz_0_0_0[i] = g_0_xxxxz_0_0_0[i] * pb_y + g_0_xxxxz_0_0_1[i] * wp_y[i];

        g_0_xxxxzz_0_0_0[i] = 3.0 * g_0_xxzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_0_1[i] * fti_ab_0 + g_0_xxxzz_0_0_0[i] * pb_x + g_0_xxxzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyy_0_0_0[i] = 2.0 * g_0_xyyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_0_1[i] * fti_ab_0 + g_0_xxyyy_0_0_0[i] * pb_x + g_0_xxyyy_0_0_1[i] * wp_x[i];

        g_0_xxxyyz_0_0_0[i] = g_0_xxxyy_0_0_0[i] * pb_z + g_0_xxxyy_0_0_1[i] * wp_z[i];

        g_0_xxxyzz_0_0_0[i] = g_0_xxxzz_0_0_0[i] * pb_y + g_0_xxxzz_0_0_1[i] * wp_y[i];

        g_0_xxxzzz_0_0_0[i] = 2.0 * g_0_xzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_0_1[i] * fti_ab_0 + g_0_xxzzz_0_0_0[i] * pb_x + g_0_xxzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyy_0_0_0[i] = g_0_yyyy_0_0_0[i] * fi_ab_0 - g_0_yyyy_0_0_1[i] * fti_ab_0 + g_0_xyyyy_0_0_0[i] * pb_x + g_0_xyyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyyz_0_0_0[i] = g_0_xxyyy_0_0_0[i] * pb_z + g_0_xxyyy_0_0_1[i] * wp_z[i];

        g_0_xxyyzz_0_0_0[i] = g_0_yyzz_0_0_0[i] * fi_ab_0 - g_0_yyzz_0_0_1[i] * fti_ab_0 + g_0_xyyzz_0_0_0[i] * pb_x + g_0_xyyzz_0_0_1[i] * wp_x[i];

        g_0_xxyzzz_0_0_0[i] = g_0_xxzzz_0_0_0[i] * pb_y + g_0_xxzzz_0_0_1[i] * wp_y[i];

        g_0_xxzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * fi_ab_0 - g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_xzzzz_0_0_0[i] * pb_x + g_0_xzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyy_0_0_0[i] = g_0_yyyyy_0_0_0[i] * pb_x + g_0_yyyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyyz_0_0_0[i] = g_0_yyyyz_0_0_0[i] * pb_x + g_0_yyyyz_0_0_1[i] * wp_x[i];

        g_0_xyyyzz_0_0_0[i] = g_0_yyyzz_0_0_0[i] * pb_x + g_0_yyyzz_0_0_1[i] * wp_x[i];

        g_0_xyyzzz_0_0_0[i] = g_0_yyzzz_0_0_0[i] * pb_x + g_0_yyzzz_0_0_1[i] * wp_x[i];

        g_0_xyzzzz_0_0_0[i] = g_0_yzzzz_0_0_0[i] * pb_x + g_0_yzzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * pb_x + g_0_zzzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyyy_0_0_0[i] = 5.0 * g_0_yyyy_0_0_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_0_1[i] * fti_ab_0 + g_0_yyyyy_0_0_0[i] * pb_y + g_0_yyyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyyz_0_0_0[i] = g_0_yyyyy_0_0_0[i] * pb_z + g_0_yyyyy_0_0_1[i] * wp_z[i];

        g_0_yyyyzz_0_0_0[i] = 3.0 * g_0_yyzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_0_1[i] * fti_ab_0 + g_0_yyyzz_0_0_0[i] * pb_y + g_0_yyyzz_0_0_1[i] * wp_y[i];

        g_0_yyyzzz_0_0_0[i] = 2.0 * g_0_yzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_0_1[i] * fti_ab_0 + g_0_yyzzz_0_0_0[i] * pb_y + g_0_yyzzz_0_0_1[i] * wp_y[i];

        g_0_yyzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * fi_ab_0 - g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_yzzzz_0_0_0[i] * pb_y + g_0_yzzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * pb_y + g_0_zzzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzzz_0_0_0[i] = 5.0 * g_0_zzzz_0_0_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_zzzzz_0_0_0[i] * pb_z + g_0_zzzzz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace


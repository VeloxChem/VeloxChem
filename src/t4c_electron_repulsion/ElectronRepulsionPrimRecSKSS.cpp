#include "ElectronRepulsionPrimRecSKSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_skss(CSimdArray<double>& prim_buffer_0_skss,
                                  const CSimdArray<double>& prim_buffer_0_shss,
                                  const CSimdArray<double>& prim_buffer_1_shss,
                                  const CSimdArray<double>& prim_buffer_0_siss,
                                  const CSimdArray<double>& prim_buffer_1_siss,
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
    const auto ndims = prim_buffer_0_skss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_shss

    auto g_0_xxxxx_0_0_0 = prim_buffer_0_shss[0];

    auto g_0_xxxyy_0_0_0 = prim_buffer_0_shss[3];

    auto g_0_xxxzz_0_0_0 = prim_buffer_0_shss[5];

    auto g_0_xxyyy_0_0_0 = prim_buffer_0_shss[6];

    auto g_0_xxzzz_0_0_0 = prim_buffer_0_shss[9];

    auto g_0_xyyyy_0_0_0 = prim_buffer_0_shss[10];

    auto g_0_xyyzz_0_0_0 = prim_buffer_0_shss[12];

    auto g_0_xzzzz_0_0_0 = prim_buffer_0_shss[14];

    auto g_0_yyyyy_0_0_0 = prim_buffer_0_shss[15];

    auto g_0_yyyzz_0_0_0 = prim_buffer_0_shss[17];

    auto g_0_yyzzz_0_0_0 = prim_buffer_0_shss[18];

    auto g_0_yzzzz_0_0_0 = prim_buffer_0_shss[19];

    auto g_0_zzzzz_0_0_0 = prim_buffer_0_shss[20];

    /// Set up components of auxilary buffer : prim_buffer_1_shss

    auto g_0_xxxxx_0_0_1 = prim_buffer_1_shss[0];

    auto g_0_xxxyy_0_0_1 = prim_buffer_1_shss[3];

    auto g_0_xxxzz_0_0_1 = prim_buffer_1_shss[5];

    auto g_0_xxyyy_0_0_1 = prim_buffer_1_shss[6];

    auto g_0_xxzzz_0_0_1 = prim_buffer_1_shss[9];

    auto g_0_xyyyy_0_0_1 = prim_buffer_1_shss[10];

    auto g_0_xyyzz_0_0_1 = prim_buffer_1_shss[12];

    auto g_0_xzzzz_0_0_1 = prim_buffer_1_shss[14];

    auto g_0_yyyyy_0_0_1 = prim_buffer_1_shss[15];

    auto g_0_yyyzz_0_0_1 = prim_buffer_1_shss[17];

    auto g_0_yyzzz_0_0_1 = prim_buffer_1_shss[18];

    auto g_0_yzzzz_0_0_1 = prim_buffer_1_shss[19];

    auto g_0_zzzzz_0_0_1 = prim_buffer_1_shss[20];

    /// Set up components of auxilary buffer : prim_buffer_0_siss

    auto g_0_xxxxxx_0_0_0 = prim_buffer_0_siss[0];

    auto g_0_xxxxxz_0_0_0 = prim_buffer_0_siss[2];

    auto g_0_xxxxyy_0_0_0 = prim_buffer_0_siss[3];

    auto g_0_xxxxzz_0_0_0 = prim_buffer_0_siss[5];

    auto g_0_xxxyyy_0_0_0 = prim_buffer_0_siss[6];

    auto g_0_xxxzzz_0_0_0 = prim_buffer_0_siss[9];

    auto g_0_xxyyyy_0_0_0 = prim_buffer_0_siss[10];

    auto g_0_xxyyzz_0_0_0 = prim_buffer_0_siss[12];

    auto g_0_xxzzzz_0_0_0 = prim_buffer_0_siss[14];

    auto g_0_xyyyyy_0_0_0 = prim_buffer_0_siss[15];

    auto g_0_xyyyzz_0_0_0 = prim_buffer_0_siss[17];

    auto g_0_xyyzzz_0_0_0 = prim_buffer_0_siss[18];

    auto g_0_xzzzzz_0_0_0 = prim_buffer_0_siss[20];

    auto g_0_yyyyyy_0_0_0 = prim_buffer_0_siss[21];

    auto g_0_yyyyyz_0_0_0 = prim_buffer_0_siss[22];

    auto g_0_yyyyzz_0_0_0 = prim_buffer_0_siss[23];

    auto g_0_yyyzzz_0_0_0 = prim_buffer_0_siss[24];

    auto g_0_yyzzzz_0_0_0 = prim_buffer_0_siss[25];

    auto g_0_yzzzzz_0_0_0 = prim_buffer_0_siss[26];

    auto g_0_zzzzzz_0_0_0 = prim_buffer_0_siss[27];

    /// Set up components of auxilary buffer : prim_buffer_1_siss

    auto g_0_xxxxxx_0_0_1 = prim_buffer_1_siss[0];

    auto g_0_xxxxxz_0_0_1 = prim_buffer_1_siss[2];

    auto g_0_xxxxyy_0_0_1 = prim_buffer_1_siss[3];

    auto g_0_xxxxzz_0_0_1 = prim_buffer_1_siss[5];

    auto g_0_xxxyyy_0_0_1 = prim_buffer_1_siss[6];

    auto g_0_xxxzzz_0_0_1 = prim_buffer_1_siss[9];

    auto g_0_xxyyyy_0_0_1 = prim_buffer_1_siss[10];

    auto g_0_xxyyzz_0_0_1 = prim_buffer_1_siss[12];

    auto g_0_xxzzzz_0_0_1 = prim_buffer_1_siss[14];

    auto g_0_xyyyyy_0_0_1 = prim_buffer_1_siss[15];

    auto g_0_xyyyzz_0_0_1 = prim_buffer_1_siss[17];

    auto g_0_xyyzzz_0_0_1 = prim_buffer_1_siss[18];

    auto g_0_xzzzzz_0_0_1 = prim_buffer_1_siss[20];

    auto g_0_yyyyyy_0_0_1 = prim_buffer_1_siss[21];

    auto g_0_yyyyyz_0_0_1 = prim_buffer_1_siss[22];

    auto g_0_yyyyzz_0_0_1 = prim_buffer_1_siss[23];

    auto g_0_yyyzzz_0_0_1 = prim_buffer_1_siss[24];

    auto g_0_yyzzzz_0_0_1 = prim_buffer_1_siss[25];

    auto g_0_yzzzzz_0_0_1 = prim_buffer_1_siss[26];

    auto g_0_zzzzzz_0_0_1 = prim_buffer_1_siss[27];

    /// Set up components of targeted buffer : prim_buffer_0_skss

    auto g_0_xxxxxxx_0_0_0 = prim_buffer_0_skss[0];

    auto g_0_xxxxxxy_0_0_0 = prim_buffer_0_skss[1];

    auto g_0_xxxxxxz_0_0_0 = prim_buffer_0_skss[2];

    auto g_0_xxxxxyy_0_0_0 = prim_buffer_0_skss[3];

    auto g_0_xxxxxyz_0_0_0 = prim_buffer_0_skss[4];

    auto g_0_xxxxxzz_0_0_0 = prim_buffer_0_skss[5];

    auto g_0_xxxxyyy_0_0_0 = prim_buffer_0_skss[6];

    auto g_0_xxxxyyz_0_0_0 = prim_buffer_0_skss[7];

    auto g_0_xxxxyzz_0_0_0 = prim_buffer_0_skss[8];

    auto g_0_xxxxzzz_0_0_0 = prim_buffer_0_skss[9];

    auto g_0_xxxyyyy_0_0_0 = prim_buffer_0_skss[10];

    auto g_0_xxxyyyz_0_0_0 = prim_buffer_0_skss[11];

    auto g_0_xxxyyzz_0_0_0 = prim_buffer_0_skss[12];

    auto g_0_xxxyzzz_0_0_0 = prim_buffer_0_skss[13];

    auto g_0_xxxzzzz_0_0_0 = prim_buffer_0_skss[14];

    auto g_0_xxyyyyy_0_0_0 = prim_buffer_0_skss[15];

    auto g_0_xxyyyyz_0_0_0 = prim_buffer_0_skss[16];

    auto g_0_xxyyyzz_0_0_0 = prim_buffer_0_skss[17];

    auto g_0_xxyyzzz_0_0_0 = prim_buffer_0_skss[18];

    auto g_0_xxyzzzz_0_0_0 = prim_buffer_0_skss[19];

    auto g_0_xxzzzzz_0_0_0 = prim_buffer_0_skss[20];

    auto g_0_xyyyyyy_0_0_0 = prim_buffer_0_skss[21];

    auto g_0_xyyyyyz_0_0_0 = prim_buffer_0_skss[22];

    auto g_0_xyyyyzz_0_0_0 = prim_buffer_0_skss[23];

    auto g_0_xyyyzzz_0_0_0 = prim_buffer_0_skss[24];

    auto g_0_xyyzzzz_0_0_0 = prim_buffer_0_skss[25];

    auto g_0_xyzzzzz_0_0_0 = prim_buffer_0_skss[26];

    auto g_0_xzzzzzz_0_0_0 = prim_buffer_0_skss[27];

    auto g_0_yyyyyyy_0_0_0 = prim_buffer_0_skss[28];

    auto g_0_yyyyyyz_0_0_0 = prim_buffer_0_skss[29];

    auto g_0_yyyyyzz_0_0_0 = prim_buffer_0_skss[30];

    auto g_0_yyyyzzz_0_0_0 = prim_buffer_0_skss[31];

    auto g_0_yyyzzzz_0_0_0 = prim_buffer_0_skss[32];

    auto g_0_yyzzzzz_0_0_0 = prim_buffer_0_skss[33];

    auto g_0_yzzzzzz_0_0_0 = prim_buffer_0_skss[34];

    auto g_0_zzzzzzz_0_0_0 = prim_buffer_0_skss[35];

    #pragma omp simd aligned(g_0_xxxxx_0_0_0, g_0_xxxxx_0_0_1, g_0_xxxxxx_0_0_0, g_0_xxxxxx_0_0_1, g_0_xxxxxxx_0_0_0, g_0_xxxxxxy_0_0_0, g_0_xxxxxxz_0_0_0, g_0_xxxxxyy_0_0_0, g_0_xxxxxyz_0_0_0, g_0_xxxxxz_0_0_0, g_0_xxxxxz_0_0_1, g_0_xxxxxzz_0_0_0, g_0_xxxxyy_0_0_0, g_0_xxxxyy_0_0_1, g_0_xxxxyyy_0_0_0, g_0_xxxxyyz_0_0_0, g_0_xxxxyzz_0_0_0, g_0_xxxxzz_0_0_0, g_0_xxxxzz_0_0_1, g_0_xxxxzzz_0_0_0, g_0_xxxyy_0_0_0, g_0_xxxyy_0_0_1, g_0_xxxyyy_0_0_0, g_0_xxxyyy_0_0_1, g_0_xxxyyyy_0_0_0, g_0_xxxyyyz_0_0_0, g_0_xxxyyzz_0_0_0, g_0_xxxyzzz_0_0_0, g_0_xxxzz_0_0_0, g_0_xxxzz_0_0_1, g_0_xxxzzz_0_0_0, g_0_xxxzzz_0_0_1, g_0_xxxzzzz_0_0_0, g_0_xxyyy_0_0_0, g_0_xxyyy_0_0_1, g_0_xxyyyy_0_0_0, g_0_xxyyyy_0_0_1, g_0_xxyyyyy_0_0_0, g_0_xxyyyyz_0_0_0, g_0_xxyyyzz_0_0_0, g_0_xxyyzz_0_0_0, g_0_xxyyzz_0_0_1, g_0_xxyyzzz_0_0_0, g_0_xxyzzzz_0_0_0, g_0_xxzzz_0_0_0, g_0_xxzzz_0_0_1, g_0_xxzzzz_0_0_0, g_0_xxzzzz_0_0_1, g_0_xxzzzzz_0_0_0, g_0_xyyyy_0_0_0, g_0_xyyyy_0_0_1, g_0_xyyyyy_0_0_0, g_0_xyyyyy_0_0_1, g_0_xyyyyyy_0_0_0, g_0_xyyyyyz_0_0_0, g_0_xyyyyzz_0_0_0, g_0_xyyyzz_0_0_0, g_0_xyyyzz_0_0_1, g_0_xyyyzzz_0_0_0, g_0_xyyzz_0_0_0, g_0_xyyzz_0_0_1, g_0_xyyzzz_0_0_0, g_0_xyyzzz_0_0_1, g_0_xyyzzzz_0_0_0, g_0_xyzzzzz_0_0_0, g_0_xzzzz_0_0_0, g_0_xzzzz_0_0_1, g_0_xzzzzz_0_0_0, g_0_xzzzzz_0_0_1, g_0_xzzzzzz_0_0_0, g_0_yyyyy_0_0_0, g_0_yyyyy_0_0_1, g_0_yyyyyy_0_0_0, g_0_yyyyyy_0_0_1, g_0_yyyyyyy_0_0_0, g_0_yyyyyyz_0_0_0, g_0_yyyyyz_0_0_0, g_0_yyyyyz_0_0_1, g_0_yyyyyzz_0_0_0, g_0_yyyyzz_0_0_0, g_0_yyyyzz_0_0_1, g_0_yyyyzzz_0_0_0, g_0_yyyzz_0_0_0, g_0_yyyzz_0_0_1, g_0_yyyzzz_0_0_0, g_0_yyyzzz_0_0_1, g_0_yyyzzzz_0_0_0, g_0_yyzzz_0_0_0, g_0_yyzzz_0_0_1, g_0_yyzzzz_0_0_0, g_0_yyzzzz_0_0_1, g_0_yyzzzzz_0_0_0, g_0_yzzzz_0_0_0, g_0_yzzzz_0_0_1, g_0_yzzzzz_0_0_0, g_0_yzzzzz_0_0_1, g_0_yzzzzzz_0_0_0, g_0_zzzzz_0_0_0, g_0_zzzzz_0_0_1, g_0_zzzzzz_0_0_0, g_0_zzzzzz_0_0_1, g_0_zzzzzzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_0_0[i] = 6.0 * g_0_xxxxx_0_0_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_0_1[i] * fti_ab_0 + g_0_xxxxxx_0_0_0[i] * pb_x + g_0_xxxxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxxxy_0_0_0[i] = g_0_xxxxxx_0_0_0[i] * pb_y + g_0_xxxxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxxxz_0_0_0[i] = g_0_xxxxxx_0_0_0[i] * pb_z + g_0_xxxxxx_0_0_1[i] * wp_z[i];

        g_0_xxxxxyy_0_0_0[i] = 4.0 * g_0_xxxyy_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_0_1[i] * fti_ab_0 + g_0_xxxxyy_0_0_0[i] * pb_x + g_0_xxxxyy_0_0_1[i] * wp_x[i];

        g_0_xxxxxyz_0_0_0[i] = g_0_xxxxxz_0_0_0[i] * pb_y + g_0_xxxxxz_0_0_1[i] * wp_y[i];

        g_0_xxxxxzz_0_0_0[i] = 4.0 * g_0_xxxzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_0_1[i] * fti_ab_0 + g_0_xxxxzz_0_0_0[i] * pb_x + g_0_xxxxzz_0_0_1[i] * wp_x[i];

        g_0_xxxxyyy_0_0_0[i] = 3.0 * g_0_xxyyy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_0_1[i] * fti_ab_0 + g_0_xxxyyy_0_0_0[i] * pb_x + g_0_xxxyyy_0_0_1[i] * wp_x[i];

        g_0_xxxxyyz_0_0_0[i] = g_0_xxxxyy_0_0_0[i] * pb_z + g_0_xxxxyy_0_0_1[i] * wp_z[i];

        g_0_xxxxyzz_0_0_0[i] = g_0_xxxxzz_0_0_0[i] * pb_y + g_0_xxxxzz_0_0_1[i] * wp_y[i];

        g_0_xxxxzzz_0_0_0[i] = 3.0 * g_0_xxzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_0_1[i] * fti_ab_0 + g_0_xxxzzz_0_0_0[i] * pb_x + g_0_xxxzzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyyy_0_0_0[i] = 2.0 * g_0_xyyyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_0_1[i] * fti_ab_0 + g_0_xxyyyy_0_0_0[i] * pb_x + g_0_xxyyyy_0_0_1[i] * wp_x[i];

        g_0_xxxyyyz_0_0_0[i] = g_0_xxxyyy_0_0_0[i] * pb_z + g_0_xxxyyy_0_0_1[i] * wp_z[i];

        g_0_xxxyyzz_0_0_0[i] = 2.0 * g_0_xyyzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_0_1[i] * fti_ab_0 + g_0_xxyyzz_0_0_0[i] * pb_x + g_0_xxyyzz_0_0_1[i] * wp_x[i];

        g_0_xxxyzzz_0_0_0[i] = g_0_xxxzzz_0_0_0[i] * pb_y + g_0_xxxzzz_0_0_1[i] * wp_y[i];

        g_0_xxxzzzz_0_0_0[i] = 2.0 * g_0_xzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_0_1[i] * fti_ab_0 + g_0_xxzzzz_0_0_0[i] * pb_x + g_0_xxzzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyyy_0_0_0[i] = g_0_yyyyy_0_0_0[i] * fi_ab_0 - g_0_yyyyy_0_0_1[i] * fti_ab_0 + g_0_xyyyyy_0_0_0[i] * pb_x + g_0_xyyyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyyyz_0_0_0[i] = g_0_xxyyyy_0_0_0[i] * pb_z + g_0_xxyyyy_0_0_1[i] * wp_z[i];

        g_0_xxyyyzz_0_0_0[i] = g_0_yyyzz_0_0_0[i] * fi_ab_0 - g_0_yyyzz_0_0_1[i] * fti_ab_0 + g_0_xyyyzz_0_0_0[i] * pb_x + g_0_xyyyzz_0_0_1[i] * wp_x[i];

        g_0_xxyyzzz_0_0_0[i] = g_0_yyzzz_0_0_0[i] * fi_ab_0 - g_0_yyzzz_0_0_1[i] * fti_ab_0 + g_0_xyyzzz_0_0_0[i] * pb_x + g_0_xyyzzz_0_0_1[i] * wp_x[i];

        g_0_xxyzzzz_0_0_0[i] = g_0_xxzzzz_0_0_0[i] * pb_y + g_0_xxzzzz_0_0_1[i] * wp_y[i];

        g_0_xxzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_xzzzzz_0_0_0[i] * pb_x + g_0_xzzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyyy_0_0_0[i] = g_0_yyyyyy_0_0_0[i] * pb_x + g_0_yyyyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyyyz_0_0_0[i] = g_0_yyyyyz_0_0_0[i] * pb_x + g_0_yyyyyz_0_0_1[i] * wp_x[i];

        g_0_xyyyyzz_0_0_0[i] = g_0_yyyyzz_0_0_0[i] * pb_x + g_0_yyyyzz_0_0_1[i] * wp_x[i];

        g_0_xyyyzzz_0_0_0[i] = g_0_yyyzzz_0_0_0[i] * pb_x + g_0_yyyzzz_0_0_1[i] * wp_x[i];

        g_0_xyyzzzz_0_0_0[i] = g_0_yyzzzz_0_0_0[i] * pb_x + g_0_yyzzzz_0_0_1[i] * wp_x[i];

        g_0_xyzzzzz_0_0_0[i] = g_0_yzzzzz_0_0_0[i] * pb_x + g_0_yzzzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * pb_x + g_0_zzzzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyyyy_0_0_0[i] = 6.0 * g_0_yyyyy_0_0_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_0_1[i] * fti_ab_0 + g_0_yyyyyy_0_0_0[i] * pb_y + g_0_yyyyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyyyz_0_0_0[i] = g_0_yyyyyy_0_0_0[i] * pb_z + g_0_yyyyyy_0_0_1[i] * wp_z[i];

        g_0_yyyyyzz_0_0_0[i] = 4.0 * g_0_yyyzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_0_1[i] * fti_ab_0 + g_0_yyyyzz_0_0_0[i] * pb_y + g_0_yyyyzz_0_0_1[i] * wp_y[i];

        g_0_yyyyzzz_0_0_0[i] = 3.0 * g_0_yyzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_0_1[i] * fti_ab_0 + g_0_yyyzzz_0_0_0[i] * pb_y + g_0_yyyzzz_0_0_1[i] * wp_y[i];

        g_0_yyyzzzz_0_0_0[i] = 2.0 * g_0_yzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_0_1[i] * fti_ab_0 + g_0_yyzzzz_0_0_0[i] * pb_y + g_0_yyzzzz_0_0_1[i] * wp_y[i];

        g_0_yyzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_yzzzzz_0_0_0[i] * pb_y + g_0_yzzzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * pb_y + g_0_zzzzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzzzz_0_0_0[i] = 6.0 * g_0_zzzzz_0_0_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_zzzzzz_0_0_0[i] * pb_z + g_0_zzzzzz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace


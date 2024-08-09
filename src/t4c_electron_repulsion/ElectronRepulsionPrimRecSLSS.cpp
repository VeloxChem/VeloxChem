#include "ElectronRepulsionPrimRecSLSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slss(CSimdArray<double>& prim_buffer_0_slss,
                                  const CSimdArray<double>& prim_buffer_0_siss,
                                  const CSimdArray<double>& prim_buffer_1_siss,
                                  const CSimdArray<double>& prim_buffer_0_skss,
                                  const CSimdArray<double>& prim_buffer_1_skss,
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
    const auto ndims = prim_buffer_0_slss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_siss

    auto g_0_xxxxxx_0_0_0 = prim_buffer_0_siss[0];

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

    auto g_0_yyyyzz_0_0_0 = prim_buffer_0_siss[23];

    auto g_0_yyyzzz_0_0_0 = prim_buffer_0_siss[24];

    auto g_0_yyzzzz_0_0_0 = prim_buffer_0_siss[25];

    auto g_0_yzzzzz_0_0_0 = prim_buffer_0_siss[26];

    auto g_0_zzzzzz_0_0_0 = prim_buffer_0_siss[27];

    /// Set up components of auxilary buffer : prim_buffer_1_siss

    auto g_0_xxxxxx_0_0_1 = prim_buffer_1_siss[0];

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

    auto g_0_yyyyzz_0_0_1 = prim_buffer_1_siss[23];

    auto g_0_yyyzzz_0_0_1 = prim_buffer_1_siss[24];

    auto g_0_yyzzzz_0_0_1 = prim_buffer_1_siss[25];

    auto g_0_yzzzzz_0_0_1 = prim_buffer_1_siss[26];

    auto g_0_zzzzzz_0_0_1 = prim_buffer_1_siss[27];

    /// Set up components of auxilary buffer : prim_buffer_0_skss

    auto g_0_xxxxxxx_0_0_0 = prim_buffer_0_skss[0];

    auto g_0_xxxxxxz_0_0_0 = prim_buffer_0_skss[2];

    auto g_0_xxxxxyy_0_0_0 = prim_buffer_0_skss[3];

    auto g_0_xxxxxzz_0_0_0 = prim_buffer_0_skss[5];

    auto g_0_xxxxyyy_0_0_0 = prim_buffer_0_skss[6];

    auto g_0_xxxxzzz_0_0_0 = prim_buffer_0_skss[9];

    auto g_0_xxxyyyy_0_0_0 = prim_buffer_0_skss[10];

    auto g_0_xxxyyzz_0_0_0 = prim_buffer_0_skss[12];

    auto g_0_xxxzzzz_0_0_0 = prim_buffer_0_skss[14];

    auto g_0_xxyyyyy_0_0_0 = prim_buffer_0_skss[15];

    auto g_0_xxyyyzz_0_0_0 = prim_buffer_0_skss[17];

    auto g_0_xxyyzzz_0_0_0 = prim_buffer_0_skss[18];

    auto g_0_xxzzzzz_0_0_0 = prim_buffer_0_skss[20];

    auto g_0_xyyyyyy_0_0_0 = prim_buffer_0_skss[21];

    auto g_0_xyyyyzz_0_0_0 = prim_buffer_0_skss[23];

    auto g_0_xyyyzzz_0_0_0 = prim_buffer_0_skss[24];

    auto g_0_xyyzzzz_0_0_0 = prim_buffer_0_skss[25];

    auto g_0_xzzzzzz_0_0_0 = prim_buffer_0_skss[27];

    auto g_0_yyyyyyy_0_0_0 = prim_buffer_0_skss[28];

    auto g_0_yyyyyyz_0_0_0 = prim_buffer_0_skss[29];

    auto g_0_yyyyyzz_0_0_0 = prim_buffer_0_skss[30];

    auto g_0_yyyyzzz_0_0_0 = prim_buffer_0_skss[31];

    auto g_0_yyyzzzz_0_0_0 = prim_buffer_0_skss[32];

    auto g_0_yyzzzzz_0_0_0 = prim_buffer_0_skss[33];

    auto g_0_yzzzzzz_0_0_0 = prim_buffer_0_skss[34];

    auto g_0_zzzzzzz_0_0_0 = prim_buffer_0_skss[35];

    /// Set up components of auxilary buffer : prim_buffer_1_skss

    auto g_0_xxxxxxx_0_0_1 = prim_buffer_1_skss[0];

    auto g_0_xxxxxxz_0_0_1 = prim_buffer_1_skss[2];

    auto g_0_xxxxxyy_0_0_1 = prim_buffer_1_skss[3];

    auto g_0_xxxxxzz_0_0_1 = prim_buffer_1_skss[5];

    auto g_0_xxxxyyy_0_0_1 = prim_buffer_1_skss[6];

    auto g_0_xxxxzzz_0_0_1 = prim_buffer_1_skss[9];

    auto g_0_xxxyyyy_0_0_1 = prim_buffer_1_skss[10];

    auto g_0_xxxyyzz_0_0_1 = prim_buffer_1_skss[12];

    auto g_0_xxxzzzz_0_0_1 = prim_buffer_1_skss[14];

    auto g_0_xxyyyyy_0_0_1 = prim_buffer_1_skss[15];

    auto g_0_xxyyyzz_0_0_1 = prim_buffer_1_skss[17];

    auto g_0_xxyyzzz_0_0_1 = prim_buffer_1_skss[18];

    auto g_0_xxzzzzz_0_0_1 = prim_buffer_1_skss[20];

    auto g_0_xyyyyyy_0_0_1 = prim_buffer_1_skss[21];

    auto g_0_xyyyyzz_0_0_1 = prim_buffer_1_skss[23];

    auto g_0_xyyyzzz_0_0_1 = prim_buffer_1_skss[24];

    auto g_0_xyyzzzz_0_0_1 = prim_buffer_1_skss[25];

    auto g_0_xzzzzzz_0_0_1 = prim_buffer_1_skss[27];

    auto g_0_yyyyyyy_0_0_1 = prim_buffer_1_skss[28];

    auto g_0_yyyyyyz_0_0_1 = prim_buffer_1_skss[29];

    auto g_0_yyyyyzz_0_0_1 = prim_buffer_1_skss[30];

    auto g_0_yyyyzzz_0_0_1 = prim_buffer_1_skss[31];

    auto g_0_yyyzzzz_0_0_1 = prim_buffer_1_skss[32];

    auto g_0_yyzzzzz_0_0_1 = prim_buffer_1_skss[33];

    auto g_0_yzzzzzz_0_0_1 = prim_buffer_1_skss[34];

    auto g_0_zzzzzzz_0_0_1 = prim_buffer_1_skss[35];

    /// Set up components of targeted buffer : prim_buffer_0_slss

    auto g_0_xxxxxxxx_0_0_0 = prim_buffer_0_slss[0];

    auto g_0_xxxxxxxy_0_0_0 = prim_buffer_0_slss[1];

    auto g_0_xxxxxxxz_0_0_0 = prim_buffer_0_slss[2];

    auto g_0_xxxxxxyy_0_0_0 = prim_buffer_0_slss[3];

    auto g_0_xxxxxxyz_0_0_0 = prim_buffer_0_slss[4];

    auto g_0_xxxxxxzz_0_0_0 = prim_buffer_0_slss[5];

    auto g_0_xxxxxyyy_0_0_0 = prim_buffer_0_slss[6];

    auto g_0_xxxxxyyz_0_0_0 = prim_buffer_0_slss[7];

    auto g_0_xxxxxyzz_0_0_0 = prim_buffer_0_slss[8];

    auto g_0_xxxxxzzz_0_0_0 = prim_buffer_0_slss[9];

    auto g_0_xxxxyyyy_0_0_0 = prim_buffer_0_slss[10];

    auto g_0_xxxxyyyz_0_0_0 = prim_buffer_0_slss[11];

    auto g_0_xxxxyyzz_0_0_0 = prim_buffer_0_slss[12];

    auto g_0_xxxxyzzz_0_0_0 = prim_buffer_0_slss[13];

    auto g_0_xxxxzzzz_0_0_0 = prim_buffer_0_slss[14];

    auto g_0_xxxyyyyy_0_0_0 = prim_buffer_0_slss[15];

    auto g_0_xxxyyyyz_0_0_0 = prim_buffer_0_slss[16];

    auto g_0_xxxyyyzz_0_0_0 = prim_buffer_0_slss[17];

    auto g_0_xxxyyzzz_0_0_0 = prim_buffer_0_slss[18];

    auto g_0_xxxyzzzz_0_0_0 = prim_buffer_0_slss[19];

    auto g_0_xxxzzzzz_0_0_0 = prim_buffer_0_slss[20];

    auto g_0_xxyyyyyy_0_0_0 = prim_buffer_0_slss[21];

    auto g_0_xxyyyyyz_0_0_0 = prim_buffer_0_slss[22];

    auto g_0_xxyyyyzz_0_0_0 = prim_buffer_0_slss[23];

    auto g_0_xxyyyzzz_0_0_0 = prim_buffer_0_slss[24];

    auto g_0_xxyyzzzz_0_0_0 = prim_buffer_0_slss[25];

    auto g_0_xxyzzzzz_0_0_0 = prim_buffer_0_slss[26];

    auto g_0_xxzzzzzz_0_0_0 = prim_buffer_0_slss[27];

    auto g_0_xyyyyyyy_0_0_0 = prim_buffer_0_slss[28];

    auto g_0_xyyyyyyz_0_0_0 = prim_buffer_0_slss[29];

    auto g_0_xyyyyyzz_0_0_0 = prim_buffer_0_slss[30];

    auto g_0_xyyyyzzz_0_0_0 = prim_buffer_0_slss[31];

    auto g_0_xyyyzzzz_0_0_0 = prim_buffer_0_slss[32];

    auto g_0_xyyzzzzz_0_0_0 = prim_buffer_0_slss[33];

    auto g_0_xyzzzzzz_0_0_0 = prim_buffer_0_slss[34];

    auto g_0_xzzzzzzz_0_0_0 = prim_buffer_0_slss[35];

    auto g_0_yyyyyyyy_0_0_0 = prim_buffer_0_slss[36];

    auto g_0_yyyyyyyz_0_0_0 = prim_buffer_0_slss[37];

    auto g_0_yyyyyyzz_0_0_0 = prim_buffer_0_slss[38];

    auto g_0_yyyyyzzz_0_0_0 = prim_buffer_0_slss[39];

    auto g_0_yyyyzzzz_0_0_0 = prim_buffer_0_slss[40];

    auto g_0_yyyzzzzz_0_0_0 = prim_buffer_0_slss[41];

    auto g_0_yyzzzzzz_0_0_0 = prim_buffer_0_slss[42];

    auto g_0_yzzzzzzz_0_0_0 = prim_buffer_0_slss[43];

    auto g_0_zzzzzzzz_0_0_0 = prim_buffer_0_slss[44];

    #pragma omp simd aligned(g_0_xxxxxx_0_0_0, g_0_xxxxxx_0_0_1, g_0_xxxxxxx_0_0_0, g_0_xxxxxxx_0_0_1, g_0_xxxxxxxx_0_0_0, g_0_xxxxxxxy_0_0_0, g_0_xxxxxxxz_0_0_0, g_0_xxxxxxyy_0_0_0, g_0_xxxxxxyz_0_0_0, g_0_xxxxxxz_0_0_0, g_0_xxxxxxz_0_0_1, g_0_xxxxxxzz_0_0_0, g_0_xxxxxyy_0_0_0, g_0_xxxxxyy_0_0_1, g_0_xxxxxyyy_0_0_0, g_0_xxxxxyyz_0_0_0, g_0_xxxxxyzz_0_0_0, g_0_xxxxxzz_0_0_0, g_0_xxxxxzz_0_0_1, g_0_xxxxxzzz_0_0_0, g_0_xxxxyy_0_0_0, g_0_xxxxyy_0_0_1, g_0_xxxxyyy_0_0_0, g_0_xxxxyyy_0_0_1, g_0_xxxxyyyy_0_0_0, g_0_xxxxyyyz_0_0_0, g_0_xxxxyyzz_0_0_0, g_0_xxxxyzzz_0_0_0, g_0_xxxxzz_0_0_0, g_0_xxxxzz_0_0_1, g_0_xxxxzzz_0_0_0, g_0_xxxxzzz_0_0_1, g_0_xxxxzzzz_0_0_0, g_0_xxxyyy_0_0_0, g_0_xxxyyy_0_0_1, g_0_xxxyyyy_0_0_0, g_0_xxxyyyy_0_0_1, g_0_xxxyyyyy_0_0_0, g_0_xxxyyyyz_0_0_0, g_0_xxxyyyzz_0_0_0, g_0_xxxyyzz_0_0_0, g_0_xxxyyzz_0_0_1, g_0_xxxyyzzz_0_0_0, g_0_xxxyzzzz_0_0_0, g_0_xxxzzz_0_0_0, g_0_xxxzzz_0_0_1, g_0_xxxzzzz_0_0_0, g_0_xxxzzzz_0_0_1, g_0_xxxzzzzz_0_0_0, g_0_xxyyyy_0_0_0, g_0_xxyyyy_0_0_1, g_0_xxyyyyy_0_0_0, g_0_xxyyyyy_0_0_1, g_0_xxyyyyyy_0_0_0, g_0_xxyyyyyz_0_0_0, g_0_xxyyyyzz_0_0_0, g_0_xxyyyzz_0_0_0, g_0_xxyyyzz_0_0_1, g_0_xxyyyzzz_0_0_0, g_0_xxyyzz_0_0_0, g_0_xxyyzz_0_0_1, g_0_xxyyzzz_0_0_0, g_0_xxyyzzz_0_0_1, g_0_xxyyzzzz_0_0_0, g_0_xxyzzzzz_0_0_0, g_0_xxzzzz_0_0_0, g_0_xxzzzz_0_0_1, g_0_xxzzzzz_0_0_0, g_0_xxzzzzz_0_0_1, g_0_xxzzzzzz_0_0_0, g_0_xyyyyy_0_0_0, g_0_xyyyyy_0_0_1, g_0_xyyyyyy_0_0_0, g_0_xyyyyyy_0_0_1, g_0_xyyyyyyy_0_0_0, g_0_xyyyyyyz_0_0_0, g_0_xyyyyyzz_0_0_0, g_0_xyyyyzz_0_0_0, g_0_xyyyyzz_0_0_1, g_0_xyyyyzzz_0_0_0, g_0_xyyyzz_0_0_0, g_0_xyyyzz_0_0_1, g_0_xyyyzzz_0_0_0, g_0_xyyyzzz_0_0_1, g_0_xyyyzzzz_0_0_0, g_0_xyyzzz_0_0_0, g_0_xyyzzz_0_0_1, g_0_xyyzzzz_0_0_0, g_0_xyyzzzz_0_0_1, g_0_xyyzzzzz_0_0_0, g_0_xyzzzzzz_0_0_0, g_0_xzzzzz_0_0_0, g_0_xzzzzz_0_0_1, g_0_xzzzzzz_0_0_0, g_0_xzzzzzz_0_0_1, g_0_xzzzzzzz_0_0_0, g_0_yyyyyy_0_0_0, g_0_yyyyyy_0_0_1, g_0_yyyyyyy_0_0_0, g_0_yyyyyyy_0_0_1, g_0_yyyyyyyy_0_0_0, g_0_yyyyyyyz_0_0_0, g_0_yyyyyyz_0_0_0, g_0_yyyyyyz_0_0_1, g_0_yyyyyyzz_0_0_0, g_0_yyyyyzz_0_0_0, g_0_yyyyyzz_0_0_1, g_0_yyyyyzzz_0_0_0, g_0_yyyyzz_0_0_0, g_0_yyyyzz_0_0_1, g_0_yyyyzzz_0_0_0, g_0_yyyyzzz_0_0_1, g_0_yyyyzzzz_0_0_0, g_0_yyyzzz_0_0_0, g_0_yyyzzz_0_0_1, g_0_yyyzzzz_0_0_0, g_0_yyyzzzz_0_0_1, g_0_yyyzzzzz_0_0_0, g_0_yyzzzz_0_0_0, g_0_yyzzzz_0_0_1, g_0_yyzzzzz_0_0_0, g_0_yyzzzzz_0_0_1, g_0_yyzzzzzz_0_0_0, g_0_yzzzzz_0_0_0, g_0_yzzzzz_0_0_1, g_0_yzzzzzz_0_0_0, g_0_yzzzzzz_0_0_1, g_0_yzzzzzzz_0_0_0, g_0_zzzzzz_0_0_0, g_0_zzzzzz_0_0_1, g_0_zzzzzzz_0_0_0, g_0_zzzzzzz_0_0_1, g_0_zzzzzzzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_0_0[i] = 7.0 * g_0_xxxxxx_0_0_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_0_1[i] * fti_ab_0 + g_0_xxxxxxx_0_0_0[i] * pb_x + g_0_xxxxxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxxxxy_0_0_0[i] = g_0_xxxxxxx_0_0_0[i] * pb_y + g_0_xxxxxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxxxxz_0_0_0[i] = g_0_xxxxxxx_0_0_0[i] * pb_z + g_0_xxxxxxx_0_0_1[i] * wp_z[i];

        g_0_xxxxxxyy_0_0_0[i] = 5.0 * g_0_xxxxyy_0_0_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_0_1[i] * fti_ab_0 + g_0_xxxxxyy_0_0_0[i] * pb_x + g_0_xxxxxyy_0_0_1[i] * wp_x[i];

        g_0_xxxxxxyz_0_0_0[i] = g_0_xxxxxxz_0_0_0[i] * pb_y + g_0_xxxxxxz_0_0_1[i] * wp_y[i];

        g_0_xxxxxxzz_0_0_0[i] = 5.0 * g_0_xxxxzz_0_0_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_0_1[i] * fti_ab_0 + g_0_xxxxxzz_0_0_0[i] * pb_x + g_0_xxxxxzz_0_0_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_0_0[i] = 4.0 * g_0_xxxyyy_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_0_1[i] * fti_ab_0 + g_0_xxxxyyy_0_0_0[i] * pb_x + g_0_xxxxyyy_0_0_1[i] * wp_x[i];

        g_0_xxxxxyyz_0_0_0[i] = g_0_xxxxxyy_0_0_0[i] * pb_z + g_0_xxxxxyy_0_0_1[i] * wp_z[i];

        g_0_xxxxxyzz_0_0_0[i] = g_0_xxxxxzz_0_0_0[i] * pb_y + g_0_xxxxxzz_0_0_1[i] * wp_y[i];

        g_0_xxxxxzzz_0_0_0[i] = 4.0 * g_0_xxxzzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_0_1[i] * fti_ab_0 + g_0_xxxxzzz_0_0_0[i] * pb_x + g_0_xxxxzzz_0_0_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_0_0[i] = 3.0 * g_0_xxyyyy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_0_1[i] * fti_ab_0 + g_0_xxxyyyy_0_0_0[i] * pb_x + g_0_xxxyyyy_0_0_1[i] * wp_x[i];

        g_0_xxxxyyyz_0_0_0[i] = g_0_xxxxyyy_0_0_0[i] * pb_z + g_0_xxxxyyy_0_0_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_0_0[i] = 3.0 * g_0_xxyyzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_0_1[i] * fti_ab_0 + g_0_xxxyyzz_0_0_0[i] * pb_x + g_0_xxxyyzz_0_0_1[i] * wp_x[i];

        g_0_xxxxyzzz_0_0_0[i] = g_0_xxxxzzz_0_0_0[i] * pb_y + g_0_xxxxzzz_0_0_1[i] * wp_y[i];

        g_0_xxxxzzzz_0_0_0[i] = 3.0 * g_0_xxzzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_0_1[i] * fti_ab_0 + g_0_xxxzzzz_0_0_0[i] * pb_x + g_0_xxxzzzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_0_0[i] = 2.0 * g_0_xyyyyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_0_1[i] * fti_ab_0 + g_0_xxyyyyy_0_0_0[i] * pb_x + g_0_xxyyyyy_0_0_1[i] * wp_x[i];

        g_0_xxxyyyyz_0_0_0[i] = g_0_xxxyyyy_0_0_0[i] * pb_z + g_0_xxxyyyy_0_0_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_0_0[i] = 2.0 * g_0_xyyyzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_0_1[i] * fti_ab_0 + g_0_xxyyyzz_0_0_0[i] * pb_x + g_0_xxyyyzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_0_0[i] = 2.0 * g_0_xyyzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_0_1[i] * fti_ab_0 + g_0_xxyyzzz_0_0_0[i] * pb_x + g_0_xxyyzzz_0_0_1[i] * wp_x[i];

        g_0_xxxyzzzz_0_0_0[i] = g_0_xxxzzzz_0_0_0[i] * pb_y + g_0_xxxzzzz_0_0_1[i] * wp_y[i];

        g_0_xxxzzzzz_0_0_0[i] = 2.0 * g_0_xzzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_0_1[i] * fti_ab_0 + g_0_xxzzzzz_0_0_0[i] * pb_x + g_0_xxzzzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_0_0[i] = g_0_yyyyyy_0_0_0[i] * fi_ab_0 - g_0_yyyyyy_0_0_1[i] * fti_ab_0 + g_0_xyyyyyy_0_0_0[i] * pb_x + g_0_xyyyyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyyyyz_0_0_0[i] = g_0_xxyyyyy_0_0_0[i] * pb_z + g_0_xxyyyyy_0_0_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_0_0[i] = g_0_yyyyzz_0_0_0[i] * fi_ab_0 - g_0_yyyyzz_0_0_1[i] * fti_ab_0 + g_0_xyyyyzz_0_0_0[i] * pb_x + g_0_xyyyyzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_0_0[i] = g_0_yyyzzz_0_0_0[i] * fi_ab_0 - g_0_yyyzzz_0_0_1[i] * fti_ab_0 + g_0_xyyyzzz_0_0_0[i] * pb_x + g_0_xyyyzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_0_0[i] = g_0_yyzzzz_0_0_0[i] * fi_ab_0 - g_0_yyzzzz_0_0_1[i] * fti_ab_0 + g_0_xyyzzzz_0_0_0[i] * pb_x + g_0_xyyzzzz_0_0_1[i] * wp_x[i];

        g_0_xxyzzzzz_0_0_0[i] = g_0_xxzzzzz_0_0_0[i] * pb_y + g_0_xxzzzzz_0_0_1[i] * wp_y[i];

        g_0_xxzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzzz_0_0_1[i] * fti_ab_0 + g_0_xzzzzzz_0_0_0[i] * pb_x + g_0_xzzzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_0_0[i] = g_0_yyyyyyy_0_0_0[i] * pb_x + g_0_yyyyyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_0_0[i] = g_0_yyyyyyz_0_0_0[i] * pb_x + g_0_yyyyyyz_0_0_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_0_0[i] = g_0_yyyyyzz_0_0_0[i] * pb_x + g_0_yyyyyzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_0_0[i] = g_0_yyyyzzz_0_0_0[i] * pb_x + g_0_yyyyzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_0_0[i] = g_0_yyyzzzz_0_0_0[i] * pb_x + g_0_yyyzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_0_0[i] = g_0_yyzzzzz_0_0_0[i] * pb_x + g_0_yyzzzzz_0_0_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_0_0[i] = g_0_yzzzzzz_0_0_0[i] * pb_x + g_0_yzzzzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_0_0[i] = g_0_zzzzzzz_0_0_0[i] * pb_x + g_0_zzzzzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyyyyy_0_0_0[i] = 7.0 * g_0_yyyyyy_0_0_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_0_1[i] * fti_ab_0 + g_0_yyyyyyy_0_0_0[i] * pb_y + g_0_yyyyyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyyyyz_0_0_0[i] = g_0_yyyyyyy_0_0_0[i] * pb_z + g_0_yyyyyyy_0_0_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_0_0[i] = 5.0 * g_0_yyyyzz_0_0_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_0_1[i] * fti_ab_0 + g_0_yyyyyzz_0_0_0[i] * pb_y + g_0_yyyyyzz_0_0_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_0_0[i] = 4.0 * g_0_yyyzzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_0_1[i] * fti_ab_0 + g_0_yyyyzzz_0_0_0[i] * pb_y + g_0_yyyyzzz_0_0_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_0_0[i] = 3.0 * g_0_yyzzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_0_1[i] * fti_ab_0 + g_0_yyyzzzz_0_0_0[i] * pb_y + g_0_yyyzzzz_0_0_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_0_0[i] = 2.0 * g_0_yzzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_0_1[i] * fti_ab_0 + g_0_yyzzzzz_0_0_0[i] * pb_y + g_0_yyzzzzz_0_0_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzzz_0_0_1[i] * fti_ab_0 + g_0_yzzzzzz_0_0_0[i] * pb_y + g_0_yzzzzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_0_0[i] = g_0_zzzzzzz_0_0_0[i] * pb_y + g_0_zzzzzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzzzzz_0_0_0[i] = 7.0 * g_0_zzzzzz_0_0_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_0_1[i] * fti_ab_0 + g_0_zzzzzzz_0_0_0[i] * pb_z + g_0_zzzzzzz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace


#include "ElectronRepulsionPrimRecSPSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsk(CSimdArray<double>& prim_buffer_0_spsk,
                                  const CSimdArray<double>& prim_buffer_1_sssi,
                                  const CSimdArray<double>& prim_buffer_0_sssk,
                                  const CSimdArray<double>& prim_buffer_1_sssk,
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
    const auto ndims = prim_buffer_0_spsk.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_1_sssi

    auto g_0_0_0_xxxxxx_1 = prim_buffer_1_sssi[0];

    auto g_0_0_0_xxxxxy_1 = prim_buffer_1_sssi[1];

    auto g_0_0_0_xxxxxz_1 = prim_buffer_1_sssi[2];

    auto g_0_0_0_xxxxyy_1 = prim_buffer_1_sssi[3];

    auto g_0_0_0_xxxxyz_1 = prim_buffer_1_sssi[4];

    auto g_0_0_0_xxxxzz_1 = prim_buffer_1_sssi[5];

    auto g_0_0_0_xxxyyy_1 = prim_buffer_1_sssi[6];

    auto g_0_0_0_xxxyyz_1 = prim_buffer_1_sssi[7];

    auto g_0_0_0_xxxyzz_1 = prim_buffer_1_sssi[8];

    auto g_0_0_0_xxxzzz_1 = prim_buffer_1_sssi[9];

    auto g_0_0_0_xxyyyy_1 = prim_buffer_1_sssi[10];

    auto g_0_0_0_xxyyyz_1 = prim_buffer_1_sssi[11];

    auto g_0_0_0_xxyyzz_1 = prim_buffer_1_sssi[12];

    auto g_0_0_0_xxyzzz_1 = prim_buffer_1_sssi[13];

    auto g_0_0_0_xxzzzz_1 = prim_buffer_1_sssi[14];

    auto g_0_0_0_xyyyyy_1 = prim_buffer_1_sssi[15];

    auto g_0_0_0_xyyyyz_1 = prim_buffer_1_sssi[16];

    auto g_0_0_0_xyyyzz_1 = prim_buffer_1_sssi[17];

    auto g_0_0_0_xyyzzz_1 = prim_buffer_1_sssi[18];

    auto g_0_0_0_xyzzzz_1 = prim_buffer_1_sssi[19];

    auto g_0_0_0_xzzzzz_1 = prim_buffer_1_sssi[20];

    auto g_0_0_0_yyyyyy_1 = prim_buffer_1_sssi[21];

    auto g_0_0_0_yyyyyz_1 = prim_buffer_1_sssi[22];

    auto g_0_0_0_yyyyzz_1 = prim_buffer_1_sssi[23];

    auto g_0_0_0_yyyzzz_1 = prim_buffer_1_sssi[24];

    auto g_0_0_0_yyzzzz_1 = prim_buffer_1_sssi[25];

    auto g_0_0_0_yzzzzz_1 = prim_buffer_1_sssi[26];

    auto g_0_0_0_zzzzzz_1 = prim_buffer_1_sssi[27];

    /// Set up components of auxilary buffer : prim_buffer_0_sssk

    auto g_0_0_0_xxxxxxx_0 = prim_buffer_0_sssk[0];

    auto g_0_0_0_xxxxxxy_0 = prim_buffer_0_sssk[1];

    auto g_0_0_0_xxxxxxz_0 = prim_buffer_0_sssk[2];

    auto g_0_0_0_xxxxxyy_0 = prim_buffer_0_sssk[3];

    auto g_0_0_0_xxxxxyz_0 = prim_buffer_0_sssk[4];

    auto g_0_0_0_xxxxxzz_0 = prim_buffer_0_sssk[5];

    auto g_0_0_0_xxxxyyy_0 = prim_buffer_0_sssk[6];

    auto g_0_0_0_xxxxyyz_0 = prim_buffer_0_sssk[7];

    auto g_0_0_0_xxxxyzz_0 = prim_buffer_0_sssk[8];

    auto g_0_0_0_xxxxzzz_0 = prim_buffer_0_sssk[9];

    auto g_0_0_0_xxxyyyy_0 = prim_buffer_0_sssk[10];

    auto g_0_0_0_xxxyyyz_0 = prim_buffer_0_sssk[11];

    auto g_0_0_0_xxxyyzz_0 = prim_buffer_0_sssk[12];

    auto g_0_0_0_xxxyzzz_0 = prim_buffer_0_sssk[13];

    auto g_0_0_0_xxxzzzz_0 = prim_buffer_0_sssk[14];

    auto g_0_0_0_xxyyyyy_0 = prim_buffer_0_sssk[15];

    auto g_0_0_0_xxyyyyz_0 = prim_buffer_0_sssk[16];

    auto g_0_0_0_xxyyyzz_0 = prim_buffer_0_sssk[17];

    auto g_0_0_0_xxyyzzz_0 = prim_buffer_0_sssk[18];

    auto g_0_0_0_xxyzzzz_0 = prim_buffer_0_sssk[19];

    auto g_0_0_0_xxzzzzz_0 = prim_buffer_0_sssk[20];

    auto g_0_0_0_xyyyyyy_0 = prim_buffer_0_sssk[21];

    auto g_0_0_0_xyyyyyz_0 = prim_buffer_0_sssk[22];

    auto g_0_0_0_xyyyyzz_0 = prim_buffer_0_sssk[23];

    auto g_0_0_0_xyyyzzz_0 = prim_buffer_0_sssk[24];

    auto g_0_0_0_xyyzzzz_0 = prim_buffer_0_sssk[25];

    auto g_0_0_0_xyzzzzz_0 = prim_buffer_0_sssk[26];

    auto g_0_0_0_xzzzzzz_0 = prim_buffer_0_sssk[27];

    auto g_0_0_0_yyyyyyy_0 = prim_buffer_0_sssk[28];

    auto g_0_0_0_yyyyyyz_0 = prim_buffer_0_sssk[29];

    auto g_0_0_0_yyyyyzz_0 = prim_buffer_0_sssk[30];

    auto g_0_0_0_yyyyzzz_0 = prim_buffer_0_sssk[31];

    auto g_0_0_0_yyyzzzz_0 = prim_buffer_0_sssk[32];

    auto g_0_0_0_yyzzzzz_0 = prim_buffer_0_sssk[33];

    auto g_0_0_0_yzzzzzz_0 = prim_buffer_0_sssk[34];

    auto g_0_0_0_zzzzzzz_0 = prim_buffer_0_sssk[35];

    /// Set up components of auxilary buffer : prim_buffer_1_sssk

    auto g_0_0_0_xxxxxxx_1 = prim_buffer_1_sssk[0];

    auto g_0_0_0_xxxxxxy_1 = prim_buffer_1_sssk[1];

    auto g_0_0_0_xxxxxxz_1 = prim_buffer_1_sssk[2];

    auto g_0_0_0_xxxxxyy_1 = prim_buffer_1_sssk[3];

    auto g_0_0_0_xxxxxyz_1 = prim_buffer_1_sssk[4];

    auto g_0_0_0_xxxxxzz_1 = prim_buffer_1_sssk[5];

    auto g_0_0_0_xxxxyyy_1 = prim_buffer_1_sssk[6];

    auto g_0_0_0_xxxxyyz_1 = prim_buffer_1_sssk[7];

    auto g_0_0_0_xxxxyzz_1 = prim_buffer_1_sssk[8];

    auto g_0_0_0_xxxxzzz_1 = prim_buffer_1_sssk[9];

    auto g_0_0_0_xxxyyyy_1 = prim_buffer_1_sssk[10];

    auto g_0_0_0_xxxyyyz_1 = prim_buffer_1_sssk[11];

    auto g_0_0_0_xxxyyzz_1 = prim_buffer_1_sssk[12];

    auto g_0_0_0_xxxyzzz_1 = prim_buffer_1_sssk[13];

    auto g_0_0_0_xxxzzzz_1 = prim_buffer_1_sssk[14];

    auto g_0_0_0_xxyyyyy_1 = prim_buffer_1_sssk[15];

    auto g_0_0_0_xxyyyyz_1 = prim_buffer_1_sssk[16];

    auto g_0_0_0_xxyyyzz_1 = prim_buffer_1_sssk[17];

    auto g_0_0_0_xxyyzzz_1 = prim_buffer_1_sssk[18];

    auto g_0_0_0_xxyzzzz_1 = prim_buffer_1_sssk[19];

    auto g_0_0_0_xxzzzzz_1 = prim_buffer_1_sssk[20];

    auto g_0_0_0_xyyyyyy_1 = prim_buffer_1_sssk[21];

    auto g_0_0_0_xyyyyyz_1 = prim_buffer_1_sssk[22];

    auto g_0_0_0_xyyyyzz_1 = prim_buffer_1_sssk[23];

    auto g_0_0_0_xyyyzzz_1 = prim_buffer_1_sssk[24];

    auto g_0_0_0_xyyzzzz_1 = prim_buffer_1_sssk[25];

    auto g_0_0_0_xyzzzzz_1 = prim_buffer_1_sssk[26];

    auto g_0_0_0_xzzzzzz_1 = prim_buffer_1_sssk[27];

    auto g_0_0_0_yyyyyyy_1 = prim_buffer_1_sssk[28];

    auto g_0_0_0_yyyyyyz_1 = prim_buffer_1_sssk[29];

    auto g_0_0_0_yyyyyzz_1 = prim_buffer_1_sssk[30];

    auto g_0_0_0_yyyyzzz_1 = prim_buffer_1_sssk[31];

    auto g_0_0_0_yyyzzzz_1 = prim_buffer_1_sssk[32];

    auto g_0_0_0_yyzzzzz_1 = prim_buffer_1_sssk[33];

    auto g_0_0_0_yzzzzzz_1 = prim_buffer_1_sssk[34];

    auto g_0_0_0_zzzzzzz_1 = prim_buffer_1_sssk[35];

    /// Set up 0-36 components of targeted buffer : prim_buffer_0_spsk

    auto g_0_x_0_xxxxxxx_0 = prim_buffer_0_spsk[0];

    auto g_0_x_0_xxxxxxy_0 = prim_buffer_0_spsk[1];

    auto g_0_x_0_xxxxxxz_0 = prim_buffer_0_spsk[2];

    auto g_0_x_0_xxxxxyy_0 = prim_buffer_0_spsk[3];

    auto g_0_x_0_xxxxxyz_0 = prim_buffer_0_spsk[4];

    auto g_0_x_0_xxxxxzz_0 = prim_buffer_0_spsk[5];

    auto g_0_x_0_xxxxyyy_0 = prim_buffer_0_spsk[6];

    auto g_0_x_0_xxxxyyz_0 = prim_buffer_0_spsk[7];

    auto g_0_x_0_xxxxyzz_0 = prim_buffer_0_spsk[8];

    auto g_0_x_0_xxxxzzz_0 = prim_buffer_0_spsk[9];

    auto g_0_x_0_xxxyyyy_0 = prim_buffer_0_spsk[10];

    auto g_0_x_0_xxxyyyz_0 = prim_buffer_0_spsk[11];

    auto g_0_x_0_xxxyyzz_0 = prim_buffer_0_spsk[12];

    auto g_0_x_0_xxxyzzz_0 = prim_buffer_0_spsk[13];

    auto g_0_x_0_xxxzzzz_0 = prim_buffer_0_spsk[14];

    auto g_0_x_0_xxyyyyy_0 = prim_buffer_0_spsk[15];

    auto g_0_x_0_xxyyyyz_0 = prim_buffer_0_spsk[16];

    auto g_0_x_0_xxyyyzz_0 = prim_buffer_0_spsk[17];

    auto g_0_x_0_xxyyzzz_0 = prim_buffer_0_spsk[18];

    auto g_0_x_0_xxyzzzz_0 = prim_buffer_0_spsk[19];

    auto g_0_x_0_xxzzzzz_0 = prim_buffer_0_spsk[20];

    auto g_0_x_0_xyyyyyy_0 = prim_buffer_0_spsk[21];

    auto g_0_x_0_xyyyyyz_0 = prim_buffer_0_spsk[22];

    auto g_0_x_0_xyyyyzz_0 = prim_buffer_0_spsk[23];

    auto g_0_x_0_xyyyzzz_0 = prim_buffer_0_spsk[24];

    auto g_0_x_0_xyyzzzz_0 = prim_buffer_0_spsk[25];

    auto g_0_x_0_xyzzzzz_0 = prim_buffer_0_spsk[26];

    auto g_0_x_0_xzzzzzz_0 = prim_buffer_0_spsk[27];

    auto g_0_x_0_yyyyyyy_0 = prim_buffer_0_spsk[28];

    auto g_0_x_0_yyyyyyz_0 = prim_buffer_0_spsk[29];

    auto g_0_x_0_yyyyyzz_0 = prim_buffer_0_spsk[30];

    auto g_0_x_0_yyyyzzz_0 = prim_buffer_0_spsk[31];

    auto g_0_x_0_yyyzzzz_0 = prim_buffer_0_spsk[32];

    auto g_0_x_0_yyzzzzz_0 = prim_buffer_0_spsk[33];

    auto g_0_x_0_yzzzzzz_0 = prim_buffer_0_spsk[34];

    auto g_0_x_0_zzzzzzz_0 = prim_buffer_0_spsk[35];

    #pragma omp simd aligned(g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_x_0_xxxxxxx_0, g_0_x_0_xxxxxxy_0, g_0_x_0_xxxxxxz_0, g_0_x_0_xxxxxyy_0, g_0_x_0_xxxxxyz_0, g_0_x_0_xxxxxzz_0, g_0_x_0_xxxxyyy_0, g_0_x_0_xxxxyyz_0, g_0_x_0_xxxxyzz_0, g_0_x_0_xxxxzzz_0, g_0_x_0_xxxyyyy_0, g_0_x_0_xxxyyyz_0, g_0_x_0_xxxyyzz_0, g_0_x_0_xxxyzzz_0, g_0_x_0_xxxzzzz_0, g_0_x_0_xxyyyyy_0, g_0_x_0_xxyyyyz_0, g_0_x_0_xxyyyzz_0, g_0_x_0_xxyyzzz_0, g_0_x_0_xxyzzzz_0, g_0_x_0_xxzzzzz_0, g_0_x_0_xyyyyyy_0, g_0_x_0_xyyyyyz_0, g_0_x_0_xyyyyzz_0, g_0_x_0_xyyyzzz_0, g_0_x_0_xyyzzzz_0, g_0_x_0_xyzzzzz_0, g_0_x_0_xzzzzzz_0, g_0_x_0_yyyyyyy_0, g_0_x_0_yyyyyyz_0, g_0_x_0_yyyyyzz_0, g_0_x_0_yyyyzzz_0, g_0_x_0_yyyzzzz_0, g_0_x_0_yyzzzzz_0, g_0_x_0_yzzzzzz_0, g_0_x_0_zzzzzzz_0, wp_x  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxxxxxx_0[i] = 7.0 * g_0_0_0_xxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxx_0[i] * pb_x + g_0_0_0_xxxxxxx_1[i] * wp_x[i];

        g_0_x_0_xxxxxxy_0[i] = 6.0 * g_0_0_0_xxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxy_0[i] * pb_x + g_0_0_0_xxxxxxy_1[i] * wp_x[i];

        g_0_x_0_xxxxxxz_0[i] = 6.0 * g_0_0_0_xxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxz_0[i] * pb_x + g_0_0_0_xxxxxxz_1[i] * wp_x[i];

        g_0_x_0_xxxxxyy_0[i] = 5.0 * g_0_0_0_xxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyy_0[i] * pb_x + g_0_0_0_xxxxxyy_1[i] * wp_x[i];

        g_0_x_0_xxxxxyz_0[i] = 5.0 * g_0_0_0_xxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyz_0[i] * pb_x + g_0_0_0_xxxxxyz_1[i] * wp_x[i];

        g_0_x_0_xxxxxzz_0[i] = 5.0 * g_0_0_0_xxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxzz_0[i] * pb_x + g_0_0_0_xxxxxzz_1[i] * wp_x[i];

        g_0_x_0_xxxxyyy_0[i] = 4.0 * g_0_0_0_xxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyy_0[i] * pb_x + g_0_0_0_xxxxyyy_1[i] * wp_x[i];

        g_0_x_0_xxxxyyz_0[i] = 4.0 * g_0_0_0_xxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyz_0[i] * pb_x + g_0_0_0_xxxxyyz_1[i] * wp_x[i];

        g_0_x_0_xxxxyzz_0[i] = 4.0 * g_0_0_0_xxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzz_0[i] * pb_x + g_0_0_0_xxxxyzz_1[i] * wp_x[i];

        g_0_x_0_xxxxzzz_0[i] = 4.0 * g_0_0_0_xxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzzz_0[i] * pb_x + g_0_0_0_xxxxzzz_1[i] * wp_x[i];

        g_0_x_0_xxxyyyy_0[i] = 3.0 * g_0_0_0_xxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyy_0[i] * pb_x + g_0_0_0_xxxyyyy_1[i] * wp_x[i];

        g_0_x_0_xxxyyyz_0[i] = 3.0 * g_0_0_0_xxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyz_0[i] * pb_x + g_0_0_0_xxxyyyz_1[i] * wp_x[i];

        g_0_x_0_xxxyyzz_0[i] = 3.0 * g_0_0_0_xxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzz_0[i] * pb_x + g_0_0_0_xxxyyzz_1[i] * wp_x[i];

        g_0_x_0_xxxyzzz_0[i] = 3.0 * g_0_0_0_xxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzz_0[i] * pb_x + g_0_0_0_xxxyzzz_1[i] * wp_x[i];

        g_0_x_0_xxxzzzz_0[i] = 3.0 * g_0_0_0_xxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzzz_0[i] * pb_x + g_0_0_0_xxxzzzz_1[i] * wp_x[i];

        g_0_x_0_xxyyyyy_0[i] = 2.0 * g_0_0_0_xyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyy_0[i] * pb_x + g_0_0_0_xxyyyyy_1[i] * wp_x[i];

        g_0_x_0_xxyyyyz_0[i] = 2.0 * g_0_0_0_xyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyz_0[i] * pb_x + g_0_0_0_xxyyyyz_1[i] * wp_x[i];

        g_0_x_0_xxyyyzz_0[i] = 2.0 * g_0_0_0_xyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzz_0[i] * pb_x + g_0_0_0_xxyyyzz_1[i] * wp_x[i];

        g_0_x_0_xxyyzzz_0[i] = 2.0 * g_0_0_0_xyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzz_0[i] * pb_x + g_0_0_0_xxyyzzz_1[i] * wp_x[i];

        g_0_x_0_xxyzzzz_0[i] = 2.0 * g_0_0_0_xyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzz_0[i] * pb_x + g_0_0_0_xxyzzzz_1[i] * wp_x[i];

        g_0_x_0_xxzzzzz_0[i] = 2.0 * g_0_0_0_xzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzzz_0[i] * pb_x + g_0_0_0_xxzzzzz_1[i] * wp_x[i];

        g_0_x_0_xyyyyyy_0[i] = g_0_0_0_yyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyy_0[i] * pb_x + g_0_0_0_xyyyyyy_1[i] * wp_x[i];

        g_0_x_0_xyyyyyz_0[i] = g_0_0_0_yyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyz_0[i] * pb_x + g_0_0_0_xyyyyyz_1[i] * wp_x[i];

        g_0_x_0_xyyyyzz_0[i] = g_0_0_0_yyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzz_0[i] * pb_x + g_0_0_0_xyyyyzz_1[i] * wp_x[i];

        g_0_x_0_xyyyzzz_0[i] = g_0_0_0_yyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzz_0[i] * pb_x + g_0_0_0_xyyyzzz_1[i] * wp_x[i];

        g_0_x_0_xyyzzzz_0[i] = g_0_0_0_yyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzz_0[i] * pb_x + g_0_0_0_xyyzzzz_1[i] * wp_x[i];

        g_0_x_0_xyzzzzz_0[i] = g_0_0_0_yzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzz_0[i] * pb_x + g_0_0_0_xyzzzzz_1[i] * wp_x[i];

        g_0_x_0_xzzzzzz_0[i] = g_0_0_0_zzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzzz_0[i] * pb_x + g_0_0_0_xzzzzzz_1[i] * wp_x[i];

        g_0_x_0_yyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * pb_x + g_0_0_0_yyyyyyy_1[i] * wp_x[i];

        g_0_x_0_yyyyyyz_0[i] = g_0_0_0_yyyyyyz_0[i] * pb_x + g_0_0_0_yyyyyyz_1[i] * wp_x[i];

        g_0_x_0_yyyyyzz_0[i] = g_0_0_0_yyyyyzz_0[i] * pb_x + g_0_0_0_yyyyyzz_1[i] * wp_x[i];

        g_0_x_0_yyyyzzz_0[i] = g_0_0_0_yyyyzzz_0[i] * pb_x + g_0_0_0_yyyyzzz_1[i] * wp_x[i];

        g_0_x_0_yyyzzzz_0[i] = g_0_0_0_yyyzzzz_0[i] * pb_x + g_0_0_0_yyyzzzz_1[i] * wp_x[i];

        g_0_x_0_yyzzzzz_0[i] = g_0_0_0_yyzzzzz_0[i] * pb_x + g_0_0_0_yyzzzzz_1[i] * wp_x[i];

        g_0_x_0_yzzzzzz_0[i] = g_0_0_0_yzzzzzz_0[i] * pb_x + g_0_0_0_yzzzzzz_1[i] * wp_x[i];

        g_0_x_0_zzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * pb_x + g_0_0_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : prim_buffer_0_spsk

    auto g_0_y_0_xxxxxxx_0 = prim_buffer_0_spsk[36];

    auto g_0_y_0_xxxxxxy_0 = prim_buffer_0_spsk[37];

    auto g_0_y_0_xxxxxxz_0 = prim_buffer_0_spsk[38];

    auto g_0_y_0_xxxxxyy_0 = prim_buffer_0_spsk[39];

    auto g_0_y_0_xxxxxyz_0 = prim_buffer_0_spsk[40];

    auto g_0_y_0_xxxxxzz_0 = prim_buffer_0_spsk[41];

    auto g_0_y_0_xxxxyyy_0 = prim_buffer_0_spsk[42];

    auto g_0_y_0_xxxxyyz_0 = prim_buffer_0_spsk[43];

    auto g_0_y_0_xxxxyzz_0 = prim_buffer_0_spsk[44];

    auto g_0_y_0_xxxxzzz_0 = prim_buffer_0_spsk[45];

    auto g_0_y_0_xxxyyyy_0 = prim_buffer_0_spsk[46];

    auto g_0_y_0_xxxyyyz_0 = prim_buffer_0_spsk[47];

    auto g_0_y_0_xxxyyzz_0 = prim_buffer_0_spsk[48];

    auto g_0_y_0_xxxyzzz_0 = prim_buffer_0_spsk[49];

    auto g_0_y_0_xxxzzzz_0 = prim_buffer_0_spsk[50];

    auto g_0_y_0_xxyyyyy_0 = prim_buffer_0_spsk[51];

    auto g_0_y_0_xxyyyyz_0 = prim_buffer_0_spsk[52];

    auto g_0_y_0_xxyyyzz_0 = prim_buffer_0_spsk[53];

    auto g_0_y_0_xxyyzzz_0 = prim_buffer_0_spsk[54];

    auto g_0_y_0_xxyzzzz_0 = prim_buffer_0_spsk[55];

    auto g_0_y_0_xxzzzzz_0 = prim_buffer_0_spsk[56];

    auto g_0_y_0_xyyyyyy_0 = prim_buffer_0_spsk[57];

    auto g_0_y_0_xyyyyyz_0 = prim_buffer_0_spsk[58];

    auto g_0_y_0_xyyyyzz_0 = prim_buffer_0_spsk[59];

    auto g_0_y_0_xyyyzzz_0 = prim_buffer_0_spsk[60];

    auto g_0_y_0_xyyzzzz_0 = prim_buffer_0_spsk[61];

    auto g_0_y_0_xyzzzzz_0 = prim_buffer_0_spsk[62];

    auto g_0_y_0_xzzzzzz_0 = prim_buffer_0_spsk[63];

    auto g_0_y_0_yyyyyyy_0 = prim_buffer_0_spsk[64];

    auto g_0_y_0_yyyyyyz_0 = prim_buffer_0_spsk[65];

    auto g_0_y_0_yyyyyzz_0 = prim_buffer_0_spsk[66];

    auto g_0_y_0_yyyyzzz_0 = prim_buffer_0_spsk[67];

    auto g_0_y_0_yyyzzzz_0 = prim_buffer_0_spsk[68];

    auto g_0_y_0_yyzzzzz_0 = prim_buffer_0_spsk[69];

    auto g_0_y_0_yzzzzzz_0 = prim_buffer_0_spsk[70];

    auto g_0_y_0_zzzzzzz_0 = prim_buffer_0_spsk[71];

    #pragma omp simd aligned(g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_y_0_xxxxxxx_0, g_0_y_0_xxxxxxy_0, g_0_y_0_xxxxxxz_0, g_0_y_0_xxxxxyy_0, g_0_y_0_xxxxxyz_0, g_0_y_0_xxxxxzz_0, g_0_y_0_xxxxyyy_0, g_0_y_0_xxxxyyz_0, g_0_y_0_xxxxyzz_0, g_0_y_0_xxxxzzz_0, g_0_y_0_xxxyyyy_0, g_0_y_0_xxxyyyz_0, g_0_y_0_xxxyyzz_0, g_0_y_0_xxxyzzz_0, g_0_y_0_xxxzzzz_0, g_0_y_0_xxyyyyy_0, g_0_y_0_xxyyyyz_0, g_0_y_0_xxyyyzz_0, g_0_y_0_xxyyzzz_0, g_0_y_0_xxyzzzz_0, g_0_y_0_xxzzzzz_0, g_0_y_0_xyyyyyy_0, g_0_y_0_xyyyyyz_0, g_0_y_0_xyyyyzz_0, g_0_y_0_xyyyzzz_0, g_0_y_0_xyyzzzz_0, g_0_y_0_xyzzzzz_0, g_0_y_0_xzzzzzz_0, g_0_y_0_yyyyyyy_0, g_0_y_0_yyyyyyz_0, g_0_y_0_yyyyyzz_0, g_0_y_0_yyyyzzz_0, g_0_y_0_yyyzzzz_0, g_0_y_0_yyzzzzz_0, g_0_y_0_yzzzzzz_0, g_0_y_0_zzzzzzz_0, wp_y  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxxxxxx_0[i] = g_0_0_0_xxxxxxx_0[i] * pb_y + g_0_0_0_xxxxxxx_1[i] * wp_y[i];

        g_0_y_0_xxxxxxy_0[i] = g_0_0_0_xxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxy_0[i] * pb_y + g_0_0_0_xxxxxxy_1[i] * wp_y[i];

        g_0_y_0_xxxxxxz_0[i] = g_0_0_0_xxxxxxz_0[i] * pb_y + g_0_0_0_xxxxxxz_1[i] * wp_y[i];

        g_0_y_0_xxxxxyy_0[i] = 2.0 * g_0_0_0_xxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyy_0[i] * pb_y + g_0_0_0_xxxxxyy_1[i] * wp_y[i];

        g_0_y_0_xxxxxyz_0[i] = g_0_0_0_xxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyz_0[i] * pb_y + g_0_0_0_xxxxxyz_1[i] * wp_y[i];

        g_0_y_0_xxxxxzz_0[i] = g_0_0_0_xxxxxzz_0[i] * pb_y + g_0_0_0_xxxxxzz_1[i] * wp_y[i];

        g_0_y_0_xxxxyyy_0[i] = 3.0 * g_0_0_0_xxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyy_0[i] * pb_y + g_0_0_0_xxxxyyy_1[i] * wp_y[i];

        g_0_y_0_xxxxyyz_0[i] = 2.0 * g_0_0_0_xxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyz_0[i] * pb_y + g_0_0_0_xxxxyyz_1[i] * wp_y[i];

        g_0_y_0_xxxxyzz_0[i] = g_0_0_0_xxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzz_0[i] * pb_y + g_0_0_0_xxxxyzz_1[i] * wp_y[i];

        g_0_y_0_xxxxzzz_0[i] = g_0_0_0_xxxxzzz_0[i] * pb_y + g_0_0_0_xxxxzzz_1[i] * wp_y[i];

        g_0_y_0_xxxyyyy_0[i] = 4.0 * g_0_0_0_xxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyy_0[i] * pb_y + g_0_0_0_xxxyyyy_1[i] * wp_y[i];

        g_0_y_0_xxxyyyz_0[i] = 3.0 * g_0_0_0_xxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyz_0[i] * pb_y + g_0_0_0_xxxyyyz_1[i] * wp_y[i];

        g_0_y_0_xxxyyzz_0[i] = 2.0 * g_0_0_0_xxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzz_0[i] * pb_y + g_0_0_0_xxxyyzz_1[i] * wp_y[i];

        g_0_y_0_xxxyzzz_0[i] = g_0_0_0_xxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzz_0[i] * pb_y + g_0_0_0_xxxyzzz_1[i] * wp_y[i];

        g_0_y_0_xxxzzzz_0[i] = g_0_0_0_xxxzzzz_0[i] * pb_y + g_0_0_0_xxxzzzz_1[i] * wp_y[i];

        g_0_y_0_xxyyyyy_0[i] = 5.0 * g_0_0_0_xxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyy_0[i] * pb_y + g_0_0_0_xxyyyyy_1[i] * wp_y[i];

        g_0_y_0_xxyyyyz_0[i] = 4.0 * g_0_0_0_xxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyz_0[i] * pb_y + g_0_0_0_xxyyyyz_1[i] * wp_y[i];

        g_0_y_0_xxyyyzz_0[i] = 3.0 * g_0_0_0_xxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzz_0[i] * pb_y + g_0_0_0_xxyyyzz_1[i] * wp_y[i];

        g_0_y_0_xxyyzzz_0[i] = 2.0 * g_0_0_0_xxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzz_0[i] * pb_y + g_0_0_0_xxyyzzz_1[i] * wp_y[i];

        g_0_y_0_xxyzzzz_0[i] = g_0_0_0_xxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzz_0[i] * pb_y + g_0_0_0_xxyzzzz_1[i] * wp_y[i];

        g_0_y_0_xxzzzzz_0[i] = g_0_0_0_xxzzzzz_0[i] * pb_y + g_0_0_0_xxzzzzz_1[i] * wp_y[i];

        g_0_y_0_xyyyyyy_0[i] = 6.0 * g_0_0_0_xyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyy_0[i] * pb_y + g_0_0_0_xyyyyyy_1[i] * wp_y[i];

        g_0_y_0_xyyyyyz_0[i] = 5.0 * g_0_0_0_xyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyz_0[i] * pb_y + g_0_0_0_xyyyyyz_1[i] * wp_y[i];

        g_0_y_0_xyyyyzz_0[i] = 4.0 * g_0_0_0_xyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzz_0[i] * pb_y + g_0_0_0_xyyyyzz_1[i] * wp_y[i];

        g_0_y_0_xyyyzzz_0[i] = 3.0 * g_0_0_0_xyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzz_0[i] * pb_y + g_0_0_0_xyyyzzz_1[i] * wp_y[i];

        g_0_y_0_xyyzzzz_0[i] = 2.0 * g_0_0_0_xyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzz_0[i] * pb_y + g_0_0_0_xyyzzzz_1[i] * wp_y[i];

        g_0_y_0_xyzzzzz_0[i] = g_0_0_0_xzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzz_0[i] * pb_y + g_0_0_0_xyzzzzz_1[i] * wp_y[i];

        g_0_y_0_xzzzzzz_0[i] = g_0_0_0_xzzzzzz_0[i] * pb_y + g_0_0_0_xzzzzzz_1[i] * wp_y[i];

        g_0_y_0_yyyyyyy_0[i] = 7.0 * g_0_0_0_yyyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyy_0[i] * pb_y + g_0_0_0_yyyyyyy_1[i] * wp_y[i];

        g_0_y_0_yyyyyyz_0[i] = 6.0 * g_0_0_0_yyyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyz_0[i] * pb_y + g_0_0_0_yyyyyyz_1[i] * wp_y[i];

        g_0_y_0_yyyyyzz_0[i] = 5.0 * g_0_0_0_yyyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyzz_0[i] * pb_y + g_0_0_0_yyyyyzz_1[i] * wp_y[i];

        g_0_y_0_yyyyzzz_0[i] = 4.0 * g_0_0_0_yyyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzzz_0[i] * pb_y + g_0_0_0_yyyyzzz_1[i] * wp_y[i];

        g_0_y_0_yyyzzzz_0[i] = 3.0 * g_0_0_0_yyzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzzz_0[i] * pb_y + g_0_0_0_yyyzzzz_1[i] * wp_y[i];

        g_0_y_0_yyzzzzz_0[i] = 2.0 * g_0_0_0_yzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzzz_0[i] * pb_y + g_0_0_0_yyzzzzz_1[i] * wp_y[i];

        g_0_y_0_yzzzzzz_0[i] = g_0_0_0_zzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzzz_0[i] * pb_y + g_0_0_0_yzzzzzz_1[i] * wp_y[i];

        g_0_y_0_zzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * pb_y + g_0_0_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : prim_buffer_0_spsk

    auto g_0_z_0_xxxxxxx_0 = prim_buffer_0_spsk[72];

    auto g_0_z_0_xxxxxxy_0 = prim_buffer_0_spsk[73];

    auto g_0_z_0_xxxxxxz_0 = prim_buffer_0_spsk[74];

    auto g_0_z_0_xxxxxyy_0 = prim_buffer_0_spsk[75];

    auto g_0_z_0_xxxxxyz_0 = prim_buffer_0_spsk[76];

    auto g_0_z_0_xxxxxzz_0 = prim_buffer_0_spsk[77];

    auto g_0_z_0_xxxxyyy_0 = prim_buffer_0_spsk[78];

    auto g_0_z_0_xxxxyyz_0 = prim_buffer_0_spsk[79];

    auto g_0_z_0_xxxxyzz_0 = prim_buffer_0_spsk[80];

    auto g_0_z_0_xxxxzzz_0 = prim_buffer_0_spsk[81];

    auto g_0_z_0_xxxyyyy_0 = prim_buffer_0_spsk[82];

    auto g_0_z_0_xxxyyyz_0 = prim_buffer_0_spsk[83];

    auto g_0_z_0_xxxyyzz_0 = prim_buffer_0_spsk[84];

    auto g_0_z_0_xxxyzzz_0 = prim_buffer_0_spsk[85];

    auto g_0_z_0_xxxzzzz_0 = prim_buffer_0_spsk[86];

    auto g_0_z_0_xxyyyyy_0 = prim_buffer_0_spsk[87];

    auto g_0_z_0_xxyyyyz_0 = prim_buffer_0_spsk[88];

    auto g_0_z_0_xxyyyzz_0 = prim_buffer_0_spsk[89];

    auto g_0_z_0_xxyyzzz_0 = prim_buffer_0_spsk[90];

    auto g_0_z_0_xxyzzzz_0 = prim_buffer_0_spsk[91];

    auto g_0_z_0_xxzzzzz_0 = prim_buffer_0_spsk[92];

    auto g_0_z_0_xyyyyyy_0 = prim_buffer_0_spsk[93];

    auto g_0_z_0_xyyyyyz_0 = prim_buffer_0_spsk[94];

    auto g_0_z_0_xyyyyzz_0 = prim_buffer_0_spsk[95];

    auto g_0_z_0_xyyyzzz_0 = prim_buffer_0_spsk[96];

    auto g_0_z_0_xyyzzzz_0 = prim_buffer_0_spsk[97];

    auto g_0_z_0_xyzzzzz_0 = prim_buffer_0_spsk[98];

    auto g_0_z_0_xzzzzzz_0 = prim_buffer_0_spsk[99];

    auto g_0_z_0_yyyyyyy_0 = prim_buffer_0_spsk[100];

    auto g_0_z_0_yyyyyyz_0 = prim_buffer_0_spsk[101];

    auto g_0_z_0_yyyyyzz_0 = prim_buffer_0_spsk[102];

    auto g_0_z_0_yyyyzzz_0 = prim_buffer_0_spsk[103];

    auto g_0_z_0_yyyzzzz_0 = prim_buffer_0_spsk[104];

    auto g_0_z_0_yyzzzzz_0 = prim_buffer_0_spsk[105];

    auto g_0_z_0_yzzzzzz_0 = prim_buffer_0_spsk[106];

    auto g_0_z_0_zzzzzzz_0 = prim_buffer_0_spsk[107];

    #pragma omp simd aligned(g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_z_0_xxxxxxx_0, g_0_z_0_xxxxxxy_0, g_0_z_0_xxxxxxz_0, g_0_z_0_xxxxxyy_0, g_0_z_0_xxxxxyz_0, g_0_z_0_xxxxxzz_0, g_0_z_0_xxxxyyy_0, g_0_z_0_xxxxyyz_0, g_0_z_0_xxxxyzz_0, g_0_z_0_xxxxzzz_0, g_0_z_0_xxxyyyy_0, g_0_z_0_xxxyyyz_0, g_0_z_0_xxxyyzz_0, g_0_z_0_xxxyzzz_0, g_0_z_0_xxxzzzz_0, g_0_z_0_xxyyyyy_0, g_0_z_0_xxyyyyz_0, g_0_z_0_xxyyyzz_0, g_0_z_0_xxyyzzz_0, g_0_z_0_xxyzzzz_0, g_0_z_0_xxzzzzz_0, g_0_z_0_xyyyyyy_0, g_0_z_0_xyyyyyz_0, g_0_z_0_xyyyyzz_0, g_0_z_0_xyyyzzz_0, g_0_z_0_xyyzzzz_0, g_0_z_0_xyzzzzz_0, g_0_z_0_xzzzzzz_0, g_0_z_0_yyyyyyy_0, g_0_z_0_yyyyyyz_0, g_0_z_0_yyyyyzz_0, g_0_z_0_yyyyzzz_0, g_0_z_0_yyyzzzz_0, g_0_z_0_yyzzzzz_0, g_0_z_0_yzzzzzz_0, g_0_z_0_zzzzzzz_0, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxxxxxx_0[i] = g_0_0_0_xxxxxxx_0[i] * pb_z + g_0_0_0_xxxxxxx_1[i] * wp_z[i];

        g_0_z_0_xxxxxxy_0[i] = g_0_0_0_xxxxxxy_0[i] * pb_z + g_0_0_0_xxxxxxy_1[i] * wp_z[i];

        g_0_z_0_xxxxxxz_0[i] = g_0_0_0_xxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxz_0[i] * pb_z + g_0_0_0_xxxxxxz_1[i] * wp_z[i];

        g_0_z_0_xxxxxyy_0[i] = g_0_0_0_xxxxxyy_0[i] * pb_z + g_0_0_0_xxxxxyy_1[i] * wp_z[i];

        g_0_z_0_xxxxxyz_0[i] = g_0_0_0_xxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyz_0[i] * pb_z + g_0_0_0_xxxxxyz_1[i] * wp_z[i];

        g_0_z_0_xxxxxzz_0[i] = 2.0 * g_0_0_0_xxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxzz_0[i] * pb_z + g_0_0_0_xxxxxzz_1[i] * wp_z[i];

        g_0_z_0_xxxxyyy_0[i] = g_0_0_0_xxxxyyy_0[i] * pb_z + g_0_0_0_xxxxyyy_1[i] * wp_z[i];

        g_0_z_0_xxxxyyz_0[i] = g_0_0_0_xxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyz_0[i] * pb_z + g_0_0_0_xxxxyyz_1[i] * wp_z[i];

        g_0_z_0_xxxxyzz_0[i] = 2.0 * g_0_0_0_xxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzz_0[i] * pb_z + g_0_0_0_xxxxyzz_1[i] * wp_z[i];

        g_0_z_0_xxxxzzz_0[i] = 3.0 * g_0_0_0_xxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzzz_0[i] * pb_z + g_0_0_0_xxxxzzz_1[i] * wp_z[i];

        g_0_z_0_xxxyyyy_0[i] = g_0_0_0_xxxyyyy_0[i] * pb_z + g_0_0_0_xxxyyyy_1[i] * wp_z[i];

        g_0_z_0_xxxyyyz_0[i] = g_0_0_0_xxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyz_0[i] * pb_z + g_0_0_0_xxxyyyz_1[i] * wp_z[i];

        g_0_z_0_xxxyyzz_0[i] = 2.0 * g_0_0_0_xxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzz_0[i] * pb_z + g_0_0_0_xxxyyzz_1[i] * wp_z[i];

        g_0_z_0_xxxyzzz_0[i] = 3.0 * g_0_0_0_xxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzz_0[i] * pb_z + g_0_0_0_xxxyzzz_1[i] * wp_z[i];

        g_0_z_0_xxxzzzz_0[i] = 4.0 * g_0_0_0_xxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzzz_0[i] * pb_z + g_0_0_0_xxxzzzz_1[i] * wp_z[i];

        g_0_z_0_xxyyyyy_0[i] = g_0_0_0_xxyyyyy_0[i] * pb_z + g_0_0_0_xxyyyyy_1[i] * wp_z[i];

        g_0_z_0_xxyyyyz_0[i] = g_0_0_0_xxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyz_0[i] * pb_z + g_0_0_0_xxyyyyz_1[i] * wp_z[i];

        g_0_z_0_xxyyyzz_0[i] = 2.0 * g_0_0_0_xxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzz_0[i] * pb_z + g_0_0_0_xxyyyzz_1[i] * wp_z[i];

        g_0_z_0_xxyyzzz_0[i] = 3.0 * g_0_0_0_xxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzz_0[i] * pb_z + g_0_0_0_xxyyzzz_1[i] * wp_z[i];

        g_0_z_0_xxyzzzz_0[i] = 4.0 * g_0_0_0_xxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzz_0[i] * pb_z + g_0_0_0_xxyzzzz_1[i] * wp_z[i];

        g_0_z_0_xxzzzzz_0[i] = 5.0 * g_0_0_0_xxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzzz_0[i] * pb_z + g_0_0_0_xxzzzzz_1[i] * wp_z[i];

        g_0_z_0_xyyyyyy_0[i] = g_0_0_0_xyyyyyy_0[i] * pb_z + g_0_0_0_xyyyyyy_1[i] * wp_z[i];

        g_0_z_0_xyyyyyz_0[i] = g_0_0_0_xyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyz_0[i] * pb_z + g_0_0_0_xyyyyyz_1[i] * wp_z[i];

        g_0_z_0_xyyyyzz_0[i] = 2.0 * g_0_0_0_xyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzz_0[i] * pb_z + g_0_0_0_xyyyyzz_1[i] * wp_z[i];

        g_0_z_0_xyyyzzz_0[i] = 3.0 * g_0_0_0_xyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzz_0[i] * pb_z + g_0_0_0_xyyyzzz_1[i] * wp_z[i];

        g_0_z_0_xyyzzzz_0[i] = 4.0 * g_0_0_0_xyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzz_0[i] * pb_z + g_0_0_0_xyyzzzz_1[i] * wp_z[i];

        g_0_z_0_xyzzzzz_0[i] = 5.0 * g_0_0_0_xyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzz_0[i] * pb_z + g_0_0_0_xyzzzzz_1[i] * wp_z[i];

        g_0_z_0_xzzzzzz_0[i] = 6.0 * g_0_0_0_xzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzzz_0[i] * pb_z + g_0_0_0_xzzzzzz_1[i] * wp_z[i];

        g_0_z_0_yyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * pb_z + g_0_0_0_yyyyyyy_1[i] * wp_z[i];

        g_0_z_0_yyyyyyz_0[i] = g_0_0_0_yyyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyz_0[i] * pb_z + g_0_0_0_yyyyyyz_1[i] * wp_z[i];

        g_0_z_0_yyyyyzz_0[i] = 2.0 * g_0_0_0_yyyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyzz_0[i] * pb_z + g_0_0_0_yyyyyzz_1[i] * wp_z[i];

        g_0_z_0_yyyyzzz_0[i] = 3.0 * g_0_0_0_yyyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzzz_0[i] * pb_z + g_0_0_0_yyyyzzz_1[i] * wp_z[i];

        g_0_z_0_yyyzzzz_0[i] = 4.0 * g_0_0_0_yyyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzzz_0[i] * pb_z + g_0_0_0_yyyzzzz_1[i] * wp_z[i];

        g_0_z_0_yyzzzzz_0[i] = 5.0 * g_0_0_0_yyzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzzz_0[i] * pb_z + g_0_0_0_yyzzzzz_1[i] * wp_z[i];

        g_0_z_0_yzzzzzz_0[i] = 6.0 * g_0_0_0_yzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzzz_0[i] * pb_z + g_0_0_0_yzzzzzz_1[i] * wp_z[i];

        g_0_z_0_zzzzzzz_0[i] = 7.0 * g_0_0_0_zzzzzz_1[i] * fi_abcd_0 + g_0_0_0_zzzzzzz_0[i] * pb_z + g_0_0_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


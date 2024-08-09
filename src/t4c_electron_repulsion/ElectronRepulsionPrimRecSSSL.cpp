#include "ElectronRepulsionPrimRecSSSL.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssl(CSimdArray<double>& prim_buffer_0_sssl,
                                  const CSimdArray<double>& prim_buffer_0_sssi,
                                  const CSimdArray<double>& prim_buffer_1_sssi,
                                  const CSimdArray<double>& prim_buffer_0_sssk,
                                  const CSimdArray<double>& prim_buffer_1_sssk,
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
    const auto ndims = prim_buffer_0_sssl.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssi

    auto g_0_0_0_xxxxxx_0 = prim_buffer_0_sssi[0];

    auto g_0_0_0_xxxxyy_0 = prim_buffer_0_sssi[3];

    auto g_0_0_0_xxxxzz_0 = prim_buffer_0_sssi[5];

    auto g_0_0_0_xxxyyy_0 = prim_buffer_0_sssi[6];

    auto g_0_0_0_xxxzzz_0 = prim_buffer_0_sssi[9];

    auto g_0_0_0_xxyyyy_0 = prim_buffer_0_sssi[10];

    auto g_0_0_0_xxyyzz_0 = prim_buffer_0_sssi[12];

    auto g_0_0_0_xxzzzz_0 = prim_buffer_0_sssi[14];

    auto g_0_0_0_xyyyyy_0 = prim_buffer_0_sssi[15];

    auto g_0_0_0_xyyyzz_0 = prim_buffer_0_sssi[17];

    auto g_0_0_0_xyyzzz_0 = prim_buffer_0_sssi[18];

    auto g_0_0_0_xzzzzz_0 = prim_buffer_0_sssi[20];

    auto g_0_0_0_yyyyyy_0 = prim_buffer_0_sssi[21];

    auto g_0_0_0_yyyyzz_0 = prim_buffer_0_sssi[23];

    auto g_0_0_0_yyyzzz_0 = prim_buffer_0_sssi[24];

    auto g_0_0_0_yyzzzz_0 = prim_buffer_0_sssi[25];

    auto g_0_0_0_yzzzzz_0 = prim_buffer_0_sssi[26];

    auto g_0_0_0_zzzzzz_0 = prim_buffer_0_sssi[27];

    /// Set up components of auxilary buffer : prim_buffer_1_sssi

    auto g_0_0_0_xxxxxx_1 = prim_buffer_1_sssi[0];

    auto g_0_0_0_xxxxyy_1 = prim_buffer_1_sssi[3];

    auto g_0_0_0_xxxxzz_1 = prim_buffer_1_sssi[5];

    auto g_0_0_0_xxxyyy_1 = prim_buffer_1_sssi[6];

    auto g_0_0_0_xxxzzz_1 = prim_buffer_1_sssi[9];

    auto g_0_0_0_xxyyyy_1 = prim_buffer_1_sssi[10];

    auto g_0_0_0_xxyyzz_1 = prim_buffer_1_sssi[12];

    auto g_0_0_0_xxzzzz_1 = prim_buffer_1_sssi[14];

    auto g_0_0_0_xyyyyy_1 = prim_buffer_1_sssi[15];

    auto g_0_0_0_xyyyzz_1 = prim_buffer_1_sssi[17];

    auto g_0_0_0_xyyzzz_1 = prim_buffer_1_sssi[18];

    auto g_0_0_0_xzzzzz_1 = prim_buffer_1_sssi[20];

    auto g_0_0_0_yyyyyy_1 = prim_buffer_1_sssi[21];

    auto g_0_0_0_yyyyzz_1 = prim_buffer_1_sssi[23];

    auto g_0_0_0_yyyzzz_1 = prim_buffer_1_sssi[24];

    auto g_0_0_0_yyzzzz_1 = prim_buffer_1_sssi[25];

    auto g_0_0_0_yzzzzz_1 = prim_buffer_1_sssi[26];

    auto g_0_0_0_zzzzzz_1 = prim_buffer_1_sssi[27];

    /// Set up components of auxilary buffer : prim_buffer_0_sssk

    auto g_0_0_0_xxxxxxx_0 = prim_buffer_0_sssk[0];

    auto g_0_0_0_xxxxxxz_0 = prim_buffer_0_sssk[2];

    auto g_0_0_0_xxxxxyy_0 = prim_buffer_0_sssk[3];

    auto g_0_0_0_xxxxxzz_0 = prim_buffer_0_sssk[5];

    auto g_0_0_0_xxxxyyy_0 = prim_buffer_0_sssk[6];

    auto g_0_0_0_xxxxzzz_0 = prim_buffer_0_sssk[9];

    auto g_0_0_0_xxxyyyy_0 = prim_buffer_0_sssk[10];

    auto g_0_0_0_xxxyyzz_0 = prim_buffer_0_sssk[12];

    auto g_0_0_0_xxxzzzz_0 = prim_buffer_0_sssk[14];

    auto g_0_0_0_xxyyyyy_0 = prim_buffer_0_sssk[15];

    auto g_0_0_0_xxyyyzz_0 = prim_buffer_0_sssk[17];

    auto g_0_0_0_xxyyzzz_0 = prim_buffer_0_sssk[18];

    auto g_0_0_0_xxzzzzz_0 = prim_buffer_0_sssk[20];

    auto g_0_0_0_xyyyyyy_0 = prim_buffer_0_sssk[21];

    auto g_0_0_0_xyyyyzz_0 = prim_buffer_0_sssk[23];

    auto g_0_0_0_xyyyzzz_0 = prim_buffer_0_sssk[24];

    auto g_0_0_0_xyyzzzz_0 = prim_buffer_0_sssk[25];

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

    auto g_0_0_0_xxxxxxz_1 = prim_buffer_1_sssk[2];

    auto g_0_0_0_xxxxxyy_1 = prim_buffer_1_sssk[3];

    auto g_0_0_0_xxxxxzz_1 = prim_buffer_1_sssk[5];

    auto g_0_0_0_xxxxyyy_1 = prim_buffer_1_sssk[6];

    auto g_0_0_0_xxxxzzz_1 = prim_buffer_1_sssk[9];

    auto g_0_0_0_xxxyyyy_1 = prim_buffer_1_sssk[10];

    auto g_0_0_0_xxxyyzz_1 = prim_buffer_1_sssk[12];

    auto g_0_0_0_xxxzzzz_1 = prim_buffer_1_sssk[14];

    auto g_0_0_0_xxyyyyy_1 = prim_buffer_1_sssk[15];

    auto g_0_0_0_xxyyyzz_1 = prim_buffer_1_sssk[17];

    auto g_0_0_0_xxyyzzz_1 = prim_buffer_1_sssk[18];

    auto g_0_0_0_xxzzzzz_1 = prim_buffer_1_sssk[20];

    auto g_0_0_0_xyyyyyy_1 = prim_buffer_1_sssk[21];

    auto g_0_0_0_xyyyyzz_1 = prim_buffer_1_sssk[23];

    auto g_0_0_0_xyyyzzz_1 = prim_buffer_1_sssk[24];

    auto g_0_0_0_xyyzzzz_1 = prim_buffer_1_sssk[25];

    auto g_0_0_0_xzzzzzz_1 = prim_buffer_1_sssk[27];

    auto g_0_0_0_yyyyyyy_1 = prim_buffer_1_sssk[28];

    auto g_0_0_0_yyyyyyz_1 = prim_buffer_1_sssk[29];

    auto g_0_0_0_yyyyyzz_1 = prim_buffer_1_sssk[30];

    auto g_0_0_0_yyyyzzz_1 = prim_buffer_1_sssk[31];

    auto g_0_0_0_yyyzzzz_1 = prim_buffer_1_sssk[32];

    auto g_0_0_0_yyzzzzz_1 = prim_buffer_1_sssk[33];

    auto g_0_0_0_yzzzzzz_1 = prim_buffer_1_sssk[34];

    auto g_0_0_0_zzzzzzz_1 = prim_buffer_1_sssk[35];

    /// Set up components of targeted buffer : prim_buffer_0_sssl

    auto g_0_0_0_xxxxxxxx_0 = prim_buffer_0_sssl[0];

    auto g_0_0_0_xxxxxxxy_0 = prim_buffer_0_sssl[1];

    auto g_0_0_0_xxxxxxxz_0 = prim_buffer_0_sssl[2];

    auto g_0_0_0_xxxxxxyy_0 = prim_buffer_0_sssl[3];

    auto g_0_0_0_xxxxxxyz_0 = prim_buffer_0_sssl[4];

    auto g_0_0_0_xxxxxxzz_0 = prim_buffer_0_sssl[5];

    auto g_0_0_0_xxxxxyyy_0 = prim_buffer_0_sssl[6];

    auto g_0_0_0_xxxxxyyz_0 = prim_buffer_0_sssl[7];

    auto g_0_0_0_xxxxxyzz_0 = prim_buffer_0_sssl[8];

    auto g_0_0_0_xxxxxzzz_0 = prim_buffer_0_sssl[9];

    auto g_0_0_0_xxxxyyyy_0 = prim_buffer_0_sssl[10];

    auto g_0_0_0_xxxxyyyz_0 = prim_buffer_0_sssl[11];

    auto g_0_0_0_xxxxyyzz_0 = prim_buffer_0_sssl[12];

    auto g_0_0_0_xxxxyzzz_0 = prim_buffer_0_sssl[13];

    auto g_0_0_0_xxxxzzzz_0 = prim_buffer_0_sssl[14];

    auto g_0_0_0_xxxyyyyy_0 = prim_buffer_0_sssl[15];

    auto g_0_0_0_xxxyyyyz_0 = prim_buffer_0_sssl[16];

    auto g_0_0_0_xxxyyyzz_0 = prim_buffer_0_sssl[17];

    auto g_0_0_0_xxxyyzzz_0 = prim_buffer_0_sssl[18];

    auto g_0_0_0_xxxyzzzz_0 = prim_buffer_0_sssl[19];

    auto g_0_0_0_xxxzzzzz_0 = prim_buffer_0_sssl[20];

    auto g_0_0_0_xxyyyyyy_0 = prim_buffer_0_sssl[21];

    auto g_0_0_0_xxyyyyyz_0 = prim_buffer_0_sssl[22];

    auto g_0_0_0_xxyyyyzz_0 = prim_buffer_0_sssl[23];

    auto g_0_0_0_xxyyyzzz_0 = prim_buffer_0_sssl[24];

    auto g_0_0_0_xxyyzzzz_0 = prim_buffer_0_sssl[25];

    auto g_0_0_0_xxyzzzzz_0 = prim_buffer_0_sssl[26];

    auto g_0_0_0_xxzzzzzz_0 = prim_buffer_0_sssl[27];

    auto g_0_0_0_xyyyyyyy_0 = prim_buffer_0_sssl[28];

    auto g_0_0_0_xyyyyyyz_0 = prim_buffer_0_sssl[29];

    auto g_0_0_0_xyyyyyzz_0 = prim_buffer_0_sssl[30];

    auto g_0_0_0_xyyyyzzz_0 = prim_buffer_0_sssl[31];

    auto g_0_0_0_xyyyzzzz_0 = prim_buffer_0_sssl[32];

    auto g_0_0_0_xyyzzzzz_0 = prim_buffer_0_sssl[33];

    auto g_0_0_0_xyzzzzzz_0 = prim_buffer_0_sssl[34];

    auto g_0_0_0_xzzzzzzz_0 = prim_buffer_0_sssl[35];

    auto g_0_0_0_yyyyyyyy_0 = prim_buffer_0_sssl[36];

    auto g_0_0_0_yyyyyyyz_0 = prim_buffer_0_sssl[37];

    auto g_0_0_0_yyyyyyzz_0 = prim_buffer_0_sssl[38];

    auto g_0_0_0_yyyyyzzz_0 = prim_buffer_0_sssl[39];

    auto g_0_0_0_yyyyzzzz_0 = prim_buffer_0_sssl[40];

    auto g_0_0_0_yyyzzzzz_0 = prim_buffer_0_sssl[41];

    auto g_0_0_0_yyzzzzzz_0 = prim_buffer_0_sssl[42];

    auto g_0_0_0_yzzzzzzz_0 = prim_buffer_0_sssl[43];

    auto g_0_0_0_zzzzzzzz_0 = prim_buffer_0_sssl[44];

    #pragma omp simd aligned(g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxzz_0, g_0_0_0_xxxxzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxzzz_0, g_0_0_0_xxxzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyzz_0, g_0_0_0_xxyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxzzzz_0, g_0_0_0_xxzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyzz_0, g_0_0_0_xyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyzzz_0, g_0_0_0_xyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyzzzzzz_0, g_0_0_0_xzzzzz_0, g_0_0_0_xzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyzz_0, g_0_0_0_yyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyzzz_0, g_0_0_0_yyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyzzzz_0, g_0_0_0_yyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yzzzzz_0, g_0_0_0_yzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_zzzzzz_0, g_0_0_0_zzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_0_0_zzzzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxxxxxx_0[i] = 7.0 * g_0_0_0_xxxxxx_0[i] * fi_cd_0 - 7.0 * g_0_0_0_xxxxxx_1[i] * fti_cd_0 + g_0_0_0_xxxxxxx_0[i] * qd_x[i] + g_0_0_0_xxxxxxx_1[i] * wq_x[i];

        g_0_0_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxx_0[i] * qd_y[i] + g_0_0_0_xxxxxxx_1[i] * wq_y[i];

        g_0_0_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxx_0[i] * qd_z[i] + g_0_0_0_xxxxxxx_1[i] * wq_z[i];

        g_0_0_0_xxxxxxyy_0[i] = 5.0 * g_0_0_0_xxxxyy_0[i] * fi_cd_0 - 5.0 * g_0_0_0_xxxxyy_1[i] * fti_cd_0 + g_0_0_0_xxxxxyy_0[i] * qd_x[i] + g_0_0_0_xxxxxyy_1[i] * wq_x[i];

        g_0_0_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxz_0[i] * qd_y[i] + g_0_0_0_xxxxxxz_1[i] * wq_y[i];

        g_0_0_0_xxxxxxzz_0[i] = 5.0 * g_0_0_0_xxxxzz_0[i] * fi_cd_0 - 5.0 * g_0_0_0_xxxxzz_1[i] * fti_cd_0 + g_0_0_0_xxxxxzz_0[i] * qd_x[i] + g_0_0_0_xxxxxzz_1[i] * wq_x[i];

        g_0_0_0_xxxxxyyy_0[i] = 4.0 * g_0_0_0_xxxyyy_0[i] * fi_cd_0 - 4.0 * g_0_0_0_xxxyyy_1[i] * fti_cd_0 + g_0_0_0_xxxxyyy_0[i] * qd_x[i] + g_0_0_0_xxxxyyy_1[i] * wq_x[i];

        g_0_0_0_xxxxxyyz_0[i] = g_0_0_0_xxxxxyy_0[i] * qd_z[i] + g_0_0_0_xxxxxyy_1[i] * wq_z[i];

        g_0_0_0_xxxxxyzz_0[i] = g_0_0_0_xxxxxzz_0[i] * qd_y[i] + g_0_0_0_xxxxxzz_1[i] * wq_y[i];

        g_0_0_0_xxxxxzzz_0[i] = 4.0 * g_0_0_0_xxxzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_0_xxxzzz_1[i] * fti_cd_0 + g_0_0_0_xxxxzzz_0[i] * qd_x[i] + g_0_0_0_xxxxzzz_1[i] * wq_x[i];

        g_0_0_0_xxxxyyyy_0[i] = 3.0 * g_0_0_0_xxyyyy_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxyyyy_1[i] * fti_cd_0 + g_0_0_0_xxxyyyy_0[i] * qd_x[i] + g_0_0_0_xxxyyyy_1[i] * wq_x[i];

        g_0_0_0_xxxxyyyz_0[i] = g_0_0_0_xxxxyyy_0[i] * qd_z[i] + g_0_0_0_xxxxyyy_1[i] * wq_z[i];

        g_0_0_0_xxxxyyzz_0[i] = 3.0 * g_0_0_0_xxyyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxyyzz_1[i] * fti_cd_0 + g_0_0_0_xxxyyzz_0[i] * qd_x[i] + g_0_0_0_xxxyyzz_1[i] * wq_x[i];

        g_0_0_0_xxxxyzzz_0[i] = g_0_0_0_xxxxzzz_0[i] * qd_y[i] + g_0_0_0_xxxxzzz_1[i] * wq_y[i];

        g_0_0_0_xxxxzzzz_0[i] = 3.0 * g_0_0_0_xxzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxzzzz_1[i] * fti_cd_0 + g_0_0_0_xxxzzzz_0[i] * qd_x[i] + g_0_0_0_xxxzzzz_1[i] * wq_x[i];

        g_0_0_0_xxxyyyyy_0[i] = 2.0 * g_0_0_0_xyyyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyyyy_1[i] * fti_cd_0 + g_0_0_0_xxyyyyy_0[i] * qd_x[i] + g_0_0_0_xxyyyyy_1[i] * wq_x[i];

        g_0_0_0_xxxyyyyz_0[i] = g_0_0_0_xxxyyyy_0[i] * qd_z[i] + g_0_0_0_xxxyyyy_1[i] * wq_z[i];

        g_0_0_0_xxxyyyzz_0[i] = 2.0 * g_0_0_0_xyyyzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyyzz_1[i] * fti_cd_0 + g_0_0_0_xxyyyzz_0[i] * qd_x[i] + g_0_0_0_xxyyyzz_1[i] * wq_x[i];

        g_0_0_0_xxxyyzzz_0[i] = 2.0 * g_0_0_0_xyyzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyzzz_1[i] * fti_cd_0 + g_0_0_0_xxyyzzz_0[i] * qd_x[i] + g_0_0_0_xxyyzzz_1[i] * wq_x[i];

        g_0_0_0_xxxyzzzz_0[i] = g_0_0_0_xxxzzzz_0[i] * qd_y[i] + g_0_0_0_xxxzzzz_1[i] * wq_y[i];

        g_0_0_0_xxxzzzzz_0[i] = 2.0 * g_0_0_0_xzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xzzzzz_1[i] * fti_cd_0 + g_0_0_0_xxzzzzz_0[i] * qd_x[i] + g_0_0_0_xxzzzzz_1[i] * wq_x[i];

        g_0_0_0_xxyyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * fi_cd_0 - g_0_0_0_yyyyyy_1[i] * fti_cd_0 + g_0_0_0_xyyyyyy_0[i] * qd_x[i] + g_0_0_0_xyyyyyy_1[i] * wq_x[i];

        g_0_0_0_xxyyyyyz_0[i] = g_0_0_0_xxyyyyy_0[i] * qd_z[i] + g_0_0_0_xxyyyyy_1[i] * wq_z[i];

        g_0_0_0_xxyyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * fi_cd_0 - g_0_0_0_yyyyzz_1[i] * fti_cd_0 + g_0_0_0_xyyyyzz_0[i] * qd_x[i] + g_0_0_0_xyyyyzz_1[i] * wq_x[i];

        g_0_0_0_xxyyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * fi_cd_0 - g_0_0_0_yyyzzz_1[i] * fti_cd_0 + g_0_0_0_xyyyzzz_0[i] * qd_x[i] + g_0_0_0_xyyyzzz_1[i] * wq_x[i];

        g_0_0_0_xxyyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * fi_cd_0 - g_0_0_0_yyzzzz_1[i] * fti_cd_0 + g_0_0_0_xyyzzzz_0[i] * qd_x[i] + g_0_0_0_xyyzzzz_1[i] * wq_x[i];

        g_0_0_0_xxyzzzzz_0[i] = g_0_0_0_xxzzzzz_0[i] * qd_y[i] + g_0_0_0_xxzzzzz_1[i] * wq_y[i];

        g_0_0_0_xxzzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * fi_cd_0 - g_0_0_0_zzzzzz_1[i] * fti_cd_0 + g_0_0_0_xzzzzzz_0[i] * qd_x[i] + g_0_0_0_xzzzzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * qd_x[i] + g_0_0_0_yyyyyyy_1[i] * wq_x[i];

        g_0_0_0_xyyyyyyz_0[i] = g_0_0_0_yyyyyyz_0[i] * qd_x[i] + g_0_0_0_yyyyyyz_1[i] * wq_x[i];

        g_0_0_0_xyyyyyzz_0[i] = g_0_0_0_yyyyyzz_0[i] * qd_x[i] + g_0_0_0_yyyyyzz_1[i] * wq_x[i];

        g_0_0_0_xyyyyzzz_0[i] = g_0_0_0_yyyyzzz_0[i] * qd_x[i] + g_0_0_0_yyyyzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyzzzz_0[i] = g_0_0_0_yyyzzzz_0[i] * qd_x[i] + g_0_0_0_yyyzzzz_1[i] * wq_x[i];

        g_0_0_0_xyyzzzzz_0[i] = g_0_0_0_yyzzzzz_0[i] * qd_x[i] + g_0_0_0_yyzzzzz_1[i] * wq_x[i];

        g_0_0_0_xyzzzzzz_0[i] = g_0_0_0_yzzzzzz_0[i] * qd_x[i] + g_0_0_0_yzzzzzz_1[i] * wq_x[i];

        g_0_0_0_xzzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * qd_x[i] + g_0_0_0_zzzzzzz_1[i] * wq_x[i];

        g_0_0_0_yyyyyyyy_0[i] = 7.0 * g_0_0_0_yyyyyy_0[i] * fi_cd_0 - 7.0 * g_0_0_0_yyyyyy_1[i] * fti_cd_0 + g_0_0_0_yyyyyyy_0[i] * qd_y[i] + g_0_0_0_yyyyyyy_1[i] * wq_y[i];

        g_0_0_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyy_0[i] * qd_z[i] + g_0_0_0_yyyyyyy_1[i] * wq_z[i];

        g_0_0_0_yyyyyyzz_0[i] = 5.0 * g_0_0_0_yyyyzz_0[i] * fi_cd_0 - 5.0 * g_0_0_0_yyyyzz_1[i] * fti_cd_0 + g_0_0_0_yyyyyzz_0[i] * qd_y[i] + g_0_0_0_yyyyyzz_1[i] * wq_y[i];

        g_0_0_0_yyyyyzzz_0[i] = 4.0 * g_0_0_0_yyyzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_0_yyyzzz_1[i] * fti_cd_0 + g_0_0_0_yyyyzzz_0[i] * qd_y[i] + g_0_0_0_yyyyzzz_1[i] * wq_y[i];

        g_0_0_0_yyyyzzzz_0[i] = 3.0 * g_0_0_0_yyzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_yyzzzz_1[i] * fti_cd_0 + g_0_0_0_yyyzzzz_0[i] * qd_y[i] + g_0_0_0_yyyzzzz_1[i] * wq_y[i];

        g_0_0_0_yyyzzzzz_0[i] = 2.0 * g_0_0_0_yzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_yzzzzz_1[i] * fti_cd_0 + g_0_0_0_yyzzzzz_0[i] * qd_y[i] + g_0_0_0_yyzzzzz_1[i] * wq_y[i];

        g_0_0_0_yyzzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * fi_cd_0 - g_0_0_0_zzzzzz_1[i] * fti_cd_0 + g_0_0_0_yzzzzzz_0[i] * qd_y[i] + g_0_0_0_yzzzzzz_1[i] * wq_y[i];

        g_0_0_0_yzzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * qd_y[i] + g_0_0_0_zzzzzzz_1[i] * wq_y[i];

        g_0_0_0_zzzzzzzz_0[i] = 7.0 * g_0_0_0_zzzzzz_0[i] * fi_cd_0 - 7.0 * g_0_0_0_zzzzzz_1[i] * fti_cd_0 + g_0_0_0_zzzzzzz_0[i] * qd_z[i] + g_0_0_0_zzzzzzz_1[i] * wq_z[i];
    }
}

} // erirec namespace


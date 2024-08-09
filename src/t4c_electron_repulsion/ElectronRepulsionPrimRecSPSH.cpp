#include "ElectronRepulsionPrimRecSPSH.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsh(CSimdArray<double>& prim_buffer_0_spsh,
                                  const CSimdArray<double>& prim_buffer_1_sssg,
                                  const CSimdArray<double>& prim_buffer_0_sssh,
                                  const CSimdArray<double>& prim_buffer_1_sssh,
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
    const auto ndims = prim_buffer_0_spsh.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_1_sssg

    auto g_0_0_0_xxxx_1 = prim_buffer_1_sssg[0];

    auto g_0_0_0_xxxy_1 = prim_buffer_1_sssg[1];

    auto g_0_0_0_xxxz_1 = prim_buffer_1_sssg[2];

    auto g_0_0_0_xxyy_1 = prim_buffer_1_sssg[3];

    auto g_0_0_0_xxyz_1 = prim_buffer_1_sssg[4];

    auto g_0_0_0_xxzz_1 = prim_buffer_1_sssg[5];

    auto g_0_0_0_xyyy_1 = prim_buffer_1_sssg[6];

    auto g_0_0_0_xyyz_1 = prim_buffer_1_sssg[7];

    auto g_0_0_0_xyzz_1 = prim_buffer_1_sssg[8];

    auto g_0_0_0_xzzz_1 = prim_buffer_1_sssg[9];

    auto g_0_0_0_yyyy_1 = prim_buffer_1_sssg[10];

    auto g_0_0_0_yyyz_1 = prim_buffer_1_sssg[11];

    auto g_0_0_0_yyzz_1 = prim_buffer_1_sssg[12];

    auto g_0_0_0_yzzz_1 = prim_buffer_1_sssg[13];

    auto g_0_0_0_zzzz_1 = prim_buffer_1_sssg[14];

    /// Set up components of auxilary buffer : prim_buffer_0_sssh

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

    /// Set up components of auxilary buffer : prim_buffer_1_sssh

    auto g_0_0_0_xxxxx_1 = prim_buffer_1_sssh[0];

    auto g_0_0_0_xxxxy_1 = prim_buffer_1_sssh[1];

    auto g_0_0_0_xxxxz_1 = prim_buffer_1_sssh[2];

    auto g_0_0_0_xxxyy_1 = prim_buffer_1_sssh[3];

    auto g_0_0_0_xxxyz_1 = prim_buffer_1_sssh[4];

    auto g_0_0_0_xxxzz_1 = prim_buffer_1_sssh[5];

    auto g_0_0_0_xxyyy_1 = prim_buffer_1_sssh[6];

    auto g_0_0_0_xxyyz_1 = prim_buffer_1_sssh[7];

    auto g_0_0_0_xxyzz_1 = prim_buffer_1_sssh[8];

    auto g_0_0_0_xxzzz_1 = prim_buffer_1_sssh[9];

    auto g_0_0_0_xyyyy_1 = prim_buffer_1_sssh[10];

    auto g_0_0_0_xyyyz_1 = prim_buffer_1_sssh[11];

    auto g_0_0_0_xyyzz_1 = prim_buffer_1_sssh[12];

    auto g_0_0_0_xyzzz_1 = prim_buffer_1_sssh[13];

    auto g_0_0_0_xzzzz_1 = prim_buffer_1_sssh[14];

    auto g_0_0_0_yyyyy_1 = prim_buffer_1_sssh[15];

    auto g_0_0_0_yyyyz_1 = prim_buffer_1_sssh[16];

    auto g_0_0_0_yyyzz_1 = prim_buffer_1_sssh[17];

    auto g_0_0_0_yyzzz_1 = prim_buffer_1_sssh[18];

    auto g_0_0_0_yzzzz_1 = prim_buffer_1_sssh[19];

    auto g_0_0_0_zzzzz_1 = prim_buffer_1_sssh[20];

    /// Set up 0-21 components of targeted buffer : prim_buffer_0_spsh

    auto g_0_x_0_xxxxx_0 = prim_buffer_0_spsh[0];

    auto g_0_x_0_xxxxy_0 = prim_buffer_0_spsh[1];

    auto g_0_x_0_xxxxz_0 = prim_buffer_0_spsh[2];

    auto g_0_x_0_xxxyy_0 = prim_buffer_0_spsh[3];

    auto g_0_x_0_xxxyz_0 = prim_buffer_0_spsh[4];

    auto g_0_x_0_xxxzz_0 = prim_buffer_0_spsh[5];

    auto g_0_x_0_xxyyy_0 = prim_buffer_0_spsh[6];

    auto g_0_x_0_xxyyz_0 = prim_buffer_0_spsh[7];

    auto g_0_x_0_xxyzz_0 = prim_buffer_0_spsh[8];

    auto g_0_x_0_xxzzz_0 = prim_buffer_0_spsh[9];

    auto g_0_x_0_xyyyy_0 = prim_buffer_0_spsh[10];

    auto g_0_x_0_xyyyz_0 = prim_buffer_0_spsh[11];

    auto g_0_x_0_xyyzz_0 = prim_buffer_0_spsh[12];

    auto g_0_x_0_xyzzz_0 = prim_buffer_0_spsh[13];

    auto g_0_x_0_xzzzz_0 = prim_buffer_0_spsh[14];

    auto g_0_x_0_yyyyy_0 = prim_buffer_0_spsh[15];

    auto g_0_x_0_yyyyz_0 = prim_buffer_0_spsh[16];

    auto g_0_x_0_yyyzz_0 = prim_buffer_0_spsh[17];

    auto g_0_x_0_yyzzz_0 = prim_buffer_0_spsh[18];

    auto g_0_x_0_yzzzz_0 = prim_buffer_0_spsh[19];

    auto g_0_x_0_zzzzz_0 = prim_buffer_0_spsh[20];

    #pragma omp simd aligned(g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxx_1, g_0_0_0_xxxxy_0, g_0_0_0_xxxxy_1, g_0_0_0_xxxxz_0, g_0_0_0_xxxxz_1, g_0_0_0_xxxy_1, g_0_0_0_xxxyy_0, g_0_0_0_xxxyy_1, g_0_0_0_xxxyz_0, g_0_0_0_xxxyz_1, g_0_0_0_xxxz_1, g_0_0_0_xxxzz_0, g_0_0_0_xxxzz_1, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyy_1, g_0_0_0_xxyyz_0, g_0_0_0_xxyyz_1, g_0_0_0_xxyz_1, g_0_0_0_xxyzz_0, g_0_0_0_xxyzz_1, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xxzzz_1, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyy_1, g_0_0_0_xyyyz_0, g_0_0_0_xyyyz_1, g_0_0_0_xyyz_1, g_0_0_0_xyyzz_0, g_0_0_0_xyyzz_1, g_0_0_0_xyzz_1, g_0_0_0_xyzzz_0, g_0_0_0_xyzzz_1, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_xzzzz_1, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyy_1, g_0_0_0_yyyyz_0, g_0_0_0_yyyyz_1, g_0_0_0_yyyz_1, g_0_0_0_yyyzz_0, g_0_0_0_yyyzz_1, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yyzzz_1, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_yzzzz_1, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, g_0_0_0_zzzzz_1, g_0_x_0_xxxxx_0, g_0_x_0_xxxxy_0, g_0_x_0_xxxxz_0, g_0_x_0_xxxyy_0, g_0_x_0_xxxyz_0, g_0_x_0_xxxzz_0, g_0_x_0_xxyyy_0, g_0_x_0_xxyyz_0, g_0_x_0_xxyzz_0, g_0_x_0_xxzzz_0, g_0_x_0_xyyyy_0, g_0_x_0_xyyyz_0, g_0_x_0_xyyzz_0, g_0_x_0_xyzzz_0, g_0_x_0_xzzzz_0, g_0_x_0_yyyyy_0, g_0_x_0_yyyyz_0, g_0_x_0_yyyzz_0, g_0_x_0_yyzzz_0, g_0_x_0_yzzzz_0, g_0_x_0_zzzzz_0, wp_x  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxxxx_0[i] = 5.0 * g_0_0_0_xxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxx_0[i] * pb_x + g_0_0_0_xxxxx_1[i] * wp_x[i];

        g_0_x_0_xxxxy_0[i] = 4.0 * g_0_0_0_xxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxy_0[i] * pb_x + g_0_0_0_xxxxy_1[i] * wp_x[i];

        g_0_x_0_xxxxz_0[i] = 4.0 * g_0_0_0_xxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxz_0[i] * pb_x + g_0_0_0_xxxxz_1[i] * wp_x[i];

        g_0_x_0_xxxyy_0[i] = 3.0 * g_0_0_0_xxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyy_0[i] * pb_x + g_0_0_0_xxxyy_1[i] * wp_x[i];

        g_0_x_0_xxxyz_0[i] = 3.0 * g_0_0_0_xxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyz_0[i] * pb_x + g_0_0_0_xxxyz_1[i] * wp_x[i];

        g_0_x_0_xxxzz_0[i] = 3.0 * g_0_0_0_xxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzz_0[i] * pb_x + g_0_0_0_xxxzz_1[i] * wp_x[i];

        g_0_x_0_xxyyy_0[i] = 2.0 * g_0_0_0_xyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyy_0[i] * pb_x + g_0_0_0_xxyyy_1[i] * wp_x[i];

        g_0_x_0_xxyyz_0[i] = 2.0 * g_0_0_0_xyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyz_0[i] * pb_x + g_0_0_0_xxyyz_1[i] * wp_x[i];

        g_0_x_0_xxyzz_0[i] = 2.0 * g_0_0_0_xyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzz_0[i] * pb_x + g_0_0_0_xxyzz_1[i] * wp_x[i];

        g_0_x_0_xxzzz_0[i] = 2.0 * g_0_0_0_xzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzz_0[i] * pb_x + g_0_0_0_xxzzz_1[i] * wp_x[i];

        g_0_x_0_xyyyy_0[i] = g_0_0_0_yyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyy_0[i] * pb_x + g_0_0_0_xyyyy_1[i] * wp_x[i];

        g_0_x_0_xyyyz_0[i] = g_0_0_0_yyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyz_0[i] * pb_x + g_0_0_0_xyyyz_1[i] * wp_x[i];

        g_0_x_0_xyyzz_0[i] = g_0_0_0_yyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzz_0[i] * pb_x + g_0_0_0_xyyzz_1[i] * wp_x[i];

        g_0_x_0_xyzzz_0[i] = g_0_0_0_yzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzz_0[i] * pb_x + g_0_0_0_xyzzz_1[i] * wp_x[i];

        g_0_x_0_xzzzz_0[i] = g_0_0_0_zzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzz_0[i] * pb_x + g_0_0_0_xzzzz_1[i] * wp_x[i];

        g_0_x_0_yyyyy_0[i] = g_0_0_0_yyyyy_0[i] * pb_x + g_0_0_0_yyyyy_1[i] * wp_x[i];

        g_0_x_0_yyyyz_0[i] = g_0_0_0_yyyyz_0[i] * pb_x + g_0_0_0_yyyyz_1[i] * wp_x[i];

        g_0_x_0_yyyzz_0[i] = g_0_0_0_yyyzz_0[i] * pb_x + g_0_0_0_yyyzz_1[i] * wp_x[i];

        g_0_x_0_yyzzz_0[i] = g_0_0_0_yyzzz_0[i] * pb_x + g_0_0_0_yyzzz_1[i] * wp_x[i];

        g_0_x_0_yzzzz_0[i] = g_0_0_0_yzzzz_0[i] * pb_x + g_0_0_0_yzzzz_1[i] * wp_x[i];

        g_0_x_0_zzzzz_0[i] = g_0_0_0_zzzzz_0[i] * pb_x + g_0_0_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : prim_buffer_0_spsh

    auto g_0_y_0_xxxxx_0 = prim_buffer_0_spsh[21];

    auto g_0_y_0_xxxxy_0 = prim_buffer_0_spsh[22];

    auto g_0_y_0_xxxxz_0 = prim_buffer_0_spsh[23];

    auto g_0_y_0_xxxyy_0 = prim_buffer_0_spsh[24];

    auto g_0_y_0_xxxyz_0 = prim_buffer_0_spsh[25];

    auto g_0_y_0_xxxzz_0 = prim_buffer_0_spsh[26];

    auto g_0_y_0_xxyyy_0 = prim_buffer_0_spsh[27];

    auto g_0_y_0_xxyyz_0 = prim_buffer_0_spsh[28];

    auto g_0_y_0_xxyzz_0 = prim_buffer_0_spsh[29];

    auto g_0_y_0_xxzzz_0 = prim_buffer_0_spsh[30];

    auto g_0_y_0_xyyyy_0 = prim_buffer_0_spsh[31];

    auto g_0_y_0_xyyyz_0 = prim_buffer_0_spsh[32];

    auto g_0_y_0_xyyzz_0 = prim_buffer_0_spsh[33];

    auto g_0_y_0_xyzzz_0 = prim_buffer_0_spsh[34];

    auto g_0_y_0_xzzzz_0 = prim_buffer_0_spsh[35];

    auto g_0_y_0_yyyyy_0 = prim_buffer_0_spsh[36];

    auto g_0_y_0_yyyyz_0 = prim_buffer_0_spsh[37];

    auto g_0_y_0_yyyzz_0 = prim_buffer_0_spsh[38];

    auto g_0_y_0_yyzzz_0 = prim_buffer_0_spsh[39];

    auto g_0_y_0_yzzzz_0 = prim_buffer_0_spsh[40];

    auto g_0_y_0_zzzzz_0 = prim_buffer_0_spsh[41];

    #pragma omp simd aligned(g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxx_1, g_0_0_0_xxxxy_0, g_0_0_0_xxxxy_1, g_0_0_0_xxxxz_0, g_0_0_0_xxxxz_1, g_0_0_0_xxxy_1, g_0_0_0_xxxyy_0, g_0_0_0_xxxyy_1, g_0_0_0_xxxyz_0, g_0_0_0_xxxyz_1, g_0_0_0_xxxz_1, g_0_0_0_xxxzz_0, g_0_0_0_xxxzz_1, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyy_1, g_0_0_0_xxyyz_0, g_0_0_0_xxyyz_1, g_0_0_0_xxyz_1, g_0_0_0_xxyzz_0, g_0_0_0_xxyzz_1, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xxzzz_1, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyy_1, g_0_0_0_xyyyz_0, g_0_0_0_xyyyz_1, g_0_0_0_xyyz_1, g_0_0_0_xyyzz_0, g_0_0_0_xyyzz_1, g_0_0_0_xyzz_1, g_0_0_0_xyzzz_0, g_0_0_0_xyzzz_1, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_xzzzz_1, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyy_1, g_0_0_0_yyyyz_0, g_0_0_0_yyyyz_1, g_0_0_0_yyyz_1, g_0_0_0_yyyzz_0, g_0_0_0_yyyzz_1, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yyzzz_1, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_yzzzz_1, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, g_0_0_0_zzzzz_1, g_0_y_0_xxxxx_0, g_0_y_0_xxxxy_0, g_0_y_0_xxxxz_0, g_0_y_0_xxxyy_0, g_0_y_0_xxxyz_0, g_0_y_0_xxxzz_0, g_0_y_0_xxyyy_0, g_0_y_0_xxyyz_0, g_0_y_0_xxyzz_0, g_0_y_0_xxzzz_0, g_0_y_0_xyyyy_0, g_0_y_0_xyyyz_0, g_0_y_0_xyyzz_0, g_0_y_0_xyzzz_0, g_0_y_0_xzzzz_0, g_0_y_0_yyyyy_0, g_0_y_0_yyyyz_0, g_0_y_0_yyyzz_0, g_0_y_0_yyzzz_0, g_0_y_0_yzzzz_0, g_0_y_0_zzzzz_0, wp_y  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxxxx_0[i] = g_0_0_0_xxxxx_0[i] * pb_y + g_0_0_0_xxxxx_1[i] * wp_y[i];

        g_0_y_0_xxxxy_0[i] = g_0_0_0_xxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxy_0[i] * pb_y + g_0_0_0_xxxxy_1[i] * wp_y[i];

        g_0_y_0_xxxxz_0[i] = g_0_0_0_xxxxz_0[i] * pb_y + g_0_0_0_xxxxz_1[i] * wp_y[i];

        g_0_y_0_xxxyy_0[i] = 2.0 * g_0_0_0_xxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxyy_0[i] * pb_y + g_0_0_0_xxxyy_1[i] * wp_y[i];

        g_0_y_0_xxxyz_0[i] = g_0_0_0_xxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxyz_0[i] * pb_y + g_0_0_0_xxxyz_1[i] * wp_y[i];

        g_0_y_0_xxxzz_0[i] = g_0_0_0_xxxzz_0[i] * pb_y + g_0_0_0_xxxzz_1[i] * wp_y[i];

        g_0_y_0_xxyyy_0[i] = 3.0 * g_0_0_0_xxyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyy_0[i] * pb_y + g_0_0_0_xxyyy_1[i] * wp_y[i];

        g_0_y_0_xxyyz_0[i] = 2.0 * g_0_0_0_xxyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyz_0[i] * pb_y + g_0_0_0_xxyyz_1[i] * wp_y[i];

        g_0_y_0_xxyzz_0[i] = g_0_0_0_xxzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzz_0[i] * pb_y + g_0_0_0_xxyzz_1[i] * wp_y[i];

        g_0_y_0_xxzzz_0[i] = g_0_0_0_xxzzz_0[i] * pb_y + g_0_0_0_xxzzz_1[i] * wp_y[i];

        g_0_y_0_xyyyy_0[i] = 4.0 * g_0_0_0_xyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyy_0[i] * pb_y + g_0_0_0_xyyyy_1[i] * wp_y[i];

        g_0_y_0_xyyyz_0[i] = 3.0 * g_0_0_0_xyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyz_0[i] * pb_y + g_0_0_0_xyyyz_1[i] * wp_y[i];

        g_0_y_0_xyyzz_0[i] = 2.0 * g_0_0_0_xyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzz_0[i] * pb_y + g_0_0_0_xyyzz_1[i] * wp_y[i];

        g_0_y_0_xyzzz_0[i] = g_0_0_0_xzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzz_0[i] * pb_y + g_0_0_0_xyzzz_1[i] * wp_y[i];

        g_0_y_0_xzzzz_0[i] = g_0_0_0_xzzzz_0[i] * pb_y + g_0_0_0_xzzzz_1[i] * wp_y[i];

        g_0_y_0_yyyyy_0[i] = 5.0 * g_0_0_0_yyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyy_0[i] * pb_y + g_0_0_0_yyyyy_1[i] * wp_y[i];

        g_0_y_0_yyyyz_0[i] = 4.0 * g_0_0_0_yyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyz_0[i] * pb_y + g_0_0_0_yyyyz_1[i] * wp_y[i];

        g_0_y_0_yyyzz_0[i] = 3.0 * g_0_0_0_yyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzz_0[i] * pb_y + g_0_0_0_yyyzz_1[i] * wp_y[i];

        g_0_y_0_yyzzz_0[i] = 2.0 * g_0_0_0_yzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzz_0[i] * pb_y + g_0_0_0_yyzzz_1[i] * wp_y[i];

        g_0_y_0_yzzzz_0[i] = g_0_0_0_zzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzz_0[i] * pb_y + g_0_0_0_yzzzz_1[i] * wp_y[i];

        g_0_y_0_zzzzz_0[i] = g_0_0_0_zzzzz_0[i] * pb_y + g_0_0_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : prim_buffer_0_spsh

    auto g_0_z_0_xxxxx_0 = prim_buffer_0_spsh[42];

    auto g_0_z_0_xxxxy_0 = prim_buffer_0_spsh[43];

    auto g_0_z_0_xxxxz_0 = prim_buffer_0_spsh[44];

    auto g_0_z_0_xxxyy_0 = prim_buffer_0_spsh[45];

    auto g_0_z_0_xxxyz_0 = prim_buffer_0_spsh[46];

    auto g_0_z_0_xxxzz_0 = prim_buffer_0_spsh[47];

    auto g_0_z_0_xxyyy_0 = prim_buffer_0_spsh[48];

    auto g_0_z_0_xxyyz_0 = prim_buffer_0_spsh[49];

    auto g_0_z_0_xxyzz_0 = prim_buffer_0_spsh[50];

    auto g_0_z_0_xxzzz_0 = prim_buffer_0_spsh[51];

    auto g_0_z_0_xyyyy_0 = prim_buffer_0_spsh[52];

    auto g_0_z_0_xyyyz_0 = prim_buffer_0_spsh[53];

    auto g_0_z_0_xyyzz_0 = prim_buffer_0_spsh[54];

    auto g_0_z_0_xyzzz_0 = prim_buffer_0_spsh[55];

    auto g_0_z_0_xzzzz_0 = prim_buffer_0_spsh[56];

    auto g_0_z_0_yyyyy_0 = prim_buffer_0_spsh[57];

    auto g_0_z_0_yyyyz_0 = prim_buffer_0_spsh[58];

    auto g_0_z_0_yyyzz_0 = prim_buffer_0_spsh[59];

    auto g_0_z_0_yyzzz_0 = prim_buffer_0_spsh[60];

    auto g_0_z_0_yzzzz_0 = prim_buffer_0_spsh[61];

    auto g_0_z_0_zzzzz_0 = prim_buffer_0_spsh[62];

    #pragma omp simd aligned(g_0_0_0_xxxx_1, g_0_0_0_xxxxx_0, g_0_0_0_xxxxx_1, g_0_0_0_xxxxy_0, g_0_0_0_xxxxy_1, g_0_0_0_xxxxz_0, g_0_0_0_xxxxz_1, g_0_0_0_xxxy_1, g_0_0_0_xxxyy_0, g_0_0_0_xxxyy_1, g_0_0_0_xxxyz_0, g_0_0_0_xxxyz_1, g_0_0_0_xxxz_1, g_0_0_0_xxxzz_0, g_0_0_0_xxxzz_1, g_0_0_0_xxyy_1, g_0_0_0_xxyyy_0, g_0_0_0_xxyyy_1, g_0_0_0_xxyyz_0, g_0_0_0_xxyyz_1, g_0_0_0_xxyz_1, g_0_0_0_xxyzz_0, g_0_0_0_xxyzz_1, g_0_0_0_xxzz_1, g_0_0_0_xxzzz_0, g_0_0_0_xxzzz_1, g_0_0_0_xyyy_1, g_0_0_0_xyyyy_0, g_0_0_0_xyyyy_1, g_0_0_0_xyyyz_0, g_0_0_0_xyyyz_1, g_0_0_0_xyyz_1, g_0_0_0_xyyzz_0, g_0_0_0_xyyzz_1, g_0_0_0_xyzz_1, g_0_0_0_xyzzz_0, g_0_0_0_xyzzz_1, g_0_0_0_xzzz_1, g_0_0_0_xzzzz_0, g_0_0_0_xzzzz_1, g_0_0_0_yyyy_1, g_0_0_0_yyyyy_0, g_0_0_0_yyyyy_1, g_0_0_0_yyyyz_0, g_0_0_0_yyyyz_1, g_0_0_0_yyyz_1, g_0_0_0_yyyzz_0, g_0_0_0_yyyzz_1, g_0_0_0_yyzz_1, g_0_0_0_yyzzz_0, g_0_0_0_yyzzz_1, g_0_0_0_yzzz_1, g_0_0_0_yzzzz_0, g_0_0_0_yzzzz_1, g_0_0_0_zzzz_1, g_0_0_0_zzzzz_0, g_0_0_0_zzzzz_1, g_0_z_0_xxxxx_0, g_0_z_0_xxxxy_0, g_0_z_0_xxxxz_0, g_0_z_0_xxxyy_0, g_0_z_0_xxxyz_0, g_0_z_0_xxxzz_0, g_0_z_0_xxyyy_0, g_0_z_0_xxyyz_0, g_0_z_0_xxyzz_0, g_0_z_0_xxzzz_0, g_0_z_0_xyyyy_0, g_0_z_0_xyyyz_0, g_0_z_0_xyyzz_0, g_0_z_0_xyzzz_0, g_0_z_0_xzzzz_0, g_0_z_0_yyyyy_0, g_0_z_0_yyyyz_0, g_0_z_0_yyyzz_0, g_0_z_0_yyzzz_0, g_0_z_0_yzzzz_0, g_0_z_0_zzzzz_0, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxxxx_0[i] = g_0_0_0_xxxxx_0[i] * pb_z + g_0_0_0_xxxxx_1[i] * wp_z[i];

        g_0_z_0_xxxxy_0[i] = g_0_0_0_xxxxy_0[i] * pb_z + g_0_0_0_xxxxy_1[i] * wp_z[i];

        g_0_z_0_xxxxz_0[i] = g_0_0_0_xxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxz_0[i] * pb_z + g_0_0_0_xxxxz_1[i] * wp_z[i];

        g_0_z_0_xxxyy_0[i] = g_0_0_0_xxxyy_0[i] * pb_z + g_0_0_0_xxxyy_1[i] * wp_z[i];

        g_0_z_0_xxxyz_0[i] = g_0_0_0_xxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxyz_0[i] * pb_z + g_0_0_0_xxxyz_1[i] * wp_z[i];

        g_0_z_0_xxxzz_0[i] = 2.0 * g_0_0_0_xxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxzz_0[i] * pb_z + g_0_0_0_xxxzz_1[i] * wp_z[i];

        g_0_z_0_xxyyy_0[i] = g_0_0_0_xxyyy_0[i] * pb_z + g_0_0_0_xxyyy_1[i] * wp_z[i];

        g_0_z_0_xxyyz_0[i] = g_0_0_0_xxyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyz_0[i] * pb_z + g_0_0_0_xxyyz_1[i] * wp_z[i];

        g_0_z_0_xxyzz_0[i] = 2.0 * g_0_0_0_xxyz_1[i] * fi_abcd_0 + g_0_0_0_xxyzz_0[i] * pb_z + g_0_0_0_xxyzz_1[i] * wp_z[i];

        g_0_z_0_xxzzz_0[i] = 3.0 * g_0_0_0_xxzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzz_0[i] * pb_z + g_0_0_0_xxzzz_1[i] * wp_z[i];

        g_0_z_0_xyyyy_0[i] = g_0_0_0_xyyyy_0[i] * pb_z + g_0_0_0_xyyyy_1[i] * wp_z[i];

        g_0_z_0_xyyyz_0[i] = g_0_0_0_xyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyz_0[i] * pb_z + g_0_0_0_xyyyz_1[i] * wp_z[i];

        g_0_z_0_xyyzz_0[i] = 2.0 * g_0_0_0_xyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyzz_0[i] * pb_z + g_0_0_0_xyyzz_1[i] * wp_z[i];

        g_0_z_0_xyzzz_0[i] = 3.0 * g_0_0_0_xyzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzz_0[i] * pb_z + g_0_0_0_xyzzz_1[i] * wp_z[i];

        g_0_z_0_xzzzz_0[i] = 4.0 * g_0_0_0_xzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzz_0[i] * pb_z + g_0_0_0_xzzzz_1[i] * wp_z[i];

        g_0_z_0_yyyyy_0[i] = g_0_0_0_yyyyy_0[i] * pb_z + g_0_0_0_yyyyy_1[i] * wp_z[i];

        g_0_z_0_yyyyz_0[i] = g_0_0_0_yyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyz_0[i] * pb_z + g_0_0_0_yyyyz_1[i] * wp_z[i];

        g_0_z_0_yyyzz_0[i] = 2.0 * g_0_0_0_yyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyzz_0[i] * pb_z + g_0_0_0_yyyzz_1[i] * wp_z[i];

        g_0_z_0_yyzzz_0[i] = 3.0 * g_0_0_0_yyzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzz_0[i] * pb_z + g_0_0_0_yyzzz_1[i] * wp_z[i];

        g_0_z_0_yzzzz_0[i] = 4.0 * g_0_0_0_yzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzz_0[i] * pb_z + g_0_0_0_yzzzz_1[i] * wp_z[i];

        g_0_z_0_zzzzz_0[i] = 5.0 * g_0_0_0_zzzz_1[i] * fi_abcd_0 + g_0_0_0_zzzzz_0[i] * pb_z + g_0_0_0_zzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


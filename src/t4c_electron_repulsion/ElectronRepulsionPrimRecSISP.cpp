#include "ElectronRepulsionPrimRecSISP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sisp(CSimdArray<double>& prim_buffer_0_sisp,
                                  const CSimdArray<double>& prim_buffer_0_sgsp,
                                  const CSimdArray<double>& prim_buffer_1_sgsp,
                                  const CSimdArray<double>& prim_buffer_1_shss,
                                  const CSimdArray<double>& prim_buffer_0_shsp,
                                  const CSimdArray<double>& prim_buffer_1_shsp,
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
    const auto ndims = prim_buffer_0_sisp.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sgsp

    auto g_0_xxxx_0_x_0 = prim_buffer_0_sgsp[0];

    auto g_0_xxxx_0_y_0 = prim_buffer_0_sgsp[1];

    auto g_0_xxxx_0_z_0 = prim_buffer_0_sgsp[2];

    auto g_0_xxxy_0_x_0 = prim_buffer_0_sgsp[3];

    auto g_0_xxxz_0_x_0 = prim_buffer_0_sgsp[6];

    auto g_0_xxyy_0_x_0 = prim_buffer_0_sgsp[9];

    auto g_0_xxyy_0_y_0 = prim_buffer_0_sgsp[10];

    auto g_0_xxyy_0_z_0 = prim_buffer_0_sgsp[11];

    auto g_0_xxzz_0_x_0 = prim_buffer_0_sgsp[15];

    auto g_0_xxzz_0_y_0 = prim_buffer_0_sgsp[16];

    auto g_0_xxzz_0_z_0 = prim_buffer_0_sgsp[17];

    auto g_0_xyyy_0_y_0 = prim_buffer_0_sgsp[19];

    auto g_0_xyyy_0_z_0 = prim_buffer_0_sgsp[20];

    auto g_0_xzzz_0_y_0 = prim_buffer_0_sgsp[28];

    auto g_0_xzzz_0_z_0 = prim_buffer_0_sgsp[29];

    auto g_0_yyyy_0_x_0 = prim_buffer_0_sgsp[30];

    auto g_0_yyyy_0_y_0 = prim_buffer_0_sgsp[31];

    auto g_0_yyyy_0_z_0 = prim_buffer_0_sgsp[32];

    auto g_0_yyyz_0_y_0 = prim_buffer_0_sgsp[34];

    auto g_0_yyzz_0_x_0 = prim_buffer_0_sgsp[36];

    auto g_0_yyzz_0_y_0 = prim_buffer_0_sgsp[37];

    auto g_0_yyzz_0_z_0 = prim_buffer_0_sgsp[38];

    auto g_0_yzzz_0_x_0 = prim_buffer_0_sgsp[39];

    auto g_0_yzzz_0_z_0 = prim_buffer_0_sgsp[41];

    auto g_0_zzzz_0_x_0 = prim_buffer_0_sgsp[42];

    auto g_0_zzzz_0_y_0 = prim_buffer_0_sgsp[43];

    auto g_0_zzzz_0_z_0 = prim_buffer_0_sgsp[44];

    /// Set up components of auxilary buffer : prim_buffer_1_sgsp

    auto g_0_xxxx_0_x_1 = prim_buffer_1_sgsp[0];

    auto g_0_xxxx_0_y_1 = prim_buffer_1_sgsp[1];

    auto g_0_xxxx_0_z_1 = prim_buffer_1_sgsp[2];

    auto g_0_xxxy_0_x_1 = prim_buffer_1_sgsp[3];

    auto g_0_xxxz_0_x_1 = prim_buffer_1_sgsp[6];

    auto g_0_xxyy_0_x_1 = prim_buffer_1_sgsp[9];

    auto g_0_xxyy_0_y_1 = prim_buffer_1_sgsp[10];

    auto g_0_xxyy_0_z_1 = prim_buffer_1_sgsp[11];

    auto g_0_xxzz_0_x_1 = prim_buffer_1_sgsp[15];

    auto g_0_xxzz_0_y_1 = prim_buffer_1_sgsp[16];

    auto g_0_xxzz_0_z_1 = prim_buffer_1_sgsp[17];

    auto g_0_xyyy_0_y_1 = prim_buffer_1_sgsp[19];

    auto g_0_xyyy_0_z_1 = prim_buffer_1_sgsp[20];

    auto g_0_xzzz_0_y_1 = prim_buffer_1_sgsp[28];

    auto g_0_xzzz_0_z_1 = prim_buffer_1_sgsp[29];

    auto g_0_yyyy_0_x_1 = prim_buffer_1_sgsp[30];

    auto g_0_yyyy_0_y_1 = prim_buffer_1_sgsp[31];

    auto g_0_yyyy_0_z_1 = prim_buffer_1_sgsp[32];

    auto g_0_yyyz_0_y_1 = prim_buffer_1_sgsp[34];

    auto g_0_yyzz_0_x_1 = prim_buffer_1_sgsp[36];

    auto g_0_yyzz_0_y_1 = prim_buffer_1_sgsp[37];

    auto g_0_yyzz_0_z_1 = prim_buffer_1_sgsp[38];

    auto g_0_yzzz_0_x_1 = prim_buffer_1_sgsp[39];

    auto g_0_yzzz_0_z_1 = prim_buffer_1_sgsp[41];

    auto g_0_zzzz_0_x_1 = prim_buffer_1_sgsp[42];

    auto g_0_zzzz_0_y_1 = prim_buffer_1_sgsp[43];

    auto g_0_zzzz_0_z_1 = prim_buffer_1_sgsp[44];

    /// Set up components of auxilary buffer : prim_buffer_1_shss

    auto g_0_xxxxx_0_0_1 = prim_buffer_1_shss[0];

    auto g_0_xxxyy_0_0_1 = prim_buffer_1_shss[3];

    auto g_0_xxxzz_0_0_1 = prim_buffer_1_shss[5];

    auto g_0_xxyyy_0_0_1 = prim_buffer_1_shss[6];

    auto g_0_xxzzz_0_0_1 = prim_buffer_1_shss[9];

    auto g_0_yyyyy_0_0_1 = prim_buffer_1_shss[15];

    auto g_0_yyyzz_0_0_1 = prim_buffer_1_shss[17];

    auto g_0_yyzzz_0_0_1 = prim_buffer_1_shss[18];

    auto g_0_zzzzz_0_0_1 = prim_buffer_1_shss[20];

    /// Set up components of auxilary buffer : prim_buffer_0_shsp

    auto g_0_xxxxx_0_x_0 = prim_buffer_0_shsp[0];

    auto g_0_xxxxx_0_y_0 = prim_buffer_0_shsp[1];

    auto g_0_xxxxx_0_z_0 = prim_buffer_0_shsp[2];

    auto g_0_xxxxy_0_x_0 = prim_buffer_0_shsp[3];

    auto g_0_xxxxy_0_y_0 = prim_buffer_0_shsp[4];

    auto g_0_xxxxz_0_x_0 = prim_buffer_0_shsp[6];

    auto g_0_xxxxz_0_z_0 = prim_buffer_0_shsp[8];

    auto g_0_xxxyy_0_x_0 = prim_buffer_0_shsp[9];

    auto g_0_xxxyy_0_y_0 = prim_buffer_0_shsp[10];

    auto g_0_xxxyy_0_z_0 = prim_buffer_0_shsp[11];

    auto g_0_xxxzz_0_x_0 = prim_buffer_0_shsp[15];

    auto g_0_xxxzz_0_y_0 = prim_buffer_0_shsp[16];

    auto g_0_xxxzz_0_z_0 = prim_buffer_0_shsp[17];

    auto g_0_xxyyy_0_x_0 = prim_buffer_0_shsp[18];

    auto g_0_xxyyy_0_y_0 = prim_buffer_0_shsp[19];

    auto g_0_xxyyy_0_z_0 = prim_buffer_0_shsp[20];

    auto g_0_xxyzz_0_x_0 = prim_buffer_0_shsp[24];

    auto g_0_xxzzz_0_x_0 = prim_buffer_0_shsp[27];

    auto g_0_xxzzz_0_y_0 = prim_buffer_0_shsp[28];

    auto g_0_xxzzz_0_z_0 = prim_buffer_0_shsp[29];

    auto g_0_xyyyy_0_x_0 = prim_buffer_0_shsp[30];

    auto g_0_xyyyy_0_y_0 = prim_buffer_0_shsp[31];

    auto g_0_xyyyy_0_z_0 = prim_buffer_0_shsp[32];

    auto g_0_xyyzz_0_y_0 = prim_buffer_0_shsp[37];

    auto g_0_xyyzz_0_z_0 = prim_buffer_0_shsp[38];

    auto g_0_xzzzz_0_x_0 = prim_buffer_0_shsp[42];

    auto g_0_xzzzz_0_y_0 = prim_buffer_0_shsp[43];

    auto g_0_xzzzz_0_z_0 = prim_buffer_0_shsp[44];

    auto g_0_yyyyy_0_x_0 = prim_buffer_0_shsp[45];

    auto g_0_yyyyy_0_y_0 = prim_buffer_0_shsp[46];

    auto g_0_yyyyy_0_z_0 = prim_buffer_0_shsp[47];

    auto g_0_yyyyz_0_y_0 = prim_buffer_0_shsp[49];

    auto g_0_yyyyz_0_z_0 = prim_buffer_0_shsp[50];

    auto g_0_yyyzz_0_x_0 = prim_buffer_0_shsp[51];

    auto g_0_yyyzz_0_y_0 = prim_buffer_0_shsp[52];

    auto g_0_yyyzz_0_z_0 = prim_buffer_0_shsp[53];

    auto g_0_yyzzz_0_x_0 = prim_buffer_0_shsp[54];

    auto g_0_yyzzz_0_y_0 = prim_buffer_0_shsp[55];

    auto g_0_yyzzz_0_z_0 = prim_buffer_0_shsp[56];

    auto g_0_yzzzz_0_x_0 = prim_buffer_0_shsp[57];

    auto g_0_yzzzz_0_y_0 = prim_buffer_0_shsp[58];

    auto g_0_yzzzz_0_z_0 = prim_buffer_0_shsp[59];

    auto g_0_zzzzz_0_x_0 = prim_buffer_0_shsp[60];

    auto g_0_zzzzz_0_y_0 = prim_buffer_0_shsp[61];

    auto g_0_zzzzz_0_z_0 = prim_buffer_0_shsp[62];

    /// Set up components of auxilary buffer : prim_buffer_1_shsp

    auto g_0_xxxxx_0_x_1 = prim_buffer_1_shsp[0];

    auto g_0_xxxxx_0_y_1 = prim_buffer_1_shsp[1];

    auto g_0_xxxxx_0_z_1 = prim_buffer_1_shsp[2];

    auto g_0_xxxxy_0_x_1 = prim_buffer_1_shsp[3];

    auto g_0_xxxxy_0_y_1 = prim_buffer_1_shsp[4];

    auto g_0_xxxxz_0_x_1 = prim_buffer_1_shsp[6];

    auto g_0_xxxxz_0_z_1 = prim_buffer_1_shsp[8];

    auto g_0_xxxyy_0_x_1 = prim_buffer_1_shsp[9];

    auto g_0_xxxyy_0_y_1 = prim_buffer_1_shsp[10];

    auto g_0_xxxyy_0_z_1 = prim_buffer_1_shsp[11];

    auto g_0_xxxzz_0_x_1 = prim_buffer_1_shsp[15];

    auto g_0_xxxzz_0_y_1 = prim_buffer_1_shsp[16];

    auto g_0_xxxzz_0_z_1 = prim_buffer_1_shsp[17];

    auto g_0_xxyyy_0_x_1 = prim_buffer_1_shsp[18];

    auto g_0_xxyyy_0_y_1 = prim_buffer_1_shsp[19];

    auto g_0_xxyyy_0_z_1 = prim_buffer_1_shsp[20];

    auto g_0_xxyzz_0_x_1 = prim_buffer_1_shsp[24];

    auto g_0_xxzzz_0_x_1 = prim_buffer_1_shsp[27];

    auto g_0_xxzzz_0_y_1 = prim_buffer_1_shsp[28];

    auto g_0_xxzzz_0_z_1 = prim_buffer_1_shsp[29];

    auto g_0_xyyyy_0_x_1 = prim_buffer_1_shsp[30];

    auto g_0_xyyyy_0_y_1 = prim_buffer_1_shsp[31];

    auto g_0_xyyyy_0_z_1 = prim_buffer_1_shsp[32];

    auto g_0_xyyzz_0_y_1 = prim_buffer_1_shsp[37];

    auto g_0_xyyzz_0_z_1 = prim_buffer_1_shsp[38];

    auto g_0_xzzzz_0_x_1 = prim_buffer_1_shsp[42];

    auto g_0_xzzzz_0_y_1 = prim_buffer_1_shsp[43];

    auto g_0_xzzzz_0_z_1 = prim_buffer_1_shsp[44];

    auto g_0_yyyyy_0_x_1 = prim_buffer_1_shsp[45];

    auto g_0_yyyyy_0_y_1 = prim_buffer_1_shsp[46];

    auto g_0_yyyyy_0_z_1 = prim_buffer_1_shsp[47];

    auto g_0_yyyyz_0_y_1 = prim_buffer_1_shsp[49];

    auto g_0_yyyyz_0_z_1 = prim_buffer_1_shsp[50];

    auto g_0_yyyzz_0_x_1 = prim_buffer_1_shsp[51];

    auto g_0_yyyzz_0_y_1 = prim_buffer_1_shsp[52];

    auto g_0_yyyzz_0_z_1 = prim_buffer_1_shsp[53];

    auto g_0_yyzzz_0_x_1 = prim_buffer_1_shsp[54];

    auto g_0_yyzzz_0_y_1 = prim_buffer_1_shsp[55];

    auto g_0_yyzzz_0_z_1 = prim_buffer_1_shsp[56];

    auto g_0_yzzzz_0_x_1 = prim_buffer_1_shsp[57];

    auto g_0_yzzzz_0_y_1 = prim_buffer_1_shsp[58];

    auto g_0_yzzzz_0_z_1 = prim_buffer_1_shsp[59];

    auto g_0_zzzzz_0_x_1 = prim_buffer_1_shsp[60];

    auto g_0_zzzzz_0_y_1 = prim_buffer_1_shsp[61];

    auto g_0_zzzzz_0_z_1 = prim_buffer_1_shsp[62];

    /// Set up 0-3 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxxx_0_x_0 = prim_buffer_0_sisp[0];

    auto g_0_xxxxxx_0_y_0 = prim_buffer_0_sisp[1];

    auto g_0_xxxxxx_0_z_0 = prim_buffer_0_sisp[2];

    #pragma omp simd aligned(g_0_xxxx_0_x_0, g_0_xxxx_0_x_1, g_0_xxxx_0_y_0, g_0_xxxx_0_y_1, g_0_xxxx_0_z_0, g_0_xxxx_0_z_1, g_0_xxxxx_0_0_1, g_0_xxxxx_0_x_0, g_0_xxxxx_0_x_1, g_0_xxxxx_0_y_0, g_0_xxxxx_0_y_1, g_0_xxxxx_0_z_0, g_0_xxxxx_0_z_1, g_0_xxxxxx_0_x_0, g_0_xxxxxx_0_y_0, g_0_xxxxxx_0_z_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_x_0[i] = 5.0 * g_0_xxxx_0_x_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxx_0_x_0[i] * pb_x + g_0_xxxxx_0_x_1[i] * wp_x[i];

        g_0_xxxxxx_0_y_0[i] = 5.0 * g_0_xxxx_0_y_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_y_1[i] * fti_ab_0 + g_0_xxxxx_0_y_0[i] * pb_x + g_0_xxxxx_0_y_1[i] * wp_x[i];

        g_0_xxxxxx_0_z_0[i] = 5.0 * g_0_xxxx_0_z_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_z_1[i] * fti_ab_0 + g_0_xxxxx_0_z_0[i] * pb_x + g_0_xxxxx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxxy_0_x_0 = prim_buffer_0_sisp[3];

    auto g_0_xxxxxy_0_y_0 = prim_buffer_0_sisp[4];

    auto g_0_xxxxxy_0_z_0 = prim_buffer_0_sisp[5];

    #pragma omp simd aligned(g_0_xxxxx_0_0_1, g_0_xxxxx_0_x_0, g_0_xxxxx_0_x_1, g_0_xxxxx_0_y_0, g_0_xxxxx_0_y_1, g_0_xxxxx_0_z_0, g_0_xxxxx_0_z_1, g_0_xxxxxy_0_x_0, g_0_xxxxxy_0_y_0, g_0_xxxxxy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_x_0[i] = g_0_xxxxx_0_x_0[i] * pb_y + g_0_xxxxx_0_x_1[i] * wp_y[i];

        g_0_xxxxxy_0_y_0[i] = g_0_xxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxx_0_y_0[i] * pb_y + g_0_xxxxx_0_y_1[i] * wp_y[i];

        g_0_xxxxxy_0_z_0[i] = g_0_xxxxx_0_z_0[i] * pb_y + g_0_xxxxx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxxz_0_x_0 = prim_buffer_0_sisp[6];

    auto g_0_xxxxxz_0_y_0 = prim_buffer_0_sisp[7];

    auto g_0_xxxxxz_0_z_0 = prim_buffer_0_sisp[8];

    #pragma omp simd aligned(g_0_xxxxx_0_0_1, g_0_xxxxx_0_x_0, g_0_xxxxx_0_x_1, g_0_xxxxx_0_y_0, g_0_xxxxx_0_y_1, g_0_xxxxx_0_z_0, g_0_xxxxx_0_z_1, g_0_xxxxxz_0_x_0, g_0_xxxxxz_0_y_0, g_0_xxxxxz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_x_0[i] = g_0_xxxxx_0_x_0[i] * pb_z + g_0_xxxxx_0_x_1[i] * wp_z[i];

        g_0_xxxxxz_0_y_0[i] = g_0_xxxxx_0_y_0[i] * pb_z + g_0_xxxxx_0_y_1[i] * wp_z[i];

        g_0_xxxxxz_0_z_0[i] = g_0_xxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxx_0_z_0[i] * pb_z + g_0_xxxxx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxyy_0_x_0 = prim_buffer_0_sisp[9];

    auto g_0_xxxxyy_0_y_0 = prim_buffer_0_sisp[10];

    auto g_0_xxxxyy_0_z_0 = prim_buffer_0_sisp[11];

    #pragma omp simd aligned(g_0_xxxx_0_x_0, g_0_xxxx_0_x_1, g_0_xxxxy_0_x_0, g_0_xxxxy_0_x_1, g_0_xxxxyy_0_x_0, g_0_xxxxyy_0_y_0, g_0_xxxxyy_0_z_0, g_0_xxxyy_0_y_0, g_0_xxxyy_0_y_1, g_0_xxxyy_0_z_0, g_0_xxxyy_0_z_1, g_0_xxyy_0_y_0, g_0_xxyy_0_y_1, g_0_xxyy_0_z_0, g_0_xxyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_x_0[i] = g_0_xxxx_0_x_0[i] * fi_ab_0 - g_0_xxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxy_0_x_0[i] * pb_y + g_0_xxxxy_0_x_1[i] * wp_y[i];

        g_0_xxxxyy_0_y_0[i] = 3.0 * g_0_xxyy_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_y_1[i] * fti_ab_0 + g_0_xxxyy_0_y_0[i] * pb_x + g_0_xxxyy_0_y_1[i] * wp_x[i];

        g_0_xxxxyy_0_z_0[i] = 3.0 * g_0_xxyy_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_z_1[i] * fti_ab_0 + g_0_xxxyy_0_z_0[i] * pb_x + g_0_xxxyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxyz_0_x_0 = prim_buffer_0_sisp[12];

    auto g_0_xxxxyz_0_y_0 = prim_buffer_0_sisp[13];

    auto g_0_xxxxyz_0_z_0 = prim_buffer_0_sisp[14];

    #pragma omp simd aligned(g_0_xxxxy_0_y_0, g_0_xxxxy_0_y_1, g_0_xxxxyz_0_x_0, g_0_xxxxyz_0_y_0, g_0_xxxxyz_0_z_0, g_0_xxxxz_0_x_0, g_0_xxxxz_0_x_1, g_0_xxxxz_0_z_0, g_0_xxxxz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_xxxxyz_0_x_0[i] = g_0_xxxxz_0_x_0[i] * pb_y + g_0_xxxxz_0_x_1[i] * wp_y[i];

        g_0_xxxxyz_0_y_0[i] = g_0_xxxxy_0_y_0[i] * pb_z + g_0_xxxxy_0_y_1[i] * wp_z[i];

        g_0_xxxxyz_0_z_0[i] = g_0_xxxxz_0_z_0[i] * pb_y + g_0_xxxxz_0_z_1[i] * wp_y[i];
    }

    /// Set up 15-18 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxxzz_0_x_0 = prim_buffer_0_sisp[15];

    auto g_0_xxxxzz_0_y_0 = prim_buffer_0_sisp[16];

    auto g_0_xxxxzz_0_z_0 = prim_buffer_0_sisp[17];

    #pragma omp simd aligned(g_0_xxxx_0_x_0, g_0_xxxx_0_x_1, g_0_xxxxz_0_x_0, g_0_xxxxz_0_x_1, g_0_xxxxzz_0_x_0, g_0_xxxxzz_0_y_0, g_0_xxxxzz_0_z_0, g_0_xxxzz_0_y_0, g_0_xxxzz_0_y_1, g_0_xxxzz_0_z_0, g_0_xxxzz_0_z_1, g_0_xxzz_0_y_0, g_0_xxzz_0_y_1, g_0_xxzz_0_z_0, g_0_xxzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_x_0[i] = g_0_xxxx_0_x_0[i] * fi_ab_0 - g_0_xxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxz_0_x_0[i] * pb_z + g_0_xxxxz_0_x_1[i] * wp_z[i];

        g_0_xxxxzz_0_y_0[i] = 3.0 * g_0_xxzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_y_1[i] * fti_ab_0 + g_0_xxxzz_0_y_0[i] * pb_x + g_0_xxxzz_0_y_1[i] * wp_x[i];

        g_0_xxxxzz_0_z_0[i] = 3.0 * g_0_xxzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_z_1[i] * fti_ab_0 + g_0_xxxzz_0_z_0[i] * pb_x + g_0_xxxzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxyyy_0_x_0 = prim_buffer_0_sisp[18];

    auto g_0_xxxyyy_0_y_0 = prim_buffer_0_sisp[19];

    auto g_0_xxxyyy_0_z_0 = prim_buffer_0_sisp[20];

    #pragma omp simd aligned(g_0_xxxy_0_x_0, g_0_xxxy_0_x_1, g_0_xxxyy_0_x_0, g_0_xxxyy_0_x_1, g_0_xxxyyy_0_x_0, g_0_xxxyyy_0_y_0, g_0_xxxyyy_0_z_0, g_0_xxyyy_0_y_0, g_0_xxyyy_0_y_1, g_0_xxyyy_0_z_0, g_0_xxyyy_0_z_1, g_0_xyyy_0_y_0, g_0_xyyy_0_y_1, g_0_xyyy_0_z_0, g_0_xyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_x_0[i] = 2.0 * g_0_xxxy_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_x_1[i] * fti_ab_0 + g_0_xxxyy_0_x_0[i] * pb_y + g_0_xxxyy_0_x_1[i] * wp_y[i];

        g_0_xxxyyy_0_y_0[i] = 2.0 * g_0_xyyy_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_y_1[i] * fti_ab_0 + g_0_xxyyy_0_y_0[i] * pb_x + g_0_xxyyy_0_y_1[i] * wp_x[i];

        g_0_xxxyyy_0_z_0[i] = 2.0 * g_0_xyyy_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_z_1[i] * fti_ab_0 + g_0_xxyyy_0_z_0[i] * pb_x + g_0_xxyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 21-24 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxyyz_0_x_0 = prim_buffer_0_sisp[21];

    auto g_0_xxxyyz_0_y_0 = prim_buffer_0_sisp[22];

    auto g_0_xxxyyz_0_z_0 = prim_buffer_0_sisp[23];

    #pragma omp simd aligned(g_0_xxxyy_0_0_1, g_0_xxxyy_0_x_0, g_0_xxxyy_0_x_1, g_0_xxxyy_0_y_0, g_0_xxxyy_0_y_1, g_0_xxxyy_0_z_0, g_0_xxxyy_0_z_1, g_0_xxxyyz_0_x_0, g_0_xxxyyz_0_y_0, g_0_xxxyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_x_0[i] = g_0_xxxyy_0_x_0[i] * pb_z + g_0_xxxyy_0_x_1[i] * wp_z[i];

        g_0_xxxyyz_0_y_0[i] = g_0_xxxyy_0_y_0[i] * pb_z + g_0_xxxyy_0_y_1[i] * wp_z[i];

        g_0_xxxyyz_0_z_0[i] = g_0_xxxyy_0_0_1[i] * fi_abcd_0 + g_0_xxxyy_0_z_0[i] * pb_z + g_0_xxxyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 24-27 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxyzz_0_x_0 = prim_buffer_0_sisp[24];

    auto g_0_xxxyzz_0_y_0 = prim_buffer_0_sisp[25];

    auto g_0_xxxyzz_0_z_0 = prim_buffer_0_sisp[26];

    #pragma omp simd aligned(g_0_xxxyzz_0_x_0, g_0_xxxyzz_0_y_0, g_0_xxxyzz_0_z_0, g_0_xxxzz_0_0_1, g_0_xxxzz_0_x_0, g_0_xxxzz_0_x_1, g_0_xxxzz_0_y_0, g_0_xxxzz_0_y_1, g_0_xxxzz_0_z_0, g_0_xxxzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_x_0[i] = g_0_xxxzz_0_x_0[i] * pb_y + g_0_xxxzz_0_x_1[i] * wp_y[i];

        g_0_xxxyzz_0_y_0[i] = g_0_xxxzz_0_0_1[i] * fi_abcd_0 + g_0_xxxzz_0_y_0[i] * pb_y + g_0_xxxzz_0_y_1[i] * wp_y[i];

        g_0_xxxyzz_0_z_0[i] = g_0_xxxzz_0_z_0[i] * pb_y + g_0_xxxzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 27-30 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxxzzz_0_x_0 = prim_buffer_0_sisp[27];

    auto g_0_xxxzzz_0_y_0 = prim_buffer_0_sisp[28];

    auto g_0_xxxzzz_0_z_0 = prim_buffer_0_sisp[29];

    #pragma omp simd aligned(g_0_xxxz_0_x_0, g_0_xxxz_0_x_1, g_0_xxxzz_0_x_0, g_0_xxxzz_0_x_1, g_0_xxxzzz_0_x_0, g_0_xxxzzz_0_y_0, g_0_xxxzzz_0_z_0, g_0_xxzzz_0_y_0, g_0_xxzzz_0_y_1, g_0_xxzzz_0_z_0, g_0_xxzzz_0_z_1, g_0_xzzz_0_y_0, g_0_xzzz_0_y_1, g_0_xzzz_0_z_0, g_0_xzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_x_0[i] = 2.0 * g_0_xxxz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_x_1[i] * fti_ab_0 + g_0_xxxzz_0_x_0[i] * pb_z + g_0_xxxzz_0_x_1[i] * wp_z[i];

        g_0_xxxzzz_0_y_0[i] = 2.0 * g_0_xzzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_y_1[i] * fti_ab_0 + g_0_xxzzz_0_y_0[i] * pb_x + g_0_xxzzz_0_y_1[i] * wp_x[i];

        g_0_xxxzzz_0_z_0[i] = 2.0 * g_0_xzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_z_1[i] * fti_ab_0 + g_0_xxzzz_0_z_0[i] * pb_x + g_0_xxzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 30-33 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxyyyy_0_x_0 = prim_buffer_0_sisp[30];

    auto g_0_xxyyyy_0_y_0 = prim_buffer_0_sisp[31];

    auto g_0_xxyyyy_0_z_0 = prim_buffer_0_sisp[32];

    #pragma omp simd aligned(g_0_xxyy_0_x_0, g_0_xxyy_0_x_1, g_0_xxyyy_0_x_0, g_0_xxyyy_0_x_1, g_0_xxyyyy_0_x_0, g_0_xxyyyy_0_y_0, g_0_xxyyyy_0_z_0, g_0_xyyyy_0_y_0, g_0_xyyyy_0_y_1, g_0_xyyyy_0_z_0, g_0_xyyyy_0_z_1, g_0_yyyy_0_y_0, g_0_yyyy_0_y_1, g_0_yyyy_0_z_0, g_0_yyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_x_0[i] = 3.0 * g_0_xxyy_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_x_1[i] * fti_ab_0 + g_0_xxyyy_0_x_0[i] * pb_y + g_0_xxyyy_0_x_1[i] * wp_y[i];

        g_0_xxyyyy_0_y_0[i] = g_0_yyyy_0_y_0[i] * fi_ab_0 - g_0_yyyy_0_y_1[i] * fti_ab_0 + g_0_xyyyy_0_y_0[i] * pb_x + g_0_xyyyy_0_y_1[i] * wp_x[i];

        g_0_xxyyyy_0_z_0[i] = g_0_yyyy_0_z_0[i] * fi_ab_0 - g_0_yyyy_0_z_1[i] * fti_ab_0 + g_0_xyyyy_0_z_0[i] * pb_x + g_0_xyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 33-36 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxyyyz_0_x_0 = prim_buffer_0_sisp[33];

    auto g_0_xxyyyz_0_y_0 = prim_buffer_0_sisp[34];

    auto g_0_xxyyyz_0_z_0 = prim_buffer_0_sisp[35];

    #pragma omp simd aligned(g_0_xxyyy_0_0_1, g_0_xxyyy_0_x_0, g_0_xxyyy_0_x_1, g_0_xxyyy_0_y_0, g_0_xxyyy_0_y_1, g_0_xxyyy_0_z_0, g_0_xxyyy_0_z_1, g_0_xxyyyz_0_x_0, g_0_xxyyyz_0_y_0, g_0_xxyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_x_0[i] = g_0_xxyyy_0_x_0[i] * pb_z + g_0_xxyyy_0_x_1[i] * wp_z[i];

        g_0_xxyyyz_0_y_0[i] = g_0_xxyyy_0_y_0[i] * pb_z + g_0_xxyyy_0_y_1[i] * wp_z[i];

        g_0_xxyyyz_0_z_0[i] = g_0_xxyyy_0_0_1[i] * fi_abcd_0 + g_0_xxyyy_0_z_0[i] * pb_z + g_0_xxyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 36-39 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxyyzz_0_x_0 = prim_buffer_0_sisp[36];

    auto g_0_xxyyzz_0_y_0 = prim_buffer_0_sisp[37];

    auto g_0_xxyyzz_0_z_0 = prim_buffer_0_sisp[38];

    #pragma omp simd aligned(g_0_xxyyzz_0_x_0, g_0_xxyyzz_0_y_0, g_0_xxyyzz_0_z_0, g_0_xxyzz_0_x_0, g_0_xxyzz_0_x_1, g_0_xxzz_0_x_0, g_0_xxzz_0_x_1, g_0_xyyzz_0_y_0, g_0_xyyzz_0_y_1, g_0_xyyzz_0_z_0, g_0_xyyzz_0_z_1, g_0_yyzz_0_y_0, g_0_yyzz_0_y_1, g_0_yyzz_0_z_0, g_0_yyzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_x_0[i] = g_0_xxzz_0_x_0[i] * fi_ab_0 - g_0_xxzz_0_x_1[i] * fti_ab_0 + g_0_xxyzz_0_x_0[i] * pb_y + g_0_xxyzz_0_x_1[i] * wp_y[i];

        g_0_xxyyzz_0_y_0[i] = g_0_yyzz_0_y_0[i] * fi_ab_0 - g_0_yyzz_0_y_1[i] * fti_ab_0 + g_0_xyyzz_0_y_0[i] * pb_x + g_0_xyyzz_0_y_1[i] * wp_x[i];

        g_0_xxyyzz_0_z_0[i] = g_0_yyzz_0_z_0[i] * fi_ab_0 - g_0_yyzz_0_z_1[i] * fti_ab_0 + g_0_xyyzz_0_z_0[i] * pb_x + g_0_xyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 39-42 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxyzzz_0_x_0 = prim_buffer_0_sisp[39];

    auto g_0_xxyzzz_0_y_0 = prim_buffer_0_sisp[40];

    auto g_0_xxyzzz_0_z_0 = prim_buffer_0_sisp[41];

    #pragma omp simd aligned(g_0_xxyzzz_0_x_0, g_0_xxyzzz_0_y_0, g_0_xxyzzz_0_z_0, g_0_xxzzz_0_0_1, g_0_xxzzz_0_x_0, g_0_xxzzz_0_x_1, g_0_xxzzz_0_y_0, g_0_xxzzz_0_y_1, g_0_xxzzz_0_z_0, g_0_xxzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_x_0[i] = g_0_xxzzz_0_x_0[i] * pb_y + g_0_xxzzz_0_x_1[i] * wp_y[i];

        g_0_xxyzzz_0_y_0[i] = g_0_xxzzz_0_0_1[i] * fi_abcd_0 + g_0_xxzzz_0_y_0[i] * pb_y + g_0_xxzzz_0_y_1[i] * wp_y[i];

        g_0_xxyzzz_0_z_0[i] = g_0_xxzzz_0_z_0[i] * pb_y + g_0_xxzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 42-45 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xxzzzz_0_x_0 = prim_buffer_0_sisp[42];

    auto g_0_xxzzzz_0_y_0 = prim_buffer_0_sisp[43];

    auto g_0_xxzzzz_0_z_0 = prim_buffer_0_sisp[44];

    #pragma omp simd aligned(g_0_xxzz_0_x_0, g_0_xxzz_0_x_1, g_0_xxzzz_0_x_0, g_0_xxzzz_0_x_1, g_0_xxzzzz_0_x_0, g_0_xxzzzz_0_y_0, g_0_xxzzzz_0_z_0, g_0_xzzzz_0_y_0, g_0_xzzzz_0_y_1, g_0_xzzzz_0_z_0, g_0_xzzzz_0_z_1, g_0_zzzz_0_y_0, g_0_zzzz_0_y_1, g_0_zzzz_0_z_0, g_0_zzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_x_0[i] = 3.0 * g_0_xxzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_x_1[i] * fti_ab_0 + g_0_xxzzz_0_x_0[i] * pb_z + g_0_xxzzz_0_x_1[i] * wp_z[i];

        g_0_xxzzzz_0_y_0[i] = g_0_zzzz_0_y_0[i] * fi_ab_0 - g_0_zzzz_0_y_1[i] * fti_ab_0 + g_0_xzzzz_0_y_0[i] * pb_x + g_0_xzzzz_0_y_1[i] * wp_x[i];

        g_0_xxzzzz_0_z_0[i] = g_0_zzzz_0_z_0[i] * fi_ab_0 - g_0_zzzz_0_z_1[i] * fti_ab_0 + g_0_xzzzz_0_z_0[i] * pb_x + g_0_xzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 45-48 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xyyyyy_0_x_0 = prim_buffer_0_sisp[45];

    auto g_0_xyyyyy_0_y_0 = prim_buffer_0_sisp[46];

    auto g_0_xyyyyy_0_z_0 = prim_buffer_0_sisp[47];

    #pragma omp simd aligned(g_0_xyyyyy_0_x_0, g_0_xyyyyy_0_y_0, g_0_xyyyyy_0_z_0, g_0_yyyyy_0_0_1, g_0_yyyyy_0_x_0, g_0_yyyyy_0_x_1, g_0_yyyyy_0_y_0, g_0_yyyyy_0_y_1, g_0_yyyyy_0_z_0, g_0_yyyyy_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_x_0[i] = g_0_yyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyy_0_x_0[i] * pb_x + g_0_yyyyy_0_x_1[i] * wp_x[i];

        g_0_xyyyyy_0_y_0[i] = g_0_yyyyy_0_y_0[i] * pb_x + g_0_yyyyy_0_y_1[i] * wp_x[i];

        g_0_xyyyyy_0_z_0[i] = g_0_yyyyy_0_z_0[i] * pb_x + g_0_yyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 48-51 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xyyyyz_0_x_0 = prim_buffer_0_sisp[48];

    auto g_0_xyyyyz_0_y_0 = prim_buffer_0_sisp[49];

    auto g_0_xyyyyz_0_z_0 = prim_buffer_0_sisp[50];

    #pragma omp simd aligned(g_0_xyyyy_0_x_0, g_0_xyyyy_0_x_1, g_0_xyyyyz_0_x_0, g_0_xyyyyz_0_y_0, g_0_xyyyyz_0_z_0, g_0_yyyyz_0_y_0, g_0_yyyyz_0_y_1, g_0_yyyyz_0_z_0, g_0_yyyyz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_xyyyyz_0_x_0[i] = g_0_xyyyy_0_x_0[i] * pb_z + g_0_xyyyy_0_x_1[i] * wp_z[i];

        g_0_xyyyyz_0_y_0[i] = g_0_yyyyz_0_y_0[i] * pb_x + g_0_yyyyz_0_y_1[i] * wp_x[i];

        g_0_xyyyyz_0_z_0[i] = g_0_yyyyz_0_z_0[i] * pb_x + g_0_yyyyz_0_z_1[i] * wp_x[i];
    }

    /// Set up 51-54 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xyyyzz_0_x_0 = prim_buffer_0_sisp[51];

    auto g_0_xyyyzz_0_y_0 = prim_buffer_0_sisp[52];

    auto g_0_xyyyzz_0_z_0 = prim_buffer_0_sisp[53];

    #pragma omp simd aligned(g_0_xyyyzz_0_x_0, g_0_xyyyzz_0_y_0, g_0_xyyyzz_0_z_0, g_0_yyyzz_0_0_1, g_0_yyyzz_0_x_0, g_0_yyyzz_0_x_1, g_0_yyyzz_0_y_0, g_0_yyyzz_0_y_1, g_0_yyyzz_0_z_0, g_0_yyyzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_x_0[i] = g_0_yyyzz_0_0_1[i] * fi_abcd_0 + g_0_yyyzz_0_x_0[i] * pb_x + g_0_yyyzz_0_x_1[i] * wp_x[i];

        g_0_xyyyzz_0_y_0[i] = g_0_yyyzz_0_y_0[i] * pb_x + g_0_yyyzz_0_y_1[i] * wp_x[i];

        g_0_xyyyzz_0_z_0[i] = g_0_yyyzz_0_z_0[i] * pb_x + g_0_yyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 54-57 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xyyzzz_0_x_0 = prim_buffer_0_sisp[54];

    auto g_0_xyyzzz_0_y_0 = prim_buffer_0_sisp[55];

    auto g_0_xyyzzz_0_z_0 = prim_buffer_0_sisp[56];

    #pragma omp simd aligned(g_0_xyyzzz_0_x_0, g_0_xyyzzz_0_y_0, g_0_xyyzzz_0_z_0, g_0_yyzzz_0_0_1, g_0_yyzzz_0_x_0, g_0_yyzzz_0_x_1, g_0_yyzzz_0_y_0, g_0_yyzzz_0_y_1, g_0_yyzzz_0_z_0, g_0_yyzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_x_0[i] = g_0_yyzzz_0_0_1[i] * fi_abcd_0 + g_0_yyzzz_0_x_0[i] * pb_x + g_0_yyzzz_0_x_1[i] * wp_x[i];

        g_0_xyyzzz_0_y_0[i] = g_0_yyzzz_0_y_0[i] * pb_x + g_0_yyzzz_0_y_1[i] * wp_x[i];

        g_0_xyyzzz_0_z_0[i] = g_0_yyzzz_0_z_0[i] * pb_x + g_0_yyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 57-60 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xyzzzz_0_x_0 = prim_buffer_0_sisp[57];

    auto g_0_xyzzzz_0_y_0 = prim_buffer_0_sisp[58];

    auto g_0_xyzzzz_0_z_0 = prim_buffer_0_sisp[59];

    #pragma omp simd aligned(g_0_xyzzzz_0_x_0, g_0_xyzzzz_0_y_0, g_0_xyzzzz_0_z_0, g_0_xzzzz_0_x_0, g_0_xzzzz_0_x_1, g_0_yzzzz_0_y_0, g_0_yzzzz_0_y_1, g_0_yzzzz_0_z_0, g_0_yzzzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_xyzzzz_0_x_0[i] = g_0_xzzzz_0_x_0[i] * pb_y + g_0_xzzzz_0_x_1[i] * wp_y[i];

        g_0_xyzzzz_0_y_0[i] = g_0_yzzzz_0_y_0[i] * pb_x + g_0_yzzzz_0_y_1[i] * wp_x[i];

        g_0_xyzzzz_0_z_0[i] = g_0_yzzzz_0_z_0[i] * pb_x + g_0_yzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 60-63 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_xzzzzz_0_x_0 = prim_buffer_0_sisp[60];

    auto g_0_xzzzzz_0_y_0 = prim_buffer_0_sisp[61];

    auto g_0_xzzzzz_0_z_0 = prim_buffer_0_sisp[62];

    #pragma omp simd aligned(g_0_xzzzzz_0_x_0, g_0_xzzzzz_0_y_0, g_0_xzzzzz_0_z_0, g_0_zzzzz_0_0_1, g_0_zzzzz_0_x_0, g_0_zzzzz_0_x_1, g_0_zzzzz_0_y_0, g_0_zzzzz_0_y_1, g_0_zzzzz_0_z_0, g_0_zzzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_x_0[i] = g_0_zzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzz_0_x_0[i] * pb_x + g_0_zzzzz_0_x_1[i] * wp_x[i];

        g_0_xzzzzz_0_y_0[i] = g_0_zzzzz_0_y_0[i] * pb_x + g_0_zzzzz_0_y_1[i] * wp_x[i];

        g_0_xzzzzz_0_z_0[i] = g_0_zzzzz_0_z_0[i] * pb_x + g_0_zzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 63-66 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yyyyyy_0_x_0 = prim_buffer_0_sisp[63];

    auto g_0_yyyyyy_0_y_0 = prim_buffer_0_sisp[64];

    auto g_0_yyyyyy_0_z_0 = prim_buffer_0_sisp[65];

    #pragma omp simd aligned(g_0_yyyy_0_x_0, g_0_yyyy_0_x_1, g_0_yyyy_0_y_0, g_0_yyyy_0_y_1, g_0_yyyy_0_z_0, g_0_yyyy_0_z_1, g_0_yyyyy_0_0_1, g_0_yyyyy_0_x_0, g_0_yyyyy_0_x_1, g_0_yyyyy_0_y_0, g_0_yyyyy_0_y_1, g_0_yyyyy_0_z_0, g_0_yyyyy_0_z_1, g_0_yyyyyy_0_x_0, g_0_yyyyyy_0_y_0, g_0_yyyyyy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_x_0[i] = 5.0 * g_0_yyyy_0_x_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_x_1[i] * fti_ab_0 + g_0_yyyyy_0_x_0[i] * pb_y + g_0_yyyyy_0_x_1[i] * wp_y[i];

        g_0_yyyyyy_0_y_0[i] = 5.0 * g_0_yyyy_0_y_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyy_0_y_0[i] * pb_y + g_0_yyyyy_0_y_1[i] * wp_y[i];

        g_0_yyyyyy_0_z_0[i] = 5.0 * g_0_yyyy_0_z_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_z_1[i] * fti_ab_0 + g_0_yyyyy_0_z_0[i] * pb_y + g_0_yyyyy_0_z_1[i] * wp_y[i];
    }

    /// Set up 66-69 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yyyyyz_0_x_0 = prim_buffer_0_sisp[66];

    auto g_0_yyyyyz_0_y_0 = prim_buffer_0_sisp[67];

    auto g_0_yyyyyz_0_z_0 = prim_buffer_0_sisp[68];

    #pragma omp simd aligned(g_0_yyyyy_0_0_1, g_0_yyyyy_0_x_0, g_0_yyyyy_0_x_1, g_0_yyyyy_0_y_0, g_0_yyyyy_0_y_1, g_0_yyyyy_0_z_0, g_0_yyyyy_0_z_1, g_0_yyyyyz_0_x_0, g_0_yyyyyz_0_y_0, g_0_yyyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_x_0[i] = g_0_yyyyy_0_x_0[i] * pb_z + g_0_yyyyy_0_x_1[i] * wp_z[i];

        g_0_yyyyyz_0_y_0[i] = g_0_yyyyy_0_y_0[i] * pb_z + g_0_yyyyy_0_y_1[i] * wp_z[i];

        g_0_yyyyyz_0_z_0[i] = g_0_yyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyy_0_z_0[i] * pb_z + g_0_yyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 69-72 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yyyyzz_0_x_0 = prim_buffer_0_sisp[69];

    auto g_0_yyyyzz_0_y_0 = prim_buffer_0_sisp[70];

    auto g_0_yyyyzz_0_z_0 = prim_buffer_0_sisp[71];

    #pragma omp simd aligned(g_0_yyyy_0_y_0, g_0_yyyy_0_y_1, g_0_yyyyz_0_y_0, g_0_yyyyz_0_y_1, g_0_yyyyzz_0_x_0, g_0_yyyyzz_0_y_0, g_0_yyyyzz_0_z_0, g_0_yyyzz_0_x_0, g_0_yyyzz_0_x_1, g_0_yyyzz_0_z_0, g_0_yyyzz_0_z_1, g_0_yyzz_0_x_0, g_0_yyzz_0_x_1, g_0_yyzz_0_z_0, g_0_yyzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_x_0[i] = 3.0 * g_0_yyzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_x_1[i] * fti_ab_0 + g_0_yyyzz_0_x_0[i] * pb_y + g_0_yyyzz_0_x_1[i] * wp_y[i];

        g_0_yyyyzz_0_y_0[i] = g_0_yyyy_0_y_0[i] * fi_ab_0 - g_0_yyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyz_0_y_0[i] * pb_z + g_0_yyyyz_0_y_1[i] * wp_z[i];

        g_0_yyyyzz_0_z_0[i] = 3.0 * g_0_yyzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_z_1[i] * fti_ab_0 + g_0_yyyzz_0_z_0[i] * pb_y + g_0_yyyzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 72-75 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yyyzzz_0_x_0 = prim_buffer_0_sisp[72];

    auto g_0_yyyzzz_0_y_0 = prim_buffer_0_sisp[73];

    auto g_0_yyyzzz_0_z_0 = prim_buffer_0_sisp[74];

    #pragma omp simd aligned(g_0_yyyz_0_y_0, g_0_yyyz_0_y_1, g_0_yyyzz_0_y_0, g_0_yyyzz_0_y_1, g_0_yyyzzz_0_x_0, g_0_yyyzzz_0_y_0, g_0_yyyzzz_0_z_0, g_0_yyzzz_0_x_0, g_0_yyzzz_0_x_1, g_0_yyzzz_0_z_0, g_0_yyzzz_0_z_1, g_0_yzzz_0_x_0, g_0_yzzz_0_x_1, g_0_yzzz_0_z_0, g_0_yzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_x_0[i] = 2.0 * g_0_yzzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_x_1[i] * fti_ab_0 + g_0_yyzzz_0_x_0[i] * pb_y + g_0_yyzzz_0_x_1[i] * wp_y[i];

        g_0_yyyzzz_0_y_0[i] = 2.0 * g_0_yyyz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_y_1[i] * fti_ab_0 + g_0_yyyzz_0_y_0[i] * pb_z + g_0_yyyzz_0_y_1[i] * wp_z[i];

        g_0_yyyzzz_0_z_0[i] = 2.0 * g_0_yzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_z_1[i] * fti_ab_0 + g_0_yyzzz_0_z_0[i] * pb_y + g_0_yyzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 75-78 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yyzzzz_0_x_0 = prim_buffer_0_sisp[75];

    auto g_0_yyzzzz_0_y_0 = prim_buffer_0_sisp[76];

    auto g_0_yyzzzz_0_z_0 = prim_buffer_0_sisp[77];

    #pragma omp simd aligned(g_0_yyzz_0_y_0, g_0_yyzz_0_y_1, g_0_yyzzz_0_y_0, g_0_yyzzz_0_y_1, g_0_yyzzzz_0_x_0, g_0_yyzzzz_0_y_0, g_0_yyzzzz_0_z_0, g_0_yzzzz_0_x_0, g_0_yzzzz_0_x_1, g_0_yzzzz_0_z_0, g_0_yzzzz_0_z_1, g_0_zzzz_0_x_0, g_0_zzzz_0_x_1, g_0_zzzz_0_z_0, g_0_zzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_x_0[i] = g_0_zzzz_0_x_0[i] * fi_ab_0 - g_0_zzzz_0_x_1[i] * fti_ab_0 + g_0_yzzzz_0_x_0[i] * pb_y + g_0_yzzzz_0_x_1[i] * wp_y[i];

        g_0_yyzzzz_0_y_0[i] = 3.0 * g_0_yyzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_y_1[i] * fti_ab_0 + g_0_yyzzz_0_y_0[i] * pb_z + g_0_yyzzz_0_y_1[i] * wp_z[i];

        g_0_yyzzzz_0_z_0[i] = g_0_zzzz_0_z_0[i] * fi_ab_0 - g_0_zzzz_0_z_1[i] * fti_ab_0 + g_0_yzzzz_0_z_0[i] * pb_y + g_0_yzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 78-81 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_yzzzzz_0_x_0 = prim_buffer_0_sisp[78];

    auto g_0_yzzzzz_0_y_0 = prim_buffer_0_sisp[79];

    auto g_0_yzzzzz_0_z_0 = prim_buffer_0_sisp[80];

    #pragma omp simd aligned(g_0_yzzzzz_0_x_0, g_0_yzzzzz_0_y_0, g_0_yzzzzz_0_z_0, g_0_zzzzz_0_0_1, g_0_zzzzz_0_x_0, g_0_zzzzz_0_x_1, g_0_zzzzz_0_y_0, g_0_zzzzz_0_y_1, g_0_zzzzz_0_z_0, g_0_zzzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_x_0[i] = g_0_zzzzz_0_x_0[i] * pb_y + g_0_zzzzz_0_x_1[i] * wp_y[i];

        g_0_yzzzzz_0_y_0[i] = g_0_zzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzz_0_y_0[i] * pb_y + g_0_zzzzz_0_y_1[i] * wp_y[i];

        g_0_yzzzzz_0_z_0[i] = g_0_zzzzz_0_z_0[i] * pb_y + g_0_zzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 81-84 components of targeted buffer : prim_buffer_0_sisp

    auto g_0_zzzzzz_0_x_0 = prim_buffer_0_sisp[81];

    auto g_0_zzzzzz_0_y_0 = prim_buffer_0_sisp[82];

    auto g_0_zzzzzz_0_z_0 = prim_buffer_0_sisp[83];

    #pragma omp simd aligned(g_0_zzzz_0_x_0, g_0_zzzz_0_x_1, g_0_zzzz_0_y_0, g_0_zzzz_0_y_1, g_0_zzzz_0_z_0, g_0_zzzz_0_z_1, g_0_zzzzz_0_0_1, g_0_zzzzz_0_x_0, g_0_zzzzz_0_x_1, g_0_zzzzz_0_y_0, g_0_zzzzz_0_y_1, g_0_zzzzz_0_z_0, g_0_zzzzz_0_z_1, g_0_zzzzzz_0_x_0, g_0_zzzzzz_0_y_0, g_0_zzzzzz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_x_0[i] = 5.0 * g_0_zzzz_0_x_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_x_1[i] * fti_ab_0 + g_0_zzzzz_0_x_0[i] * pb_z + g_0_zzzzz_0_x_1[i] * wp_z[i];

        g_0_zzzzzz_0_y_0[i] = 5.0 * g_0_zzzz_0_y_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_y_1[i] * fti_ab_0 + g_0_zzzzz_0_y_0[i] * pb_z + g_0_zzzzz_0_y_1[i] * wp_z[i];

        g_0_zzzzzz_0_z_0[i] = 5.0 * g_0_zzzz_0_z_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_z_1[i] * fti_ab_0 + g_0_zzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzz_0_z_0[i] * pb_z + g_0_zzzzz_0_z_1[i] * wp_z[i];
    }
}

} // erirec namespace


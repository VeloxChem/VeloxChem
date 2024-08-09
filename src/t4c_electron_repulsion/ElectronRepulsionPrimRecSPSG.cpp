#include "ElectronRepulsionPrimRecSPSG.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsg(CSimdArray<double>& prim_buffer_0_spsg,
                                  const CSimdArray<double>& prim_buffer_1_sssf,
                                  const CSimdArray<double>& prim_buffer_0_sssg,
                                  const CSimdArray<double>& prim_buffer_1_sssg,
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
    const auto ndims = prim_buffer_0_spsg.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_1_sssf

    auto g_0_0_0_xxx_1 = prim_buffer_1_sssf[0];

    auto g_0_0_0_xxy_1 = prim_buffer_1_sssf[1];

    auto g_0_0_0_xxz_1 = prim_buffer_1_sssf[2];

    auto g_0_0_0_xyy_1 = prim_buffer_1_sssf[3];

    auto g_0_0_0_xyz_1 = prim_buffer_1_sssf[4];

    auto g_0_0_0_xzz_1 = prim_buffer_1_sssf[5];

    auto g_0_0_0_yyy_1 = prim_buffer_1_sssf[6];

    auto g_0_0_0_yyz_1 = prim_buffer_1_sssf[7];

    auto g_0_0_0_yzz_1 = prim_buffer_1_sssf[8];

    auto g_0_0_0_zzz_1 = prim_buffer_1_sssf[9];

    /// Set up components of auxilary buffer : prim_buffer_0_sssg

    auto g_0_0_0_xxxx_0 = prim_buffer_0_sssg[0];

    auto g_0_0_0_xxxy_0 = prim_buffer_0_sssg[1];

    auto g_0_0_0_xxxz_0 = prim_buffer_0_sssg[2];

    auto g_0_0_0_xxyy_0 = prim_buffer_0_sssg[3];

    auto g_0_0_0_xxyz_0 = prim_buffer_0_sssg[4];

    auto g_0_0_0_xxzz_0 = prim_buffer_0_sssg[5];

    auto g_0_0_0_xyyy_0 = prim_buffer_0_sssg[6];

    auto g_0_0_0_xyyz_0 = prim_buffer_0_sssg[7];

    auto g_0_0_0_xyzz_0 = prim_buffer_0_sssg[8];

    auto g_0_0_0_xzzz_0 = prim_buffer_0_sssg[9];

    auto g_0_0_0_yyyy_0 = prim_buffer_0_sssg[10];

    auto g_0_0_0_yyyz_0 = prim_buffer_0_sssg[11];

    auto g_0_0_0_yyzz_0 = prim_buffer_0_sssg[12];

    auto g_0_0_0_yzzz_0 = prim_buffer_0_sssg[13];

    auto g_0_0_0_zzzz_0 = prim_buffer_0_sssg[14];

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

    /// Set up 0-15 components of targeted buffer : prim_buffer_0_spsg

    auto g_0_x_0_xxxx_0 = prim_buffer_0_spsg[0];

    auto g_0_x_0_xxxy_0 = prim_buffer_0_spsg[1];

    auto g_0_x_0_xxxz_0 = prim_buffer_0_spsg[2];

    auto g_0_x_0_xxyy_0 = prim_buffer_0_spsg[3];

    auto g_0_x_0_xxyz_0 = prim_buffer_0_spsg[4];

    auto g_0_x_0_xxzz_0 = prim_buffer_0_spsg[5];

    auto g_0_x_0_xyyy_0 = prim_buffer_0_spsg[6];

    auto g_0_x_0_xyyz_0 = prim_buffer_0_spsg[7];

    auto g_0_x_0_xyzz_0 = prim_buffer_0_spsg[8];

    auto g_0_x_0_xzzz_0 = prim_buffer_0_spsg[9];

    auto g_0_x_0_yyyy_0 = prim_buffer_0_spsg[10];

    auto g_0_x_0_yyyz_0 = prim_buffer_0_spsg[11];

    auto g_0_x_0_yyzz_0 = prim_buffer_0_spsg[12];

    auto g_0_x_0_yzzz_0 = prim_buffer_0_spsg[13];

    auto g_0_x_0_zzzz_0 = prim_buffer_0_spsg[14];

    #pragma omp simd aligned(g_0_0_0_xxx_1, g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxy_0, g_0_0_0_xxxy_1, g_0_0_0_xxxz_0, g_0_0_0_xxxz_1, g_0_0_0_xxy_1, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyz_0, g_0_0_0_xxyz_1, g_0_0_0_xxz_1, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xyy_1, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyz_0, g_0_0_0_xyyz_1, g_0_0_0_xyz_1, g_0_0_0_xyzz_0, g_0_0_0_xyzz_1, g_0_0_0_xzz_1, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_yyy_1, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyz_0, g_0_0_0_yyyz_1, g_0_0_0_yyz_1, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yzz_1, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_zzz_1, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_x_0_xxxx_0, g_0_x_0_xxxy_0, g_0_x_0_xxxz_0, g_0_x_0_xxyy_0, g_0_x_0_xxyz_0, g_0_x_0_xxzz_0, g_0_x_0_xyyy_0, g_0_x_0_xyyz_0, g_0_x_0_xyzz_0, g_0_x_0_xzzz_0, g_0_x_0_yyyy_0, g_0_x_0_yyyz_0, g_0_x_0_yyzz_0, g_0_x_0_yzzz_0, g_0_x_0_zzzz_0, wp_x  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxxx_0[i] = 4.0 * g_0_0_0_xxx_1[i] * fi_abcd_0 + g_0_0_0_xxxx_0[i] * pb_x + g_0_0_0_xxxx_1[i] * wp_x[i];

        g_0_x_0_xxxy_0[i] = 3.0 * g_0_0_0_xxy_1[i] * fi_abcd_0 + g_0_0_0_xxxy_0[i] * pb_x + g_0_0_0_xxxy_1[i] * wp_x[i];

        g_0_x_0_xxxz_0[i] = 3.0 * g_0_0_0_xxz_1[i] * fi_abcd_0 + g_0_0_0_xxxz_0[i] * pb_x + g_0_0_0_xxxz_1[i] * wp_x[i];

        g_0_x_0_xxyy_0[i] = 2.0 * g_0_0_0_xyy_1[i] * fi_abcd_0 + g_0_0_0_xxyy_0[i] * pb_x + g_0_0_0_xxyy_1[i] * wp_x[i];

        g_0_x_0_xxyz_0[i] = 2.0 * g_0_0_0_xyz_1[i] * fi_abcd_0 + g_0_0_0_xxyz_0[i] * pb_x + g_0_0_0_xxyz_1[i] * wp_x[i];

        g_0_x_0_xxzz_0[i] = 2.0 * g_0_0_0_xzz_1[i] * fi_abcd_0 + g_0_0_0_xxzz_0[i] * pb_x + g_0_0_0_xxzz_1[i] * wp_x[i];

        g_0_x_0_xyyy_0[i] = g_0_0_0_yyy_1[i] * fi_abcd_0 + g_0_0_0_xyyy_0[i] * pb_x + g_0_0_0_xyyy_1[i] * wp_x[i];

        g_0_x_0_xyyz_0[i] = g_0_0_0_yyz_1[i] * fi_abcd_0 + g_0_0_0_xyyz_0[i] * pb_x + g_0_0_0_xyyz_1[i] * wp_x[i];

        g_0_x_0_xyzz_0[i] = g_0_0_0_yzz_1[i] * fi_abcd_0 + g_0_0_0_xyzz_0[i] * pb_x + g_0_0_0_xyzz_1[i] * wp_x[i];

        g_0_x_0_xzzz_0[i] = g_0_0_0_zzz_1[i] * fi_abcd_0 + g_0_0_0_xzzz_0[i] * pb_x + g_0_0_0_xzzz_1[i] * wp_x[i];

        g_0_x_0_yyyy_0[i] = g_0_0_0_yyyy_0[i] * pb_x + g_0_0_0_yyyy_1[i] * wp_x[i];

        g_0_x_0_yyyz_0[i] = g_0_0_0_yyyz_0[i] * pb_x + g_0_0_0_yyyz_1[i] * wp_x[i];

        g_0_x_0_yyzz_0[i] = g_0_0_0_yyzz_0[i] * pb_x + g_0_0_0_yyzz_1[i] * wp_x[i];

        g_0_x_0_yzzz_0[i] = g_0_0_0_yzzz_0[i] * pb_x + g_0_0_0_yzzz_1[i] * wp_x[i];

        g_0_x_0_zzzz_0[i] = g_0_0_0_zzzz_0[i] * pb_x + g_0_0_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : prim_buffer_0_spsg

    auto g_0_y_0_xxxx_0 = prim_buffer_0_spsg[15];

    auto g_0_y_0_xxxy_0 = prim_buffer_0_spsg[16];

    auto g_0_y_0_xxxz_0 = prim_buffer_0_spsg[17];

    auto g_0_y_0_xxyy_0 = prim_buffer_0_spsg[18];

    auto g_0_y_0_xxyz_0 = prim_buffer_0_spsg[19];

    auto g_0_y_0_xxzz_0 = prim_buffer_0_spsg[20];

    auto g_0_y_0_xyyy_0 = prim_buffer_0_spsg[21];

    auto g_0_y_0_xyyz_0 = prim_buffer_0_spsg[22];

    auto g_0_y_0_xyzz_0 = prim_buffer_0_spsg[23];

    auto g_0_y_0_xzzz_0 = prim_buffer_0_spsg[24];

    auto g_0_y_0_yyyy_0 = prim_buffer_0_spsg[25];

    auto g_0_y_0_yyyz_0 = prim_buffer_0_spsg[26];

    auto g_0_y_0_yyzz_0 = prim_buffer_0_spsg[27];

    auto g_0_y_0_yzzz_0 = prim_buffer_0_spsg[28];

    auto g_0_y_0_zzzz_0 = prim_buffer_0_spsg[29];

    #pragma omp simd aligned(g_0_0_0_xxx_1, g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxy_0, g_0_0_0_xxxy_1, g_0_0_0_xxxz_0, g_0_0_0_xxxz_1, g_0_0_0_xxy_1, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyz_0, g_0_0_0_xxyz_1, g_0_0_0_xxz_1, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xyy_1, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyz_0, g_0_0_0_xyyz_1, g_0_0_0_xyz_1, g_0_0_0_xyzz_0, g_0_0_0_xyzz_1, g_0_0_0_xzz_1, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_yyy_1, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyz_0, g_0_0_0_yyyz_1, g_0_0_0_yyz_1, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yzz_1, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_zzz_1, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_y_0_xxxx_0, g_0_y_0_xxxy_0, g_0_y_0_xxxz_0, g_0_y_0_xxyy_0, g_0_y_0_xxyz_0, g_0_y_0_xxzz_0, g_0_y_0_xyyy_0, g_0_y_0_xyyz_0, g_0_y_0_xyzz_0, g_0_y_0_xzzz_0, g_0_y_0_yyyy_0, g_0_y_0_yyyz_0, g_0_y_0_yyzz_0, g_0_y_0_yzzz_0, g_0_y_0_zzzz_0, wp_y  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxxx_0[i] = g_0_0_0_xxxx_0[i] * pb_y + g_0_0_0_xxxx_1[i] * wp_y[i];

        g_0_y_0_xxxy_0[i] = g_0_0_0_xxx_1[i] * fi_abcd_0 + g_0_0_0_xxxy_0[i] * pb_y + g_0_0_0_xxxy_1[i] * wp_y[i];

        g_0_y_0_xxxz_0[i] = g_0_0_0_xxxz_0[i] * pb_y + g_0_0_0_xxxz_1[i] * wp_y[i];

        g_0_y_0_xxyy_0[i] = 2.0 * g_0_0_0_xxy_1[i] * fi_abcd_0 + g_0_0_0_xxyy_0[i] * pb_y + g_0_0_0_xxyy_1[i] * wp_y[i];

        g_0_y_0_xxyz_0[i] = g_0_0_0_xxz_1[i] * fi_abcd_0 + g_0_0_0_xxyz_0[i] * pb_y + g_0_0_0_xxyz_1[i] * wp_y[i];

        g_0_y_0_xxzz_0[i] = g_0_0_0_xxzz_0[i] * pb_y + g_0_0_0_xxzz_1[i] * wp_y[i];

        g_0_y_0_xyyy_0[i] = 3.0 * g_0_0_0_xyy_1[i] * fi_abcd_0 + g_0_0_0_xyyy_0[i] * pb_y + g_0_0_0_xyyy_1[i] * wp_y[i];

        g_0_y_0_xyyz_0[i] = 2.0 * g_0_0_0_xyz_1[i] * fi_abcd_0 + g_0_0_0_xyyz_0[i] * pb_y + g_0_0_0_xyyz_1[i] * wp_y[i];

        g_0_y_0_xyzz_0[i] = g_0_0_0_xzz_1[i] * fi_abcd_0 + g_0_0_0_xyzz_0[i] * pb_y + g_0_0_0_xyzz_1[i] * wp_y[i];

        g_0_y_0_xzzz_0[i] = g_0_0_0_xzzz_0[i] * pb_y + g_0_0_0_xzzz_1[i] * wp_y[i];

        g_0_y_0_yyyy_0[i] = 4.0 * g_0_0_0_yyy_1[i] * fi_abcd_0 + g_0_0_0_yyyy_0[i] * pb_y + g_0_0_0_yyyy_1[i] * wp_y[i];

        g_0_y_0_yyyz_0[i] = 3.0 * g_0_0_0_yyz_1[i] * fi_abcd_0 + g_0_0_0_yyyz_0[i] * pb_y + g_0_0_0_yyyz_1[i] * wp_y[i];

        g_0_y_0_yyzz_0[i] = 2.0 * g_0_0_0_yzz_1[i] * fi_abcd_0 + g_0_0_0_yyzz_0[i] * pb_y + g_0_0_0_yyzz_1[i] * wp_y[i];

        g_0_y_0_yzzz_0[i] = g_0_0_0_zzz_1[i] * fi_abcd_0 + g_0_0_0_yzzz_0[i] * pb_y + g_0_0_0_yzzz_1[i] * wp_y[i];

        g_0_y_0_zzzz_0[i] = g_0_0_0_zzzz_0[i] * pb_y + g_0_0_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : prim_buffer_0_spsg

    auto g_0_z_0_xxxx_0 = prim_buffer_0_spsg[30];

    auto g_0_z_0_xxxy_0 = prim_buffer_0_spsg[31];

    auto g_0_z_0_xxxz_0 = prim_buffer_0_spsg[32];

    auto g_0_z_0_xxyy_0 = prim_buffer_0_spsg[33];

    auto g_0_z_0_xxyz_0 = prim_buffer_0_spsg[34];

    auto g_0_z_0_xxzz_0 = prim_buffer_0_spsg[35];

    auto g_0_z_0_xyyy_0 = prim_buffer_0_spsg[36];

    auto g_0_z_0_xyyz_0 = prim_buffer_0_spsg[37];

    auto g_0_z_0_xyzz_0 = prim_buffer_0_spsg[38];

    auto g_0_z_0_xzzz_0 = prim_buffer_0_spsg[39];

    auto g_0_z_0_yyyy_0 = prim_buffer_0_spsg[40];

    auto g_0_z_0_yyyz_0 = prim_buffer_0_spsg[41];

    auto g_0_z_0_yyzz_0 = prim_buffer_0_spsg[42];

    auto g_0_z_0_yzzz_0 = prim_buffer_0_spsg[43];

    auto g_0_z_0_zzzz_0 = prim_buffer_0_spsg[44];

    #pragma omp simd aligned(g_0_0_0_xxx_1, g_0_0_0_xxxx_0, g_0_0_0_xxxx_1, g_0_0_0_xxxy_0, g_0_0_0_xxxy_1, g_0_0_0_xxxz_0, g_0_0_0_xxxz_1, g_0_0_0_xxy_1, g_0_0_0_xxyy_0, g_0_0_0_xxyy_1, g_0_0_0_xxyz_0, g_0_0_0_xxyz_1, g_0_0_0_xxz_1, g_0_0_0_xxzz_0, g_0_0_0_xxzz_1, g_0_0_0_xyy_1, g_0_0_0_xyyy_0, g_0_0_0_xyyy_1, g_0_0_0_xyyz_0, g_0_0_0_xyyz_1, g_0_0_0_xyz_1, g_0_0_0_xyzz_0, g_0_0_0_xyzz_1, g_0_0_0_xzz_1, g_0_0_0_xzzz_0, g_0_0_0_xzzz_1, g_0_0_0_yyy_1, g_0_0_0_yyyy_0, g_0_0_0_yyyy_1, g_0_0_0_yyyz_0, g_0_0_0_yyyz_1, g_0_0_0_yyz_1, g_0_0_0_yyzz_0, g_0_0_0_yyzz_1, g_0_0_0_yzz_1, g_0_0_0_yzzz_0, g_0_0_0_yzzz_1, g_0_0_0_zzz_1, g_0_0_0_zzzz_0, g_0_0_0_zzzz_1, g_0_z_0_xxxx_0, g_0_z_0_xxxy_0, g_0_z_0_xxxz_0, g_0_z_0_xxyy_0, g_0_z_0_xxyz_0, g_0_z_0_xxzz_0, g_0_z_0_xyyy_0, g_0_z_0_xyyz_0, g_0_z_0_xyzz_0, g_0_z_0_xzzz_0, g_0_z_0_yyyy_0, g_0_z_0_yyyz_0, g_0_z_0_yyzz_0, g_0_z_0_yzzz_0, g_0_z_0_zzzz_0, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxxx_0[i] = g_0_0_0_xxxx_0[i] * pb_z + g_0_0_0_xxxx_1[i] * wp_z[i];

        g_0_z_0_xxxy_0[i] = g_0_0_0_xxxy_0[i] * pb_z + g_0_0_0_xxxy_1[i] * wp_z[i];

        g_0_z_0_xxxz_0[i] = g_0_0_0_xxx_1[i] * fi_abcd_0 + g_0_0_0_xxxz_0[i] * pb_z + g_0_0_0_xxxz_1[i] * wp_z[i];

        g_0_z_0_xxyy_0[i] = g_0_0_0_xxyy_0[i] * pb_z + g_0_0_0_xxyy_1[i] * wp_z[i];

        g_0_z_0_xxyz_0[i] = g_0_0_0_xxy_1[i] * fi_abcd_0 + g_0_0_0_xxyz_0[i] * pb_z + g_0_0_0_xxyz_1[i] * wp_z[i];

        g_0_z_0_xxzz_0[i] = 2.0 * g_0_0_0_xxz_1[i] * fi_abcd_0 + g_0_0_0_xxzz_0[i] * pb_z + g_0_0_0_xxzz_1[i] * wp_z[i];

        g_0_z_0_xyyy_0[i] = g_0_0_0_xyyy_0[i] * pb_z + g_0_0_0_xyyy_1[i] * wp_z[i];

        g_0_z_0_xyyz_0[i] = g_0_0_0_xyy_1[i] * fi_abcd_0 + g_0_0_0_xyyz_0[i] * pb_z + g_0_0_0_xyyz_1[i] * wp_z[i];

        g_0_z_0_xyzz_0[i] = 2.0 * g_0_0_0_xyz_1[i] * fi_abcd_0 + g_0_0_0_xyzz_0[i] * pb_z + g_0_0_0_xyzz_1[i] * wp_z[i];

        g_0_z_0_xzzz_0[i] = 3.0 * g_0_0_0_xzz_1[i] * fi_abcd_0 + g_0_0_0_xzzz_0[i] * pb_z + g_0_0_0_xzzz_1[i] * wp_z[i];

        g_0_z_0_yyyy_0[i] = g_0_0_0_yyyy_0[i] * pb_z + g_0_0_0_yyyy_1[i] * wp_z[i];

        g_0_z_0_yyyz_0[i] = g_0_0_0_yyy_1[i] * fi_abcd_0 + g_0_0_0_yyyz_0[i] * pb_z + g_0_0_0_yyyz_1[i] * wp_z[i];

        g_0_z_0_yyzz_0[i] = 2.0 * g_0_0_0_yyz_1[i] * fi_abcd_0 + g_0_0_0_yyzz_0[i] * pb_z + g_0_0_0_yyzz_1[i] * wp_z[i];

        g_0_z_0_yzzz_0[i] = 3.0 * g_0_0_0_yzz_1[i] * fi_abcd_0 + g_0_0_0_yzzz_0[i] * pb_z + g_0_0_0_yzzz_1[i] * wp_z[i];

        g_0_z_0_zzzz_0[i] = 4.0 * g_0_0_0_zzz_1[i] * fi_abcd_0 + g_0_0_0_zzzz_0[i] * pb_z + g_0_0_0_zzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


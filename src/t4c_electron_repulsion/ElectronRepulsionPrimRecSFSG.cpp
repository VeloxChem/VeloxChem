#include "ElectronRepulsionPrimRecSFSG.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfsg(CSimdArray<double>& prim_buffer_0_sfsg,
                                  const CSimdArray<double>& prim_buffer_0_spsg,
                                  const CSimdArray<double>& prim_buffer_1_spsg,
                                  const CSimdArray<double>& prim_buffer_1_sdsf,
                                  const CSimdArray<double>& prim_buffer_0_sdsg,
                                  const CSimdArray<double>& prim_buffer_1_sdsg,
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
    const auto ndims = prim_buffer_0_sfsg.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_spsg

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

    /// Set up components of auxilary buffer : prim_buffer_1_spsg

    auto g_0_x_0_xxxx_1 = prim_buffer_1_spsg[0];

    auto g_0_x_0_xxxy_1 = prim_buffer_1_spsg[1];

    auto g_0_x_0_xxxz_1 = prim_buffer_1_spsg[2];

    auto g_0_x_0_xxyy_1 = prim_buffer_1_spsg[3];

    auto g_0_x_0_xxyz_1 = prim_buffer_1_spsg[4];

    auto g_0_x_0_xxzz_1 = prim_buffer_1_spsg[5];

    auto g_0_x_0_xyyy_1 = prim_buffer_1_spsg[6];

    auto g_0_x_0_xyyz_1 = prim_buffer_1_spsg[7];

    auto g_0_x_0_xyzz_1 = prim_buffer_1_spsg[8];

    auto g_0_x_0_xzzz_1 = prim_buffer_1_spsg[9];

    auto g_0_x_0_yyyy_1 = prim_buffer_1_spsg[10];

    auto g_0_x_0_yyyz_1 = prim_buffer_1_spsg[11];

    auto g_0_x_0_yyzz_1 = prim_buffer_1_spsg[12];

    auto g_0_x_0_yzzz_1 = prim_buffer_1_spsg[13];

    auto g_0_x_0_zzzz_1 = prim_buffer_1_spsg[14];

    auto g_0_y_0_xxxx_1 = prim_buffer_1_spsg[15];

    auto g_0_y_0_xxxy_1 = prim_buffer_1_spsg[16];

    auto g_0_y_0_xxxz_1 = prim_buffer_1_spsg[17];

    auto g_0_y_0_xxyy_1 = prim_buffer_1_spsg[18];

    auto g_0_y_0_xxyz_1 = prim_buffer_1_spsg[19];

    auto g_0_y_0_xxzz_1 = prim_buffer_1_spsg[20];

    auto g_0_y_0_xyyy_1 = prim_buffer_1_spsg[21];

    auto g_0_y_0_xyyz_1 = prim_buffer_1_spsg[22];

    auto g_0_y_0_xyzz_1 = prim_buffer_1_spsg[23];

    auto g_0_y_0_xzzz_1 = prim_buffer_1_spsg[24];

    auto g_0_y_0_yyyy_1 = prim_buffer_1_spsg[25];

    auto g_0_y_0_yyyz_1 = prim_buffer_1_spsg[26];

    auto g_0_y_0_yyzz_1 = prim_buffer_1_spsg[27];

    auto g_0_y_0_yzzz_1 = prim_buffer_1_spsg[28];

    auto g_0_y_0_zzzz_1 = prim_buffer_1_spsg[29];

    auto g_0_z_0_xxxx_1 = prim_buffer_1_spsg[30];

    auto g_0_z_0_xxxy_1 = prim_buffer_1_spsg[31];

    auto g_0_z_0_xxxz_1 = prim_buffer_1_spsg[32];

    auto g_0_z_0_xxyy_1 = prim_buffer_1_spsg[33];

    auto g_0_z_0_xxyz_1 = prim_buffer_1_spsg[34];

    auto g_0_z_0_xxzz_1 = prim_buffer_1_spsg[35];

    auto g_0_z_0_xyyy_1 = prim_buffer_1_spsg[36];

    auto g_0_z_0_xyyz_1 = prim_buffer_1_spsg[37];

    auto g_0_z_0_xyzz_1 = prim_buffer_1_spsg[38];

    auto g_0_z_0_xzzz_1 = prim_buffer_1_spsg[39];

    auto g_0_z_0_yyyy_1 = prim_buffer_1_spsg[40];

    auto g_0_z_0_yyyz_1 = prim_buffer_1_spsg[41];

    auto g_0_z_0_yyzz_1 = prim_buffer_1_spsg[42];

    auto g_0_z_0_yzzz_1 = prim_buffer_1_spsg[43];

    auto g_0_z_0_zzzz_1 = prim_buffer_1_spsg[44];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsf

    auto g_0_xx_0_xxx_1 = prim_buffer_1_sdsf[0];

    auto g_0_xx_0_xxy_1 = prim_buffer_1_sdsf[1];

    auto g_0_xx_0_xxz_1 = prim_buffer_1_sdsf[2];

    auto g_0_xx_0_xyy_1 = prim_buffer_1_sdsf[3];

    auto g_0_xx_0_xyz_1 = prim_buffer_1_sdsf[4];

    auto g_0_xx_0_xzz_1 = prim_buffer_1_sdsf[5];

    auto g_0_xx_0_yyy_1 = prim_buffer_1_sdsf[6];

    auto g_0_xx_0_yyz_1 = prim_buffer_1_sdsf[7];

    auto g_0_xx_0_yzz_1 = prim_buffer_1_sdsf[8];

    auto g_0_xx_0_zzz_1 = prim_buffer_1_sdsf[9];

    auto g_0_yy_0_xxx_1 = prim_buffer_1_sdsf[30];

    auto g_0_yy_0_xxy_1 = prim_buffer_1_sdsf[31];

    auto g_0_yy_0_xxz_1 = prim_buffer_1_sdsf[32];

    auto g_0_yy_0_xyy_1 = prim_buffer_1_sdsf[33];

    auto g_0_yy_0_xyz_1 = prim_buffer_1_sdsf[34];

    auto g_0_yy_0_xzz_1 = prim_buffer_1_sdsf[35];

    auto g_0_yy_0_yyy_1 = prim_buffer_1_sdsf[36];

    auto g_0_yy_0_yyz_1 = prim_buffer_1_sdsf[37];

    auto g_0_yy_0_yzz_1 = prim_buffer_1_sdsf[38];

    auto g_0_yy_0_zzz_1 = prim_buffer_1_sdsf[39];

    auto g_0_yz_0_xyz_1 = prim_buffer_1_sdsf[44];

    auto g_0_yz_0_yyz_1 = prim_buffer_1_sdsf[47];

    auto g_0_yz_0_yzz_1 = prim_buffer_1_sdsf[48];

    auto g_0_zz_0_xxx_1 = prim_buffer_1_sdsf[50];

    auto g_0_zz_0_xxy_1 = prim_buffer_1_sdsf[51];

    auto g_0_zz_0_xxz_1 = prim_buffer_1_sdsf[52];

    auto g_0_zz_0_xyy_1 = prim_buffer_1_sdsf[53];

    auto g_0_zz_0_xyz_1 = prim_buffer_1_sdsf[54];

    auto g_0_zz_0_xzz_1 = prim_buffer_1_sdsf[55];

    auto g_0_zz_0_yyy_1 = prim_buffer_1_sdsf[56];

    auto g_0_zz_0_yyz_1 = prim_buffer_1_sdsf[57];

    auto g_0_zz_0_yzz_1 = prim_buffer_1_sdsf[58];

    auto g_0_zz_0_zzz_1 = prim_buffer_1_sdsf[59];

    /// Set up components of auxilary buffer : prim_buffer_0_sdsg

    auto g_0_xx_0_xxxx_0 = prim_buffer_0_sdsg[0];

    auto g_0_xx_0_xxxy_0 = prim_buffer_0_sdsg[1];

    auto g_0_xx_0_xxxz_0 = prim_buffer_0_sdsg[2];

    auto g_0_xx_0_xxyy_0 = prim_buffer_0_sdsg[3];

    auto g_0_xx_0_xxyz_0 = prim_buffer_0_sdsg[4];

    auto g_0_xx_0_xxzz_0 = prim_buffer_0_sdsg[5];

    auto g_0_xx_0_xyyy_0 = prim_buffer_0_sdsg[6];

    auto g_0_xx_0_xyyz_0 = prim_buffer_0_sdsg[7];

    auto g_0_xx_0_xyzz_0 = prim_buffer_0_sdsg[8];

    auto g_0_xx_0_xzzz_0 = prim_buffer_0_sdsg[9];

    auto g_0_xx_0_yyyy_0 = prim_buffer_0_sdsg[10];

    auto g_0_xx_0_yyyz_0 = prim_buffer_0_sdsg[11];

    auto g_0_xx_0_yyzz_0 = prim_buffer_0_sdsg[12];

    auto g_0_xx_0_yzzz_0 = prim_buffer_0_sdsg[13];

    auto g_0_xx_0_zzzz_0 = prim_buffer_0_sdsg[14];

    auto g_0_xy_0_xxxy_0 = prim_buffer_0_sdsg[16];

    auto g_0_xy_0_xxyy_0 = prim_buffer_0_sdsg[18];

    auto g_0_xy_0_xyyy_0 = prim_buffer_0_sdsg[21];

    auto g_0_xz_0_xxxx_0 = prim_buffer_0_sdsg[30];

    auto g_0_xz_0_xxxz_0 = prim_buffer_0_sdsg[32];

    auto g_0_xz_0_xxzz_0 = prim_buffer_0_sdsg[35];

    auto g_0_xz_0_xzzz_0 = prim_buffer_0_sdsg[39];

    auto g_0_yy_0_xxxx_0 = prim_buffer_0_sdsg[45];

    auto g_0_yy_0_xxxy_0 = prim_buffer_0_sdsg[46];

    auto g_0_yy_0_xxxz_0 = prim_buffer_0_sdsg[47];

    auto g_0_yy_0_xxyy_0 = prim_buffer_0_sdsg[48];

    auto g_0_yy_0_xxyz_0 = prim_buffer_0_sdsg[49];

    auto g_0_yy_0_xxzz_0 = prim_buffer_0_sdsg[50];

    auto g_0_yy_0_xyyy_0 = prim_buffer_0_sdsg[51];

    auto g_0_yy_0_xyyz_0 = prim_buffer_0_sdsg[52];

    auto g_0_yy_0_xyzz_0 = prim_buffer_0_sdsg[53];

    auto g_0_yy_0_xzzz_0 = prim_buffer_0_sdsg[54];

    auto g_0_yy_0_yyyy_0 = prim_buffer_0_sdsg[55];

    auto g_0_yy_0_yyyz_0 = prim_buffer_0_sdsg[56];

    auto g_0_yy_0_yyzz_0 = prim_buffer_0_sdsg[57];

    auto g_0_yy_0_yzzz_0 = prim_buffer_0_sdsg[58];

    auto g_0_yy_0_zzzz_0 = prim_buffer_0_sdsg[59];

    auto g_0_yz_0_xxyz_0 = prim_buffer_0_sdsg[64];

    auto g_0_yz_0_xyyz_0 = prim_buffer_0_sdsg[67];

    auto g_0_yz_0_xyzz_0 = prim_buffer_0_sdsg[68];

    auto g_0_yz_0_yyyy_0 = prim_buffer_0_sdsg[70];

    auto g_0_yz_0_yyyz_0 = prim_buffer_0_sdsg[71];

    auto g_0_yz_0_yyzz_0 = prim_buffer_0_sdsg[72];

    auto g_0_yz_0_yzzz_0 = prim_buffer_0_sdsg[73];

    auto g_0_yz_0_zzzz_0 = prim_buffer_0_sdsg[74];

    auto g_0_zz_0_xxxx_0 = prim_buffer_0_sdsg[75];

    auto g_0_zz_0_xxxy_0 = prim_buffer_0_sdsg[76];

    auto g_0_zz_0_xxxz_0 = prim_buffer_0_sdsg[77];

    auto g_0_zz_0_xxyy_0 = prim_buffer_0_sdsg[78];

    auto g_0_zz_0_xxyz_0 = prim_buffer_0_sdsg[79];

    auto g_0_zz_0_xxzz_0 = prim_buffer_0_sdsg[80];

    auto g_0_zz_0_xyyy_0 = prim_buffer_0_sdsg[81];

    auto g_0_zz_0_xyyz_0 = prim_buffer_0_sdsg[82];

    auto g_0_zz_0_xyzz_0 = prim_buffer_0_sdsg[83];

    auto g_0_zz_0_xzzz_0 = prim_buffer_0_sdsg[84];

    auto g_0_zz_0_yyyy_0 = prim_buffer_0_sdsg[85];

    auto g_0_zz_0_yyyz_0 = prim_buffer_0_sdsg[86];

    auto g_0_zz_0_yyzz_0 = prim_buffer_0_sdsg[87];

    auto g_0_zz_0_yzzz_0 = prim_buffer_0_sdsg[88];

    auto g_0_zz_0_zzzz_0 = prim_buffer_0_sdsg[89];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsg

    auto g_0_xx_0_xxxx_1 = prim_buffer_1_sdsg[0];

    auto g_0_xx_0_xxxy_1 = prim_buffer_1_sdsg[1];

    auto g_0_xx_0_xxxz_1 = prim_buffer_1_sdsg[2];

    auto g_0_xx_0_xxyy_1 = prim_buffer_1_sdsg[3];

    auto g_0_xx_0_xxyz_1 = prim_buffer_1_sdsg[4];

    auto g_0_xx_0_xxzz_1 = prim_buffer_1_sdsg[5];

    auto g_0_xx_0_xyyy_1 = prim_buffer_1_sdsg[6];

    auto g_0_xx_0_xyyz_1 = prim_buffer_1_sdsg[7];

    auto g_0_xx_0_xyzz_1 = prim_buffer_1_sdsg[8];

    auto g_0_xx_0_xzzz_1 = prim_buffer_1_sdsg[9];

    auto g_0_xx_0_yyyy_1 = prim_buffer_1_sdsg[10];

    auto g_0_xx_0_yyyz_1 = prim_buffer_1_sdsg[11];

    auto g_0_xx_0_yyzz_1 = prim_buffer_1_sdsg[12];

    auto g_0_xx_0_yzzz_1 = prim_buffer_1_sdsg[13];

    auto g_0_xx_0_zzzz_1 = prim_buffer_1_sdsg[14];

    auto g_0_xy_0_xxxy_1 = prim_buffer_1_sdsg[16];

    auto g_0_xy_0_xxyy_1 = prim_buffer_1_sdsg[18];

    auto g_0_xy_0_xyyy_1 = prim_buffer_1_sdsg[21];

    auto g_0_xz_0_xxxx_1 = prim_buffer_1_sdsg[30];

    auto g_0_xz_0_xxxz_1 = prim_buffer_1_sdsg[32];

    auto g_0_xz_0_xxzz_1 = prim_buffer_1_sdsg[35];

    auto g_0_xz_0_xzzz_1 = prim_buffer_1_sdsg[39];

    auto g_0_yy_0_xxxx_1 = prim_buffer_1_sdsg[45];

    auto g_0_yy_0_xxxy_1 = prim_buffer_1_sdsg[46];

    auto g_0_yy_0_xxxz_1 = prim_buffer_1_sdsg[47];

    auto g_0_yy_0_xxyy_1 = prim_buffer_1_sdsg[48];

    auto g_0_yy_0_xxyz_1 = prim_buffer_1_sdsg[49];

    auto g_0_yy_0_xxzz_1 = prim_buffer_1_sdsg[50];

    auto g_0_yy_0_xyyy_1 = prim_buffer_1_sdsg[51];

    auto g_0_yy_0_xyyz_1 = prim_buffer_1_sdsg[52];

    auto g_0_yy_0_xyzz_1 = prim_buffer_1_sdsg[53];

    auto g_0_yy_0_xzzz_1 = prim_buffer_1_sdsg[54];

    auto g_0_yy_0_yyyy_1 = prim_buffer_1_sdsg[55];

    auto g_0_yy_0_yyyz_1 = prim_buffer_1_sdsg[56];

    auto g_0_yy_0_yyzz_1 = prim_buffer_1_sdsg[57];

    auto g_0_yy_0_yzzz_1 = prim_buffer_1_sdsg[58];

    auto g_0_yy_0_zzzz_1 = prim_buffer_1_sdsg[59];

    auto g_0_yz_0_xxyz_1 = prim_buffer_1_sdsg[64];

    auto g_0_yz_0_xyyz_1 = prim_buffer_1_sdsg[67];

    auto g_0_yz_0_xyzz_1 = prim_buffer_1_sdsg[68];

    auto g_0_yz_0_yyyy_1 = prim_buffer_1_sdsg[70];

    auto g_0_yz_0_yyyz_1 = prim_buffer_1_sdsg[71];

    auto g_0_yz_0_yyzz_1 = prim_buffer_1_sdsg[72];

    auto g_0_yz_0_yzzz_1 = prim_buffer_1_sdsg[73];

    auto g_0_yz_0_zzzz_1 = prim_buffer_1_sdsg[74];

    auto g_0_zz_0_xxxx_1 = prim_buffer_1_sdsg[75];

    auto g_0_zz_0_xxxy_1 = prim_buffer_1_sdsg[76];

    auto g_0_zz_0_xxxz_1 = prim_buffer_1_sdsg[77];

    auto g_0_zz_0_xxyy_1 = prim_buffer_1_sdsg[78];

    auto g_0_zz_0_xxyz_1 = prim_buffer_1_sdsg[79];

    auto g_0_zz_0_xxzz_1 = prim_buffer_1_sdsg[80];

    auto g_0_zz_0_xyyy_1 = prim_buffer_1_sdsg[81];

    auto g_0_zz_0_xyyz_1 = prim_buffer_1_sdsg[82];

    auto g_0_zz_0_xyzz_1 = prim_buffer_1_sdsg[83];

    auto g_0_zz_0_xzzz_1 = prim_buffer_1_sdsg[84];

    auto g_0_zz_0_yyyy_1 = prim_buffer_1_sdsg[85];

    auto g_0_zz_0_yyyz_1 = prim_buffer_1_sdsg[86];

    auto g_0_zz_0_yyzz_1 = prim_buffer_1_sdsg[87];

    auto g_0_zz_0_yzzz_1 = prim_buffer_1_sdsg[88];

    auto g_0_zz_0_zzzz_1 = prim_buffer_1_sdsg[89];

    /// Set up 0-15 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xxx_0_xxxx_0 = prim_buffer_0_sfsg[0];

    auto g_0_xxx_0_xxxy_0 = prim_buffer_0_sfsg[1];

    auto g_0_xxx_0_xxxz_0 = prim_buffer_0_sfsg[2];

    auto g_0_xxx_0_xxyy_0 = prim_buffer_0_sfsg[3];

    auto g_0_xxx_0_xxyz_0 = prim_buffer_0_sfsg[4];

    auto g_0_xxx_0_xxzz_0 = prim_buffer_0_sfsg[5];

    auto g_0_xxx_0_xyyy_0 = prim_buffer_0_sfsg[6];

    auto g_0_xxx_0_xyyz_0 = prim_buffer_0_sfsg[7];

    auto g_0_xxx_0_xyzz_0 = prim_buffer_0_sfsg[8];

    auto g_0_xxx_0_xzzz_0 = prim_buffer_0_sfsg[9];

    auto g_0_xxx_0_yyyy_0 = prim_buffer_0_sfsg[10];

    auto g_0_xxx_0_yyyz_0 = prim_buffer_0_sfsg[11];

    auto g_0_xxx_0_yyzz_0 = prim_buffer_0_sfsg[12];

    auto g_0_xxx_0_yzzz_0 = prim_buffer_0_sfsg[13];

    auto g_0_xxx_0_zzzz_0 = prim_buffer_0_sfsg[14];

    #pragma omp simd aligned(g_0_x_0_xxxx_0, g_0_x_0_xxxx_1, g_0_x_0_xxxy_0, g_0_x_0_xxxy_1, g_0_x_0_xxxz_0, g_0_x_0_xxxz_1, g_0_x_0_xxyy_0, g_0_x_0_xxyy_1, g_0_x_0_xxyz_0, g_0_x_0_xxyz_1, g_0_x_0_xxzz_0, g_0_x_0_xxzz_1, g_0_x_0_xyyy_0, g_0_x_0_xyyy_1, g_0_x_0_xyyz_0, g_0_x_0_xyyz_1, g_0_x_0_xyzz_0, g_0_x_0_xyzz_1, g_0_x_0_xzzz_0, g_0_x_0_xzzz_1, g_0_x_0_yyyy_0, g_0_x_0_yyyy_1, g_0_x_0_yyyz_0, g_0_x_0_yyyz_1, g_0_x_0_yyzz_0, g_0_x_0_yyzz_1, g_0_x_0_yzzz_0, g_0_x_0_yzzz_1, g_0_x_0_zzzz_0, g_0_x_0_zzzz_1, g_0_xx_0_xxx_1, g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxy_0, g_0_xx_0_xxxy_1, g_0_xx_0_xxxz_0, g_0_xx_0_xxxz_1, g_0_xx_0_xxy_1, g_0_xx_0_xxyy_0, g_0_xx_0_xxyy_1, g_0_xx_0_xxyz_0, g_0_xx_0_xxyz_1, g_0_xx_0_xxz_1, g_0_xx_0_xxzz_0, g_0_xx_0_xxzz_1, g_0_xx_0_xyy_1, g_0_xx_0_xyyy_0, g_0_xx_0_xyyy_1, g_0_xx_0_xyyz_0, g_0_xx_0_xyyz_1, g_0_xx_0_xyz_1, g_0_xx_0_xyzz_0, g_0_xx_0_xyzz_1, g_0_xx_0_xzz_1, g_0_xx_0_xzzz_0, g_0_xx_0_xzzz_1, g_0_xx_0_yyy_1, g_0_xx_0_yyyy_0, g_0_xx_0_yyyy_1, g_0_xx_0_yyyz_0, g_0_xx_0_yyyz_1, g_0_xx_0_yyz_1, g_0_xx_0_yyzz_0, g_0_xx_0_yyzz_1, g_0_xx_0_yzz_1, g_0_xx_0_yzzz_0, g_0_xx_0_yzzz_1, g_0_xx_0_zzz_1, g_0_xx_0_zzzz_0, g_0_xx_0_zzzz_1, g_0_xxx_0_xxxx_0, g_0_xxx_0_xxxy_0, g_0_xxx_0_xxxz_0, g_0_xxx_0_xxyy_0, g_0_xxx_0_xxyz_0, g_0_xxx_0_xxzz_0, g_0_xxx_0_xyyy_0, g_0_xxx_0_xyyz_0, g_0_xxx_0_xyzz_0, g_0_xxx_0_xzzz_0, g_0_xxx_0_yyyy_0, g_0_xxx_0_yyyz_0, g_0_xxx_0_yyzz_0, g_0_xxx_0_yzzz_0, g_0_xxx_0_zzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxxx_0[i] = 2.0 * g_0_x_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxx_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxx_0[i] * pb_x + g_0_xx_0_xxxx_1[i] * wp_x[i];

        g_0_xxx_0_xxxy_0[i] = 2.0 * g_0_x_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxy_0[i] * pb_x + g_0_xx_0_xxxy_1[i] * wp_x[i];

        g_0_xxx_0_xxxz_0[i] = 2.0 * g_0_x_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxz_0[i] * pb_x + g_0_xx_0_xxxz_1[i] * wp_x[i];

        g_0_xxx_0_xxyy_0[i] = 2.0 * g_0_x_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyy_0[i] * pb_x + g_0_xx_0_xxyy_1[i] * wp_x[i];

        g_0_xxx_0_xxyz_0[i] = 2.0 * g_0_x_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyz_0[i] * pb_x + g_0_xx_0_xxyz_1[i] * wp_x[i];

        g_0_xxx_0_xxzz_0[i] = 2.0 * g_0_x_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzz_0[i] * pb_x + g_0_xx_0_xxzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyy_0[i] = 2.0 * g_0_x_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyy_1[i] * fti_ab_0 + g_0_xx_0_yyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyy_0[i] * pb_x + g_0_xx_0_xyyy_1[i] * wp_x[i];

        g_0_xxx_0_xyyz_0[i] = 2.0 * g_0_x_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyz_1[i] * fti_ab_0 + g_0_xx_0_yyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyz_0[i] * pb_x + g_0_xx_0_xyyz_1[i] * wp_x[i];

        g_0_xxx_0_xyzz_0[i] = 2.0 * g_0_x_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyzz_1[i] * fti_ab_0 + g_0_xx_0_yzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzz_0[i] * pb_x + g_0_xx_0_xyzz_1[i] * wp_x[i];

        g_0_xxx_0_xzzz_0[i] = 2.0 * g_0_x_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzzz_1[i] * fti_ab_0 + g_0_xx_0_zzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzz_0[i] * pb_x + g_0_xx_0_xzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyy_0[i] = 2.0 * g_0_x_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyy_0[i] * pb_x + g_0_xx_0_yyyy_1[i] * wp_x[i];

        g_0_xxx_0_yyyz_0[i] = 2.0 * g_0_x_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyz_0[i] * pb_x + g_0_xx_0_yyyz_1[i] * wp_x[i];

        g_0_xxx_0_yyzz_0[i] = 2.0 * g_0_x_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyzz_1[i] * fti_ab_0 + g_0_xx_0_yyzz_0[i] * pb_x + g_0_xx_0_yyzz_1[i] * wp_x[i];

        g_0_xxx_0_yzzz_0[i] = 2.0 * g_0_x_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzz_0[i] * pb_x + g_0_xx_0_yzzz_1[i] * wp_x[i];

        g_0_xxx_0_zzzz_0[i] = 2.0 * g_0_x_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzz_0[i] * pb_x + g_0_xx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xxy_0_xxxx_0 = prim_buffer_0_sfsg[15];

    auto g_0_xxy_0_xxxy_0 = prim_buffer_0_sfsg[16];

    auto g_0_xxy_0_xxxz_0 = prim_buffer_0_sfsg[17];

    auto g_0_xxy_0_xxyy_0 = prim_buffer_0_sfsg[18];

    auto g_0_xxy_0_xxyz_0 = prim_buffer_0_sfsg[19];

    auto g_0_xxy_0_xxzz_0 = prim_buffer_0_sfsg[20];

    auto g_0_xxy_0_xyyy_0 = prim_buffer_0_sfsg[21];

    auto g_0_xxy_0_xyyz_0 = prim_buffer_0_sfsg[22];

    auto g_0_xxy_0_xyzz_0 = prim_buffer_0_sfsg[23];

    auto g_0_xxy_0_xzzz_0 = prim_buffer_0_sfsg[24];

    auto g_0_xxy_0_yyyy_0 = prim_buffer_0_sfsg[25];

    auto g_0_xxy_0_yyyz_0 = prim_buffer_0_sfsg[26];

    auto g_0_xxy_0_yyzz_0 = prim_buffer_0_sfsg[27];

    auto g_0_xxy_0_yzzz_0 = prim_buffer_0_sfsg[28];

    auto g_0_xxy_0_zzzz_0 = prim_buffer_0_sfsg[29];

    #pragma omp simd aligned(g_0_xx_0_xxx_1, g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxy_0, g_0_xx_0_xxxy_1, g_0_xx_0_xxxz_0, g_0_xx_0_xxxz_1, g_0_xx_0_xxy_1, g_0_xx_0_xxyy_0, g_0_xx_0_xxyy_1, g_0_xx_0_xxyz_0, g_0_xx_0_xxyz_1, g_0_xx_0_xxz_1, g_0_xx_0_xxzz_0, g_0_xx_0_xxzz_1, g_0_xx_0_xyy_1, g_0_xx_0_xyyy_0, g_0_xx_0_xyyy_1, g_0_xx_0_xyyz_0, g_0_xx_0_xyyz_1, g_0_xx_0_xyz_1, g_0_xx_0_xyzz_0, g_0_xx_0_xyzz_1, g_0_xx_0_xzz_1, g_0_xx_0_xzzz_0, g_0_xx_0_xzzz_1, g_0_xx_0_yyy_1, g_0_xx_0_yyyy_0, g_0_xx_0_yyyy_1, g_0_xx_0_yyyz_0, g_0_xx_0_yyyz_1, g_0_xx_0_yyz_1, g_0_xx_0_yyzz_0, g_0_xx_0_yyzz_1, g_0_xx_0_yzz_1, g_0_xx_0_yzzz_0, g_0_xx_0_yzzz_1, g_0_xx_0_zzz_1, g_0_xx_0_zzzz_0, g_0_xx_0_zzzz_1, g_0_xxy_0_xxxx_0, g_0_xxy_0_xxxy_0, g_0_xxy_0_xxxz_0, g_0_xxy_0_xxyy_0, g_0_xxy_0_xxyz_0, g_0_xxy_0_xxzz_0, g_0_xxy_0_xyyy_0, g_0_xxy_0_xyyz_0, g_0_xxy_0_xyzz_0, g_0_xxy_0_xzzz_0, g_0_xxy_0_yyyy_0, g_0_xxy_0_yyyz_0, g_0_xxy_0_yyzz_0, g_0_xxy_0_yzzz_0, g_0_xxy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxxx_0[i] = g_0_xx_0_xxxx_0[i] * pb_y + g_0_xx_0_xxxx_1[i] * wp_y[i];

        g_0_xxy_0_xxxy_0[i] = g_0_xx_0_xxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxy_0[i] * pb_y + g_0_xx_0_xxxy_1[i] * wp_y[i];

        g_0_xxy_0_xxxz_0[i] = g_0_xx_0_xxxz_0[i] * pb_y + g_0_xx_0_xxxz_1[i] * wp_y[i];

        g_0_xxy_0_xxyy_0[i] = 2.0 * g_0_xx_0_xxy_1[i] * fi_abcd_0 + g_0_xx_0_xxyy_0[i] * pb_y + g_0_xx_0_xxyy_1[i] * wp_y[i];

        g_0_xxy_0_xxyz_0[i] = g_0_xx_0_xxz_1[i] * fi_abcd_0 + g_0_xx_0_xxyz_0[i] * pb_y + g_0_xx_0_xxyz_1[i] * wp_y[i];

        g_0_xxy_0_xxzz_0[i] = g_0_xx_0_xxzz_0[i] * pb_y + g_0_xx_0_xxzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyy_0[i] = 3.0 * g_0_xx_0_xyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyy_0[i] * pb_y + g_0_xx_0_xyyy_1[i] * wp_y[i];

        g_0_xxy_0_xyyz_0[i] = 2.0 * g_0_xx_0_xyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyz_0[i] * pb_y + g_0_xx_0_xyyz_1[i] * wp_y[i];

        g_0_xxy_0_xyzz_0[i] = g_0_xx_0_xzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzz_0[i] * pb_y + g_0_xx_0_xyzz_1[i] * wp_y[i];

        g_0_xxy_0_xzzz_0[i] = g_0_xx_0_xzzz_0[i] * pb_y + g_0_xx_0_xzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyy_0[i] = 4.0 * g_0_xx_0_yyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyy_0[i] * pb_y + g_0_xx_0_yyyy_1[i] * wp_y[i];

        g_0_xxy_0_yyyz_0[i] = 3.0 * g_0_xx_0_yyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyz_0[i] * pb_y + g_0_xx_0_yyyz_1[i] * wp_y[i];

        g_0_xxy_0_yyzz_0[i] = 2.0 * g_0_xx_0_yzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzz_0[i] * pb_y + g_0_xx_0_yyzz_1[i] * wp_y[i];

        g_0_xxy_0_yzzz_0[i] = g_0_xx_0_zzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzz_0[i] * pb_y + g_0_xx_0_yzzz_1[i] * wp_y[i];

        g_0_xxy_0_zzzz_0[i] = g_0_xx_0_zzzz_0[i] * pb_y + g_0_xx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xxz_0_xxxx_0 = prim_buffer_0_sfsg[30];

    auto g_0_xxz_0_xxxy_0 = prim_buffer_0_sfsg[31];

    auto g_0_xxz_0_xxxz_0 = prim_buffer_0_sfsg[32];

    auto g_0_xxz_0_xxyy_0 = prim_buffer_0_sfsg[33];

    auto g_0_xxz_0_xxyz_0 = prim_buffer_0_sfsg[34];

    auto g_0_xxz_0_xxzz_0 = prim_buffer_0_sfsg[35];

    auto g_0_xxz_0_xyyy_0 = prim_buffer_0_sfsg[36];

    auto g_0_xxz_0_xyyz_0 = prim_buffer_0_sfsg[37];

    auto g_0_xxz_0_xyzz_0 = prim_buffer_0_sfsg[38];

    auto g_0_xxz_0_xzzz_0 = prim_buffer_0_sfsg[39];

    auto g_0_xxz_0_yyyy_0 = prim_buffer_0_sfsg[40];

    auto g_0_xxz_0_yyyz_0 = prim_buffer_0_sfsg[41];

    auto g_0_xxz_0_yyzz_0 = prim_buffer_0_sfsg[42];

    auto g_0_xxz_0_yzzz_0 = prim_buffer_0_sfsg[43];

    auto g_0_xxz_0_zzzz_0 = prim_buffer_0_sfsg[44];

    #pragma omp simd aligned(g_0_xx_0_xxx_1, g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxy_0, g_0_xx_0_xxxy_1, g_0_xx_0_xxxz_0, g_0_xx_0_xxxz_1, g_0_xx_0_xxy_1, g_0_xx_0_xxyy_0, g_0_xx_0_xxyy_1, g_0_xx_0_xxyz_0, g_0_xx_0_xxyz_1, g_0_xx_0_xxz_1, g_0_xx_0_xxzz_0, g_0_xx_0_xxzz_1, g_0_xx_0_xyy_1, g_0_xx_0_xyyy_0, g_0_xx_0_xyyy_1, g_0_xx_0_xyyz_0, g_0_xx_0_xyyz_1, g_0_xx_0_xyz_1, g_0_xx_0_xyzz_0, g_0_xx_0_xyzz_1, g_0_xx_0_xzz_1, g_0_xx_0_xzzz_0, g_0_xx_0_xzzz_1, g_0_xx_0_yyy_1, g_0_xx_0_yyyy_0, g_0_xx_0_yyyy_1, g_0_xx_0_yyyz_0, g_0_xx_0_yyyz_1, g_0_xx_0_yyz_1, g_0_xx_0_yyzz_0, g_0_xx_0_yyzz_1, g_0_xx_0_yzz_1, g_0_xx_0_yzzz_0, g_0_xx_0_yzzz_1, g_0_xx_0_zzz_1, g_0_xx_0_zzzz_0, g_0_xx_0_zzzz_1, g_0_xxz_0_xxxx_0, g_0_xxz_0_xxxy_0, g_0_xxz_0_xxxz_0, g_0_xxz_0_xxyy_0, g_0_xxz_0_xxyz_0, g_0_xxz_0_xxzz_0, g_0_xxz_0_xyyy_0, g_0_xxz_0_xyyz_0, g_0_xxz_0_xyzz_0, g_0_xxz_0_xzzz_0, g_0_xxz_0_yyyy_0, g_0_xxz_0_yyyz_0, g_0_xxz_0_yyzz_0, g_0_xxz_0_yzzz_0, g_0_xxz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxxx_0[i] = g_0_xx_0_xxxx_0[i] * pb_z + g_0_xx_0_xxxx_1[i] * wp_z[i];

        g_0_xxz_0_xxxy_0[i] = g_0_xx_0_xxxy_0[i] * pb_z + g_0_xx_0_xxxy_1[i] * wp_z[i];

        g_0_xxz_0_xxxz_0[i] = g_0_xx_0_xxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxz_0[i] * pb_z + g_0_xx_0_xxxz_1[i] * wp_z[i];

        g_0_xxz_0_xxyy_0[i] = g_0_xx_0_xxyy_0[i] * pb_z + g_0_xx_0_xxyy_1[i] * wp_z[i];

        g_0_xxz_0_xxyz_0[i] = g_0_xx_0_xxy_1[i] * fi_abcd_0 + g_0_xx_0_xxyz_0[i] * pb_z + g_0_xx_0_xxyz_1[i] * wp_z[i];

        g_0_xxz_0_xxzz_0[i] = 2.0 * g_0_xx_0_xxz_1[i] * fi_abcd_0 + g_0_xx_0_xxzz_0[i] * pb_z + g_0_xx_0_xxzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyy_0[i] = g_0_xx_0_xyyy_0[i] * pb_z + g_0_xx_0_xyyy_1[i] * wp_z[i];

        g_0_xxz_0_xyyz_0[i] = g_0_xx_0_xyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyz_0[i] * pb_z + g_0_xx_0_xyyz_1[i] * wp_z[i];

        g_0_xxz_0_xyzz_0[i] = 2.0 * g_0_xx_0_xyz_1[i] * fi_abcd_0 + g_0_xx_0_xyzz_0[i] * pb_z + g_0_xx_0_xyzz_1[i] * wp_z[i];

        g_0_xxz_0_xzzz_0[i] = 3.0 * g_0_xx_0_xzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzz_0[i] * pb_z + g_0_xx_0_xzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyy_0[i] = g_0_xx_0_yyyy_0[i] * pb_z + g_0_xx_0_yyyy_1[i] * wp_z[i];

        g_0_xxz_0_yyyz_0[i] = g_0_xx_0_yyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyz_0[i] * pb_z + g_0_xx_0_yyyz_1[i] * wp_z[i];

        g_0_xxz_0_yyzz_0[i] = 2.0 * g_0_xx_0_yyz_1[i] * fi_abcd_0 + g_0_xx_0_yyzz_0[i] * pb_z + g_0_xx_0_yyzz_1[i] * wp_z[i];

        g_0_xxz_0_yzzz_0[i] = 3.0 * g_0_xx_0_yzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzz_0[i] * pb_z + g_0_xx_0_yzzz_1[i] * wp_z[i];

        g_0_xxz_0_zzzz_0[i] = 4.0 * g_0_xx_0_zzz_1[i] * fi_abcd_0 + g_0_xx_0_zzzz_0[i] * pb_z + g_0_xx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xyy_0_xxxx_0 = prim_buffer_0_sfsg[45];

    auto g_0_xyy_0_xxxy_0 = prim_buffer_0_sfsg[46];

    auto g_0_xyy_0_xxxz_0 = prim_buffer_0_sfsg[47];

    auto g_0_xyy_0_xxyy_0 = prim_buffer_0_sfsg[48];

    auto g_0_xyy_0_xxyz_0 = prim_buffer_0_sfsg[49];

    auto g_0_xyy_0_xxzz_0 = prim_buffer_0_sfsg[50];

    auto g_0_xyy_0_xyyy_0 = prim_buffer_0_sfsg[51];

    auto g_0_xyy_0_xyyz_0 = prim_buffer_0_sfsg[52];

    auto g_0_xyy_0_xyzz_0 = prim_buffer_0_sfsg[53];

    auto g_0_xyy_0_xzzz_0 = prim_buffer_0_sfsg[54];

    auto g_0_xyy_0_yyyy_0 = prim_buffer_0_sfsg[55];

    auto g_0_xyy_0_yyyz_0 = prim_buffer_0_sfsg[56];

    auto g_0_xyy_0_yyzz_0 = prim_buffer_0_sfsg[57];

    auto g_0_xyy_0_yzzz_0 = prim_buffer_0_sfsg[58];

    auto g_0_xyy_0_zzzz_0 = prim_buffer_0_sfsg[59];

    #pragma omp simd aligned(g_0_xyy_0_xxxx_0, g_0_xyy_0_xxxy_0, g_0_xyy_0_xxxz_0, g_0_xyy_0_xxyy_0, g_0_xyy_0_xxyz_0, g_0_xyy_0_xxzz_0, g_0_xyy_0_xyyy_0, g_0_xyy_0_xyyz_0, g_0_xyy_0_xyzz_0, g_0_xyy_0_xzzz_0, g_0_xyy_0_yyyy_0, g_0_xyy_0_yyyz_0, g_0_xyy_0_yyzz_0, g_0_xyy_0_yzzz_0, g_0_xyy_0_zzzz_0, g_0_yy_0_xxx_1, g_0_yy_0_xxxx_0, g_0_yy_0_xxxx_1, g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxxz_0, g_0_yy_0_xxxz_1, g_0_yy_0_xxy_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xxyz_0, g_0_yy_0_xxyz_1, g_0_yy_0_xxz_1, g_0_yy_0_xxzz_0, g_0_yy_0_xxzz_1, g_0_yy_0_xyy_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_xyyz_0, g_0_yy_0_xyyz_1, g_0_yy_0_xyz_1, g_0_yy_0_xyzz_0, g_0_yy_0_xyzz_1, g_0_yy_0_xzz_1, g_0_yy_0_xzzz_0, g_0_yy_0_xzzz_1, g_0_yy_0_yyy_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yy_0_yyyz_0, g_0_yy_0_yyyz_1, g_0_yy_0_yyz_1, g_0_yy_0_yyzz_0, g_0_yy_0_yyzz_1, g_0_yy_0_yzz_1, g_0_yy_0_yzzz_0, g_0_yy_0_yzzz_1, g_0_yy_0_zzz_1, g_0_yy_0_zzzz_0, g_0_yy_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxxx_0[i] = 4.0 * g_0_yy_0_xxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxx_0[i] * pb_x + g_0_yy_0_xxxx_1[i] * wp_x[i];

        g_0_xyy_0_xxxy_0[i] = 3.0 * g_0_yy_0_xxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxy_0[i] * pb_x + g_0_yy_0_xxxy_1[i] * wp_x[i];

        g_0_xyy_0_xxxz_0[i] = 3.0 * g_0_yy_0_xxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxz_0[i] * pb_x + g_0_yy_0_xxxz_1[i] * wp_x[i];

        g_0_xyy_0_xxyy_0[i] = 2.0 * g_0_yy_0_xyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyy_0[i] * pb_x + g_0_yy_0_xxyy_1[i] * wp_x[i];

        g_0_xyy_0_xxyz_0[i] = 2.0 * g_0_yy_0_xyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyz_0[i] * pb_x + g_0_yy_0_xxyz_1[i] * wp_x[i];

        g_0_xyy_0_xxzz_0[i] = 2.0 * g_0_yy_0_xzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzz_0[i] * pb_x + g_0_yy_0_xxzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyy_0[i] = g_0_yy_0_yyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyy_0[i] * pb_x + g_0_yy_0_xyyy_1[i] * wp_x[i];

        g_0_xyy_0_xyyz_0[i] = g_0_yy_0_yyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyz_0[i] * pb_x + g_0_yy_0_xyyz_1[i] * wp_x[i];

        g_0_xyy_0_xyzz_0[i] = g_0_yy_0_yzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzz_0[i] * pb_x + g_0_yy_0_xyzz_1[i] * wp_x[i];

        g_0_xyy_0_xzzz_0[i] = g_0_yy_0_zzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzz_0[i] * pb_x + g_0_yy_0_xzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyy_0[i] = g_0_yy_0_yyyy_0[i] * pb_x + g_0_yy_0_yyyy_1[i] * wp_x[i];

        g_0_xyy_0_yyyz_0[i] = g_0_yy_0_yyyz_0[i] * pb_x + g_0_yy_0_yyyz_1[i] * wp_x[i];

        g_0_xyy_0_yyzz_0[i] = g_0_yy_0_yyzz_0[i] * pb_x + g_0_yy_0_yyzz_1[i] * wp_x[i];

        g_0_xyy_0_yzzz_0[i] = g_0_yy_0_yzzz_0[i] * pb_x + g_0_yy_0_yzzz_1[i] * wp_x[i];

        g_0_xyy_0_zzzz_0[i] = g_0_yy_0_zzzz_0[i] * pb_x + g_0_yy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xyz_0_xxxx_0 = prim_buffer_0_sfsg[60];

    auto g_0_xyz_0_xxxy_0 = prim_buffer_0_sfsg[61];

    auto g_0_xyz_0_xxxz_0 = prim_buffer_0_sfsg[62];

    auto g_0_xyz_0_xxyy_0 = prim_buffer_0_sfsg[63];

    auto g_0_xyz_0_xxyz_0 = prim_buffer_0_sfsg[64];

    auto g_0_xyz_0_xxzz_0 = prim_buffer_0_sfsg[65];

    auto g_0_xyz_0_xyyy_0 = prim_buffer_0_sfsg[66];

    auto g_0_xyz_0_xyyz_0 = prim_buffer_0_sfsg[67];

    auto g_0_xyz_0_xyzz_0 = prim_buffer_0_sfsg[68];

    auto g_0_xyz_0_xzzz_0 = prim_buffer_0_sfsg[69];

    auto g_0_xyz_0_yyyy_0 = prim_buffer_0_sfsg[70];

    auto g_0_xyz_0_yyyz_0 = prim_buffer_0_sfsg[71];

    auto g_0_xyz_0_yyzz_0 = prim_buffer_0_sfsg[72];

    auto g_0_xyz_0_yzzz_0 = prim_buffer_0_sfsg[73];

    auto g_0_xyz_0_zzzz_0 = prim_buffer_0_sfsg[74];

    #pragma omp simd aligned(g_0_xy_0_xxxy_0, g_0_xy_0_xxxy_1, g_0_xy_0_xxyy_0, g_0_xy_0_xxyy_1, g_0_xy_0_xyyy_0, g_0_xy_0_xyyy_1, g_0_xyz_0_xxxx_0, g_0_xyz_0_xxxy_0, g_0_xyz_0_xxxz_0, g_0_xyz_0_xxyy_0, g_0_xyz_0_xxyz_0, g_0_xyz_0_xxzz_0, g_0_xyz_0_xyyy_0, g_0_xyz_0_xyyz_0, g_0_xyz_0_xyzz_0, g_0_xyz_0_xzzz_0, g_0_xyz_0_yyyy_0, g_0_xyz_0_yyyz_0, g_0_xyz_0_yyzz_0, g_0_xyz_0_yzzz_0, g_0_xyz_0_zzzz_0, g_0_xz_0_xxxx_0, g_0_xz_0_xxxx_1, g_0_xz_0_xxxz_0, g_0_xz_0_xxxz_1, g_0_xz_0_xxzz_0, g_0_xz_0_xxzz_1, g_0_xz_0_xzzz_0, g_0_xz_0_xzzz_1, g_0_yz_0_xxyz_0, g_0_yz_0_xxyz_1, g_0_yz_0_xyyz_0, g_0_yz_0_xyyz_1, g_0_yz_0_xyz_1, g_0_yz_0_xyzz_0, g_0_yz_0_xyzz_1, g_0_yz_0_yyyy_0, g_0_yz_0_yyyy_1, g_0_yz_0_yyyz_0, g_0_yz_0_yyyz_1, g_0_yz_0_yyz_1, g_0_yz_0_yyzz_0, g_0_yz_0_yyzz_1, g_0_yz_0_yzz_1, g_0_yz_0_yzzz_0, g_0_yz_0_yzzz_1, g_0_yz_0_zzzz_0, g_0_yz_0_zzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxxx_0[i] = g_0_xz_0_xxxx_0[i] * pb_y + g_0_xz_0_xxxx_1[i] * wp_y[i];

        g_0_xyz_0_xxxy_0[i] = g_0_xy_0_xxxy_0[i] * pb_z + g_0_xy_0_xxxy_1[i] * wp_z[i];

        g_0_xyz_0_xxxz_0[i] = g_0_xz_0_xxxz_0[i] * pb_y + g_0_xz_0_xxxz_1[i] * wp_y[i];

        g_0_xyz_0_xxyy_0[i] = g_0_xy_0_xxyy_0[i] * pb_z + g_0_xy_0_xxyy_1[i] * wp_z[i];

        g_0_xyz_0_xxyz_0[i] = 2.0 * g_0_yz_0_xyz_1[i] * fi_abcd_0 + g_0_yz_0_xxyz_0[i] * pb_x + g_0_yz_0_xxyz_1[i] * wp_x[i];

        g_0_xyz_0_xxzz_0[i] = g_0_xz_0_xxzz_0[i] * pb_y + g_0_xz_0_xxzz_1[i] * wp_y[i];

        g_0_xyz_0_xyyy_0[i] = g_0_xy_0_xyyy_0[i] * pb_z + g_0_xy_0_xyyy_1[i] * wp_z[i];

        g_0_xyz_0_xyyz_0[i] = g_0_yz_0_yyz_1[i] * fi_abcd_0 + g_0_yz_0_xyyz_0[i] * pb_x + g_0_yz_0_xyyz_1[i] * wp_x[i];

        g_0_xyz_0_xyzz_0[i] = g_0_yz_0_yzz_1[i] * fi_abcd_0 + g_0_yz_0_xyzz_0[i] * pb_x + g_0_yz_0_xyzz_1[i] * wp_x[i];

        g_0_xyz_0_xzzz_0[i] = g_0_xz_0_xzzz_0[i] * pb_y + g_0_xz_0_xzzz_1[i] * wp_y[i];

        g_0_xyz_0_yyyy_0[i] = g_0_yz_0_yyyy_0[i] * pb_x + g_0_yz_0_yyyy_1[i] * wp_x[i];

        g_0_xyz_0_yyyz_0[i] = g_0_yz_0_yyyz_0[i] * pb_x + g_0_yz_0_yyyz_1[i] * wp_x[i];

        g_0_xyz_0_yyzz_0[i] = g_0_yz_0_yyzz_0[i] * pb_x + g_0_yz_0_yyzz_1[i] * wp_x[i];

        g_0_xyz_0_yzzz_0[i] = g_0_yz_0_yzzz_0[i] * pb_x + g_0_yz_0_yzzz_1[i] * wp_x[i];

        g_0_xyz_0_zzzz_0[i] = g_0_yz_0_zzzz_0[i] * pb_x + g_0_yz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 75-90 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_xzz_0_xxxx_0 = prim_buffer_0_sfsg[75];

    auto g_0_xzz_0_xxxy_0 = prim_buffer_0_sfsg[76];

    auto g_0_xzz_0_xxxz_0 = prim_buffer_0_sfsg[77];

    auto g_0_xzz_0_xxyy_0 = prim_buffer_0_sfsg[78];

    auto g_0_xzz_0_xxyz_0 = prim_buffer_0_sfsg[79];

    auto g_0_xzz_0_xxzz_0 = prim_buffer_0_sfsg[80];

    auto g_0_xzz_0_xyyy_0 = prim_buffer_0_sfsg[81];

    auto g_0_xzz_0_xyyz_0 = prim_buffer_0_sfsg[82];

    auto g_0_xzz_0_xyzz_0 = prim_buffer_0_sfsg[83];

    auto g_0_xzz_0_xzzz_0 = prim_buffer_0_sfsg[84];

    auto g_0_xzz_0_yyyy_0 = prim_buffer_0_sfsg[85];

    auto g_0_xzz_0_yyyz_0 = prim_buffer_0_sfsg[86];

    auto g_0_xzz_0_yyzz_0 = prim_buffer_0_sfsg[87];

    auto g_0_xzz_0_yzzz_0 = prim_buffer_0_sfsg[88];

    auto g_0_xzz_0_zzzz_0 = prim_buffer_0_sfsg[89];

    #pragma omp simd aligned(g_0_xzz_0_xxxx_0, g_0_xzz_0_xxxy_0, g_0_xzz_0_xxxz_0, g_0_xzz_0_xxyy_0, g_0_xzz_0_xxyz_0, g_0_xzz_0_xxzz_0, g_0_xzz_0_xyyy_0, g_0_xzz_0_xyyz_0, g_0_xzz_0_xyzz_0, g_0_xzz_0_xzzz_0, g_0_xzz_0_yyyy_0, g_0_xzz_0_yyyz_0, g_0_xzz_0_yyzz_0, g_0_xzz_0_yzzz_0, g_0_xzz_0_zzzz_0, g_0_zz_0_xxx_1, g_0_zz_0_xxxx_0, g_0_zz_0_xxxx_1, g_0_zz_0_xxxy_0, g_0_zz_0_xxxy_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxy_1, g_0_zz_0_xxyy_0, g_0_zz_0_xxyy_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyy_1, g_0_zz_0_xyyy_0, g_0_zz_0_xyyy_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyy_1, g_0_zz_0_yyyy_0, g_0_zz_0_yyyy_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxxx_0[i] = 4.0 * g_0_zz_0_xxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxx_0[i] * pb_x + g_0_zz_0_xxxx_1[i] * wp_x[i];

        g_0_xzz_0_xxxy_0[i] = 3.0 * g_0_zz_0_xxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxy_0[i] * pb_x + g_0_zz_0_xxxy_1[i] * wp_x[i];

        g_0_xzz_0_xxxz_0[i] = 3.0 * g_0_zz_0_xxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxz_0[i] * pb_x + g_0_zz_0_xxxz_1[i] * wp_x[i];

        g_0_xzz_0_xxyy_0[i] = 2.0 * g_0_zz_0_xyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyy_0[i] * pb_x + g_0_zz_0_xxyy_1[i] * wp_x[i];

        g_0_xzz_0_xxyz_0[i] = 2.0 * g_0_zz_0_xyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyz_0[i] * pb_x + g_0_zz_0_xxyz_1[i] * wp_x[i];

        g_0_xzz_0_xxzz_0[i] = 2.0 * g_0_zz_0_xzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzz_0[i] * pb_x + g_0_zz_0_xxzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyy_0[i] = g_0_zz_0_yyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyy_0[i] * pb_x + g_0_zz_0_xyyy_1[i] * wp_x[i];

        g_0_xzz_0_xyyz_0[i] = g_0_zz_0_yyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyz_0[i] * pb_x + g_0_zz_0_xyyz_1[i] * wp_x[i];

        g_0_xzz_0_xyzz_0[i] = g_0_zz_0_yzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzz_0[i] * pb_x + g_0_zz_0_xyzz_1[i] * wp_x[i];

        g_0_xzz_0_xzzz_0[i] = g_0_zz_0_zzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzz_0[i] * pb_x + g_0_zz_0_xzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyy_0[i] = g_0_zz_0_yyyy_0[i] * pb_x + g_0_zz_0_yyyy_1[i] * wp_x[i];

        g_0_xzz_0_yyyz_0[i] = g_0_zz_0_yyyz_0[i] * pb_x + g_0_zz_0_yyyz_1[i] * wp_x[i];

        g_0_xzz_0_yyzz_0[i] = g_0_zz_0_yyzz_0[i] * pb_x + g_0_zz_0_yyzz_1[i] * wp_x[i];

        g_0_xzz_0_yzzz_0[i] = g_0_zz_0_yzzz_0[i] * pb_x + g_0_zz_0_yzzz_1[i] * wp_x[i];

        g_0_xzz_0_zzzz_0[i] = g_0_zz_0_zzzz_0[i] * pb_x + g_0_zz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_yyy_0_xxxx_0 = prim_buffer_0_sfsg[90];

    auto g_0_yyy_0_xxxy_0 = prim_buffer_0_sfsg[91];

    auto g_0_yyy_0_xxxz_0 = prim_buffer_0_sfsg[92];

    auto g_0_yyy_0_xxyy_0 = prim_buffer_0_sfsg[93];

    auto g_0_yyy_0_xxyz_0 = prim_buffer_0_sfsg[94];

    auto g_0_yyy_0_xxzz_0 = prim_buffer_0_sfsg[95];

    auto g_0_yyy_0_xyyy_0 = prim_buffer_0_sfsg[96];

    auto g_0_yyy_0_xyyz_0 = prim_buffer_0_sfsg[97];

    auto g_0_yyy_0_xyzz_0 = prim_buffer_0_sfsg[98];

    auto g_0_yyy_0_xzzz_0 = prim_buffer_0_sfsg[99];

    auto g_0_yyy_0_yyyy_0 = prim_buffer_0_sfsg[100];

    auto g_0_yyy_0_yyyz_0 = prim_buffer_0_sfsg[101];

    auto g_0_yyy_0_yyzz_0 = prim_buffer_0_sfsg[102];

    auto g_0_yyy_0_yzzz_0 = prim_buffer_0_sfsg[103];

    auto g_0_yyy_0_zzzz_0 = prim_buffer_0_sfsg[104];

    #pragma omp simd aligned(g_0_y_0_xxxx_0, g_0_y_0_xxxx_1, g_0_y_0_xxxy_0, g_0_y_0_xxxy_1, g_0_y_0_xxxz_0, g_0_y_0_xxxz_1, g_0_y_0_xxyy_0, g_0_y_0_xxyy_1, g_0_y_0_xxyz_0, g_0_y_0_xxyz_1, g_0_y_0_xxzz_0, g_0_y_0_xxzz_1, g_0_y_0_xyyy_0, g_0_y_0_xyyy_1, g_0_y_0_xyyz_0, g_0_y_0_xyyz_1, g_0_y_0_xyzz_0, g_0_y_0_xyzz_1, g_0_y_0_xzzz_0, g_0_y_0_xzzz_1, g_0_y_0_yyyy_0, g_0_y_0_yyyy_1, g_0_y_0_yyyz_0, g_0_y_0_yyyz_1, g_0_y_0_yyzz_0, g_0_y_0_yyzz_1, g_0_y_0_yzzz_0, g_0_y_0_yzzz_1, g_0_y_0_zzzz_0, g_0_y_0_zzzz_1, g_0_yy_0_xxx_1, g_0_yy_0_xxxx_0, g_0_yy_0_xxxx_1, g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxxz_0, g_0_yy_0_xxxz_1, g_0_yy_0_xxy_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xxyz_0, g_0_yy_0_xxyz_1, g_0_yy_0_xxz_1, g_0_yy_0_xxzz_0, g_0_yy_0_xxzz_1, g_0_yy_0_xyy_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_xyyz_0, g_0_yy_0_xyyz_1, g_0_yy_0_xyz_1, g_0_yy_0_xyzz_0, g_0_yy_0_xyzz_1, g_0_yy_0_xzz_1, g_0_yy_0_xzzz_0, g_0_yy_0_xzzz_1, g_0_yy_0_yyy_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yy_0_yyyz_0, g_0_yy_0_yyyz_1, g_0_yy_0_yyz_1, g_0_yy_0_yyzz_0, g_0_yy_0_yyzz_1, g_0_yy_0_yzz_1, g_0_yy_0_yzzz_0, g_0_yy_0_yzzz_1, g_0_yy_0_zzz_1, g_0_yy_0_zzzz_0, g_0_yy_0_zzzz_1, g_0_yyy_0_xxxx_0, g_0_yyy_0_xxxy_0, g_0_yyy_0_xxxz_0, g_0_yyy_0_xxyy_0, g_0_yyy_0_xxyz_0, g_0_yyy_0_xxzz_0, g_0_yyy_0_xyyy_0, g_0_yyy_0_xyyz_0, g_0_yyy_0_xyzz_0, g_0_yyy_0_xzzz_0, g_0_yyy_0_yyyy_0, g_0_yyy_0_yyyz_0, g_0_yyy_0_yyzz_0, g_0_yyy_0_yzzz_0, g_0_yyy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxxx_0[i] = 2.0 * g_0_y_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxx_1[i] * fti_ab_0 + g_0_yy_0_xxxx_0[i] * pb_y + g_0_yy_0_xxxx_1[i] * wp_y[i];

        g_0_yyy_0_xxxy_0[i] = 2.0 * g_0_y_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxy_1[i] * fti_ab_0 + g_0_yy_0_xxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxy_0[i] * pb_y + g_0_yy_0_xxxy_1[i] * wp_y[i];

        g_0_yyy_0_xxxz_0[i] = 2.0 * g_0_y_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxz_1[i] * fti_ab_0 + g_0_yy_0_xxxz_0[i] * pb_y + g_0_yy_0_xxxz_1[i] * wp_y[i];

        g_0_yyy_0_xxyy_0[i] = 2.0 * g_0_y_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxy_1[i] * fi_abcd_0 + g_0_yy_0_xxyy_0[i] * pb_y + g_0_yy_0_xxyy_1[i] * wp_y[i];

        g_0_yyy_0_xxyz_0[i] = 2.0 * g_0_y_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyz_1[i] * fti_ab_0 + g_0_yy_0_xxz_1[i] * fi_abcd_0 + g_0_yy_0_xxyz_0[i] * pb_y + g_0_yy_0_xxyz_1[i] * wp_y[i];

        g_0_yyy_0_xxzz_0[i] = 2.0 * g_0_y_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxzz_1[i] * fti_ab_0 + g_0_yy_0_xxzz_0[i] * pb_y + g_0_yy_0_xxzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyy_0[i] = 2.0 * g_0_y_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyy_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyy_0[i] * pb_y + g_0_yy_0_xyyy_1[i] * wp_y[i];

        g_0_yyy_0_xyyz_0[i] = 2.0 * g_0_y_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyz_0[i] * pb_y + g_0_yy_0_xyyz_1[i] * wp_y[i];

        g_0_yyy_0_xyzz_0[i] = 2.0 * g_0_y_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyzz_1[i] * fti_ab_0 + g_0_yy_0_xzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzz_0[i] * pb_y + g_0_yy_0_xyzz_1[i] * wp_y[i];

        g_0_yyy_0_xzzz_0[i] = 2.0 * g_0_y_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzz_0[i] * pb_y + g_0_yy_0_xzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyy_0[i] = 2.0 * g_0_y_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyy_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_yyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyy_0[i] * pb_y + g_0_yy_0_yyyy_1[i] * wp_y[i];

        g_0_yyy_0_yyyz_0[i] = 2.0 * g_0_y_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_yyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyz_0[i] * pb_y + g_0_yy_0_yyyz_1[i] * wp_y[i];

        g_0_yyy_0_yyzz_0[i] = 2.0 * g_0_y_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_yzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzz_0[i] * pb_y + g_0_yy_0_yyzz_1[i] * wp_y[i];

        g_0_yyy_0_yzzz_0[i] = 2.0 * g_0_y_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzzz_1[i] * fti_ab_0 + g_0_yy_0_zzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzz_0[i] * pb_y + g_0_yy_0_yzzz_1[i] * wp_y[i];

        g_0_yyy_0_zzzz_0[i] = 2.0 * g_0_y_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzz_0[i] * pb_y + g_0_yy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 105-120 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_yyz_0_xxxx_0 = prim_buffer_0_sfsg[105];

    auto g_0_yyz_0_xxxy_0 = prim_buffer_0_sfsg[106];

    auto g_0_yyz_0_xxxz_0 = prim_buffer_0_sfsg[107];

    auto g_0_yyz_0_xxyy_0 = prim_buffer_0_sfsg[108];

    auto g_0_yyz_0_xxyz_0 = prim_buffer_0_sfsg[109];

    auto g_0_yyz_0_xxzz_0 = prim_buffer_0_sfsg[110];

    auto g_0_yyz_0_xyyy_0 = prim_buffer_0_sfsg[111];

    auto g_0_yyz_0_xyyz_0 = prim_buffer_0_sfsg[112];

    auto g_0_yyz_0_xyzz_0 = prim_buffer_0_sfsg[113];

    auto g_0_yyz_0_xzzz_0 = prim_buffer_0_sfsg[114];

    auto g_0_yyz_0_yyyy_0 = prim_buffer_0_sfsg[115];

    auto g_0_yyz_0_yyyz_0 = prim_buffer_0_sfsg[116];

    auto g_0_yyz_0_yyzz_0 = prim_buffer_0_sfsg[117];

    auto g_0_yyz_0_yzzz_0 = prim_buffer_0_sfsg[118];

    auto g_0_yyz_0_zzzz_0 = prim_buffer_0_sfsg[119];

    #pragma omp simd aligned(g_0_yy_0_xxx_1, g_0_yy_0_xxxx_0, g_0_yy_0_xxxx_1, g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxxz_0, g_0_yy_0_xxxz_1, g_0_yy_0_xxy_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xxyz_0, g_0_yy_0_xxyz_1, g_0_yy_0_xxz_1, g_0_yy_0_xxzz_0, g_0_yy_0_xxzz_1, g_0_yy_0_xyy_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_xyyz_0, g_0_yy_0_xyyz_1, g_0_yy_0_xyz_1, g_0_yy_0_xyzz_0, g_0_yy_0_xyzz_1, g_0_yy_0_xzz_1, g_0_yy_0_xzzz_0, g_0_yy_0_xzzz_1, g_0_yy_0_yyy_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yy_0_yyyz_0, g_0_yy_0_yyyz_1, g_0_yy_0_yyz_1, g_0_yy_0_yyzz_0, g_0_yy_0_yyzz_1, g_0_yy_0_yzz_1, g_0_yy_0_yzzz_0, g_0_yy_0_yzzz_1, g_0_yy_0_zzz_1, g_0_yy_0_zzzz_0, g_0_yy_0_zzzz_1, g_0_yyz_0_xxxx_0, g_0_yyz_0_xxxy_0, g_0_yyz_0_xxxz_0, g_0_yyz_0_xxyy_0, g_0_yyz_0_xxyz_0, g_0_yyz_0_xxzz_0, g_0_yyz_0_xyyy_0, g_0_yyz_0_xyyz_0, g_0_yyz_0_xyzz_0, g_0_yyz_0_xzzz_0, g_0_yyz_0_yyyy_0, g_0_yyz_0_yyyz_0, g_0_yyz_0_yyzz_0, g_0_yyz_0_yzzz_0, g_0_yyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxxx_0[i] = g_0_yy_0_xxxx_0[i] * pb_z + g_0_yy_0_xxxx_1[i] * wp_z[i];

        g_0_yyz_0_xxxy_0[i] = g_0_yy_0_xxxy_0[i] * pb_z + g_0_yy_0_xxxy_1[i] * wp_z[i];

        g_0_yyz_0_xxxz_0[i] = g_0_yy_0_xxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxz_0[i] * pb_z + g_0_yy_0_xxxz_1[i] * wp_z[i];

        g_0_yyz_0_xxyy_0[i] = g_0_yy_0_xxyy_0[i] * pb_z + g_0_yy_0_xxyy_1[i] * wp_z[i];

        g_0_yyz_0_xxyz_0[i] = g_0_yy_0_xxy_1[i] * fi_abcd_0 + g_0_yy_0_xxyz_0[i] * pb_z + g_0_yy_0_xxyz_1[i] * wp_z[i];

        g_0_yyz_0_xxzz_0[i] = 2.0 * g_0_yy_0_xxz_1[i] * fi_abcd_0 + g_0_yy_0_xxzz_0[i] * pb_z + g_0_yy_0_xxzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyy_0[i] = g_0_yy_0_xyyy_0[i] * pb_z + g_0_yy_0_xyyy_1[i] * wp_z[i];

        g_0_yyz_0_xyyz_0[i] = g_0_yy_0_xyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyz_0[i] * pb_z + g_0_yy_0_xyyz_1[i] * wp_z[i];

        g_0_yyz_0_xyzz_0[i] = 2.0 * g_0_yy_0_xyz_1[i] * fi_abcd_0 + g_0_yy_0_xyzz_0[i] * pb_z + g_0_yy_0_xyzz_1[i] * wp_z[i];

        g_0_yyz_0_xzzz_0[i] = 3.0 * g_0_yy_0_xzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzz_0[i] * pb_z + g_0_yy_0_xzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyy_0[i] = g_0_yy_0_yyyy_0[i] * pb_z + g_0_yy_0_yyyy_1[i] * wp_z[i];

        g_0_yyz_0_yyyz_0[i] = g_0_yy_0_yyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyz_0[i] * pb_z + g_0_yy_0_yyyz_1[i] * wp_z[i];

        g_0_yyz_0_yyzz_0[i] = 2.0 * g_0_yy_0_yyz_1[i] * fi_abcd_0 + g_0_yy_0_yyzz_0[i] * pb_z + g_0_yy_0_yyzz_1[i] * wp_z[i];

        g_0_yyz_0_yzzz_0[i] = 3.0 * g_0_yy_0_yzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzz_0[i] * pb_z + g_0_yy_0_yzzz_1[i] * wp_z[i];

        g_0_yyz_0_zzzz_0[i] = 4.0 * g_0_yy_0_zzz_1[i] * fi_abcd_0 + g_0_yy_0_zzzz_0[i] * pb_z + g_0_yy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 120-135 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_yzz_0_xxxx_0 = prim_buffer_0_sfsg[120];

    auto g_0_yzz_0_xxxy_0 = prim_buffer_0_sfsg[121];

    auto g_0_yzz_0_xxxz_0 = prim_buffer_0_sfsg[122];

    auto g_0_yzz_0_xxyy_0 = prim_buffer_0_sfsg[123];

    auto g_0_yzz_0_xxyz_0 = prim_buffer_0_sfsg[124];

    auto g_0_yzz_0_xxzz_0 = prim_buffer_0_sfsg[125];

    auto g_0_yzz_0_xyyy_0 = prim_buffer_0_sfsg[126];

    auto g_0_yzz_0_xyyz_0 = prim_buffer_0_sfsg[127];

    auto g_0_yzz_0_xyzz_0 = prim_buffer_0_sfsg[128];

    auto g_0_yzz_0_xzzz_0 = prim_buffer_0_sfsg[129];

    auto g_0_yzz_0_yyyy_0 = prim_buffer_0_sfsg[130];

    auto g_0_yzz_0_yyyz_0 = prim_buffer_0_sfsg[131];

    auto g_0_yzz_0_yyzz_0 = prim_buffer_0_sfsg[132];

    auto g_0_yzz_0_yzzz_0 = prim_buffer_0_sfsg[133];

    auto g_0_yzz_0_zzzz_0 = prim_buffer_0_sfsg[134];

    #pragma omp simd aligned(g_0_yzz_0_xxxx_0, g_0_yzz_0_xxxy_0, g_0_yzz_0_xxxz_0, g_0_yzz_0_xxyy_0, g_0_yzz_0_xxyz_0, g_0_yzz_0_xxzz_0, g_0_yzz_0_xyyy_0, g_0_yzz_0_xyyz_0, g_0_yzz_0_xyzz_0, g_0_yzz_0_xzzz_0, g_0_yzz_0_yyyy_0, g_0_yzz_0_yyyz_0, g_0_yzz_0_yyzz_0, g_0_yzz_0_yzzz_0, g_0_yzz_0_zzzz_0, g_0_zz_0_xxx_1, g_0_zz_0_xxxx_0, g_0_zz_0_xxxx_1, g_0_zz_0_xxxy_0, g_0_zz_0_xxxy_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxy_1, g_0_zz_0_xxyy_0, g_0_zz_0_xxyy_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyy_1, g_0_zz_0_xyyy_0, g_0_zz_0_xyyy_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyy_1, g_0_zz_0_yyyy_0, g_0_zz_0_yyyy_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxxx_0[i] = g_0_zz_0_xxxx_0[i] * pb_y + g_0_zz_0_xxxx_1[i] * wp_y[i];

        g_0_yzz_0_xxxy_0[i] = g_0_zz_0_xxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxy_0[i] * pb_y + g_0_zz_0_xxxy_1[i] * wp_y[i];

        g_0_yzz_0_xxxz_0[i] = g_0_zz_0_xxxz_0[i] * pb_y + g_0_zz_0_xxxz_1[i] * wp_y[i];

        g_0_yzz_0_xxyy_0[i] = 2.0 * g_0_zz_0_xxy_1[i] * fi_abcd_0 + g_0_zz_0_xxyy_0[i] * pb_y + g_0_zz_0_xxyy_1[i] * wp_y[i];

        g_0_yzz_0_xxyz_0[i] = g_0_zz_0_xxz_1[i] * fi_abcd_0 + g_0_zz_0_xxyz_0[i] * pb_y + g_0_zz_0_xxyz_1[i] * wp_y[i];

        g_0_yzz_0_xxzz_0[i] = g_0_zz_0_xxzz_0[i] * pb_y + g_0_zz_0_xxzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyy_0[i] = 3.0 * g_0_zz_0_xyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyy_0[i] * pb_y + g_0_zz_0_xyyy_1[i] * wp_y[i];

        g_0_yzz_0_xyyz_0[i] = 2.0 * g_0_zz_0_xyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyz_0[i] * pb_y + g_0_zz_0_xyyz_1[i] * wp_y[i];

        g_0_yzz_0_xyzz_0[i] = g_0_zz_0_xzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzz_0[i] * pb_y + g_0_zz_0_xyzz_1[i] * wp_y[i];

        g_0_yzz_0_xzzz_0[i] = g_0_zz_0_xzzz_0[i] * pb_y + g_0_zz_0_xzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyy_0[i] = 4.0 * g_0_zz_0_yyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyy_0[i] * pb_y + g_0_zz_0_yyyy_1[i] * wp_y[i];

        g_0_yzz_0_yyyz_0[i] = 3.0 * g_0_zz_0_yyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyz_0[i] * pb_y + g_0_zz_0_yyyz_1[i] * wp_y[i];

        g_0_yzz_0_yyzz_0[i] = 2.0 * g_0_zz_0_yzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzz_0[i] * pb_y + g_0_zz_0_yyzz_1[i] * wp_y[i];

        g_0_yzz_0_yzzz_0[i] = g_0_zz_0_zzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzz_0[i] * pb_y + g_0_zz_0_yzzz_1[i] * wp_y[i];

        g_0_yzz_0_zzzz_0[i] = g_0_zz_0_zzzz_0[i] * pb_y + g_0_zz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 135-150 components of targeted buffer : prim_buffer_0_sfsg

    auto g_0_zzz_0_xxxx_0 = prim_buffer_0_sfsg[135];

    auto g_0_zzz_0_xxxy_0 = prim_buffer_0_sfsg[136];

    auto g_0_zzz_0_xxxz_0 = prim_buffer_0_sfsg[137];

    auto g_0_zzz_0_xxyy_0 = prim_buffer_0_sfsg[138];

    auto g_0_zzz_0_xxyz_0 = prim_buffer_0_sfsg[139];

    auto g_0_zzz_0_xxzz_0 = prim_buffer_0_sfsg[140];

    auto g_0_zzz_0_xyyy_0 = prim_buffer_0_sfsg[141];

    auto g_0_zzz_0_xyyz_0 = prim_buffer_0_sfsg[142];

    auto g_0_zzz_0_xyzz_0 = prim_buffer_0_sfsg[143];

    auto g_0_zzz_0_xzzz_0 = prim_buffer_0_sfsg[144];

    auto g_0_zzz_0_yyyy_0 = prim_buffer_0_sfsg[145];

    auto g_0_zzz_0_yyyz_0 = prim_buffer_0_sfsg[146];

    auto g_0_zzz_0_yyzz_0 = prim_buffer_0_sfsg[147];

    auto g_0_zzz_0_yzzz_0 = prim_buffer_0_sfsg[148];

    auto g_0_zzz_0_zzzz_0 = prim_buffer_0_sfsg[149];

    #pragma omp simd aligned(g_0_z_0_xxxx_0, g_0_z_0_xxxx_1, g_0_z_0_xxxy_0, g_0_z_0_xxxy_1, g_0_z_0_xxxz_0, g_0_z_0_xxxz_1, g_0_z_0_xxyy_0, g_0_z_0_xxyy_1, g_0_z_0_xxyz_0, g_0_z_0_xxyz_1, g_0_z_0_xxzz_0, g_0_z_0_xxzz_1, g_0_z_0_xyyy_0, g_0_z_0_xyyy_1, g_0_z_0_xyyz_0, g_0_z_0_xyyz_1, g_0_z_0_xyzz_0, g_0_z_0_xyzz_1, g_0_z_0_xzzz_0, g_0_z_0_xzzz_1, g_0_z_0_yyyy_0, g_0_z_0_yyyy_1, g_0_z_0_yyyz_0, g_0_z_0_yyyz_1, g_0_z_0_yyzz_0, g_0_z_0_yyzz_1, g_0_z_0_yzzz_0, g_0_z_0_yzzz_1, g_0_z_0_zzzz_0, g_0_z_0_zzzz_1, g_0_zz_0_xxx_1, g_0_zz_0_xxxx_0, g_0_zz_0_xxxx_1, g_0_zz_0_xxxy_0, g_0_zz_0_xxxy_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxy_1, g_0_zz_0_xxyy_0, g_0_zz_0_xxyy_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyy_1, g_0_zz_0_xyyy_0, g_0_zz_0_xyyy_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyy_1, g_0_zz_0_yyyy_0, g_0_zz_0_yyyy_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, g_0_zzz_0_xxxx_0, g_0_zzz_0_xxxy_0, g_0_zzz_0_xxxz_0, g_0_zzz_0_xxyy_0, g_0_zzz_0_xxyz_0, g_0_zzz_0_xxzz_0, g_0_zzz_0_xyyy_0, g_0_zzz_0_xyyz_0, g_0_zzz_0_xyzz_0, g_0_zzz_0_xzzz_0, g_0_zzz_0_yyyy_0, g_0_zzz_0_yyyz_0, g_0_zzz_0_yyzz_0, g_0_zzz_0_yzzz_0, g_0_zzz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxxx_0[i] = 2.0 * g_0_z_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxx_1[i] * fti_ab_0 + g_0_zz_0_xxxx_0[i] * pb_z + g_0_zz_0_xxxx_1[i] * wp_z[i];

        g_0_zzz_0_xxxy_0[i] = 2.0 * g_0_z_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxy_1[i] * fti_ab_0 + g_0_zz_0_xxxy_0[i] * pb_z + g_0_zz_0_xxxy_1[i] * wp_z[i];

        g_0_zzz_0_xxxz_0[i] = 2.0 * g_0_z_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxz_1[i] * fti_ab_0 + g_0_zz_0_xxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxz_0[i] * pb_z + g_0_zz_0_xxxz_1[i] * wp_z[i];

        g_0_zzz_0_xxyy_0[i] = 2.0 * g_0_z_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyy_1[i] * fti_ab_0 + g_0_zz_0_xxyy_0[i] * pb_z + g_0_zz_0_xxyy_1[i] * wp_z[i];

        g_0_zzz_0_xxyz_0[i] = 2.0 * g_0_z_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyz_1[i] * fti_ab_0 + g_0_zz_0_xxy_1[i] * fi_abcd_0 + g_0_zz_0_xxyz_0[i] * pb_z + g_0_zz_0_xxyz_1[i] * wp_z[i];

        g_0_zzz_0_xxzz_0[i] = 2.0 * g_0_z_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxz_1[i] * fi_abcd_0 + g_0_zz_0_xxzz_0[i] * pb_z + g_0_zz_0_xxzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyy_0[i] = 2.0 * g_0_z_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyy_1[i] * fti_ab_0 + g_0_zz_0_xyyy_0[i] * pb_z + g_0_zz_0_xyyy_1[i] * wp_z[i];

        g_0_zzz_0_xyyz_0[i] = 2.0 * g_0_z_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyz_1[i] * fti_ab_0 + g_0_zz_0_xyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyz_0[i] * pb_z + g_0_zz_0_xyyz_1[i] * wp_z[i];

        g_0_zzz_0_xyzz_0[i] = 2.0 * g_0_z_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xyz_1[i] * fi_abcd_0 + g_0_zz_0_xyzz_0[i] * pb_z + g_0_zz_0_xyzz_1[i] * wp_z[i];

        g_0_zzz_0_xzzz_0[i] = 2.0 * g_0_z_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzz_0[i] * pb_z + g_0_zz_0_xzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyy_0[i] = 2.0 * g_0_z_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyy_1[i] * fti_ab_0 + g_0_zz_0_yyyy_0[i] * pb_z + g_0_zz_0_yyyy_1[i] * wp_z[i];

        g_0_zzz_0_yyyz_0[i] = 2.0 * g_0_z_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyz_1[i] * fti_ab_0 + g_0_zz_0_yyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyz_0[i] * pb_z + g_0_zz_0_yyyz_1[i] * wp_z[i];

        g_0_zzz_0_yyzz_0[i] = 2.0 * g_0_z_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_yyz_1[i] * fi_abcd_0 + g_0_zz_0_yyzz_0[i] * pb_z + g_0_zz_0_yyzz_1[i] * wp_z[i];

        g_0_zzz_0_yzzz_0[i] = 2.0 * g_0_z_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_yzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzz_0[i] * pb_z + g_0_zz_0_yzzz_1[i] * wp_z[i];

        g_0_zzz_0_zzzz_0[i] = 2.0 * g_0_z_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_zzz_1[i] * fi_abcd_0 + g_0_zz_0_zzzz_0[i] * pb_z + g_0_zz_0_zzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


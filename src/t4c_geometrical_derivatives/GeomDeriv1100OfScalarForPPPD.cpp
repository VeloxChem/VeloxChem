#include "GeomDeriv1100OfScalarForPPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pppd_0(CSimdArray<double>& buffer_1100_pppd,
                     const CSimdArray<double>& buffer_sspd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_dspd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sspd

    auto g_0_0_x_xx = buffer_sspd[0];

    auto g_0_0_x_xy = buffer_sspd[1];

    auto g_0_0_x_xz = buffer_sspd[2];

    auto g_0_0_x_yy = buffer_sspd[3];

    auto g_0_0_x_yz = buffer_sspd[4];

    auto g_0_0_x_zz = buffer_sspd[5];

    auto g_0_0_y_xx = buffer_sspd[6];

    auto g_0_0_y_xy = buffer_sspd[7];

    auto g_0_0_y_xz = buffer_sspd[8];

    auto g_0_0_y_yy = buffer_sspd[9];

    auto g_0_0_y_yz = buffer_sspd[10];

    auto g_0_0_y_zz = buffer_sspd[11];

    auto g_0_0_z_xx = buffer_sspd[12];

    auto g_0_0_z_xy = buffer_sspd[13];

    auto g_0_0_z_xz = buffer_sspd[14];

    auto g_0_0_z_yy = buffer_sspd[15];

    auto g_0_0_z_yz = buffer_sspd[16];

    auto g_0_0_z_zz = buffer_sspd[17];

    /// Set up components of auxilary buffer : buffer_sdpd

    auto g_0_xx_x_xx = buffer_sdpd[0];

    auto g_0_xx_x_xy = buffer_sdpd[1];

    auto g_0_xx_x_xz = buffer_sdpd[2];

    auto g_0_xx_x_yy = buffer_sdpd[3];

    auto g_0_xx_x_yz = buffer_sdpd[4];

    auto g_0_xx_x_zz = buffer_sdpd[5];

    auto g_0_xx_y_xx = buffer_sdpd[6];

    auto g_0_xx_y_xy = buffer_sdpd[7];

    auto g_0_xx_y_xz = buffer_sdpd[8];

    auto g_0_xx_y_yy = buffer_sdpd[9];

    auto g_0_xx_y_yz = buffer_sdpd[10];

    auto g_0_xx_y_zz = buffer_sdpd[11];

    auto g_0_xx_z_xx = buffer_sdpd[12];

    auto g_0_xx_z_xy = buffer_sdpd[13];

    auto g_0_xx_z_xz = buffer_sdpd[14];

    auto g_0_xx_z_yy = buffer_sdpd[15];

    auto g_0_xx_z_yz = buffer_sdpd[16];

    auto g_0_xx_z_zz = buffer_sdpd[17];

    auto g_0_xy_x_xx = buffer_sdpd[18];

    auto g_0_xy_x_xy = buffer_sdpd[19];

    auto g_0_xy_x_xz = buffer_sdpd[20];

    auto g_0_xy_x_yy = buffer_sdpd[21];

    auto g_0_xy_x_yz = buffer_sdpd[22];

    auto g_0_xy_x_zz = buffer_sdpd[23];

    auto g_0_xy_y_xx = buffer_sdpd[24];

    auto g_0_xy_y_xy = buffer_sdpd[25];

    auto g_0_xy_y_xz = buffer_sdpd[26];

    auto g_0_xy_y_yy = buffer_sdpd[27];

    auto g_0_xy_y_yz = buffer_sdpd[28];

    auto g_0_xy_y_zz = buffer_sdpd[29];

    auto g_0_xy_z_xx = buffer_sdpd[30];

    auto g_0_xy_z_xy = buffer_sdpd[31];

    auto g_0_xy_z_xz = buffer_sdpd[32];

    auto g_0_xy_z_yy = buffer_sdpd[33];

    auto g_0_xy_z_yz = buffer_sdpd[34];

    auto g_0_xy_z_zz = buffer_sdpd[35];

    auto g_0_xz_x_xx = buffer_sdpd[36];

    auto g_0_xz_x_xy = buffer_sdpd[37];

    auto g_0_xz_x_xz = buffer_sdpd[38];

    auto g_0_xz_x_yy = buffer_sdpd[39];

    auto g_0_xz_x_yz = buffer_sdpd[40];

    auto g_0_xz_x_zz = buffer_sdpd[41];

    auto g_0_xz_y_xx = buffer_sdpd[42];

    auto g_0_xz_y_xy = buffer_sdpd[43];

    auto g_0_xz_y_xz = buffer_sdpd[44];

    auto g_0_xz_y_yy = buffer_sdpd[45];

    auto g_0_xz_y_yz = buffer_sdpd[46];

    auto g_0_xz_y_zz = buffer_sdpd[47];

    auto g_0_xz_z_xx = buffer_sdpd[48];

    auto g_0_xz_z_xy = buffer_sdpd[49];

    auto g_0_xz_z_xz = buffer_sdpd[50];

    auto g_0_xz_z_yy = buffer_sdpd[51];

    auto g_0_xz_z_yz = buffer_sdpd[52];

    auto g_0_xz_z_zz = buffer_sdpd[53];

    auto g_0_yy_x_xx = buffer_sdpd[54];

    auto g_0_yy_x_xy = buffer_sdpd[55];

    auto g_0_yy_x_xz = buffer_sdpd[56];

    auto g_0_yy_x_yy = buffer_sdpd[57];

    auto g_0_yy_x_yz = buffer_sdpd[58];

    auto g_0_yy_x_zz = buffer_sdpd[59];

    auto g_0_yy_y_xx = buffer_sdpd[60];

    auto g_0_yy_y_xy = buffer_sdpd[61];

    auto g_0_yy_y_xz = buffer_sdpd[62];

    auto g_0_yy_y_yy = buffer_sdpd[63];

    auto g_0_yy_y_yz = buffer_sdpd[64];

    auto g_0_yy_y_zz = buffer_sdpd[65];

    auto g_0_yy_z_xx = buffer_sdpd[66];

    auto g_0_yy_z_xy = buffer_sdpd[67];

    auto g_0_yy_z_xz = buffer_sdpd[68];

    auto g_0_yy_z_yy = buffer_sdpd[69];

    auto g_0_yy_z_yz = buffer_sdpd[70];

    auto g_0_yy_z_zz = buffer_sdpd[71];

    auto g_0_yz_x_xx = buffer_sdpd[72];

    auto g_0_yz_x_xy = buffer_sdpd[73];

    auto g_0_yz_x_xz = buffer_sdpd[74];

    auto g_0_yz_x_yy = buffer_sdpd[75];

    auto g_0_yz_x_yz = buffer_sdpd[76];

    auto g_0_yz_x_zz = buffer_sdpd[77];

    auto g_0_yz_y_xx = buffer_sdpd[78];

    auto g_0_yz_y_xy = buffer_sdpd[79];

    auto g_0_yz_y_xz = buffer_sdpd[80];

    auto g_0_yz_y_yy = buffer_sdpd[81];

    auto g_0_yz_y_yz = buffer_sdpd[82];

    auto g_0_yz_y_zz = buffer_sdpd[83];

    auto g_0_yz_z_xx = buffer_sdpd[84];

    auto g_0_yz_z_xy = buffer_sdpd[85];

    auto g_0_yz_z_xz = buffer_sdpd[86];

    auto g_0_yz_z_yy = buffer_sdpd[87];

    auto g_0_yz_z_yz = buffer_sdpd[88];

    auto g_0_yz_z_zz = buffer_sdpd[89];

    auto g_0_zz_x_xx = buffer_sdpd[90];

    auto g_0_zz_x_xy = buffer_sdpd[91];

    auto g_0_zz_x_xz = buffer_sdpd[92];

    auto g_0_zz_x_yy = buffer_sdpd[93];

    auto g_0_zz_x_yz = buffer_sdpd[94];

    auto g_0_zz_x_zz = buffer_sdpd[95];

    auto g_0_zz_y_xx = buffer_sdpd[96];

    auto g_0_zz_y_xy = buffer_sdpd[97];

    auto g_0_zz_y_xz = buffer_sdpd[98];

    auto g_0_zz_y_yy = buffer_sdpd[99];

    auto g_0_zz_y_yz = buffer_sdpd[100];

    auto g_0_zz_y_zz = buffer_sdpd[101];

    auto g_0_zz_z_xx = buffer_sdpd[102];

    auto g_0_zz_z_xy = buffer_sdpd[103];

    auto g_0_zz_z_xz = buffer_sdpd[104];

    auto g_0_zz_z_yy = buffer_sdpd[105];

    auto g_0_zz_z_yz = buffer_sdpd[106];

    auto g_0_zz_z_zz = buffer_sdpd[107];

    /// Set up components of auxilary buffer : buffer_dspd

    auto g_xx_0_x_xx = buffer_dspd[0];

    auto g_xx_0_x_xy = buffer_dspd[1];

    auto g_xx_0_x_xz = buffer_dspd[2];

    auto g_xx_0_x_yy = buffer_dspd[3];

    auto g_xx_0_x_yz = buffer_dspd[4];

    auto g_xx_0_x_zz = buffer_dspd[5];

    auto g_xx_0_y_xx = buffer_dspd[6];

    auto g_xx_0_y_xy = buffer_dspd[7];

    auto g_xx_0_y_xz = buffer_dspd[8];

    auto g_xx_0_y_yy = buffer_dspd[9];

    auto g_xx_0_y_yz = buffer_dspd[10];

    auto g_xx_0_y_zz = buffer_dspd[11];

    auto g_xx_0_z_xx = buffer_dspd[12];

    auto g_xx_0_z_xy = buffer_dspd[13];

    auto g_xx_0_z_xz = buffer_dspd[14];

    auto g_xx_0_z_yy = buffer_dspd[15];

    auto g_xx_0_z_yz = buffer_dspd[16];

    auto g_xx_0_z_zz = buffer_dspd[17];

    auto g_xy_0_x_xx = buffer_dspd[18];

    auto g_xy_0_x_xy = buffer_dspd[19];

    auto g_xy_0_x_xz = buffer_dspd[20];

    auto g_xy_0_x_yy = buffer_dspd[21];

    auto g_xy_0_x_yz = buffer_dspd[22];

    auto g_xy_0_x_zz = buffer_dspd[23];

    auto g_xy_0_y_xx = buffer_dspd[24];

    auto g_xy_0_y_xy = buffer_dspd[25];

    auto g_xy_0_y_xz = buffer_dspd[26];

    auto g_xy_0_y_yy = buffer_dspd[27];

    auto g_xy_0_y_yz = buffer_dspd[28];

    auto g_xy_0_y_zz = buffer_dspd[29];

    auto g_xy_0_z_xx = buffer_dspd[30];

    auto g_xy_0_z_xy = buffer_dspd[31];

    auto g_xy_0_z_xz = buffer_dspd[32];

    auto g_xy_0_z_yy = buffer_dspd[33];

    auto g_xy_0_z_yz = buffer_dspd[34];

    auto g_xy_0_z_zz = buffer_dspd[35];

    auto g_xz_0_x_xx = buffer_dspd[36];

    auto g_xz_0_x_xy = buffer_dspd[37];

    auto g_xz_0_x_xz = buffer_dspd[38];

    auto g_xz_0_x_yy = buffer_dspd[39];

    auto g_xz_0_x_yz = buffer_dspd[40];

    auto g_xz_0_x_zz = buffer_dspd[41];

    auto g_xz_0_y_xx = buffer_dspd[42];

    auto g_xz_0_y_xy = buffer_dspd[43];

    auto g_xz_0_y_xz = buffer_dspd[44];

    auto g_xz_0_y_yy = buffer_dspd[45];

    auto g_xz_0_y_yz = buffer_dspd[46];

    auto g_xz_0_y_zz = buffer_dspd[47];

    auto g_xz_0_z_xx = buffer_dspd[48];

    auto g_xz_0_z_xy = buffer_dspd[49];

    auto g_xz_0_z_xz = buffer_dspd[50];

    auto g_xz_0_z_yy = buffer_dspd[51];

    auto g_xz_0_z_yz = buffer_dspd[52];

    auto g_xz_0_z_zz = buffer_dspd[53];

    auto g_yy_0_x_xx = buffer_dspd[54];

    auto g_yy_0_x_xy = buffer_dspd[55];

    auto g_yy_0_x_xz = buffer_dspd[56];

    auto g_yy_0_x_yy = buffer_dspd[57];

    auto g_yy_0_x_yz = buffer_dspd[58];

    auto g_yy_0_x_zz = buffer_dspd[59];

    auto g_yy_0_y_xx = buffer_dspd[60];

    auto g_yy_0_y_xy = buffer_dspd[61];

    auto g_yy_0_y_xz = buffer_dspd[62];

    auto g_yy_0_y_yy = buffer_dspd[63];

    auto g_yy_0_y_yz = buffer_dspd[64];

    auto g_yy_0_y_zz = buffer_dspd[65];

    auto g_yy_0_z_xx = buffer_dspd[66];

    auto g_yy_0_z_xy = buffer_dspd[67];

    auto g_yy_0_z_xz = buffer_dspd[68];

    auto g_yy_0_z_yy = buffer_dspd[69];

    auto g_yy_0_z_yz = buffer_dspd[70];

    auto g_yy_0_z_zz = buffer_dspd[71];

    auto g_yz_0_x_xx = buffer_dspd[72];

    auto g_yz_0_x_xy = buffer_dspd[73];

    auto g_yz_0_x_xz = buffer_dspd[74];

    auto g_yz_0_x_yy = buffer_dspd[75];

    auto g_yz_0_x_yz = buffer_dspd[76];

    auto g_yz_0_x_zz = buffer_dspd[77];

    auto g_yz_0_y_xx = buffer_dspd[78];

    auto g_yz_0_y_xy = buffer_dspd[79];

    auto g_yz_0_y_xz = buffer_dspd[80];

    auto g_yz_0_y_yy = buffer_dspd[81];

    auto g_yz_0_y_yz = buffer_dspd[82];

    auto g_yz_0_y_zz = buffer_dspd[83];

    auto g_yz_0_z_xx = buffer_dspd[84];

    auto g_yz_0_z_xy = buffer_dspd[85];

    auto g_yz_0_z_xz = buffer_dspd[86];

    auto g_yz_0_z_yy = buffer_dspd[87];

    auto g_yz_0_z_yz = buffer_dspd[88];

    auto g_yz_0_z_zz = buffer_dspd[89];

    auto g_zz_0_x_xx = buffer_dspd[90];

    auto g_zz_0_x_xy = buffer_dspd[91];

    auto g_zz_0_x_xz = buffer_dspd[92];

    auto g_zz_0_x_yy = buffer_dspd[93];

    auto g_zz_0_x_yz = buffer_dspd[94];

    auto g_zz_0_x_zz = buffer_dspd[95];

    auto g_zz_0_y_xx = buffer_dspd[96];

    auto g_zz_0_y_xy = buffer_dspd[97];

    auto g_zz_0_y_xz = buffer_dspd[98];

    auto g_zz_0_y_yy = buffer_dspd[99];

    auto g_zz_0_y_yz = buffer_dspd[100];

    auto g_zz_0_y_zz = buffer_dspd[101];

    auto g_zz_0_z_xx = buffer_dspd[102];

    auto g_zz_0_z_xy = buffer_dspd[103];

    auto g_zz_0_z_xz = buffer_dspd[104];

    auto g_zz_0_z_yy = buffer_dspd[105];

    auto g_zz_0_z_yz = buffer_dspd[106];

    auto g_zz_0_z_zz = buffer_dspd[107];

    /// Set up components of auxilary buffer : buffer_ddpd

    auto g_xx_xx_x_xx = buffer_ddpd[0];

    auto g_xx_xx_x_xy = buffer_ddpd[1];

    auto g_xx_xx_x_xz = buffer_ddpd[2];

    auto g_xx_xx_x_yy = buffer_ddpd[3];

    auto g_xx_xx_x_yz = buffer_ddpd[4];

    auto g_xx_xx_x_zz = buffer_ddpd[5];

    auto g_xx_xx_y_xx = buffer_ddpd[6];

    auto g_xx_xx_y_xy = buffer_ddpd[7];

    auto g_xx_xx_y_xz = buffer_ddpd[8];

    auto g_xx_xx_y_yy = buffer_ddpd[9];

    auto g_xx_xx_y_yz = buffer_ddpd[10];

    auto g_xx_xx_y_zz = buffer_ddpd[11];

    auto g_xx_xx_z_xx = buffer_ddpd[12];

    auto g_xx_xx_z_xy = buffer_ddpd[13];

    auto g_xx_xx_z_xz = buffer_ddpd[14];

    auto g_xx_xx_z_yy = buffer_ddpd[15];

    auto g_xx_xx_z_yz = buffer_ddpd[16];

    auto g_xx_xx_z_zz = buffer_ddpd[17];

    auto g_xx_xy_x_xx = buffer_ddpd[18];

    auto g_xx_xy_x_xy = buffer_ddpd[19];

    auto g_xx_xy_x_xz = buffer_ddpd[20];

    auto g_xx_xy_x_yy = buffer_ddpd[21];

    auto g_xx_xy_x_yz = buffer_ddpd[22];

    auto g_xx_xy_x_zz = buffer_ddpd[23];

    auto g_xx_xy_y_xx = buffer_ddpd[24];

    auto g_xx_xy_y_xy = buffer_ddpd[25];

    auto g_xx_xy_y_xz = buffer_ddpd[26];

    auto g_xx_xy_y_yy = buffer_ddpd[27];

    auto g_xx_xy_y_yz = buffer_ddpd[28];

    auto g_xx_xy_y_zz = buffer_ddpd[29];

    auto g_xx_xy_z_xx = buffer_ddpd[30];

    auto g_xx_xy_z_xy = buffer_ddpd[31];

    auto g_xx_xy_z_xz = buffer_ddpd[32];

    auto g_xx_xy_z_yy = buffer_ddpd[33];

    auto g_xx_xy_z_yz = buffer_ddpd[34];

    auto g_xx_xy_z_zz = buffer_ddpd[35];

    auto g_xx_xz_x_xx = buffer_ddpd[36];

    auto g_xx_xz_x_xy = buffer_ddpd[37];

    auto g_xx_xz_x_xz = buffer_ddpd[38];

    auto g_xx_xz_x_yy = buffer_ddpd[39];

    auto g_xx_xz_x_yz = buffer_ddpd[40];

    auto g_xx_xz_x_zz = buffer_ddpd[41];

    auto g_xx_xz_y_xx = buffer_ddpd[42];

    auto g_xx_xz_y_xy = buffer_ddpd[43];

    auto g_xx_xz_y_xz = buffer_ddpd[44];

    auto g_xx_xz_y_yy = buffer_ddpd[45];

    auto g_xx_xz_y_yz = buffer_ddpd[46];

    auto g_xx_xz_y_zz = buffer_ddpd[47];

    auto g_xx_xz_z_xx = buffer_ddpd[48];

    auto g_xx_xz_z_xy = buffer_ddpd[49];

    auto g_xx_xz_z_xz = buffer_ddpd[50];

    auto g_xx_xz_z_yy = buffer_ddpd[51];

    auto g_xx_xz_z_yz = buffer_ddpd[52];

    auto g_xx_xz_z_zz = buffer_ddpd[53];

    auto g_xx_yy_x_xx = buffer_ddpd[54];

    auto g_xx_yy_x_xy = buffer_ddpd[55];

    auto g_xx_yy_x_xz = buffer_ddpd[56];

    auto g_xx_yy_x_yy = buffer_ddpd[57];

    auto g_xx_yy_x_yz = buffer_ddpd[58];

    auto g_xx_yy_x_zz = buffer_ddpd[59];

    auto g_xx_yy_y_xx = buffer_ddpd[60];

    auto g_xx_yy_y_xy = buffer_ddpd[61];

    auto g_xx_yy_y_xz = buffer_ddpd[62];

    auto g_xx_yy_y_yy = buffer_ddpd[63];

    auto g_xx_yy_y_yz = buffer_ddpd[64];

    auto g_xx_yy_y_zz = buffer_ddpd[65];

    auto g_xx_yy_z_xx = buffer_ddpd[66];

    auto g_xx_yy_z_xy = buffer_ddpd[67];

    auto g_xx_yy_z_xz = buffer_ddpd[68];

    auto g_xx_yy_z_yy = buffer_ddpd[69];

    auto g_xx_yy_z_yz = buffer_ddpd[70];

    auto g_xx_yy_z_zz = buffer_ddpd[71];

    auto g_xx_yz_x_xx = buffer_ddpd[72];

    auto g_xx_yz_x_xy = buffer_ddpd[73];

    auto g_xx_yz_x_xz = buffer_ddpd[74];

    auto g_xx_yz_x_yy = buffer_ddpd[75];

    auto g_xx_yz_x_yz = buffer_ddpd[76];

    auto g_xx_yz_x_zz = buffer_ddpd[77];

    auto g_xx_yz_y_xx = buffer_ddpd[78];

    auto g_xx_yz_y_xy = buffer_ddpd[79];

    auto g_xx_yz_y_xz = buffer_ddpd[80];

    auto g_xx_yz_y_yy = buffer_ddpd[81];

    auto g_xx_yz_y_yz = buffer_ddpd[82];

    auto g_xx_yz_y_zz = buffer_ddpd[83];

    auto g_xx_yz_z_xx = buffer_ddpd[84];

    auto g_xx_yz_z_xy = buffer_ddpd[85];

    auto g_xx_yz_z_xz = buffer_ddpd[86];

    auto g_xx_yz_z_yy = buffer_ddpd[87];

    auto g_xx_yz_z_yz = buffer_ddpd[88];

    auto g_xx_yz_z_zz = buffer_ddpd[89];

    auto g_xx_zz_x_xx = buffer_ddpd[90];

    auto g_xx_zz_x_xy = buffer_ddpd[91];

    auto g_xx_zz_x_xz = buffer_ddpd[92];

    auto g_xx_zz_x_yy = buffer_ddpd[93];

    auto g_xx_zz_x_yz = buffer_ddpd[94];

    auto g_xx_zz_x_zz = buffer_ddpd[95];

    auto g_xx_zz_y_xx = buffer_ddpd[96];

    auto g_xx_zz_y_xy = buffer_ddpd[97];

    auto g_xx_zz_y_xz = buffer_ddpd[98];

    auto g_xx_zz_y_yy = buffer_ddpd[99];

    auto g_xx_zz_y_yz = buffer_ddpd[100];

    auto g_xx_zz_y_zz = buffer_ddpd[101];

    auto g_xx_zz_z_xx = buffer_ddpd[102];

    auto g_xx_zz_z_xy = buffer_ddpd[103];

    auto g_xx_zz_z_xz = buffer_ddpd[104];

    auto g_xx_zz_z_yy = buffer_ddpd[105];

    auto g_xx_zz_z_yz = buffer_ddpd[106];

    auto g_xx_zz_z_zz = buffer_ddpd[107];

    auto g_xy_xx_x_xx = buffer_ddpd[108];

    auto g_xy_xx_x_xy = buffer_ddpd[109];

    auto g_xy_xx_x_xz = buffer_ddpd[110];

    auto g_xy_xx_x_yy = buffer_ddpd[111];

    auto g_xy_xx_x_yz = buffer_ddpd[112];

    auto g_xy_xx_x_zz = buffer_ddpd[113];

    auto g_xy_xx_y_xx = buffer_ddpd[114];

    auto g_xy_xx_y_xy = buffer_ddpd[115];

    auto g_xy_xx_y_xz = buffer_ddpd[116];

    auto g_xy_xx_y_yy = buffer_ddpd[117];

    auto g_xy_xx_y_yz = buffer_ddpd[118];

    auto g_xy_xx_y_zz = buffer_ddpd[119];

    auto g_xy_xx_z_xx = buffer_ddpd[120];

    auto g_xy_xx_z_xy = buffer_ddpd[121];

    auto g_xy_xx_z_xz = buffer_ddpd[122];

    auto g_xy_xx_z_yy = buffer_ddpd[123];

    auto g_xy_xx_z_yz = buffer_ddpd[124];

    auto g_xy_xx_z_zz = buffer_ddpd[125];

    auto g_xy_xy_x_xx = buffer_ddpd[126];

    auto g_xy_xy_x_xy = buffer_ddpd[127];

    auto g_xy_xy_x_xz = buffer_ddpd[128];

    auto g_xy_xy_x_yy = buffer_ddpd[129];

    auto g_xy_xy_x_yz = buffer_ddpd[130];

    auto g_xy_xy_x_zz = buffer_ddpd[131];

    auto g_xy_xy_y_xx = buffer_ddpd[132];

    auto g_xy_xy_y_xy = buffer_ddpd[133];

    auto g_xy_xy_y_xz = buffer_ddpd[134];

    auto g_xy_xy_y_yy = buffer_ddpd[135];

    auto g_xy_xy_y_yz = buffer_ddpd[136];

    auto g_xy_xy_y_zz = buffer_ddpd[137];

    auto g_xy_xy_z_xx = buffer_ddpd[138];

    auto g_xy_xy_z_xy = buffer_ddpd[139];

    auto g_xy_xy_z_xz = buffer_ddpd[140];

    auto g_xy_xy_z_yy = buffer_ddpd[141];

    auto g_xy_xy_z_yz = buffer_ddpd[142];

    auto g_xy_xy_z_zz = buffer_ddpd[143];

    auto g_xy_xz_x_xx = buffer_ddpd[144];

    auto g_xy_xz_x_xy = buffer_ddpd[145];

    auto g_xy_xz_x_xz = buffer_ddpd[146];

    auto g_xy_xz_x_yy = buffer_ddpd[147];

    auto g_xy_xz_x_yz = buffer_ddpd[148];

    auto g_xy_xz_x_zz = buffer_ddpd[149];

    auto g_xy_xz_y_xx = buffer_ddpd[150];

    auto g_xy_xz_y_xy = buffer_ddpd[151];

    auto g_xy_xz_y_xz = buffer_ddpd[152];

    auto g_xy_xz_y_yy = buffer_ddpd[153];

    auto g_xy_xz_y_yz = buffer_ddpd[154];

    auto g_xy_xz_y_zz = buffer_ddpd[155];

    auto g_xy_xz_z_xx = buffer_ddpd[156];

    auto g_xy_xz_z_xy = buffer_ddpd[157];

    auto g_xy_xz_z_xz = buffer_ddpd[158];

    auto g_xy_xz_z_yy = buffer_ddpd[159];

    auto g_xy_xz_z_yz = buffer_ddpd[160];

    auto g_xy_xz_z_zz = buffer_ddpd[161];

    auto g_xy_yy_x_xx = buffer_ddpd[162];

    auto g_xy_yy_x_xy = buffer_ddpd[163];

    auto g_xy_yy_x_xz = buffer_ddpd[164];

    auto g_xy_yy_x_yy = buffer_ddpd[165];

    auto g_xy_yy_x_yz = buffer_ddpd[166];

    auto g_xy_yy_x_zz = buffer_ddpd[167];

    auto g_xy_yy_y_xx = buffer_ddpd[168];

    auto g_xy_yy_y_xy = buffer_ddpd[169];

    auto g_xy_yy_y_xz = buffer_ddpd[170];

    auto g_xy_yy_y_yy = buffer_ddpd[171];

    auto g_xy_yy_y_yz = buffer_ddpd[172];

    auto g_xy_yy_y_zz = buffer_ddpd[173];

    auto g_xy_yy_z_xx = buffer_ddpd[174];

    auto g_xy_yy_z_xy = buffer_ddpd[175];

    auto g_xy_yy_z_xz = buffer_ddpd[176];

    auto g_xy_yy_z_yy = buffer_ddpd[177];

    auto g_xy_yy_z_yz = buffer_ddpd[178];

    auto g_xy_yy_z_zz = buffer_ddpd[179];

    auto g_xy_yz_x_xx = buffer_ddpd[180];

    auto g_xy_yz_x_xy = buffer_ddpd[181];

    auto g_xy_yz_x_xz = buffer_ddpd[182];

    auto g_xy_yz_x_yy = buffer_ddpd[183];

    auto g_xy_yz_x_yz = buffer_ddpd[184];

    auto g_xy_yz_x_zz = buffer_ddpd[185];

    auto g_xy_yz_y_xx = buffer_ddpd[186];

    auto g_xy_yz_y_xy = buffer_ddpd[187];

    auto g_xy_yz_y_xz = buffer_ddpd[188];

    auto g_xy_yz_y_yy = buffer_ddpd[189];

    auto g_xy_yz_y_yz = buffer_ddpd[190];

    auto g_xy_yz_y_zz = buffer_ddpd[191];

    auto g_xy_yz_z_xx = buffer_ddpd[192];

    auto g_xy_yz_z_xy = buffer_ddpd[193];

    auto g_xy_yz_z_xz = buffer_ddpd[194];

    auto g_xy_yz_z_yy = buffer_ddpd[195];

    auto g_xy_yz_z_yz = buffer_ddpd[196];

    auto g_xy_yz_z_zz = buffer_ddpd[197];

    auto g_xy_zz_x_xx = buffer_ddpd[198];

    auto g_xy_zz_x_xy = buffer_ddpd[199];

    auto g_xy_zz_x_xz = buffer_ddpd[200];

    auto g_xy_zz_x_yy = buffer_ddpd[201];

    auto g_xy_zz_x_yz = buffer_ddpd[202];

    auto g_xy_zz_x_zz = buffer_ddpd[203];

    auto g_xy_zz_y_xx = buffer_ddpd[204];

    auto g_xy_zz_y_xy = buffer_ddpd[205];

    auto g_xy_zz_y_xz = buffer_ddpd[206];

    auto g_xy_zz_y_yy = buffer_ddpd[207];

    auto g_xy_zz_y_yz = buffer_ddpd[208];

    auto g_xy_zz_y_zz = buffer_ddpd[209];

    auto g_xy_zz_z_xx = buffer_ddpd[210];

    auto g_xy_zz_z_xy = buffer_ddpd[211];

    auto g_xy_zz_z_xz = buffer_ddpd[212];

    auto g_xy_zz_z_yy = buffer_ddpd[213];

    auto g_xy_zz_z_yz = buffer_ddpd[214];

    auto g_xy_zz_z_zz = buffer_ddpd[215];

    auto g_xz_xx_x_xx = buffer_ddpd[216];

    auto g_xz_xx_x_xy = buffer_ddpd[217];

    auto g_xz_xx_x_xz = buffer_ddpd[218];

    auto g_xz_xx_x_yy = buffer_ddpd[219];

    auto g_xz_xx_x_yz = buffer_ddpd[220];

    auto g_xz_xx_x_zz = buffer_ddpd[221];

    auto g_xz_xx_y_xx = buffer_ddpd[222];

    auto g_xz_xx_y_xy = buffer_ddpd[223];

    auto g_xz_xx_y_xz = buffer_ddpd[224];

    auto g_xz_xx_y_yy = buffer_ddpd[225];

    auto g_xz_xx_y_yz = buffer_ddpd[226];

    auto g_xz_xx_y_zz = buffer_ddpd[227];

    auto g_xz_xx_z_xx = buffer_ddpd[228];

    auto g_xz_xx_z_xy = buffer_ddpd[229];

    auto g_xz_xx_z_xz = buffer_ddpd[230];

    auto g_xz_xx_z_yy = buffer_ddpd[231];

    auto g_xz_xx_z_yz = buffer_ddpd[232];

    auto g_xz_xx_z_zz = buffer_ddpd[233];

    auto g_xz_xy_x_xx = buffer_ddpd[234];

    auto g_xz_xy_x_xy = buffer_ddpd[235];

    auto g_xz_xy_x_xz = buffer_ddpd[236];

    auto g_xz_xy_x_yy = buffer_ddpd[237];

    auto g_xz_xy_x_yz = buffer_ddpd[238];

    auto g_xz_xy_x_zz = buffer_ddpd[239];

    auto g_xz_xy_y_xx = buffer_ddpd[240];

    auto g_xz_xy_y_xy = buffer_ddpd[241];

    auto g_xz_xy_y_xz = buffer_ddpd[242];

    auto g_xz_xy_y_yy = buffer_ddpd[243];

    auto g_xz_xy_y_yz = buffer_ddpd[244];

    auto g_xz_xy_y_zz = buffer_ddpd[245];

    auto g_xz_xy_z_xx = buffer_ddpd[246];

    auto g_xz_xy_z_xy = buffer_ddpd[247];

    auto g_xz_xy_z_xz = buffer_ddpd[248];

    auto g_xz_xy_z_yy = buffer_ddpd[249];

    auto g_xz_xy_z_yz = buffer_ddpd[250];

    auto g_xz_xy_z_zz = buffer_ddpd[251];

    auto g_xz_xz_x_xx = buffer_ddpd[252];

    auto g_xz_xz_x_xy = buffer_ddpd[253];

    auto g_xz_xz_x_xz = buffer_ddpd[254];

    auto g_xz_xz_x_yy = buffer_ddpd[255];

    auto g_xz_xz_x_yz = buffer_ddpd[256];

    auto g_xz_xz_x_zz = buffer_ddpd[257];

    auto g_xz_xz_y_xx = buffer_ddpd[258];

    auto g_xz_xz_y_xy = buffer_ddpd[259];

    auto g_xz_xz_y_xz = buffer_ddpd[260];

    auto g_xz_xz_y_yy = buffer_ddpd[261];

    auto g_xz_xz_y_yz = buffer_ddpd[262];

    auto g_xz_xz_y_zz = buffer_ddpd[263];

    auto g_xz_xz_z_xx = buffer_ddpd[264];

    auto g_xz_xz_z_xy = buffer_ddpd[265];

    auto g_xz_xz_z_xz = buffer_ddpd[266];

    auto g_xz_xz_z_yy = buffer_ddpd[267];

    auto g_xz_xz_z_yz = buffer_ddpd[268];

    auto g_xz_xz_z_zz = buffer_ddpd[269];

    auto g_xz_yy_x_xx = buffer_ddpd[270];

    auto g_xz_yy_x_xy = buffer_ddpd[271];

    auto g_xz_yy_x_xz = buffer_ddpd[272];

    auto g_xz_yy_x_yy = buffer_ddpd[273];

    auto g_xz_yy_x_yz = buffer_ddpd[274];

    auto g_xz_yy_x_zz = buffer_ddpd[275];

    auto g_xz_yy_y_xx = buffer_ddpd[276];

    auto g_xz_yy_y_xy = buffer_ddpd[277];

    auto g_xz_yy_y_xz = buffer_ddpd[278];

    auto g_xz_yy_y_yy = buffer_ddpd[279];

    auto g_xz_yy_y_yz = buffer_ddpd[280];

    auto g_xz_yy_y_zz = buffer_ddpd[281];

    auto g_xz_yy_z_xx = buffer_ddpd[282];

    auto g_xz_yy_z_xy = buffer_ddpd[283];

    auto g_xz_yy_z_xz = buffer_ddpd[284];

    auto g_xz_yy_z_yy = buffer_ddpd[285];

    auto g_xz_yy_z_yz = buffer_ddpd[286];

    auto g_xz_yy_z_zz = buffer_ddpd[287];

    auto g_xz_yz_x_xx = buffer_ddpd[288];

    auto g_xz_yz_x_xy = buffer_ddpd[289];

    auto g_xz_yz_x_xz = buffer_ddpd[290];

    auto g_xz_yz_x_yy = buffer_ddpd[291];

    auto g_xz_yz_x_yz = buffer_ddpd[292];

    auto g_xz_yz_x_zz = buffer_ddpd[293];

    auto g_xz_yz_y_xx = buffer_ddpd[294];

    auto g_xz_yz_y_xy = buffer_ddpd[295];

    auto g_xz_yz_y_xz = buffer_ddpd[296];

    auto g_xz_yz_y_yy = buffer_ddpd[297];

    auto g_xz_yz_y_yz = buffer_ddpd[298];

    auto g_xz_yz_y_zz = buffer_ddpd[299];

    auto g_xz_yz_z_xx = buffer_ddpd[300];

    auto g_xz_yz_z_xy = buffer_ddpd[301];

    auto g_xz_yz_z_xz = buffer_ddpd[302];

    auto g_xz_yz_z_yy = buffer_ddpd[303];

    auto g_xz_yz_z_yz = buffer_ddpd[304];

    auto g_xz_yz_z_zz = buffer_ddpd[305];

    auto g_xz_zz_x_xx = buffer_ddpd[306];

    auto g_xz_zz_x_xy = buffer_ddpd[307];

    auto g_xz_zz_x_xz = buffer_ddpd[308];

    auto g_xz_zz_x_yy = buffer_ddpd[309];

    auto g_xz_zz_x_yz = buffer_ddpd[310];

    auto g_xz_zz_x_zz = buffer_ddpd[311];

    auto g_xz_zz_y_xx = buffer_ddpd[312];

    auto g_xz_zz_y_xy = buffer_ddpd[313];

    auto g_xz_zz_y_xz = buffer_ddpd[314];

    auto g_xz_zz_y_yy = buffer_ddpd[315];

    auto g_xz_zz_y_yz = buffer_ddpd[316];

    auto g_xz_zz_y_zz = buffer_ddpd[317];

    auto g_xz_zz_z_xx = buffer_ddpd[318];

    auto g_xz_zz_z_xy = buffer_ddpd[319];

    auto g_xz_zz_z_xz = buffer_ddpd[320];

    auto g_xz_zz_z_yy = buffer_ddpd[321];

    auto g_xz_zz_z_yz = buffer_ddpd[322];

    auto g_xz_zz_z_zz = buffer_ddpd[323];

    auto g_yy_xx_x_xx = buffer_ddpd[324];

    auto g_yy_xx_x_xy = buffer_ddpd[325];

    auto g_yy_xx_x_xz = buffer_ddpd[326];

    auto g_yy_xx_x_yy = buffer_ddpd[327];

    auto g_yy_xx_x_yz = buffer_ddpd[328];

    auto g_yy_xx_x_zz = buffer_ddpd[329];

    auto g_yy_xx_y_xx = buffer_ddpd[330];

    auto g_yy_xx_y_xy = buffer_ddpd[331];

    auto g_yy_xx_y_xz = buffer_ddpd[332];

    auto g_yy_xx_y_yy = buffer_ddpd[333];

    auto g_yy_xx_y_yz = buffer_ddpd[334];

    auto g_yy_xx_y_zz = buffer_ddpd[335];

    auto g_yy_xx_z_xx = buffer_ddpd[336];

    auto g_yy_xx_z_xy = buffer_ddpd[337];

    auto g_yy_xx_z_xz = buffer_ddpd[338];

    auto g_yy_xx_z_yy = buffer_ddpd[339];

    auto g_yy_xx_z_yz = buffer_ddpd[340];

    auto g_yy_xx_z_zz = buffer_ddpd[341];

    auto g_yy_xy_x_xx = buffer_ddpd[342];

    auto g_yy_xy_x_xy = buffer_ddpd[343];

    auto g_yy_xy_x_xz = buffer_ddpd[344];

    auto g_yy_xy_x_yy = buffer_ddpd[345];

    auto g_yy_xy_x_yz = buffer_ddpd[346];

    auto g_yy_xy_x_zz = buffer_ddpd[347];

    auto g_yy_xy_y_xx = buffer_ddpd[348];

    auto g_yy_xy_y_xy = buffer_ddpd[349];

    auto g_yy_xy_y_xz = buffer_ddpd[350];

    auto g_yy_xy_y_yy = buffer_ddpd[351];

    auto g_yy_xy_y_yz = buffer_ddpd[352];

    auto g_yy_xy_y_zz = buffer_ddpd[353];

    auto g_yy_xy_z_xx = buffer_ddpd[354];

    auto g_yy_xy_z_xy = buffer_ddpd[355];

    auto g_yy_xy_z_xz = buffer_ddpd[356];

    auto g_yy_xy_z_yy = buffer_ddpd[357];

    auto g_yy_xy_z_yz = buffer_ddpd[358];

    auto g_yy_xy_z_zz = buffer_ddpd[359];

    auto g_yy_xz_x_xx = buffer_ddpd[360];

    auto g_yy_xz_x_xy = buffer_ddpd[361];

    auto g_yy_xz_x_xz = buffer_ddpd[362];

    auto g_yy_xz_x_yy = buffer_ddpd[363];

    auto g_yy_xz_x_yz = buffer_ddpd[364];

    auto g_yy_xz_x_zz = buffer_ddpd[365];

    auto g_yy_xz_y_xx = buffer_ddpd[366];

    auto g_yy_xz_y_xy = buffer_ddpd[367];

    auto g_yy_xz_y_xz = buffer_ddpd[368];

    auto g_yy_xz_y_yy = buffer_ddpd[369];

    auto g_yy_xz_y_yz = buffer_ddpd[370];

    auto g_yy_xz_y_zz = buffer_ddpd[371];

    auto g_yy_xz_z_xx = buffer_ddpd[372];

    auto g_yy_xz_z_xy = buffer_ddpd[373];

    auto g_yy_xz_z_xz = buffer_ddpd[374];

    auto g_yy_xz_z_yy = buffer_ddpd[375];

    auto g_yy_xz_z_yz = buffer_ddpd[376];

    auto g_yy_xz_z_zz = buffer_ddpd[377];

    auto g_yy_yy_x_xx = buffer_ddpd[378];

    auto g_yy_yy_x_xy = buffer_ddpd[379];

    auto g_yy_yy_x_xz = buffer_ddpd[380];

    auto g_yy_yy_x_yy = buffer_ddpd[381];

    auto g_yy_yy_x_yz = buffer_ddpd[382];

    auto g_yy_yy_x_zz = buffer_ddpd[383];

    auto g_yy_yy_y_xx = buffer_ddpd[384];

    auto g_yy_yy_y_xy = buffer_ddpd[385];

    auto g_yy_yy_y_xz = buffer_ddpd[386];

    auto g_yy_yy_y_yy = buffer_ddpd[387];

    auto g_yy_yy_y_yz = buffer_ddpd[388];

    auto g_yy_yy_y_zz = buffer_ddpd[389];

    auto g_yy_yy_z_xx = buffer_ddpd[390];

    auto g_yy_yy_z_xy = buffer_ddpd[391];

    auto g_yy_yy_z_xz = buffer_ddpd[392];

    auto g_yy_yy_z_yy = buffer_ddpd[393];

    auto g_yy_yy_z_yz = buffer_ddpd[394];

    auto g_yy_yy_z_zz = buffer_ddpd[395];

    auto g_yy_yz_x_xx = buffer_ddpd[396];

    auto g_yy_yz_x_xy = buffer_ddpd[397];

    auto g_yy_yz_x_xz = buffer_ddpd[398];

    auto g_yy_yz_x_yy = buffer_ddpd[399];

    auto g_yy_yz_x_yz = buffer_ddpd[400];

    auto g_yy_yz_x_zz = buffer_ddpd[401];

    auto g_yy_yz_y_xx = buffer_ddpd[402];

    auto g_yy_yz_y_xy = buffer_ddpd[403];

    auto g_yy_yz_y_xz = buffer_ddpd[404];

    auto g_yy_yz_y_yy = buffer_ddpd[405];

    auto g_yy_yz_y_yz = buffer_ddpd[406];

    auto g_yy_yz_y_zz = buffer_ddpd[407];

    auto g_yy_yz_z_xx = buffer_ddpd[408];

    auto g_yy_yz_z_xy = buffer_ddpd[409];

    auto g_yy_yz_z_xz = buffer_ddpd[410];

    auto g_yy_yz_z_yy = buffer_ddpd[411];

    auto g_yy_yz_z_yz = buffer_ddpd[412];

    auto g_yy_yz_z_zz = buffer_ddpd[413];

    auto g_yy_zz_x_xx = buffer_ddpd[414];

    auto g_yy_zz_x_xy = buffer_ddpd[415];

    auto g_yy_zz_x_xz = buffer_ddpd[416];

    auto g_yy_zz_x_yy = buffer_ddpd[417];

    auto g_yy_zz_x_yz = buffer_ddpd[418];

    auto g_yy_zz_x_zz = buffer_ddpd[419];

    auto g_yy_zz_y_xx = buffer_ddpd[420];

    auto g_yy_zz_y_xy = buffer_ddpd[421];

    auto g_yy_zz_y_xz = buffer_ddpd[422];

    auto g_yy_zz_y_yy = buffer_ddpd[423];

    auto g_yy_zz_y_yz = buffer_ddpd[424];

    auto g_yy_zz_y_zz = buffer_ddpd[425];

    auto g_yy_zz_z_xx = buffer_ddpd[426];

    auto g_yy_zz_z_xy = buffer_ddpd[427];

    auto g_yy_zz_z_xz = buffer_ddpd[428];

    auto g_yy_zz_z_yy = buffer_ddpd[429];

    auto g_yy_zz_z_yz = buffer_ddpd[430];

    auto g_yy_zz_z_zz = buffer_ddpd[431];

    auto g_yz_xx_x_xx = buffer_ddpd[432];

    auto g_yz_xx_x_xy = buffer_ddpd[433];

    auto g_yz_xx_x_xz = buffer_ddpd[434];

    auto g_yz_xx_x_yy = buffer_ddpd[435];

    auto g_yz_xx_x_yz = buffer_ddpd[436];

    auto g_yz_xx_x_zz = buffer_ddpd[437];

    auto g_yz_xx_y_xx = buffer_ddpd[438];

    auto g_yz_xx_y_xy = buffer_ddpd[439];

    auto g_yz_xx_y_xz = buffer_ddpd[440];

    auto g_yz_xx_y_yy = buffer_ddpd[441];

    auto g_yz_xx_y_yz = buffer_ddpd[442];

    auto g_yz_xx_y_zz = buffer_ddpd[443];

    auto g_yz_xx_z_xx = buffer_ddpd[444];

    auto g_yz_xx_z_xy = buffer_ddpd[445];

    auto g_yz_xx_z_xz = buffer_ddpd[446];

    auto g_yz_xx_z_yy = buffer_ddpd[447];

    auto g_yz_xx_z_yz = buffer_ddpd[448];

    auto g_yz_xx_z_zz = buffer_ddpd[449];

    auto g_yz_xy_x_xx = buffer_ddpd[450];

    auto g_yz_xy_x_xy = buffer_ddpd[451];

    auto g_yz_xy_x_xz = buffer_ddpd[452];

    auto g_yz_xy_x_yy = buffer_ddpd[453];

    auto g_yz_xy_x_yz = buffer_ddpd[454];

    auto g_yz_xy_x_zz = buffer_ddpd[455];

    auto g_yz_xy_y_xx = buffer_ddpd[456];

    auto g_yz_xy_y_xy = buffer_ddpd[457];

    auto g_yz_xy_y_xz = buffer_ddpd[458];

    auto g_yz_xy_y_yy = buffer_ddpd[459];

    auto g_yz_xy_y_yz = buffer_ddpd[460];

    auto g_yz_xy_y_zz = buffer_ddpd[461];

    auto g_yz_xy_z_xx = buffer_ddpd[462];

    auto g_yz_xy_z_xy = buffer_ddpd[463];

    auto g_yz_xy_z_xz = buffer_ddpd[464];

    auto g_yz_xy_z_yy = buffer_ddpd[465];

    auto g_yz_xy_z_yz = buffer_ddpd[466];

    auto g_yz_xy_z_zz = buffer_ddpd[467];

    auto g_yz_xz_x_xx = buffer_ddpd[468];

    auto g_yz_xz_x_xy = buffer_ddpd[469];

    auto g_yz_xz_x_xz = buffer_ddpd[470];

    auto g_yz_xz_x_yy = buffer_ddpd[471];

    auto g_yz_xz_x_yz = buffer_ddpd[472];

    auto g_yz_xz_x_zz = buffer_ddpd[473];

    auto g_yz_xz_y_xx = buffer_ddpd[474];

    auto g_yz_xz_y_xy = buffer_ddpd[475];

    auto g_yz_xz_y_xz = buffer_ddpd[476];

    auto g_yz_xz_y_yy = buffer_ddpd[477];

    auto g_yz_xz_y_yz = buffer_ddpd[478];

    auto g_yz_xz_y_zz = buffer_ddpd[479];

    auto g_yz_xz_z_xx = buffer_ddpd[480];

    auto g_yz_xz_z_xy = buffer_ddpd[481];

    auto g_yz_xz_z_xz = buffer_ddpd[482];

    auto g_yz_xz_z_yy = buffer_ddpd[483];

    auto g_yz_xz_z_yz = buffer_ddpd[484];

    auto g_yz_xz_z_zz = buffer_ddpd[485];

    auto g_yz_yy_x_xx = buffer_ddpd[486];

    auto g_yz_yy_x_xy = buffer_ddpd[487];

    auto g_yz_yy_x_xz = buffer_ddpd[488];

    auto g_yz_yy_x_yy = buffer_ddpd[489];

    auto g_yz_yy_x_yz = buffer_ddpd[490];

    auto g_yz_yy_x_zz = buffer_ddpd[491];

    auto g_yz_yy_y_xx = buffer_ddpd[492];

    auto g_yz_yy_y_xy = buffer_ddpd[493];

    auto g_yz_yy_y_xz = buffer_ddpd[494];

    auto g_yz_yy_y_yy = buffer_ddpd[495];

    auto g_yz_yy_y_yz = buffer_ddpd[496];

    auto g_yz_yy_y_zz = buffer_ddpd[497];

    auto g_yz_yy_z_xx = buffer_ddpd[498];

    auto g_yz_yy_z_xy = buffer_ddpd[499];

    auto g_yz_yy_z_xz = buffer_ddpd[500];

    auto g_yz_yy_z_yy = buffer_ddpd[501];

    auto g_yz_yy_z_yz = buffer_ddpd[502];

    auto g_yz_yy_z_zz = buffer_ddpd[503];

    auto g_yz_yz_x_xx = buffer_ddpd[504];

    auto g_yz_yz_x_xy = buffer_ddpd[505];

    auto g_yz_yz_x_xz = buffer_ddpd[506];

    auto g_yz_yz_x_yy = buffer_ddpd[507];

    auto g_yz_yz_x_yz = buffer_ddpd[508];

    auto g_yz_yz_x_zz = buffer_ddpd[509];

    auto g_yz_yz_y_xx = buffer_ddpd[510];

    auto g_yz_yz_y_xy = buffer_ddpd[511];

    auto g_yz_yz_y_xz = buffer_ddpd[512];

    auto g_yz_yz_y_yy = buffer_ddpd[513];

    auto g_yz_yz_y_yz = buffer_ddpd[514];

    auto g_yz_yz_y_zz = buffer_ddpd[515];

    auto g_yz_yz_z_xx = buffer_ddpd[516];

    auto g_yz_yz_z_xy = buffer_ddpd[517];

    auto g_yz_yz_z_xz = buffer_ddpd[518];

    auto g_yz_yz_z_yy = buffer_ddpd[519];

    auto g_yz_yz_z_yz = buffer_ddpd[520];

    auto g_yz_yz_z_zz = buffer_ddpd[521];

    auto g_yz_zz_x_xx = buffer_ddpd[522];

    auto g_yz_zz_x_xy = buffer_ddpd[523];

    auto g_yz_zz_x_xz = buffer_ddpd[524];

    auto g_yz_zz_x_yy = buffer_ddpd[525];

    auto g_yz_zz_x_yz = buffer_ddpd[526];

    auto g_yz_zz_x_zz = buffer_ddpd[527];

    auto g_yz_zz_y_xx = buffer_ddpd[528];

    auto g_yz_zz_y_xy = buffer_ddpd[529];

    auto g_yz_zz_y_xz = buffer_ddpd[530];

    auto g_yz_zz_y_yy = buffer_ddpd[531];

    auto g_yz_zz_y_yz = buffer_ddpd[532];

    auto g_yz_zz_y_zz = buffer_ddpd[533];

    auto g_yz_zz_z_xx = buffer_ddpd[534];

    auto g_yz_zz_z_xy = buffer_ddpd[535];

    auto g_yz_zz_z_xz = buffer_ddpd[536];

    auto g_yz_zz_z_yy = buffer_ddpd[537];

    auto g_yz_zz_z_yz = buffer_ddpd[538];

    auto g_yz_zz_z_zz = buffer_ddpd[539];

    auto g_zz_xx_x_xx = buffer_ddpd[540];

    auto g_zz_xx_x_xy = buffer_ddpd[541];

    auto g_zz_xx_x_xz = buffer_ddpd[542];

    auto g_zz_xx_x_yy = buffer_ddpd[543];

    auto g_zz_xx_x_yz = buffer_ddpd[544];

    auto g_zz_xx_x_zz = buffer_ddpd[545];

    auto g_zz_xx_y_xx = buffer_ddpd[546];

    auto g_zz_xx_y_xy = buffer_ddpd[547];

    auto g_zz_xx_y_xz = buffer_ddpd[548];

    auto g_zz_xx_y_yy = buffer_ddpd[549];

    auto g_zz_xx_y_yz = buffer_ddpd[550];

    auto g_zz_xx_y_zz = buffer_ddpd[551];

    auto g_zz_xx_z_xx = buffer_ddpd[552];

    auto g_zz_xx_z_xy = buffer_ddpd[553];

    auto g_zz_xx_z_xz = buffer_ddpd[554];

    auto g_zz_xx_z_yy = buffer_ddpd[555];

    auto g_zz_xx_z_yz = buffer_ddpd[556];

    auto g_zz_xx_z_zz = buffer_ddpd[557];

    auto g_zz_xy_x_xx = buffer_ddpd[558];

    auto g_zz_xy_x_xy = buffer_ddpd[559];

    auto g_zz_xy_x_xz = buffer_ddpd[560];

    auto g_zz_xy_x_yy = buffer_ddpd[561];

    auto g_zz_xy_x_yz = buffer_ddpd[562];

    auto g_zz_xy_x_zz = buffer_ddpd[563];

    auto g_zz_xy_y_xx = buffer_ddpd[564];

    auto g_zz_xy_y_xy = buffer_ddpd[565];

    auto g_zz_xy_y_xz = buffer_ddpd[566];

    auto g_zz_xy_y_yy = buffer_ddpd[567];

    auto g_zz_xy_y_yz = buffer_ddpd[568];

    auto g_zz_xy_y_zz = buffer_ddpd[569];

    auto g_zz_xy_z_xx = buffer_ddpd[570];

    auto g_zz_xy_z_xy = buffer_ddpd[571];

    auto g_zz_xy_z_xz = buffer_ddpd[572];

    auto g_zz_xy_z_yy = buffer_ddpd[573];

    auto g_zz_xy_z_yz = buffer_ddpd[574];

    auto g_zz_xy_z_zz = buffer_ddpd[575];

    auto g_zz_xz_x_xx = buffer_ddpd[576];

    auto g_zz_xz_x_xy = buffer_ddpd[577];

    auto g_zz_xz_x_xz = buffer_ddpd[578];

    auto g_zz_xz_x_yy = buffer_ddpd[579];

    auto g_zz_xz_x_yz = buffer_ddpd[580];

    auto g_zz_xz_x_zz = buffer_ddpd[581];

    auto g_zz_xz_y_xx = buffer_ddpd[582];

    auto g_zz_xz_y_xy = buffer_ddpd[583];

    auto g_zz_xz_y_xz = buffer_ddpd[584];

    auto g_zz_xz_y_yy = buffer_ddpd[585];

    auto g_zz_xz_y_yz = buffer_ddpd[586];

    auto g_zz_xz_y_zz = buffer_ddpd[587];

    auto g_zz_xz_z_xx = buffer_ddpd[588];

    auto g_zz_xz_z_xy = buffer_ddpd[589];

    auto g_zz_xz_z_xz = buffer_ddpd[590];

    auto g_zz_xz_z_yy = buffer_ddpd[591];

    auto g_zz_xz_z_yz = buffer_ddpd[592];

    auto g_zz_xz_z_zz = buffer_ddpd[593];

    auto g_zz_yy_x_xx = buffer_ddpd[594];

    auto g_zz_yy_x_xy = buffer_ddpd[595];

    auto g_zz_yy_x_xz = buffer_ddpd[596];

    auto g_zz_yy_x_yy = buffer_ddpd[597];

    auto g_zz_yy_x_yz = buffer_ddpd[598];

    auto g_zz_yy_x_zz = buffer_ddpd[599];

    auto g_zz_yy_y_xx = buffer_ddpd[600];

    auto g_zz_yy_y_xy = buffer_ddpd[601];

    auto g_zz_yy_y_xz = buffer_ddpd[602];

    auto g_zz_yy_y_yy = buffer_ddpd[603];

    auto g_zz_yy_y_yz = buffer_ddpd[604];

    auto g_zz_yy_y_zz = buffer_ddpd[605];

    auto g_zz_yy_z_xx = buffer_ddpd[606];

    auto g_zz_yy_z_xy = buffer_ddpd[607];

    auto g_zz_yy_z_xz = buffer_ddpd[608];

    auto g_zz_yy_z_yy = buffer_ddpd[609];

    auto g_zz_yy_z_yz = buffer_ddpd[610];

    auto g_zz_yy_z_zz = buffer_ddpd[611];

    auto g_zz_yz_x_xx = buffer_ddpd[612];

    auto g_zz_yz_x_xy = buffer_ddpd[613];

    auto g_zz_yz_x_xz = buffer_ddpd[614];

    auto g_zz_yz_x_yy = buffer_ddpd[615];

    auto g_zz_yz_x_yz = buffer_ddpd[616];

    auto g_zz_yz_x_zz = buffer_ddpd[617];

    auto g_zz_yz_y_xx = buffer_ddpd[618];

    auto g_zz_yz_y_xy = buffer_ddpd[619];

    auto g_zz_yz_y_xz = buffer_ddpd[620];

    auto g_zz_yz_y_yy = buffer_ddpd[621];

    auto g_zz_yz_y_yz = buffer_ddpd[622];

    auto g_zz_yz_y_zz = buffer_ddpd[623];

    auto g_zz_yz_z_xx = buffer_ddpd[624];

    auto g_zz_yz_z_xy = buffer_ddpd[625];

    auto g_zz_yz_z_xz = buffer_ddpd[626];

    auto g_zz_yz_z_yy = buffer_ddpd[627];

    auto g_zz_yz_z_yz = buffer_ddpd[628];

    auto g_zz_yz_z_zz = buffer_ddpd[629];

    auto g_zz_zz_x_xx = buffer_ddpd[630];

    auto g_zz_zz_x_xy = buffer_ddpd[631];

    auto g_zz_zz_x_xz = buffer_ddpd[632];

    auto g_zz_zz_x_yy = buffer_ddpd[633];

    auto g_zz_zz_x_yz = buffer_ddpd[634];

    auto g_zz_zz_x_zz = buffer_ddpd[635];

    auto g_zz_zz_y_xx = buffer_ddpd[636];

    auto g_zz_zz_y_xy = buffer_ddpd[637];

    auto g_zz_zz_y_xz = buffer_ddpd[638];

    auto g_zz_zz_y_yy = buffer_ddpd[639];

    auto g_zz_zz_y_yz = buffer_ddpd[640];

    auto g_zz_zz_y_zz = buffer_ddpd[641];

    auto g_zz_zz_z_xx = buffer_ddpd[642];

    auto g_zz_zz_z_xy = buffer_ddpd[643];

    auto g_zz_zz_z_xz = buffer_ddpd[644];

    auto g_zz_zz_z_yy = buffer_ddpd[645];

    auto g_zz_zz_z_yz = buffer_ddpd[646];

    auto g_zz_zz_z_zz = buffer_ddpd[647];

    /// Set up components of integrals buffer : buffer_1100_pppd

    auto g_x_x_0_0_x_x_x_xx = buffer_1100_pppd[0];

    auto g_x_x_0_0_x_x_x_xy = buffer_1100_pppd[1];

    auto g_x_x_0_0_x_x_x_xz = buffer_1100_pppd[2];

    auto g_x_x_0_0_x_x_x_yy = buffer_1100_pppd[3];

    auto g_x_x_0_0_x_x_x_yz = buffer_1100_pppd[4];

    auto g_x_x_0_0_x_x_x_zz = buffer_1100_pppd[5];

    auto g_x_x_0_0_x_x_y_xx = buffer_1100_pppd[6];

    auto g_x_x_0_0_x_x_y_xy = buffer_1100_pppd[7];

    auto g_x_x_0_0_x_x_y_xz = buffer_1100_pppd[8];

    auto g_x_x_0_0_x_x_y_yy = buffer_1100_pppd[9];

    auto g_x_x_0_0_x_x_y_yz = buffer_1100_pppd[10];

    auto g_x_x_0_0_x_x_y_zz = buffer_1100_pppd[11];

    auto g_x_x_0_0_x_x_z_xx = buffer_1100_pppd[12];

    auto g_x_x_0_0_x_x_z_xy = buffer_1100_pppd[13];

    auto g_x_x_0_0_x_x_z_xz = buffer_1100_pppd[14];

    auto g_x_x_0_0_x_x_z_yy = buffer_1100_pppd[15];

    auto g_x_x_0_0_x_x_z_yz = buffer_1100_pppd[16];

    auto g_x_x_0_0_x_x_z_zz = buffer_1100_pppd[17];

    auto g_x_x_0_0_x_y_x_xx = buffer_1100_pppd[18];

    auto g_x_x_0_0_x_y_x_xy = buffer_1100_pppd[19];

    auto g_x_x_0_0_x_y_x_xz = buffer_1100_pppd[20];

    auto g_x_x_0_0_x_y_x_yy = buffer_1100_pppd[21];

    auto g_x_x_0_0_x_y_x_yz = buffer_1100_pppd[22];

    auto g_x_x_0_0_x_y_x_zz = buffer_1100_pppd[23];

    auto g_x_x_0_0_x_y_y_xx = buffer_1100_pppd[24];

    auto g_x_x_0_0_x_y_y_xy = buffer_1100_pppd[25];

    auto g_x_x_0_0_x_y_y_xz = buffer_1100_pppd[26];

    auto g_x_x_0_0_x_y_y_yy = buffer_1100_pppd[27];

    auto g_x_x_0_0_x_y_y_yz = buffer_1100_pppd[28];

    auto g_x_x_0_0_x_y_y_zz = buffer_1100_pppd[29];

    auto g_x_x_0_0_x_y_z_xx = buffer_1100_pppd[30];

    auto g_x_x_0_0_x_y_z_xy = buffer_1100_pppd[31];

    auto g_x_x_0_0_x_y_z_xz = buffer_1100_pppd[32];

    auto g_x_x_0_0_x_y_z_yy = buffer_1100_pppd[33];

    auto g_x_x_0_0_x_y_z_yz = buffer_1100_pppd[34];

    auto g_x_x_0_0_x_y_z_zz = buffer_1100_pppd[35];

    auto g_x_x_0_0_x_z_x_xx = buffer_1100_pppd[36];

    auto g_x_x_0_0_x_z_x_xy = buffer_1100_pppd[37];

    auto g_x_x_0_0_x_z_x_xz = buffer_1100_pppd[38];

    auto g_x_x_0_0_x_z_x_yy = buffer_1100_pppd[39];

    auto g_x_x_0_0_x_z_x_yz = buffer_1100_pppd[40];

    auto g_x_x_0_0_x_z_x_zz = buffer_1100_pppd[41];

    auto g_x_x_0_0_x_z_y_xx = buffer_1100_pppd[42];

    auto g_x_x_0_0_x_z_y_xy = buffer_1100_pppd[43];

    auto g_x_x_0_0_x_z_y_xz = buffer_1100_pppd[44];

    auto g_x_x_0_0_x_z_y_yy = buffer_1100_pppd[45];

    auto g_x_x_0_0_x_z_y_yz = buffer_1100_pppd[46];

    auto g_x_x_0_0_x_z_y_zz = buffer_1100_pppd[47];

    auto g_x_x_0_0_x_z_z_xx = buffer_1100_pppd[48];

    auto g_x_x_0_0_x_z_z_xy = buffer_1100_pppd[49];

    auto g_x_x_0_0_x_z_z_xz = buffer_1100_pppd[50];

    auto g_x_x_0_0_x_z_z_yy = buffer_1100_pppd[51];

    auto g_x_x_0_0_x_z_z_yz = buffer_1100_pppd[52];

    auto g_x_x_0_0_x_z_z_zz = buffer_1100_pppd[53];

    auto g_x_x_0_0_y_x_x_xx = buffer_1100_pppd[54];

    auto g_x_x_0_0_y_x_x_xy = buffer_1100_pppd[55];

    auto g_x_x_0_0_y_x_x_xz = buffer_1100_pppd[56];

    auto g_x_x_0_0_y_x_x_yy = buffer_1100_pppd[57];

    auto g_x_x_0_0_y_x_x_yz = buffer_1100_pppd[58];

    auto g_x_x_0_0_y_x_x_zz = buffer_1100_pppd[59];

    auto g_x_x_0_0_y_x_y_xx = buffer_1100_pppd[60];

    auto g_x_x_0_0_y_x_y_xy = buffer_1100_pppd[61];

    auto g_x_x_0_0_y_x_y_xz = buffer_1100_pppd[62];

    auto g_x_x_0_0_y_x_y_yy = buffer_1100_pppd[63];

    auto g_x_x_0_0_y_x_y_yz = buffer_1100_pppd[64];

    auto g_x_x_0_0_y_x_y_zz = buffer_1100_pppd[65];

    auto g_x_x_0_0_y_x_z_xx = buffer_1100_pppd[66];

    auto g_x_x_0_0_y_x_z_xy = buffer_1100_pppd[67];

    auto g_x_x_0_0_y_x_z_xz = buffer_1100_pppd[68];

    auto g_x_x_0_0_y_x_z_yy = buffer_1100_pppd[69];

    auto g_x_x_0_0_y_x_z_yz = buffer_1100_pppd[70];

    auto g_x_x_0_0_y_x_z_zz = buffer_1100_pppd[71];

    auto g_x_x_0_0_y_y_x_xx = buffer_1100_pppd[72];

    auto g_x_x_0_0_y_y_x_xy = buffer_1100_pppd[73];

    auto g_x_x_0_0_y_y_x_xz = buffer_1100_pppd[74];

    auto g_x_x_0_0_y_y_x_yy = buffer_1100_pppd[75];

    auto g_x_x_0_0_y_y_x_yz = buffer_1100_pppd[76];

    auto g_x_x_0_0_y_y_x_zz = buffer_1100_pppd[77];

    auto g_x_x_0_0_y_y_y_xx = buffer_1100_pppd[78];

    auto g_x_x_0_0_y_y_y_xy = buffer_1100_pppd[79];

    auto g_x_x_0_0_y_y_y_xz = buffer_1100_pppd[80];

    auto g_x_x_0_0_y_y_y_yy = buffer_1100_pppd[81];

    auto g_x_x_0_0_y_y_y_yz = buffer_1100_pppd[82];

    auto g_x_x_0_0_y_y_y_zz = buffer_1100_pppd[83];

    auto g_x_x_0_0_y_y_z_xx = buffer_1100_pppd[84];

    auto g_x_x_0_0_y_y_z_xy = buffer_1100_pppd[85];

    auto g_x_x_0_0_y_y_z_xz = buffer_1100_pppd[86];

    auto g_x_x_0_0_y_y_z_yy = buffer_1100_pppd[87];

    auto g_x_x_0_0_y_y_z_yz = buffer_1100_pppd[88];

    auto g_x_x_0_0_y_y_z_zz = buffer_1100_pppd[89];

    auto g_x_x_0_0_y_z_x_xx = buffer_1100_pppd[90];

    auto g_x_x_0_0_y_z_x_xy = buffer_1100_pppd[91];

    auto g_x_x_0_0_y_z_x_xz = buffer_1100_pppd[92];

    auto g_x_x_0_0_y_z_x_yy = buffer_1100_pppd[93];

    auto g_x_x_0_0_y_z_x_yz = buffer_1100_pppd[94];

    auto g_x_x_0_0_y_z_x_zz = buffer_1100_pppd[95];

    auto g_x_x_0_0_y_z_y_xx = buffer_1100_pppd[96];

    auto g_x_x_0_0_y_z_y_xy = buffer_1100_pppd[97];

    auto g_x_x_0_0_y_z_y_xz = buffer_1100_pppd[98];

    auto g_x_x_0_0_y_z_y_yy = buffer_1100_pppd[99];

    auto g_x_x_0_0_y_z_y_yz = buffer_1100_pppd[100];

    auto g_x_x_0_0_y_z_y_zz = buffer_1100_pppd[101];

    auto g_x_x_0_0_y_z_z_xx = buffer_1100_pppd[102];

    auto g_x_x_0_0_y_z_z_xy = buffer_1100_pppd[103];

    auto g_x_x_0_0_y_z_z_xz = buffer_1100_pppd[104];

    auto g_x_x_0_0_y_z_z_yy = buffer_1100_pppd[105];

    auto g_x_x_0_0_y_z_z_yz = buffer_1100_pppd[106];

    auto g_x_x_0_0_y_z_z_zz = buffer_1100_pppd[107];

    auto g_x_x_0_0_z_x_x_xx = buffer_1100_pppd[108];

    auto g_x_x_0_0_z_x_x_xy = buffer_1100_pppd[109];

    auto g_x_x_0_0_z_x_x_xz = buffer_1100_pppd[110];

    auto g_x_x_0_0_z_x_x_yy = buffer_1100_pppd[111];

    auto g_x_x_0_0_z_x_x_yz = buffer_1100_pppd[112];

    auto g_x_x_0_0_z_x_x_zz = buffer_1100_pppd[113];

    auto g_x_x_0_0_z_x_y_xx = buffer_1100_pppd[114];

    auto g_x_x_0_0_z_x_y_xy = buffer_1100_pppd[115];

    auto g_x_x_0_0_z_x_y_xz = buffer_1100_pppd[116];

    auto g_x_x_0_0_z_x_y_yy = buffer_1100_pppd[117];

    auto g_x_x_0_0_z_x_y_yz = buffer_1100_pppd[118];

    auto g_x_x_0_0_z_x_y_zz = buffer_1100_pppd[119];

    auto g_x_x_0_0_z_x_z_xx = buffer_1100_pppd[120];

    auto g_x_x_0_0_z_x_z_xy = buffer_1100_pppd[121];

    auto g_x_x_0_0_z_x_z_xz = buffer_1100_pppd[122];

    auto g_x_x_0_0_z_x_z_yy = buffer_1100_pppd[123];

    auto g_x_x_0_0_z_x_z_yz = buffer_1100_pppd[124];

    auto g_x_x_0_0_z_x_z_zz = buffer_1100_pppd[125];

    auto g_x_x_0_0_z_y_x_xx = buffer_1100_pppd[126];

    auto g_x_x_0_0_z_y_x_xy = buffer_1100_pppd[127];

    auto g_x_x_0_0_z_y_x_xz = buffer_1100_pppd[128];

    auto g_x_x_0_0_z_y_x_yy = buffer_1100_pppd[129];

    auto g_x_x_0_0_z_y_x_yz = buffer_1100_pppd[130];

    auto g_x_x_0_0_z_y_x_zz = buffer_1100_pppd[131];

    auto g_x_x_0_0_z_y_y_xx = buffer_1100_pppd[132];

    auto g_x_x_0_0_z_y_y_xy = buffer_1100_pppd[133];

    auto g_x_x_0_0_z_y_y_xz = buffer_1100_pppd[134];

    auto g_x_x_0_0_z_y_y_yy = buffer_1100_pppd[135];

    auto g_x_x_0_0_z_y_y_yz = buffer_1100_pppd[136];

    auto g_x_x_0_0_z_y_y_zz = buffer_1100_pppd[137];

    auto g_x_x_0_0_z_y_z_xx = buffer_1100_pppd[138];

    auto g_x_x_0_0_z_y_z_xy = buffer_1100_pppd[139];

    auto g_x_x_0_0_z_y_z_xz = buffer_1100_pppd[140];

    auto g_x_x_0_0_z_y_z_yy = buffer_1100_pppd[141];

    auto g_x_x_0_0_z_y_z_yz = buffer_1100_pppd[142];

    auto g_x_x_0_0_z_y_z_zz = buffer_1100_pppd[143];

    auto g_x_x_0_0_z_z_x_xx = buffer_1100_pppd[144];

    auto g_x_x_0_0_z_z_x_xy = buffer_1100_pppd[145];

    auto g_x_x_0_0_z_z_x_xz = buffer_1100_pppd[146];

    auto g_x_x_0_0_z_z_x_yy = buffer_1100_pppd[147];

    auto g_x_x_0_0_z_z_x_yz = buffer_1100_pppd[148];

    auto g_x_x_0_0_z_z_x_zz = buffer_1100_pppd[149];

    auto g_x_x_0_0_z_z_y_xx = buffer_1100_pppd[150];

    auto g_x_x_0_0_z_z_y_xy = buffer_1100_pppd[151];

    auto g_x_x_0_0_z_z_y_xz = buffer_1100_pppd[152];

    auto g_x_x_0_0_z_z_y_yy = buffer_1100_pppd[153];

    auto g_x_x_0_0_z_z_y_yz = buffer_1100_pppd[154];

    auto g_x_x_0_0_z_z_y_zz = buffer_1100_pppd[155];

    auto g_x_x_0_0_z_z_z_xx = buffer_1100_pppd[156];

    auto g_x_x_0_0_z_z_z_xy = buffer_1100_pppd[157];

    auto g_x_x_0_0_z_z_z_xz = buffer_1100_pppd[158];

    auto g_x_x_0_0_z_z_z_yy = buffer_1100_pppd[159];

    auto g_x_x_0_0_z_z_z_yz = buffer_1100_pppd[160];

    auto g_x_x_0_0_z_z_z_zz = buffer_1100_pppd[161];

    auto g_x_y_0_0_x_x_x_xx = buffer_1100_pppd[162];

    auto g_x_y_0_0_x_x_x_xy = buffer_1100_pppd[163];

    auto g_x_y_0_0_x_x_x_xz = buffer_1100_pppd[164];

    auto g_x_y_0_0_x_x_x_yy = buffer_1100_pppd[165];

    auto g_x_y_0_0_x_x_x_yz = buffer_1100_pppd[166];

    auto g_x_y_0_0_x_x_x_zz = buffer_1100_pppd[167];

    auto g_x_y_0_0_x_x_y_xx = buffer_1100_pppd[168];

    auto g_x_y_0_0_x_x_y_xy = buffer_1100_pppd[169];

    auto g_x_y_0_0_x_x_y_xz = buffer_1100_pppd[170];

    auto g_x_y_0_0_x_x_y_yy = buffer_1100_pppd[171];

    auto g_x_y_0_0_x_x_y_yz = buffer_1100_pppd[172];

    auto g_x_y_0_0_x_x_y_zz = buffer_1100_pppd[173];

    auto g_x_y_0_0_x_x_z_xx = buffer_1100_pppd[174];

    auto g_x_y_0_0_x_x_z_xy = buffer_1100_pppd[175];

    auto g_x_y_0_0_x_x_z_xz = buffer_1100_pppd[176];

    auto g_x_y_0_0_x_x_z_yy = buffer_1100_pppd[177];

    auto g_x_y_0_0_x_x_z_yz = buffer_1100_pppd[178];

    auto g_x_y_0_0_x_x_z_zz = buffer_1100_pppd[179];

    auto g_x_y_0_0_x_y_x_xx = buffer_1100_pppd[180];

    auto g_x_y_0_0_x_y_x_xy = buffer_1100_pppd[181];

    auto g_x_y_0_0_x_y_x_xz = buffer_1100_pppd[182];

    auto g_x_y_0_0_x_y_x_yy = buffer_1100_pppd[183];

    auto g_x_y_0_0_x_y_x_yz = buffer_1100_pppd[184];

    auto g_x_y_0_0_x_y_x_zz = buffer_1100_pppd[185];

    auto g_x_y_0_0_x_y_y_xx = buffer_1100_pppd[186];

    auto g_x_y_0_0_x_y_y_xy = buffer_1100_pppd[187];

    auto g_x_y_0_0_x_y_y_xz = buffer_1100_pppd[188];

    auto g_x_y_0_0_x_y_y_yy = buffer_1100_pppd[189];

    auto g_x_y_0_0_x_y_y_yz = buffer_1100_pppd[190];

    auto g_x_y_0_0_x_y_y_zz = buffer_1100_pppd[191];

    auto g_x_y_0_0_x_y_z_xx = buffer_1100_pppd[192];

    auto g_x_y_0_0_x_y_z_xy = buffer_1100_pppd[193];

    auto g_x_y_0_0_x_y_z_xz = buffer_1100_pppd[194];

    auto g_x_y_0_0_x_y_z_yy = buffer_1100_pppd[195];

    auto g_x_y_0_0_x_y_z_yz = buffer_1100_pppd[196];

    auto g_x_y_0_0_x_y_z_zz = buffer_1100_pppd[197];

    auto g_x_y_0_0_x_z_x_xx = buffer_1100_pppd[198];

    auto g_x_y_0_0_x_z_x_xy = buffer_1100_pppd[199];

    auto g_x_y_0_0_x_z_x_xz = buffer_1100_pppd[200];

    auto g_x_y_0_0_x_z_x_yy = buffer_1100_pppd[201];

    auto g_x_y_0_0_x_z_x_yz = buffer_1100_pppd[202];

    auto g_x_y_0_0_x_z_x_zz = buffer_1100_pppd[203];

    auto g_x_y_0_0_x_z_y_xx = buffer_1100_pppd[204];

    auto g_x_y_0_0_x_z_y_xy = buffer_1100_pppd[205];

    auto g_x_y_0_0_x_z_y_xz = buffer_1100_pppd[206];

    auto g_x_y_0_0_x_z_y_yy = buffer_1100_pppd[207];

    auto g_x_y_0_0_x_z_y_yz = buffer_1100_pppd[208];

    auto g_x_y_0_0_x_z_y_zz = buffer_1100_pppd[209];

    auto g_x_y_0_0_x_z_z_xx = buffer_1100_pppd[210];

    auto g_x_y_0_0_x_z_z_xy = buffer_1100_pppd[211];

    auto g_x_y_0_0_x_z_z_xz = buffer_1100_pppd[212];

    auto g_x_y_0_0_x_z_z_yy = buffer_1100_pppd[213];

    auto g_x_y_0_0_x_z_z_yz = buffer_1100_pppd[214];

    auto g_x_y_0_0_x_z_z_zz = buffer_1100_pppd[215];

    auto g_x_y_0_0_y_x_x_xx = buffer_1100_pppd[216];

    auto g_x_y_0_0_y_x_x_xy = buffer_1100_pppd[217];

    auto g_x_y_0_0_y_x_x_xz = buffer_1100_pppd[218];

    auto g_x_y_0_0_y_x_x_yy = buffer_1100_pppd[219];

    auto g_x_y_0_0_y_x_x_yz = buffer_1100_pppd[220];

    auto g_x_y_0_0_y_x_x_zz = buffer_1100_pppd[221];

    auto g_x_y_0_0_y_x_y_xx = buffer_1100_pppd[222];

    auto g_x_y_0_0_y_x_y_xy = buffer_1100_pppd[223];

    auto g_x_y_0_0_y_x_y_xz = buffer_1100_pppd[224];

    auto g_x_y_0_0_y_x_y_yy = buffer_1100_pppd[225];

    auto g_x_y_0_0_y_x_y_yz = buffer_1100_pppd[226];

    auto g_x_y_0_0_y_x_y_zz = buffer_1100_pppd[227];

    auto g_x_y_0_0_y_x_z_xx = buffer_1100_pppd[228];

    auto g_x_y_0_0_y_x_z_xy = buffer_1100_pppd[229];

    auto g_x_y_0_0_y_x_z_xz = buffer_1100_pppd[230];

    auto g_x_y_0_0_y_x_z_yy = buffer_1100_pppd[231];

    auto g_x_y_0_0_y_x_z_yz = buffer_1100_pppd[232];

    auto g_x_y_0_0_y_x_z_zz = buffer_1100_pppd[233];

    auto g_x_y_0_0_y_y_x_xx = buffer_1100_pppd[234];

    auto g_x_y_0_0_y_y_x_xy = buffer_1100_pppd[235];

    auto g_x_y_0_0_y_y_x_xz = buffer_1100_pppd[236];

    auto g_x_y_0_0_y_y_x_yy = buffer_1100_pppd[237];

    auto g_x_y_0_0_y_y_x_yz = buffer_1100_pppd[238];

    auto g_x_y_0_0_y_y_x_zz = buffer_1100_pppd[239];

    auto g_x_y_0_0_y_y_y_xx = buffer_1100_pppd[240];

    auto g_x_y_0_0_y_y_y_xy = buffer_1100_pppd[241];

    auto g_x_y_0_0_y_y_y_xz = buffer_1100_pppd[242];

    auto g_x_y_0_0_y_y_y_yy = buffer_1100_pppd[243];

    auto g_x_y_0_0_y_y_y_yz = buffer_1100_pppd[244];

    auto g_x_y_0_0_y_y_y_zz = buffer_1100_pppd[245];

    auto g_x_y_0_0_y_y_z_xx = buffer_1100_pppd[246];

    auto g_x_y_0_0_y_y_z_xy = buffer_1100_pppd[247];

    auto g_x_y_0_0_y_y_z_xz = buffer_1100_pppd[248];

    auto g_x_y_0_0_y_y_z_yy = buffer_1100_pppd[249];

    auto g_x_y_0_0_y_y_z_yz = buffer_1100_pppd[250];

    auto g_x_y_0_0_y_y_z_zz = buffer_1100_pppd[251];

    auto g_x_y_0_0_y_z_x_xx = buffer_1100_pppd[252];

    auto g_x_y_0_0_y_z_x_xy = buffer_1100_pppd[253];

    auto g_x_y_0_0_y_z_x_xz = buffer_1100_pppd[254];

    auto g_x_y_0_0_y_z_x_yy = buffer_1100_pppd[255];

    auto g_x_y_0_0_y_z_x_yz = buffer_1100_pppd[256];

    auto g_x_y_0_0_y_z_x_zz = buffer_1100_pppd[257];

    auto g_x_y_0_0_y_z_y_xx = buffer_1100_pppd[258];

    auto g_x_y_0_0_y_z_y_xy = buffer_1100_pppd[259];

    auto g_x_y_0_0_y_z_y_xz = buffer_1100_pppd[260];

    auto g_x_y_0_0_y_z_y_yy = buffer_1100_pppd[261];

    auto g_x_y_0_0_y_z_y_yz = buffer_1100_pppd[262];

    auto g_x_y_0_0_y_z_y_zz = buffer_1100_pppd[263];

    auto g_x_y_0_0_y_z_z_xx = buffer_1100_pppd[264];

    auto g_x_y_0_0_y_z_z_xy = buffer_1100_pppd[265];

    auto g_x_y_0_0_y_z_z_xz = buffer_1100_pppd[266];

    auto g_x_y_0_0_y_z_z_yy = buffer_1100_pppd[267];

    auto g_x_y_0_0_y_z_z_yz = buffer_1100_pppd[268];

    auto g_x_y_0_0_y_z_z_zz = buffer_1100_pppd[269];

    auto g_x_y_0_0_z_x_x_xx = buffer_1100_pppd[270];

    auto g_x_y_0_0_z_x_x_xy = buffer_1100_pppd[271];

    auto g_x_y_0_0_z_x_x_xz = buffer_1100_pppd[272];

    auto g_x_y_0_0_z_x_x_yy = buffer_1100_pppd[273];

    auto g_x_y_0_0_z_x_x_yz = buffer_1100_pppd[274];

    auto g_x_y_0_0_z_x_x_zz = buffer_1100_pppd[275];

    auto g_x_y_0_0_z_x_y_xx = buffer_1100_pppd[276];

    auto g_x_y_0_0_z_x_y_xy = buffer_1100_pppd[277];

    auto g_x_y_0_0_z_x_y_xz = buffer_1100_pppd[278];

    auto g_x_y_0_0_z_x_y_yy = buffer_1100_pppd[279];

    auto g_x_y_0_0_z_x_y_yz = buffer_1100_pppd[280];

    auto g_x_y_0_0_z_x_y_zz = buffer_1100_pppd[281];

    auto g_x_y_0_0_z_x_z_xx = buffer_1100_pppd[282];

    auto g_x_y_0_0_z_x_z_xy = buffer_1100_pppd[283];

    auto g_x_y_0_0_z_x_z_xz = buffer_1100_pppd[284];

    auto g_x_y_0_0_z_x_z_yy = buffer_1100_pppd[285];

    auto g_x_y_0_0_z_x_z_yz = buffer_1100_pppd[286];

    auto g_x_y_0_0_z_x_z_zz = buffer_1100_pppd[287];

    auto g_x_y_0_0_z_y_x_xx = buffer_1100_pppd[288];

    auto g_x_y_0_0_z_y_x_xy = buffer_1100_pppd[289];

    auto g_x_y_0_0_z_y_x_xz = buffer_1100_pppd[290];

    auto g_x_y_0_0_z_y_x_yy = buffer_1100_pppd[291];

    auto g_x_y_0_0_z_y_x_yz = buffer_1100_pppd[292];

    auto g_x_y_0_0_z_y_x_zz = buffer_1100_pppd[293];

    auto g_x_y_0_0_z_y_y_xx = buffer_1100_pppd[294];

    auto g_x_y_0_0_z_y_y_xy = buffer_1100_pppd[295];

    auto g_x_y_0_0_z_y_y_xz = buffer_1100_pppd[296];

    auto g_x_y_0_0_z_y_y_yy = buffer_1100_pppd[297];

    auto g_x_y_0_0_z_y_y_yz = buffer_1100_pppd[298];

    auto g_x_y_0_0_z_y_y_zz = buffer_1100_pppd[299];

    auto g_x_y_0_0_z_y_z_xx = buffer_1100_pppd[300];

    auto g_x_y_0_0_z_y_z_xy = buffer_1100_pppd[301];

    auto g_x_y_0_0_z_y_z_xz = buffer_1100_pppd[302];

    auto g_x_y_0_0_z_y_z_yy = buffer_1100_pppd[303];

    auto g_x_y_0_0_z_y_z_yz = buffer_1100_pppd[304];

    auto g_x_y_0_0_z_y_z_zz = buffer_1100_pppd[305];

    auto g_x_y_0_0_z_z_x_xx = buffer_1100_pppd[306];

    auto g_x_y_0_0_z_z_x_xy = buffer_1100_pppd[307];

    auto g_x_y_0_0_z_z_x_xz = buffer_1100_pppd[308];

    auto g_x_y_0_0_z_z_x_yy = buffer_1100_pppd[309];

    auto g_x_y_0_0_z_z_x_yz = buffer_1100_pppd[310];

    auto g_x_y_0_0_z_z_x_zz = buffer_1100_pppd[311];

    auto g_x_y_0_0_z_z_y_xx = buffer_1100_pppd[312];

    auto g_x_y_0_0_z_z_y_xy = buffer_1100_pppd[313];

    auto g_x_y_0_0_z_z_y_xz = buffer_1100_pppd[314];

    auto g_x_y_0_0_z_z_y_yy = buffer_1100_pppd[315];

    auto g_x_y_0_0_z_z_y_yz = buffer_1100_pppd[316];

    auto g_x_y_0_0_z_z_y_zz = buffer_1100_pppd[317];

    auto g_x_y_0_0_z_z_z_xx = buffer_1100_pppd[318];

    auto g_x_y_0_0_z_z_z_xy = buffer_1100_pppd[319];

    auto g_x_y_0_0_z_z_z_xz = buffer_1100_pppd[320];

    auto g_x_y_0_0_z_z_z_yy = buffer_1100_pppd[321];

    auto g_x_y_0_0_z_z_z_yz = buffer_1100_pppd[322];

    auto g_x_y_0_0_z_z_z_zz = buffer_1100_pppd[323];

    auto g_x_z_0_0_x_x_x_xx = buffer_1100_pppd[324];

    auto g_x_z_0_0_x_x_x_xy = buffer_1100_pppd[325];

    auto g_x_z_0_0_x_x_x_xz = buffer_1100_pppd[326];

    auto g_x_z_0_0_x_x_x_yy = buffer_1100_pppd[327];

    auto g_x_z_0_0_x_x_x_yz = buffer_1100_pppd[328];

    auto g_x_z_0_0_x_x_x_zz = buffer_1100_pppd[329];

    auto g_x_z_0_0_x_x_y_xx = buffer_1100_pppd[330];

    auto g_x_z_0_0_x_x_y_xy = buffer_1100_pppd[331];

    auto g_x_z_0_0_x_x_y_xz = buffer_1100_pppd[332];

    auto g_x_z_0_0_x_x_y_yy = buffer_1100_pppd[333];

    auto g_x_z_0_0_x_x_y_yz = buffer_1100_pppd[334];

    auto g_x_z_0_0_x_x_y_zz = buffer_1100_pppd[335];

    auto g_x_z_0_0_x_x_z_xx = buffer_1100_pppd[336];

    auto g_x_z_0_0_x_x_z_xy = buffer_1100_pppd[337];

    auto g_x_z_0_0_x_x_z_xz = buffer_1100_pppd[338];

    auto g_x_z_0_0_x_x_z_yy = buffer_1100_pppd[339];

    auto g_x_z_0_0_x_x_z_yz = buffer_1100_pppd[340];

    auto g_x_z_0_0_x_x_z_zz = buffer_1100_pppd[341];

    auto g_x_z_0_0_x_y_x_xx = buffer_1100_pppd[342];

    auto g_x_z_0_0_x_y_x_xy = buffer_1100_pppd[343];

    auto g_x_z_0_0_x_y_x_xz = buffer_1100_pppd[344];

    auto g_x_z_0_0_x_y_x_yy = buffer_1100_pppd[345];

    auto g_x_z_0_0_x_y_x_yz = buffer_1100_pppd[346];

    auto g_x_z_0_0_x_y_x_zz = buffer_1100_pppd[347];

    auto g_x_z_0_0_x_y_y_xx = buffer_1100_pppd[348];

    auto g_x_z_0_0_x_y_y_xy = buffer_1100_pppd[349];

    auto g_x_z_0_0_x_y_y_xz = buffer_1100_pppd[350];

    auto g_x_z_0_0_x_y_y_yy = buffer_1100_pppd[351];

    auto g_x_z_0_0_x_y_y_yz = buffer_1100_pppd[352];

    auto g_x_z_0_0_x_y_y_zz = buffer_1100_pppd[353];

    auto g_x_z_0_0_x_y_z_xx = buffer_1100_pppd[354];

    auto g_x_z_0_0_x_y_z_xy = buffer_1100_pppd[355];

    auto g_x_z_0_0_x_y_z_xz = buffer_1100_pppd[356];

    auto g_x_z_0_0_x_y_z_yy = buffer_1100_pppd[357];

    auto g_x_z_0_0_x_y_z_yz = buffer_1100_pppd[358];

    auto g_x_z_0_0_x_y_z_zz = buffer_1100_pppd[359];

    auto g_x_z_0_0_x_z_x_xx = buffer_1100_pppd[360];

    auto g_x_z_0_0_x_z_x_xy = buffer_1100_pppd[361];

    auto g_x_z_0_0_x_z_x_xz = buffer_1100_pppd[362];

    auto g_x_z_0_0_x_z_x_yy = buffer_1100_pppd[363];

    auto g_x_z_0_0_x_z_x_yz = buffer_1100_pppd[364];

    auto g_x_z_0_0_x_z_x_zz = buffer_1100_pppd[365];

    auto g_x_z_0_0_x_z_y_xx = buffer_1100_pppd[366];

    auto g_x_z_0_0_x_z_y_xy = buffer_1100_pppd[367];

    auto g_x_z_0_0_x_z_y_xz = buffer_1100_pppd[368];

    auto g_x_z_0_0_x_z_y_yy = buffer_1100_pppd[369];

    auto g_x_z_0_0_x_z_y_yz = buffer_1100_pppd[370];

    auto g_x_z_0_0_x_z_y_zz = buffer_1100_pppd[371];

    auto g_x_z_0_0_x_z_z_xx = buffer_1100_pppd[372];

    auto g_x_z_0_0_x_z_z_xy = buffer_1100_pppd[373];

    auto g_x_z_0_0_x_z_z_xz = buffer_1100_pppd[374];

    auto g_x_z_0_0_x_z_z_yy = buffer_1100_pppd[375];

    auto g_x_z_0_0_x_z_z_yz = buffer_1100_pppd[376];

    auto g_x_z_0_0_x_z_z_zz = buffer_1100_pppd[377];

    auto g_x_z_0_0_y_x_x_xx = buffer_1100_pppd[378];

    auto g_x_z_0_0_y_x_x_xy = buffer_1100_pppd[379];

    auto g_x_z_0_0_y_x_x_xz = buffer_1100_pppd[380];

    auto g_x_z_0_0_y_x_x_yy = buffer_1100_pppd[381];

    auto g_x_z_0_0_y_x_x_yz = buffer_1100_pppd[382];

    auto g_x_z_0_0_y_x_x_zz = buffer_1100_pppd[383];

    auto g_x_z_0_0_y_x_y_xx = buffer_1100_pppd[384];

    auto g_x_z_0_0_y_x_y_xy = buffer_1100_pppd[385];

    auto g_x_z_0_0_y_x_y_xz = buffer_1100_pppd[386];

    auto g_x_z_0_0_y_x_y_yy = buffer_1100_pppd[387];

    auto g_x_z_0_0_y_x_y_yz = buffer_1100_pppd[388];

    auto g_x_z_0_0_y_x_y_zz = buffer_1100_pppd[389];

    auto g_x_z_0_0_y_x_z_xx = buffer_1100_pppd[390];

    auto g_x_z_0_0_y_x_z_xy = buffer_1100_pppd[391];

    auto g_x_z_0_0_y_x_z_xz = buffer_1100_pppd[392];

    auto g_x_z_0_0_y_x_z_yy = buffer_1100_pppd[393];

    auto g_x_z_0_0_y_x_z_yz = buffer_1100_pppd[394];

    auto g_x_z_0_0_y_x_z_zz = buffer_1100_pppd[395];

    auto g_x_z_0_0_y_y_x_xx = buffer_1100_pppd[396];

    auto g_x_z_0_0_y_y_x_xy = buffer_1100_pppd[397];

    auto g_x_z_0_0_y_y_x_xz = buffer_1100_pppd[398];

    auto g_x_z_0_0_y_y_x_yy = buffer_1100_pppd[399];

    auto g_x_z_0_0_y_y_x_yz = buffer_1100_pppd[400];

    auto g_x_z_0_0_y_y_x_zz = buffer_1100_pppd[401];

    auto g_x_z_0_0_y_y_y_xx = buffer_1100_pppd[402];

    auto g_x_z_0_0_y_y_y_xy = buffer_1100_pppd[403];

    auto g_x_z_0_0_y_y_y_xz = buffer_1100_pppd[404];

    auto g_x_z_0_0_y_y_y_yy = buffer_1100_pppd[405];

    auto g_x_z_0_0_y_y_y_yz = buffer_1100_pppd[406];

    auto g_x_z_0_0_y_y_y_zz = buffer_1100_pppd[407];

    auto g_x_z_0_0_y_y_z_xx = buffer_1100_pppd[408];

    auto g_x_z_0_0_y_y_z_xy = buffer_1100_pppd[409];

    auto g_x_z_0_0_y_y_z_xz = buffer_1100_pppd[410];

    auto g_x_z_0_0_y_y_z_yy = buffer_1100_pppd[411];

    auto g_x_z_0_0_y_y_z_yz = buffer_1100_pppd[412];

    auto g_x_z_0_0_y_y_z_zz = buffer_1100_pppd[413];

    auto g_x_z_0_0_y_z_x_xx = buffer_1100_pppd[414];

    auto g_x_z_0_0_y_z_x_xy = buffer_1100_pppd[415];

    auto g_x_z_0_0_y_z_x_xz = buffer_1100_pppd[416];

    auto g_x_z_0_0_y_z_x_yy = buffer_1100_pppd[417];

    auto g_x_z_0_0_y_z_x_yz = buffer_1100_pppd[418];

    auto g_x_z_0_0_y_z_x_zz = buffer_1100_pppd[419];

    auto g_x_z_0_0_y_z_y_xx = buffer_1100_pppd[420];

    auto g_x_z_0_0_y_z_y_xy = buffer_1100_pppd[421];

    auto g_x_z_0_0_y_z_y_xz = buffer_1100_pppd[422];

    auto g_x_z_0_0_y_z_y_yy = buffer_1100_pppd[423];

    auto g_x_z_0_0_y_z_y_yz = buffer_1100_pppd[424];

    auto g_x_z_0_0_y_z_y_zz = buffer_1100_pppd[425];

    auto g_x_z_0_0_y_z_z_xx = buffer_1100_pppd[426];

    auto g_x_z_0_0_y_z_z_xy = buffer_1100_pppd[427];

    auto g_x_z_0_0_y_z_z_xz = buffer_1100_pppd[428];

    auto g_x_z_0_0_y_z_z_yy = buffer_1100_pppd[429];

    auto g_x_z_0_0_y_z_z_yz = buffer_1100_pppd[430];

    auto g_x_z_0_0_y_z_z_zz = buffer_1100_pppd[431];

    auto g_x_z_0_0_z_x_x_xx = buffer_1100_pppd[432];

    auto g_x_z_0_0_z_x_x_xy = buffer_1100_pppd[433];

    auto g_x_z_0_0_z_x_x_xz = buffer_1100_pppd[434];

    auto g_x_z_0_0_z_x_x_yy = buffer_1100_pppd[435];

    auto g_x_z_0_0_z_x_x_yz = buffer_1100_pppd[436];

    auto g_x_z_0_0_z_x_x_zz = buffer_1100_pppd[437];

    auto g_x_z_0_0_z_x_y_xx = buffer_1100_pppd[438];

    auto g_x_z_0_0_z_x_y_xy = buffer_1100_pppd[439];

    auto g_x_z_0_0_z_x_y_xz = buffer_1100_pppd[440];

    auto g_x_z_0_0_z_x_y_yy = buffer_1100_pppd[441];

    auto g_x_z_0_0_z_x_y_yz = buffer_1100_pppd[442];

    auto g_x_z_0_0_z_x_y_zz = buffer_1100_pppd[443];

    auto g_x_z_0_0_z_x_z_xx = buffer_1100_pppd[444];

    auto g_x_z_0_0_z_x_z_xy = buffer_1100_pppd[445];

    auto g_x_z_0_0_z_x_z_xz = buffer_1100_pppd[446];

    auto g_x_z_0_0_z_x_z_yy = buffer_1100_pppd[447];

    auto g_x_z_0_0_z_x_z_yz = buffer_1100_pppd[448];

    auto g_x_z_0_0_z_x_z_zz = buffer_1100_pppd[449];

    auto g_x_z_0_0_z_y_x_xx = buffer_1100_pppd[450];

    auto g_x_z_0_0_z_y_x_xy = buffer_1100_pppd[451];

    auto g_x_z_0_0_z_y_x_xz = buffer_1100_pppd[452];

    auto g_x_z_0_0_z_y_x_yy = buffer_1100_pppd[453];

    auto g_x_z_0_0_z_y_x_yz = buffer_1100_pppd[454];

    auto g_x_z_0_0_z_y_x_zz = buffer_1100_pppd[455];

    auto g_x_z_0_0_z_y_y_xx = buffer_1100_pppd[456];

    auto g_x_z_0_0_z_y_y_xy = buffer_1100_pppd[457];

    auto g_x_z_0_0_z_y_y_xz = buffer_1100_pppd[458];

    auto g_x_z_0_0_z_y_y_yy = buffer_1100_pppd[459];

    auto g_x_z_0_0_z_y_y_yz = buffer_1100_pppd[460];

    auto g_x_z_0_0_z_y_y_zz = buffer_1100_pppd[461];

    auto g_x_z_0_0_z_y_z_xx = buffer_1100_pppd[462];

    auto g_x_z_0_0_z_y_z_xy = buffer_1100_pppd[463];

    auto g_x_z_0_0_z_y_z_xz = buffer_1100_pppd[464];

    auto g_x_z_0_0_z_y_z_yy = buffer_1100_pppd[465];

    auto g_x_z_0_0_z_y_z_yz = buffer_1100_pppd[466];

    auto g_x_z_0_0_z_y_z_zz = buffer_1100_pppd[467];

    auto g_x_z_0_0_z_z_x_xx = buffer_1100_pppd[468];

    auto g_x_z_0_0_z_z_x_xy = buffer_1100_pppd[469];

    auto g_x_z_0_0_z_z_x_xz = buffer_1100_pppd[470];

    auto g_x_z_0_0_z_z_x_yy = buffer_1100_pppd[471];

    auto g_x_z_0_0_z_z_x_yz = buffer_1100_pppd[472];

    auto g_x_z_0_0_z_z_x_zz = buffer_1100_pppd[473];

    auto g_x_z_0_0_z_z_y_xx = buffer_1100_pppd[474];

    auto g_x_z_0_0_z_z_y_xy = buffer_1100_pppd[475];

    auto g_x_z_0_0_z_z_y_xz = buffer_1100_pppd[476];

    auto g_x_z_0_0_z_z_y_yy = buffer_1100_pppd[477];

    auto g_x_z_0_0_z_z_y_yz = buffer_1100_pppd[478];

    auto g_x_z_0_0_z_z_y_zz = buffer_1100_pppd[479];

    auto g_x_z_0_0_z_z_z_xx = buffer_1100_pppd[480];

    auto g_x_z_0_0_z_z_z_xy = buffer_1100_pppd[481];

    auto g_x_z_0_0_z_z_z_xz = buffer_1100_pppd[482];

    auto g_x_z_0_0_z_z_z_yy = buffer_1100_pppd[483];

    auto g_x_z_0_0_z_z_z_yz = buffer_1100_pppd[484];

    auto g_x_z_0_0_z_z_z_zz = buffer_1100_pppd[485];

    auto g_y_x_0_0_x_x_x_xx = buffer_1100_pppd[486];

    auto g_y_x_0_0_x_x_x_xy = buffer_1100_pppd[487];

    auto g_y_x_0_0_x_x_x_xz = buffer_1100_pppd[488];

    auto g_y_x_0_0_x_x_x_yy = buffer_1100_pppd[489];

    auto g_y_x_0_0_x_x_x_yz = buffer_1100_pppd[490];

    auto g_y_x_0_0_x_x_x_zz = buffer_1100_pppd[491];

    auto g_y_x_0_0_x_x_y_xx = buffer_1100_pppd[492];

    auto g_y_x_0_0_x_x_y_xy = buffer_1100_pppd[493];

    auto g_y_x_0_0_x_x_y_xz = buffer_1100_pppd[494];

    auto g_y_x_0_0_x_x_y_yy = buffer_1100_pppd[495];

    auto g_y_x_0_0_x_x_y_yz = buffer_1100_pppd[496];

    auto g_y_x_0_0_x_x_y_zz = buffer_1100_pppd[497];

    auto g_y_x_0_0_x_x_z_xx = buffer_1100_pppd[498];

    auto g_y_x_0_0_x_x_z_xy = buffer_1100_pppd[499];

    auto g_y_x_0_0_x_x_z_xz = buffer_1100_pppd[500];

    auto g_y_x_0_0_x_x_z_yy = buffer_1100_pppd[501];

    auto g_y_x_0_0_x_x_z_yz = buffer_1100_pppd[502];

    auto g_y_x_0_0_x_x_z_zz = buffer_1100_pppd[503];

    auto g_y_x_0_0_x_y_x_xx = buffer_1100_pppd[504];

    auto g_y_x_0_0_x_y_x_xy = buffer_1100_pppd[505];

    auto g_y_x_0_0_x_y_x_xz = buffer_1100_pppd[506];

    auto g_y_x_0_0_x_y_x_yy = buffer_1100_pppd[507];

    auto g_y_x_0_0_x_y_x_yz = buffer_1100_pppd[508];

    auto g_y_x_0_0_x_y_x_zz = buffer_1100_pppd[509];

    auto g_y_x_0_0_x_y_y_xx = buffer_1100_pppd[510];

    auto g_y_x_0_0_x_y_y_xy = buffer_1100_pppd[511];

    auto g_y_x_0_0_x_y_y_xz = buffer_1100_pppd[512];

    auto g_y_x_0_0_x_y_y_yy = buffer_1100_pppd[513];

    auto g_y_x_0_0_x_y_y_yz = buffer_1100_pppd[514];

    auto g_y_x_0_0_x_y_y_zz = buffer_1100_pppd[515];

    auto g_y_x_0_0_x_y_z_xx = buffer_1100_pppd[516];

    auto g_y_x_0_0_x_y_z_xy = buffer_1100_pppd[517];

    auto g_y_x_0_0_x_y_z_xz = buffer_1100_pppd[518];

    auto g_y_x_0_0_x_y_z_yy = buffer_1100_pppd[519];

    auto g_y_x_0_0_x_y_z_yz = buffer_1100_pppd[520];

    auto g_y_x_0_0_x_y_z_zz = buffer_1100_pppd[521];

    auto g_y_x_0_0_x_z_x_xx = buffer_1100_pppd[522];

    auto g_y_x_0_0_x_z_x_xy = buffer_1100_pppd[523];

    auto g_y_x_0_0_x_z_x_xz = buffer_1100_pppd[524];

    auto g_y_x_0_0_x_z_x_yy = buffer_1100_pppd[525];

    auto g_y_x_0_0_x_z_x_yz = buffer_1100_pppd[526];

    auto g_y_x_0_0_x_z_x_zz = buffer_1100_pppd[527];

    auto g_y_x_0_0_x_z_y_xx = buffer_1100_pppd[528];

    auto g_y_x_0_0_x_z_y_xy = buffer_1100_pppd[529];

    auto g_y_x_0_0_x_z_y_xz = buffer_1100_pppd[530];

    auto g_y_x_0_0_x_z_y_yy = buffer_1100_pppd[531];

    auto g_y_x_0_0_x_z_y_yz = buffer_1100_pppd[532];

    auto g_y_x_0_0_x_z_y_zz = buffer_1100_pppd[533];

    auto g_y_x_0_0_x_z_z_xx = buffer_1100_pppd[534];

    auto g_y_x_0_0_x_z_z_xy = buffer_1100_pppd[535];

    auto g_y_x_0_0_x_z_z_xz = buffer_1100_pppd[536];

    auto g_y_x_0_0_x_z_z_yy = buffer_1100_pppd[537];

    auto g_y_x_0_0_x_z_z_yz = buffer_1100_pppd[538];

    auto g_y_x_0_0_x_z_z_zz = buffer_1100_pppd[539];

    auto g_y_x_0_0_y_x_x_xx = buffer_1100_pppd[540];

    auto g_y_x_0_0_y_x_x_xy = buffer_1100_pppd[541];

    auto g_y_x_0_0_y_x_x_xz = buffer_1100_pppd[542];

    auto g_y_x_0_0_y_x_x_yy = buffer_1100_pppd[543];

    auto g_y_x_0_0_y_x_x_yz = buffer_1100_pppd[544];

    auto g_y_x_0_0_y_x_x_zz = buffer_1100_pppd[545];

    auto g_y_x_0_0_y_x_y_xx = buffer_1100_pppd[546];

    auto g_y_x_0_0_y_x_y_xy = buffer_1100_pppd[547];

    auto g_y_x_0_0_y_x_y_xz = buffer_1100_pppd[548];

    auto g_y_x_0_0_y_x_y_yy = buffer_1100_pppd[549];

    auto g_y_x_0_0_y_x_y_yz = buffer_1100_pppd[550];

    auto g_y_x_0_0_y_x_y_zz = buffer_1100_pppd[551];

    auto g_y_x_0_0_y_x_z_xx = buffer_1100_pppd[552];

    auto g_y_x_0_0_y_x_z_xy = buffer_1100_pppd[553];

    auto g_y_x_0_0_y_x_z_xz = buffer_1100_pppd[554];

    auto g_y_x_0_0_y_x_z_yy = buffer_1100_pppd[555];

    auto g_y_x_0_0_y_x_z_yz = buffer_1100_pppd[556];

    auto g_y_x_0_0_y_x_z_zz = buffer_1100_pppd[557];

    auto g_y_x_0_0_y_y_x_xx = buffer_1100_pppd[558];

    auto g_y_x_0_0_y_y_x_xy = buffer_1100_pppd[559];

    auto g_y_x_0_0_y_y_x_xz = buffer_1100_pppd[560];

    auto g_y_x_0_0_y_y_x_yy = buffer_1100_pppd[561];

    auto g_y_x_0_0_y_y_x_yz = buffer_1100_pppd[562];

    auto g_y_x_0_0_y_y_x_zz = buffer_1100_pppd[563];

    auto g_y_x_0_0_y_y_y_xx = buffer_1100_pppd[564];

    auto g_y_x_0_0_y_y_y_xy = buffer_1100_pppd[565];

    auto g_y_x_0_0_y_y_y_xz = buffer_1100_pppd[566];

    auto g_y_x_0_0_y_y_y_yy = buffer_1100_pppd[567];

    auto g_y_x_0_0_y_y_y_yz = buffer_1100_pppd[568];

    auto g_y_x_0_0_y_y_y_zz = buffer_1100_pppd[569];

    auto g_y_x_0_0_y_y_z_xx = buffer_1100_pppd[570];

    auto g_y_x_0_0_y_y_z_xy = buffer_1100_pppd[571];

    auto g_y_x_0_0_y_y_z_xz = buffer_1100_pppd[572];

    auto g_y_x_0_0_y_y_z_yy = buffer_1100_pppd[573];

    auto g_y_x_0_0_y_y_z_yz = buffer_1100_pppd[574];

    auto g_y_x_0_0_y_y_z_zz = buffer_1100_pppd[575];

    auto g_y_x_0_0_y_z_x_xx = buffer_1100_pppd[576];

    auto g_y_x_0_0_y_z_x_xy = buffer_1100_pppd[577];

    auto g_y_x_0_0_y_z_x_xz = buffer_1100_pppd[578];

    auto g_y_x_0_0_y_z_x_yy = buffer_1100_pppd[579];

    auto g_y_x_0_0_y_z_x_yz = buffer_1100_pppd[580];

    auto g_y_x_0_0_y_z_x_zz = buffer_1100_pppd[581];

    auto g_y_x_0_0_y_z_y_xx = buffer_1100_pppd[582];

    auto g_y_x_0_0_y_z_y_xy = buffer_1100_pppd[583];

    auto g_y_x_0_0_y_z_y_xz = buffer_1100_pppd[584];

    auto g_y_x_0_0_y_z_y_yy = buffer_1100_pppd[585];

    auto g_y_x_0_0_y_z_y_yz = buffer_1100_pppd[586];

    auto g_y_x_0_0_y_z_y_zz = buffer_1100_pppd[587];

    auto g_y_x_0_0_y_z_z_xx = buffer_1100_pppd[588];

    auto g_y_x_0_0_y_z_z_xy = buffer_1100_pppd[589];

    auto g_y_x_0_0_y_z_z_xz = buffer_1100_pppd[590];

    auto g_y_x_0_0_y_z_z_yy = buffer_1100_pppd[591];

    auto g_y_x_0_0_y_z_z_yz = buffer_1100_pppd[592];

    auto g_y_x_0_0_y_z_z_zz = buffer_1100_pppd[593];

    auto g_y_x_0_0_z_x_x_xx = buffer_1100_pppd[594];

    auto g_y_x_0_0_z_x_x_xy = buffer_1100_pppd[595];

    auto g_y_x_0_0_z_x_x_xz = buffer_1100_pppd[596];

    auto g_y_x_0_0_z_x_x_yy = buffer_1100_pppd[597];

    auto g_y_x_0_0_z_x_x_yz = buffer_1100_pppd[598];

    auto g_y_x_0_0_z_x_x_zz = buffer_1100_pppd[599];

    auto g_y_x_0_0_z_x_y_xx = buffer_1100_pppd[600];

    auto g_y_x_0_0_z_x_y_xy = buffer_1100_pppd[601];

    auto g_y_x_0_0_z_x_y_xz = buffer_1100_pppd[602];

    auto g_y_x_0_0_z_x_y_yy = buffer_1100_pppd[603];

    auto g_y_x_0_0_z_x_y_yz = buffer_1100_pppd[604];

    auto g_y_x_0_0_z_x_y_zz = buffer_1100_pppd[605];

    auto g_y_x_0_0_z_x_z_xx = buffer_1100_pppd[606];

    auto g_y_x_0_0_z_x_z_xy = buffer_1100_pppd[607];

    auto g_y_x_0_0_z_x_z_xz = buffer_1100_pppd[608];

    auto g_y_x_0_0_z_x_z_yy = buffer_1100_pppd[609];

    auto g_y_x_0_0_z_x_z_yz = buffer_1100_pppd[610];

    auto g_y_x_0_0_z_x_z_zz = buffer_1100_pppd[611];

    auto g_y_x_0_0_z_y_x_xx = buffer_1100_pppd[612];

    auto g_y_x_0_0_z_y_x_xy = buffer_1100_pppd[613];

    auto g_y_x_0_0_z_y_x_xz = buffer_1100_pppd[614];

    auto g_y_x_0_0_z_y_x_yy = buffer_1100_pppd[615];

    auto g_y_x_0_0_z_y_x_yz = buffer_1100_pppd[616];

    auto g_y_x_0_0_z_y_x_zz = buffer_1100_pppd[617];

    auto g_y_x_0_0_z_y_y_xx = buffer_1100_pppd[618];

    auto g_y_x_0_0_z_y_y_xy = buffer_1100_pppd[619];

    auto g_y_x_0_0_z_y_y_xz = buffer_1100_pppd[620];

    auto g_y_x_0_0_z_y_y_yy = buffer_1100_pppd[621];

    auto g_y_x_0_0_z_y_y_yz = buffer_1100_pppd[622];

    auto g_y_x_0_0_z_y_y_zz = buffer_1100_pppd[623];

    auto g_y_x_0_0_z_y_z_xx = buffer_1100_pppd[624];

    auto g_y_x_0_0_z_y_z_xy = buffer_1100_pppd[625];

    auto g_y_x_0_0_z_y_z_xz = buffer_1100_pppd[626];

    auto g_y_x_0_0_z_y_z_yy = buffer_1100_pppd[627];

    auto g_y_x_0_0_z_y_z_yz = buffer_1100_pppd[628];

    auto g_y_x_0_0_z_y_z_zz = buffer_1100_pppd[629];

    auto g_y_x_0_0_z_z_x_xx = buffer_1100_pppd[630];

    auto g_y_x_0_0_z_z_x_xy = buffer_1100_pppd[631];

    auto g_y_x_0_0_z_z_x_xz = buffer_1100_pppd[632];

    auto g_y_x_0_0_z_z_x_yy = buffer_1100_pppd[633];

    auto g_y_x_0_0_z_z_x_yz = buffer_1100_pppd[634];

    auto g_y_x_0_0_z_z_x_zz = buffer_1100_pppd[635];

    auto g_y_x_0_0_z_z_y_xx = buffer_1100_pppd[636];

    auto g_y_x_0_0_z_z_y_xy = buffer_1100_pppd[637];

    auto g_y_x_0_0_z_z_y_xz = buffer_1100_pppd[638];

    auto g_y_x_0_0_z_z_y_yy = buffer_1100_pppd[639];

    auto g_y_x_0_0_z_z_y_yz = buffer_1100_pppd[640];

    auto g_y_x_0_0_z_z_y_zz = buffer_1100_pppd[641];

    auto g_y_x_0_0_z_z_z_xx = buffer_1100_pppd[642];

    auto g_y_x_0_0_z_z_z_xy = buffer_1100_pppd[643];

    auto g_y_x_0_0_z_z_z_xz = buffer_1100_pppd[644];

    auto g_y_x_0_0_z_z_z_yy = buffer_1100_pppd[645];

    auto g_y_x_0_0_z_z_z_yz = buffer_1100_pppd[646];

    auto g_y_x_0_0_z_z_z_zz = buffer_1100_pppd[647];

    auto g_y_y_0_0_x_x_x_xx = buffer_1100_pppd[648];

    auto g_y_y_0_0_x_x_x_xy = buffer_1100_pppd[649];

    auto g_y_y_0_0_x_x_x_xz = buffer_1100_pppd[650];

    auto g_y_y_0_0_x_x_x_yy = buffer_1100_pppd[651];

    auto g_y_y_0_0_x_x_x_yz = buffer_1100_pppd[652];

    auto g_y_y_0_0_x_x_x_zz = buffer_1100_pppd[653];

    auto g_y_y_0_0_x_x_y_xx = buffer_1100_pppd[654];

    auto g_y_y_0_0_x_x_y_xy = buffer_1100_pppd[655];

    auto g_y_y_0_0_x_x_y_xz = buffer_1100_pppd[656];

    auto g_y_y_0_0_x_x_y_yy = buffer_1100_pppd[657];

    auto g_y_y_0_0_x_x_y_yz = buffer_1100_pppd[658];

    auto g_y_y_0_0_x_x_y_zz = buffer_1100_pppd[659];

    auto g_y_y_0_0_x_x_z_xx = buffer_1100_pppd[660];

    auto g_y_y_0_0_x_x_z_xy = buffer_1100_pppd[661];

    auto g_y_y_0_0_x_x_z_xz = buffer_1100_pppd[662];

    auto g_y_y_0_0_x_x_z_yy = buffer_1100_pppd[663];

    auto g_y_y_0_0_x_x_z_yz = buffer_1100_pppd[664];

    auto g_y_y_0_0_x_x_z_zz = buffer_1100_pppd[665];

    auto g_y_y_0_0_x_y_x_xx = buffer_1100_pppd[666];

    auto g_y_y_0_0_x_y_x_xy = buffer_1100_pppd[667];

    auto g_y_y_0_0_x_y_x_xz = buffer_1100_pppd[668];

    auto g_y_y_0_0_x_y_x_yy = buffer_1100_pppd[669];

    auto g_y_y_0_0_x_y_x_yz = buffer_1100_pppd[670];

    auto g_y_y_0_0_x_y_x_zz = buffer_1100_pppd[671];

    auto g_y_y_0_0_x_y_y_xx = buffer_1100_pppd[672];

    auto g_y_y_0_0_x_y_y_xy = buffer_1100_pppd[673];

    auto g_y_y_0_0_x_y_y_xz = buffer_1100_pppd[674];

    auto g_y_y_0_0_x_y_y_yy = buffer_1100_pppd[675];

    auto g_y_y_0_0_x_y_y_yz = buffer_1100_pppd[676];

    auto g_y_y_0_0_x_y_y_zz = buffer_1100_pppd[677];

    auto g_y_y_0_0_x_y_z_xx = buffer_1100_pppd[678];

    auto g_y_y_0_0_x_y_z_xy = buffer_1100_pppd[679];

    auto g_y_y_0_0_x_y_z_xz = buffer_1100_pppd[680];

    auto g_y_y_0_0_x_y_z_yy = buffer_1100_pppd[681];

    auto g_y_y_0_0_x_y_z_yz = buffer_1100_pppd[682];

    auto g_y_y_0_0_x_y_z_zz = buffer_1100_pppd[683];

    auto g_y_y_0_0_x_z_x_xx = buffer_1100_pppd[684];

    auto g_y_y_0_0_x_z_x_xy = buffer_1100_pppd[685];

    auto g_y_y_0_0_x_z_x_xz = buffer_1100_pppd[686];

    auto g_y_y_0_0_x_z_x_yy = buffer_1100_pppd[687];

    auto g_y_y_0_0_x_z_x_yz = buffer_1100_pppd[688];

    auto g_y_y_0_0_x_z_x_zz = buffer_1100_pppd[689];

    auto g_y_y_0_0_x_z_y_xx = buffer_1100_pppd[690];

    auto g_y_y_0_0_x_z_y_xy = buffer_1100_pppd[691];

    auto g_y_y_0_0_x_z_y_xz = buffer_1100_pppd[692];

    auto g_y_y_0_0_x_z_y_yy = buffer_1100_pppd[693];

    auto g_y_y_0_0_x_z_y_yz = buffer_1100_pppd[694];

    auto g_y_y_0_0_x_z_y_zz = buffer_1100_pppd[695];

    auto g_y_y_0_0_x_z_z_xx = buffer_1100_pppd[696];

    auto g_y_y_0_0_x_z_z_xy = buffer_1100_pppd[697];

    auto g_y_y_0_0_x_z_z_xz = buffer_1100_pppd[698];

    auto g_y_y_0_0_x_z_z_yy = buffer_1100_pppd[699];

    auto g_y_y_0_0_x_z_z_yz = buffer_1100_pppd[700];

    auto g_y_y_0_0_x_z_z_zz = buffer_1100_pppd[701];

    auto g_y_y_0_0_y_x_x_xx = buffer_1100_pppd[702];

    auto g_y_y_0_0_y_x_x_xy = buffer_1100_pppd[703];

    auto g_y_y_0_0_y_x_x_xz = buffer_1100_pppd[704];

    auto g_y_y_0_0_y_x_x_yy = buffer_1100_pppd[705];

    auto g_y_y_0_0_y_x_x_yz = buffer_1100_pppd[706];

    auto g_y_y_0_0_y_x_x_zz = buffer_1100_pppd[707];

    auto g_y_y_0_0_y_x_y_xx = buffer_1100_pppd[708];

    auto g_y_y_0_0_y_x_y_xy = buffer_1100_pppd[709];

    auto g_y_y_0_0_y_x_y_xz = buffer_1100_pppd[710];

    auto g_y_y_0_0_y_x_y_yy = buffer_1100_pppd[711];

    auto g_y_y_0_0_y_x_y_yz = buffer_1100_pppd[712];

    auto g_y_y_0_0_y_x_y_zz = buffer_1100_pppd[713];

    auto g_y_y_0_0_y_x_z_xx = buffer_1100_pppd[714];

    auto g_y_y_0_0_y_x_z_xy = buffer_1100_pppd[715];

    auto g_y_y_0_0_y_x_z_xz = buffer_1100_pppd[716];

    auto g_y_y_0_0_y_x_z_yy = buffer_1100_pppd[717];

    auto g_y_y_0_0_y_x_z_yz = buffer_1100_pppd[718];

    auto g_y_y_0_0_y_x_z_zz = buffer_1100_pppd[719];

    auto g_y_y_0_0_y_y_x_xx = buffer_1100_pppd[720];

    auto g_y_y_0_0_y_y_x_xy = buffer_1100_pppd[721];

    auto g_y_y_0_0_y_y_x_xz = buffer_1100_pppd[722];

    auto g_y_y_0_0_y_y_x_yy = buffer_1100_pppd[723];

    auto g_y_y_0_0_y_y_x_yz = buffer_1100_pppd[724];

    auto g_y_y_0_0_y_y_x_zz = buffer_1100_pppd[725];

    auto g_y_y_0_0_y_y_y_xx = buffer_1100_pppd[726];

    auto g_y_y_0_0_y_y_y_xy = buffer_1100_pppd[727];

    auto g_y_y_0_0_y_y_y_xz = buffer_1100_pppd[728];

    auto g_y_y_0_0_y_y_y_yy = buffer_1100_pppd[729];

    auto g_y_y_0_0_y_y_y_yz = buffer_1100_pppd[730];

    auto g_y_y_0_0_y_y_y_zz = buffer_1100_pppd[731];

    auto g_y_y_0_0_y_y_z_xx = buffer_1100_pppd[732];

    auto g_y_y_0_0_y_y_z_xy = buffer_1100_pppd[733];

    auto g_y_y_0_0_y_y_z_xz = buffer_1100_pppd[734];

    auto g_y_y_0_0_y_y_z_yy = buffer_1100_pppd[735];

    auto g_y_y_0_0_y_y_z_yz = buffer_1100_pppd[736];

    auto g_y_y_0_0_y_y_z_zz = buffer_1100_pppd[737];

    auto g_y_y_0_0_y_z_x_xx = buffer_1100_pppd[738];

    auto g_y_y_0_0_y_z_x_xy = buffer_1100_pppd[739];

    auto g_y_y_0_0_y_z_x_xz = buffer_1100_pppd[740];

    auto g_y_y_0_0_y_z_x_yy = buffer_1100_pppd[741];

    auto g_y_y_0_0_y_z_x_yz = buffer_1100_pppd[742];

    auto g_y_y_0_0_y_z_x_zz = buffer_1100_pppd[743];

    auto g_y_y_0_0_y_z_y_xx = buffer_1100_pppd[744];

    auto g_y_y_0_0_y_z_y_xy = buffer_1100_pppd[745];

    auto g_y_y_0_0_y_z_y_xz = buffer_1100_pppd[746];

    auto g_y_y_0_0_y_z_y_yy = buffer_1100_pppd[747];

    auto g_y_y_0_0_y_z_y_yz = buffer_1100_pppd[748];

    auto g_y_y_0_0_y_z_y_zz = buffer_1100_pppd[749];

    auto g_y_y_0_0_y_z_z_xx = buffer_1100_pppd[750];

    auto g_y_y_0_0_y_z_z_xy = buffer_1100_pppd[751];

    auto g_y_y_0_0_y_z_z_xz = buffer_1100_pppd[752];

    auto g_y_y_0_0_y_z_z_yy = buffer_1100_pppd[753];

    auto g_y_y_0_0_y_z_z_yz = buffer_1100_pppd[754];

    auto g_y_y_0_0_y_z_z_zz = buffer_1100_pppd[755];

    auto g_y_y_0_0_z_x_x_xx = buffer_1100_pppd[756];

    auto g_y_y_0_0_z_x_x_xy = buffer_1100_pppd[757];

    auto g_y_y_0_0_z_x_x_xz = buffer_1100_pppd[758];

    auto g_y_y_0_0_z_x_x_yy = buffer_1100_pppd[759];

    auto g_y_y_0_0_z_x_x_yz = buffer_1100_pppd[760];

    auto g_y_y_0_0_z_x_x_zz = buffer_1100_pppd[761];

    auto g_y_y_0_0_z_x_y_xx = buffer_1100_pppd[762];

    auto g_y_y_0_0_z_x_y_xy = buffer_1100_pppd[763];

    auto g_y_y_0_0_z_x_y_xz = buffer_1100_pppd[764];

    auto g_y_y_0_0_z_x_y_yy = buffer_1100_pppd[765];

    auto g_y_y_0_0_z_x_y_yz = buffer_1100_pppd[766];

    auto g_y_y_0_0_z_x_y_zz = buffer_1100_pppd[767];

    auto g_y_y_0_0_z_x_z_xx = buffer_1100_pppd[768];

    auto g_y_y_0_0_z_x_z_xy = buffer_1100_pppd[769];

    auto g_y_y_0_0_z_x_z_xz = buffer_1100_pppd[770];

    auto g_y_y_0_0_z_x_z_yy = buffer_1100_pppd[771];

    auto g_y_y_0_0_z_x_z_yz = buffer_1100_pppd[772];

    auto g_y_y_0_0_z_x_z_zz = buffer_1100_pppd[773];

    auto g_y_y_0_0_z_y_x_xx = buffer_1100_pppd[774];

    auto g_y_y_0_0_z_y_x_xy = buffer_1100_pppd[775];

    auto g_y_y_0_0_z_y_x_xz = buffer_1100_pppd[776];

    auto g_y_y_0_0_z_y_x_yy = buffer_1100_pppd[777];

    auto g_y_y_0_0_z_y_x_yz = buffer_1100_pppd[778];

    auto g_y_y_0_0_z_y_x_zz = buffer_1100_pppd[779];

    auto g_y_y_0_0_z_y_y_xx = buffer_1100_pppd[780];

    auto g_y_y_0_0_z_y_y_xy = buffer_1100_pppd[781];

    auto g_y_y_0_0_z_y_y_xz = buffer_1100_pppd[782];

    auto g_y_y_0_0_z_y_y_yy = buffer_1100_pppd[783];

    auto g_y_y_0_0_z_y_y_yz = buffer_1100_pppd[784];

    auto g_y_y_0_0_z_y_y_zz = buffer_1100_pppd[785];

    auto g_y_y_0_0_z_y_z_xx = buffer_1100_pppd[786];

    auto g_y_y_0_0_z_y_z_xy = buffer_1100_pppd[787];

    auto g_y_y_0_0_z_y_z_xz = buffer_1100_pppd[788];

    auto g_y_y_0_0_z_y_z_yy = buffer_1100_pppd[789];

    auto g_y_y_0_0_z_y_z_yz = buffer_1100_pppd[790];

    auto g_y_y_0_0_z_y_z_zz = buffer_1100_pppd[791];

    auto g_y_y_0_0_z_z_x_xx = buffer_1100_pppd[792];

    auto g_y_y_0_0_z_z_x_xy = buffer_1100_pppd[793];

    auto g_y_y_0_0_z_z_x_xz = buffer_1100_pppd[794];

    auto g_y_y_0_0_z_z_x_yy = buffer_1100_pppd[795];

    auto g_y_y_0_0_z_z_x_yz = buffer_1100_pppd[796];

    auto g_y_y_0_0_z_z_x_zz = buffer_1100_pppd[797];

    auto g_y_y_0_0_z_z_y_xx = buffer_1100_pppd[798];

    auto g_y_y_0_0_z_z_y_xy = buffer_1100_pppd[799];

    auto g_y_y_0_0_z_z_y_xz = buffer_1100_pppd[800];

    auto g_y_y_0_0_z_z_y_yy = buffer_1100_pppd[801];

    auto g_y_y_0_0_z_z_y_yz = buffer_1100_pppd[802];

    auto g_y_y_0_0_z_z_y_zz = buffer_1100_pppd[803];

    auto g_y_y_0_0_z_z_z_xx = buffer_1100_pppd[804];

    auto g_y_y_0_0_z_z_z_xy = buffer_1100_pppd[805];

    auto g_y_y_0_0_z_z_z_xz = buffer_1100_pppd[806];

    auto g_y_y_0_0_z_z_z_yy = buffer_1100_pppd[807];

    auto g_y_y_0_0_z_z_z_yz = buffer_1100_pppd[808];

    auto g_y_y_0_0_z_z_z_zz = buffer_1100_pppd[809];

    auto g_y_z_0_0_x_x_x_xx = buffer_1100_pppd[810];

    auto g_y_z_0_0_x_x_x_xy = buffer_1100_pppd[811];

    auto g_y_z_0_0_x_x_x_xz = buffer_1100_pppd[812];

    auto g_y_z_0_0_x_x_x_yy = buffer_1100_pppd[813];

    auto g_y_z_0_0_x_x_x_yz = buffer_1100_pppd[814];

    auto g_y_z_0_0_x_x_x_zz = buffer_1100_pppd[815];

    auto g_y_z_0_0_x_x_y_xx = buffer_1100_pppd[816];

    auto g_y_z_0_0_x_x_y_xy = buffer_1100_pppd[817];

    auto g_y_z_0_0_x_x_y_xz = buffer_1100_pppd[818];

    auto g_y_z_0_0_x_x_y_yy = buffer_1100_pppd[819];

    auto g_y_z_0_0_x_x_y_yz = buffer_1100_pppd[820];

    auto g_y_z_0_0_x_x_y_zz = buffer_1100_pppd[821];

    auto g_y_z_0_0_x_x_z_xx = buffer_1100_pppd[822];

    auto g_y_z_0_0_x_x_z_xy = buffer_1100_pppd[823];

    auto g_y_z_0_0_x_x_z_xz = buffer_1100_pppd[824];

    auto g_y_z_0_0_x_x_z_yy = buffer_1100_pppd[825];

    auto g_y_z_0_0_x_x_z_yz = buffer_1100_pppd[826];

    auto g_y_z_0_0_x_x_z_zz = buffer_1100_pppd[827];

    auto g_y_z_0_0_x_y_x_xx = buffer_1100_pppd[828];

    auto g_y_z_0_0_x_y_x_xy = buffer_1100_pppd[829];

    auto g_y_z_0_0_x_y_x_xz = buffer_1100_pppd[830];

    auto g_y_z_0_0_x_y_x_yy = buffer_1100_pppd[831];

    auto g_y_z_0_0_x_y_x_yz = buffer_1100_pppd[832];

    auto g_y_z_0_0_x_y_x_zz = buffer_1100_pppd[833];

    auto g_y_z_0_0_x_y_y_xx = buffer_1100_pppd[834];

    auto g_y_z_0_0_x_y_y_xy = buffer_1100_pppd[835];

    auto g_y_z_0_0_x_y_y_xz = buffer_1100_pppd[836];

    auto g_y_z_0_0_x_y_y_yy = buffer_1100_pppd[837];

    auto g_y_z_0_0_x_y_y_yz = buffer_1100_pppd[838];

    auto g_y_z_0_0_x_y_y_zz = buffer_1100_pppd[839];

    auto g_y_z_0_0_x_y_z_xx = buffer_1100_pppd[840];

    auto g_y_z_0_0_x_y_z_xy = buffer_1100_pppd[841];

    auto g_y_z_0_0_x_y_z_xz = buffer_1100_pppd[842];

    auto g_y_z_0_0_x_y_z_yy = buffer_1100_pppd[843];

    auto g_y_z_0_0_x_y_z_yz = buffer_1100_pppd[844];

    auto g_y_z_0_0_x_y_z_zz = buffer_1100_pppd[845];

    auto g_y_z_0_0_x_z_x_xx = buffer_1100_pppd[846];

    auto g_y_z_0_0_x_z_x_xy = buffer_1100_pppd[847];

    auto g_y_z_0_0_x_z_x_xz = buffer_1100_pppd[848];

    auto g_y_z_0_0_x_z_x_yy = buffer_1100_pppd[849];

    auto g_y_z_0_0_x_z_x_yz = buffer_1100_pppd[850];

    auto g_y_z_0_0_x_z_x_zz = buffer_1100_pppd[851];

    auto g_y_z_0_0_x_z_y_xx = buffer_1100_pppd[852];

    auto g_y_z_0_0_x_z_y_xy = buffer_1100_pppd[853];

    auto g_y_z_0_0_x_z_y_xz = buffer_1100_pppd[854];

    auto g_y_z_0_0_x_z_y_yy = buffer_1100_pppd[855];

    auto g_y_z_0_0_x_z_y_yz = buffer_1100_pppd[856];

    auto g_y_z_0_0_x_z_y_zz = buffer_1100_pppd[857];

    auto g_y_z_0_0_x_z_z_xx = buffer_1100_pppd[858];

    auto g_y_z_0_0_x_z_z_xy = buffer_1100_pppd[859];

    auto g_y_z_0_0_x_z_z_xz = buffer_1100_pppd[860];

    auto g_y_z_0_0_x_z_z_yy = buffer_1100_pppd[861];

    auto g_y_z_0_0_x_z_z_yz = buffer_1100_pppd[862];

    auto g_y_z_0_0_x_z_z_zz = buffer_1100_pppd[863];

    auto g_y_z_0_0_y_x_x_xx = buffer_1100_pppd[864];

    auto g_y_z_0_0_y_x_x_xy = buffer_1100_pppd[865];

    auto g_y_z_0_0_y_x_x_xz = buffer_1100_pppd[866];

    auto g_y_z_0_0_y_x_x_yy = buffer_1100_pppd[867];

    auto g_y_z_0_0_y_x_x_yz = buffer_1100_pppd[868];

    auto g_y_z_0_0_y_x_x_zz = buffer_1100_pppd[869];

    auto g_y_z_0_0_y_x_y_xx = buffer_1100_pppd[870];

    auto g_y_z_0_0_y_x_y_xy = buffer_1100_pppd[871];

    auto g_y_z_0_0_y_x_y_xz = buffer_1100_pppd[872];

    auto g_y_z_0_0_y_x_y_yy = buffer_1100_pppd[873];

    auto g_y_z_0_0_y_x_y_yz = buffer_1100_pppd[874];

    auto g_y_z_0_0_y_x_y_zz = buffer_1100_pppd[875];

    auto g_y_z_0_0_y_x_z_xx = buffer_1100_pppd[876];

    auto g_y_z_0_0_y_x_z_xy = buffer_1100_pppd[877];

    auto g_y_z_0_0_y_x_z_xz = buffer_1100_pppd[878];

    auto g_y_z_0_0_y_x_z_yy = buffer_1100_pppd[879];

    auto g_y_z_0_0_y_x_z_yz = buffer_1100_pppd[880];

    auto g_y_z_0_0_y_x_z_zz = buffer_1100_pppd[881];

    auto g_y_z_0_0_y_y_x_xx = buffer_1100_pppd[882];

    auto g_y_z_0_0_y_y_x_xy = buffer_1100_pppd[883];

    auto g_y_z_0_0_y_y_x_xz = buffer_1100_pppd[884];

    auto g_y_z_0_0_y_y_x_yy = buffer_1100_pppd[885];

    auto g_y_z_0_0_y_y_x_yz = buffer_1100_pppd[886];

    auto g_y_z_0_0_y_y_x_zz = buffer_1100_pppd[887];

    auto g_y_z_0_0_y_y_y_xx = buffer_1100_pppd[888];

    auto g_y_z_0_0_y_y_y_xy = buffer_1100_pppd[889];

    auto g_y_z_0_0_y_y_y_xz = buffer_1100_pppd[890];

    auto g_y_z_0_0_y_y_y_yy = buffer_1100_pppd[891];

    auto g_y_z_0_0_y_y_y_yz = buffer_1100_pppd[892];

    auto g_y_z_0_0_y_y_y_zz = buffer_1100_pppd[893];

    auto g_y_z_0_0_y_y_z_xx = buffer_1100_pppd[894];

    auto g_y_z_0_0_y_y_z_xy = buffer_1100_pppd[895];

    auto g_y_z_0_0_y_y_z_xz = buffer_1100_pppd[896];

    auto g_y_z_0_0_y_y_z_yy = buffer_1100_pppd[897];

    auto g_y_z_0_0_y_y_z_yz = buffer_1100_pppd[898];

    auto g_y_z_0_0_y_y_z_zz = buffer_1100_pppd[899];

    auto g_y_z_0_0_y_z_x_xx = buffer_1100_pppd[900];

    auto g_y_z_0_0_y_z_x_xy = buffer_1100_pppd[901];

    auto g_y_z_0_0_y_z_x_xz = buffer_1100_pppd[902];

    auto g_y_z_0_0_y_z_x_yy = buffer_1100_pppd[903];

    auto g_y_z_0_0_y_z_x_yz = buffer_1100_pppd[904];

    auto g_y_z_0_0_y_z_x_zz = buffer_1100_pppd[905];

    auto g_y_z_0_0_y_z_y_xx = buffer_1100_pppd[906];

    auto g_y_z_0_0_y_z_y_xy = buffer_1100_pppd[907];

    auto g_y_z_0_0_y_z_y_xz = buffer_1100_pppd[908];

    auto g_y_z_0_0_y_z_y_yy = buffer_1100_pppd[909];

    auto g_y_z_0_0_y_z_y_yz = buffer_1100_pppd[910];

    auto g_y_z_0_0_y_z_y_zz = buffer_1100_pppd[911];

    auto g_y_z_0_0_y_z_z_xx = buffer_1100_pppd[912];

    auto g_y_z_0_0_y_z_z_xy = buffer_1100_pppd[913];

    auto g_y_z_0_0_y_z_z_xz = buffer_1100_pppd[914];

    auto g_y_z_0_0_y_z_z_yy = buffer_1100_pppd[915];

    auto g_y_z_0_0_y_z_z_yz = buffer_1100_pppd[916];

    auto g_y_z_0_0_y_z_z_zz = buffer_1100_pppd[917];

    auto g_y_z_0_0_z_x_x_xx = buffer_1100_pppd[918];

    auto g_y_z_0_0_z_x_x_xy = buffer_1100_pppd[919];

    auto g_y_z_0_0_z_x_x_xz = buffer_1100_pppd[920];

    auto g_y_z_0_0_z_x_x_yy = buffer_1100_pppd[921];

    auto g_y_z_0_0_z_x_x_yz = buffer_1100_pppd[922];

    auto g_y_z_0_0_z_x_x_zz = buffer_1100_pppd[923];

    auto g_y_z_0_0_z_x_y_xx = buffer_1100_pppd[924];

    auto g_y_z_0_0_z_x_y_xy = buffer_1100_pppd[925];

    auto g_y_z_0_0_z_x_y_xz = buffer_1100_pppd[926];

    auto g_y_z_0_0_z_x_y_yy = buffer_1100_pppd[927];

    auto g_y_z_0_0_z_x_y_yz = buffer_1100_pppd[928];

    auto g_y_z_0_0_z_x_y_zz = buffer_1100_pppd[929];

    auto g_y_z_0_0_z_x_z_xx = buffer_1100_pppd[930];

    auto g_y_z_0_0_z_x_z_xy = buffer_1100_pppd[931];

    auto g_y_z_0_0_z_x_z_xz = buffer_1100_pppd[932];

    auto g_y_z_0_0_z_x_z_yy = buffer_1100_pppd[933];

    auto g_y_z_0_0_z_x_z_yz = buffer_1100_pppd[934];

    auto g_y_z_0_0_z_x_z_zz = buffer_1100_pppd[935];

    auto g_y_z_0_0_z_y_x_xx = buffer_1100_pppd[936];

    auto g_y_z_0_0_z_y_x_xy = buffer_1100_pppd[937];

    auto g_y_z_0_0_z_y_x_xz = buffer_1100_pppd[938];

    auto g_y_z_0_0_z_y_x_yy = buffer_1100_pppd[939];

    auto g_y_z_0_0_z_y_x_yz = buffer_1100_pppd[940];

    auto g_y_z_0_0_z_y_x_zz = buffer_1100_pppd[941];

    auto g_y_z_0_0_z_y_y_xx = buffer_1100_pppd[942];

    auto g_y_z_0_0_z_y_y_xy = buffer_1100_pppd[943];

    auto g_y_z_0_0_z_y_y_xz = buffer_1100_pppd[944];

    auto g_y_z_0_0_z_y_y_yy = buffer_1100_pppd[945];

    auto g_y_z_0_0_z_y_y_yz = buffer_1100_pppd[946];

    auto g_y_z_0_0_z_y_y_zz = buffer_1100_pppd[947];

    auto g_y_z_0_0_z_y_z_xx = buffer_1100_pppd[948];

    auto g_y_z_0_0_z_y_z_xy = buffer_1100_pppd[949];

    auto g_y_z_0_0_z_y_z_xz = buffer_1100_pppd[950];

    auto g_y_z_0_0_z_y_z_yy = buffer_1100_pppd[951];

    auto g_y_z_0_0_z_y_z_yz = buffer_1100_pppd[952];

    auto g_y_z_0_0_z_y_z_zz = buffer_1100_pppd[953];

    auto g_y_z_0_0_z_z_x_xx = buffer_1100_pppd[954];

    auto g_y_z_0_0_z_z_x_xy = buffer_1100_pppd[955];

    auto g_y_z_0_0_z_z_x_xz = buffer_1100_pppd[956];

    auto g_y_z_0_0_z_z_x_yy = buffer_1100_pppd[957];

    auto g_y_z_0_0_z_z_x_yz = buffer_1100_pppd[958];

    auto g_y_z_0_0_z_z_x_zz = buffer_1100_pppd[959];

    auto g_y_z_0_0_z_z_y_xx = buffer_1100_pppd[960];

    auto g_y_z_0_0_z_z_y_xy = buffer_1100_pppd[961];

    auto g_y_z_0_0_z_z_y_xz = buffer_1100_pppd[962];

    auto g_y_z_0_0_z_z_y_yy = buffer_1100_pppd[963];

    auto g_y_z_0_0_z_z_y_yz = buffer_1100_pppd[964];

    auto g_y_z_0_0_z_z_y_zz = buffer_1100_pppd[965];

    auto g_y_z_0_0_z_z_z_xx = buffer_1100_pppd[966];

    auto g_y_z_0_0_z_z_z_xy = buffer_1100_pppd[967];

    auto g_y_z_0_0_z_z_z_xz = buffer_1100_pppd[968];

    auto g_y_z_0_0_z_z_z_yy = buffer_1100_pppd[969];

    auto g_y_z_0_0_z_z_z_yz = buffer_1100_pppd[970];

    auto g_y_z_0_0_z_z_z_zz = buffer_1100_pppd[971];

    auto g_z_x_0_0_x_x_x_xx = buffer_1100_pppd[972];

    auto g_z_x_0_0_x_x_x_xy = buffer_1100_pppd[973];

    auto g_z_x_0_0_x_x_x_xz = buffer_1100_pppd[974];

    auto g_z_x_0_0_x_x_x_yy = buffer_1100_pppd[975];

    auto g_z_x_0_0_x_x_x_yz = buffer_1100_pppd[976];

    auto g_z_x_0_0_x_x_x_zz = buffer_1100_pppd[977];

    auto g_z_x_0_0_x_x_y_xx = buffer_1100_pppd[978];

    auto g_z_x_0_0_x_x_y_xy = buffer_1100_pppd[979];

    auto g_z_x_0_0_x_x_y_xz = buffer_1100_pppd[980];

    auto g_z_x_0_0_x_x_y_yy = buffer_1100_pppd[981];

    auto g_z_x_0_0_x_x_y_yz = buffer_1100_pppd[982];

    auto g_z_x_0_0_x_x_y_zz = buffer_1100_pppd[983];

    auto g_z_x_0_0_x_x_z_xx = buffer_1100_pppd[984];

    auto g_z_x_0_0_x_x_z_xy = buffer_1100_pppd[985];

    auto g_z_x_0_0_x_x_z_xz = buffer_1100_pppd[986];

    auto g_z_x_0_0_x_x_z_yy = buffer_1100_pppd[987];

    auto g_z_x_0_0_x_x_z_yz = buffer_1100_pppd[988];

    auto g_z_x_0_0_x_x_z_zz = buffer_1100_pppd[989];

    auto g_z_x_0_0_x_y_x_xx = buffer_1100_pppd[990];

    auto g_z_x_0_0_x_y_x_xy = buffer_1100_pppd[991];

    auto g_z_x_0_0_x_y_x_xz = buffer_1100_pppd[992];

    auto g_z_x_0_0_x_y_x_yy = buffer_1100_pppd[993];

    auto g_z_x_0_0_x_y_x_yz = buffer_1100_pppd[994];

    auto g_z_x_0_0_x_y_x_zz = buffer_1100_pppd[995];

    auto g_z_x_0_0_x_y_y_xx = buffer_1100_pppd[996];

    auto g_z_x_0_0_x_y_y_xy = buffer_1100_pppd[997];

    auto g_z_x_0_0_x_y_y_xz = buffer_1100_pppd[998];

    auto g_z_x_0_0_x_y_y_yy = buffer_1100_pppd[999];

    auto g_z_x_0_0_x_y_y_yz = buffer_1100_pppd[1000];

    auto g_z_x_0_0_x_y_y_zz = buffer_1100_pppd[1001];

    auto g_z_x_0_0_x_y_z_xx = buffer_1100_pppd[1002];

    auto g_z_x_0_0_x_y_z_xy = buffer_1100_pppd[1003];

    auto g_z_x_0_0_x_y_z_xz = buffer_1100_pppd[1004];

    auto g_z_x_0_0_x_y_z_yy = buffer_1100_pppd[1005];

    auto g_z_x_0_0_x_y_z_yz = buffer_1100_pppd[1006];

    auto g_z_x_0_0_x_y_z_zz = buffer_1100_pppd[1007];

    auto g_z_x_0_0_x_z_x_xx = buffer_1100_pppd[1008];

    auto g_z_x_0_0_x_z_x_xy = buffer_1100_pppd[1009];

    auto g_z_x_0_0_x_z_x_xz = buffer_1100_pppd[1010];

    auto g_z_x_0_0_x_z_x_yy = buffer_1100_pppd[1011];

    auto g_z_x_0_0_x_z_x_yz = buffer_1100_pppd[1012];

    auto g_z_x_0_0_x_z_x_zz = buffer_1100_pppd[1013];

    auto g_z_x_0_0_x_z_y_xx = buffer_1100_pppd[1014];

    auto g_z_x_0_0_x_z_y_xy = buffer_1100_pppd[1015];

    auto g_z_x_0_0_x_z_y_xz = buffer_1100_pppd[1016];

    auto g_z_x_0_0_x_z_y_yy = buffer_1100_pppd[1017];

    auto g_z_x_0_0_x_z_y_yz = buffer_1100_pppd[1018];

    auto g_z_x_0_0_x_z_y_zz = buffer_1100_pppd[1019];

    auto g_z_x_0_0_x_z_z_xx = buffer_1100_pppd[1020];

    auto g_z_x_0_0_x_z_z_xy = buffer_1100_pppd[1021];

    auto g_z_x_0_0_x_z_z_xz = buffer_1100_pppd[1022];

    auto g_z_x_0_0_x_z_z_yy = buffer_1100_pppd[1023];

    auto g_z_x_0_0_x_z_z_yz = buffer_1100_pppd[1024];

    auto g_z_x_0_0_x_z_z_zz = buffer_1100_pppd[1025];

    auto g_z_x_0_0_y_x_x_xx = buffer_1100_pppd[1026];

    auto g_z_x_0_0_y_x_x_xy = buffer_1100_pppd[1027];

    auto g_z_x_0_0_y_x_x_xz = buffer_1100_pppd[1028];

    auto g_z_x_0_0_y_x_x_yy = buffer_1100_pppd[1029];

    auto g_z_x_0_0_y_x_x_yz = buffer_1100_pppd[1030];

    auto g_z_x_0_0_y_x_x_zz = buffer_1100_pppd[1031];

    auto g_z_x_0_0_y_x_y_xx = buffer_1100_pppd[1032];

    auto g_z_x_0_0_y_x_y_xy = buffer_1100_pppd[1033];

    auto g_z_x_0_0_y_x_y_xz = buffer_1100_pppd[1034];

    auto g_z_x_0_0_y_x_y_yy = buffer_1100_pppd[1035];

    auto g_z_x_0_0_y_x_y_yz = buffer_1100_pppd[1036];

    auto g_z_x_0_0_y_x_y_zz = buffer_1100_pppd[1037];

    auto g_z_x_0_0_y_x_z_xx = buffer_1100_pppd[1038];

    auto g_z_x_0_0_y_x_z_xy = buffer_1100_pppd[1039];

    auto g_z_x_0_0_y_x_z_xz = buffer_1100_pppd[1040];

    auto g_z_x_0_0_y_x_z_yy = buffer_1100_pppd[1041];

    auto g_z_x_0_0_y_x_z_yz = buffer_1100_pppd[1042];

    auto g_z_x_0_0_y_x_z_zz = buffer_1100_pppd[1043];

    auto g_z_x_0_0_y_y_x_xx = buffer_1100_pppd[1044];

    auto g_z_x_0_0_y_y_x_xy = buffer_1100_pppd[1045];

    auto g_z_x_0_0_y_y_x_xz = buffer_1100_pppd[1046];

    auto g_z_x_0_0_y_y_x_yy = buffer_1100_pppd[1047];

    auto g_z_x_0_0_y_y_x_yz = buffer_1100_pppd[1048];

    auto g_z_x_0_0_y_y_x_zz = buffer_1100_pppd[1049];

    auto g_z_x_0_0_y_y_y_xx = buffer_1100_pppd[1050];

    auto g_z_x_0_0_y_y_y_xy = buffer_1100_pppd[1051];

    auto g_z_x_0_0_y_y_y_xz = buffer_1100_pppd[1052];

    auto g_z_x_0_0_y_y_y_yy = buffer_1100_pppd[1053];

    auto g_z_x_0_0_y_y_y_yz = buffer_1100_pppd[1054];

    auto g_z_x_0_0_y_y_y_zz = buffer_1100_pppd[1055];

    auto g_z_x_0_0_y_y_z_xx = buffer_1100_pppd[1056];

    auto g_z_x_0_0_y_y_z_xy = buffer_1100_pppd[1057];

    auto g_z_x_0_0_y_y_z_xz = buffer_1100_pppd[1058];

    auto g_z_x_0_0_y_y_z_yy = buffer_1100_pppd[1059];

    auto g_z_x_0_0_y_y_z_yz = buffer_1100_pppd[1060];

    auto g_z_x_0_0_y_y_z_zz = buffer_1100_pppd[1061];

    auto g_z_x_0_0_y_z_x_xx = buffer_1100_pppd[1062];

    auto g_z_x_0_0_y_z_x_xy = buffer_1100_pppd[1063];

    auto g_z_x_0_0_y_z_x_xz = buffer_1100_pppd[1064];

    auto g_z_x_0_0_y_z_x_yy = buffer_1100_pppd[1065];

    auto g_z_x_0_0_y_z_x_yz = buffer_1100_pppd[1066];

    auto g_z_x_0_0_y_z_x_zz = buffer_1100_pppd[1067];

    auto g_z_x_0_0_y_z_y_xx = buffer_1100_pppd[1068];

    auto g_z_x_0_0_y_z_y_xy = buffer_1100_pppd[1069];

    auto g_z_x_0_0_y_z_y_xz = buffer_1100_pppd[1070];

    auto g_z_x_0_0_y_z_y_yy = buffer_1100_pppd[1071];

    auto g_z_x_0_0_y_z_y_yz = buffer_1100_pppd[1072];

    auto g_z_x_0_0_y_z_y_zz = buffer_1100_pppd[1073];

    auto g_z_x_0_0_y_z_z_xx = buffer_1100_pppd[1074];

    auto g_z_x_0_0_y_z_z_xy = buffer_1100_pppd[1075];

    auto g_z_x_0_0_y_z_z_xz = buffer_1100_pppd[1076];

    auto g_z_x_0_0_y_z_z_yy = buffer_1100_pppd[1077];

    auto g_z_x_0_0_y_z_z_yz = buffer_1100_pppd[1078];

    auto g_z_x_0_0_y_z_z_zz = buffer_1100_pppd[1079];

    auto g_z_x_0_0_z_x_x_xx = buffer_1100_pppd[1080];

    auto g_z_x_0_0_z_x_x_xy = buffer_1100_pppd[1081];

    auto g_z_x_0_0_z_x_x_xz = buffer_1100_pppd[1082];

    auto g_z_x_0_0_z_x_x_yy = buffer_1100_pppd[1083];

    auto g_z_x_0_0_z_x_x_yz = buffer_1100_pppd[1084];

    auto g_z_x_0_0_z_x_x_zz = buffer_1100_pppd[1085];

    auto g_z_x_0_0_z_x_y_xx = buffer_1100_pppd[1086];

    auto g_z_x_0_0_z_x_y_xy = buffer_1100_pppd[1087];

    auto g_z_x_0_0_z_x_y_xz = buffer_1100_pppd[1088];

    auto g_z_x_0_0_z_x_y_yy = buffer_1100_pppd[1089];

    auto g_z_x_0_0_z_x_y_yz = buffer_1100_pppd[1090];

    auto g_z_x_0_0_z_x_y_zz = buffer_1100_pppd[1091];

    auto g_z_x_0_0_z_x_z_xx = buffer_1100_pppd[1092];

    auto g_z_x_0_0_z_x_z_xy = buffer_1100_pppd[1093];

    auto g_z_x_0_0_z_x_z_xz = buffer_1100_pppd[1094];

    auto g_z_x_0_0_z_x_z_yy = buffer_1100_pppd[1095];

    auto g_z_x_0_0_z_x_z_yz = buffer_1100_pppd[1096];

    auto g_z_x_0_0_z_x_z_zz = buffer_1100_pppd[1097];

    auto g_z_x_0_0_z_y_x_xx = buffer_1100_pppd[1098];

    auto g_z_x_0_0_z_y_x_xy = buffer_1100_pppd[1099];

    auto g_z_x_0_0_z_y_x_xz = buffer_1100_pppd[1100];

    auto g_z_x_0_0_z_y_x_yy = buffer_1100_pppd[1101];

    auto g_z_x_0_0_z_y_x_yz = buffer_1100_pppd[1102];

    auto g_z_x_0_0_z_y_x_zz = buffer_1100_pppd[1103];

    auto g_z_x_0_0_z_y_y_xx = buffer_1100_pppd[1104];

    auto g_z_x_0_0_z_y_y_xy = buffer_1100_pppd[1105];

    auto g_z_x_0_0_z_y_y_xz = buffer_1100_pppd[1106];

    auto g_z_x_0_0_z_y_y_yy = buffer_1100_pppd[1107];

    auto g_z_x_0_0_z_y_y_yz = buffer_1100_pppd[1108];

    auto g_z_x_0_0_z_y_y_zz = buffer_1100_pppd[1109];

    auto g_z_x_0_0_z_y_z_xx = buffer_1100_pppd[1110];

    auto g_z_x_0_0_z_y_z_xy = buffer_1100_pppd[1111];

    auto g_z_x_0_0_z_y_z_xz = buffer_1100_pppd[1112];

    auto g_z_x_0_0_z_y_z_yy = buffer_1100_pppd[1113];

    auto g_z_x_0_0_z_y_z_yz = buffer_1100_pppd[1114];

    auto g_z_x_0_0_z_y_z_zz = buffer_1100_pppd[1115];

    auto g_z_x_0_0_z_z_x_xx = buffer_1100_pppd[1116];

    auto g_z_x_0_0_z_z_x_xy = buffer_1100_pppd[1117];

    auto g_z_x_0_0_z_z_x_xz = buffer_1100_pppd[1118];

    auto g_z_x_0_0_z_z_x_yy = buffer_1100_pppd[1119];

    auto g_z_x_0_0_z_z_x_yz = buffer_1100_pppd[1120];

    auto g_z_x_0_0_z_z_x_zz = buffer_1100_pppd[1121];

    auto g_z_x_0_0_z_z_y_xx = buffer_1100_pppd[1122];

    auto g_z_x_0_0_z_z_y_xy = buffer_1100_pppd[1123];

    auto g_z_x_0_0_z_z_y_xz = buffer_1100_pppd[1124];

    auto g_z_x_0_0_z_z_y_yy = buffer_1100_pppd[1125];

    auto g_z_x_0_0_z_z_y_yz = buffer_1100_pppd[1126];

    auto g_z_x_0_0_z_z_y_zz = buffer_1100_pppd[1127];

    auto g_z_x_0_0_z_z_z_xx = buffer_1100_pppd[1128];

    auto g_z_x_0_0_z_z_z_xy = buffer_1100_pppd[1129];

    auto g_z_x_0_0_z_z_z_xz = buffer_1100_pppd[1130];

    auto g_z_x_0_0_z_z_z_yy = buffer_1100_pppd[1131];

    auto g_z_x_0_0_z_z_z_yz = buffer_1100_pppd[1132];

    auto g_z_x_0_0_z_z_z_zz = buffer_1100_pppd[1133];

    auto g_z_y_0_0_x_x_x_xx = buffer_1100_pppd[1134];

    auto g_z_y_0_0_x_x_x_xy = buffer_1100_pppd[1135];

    auto g_z_y_0_0_x_x_x_xz = buffer_1100_pppd[1136];

    auto g_z_y_0_0_x_x_x_yy = buffer_1100_pppd[1137];

    auto g_z_y_0_0_x_x_x_yz = buffer_1100_pppd[1138];

    auto g_z_y_0_0_x_x_x_zz = buffer_1100_pppd[1139];

    auto g_z_y_0_0_x_x_y_xx = buffer_1100_pppd[1140];

    auto g_z_y_0_0_x_x_y_xy = buffer_1100_pppd[1141];

    auto g_z_y_0_0_x_x_y_xz = buffer_1100_pppd[1142];

    auto g_z_y_0_0_x_x_y_yy = buffer_1100_pppd[1143];

    auto g_z_y_0_0_x_x_y_yz = buffer_1100_pppd[1144];

    auto g_z_y_0_0_x_x_y_zz = buffer_1100_pppd[1145];

    auto g_z_y_0_0_x_x_z_xx = buffer_1100_pppd[1146];

    auto g_z_y_0_0_x_x_z_xy = buffer_1100_pppd[1147];

    auto g_z_y_0_0_x_x_z_xz = buffer_1100_pppd[1148];

    auto g_z_y_0_0_x_x_z_yy = buffer_1100_pppd[1149];

    auto g_z_y_0_0_x_x_z_yz = buffer_1100_pppd[1150];

    auto g_z_y_0_0_x_x_z_zz = buffer_1100_pppd[1151];

    auto g_z_y_0_0_x_y_x_xx = buffer_1100_pppd[1152];

    auto g_z_y_0_0_x_y_x_xy = buffer_1100_pppd[1153];

    auto g_z_y_0_0_x_y_x_xz = buffer_1100_pppd[1154];

    auto g_z_y_0_0_x_y_x_yy = buffer_1100_pppd[1155];

    auto g_z_y_0_0_x_y_x_yz = buffer_1100_pppd[1156];

    auto g_z_y_0_0_x_y_x_zz = buffer_1100_pppd[1157];

    auto g_z_y_0_0_x_y_y_xx = buffer_1100_pppd[1158];

    auto g_z_y_0_0_x_y_y_xy = buffer_1100_pppd[1159];

    auto g_z_y_0_0_x_y_y_xz = buffer_1100_pppd[1160];

    auto g_z_y_0_0_x_y_y_yy = buffer_1100_pppd[1161];

    auto g_z_y_0_0_x_y_y_yz = buffer_1100_pppd[1162];

    auto g_z_y_0_0_x_y_y_zz = buffer_1100_pppd[1163];

    auto g_z_y_0_0_x_y_z_xx = buffer_1100_pppd[1164];

    auto g_z_y_0_0_x_y_z_xy = buffer_1100_pppd[1165];

    auto g_z_y_0_0_x_y_z_xz = buffer_1100_pppd[1166];

    auto g_z_y_0_0_x_y_z_yy = buffer_1100_pppd[1167];

    auto g_z_y_0_0_x_y_z_yz = buffer_1100_pppd[1168];

    auto g_z_y_0_0_x_y_z_zz = buffer_1100_pppd[1169];

    auto g_z_y_0_0_x_z_x_xx = buffer_1100_pppd[1170];

    auto g_z_y_0_0_x_z_x_xy = buffer_1100_pppd[1171];

    auto g_z_y_0_0_x_z_x_xz = buffer_1100_pppd[1172];

    auto g_z_y_0_0_x_z_x_yy = buffer_1100_pppd[1173];

    auto g_z_y_0_0_x_z_x_yz = buffer_1100_pppd[1174];

    auto g_z_y_0_0_x_z_x_zz = buffer_1100_pppd[1175];

    auto g_z_y_0_0_x_z_y_xx = buffer_1100_pppd[1176];

    auto g_z_y_0_0_x_z_y_xy = buffer_1100_pppd[1177];

    auto g_z_y_0_0_x_z_y_xz = buffer_1100_pppd[1178];

    auto g_z_y_0_0_x_z_y_yy = buffer_1100_pppd[1179];

    auto g_z_y_0_0_x_z_y_yz = buffer_1100_pppd[1180];

    auto g_z_y_0_0_x_z_y_zz = buffer_1100_pppd[1181];

    auto g_z_y_0_0_x_z_z_xx = buffer_1100_pppd[1182];

    auto g_z_y_0_0_x_z_z_xy = buffer_1100_pppd[1183];

    auto g_z_y_0_0_x_z_z_xz = buffer_1100_pppd[1184];

    auto g_z_y_0_0_x_z_z_yy = buffer_1100_pppd[1185];

    auto g_z_y_0_0_x_z_z_yz = buffer_1100_pppd[1186];

    auto g_z_y_0_0_x_z_z_zz = buffer_1100_pppd[1187];

    auto g_z_y_0_0_y_x_x_xx = buffer_1100_pppd[1188];

    auto g_z_y_0_0_y_x_x_xy = buffer_1100_pppd[1189];

    auto g_z_y_0_0_y_x_x_xz = buffer_1100_pppd[1190];

    auto g_z_y_0_0_y_x_x_yy = buffer_1100_pppd[1191];

    auto g_z_y_0_0_y_x_x_yz = buffer_1100_pppd[1192];

    auto g_z_y_0_0_y_x_x_zz = buffer_1100_pppd[1193];

    auto g_z_y_0_0_y_x_y_xx = buffer_1100_pppd[1194];

    auto g_z_y_0_0_y_x_y_xy = buffer_1100_pppd[1195];

    auto g_z_y_0_0_y_x_y_xz = buffer_1100_pppd[1196];

    auto g_z_y_0_0_y_x_y_yy = buffer_1100_pppd[1197];

    auto g_z_y_0_0_y_x_y_yz = buffer_1100_pppd[1198];

    auto g_z_y_0_0_y_x_y_zz = buffer_1100_pppd[1199];

    auto g_z_y_0_0_y_x_z_xx = buffer_1100_pppd[1200];

    auto g_z_y_0_0_y_x_z_xy = buffer_1100_pppd[1201];

    auto g_z_y_0_0_y_x_z_xz = buffer_1100_pppd[1202];

    auto g_z_y_0_0_y_x_z_yy = buffer_1100_pppd[1203];

    auto g_z_y_0_0_y_x_z_yz = buffer_1100_pppd[1204];

    auto g_z_y_0_0_y_x_z_zz = buffer_1100_pppd[1205];

    auto g_z_y_0_0_y_y_x_xx = buffer_1100_pppd[1206];

    auto g_z_y_0_0_y_y_x_xy = buffer_1100_pppd[1207];

    auto g_z_y_0_0_y_y_x_xz = buffer_1100_pppd[1208];

    auto g_z_y_0_0_y_y_x_yy = buffer_1100_pppd[1209];

    auto g_z_y_0_0_y_y_x_yz = buffer_1100_pppd[1210];

    auto g_z_y_0_0_y_y_x_zz = buffer_1100_pppd[1211];

    auto g_z_y_0_0_y_y_y_xx = buffer_1100_pppd[1212];

    auto g_z_y_0_0_y_y_y_xy = buffer_1100_pppd[1213];

    auto g_z_y_0_0_y_y_y_xz = buffer_1100_pppd[1214];

    auto g_z_y_0_0_y_y_y_yy = buffer_1100_pppd[1215];

    auto g_z_y_0_0_y_y_y_yz = buffer_1100_pppd[1216];

    auto g_z_y_0_0_y_y_y_zz = buffer_1100_pppd[1217];

    auto g_z_y_0_0_y_y_z_xx = buffer_1100_pppd[1218];

    auto g_z_y_0_0_y_y_z_xy = buffer_1100_pppd[1219];

    auto g_z_y_0_0_y_y_z_xz = buffer_1100_pppd[1220];

    auto g_z_y_0_0_y_y_z_yy = buffer_1100_pppd[1221];

    auto g_z_y_0_0_y_y_z_yz = buffer_1100_pppd[1222];

    auto g_z_y_0_0_y_y_z_zz = buffer_1100_pppd[1223];

    auto g_z_y_0_0_y_z_x_xx = buffer_1100_pppd[1224];

    auto g_z_y_0_0_y_z_x_xy = buffer_1100_pppd[1225];

    auto g_z_y_0_0_y_z_x_xz = buffer_1100_pppd[1226];

    auto g_z_y_0_0_y_z_x_yy = buffer_1100_pppd[1227];

    auto g_z_y_0_0_y_z_x_yz = buffer_1100_pppd[1228];

    auto g_z_y_0_0_y_z_x_zz = buffer_1100_pppd[1229];

    auto g_z_y_0_0_y_z_y_xx = buffer_1100_pppd[1230];

    auto g_z_y_0_0_y_z_y_xy = buffer_1100_pppd[1231];

    auto g_z_y_0_0_y_z_y_xz = buffer_1100_pppd[1232];

    auto g_z_y_0_0_y_z_y_yy = buffer_1100_pppd[1233];

    auto g_z_y_0_0_y_z_y_yz = buffer_1100_pppd[1234];

    auto g_z_y_0_0_y_z_y_zz = buffer_1100_pppd[1235];

    auto g_z_y_0_0_y_z_z_xx = buffer_1100_pppd[1236];

    auto g_z_y_0_0_y_z_z_xy = buffer_1100_pppd[1237];

    auto g_z_y_0_0_y_z_z_xz = buffer_1100_pppd[1238];

    auto g_z_y_0_0_y_z_z_yy = buffer_1100_pppd[1239];

    auto g_z_y_0_0_y_z_z_yz = buffer_1100_pppd[1240];

    auto g_z_y_0_0_y_z_z_zz = buffer_1100_pppd[1241];

    auto g_z_y_0_0_z_x_x_xx = buffer_1100_pppd[1242];

    auto g_z_y_0_0_z_x_x_xy = buffer_1100_pppd[1243];

    auto g_z_y_0_0_z_x_x_xz = buffer_1100_pppd[1244];

    auto g_z_y_0_0_z_x_x_yy = buffer_1100_pppd[1245];

    auto g_z_y_0_0_z_x_x_yz = buffer_1100_pppd[1246];

    auto g_z_y_0_0_z_x_x_zz = buffer_1100_pppd[1247];

    auto g_z_y_0_0_z_x_y_xx = buffer_1100_pppd[1248];

    auto g_z_y_0_0_z_x_y_xy = buffer_1100_pppd[1249];

    auto g_z_y_0_0_z_x_y_xz = buffer_1100_pppd[1250];

    auto g_z_y_0_0_z_x_y_yy = buffer_1100_pppd[1251];

    auto g_z_y_0_0_z_x_y_yz = buffer_1100_pppd[1252];

    auto g_z_y_0_0_z_x_y_zz = buffer_1100_pppd[1253];

    auto g_z_y_0_0_z_x_z_xx = buffer_1100_pppd[1254];

    auto g_z_y_0_0_z_x_z_xy = buffer_1100_pppd[1255];

    auto g_z_y_0_0_z_x_z_xz = buffer_1100_pppd[1256];

    auto g_z_y_0_0_z_x_z_yy = buffer_1100_pppd[1257];

    auto g_z_y_0_0_z_x_z_yz = buffer_1100_pppd[1258];

    auto g_z_y_0_0_z_x_z_zz = buffer_1100_pppd[1259];

    auto g_z_y_0_0_z_y_x_xx = buffer_1100_pppd[1260];

    auto g_z_y_0_0_z_y_x_xy = buffer_1100_pppd[1261];

    auto g_z_y_0_0_z_y_x_xz = buffer_1100_pppd[1262];

    auto g_z_y_0_0_z_y_x_yy = buffer_1100_pppd[1263];

    auto g_z_y_0_0_z_y_x_yz = buffer_1100_pppd[1264];

    auto g_z_y_0_0_z_y_x_zz = buffer_1100_pppd[1265];

    auto g_z_y_0_0_z_y_y_xx = buffer_1100_pppd[1266];

    auto g_z_y_0_0_z_y_y_xy = buffer_1100_pppd[1267];

    auto g_z_y_0_0_z_y_y_xz = buffer_1100_pppd[1268];

    auto g_z_y_0_0_z_y_y_yy = buffer_1100_pppd[1269];

    auto g_z_y_0_0_z_y_y_yz = buffer_1100_pppd[1270];

    auto g_z_y_0_0_z_y_y_zz = buffer_1100_pppd[1271];

    auto g_z_y_0_0_z_y_z_xx = buffer_1100_pppd[1272];

    auto g_z_y_0_0_z_y_z_xy = buffer_1100_pppd[1273];

    auto g_z_y_0_0_z_y_z_xz = buffer_1100_pppd[1274];

    auto g_z_y_0_0_z_y_z_yy = buffer_1100_pppd[1275];

    auto g_z_y_0_0_z_y_z_yz = buffer_1100_pppd[1276];

    auto g_z_y_0_0_z_y_z_zz = buffer_1100_pppd[1277];

    auto g_z_y_0_0_z_z_x_xx = buffer_1100_pppd[1278];

    auto g_z_y_0_0_z_z_x_xy = buffer_1100_pppd[1279];

    auto g_z_y_0_0_z_z_x_xz = buffer_1100_pppd[1280];

    auto g_z_y_0_0_z_z_x_yy = buffer_1100_pppd[1281];

    auto g_z_y_0_0_z_z_x_yz = buffer_1100_pppd[1282];

    auto g_z_y_0_0_z_z_x_zz = buffer_1100_pppd[1283];

    auto g_z_y_0_0_z_z_y_xx = buffer_1100_pppd[1284];

    auto g_z_y_0_0_z_z_y_xy = buffer_1100_pppd[1285];

    auto g_z_y_0_0_z_z_y_xz = buffer_1100_pppd[1286];

    auto g_z_y_0_0_z_z_y_yy = buffer_1100_pppd[1287];

    auto g_z_y_0_0_z_z_y_yz = buffer_1100_pppd[1288];

    auto g_z_y_0_0_z_z_y_zz = buffer_1100_pppd[1289];

    auto g_z_y_0_0_z_z_z_xx = buffer_1100_pppd[1290];

    auto g_z_y_0_0_z_z_z_xy = buffer_1100_pppd[1291];

    auto g_z_y_0_0_z_z_z_xz = buffer_1100_pppd[1292];

    auto g_z_y_0_0_z_z_z_yy = buffer_1100_pppd[1293];

    auto g_z_y_0_0_z_z_z_yz = buffer_1100_pppd[1294];

    auto g_z_y_0_0_z_z_z_zz = buffer_1100_pppd[1295];

    auto g_z_z_0_0_x_x_x_xx = buffer_1100_pppd[1296];

    auto g_z_z_0_0_x_x_x_xy = buffer_1100_pppd[1297];

    auto g_z_z_0_0_x_x_x_xz = buffer_1100_pppd[1298];

    auto g_z_z_0_0_x_x_x_yy = buffer_1100_pppd[1299];

    auto g_z_z_0_0_x_x_x_yz = buffer_1100_pppd[1300];

    auto g_z_z_0_0_x_x_x_zz = buffer_1100_pppd[1301];

    auto g_z_z_0_0_x_x_y_xx = buffer_1100_pppd[1302];

    auto g_z_z_0_0_x_x_y_xy = buffer_1100_pppd[1303];

    auto g_z_z_0_0_x_x_y_xz = buffer_1100_pppd[1304];

    auto g_z_z_0_0_x_x_y_yy = buffer_1100_pppd[1305];

    auto g_z_z_0_0_x_x_y_yz = buffer_1100_pppd[1306];

    auto g_z_z_0_0_x_x_y_zz = buffer_1100_pppd[1307];

    auto g_z_z_0_0_x_x_z_xx = buffer_1100_pppd[1308];

    auto g_z_z_0_0_x_x_z_xy = buffer_1100_pppd[1309];

    auto g_z_z_0_0_x_x_z_xz = buffer_1100_pppd[1310];

    auto g_z_z_0_0_x_x_z_yy = buffer_1100_pppd[1311];

    auto g_z_z_0_0_x_x_z_yz = buffer_1100_pppd[1312];

    auto g_z_z_0_0_x_x_z_zz = buffer_1100_pppd[1313];

    auto g_z_z_0_0_x_y_x_xx = buffer_1100_pppd[1314];

    auto g_z_z_0_0_x_y_x_xy = buffer_1100_pppd[1315];

    auto g_z_z_0_0_x_y_x_xz = buffer_1100_pppd[1316];

    auto g_z_z_0_0_x_y_x_yy = buffer_1100_pppd[1317];

    auto g_z_z_0_0_x_y_x_yz = buffer_1100_pppd[1318];

    auto g_z_z_0_0_x_y_x_zz = buffer_1100_pppd[1319];

    auto g_z_z_0_0_x_y_y_xx = buffer_1100_pppd[1320];

    auto g_z_z_0_0_x_y_y_xy = buffer_1100_pppd[1321];

    auto g_z_z_0_0_x_y_y_xz = buffer_1100_pppd[1322];

    auto g_z_z_0_0_x_y_y_yy = buffer_1100_pppd[1323];

    auto g_z_z_0_0_x_y_y_yz = buffer_1100_pppd[1324];

    auto g_z_z_0_0_x_y_y_zz = buffer_1100_pppd[1325];

    auto g_z_z_0_0_x_y_z_xx = buffer_1100_pppd[1326];

    auto g_z_z_0_0_x_y_z_xy = buffer_1100_pppd[1327];

    auto g_z_z_0_0_x_y_z_xz = buffer_1100_pppd[1328];

    auto g_z_z_0_0_x_y_z_yy = buffer_1100_pppd[1329];

    auto g_z_z_0_0_x_y_z_yz = buffer_1100_pppd[1330];

    auto g_z_z_0_0_x_y_z_zz = buffer_1100_pppd[1331];

    auto g_z_z_0_0_x_z_x_xx = buffer_1100_pppd[1332];

    auto g_z_z_0_0_x_z_x_xy = buffer_1100_pppd[1333];

    auto g_z_z_0_0_x_z_x_xz = buffer_1100_pppd[1334];

    auto g_z_z_0_0_x_z_x_yy = buffer_1100_pppd[1335];

    auto g_z_z_0_0_x_z_x_yz = buffer_1100_pppd[1336];

    auto g_z_z_0_0_x_z_x_zz = buffer_1100_pppd[1337];

    auto g_z_z_0_0_x_z_y_xx = buffer_1100_pppd[1338];

    auto g_z_z_0_0_x_z_y_xy = buffer_1100_pppd[1339];

    auto g_z_z_0_0_x_z_y_xz = buffer_1100_pppd[1340];

    auto g_z_z_0_0_x_z_y_yy = buffer_1100_pppd[1341];

    auto g_z_z_0_0_x_z_y_yz = buffer_1100_pppd[1342];

    auto g_z_z_0_0_x_z_y_zz = buffer_1100_pppd[1343];

    auto g_z_z_0_0_x_z_z_xx = buffer_1100_pppd[1344];

    auto g_z_z_0_0_x_z_z_xy = buffer_1100_pppd[1345];

    auto g_z_z_0_0_x_z_z_xz = buffer_1100_pppd[1346];

    auto g_z_z_0_0_x_z_z_yy = buffer_1100_pppd[1347];

    auto g_z_z_0_0_x_z_z_yz = buffer_1100_pppd[1348];

    auto g_z_z_0_0_x_z_z_zz = buffer_1100_pppd[1349];

    auto g_z_z_0_0_y_x_x_xx = buffer_1100_pppd[1350];

    auto g_z_z_0_0_y_x_x_xy = buffer_1100_pppd[1351];

    auto g_z_z_0_0_y_x_x_xz = buffer_1100_pppd[1352];

    auto g_z_z_0_0_y_x_x_yy = buffer_1100_pppd[1353];

    auto g_z_z_0_0_y_x_x_yz = buffer_1100_pppd[1354];

    auto g_z_z_0_0_y_x_x_zz = buffer_1100_pppd[1355];

    auto g_z_z_0_0_y_x_y_xx = buffer_1100_pppd[1356];

    auto g_z_z_0_0_y_x_y_xy = buffer_1100_pppd[1357];

    auto g_z_z_0_0_y_x_y_xz = buffer_1100_pppd[1358];

    auto g_z_z_0_0_y_x_y_yy = buffer_1100_pppd[1359];

    auto g_z_z_0_0_y_x_y_yz = buffer_1100_pppd[1360];

    auto g_z_z_0_0_y_x_y_zz = buffer_1100_pppd[1361];

    auto g_z_z_0_0_y_x_z_xx = buffer_1100_pppd[1362];

    auto g_z_z_0_0_y_x_z_xy = buffer_1100_pppd[1363];

    auto g_z_z_0_0_y_x_z_xz = buffer_1100_pppd[1364];

    auto g_z_z_0_0_y_x_z_yy = buffer_1100_pppd[1365];

    auto g_z_z_0_0_y_x_z_yz = buffer_1100_pppd[1366];

    auto g_z_z_0_0_y_x_z_zz = buffer_1100_pppd[1367];

    auto g_z_z_0_0_y_y_x_xx = buffer_1100_pppd[1368];

    auto g_z_z_0_0_y_y_x_xy = buffer_1100_pppd[1369];

    auto g_z_z_0_0_y_y_x_xz = buffer_1100_pppd[1370];

    auto g_z_z_0_0_y_y_x_yy = buffer_1100_pppd[1371];

    auto g_z_z_0_0_y_y_x_yz = buffer_1100_pppd[1372];

    auto g_z_z_0_0_y_y_x_zz = buffer_1100_pppd[1373];

    auto g_z_z_0_0_y_y_y_xx = buffer_1100_pppd[1374];

    auto g_z_z_0_0_y_y_y_xy = buffer_1100_pppd[1375];

    auto g_z_z_0_0_y_y_y_xz = buffer_1100_pppd[1376];

    auto g_z_z_0_0_y_y_y_yy = buffer_1100_pppd[1377];

    auto g_z_z_0_0_y_y_y_yz = buffer_1100_pppd[1378];

    auto g_z_z_0_0_y_y_y_zz = buffer_1100_pppd[1379];

    auto g_z_z_0_0_y_y_z_xx = buffer_1100_pppd[1380];

    auto g_z_z_0_0_y_y_z_xy = buffer_1100_pppd[1381];

    auto g_z_z_0_0_y_y_z_xz = buffer_1100_pppd[1382];

    auto g_z_z_0_0_y_y_z_yy = buffer_1100_pppd[1383];

    auto g_z_z_0_0_y_y_z_yz = buffer_1100_pppd[1384];

    auto g_z_z_0_0_y_y_z_zz = buffer_1100_pppd[1385];

    auto g_z_z_0_0_y_z_x_xx = buffer_1100_pppd[1386];

    auto g_z_z_0_0_y_z_x_xy = buffer_1100_pppd[1387];

    auto g_z_z_0_0_y_z_x_xz = buffer_1100_pppd[1388];

    auto g_z_z_0_0_y_z_x_yy = buffer_1100_pppd[1389];

    auto g_z_z_0_0_y_z_x_yz = buffer_1100_pppd[1390];

    auto g_z_z_0_0_y_z_x_zz = buffer_1100_pppd[1391];

    auto g_z_z_0_0_y_z_y_xx = buffer_1100_pppd[1392];

    auto g_z_z_0_0_y_z_y_xy = buffer_1100_pppd[1393];

    auto g_z_z_0_0_y_z_y_xz = buffer_1100_pppd[1394];

    auto g_z_z_0_0_y_z_y_yy = buffer_1100_pppd[1395];

    auto g_z_z_0_0_y_z_y_yz = buffer_1100_pppd[1396];

    auto g_z_z_0_0_y_z_y_zz = buffer_1100_pppd[1397];

    auto g_z_z_0_0_y_z_z_xx = buffer_1100_pppd[1398];

    auto g_z_z_0_0_y_z_z_xy = buffer_1100_pppd[1399];

    auto g_z_z_0_0_y_z_z_xz = buffer_1100_pppd[1400];

    auto g_z_z_0_0_y_z_z_yy = buffer_1100_pppd[1401];

    auto g_z_z_0_0_y_z_z_yz = buffer_1100_pppd[1402];

    auto g_z_z_0_0_y_z_z_zz = buffer_1100_pppd[1403];

    auto g_z_z_0_0_z_x_x_xx = buffer_1100_pppd[1404];

    auto g_z_z_0_0_z_x_x_xy = buffer_1100_pppd[1405];

    auto g_z_z_0_0_z_x_x_xz = buffer_1100_pppd[1406];

    auto g_z_z_0_0_z_x_x_yy = buffer_1100_pppd[1407];

    auto g_z_z_0_0_z_x_x_yz = buffer_1100_pppd[1408];

    auto g_z_z_0_0_z_x_x_zz = buffer_1100_pppd[1409];

    auto g_z_z_0_0_z_x_y_xx = buffer_1100_pppd[1410];

    auto g_z_z_0_0_z_x_y_xy = buffer_1100_pppd[1411];

    auto g_z_z_0_0_z_x_y_xz = buffer_1100_pppd[1412];

    auto g_z_z_0_0_z_x_y_yy = buffer_1100_pppd[1413];

    auto g_z_z_0_0_z_x_y_yz = buffer_1100_pppd[1414];

    auto g_z_z_0_0_z_x_y_zz = buffer_1100_pppd[1415];

    auto g_z_z_0_0_z_x_z_xx = buffer_1100_pppd[1416];

    auto g_z_z_0_0_z_x_z_xy = buffer_1100_pppd[1417];

    auto g_z_z_0_0_z_x_z_xz = buffer_1100_pppd[1418];

    auto g_z_z_0_0_z_x_z_yy = buffer_1100_pppd[1419];

    auto g_z_z_0_0_z_x_z_yz = buffer_1100_pppd[1420];

    auto g_z_z_0_0_z_x_z_zz = buffer_1100_pppd[1421];

    auto g_z_z_0_0_z_y_x_xx = buffer_1100_pppd[1422];

    auto g_z_z_0_0_z_y_x_xy = buffer_1100_pppd[1423];

    auto g_z_z_0_0_z_y_x_xz = buffer_1100_pppd[1424];

    auto g_z_z_0_0_z_y_x_yy = buffer_1100_pppd[1425];

    auto g_z_z_0_0_z_y_x_yz = buffer_1100_pppd[1426];

    auto g_z_z_0_0_z_y_x_zz = buffer_1100_pppd[1427];

    auto g_z_z_0_0_z_y_y_xx = buffer_1100_pppd[1428];

    auto g_z_z_0_0_z_y_y_xy = buffer_1100_pppd[1429];

    auto g_z_z_0_0_z_y_y_xz = buffer_1100_pppd[1430];

    auto g_z_z_0_0_z_y_y_yy = buffer_1100_pppd[1431];

    auto g_z_z_0_0_z_y_y_yz = buffer_1100_pppd[1432];

    auto g_z_z_0_0_z_y_y_zz = buffer_1100_pppd[1433];

    auto g_z_z_0_0_z_y_z_xx = buffer_1100_pppd[1434];

    auto g_z_z_0_0_z_y_z_xy = buffer_1100_pppd[1435];

    auto g_z_z_0_0_z_y_z_xz = buffer_1100_pppd[1436];

    auto g_z_z_0_0_z_y_z_yy = buffer_1100_pppd[1437];

    auto g_z_z_0_0_z_y_z_yz = buffer_1100_pppd[1438];

    auto g_z_z_0_0_z_y_z_zz = buffer_1100_pppd[1439];

    auto g_z_z_0_0_z_z_x_xx = buffer_1100_pppd[1440];

    auto g_z_z_0_0_z_z_x_xy = buffer_1100_pppd[1441];

    auto g_z_z_0_0_z_z_x_xz = buffer_1100_pppd[1442];

    auto g_z_z_0_0_z_z_x_yy = buffer_1100_pppd[1443];

    auto g_z_z_0_0_z_z_x_yz = buffer_1100_pppd[1444];

    auto g_z_z_0_0_z_z_x_zz = buffer_1100_pppd[1445];

    auto g_z_z_0_0_z_z_y_xx = buffer_1100_pppd[1446];

    auto g_z_z_0_0_z_z_y_xy = buffer_1100_pppd[1447];

    auto g_z_z_0_0_z_z_y_xz = buffer_1100_pppd[1448];

    auto g_z_z_0_0_z_z_y_yy = buffer_1100_pppd[1449];

    auto g_z_z_0_0_z_z_y_yz = buffer_1100_pppd[1450];

    auto g_z_z_0_0_z_z_y_zz = buffer_1100_pppd[1451];

    auto g_z_z_0_0_z_z_z_xx = buffer_1100_pppd[1452];

    auto g_z_z_0_0_z_z_z_xy = buffer_1100_pppd[1453];

    auto g_z_z_0_0_z_z_z_xz = buffer_1100_pppd[1454];

    auto g_z_z_0_0_z_z_z_yy = buffer_1100_pppd[1455];

    auto g_z_z_0_0_z_z_z_yz = buffer_1100_pppd[1456];

    auto g_z_z_0_0_z_z_z_zz = buffer_1100_pppd[1457];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_x_x_0_0_x_x_x_xx, g_x_x_0_0_x_x_x_xy, g_x_x_0_0_x_x_x_xz, g_x_x_0_0_x_x_x_yy, g_x_x_0_0_x_x_x_yz, g_x_x_0_0_x_x_x_zz, g_xx_0_x_xx, g_xx_0_x_xy, g_xx_0_x_xz, g_xx_0_x_yy, g_xx_0_x_yz, g_xx_0_x_zz, g_xx_xx_x_xx, g_xx_xx_x_xy, g_xx_xx_x_xz, g_xx_xx_x_yy, g_xx_xx_x_yz, g_xx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_xx_x_xx[i] * b_exp - 2.0 * g_xx_0_x_xx[i] * a_exp + 4.0 * g_xx_xx_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_xx_x_xy[i] * b_exp - 2.0 * g_xx_0_x_xy[i] * a_exp + 4.0 * g_xx_xx_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_xx_x_xz[i] * b_exp - 2.0 * g_xx_0_x_xz[i] * a_exp + 4.0 * g_xx_xx_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_xx_x_yy[i] * b_exp - 2.0 * g_xx_0_x_yy[i] * a_exp + 4.0 * g_xx_xx_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_xx_x_yz[i] * b_exp - 2.0 * g_xx_0_x_yz[i] * a_exp + 4.0 * g_xx_xx_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_xx_x_zz[i] * b_exp - 2.0 * g_xx_0_x_zz[i] * a_exp + 4.0 * g_xx_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_x_x_0_0_x_x_y_xx, g_x_x_0_0_x_x_y_xy, g_x_x_0_0_x_x_y_xz, g_x_x_0_0_x_x_y_yy, g_x_x_0_0_x_x_y_yz, g_x_x_0_0_x_x_y_zz, g_xx_0_y_xx, g_xx_0_y_xy, g_xx_0_y_xz, g_xx_0_y_yy, g_xx_0_y_yz, g_xx_0_y_zz, g_xx_xx_y_xx, g_xx_xx_y_xy, g_xx_xx_y_xz, g_xx_xx_y_yy, g_xx_xx_y_yz, g_xx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_xx_y_xx[i] * b_exp - 2.0 * g_xx_0_y_xx[i] * a_exp + 4.0 * g_xx_xx_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_xx_y_xy[i] * b_exp - 2.0 * g_xx_0_y_xy[i] * a_exp + 4.0 * g_xx_xx_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_xx_y_xz[i] * b_exp - 2.0 * g_xx_0_y_xz[i] * a_exp + 4.0 * g_xx_xx_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_xx_y_yy[i] * b_exp - 2.0 * g_xx_0_y_yy[i] * a_exp + 4.0 * g_xx_xx_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_xx_y_yz[i] * b_exp - 2.0 * g_xx_0_y_yz[i] * a_exp + 4.0 * g_xx_xx_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_xx_y_zz[i] * b_exp - 2.0 * g_xx_0_y_zz[i] * a_exp + 4.0 * g_xx_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_x_x_0_0_x_x_z_xx, g_x_x_0_0_x_x_z_xy, g_x_x_0_0_x_x_z_xz, g_x_x_0_0_x_x_z_yy, g_x_x_0_0_x_x_z_yz, g_x_x_0_0_x_x_z_zz, g_xx_0_z_xx, g_xx_0_z_xy, g_xx_0_z_xz, g_xx_0_z_yy, g_xx_0_z_yz, g_xx_0_z_zz, g_xx_xx_z_xx, g_xx_xx_z_xy, g_xx_xx_z_xz, g_xx_xx_z_yy, g_xx_xx_z_yz, g_xx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_xx_z_xx[i] * b_exp - 2.0 * g_xx_0_z_xx[i] * a_exp + 4.0 * g_xx_xx_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_xx_z_xy[i] * b_exp - 2.0 * g_xx_0_z_xy[i] * a_exp + 4.0 * g_xx_xx_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_xx_z_xz[i] * b_exp - 2.0 * g_xx_0_z_xz[i] * a_exp + 4.0 * g_xx_xx_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_xx_z_yy[i] * b_exp - 2.0 * g_xx_0_z_yy[i] * a_exp + 4.0 * g_xx_xx_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_xx_z_yz[i] * b_exp - 2.0 * g_xx_0_z_yz[i] * a_exp + 4.0 * g_xx_xx_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_xx_z_zz[i] * b_exp - 2.0 * g_xx_0_z_zz[i] * a_exp + 4.0 * g_xx_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_x_x_0_0_x_y_x_xx, g_x_x_0_0_x_y_x_xy, g_x_x_0_0_x_y_x_xz, g_x_x_0_0_x_y_x_yy, g_x_x_0_0_x_y_x_yz, g_x_x_0_0_x_y_x_zz, g_xx_xy_x_xx, g_xx_xy_x_xy, g_xx_xy_x_xz, g_xx_xy_x_yy, g_xx_xy_x_yz, g_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_xx_xy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_xx_xy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_xx_xy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_xx_xy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_xx_xy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_xx_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_x_x_0_0_x_y_y_xx, g_x_x_0_0_x_y_y_xy, g_x_x_0_0_x_y_y_xz, g_x_x_0_0_x_y_y_yy, g_x_x_0_0_x_y_y_yz, g_x_x_0_0_x_y_y_zz, g_xx_xy_y_xx, g_xx_xy_y_xy, g_xx_xy_y_xz, g_xx_xy_y_yy, g_xx_xy_y_yz, g_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_xx_xy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_xx_xy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_xx_xy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_xx_xy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_xx_xy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_xx_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_x_x_0_0_x_y_z_xx, g_x_x_0_0_x_y_z_xy, g_x_x_0_0_x_y_z_xz, g_x_x_0_0_x_y_z_yy, g_x_x_0_0_x_y_z_yz, g_x_x_0_0_x_y_z_zz, g_xx_xy_z_xx, g_xx_xy_z_xy, g_xx_xy_z_xz, g_xx_xy_z_yy, g_xx_xy_z_yz, g_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_xx_xy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_xx_xy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_xx_xy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_xx_xy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_xx_xy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_xx_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_x_x_0_0_x_z_x_xx, g_x_x_0_0_x_z_x_xy, g_x_x_0_0_x_z_x_xz, g_x_x_0_0_x_z_x_yy, g_x_x_0_0_x_z_x_yz, g_x_x_0_0_x_z_x_zz, g_xx_xz_x_xx, g_xx_xz_x_xy, g_xx_xz_x_xz, g_xx_xz_x_yy, g_xx_xz_x_yz, g_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_xx_xz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_xx_xz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_xx_xz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_xx_xz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_xx_xz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_xx_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_x_x_0_0_x_z_y_xx, g_x_x_0_0_x_z_y_xy, g_x_x_0_0_x_z_y_xz, g_x_x_0_0_x_z_y_yy, g_x_x_0_0_x_z_y_yz, g_x_x_0_0_x_z_y_zz, g_xx_xz_y_xx, g_xx_xz_y_xy, g_xx_xz_y_xz, g_xx_xz_y_yy, g_xx_xz_y_yz, g_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_xx_xz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_xx_xz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_xx_xz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_xx_xz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_xx_xz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_xx_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_x_x_0_0_x_z_z_xx, g_x_x_0_0_x_z_z_xy, g_x_x_0_0_x_z_z_xz, g_x_x_0_0_x_z_z_yy, g_x_x_0_0_x_z_z_yz, g_x_x_0_0_x_z_z_zz, g_xx_xz_z_xx, g_xx_xz_z_xy, g_xx_xz_z_xz, g_xx_xz_z_yy, g_xx_xz_z_yz, g_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_xx_xz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_xx_xz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_xx_xz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_xx_xz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_xx_xz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_xx_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_y_x_x_xx, g_x_x_0_0_y_x_x_xy, g_x_x_0_0_y_x_x_xz, g_x_x_0_0_y_x_x_yy, g_x_x_0_0_y_x_x_yz, g_x_x_0_0_y_x_x_zz, g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_xx_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_xx_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_xx_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_xx_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_xx_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_y_x_y_xx, g_x_x_0_0_y_x_y_xy, g_x_x_0_0_y_x_y_xz, g_x_x_0_0_y_x_y_yy, g_x_x_0_0_y_x_y_yz, g_x_x_0_0_y_x_y_zz, g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_xx_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_xx_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_xx_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_xx_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_xx_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_y_x_z_xx, g_x_x_0_0_y_x_z_xy, g_x_x_0_0_y_x_z_xz, g_x_x_0_0_y_x_z_yy, g_x_x_0_0_y_x_z_yz, g_x_x_0_0_y_x_z_zz, g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_xx_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_xx_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_xx_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_xx_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_xx_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_y_y_x_xx, g_x_x_0_0_y_y_x_xy, g_x_x_0_0_y_y_x_xz, g_x_x_0_0_y_y_x_yy, g_x_x_0_0_y_y_x_yz, g_x_x_0_0_y_y_x_zz, g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_x_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_y_y_y_xx, g_x_x_0_0_y_y_y_xy, g_x_x_0_0_y_y_y_xz, g_x_x_0_0_y_y_y_yy, g_x_x_0_0_y_y_y_yz, g_x_x_0_0_y_y_y_zz, g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_y_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_y_y_z_xx, g_x_x_0_0_y_y_z_xy, g_x_x_0_0_y_y_z_xz, g_x_x_0_0_y_y_z_yy, g_x_x_0_0_y_y_z_yz, g_x_x_0_0_y_y_z_zz, g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_z_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_y_z_x_xx, g_x_x_0_0_y_z_x_xy, g_x_x_0_0_y_z_x_xz, g_x_x_0_0_y_z_x_yy, g_x_x_0_0_y_z_x_yz, g_x_x_0_0_y_z_x_zz, g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_x_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_y_z_y_xx, g_x_x_0_0_y_z_y_xy, g_x_x_0_0_y_z_y_xz, g_x_x_0_0_y_z_y_yy, g_x_x_0_0_y_z_y_yz, g_x_x_0_0_y_z_y_zz, g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_y_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_y_z_z_xx, g_x_x_0_0_y_z_z_xy, g_x_x_0_0_y_z_z_xz, g_x_x_0_0_y_z_z_yy, g_x_x_0_0_y_z_z_yz, g_x_x_0_0_y_z_z_zz, g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_z_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_x_0_0_z_x_x_xx, g_x_x_0_0_z_x_x_xy, g_x_x_0_0_z_x_x_xz, g_x_x_0_0_z_x_x_yy, g_x_x_0_0_z_x_x_yz, g_x_x_0_0_z_x_x_zz, g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_xx_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_xx_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_xx_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_xx_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_xx_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_x_0_0_z_x_y_xx, g_x_x_0_0_z_x_y_xy, g_x_x_0_0_z_x_y_xz, g_x_x_0_0_z_x_y_yy, g_x_x_0_0_z_x_y_yz, g_x_x_0_0_z_x_y_zz, g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_xx_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_xx_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_xx_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_xx_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_xx_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_x_0_0_z_x_z_xx, g_x_x_0_0_z_x_z_xy, g_x_x_0_0_z_x_z_xz, g_x_x_0_0_z_x_z_yy, g_x_x_0_0_z_x_z_yz, g_x_x_0_0_z_x_z_zz, g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_xx_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_xx_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_xx_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_xx_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_xx_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_x_0_0_z_y_x_xx, g_x_x_0_0_z_y_x_xy, g_x_x_0_0_z_y_x_xz, g_x_x_0_0_z_y_x_yy, g_x_x_0_0_z_y_x_yz, g_x_x_0_0_z_y_x_zz, g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_x_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_x_0_0_z_y_y_xx, g_x_x_0_0_z_y_y_xy, g_x_x_0_0_z_y_y_xz, g_x_x_0_0_z_y_y_yy, g_x_x_0_0_z_y_y_yz, g_x_x_0_0_z_y_y_zz, g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_y_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_x_0_0_z_y_z_xx, g_x_x_0_0_z_y_z_xy, g_x_x_0_0_z_y_z_xz, g_x_x_0_0_z_y_z_yy, g_x_x_0_0_z_y_z_yz, g_x_x_0_0_z_y_z_zz, g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_z_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_x_0_0_z_z_x_xx, g_x_x_0_0_z_z_x_xy, g_x_x_0_0_z_z_x_xz, g_x_x_0_0_z_z_x_yy, g_x_x_0_0_z_z_x_yz, g_x_x_0_0_z_z_x_zz, g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_x_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_x_0_0_z_z_y_xx, g_x_x_0_0_z_z_y_xy, g_x_x_0_0_z_z_y_xz, g_x_x_0_0_z_z_y_yy, g_x_x_0_0_z_z_y_yz, g_x_x_0_0_z_z_y_zz, g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_y_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_x_0_0_z_z_z_xx, g_x_x_0_0_z_z_z_xy, g_x_x_0_0_z_z_z_xz, g_x_x_0_0_z_z_z_yy, g_x_x_0_0_z_z_z_yz, g_x_x_0_0_z_z_z_zz, g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_z_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_x_y_0_0_x_x_x_xx, g_x_y_0_0_x_x_x_xy, g_x_y_0_0_x_x_x_xz, g_x_y_0_0_x_x_x_yy, g_x_y_0_0_x_x_x_yz, g_x_y_0_0_x_x_x_zz, g_xx_xy_x_xx, g_xx_xy_x_xy, g_xx_xy_x_xz, g_xx_xy_x_yy, g_xx_xy_x_yz, g_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_xx_xy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_xx_xy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_xx_xy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_xx_xy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_xx_xy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_xx_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_x_y_0_0_x_x_y_xx, g_x_y_0_0_x_x_y_xy, g_x_y_0_0_x_x_y_xz, g_x_y_0_0_x_x_y_yy, g_x_y_0_0_x_x_y_yz, g_x_y_0_0_x_x_y_zz, g_xx_xy_y_xx, g_xx_xy_y_xy, g_xx_xy_y_xz, g_xx_xy_y_yy, g_xx_xy_y_yz, g_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_xx_xy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_xx_xy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_xx_xy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_xx_xy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_xx_xy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_xx_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_x_y_0_0_x_x_z_xx, g_x_y_0_0_x_x_z_xy, g_x_y_0_0_x_x_z_xz, g_x_y_0_0_x_x_z_yy, g_x_y_0_0_x_x_z_yz, g_x_y_0_0_x_x_z_zz, g_xx_xy_z_xx, g_xx_xy_z_xy, g_xx_xy_z_xz, g_xx_xy_z_yy, g_xx_xy_z_yz, g_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_xx_xy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_xx_xy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_xx_xy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_xx_xy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_xx_xy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_xx_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_x_y_0_0_x_y_x_xx, g_x_y_0_0_x_y_x_xy, g_x_y_0_0_x_y_x_xz, g_x_y_0_0_x_y_x_yy, g_x_y_0_0_x_y_x_yz, g_x_y_0_0_x_y_x_zz, g_xx_0_x_xx, g_xx_0_x_xy, g_xx_0_x_xz, g_xx_0_x_yy, g_xx_0_x_yz, g_xx_0_x_zz, g_xx_yy_x_xx, g_xx_yy_x_xy, g_xx_yy_x_xz, g_xx_yy_x_yy, g_xx_yy_x_yz, g_xx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_yy_x_xx[i] * b_exp - 2.0 * g_xx_0_x_xx[i] * a_exp + 4.0 * g_xx_yy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_yy_x_xy[i] * b_exp - 2.0 * g_xx_0_x_xy[i] * a_exp + 4.0 * g_xx_yy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_yy_x_xz[i] * b_exp - 2.0 * g_xx_0_x_xz[i] * a_exp + 4.0 * g_xx_yy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_yy_x_yy[i] * b_exp - 2.0 * g_xx_0_x_yy[i] * a_exp + 4.0 * g_xx_yy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_yy_x_yz[i] * b_exp - 2.0 * g_xx_0_x_yz[i] * a_exp + 4.0 * g_xx_yy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_yy_x_zz[i] * b_exp - 2.0 * g_xx_0_x_zz[i] * a_exp + 4.0 * g_xx_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_x_y_0_0_x_y_y_xx, g_x_y_0_0_x_y_y_xy, g_x_y_0_0_x_y_y_xz, g_x_y_0_0_x_y_y_yy, g_x_y_0_0_x_y_y_yz, g_x_y_0_0_x_y_y_zz, g_xx_0_y_xx, g_xx_0_y_xy, g_xx_0_y_xz, g_xx_0_y_yy, g_xx_0_y_yz, g_xx_0_y_zz, g_xx_yy_y_xx, g_xx_yy_y_xy, g_xx_yy_y_xz, g_xx_yy_y_yy, g_xx_yy_y_yz, g_xx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_yy_y_xx[i] * b_exp - 2.0 * g_xx_0_y_xx[i] * a_exp + 4.0 * g_xx_yy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_yy_y_xy[i] * b_exp - 2.0 * g_xx_0_y_xy[i] * a_exp + 4.0 * g_xx_yy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_yy_y_xz[i] * b_exp - 2.0 * g_xx_0_y_xz[i] * a_exp + 4.0 * g_xx_yy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_yy_y_yy[i] * b_exp - 2.0 * g_xx_0_y_yy[i] * a_exp + 4.0 * g_xx_yy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_yy_y_yz[i] * b_exp - 2.0 * g_xx_0_y_yz[i] * a_exp + 4.0 * g_xx_yy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_yy_y_zz[i] * b_exp - 2.0 * g_xx_0_y_zz[i] * a_exp + 4.0 * g_xx_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_x_y_0_0_x_y_z_xx, g_x_y_0_0_x_y_z_xy, g_x_y_0_0_x_y_z_xz, g_x_y_0_0_x_y_z_yy, g_x_y_0_0_x_y_z_yz, g_x_y_0_0_x_y_z_zz, g_xx_0_z_xx, g_xx_0_z_xy, g_xx_0_z_xz, g_xx_0_z_yy, g_xx_0_z_yz, g_xx_0_z_zz, g_xx_yy_z_xx, g_xx_yy_z_xy, g_xx_yy_z_xz, g_xx_yy_z_yy, g_xx_yy_z_yz, g_xx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_yy_z_xx[i] * b_exp - 2.0 * g_xx_0_z_xx[i] * a_exp + 4.0 * g_xx_yy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_yy_z_xy[i] * b_exp - 2.0 * g_xx_0_z_xy[i] * a_exp + 4.0 * g_xx_yy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_yy_z_xz[i] * b_exp - 2.0 * g_xx_0_z_xz[i] * a_exp + 4.0 * g_xx_yy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_yy_z_yy[i] * b_exp - 2.0 * g_xx_0_z_yy[i] * a_exp + 4.0 * g_xx_yy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_yy_z_yz[i] * b_exp - 2.0 * g_xx_0_z_yz[i] * a_exp + 4.0 * g_xx_yy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_yy_z_zz[i] * b_exp - 2.0 * g_xx_0_z_zz[i] * a_exp + 4.0 * g_xx_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_x_y_0_0_x_z_x_xx, g_x_y_0_0_x_z_x_xy, g_x_y_0_0_x_z_x_xz, g_x_y_0_0_x_z_x_yy, g_x_y_0_0_x_z_x_yz, g_x_y_0_0_x_z_x_zz, g_xx_yz_x_xx, g_xx_yz_x_xy, g_xx_yz_x_xz, g_xx_yz_x_yy, g_xx_yz_x_yz, g_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_xx_yz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_xx_yz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_xx_yz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_xx_yz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_xx_yz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_xx_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_x_y_0_0_x_z_y_xx, g_x_y_0_0_x_z_y_xy, g_x_y_0_0_x_z_y_xz, g_x_y_0_0_x_z_y_yy, g_x_y_0_0_x_z_y_yz, g_x_y_0_0_x_z_y_zz, g_xx_yz_y_xx, g_xx_yz_y_xy, g_xx_yz_y_xz, g_xx_yz_y_yy, g_xx_yz_y_yz, g_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_xx_yz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_xx_yz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_xx_yz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_xx_yz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_xx_yz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_xx_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_x_y_0_0_x_z_z_xx, g_x_y_0_0_x_z_z_xy, g_x_y_0_0_x_z_z_xz, g_x_y_0_0_x_z_z_yy, g_x_y_0_0_x_z_z_yz, g_x_y_0_0_x_z_z_zz, g_xx_yz_z_xx, g_xx_yz_z_xy, g_xx_yz_z_xz, g_xx_yz_z_yy, g_xx_yz_z_yz, g_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_xx_yz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_xx_yz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_xx_yz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_xx_yz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_xx_yz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_xx_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_y_0_0_y_x_x_xx, g_x_y_0_0_y_x_x_xy, g_x_y_0_0_y_x_x_xz, g_x_y_0_0_y_x_x_yy, g_x_y_0_0_y_x_x_yz, g_x_y_0_0_y_x_x_zz, g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_x_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_y_0_0_y_x_y_xx, g_x_y_0_0_y_x_y_xy, g_x_y_0_0_y_x_y_xz, g_x_y_0_0_y_x_y_yy, g_x_y_0_0_y_x_y_yz, g_x_y_0_0_y_x_y_zz, g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_y_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_y_0_0_y_x_z_xx, g_x_y_0_0_y_x_z_xy, g_x_y_0_0_y_x_z_xz, g_x_y_0_0_y_x_z_yy, g_x_y_0_0_y_x_z_yz, g_x_y_0_0_y_x_z_zz, g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_z_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_y_0_0_y_y_x_xx, g_x_y_0_0_y_y_x_xy, g_x_y_0_0_y_y_x_xz, g_x_y_0_0_y_y_x_yy, g_x_y_0_0_y_y_x_yz, g_x_y_0_0_y_y_x_zz, g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_yy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_yy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_yy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_yy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_yy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_y_0_0_y_y_y_xx, g_x_y_0_0_y_y_y_xy, g_x_y_0_0_y_y_y_xz, g_x_y_0_0_y_y_y_yy, g_x_y_0_0_y_y_y_yz, g_x_y_0_0_y_y_y_zz, g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_yy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_yy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_yy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_yy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_yy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_y_0_0_y_y_z_xx, g_x_y_0_0_y_y_z_xy, g_x_y_0_0_y_y_z_xz, g_x_y_0_0_y_y_z_yy, g_x_y_0_0_y_y_z_yz, g_x_y_0_0_y_y_z_zz, g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_yy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_yy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_yy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_yy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_yy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_y_0_0_y_z_x_xx, g_x_y_0_0_y_z_x_xy, g_x_y_0_0_y_z_x_xz, g_x_y_0_0_y_z_x_yy, g_x_y_0_0_y_z_x_yz, g_x_y_0_0_y_z_x_zz, g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_x_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_y_0_0_y_z_y_xx, g_x_y_0_0_y_z_y_xy, g_x_y_0_0_y_z_y_xz, g_x_y_0_0_y_z_y_yy, g_x_y_0_0_y_z_y_yz, g_x_y_0_0_y_z_y_zz, g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_y_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_y_0_0_y_z_z_xx, g_x_y_0_0_y_z_z_xy, g_x_y_0_0_y_z_z_xz, g_x_y_0_0_y_z_z_yy, g_x_y_0_0_y_z_z_yz, g_x_y_0_0_y_z_z_zz, g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_z_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_y_0_0_z_x_x_xx, g_x_y_0_0_z_x_x_xy, g_x_y_0_0_z_x_x_xz, g_x_y_0_0_z_x_x_yy, g_x_y_0_0_z_x_x_yz, g_x_y_0_0_z_x_x_zz, g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_x_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_y_0_0_z_x_y_xx, g_x_y_0_0_z_x_y_xy, g_x_y_0_0_z_x_y_xz, g_x_y_0_0_z_x_y_yy, g_x_y_0_0_z_x_y_yz, g_x_y_0_0_z_x_y_zz, g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_y_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_y_0_0_z_x_z_xx, g_x_y_0_0_z_x_z_xy, g_x_y_0_0_z_x_z_xz, g_x_y_0_0_z_x_z_yy, g_x_y_0_0_z_x_z_yz, g_x_y_0_0_z_x_z_zz, g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_z_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_y_0_0_z_y_x_xx, g_x_y_0_0_z_y_x_xy, g_x_y_0_0_z_y_x_xz, g_x_y_0_0_z_y_x_yy, g_x_y_0_0_z_y_x_yz, g_x_y_0_0_z_y_x_zz, g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_yy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_yy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_yy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_yy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_yy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_y_0_0_z_y_y_xx, g_x_y_0_0_z_y_y_xy, g_x_y_0_0_z_y_y_xz, g_x_y_0_0_z_y_y_yy, g_x_y_0_0_z_y_y_yz, g_x_y_0_0_z_y_y_zz, g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_yy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_yy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_yy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_yy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_yy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_y_0_0_z_y_z_xx, g_x_y_0_0_z_y_z_xy, g_x_y_0_0_z_y_z_xz, g_x_y_0_0_z_y_z_yy, g_x_y_0_0_z_y_z_yz, g_x_y_0_0_z_y_z_zz, g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_yy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_yy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_yy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_yy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_yy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_y_0_0_z_z_x_xx, g_x_y_0_0_z_z_x_xy, g_x_y_0_0_z_z_x_xz, g_x_y_0_0_z_z_x_yy, g_x_y_0_0_z_z_x_yz, g_x_y_0_0_z_z_x_zz, g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_x_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_y_0_0_z_z_y_xx, g_x_y_0_0_z_z_y_xy, g_x_y_0_0_z_z_y_xz, g_x_y_0_0_z_z_y_yy, g_x_y_0_0_z_z_y_yz, g_x_y_0_0_z_z_y_zz, g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_y_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_y_0_0_z_z_z_xx, g_x_y_0_0_z_z_z_xy, g_x_y_0_0_z_z_z_xz, g_x_y_0_0_z_z_z_yy, g_x_y_0_0_z_z_z_yz, g_x_y_0_0_z_z_z_zz, g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_z_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_x_z_0_0_x_x_x_xx, g_x_z_0_0_x_x_x_xy, g_x_z_0_0_x_x_x_xz, g_x_z_0_0_x_x_x_yy, g_x_z_0_0_x_x_x_yz, g_x_z_0_0_x_x_x_zz, g_xx_xz_x_xx, g_xx_xz_x_xy, g_xx_xz_x_xz, g_xx_xz_x_yy, g_xx_xz_x_yz, g_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_xx_xz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_xx_xz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_xx_xz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_xx_xz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_xx_xz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_xx_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_x_z_0_0_x_x_y_xx, g_x_z_0_0_x_x_y_xy, g_x_z_0_0_x_x_y_xz, g_x_z_0_0_x_x_y_yy, g_x_z_0_0_x_x_y_yz, g_x_z_0_0_x_x_y_zz, g_xx_xz_y_xx, g_xx_xz_y_xy, g_xx_xz_y_xz, g_xx_xz_y_yy, g_xx_xz_y_yz, g_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_xx_xz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_xx_xz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_xx_xz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_xx_xz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_xx_xz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_xx_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_x_z_0_0_x_x_z_xx, g_x_z_0_0_x_x_z_xy, g_x_z_0_0_x_x_z_xz, g_x_z_0_0_x_x_z_yy, g_x_z_0_0_x_x_z_yz, g_x_z_0_0_x_x_z_zz, g_xx_xz_z_xx, g_xx_xz_z_xy, g_xx_xz_z_xz, g_xx_xz_z_yy, g_xx_xz_z_yz, g_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_xx_xz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_xx_xz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_xx_xz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_xx_xz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_xx_xz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_xx_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_x_z_0_0_x_y_x_xx, g_x_z_0_0_x_y_x_xy, g_x_z_0_0_x_y_x_xz, g_x_z_0_0_x_y_x_yy, g_x_z_0_0_x_y_x_yz, g_x_z_0_0_x_y_x_zz, g_xx_yz_x_xx, g_xx_yz_x_xy, g_xx_yz_x_xz, g_xx_yz_x_yy, g_xx_yz_x_yz, g_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_xx_yz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_xx_yz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_xx_yz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_xx_yz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_xx_yz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_xx_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_x_z_0_0_x_y_y_xx, g_x_z_0_0_x_y_y_xy, g_x_z_0_0_x_y_y_xz, g_x_z_0_0_x_y_y_yy, g_x_z_0_0_x_y_y_yz, g_x_z_0_0_x_y_y_zz, g_xx_yz_y_xx, g_xx_yz_y_xy, g_xx_yz_y_xz, g_xx_yz_y_yy, g_xx_yz_y_yz, g_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_xx_yz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_xx_yz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_xx_yz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_xx_yz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_xx_yz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_xx_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_x_z_0_0_x_y_z_xx, g_x_z_0_0_x_y_z_xy, g_x_z_0_0_x_y_z_xz, g_x_z_0_0_x_y_z_yy, g_x_z_0_0_x_y_z_yz, g_x_z_0_0_x_y_z_zz, g_xx_yz_z_xx, g_xx_yz_z_xy, g_xx_yz_z_xz, g_xx_yz_z_yy, g_xx_yz_z_yz, g_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_xx_yz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_xx_yz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_xx_yz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_xx_yz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_xx_yz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_xx_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_x_z_0_0_x_z_x_xx, g_x_z_0_0_x_z_x_xy, g_x_z_0_0_x_z_x_xz, g_x_z_0_0_x_z_x_yy, g_x_z_0_0_x_z_x_yz, g_x_z_0_0_x_z_x_zz, g_xx_0_x_xx, g_xx_0_x_xy, g_xx_0_x_xz, g_xx_0_x_yy, g_xx_0_x_yz, g_xx_0_x_zz, g_xx_zz_x_xx, g_xx_zz_x_xy, g_xx_zz_x_xz, g_xx_zz_x_yy, g_xx_zz_x_yz, g_xx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_zz_x_xx[i] * b_exp - 2.0 * g_xx_0_x_xx[i] * a_exp + 4.0 * g_xx_zz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_zz_x_xy[i] * b_exp - 2.0 * g_xx_0_x_xy[i] * a_exp + 4.0 * g_xx_zz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_zz_x_xz[i] * b_exp - 2.0 * g_xx_0_x_xz[i] * a_exp + 4.0 * g_xx_zz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_zz_x_yy[i] * b_exp - 2.0 * g_xx_0_x_yy[i] * a_exp + 4.0 * g_xx_zz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_zz_x_yz[i] * b_exp - 2.0 * g_xx_0_x_yz[i] * a_exp + 4.0 * g_xx_zz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_zz_x_zz[i] * b_exp - 2.0 * g_xx_0_x_zz[i] * a_exp + 4.0 * g_xx_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_x_z_0_0_x_z_y_xx, g_x_z_0_0_x_z_y_xy, g_x_z_0_0_x_z_y_xz, g_x_z_0_0_x_z_y_yy, g_x_z_0_0_x_z_y_yz, g_x_z_0_0_x_z_y_zz, g_xx_0_y_xx, g_xx_0_y_xy, g_xx_0_y_xz, g_xx_0_y_yy, g_xx_0_y_yz, g_xx_0_y_zz, g_xx_zz_y_xx, g_xx_zz_y_xy, g_xx_zz_y_xz, g_xx_zz_y_yy, g_xx_zz_y_yz, g_xx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_zz_y_xx[i] * b_exp - 2.0 * g_xx_0_y_xx[i] * a_exp + 4.0 * g_xx_zz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_zz_y_xy[i] * b_exp - 2.0 * g_xx_0_y_xy[i] * a_exp + 4.0 * g_xx_zz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_zz_y_xz[i] * b_exp - 2.0 * g_xx_0_y_xz[i] * a_exp + 4.0 * g_xx_zz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_zz_y_yy[i] * b_exp - 2.0 * g_xx_0_y_yy[i] * a_exp + 4.0 * g_xx_zz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_zz_y_yz[i] * b_exp - 2.0 * g_xx_0_y_yz[i] * a_exp + 4.0 * g_xx_zz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_zz_y_zz[i] * b_exp - 2.0 * g_xx_0_y_zz[i] * a_exp + 4.0 * g_xx_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_x_z_0_0_x_z_z_xx, g_x_z_0_0_x_z_z_xy, g_x_z_0_0_x_z_z_xz, g_x_z_0_0_x_z_z_yy, g_x_z_0_0_x_z_z_yz, g_x_z_0_0_x_z_z_zz, g_xx_0_z_xx, g_xx_0_z_xy, g_xx_0_z_xz, g_xx_0_z_yy, g_xx_0_z_yz, g_xx_0_z_zz, g_xx_zz_z_xx, g_xx_zz_z_xy, g_xx_zz_z_xz, g_xx_zz_z_yy, g_xx_zz_z_yz, g_xx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_zz_z_xx[i] * b_exp - 2.0 * g_xx_0_z_xx[i] * a_exp + 4.0 * g_xx_zz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_zz_z_xy[i] * b_exp - 2.0 * g_xx_0_z_xy[i] * a_exp + 4.0 * g_xx_zz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_zz_z_xz[i] * b_exp - 2.0 * g_xx_0_z_xz[i] * a_exp + 4.0 * g_xx_zz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_zz_z_yy[i] * b_exp - 2.0 * g_xx_0_z_yy[i] * a_exp + 4.0 * g_xx_zz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_zz_z_yz[i] * b_exp - 2.0 * g_xx_0_z_yz[i] * a_exp + 4.0 * g_xx_zz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_zz_z_zz[i] * b_exp - 2.0 * g_xx_0_z_zz[i] * a_exp + 4.0 * g_xx_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_z_0_0_y_x_x_xx, g_x_z_0_0_y_x_x_xy, g_x_z_0_0_y_x_x_xz, g_x_z_0_0_y_x_x_yy, g_x_z_0_0_y_x_x_yz, g_x_z_0_0_y_x_x_zz, g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_x_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_z_0_0_y_x_y_xx, g_x_z_0_0_y_x_y_xy, g_x_z_0_0_y_x_y_xz, g_x_z_0_0_y_x_y_yy, g_x_z_0_0_y_x_y_yz, g_x_z_0_0_y_x_y_zz, g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_y_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_z_0_0_y_x_z_xx, g_x_z_0_0_y_x_z_xy, g_x_z_0_0_y_x_z_xz, g_x_z_0_0_y_x_z_yy, g_x_z_0_0_y_x_z_yz, g_x_z_0_0_y_x_z_zz, g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_z_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_z_0_0_y_y_x_xx, g_x_z_0_0_y_y_x_xy, g_x_z_0_0_y_y_x_xz, g_x_z_0_0_y_y_x_yy, g_x_z_0_0_y_y_x_yz, g_x_z_0_0_y_y_x_zz, g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_x_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_z_0_0_y_y_y_xx, g_x_z_0_0_y_y_y_xy, g_x_z_0_0_y_y_y_xz, g_x_z_0_0_y_y_y_yy, g_x_z_0_0_y_y_y_yz, g_x_z_0_0_y_y_y_zz, g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_y_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_z_0_0_y_y_z_xx, g_x_z_0_0_y_y_z_xy, g_x_z_0_0_y_y_z_xz, g_x_z_0_0_y_y_z_yy, g_x_z_0_0_y_y_z_yz, g_x_z_0_0_y_y_z_zz, g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_z_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_z_0_0_y_z_x_xx, g_x_z_0_0_y_z_x_xy, g_x_z_0_0_y_z_x_xz, g_x_z_0_0_y_z_x_yy, g_x_z_0_0_y_z_x_yz, g_x_z_0_0_y_z_x_zz, g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_zz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_zz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_zz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_zz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_zz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_z_0_0_y_z_y_xx, g_x_z_0_0_y_z_y_xy, g_x_z_0_0_y_z_y_xz, g_x_z_0_0_y_z_y_yy, g_x_z_0_0_y_z_y_yz, g_x_z_0_0_y_z_y_zz, g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_zz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_zz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_zz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_zz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_zz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_z_0_0_y_z_z_xx, g_x_z_0_0_y_z_z_xy, g_x_z_0_0_y_z_z_xz, g_x_z_0_0_y_z_z_yy, g_x_z_0_0_y_z_z_yz, g_x_z_0_0_y_z_z_zz, g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_zz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_zz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_zz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_zz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_zz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_z_0_0_z_x_x_xx, g_x_z_0_0_z_x_x_xy, g_x_z_0_0_z_x_x_xz, g_x_z_0_0_z_x_x_yy, g_x_z_0_0_z_x_x_yz, g_x_z_0_0_z_x_x_zz, g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_x_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_z_0_0_z_x_y_xx, g_x_z_0_0_z_x_y_xy, g_x_z_0_0_z_x_y_xz, g_x_z_0_0_z_x_y_yy, g_x_z_0_0_z_x_y_yz, g_x_z_0_0_z_x_y_zz, g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_y_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_z_0_0_z_x_z_xx, g_x_z_0_0_z_x_z_xy, g_x_z_0_0_z_x_z_xz, g_x_z_0_0_z_x_z_yy, g_x_z_0_0_z_x_z_yz, g_x_z_0_0_z_x_z_zz, g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_z_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_z_0_0_z_y_x_xx, g_x_z_0_0_z_y_x_xy, g_x_z_0_0_z_y_x_xz, g_x_z_0_0_z_y_x_yy, g_x_z_0_0_z_y_x_yz, g_x_z_0_0_z_y_x_zz, g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_x_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_z_0_0_z_y_y_xx, g_x_z_0_0_z_y_y_xy, g_x_z_0_0_z_y_y_xz, g_x_z_0_0_z_y_y_yy, g_x_z_0_0_z_y_y_yz, g_x_z_0_0_z_y_y_zz, g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_y_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_z_0_0_z_y_z_xx, g_x_z_0_0_z_y_z_xy, g_x_z_0_0_z_y_z_xz, g_x_z_0_0_z_y_z_yy, g_x_z_0_0_z_y_z_yz, g_x_z_0_0_z_y_z_zz, g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_z_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_z_0_0_z_z_x_xx, g_x_z_0_0_z_z_x_xy, g_x_z_0_0_z_z_x_xz, g_x_z_0_0_z_z_x_yy, g_x_z_0_0_z_z_x_yz, g_x_z_0_0_z_z_x_zz, g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_zz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_zz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_zz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_zz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_zz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_z_0_0_z_z_y_xx, g_x_z_0_0_z_z_y_xy, g_x_z_0_0_z_z_y_xz, g_x_z_0_0_z_z_y_yy, g_x_z_0_0_z_z_y_yz, g_x_z_0_0_z_z_y_zz, g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_zz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_zz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_zz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_zz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_zz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_z_0_0_z_z_z_xx, g_x_z_0_0_z_z_z_xy, g_x_z_0_0_z_z_z_xz, g_x_z_0_0_z_z_z_yy, g_x_z_0_0_z_z_z_yz, g_x_z_0_0_z_z_z_zz, g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_zz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_zz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_zz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_zz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_zz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz, g_y_x_0_0_x_x_x_xx, g_y_x_0_0_x_x_x_xy, g_y_x_0_0_x_x_x_xz, g_y_x_0_0_x_x_x_yy, g_y_x_0_0_x_x_x_yz, g_y_x_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_xx_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_xx_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_xx_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_xx_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_xx_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz, g_y_x_0_0_x_x_y_xx, g_y_x_0_0_x_x_y_xy, g_y_x_0_0_x_x_y_xz, g_y_x_0_0_x_x_y_yy, g_y_x_0_0_x_x_y_yz, g_y_x_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_xx_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_xx_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_xx_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_xx_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_xx_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz, g_y_x_0_0_x_x_z_xx, g_y_x_0_0_x_x_z_xy, g_y_x_0_0_x_x_z_xz, g_y_x_0_0_x_x_z_yy, g_y_x_0_0_x_x_z_yz, g_y_x_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_xx_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_xx_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_xx_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_xx_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_xx_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz, g_y_x_0_0_x_y_x_xx, g_y_x_0_0_x_y_x_xy, g_y_x_0_0_x_y_x_xz, g_y_x_0_0_x_y_x_yy, g_y_x_0_0_x_y_x_yz, g_y_x_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_x_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz, g_y_x_0_0_x_y_y_xx, g_y_x_0_0_x_y_y_xy, g_y_x_0_0_x_y_y_xz, g_y_x_0_0_x_y_y_yy, g_y_x_0_0_x_y_y_yz, g_y_x_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_y_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz, g_y_x_0_0_x_y_z_xx, g_y_x_0_0_x_y_z_xy, g_y_x_0_0_x_y_z_xz, g_y_x_0_0_x_y_z_yy, g_y_x_0_0_x_y_z_yz, g_y_x_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_z_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz, g_y_x_0_0_x_z_x_xx, g_y_x_0_0_x_z_x_xy, g_y_x_0_0_x_z_x_xz, g_y_x_0_0_x_z_x_yy, g_y_x_0_0_x_z_x_yz, g_y_x_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_x_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz, g_y_x_0_0_x_z_y_xx, g_y_x_0_0_x_z_y_xy, g_y_x_0_0_x_z_y_xz, g_y_x_0_0_x_z_y_yy, g_y_x_0_0_x_z_y_yz, g_y_x_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_y_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz, g_y_x_0_0_x_z_z_xx, g_y_x_0_0_x_z_z_xy, g_y_x_0_0_x_z_z_xz, g_y_x_0_0_x_z_z_yy, g_y_x_0_0_x_z_z_yz, g_y_x_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_z_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_y_x_0_0_y_x_x_xx, g_y_x_0_0_y_x_x_xy, g_y_x_0_0_y_x_x_xz, g_y_x_0_0_y_x_x_yy, g_y_x_0_0_y_x_x_yz, g_y_x_0_0_y_x_x_zz, g_yy_0_x_xx, g_yy_0_x_xy, g_yy_0_x_xz, g_yy_0_x_yy, g_yy_0_x_yz, g_yy_0_x_zz, g_yy_xx_x_xx, g_yy_xx_x_xy, g_yy_xx_x_xz, g_yy_xx_x_yy, g_yy_xx_x_yz, g_yy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_xx_x_xx[i] * b_exp - 2.0 * g_yy_0_x_xx[i] * a_exp + 4.0 * g_yy_xx_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_xx_x_xy[i] * b_exp - 2.0 * g_yy_0_x_xy[i] * a_exp + 4.0 * g_yy_xx_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_xx_x_xz[i] * b_exp - 2.0 * g_yy_0_x_xz[i] * a_exp + 4.0 * g_yy_xx_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_xx_x_yy[i] * b_exp - 2.0 * g_yy_0_x_yy[i] * a_exp + 4.0 * g_yy_xx_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_xx_x_yz[i] * b_exp - 2.0 * g_yy_0_x_yz[i] * a_exp + 4.0 * g_yy_xx_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_xx_x_zz[i] * b_exp - 2.0 * g_yy_0_x_zz[i] * a_exp + 4.0 * g_yy_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_y_x_0_0_y_x_y_xx, g_y_x_0_0_y_x_y_xy, g_y_x_0_0_y_x_y_xz, g_y_x_0_0_y_x_y_yy, g_y_x_0_0_y_x_y_yz, g_y_x_0_0_y_x_y_zz, g_yy_0_y_xx, g_yy_0_y_xy, g_yy_0_y_xz, g_yy_0_y_yy, g_yy_0_y_yz, g_yy_0_y_zz, g_yy_xx_y_xx, g_yy_xx_y_xy, g_yy_xx_y_xz, g_yy_xx_y_yy, g_yy_xx_y_yz, g_yy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_xx_y_xx[i] * b_exp - 2.0 * g_yy_0_y_xx[i] * a_exp + 4.0 * g_yy_xx_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_xx_y_xy[i] * b_exp - 2.0 * g_yy_0_y_xy[i] * a_exp + 4.0 * g_yy_xx_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_xx_y_xz[i] * b_exp - 2.0 * g_yy_0_y_xz[i] * a_exp + 4.0 * g_yy_xx_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_xx_y_yy[i] * b_exp - 2.0 * g_yy_0_y_yy[i] * a_exp + 4.0 * g_yy_xx_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_xx_y_yz[i] * b_exp - 2.0 * g_yy_0_y_yz[i] * a_exp + 4.0 * g_yy_xx_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_xx_y_zz[i] * b_exp - 2.0 * g_yy_0_y_zz[i] * a_exp + 4.0 * g_yy_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_y_x_0_0_y_x_z_xx, g_y_x_0_0_y_x_z_xy, g_y_x_0_0_y_x_z_xz, g_y_x_0_0_y_x_z_yy, g_y_x_0_0_y_x_z_yz, g_y_x_0_0_y_x_z_zz, g_yy_0_z_xx, g_yy_0_z_xy, g_yy_0_z_xz, g_yy_0_z_yy, g_yy_0_z_yz, g_yy_0_z_zz, g_yy_xx_z_xx, g_yy_xx_z_xy, g_yy_xx_z_xz, g_yy_xx_z_yy, g_yy_xx_z_yz, g_yy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_xx_z_xx[i] * b_exp - 2.0 * g_yy_0_z_xx[i] * a_exp + 4.0 * g_yy_xx_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_xx_z_xy[i] * b_exp - 2.0 * g_yy_0_z_xy[i] * a_exp + 4.0 * g_yy_xx_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_xx_z_xz[i] * b_exp - 2.0 * g_yy_0_z_xz[i] * a_exp + 4.0 * g_yy_xx_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_xx_z_yy[i] * b_exp - 2.0 * g_yy_0_z_yy[i] * a_exp + 4.0 * g_yy_xx_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_xx_z_yz[i] * b_exp - 2.0 * g_yy_0_z_yz[i] * a_exp + 4.0 * g_yy_xx_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_xx_z_zz[i] * b_exp - 2.0 * g_yy_0_z_zz[i] * a_exp + 4.0 * g_yy_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_y_x_0_0_y_y_x_xx, g_y_x_0_0_y_y_x_xy, g_y_x_0_0_y_y_x_xz, g_y_x_0_0_y_y_x_yy, g_y_x_0_0_y_y_x_yz, g_y_x_0_0_y_y_x_zz, g_yy_xy_x_xx, g_yy_xy_x_xy, g_yy_xy_x_xz, g_yy_xy_x_yy, g_yy_xy_x_yz, g_yy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_yy_xy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_yy_xy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_yy_xy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_yy_xy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_yy_xy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_yy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_y_x_0_0_y_y_y_xx, g_y_x_0_0_y_y_y_xy, g_y_x_0_0_y_y_y_xz, g_y_x_0_0_y_y_y_yy, g_y_x_0_0_y_y_y_yz, g_y_x_0_0_y_y_y_zz, g_yy_xy_y_xx, g_yy_xy_y_xy, g_yy_xy_y_xz, g_yy_xy_y_yy, g_yy_xy_y_yz, g_yy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_yy_xy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_yy_xy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_yy_xy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_yy_xy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_yy_xy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_yy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_y_x_0_0_y_y_z_xx, g_y_x_0_0_y_y_z_xy, g_y_x_0_0_y_y_z_xz, g_y_x_0_0_y_y_z_yy, g_y_x_0_0_y_y_z_yz, g_y_x_0_0_y_y_z_zz, g_yy_xy_z_xx, g_yy_xy_z_xy, g_yy_xy_z_xz, g_yy_xy_z_yy, g_yy_xy_z_yz, g_yy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_yy_xy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_yy_xy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_yy_xy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_yy_xy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_yy_xy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_yy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_y_x_0_0_y_z_x_xx, g_y_x_0_0_y_z_x_xy, g_y_x_0_0_y_z_x_xz, g_y_x_0_0_y_z_x_yy, g_y_x_0_0_y_z_x_yz, g_y_x_0_0_y_z_x_zz, g_yy_xz_x_xx, g_yy_xz_x_xy, g_yy_xz_x_xz, g_yy_xz_x_yy, g_yy_xz_x_yz, g_yy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_yy_xz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_yy_xz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_yy_xz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_yy_xz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_yy_xz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_yy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_y_x_0_0_y_z_y_xx, g_y_x_0_0_y_z_y_xy, g_y_x_0_0_y_z_y_xz, g_y_x_0_0_y_z_y_yy, g_y_x_0_0_y_z_y_yz, g_y_x_0_0_y_z_y_zz, g_yy_xz_y_xx, g_yy_xz_y_xy, g_yy_xz_y_xz, g_yy_xz_y_yy, g_yy_xz_y_yz, g_yy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_yy_xz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_yy_xz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_yy_xz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_yy_xz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_yy_xz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_yy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_y_x_0_0_y_z_z_xx, g_y_x_0_0_y_z_z_xy, g_y_x_0_0_y_z_z_xz, g_y_x_0_0_y_z_z_yy, g_y_x_0_0_y_z_z_yz, g_y_x_0_0_y_z_z_zz, g_yy_xz_z_xx, g_yy_xz_z_xy, g_yy_xz_z_xz, g_yy_xz_z_yy, g_yy_xz_z_yz, g_yy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_yy_xz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_yy_xz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_yy_xz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_yy_xz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_yy_xz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_yy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_x_0_0_z_x_x_xx, g_y_x_0_0_z_x_x_xy, g_y_x_0_0_z_x_x_xz, g_y_x_0_0_z_x_x_yy, g_y_x_0_0_z_x_x_yz, g_y_x_0_0_z_x_x_zz, g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_xx_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_xx_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_xx_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_xx_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_xx_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_x_0_0_z_x_y_xx, g_y_x_0_0_z_x_y_xy, g_y_x_0_0_z_x_y_xz, g_y_x_0_0_z_x_y_yy, g_y_x_0_0_z_x_y_yz, g_y_x_0_0_z_x_y_zz, g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_xx_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_xx_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_xx_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_xx_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_xx_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_x_0_0_z_x_z_xx, g_y_x_0_0_z_x_z_xy, g_y_x_0_0_z_x_z_xz, g_y_x_0_0_z_x_z_yy, g_y_x_0_0_z_x_z_yz, g_y_x_0_0_z_x_z_zz, g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_xx_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_xx_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_xx_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_xx_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_xx_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_x_0_0_z_y_x_xx, g_y_x_0_0_z_y_x_xy, g_y_x_0_0_z_y_x_xz, g_y_x_0_0_z_y_x_yy, g_y_x_0_0_z_y_x_yz, g_y_x_0_0_z_y_x_zz, g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_x_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_x_0_0_z_y_y_xx, g_y_x_0_0_z_y_y_xy, g_y_x_0_0_z_y_y_xz, g_y_x_0_0_z_y_y_yy, g_y_x_0_0_z_y_y_yz, g_y_x_0_0_z_y_y_zz, g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_y_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_x_0_0_z_y_z_xx, g_y_x_0_0_z_y_z_xy, g_y_x_0_0_z_y_z_xz, g_y_x_0_0_z_y_z_yy, g_y_x_0_0_z_y_z_yz, g_y_x_0_0_z_y_z_zz, g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_z_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_x_0_0_z_z_x_xx, g_y_x_0_0_z_z_x_xy, g_y_x_0_0_z_z_x_xz, g_y_x_0_0_z_z_x_yy, g_y_x_0_0_z_z_x_yz, g_y_x_0_0_z_z_x_zz, g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_x_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_x_0_0_z_z_y_xx, g_y_x_0_0_z_z_y_xy, g_y_x_0_0_z_z_y_xz, g_y_x_0_0_z_z_y_yy, g_y_x_0_0_z_z_y_yz, g_y_x_0_0_z_z_y_zz, g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_y_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_x_0_0_z_z_z_xx, g_y_x_0_0_z_z_z_xy, g_y_x_0_0_z_z_z_xz, g_y_x_0_0_z_z_z_yy, g_y_x_0_0_z_z_z_yz, g_y_x_0_0_z_z_z_zz, g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_z_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz, g_y_y_0_0_x_x_x_xx, g_y_y_0_0_x_x_x_xy, g_y_y_0_0_x_x_x_xz, g_y_y_0_0_x_x_x_yy, g_y_y_0_0_x_x_x_yz, g_y_y_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_x_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz, g_y_y_0_0_x_x_y_xx, g_y_y_0_0_x_x_y_xy, g_y_y_0_0_x_x_y_xz, g_y_y_0_0_x_x_y_yy, g_y_y_0_0_x_x_y_yz, g_y_y_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_y_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz, g_y_y_0_0_x_x_z_xx, g_y_y_0_0_x_x_z_xy, g_y_y_0_0_x_x_z_xz, g_y_y_0_0_x_x_z_yy, g_y_y_0_0_x_x_z_yz, g_y_y_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_z_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz, g_y_y_0_0_x_y_x_xx, g_y_y_0_0_x_y_x_xy, g_y_y_0_0_x_y_x_xz, g_y_y_0_0_x_y_x_yy, g_y_y_0_0_x_y_x_yz, g_y_y_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_yy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_yy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_yy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_yy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_yy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz, g_y_y_0_0_x_y_y_xx, g_y_y_0_0_x_y_y_xy, g_y_y_0_0_x_y_y_xz, g_y_y_0_0_x_y_y_yy, g_y_y_0_0_x_y_y_yz, g_y_y_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_yy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_yy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_yy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_yy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_yy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz, g_y_y_0_0_x_y_z_xx, g_y_y_0_0_x_y_z_xy, g_y_y_0_0_x_y_z_xz, g_y_y_0_0_x_y_z_yy, g_y_y_0_0_x_y_z_yz, g_y_y_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_yy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_yy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_yy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_yy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_yy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz, g_y_y_0_0_x_z_x_xx, g_y_y_0_0_x_z_x_xy, g_y_y_0_0_x_z_x_xz, g_y_y_0_0_x_z_x_yy, g_y_y_0_0_x_z_x_yz, g_y_y_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_x_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz, g_y_y_0_0_x_z_y_xx, g_y_y_0_0_x_z_y_xy, g_y_y_0_0_x_z_y_xz, g_y_y_0_0_x_z_y_yy, g_y_y_0_0_x_z_y_yz, g_y_y_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_y_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz, g_y_y_0_0_x_z_z_xx, g_y_y_0_0_x_z_z_xy, g_y_y_0_0_x_z_z_xz, g_y_y_0_0_x_z_z_yy, g_y_y_0_0_x_z_z_yz, g_y_y_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_z_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_y_y_0_0_y_x_x_xx, g_y_y_0_0_y_x_x_xy, g_y_y_0_0_y_x_x_xz, g_y_y_0_0_y_x_x_yy, g_y_y_0_0_y_x_x_yz, g_y_y_0_0_y_x_x_zz, g_yy_xy_x_xx, g_yy_xy_x_xy, g_yy_xy_x_xz, g_yy_xy_x_yy, g_yy_xy_x_yz, g_yy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_yy_xy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_yy_xy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_yy_xy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_yy_xy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_yy_xy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_yy_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_y_y_0_0_y_x_y_xx, g_y_y_0_0_y_x_y_xy, g_y_y_0_0_y_x_y_xz, g_y_y_0_0_y_x_y_yy, g_y_y_0_0_y_x_y_yz, g_y_y_0_0_y_x_y_zz, g_yy_xy_y_xx, g_yy_xy_y_xy, g_yy_xy_y_xz, g_yy_xy_y_yy, g_yy_xy_y_yz, g_yy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_yy_xy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_yy_xy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_yy_xy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_yy_xy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_yy_xy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_yy_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_y_y_0_0_y_x_z_xx, g_y_y_0_0_y_x_z_xy, g_y_y_0_0_y_x_z_xz, g_y_y_0_0_y_x_z_yy, g_y_y_0_0_y_x_z_yz, g_y_y_0_0_y_x_z_zz, g_yy_xy_z_xx, g_yy_xy_z_xy, g_yy_xy_z_xz, g_yy_xy_z_yy, g_yy_xy_z_yz, g_yy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_yy_xy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_yy_xy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_yy_xy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_yy_xy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_yy_xy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_yy_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_y_y_0_0_y_y_x_xx, g_y_y_0_0_y_y_x_xy, g_y_y_0_0_y_y_x_xz, g_y_y_0_0_y_y_x_yy, g_y_y_0_0_y_y_x_yz, g_y_y_0_0_y_y_x_zz, g_yy_0_x_xx, g_yy_0_x_xy, g_yy_0_x_xz, g_yy_0_x_yy, g_yy_0_x_yz, g_yy_0_x_zz, g_yy_yy_x_xx, g_yy_yy_x_xy, g_yy_yy_x_xz, g_yy_yy_x_yy, g_yy_yy_x_yz, g_yy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_yy_x_xx[i] * b_exp - 2.0 * g_yy_0_x_xx[i] * a_exp + 4.0 * g_yy_yy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_yy_x_xy[i] * b_exp - 2.0 * g_yy_0_x_xy[i] * a_exp + 4.0 * g_yy_yy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_yy_x_xz[i] * b_exp - 2.0 * g_yy_0_x_xz[i] * a_exp + 4.0 * g_yy_yy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_yy_x_yy[i] * b_exp - 2.0 * g_yy_0_x_yy[i] * a_exp + 4.0 * g_yy_yy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_yy_x_yz[i] * b_exp - 2.0 * g_yy_0_x_yz[i] * a_exp + 4.0 * g_yy_yy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_yy_x_zz[i] * b_exp - 2.0 * g_yy_0_x_zz[i] * a_exp + 4.0 * g_yy_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_y_y_0_0_y_y_y_xx, g_y_y_0_0_y_y_y_xy, g_y_y_0_0_y_y_y_xz, g_y_y_0_0_y_y_y_yy, g_y_y_0_0_y_y_y_yz, g_y_y_0_0_y_y_y_zz, g_yy_0_y_xx, g_yy_0_y_xy, g_yy_0_y_xz, g_yy_0_y_yy, g_yy_0_y_yz, g_yy_0_y_zz, g_yy_yy_y_xx, g_yy_yy_y_xy, g_yy_yy_y_xz, g_yy_yy_y_yy, g_yy_yy_y_yz, g_yy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_yy_y_xx[i] * b_exp - 2.0 * g_yy_0_y_xx[i] * a_exp + 4.0 * g_yy_yy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_yy_y_xy[i] * b_exp - 2.0 * g_yy_0_y_xy[i] * a_exp + 4.0 * g_yy_yy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_yy_y_xz[i] * b_exp - 2.0 * g_yy_0_y_xz[i] * a_exp + 4.0 * g_yy_yy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_yy_y_yy[i] * b_exp - 2.0 * g_yy_0_y_yy[i] * a_exp + 4.0 * g_yy_yy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_yy_y_yz[i] * b_exp - 2.0 * g_yy_0_y_yz[i] * a_exp + 4.0 * g_yy_yy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_yy_y_zz[i] * b_exp - 2.0 * g_yy_0_y_zz[i] * a_exp + 4.0 * g_yy_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_y_y_0_0_y_y_z_xx, g_y_y_0_0_y_y_z_xy, g_y_y_0_0_y_y_z_xz, g_y_y_0_0_y_y_z_yy, g_y_y_0_0_y_y_z_yz, g_y_y_0_0_y_y_z_zz, g_yy_0_z_xx, g_yy_0_z_xy, g_yy_0_z_xz, g_yy_0_z_yy, g_yy_0_z_yz, g_yy_0_z_zz, g_yy_yy_z_xx, g_yy_yy_z_xy, g_yy_yy_z_xz, g_yy_yy_z_yy, g_yy_yy_z_yz, g_yy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_yy_z_xx[i] * b_exp - 2.0 * g_yy_0_z_xx[i] * a_exp + 4.0 * g_yy_yy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_yy_z_xy[i] * b_exp - 2.0 * g_yy_0_z_xy[i] * a_exp + 4.0 * g_yy_yy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_yy_z_xz[i] * b_exp - 2.0 * g_yy_0_z_xz[i] * a_exp + 4.0 * g_yy_yy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_yy_z_yy[i] * b_exp - 2.0 * g_yy_0_z_yy[i] * a_exp + 4.0 * g_yy_yy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_yy_z_yz[i] * b_exp - 2.0 * g_yy_0_z_yz[i] * a_exp + 4.0 * g_yy_yy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_yy_z_zz[i] * b_exp - 2.0 * g_yy_0_z_zz[i] * a_exp + 4.0 * g_yy_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_y_y_0_0_y_z_x_xx, g_y_y_0_0_y_z_x_xy, g_y_y_0_0_y_z_x_xz, g_y_y_0_0_y_z_x_yy, g_y_y_0_0_y_z_x_yz, g_y_y_0_0_y_z_x_zz, g_yy_yz_x_xx, g_yy_yz_x_xy, g_yy_yz_x_xz, g_yy_yz_x_yy, g_yy_yz_x_yz, g_yy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_yy_yz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_yy_yz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_yy_yz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_yy_yz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_yy_yz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_yy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_y_y_0_0_y_z_y_xx, g_y_y_0_0_y_z_y_xy, g_y_y_0_0_y_z_y_xz, g_y_y_0_0_y_z_y_yy, g_y_y_0_0_y_z_y_yz, g_y_y_0_0_y_z_y_zz, g_yy_yz_y_xx, g_yy_yz_y_xy, g_yy_yz_y_xz, g_yy_yz_y_yy, g_yy_yz_y_yz, g_yy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_yy_yz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_yy_yz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_yy_yz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_yy_yz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_yy_yz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_yy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_y_y_0_0_y_z_z_xx, g_y_y_0_0_y_z_z_xy, g_y_y_0_0_y_z_z_xz, g_y_y_0_0_y_z_z_yy, g_y_y_0_0_y_z_z_yz, g_y_y_0_0_y_z_z_zz, g_yy_yz_z_xx, g_yy_yz_z_xy, g_yy_yz_z_xz, g_yy_yz_z_yy, g_yy_yz_z_yz, g_yy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_yy_yz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_yy_yz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_yy_yz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_yy_yz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_yy_yz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_yy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_y_0_0_z_x_x_xx, g_y_y_0_0_z_x_x_xy, g_y_y_0_0_z_x_x_xz, g_y_y_0_0_z_x_x_yy, g_y_y_0_0_z_x_x_yz, g_y_y_0_0_z_x_x_zz, g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_x_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_y_0_0_z_x_y_xx, g_y_y_0_0_z_x_y_xy, g_y_y_0_0_z_x_y_xz, g_y_y_0_0_z_x_y_yy, g_y_y_0_0_z_x_y_yz, g_y_y_0_0_z_x_y_zz, g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_y_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_y_0_0_z_x_z_xx, g_y_y_0_0_z_x_z_xy, g_y_y_0_0_z_x_z_xz, g_y_y_0_0_z_x_z_yy, g_y_y_0_0_z_x_z_yz, g_y_y_0_0_z_x_z_zz, g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_z_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_y_0_0_z_y_x_xx, g_y_y_0_0_z_y_x_xy, g_y_y_0_0_z_y_x_xz, g_y_y_0_0_z_y_x_yy, g_y_y_0_0_z_y_x_yz, g_y_y_0_0_z_y_x_zz, g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_yy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_yy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_yy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_yy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_yy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_y_0_0_z_y_y_xx, g_y_y_0_0_z_y_y_xy, g_y_y_0_0_z_y_y_xz, g_y_y_0_0_z_y_y_yy, g_y_y_0_0_z_y_y_yz, g_y_y_0_0_z_y_y_zz, g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_yy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_yy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_yy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_yy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_yy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_y_0_0_z_y_z_xx, g_y_y_0_0_z_y_z_xy, g_y_y_0_0_z_y_z_xz, g_y_y_0_0_z_y_z_yy, g_y_y_0_0_z_y_z_yz, g_y_y_0_0_z_y_z_zz, g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_yy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_yy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_yy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_yy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_yy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_y_0_0_z_z_x_xx, g_y_y_0_0_z_z_x_xy, g_y_y_0_0_z_z_x_xz, g_y_y_0_0_z_z_x_yy, g_y_y_0_0_z_z_x_yz, g_y_y_0_0_z_z_x_zz, g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_x_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_y_0_0_z_z_y_xx, g_y_y_0_0_z_z_y_xy, g_y_y_0_0_z_z_y_xz, g_y_y_0_0_z_z_y_yy, g_y_y_0_0_z_z_y_yz, g_y_y_0_0_z_z_y_zz, g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_y_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_y_0_0_z_z_z_xx, g_y_y_0_0_z_z_z_xy, g_y_y_0_0_z_z_z_xz, g_y_y_0_0_z_z_z_yy, g_y_y_0_0_z_z_z_yz, g_y_y_0_0_z_z_z_zz, g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_z_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz, g_y_z_0_0_x_x_x_xx, g_y_z_0_0_x_x_x_xy, g_y_z_0_0_x_x_x_xz, g_y_z_0_0_x_x_x_yy, g_y_z_0_0_x_x_x_yz, g_y_z_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_x_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz, g_y_z_0_0_x_x_y_xx, g_y_z_0_0_x_x_y_xy, g_y_z_0_0_x_x_y_xz, g_y_z_0_0_x_x_y_yy, g_y_z_0_0_x_x_y_yz, g_y_z_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_y_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz, g_y_z_0_0_x_x_z_xx, g_y_z_0_0_x_x_z_xy, g_y_z_0_0_x_x_z_xz, g_y_z_0_0_x_x_z_yy, g_y_z_0_0_x_x_z_yz, g_y_z_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_z_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz, g_y_z_0_0_x_y_x_xx, g_y_z_0_0_x_y_x_xy, g_y_z_0_0_x_y_x_xz, g_y_z_0_0_x_y_x_yy, g_y_z_0_0_x_y_x_yz, g_y_z_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_x_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz, g_y_z_0_0_x_y_y_xx, g_y_z_0_0_x_y_y_xy, g_y_z_0_0_x_y_y_xz, g_y_z_0_0_x_y_y_yy, g_y_z_0_0_x_y_y_yz, g_y_z_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_y_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz, g_y_z_0_0_x_y_z_xx, g_y_z_0_0_x_y_z_xy, g_y_z_0_0_x_y_z_xz, g_y_z_0_0_x_y_z_yy, g_y_z_0_0_x_y_z_yz, g_y_z_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_z_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz, g_y_z_0_0_x_z_x_xx, g_y_z_0_0_x_z_x_xy, g_y_z_0_0_x_z_x_xz, g_y_z_0_0_x_z_x_yy, g_y_z_0_0_x_z_x_yz, g_y_z_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_x_xx[i] = -2.0 * g_xy_0_x_xx[i] * a_exp + 4.0 * g_xy_zz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_xy[i] = -2.0 * g_xy_0_x_xy[i] * a_exp + 4.0 * g_xy_zz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_xz[i] = -2.0 * g_xy_0_x_xz[i] * a_exp + 4.0 * g_xy_zz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_yy[i] = -2.0 * g_xy_0_x_yy[i] * a_exp + 4.0 * g_xy_zz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_yz[i] = -2.0 * g_xy_0_x_yz[i] * a_exp + 4.0 * g_xy_zz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_zz[i] = -2.0 * g_xy_0_x_zz[i] * a_exp + 4.0 * g_xy_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz, g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz, g_y_z_0_0_x_z_y_xx, g_y_z_0_0_x_z_y_xy, g_y_z_0_0_x_z_y_xz, g_y_z_0_0_x_z_y_yy, g_y_z_0_0_x_z_y_yz, g_y_z_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_y_xx[i] = -2.0 * g_xy_0_y_xx[i] * a_exp + 4.0 * g_xy_zz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_xy[i] = -2.0 * g_xy_0_y_xy[i] * a_exp + 4.0 * g_xy_zz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_xz[i] = -2.0 * g_xy_0_y_xz[i] * a_exp + 4.0 * g_xy_zz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_yy[i] = -2.0 * g_xy_0_y_yy[i] * a_exp + 4.0 * g_xy_zz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_yz[i] = -2.0 * g_xy_0_y_yz[i] * a_exp + 4.0 * g_xy_zz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_zz[i] = -2.0 * g_xy_0_y_zz[i] * a_exp + 4.0 * g_xy_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz, g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz, g_y_z_0_0_x_z_z_xx, g_y_z_0_0_x_z_z_xy, g_y_z_0_0_x_z_z_xz, g_y_z_0_0_x_z_z_yy, g_y_z_0_0_x_z_z_yz, g_y_z_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_z_xx[i] = -2.0 * g_xy_0_z_xx[i] * a_exp + 4.0 * g_xy_zz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_xy[i] = -2.0 * g_xy_0_z_xy[i] * a_exp + 4.0 * g_xy_zz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_xz[i] = -2.0 * g_xy_0_z_xz[i] * a_exp + 4.0 * g_xy_zz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_yy[i] = -2.0 * g_xy_0_z_yy[i] * a_exp + 4.0 * g_xy_zz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_yz[i] = -2.0 * g_xy_0_z_yz[i] * a_exp + 4.0 * g_xy_zz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_zz[i] = -2.0 * g_xy_0_z_zz[i] * a_exp + 4.0 * g_xy_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_y_z_0_0_y_x_x_xx, g_y_z_0_0_y_x_x_xy, g_y_z_0_0_y_x_x_xz, g_y_z_0_0_y_x_x_yy, g_y_z_0_0_y_x_x_yz, g_y_z_0_0_y_x_x_zz, g_yy_xz_x_xx, g_yy_xz_x_xy, g_yy_xz_x_xz, g_yy_xz_x_yy, g_yy_xz_x_yz, g_yy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_yy_xz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_yy_xz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_yy_xz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_yy_xz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_yy_xz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_yy_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_y_z_0_0_y_x_y_xx, g_y_z_0_0_y_x_y_xy, g_y_z_0_0_y_x_y_xz, g_y_z_0_0_y_x_y_yy, g_y_z_0_0_y_x_y_yz, g_y_z_0_0_y_x_y_zz, g_yy_xz_y_xx, g_yy_xz_y_xy, g_yy_xz_y_xz, g_yy_xz_y_yy, g_yy_xz_y_yz, g_yy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_yy_xz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_yy_xz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_yy_xz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_yy_xz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_yy_xz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_yy_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_y_z_0_0_y_x_z_xx, g_y_z_0_0_y_x_z_xy, g_y_z_0_0_y_x_z_xz, g_y_z_0_0_y_x_z_yy, g_y_z_0_0_y_x_z_yz, g_y_z_0_0_y_x_z_zz, g_yy_xz_z_xx, g_yy_xz_z_xy, g_yy_xz_z_xz, g_yy_xz_z_yy, g_yy_xz_z_yz, g_yy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_yy_xz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_yy_xz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_yy_xz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_yy_xz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_yy_xz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_yy_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_y_z_0_0_y_y_x_xx, g_y_z_0_0_y_y_x_xy, g_y_z_0_0_y_y_x_xz, g_y_z_0_0_y_y_x_yy, g_y_z_0_0_y_y_x_yz, g_y_z_0_0_y_y_x_zz, g_yy_yz_x_xx, g_yy_yz_x_xy, g_yy_yz_x_xz, g_yy_yz_x_yy, g_yy_yz_x_yz, g_yy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_yy_yz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_yy_yz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_yy_yz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_yy_yz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_yy_yz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_yy_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_y_z_0_0_y_y_y_xx, g_y_z_0_0_y_y_y_xy, g_y_z_0_0_y_y_y_xz, g_y_z_0_0_y_y_y_yy, g_y_z_0_0_y_y_y_yz, g_y_z_0_0_y_y_y_zz, g_yy_yz_y_xx, g_yy_yz_y_xy, g_yy_yz_y_xz, g_yy_yz_y_yy, g_yy_yz_y_yz, g_yy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_yy_yz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_yy_yz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_yy_yz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_yy_yz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_yy_yz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_yy_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_y_z_0_0_y_y_z_xx, g_y_z_0_0_y_y_z_xy, g_y_z_0_0_y_y_z_xz, g_y_z_0_0_y_y_z_yy, g_y_z_0_0_y_y_z_yz, g_y_z_0_0_y_y_z_zz, g_yy_yz_z_xx, g_yy_yz_z_xy, g_yy_yz_z_xz, g_yy_yz_z_yy, g_yy_yz_z_yz, g_yy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_yy_yz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_yy_yz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_yy_yz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_yy_yz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_yy_yz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_yy_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_y_z_0_0_y_z_x_xx, g_y_z_0_0_y_z_x_xy, g_y_z_0_0_y_z_x_xz, g_y_z_0_0_y_z_x_yy, g_y_z_0_0_y_z_x_yz, g_y_z_0_0_y_z_x_zz, g_yy_0_x_xx, g_yy_0_x_xy, g_yy_0_x_xz, g_yy_0_x_yy, g_yy_0_x_yz, g_yy_0_x_zz, g_yy_zz_x_xx, g_yy_zz_x_xy, g_yy_zz_x_xz, g_yy_zz_x_yy, g_yy_zz_x_yz, g_yy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_zz_x_xx[i] * b_exp - 2.0 * g_yy_0_x_xx[i] * a_exp + 4.0 * g_yy_zz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_zz_x_xy[i] * b_exp - 2.0 * g_yy_0_x_xy[i] * a_exp + 4.0 * g_yy_zz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_zz_x_xz[i] * b_exp - 2.0 * g_yy_0_x_xz[i] * a_exp + 4.0 * g_yy_zz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_zz_x_yy[i] * b_exp - 2.0 * g_yy_0_x_yy[i] * a_exp + 4.0 * g_yy_zz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_zz_x_yz[i] * b_exp - 2.0 * g_yy_0_x_yz[i] * a_exp + 4.0 * g_yy_zz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_zz_x_zz[i] * b_exp - 2.0 * g_yy_0_x_zz[i] * a_exp + 4.0 * g_yy_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_y_z_0_0_y_z_y_xx, g_y_z_0_0_y_z_y_xy, g_y_z_0_0_y_z_y_xz, g_y_z_0_0_y_z_y_yy, g_y_z_0_0_y_z_y_yz, g_y_z_0_0_y_z_y_zz, g_yy_0_y_xx, g_yy_0_y_xy, g_yy_0_y_xz, g_yy_0_y_yy, g_yy_0_y_yz, g_yy_0_y_zz, g_yy_zz_y_xx, g_yy_zz_y_xy, g_yy_zz_y_xz, g_yy_zz_y_yy, g_yy_zz_y_yz, g_yy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_zz_y_xx[i] * b_exp - 2.0 * g_yy_0_y_xx[i] * a_exp + 4.0 * g_yy_zz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_zz_y_xy[i] * b_exp - 2.0 * g_yy_0_y_xy[i] * a_exp + 4.0 * g_yy_zz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_zz_y_xz[i] * b_exp - 2.0 * g_yy_0_y_xz[i] * a_exp + 4.0 * g_yy_zz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_zz_y_yy[i] * b_exp - 2.0 * g_yy_0_y_yy[i] * a_exp + 4.0 * g_yy_zz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_zz_y_yz[i] * b_exp - 2.0 * g_yy_0_y_yz[i] * a_exp + 4.0 * g_yy_zz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_zz_y_zz[i] * b_exp - 2.0 * g_yy_0_y_zz[i] * a_exp + 4.0 * g_yy_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_y_z_0_0_y_z_z_xx, g_y_z_0_0_y_z_z_xy, g_y_z_0_0_y_z_z_xz, g_y_z_0_0_y_z_z_yy, g_y_z_0_0_y_z_z_yz, g_y_z_0_0_y_z_z_zz, g_yy_0_z_xx, g_yy_0_z_xy, g_yy_0_z_xz, g_yy_0_z_yy, g_yy_0_z_yz, g_yy_0_z_zz, g_yy_zz_z_xx, g_yy_zz_z_xy, g_yy_zz_z_xz, g_yy_zz_z_yy, g_yy_zz_z_yz, g_yy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_zz_z_xx[i] * b_exp - 2.0 * g_yy_0_z_xx[i] * a_exp + 4.0 * g_yy_zz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_zz_z_xy[i] * b_exp - 2.0 * g_yy_0_z_xy[i] * a_exp + 4.0 * g_yy_zz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_zz_z_xz[i] * b_exp - 2.0 * g_yy_0_z_xz[i] * a_exp + 4.0 * g_yy_zz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_zz_z_yy[i] * b_exp - 2.0 * g_yy_0_z_yy[i] * a_exp + 4.0 * g_yy_zz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_zz_z_yz[i] * b_exp - 2.0 * g_yy_0_z_yz[i] * a_exp + 4.0 * g_yy_zz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_zz_z_zz[i] * b_exp - 2.0 * g_yy_0_z_zz[i] * a_exp + 4.0 * g_yy_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_y_z_0_0_z_x_x_xx, g_y_z_0_0_z_x_x_xy, g_y_z_0_0_z_x_x_xz, g_y_z_0_0_z_x_x_yy, g_y_z_0_0_z_x_x_yz, g_y_z_0_0_z_x_x_zz, g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_x_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_y_z_0_0_z_x_y_xx, g_y_z_0_0_z_x_y_xy, g_y_z_0_0_z_x_y_xz, g_y_z_0_0_z_x_y_yy, g_y_z_0_0_z_x_y_yz, g_y_z_0_0_z_x_y_zz, g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_y_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_y_z_0_0_z_x_z_xx, g_y_z_0_0_z_x_z_xy, g_y_z_0_0_z_x_z_xz, g_y_z_0_0_z_x_z_yy, g_y_z_0_0_z_x_z_yz, g_y_z_0_0_z_x_z_zz, g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_z_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_y_z_0_0_z_y_x_xx, g_y_z_0_0_z_y_x_xy, g_y_z_0_0_z_y_x_xz, g_y_z_0_0_z_y_x_yy, g_y_z_0_0_z_y_x_yz, g_y_z_0_0_z_y_x_zz, g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_x_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_y_z_0_0_z_y_y_xx, g_y_z_0_0_z_y_y_xy, g_y_z_0_0_z_y_y_xz, g_y_z_0_0_z_y_y_yy, g_y_z_0_0_z_y_y_yz, g_y_z_0_0_z_y_y_zz, g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_y_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_y_z_0_0_z_y_z_xx, g_y_z_0_0_z_y_z_xy, g_y_z_0_0_z_y_z_xz, g_y_z_0_0_z_y_z_yy, g_y_z_0_0_z_y_z_yz, g_y_z_0_0_z_y_z_zz, g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_z_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_y_z_0_0_z_z_x_xx, g_y_z_0_0_z_z_x_xy, g_y_z_0_0_z_z_x_xz, g_y_z_0_0_z_z_x_yy, g_y_z_0_0_z_z_x_yz, g_y_z_0_0_z_z_x_zz, g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_zz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_zz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_zz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_zz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_zz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_y_z_0_0_z_z_y_xx, g_y_z_0_0_z_z_y_xy, g_y_z_0_0_z_z_y_xz, g_y_z_0_0_z_z_y_yy, g_y_z_0_0_z_z_y_yz, g_y_z_0_0_z_z_y_zz, g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_zz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_zz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_zz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_zz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_zz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_y_z_0_0_z_z_z_xx, g_y_z_0_0_z_z_z_xy, g_y_z_0_0_z_z_z_xz, g_y_z_0_0_z_z_z_yy, g_y_z_0_0_z_z_z_yz, g_y_z_0_0_z_z_z_zz, g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_zz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_zz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_zz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_zz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_zz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz, g_z_x_0_0_x_x_x_xx, g_z_x_0_0_x_x_x_xy, g_z_x_0_0_x_x_x_xz, g_z_x_0_0_x_x_x_yy, g_z_x_0_0_x_x_x_yz, g_z_x_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_xx_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_xx_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_xx_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_xx_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_xx_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz, g_z_x_0_0_x_x_y_xx, g_z_x_0_0_x_x_y_xy, g_z_x_0_0_x_x_y_xz, g_z_x_0_0_x_x_y_yy, g_z_x_0_0_x_x_y_yz, g_z_x_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_xx_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_xx_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_xx_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_xx_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_xx_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz, g_z_x_0_0_x_x_z_xx, g_z_x_0_0_x_x_z_xy, g_z_x_0_0_x_x_z_xz, g_z_x_0_0_x_x_z_yy, g_z_x_0_0_x_x_z_yz, g_z_x_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_xx_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_xx_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_xx_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_xx_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_xx_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz, g_z_x_0_0_x_y_x_xx, g_z_x_0_0_x_y_x_xy, g_z_x_0_0_x_y_x_xz, g_z_x_0_0_x_y_x_yy, g_z_x_0_0_x_y_x_yz, g_z_x_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_x_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz, g_z_x_0_0_x_y_y_xx, g_z_x_0_0_x_y_y_xy, g_z_x_0_0_x_y_y_xz, g_z_x_0_0_x_y_y_yy, g_z_x_0_0_x_y_y_yz, g_z_x_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_y_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz, g_z_x_0_0_x_y_z_xx, g_z_x_0_0_x_y_z_xy, g_z_x_0_0_x_y_z_xz, g_z_x_0_0_x_y_z_yy, g_z_x_0_0_x_y_z_yz, g_z_x_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_z_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz, g_z_x_0_0_x_z_x_xx, g_z_x_0_0_x_z_x_xy, g_z_x_0_0_x_z_x_xz, g_z_x_0_0_x_z_x_yy, g_z_x_0_0_x_z_x_yz, g_z_x_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_x_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz, g_z_x_0_0_x_z_y_xx, g_z_x_0_0_x_z_y_xy, g_z_x_0_0_x_z_y_xz, g_z_x_0_0_x_z_y_yy, g_z_x_0_0_x_z_y_yz, g_z_x_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_y_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz, g_z_x_0_0_x_z_z_xx, g_z_x_0_0_x_z_z_xy, g_z_x_0_0_x_z_z_xz, g_z_x_0_0_x_z_z_yy, g_z_x_0_0_x_z_z_yz, g_z_x_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_z_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz, g_z_x_0_0_y_x_x_xx, g_z_x_0_0_y_x_x_xy, g_z_x_0_0_y_x_x_xz, g_z_x_0_0_y_x_x_yy, g_z_x_0_0_y_x_x_yz, g_z_x_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_xx_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_xx_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_xx_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_xx_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_xx_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz, g_z_x_0_0_y_x_y_xx, g_z_x_0_0_y_x_y_xy, g_z_x_0_0_y_x_y_xz, g_z_x_0_0_y_x_y_yy, g_z_x_0_0_y_x_y_yz, g_z_x_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_xx_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_xx_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_xx_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_xx_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_xx_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz, g_z_x_0_0_y_x_z_xx, g_z_x_0_0_y_x_z_xy, g_z_x_0_0_y_x_z_xz, g_z_x_0_0_y_x_z_yy, g_z_x_0_0_y_x_z_yz, g_z_x_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_xx_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_xx_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_xx_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_xx_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_xx_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz, g_z_x_0_0_y_y_x_xx, g_z_x_0_0_y_y_x_xy, g_z_x_0_0_y_y_x_xz, g_z_x_0_0_y_y_x_yy, g_z_x_0_0_y_y_x_yz, g_z_x_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_x_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz, g_z_x_0_0_y_y_y_xx, g_z_x_0_0_y_y_y_xy, g_z_x_0_0_y_y_y_xz, g_z_x_0_0_y_y_y_yy, g_z_x_0_0_y_y_y_yz, g_z_x_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_y_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz, g_z_x_0_0_y_y_z_xx, g_z_x_0_0_y_y_z_xy, g_z_x_0_0_y_y_z_xz, g_z_x_0_0_y_y_z_yy, g_z_x_0_0_y_y_z_yz, g_z_x_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_z_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz, g_z_x_0_0_y_z_x_xx, g_z_x_0_0_y_z_x_xy, g_z_x_0_0_y_z_x_xz, g_z_x_0_0_y_z_x_yy, g_z_x_0_0_y_z_x_yz, g_z_x_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_x_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz, g_z_x_0_0_y_z_y_xx, g_z_x_0_0_y_z_y_xy, g_z_x_0_0_y_z_y_xz, g_z_x_0_0_y_z_y_yy, g_z_x_0_0_y_z_y_yz, g_z_x_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_y_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz, g_z_x_0_0_y_z_z_xx, g_z_x_0_0_y_z_z_xy, g_z_x_0_0_y_z_z_xz, g_z_x_0_0_y_z_z_yy, g_z_x_0_0_y_z_z_yz, g_z_x_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_z_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_z_x_0_0_z_x_x_xx, g_z_x_0_0_z_x_x_xy, g_z_x_0_0_z_x_x_xz, g_z_x_0_0_z_x_x_yy, g_z_x_0_0_z_x_x_yz, g_z_x_0_0_z_x_x_zz, g_zz_0_x_xx, g_zz_0_x_xy, g_zz_0_x_xz, g_zz_0_x_yy, g_zz_0_x_yz, g_zz_0_x_zz, g_zz_xx_x_xx, g_zz_xx_x_xy, g_zz_xx_x_xz, g_zz_xx_x_yy, g_zz_xx_x_yz, g_zz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_xx_x_xx[i] * b_exp - 2.0 * g_zz_0_x_xx[i] * a_exp + 4.0 * g_zz_xx_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_xx_x_xy[i] * b_exp - 2.0 * g_zz_0_x_xy[i] * a_exp + 4.0 * g_zz_xx_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_xx_x_xz[i] * b_exp - 2.0 * g_zz_0_x_xz[i] * a_exp + 4.0 * g_zz_xx_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_xx_x_yy[i] * b_exp - 2.0 * g_zz_0_x_yy[i] * a_exp + 4.0 * g_zz_xx_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_xx_x_yz[i] * b_exp - 2.0 * g_zz_0_x_yz[i] * a_exp + 4.0 * g_zz_xx_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_xx_x_zz[i] * b_exp - 2.0 * g_zz_0_x_zz[i] * a_exp + 4.0 * g_zz_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_z_x_0_0_z_x_y_xx, g_z_x_0_0_z_x_y_xy, g_z_x_0_0_z_x_y_xz, g_z_x_0_0_z_x_y_yy, g_z_x_0_0_z_x_y_yz, g_z_x_0_0_z_x_y_zz, g_zz_0_y_xx, g_zz_0_y_xy, g_zz_0_y_xz, g_zz_0_y_yy, g_zz_0_y_yz, g_zz_0_y_zz, g_zz_xx_y_xx, g_zz_xx_y_xy, g_zz_xx_y_xz, g_zz_xx_y_yy, g_zz_xx_y_yz, g_zz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_xx_y_xx[i] * b_exp - 2.0 * g_zz_0_y_xx[i] * a_exp + 4.0 * g_zz_xx_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_xx_y_xy[i] * b_exp - 2.0 * g_zz_0_y_xy[i] * a_exp + 4.0 * g_zz_xx_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_xx_y_xz[i] * b_exp - 2.0 * g_zz_0_y_xz[i] * a_exp + 4.0 * g_zz_xx_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_xx_y_yy[i] * b_exp - 2.0 * g_zz_0_y_yy[i] * a_exp + 4.0 * g_zz_xx_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_xx_y_yz[i] * b_exp - 2.0 * g_zz_0_y_yz[i] * a_exp + 4.0 * g_zz_xx_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_xx_y_zz[i] * b_exp - 2.0 * g_zz_0_y_zz[i] * a_exp + 4.0 * g_zz_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_z_x_0_0_z_x_z_xx, g_z_x_0_0_z_x_z_xy, g_z_x_0_0_z_x_z_xz, g_z_x_0_0_z_x_z_yy, g_z_x_0_0_z_x_z_yz, g_z_x_0_0_z_x_z_zz, g_zz_0_z_xx, g_zz_0_z_xy, g_zz_0_z_xz, g_zz_0_z_yy, g_zz_0_z_yz, g_zz_0_z_zz, g_zz_xx_z_xx, g_zz_xx_z_xy, g_zz_xx_z_xz, g_zz_xx_z_yy, g_zz_xx_z_yz, g_zz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_xx_z_xx[i] * b_exp - 2.0 * g_zz_0_z_xx[i] * a_exp + 4.0 * g_zz_xx_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_xx_z_xy[i] * b_exp - 2.0 * g_zz_0_z_xy[i] * a_exp + 4.0 * g_zz_xx_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_xx_z_xz[i] * b_exp - 2.0 * g_zz_0_z_xz[i] * a_exp + 4.0 * g_zz_xx_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_xx_z_yy[i] * b_exp - 2.0 * g_zz_0_z_yy[i] * a_exp + 4.0 * g_zz_xx_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_xx_z_yz[i] * b_exp - 2.0 * g_zz_0_z_yz[i] * a_exp + 4.0 * g_zz_xx_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_xx_z_zz[i] * b_exp - 2.0 * g_zz_0_z_zz[i] * a_exp + 4.0 * g_zz_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_z_x_0_0_z_y_x_xx, g_z_x_0_0_z_y_x_xy, g_z_x_0_0_z_y_x_xz, g_z_x_0_0_z_y_x_yy, g_z_x_0_0_z_y_x_yz, g_z_x_0_0_z_y_x_zz, g_zz_xy_x_xx, g_zz_xy_x_xy, g_zz_xy_x_xz, g_zz_xy_x_yy, g_zz_xy_x_yz, g_zz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_zz_xy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_zz_xy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_zz_xy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_zz_xy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_zz_xy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_zz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_z_x_0_0_z_y_y_xx, g_z_x_0_0_z_y_y_xy, g_z_x_0_0_z_y_y_xz, g_z_x_0_0_z_y_y_yy, g_z_x_0_0_z_y_y_yz, g_z_x_0_0_z_y_y_zz, g_zz_xy_y_xx, g_zz_xy_y_xy, g_zz_xy_y_xz, g_zz_xy_y_yy, g_zz_xy_y_yz, g_zz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_zz_xy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_zz_xy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_zz_xy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_zz_xy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_zz_xy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_zz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_z_x_0_0_z_y_z_xx, g_z_x_0_0_z_y_z_xy, g_z_x_0_0_z_y_z_xz, g_z_x_0_0_z_y_z_yy, g_z_x_0_0_z_y_z_yz, g_z_x_0_0_z_y_z_zz, g_zz_xy_z_xx, g_zz_xy_z_xy, g_zz_xy_z_xz, g_zz_xy_z_yy, g_zz_xy_z_yz, g_zz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_zz_xy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_zz_xy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_zz_xy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_zz_xy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_zz_xy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_zz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_z_x_0_0_z_z_x_xx, g_z_x_0_0_z_z_x_xy, g_z_x_0_0_z_z_x_xz, g_z_x_0_0_z_z_x_yy, g_z_x_0_0_z_z_x_yz, g_z_x_0_0_z_z_x_zz, g_zz_xz_x_xx, g_zz_xz_x_xy, g_zz_xz_x_xz, g_zz_xz_x_yy, g_zz_xz_x_yz, g_zz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_zz_xz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_zz_xz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_zz_xz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_zz_xz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_zz_xz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_zz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_z_x_0_0_z_z_y_xx, g_z_x_0_0_z_z_y_xy, g_z_x_0_0_z_z_y_xz, g_z_x_0_0_z_z_y_yy, g_z_x_0_0_z_z_y_yz, g_z_x_0_0_z_z_y_zz, g_zz_xz_y_xx, g_zz_xz_y_xy, g_zz_xz_y_xz, g_zz_xz_y_yy, g_zz_xz_y_yz, g_zz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_zz_xz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_zz_xz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_zz_xz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_zz_xz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_zz_xz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_zz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_z_x_0_0_z_z_z_xx, g_z_x_0_0_z_z_z_xy, g_z_x_0_0_z_z_z_xz, g_z_x_0_0_z_z_z_yy, g_z_x_0_0_z_z_z_yz, g_z_x_0_0_z_z_z_zz, g_zz_xz_z_xx, g_zz_xz_z_xy, g_zz_xz_z_xz, g_zz_xz_z_yy, g_zz_xz_z_yz, g_zz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_zz_xz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_zz_xz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_zz_xz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_zz_xz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_zz_xz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_zz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz, g_z_y_0_0_x_x_x_xx, g_z_y_0_0_x_x_x_xy, g_z_y_0_0_x_x_x_xz, g_z_y_0_0_x_x_x_yy, g_z_y_0_0_x_x_x_yz, g_z_y_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_x_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz, g_z_y_0_0_x_x_y_xx, g_z_y_0_0_x_x_y_xy, g_z_y_0_0_x_x_y_xz, g_z_y_0_0_x_x_y_yy, g_z_y_0_0_x_x_y_yz, g_z_y_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_y_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz, g_z_y_0_0_x_x_z_xx, g_z_y_0_0_x_x_z_xy, g_z_y_0_0_x_x_z_xz, g_z_y_0_0_x_x_z_yy, g_z_y_0_0_x_x_z_yz, g_z_y_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_z_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz, g_z_y_0_0_x_y_x_xx, g_z_y_0_0_x_y_x_xy, g_z_y_0_0_x_y_x_xz, g_z_y_0_0_x_y_x_yy, g_z_y_0_0_x_y_x_yz, g_z_y_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_yy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_yy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_yy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_yy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_yy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz, g_z_y_0_0_x_y_y_xx, g_z_y_0_0_x_y_y_xy, g_z_y_0_0_x_y_y_xz, g_z_y_0_0_x_y_y_yy, g_z_y_0_0_x_y_y_yz, g_z_y_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_yy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_yy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_yy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_yy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_yy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz, g_z_y_0_0_x_y_z_xx, g_z_y_0_0_x_y_z_xy, g_z_y_0_0_x_y_z_xz, g_z_y_0_0_x_y_z_yy, g_z_y_0_0_x_y_z_yz, g_z_y_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_yy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_yy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_yy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_yy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_yy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz, g_z_y_0_0_x_z_x_xx, g_z_y_0_0_x_z_x_xy, g_z_y_0_0_x_z_x_xz, g_z_y_0_0_x_z_x_yy, g_z_y_0_0_x_z_x_yz, g_z_y_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_x_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz, g_z_y_0_0_x_z_y_xx, g_z_y_0_0_x_z_y_xy, g_z_y_0_0_x_z_y_xz, g_z_y_0_0_x_z_y_yy, g_z_y_0_0_x_z_y_yz, g_z_y_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_y_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz, g_z_y_0_0_x_z_z_xx, g_z_y_0_0_x_z_z_xy, g_z_y_0_0_x_z_z_xz, g_z_y_0_0_x_z_z_yy, g_z_y_0_0_x_z_z_yz, g_z_y_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_z_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz, g_z_y_0_0_y_x_x_xx, g_z_y_0_0_y_x_x_xy, g_z_y_0_0_y_x_x_xz, g_z_y_0_0_y_x_x_yy, g_z_y_0_0_y_x_x_yz, g_z_y_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_x_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz, g_z_y_0_0_y_x_y_xx, g_z_y_0_0_y_x_y_xy, g_z_y_0_0_y_x_y_xz, g_z_y_0_0_y_x_y_yy, g_z_y_0_0_y_x_y_yz, g_z_y_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_y_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz, g_z_y_0_0_y_x_z_xx, g_z_y_0_0_y_x_z_xy, g_z_y_0_0_y_x_z_xz, g_z_y_0_0_y_x_z_yy, g_z_y_0_0_y_x_z_yz, g_z_y_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_z_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz, g_z_y_0_0_y_y_x_xx, g_z_y_0_0_y_y_x_xy, g_z_y_0_0_y_y_x_xz, g_z_y_0_0_y_y_x_yy, g_z_y_0_0_y_y_x_yz, g_z_y_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_yy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_yy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_yy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_yy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_yy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz, g_z_y_0_0_y_y_y_xx, g_z_y_0_0_y_y_y_xy, g_z_y_0_0_y_y_y_xz, g_z_y_0_0_y_y_y_yy, g_z_y_0_0_y_y_y_yz, g_z_y_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_yy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_yy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_yy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_yy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_yy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz, g_z_y_0_0_y_y_z_xx, g_z_y_0_0_y_y_z_xy, g_z_y_0_0_y_y_z_xz, g_z_y_0_0_y_y_z_yy, g_z_y_0_0_y_y_z_yz, g_z_y_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_yy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_yy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_yy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_yy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_yy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz, g_z_y_0_0_y_z_x_xx, g_z_y_0_0_y_z_x_xy, g_z_y_0_0_y_z_x_xz, g_z_y_0_0_y_z_x_yy, g_z_y_0_0_y_z_x_yz, g_z_y_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_x_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz, g_z_y_0_0_y_z_y_xx, g_z_y_0_0_y_z_y_xy, g_z_y_0_0_y_z_y_xz, g_z_y_0_0_y_z_y_yy, g_z_y_0_0_y_z_y_yz, g_z_y_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_y_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz, g_z_y_0_0_y_z_z_xx, g_z_y_0_0_y_z_z_xy, g_z_y_0_0_y_z_z_xz, g_z_y_0_0_y_z_z_yy, g_z_y_0_0_y_z_z_yz, g_z_y_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_z_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_z_y_0_0_z_x_x_xx, g_z_y_0_0_z_x_x_xy, g_z_y_0_0_z_x_x_xz, g_z_y_0_0_z_x_x_yy, g_z_y_0_0_z_x_x_yz, g_z_y_0_0_z_x_x_zz, g_zz_xy_x_xx, g_zz_xy_x_xy, g_zz_xy_x_xz, g_zz_xy_x_yy, g_zz_xy_x_yz, g_zz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * b_exp + 4.0 * g_zz_xy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * b_exp + 4.0 * g_zz_xy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * b_exp + 4.0 * g_zz_xy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * b_exp + 4.0 * g_zz_xy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * b_exp + 4.0 * g_zz_xy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * b_exp + 4.0 * g_zz_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_z_y_0_0_z_x_y_xx, g_z_y_0_0_z_x_y_xy, g_z_y_0_0_z_x_y_xz, g_z_y_0_0_z_x_y_yy, g_z_y_0_0_z_x_y_yz, g_z_y_0_0_z_x_y_zz, g_zz_xy_y_xx, g_zz_xy_y_xy, g_zz_xy_y_xz, g_zz_xy_y_yy, g_zz_xy_y_yz, g_zz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * b_exp + 4.0 * g_zz_xy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * b_exp + 4.0 * g_zz_xy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * b_exp + 4.0 * g_zz_xy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * b_exp + 4.0 * g_zz_xy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * b_exp + 4.0 * g_zz_xy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * b_exp + 4.0 * g_zz_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_z_y_0_0_z_x_z_xx, g_z_y_0_0_z_x_z_xy, g_z_y_0_0_z_x_z_xz, g_z_y_0_0_z_x_z_yy, g_z_y_0_0_z_x_z_yz, g_z_y_0_0_z_x_z_zz, g_zz_xy_z_xx, g_zz_xy_z_xy, g_zz_xy_z_xz, g_zz_xy_z_yy, g_zz_xy_z_yz, g_zz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * b_exp + 4.0 * g_zz_xy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * b_exp + 4.0 * g_zz_xy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * b_exp + 4.0 * g_zz_xy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * b_exp + 4.0 * g_zz_xy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * b_exp + 4.0 * g_zz_xy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * b_exp + 4.0 * g_zz_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_z_y_0_0_z_y_x_xx, g_z_y_0_0_z_y_x_xy, g_z_y_0_0_z_y_x_xz, g_z_y_0_0_z_y_x_yy, g_z_y_0_0_z_y_x_yz, g_z_y_0_0_z_y_x_zz, g_zz_0_x_xx, g_zz_0_x_xy, g_zz_0_x_xz, g_zz_0_x_yy, g_zz_0_x_yz, g_zz_0_x_zz, g_zz_yy_x_xx, g_zz_yy_x_xy, g_zz_yy_x_xz, g_zz_yy_x_yy, g_zz_yy_x_yz, g_zz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_yy_x_xx[i] * b_exp - 2.0 * g_zz_0_x_xx[i] * a_exp + 4.0 * g_zz_yy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_yy_x_xy[i] * b_exp - 2.0 * g_zz_0_x_xy[i] * a_exp + 4.0 * g_zz_yy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_yy_x_xz[i] * b_exp - 2.0 * g_zz_0_x_xz[i] * a_exp + 4.0 * g_zz_yy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_yy_x_yy[i] * b_exp - 2.0 * g_zz_0_x_yy[i] * a_exp + 4.0 * g_zz_yy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_yy_x_yz[i] * b_exp - 2.0 * g_zz_0_x_yz[i] * a_exp + 4.0 * g_zz_yy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_yy_x_zz[i] * b_exp - 2.0 * g_zz_0_x_zz[i] * a_exp + 4.0 * g_zz_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_z_y_0_0_z_y_y_xx, g_z_y_0_0_z_y_y_xy, g_z_y_0_0_z_y_y_xz, g_z_y_0_0_z_y_y_yy, g_z_y_0_0_z_y_y_yz, g_z_y_0_0_z_y_y_zz, g_zz_0_y_xx, g_zz_0_y_xy, g_zz_0_y_xz, g_zz_0_y_yy, g_zz_0_y_yz, g_zz_0_y_zz, g_zz_yy_y_xx, g_zz_yy_y_xy, g_zz_yy_y_xz, g_zz_yy_y_yy, g_zz_yy_y_yz, g_zz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_yy_y_xx[i] * b_exp - 2.0 * g_zz_0_y_xx[i] * a_exp + 4.0 * g_zz_yy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_yy_y_xy[i] * b_exp - 2.0 * g_zz_0_y_xy[i] * a_exp + 4.0 * g_zz_yy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_yy_y_xz[i] * b_exp - 2.0 * g_zz_0_y_xz[i] * a_exp + 4.0 * g_zz_yy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_yy_y_yy[i] * b_exp - 2.0 * g_zz_0_y_yy[i] * a_exp + 4.0 * g_zz_yy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_yy_y_yz[i] * b_exp - 2.0 * g_zz_0_y_yz[i] * a_exp + 4.0 * g_zz_yy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_yy_y_zz[i] * b_exp - 2.0 * g_zz_0_y_zz[i] * a_exp + 4.0 * g_zz_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_z_y_0_0_z_y_z_xx, g_z_y_0_0_z_y_z_xy, g_z_y_0_0_z_y_z_xz, g_z_y_0_0_z_y_z_yy, g_z_y_0_0_z_y_z_yz, g_z_y_0_0_z_y_z_zz, g_zz_0_z_xx, g_zz_0_z_xy, g_zz_0_z_xz, g_zz_0_z_yy, g_zz_0_z_yz, g_zz_0_z_zz, g_zz_yy_z_xx, g_zz_yy_z_xy, g_zz_yy_z_xz, g_zz_yy_z_yy, g_zz_yy_z_yz, g_zz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_yy_z_xx[i] * b_exp - 2.0 * g_zz_0_z_xx[i] * a_exp + 4.0 * g_zz_yy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_yy_z_xy[i] * b_exp - 2.0 * g_zz_0_z_xy[i] * a_exp + 4.0 * g_zz_yy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_yy_z_xz[i] * b_exp - 2.0 * g_zz_0_z_xz[i] * a_exp + 4.0 * g_zz_yy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_yy_z_yy[i] * b_exp - 2.0 * g_zz_0_z_yy[i] * a_exp + 4.0 * g_zz_yy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_yy_z_yz[i] * b_exp - 2.0 * g_zz_0_z_yz[i] * a_exp + 4.0 * g_zz_yy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_yy_z_zz[i] * b_exp - 2.0 * g_zz_0_z_zz[i] * a_exp + 4.0 * g_zz_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_z_y_0_0_z_z_x_xx, g_z_y_0_0_z_z_x_xy, g_z_y_0_0_z_z_x_xz, g_z_y_0_0_z_z_x_yy, g_z_y_0_0_z_z_x_yz, g_z_y_0_0_z_z_x_zz, g_zz_yz_x_xx, g_zz_yz_x_xy, g_zz_yz_x_xz, g_zz_yz_x_yy, g_zz_yz_x_yz, g_zz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_zz_yz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_zz_yz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_zz_yz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_zz_yz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_zz_yz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_zz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_z_y_0_0_z_z_y_xx, g_z_y_0_0_z_z_y_xy, g_z_y_0_0_z_z_y_xz, g_z_y_0_0_z_z_y_yy, g_z_y_0_0_z_z_y_yz, g_z_y_0_0_z_z_y_zz, g_zz_yz_y_xx, g_zz_yz_y_xy, g_zz_yz_y_xz, g_zz_yz_y_yy, g_zz_yz_y_yz, g_zz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_zz_yz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_zz_yz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_zz_yz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_zz_yz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_zz_yz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_zz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_z_y_0_0_z_z_z_xx, g_z_y_0_0_z_z_z_xy, g_z_y_0_0_z_z_z_xz, g_z_y_0_0_z_z_z_yy, g_z_y_0_0_z_z_z_yz, g_z_y_0_0_z_z_z_zz, g_zz_yz_z_xx, g_zz_yz_z_xy, g_zz_yz_z_xz, g_zz_yz_z_yy, g_zz_yz_z_yz, g_zz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_zz_yz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_zz_yz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_zz_yz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_zz_yz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_zz_yz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_zz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz, g_z_z_0_0_x_x_x_xx, g_z_z_0_0_x_x_x_xy, g_z_z_0_0_x_x_x_xz, g_z_z_0_0_x_x_x_yy, g_z_z_0_0_x_x_x_yz, g_z_z_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_x_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz, g_z_z_0_0_x_x_y_xx, g_z_z_0_0_x_x_y_xy, g_z_z_0_0_x_x_y_xz, g_z_z_0_0_x_x_y_yy, g_z_z_0_0_x_x_y_yz, g_z_z_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_y_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz, g_z_z_0_0_x_x_z_xx, g_z_z_0_0_x_x_z_xy, g_z_z_0_0_x_x_z_xz, g_z_z_0_0_x_x_z_yy, g_z_z_0_0_x_x_z_yz, g_z_z_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_z_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz, g_z_z_0_0_x_y_x_xx, g_z_z_0_0_x_y_x_xy, g_z_z_0_0_x_y_x_xz, g_z_z_0_0_x_y_x_yy, g_z_z_0_0_x_y_x_yz, g_z_z_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_x_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz, g_z_z_0_0_x_y_y_xx, g_z_z_0_0_x_y_y_xy, g_z_z_0_0_x_y_y_xz, g_z_z_0_0_x_y_y_yy, g_z_z_0_0_x_y_y_yz, g_z_z_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_y_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz, g_z_z_0_0_x_y_z_xx, g_z_z_0_0_x_y_z_xy, g_z_z_0_0_x_y_z_xz, g_z_z_0_0_x_y_z_yy, g_z_z_0_0_x_y_z_yz, g_z_z_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_z_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz, g_z_z_0_0_x_z_x_xx, g_z_z_0_0_x_z_x_xy, g_z_z_0_0_x_z_x_xz, g_z_z_0_0_x_z_x_yy, g_z_z_0_0_x_z_x_yz, g_z_z_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_x_xx[i] = -2.0 * g_xz_0_x_xx[i] * a_exp + 4.0 * g_xz_zz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_xy[i] = -2.0 * g_xz_0_x_xy[i] * a_exp + 4.0 * g_xz_zz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_xz[i] = -2.0 * g_xz_0_x_xz[i] * a_exp + 4.0 * g_xz_zz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_yy[i] = -2.0 * g_xz_0_x_yy[i] * a_exp + 4.0 * g_xz_zz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_yz[i] = -2.0 * g_xz_0_x_yz[i] * a_exp + 4.0 * g_xz_zz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_zz[i] = -2.0 * g_xz_0_x_zz[i] * a_exp + 4.0 * g_xz_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz, g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz, g_z_z_0_0_x_z_y_xx, g_z_z_0_0_x_z_y_xy, g_z_z_0_0_x_z_y_xz, g_z_z_0_0_x_z_y_yy, g_z_z_0_0_x_z_y_yz, g_z_z_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_y_xx[i] = -2.0 * g_xz_0_y_xx[i] * a_exp + 4.0 * g_xz_zz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_xy[i] = -2.0 * g_xz_0_y_xy[i] * a_exp + 4.0 * g_xz_zz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_xz[i] = -2.0 * g_xz_0_y_xz[i] * a_exp + 4.0 * g_xz_zz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_yy[i] = -2.0 * g_xz_0_y_yy[i] * a_exp + 4.0 * g_xz_zz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_yz[i] = -2.0 * g_xz_0_y_yz[i] * a_exp + 4.0 * g_xz_zz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_zz[i] = -2.0 * g_xz_0_y_zz[i] * a_exp + 4.0 * g_xz_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz, g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz, g_z_z_0_0_x_z_z_xx, g_z_z_0_0_x_z_z_xy, g_z_z_0_0_x_z_z_xz, g_z_z_0_0_x_z_z_yy, g_z_z_0_0_x_z_z_yz, g_z_z_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_z_xx[i] = -2.0 * g_xz_0_z_xx[i] * a_exp + 4.0 * g_xz_zz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_xy[i] = -2.0 * g_xz_0_z_xy[i] * a_exp + 4.0 * g_xz_zz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_xz[i] = -2.0 * g_xz_0_z_xz[i] * a_exp + 4.0 * g_xz_zz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_yy[i] = -2.0 * g_xz_0_z_yy[i] * a_exp + 4.0 * g_xz_zz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_yz[i] = -2.0 * g_xz_0_z_yz[i] * a_exp + 4.0 * g_xz_zz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_zz[i] = -2.0 * g_xz_0_z_zz[i] * a_exp + 4.0 * g_xz_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz, g_z_z_0_0_y_x_x_xx, g_z_z_0_0_y_x_x_xy, g_z_z_0_0_y_x_x_xz, g_z_z_0_0_y_x_x_yy, g_z_z_0_0_y_x_x_yz, g_z_z_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_x_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz, g_z_z_0_0_y_x_y_xx, g_z_z_0_0_y_x_y_xy, g_z_z_0_0_y_x_y_xz, g_z_z_0_0_y_x_y_yy, g_z_z_0_0_y_x_y_yz, g_z_z_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_y_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz, g_z_z_0_0_y_x_z_xx, g_z_z_0_0_y_x_z_xy, g_z_z_0_0_y_x_z_xz, g_z_z_0_0_y_x_z_yy, g_z_z_0_0_y_x_z_yz, g_z_z_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_z_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz, g_z_z_0_0_y_y_x_xx, g_z_z_0_0_y_y_x_xy, g_z_z_0_0_y_y_x_xz, g_z_z_0_0_y_y_x_yy, g_z_z_0_0_y_y_x_yz, g_z_z_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_x_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz, g_z_z_0_0_y_y_y_xx, g_z_z_0_0_y_y_y_xy, g_z_z_0_0_y_y_y_xz, g_z_z_0_0_y_y_y_yy, g_z_z_0_0_y_y_y_yz, g_z_z_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_y_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz, g_z_z_0_0_y_y_z_xx, g_z_z_0_0_y_y_z_xy, g_z_z_0_0_y_y_z_xz, g_z_z_0_0_y_y_z_yy, g_z_z_0_0_y_y_z_yz, g_z_z_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_z_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz, g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz, g_z_z_0_0_y_z_x_xx, g_z_z_0_0_y_z_x_xy, g_z_z_0_0_y_z_x_xz, g_z_z_0_0_y_z_x_yy, g_z_z_0_0_y_z_x_yz, g_z_z_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_x_xx[i] = -2.0 * g_yz_0_x_xx[i] * a_exp + 4.0 * g_yz_zz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_xy[i] = -2.0 * g_yz_0_x_xy[i] * a_exp + 4.0 * g_yz_zz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_xz[i] = -2.0 * g_yz_0_x_xz[i] * a_exp + 4.0 * g_yz_zz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_yy[i] = -2.0 * g_yz_0_x_yy[i] * a_exp + 4.0 * g_yz_zz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_yz[i] = -2.0 * g_yz_0_x_yz[i] * a_exp + 4.0 * g_yz_zz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_zz[i] = -2.0 * g_yz_0_x_zz[i] * a_exp + 4.0 * g_yz_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz, g_z_z_0_0_y_z_y_xx, g_z_z_0_0_y_z_y_xy, g_z_z_0_0_y_z_y_xz, g_z_z_0_0_y_z_y_yy, g_z_z_0_0_y_z_y_yz, g_z_z_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_y_xx[i] = -2.0 * g_yz_0_y_xx[i] * a_exp + 4.0 * g_yz_zz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_xy[i] = -2.0 * g_yz_0_y_xy[i] * a_exp + 4.0 * g_yz_zz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_xz[i] = -2.0 * g_yz_0_y_xz[i] * a_exp + 4.0 * g_yz_zz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_yy[i] = -2.0 * g_yz_0_y_yy[i] * a_exp + 4.0 * g_yz_zz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_yz[i] = -2.0 * g_yz_0_y_yz[i] * a_exp + 4.0 * g_yz_zz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_zz[i] = -2.0 * g_yz_0_y_zz[i] * a_exp + 4.0 * g_yz_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz, g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz, g_z_z_0_0_y_z_z_xx, g_z_z_0_0_y_z_z_xy, g_z_z_0_0_y_z_z_xz, g_z_z_0_0_y_z_z_yy, g_z_z_0_0_y_z_z_yz, g_z_z_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_z_xx[i] = -2.0 * g_yz_0_z_xx[i] * a_exp + 4.0 * g_yz_zz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_xy[i] = -2.0 * g_yz_0_z_xy[i] * a_exp + 4.0 * g_yz_zz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_xz[i] = -2.0 * g_yz_0_z_xz[i] * a_exp + 4.0 * g_yz_zz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_yy[i] = -2.0 * g_yz_0_z_yy[i] * a_exp + 4.0 * g_yz_zz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_yz[i] = -2.0 * g_yz_0_z_yz[i] * a_exp + 4.0 * g_yz_zz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_zz[i] = -2.0 * g_yz_0_z_zz[i] * a_exp + 4.0 * g_yz_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_z_z_0_0_z_x_x_xx, g_z_z_0_0_z_x_x_xy, g_z_z_0_0_z_x_x_xz, g_z_z_0_0_z_x_x_yy, g_z_z_0_0_z_x_x_yz, g_z_z_0_0_z_x_x_zz, g_zz_xz_x_xx, g_zz_xz_x_xy, g_zz_xz_x_xz, g_zz_xz_x_yy, g_zz_xz_x_yz, g_zz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * b_exp + 4.0 * g_zz_xz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * b_exp + 4.0 * g_zz_xz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * b_exp + 4.0 * g_zz_xz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * b_exp + 4.0 * g_zz_xz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * b_exp + 4.0 * g_zz_xz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * b_exp + 4.0 * g_zz_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_z_z_0_0_z_x_y_xx, g_z_z_0_0_z_x_y_xy, g_z_z_0_0_z_x_y_xz, g_z_z_0_0_z_x_y_yy, g_z_z_0_0_z_x_y_yz, g_z_z_0_0_z_x_y_zz, g_zz_xz_y_xx, g_zz_xz_y_xy, g_zz_xz_y_xz, g_zz_xz_y_yy, g_zz_xz_y_yz, g_zz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * b_exp + 4.0 * g_zz_xz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * b_exp + 4.0 * g_zz_xz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * b_exp + 4.0 * g_zz_xz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * b_exp + 4.0 * g_zz_xz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * b_exp + 4.0 * g_zz_xz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * b_exp + 4.0 * g_zz_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_z_z_0_0_z_x_z_xx, g_z_z_0_0_z_x_z_xy, g_z_z_0_0_z_x_z_xz, g_z_z_0_0_z_x_z_yy, g_z_z_0_0_z_x_z_yz, g_z_z_0_0_z_x_z_zz, g_zz_xz_z_xx, g_zz_xz_z_xy, g_zz_xz_z_xz, g_zz_xz_z_yy, g_zz_xz_z_yz, g_zz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * b_exp + 4.0 * g_zz_xz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * b_exp + 4.0 * g_zz_xz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * b_exp + 4.0 * g_zz_xz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * b_exp + 4.0 * g_zz_xz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * b_exp + 4.0 * g_zz_xz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * b_exp + 4.0 * g_zz_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_z_z_0_0_z_y_x_xx, g_z_z_0_0_z_y_x_xy, g_z_z_0_0_z_y_x_xz, g_z_z_0_0_z_y_x_yy, g_z_z_0_0_z_y_x_yz, g_z_z_0_0_z_y_x_zz, g_zz_yz_x_xx, g_zz_yz_x_xy, g_zz_yz_x_xz, g_zz_yz_x_yy, g_zz_yz_x_yz, g_zz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * b_exp + 4.0 * g_zz_yz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * b_exp + 4.0 * g_zz_yz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * b_exp + 4.0 * g_zz_yz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * b_exp + 4.0 * g_zz_yz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * b_exp + 4.0 * g_zz_yz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * b_exp + 4.0 * g_zz_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_z_z_0_0_z_y_y_xx, g_z_z_0_0_z_y_y_xy, g_z_z_0_0_z_y_y_xz, g_z_z_0_0_z_y_y_yy, g_z_z_0_0_z_y_y_yz, g_z_z_0_0_z_y_y_zz, g_zz_yz_y_xx, g_zz_yz_y_xy, g_zz_yz_y_xz, g_zz_yz_y_yy, g_zz_yz_y_yz, g_zz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * b_exp + 4.0 * g_zz_yz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * b_exp + 4.0 * g_zz_yz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * b_exp + 4.0 * g_zz_yz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * b_exp + 4.0 * g_zz_yz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * b_exp + 4.0 * g_zz_yz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * b_exp + 4.0 * g_zz_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_z_z_0_0_z_y_z_xx, g_z_z_0_0_z_y_z_xy, g_z_z_0_0_z_y_z_xz, g_z_z_0_0_z_y_z_yy, g_z_z_0_0_z_y_z_yz, g_z_z_0_0_z_y_z_zz, g_zz_yz_z_xx, g_zz_yz_z_xy, g_zz_yz_z_xz, g_zz_yz_z_yy, g_zz_yz_z_yz, g_zz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * b_exp + 4.0 * g_zz_yz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * b_exp + 4.0 * g_zz_yz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * b_exp + 4.0 * g_zz_yz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * b_exp + 4.0 * g_zz_yz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * b_exp + 4.0 * g_zz_yz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * b_exp + 4.0 * g_zz_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_z_z_0_0_z_z_x_xx, g_z_z_0_0_z_z_x_xy, g_z_z_0_0_z_z_x_xz, g_z_z_0_0_z_z_x_yy, g_z_z_0_0_z_z_x_yz, g_z_z_0_0_z_z_x_zz, g_zz_0_x_xx, g_zz_0_x_xy, g_zz_0_x_xz, g_zz_0_x_yy, g_zz_0_x_yz, g_zz_0_x_zz, g_zz_zz_x_xx, g_zz_zz_x_xy, g_zz_zz_x_xz, g_zz_zz_x_yy, g_zz_zz_x_yz, g_zz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_x_xx[i] = g_0_0_x_xx[i] - 2.0 * g_0_zz_x_xx[i] * b_exp - 2.0 * g_zz_0_x_xx[i] * a_exp + 4.0 * g_zz_zz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_xy[i] = g_0_0_x_xy[i] - 2.0 * g_0_zz_x_xy[i] * b_exp - 2.0 * g_zz_0_x_xy[i] * a_exp + 4.0 * g_zz_zz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_xz[i] = g_0_0_x_xz[i] - 2.0 * g_0_zz_x_xz[i] * b_exp - 2.0 * g_zz_0_x_xz[i] * a_exp + 4.0 * g_zz_zz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_yy[i] = g_0_0_x_yy[i] - 2.0 * g_0_zz_x_yy[i] * b_exp - 2.0 * g_zz_0_x_yy[i] * a_exp + 4.0 * g_zz_zz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_yz[i] = g_0_0_x_yz[i] - 2.0 * g_0_zz_x_yz[i] * b_exp - 2.0 * g_zz_0_x_yz[i] * a_exp + 4.0 * g_zz_zz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_zz[i] = g_0_0_x_zz[i] - 2.0 * g_0_zz_x_zz[i] * b_exp - 2.0 * g_zz_0_x_zz[i] * a_exp + 4.0 * g_zz_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_z_z_0_0_z_z_y_xx, g_z_z_0_0_z_z_y_xy, g_z_z_0_0_z_z_y_xz, g_z_z_0_0_z_z_y_yy, g_z_z_0_0_z_z_y_yz, g_z_z_0_0_z_z_y_zz, g_zz_0_y_xx, g_zz_0_y_xy, g_zz_0_y_xz, g_zz_0_y_yy, g_zz_0_y_yz, g_zz_0_y_zz, g_zz_zz_y_xx, g_zz_zz_y_xy, g_zz_zz_y_xz, g_zz_zz_y_yy, g_zz_zz_y_yz, g_zz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_y_xx[i] = g_0_0_y_xx[i] - 2.0 * g_0_zz_y_xx[i] * b_exp - 2.0 * g_zz_0_y_xx[i] * a_exp + 4.0 * g_zz_zz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_xy[i] = g_0_0_y_xy[i] - 2.0 * g_0_zz_y_xy[i] * b_exp - 2.0 * g_zz_0_y_xy[i] * a_exp + 4.0 * g_zz_zz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_xz[i] = g_0_0_y_xz[i] - 2.0 * g_0_zz_y_xz[i] * b_exp - 2.0 * g_zz_0_y_xz[i] * a_exp + 4.0 * g_zz_zz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_yy[i] = g_0_0_y_yy[i] - 2.0 * g_0_zz_y_yy[i] * b_exp - 2.0 * g_zz_0_y_yy[i] * a_exp + 4.0 * g_zz_zz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_yz[i] = g_0_0_y_yz[i] - 2.0 * g_0_zz_y_yz[i] * b_exp - 2.0 * g_zz_0_y_yz[i] * a_exp + 4.0 * g_zz_zz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_zz[i] = g_0_0_y_zz[i] - 2.0 * g_0_zz_y_zz[i] * b_exp - 2.0 * g_zz_0_y_zz[i] * a_exp + 4.0 * g_zz_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_z_z_0_0_z_z_z_xx, g_z_z_0_0_z_z_z_xy, g_z_z_0_0_z_z_z_xz, g_z_z_0_0_z_z_z_yy, g_z_z_0_0_z_z_z_yz, g_z_z_0_0_z_z_z_zz, g_zz_0_z_xx, g_zz_0_z_xy, g_zz_0_z_xz, g_zz_0_z_yy, g_zz_0_z_yz, g_zz_0_z_zz, g_zz_zz_z_xx, g_zz_zz_z_xy, g_zz_zz_z_xz, g_zz_zz_z_yy, g_zz_zz_z_yz, g_zz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_z_xx[i] = g_0_0_z_xx[i] - 2.0 * g_0_zz_z_xx[i] * b_exp - 2.0 * g_zz_0_z_xx[i] * a_exp + 4.0 * g_zz_zz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_xy[i] = g_0_0_z_xy[i] - 2.0 * g_0_zz_z_xy[i] * b_exp - 2.0 * g_zz_0_z_xy[i] * a_exp + 4.0 * g_zz_zz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_xz[i] = g_0_0_z_xz[i] - 2.0 * g_0_zz_z_xz[i] * b_exp - 2.0 * g_zz_0_z_xz[i] * a_exp + 4.0 * g_zz_zz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_yy[i] = g_0_0_z_yy[i] - 2.0 * g_0_zz_z_yy[i] * b_exp - 2.0 * g_zz_0_z_yy[i] * a_exp + 4.0 * g_zz_zz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_yz[i] = g_0_0_z_yz[i] - 2.0 * g_0_zz_z_yz[i] * b_exp - 2.0 * g_zz_0_z_yz[i] * a_exp + 4.0 * g_zz_zz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_zz[i] = g_0_0_z_zz[i] - 2.0 * g_0_zz_z_zz[i] * b_exp - 2.0 * g_zz_0_z_zz[i] * a_exp + 4.0 * g_zz_zz_z_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


#include "GeomDeriv1100OfScalarForPDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pdpp_0(CSimdArray<double>& buffer_1100_pdpp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_sfpp,
                     const CSimdArray<double>& buffer_dppp,
                     const CSimdArray<double>& buffer_dfpp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pdpp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sppp

    auto g_0_x_x_x = buffer_sppp[0];

    auto g_0_x_x_y = buffer_sppp[1];

    auto g_0_x_x_z = buffer_sppp[2];

    auto g_0_x_y_x = buffer_sppp[3];

    auto g_0_x_y_y = buffer_sppp[4];

    auto g_0_x_y_z = buffer_sppp[5];

    auto g_0_x_z_x = buffer_sppp[6];

    auto g_0_x_z_y = buffer_sppp[7];

    auto g_0_x_z_z = buffer_sppp[8];

    auto g_0_y_x_x = buffer_sppp[9];

    auto g_0_y_x_y = buffer_sppp[10];

    auto g_0_y_x_z = buffer_sppp[11];

    auto g_0_y_y_x = buffer_sppp[12];

    auto g_0_y_y_y = buffer_sppp[13];

    auto g_0_y_y_z = buffer_sppp[14];

    auto g_0_y_z_x = buffer_sppp[15];

    auto g_0_y_z_y = buffer_sppp[16];

    auto g_0_y_z_z = buffer_sppp[17];

    auto g_0_z_x_x = buffer_sppp[18];

    auto g_0_z_x_y = buffer_sppp[19];

    auto g_0_z_x_z = buffer_sppp[20];

    auto g_0_z_y_x = buffer_sppp[21];

    auto g_0_z_y_y = buffer_sppp[22];

    auto g_0_z_y_z = buffer_sppp[23];

    auto g_0_z_z_x = buffer_sppp[24];

    auto g_0_z_z_y = buffer_sppp[25];

    auto g_0_z_z_z = buffer_sppp[26];

    /// Set up components of auxilary buffer : buffer_sfpp

    auto g_0_xxx_x_x = buffer_sfpp[0];

    auto g_0_xxx_x_y = buffer_sfpp[1];

    auto g_0_xxx_x_z = buffer_sfpp[2];

    auto g_0_xxx_y_x = buffer_sfpp[3];

    auto g_0_xxx_y_y = buffer_sfpp[4];

    auto g_0_xxx_y_z = buffer_sfpp[5];

    auto g_0_xxx_z_x = buffer_sfpp[6];

    auto g_0_xxx_z_y = buffer_sfpp[7];

    auto g_0_xxx_z_z = buffer_sfpp[8];

    auto g_0_xxy_x_x = buffer_sfpp[9];

    auto g_0_xxy_x_y = buffer_sfpp[10];

    auto g_0_xxy_x_z = buffer_sfpp[11];

    auto g_0_xxy_y_x = buffer_sfpp[12];

    auto g_0_xxy_y_y = buffer_sfpp[13];

    auto g_0_xxy_y_z = buffer_sfpp[14];

    auto g_0_xxy_z_x = buffer_sfpp[15];

    auto g_0_xxy_z_y = buffer_sfpp[16];

    auto g_0_xxy_z_z = buffer_sfpp[17];

    auto g_0_xxz_x_x = buffer_sfpp[18];

    auto g_0_xxz_x_y = buffer_sfpp[19];

    auto g_0_xxz_x_z = buffer_sfpp[20];

    auto g_0_xxz_y_x = buffer_sfpp[21];

    auto g_0_xxz_y_y = buffer_sfpp[22];

    auto g_0_xxz_y_z = buffer_sfpp[23];

    auto g_0_xxz_z_x = buffer_sfpp[24];

    auto g_0_xxz_z_y = buffer_sfpp[25];

    auto g_0_xxz_z_z = buffer_sfpp[26];

    auto g_0_xyy_x_x = buffer_sfpp[27];

    auto g_0_xyy_x_y = buffer_sfpp[28];

    auto g_0_xyy_x_z = buffer_sfpp[29];

    auto g_0_xyy_y_x = buffer_sfpp[30];

    auto g_0_xyy_y_y = buffer_sfpp[31];

    auto g_0_xyy_y_z = buffer_sfpp[32];

    auto g_0_xyy_z_x = buffer_sfpp[33];

    auto g_0_xyy_z_y = buffer_sfpp[34];

    auto g_0_xyy_z_z = buffer_sfpp[35];

    auto g_0_xyz_x_x = buffer_sfpp[36];

    auto g_0_xyz_x_y = buffer_sfpp[37];

    auto g_0_xyz_x_z = buffer_sfpp[38];

    auto g_0_xyz_y_x = buffer_sfpp[39];

    auto g_0_xyz_y_y = buffer_sfpp[40];

    auto g_0_xyz_y_z = buffer_sfpp[41];

    auto g_0_xyz_z_x = buffer_sfpp[42];

    auto g_0_xyz_z_y = buffer_sfpp[43];

    auto g_0_xyz_z_z = buffer_sfpp[44];

    auto g_0_xzz_x_x = buffer_sfpp[45];

    auto g_0_xzz_x_y = buffer_sfpp[46];

    auto g_0_xzz_x_z = buffer_sfpp[47];

    auto g_0_xzz_y_x = buffer_sfpp[48];

    auto g_0_xzz_y_y = buffer_sfpp[49];

    auto g_0_xzz_y_z = buffer_sfpp[50];

    auto g_0_xzz_z_x = buffer_sfpp[51];

    auto g_0_xzz_z_y = buffer_sfpp[52];

    auto g_0_xzz_z_z = buffer_sfpp[53];

    auto g_0_yyy_x_x = buffer_sfpp[54];

    auto g_0_yyy_x_y = buffer_sfpp[55];

    auto g_0_yyy_x_z = buffer_sfpp[56];

    auto g_0_yyy_y_x = buffer_sfpp[57];

    auto g_0_yyy_y_y = buffer_sfpp[58];

    auto g_0_yyy_y_z = buffer_sfpp[59];

    auto g_0_yyy_z_x = buffer_sfpp[60];

    auto g_0_yyy_z_y = buffer_sfpp[61];

    auto g_0_yyy_z_z = buffer_sfpp[62];

    auto g_0_yyz_x_x = buffer_sfpp[63];

    auto g_0_yyz_x_y = buffer_sfpp[64];

    auto g_0_yyz_x_z = buffer_sfpp[65];

    auto g_0_yyz_y_x = buffer_sfpp[66];

    auto g_0_yyz_y_y = buffer_sfpp[67];

    auto g_0_yyz_y_z = buffer_sfpp[68];

    auto g_0_yyz_z_x = buffer_sfpp[69];

    auto g_0_yyz_z_y = buffer_sfpp[70];

    auto g_0_yyz_z_z = buffer_sfpp[71];

    auto g_0_yzz_x_x = buffer_sfpp[72];

    auto g_0_yzz_x_y = buffer_sfpp[73];

    auto g_0_yzz_x_z = buffer_sfpp[74];

    auto g_0_yzz_y_x = buffer_sfpp[75];

    auto g_0_yzz_y_y = buffer_sfpp[76];

    auto g_0_yzz_y_z = buffer_sfpp[77];

    auto g_0_yzz_z_x = buffer_sfpp[78];

    auto g_0_yzz_z_y = buffer_sfpp[79];

    auto g_0_yzz_z_z = buffer_sfpp[80];

    auto g_0_zzz_x_x = buffer_sfpp[81];

    auto g_0_zzz_x_y = buffer_sfpp[82];

    auto g_0_zzz_x_z = buffer_sfpp[83];

    auto g_0_zzz_y_x = buffer_sfpp[84];

    auto g_0_zzz_y_y = buffer_sfpp[85];

    auto g_0_zzz_y_z = buffer_sfpp[86];

    auto g_0_zzz_z_x = buffer_sfpp[87];

    auto g_0_zzz_z_y = buffer_sfpp[88];

    auto g_0_zzz_z_z = buffer_sfpp[89];

    /// Set up components of auxilary buffer : buffer_dppp

    auto g_xx_x_x_x = buffer_dppp[0];

    auto g_xx_x_x_y = buffer_dppp[1];

    auto g_xx_x_x_z = buffer_dppp[2];

    auto g_xx_x_y_x = buffer_dppp[3];

    auto g_xx_x_y_y = buffer_dppp[4];

    auto g_xx_x_y_z = buffer_dppp[5];

    auto g_xx_x_z_x = buffer_dppp[6];

    auto g_xx_x_z_y = buffer_dppp[7];

    auto g_xx_x_z_z = buffer_dppp[8];

    auto g_xx_y_x_x = buffer_dppp[9];

    auto g_xx_y_x_y = buffer_dppp[10];

    auto g_xx_y_x_z = buffer_dppp[11];

    auto g_xx_y_y_x = buffer_dppp[12];

    auto g_xx_y_y_y = buffer_dppp[13];

    auto g_xx_y_y_z = buffer_dppp[14];

    auto g_xx_y_z_x = buffer_dppp[15];

    auto g_xx_y_z_y = buffer_dppp[16];

    auto g_xx_y_z_z = buffer_dppp[17];

    auto g_xx_z_x_x = buffer_dppp[18];

    auto g_xx_z_x_y = buffer_dppp[19];

    auto g_xx_z_x_z = buffer_dppp[20];

    auto g_xx_z_y_x = buffer_dppp[21];

    auto g_xx_z_y_y = buffer_dppp[22];

    auto g_xx_z_y_z = buffer_dppp[23];

    auto g_xx_z_z_x = buffer_dppp[24];

    auto g_xx_z_z_y = buffer_dppp[25];

    auto g_xx_z_z_z = buffer_dppp[26];

    auto g_xy_x_x_x = buffer_dppp[27];

    auto g_xy_x_x_y = buffer_dppp[28];

    auto g_xy_x_x_z = buffer_dppp[29];

    auto g_xy_x_y_x = buffer_dppp[30];

    auto g_xy_x_y_y = buffer_dppp[31];

    auto g_xy_x_y_z = buffer_dppp[32];

    auto g_xy_x_z_x = buffer_dppp[33];

    auto g_xy_x_z_y = buffer_dppp[34];

    auto g_xy_x_z_z = buffer_dppp[35];

    auto g_xy_y_x_x = buffer_dppp[36];

    auto g_xy_y_x_y = buffer_dppp[37];

    auto g_xy_y_x_z = buffer_dppp[38];

    auto g_xy_y_y_x = buffer_dppp[39];

    auto g_xy_y_y_y = buffer_dppp[40];

    auto g_xy_y_y_z = buffer_dppp[41];

    auto g_xy_y_z_x = buffer_dppp[42];

    auto g_xy_y_z_y = buffer_dppp[43];

    auto g_xy_y_z_z = buffer_dppp[44];

    auto g_xy_z_x_x = buffer_dppp[45];

    auto g_xy_z_x_y = buffer_dppp[46];

    auto g_xy_z_x_z = buffer_dppp[47];

    auto g_xy_z_y_x = buffer_dppp[48];

    auto g_xy_z_y_y = buffer_dppp[49];

    auto g_xy_z_y_z = buffer_dppp[50];

    auto g_xy_z_z_x = buffer_dppp[51];

    auto g_xy_z_z_y = buffer_dppp[52];

    auto g_xy_z_z_z = buffer_dppp[53];

    auto g_xz_x_x_x = buffer_dppp[54];

    auto g_xz_x_x_y = buffer_dppp[55];

    auto g_xz_x_x_z = buffer_dppp[56];

    auto g_xz_x_y_x = buffer_dppp[57];

    auto g_xz_x_y_y = buffer_dppp[58];

    auto g_xz_x_y_z = buffer_dppp[59];

    auto g_xz_x_z_x = buffer_dppp[60];

    auto g_xz_x_z_y = buffer_dppp[61];

    auto g_xz_x_z_z = buffer_dppp[62];

    auto g_xz_y_x_x = buffer_dppp[63];

    auto g_xz_y_x_y = buffer_dppp[64];

    auto g_xz_y_x_z = buffer_dppp[65];

    auto g_xz_y_y_x = buffer_dppp[66];

    auto g_xz_y_y_y = buffer_dppp[67];

    auto g_xz_y_y_z = buffer_dppp[68];

    auto g_xz_y_z_x = buffer_dppp[69];

    auto g_xz_y_z_y = buffer_dppp[70];

    auto g_xz_y_z_z = buffer_dppp[71];

    auto g_xz_z_x_x = buffer_dppp[72];

    auto g_xz_z_x_y = buffer_dppp[73];

    auto g_xz_z_x_z = buffer_dppp[74];

    auto g_xz_z_y_x = buffer_dppp[75];

    auto g_xz_z_y_y = buffer_dppp[76];

    auto g_xz_z_y_z = buffer_dppp[77];

    auto g_xz_z_z_x = buffer_dppp[78];

    auto g_xz_z_z_y = buffer_dppp[79];

    auto g_xz_z_z_z = buffer_dppp[80];

    auto g_yy_x_x_x = buffer_dppp[81];

    auto g_yy_x_x_y = buffer_dppp[82];

    auto g_yy_x_x_z = buffer_dppp[83];

    auto g_yy_x_y_x = buffer_dppp[84];

    auto g_yy_x_y_y = buffer_dppp[85];

    auto g_yy_x_y_z = buffer_dppp[86];

    auto g_yy_x_z_x = buffer_dppp[87];

    auto g_yy_x_z_y = buffer_dppp[88];

    auto g_yy_x_z_z = buffer_dppp[89];

    auto g_yy_y_x_x = buffer_dppp[90];

    auto g_yy_y_x_y = buffer_dppp[91];

    auto g_yy_y_x_z = buffer_dppp[92];

    auto g_yy_y_y_x = buffer_dppp[93];

    auto g_yy_y_y_y = buffer_dppp[94];

    auto g_yy_y_y_z = buffer_dppp[95];

    auto g_yy_y_z_x = buffer_dppp[96];

    auto g_yy_y_z_y = buffer_dppp[97];

    auto g_yy_y_z_z = buffer_dppp[98];

    auto g_yy_z_x_x = buffer_dppp[99];

    auto g_yy_z_x_y = buffer_dppp[100];

    auto g_yy_z_x_z = buffer_dppp[101];

    auto g_yy_z_y_x = buffer_dppp[102];

    auto g_yy_z_y_y = buffer_dppp[103];

    auto g_yy_z_y_z = buffer_dppp[104];

    auto g_yy_z_z_x = buffer_dppp[105];

    auto g_yy_z_z_y = buffer_dppp[106];

    auto g_yy_z_z_z = buffer_dppp[107];

    auto g_yz_x_x_x = buffer_dppp[108];

    auto g_yz_x_x_y = buffer_dppp[109];

    auto g_yz_x_x_z = buffer_dppp[110];

    auto g_yz_x_y_x = buffer_dppp[111];

    auto g_yz_x_y_y = buffer_dppp[112];

    auto g_yz_x_y_z = buffer_dppp[113];

    auto g_yz_x_z_x = buffer_dppp[114];

    auto g_yz_x_z_y = buffer_dppp[115];

    auto g_yz_x_z_z = buffer_dppp[116];

    auto g_yz_y_x_x = buffer_dppp[117];

    auto g_yz_y_x_y = buffer_dppp[118];

    auto g_yz_y_x_z = buffer_dppp[119];

    auto g_yz_y_y_x = buffer_dppp[120];

    auto g_yz_y_y_y = buffer_dppp[121];

    auto g_yz_y_y_z = buffer_dppp[122];

    auto g_yz_y_z_x = buffer_dppp[123];

    auto g_yz_y_z_y = buffer_dppp[124];

    auto g_yz_y_z_z = buffer_dppp[125];

    auto g_yz_z_x_x = buffer_dppp[126];

    auto g_yz_z_x_y = buffer_dppp[127];

    auto g_yz_z_x_z = buffer_dppp[128];

    auto g_yz_z_y_x = buffer_dppp[129];

    auto g_yz_z_y_y = buffer_dppp[130];

    auto g_yz_z_y_z = buffer_dppp[131];

    auto g_yz_z_z_x = buffer_dppp[132];

    auto g_yz_z_z_y = buffer_dppp[133];

    auto g_yz_z_z_z = buffer_dppp[134];

    auto g_zz_x_x_x = buffer_dppp[135];

    auto g_zz_x_x_y = buffer_dppp[136];

    auto g_zz_x_x_z = buffer_dppp[137];

    auto g_zz_x_y_x = buffer_dppp[138];

    auto g_zz_x_y_y = buffer_dppp[139];

    auto g_zz_x_y_z = buffer_dppp[140];

    auto g_zz_x_z_x = buffer_dppp[141];

    auto g_zz_x_z_y = buffer_dppp[142];

    auto g_zz_x_z_z = buffer_dppp[143];

    auto g_zz_y_x_x = buffer_dppp[144];

    auto g_zz_y_x_y = buffer_dppp[145];

    auto g_zz_y_x_z = buffer_dppp[146];

    auto g_zz_y_y_x = buffer_dppp[147];

    auto g_zz_y_y_y = buffer_dppp[148];

    auto g_zz_y_y_z = buffer_dppp[149];

    auto g_zz_y_z_x = buffer_dppp[150];

    auto g_zz_y_z_y = buffer_dppp[151];

    auto g_zz_y_z_z = buffer_dppp[152];

    auto g_zz_z_x_x = buffer_dppp[153];

    auto g_zz_z_x_y = buffer_dppp[154];

    auto g_zz_z_x_z = buffer_dppp[155];

    auto g_zz_z_y_x = buffer_dppp[156];

    auto g_zz_z_y_y = buffer_dppp[157];

    auto g_zz_z_y_z = buffer_dppp[158];

    auto g_zz_z_z_x = buffer_dppp[159];

    auto g_zz_z_z_y = buffer_dppp[160];

    auto g_zz_z_z_z = buffer_dppp[161];

    /// Set up components of auxilary buffer : buffer_dfpp

    auto g_xx_xxx_x_x = buffer_dfpp[0];

    auto g_xx_xxx_x_y = buffer_dfpp[1];

    auto g_xx_xxx_x_z = buffer_dfpp[2];

    auto g_xx_xxx_y_x = buffer_dfpp[3];

    auto g_xx_xxx_y_y = buffer_dfpp[4];

    auto g_xx_xxx_y_z = buffer_dfpp[5];

    auto g_xx_xxx_z_x = buffer_dfpp[6];

    auto g_xx_xxx_z_y = buffer_dfpp[7];

    auto g_xx_xxx_z_z = buffer_dfpp[8];

    auto g_xx_xxy_x_x = buffer_dfpp[9];

    auto g_xx_xxy_x_y = buffer_dfpp[10];

    auto g_xx_xxy_x_z = buffer_dfpp[11];

    auto g_xx_xxy_y_x = buffer_dfpp[12];

    auto g_xx_xxy_y_y = buffer_dfpp[13];

    auto g_xx_xxy_y_z = buffer_dfpp[14];

    auto g_xx_xxy_z_x = buffer_dfpp[15];

    auto g_xx_xxy_z_y = buffer_dfpp[16];

    auto g_xx_xxy_z_z = buffer_dfpp[17];

    auto g_xx_xxz_x_x = buffer_dfpp[18];

    auto g_xx_xxz_x_y = buffer_dfpp[19];

    auto g_xx_xxz_x_z = buffer_dfpp[20];

    auto g_xx_xxz_y_x = buffer_dfpp[21];

    auto g_xx_xxz_y_y = buffer_dfpp[22];

    auto g_xx_xxz_y_z = buffer_dfpp[23];

    auto g_xx_xxz_z_x = buffer_dfpp[24];

    auto g_xx_xxz_z_y = buffer_dfpp[25];

    auto g_xx_xxz_z_z = buffer_dfpp[26];

    auto g_xx_xyy_x_x = buffer_dfpp[27];

    auto g_xx_xyy_x_y = buffer_dfpp[28];

    auto g_xx_xyy_x_z = buffer_dfpp[29];

    auto g_xx_xyy_y_x = buffer_dfpp[30];

    auto g_xx_xyy_y_y = buffer_dfpp[31];

    auto g_xx_xyy_y_z = buffer_dfpp[32];

    auto g_xx_xyy_z_x = buffer_dfpp[33];

    auto g_xx_xyy_z_y = buffer_dfpp[34];

    auto g_xx_xyy_z_z = buffer_dfpp[35];

    auto g_xx_xyz_x_x = buffer_dfpp[36];

    auto g_xx_xyz_x_y = buffer_dfpp[37];

    auto g_xx_xyz_x_z = buffer_dfpp[38];

    auto g_xx_xyz_y_x = buffer_dfpp[39];

    auto g_xx_xyz_y_y = buffer_dfpp[40];

    auto g_xx_xyz_y_z = buffer_dfpp[41];

    auto g_xx_xyz_z_x = buffer_dfpp[42];

    auto g_xx_xyz_z_y = buffer_dfpp[43];

    auto g_xx_xyz_z_z = buffer_dfpp[44];

    auto g_xx_xzz_x_x = buffer_dfpp[45];

    auto g_xx_xzz_x_y = buffer_dfpp[46];

    auto g_xx_xzz_x_z = buffer_dfpp[47];

    auto g_xx_xzz_y_x = buffer_dfpp[48];

    auto g_xx_xzz_y_y = buffer_dfpp[49];

    auto g_xx_xzz_y_z = buffer_dfpp[50];

    auto g_xx_xzz_z_x = buffer_dfpp[51];

    auto g_xx_xzz_z_y = buffer_dfpp[52];

    auto g_xx_xzz_z_z = buffer_dfpp[53];

    auto g_xx_yyy_x_x = buffer_dfpp[54];

    auto g_xx_yyy_x_y = buffer_dfpp[55];

    auto g_xx_yyy_x_z = buffer_dfpp[56];

    auto g_xx_yyy_y_x = buffer_dfpp[57];

    auto g_xx_yyy_y_y = buffer_dfpp[58];

    auto g_xx_yyy_y_z = buffer_dfpp[59];

    auto g_xx_yyy_z_x = buffer_dfpp[60];

    auto g_xx_yyy_z_y = buffer_dfpp[61];

    auto g_xx_yyy_z_z = buffer_dfpp[62];

    auto g_xx_yyz_x_x = buffer_dfpp[63];

    auto g_xx_yyz_x_y = buffer_dfpp[64];

    auto g_xx_yyz_x_z = buffer_dfpp[65];

    auto g_xx_yyz_y_x = buffer_dfpp[66];

    auto g_xx_yyz_y_y = buffer_dfpp[67];

    auto g_xx_yyz_y_z = buffer_dfpp[68];

    auto g_xx_yyz_z_x = buffer_dfpp[69];

    auto g_xx_yyz_z_y = buffer_dfpp[70];

    auto g_xx_yyz_z_z = buffer_dfpp[71];

    auto g_xx_yzz_x_x = buffer_dfpp[72];

    auto g_xx_yzz_x_y = buffer_dfpp[73];

    auto g_xx_yzz_x_z = buffer_dfpp[74];

    auto g_xx_yzz_y_x = buffer_dfpp[75];

    auto g_xx_yzz_y_y = buffer_dfpp[76];

    auto g_xx_yzz_y_z = buffer_dfpp[77];

    auto g_xx_yzz_z_x = buffer_dfpp[78];

    auto g_xx_yzz_z_y = buffer_dfpp[79];

    auto g_xx_yzz_z_z = buffer_dfpp[80];

    auto g_xx_zzz_x_x = buffer_dfpp[81];

    auto g_xx_zzz_x_y = buffer_dfpp[82];

    auto g_xx_zzz_x_z = buffer_dfpp[83];

    auto g_xx_zzz_y_x = buffer_dfpp[84];

    auto g_xx_zzz_y_y = buffer_dfpp[85];

    auto g_xx_zzz_y_z = buffer_dfpp[86];

    auto g_xx_zzz_z_x = buffer_dfpp[87];

    auto g_xx_zzz_z_y = buffer_dfpp[88];

    auto g_xx_zzz_z_z = buffer_dfpp[89];

    auto g_xy_xxx_x_x = buffer_dfpp[90];

    auto g_xy_xxx_x_y = buffer_dfpp[91];

    auto g_xy_xxx_x_z = buffer_dfpp[92];

    auto g_xy_xxx_y_x = buffer_dfpp[93];

    auto g_xy_xxx_y_y = buffer_dfpp[94];

    auto g_xy_xxx_y_z = buffer_dfpp[95];

    auto g_xy_xxx_z_x = buffer_dfpp[96];

    auto g_xy_xxx_z_y = buffer_dfpp[97];

    auto g_xy_xxx_z_z = buffer_dfpp[98];

    auto g_xy_xxy_x_x = buffer_dfpp[99];

    auto g_xy_xxy_x_y = buffer_dfpp[100];

    auto g_xy_xxy_x_z = buffer_dfpp[101];

    auto g_xy_xxy_y_x = buffer_dfpp[102];

    auto g_xy_xxy_y_y = buffer_dfpp[103];

    auto g_xy_xxy_y_z = buffer_dfpp[104];

    auto g_xy_xxy_z_x = buffer_dfpp[105];

    auto g_xy_xxy_z_y = buffer_dfpp[106];

    auto g_xy_xxy_z_z = buffer_dfpp[107];

    auto g_xy_xxz_x_x = buffer_dfpp[108];

    auto g_xy_xxz_x_y = buffer_dfpp[109];

    auto g_xy_xxz_x_z = buffer_dfpp[110];

    auto g_xy_xxz_y_x = buffer_dfpp[111];

    auto g_xy_xxz_y_y = buffer_dfpp[112];

    auto g_xy_xxz_y_z = buffer_dfpp[113];

    auto g_xy_xxz_z_x = buffer_dfpp[114];

    auto g_xy_xxz_z_y = buffer_dfpp[115];

    auto g_xy_xxz_z_z = buffer_dfpp[116];

    auto g_xy_xyy_x_x = buffer_dfpp[117];

    auto g_xy_xyy_x_y = buffer_dfpp[118];

    auto g_xy_xyy_x_z = buffer_dfpp[119];

    auto g_xy_xyy_y_x = buffer_dfpp[120];

    auto g_xy_xyy_y_y = buffer_dfpp[121];

    auto g_xy_xyy_y_z = buffer_dfpp[122];

    auto g_xy_xyy_z_x = buffer_dfpp[123];

    auto g_xy_xyy_z_y = buffer_dfpp[124];

    auto g_xy_xyy_z_z = buffer_dfpp[125];

    auto g_xy_xyz_x_x = buffer_dfpp[126];

    auto g_xy_xyz_x_y = buffer_dfpp[127];

    auto g_xy_xyz_x_z = buffer_dfpp[128];

    auto g_xy_xyz_y_x = buffer_dfpp[129];

    auto g_xy_xyz_y_y = buffer_dfpp[130];

    auto g_xy_xyz_y_z = buffer_dfpp[131];

    auto g_xy_xyz_z_x = buffer_dfpp[132];

    auto g_xy_xyz_z_y = buffer_dfpp[133];

    auto g_xy_xyz_z_z = buffer_dfpp[134];

    auto g_xy_xzz_x_x = buffer_dfpp[135];

    auto g_xy_xzz_x_y = buffer_dfpp[136];

    auto g_xy_xzz_x_z = buffer_dfpp[137];

    auto g_xy_xzz_y_x = buffer_dfpp[138];

    auto g_xy_xzz_y_y = buffer_dfpp[139];

    auto g_xy_xzz_y_z = buffer_dfpp[140];

    auto g_xy_xzz_z_x = buffer_dfpp[141];

    auto g_xy_xzz_z_y = buffer_dfpp[142];

    auto g_xy_xzz_z_z = buffer_dfpp[143];

    auto g_xy_yyy_x_x = buffer_dfpp[144];

    auto g_xy_yyy_x_y = buffer_dfpp[145];

    auto g_xy_yyy_x_z = buffer_dfpp[146];

    auto g_xy_yyy_y_x = buffer_dfpp[147];

    auto g_xy_yyy_y_y = buffer_dfpp[148];

    auto g_xy_yyy_y_z = buffer_dfpp[149];

    auto g_xy_yyy_z_x = buffer_dfpp[150];

    auto g_xy_yyy_z_y = buffer_dfpp[151];

    auto g_xy_yyy_z_z = buffer_dfpp[152];

    auto g_xy_yyz_x_x = buffer_dfpp[153];

    auto g_xy_yyz_x_y = buffer_dfpp[154];

    auto g_xy_yyz_x_z = buffer_dfpp[155];

    auto g_xy_yyz_y_x = buffer_dfpp[156];

    auto g_xy_yyz_y_y = buffer_dfpp[157];

    auto g_xy_yyz_y_z = buffer_dfpp[158];

    auto g_xy_yyz_z_x = buffer_dfpp[159];

    auto g_xy_yyz_z_y = buffer_dfpp[160];

    auto g_xy_yyz_z_z = buffer_dfpp[161];

    auto g_xy_yzz_x_x = buffer_dfpp[162];

    auto g_xy_yzz_x_y = buffer_dfpp[163];

    auto g_xy_yzz_x_z = buffer_dfpp[164];

    auto g_xy_yzz_y_x = buffer_dfpp[165];

    auto g_xy_yzz_y_y = buffer_dfpp[166];

    auto g_xy_yzz_y_z = buffer_dfpp[167];

    auto g_xy_yzz_z_x = buffer_dfpp[168];

    auto g_xy_yzz_z_y = buffer_dfpp[169];

    auto g_xy_yzz_z_z = buffer_dfpp[170];

    auto g_xy_zzz_x_x = buffer_dfpp[171];

    auto g_xy_zzz_x_y = buffer_dfpp[172];

    auto g_xy_zzz_x_z = buffer_dfpp[173];

    auto g_xy_zzz_y_x = buffer_dfpp[174];

    auto g_xy_zzz_y_y = buffer_dfpp[175];

    auto g_xy_zzz_y_z = buffer_dfpp[176];

    auto g_xy_zzz_z_x = buffer_dfpp[177];

    auto g_xy_zzz_z_y = buffer_dfpp[178];

    auto g_xy_zzz_z_z = buffer_dfpp[179];

    auto g_xz_xxx_x_x = buffer_dfpp[180];

    auto g_xz_xxx_x_y = buffer_dfpp[181];

    auto g_xz_xxx_x_z = buffer_dfpp[182];

    auto g_xz_xxx_y_x = buffer_dfpp[183];

    auto g_xz_xxx_y_y = buffer_dfpp[184];

    auto g_xz_xxx_y_z = buffer_dfpp[185];

    auto g_xz_xxx_z_x = buffer_dfpp[186];

    auto g_xz_xxx_z_y = buffer_dfpp[187];

    auto g_xz_xxx_z_z = buffer_dfpp[188];

    auto g_xz_xxy_x_x = buffer_dfpp[189];

    auto g_xz_xxy_x_y = buffer_dfpp[190];

    auto g_xz_xxy_x_z = buffer_dfpp[191];

    auto g_xz_xxy_y_x = buffer_dfpp[192];

    auto g_xz_xxy_y_y = buffer_dfpp[193];

    auto g_xz_xxy_y_z = buffer_dfpp[194];

    auto g_xz_xxy_z_x = buffer_dfpp[195];

    auto g_xz_xxy_z_y = buffer_dfpp[196];

    auto g_xz_xxy_z_z = buffer_dfpp[197];

    auto g_xz_xxz_x_x = buffer_dfpp[198];

    auto g_xz_xxz_x_y = buffer_dfpp[199];

    auto g_xz_xxz_x_z = buffer_dfpp[200];

    auto g_xz_xxz_y_x = buffer_dfpp[201];

    auto g_xz_xxz_y_y = buffer_dfpp[202];

    auto g_xz_xxz_y_z = buffer_dfpp[203];

    auto g_xz_xxz_z_x = buffer_dfpp[204];

    auto g_xz_xxz_z_y = buffer_dfpp[205];

    auto g_xz_xxz_z_z = buffer_dfpp[206];

    auto g_xz_xyy_x_x = buffer_dfpp[207];

    auto g_xz_xyy_x_y = buffer_dfpp[208];

    auto g_xz_xyy_x_z = buffer_dfpp[209];

    auto g_xz_xyy_y_x = buffer_dfpp[210];

    auto g_xz_xyy_y_y = buffer_dfpp[211];

    auto g_xz_xyy_y_z = buffer_dfpp[212];

    auto g_xz_xyy_z_x = buffer_dfpp[213];

    auto g_xz_xyy_z_y = buffer_dfpp[214];

    auto g_xz_xyy_z_z = buffer_dfpp[215];

    auto g_xz_xyz_x_x = buffer_dfpp[216];

    auto g_xz_xyz_x_y = buffer_dfpp[217];

    auto g_xz_xyz_x_z = buffer_dfpp[218];

    auto g_xz_xyz_y_x = buffer_dfpp[219];

    auto g_xz_xyz_y_y = buffer_dfpp[220];

    auto g_xz_xyz_y_z = buffer_dfpp[221];

    auto g_xz_xyz_z_x = buffer_dfpp[222];

    auto g_xz_xyz_z_y = buffer_dfpp[223];

    auto g_xz_xyz_z_z = buffer_dfpp[224];

    auto g_xz_xzz_x_x = buffer_dfpp[225];

    auto g_xz_xzz_x_y = buffer_dfpp[226];

    auto g_xz_xzz_x_z = buffer_dfpp[227];

    auto g_xz_xzz_y_x = buffer_dfpp[228];

    auto g_xz_xzz_y_y = buffer_dfpp[229];

    auto g_xz_xzz_y_z = buffer_dfpp[230];

    auto g_xz_xzz_z_x = buffer_dfpp[231];

    auto g_xz_xzz_z_y = buffer_dfpp[232];

    auto g_xz_xzz_z_z = buffer_dfpp[233];

    auto g_xz_yyy_x_x = buffer_dfpp[234];

    auto g_xz_yyy_x_y = buffer_dfpp[235];

    auto g_xz_yyy_x_z = buffer_dfpp[236];

    auto g_xz_yyy_y_x = buffer_dfpp[237];

    auto g_xz_yyy_y_y = buffer_dfpp[238];

    auto g_xz_yyy_y_z = buffer_dfpp[239];

    auto g_xz_yyy_z_x = buffer_dfpp[240];

    auto g_xz_yyy_z_y = buffer_dfpp[241];

    auto g_xz_yyy_z_z = buffer_dfpp[242];

    auto g_xz_yyz_x_x = buffer_dfpp[243];

    auto g_xz_yyz_x_y = buffer_dfpp[244];

    auto g_xz_yyz_x_z = buffer_dfpp[245];

    auto g_xz_yyz_y_x = buffer_dfpp[246];

    auto g_xz_yyz_y_y = buffer_dfpp[247];

    auto g_xz_yyz_y_z = buffer_dfpp[248];

    auto g_xz_yyz_z_x = buffer_dfpp[249];

    auto g_xz_yyz_z_y = buffer_dfpp[250];

    auto g_xz_yyz_z_z = buffer_dfpp[251];

    auto g_xz_yzz_x_x = buffer_dfpp[252];

    auto g_xz_yzz_x_y = buffer_dfpp[253];

    auto g_xz_yzz_x_z = buffer_dfpp[254];

    auto g_xz_yzz_y_x = buffer_dfpp[255];

    auto g_xz_yzz_y_y = buffer_dfpp[256];

    auto g_xz_yzz_y_z = buffer_dfpp[257];

    auto g_xz_yzz_z_x = buffer_dfpp[258];

    auto g_xz_yzz_z_y = buffer_dfpp[259];

    auto g_xz_yzz_z_z = buffer_dfpp[260];

    auto g_xz_zzz_x_x = buffer_dfpp[261];

    auto g_xz_zzz_x_y = buffer_dfpp[262];

    auto g_xz_zzz_x_z = buffer_dfpp[263];

    auto g_xz_zzz_y_x = buffer_dfpp[264];

    auto g_xz_zzz_y_y = buffer_dfpp[265];

    auto g_xz_zzz_y_z = buffer_dfpp[266];

    auto g_xz_zzz_z_x = buffer_dfpp[267];

    auto g_xz_zzz_z_y = buffer_dfpp[268];

    auto g_xz_zzz_z_z = buffer_dfpp[269];

    auto g_yy_xxx_x_x = buffer_dfpp[270];

    auto g_yy_xxx_x_y = buffer_dfpp[271];

    auto g_yy_xxx_x_z = buffer_dfpp[272];

    auto g_yy_xxx_y_x = buffer_dfpp[273];

    auto g_yy_xxx_y_y = buffer_dfpp[274];

    auto g_yy_xxx_y_z = buffer_dfpp[275];

    auto g_yy_xxx_z_x = buffer_dfpp[276];

    auto g_yy_xxx_z_y = buffer_dfpp[277];

    auto g_yy_xxx_z_z = buffer_dfpp[278];

    auto g_yy_xxy_x_x = buffer_dfpp[279];

    auto g_yy_xxy_x_y = buffer_dfpp[280];

    auto g_yy_xxy_x_z = buffer_dfpp[281];

    auto g_yy_xxy_y_x = buffer_dfpp[282];

    auto g_yy_xxy_y_y = buffer_dfpp[283];

    auto g_yy_xxy_y_z = buffer_dfpp[284];

    auto g_yy_xxy_z_x = buffer_dfpp[285];

    auto g_yy_xxy_z_y = buffer_dfpp[286];

    auto g_yy_xxy_z_z = buffer_dfpp[287];

    auto g_yy_xxz_x_x = buffer_dfpp[288];

    auto g_yy_xxz_x_y = buffer_dfpp[289];

    auto g_yy_xxz_x_z = buffer_dfpp[290];

    auto g_yy_xxz_y_x = buffer_dfpp[291];

    auto g_yy_xxz_y_y = buffer_dfpp[292];

    auto g_yy_xxz_y_z = buffer_dfpp[293];

    auto g_yy_xxz_z_x = buffer_dfpp[294];

    auto g_yy_xxz_z_y = buffer_dfpp[295];

    auto g_yy_xxz_z_z = buffer_dfpp[296];

    auto g_yy_xyy_x_x = buffer_dfpp[297];

    auto g_yy_xyy_x_y = buffer_dfpp[298];

    auto g_yy_xyy_x_z = buffer_dfpp[299];

    auto g_yy_xyy_y_x = buffer_dfpp[300];

    auto g_yy_xyy_y_y = buffer_dfpp[301];

    auto g_yy_xyy_y_z = buffer_dfpp[302];

    auto g_yy_xyy_z_x = buffer_dfpp[303];

    auto g_yy_xyy_z_y = buffer_dfpp[304];

    auto g_yy_xyy_z_z = buffer_dfpp[305];

    auto g_yy_xyz_x_x = buffer_dfpp[306];

    auto g_yy_xyz_x_y = buffer_dfpp[307];

    auto g_yy_xyz_x_z = buffer_dfpp[308];

    auto g_yy_xyz_y_x = buffer_dfpp[309];

    auto g_yy_xyz_y_y = buffer_dfpp[310];

    auto g_yy_xyz_y_z = buffer_dfpp[311];

    auto g_yy_xyz_z_x = buffer_dfpp[312];

    auto g_yy_xyz_z_y = buffer_dfpp[313];

    auto g_yy_xyz_z_z = buffer_dfpp[314];

    auto g_yy_xzz_x_x = buffer_dfpp[315];

    auto g_yy_xzz_x_y = buffer_dfpp[316];

    auto g_yy_xzz_x_z = buffer_dfpp[317];

    auto g_yy_xzz_y_x = buffer_dfpp[318];

    auto g_yy_xzz_y_y = buffer_dfpp[319];

    auto g_yy_xzz_y_z = buffer_dfpp[320];

    auto g_yy_xzz_z_x = buffer_dfpp[321];

    auto g_yy_xzz_z_y = buffer_dfpp[322];

    auto g_yy_xzz_z_z = buffer_dfpp[323];

    auto g_yy_yyy_x_x = buffer_dfpp[324];

    auto g_yy_yyy_x_y = buffer_dfpp[325];

    auto g_yy_yyy_x_z = buffer_dfpp[326];

    auto g_yy_yyy_y_x = buffer_dfpp[327];

    auto g_yy_yyy_y_y = buffer_dfpp[328];

    auto g_yy_yyy_y_z = buffer_dfpp[329];

    auto g_yy_yyy_z_x = buffer_dfpp[330];

    auto g_yy_yyy_z_y = buffer_dfpp[331];

    auto g_yy_yyy_z_z = buffer_dfpp[332];

    auto g_yy_yyz_x_x = buffer_dfpp[333];

    auto g_yy_yyz_x_y = buffer_dfpp[334];

    auto g_yy_yyz_x_z = buffer_dfpp[335];

    auto g_yy_yyz_y_x = buffer_dfpp[336];

    auto g_yy_yyz_y_y = buffer_dfpp[337];

    auto g_yy_yyz_y_z = buffer_dfpp[338];

    auto g_yy_yyz_z_x = buffer_dfpp[339];

    auto g_yy_yyz_z_y = buffer_dfpp[340];

    auto g_yy_yyz_z_z = buffer_dfpp[341];

    auto g_yy_yzz_x_x = buffer_dfpp[342];

    auto g_yy_yzz_x_y = buffer_dfpp[343];

    auto g_yy_yzz_x_z = buffer_dfpp[344];

    auto g_yy_yzz_y_x = buffer_dfpp[345];

    auto g_yy_yzz_y_y = buffer_dfpp[346];

    auto g_yy_yzz_y_z = buffer_dfpp[347];

    auto g_yy_yzz_z_x = buffer_dfpp[348];

    auto g_yy_yzz_z_y = buffer_dfpp[349];

    auto g_yy_yzz_z_z = buffer_dfpp[350];

    auto g_yy_zzz_x_x = buffer_dfpp[351];

    auto g_yy_zzz_x_y = buffer_dfpp[352];

    auto g_yy_zzz_x_z = buffer_dfpp[353];

    auto g_yy_zzz_y_x = buffer_dfpp[354];

    auto g_yy_zzz_y_y = buffer_dfpp[355];

    auto g_yy_zzz_y_z = buffer_dfpp[356];

    auto g_yy_zzz_z_x = buffer_dfpp[357];

    auto g_yy_zzz_z_y = buffer_dfpp[358];

    auto g_yy_zzz_z_z = buffer_dfpp[359];

    auto g_yz_xxx_x_x = buffer_dfpp[360];

    auto g_yz_xxx_x_y = buffer_dfpp[361];

    auto g_yz_xxx_x_z = buffer_dfpp[362];

    auto g_yz_xxx_y_x = buffer_dfpp[363];

    auto g_yz_xxx_y_y = buffer_dfpp[364];

    auto g_yz_xxx_y_z = buffer_dfpp[365];

    auto g_yz_xxx_z_x = buffer_dfpp[366];

    auto g_yz_xxx_z_y = buffer_dfpp[367];

    auto g_yz_xxx_z_z = buffer_dfpp[368];

    auto g_yz_xxy_x_x = buffer_dfpp[369];

    auto g_yz_xxy_x_y = buffer_dfpp[370];

    auto g_yz_xxy_x_z = buffer_dfpp[371];

    auto g_yz_xxy_y_x = buffer_dfpp[372];

    auto g_yz_xxy_y_y = buffer_dfpp[373];

    auto g_yz_xxy_y_z = buffer_dfpp[374];

    auto g_yz_xxy_z_x = buffer_dfpp[375];

    auto g_yz_xxy_z_y = buffer_dfpp[376];

    auto g_yz_xxy_z_z = buffer_dfpp[377];

    auto g_yz_xxz_x_x = buffer_dfpp[378];

    auto g_yz_xxz_x_y = buffer_dfpp[379];

    auto g_yz_xxz_x_z = buffer_dfpp[380];

    auto g_yz_xxz_y_x = buffer_dfpp[381];

    auto g_yz_xxz_y_y = buffer_dfpp[382];

    auto g_yz_xxz_y_z = buffer_dfpp[383];

    auto g_yz_xxz_z_x = buffer_dfpp[384];

    auto g_yz_xxz_z_y = buffer_dfpp[385];

    auto g_yz_xxz_z_z = buffer_dfpp[386];

    auto g_yz_xyy_x_x = buffer_dfpp[387];

    auto g_yz_xyy_x_y = buffer_dfpp[388];

    auto g_yz_xyy_x_z = buffer_dfpp[389];

    auto g_yz_xyy_y_x = buffer_dfpp[390];

    auto g_yz_xyy_y_y = buffer_dfpp[391];

    auto g_yz_xyy_y_z = buffer_dfpp[392];

    auto g_yz_xyy_z_x = buffer_dfpp[393];

    auto g_yz_xyy_z_y = buffer_dfpp[394];

    auto g_yz_xyy_z_z = buffer_dfpp[395];

    auto g_yz_xyz_x_x = buffer_dfpp[396];

    auto g_yz_xyz_x_y = buffer_dfpp[397];

    auto g_yz_xyz_x_z = buffer_dfpp[398];

    auto g_yz_xyz_y_x = buffer_dfpp[399];

    auto g_yz_xyz_y_y = buffer_dfpp[400];

    auto g_yz_xyz_y_z = buffer_dfpp[401];

    auto g_yz_xyz_z_x = buffer_dfpp[402];

    auto g_yz_xyz_z_y = buffer_dfpp[403];

    auto g_yz_xyz_z_z = buffer_dfpp[404];

    auto g_yz_xzz_x_x = buffer_dfpp[405];

    auto g_yz_xzz_x_y = buffer_dfpp[406];

    auto g_yz_xzz_x_z = buffer_dfpp[407];

    auto g_yz_xzz_y_x = buffer_dfpp[408];

    auto g_yz_xzz_y_y = buffer_dfpp[409];

    auto g_yz_xzz_y_z = buffer_dfpp[410];

    auto g_yz_xzz_z_x = buffer_dfpp[411];

    auto g_yz_xzz_z_y = buffer_dfpp[412];

    auto g_yz_xzz_z_z = buffer_dfpp[413];

    auto g_yz_yyy_x_x = buffer_dfpp[414];

    auto g_yz_yyy_x_y = buffer_dfpp[415];

    auto g_yz_yyy_x_z = buffer_dfpp[416];

    auto g_yz_yyy_y_x = buffer_dfpp[417];

    auto g_yz_yyy_y_y = buffer_dfpp[418];

    auto g_yz_yyy_y_z = buffer_dfpp[419];

    auto g_yz_yyy_z_x = buffer_dfpp[420];

    auto g_yz_yyy_z_y = buffer_dfpp[421];

    auto g_yz_yyy_z_z = buffer_dfpp[422];

    auto g_yz_yyz_x_x = buffer_dfpp[423];

    auto g_yz_yyz_x_y = buffer_dfpp[424];

    auto g_yz_yyz_x_z = buffer_dfpp[425];

    auto g_yz_yyz_y_x = buffer_dfpp[426];

    auto g_yz_yyz_y_y = buffer_dfpp[427];

    auto g_yz_yyz_y_z = buffer_dfpp[428];

    auto g_yz_yyz_z_x = buffer_dfpp[429];

    auto g_yz_yyz_z_y = buffer_dfpp[430];

    auto g_yz_yyz_z_z = buffer_dfpp[431];

    auto g_yz_yzz_x_x = buffer_dfpp[432];

    auto g_yz_yzz_x_y = buffer_dfpp[433];

    auto g_yz_yzz_x_z = buffer_dfpp[434];

    auto g_yz_yzz_y_x = buffer_dfpp[435];

    auto g_yz_yzz_y_y = buffer_dfpp[436];

    auto g_yz_yzz_y_z = buffer_dfpp[437];

    auto g_yz_yzz_z_x = buffer_dfpp[438];

    auto g_yz_yzz_z_y = buffer_dfpp[439];

    auto g_yz_yzz_z_z = buffer_dfpp[440];

    auto g_yz_zzz_x_x = buffer_dfpp[441];

    auto g_yz_zzz_x_y = buffer_dfpp[442];

    auto g_yz_zzz_x_z = buffer_dfpp[443];

    auto g_yz_zzz_y_x = buffer_dfpp[444];

    auto g_yz_zzz_y_y = buffer_dfpp[445];

    auto g_yz_zzz_y_z = buffer_dfpp[446];

    auto g_yz_zzz_z_x = buffer_dfpp[447];

    auto g_yz_zzz_z_y = buffer_dfpp[448];

    auto g_yz_zzz_z_z = buffer_dfpp[449];

    auto g_zz_xxx_x_x = buffer_dfpp[450];

    auto g_zz_xxx_x_y = buffer_dfpp[451];

    auto g_zz_xxx_x_z = buffer_dfpp[452];

    auto g_zz_xxx_y_x = buffer_dfpp[453];

    auto g_zz_xxx_y_y = buffer_dfpp[454];

    auto g_zz_xxx_y_z = buffer_dfpp[455];

    auto g_zz_xxx_z_x = buffer_dfpp[456];

    auto g_zz_xxx_z_y = buffer_dfpp[457];

    auto g_zz_xxx_z_z = buffer_dfpp[458];

    auto g_zz_xxy_x_x = buffer_dfpp[459];

    auto g_zz_xxy_x_y = buffer_dfpp[460];

    auto g_zz_xxy_x_z = buffer_dfpp[461];

    auto g_zz_xxy_y_x = buffer_dfpp[462];

    auto g_zz_xxy_y_y = buffer_dfpp[463];

    auto g_zz_xxy_y_z = buffer_dfpp[464];

    auto g_zz_xxy_z_x = buffer_dfpp[465];

    auto g_zz_xxy_z_y = buffer_dfpp[466];

    auto g_zz_xxy_z_z = buffer_dfpp[467];

    auto g_zz_xxz_x_x = buffer_dfpp[468];

    auto g_zz_xxz_x_y = buffer_dfpp[469];

    auto g_zz_xxz_x_z = buffer_dfpp[470];

    auto g_zz_xxz_y_x = buffer_dfpp[471];

    auto g_zz_xxz_y_y = buffer_dfpp[472];

    auto g_zz_xxz_y_z = buffer_dfpp[473];

    auto g_zz_xxz_z_x = buffer_dfpp[474];

    auto g_zz_xxz_z_y = buffer_dfpp[475];

    auto g_zz_xxz_z_z = buffer_dfpp[476];

    auto g_zz_xyy_x_x = buffer_dfpp[477];

    auto g_zz_xyy_x_y = buffer_dfpp[478];

    auto g_zz_xyy_x_z = buffer_dfpp[479];

    auto g_zz_xyy_y_x = buffer_dfpp[480];

    auto g_zz_xyy_y_y = buffer_dfpp[481];

    auto g_zz_xyy_y_z = buffer_dfpp[482];

    auto g_zz_xyy_z_x = buffer_dfpp[483];

    auto g_zz_xyy_z_y = buffer_dfpp[484];

    auto g_zz_xyy_z_z = buffer_dfpp[485];

    auto g_zz_xyz_x_x = buffer_dfpp[486];

    auto g_zz_xyz_x_y = buffer_dfpp[487];

    auto g_zz_xyz_x_z = buffer_dfpp[488];

    auto g_zz_xyz_y_x = buffer_dfpp[489];

    auto g_zz_xyz_y_y = buffer_dfpp[490];

    auto g_zz_xyz_y_z = buffer_dfpp[491];

    auto g_zz_xyz_z_x = buffer_dfpp[492];

    auto g_zz_xyz_z_y = buffer_dfpp[493];

    auto g_zz_xyz_z_z = buffer_dfpp[494];

    auto g_zz_xzz_x_x = buffer_dfpp[495];

    auto g_zz_xzz_x_y = buffer_dfpp[496];

    auto g_zz_xzz_x_z = buffer_dfpp[497];

    auto g_zz_xzz_y_x = buffer_dfpp[498];

    auto g_zz_xzz_y_y = buffer_dfpp[499];

    auto g_zz_xzz_y_z = buffer_dfpp[500];

    auto g_zz_xzz_z_x = buffer_dfpp[501];

    auto g_zz_xzz_z_y = buffer_dfpp[502];

    auto g_zz_xzz_z_z = buffer_dfpp[503];

    auto g_zz_yyy_x_x = buffer_dfpp[504];

    auto g_zz_yyy_x_y = buffer_dfpp[505];

    auto g_zz_yyy_x_z = buffer_dfpp[506];

    auto g_zz_yyy_y_x = buffer_dfpp[507];

    auto g_zz_yyy_y_y = buffer_dfpp[508];

    auto g_zz_yyy_y_z = buffer_dfpp[509];

    auto g_zz_yyy_z_x = buffer_dfpp[510];

    auto g_zz_yyy_z_y = buffer_dfpp[511];

    auto g_zz_yyy_z_z = buffer_dfpp[512];

    auto g_zz_yyz_x_x = buffer_dfpp[513];

    auto g_zz_yyz_x_y = buffer_dfpp[514];

    auto g_zz_yyz_x_z = buffer_dfpp[515];

    auto g_zz_yyz_y_x = buffer_dfpp[516];

    auto g_zz_yyz_y_y = buffer_dfpp[517];

    auto g_zz_yyz_y_z = buffer_dfpp[518];

    auto g_zz_yyz_z_x = buffer_dfpp[519];

    auto g_zz_yyz_z_y = buffer_dfpp[520];

    auto g_zz_yyz_z_z = buffer_dfpp[521];

    auto g_zz_yzz_x_x = buffer_dfpp[522];

    auto g_zz_yzz_x_y = buffer_dfpp[523];

    auto g_zz_yzz_x_z = buffer_dfpp[524];

    auto g_zz_yzz_y_x = buffer_dfpp[525];

    auto g_zz_yzz_y_y = buffer_dfpp[526];

    auto g_zz_yzz_y_z = buffer_dfpp[527];

    auto g_zz_yzz_z_x = buffer_dfpp[528];

    auto g_zz_yzz_z_y = buffer_dfpp[529];

    auto g_zz_yzz_z_z = buffer_dfpp[530];

    auto g_zz_zzz_x_x = buffer_dfpp[531];

    auto g_zz_zzz_x_y = buffer_dfpp[532];

    auto g_zz_zzz_x_z = buffer_dfpp[533];

    auto g_zz_zzz_y_x = buffer_dfpp[534];

    auto g_zz_zzz_y_y = buffer_dfpp[535];

    auto g_zz_zzz_y_z = buffer_dfpp[536];

    auto g_zz_zzz_z_x = buffer_dfpp[537];

    auto g_zz_zzz_z_y = buffer_dfpp[538];

    auto g_zz_zzz_z_z = buffer_dfpp[539];

    /// Set up components of integrals buffer : buffer_1100_pdpp

    auto g_x_x_0_0_x_xx_x_x = buffer_1100_pdpp[0];

    auto g_x_x_0_0_x_xx_x_y = buffer_1100_pdpp[1];

    auto g_x_x_0_0_x_xx_x_z = buffer_1100_pdpp[2];

    auto g_x_x_0_0_x_xx_y_x = buffer_1100_pdpp[3];

    auto g_x_x_0_0_x_xx_y_y = buffer_1100_pdpp[4];

    auto g_x_x_0_0_x_xx_y_z = buffer_1100_pdpp[5];

    auto g_x_x_0_0_x_xx_z_x = buffer_1100_pdpp[6];

    auto g_x_x_0_0_x_xx_z_y = buffer_1100_pdpp[7];

    auto g_x_x_0_0_x_xx_z_z = buffer_1100_pdpp[8];

    auto g_x_x_0_0_x_xy_x_x = buffer_1100_pdpp[9];

    auto g_x_x_0_0_x_xy_x_y = buffer_1100_pdpp[10];

    auto g_x_x_0_0_x_xy_x_z = buffer_1100_pdpp[11];

    auto g_x_x_0_0_x_xy_y_x = buffer_1100_pdpp[12];

    auto g_x_x_0_0_x_xy_y_y = buffer_1100_pdpp[13];

    auto g_x_x_0_0_x_xy_y_z = buffer_1100_pdpp[14];

    auto g_x_x_0_0_x_xy_z_x = buffer_1100_pdpp[15];

    auto g_x_x_0_0_x_xy_z_y = buffer_1100_pdpp[16];

    auto g_x_x_0_0_x_xy_z_z = buffer_1100_pdpp[17];

    auto g_x_x_0_0_x_xz_x_x = buffer_1100_pdpp[18];

    auto g_x_x_0_0_x_xz_x_y = buffer_1100_pdpp[19];

    auto g_x_x_0_0_x_xz_x_z = buffer_1100_pdpp[20];

    auto g_x_x_0_0_x_xz_y_x = buffer_1100_pdpp[21];

    auto g_x_x_0_0_x_xz_y_y = buffer_1100_pdpp[22];

    auto g_x_x_0_0_x_xz_y_z = buffer_1100_pdpp[23];

    auto g_x_x_0_0_x_xz_z_x = buffer_1100_pdpp[24];

    auto g_x_x_0_0_x_xz_z_y = buffer_1100_pdpp[25];

    auto g_x_x_0_0_x_xz_z_z = buffer_1100_pdpp[26];

    auto g_x_x_0_0_x_yy_x_x = buffer_1100_pdpp[27];

    auto g_x_x_0_0_x_yy_x_y = buffer_1100_pdpp[28];

    auto g_x_x_0_0_x_yy_x_z = buffer_1100_pdpp[29];

    auto g_x_x_0_0_x_yy_y_x = buffer_1100_pdpp[30];

    auto g_x_x_0_0_x_yy_y_y = buffer_1100_pdpp[31];

    auto g_x_x_0_0_x_yy_y_z = buffer_1100_pdpp[32];

    auto g_x_x_0_0_x_yy_z_x = buffer_1100_pdpp[33];

    auto g_x_x_0_0_x_yy_z_y = buffer_1100_pdpp[34];

    auto g_x_x_0_0_x_yy_z_z = buffer_1100_pdpp[35];

    auto g_x_x_0_0_x_yz_x_x = buffer_1100_pdpp[36];

    auto g_x_x_0_0_x_yz_x_y = buffer_1100_pdpp[37];

    auto g_x_x_0_0_x_yz_x_z = buffer_1100_pdpp[38];

    auto g_x_x_0_0_x_yz_y_x = buffer_1100_pdpp[39];

    auto g_x_x_0_0_x_yz_y_y = buffer_1100_pdpp[40];

    auto g_x_x_0_0_x_yz_y_z = buffer_1100_pdpp[41];

    auto g_x_x_0_0_x_yz_z_x = buffer_1100_pdpp[42];

    auto g_x_x_0_0_x_yz_z_y = buffer_1100_pdpp[43];

    auto g_x_x_0_0_x_yz_z_z = buffer_1100_pdpp[44];

    auto g_x_x_0_0_x_zz_x_x = buffer_1100_pdpp[45];

    auto g_x_x_0_0_x_zz_x_y = buffer_1100_pdpp[46];

    auto g_x_x_0_0_x_zz_x_z = buffer_1100_pdpp[47];

    auto g_x_x_0_0_x_zz_y_x = buffer_1100_pdpp[48];

    auto g_x_x_0_0_x_zz_y_y = buffer_1100_pdpp[49];

    auto g_x_x_0_0_x_zz_y_z = buffer_1100_pdpp[50];

    auto g_x_x_0_0_x_zz_z_x = buffer_1100_pdpp[51];

    auto g_x_x_0_0_x_zz_z_y = buffer_1100_pdpp[52];

    auto g_x_x_0_0_x_zz_z_z = buffer_1100_pdpp[53];

    auto g_x_x_0_0_y_xx_x_x = buffer_1100_pdpp[54];

    auto g_x_x_0_0_y_xx_x_y = buffer_1100_pdpp[55];

    auto g_x_x_0_0_y_xx_x_z = buffer_1100_pdpp[56];

    auto g_x_x_0_0_y_xx_y_x = buffer_1100_pdpp[57];

    auto g_x_x_0_0_y_xx_y_y = buffer_1100_pdpp[58];

    auto g_x_x_0_0_y_xx_y_z = buffer_1100_pdpp[59];

    auto g_x_x_0_0_y_xx_z_x = buffer_1100_pdpp[60];

    auto g_x_x_0_0_y_xx_z_y = buffer_1100_pdpp[61];

    auto g_x_x_0_0_y_xx_z_z = buffer_1100_pdpp[62];

    auto g_x_x_0_0_y_xy_x_x = buffer_1100_pdpp[63];

    auto g_x_x_0_0_y_xy_x_y = buffer_1100_pdpp[64];

    auto g_x_x_0_0_y_xy_x_z = buffer_1100_pdpp[65];

    auto g_x_x_0_0_y_xy_y_x = buffer_1100_pdpp[66];

    auto g_x_x_0_0_y_xy_y_y = buffer_1100_pdpp[67];

    auto g_x_x_0_0_y_xy_y_z = buffer_1100_pdpp[68];

    auto g_x_x_0_0_y_xy_z_x = buffer_1100_pdpp[69];

    auto g_x_x_0_0_y_xy_z_y = buffer_1100_pdpp[70];

    auto g_x_x_0_0_y_xy_z_z = buffer_1100_pdpp[71];

    auto g_x_x_0_0_y_xz_x_x = buffer_1100_pdpp[72];

    auto g_x_x_0_0_y_xz_x_y = buffer_1100_pdpp[73];

    auto g_x_x_0_0_y_xz_x_z = buffer_1100_pdpp[74];

    auto g_x_x_0_0_y_xz_y_x = buffer_1100_pdpp[75];

    auto g_x_x_0_0_y_xz_y_y = buffer_1100_pdpp[76];

    auto g_x_x_0_0_y_xz_y_z = buffer_1100_pdpp[77];

    auto g_x_x_0_0_y_xz_z_x = buffer_1100_pdpp[78];

    auto g_x_x_0_0_y_xz_z_y = buffer_1100_pdpp[79];

    auto g_x_x_0_0_y_xz_z_z = buffer_1100_pdpp[80];

    auto g_x_x_0_0_y_yy_x_x = buffer_1100_pdpp[81];

    auto g_x_x_0_0_y_yy_x_y = buffer_1100_pdpp[82];

    auto g_x_x_0_0_y_yy_x_z = buffer_1100_pdpp[83];

    auto g_x_x_0_0_y_yy_y_x = buffer_1100_pdpp[84];

    auto g_x_x_0_0_y_yy_y_y = buffer_1100_pdpp[85];

    auto g_x_x_0_0_y_yy_y_z = buffer_1100_pdpp[86];

    auto g_x_x_0_0_y_yy_z_x = buffer_1100_pdpp[87];

    auto g_x_x_0_0_y_yy_z_y = buffer_1100_pdpp[88];

    auto g_x_x_0_0_y_yy_z_z = buffer_1100_pdpp[89];

    auto g_x_x_0_0_y_yz_x_x = buffer_1100_pdpp[90];

    auto g_x_x_0_0_y_yz_x_y = buffer_1100_pdpp[91];

    auto g_x_x_0_0_y_yz_x_z = buffer_1100_pdpp[92];

    auto g_x_x_0_0_y_yz_y_x = buffer_1100_pdpp[93];

    auto g_x_x_0_0_y_yz_y_y = buffer_1100_pdpp[94];

    auto g_x_x_0_0_y_yz_y_z = buffer_1100_pdpp[95];

    auto g_x_x_0_0_y_yz_z_x = buffer_1100_pdpp[96];

    auto g_x_x_0_0_y_yz_z_y = buffer_1100_pdpp[97];

    auto g_x_x_0_0_y_yz_z_z = buffer_1100_pdpp[98];

    auto g_x_x_0_0_y_zz_x_x = buffer_1100_pdpp[99];

    auto g_x_x_0_0_y_zz_x_y = buffer_1100_pdpp[100];

    auto g_x_x_0_0_y_zz_x_z = buffer_1100_pdpp[101];

    auto g_x_x_0_0_y_zz_y_x = buffer_1100_pdpp[102];

    auto g_x_x_0_0_y_zz_y_y = buffer_1100_pdpp[103];

    auto g_x_x_0_0_y_zz_y_z = buffer_1100_pdpp[104];

    auto g_x_x_0_0_y_zz_z_x = buffer_1100_pdpp[105];

    auto g_x_x_0_0_y_zz_z_y = buffer_1100_pdpp[106];

    auto g_x_x_0_0_y_zz_z_z = buffer_1100_pdpp[107];

    auto g_x_x_0_0_z_xx_x_x = buffer_1100_pdpp[108];

    auto g_x_x_0_0_z_xx_x_y = buffer_1100_pdpp[109];

    auto g_x_x_0_0_z_xx_x_z = buffer_1100_pdpp[110];

    auto g_x_x_0_0_z_xx_y_x = buffer_1100_pdpp[111];

    auto g_x_x_0_0_z_xx_y_y = buffer_1100_pdpp[112];

    auto g_x_x_0_0_z_xx_y_z = buffer_1100_pdpp[113];

    auto g_x_x_0_0_z_xx_z_x = buffer_1100_pdpp[114];

    auto g_x_x_0_0_z_xx_z_y = buffer_1100_pdpp[115];

    auto g_x_x_0_0_z_xx_z_z = buffer_1100_pdpp[116];

    auto g_x_x_0_0_z_xy_x_x = buffer_1100_pdpp[117];

    auto g_x_x_0_0_z_xy_x_y = buffer_1100_pdpp[118];

    auto g_x_x_0_0_z_xy_x_z = buffer_1100_pdpp[119];

    auto g_x_x_0_0_z_xy_y_x = buffer_1100_pdpp[120];

    auto g_x_x_0_0_z_xy_y_y = buffer_1100_pdpp[121];

    auto g_x_x_0_0_z_xy_y_z = buffer_1100_pdpp[122];

    auto g_x_x_0_0_z_xy_z_x = buffer_1100_pdpp[123];

    auto g_x_x_0_0_z_xy_z_y = buffer_1100_pdpp[124];

    auto g_x_x_0_0_z_xy_z_z = buffer_1100_pdpp[125];

    auto g_x_x_0_0_z_xz_x_x = buffer_1100_pdpp[126];

    auto g_x_x_0_0_z_xz_x_y = buffer_1100_pdpp[127];

    auto g_x_x_0_0_z_xz_x_z = buffer_1100_pdpp[128];

    auto g_x_x_0_0_z_xz_y_x = buffer_1100_pdpp[129];

    auto g_x_x_0_0_z_xz_y_y = buffer_1100_pdpp[130];

    auto g_x_x_0_0_z_xz_y_z = buffer_1100_pdpp[131];

    auto g_x_x_0_0_z_xz_z_x = buffer_1100_pdpp[132];

    auto g_x_x_0_0_z_xz_z_y = buffer_1100_pdpp[133];

    auto g_x_x_0_0_z_xz_z_z = buffer_1100_pdpp[134];

    auto g_x_x_0_0_z_yy_x_x = buffer_1100_pdpp[135];

    auto g_x_x_0_0_z_yy_x_y = buffer_1100_pdpp[136];

    auto g_x_x_0_0_z_yy_x_z = buffer_1100_pdpp[137];

    auto g_x_x_0_0_z_yy_y_x = buffer_1100_pdpp[138];

    auto g_x_x_0_0_z_yy_y_y = buffer_1100_pdpp[139];

    auto g_x_x_0_0_z_yy_y_z = buffer_1100_pdpp[140];

    auto g_x_x_0_0_z_yy_z_x = buffer_1100_pdpp[141];

    auto g_x_x_0_0_z_yy_z_y = buffer_1100_pdpp[142];

    auto g_x_x_0_0_z_yy_z_z = buffer_1100_pdpp[143];

    auto g_x_x_0_0_z_yz_x_x = buffer_1100_pdpp[144];

    auto g_x_x_0_0_z_yz_x_y = buffer_1100_pdpp[145];

    auto g_x_x_0_0_z_yz_x_z = buffer_1100_pdpp[146];

    auto g_x_x_0_0_z_yz_y_x = buffer_1100_pdpp[147];

    auto g_x_x_0_0_z_yz_y_y = buffer_1100_pdpp[148];

    auto g_x_x_0_0_z_yz_y_z = buffer_1100_pdpp[149];

    auto g_x_x_0_0_z_yz_z_x = buffer_1100_pdpp[150];

    auto g_x_x_0_0_z_yz_z_y = buffer_1100_pdpp[151];

    auto g_x_x_0_0_z_yz_z_z = buffer_1100_pdpp[152];

    auto g_x_x_0_0_z_zz_x_x = buffer_1100_pdpp[153];

    auto g_x_x_0_0_z_zz_x_y = buffer_1100_pdpp[154];

    auto g_x_x_0_0_z_zz_x_z = buffer_1100_pdpp[155];

    auto g_x_x_0_0_z_zz_y_x = buffer_1100_pdpp[156];

    auto g_x_x_0_0_z_zz_y_y = buffer_1100_pdpp[157];

    auto g_x_x_0_0_z_zz_y_z = buffer_1100_pdpp[158];

    auto g_x_x_0_0_z_zz_z_x = buffer_1100_pdpp[159];

    auto g_x_x_0_0_z_zz_z_y = buffer_1100_pdpp[160];

    auto g_x_x_0_0_z_zz_z_z = buffer_1100_pdpp[161];

    auto g_x_y_0_0_x_xx_x_x = buffer_1100_pdpp[162];

    auto g_x_y_0_0_x_xx_x_y = buffer_1100_pdpp[163];

    auto g_x_y_0_0_x_xx_x_z = buffer_1100_pdpp[164];

    auto g_x_y_0_0_x_xx_y_x = buffer_1100_pdpp[165];

    auto g_x_y_0_0_x_xx_y_y = buffer_1100_pdpp[166];

    auto g_x_y_0_0_x_xx_y_z = buffer_1100_pdpp[167];

    auto g_x_y_0_0_x_xx_z_x = buffer_1100_pdpp[168];

    auto g_x_y_0_0_x_xx_z_y = buffer_1100_pdpp[169];

    auto g_x_y_0_0_x_xx_z_z = buffer_1100_pdpp[170];

    auto g_x_y_0_0_x_xy_x_x = buffer_1100_pdpp[171];

    auto g_x_y_0_0_x_xy_x_y = buffer_1100_pdpp[172];

    auto g_x_y_0_0_x_xy_x_z = buffer_1100_pdpp[173];

    auto g_x_y_0_0_x_xy_y_x = buffer_1100_pdpp[174];

    auto g_x_y_0_0_x_xy_y_y = buffer_1100_pdpp[175];

    auto g_x_y_0_0_x_xy_y_z = buffer_1100_pdpp[176];

    auto g_x_y_0_0_x_xy_z_x = buffer_1100_pdpp[177];

    auto g_x_y_0_0_x_xy_z_y = buffer_1100_pdpp[178];

    auto g_x_y_0_0_x_xy_z_z = buffer_1100_pdpp[179];

    auto g_x_y_0_0_x_xz_x_x = buffer_1100_pdpp[180];

    auto g_x_y_0_0_x_xz_x_y = buffer_1100_pdpp[181];

    auto g_x_y_0_0_x_xz_x_z = buffer_1100_pdpp[182];

    auto g_x_y_0_0_x_xz_y_x = buffer_1100_pdpp[183];

    auto g_x_y_0_0_x_xz_y_y = buffer_1100_pdpp[184];

    auto g_x_y_0_0_x_xz_y_z = buffer_1100_pdpp[185];

    auto g_x_y_0_0_x_xz_z_x = buffer_1100_pdpp[186];

    auto g_x_y_0_0_x_xz_z_y = buffer_1100_pdpp[187];

    auto g_x_y_0_0_x_xz_z_z = buffer_1100_pdpp[188];

    auto g_x_y_0_0_x_yy_x_x = buffer_1100_pdpp[189];

    auto g_x_y_0_0_x_yy_x_y = buffer_1100_pdpp[190];

    auto g_x_y_0_0_x_yy_x_z = buffer_1100_pdpp[191];

    auto g_x_y_0_0_x_yy_y_x = buffer_1100_pdpp[192];

    auto g_x_y_0_0_x_yy_y_y = buffer_1100_pdpp[193];

    auto g_x_y_0_0_x_yy_y_z = buffer_1100_pdpp[194];

    auto g_x_y_0_0_x_yy_z_x = buffer_1100_pdpp[195];

    auto g_x_y_0_0_x_yy_z_y = buffer_1100_pdpp[196];

    auto g_x_y_0_0_x_yy_z_z = buffer_1100_pdpp[197];

    auto g_x_y_0_0_x_yz_x_x = buffer_1100_pdpp[198];

    auto g_x_y_0_0_x_yz_x_y = buffer_1100_pdpp[199];

    auto g_x_y_0_0_x_yz_x_z = buffer_1100_pdpp[200];

    auto g_x_y_0_0_x_yz_y_x = buffer_1100_pdpp[201];

    auto g_x_y_0_0_x_yz_y_y = buffer_1100_pdpp[202];

    auto g_x_y_0_0_x_yz_y_z = buffer_1100_pdpp[203];

    auto g_x_y_0_0_x_yz_z_x = buffer_1100_pdpp[204];

    auto g_x_y_0_0_x_yz_z_y = buffer_1100_pdpp[205];

    auto g_x_y_0_0_x_yz_z_z = buffer_1100_pdpp[206];

    auto g_x_y_0_0_x_zz_x_x = buffer_1100_pdpp[207];

    auto g_x_y_0_0_x_zz_x_y = buffer_1100_pdpp[208];

    auto g_x_y_0_0_x_zz_x_z = buffer_1100_pdpp[209];

    auto g_x_y_0_0_x_zz_y_x = buffer_1100_pdpp[210];

    auto g_x_y_0_0_x_zz_y_y = buffer_1100_pdpp[211];

    auto g_x_y_0_0_x_zz_y_z = buffer_1100_pdpp[212];

    auto g_x_y_0_0_x_zz_z_x = buffer_1100_pdpp[213];

    auto g_x_y_0_0_x_zz_z_y = buffer_1100_pdpp[214];

    auto g_x_y_0_0_x_zz_z_z = buffer_1100_pdpp[215];

    auto g_x_y_0_0_y_xx_x_x = buffer_1100_pdpp[216];

    auto g_x_y_0_0_y_xx_x_y = buffer_1100_pdpp[217];

    auto g_x_y_0_0_y_xx_x_z = buffer_1100_pdpp[218];

    auto g_x_y_0_0_y_xx_y_x = buffer_1100_pdpp[219];

    auto g_x_y_0_0_y_xx_y_y = buffer_1100_pdpp[220];

    auto g_x_y_0_0_y_xx_y_z = buffer_1100_pdpp[221];

    auto g_x_y_0_0_y_xx_z_x = buffer_1100_pdpp[222];

    auto g_x_y_0_0_y_xx_z_y = buffer_1100_pdpp[223];

    auto g_x_y_0_0_y_xx_z_z = buffer_1100_pdpp[224];

    auto g_x_y_0_0_y_xy_x_x = buffer_1100_pdpp[225];

    auto g_x_y_0_0_y_xy_x_y = buffer_1100_pdpp[226];

    auto g_x_y_0_0_y_xy_x_z = buffer_1100_pdpp[227];

    auto g_x_y_0_0_y_xy_y_x = buffer_1100_pdpp[228];

    auto g_x_y_0_0_y_xy_y_y = buffer_1100_pdpp[229];

    auto g_x_y_0_0_y_xy_y_z = buffer_1100_pdpp[230];

    auto g_x_y_0_0_y_xy_z_x = buffer_1100_pdpp[231];

    auto g_x_y_0_0_y_xy_z_y = buffer_1100_pdpp[232];

    auto g_x_y_0_0_y_xy_z_z = buffer_1100_pdpp[233];

    auto g_x_y_0_0_y_xz_x_x = buffer_1100_pdpp[234];

    auto g_x_y_0_0_y_xz_x_y = buffer_1100_pdpp[235];

    auto g_x_y_0_0_y_xz_x_z = buffer_1100_pdpp[236];

    auto g_x_y_0_0_y_xz_y_x = buffer_1100_pdpp[237];

    auto g_x_y_0_0_y_xz_y_y = buffer_1100_pdpp[238];

    auto g_x_y_0_0_y_xz_y_z = buffer_1100_pdpp[239];

    auto g_x_y_0_0_y_xz_z_x = buffer_1100_pdpp[240];

    auto g_x_y_0_0_y_xz_z_y = buffer_1100_pdpp[241];

    auto g_x_y_0_0_y_xz_z_z = buffer_1100_pdpp[242];

    auto g_x_y_0_0_y_yy_x_x = buffer_1100_pdpp[243];

    auto g_x_y_0_0_y_yy_x_y = buffer_1100_pdpp[244];

    auto g_x_y_0_0_y_yy_x_z = buffer_1100_pdpp[245];

    auto g_x_y_0_0_y_yy_y_x = buffer_1100_pdpp[246];

    auto g_x_y_0_0_y_yy_y_y = buffer_1100_pdpp[247];

    auto g_x_y_0_0_y_yy_y_z = buffer_1100_pdpp[248];

    auto g_x_y_0_0_y_yy_z_x = buffer_1100_pdpp[249];

    auto g_x_y_0_0_y_yy_z_y = buffer_1100_pdpp[250];

    auto g_x_y_0_0_y_yy_z_z = buffer_1100_pdpp[251];

    auto g_x_y_0_0_y_yz_x_x = buffer_1100_pdpp[252];

    auto g_x_y_0_0_y_yz_x_y = buffer_1100_pdpp[253];

    auto g_x_y_0_0_y_yz_x_z = buffer_1100_pdpp[254];

    auto g_x_y_0_0_y_yz_y_x = buffer_1100_pdpp[255];

    auto g_x_y_0_0_y_yz_y_y = buffer_1100_pdpp[256];

    auto g_x_y_0_0_y_yz_y_z = buffer_1100_pdpp[257];

    auto g_x_y_0_0_y_yz_z_x = buffer_1100_pdpp[258];

    auto g_x_y_0_0_y_yz_z_y = buffer_1100_pdpp[259];

    auto g_x_y_0_0_y_yz_z_z = buffer_1100_pdpp[260];

    auto g_x_y_0_0_y_zz_x_x = buffer_1100_pdpp[261];

    auto g_x_y_0_0_y_zz_x_y = buffer_1100_pdpp[262];

    auto g_x_y_0_0_y_zz_x_z = buffer_1100_pdpp[263];

    auto g_x_y_0_0_y_zz_y_x = buffer_1100_pdpp[264];

    auto g_x_y_0_0_y_zz_y_y = buffer_1100_pdpp[265];

    auto g_x_y_0_0_y_zz_y_z = buffer_1100_pdpp[266];

    auto g_x_y_0_0_y_zz_z_x = buffer_1100_pdpp[267];

    auto g_x_y_0_0_y_zz_z_y = buffer_1100_pdpp[268];

    auto g_x_y_0_0_y_zz_z_z = buffer_1100_pdpp[269];

    auto g_x_y_0_0_z_xx_x_x = buffer_1100_pdpp[270];

    auto g_x_y_0_0_z_xx_x_y = buffer_1100_pdpp[271];

    auto g_x_y_0_0_z_xx_x_z = buffer_1100_pdpp[272];

    auto g_x_y_0_0_z_xx_y_x = buffer_1100_pdpp[273];

    auto g_x_y_0_0_z_xx_y_y = buffer_1100_pdpp[274];

    auto g_x_y_0_0_z_xx_y_z = buffer_1100_pdpp[275];

    auto g_x_y_0_0_z_xx_z_x = buffer_1100_pdpp[276];

    auto g_x_y_0_0_z_xx_z_y = buffer_1100_pdpp[277];

    auto g_x_y_0_0_z_xx_z_z = buffer_1100_pdpp[278];

    auto g_x_y_0_0_z_xy_x_x = buffer_1100_pdpp[279];

    auto g_x_y_0_0_z_xy_x_y = buffer_1100_pdpp[280];

    auto g_x_y_0_0_z_xy_x_z = buffer_1100_pdpp[281];

    auto g_x_y_0_0_z_xy_y_x = buffer_1100_pdpp[282];

    auto g_x_y_0_0_z_xy_y_y = buffer_1100_pdpp[283];

    auto g_x_y_0_0_z_xy_y_z = buffer_1100_pdpp[284];

    auto g_x_y_0_0_z_xy_z_x = buffer_1100_pdpp[285];

    auto g_x_y_0_0_z_xy_z_y = buffer_1100_pdpp[286];

    auto g_x_y_0_0_z_xy_z_z = buffer_1100_pdpp[287];

    auto g_x_y_0_0_z_xz_x_x = buffer_1100_pdpp[288];

    auto g_x_y_0_0_z_xz_x_y = buffer_1100_pdpp[289];

    auto g_x_y_0_0_z_xz_x_z = buffer_1100_pdpp[290];

    auto g_x_y_0_0_z_xz_y_x = buffer_1100_pdpp[291];

    auto g_x_y_0_0_z_xz_y_y = buffer_1100_pdpp[292];

    auto g_x_y_0_0_z_xz_y_z = buffer_1100_pdpp[293];

    auto g_x_y_0_0_z_xz_z_x = buffer_1100_pdpp[294];

    auto g_x_y_0_0_z_xz_z_y = buffer_1100_pdpp[295];

    auto g_x_y_0_0_z_xz_z_z = buffer_1100_pdpp[296];

    auto g_x_y_0_0_z_yy_x_x = buffer_1100_pdpp[297];

    auto g_x_y_0_0_z_yy_x_y = buffer_1100_pdpp[298];

    auto g_x_y_0_0_z_yy_x_z = buffer_1100_pdpp[299];

    auto g_x_y_0_0_z_yy_y_x = buffer_1100_pdpp[300];

    auto g_x_y_0_0_z_yy_y_y = buffer_1100_pdpp[301];

    auto g_x_y_0_0_z_yy_y_z = buffer_1100_pdpp[302];

    auto g_x_y_0_0_z_yy_z_x = buffer_1100_pdpp[303];

    auto g_x_y_0_0_z_yy_z_y = buffer_1100_pdpp[304];

    auto g_x_y_0_0_z_yy_z_z = buffer_1100_pdpp[305];

    auto g_x_y_0_0_z_yz_x_x = buffer_1100_pdpp[306];

    auto g_x_y_0_0_z_yz_x_y = buffer_1100_pdpp[307];

    auto g_x_y_0_0_z_yz_x_z = buffer_1100_pdpp[308];

    auto g_x_y_0_0_z_yz_y_x = buffer_1100_pdpp[309];

    auto g_x_y_0_0_z_yz_y_y = buffer_1100_pdpp[310];

    auto g_x_y_0_0_z_yz_y_z = buffer_1100_pdpp[311];

    auto g_x_y_0_0_z_yz_z_x = buffer_1100_pdpp[312];

    auto g_x_y_0_0_z_yz_z_y = buffer_1100_pdpp[313];

    auto g_x_y_0_0_z_yz_z_z = buffer_1100_pdpp[314];

    auto g_x_y_0_0_z_zz_x_x = buffer_1100_pdpp[315];

    auto g_x_y_0_0_z_zz_x_y = buffer_1100_pdpp[316];

    auto g_x_y_0_0_z_zz_x_z = buffer_1100_pdpp[317];

    auto g_x_y_0_0_z_zz_y_x = buffer_1100_pdpp[318];

    auto g_x_y_0_0_z_zz_y_y = buffer_1100_pdpp[319];

    auto g_x_y_0_0_z_zz_y_z = buffer_1100_pdpp[320];

    auto g_x_y_0_0_z_zz_z_x = buffer_1100_pdpp[321];

    auto g_x_y_0_0_z_zz_z_y = buffer_1100_pdpp[322];

    auto g_x_y_0_0_z_zz_z_z = buffer_1100_pdpp[323];

    auto g_x_z_0_0_x_xx_x_x = buffer_1100_pdpp[324];

    auto g_x_z_0_0_x_xx_x_y = buffer_1100_pdpp[325];

    auto g_x_z_0_0_x_xx_x_z = buffer_1100_pdpp[326];

    auto g_x_z_0_0_x_xx_y_x = buffer_1100_pdpp[327];

    auto g_x_z_0_0_x_xx_y_y = buffer_1100_pdpp[328];

    auto g_x_z_0_0_x_xx_y_z = buffer_1100_pdpp[329];

    auto g_x_z_0_0_x_xx_z_x = buffer_1100_pdpp[330];

    auto g_x_z_0_0_x_xx_z_y = buffer_1100_pdpp[331];

    auto g_x_z_0_0_x_xx_z_z = buffer_1100_pdpp[332];

    auto g_x_z_0_0_x_xy_x_x = buffer_1100_pdpp[333];

    auto g_x_z_0_0_x_xy_x_y = buffer_1100_pdpp[334];

    auto g_x_z_0_0_x_xy_x_z = buffer_1100_pdpp[335];

    auto g_x_z_0_0_x_xy_y_x = buffer_1100_pdpp[336];

    auto g_x_z_0_0_x_xy_y_y = buffer_1100_pdpp[337];

    auto g_x_z_0_0_x_xy_y_z = buffer_1100_pdpp[338];

    auto g_x_z_0_0_x_xy_z_x = buffer_1100_pdpp[339];

    auto g_x_z_0_0_x_xy_z_y = buffer_1100_pdpp[340];

    auto g_x_z_0_0_x_xy_z_z = buffer_1100_pdpp[341];

    auto g_x_z_0_0_x_xz_x_x = buffer_1100_pdpp[342];

    auto g_x_z_0_0_x_xz_x_y = buffer_1100_pdpp[343];

    auto g_x_z_0_0_x_xz_x_z = buffer_1100_pdpp[344];

    auto g_x_z_0_0_x_xz_y_x = buffer_1100_pdpp[345];

    auto g_x_z_0_0_x_xz_y_y = buffer_1100_pdpp[346];

    auto g_x_z_0_0_x_xz_y_z = buffer_1100_pdpp[347];

    auto g_x_z_0_0_x_xz_z_x = buffer_1100_pdpp[348];

    auto g_x_z_0_0_x_xz_z_y = buffer_1100_pdpp[349];

    auto g_x_z_0_0_x_xz_z_z = buffer_1100_pdpp[350];

    auto g_x_z_0_0_x_yy_x_x = buffer_1100_pdpp[351];

    auto g_x_z_0_0_x_yy_x_y = buffer_1100_pdpp[352];

    auto g_x_z_0_0_x_yy_x_z = buffer_1100_pdpp[353];

    auto g_x_z_0_0_x_yy_y_x = buffer_1100_pdpp[354];

    auto g_x_z_0_0_x_yy_y_y = buffer_1100_pdpp[355];

    auto g_x_z_0_0_x_yy_y_z = buffer_1100_pdpp[356];

    auto g_x_z_0_0_x_yy_z_x = buffer_1100_pdpp[357];

    auto g_x_z_0_0_x_yy_z_y = buffer_1100_pdpp[358];

    auto g_x_z_0_0_x_yy_z_z = buffer_1100_pdpp[359];

    auto g_x_z_0_0_x_yz_x_x = buffer_1100_pdpp[360];

    auto g_x_z_0_0_x_yz_x_y = buffer_1100_pdpp[361];

    auto g_x_z_0_0_x_yz_x_z = buffer_1100_pdpp[362];

    auto g_x_z_0_0_x_yz_y_x = buffer_1100_pdpp[363];

    auto g_x_z_0_0_x_yz_y_y = buffer_1100_pdpp[364];

    auto g_x_z_0_0_x_yz_y_z = buffer_1100_pdpp[365];

    auto g_x_z_0_0_x_yz_z_x = buffer_1100_pdpp[366];

    auto g_x_z_0_0_x_yz_z_y = buffer_1100_pdpp[367];

    auto g_x_z_0_0_x_yz_z_z = buffer_1100_pdpp[368];

    auto g_x_z_0_0_x_zz_x_x = buffer_1100_pdpp[369];

    auto g_x_z_0_0_x_zz_x_y = buffer_1100_pdpp[370];

    auto g_x_z_0_0_x_zz_x_z = buffer_1100_pdpp[371];

    auto g_x_z_0_0_x_zz_y_x = buffer_1100_pdpp[372];

    auto g_x_z_0_0_x_zz_y_y = buffer_1100_pdpp[373];

    auto g_x_z_0_0_x_zz_y_z = buffer_1100_pdpp[374];

    auto g_x_z_0_0_x_zz_z_x = buffer_1100_pdpp[375];

    auto g_x_z_0_0_x_zz_z_y = buffer_1100_pdpp[376];

    auto g_x_z_0_0_x_zz_z_z = buffer_1100_pdpp[377];

    auto g_x_z_0_0_y_xx_x_x = buffer_1100_pdpp[378];

    auto g_x_z_0_0_y_xx_x_y = buffer_1100_pdpp[379];

    auto g_x_z_0_0_y_xx_x_z = buffer_1100_pdpp[380];

    auto g_x_z_0_0_y_xx_y_x = buffer_1100_pdpp[381];

    auto g_x_z_0_0_y_xx_y_y = buffer_1100_pdpp[382];

    auto g_x_z_0_0_y_xx_y_z = buffer_1100_pdpp[383];

    auto g_x_z_0_0_y_xx_z_x = buffer_1100_pdpp[384];

    auto g_x_z_0_0_y_xx_z_y = buffer_1100_pdpp[385];

    auto g_x_z_0_0_y_xx_z_z = buffer_1100_pdpp[386];

    auto g_x_z_0_0_y_xy_x_x = buffer_1100_pdpp[387];

    auto g_x_z_0_0_y_xy_x_y = buffer_1100_pdpp[388];

    auto g_x_z_0_0_y_xy_x_z = buffer_1100_pdpp[389];

    auto g_x_z_0_0_y_xy_y_x = buffer_1100_pdpp[390];

    auto g_x_z_0_0_y_xy_y_y = buffer_1100_pdpp[391];

    auto g_x_z_0_0_y_xy_y_z = buffer_1100_pdpp[392];

    auto g_x_z_0_0_y_xy_z_x = buffer_1100_pdpp[393];

    auto g_x_z_0_0_y_xy_z_y = buffer_1100_pdpp[394];

    auto g_x_z_0_0_y_xy_z_z = buffer_1100_pdpp[395];

    auto g_x_z_0_0_y_xz_x_x = buffer_1100_pdpp[396];

    auto g_x_z_0_0_y_xz_x_y = buffer_1100_pdpp[397];

    auto g_x_z_0_0_y_xz_x_z = buffer_1100_pdpp[398];

    auto g_x_z_0_0_y_xz_y_x = buffer_1100_pdpp[399];

    auto g_x_z_0_0_y_xz_y_y = buffer_1100_pdpp[400];

    auto g_x_z_0_0_y_xz_y_z = buffer_1100_pdpp[401];

    auto g_x_z_0_0_y_xz_z_x = buffer_1100_pdpp[402];

    auto g_x_z_0_0_y_xz_z_y = buffer_1100_pdpp[403];

    auto g_x_z_0_0_y_xz_z_z = buffer_1100_pdpp[404];

    auto g_x_z_0_0_y_yy_x_x = buffer_1100_pdpp[405];

    auto g_x_z_0_0_y_yy_x_y = buffer_1100_pdpp[406];

    auto g_x_z_0_0_y_yy_x_z = buffer_1100_pdpp[407];

    auto g_x_z_0_0_y_yy_y_x = buffer_1100_pdpp[408];

    auto g_x_z_0_0_y_yy_y_y = buffer_1100_pdpp[409];

    auto g_x_z_0_0_y_yy_y_z = buffer_1100_pdpp[410];

    auto g_x_z_0_0_y_yy_z_x = buffer_1100_pdpp[411];

    auto g_x_z_0_0_y_yy_z_y = buffer_1100_pdpp[412];

    auto g_x_z_0_0_y_yy_z_z = buffer_1100_pdpp[413];

    auto g_x_z_0_0_y_yz_x_x = buffer_1100_pdpp[414];

    auto g_x_z_0_0_y_yz_x_y = buffer_1100_pdpp[415];

    auto g_x_z_0_0_y_yz_x_z = buffer_1100_pdpp[416];

    auto g_x_z_0_0_y_yz_y_x = buffer_1100_pdpp[417];

    auto g_x_z_0_0_y_yz_y_y = buffer_1100_pdpp[418];

    auto g_x_z_0_0_y_yz_y_z = buffer_1100_pdpp[419];

    auto g_x_z_0_0_y_yz_z_x = buffer_1100_pdpp[420];

    auto g_x_z_0_0_y_yz_z_y = buffer_1100_pdpp[421];

    auto g_x_z_0_0_y_yz_z_z = buffer_1100_pdpp[422];

    auto g_x_z_0_0_y_zz_x_x = buffer_1100_pdpp[423];

    auto g_x_z_0_0_y_zz_x_y = buffer_1100_pdpp[424];

    auto g_x_z_0_0_y_zz_x_z = buffer_1100_pdpp[425];

    auto g_x_z_0_0_y_zz_y_x = buffer_1100_pdpp[426];

    auto g_x_z_0_0_y_zz_y_y = buffer_1100_pdpp[427];

    auto g_x_z_0_0_y_zz_y_z = buffer_1100_pdpp[428];

    auto g_x_z_0_0_y_zz_z_x = buffer_1100_pdpp[429];

    auto g_x_z_0_0_y_zz_z_y = buffer_1100_pdpp[430];

    auto g_x_z_0_0_y_zz_z_z = buffer_1100_pdpp[431];

    auto g_x_z_0_0_z_xx_x_x = buffer_1100_pdpp[432];

    auto g_x_z_0_0_z_xx_x_y = buffer_1100_pdpp[433];

    auto g_x_z_0_0_z_xx_x_z = buffer_1100_pdpp[434];

    auto g_x_z_0_0_z_xx_y_x = buffer_1100_pdpp[435];

    auto g_x_z_0_0_z_xx_y_y = buffer_1100_pdpp[436];

    auto g_x_z_0_0_z_xx_y_z = buffer_1100_pdpp[437];

    auto g_x_z_0_0_z_xx_z_x = buffer_1100_pdpp[438];

    auto g_x_z_0_0_z_xx_z_y = buffer_1100_pdpp[439];

    auto g_x_z_0_0_z_xx_z_z = buffer_1100_pdpp[440];

    auto g_x_z_0_0_z_xy_x_x = buffer_1100_pdpp[441];

    auto g_x_z_0_0_z_xy_x_y = buffer_1100_pdpp[442];

    auto g_x_z_0_0_z_xy_x_z = buffer_1100_pdpp[443];

    auto g_x_z_0_0_z_xy_y_x = buffer_1100_pdpp[444];

    auto g_x_z_0_0_z_xy_y_y = buffer_1100_pdpp[445];

    auto g_x_z_0_0_z_xy_y_z = buffer_1100_pdpp[446];

    auto g_x_z_0_0_z_xy_z_x = buffer_1100_pdpp[447];

    auto g_x_z_0_0_z_xy_z_y = buffer_1100_pdpp[448];

    auto g_x_z_0_0_z_xy_z_z = buffer_1100_pdpp[449];

    auto g_x_z_0_0_z_xz_x_x = buffer_1100_pdpp[450];

    auto g_x_z_0_0_z_xz_x_y = buffer_1100_pdpp[451];

    auto g_x_z_0_0_z_xz_x_z = buffer_1100_pdpp[452];

    auto g_x_z_0_0_z_xz_y_x = buffer_1100_pdpp[453];

    auto g_x_z_0_0_z_xz_y_y = buffer_1100_pdpp[454];

    auto g_x_z_0_0_z_xz_y_z = buffer_1100_pdpp[455];

    auto g_x_z_0_0_z_xz_z_x = buffer_1100_pdpp[456];

    auto g_x_z_0_0_z_xz_z_y = buffer_1100_pdpp[457];

    auto g_x_z_0_0_z_xz_z_z = buffer_1100_pdpp[458];

    auto g_x_z_0_0_z_yy_x_x = buffer_1100_pdpp[459];

    auto g_x_z_0_0_z_yy_x_y = buffer_1100_pdpp[460];

    auto g_x_z_0_0_z_yy_x_z = buffer_1100_pdpp[461];

    auto g_x_z_0_0_z_yy_y_x = buffer_1100_pdpp[462];

    auto g_x_z_0_0_z_yy_y_y = buffer_1100_pdpp[463];

    auto g_x_z_0_0_z_yy_y_z = buffer_1100_pdpp[464];

    auto g_x_z_0_0_z_yy_z_x = buffer_1100_pdpp[465];

    auto g_x_z_0_0_z_yy_z_y = buffer_1100_pdpp[466];

    auto g_x_z_0_0_z_yy_z_z = buffer_1100_pdpp[467];

    auto g_x_z_0_0_z_yz_x_x = buffer_1100_pdpp[468];

    auto g_x_z_0_0_z_yz_x_y = buffer_1100_pdpp[469];

    auto g_x_z_0_0_z_yz_x_z = buffer_1100_pdpp[470];

    auto g_x_z_0_0_z_yz_y_x = buffer_1100_pdpp[471];

    auto g_x_z_0_0_z_yz_y_y = buffer_1100_pdpp[472];

    auto g_x_z_0_0_z_yz_y_z = buffer_1100_pdpp[473];

    auto g_x_z_0_0_z_yz_z_x = buffer_1100_pdpp[474];

    auto g_x_z_0_0_z_yz_z_y = buffer_1100_pdpp[475];

    auto g_x_z_0_0_z_yz_z_z = buffer_1100_pdpp[476];

    auto g_x_z_0_0_z_zz_x_x = buffer_1100_pdpp[477];

    auto g_x_z_0_0_z_zz_x_y = buffer_1100_pdpp[478];

    auto g_x_z_0_0_z_zz_x_z = buffer_1100_pdpp[479];

    auto g_x_z_0_0_z_zz_y_x = buffer_1100_pdpp[480];

    auto g_x_z_0_0_z_zz_y_y = buffer_1100_pdpp[481];

    auto g_x_z_0_0_z_zz_y_z = buffer_1100_pdpp[482];

    auto g_x_z_0_0_z_zz_z_x = buffer_1100_pdpp[483];

    auto g_x_z_0_0_z_zz_z_y = buffer_1100_pdpp[484];

    auto g_x_z_0_0_z_zz_z_z = buffer_1100_pdpp[485];

    auto g_y_x_0_0_x_xx_x_x = buffer_1100_pdpp[486];

    auto g_y_x_0_0_x_xx_x_y = buffer_1100_pdpp[487];

    auto g_y_x_0_0_x_xx_x_z = buffer_1100_pdpp[488];

    auto g_y_x_0_0_x_xx_y_x = buffer_1100_pdpp[489];

    auto g_y_x_0_0_x_xx_y_y = buffer_1100_pdpp[490];

    auto g_y_x_0_0_x_xx_y_z = buffer_1100_pdpp[491];

    auto g_y_x_0_0_x_xx_z_x = buffer_1100_pdpp[492];

    auto g_y_x_0_0_x_xx_z_y = buffer_1100_pdpp[493];

    auto g_y_x_0_0_x_xx_z_z = buffer_1100_pdpp[494];

    auto g_y_x_0_0_x_xy_x_x = buffer_1100_pdpp[495];

    auto g_y_x_0_0_x_xy_x_y = buffer_1100_pdpp[496];

    auto g_y_x_0_0_x_xy_x_z = buffer_1100_pdpp[497];

    auto g_y_x_0_0_x_xy_y_x = buffer_1100_pdpp[498];

    auto g_y_x_0_0_x_xy_y_y = buffer_1100_pdpp[499];

    auto g_y_x_0_0_x_xy_y_z = buffer_1100_pdpp[500];

    auto g_y_x_0_0_x_xy_z_x = buffer_1100_pdpp[501];

    auto g_y_x_0_0_x_xy_z_y = buffer_1100_pdpp[502];

    auto g_y_x_0_0_x_xy_z_z = buffer_1100_pdpp[503];

    auto g_y_x_0_0_x_xz_x_x = buffer_1100_pdpp[504];

    auto g_y_x_0_0_x_xz_x_y = buffer_1100_pdpp[505];

    auto g_y_x_0_0_x_xz_x_z = buffer_1100_pdpp[506];

    auto g_y_x_0_0_x_xz_y_x = buffer_1100_pdpp[507];

    auto g_y_x_0_0_x_xz_y_y = buffer_1100_pdpp[508];

    auto g_y_x_0_0_x_xz_y_z = buffer_1100_pdpp[509];

    auto g_y_x_0_0_x_xz_z_x = buffer_1100_pdpp[510];

    auto g_y_x_0_0_x_xz_z_y = buffer_1100_pdpp[511];

    auto g_y_x_0_0_x_xz_z_z = buffer_1100_pdpp[512];

    auto g_y_x_0_0_x_yy_x_x = buffer_1100_pdpp[513];

    auto g_y_x_0_0_x_yy_x_y = buffer_1100_pdpp[514];

    auto g_y_x_0_0_x_yy_x_z = buffer_1100_pdpp[515];

    auto g_y_x_0_0_x_yy_y_x = buffer_1100_pdpp[516];

    auto g_y_x_0_0_x_yy_y_y = buffer_1100_pdpp[517];

    auto g_y_x_0_0_x_yy_y_z = buffer_1100_pdpp[518];

    auto g_y_x_0_0_x_yy_z_x = buffer_1100_pdpp[519];

    auto g_y_x_0_0_x_yy_z_y = buffer_1100_pdpp[520];

    auto g_y_x_0_0_x_yy_z_z = buffer_1100_pdpp[521];

    auto g_y_x_0_0_x_yz_x_x = buffer_1100_pdpp[522];

    auto g_y_x_0_0_x_yz_x_y = buffer_1100_pdpp[523];

    auto g_y_x_0_0_x_yz_x_z = buffer_1100_pdpp[524];

    auto g_y_x_0_0_x_yz_y_x = buffer_1100_pdpp[525];

    auto g_y_x_0_0_x_yz_y_y = buffer_1100_pdpp[526];

    auto g_y_x_0_0_x_yz_y_z = buffer_1100_pdpp[527];

    auto g_y_x_0_0_x_yz_z_x = buffer_1100_pdpp[528];

    auto g_y_x_0_0_x_yz_z_y = buffer_1100_pdpp[529];

    auto g_y_x_0_0_x_yz_z_z = buffer_1100_pdpp[530];

    auto g_y_x_0_0_x_zz_x_x = buffer_1100_pdpp[531];

    auto g_y_x_0_0_x_zz_x_y = buffer_1100_pdpp[532];

    auto g_y_x_0_0_x_zz_x_z = buffer_1100_pdpp[533];

    auto g_y_x_0_0_x_zz_y_x = buffer_1100_pdpp[534];

    auto g_y_x_0_0_x_zz_y_y = buffer_1100_pdpp[535];

    auto g_y_x_0_0_x_zz_y_z = buffer_1100_pdpp[536];

    auto g_y_x_0_0_x_zz_z_x = buffer_1100_pdpp[537];

    auto g_y_x_0_0_x_zz_z_y = buffer_1100_pdpp[538];

    auto g_y_x_0_0_x_zz_z_z = buffer_1100_pdpp[539];

    auto g_y_x_0_0_y_xx_x_x = buffer_1100_pdpp[540];

    auto g_y_x_0_0_y_xx_x_y = buffer_1100_pdpp[541];

    auto g_y_x_0_0_y_xx_x_z = buffer_1100_pdpp[542];

    auto g_y_x_0_0_y_xx_y_x = buffer_1100_pdpp[543];

    auto g_y_x_0_0_y_xx_y_y = buffer_1100_pdpp[544];

    auto g_y_x_0_0_y_xx_y_z = buffer_1100_pdpp[545];

    auto g_y_x_0_0_y_xx_z_x = buffer_1100_pdpp[546];

    auto g_y_x_0_0_y_xx_z_y = buffer_1100_pdpp[547];

    auto g_y_x_0_0_y_xx_z_z = buffer_1100_pdpp[548];

    auto g_y_x_0_0_y_xy_x_x = buffer_1100_pdpp[549];

    auto g_y_x_0_0_y_xy_x_y = buffer_1100_pdpp[550];

    auto g_y_x_0_0_y_xy_x_z = buffer_1100_pdpp[551];

    auto g_y_x_0_0_y_xy_y_x = buffer_1100_pdpp[552];

    auto g_y_x_0_0_y_xy_y_y = buffer_1100_pdpp[553];

    auto g_y_x_0_0_y_xy_y_z = buffer_1100_pdpp[554];

    auto g_y_x_0_0_y_xy_z_x = buffer_1100_pdpp[555];

    auto g_y_x_0_0_y_xy_z_y = buffer_1100_pdpp[556];

    auto g_y_x_0_0_y_xy_z_z = buffer_1100_pdpp[557];

    auto g_y_x_0_0_y_xz_x_x = buffer_1100_pdpp[558];

    auto g_y_x_0_0_y_xz_x_y = buffer_1100_pdpp[559];

    auto g_y_x_0_0_y_xz_x_z = buffer_1100_pdpp[560];

    auto g_y_x_0_0_y_xz_y_x = buffer_1100_pdpp[561];

    auto g_y_x_0_0_y_xz_y_y = buffer_1100_pdpp[562];

    auto g_y_x_0_0_y_xz_y_z = buffer_1100_pdpp[563];

    auto g_y_x_0_0_y_xz_z_x = buffer_1100_pdpp[564];

    auto g_y_x_0_0_y_xz_z_y = buffer_1100_pdpp[565];

    auto g_y_x_0_0_y_xz_z_z = buffer_1100_pdpp[566];

    auto g_y_x_0_0_y_yy_x_x = buffer_1100_pdpp[567];

    auto g_y_x_0_0_y_yy_x_y = buffer_1100_pdpp[568];

    auto g_y_x_0_0_y_yy_x_z = buffer_1100_pdpp[569];

    auto g_y_x_0_0_y_yy_y_x = buffer_1100_pdpp[570];

    auto g_y_x_0_0_y_yy_y_y = buffer_1100_pdpp[571];

    auto g_y_x_0_0_y_yy_y_z = buffer_1100_pdpp[572];

    auto g_y_x_0_0_y_yy_z_x = buffer_1100_pdpp[573];

    auto g_y_x_0_0_y_yy_z_y = buffer_1100_pdpp[574];

    auto g_y_x_0_0_y_yy_z_z = buffer_1100_pdpp[575];

    auto g_y_x_0_0_y_yz_x_x = buffer_1100_pdpp[576];

    auto g_y_x_0_0_y_yz_x_y = buffer_1100_pdpp[577];

    auto g_y_x_0_0_y_yz_x_z = buffer_1100_pdpp[578];

    auto g_y_x_0_0_y_yz_y_x = buffer_1100_pdpp[579];

    auto g_y_x_0_0_y_yz_y_y = buffer_1100_pdpp[580];

    auto g_y_x_0_0_y_yz_y_z = buffer_1100_pdpp[581];

    auto g_y_x_0_0_y_yz_z_x = buffer_1100_pdpp[582];

    auto g_y_x_0_0_y_yz_z_y = buffer_1100_pdpp[583];

    auto g_y_x_0_0_y_yz_z_z = buffer_1100_pdpp[584];

    auto g_y_x_0_0_y_zz_x_x = buffer_1100_pdpp[585];

    auto g_y_x_0_0_y_zz_x_y = buffer_1100_pdpp[586];

    auto g_y_x_0_0_y_zz_x_z = buffer_1100_pdpp[587];

    auto g_y_x_0_0_y_zz_y_x = buffer_1100_pdpp[588];

    auto g_y_x_0_0_y_zz_y_y = buffer_1100_pdpp[589];

    auto g_y_x_0_0_y_zz_y_z = buffer_1100_pdpp[590];

    auto g_y_x_0_0_y_zz_z_x = buffer_1100_pdpp[591];

    auto g_y_x_0_0_y_zz_z_y = buffer_1100_pdpp[592];

    auto g_y_x_0_0_y_zz_z_z = buffer_1100_pdpp[593];

    auto g_y_x_0_0_z_xx_x_x = buffer_1100_pdpp[594];

    auto g_y_x_0_0_z_xx_x_y = buffer_1100_pdpp[595];

    auto g_y_x_0_0_z_xx_x_z = buffer_1100_pdpp[596];

    auto g_y_x_0_0_z_xx_y_x = buffer_1100_pdpp[597];

    auto g_y_x_0_0_z_xx_y_y = buffer_1100_pdpp[598];

    auto g_y_x_0_0_z_xx_y_z = buffer_1100_pdpp[599];

    auto g_y_x_0_0_z_xx_z_x = buffer_1100_pdpp[600];

    auto g_y_x_0_0_z_xx_z_y = buffer_1100_pdpp[601];

    auto g_y_x_0_0_z_xx_z_z = buffer_1100_pdpp[602];

    auto g_y_x_0_0_z_xy_x_x = buffer_1100_pdpp[603];

    auto g_y_x_0_0_z_xy_x_y = buffer_1100_pdpp[604];

    auto g_y_x_0_0_z_xy_x_z = buffer_1100_pdpp[605];

    auto g_y_x_0_0_z_xy_y_x = buffer_1100_pdpp[606];

    auto g_y_x_0_0_z_xy_y_y = buffer_1100_pdpp[607];

    auto g_y_x_0_0_z_xy_y_z = buffer_1100_pdpp[608];

    auto g_y_x_0_0_z_xy_z_x = buffer_1100_pdpp[609];

    auto g_y_x_0_0_z_xy_z_y = buffer_1100_pdpp[610];

    auto g_y_x_0_0_z_xy_z_z = buffer_1100_pdpp[611];

    auto g_y_x_0_0_z_xz_x_x = buffer_1100_pdpp[612];

    auto g_y_x_0_0_z_xz_x_y = buffer_1100_pdpp[613];

    auto g_y_x_0_0_z_xz_x_z = buffer_1100_pdpp[614];

    auto g_y_x_0_0_z_xz_y_x = buffer_1100_pdpp[615];

    auto g_y_x_0_0_z_xz_y_y = buffer_1100_pdpp[616];

    auto g_y_x_0_0_z_xz_y_z = buffer_1100_pdpp[617];

    auto g_y_x_0_0_z_xz_z_x = buffer_1100_pdpp[618];

    auto g_y_x_0_0_z_xz_z_y = buffer_1100_pdpp[619];

    auto g_y_x_0_0_z_xz_z_z = buffer_1100_pdpp[620];

    auto g_y_x_0_0_z_yy_x_x = buffer_1100_pdpp[621];

    auto g_y_x_0_0_z_yy_x_y = buffer_1100_pdpp[622];

    auto g_y_x_0_0_z_yy_x_z = buffer_1100_pdpp[623];

    auto g_y_x_0_0_z_yy_y_x = buffer_1100_pdpp[624];

    auto g_y_x_0_0_z_yy_y_y = buffer_1100_pdpp[625];

    auto g_y_x_0_0_z_yy_y_z = buffer_1100_pdpp[626];

    auto g_y_x_0_0_z_yy_z_x = buffer_1100_pdpp[627];

    auto g_y_x_0_0_z_yy_z_y = buffer_1100_pdpp[628];

    auto g_y_x_0_0_z_yy_z_z = buffer_1100_pdpp[629];

    auto g_y_x_0_0_z_yz_x_x = buffer_1100_pdpp[630];

    auto g_y_x_0_0_z_yz_x_y = buffer_1100_pdpp[631];

    auto g_y_x_0_0_z_yz_x_z = buffer_1100_pdpp[632];

    auto g_y_x_0_0_z_yz_y_x = buffer_1100_pdpp[633];

    auto g_y_x_0_0_z_yz_y_y = buffer_1100_pdpp[634];

    auto g_y_x_0_0_z_yz_y_z = buffer_1100_pdpp[635];

    auto g_y_x_0_0_z_yz_z_x = buffer_1100_pdpp[636];

    auto g_y_x_0_0_z_yz_z_y = buffer_1100_pdpp[637];

    auto g_y_x_0_0_z_yz_z_z = buffer_1100_pdpp[638];

    auto g_y_x_0_0_z_zz_x_x = buffer_1100_pdpp[639];

    auto g_y_x_0_0_z_zz_x_y = buffer_1100_pdpp[640];

    auto g_y_x_0_0_z_zz_x_z = buffer_1100_pdpp[641];

    auto g_y_x_0_0_z_zz_y_x = buffer_1100_pdpp[642];

    auto g_y_x_0_0_z_zz_y_y = buffer_1100_pdpp[643];

    auto g_y_x_0_0_z_zz_y_z = buffer_1100_pdpp[644];

    auto g_y_x_0_0_z_zz_z_x = buffer_1100_pdpp[645];

    auto g_y_x_0_0_z_zz_z_y = buffer_1100_pdpp[646];

    auto g_y_x_0_0_z_zz_z_z = buffer_1100_pdpp[647];

    auto g_y_y_0_0_x_xx_x_x = buffer_1100_pdpp[648];

    auto g_y_y_0_0_x_xx_x_y = buffer_1100_pdpp[649];

    auto g_y_y_0_0_x_xx_x_z = buffer_1100_pdpp[650];

    auto g_y_y_0_0_x_xx_y_x = buffer_1100_pdpp[651];

    auto g_y_y_0_0_x_xx_y_y = buffer_1100_pdpp[652];

    auto g_y_y_0_0_x_xx_y_z = buffer_1100_pdpp[653];

    auto g_y_y_0_0_x_xx_z_x = buffer_1100_pdpp[654];

    auto g_y_y_0_0_x_xx_z_y = buffer_1100_pdpp[655];

    auto g_y_y_0_0_x_xx_z_z = buffer_1100_pdpp[656];

    auto g_y_y_0_0_x_xy_x_x = buffer_1100_pdpp[657];

    auto g_y_y_0_0_x_xy_x_y = buffer_1100_pdpp[658];

    auto g_y_y_0_0_x_xy_x_z = buffer_1100_pdpp[659];

    auto g_y_y_0_0_x_xy_y_x = buffer_1100_pdpp[660];

    auto g_y_y_0_0_x_xy_y_y = buffer_1100_pdpp[661];

    auto g_y_y_0_0_x_xy_y_z = buffer_1100_pdpp[662];

    auto g_y_y_0_0_x_xy_z_x = buffer_1100_pdpp[663];

    auto g_y_y_0_0_x_xy_z_y = buffer_1100_pdpp[664];

    auto g_y_y_0_0_x_xy_z_z = buffer_1100_pdpp[665];

    auto g_y_y_0_0_x_xz_x_x = buffer_1100_pdpp[666];

    auto g_y_y_0_0_x_xz_x_y = buffer_1100_pdpp[667];

    auto g_y_y_0_0_x_xz_x_z = buffer_1100_pdpp[668];

    auto g_y_y_0_0_x_xz_y_x = buffer_1100_pdpp[669];

    auto g_y_y_0_0_x_xz_y_y = buffer_1100_pdpp[670];

    auto g_y_y_0_0_x_xz_y_z = buffer_1100_pdpp[671];

    auto g_y_y_0_0_x_xz_z_x = buffer_1100_pdpp[672];

    auto g_y_y_0_0_x_xz_z_y = buffer_1100_pdpp[673];

    auto g_y_y_0_0_x_xz_z_z = buffer_1100_pdpp[674];

    auto g_y_y_0_0_x_yy_x_x = buffer_1100_pdpp[675];

    auto g_y_y_0_0_x_yy_x_y = buffer_1100_pdpp[676];

    auto g_y_y_0_0_x_yy_x_z = buffer_1100_pdpp[677];

    auto g_y_y_0_0_x_yy_y_x = buffer_1100_pdpp[678];

    auto g_y_y_0_0_x_yy_y_y = buffer_1100_pdpp[679];

    auto g_y_y_0_0_x_yy_y_z = buffer_1100_pdpp[680];

    auto g_y_y_0_0_x_yy_z_x = buffer_1100_pdpp[681];

    auto g_y_y_0_0_x_yy_z_y = buffer_1100_pdpp[682];

    auto g_y_y_0_0_x_yy_z_z = buffer_1100_pdpp[683];

    auto g_y_y_0_0_x_yz_x_x = buffer_1100_pdpp[684];

    auto g_y_y_0_0_x_yz_x_y = buffer_1100_pdpp[685];

    auto g_y_y_0_0_x_yz_x_z = buffer_1100_pdpp[686];

    auto g_y_y_0_0_x_yz_y_x = buffer_1100_pdpp[687];

    auto g_y_y_0_0_x_yz_y_y = buffer_1100_pdpp[688];

    auto g_y_y_0_0_x_yz_y_z = buffer_1100_pdpp[689];

    auto g_y_y_0_0_x_yz_z_x = buffer_1100_pdpp[690];

    auto g_y_y_0_0_x_yz_z_y = buffer_1100_pdpp[691];

    auto g_y_y_0_0_x_yz_z_z = buffer_1100_pdpp[692];

    auto g_y_y_0_0_x_zz_x_x = buffer_1100_pdpp[693];

    auto g_y_y_0_0_x_zz_x_y = buffer_1100_pdpp[694];

    auto g_y_y_0_0_x_zz_x_z = buffer_1100_pdpp[695];

    auto g_y_y_0_0_x_zz_y_x = buffer_1100_pdpp[696];

    auto g_y_y_0_0_x_zz_y_y = buffer_1100_pdpp[697];

    auto g_y_y_0_0_x_zz_y_z = buffer_1100_pdpp[698];

    auto g_y_y_0_0_x_zz_z_x = buffer_1100_pdpp[699];

    auto g_y_y_0_0_x_zz_z_y = buffer_1100_pdpp[700];

    auto g_y_y_0_0_x_zz_z_z = buffer_1100_pdpp[701];

    auto g_y_y_0_0_y_xx_x_x = buffer_1100_pdpp[702];

    auto g_y_y_0_0_y_xx_x_y = buffer_1100_pdpp[703];

    auto g_y_y_0_0_y_xx_x_z = buffer_1100_pdpp[704];

    auto g_y_y_0_0_y_xx_y_x = buffer_1100_pdpp[705];

    auto g_y_y_0_0_y_xx_y_y = buffer_1100_pdpp[706];

    auto g_y_y_0_0_y_xx_y_z = buffer_1100_pdpp[707];

    auto g_y_y_0_0_y_xx_z_x = buffer_1100_pdpp[708];

    auto g_y_y_0_0_y_xx_z_y = buffer_1100_pdpp[709];

    auto g_y_y_0_0_y_xx_z_z = buffer_1100_pdpp[710];

    auto g_y_y_0_0_y_xy_x_x = buffer_1100_pdpp[711];

    auto g_y_y_0_0_y_xy_x_y = buffer_1100_pdpp[712];

    auto g_y_y_0_0_y_xy_x_z = buffer_1100_pdpp[713];

    auto g_y_y_0_0_y_xy_y_x = buffer_1100_pdpp[714];

    auto g_y_y_0_0_y_xy_y_y = buffer_1100_pdpp[715];

    auto g_y_y_0_0_y_xy_y_z = buffer_1100_pdpp[716];

    auto g_y_y_0_0_y_xy_z_x = buffer_1100_pdpp[717];

    auto g_y_y_0_0_y_xy_z_y = buffer_1100_pdpp[718];

    auto g_y_y_0_0_y_xy_z_z = buffer_1100_pdpp[719];

    auto g_y_y_0_0_y_xz_x_x = buffer_1100_pdpp[720];

    auto g_y_y_0_0_y_xz_x_y = buffer_1100_pdpp[721];

    auto g_y_y_0_0_y_xz_x_z = buffer_1100_pdpp[722];

    auto g_y_y_0_0_y_xz_y_x = buffer_1100_pdpp[723];

    auto g_y_y_0_0_y_xz_y_y = buffer_1100_pdpp[724];

    auto g_y_y_0_0_y_xz_y_z = buffer_1100_pdpp[725];

    auto g_y_y_0_0_y_xz_z_x = buffer_1100_pdpp[726];

    auto g_y_y_0_0_y_xz_z_y = buffer_1100_pdpp[727];

    auto g_y_y_0_0_y_xz_z_z = buffer_1100_pdpp[728];

    auto g_y_y_0_0_y_yy_x_x = buffer_1100_pdpp[729];

    auto g_y_y_0_0_y_yy_x_y = buffer_1100_pdpp[730];

    auto g_y_y_0_0_y_yy_x_z = buffer_1100_pdpp[731];

    auto g_y_y_0_0_y_yy_y_x = buffer_1100_pdpp[732];

    auto g_y_y_0_0_y_yy_y_y = buffer_1100_pdpp[733];

    auto g_y_y_0_0_y_yy_y_z = buffer_1100_pdpp[734];

    auto g_y_y_0_0_y_yy_z_x = buffer_1100_pdpp[735];

    auto g_y_y_0_0_y_yy_z_y = buffer_1100_pdpp[736];

    auto g_y_y_0_0_y_yy_z_z = buffer_1100_pdpp[737];

    auto g_y_y_0_0_y_yz_x_x = buffer_1100_pdpp[738];

    auto g_y_y_0_0_y_yz_x_y = buffer_1100_pdpp[739];

    auto g_y_y_0_0_y_yz_x_z = buffer_1100_pdpp[740];

    auto g_y_y_0_0_y_yz_y_x = buffer_1100_pdpp[741];

    auto g_y_y_0_0_y_yz_y_y = buffer_1100_pdpp[742];

    auto g_y_y_0_0_y_yz_y_z = buffer_1100_pdpp[743];

    auto g_y_y_0_0_y_yz_z_x = buffer_1100_pdpp[744];

    auto g_y_y_0_0_y_yz_z_y = buffer_1100_pdpp[745];

    auto g_y_y_0_0_y_yz_z_z = buffer_1100_pdpp[746];

    auto g_y_y_0_0_y_zz_x_x = buffer_1100_pdpp[747];

    auto g_y_y_0_0_y_zz_x_y = buffer_1100_pdpp[748];

    auto g_y_y_0_0_y_zz_x_z = buffer_1100_pdpp[749];

    auto g_y_y_0_0_y_zz_y_x = buffer_1100_pdpp[750];

    auto g_y_y_0_0_y_zz_y_y = buffer_1100_pdpp[751];

    auto g_y_y_0_0_y_zz_y_z = buffer_1100_pdpp[752];

    auto g_y_y_0_0_y_zz_z_x = buffer_1100_pdpp[753];

    auto g_y_y_0_0_y_zz_z_y = buffer_1100_pdpp[754];

    auto g_y_y_0_0_y_zz_z_z = buffer_1100_pdpp[755];

    auto g_y_y_0_0_z_xx_x_x = buffer_1100_pdpp[756];

    auto g_y_y_0_0_z_xx_x_y = buffer_1100_pdpp[757];

    auto g_y_y_0_0_z_xx_x_z = buffer_1100_pdpp[758];

    auto g_y_y_0_0_z_xx_y_x = buffer_1100_pdpp[759];

    auto g_y_y_0_0_z_xx_y_y = buffer_1100_pdpp[760];

    auto g_y_y_0_0_z_xx_y_z = buffer_1100_pdpp[761];

    auto g_y_y_0_0_z_xx_z_x = buffer_1100_pdpp[762];

    auto g_y_y_0_0_z_xx_z_y = buffer_1100_pdpp[763];

    auto g_y_y_0_0_z_xx_z_z = buffer_1100_pdpp[764];

    auto g_y_y_0_0_z_xy_x_x = buffer_1100_pdpp[765];

    auto g_y_y_0_0_z_xy_x_y = buffer_1100_pdpp[766];

    auto g_y_y_0_0_z_xy_x_z = buffer_1100_pdpp[767];

    auto g_y_y_0_0_z_xy_y_x = buffer_1100_pdpp[768];

    auto g_y_y_0_0_z_xy_y_y = buffer_1100_pdpp[769];

    auto g_y_y_0_0_z_xy_y_z = buffer_1100_pdpp[770];

    auto g_y_y_0_0_z_xy_z_x = buffer_1100_pdpp[771];

    auto g_y_y_0_0_z_xy_z_y = buffer_1100_pdpp[772];

    auto g_y_y_0_0_z_xy_z_z = buffer_1100_pdpp[773];

    auto g_y_y_0_0_z_xz_x_x = buffer_1100_pdpp[774];

    auto g_y_y_0_0_z_xz_x_y = buffer_1100_pdpp[775];

    auto g_y_y_0_0_z_xz_x_z = buffer_1100_pdpp[776];

    auto g_y_y_0_0_z_xz_y_x = buffer_1100_pdpp[777];

    auto g_y_y_0_0_z_xz_y_y = buffer_1100_pdpp[778];

    auto g_y_y_0_0_z_xz_y_z = buffer_1100_pdpp[779];

    auto g_y_y_0_0_z_xz_z_x = buffer_1100_pdpp[780];

    auto g_y_y_0_0_z_xz_z_y = buffer_1100_pdpp[781];

    auto g_y_y_0_0_z_xz_z_z = buffer_1100_pdpp[782];

    auto g_y_y_0_0_z_yy_x_x = buffer_1100_pdpp[783];

    auto g_y_y_0_0_z_yy_x_y = buffer_1100_pdpp[784];

    auto g_y_y_0_0_z_yy_x_z = buffer_1100_pdpp[785];

    auto g_y_y_0_0_z_yy_y_x = buffer_1100_pdpp[786];

    auto g_y_y_0_0_z_yy_y_y = buffer_1100_pdpp[787];

    auto g_y_y_0_0_z_yy_y_z = buffer_1100_pdpp[788];

    auto g_y_y_0_0_z_yy_z_x = buffer_1100_pdpp[789];

    auto g_y_y_0_0_z_yy_z_y = buffer_1100_pdpp[790];

    auto g_y_y_0_0_z_yy_z_z = buffer_1100_pdpp[791];

    auto g_y_y_0_0_z_yz_x_x = buffer_1100_pdpp[792];

    auto g_y_y_0_0_z_yz_x_y = buffer_1100_pdpp[793];

    auto g_y_y_0_0_z_yz_x_z = buffer_1100_pdpp[794];

    auto g_y_y_0_0_z_yz_y_x = buffer_1100_pdpp[795];

    auto g_y_y_0_0_z_yz_y_y = buffer_1100_pdpp[796];

    auto g_y_y_0_0_z_yz_y_z = buffer_1100_pdpp[797];

    auto g_y_y_0_0_z_yz_z_x = buffer_1100_pdpp[798];

    auto g_y_y_0_0_z_yz_z_y = buffer_1100_pdpp[799];

    auto g_y_y_0_0_z_yz_z_z = buffer_1100_pdpp[800];

    auto g_y_y_0_0_z_zz_x_x = buffer_1100_pdpp[801];

    auto g_y_y_0_0_z_zz_x_y = buffer_1100_pdpp[802];

    auto g_y_y_0_0_z_zz_x_z = buffer_1100_pdpp[803];

    auto g_y_y_0_0_z_zz_y_x = buffer_1100_pdpp[804];

    auto g_y_y_0_0_z_zz_y_y = buffer_1100_pdpp[805];

    auto g_y_y_0_0_z_zz_y_z = buffer_1100_pdpp[806];

    auto g_y_y_0_0_z_zz_z_x = buffer_1100_pdpp[807];

    auto g_y_y_0_0_z_zz_z_y = buffer_1100_pdpp[808];

    auto g_y_y_0_0_z_zz_z_z = buffer_1100_pdpp[809];

    auto g_y_z_0_0_x_xx_x_x = buffer_1100_pdpp[810];

    auto g_y_z_0_0_x_xx_x_y = buffer_1100_pdpp[811];

    auto g_y_z_0_0_x_xx_x_z = buffer_1100_pdpp[812];

    auto g_y_z_0_0_x_xx_y_x = buffer_1100_pdpp[813];

    auto g_y_z_0_0_x_xx_y_y = buffer_1100_pdpp[814];

    auto g_y_z_0_0_x_xx_y_z = buffer_1100_pdpp[815];

    auto g_y_z_0_0_x_xx_z_x = buffer_1100_pdpp[816];

    auto g_y_z_0_0_x_xx_z_y = buffer_1100_pdpp[817];

    auto g_y_z_0_0_x_xx_z_z = buffer_1100_pdpp[818];

    auto g_y_z_0_0_x_xy_x_x = buffer_1100_pdpp[819];

    auto g_y_z_0_0_x_xy_x_y = buffer_1100_pdpp[820];

    auto g_y_z_0_0_x_xy_x_z = buffer_1100_pdpp[821];

    auto g_y_z_0_0_x_xy_y_x = buffer_1100_pdpp[822];

    auto g_y_z_0_0_x_xy_y_y = buffer_1100_pdpp[823];

    auto g_y_z_0_0_x_xy_y_z = buffer_1100_pdpp[824];

    auto g_y_z_0_0_x_xy_z_x = buffer_1100_pdpp[825];

    auto g_y_z_0_0_x_xy_z_y = buffer_1100_pdpp[826];

    auto g_y_z_0_0_x_xy_z_z = buffer_1100_pdpp[827];

    auto g_y_z_0_0_x_xz_x_x = buffer_1100_pdpp[828];

    auto g_y_z_0_0_x_xz_x_y = buffer_1100_pdpp[829];

    auto g_y_z_0_0_x_xz_x_z = buffer_1100_pdpp[830];

    auto g_y_z_0_0_x_xz_y_x = buffer_1100_pdpp[831];

    auto g_y_z_0_0_x_xz_y_y = buffer_1100_pdpp[832];

    auto g_y_z_0_0_x_xz_y_z = buffer_1100_pdpp[833];

    auto g_y_z_0_0_x_xz_z_x = buffer_1100_pdpp[834];

    auto g_y_z_0_0_x_xz_z_y = buffer_1100_pdpp[835];

    auto g_y_z_0_0_x_xz_z_z = buffer_1100_pdpp[836];

    auto g_y_z_0_0_x_yy_x_x = buffer_1100_pdpp[837];

    auto g_y_z_0_0_x_yy_x_y = buffer_1100_pdpp[838];

    auto g_y_z_0_0_x_yy_x_z = buffer_1100_pdpp[839];

    auto g_y_z_0_0_x_yy_y_x = buffer_1100_pdpp[840];

    auto g_y_z_0_0_x_yy_y_y = buffer_1100_pdpp[841];

    auto g_y_z_0_0_x_yy_y_z = buffer_1100_pdpp[842];

    auto g_y_z_0_0_x_yy_z_x = buffer_1100_pdpp[843];

    auto g_y_z_0_0_x_yy_z_y = buffer_1100_pdpp[844];

    auto g_y_z_0_0_x_yy_z_z = buffer_1100_pdpp[845];

    auto g_y_z_0_0_x_yz_x_x = buffer_1100_pdpp[846];

    auto g_y_z_0_0_x_yz_x_y = buffer_1100_pdpp[847];

    auto g_y_z_0_0_x_yz_x_z = buffer_1100_pdpp[848];

    auto g_y_z_0_0_x_yz_y_x = buffer_1100_pdpp[849];

    auto g_y_z_0_0_x_yz_y_y = buffer_1100_pdpp[850];

    auto g_y_z_0_0_x_yz_y_z = buffer_1100_pdpp[851];

    auto g_y_z_0_0_x_yz_z_x = buffer_1100_pdpp[852];

    auto g_y_z_0_0_x_yz_z_y = buffer_1100_pdpp[853];

    auto g_y_z_0_0_x_yz_z_z = buffer_1100_pdpp[854];

    auto g_y_z_0_0_x_zz_x_x = buffer_1100_pdpp[855];

    auto g_y_z_0_0_x_zz_x_y = buffer_1100_pdpp[856];

    auto g_y_z_0_0_x_zz_x_z = buffer_1100_pdpp[857];

    auto g_y_z_0_0_x_zz_y_x = buffer_1100_pdpp[858];

    auto g_y_z_0_0_x_zz_y_y = buffer_1100_pdpp[859];

    auto g_y_z_0_0_x_zz_y_z = buffer_1100_pdpp[860];

    auto g_y_z_0_0_x_zz_z_x = buffer_1100_pdpp[861];

    auto g_y_z_0_0_x_zz_z_y = buffer_1100_pdpp[862];

    auto g_y_z_0_0_x_zz_z_z = buffer_1100_pdpp[863];

    auto g_y_z_0_0_y_xx_x_x = buffer_1100_pdpp[864];

    auto g_y_z_0_0_y_xx_x_y = buffer_1100_pdpp[865];

    auto g_y_z_0_0_y_xx_x_z = buffer_1100_pdpp[866];

    auto g_y_z_0_0_y_xx_y_x = buffer_1100_pdpp[867];

    auto g_y_z_0_0_y_xx_y_y = buffer_1100_pdpp[868];

    auto g_y_z_0_0_y_xx_y_z = buffer_1100_pdpp[869];

    auto g_y_z_0_0_y_xx_z_x = buffer_1100_pdpp[870];

    auto g_y_z_0_0_y_xx_z_y = buffer_1100_pdpp[871];

    auto g_y_z_0_0_y_xx_z_z = buffer_1100_pdpp[872];

    auto g_y_z_0_0_y_xy_x_x = buffer_1100_pdpp[873];

    auto g_y_z_0_0_y_xy_x_y = buffer_1100_pdpp[874];

    auto g_y_z_0_0_y_xy_x_z = buffer_1100_pdpp[875];

    auto g_y_z_0_0_y_xy_y_x = buffer_1100_pdpp[876];

    auto g_y_z_0_0_y_xy_y_y = buffer_1100_pdpp[877];

    auto g_y_z_0_0_y_xy_y_z = buffer_1100_pdpp[878];

    auto g_y_z_0_0_y_xy_z_x = buffer_1100_pdpp[879];

    auto g_y_z_0_0_y_xy_z_y = buffer_1100_pdpp[880];

    auto g_y_z_0_0_y_xy_z_z = buffer_1100_pdpp[881];

    auto g_y_z_0_0_y_xz_x_x = buffer_1100_pdpp[882];

    auto g_y_z_0_0_y_xz_x_y = buffer_1100_pdpp[883];

    auto g_y_z_0_0_y_xz_x_z = buffer_1100_pdpp[884];

    auto g_y_z_0_0_y_xz_y_x = buffer_1100_pdpp[885];

    auto g_y_z_0_0_y_xz_y_y = buffer_1100_pdpp[886];

    auto g_y_z_0_0_y_xz_y_z = buffer_1100_pdpp[887];

    auto g_y_z_0_0_y_xz_z_x = buffer_1100_pdpp[888];

    auto g_y_z_0_0_y_xz_z_y = buffer_1100_pdpp[889];

    auto g_y_z_0_0_y_xz_z_z = buffer_1100_pdpp[890];

    auto g_y_z_0_0_y_yy_x_x = buffer_1100_pdpp[891];

    auto g_y_z_0_0_y_yy_x_y = buffer_1100_pdpp[892];

    auto g_y_z_0_0_y_yy_x_z = buffer_1100_pdpp[893];

    auto g_y_z_0_0_y_yy_y_x = buffer_1100_pdpp[894];

    auto g_y_z_0_0_y_yy_y_y = buffer_1100_pdpp[895];

    auto g_y_z_0_0_y_yy_y_z = buffer_1100_pdpp[896];

    auto g_y_z_0_0_y_yy_z_x = buffer_1100_pdpp[897];

    auto g_y_z_0_0_y_yy_z_y = buffer_1100_pdpp[898];

    auto g_y_z_0_0_y_yy_z_z = buffer_1100_pdpp[899];

    auto g_y_z_0_0_y_yz_x_x = buffer_1100_pdpp[900];

    auto g_y_z_0_0_y_yz_x_y = buffer_1100_pdpp[901];

    auto g_y_z_0_0_y_yz_x_z = buffer_1100_pdpp[902];

    auto g_y_z_0_0_y_yz_y_x = buffer_1100_pdpp[903];

    auto g_y_z_0_0_y_yz_y_y = buffer_1100_pdpp[904];

    auto g_y_z_0_0_y_yz_y_z = buffer_1100_pdpp[905];

    auto g_y_z_0_0_y_yz_z_x = buffer_1100_pdpp[906];

    auto g_y_z_0_0_y_yz_z_y = buffer_1100_pdpp[907];

    auto g_y_z_0_0_y_yz_z_z = buffer_1100_pdpp[908];

    auto g_y_z_0_0_y_zz_x_x = buffer_1100_pdpp[909];

    auto g_y_z_0_0_y_zz_x_y = buffer_1100_pdpp[910];

    auto g_y_z_0_0_y_zz_x_z = buffer_1100_pdpp[911];

    auto g_y_z_0_0_y_zz_y_x = buffer_1100_pdpp[912];

    auto g_y_z_0_0_y_zz_y_y = buffer_1100_pdpp[913];

    auto g_y_z_0_0_y_zz_y_z = buffer_1100_pdpp[914];

    auto g_y_z_0_0_y_zz_z_x = buffer_1100_pdpp[915];

    auto g_y_z_0_0_y_zz_z_y = buffer_1100_pdpp[916];

    auto g_y_z_0_0_y_zz_z_z = buffer_1100_pdpp[917];

    auto g_y_z_0_0_z_xx_x_x = buffer_1100_pdpp[918];

    auto g_y_z_0_0_z_xx_x_y = buffer_1100_pdpp[919];

    auto g_y_z_0_0_z_xx_x_z = buffer_1100_pdpp[920];

    auto g_y_z_0_0_z_xx_y_x = buffer_1100_pdpp[921];

    auto g_y_z_0_0_z_xx_y_y = buffer_1100_pdpp[922];

    auto g_y_z_0_0_z_xx_y_z = buffer_1100_pdpp[923];

    auto g_y_z_0_0_z_xx_z_x = buffer_1100_pdpp[924];

    auto g_y_z_0_0_z_xx_z_y = buffer_1100_pdpp[925];

    auto g_y_z_0_0_z_xx_z_z = buffer_1100_pdpp[926];

    auto g_y_z_0_0_z_xy_x_x = buffer_1100_pdpp[927];

    auto g_y_z_0_0_z_xy_x_y = buffer_1100_pdpp[928];

    auto g_y_z_0_0_z_xy_x_z = buffer_1100_pdpp[929];

    auto g_y_z_0_0_z_xy_y_x = buffer_1100_pdpp[930];

    auto g_y_z_0_0_z_xy_y_y = buffer_1100_pdpp[931];

    auto g_y_z_0_0_z_xy_y_z = buffer_1100_pdpp[932];

    auto g_y_z_0_0_z_xy_z_x = buffer_1100_pdpp[933];

    auto g_y_z_0_0_z_xy_z_y = buffer_1100_pdpp[934];

    auto g_y_z_0_0_z_xy_z_z = buffer_1100_pdpp[935];

    auto g_y_z_0_0_z_xz_x_x = buffer_1100_pdpp[936];

    auto g_y_z_0_0_z_xz_x_y = buffer_1100_pdpp[937];

    auto g_y_z_0_0_z_xz_x_z = buffer_1100_pdpp[938];

    auto g_y_z_0_0_z_xz_y_x = buffer_1100_pdpp[939];

    auto g_y_z_0_0_z_xz_y_y = buffer_1100_pdpp[940];

    auto g_y_z_0_0_z_xz_y_z = buffer_1100_pdpp[941];

    auto g_y_z_0_0_z_xz_z_x = buffer_1100_pdpp[942];

    auto g_y_z_0_0_z_xz_z_y = buffer_1100_pdpp[943];

    auto g_y_z_0_0_z_xz_z_z = buffer_1100_pdpp[944];

    auto g_y_z_0_0_z_yy_x_x = buffer_1100_pdpp[945];

    auto g_y_z_0_0_z_yy_x_y = buffer_1100_pdpp[946];

    auto g_y_z_0_0_z_yy_x_z = buffer_1100_pdpp[947];

    auto g_y_z_0_0_z_yy_y_x = buffer_1100_pdpp[948];

    auto g_y_z_0_0_z_yy_y_y = buffer_1100_pdpp[949];

    auto g_y_z_0_0_z_yy_y_z = buffer_1100_pdpp[950];

    auto g_y_z_0_0_z_yy_z_x = buffer_1100_pdpp[951];

    auto g_y_z_0_0_z_yy_z_y = buffer_1100_pdpp[952];

    auto g_y_z_0_0_z_yy_z_z = buffer_1100_pdpp[953];

    auto g_y_z_0_0_z_yz_x_x = buffer_1100_pdpp[954];

    auto g_y_z_0_0_z_yz_x_y = buffer_1100_pdpp[955];

    auto g_y_z_0_0_z_yz_x_z = buffer_1100_pdpp[956];

    auto g_y_z_0_0_z_yz_y_x = buffer_1100_pdpp[957];

    auto g_y_z_0_0_z_yz_y_y = buffer_1100_pdpp[958];

    auto g_y_z_0_0_z_yz_y_z = buffer_1100_pdpp[959];

    auto g_y_z_0_0_z_yz_z_x = buffer_1100_pdpp[960];

    auto g_y_z_0_0_z_yz_z_y = buffer_1100_pdpp[961];

    auto g_y_z_0_0_z_yz_z_z = buffer_1100_pdpp[962];

    auto g_y_z_0_0_z_zz_x_x = buffer_1100_pdpp[963];

    auto g_y_z_0_0_z_zz_x_y = buffer_1100_pdpp[964];

    auto g_y_z_0_0_z_zz_x_z = buffer_1100_pdpp[965];

    auto g_y_z_0_0_z_zz_y_x = buffer_1100_pdpp[966];

    auto g_y_z_0_0_z_zz_y_y = buffer_1100_pdpp[967];

    auto g_y_z_0_0_z_zz_y_z = buffer_1100_pdpp[968];

    auto g_y_z_0_0_z_zz_z_x = buffer_1100_pdpp[969];

    auto g_y_z_0_0_z_zz_z_y = buffer_1100_pdpp[970];

    auto g_y_z_0_0_z_zz_z_z = buffer_1100_pdpp[971];

    auto g_z_x_0_0_x_xx_x_x = buffer_1100_pdpp[972];

    auto g_z_x_0_0_x_xx_x_y = buffer_1100_pdpp[973];

    auto g_z_x_0_0_x_xx_x_z = buffer_1100_pdpp[974];

    auto g_z_x_0_0_x_xx_y_x = buffer_1100_pdpp[975];

    auto g_z_x_0_0_x_xx_y_y = buffer_1100_pdpp[976];

    auto g_z_x_0_0_x_xx_y_z = buffer_1100_pdpp[977];

    auto g_z_x_0_0_x_xx_z_x = buffer_1100_pdpp[978];

    auto g_z_x_0_0_x_xx_z_y = buffer_1100_pdpp[979];

    auto g_z_x_0_0_x_xx_z_z = buffer_1100_pdpp[980];

    auto g_z_x_0_0_x_xy_x_x = buffer_1100_pdpp[981];

    auto g_z_x_0_0_x_xy_x_y = buffer_1100_pdpp[982];

    auto g_z_x_0_0_x_xy_x_z = buffer_1100_pdpp[983];

    auto g_z_x_0_0_x_xy_y_x = buffer_1100_pdpp[984];

    auto g_z_x_0_0_x_xy_y_y = buffer_1100_pdpp[985];

    auto g_z_x_0_0_x_xy_y_z = buffer_1100_pdpp[986];

    auto g_z_x_0_0_x_xy_z_x = buffer_1100_pdpp[987];

    auto g_z_x_0_0_x_xy_z_y = buffer_1100_pdpp[988];

    auto g_z_x_0_0_x_xy_z_z = buffer_1100_pdpp[989];

    auto g_z_x_0_0_x_xz_x_x = buffer_1100_pdpp[990];

    auto g_z_x_0_0_x_xz_x_y = buffer_1100_pdpp[991];

    auto g_z_x_0_0_x_xz_x_z = buffer_1100_pdpp[992];

    auto g_z_x_0_0_x_xz_y_x = buffer_1100_pdpp[993];

    auto g_z_x_0_0_x_xz_y_y = buffer_1100_pdpp[994];

    auto g_z_x_0_0_x_xz_y_z = buffer_1100_pdpp[995];

    auto g_z_x_0_0_x_xz_z_x = buffer_1100_pdpp[996];

    auto g_z_x_0_0_x_xz_z_y = buffer_1100_pdpp[997];

    auto g_z_x_0_0_x_xz_z_z = buffer_1100_pdpp[998];

    auto g_z_x_0_0_x_yy_x_x = buffer_1100_pdpp[999];

    auto g_z_x_0_0_x_yy_x_y = buffer_1100_pdpp[1000];

    auto g_z_x_0_0_x_yy_x_z = buffer_1100_pdpp[1001];

    auto g_z_x_0_0_x_yy_y_x = buffer_1100_pdpp[1002];

    auto g_z_x_0_0_x_yy_y_y = buffer_1100_pdpp[1003];

    auto g_z_x_0_0_x_yy_y_z = buffer_1100_pdpp[1004];

    auto g_z_x_0_0_x_yy_z_x = buffer_1100_pdpp[1005];

    auto g_z_x_0_0_x_yy_z_y = buffer_1100_pdpp[1006];

    auto g_z_x_0_0_x_yy_z_z = buffer_1100_pdpp[1007];

    auto g_z_x_0_0_x_yz_x_x = buffer_1100_pdpp[1008];

    auto g_z_x_0_0_x_yz_x_y = buffer_1100_pdpp[1009];

    auto g_z_x_0_0_x_yz_x_z = buffer_1100_pdpp[1010];

    auto g_z_x_0_0_x_yz_y_x = buffer_1100_pdpp[1011];

    auto g_z_x_0_0_x_yz_y_y = buffer_1100_pdpp[1012];

    auto g_z_x_0_0_x_yz_y_z = buffer_1100_pdpp[1013];

    auto g_z_x_0_0_x_yz_z_x = buffer_1100_pdpp[1014];

    auto g_z_x_0_0_x_yz_z_y = buffer_1100_pdpp[1015];

    auto g_z_x_0_0_x_yz_z_z = buffer_1100_pdpp[1016];

    auto g_z_x_0_0_x_zz_x_x = buffer_1100_pdpp[1017];

    auto g_z_x_0_0_x_zz_x_y = buffer_1100_pdpp[1018];

    auto g_z_x_0_0_x_zz_x_z = buffer_1100_pdpp[1019];

    auto g_z_x_0_0_x_zz_y_x = buffer_1100_pdpp[1020];

    auto g_z_x_0_0_x_zz_y_y = buffer_1100_pdpp[1021];

    auto g_z_x_0_0_x_zz_y_z = buffer_1100_pdpp[1022];

    auto g_z_x_0_0_x_zz_z_x = buffer_1100_pdpp[1023];

    auto g_z_x_0_0_x_zz_z_y = buffer_1100_pdpp[1024];

    auto g_z_x_0_0_x_zz_z_z = buffer_1100_pdpp[1025];

    auto g_z_x_0_0_y_xx_x_x = buffer_1100_pdpp[1026];

    auto g_z_x_0_0_y_xx_x_y = buffer_1100_pdpp[1027];

    auto g_z_x_0_0_y_xx_x_z = buffer_1100_pdpp[1028];

    auto g_z_x_0_0_y_xx_y_x = buffer_1100_pdpp[1029];

    auto g_z_x_0_0_y_xx_y_y = buffer_1100_pdpp[1030];

    auto g_z_x_0_0_y_xx_y_z = buffer_1100_pdpp[1031];

    auto g_z_x_0_0_y_xx_z_x = buffer_1100_pdpp[1032];

    auto g_z_x_0_0_y_xx_z_y = buffer_1100_pdpp[1033];

    auto g_z_x_0_0_y_xx_z_z = buffer_1100_pdpp[1034];

    auto g_z_x_0_0_y_xy_x_x = buffer_1100_pdpp[1035];

    auto g_z_x_0_0_y_xy_x_y = buffer_1100_pdpp[1036];

    auto g_z_x_0_0_y_xy_x_z = buffer_1100_pdpp[1037];

    auto g_z_x_0_0_y_xy_y_x = buffer_1100_pdpp[1038];

    auto g_z_x_0_0_y_xy_y_y = buffer_1100_pdpp[1039];

    auto g_z_x_0_0_y_xy_y_z = buffer_1100_pdpp[1040];

    auto g_z_x_0_0_y_xy_z_x = buffer_1100_pdpp[1041];

    auto g_z_x_0_0_y_xy_z_y = buffer_1100_pdpp[1042];

    auto g_z_x_0_0_y_xy_z_z = buffer_1100_pdpp[1043];

    auto g_z_x_0_0_y_xz_x_x = buffer_1100_pdpp[1044];

    auto g_z_x_0_0_y_xz_x_y = buffer_1100_pdpp[1045];

    auto g_z_x_0_0_y_xz_x_z = buffer_1100_pdpp[1046];

    auto g_z_x_0_0_y_xz_y_x = buffer_1100_pdpp[1047];

    auto g_z_x_0_0_y_xz_y_y = buffer_1100_pdpp[1048];

    auto g_z_x_0_0_y_xz_y_z = buffer_1100_pdpp[1049];

    auto g_z_x_0_0_y_xz_z_x = buffer_1100_pdpp[1050];

    auto g_z_x_0_0_y_xz_z_y = buffer_1100_pdpp[1051];

    auto g_z_x_0_0_y_xz_z_z = buffer_1100_pdpp[1052];

    auto g_z_x_0_0_y_yy_x_x = buffer_1100_pdpp[1053];

    auto g_z_x_0_0_y_yy_x_y = buffer_1100_pdpp[1054];

    auto g_z_x_0_0_y_yy_x_z = buffer_1100_pdpp[1055];

    auto g_z_x_0_0_y_yy_y_x = buffer_1100_pdpp[1056];

    auto g_z_x_0_0_y_yy_y_y = buffer_1100_pdpp[1057];

    auto g_z_x_0_0_y_yy_y_z = buffer_1100_pdpp[1058];

    auto g_z_x_0_0_y_yy_z_x = buffer_1100_pdpp[1059];

    auto g_z_x_0_0_y_yy_z_y = buffer_1100_pdpp[1060];

    auto g_z_x_0_0_y_yy_z_z = buffer_1100_pdpp[1061];

    auto g_z_x_0_0_y_yz_x_x = buffer_1100_pdpp[1062];

    auto g_z_x_0_0_y_yz_x_y = buffer_1100_pdpp[1063];

    auto g_z_x_0_0_y_yz_x_z = buffer_1100_pdpp[1064];

    auto g_z_x_0_0_y_yz_y_x = buffer_1100_pdpp[1065];

    auto g_z_x_0_0_y_yz_y_y = buffer_1100_pdpp[1066];

    auto g_z_x_0_0_y_yz_y_z = buffer_1100_pdpp[1067];

    auto g_z_x_0_0_y_yz_z_x = buffer_1100_pdpp[1068];

    auto g_z_x_0_0_y_yz_z_y = buffer_1100_pdpp[1069];

    auto g_z_x_0_0_y_yz_z_z = buffer_1100_pdpp[1070];

    auto g_z_x_0_0_y_zz_x_x = buffer_1100_pdpp[1071];

    auto g_z_x_0_0_y_zz_x_y = buffer_1100_pdpp[1072];

    auto g_z_x_0_0_y_zz_x_z = buffer_1100_pdpp[1073];

    auto g_z_x_0_0_y_zz_y_x = buffer_1100_pdpp[1074];

    auto g_z_x_0_0_y_zz_y_y = buffer_1100_pdpp[1075];

    auto g_z_x_0_0_y_zz_y_z = buffer_1100_pdpp[1076];

    auto g_z_x_0_0_y_zz_z_x = buffer_1100_pdpp[1077];

    auto g_z_x_0_0_y_zz_z_y = buffer_1100_pdpp[1078];

    auto g_z_x_0_0_y_zz_z_z = buffer_1100_pdpp[1079];

    auto g_z_x_0_0_z_xx_x_x = buffer_1100_pdpp[1080];

    auto g_z_x_0_0_z_xx_x_y = buffer_1100_pdpp[1081];

    auto g_z_x_0_0_z_xx_x_z = buffer_1100_pdpp[1082];

    auto g_z_x_0_0_z_xx_y_x = buffer_1100_pdpp[1083];

    auto g_z_x_0_0_z_xx_y_y = buffer_1100_pdpp[1084];

    auto g_z_x_0_0_z_xx_y_z = buffer_1100_pdpp[1085];

    auto g_z_x_0_0_z_xx_z_x = buffer_1100_pdpp[1086];

    auto g_z_x_0_0_z_xx_z_y = buffer_1100_pdpp[1087];

    auto g_z_x_0_0_z_xx_z_z = buffer_1100_pdpp[1088];

    auto g_z_x_0_0_z_xy_x_x = buffer_1100_pdpp[1089];

    auto g_z_x_0_0_z_xy_x_y = buffer_1100_pdpp[1090];

    auto g_z_x_0_0_z_xy_x_z = buffer_1100_pdpp[1091];

    auto g_z_x_0_0_z_xy_y_x = buffer_1100_pdpp[1092];

    auto g_z_x_0_0_z_xy_y_y = buffer_1100_pdpp[1093];

    auto g_z_x_0_0_z_xy_y_z = buffer_1100_pdpp[1094];

    auto g_z_x_0_0_z_xy_z_x = buffer_1100_pdpp[1095];

    auto g_z_x_0_0_z_xy_z_y = buffer_1100_pdpp[1096];

    auto g_z_x_0_0_z_xy_z_z = buffer_1100_pdpp[1097];

    auto g_z_x_0_0_z_xz_x_x = buffer_1100_pdpp[1098];

    auto g_z_x_0_0_z_xz_x_y = buffer_1100_pdpp[1099];

    auto g_z_x_0_0_z_xz_x_z = buffer_1100_pdpp[1100];

    auto g_z_x_0_0_z_xz_y_x = buffer_1100_pdpp[1101];

    auto g_z_x_0_0_z_xz_y_y = buffer_1100_pdpp[1102];

    auto g_z_x_0_0_z_xz_y_z = buffer_1100_pdpp[1103];

    auto g_z_x_0_0_z_xz_z_x = buffer_1100_pdpp[1104];

    auto g_z_x_0_0_z_xz_z_y = buffer_1100_pdpp[1105];

    auto g_z_x_0_0_z_xz_z_z = buffer_1100_pdpp[1106];

    auto g_z_x_0_0_z_yy_x_x = buffer_1100_pdpp[1107];

    auto g_z_x_0_0_z_yy_x_y = buffer_1100_pdpp[1108];

    auto g_z_x_0_0_z_yy_x_z = buffer_1100_pdpp[1109];

    auto g_z_x_0_0_z_yy_y_x = buffer_1100_pdpp[1110];

    auto g_z_x_0_0_z_yy_y_y = buffer_1100_pdpp[1111];

    auto g_z_x_0_0_z_yy_y_z = buffer_1100_pdpp[1112];

    auto g_z_x_0_0_z_yy_z_x = buffer_1100_pdpp[1113];

    auto g_z_x_0_0_z_yy_z_y = buffer_1100_pdpp[1114];

    auto g_z_x_0_0_z_yy_z_z = buffer_1100_pdpp[1115];

    auto g_z_x_0_0_z_yz_x_x = buffer_1100_pdpp[1116];

    auto g_z_x_0_0_z_yz_x_y = buffer_1100_pdpp[1117];

    auto g_z_x_0_0_z_yz_x_z = buffer_1100_pdpp[1118];

    auto g_z_x_0_0_z_yz_y_x = buffer_1100_pdpp[1119];

    auto g_z_x_0_0_z_yz_y_y = buffer_1100_pdpp[1120];

    auto g_z_x_0_0_z_yz_y_z = buffer_1100_pdpp[1121];

    auto g_z_x_0_0_z_yz_z_x = buffer_1100_pdpp[1122];

    auto g_z_x_0_0_z_yz_z_y = buffer_1100_pdpp[1123];

    auto g_z_x_0_0_z_yz_z_z = buffer_1100_pdpp[1124];

    auto g_z_x_0_0_z_zz_x_x = buffer_1100_pdpp[1125];

    auto g_z_x_0_0_z_zz_x_y = buffer_1100_pdpp[1126];

    auto g_z_x_0_0_z_zz_x_z = buffer_1100_pdpp[1127];

    auto g_z_x_0_0_z_zz_y_x = buffer_1100_pdpp[1128];

    auto g_z_x_0_0_z_zz_y_y = buffer_1100_pdpp[1129];

    auto g_z_x_0_0_z_zz_y_z = buffer_1100_pdpp[1130];

    auto g_z_x_0_0_z_zz_z_x = buffer_1100_pdpp[1131];

    auto g_z_x_0_0_z_zz_z_y = buffer_1100_pdpp[1132];

    auto g_z_x_0_0_z_zz_z_z = buffer_1100_pdpp[1133];

    auto g_z_y_0_0_x_xx_x_x = buffer_1100_pdpp[1134];

    auto g_z_y_0_0_x_xx_x_y = buffer_1100_pdpp[1135];

    auto g_z_y_0_0_x_xx_x_z = buffer_1100_pdpp[1136];

    auto g_z_y_0_0_x_xx_y_x = buffer_1100_pdpp[1137];

    auto g_z_y_0_0_x_xx_y_y = buffer_1100_pdpp[1138];

    auto g_z_y_0_0_x_xx_y_z = buffer_1100_pdpp[1139];

    auto g_z_y_0_0_x_xx_z_x = buffer_1100_pdpp[1140];

    auto g_z_y_0_0_x_xx_z_y = buffer_1100_pdpp[1141];

    auto g_z_y_0_0_x_xx_z_z = buffer_1100_pdpp[1142];

    auto g_z_y_0_0_x_xy_x_x = buffer_1100_pdpp[1143];

    auto g_z_y_0_0_x_xy_x_y = buffer_1100_pdpp[1144];

    auto g_z_y_0_0_x_xy_x_z = buffer_1100_pdpp[1145];

    auto g_z_y_0_0_x_xy_y_x = buffer_1100_pdpp[1146];

    auto g_z_y_0_0_x_xy_y_y = buffer_1100_pdpp[1147];

    auto g_z_y_0_0_x_xy_y_z = buffer_1100_pdpp[1148];

    auto g_z_y_0_0_x_xy_z_x = buffer_1100_pdpp[1149];

    auto g_z_y_0_0_x_xy_z_y = buffer_1100_pdpp[1150];

    auto g_z_y_0_0_x_xy_z_z = buffer_1100_pdpp[1151];

    auto g_z_y_0_0_x_xz_x_x = buffer_1100_pdpp[1152];

    auto g_z_y_0_0_x_xz_x_y = buffer_1100_pdpp[1153];

    auto g_z_y_0_0_x_xz_x_z = buffer_1100_pdpp[1154];

    auto g_z_y_0_0_x_xz_y_x = buffer_1100_pdpp[1155];

    auto g_z_y_0_0_x_xz_y_y = buffer_1100_pdpp[1156];

    auto g_z_y_0_0_x_xz_y_z = buffer_1100_pdpp[1157];

    auto g_z_y_0_0_x_xz_z_x = buffer_1100_pdpp[1158];

    auto g_z_y_0_0_x_xz_z_y = buffer_1100_pdpp[1159];

    auto g_z_y_0_0_x_xz_z_z = buffer_1100_pdpp[1160];

    auto g_z_y_0_0_x_yy_x_x = buffer_1100_pdpp[1161];

    auto g_z_y_0_0_x_yy_x_y = buffer_1100_pdpp[1162];

    auto g_z_y_0_0_x_yy_x_z = buffer_1100_pdpp[1163];

    auto g_z_y_0_0_x_yy_y_x = buffer_1100_pdpp[1164];

    auto g_z_y_0_0_x_yy_y_y = buffer_1100_pdpp[1165];

    auto g_z_y_0_0_x_yy_y_z = buffer_1100_pdpp[1166];

    auto g_z_y_0_0_x_yy_z_x = buffer_1100_pdpp[1167];

    auto g_z_y_0_0_x_yy_z_y = buffer_1100_pdpp[1168];

    auto g_z_y_0_0_x_yy_z_z = buffer_1100_pdpp[1169];

    auto g_z_y_0_0_x_yz_x_x = buffer_1100_pdpp[1170];

    auto g_z_y_0_0_x_yz_x_y = buffer_1100_pdpp[1171];

    auto g_z_y_0_0_x_yz_x_z = buffer_1100_pdpp[1172];

    auto g_z_y_0_0_x_yz_y_x = buffer_1100_pdpp[1173];

    auto g_z_y_0_0_x_yz_y_y = buffer_1100_pdpp[1174];

    auto g_z_y_0_0_x_yz_y_z = buffer_1100_pdpp[1175];

    auto g_z_y_0_0_x_yz_z_x = buffer_1100_pdpp[1176];

    auto g_z_y_0_0_x_yz_z_y = buffer_1100_pdpp[1177];

    auto g_z_y_0_0_x_yz_z_z = buffer_1100_pdpp[1178];

    auto g_z_y_0_0_x_zz_x_x = buffer_1100_pdpp[1179];

    auto g_z_y_0_0_x_zz_x_y = buffer_1100_pdpp[1180];

    auto g_z_y_0_0_x_zz_x_z = buffer_1100_pdpp[1181];

    auto g_z_y_0_0_x_zz_y_x = buffer_1100_pdpp[1182];

    auto g_z_y_0_0_x_zz_y_y = buffer_1100_pdpp[1183];

    auto g_z_y_0_0_x_zz_y_z = buffer_1100_pdpp[1184];

    auto g_z_y_0_0_x_zz_z_x = buffer_1100_pdpp[1185];

    auto g_z_y_0_0_x_zz_z_y = buffer_1100_pdpp[1186];

    auto g_z_y_0_0_x_zz_z_z = buffer_1100_pdpp[1187];

    auto g_z_y_0_0_y_xx_x_x = buffer_1100_pdpp[1188];

    auto g_z_y_0_0_y_xx_x_y = buffer_1100_pdpp[1189];

    auto g_z_y_0_0_y_xx_x_z = buffer_1100_pdpp[1190];

    auto g_z_y_0_0_y_xx_y_x = buffer_1100_pdpp[1191];

    auto g_z_y_0_0_y_xx_y_y = buffer_1100_pdpp[1192];

    auto g_z_y_0_0_y_xx_y_z = buffer_1100_pdpp[1193];

    auto g_z_y_0_0_y_xx_z_x = buffer_1100_pdpp[1194];

    auto g_z_y_0_0_y_xx_z_y = buffer_1100_pdpp[1195];

    auto g_z_y_0_0_y_xx_z_z = buffer_1100_pdpp[1196];

    auto g_z_y_0_0_y_xy_x_x = buffer_1100_pdpp[1197];

    auto g_z_y_0_0_y_xy_x_y = buffer_1100_pdpp[1198];

    auto g_z_y_0_0_y_xy_x_z = buffer_1100_pdpp[1199];

    auto g_z_y_0_0_y_xy_y_x = buffer_1100_pdpp[1200];

    auto g_z_y_0_0_y_xy_y_y = buffer_1100_pdpp[1201];

    auto g_z_y_0_0_y_xy_y_z = buffer_1100_pdpp[1202];

    auto g_z_y_0_0_y_xy_z_x = buffer_1100_pdpp[1203];

    auto g_z_y_0_0_y_xy_z_y = buffer_1100_pdpp[1204];

    auto g_z_y_0_0_y_xy_z_z = buffer_1100_pdpp[1205];

    auto g_z_y_0_0_y_xz_x_x = buffer_1100_pdpp[1206];

    auto g_z_y_0_0_y_xz_x_y = buffer_1100_pdpp[1207];

    auto g_z_y_0_0_y_xz_x_z = buffer_1100_pdpp[1208];

    auto g_z_y_0_0_y_xz_y_x = buffer_1100_pdpp[1209];

    auto g_z_y_0_0_y_xz_y_y = buffer_1100_pdpp[1210];

    auto g_z_y_0_0_y_xz_y_z = buffer_1100_pdpp[1211];

    auto g_z_y_0_0_y_xz_z_x = buffer_1100_pdpp[1212];

    auto g_z_y_0_0_y_xz_z_y = buffer_1100_pdpp[1213];

    auto g_z_y_0_0_y_xz_z_z = buffer_1100_pdpp[1214];

    auto g_z_y_0_0_y_yy_x_x = buffer_1100_pdpp[1215];

    auto g_z_y_0_0_y_yy_x_y = buffer_1100_pdpp[1216];

    auto g_z_y_0_0_y_yy_x_z = buffer_1100_pdpp[1217];

    auto g_z_y_0_0_y_yy_y_x = buffer_1100_pdpp[1218];

    auto g_z_y_0_0_y_yy_y_y = buffer_1100_pdpp[1219];

    auto g_z_y_0_0_y_yy_y_z = buffer_1100_pdpp[1220];

    auto g_z_y_0_0_y_yy_z_x = buffer_1100_pdpp[1221];

    auto g_z_y_0_0_y_yy_z_y = buffer_1100_pdpp[1222];

    auto g_z_y_0_0_y_yy_z_z = buffer_1100_pdpp[1223];

    auto g_z_y_0_0_y_yz_x_x = buffer_1100_pdpp[1224];

    auto g_z_y_0_0_y_yz_x_y = buffer_1100_pdpp[1225];

    auto g_z_y_0_0_y_yz_x_z = buffer_1100_pdpp[1226];

    auto g_z_y_0_0_y_yz_y_x = buffer_1100_pdpp[1227];

    auto g_z_y_0_0_y_yz_y_y = buffer_1100_pdpp[1228];

    auto g_z_y_0_0_y_yz_y_z = buffer_1100_pdpp[1229];

    auto g_z_y_0_0_y_yz_z_x = buffer_1100_pdpp[1230];

    auto g_z_y_0_0_y_yz_z_y = buffer_1100_pdpp[1231];

    auto g_z_y_0_0_y_yz_z_z = buffer_1100_pdpp[1232];

    auto g_z_y_0_0_y_zz_x_x = buffer_1100_pdpp[1233];

    auto g_z_y_0_0_y_zz_x_y = buffer_1100_pdpp[1234];

    auto g_z_y_0_0_y_zz_x_z = buffer_1100_pdpp[1235];

    auto g_z_y_0_0_y_zz_y_x = buffer_1100_pdpp[1236];

    auto g_z_y_0_0_y_zz_y_y = buffer_1100_pdpp[1237];

    auto g_z_y_0_0_y_zz_y_z = buffer_1100_pdpp[1238];

    auto g_z_y_0_0_y_zz_z_x = buffer_1100_pdpp[1239];

    auto g_z_y_0_0_y_zz_z_y = buffer_1100_pdpp[1240];

    auto g_z_y_0_0_y_zz_z_z = buffer_1100_pdpp[1241];

    auto g_z_y_0_0_z_xx_x_x = buffer_1100_pdpp[1242];

    auto g_z_y_0_0_z_xx_x_y = buffer_1100_pdpp[1243];

    auto g_z_y_0_0_z_xx_x_z = buffer_1100_pdpp[1244];

    auto g_z_y_0_0_z_xx_y_x = buffer_1100_pdpp[1245];

    auto g_z_y_0_0_z_xx_y_y = buffer_1100_pdpp[1246];

    auto g_z_y_0_0_z_xx_y_z = buffer_1100_pdpp[1247];

    auto g_z_y_0_0_z_xx_z_x = buffer_1100_pdpp[1248];

    auto g_z_y_0_0_z_xx_z_y = buffer_1100_pdpp[1249];

    auto g_z_y_0_0_z_xx_z_z = buffer_1100_pdpp[1250];

    auto g_z_y_0_0_z_xy_x_x = buffer_1100_pdpp[1251];

    auto g_z_y_0_0_z_xy_x_y = buffer_1100_pdpp[1252];

    auto g_z_y_0_0_z_xy_x_z = buffer_1100_pdpp[1253];

    auto g_z_y_0_0_z_xy_y_x = buffer_1100_pdpp[1254];

    auto g_z_y_0_0_z_xy_y_y = buffer_1100_pdpp[1255];

    auto g_z_y_0_0_z_xy_y_z = buffer_1100_pdpp[1256];

    auto g_z_y_0_0_z_xy_z_x = buffer_1100_pdpp[1257];

    auto g_z_y_0_0_z_xy_z_y = buffer_1100_pdpp[1258];

    auto g_z_y_0_0_z_xy_z_z = buffer_1100_pdpp[1259];

    auto g_z_y_0_0_z_xz_x_x = buffer_1100_pdpp[1260];

    auto g_z_y_0_0_z_xz_x_y = buffer_1100_pdpp[1261];

    auto g_z_y_0_0_z_xz_x_z = buffer_1100_pdpp[1262];

    auto g_z_y_0_0_z_xz_y_x = buffer_1100_pdpp[1263];

    auto g_z_y_0_0_z_xz_y_y = buffer_1100_pdpp[1264];

    auto g_z_y_0_0_z_xz_y_z = buffer_1100_pdpp[1265];

    auto g_z_y_0_0_z_xz_z_x = buffer_1100_pdpp[1266];

    auto g_z_y_0_0_z_xz_z_y = buffer_1100_pdpp[1267];

    auto g_z_y_0_0_z_xz_z_z = buffer_1100_pdpp[1268];

    auto g_z_y_0_0_z_yy_x_x = buffer_1100_pdpp[1269];

    auto g_z_y_0_0_z_yy_x_y = buffer_1100_pdpp[1270];

    auto g_z_y_0_0_z_yy_x_z = buffer_1100_pdpp[1271];

    auto g_z_y_0_0_z_yy_y_x = buffer_1100_pdpp[1272];

    auto g_z_y_0_0_z_yy_y_y = buffer_1100_pdpp[1273];

    auto g_z_y_0_0_z_yy_y_z = buffer_1100_pdpp[1274];

    auto g_z_y_0_0_z_yy_z_x = buffer_1100_pdpp[1275];

    auto g_z_y_0_0_z_yy_z_y = buffer_1100_pdpp[1276];

    auto g_z_y_0_0_z_yy_z_z = buffer_1100_pdpp[1277];

    auto g_z_y_0_0_z_yz_x_x = buffer_1100_pdpp[1278];

    auto g_z_y_0_0_z_yz_x_y = buffer_1100_pdpp[1279];

    auto g_z_y_0_0_z_yz_x_z = buffer_1100_pdpp[1280];

    auto g_z_y_0_0_z_yz_y_x = buffer_1100_pdpp[1281];

    auto g_z_y_0_0_z_yz_y_y = buffer_1100_pdpp[1282];

    auto g_z_y_0_0_z_yz_y_z = buffer_1100_pdpp[1283];

    auto g_z_y_0_0_z_yz_z_x = buffer_1100_pdpp[1284];

    auto g_z_y_0_0_z_yz_z_y = buffer_1100_pdpp[1285];

    auto g_z_y_0_0_z_yz_z_z = buffer_1100_pdpp[1286];

    auto g_z_y_0_0_z_zz_x_x = buffer_1100_pdpp[1287];

    auto g_z_y_0_0_z_zz_x_y = buffer_1100_pdpp[1288];

    auto g_z_y_0_0_z_zz_x_z = buffer_1100_pdpp[1289];

    auto g_z_y_0_0_z_zz_y_x = buffer_1100_pdpp[1290];

    auto g_z_y_0_0_z_zz_y_y = buffer_1100_pdpp[1291];

    auto g_z_y_0_0_z_zz_y_z = buffer_1100_pdpp[1292];

    auto g_z_y_0_0_z_zz_z_x = buffer_1100_pdpp[1293];

    auto g_z_y_0_0_z_zz_z_y = buffer_1100_pdpp[1294];

    auto g_z_y_0_0_z_zz_z_z = buffer_1100_pdpp[1295];

    auto g_z_z_0_0_x_xx_x_x = buffer_1100_pdpp[1296];

    auto g_z_z_0_0_x_xx_x_y = buffer_1100_pdpp[1297];

    auto g_z_z_0_0_x_xx_x_z = buffer_1100_pdpp[1298];

    auto g_z_z_0_0_x_xx_y_x = buffer_1100_pdpp[1299];

    auto g_z_z_0_0_x_xx_y_y = buffer_1100_pdpp[1300];

    auto g_z_z_0_0_x_xx_y_z = buffer_1100_pdpp[1301];

    auto g_z_z_0_0_x_xx_z_x = buffer_1100_pdpp[1302];

    auto g_z_z_0_0_x_xx_z_y = buffer_1100_pdpp[1303];

    auto g_z_z_0_0_x_xx_z_z = buffer_1100_pdpp[1304];

    auto g_z_z_0_0_x_xy_x_x = buffer_1100_pdpp[1305];

    auto g_z_z_0_0_x_xy_x_y = buffer_1100_pdpp[1306];

    auto g_z_z_0_0_x_xy_x_z = buffer_1100_pdpp[1307];

    auto g_z_z_0_0_x_xy_y_x = buffer_1100_pdpp[1308];

    auto g_z_z_0_0_x_xy_y_y = buffer_1100_pdpp[1309];

    auto g_z_z_0_0_x_xy_y_z = buffer_1100_pdpp[1310];

    auto g_z_z_0_0_x_xy_z_x = buffer_1100_pdpp[1311];

    auto g_z_z_0_0_x_xy_z_y = buffer_1100_pdpp[1312];

    auto g_z_z_0_0_x_xy_z_z = buffer_1100_pdpp[1313];

    auto g_z_z_0_0_x_xz_x_x = buffer_1100_pdpp[1314];

    auto g_z_z_0_0_x_xz_x_y = buffer_1100_pdpp[1315];

    auto g_z_z_0_0_x_xz_x_z = buffer_1100_pdpp[1316];

    auto g_z_z_0_0_x_xz_y_x = buffer_1100_pdpp[1317];

    auto g_z_z_0_0_x_xz_y_y = buffer_1100_pdpp[1318];

    auto g_z_z_0_0_x_xz_y_z = buffer_1100_pdpp[1319];

    auto g_z_z_0_0_x_xz_z_x = buffer_1100_pdpp[1320];

    auto g_z_z_0_0_x_xz_z_y = buffer_1100_pdpp[1321];

    auto g_z_z_0_0_x_xz_z_z = buffer_1100_pdpp[1322];

    auto g_z_z_0_0_x_yy_x_x = buffer_1100_pdpp[1323];

    auto g_z_z_0_0_x_yy_x_y = buffer_1100_pdpp[1324];

    auto g_z_z_0_0_x_yy_x_z = buffer_1100_pdpp[1325];

    auto g_z_z_0_0_x_yy_y_x = buffer_1100_pdpp[1326];

    auto g_z_z_0_0_x_yy_y_y = buffer_1100_pdpp[1327];

    auto g_z_z_0_0_x_yy_y_z = buffer_1100_pdpp[1328];

    auto g_z_z_0_0_x_yy_z_x = buffer_1100_pdpp[1329];

    auto g_z_z_0_0_x_yy_z_y = buffer_1100_pdpp[1330];

    auto g_z_z_0_0_x_yy_z_z = buffer_1100_pdpp[1331];

    auto g_z_z_0_0_x_yz_x_x = buffer_1100_pdpp[1332];

    auto g_z_z_0_0_x_yz_x_y = buffer_1100_pdpp[1333];

    auto g_z_z_0_0_x_yz_x_z = buffer_1100_pdpp[1334];

    auto g_z_z_0_0_x_yz_y_x = buffer_1100_pdpp[1335];

    auto g_z_z_0_0_x_yz_y_y = buffer_1100_pdpp[1336];

    auto g_z_z_0_0_x_yz_y_z = buffer_1100_pdpp[1337];

    auto g_z_z_0_0_x_yz_z_x = buffer_1100_pdpp[1338];

    auto g_z_z_0_0_x_yz_z_y = buffer_1100_pdpp[1339];

    auto g_z_z_0_0_x_yz_z_z = buffer_1100_pdpp[1340];

    auto g_z_z_0_0_x_zz_x_x = buffer_1100_pdpp[1341];

    auto g_z_z_0_0_x_zz_x_y = buffer_1100_pdpp[1342];

    auto g_z_z_0_0_x_zz_x_z = buffer_1100_pdpp[1343];

    auto g_z_z_0_0_x_zz_y_x = buffer_1100_pdpp[1344];

    auto g_z_z_0_0_x_zz_y_y = buffer_1100_pdpp[1345];

    auto g_z_z_0_0_x_zz_y_z = buffer_1100_pdpp[1346];

    auto g_z_z_0_0_x_zz_z_x = buffer_1100_pdpp[1347];

    auto g_z_z_0_0_x_zz_z_y = buffer_1100_pdpp[1348];

    auto g_z_z_0_0_x_zz_z_z = buffer_1100_pdpp[1349];

    auto g_z_z_0_0_y_xx_x_x = buffer_1100_pdpp[1350];

    auto g_z_z_0_0_y_xx_x_y = buffer_1100_pdpp[1351];

    auto g_z_z_0_0_y_xx_x_z = buffer_1100_pdpp[1352];

    auto g_z_z_0_0_y_xx_y_x = buffer_1100_pdpp[1353];

    auto g_z_z_0_0_y_xx_y_y = buffer_1100_pdpp[1354];

    auto g_z_z_0_0_y_xx_y_z = buffer_1100_pdpp[1355];

    auto g_z_z_0_0_y_xx_z_x = buffer_1100_pdpp[1356];

    auto g_z_z_0_0_y_xx_z_y = buffer_1100_pdpp[1357];

    auto g_z_z_0_0_y_xx_z_z = buffer_1100_pdpp[1358];

    auto g_z_z_0_0_y_xy_x_x = buffer_1100_pdpp[1359];

    auto g_z_z_0_0_y_xy_x_y = buffer_1100_pdpp[1360];

    auto g_z_z_0_0_y_xy_x_z = buffer_1100_pdpp[1361];

    auto g_z_z_0_0_y_xy_y_x = buffer_1100_pdpp[1362];

    auto g_z_z_0_0_y_xy_y_y = buffer_1100_pdpp[1363];

    auto g_z_z_0_0_y_xy_y_z = buffer_1100_pdpp[1364];

    auto g_z_z_0_0_y_xy_z_x = buffer_1100_pdpp[1365];

    auto g_z_z_0_0_y_xy_z_y = buffer_1100_pdpp[1366];

    auto g_z_z_0_0_y_xy_z_z = buffer_1100_pdpp[1367];

    auto g_z_z_0_0_y_xz_x_x = buffer_1100_pdpp[1368];

    auto g_z_z_0_0_y_xz_x_y = buffer_1100_pdpp[1369];

    auto g_z_z_0_0_y_xz_x_z = buffer_1100_pdpp[1370];

    auto g_z_z_0_0_y_xz_y_x = buffer_1100_pdpp[1371];

    auto g_z_z_0_0_y_xz_y_y = buffer_1100_pdpp[1372];

    auto g_z_z_0_0_y_xz_y_z = buffer_1100_pdpp[1373];

    auto g_z_z_0_0_y_xz_z_x = buffer_1100_pdpp[1374];

    auto g_z_z_0_0_y_xz_z_y = buffer_1100_pdpp[1375];

    auto g_z_z_0_0_y_xz_z_z = buffer_1100_pdpp[1376];

    auto g_z_z_0_0_y_yy_x_x = buffer_1100_pdpp[1377];

    auto g_z_z_0_0_y_yy_x_y = buffer_1100_pdpp[1378];

    auto g_z_z_0_0_y_yy_x_z = buffer_1100_pdpp[1379];

    auto g_z_z_0_0_y_yy_y_x = buffer_1100_pdpp[1380];

    auto g_z_z_0_0_y_yy_y_y = buffer_1100_pdpp[1381];

    auto g_z_z_0_0_y_yy_y_z = buffer_1100_pdpp[1382];

    auto g_z_z_0_0_y_yy_z_x = buffer_1100_pdpp[1383];

    auto g_z_z_0_0_y_yy_z_y = buffer_1100_pdpp[1384];

    auto g_z_z_0_0_y_yy_z_z = buffer_1100_pdpp[1385];

    auto g_z_z_0_0_y_yz_x_x = buffer_1100_pdpp[1386];

    auto g_z_z_0_0_y_yz_x_y = buffer_1100_pdpp[1387];

    auto g_z_z_0_0_y_yz_x_z = buffer_1100_pdpp[1388];

    auto g_z_z_0_0_y_yz_y_x = buffer_1100_pdpp[1389];

    auto g_z_z_0_0_y_yz_y_y = buffer_1100_pdpp[1390];

    auto g_z_z_0_0_y_yz_y_z = buffer_1100_pdpp[1391];

    auto g_z_z_0_0_y_yz_z_x = buffer_1100_pdpp[1392];

    auto g_z_z_0_0_y_yz_z_y = buffer_1100_pdpp[1393];

    auto g_z_z_0_0_y_yz_z_z = buffer_1100_pdpp[1394];

    auto g_z_z_0_0_y_zz_x_x = buffer_1100_pdpp[1395];

    auto g_z_z_0_0_y_zz_x_y = buffer_1100_pdpp[1396];

    auto g_z_z_0_0_y_zz_x_z = buffer_1100_pdpp[1397];

    auto g_z_z_0_0_y_zz_y_x = buffer_1100_pdpp[1398];

    auto g_z_z_0_0_y_zz_y_y = buffer_1100_pdpp[1399];

    auto g_z_z_0_0_y_zz_y_z = buffer_1100_pdpp[1400];

    auto g_z_z_0_0_y_zz_z_x = buffer_1100_pdpp[1401];

    auto g_z_z_0_0_y_zz_z_y = buffer_1100_pdpp[1402];

    auto g_z_z_0_0_y_zz_z_z = buffer_1100_pdpp[1403];

    auto g_z_z_0_0_z_xx_x_x = buffer_1100_pdpp[1404];

    auto g_z_z_0_0_z_xx_x_y = buffer_1100_pdpp[1405];

    auto g_z_z_0_0_z_xx_x_z = buffer_1100_pdpp[1406];

    auto g_z_z_0_0_z_xx_y_x = buffer_1100_pdpp[1407];

    auto g_z_z_0_0_z_xx_y_y = buffer_1100_pdpp[1408];

    auto g_z_z_0_0_z_xx_y_z = buffer_1100_pdpp[1409];

    auto g_z_z_0_0_z_xx_z_x = buffer_1100_pdpp[1410];

    auto g_z_z_0_0_z_xx_z_y = buffer_1100_pdpp[1411];

    auto g_z_z_0_0_z_xx_z_z = buffer_1100_pdpp[1412];

    auto g_z_z_0_0_z_xy_x_x = buffer_1100_pdpp[1413];

    auto g_z_z_0_0_z_xy_x_y = buffer_1100_pdpp[1414];

    auto g_z_z_0_0_z_xy_x_z = buffer_1100_pdpp[1415];

    auto g_z_z_0_0_z_xy_y_x = buffer_1100_pdpp[1416];

    auto g_z_z_0_0_z_xy_y_y = buffer_1100_pdpp[1417];

    auto g_z_z_0_0_z_xy_y_z = buffer_1100_pdpp[1418];

    auto g_z_z_0_0_z_xy_z_x = buffer_1100_pdpp[1419];

    auto g_z_z_0_0_z_xy_z_y = buffer_1100_pdpp[1420];

    auto g_z_z_0_0_z_xy_z_z = buffer_1100_pdpp[1421];

    auto g_z_z_0_0_z_xz_x_x = buffer_1100_pdpp[1422];

    auto g_z_z_0_0_z_xz_x_y = buffer_1100_pdpp[1423];

    auto g_z_z_0_0_z_xz_x_z = buffer_1100_pdpp[1424];

    auto g_z_z_0_0_z_xz_y_x = buffer_1100_pdpp[1425];

    auto g_z_z_0_0_z_xz_y_y = buffer_1100_pdpp[1426];

    auto g_z_z_0_0_z_xz_y_z = buffer_1100_pdpp[1427];

    auto g_z_z_0_0_z_xz_z_x = buffer_1100_pdpp[1428];

    auto g_z_z_0_0_z_xz_z_y = buffer_1100_pdpp[1429];

    auto g_z_z_0_0_z_xz_z_z = buffer_1100_pdpp[1430];

    auto g_z_z_0_0_z_yy_x_x = buffer_1100_pdpp[1431];

    auto g_z_z_0_0_z_yy_x_y = buffer_1100_pdpp[1432];

    auto g_z_z_0_0_z_yy_x_z = buffer_1100_pdpp[1433];

    auto g_z_z_0_0_z_yy_y_x = buffer_1100_pdpp[1434];

    auto g_z_z_0_0_z_yy_y_y = buffer_1100_pdpp[1435];

    auto g_z_z_0_0_z_yy_y_z = buffer_1100_pdpp[1436];

    auto g_z_z_0_0_z_yy_z_x = buffer_1100_pdpp[1437];

    auto g_z_z_0_0_z_yy_z_y = buffer_1100_pdpp[1438];

    auto g_z_z_0_0_z_yy_z_z = buffer_1100_pdpp[1439];

    auto g_z_z_0_0_z_yz_x_x = buffer_1100_pdpp[1440];

    auto g_z_z_0_0_z_yz_x_y = buffer_1100_pdpp[1441];

    auto g_z_z_0_0_z_yz_x_z = buffer_1100_pdpp[1442];

    auto g_z_z_0_0_z_yz_y_x = buffer_1100_pdpp[1443];

    auto g_z_z_0_0_z_yz_y_y = buffer_1100_pdpp[1444];

    auto g_z_z_0_0_z_yz_y_z = buffer_1100_pdpp[1445];

    auto g_z_z_0_0_z_yz_z_x = buffer_1100_pdpp[1446];

    auto g_z_z_0_0_z_yz_z_y = buffer_1100_pdpp[1447];

    auto g_z_z_0_0_z_yz_z_z = buffer_1100_pdpp[1448];

    auto g_z_z_0_0_z_zz_x_x = buffer_1100_pdpp[1449];

    auto g_z_z_0_0_z_zz_x_y = buffer_1100_pdpp[1450];

    auto g_z_z_0_0_z_zz_x_z = buffer_1100_pdpp[1451];

    auto g_z_z_0_0_z_zz_y_x = buffer_1100_pdpp[1452];

    auto g_z_z_0_0_z_zz_y_y = buffer_1100_pdpp[1453];

    auto g_z_z_0_0_z_zz_y_z = buffer_1100_pdpp[1454];

    auto g_z_z_0_0_z_zz_z_x = buffer_1100_pdpp[1455];

    auto g_z_z_0_0_z_zz_z_y = buffer_1100_pdpp[1456];

    auto g_z_z_0_0_z_zz_z_z = buffer_1100_pdpp[1457];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xxx_x_x, g_0_xxx_x_y, g_0_xxx_x_z, g_x_x_0_0_x_xx_x_x, g_x_x_0_0_x_xx_x_y, g_x_x_0_0_x_xx_x_z, g_xx_x_x_x, g_xx_x_x_y, g_xx_x_x_z, g_xx_xxx_x_x, g_xx_xxx_x_y, g_xx_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_x_x[i] = 2.0 * g_0_x_x_x[i] - 2.0 * g_0_xxx_x_x[i] * b_exp - 4.0 * g_xx_x_x_x[i] * a_exp + 4.0 * g_xx_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_x_y[i] = 2.0 * g_0_x_x_y[i] - 2.0 * g_0_xxx_x_y[i] * b_exp - 4.0 * g_xx_x_x_y[i] * a_exp + 4.0 * g_xx_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_x_z[i] = 2.0 * g_0_x_x_z[i] - 2.0 * g_0_xxx_x_z[i] * b_exp - 4.0 * g_xx_x_x_z[i] * a_exp + 4.0 * g_xx_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xxx_y_x, g_0_xxx_y_y, g_0_xxx_y_z, g_x_x_0_0_x_xx_y_x, g_x_x_0_0_x_xx_y_y, g_x_x_0_0_x_xx_y_z, g_xx_x_y_x, g_xx_x_y_y, g_xx_x_y_z, g_xx_xxx_y_x, g_xx_xxx_y_y, g_xx_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_y_x[i] = 2.0 * g_0_x_y_x[i] - 2.0 * g_0_xxx_y_x[i] * b_exp - 4.0 * g_xx_x_y_x[i] * a_exp + 4.0 * g_xx_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_y_y[i] = 2.0 * g_0_x_y_y[i] - 2.0 * g_0_xxx_y_y[i] * b_exp - 4.0 * g_xx_x_y_y[i] * a_exp + 4.0 * g_xx_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_y_z[i] = 2.0 * g_0_x_y_z[i] - 2.0 * g_0_xxx_y_z[i] * b_exp - 4.0 * g_xx_x_y_z[i] * a_exp + 4.0 * g_xx_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xxx_z_x, g_0_xxx_z_y, g_0_xxx_z_z, g_x_x_0_0_x_xx_z_x, g_x_x_0_0_x_xx_z_y, g_x_x_0_0_x_xx_z_z, g_xx_x_z_x, g_xx_x_z_y, g_xx_x_z_z, g_xx_xxx_z_x, g_xx_xxx_z_y, g_xx_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_z_x[i] = 2.0 * g_0_x_z_x[i] - 2.0 * g_0_xxx_z_x[i] * b_exp - 4.0 * g_xx_x_z_x[i] * a_exp + 4.0 * g_xx_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_z_y[i] = 2.0 * g_0_x_z_y[i] - 2.0 * g_0_xxx_z_y[i] * b_exp - 4.0 * g_xx_x_z_y[i] * a_exp + 4.0 * g_xx_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_z_z[i] = 2.0 * g_0_x_z_z[i] - 2.0 * g_0_xxx_z_z[i] * b_exp - 4.0 * g_xx_x_z_z[i] * a_exp + 4.0 * g_xx_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_x_x_0_0_x_xy_x_x, g_x_x_0_0_x_xy_x_y, g_x_x_0_0_x_xy_x_z, g_xx_xxy_x_x, g_xx_xxy_x_y, g_xx_xxy_x_z, g_xx_y_x_x, g_xx_y_x_y, g_xx_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xy_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_xxy_x_x[i] * b_exp - 2.0 * g_xx_y_x_x[i] * a_exp + 4.0 * g_xx_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_xxy_x_y[i] * b_exp - 2.0 * g_xx_y_x_y[i] * a_exp + 4.0 * g_xx_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_xxy_x_z[i] * b_exp - 2.0 * g_xx_y_x_z[i] * a_exp + 4.0 * g_xx_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_x_x_0_0_x_xy_y_x, g_x_x_0_0_x_xy_y_y, g_x_x_0_0_x_xy_y_z, g_xx_xxy_y_x, g_xx_xxy_y_y, g_xx_xxy_y_z, g_xx_y_y_x, g_xx_y_y_y, g_xx_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xy_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_xxy_y_x[i] * b_exp - 2.0 * g_xx_y_y_x[i] * a_exp + 4.0 * g_xx_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_xxy_y_y[i] * b_exp - 2.0 * g_xx_y_y_y[i] * a_exp + 4.0 * g_xx_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_xxy_y_z[i] * b_exp - 2.0 * g_xx_y_y_z[i] * a_exp + 4.0 * g_xx_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_x_x_0_0_x_xy_z_x, g_x_x_0_0_x_xy_z_y, g_x_x_0_0_x_xy_z_z, g_xx_xxy_z_x, g_xx_xxy_z_y, g_xx_xxy_z_z, g_xx_y_z_x, g_xx_y_z_y, g_xx_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xy_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_xxy_z_x[i] * b_exp - 2.0 * g_xx_y_z_x[i] * a_exp + 4.0 * g_xx_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_xxy_z_y[i] * b_exp - 2.0 * g_xx_y_z_y[i] * a_exp + 4.0 * g_xx_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_xxy_z_z[i] * b_exp - 2.0 * g_xx_y_z_z[i] * a_exp + 4.0 * g_xx_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_x_x_0_0_x_xz_x_x, g_x_x_0_0_x_xz_x_y, g_x_x_0_0_x_xz_x_z, g_xx_xxz_x_x, g_xx_xxz_x_y, g_xx_xxz_x_z, g_xx_z_x_x, g_xx_z_x_y, g_xx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_xxz_x_x[i] * b_exp - 2.0 * g_xx_z_x_x[i] * a_exp + 4.0 * g_xx_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_xxz_x_y[i] * b_exp - 2.0 * g_xx_z_x_y[i] * a_exp + 4.0 * g_xx_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_xxz_x_z[i] * b_exp - 2.0 * g_xx_z_x_z[i] * a_exp + 4.0 * g_xx_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_x_x_0_0_x_xz_y_x, g_x_x_0_0_x_xz_y_y, g_x_x_0_0_x_xz_y_z, g_xx_xxz_y_x, g_xx_xxz_y_y, g_xx_xxz_y_z, g_xx_z_y_x, g_xx_z_y_y, g_xx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_xxz_y_x[i] * b_exp - 2.0 * g_xx_z_y_x[i] * a_exp + 4.0 * g_xx_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_xxz_y_y[i] * b_exp - 2.0 * g_xx_z_y_y[i] * a_exp + 4.0 * g_xx_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_xxz_y_z[i] * b_exp - 2.0 * g_xx_z_y_z[i] * a_exp + 4.0 * g_xx_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_x_x_0_0_x_xz_z_x, g_x_x_0_0_x_xz_z_y, g_x_x_0_0_x_xz_z_z, g_xx_xxz_z_x, g_xx_xxz_z_y, g_xx_xxz_z_z, g_xx_z_z_x, g_xx_z_z_y, g_xx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_xxz_z_x[i] * b_exp - 2.0 * g_xx_z_z_x[i] * a_exp + 4.0 * g_xx_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_xxz_z_y[i] * b_exp - 2.0 * g_xx_z_z_y[i] * a_exp + 4.0 * g_xx_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_xxz_z_z[i] * b_exp - 2.0 * g_xx_z_z_z[i] * a_exp + 4.0 * g_xx_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_x_x_0_0_x_yy_x_x, g_x_x_0_0_x_yy_x_y, g_x_x_0_0_x_yy_x_z, g_xx_xyy_x_x, g_xx_xyy_x_y, g_xx_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yy_x_x[i] = -2.0 * g_0_xyy_x_x[i] * b_exp + 4.0 * g_xx_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_x_y[i] = -2.0 * g_0_xyy_x_y[i] * b_exp + 4.0 * g_xx_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_x_z[i] = -2.0 * g_0_xyy_x_z[i] * b_exp + 4.0 * g_xx_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_x_x_0_0_x_yy_y_x, g_x_x_0_0_x_yy_y_y, g_x_x_0_0_x_yy_y_z, g_xx_xyy_y_x, g_xx_xyy_y_y, g_xx_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yy_y_x[i] = -2.0 * g_0_xyy_y_x[i] * b_exp + 4.0 * g_xx_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_y_y[i] = -2.0 * g_0_xyy_y_y[i] * b_exp + 4.0 * g_xx_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_y_z[i] = -2.0 * g_0_xyy_y_z[i] * b_exp + 4.0 * g_xx_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_x_x_0_0_x_yy_z_x, g_x_x_0_0_x_yy_z_y, g_x_x_0_0_x_yy_z_z, g_xx_xyy_z_x, g_xx_xyy_z_y, g_xx_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yy_z_x[i] = -2.0 * g_0_xyy_z_x[i] * b_exp + 4.0 * g_xx_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_z_y[i] = -2.0 * g_0_xyy_z_y[i] * b_exp + 4.0 * g_xx_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_z_z[i] = -2.0 * g_0_xyy_z_z[i] * b_exp + 4.0 * g_xx_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_x_x_0_0_x_yz_x_x, g_x_x_0_0_x_yz_x_y, g_x_x_0_0_x_yz_x_z, g_xx_xyz_x_x, g_xx_xyz_x_y, g_xx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_xx_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_xx_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_xx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_x_x_0_0_x_yz_y_x, g_x_x_0_0_x_yz_y_y, g_x_x_0_0_x_yz_y_z, g_xx_xyz_y_x, g_xx_xyz_y_y, g_xx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_xx_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_xx_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_xx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_x_x_0_0_x_yz_z_x, g_x_x_0_0_x_yz_z_y, g_x_x_0_0_x_yz_z_z, g_xx_xyz_z_x, g_xx_xyz_z_y, g_xx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_xx_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_xx_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_xx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_x_x_0_0_x_zz_x_x, g_x_x_0_0_x_zz_x_y, g_x_x_0_0_x_zz_x_z, g_xx_xzz_x_x, g_xx_xzz_x_y, g_xx_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_zz_x_x[i] = -2.0 * g_0_xzz_x_x[i] * b_exp + 4.0 * g_xx_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_x_y[i] = -2.0 * g_0_xzz_x_y[i] * b_exp + 4.0 * g_xx_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_x_z[i] = -2.0 * g_0_xzz_x_z[i] * b_exp + 4.0 * g_xx_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_x_x_0_0_x_zz_y_x, g_x_x_0_0_x_zz_y_y, g_x_x_0_0_x_zz_y_z, g_xx_xzz_y_x, g_xx_xzz_y_y, g_xx_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_zz_y_x[i] = -2.0 * g_0_xzz_y_x[i] * b_exp + 4.0 * g_xx_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_y_y[i] = -2.0 * g_0_xzz_y_y[i] * b_exp + 4.0 * g_xx_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_y_z[i] = -2.0 * g_0_xzz_y_z[i] * b_exp + 4.0 * g_xx_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_x_x_0_0_x_zz_z_x, g_x_x_0_0_x_zz_z_y, g_x_x_0_0_x_zz_z_z, g_xx_xzz_z_x, g_xx_xzz_z_y, g_xx_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_zz_z_x[i] = -2.0 * g_0_xzz_z_x[i] * b_exp + 4.0 * g_xx_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_z_y[i] = -2.0 * g_0_xzz_z_y[i] * b_exp + 4.0 * g_xx_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_z_z[i] = -2.0 * g_0_xzz_z_z[i] * b_exp + 4.0 * g_xx_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_x_x, g_x_x_0_0_y_xx_x_y, g_x_x_0_0_y_xx_x_z, g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xxx_x_x, g_xy_xxx_x_y, g_xy_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_x_x[i] = -4.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_x_y[i] = -4.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_x_z[i] = -4.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_y_x, g_x_x_0_0_y_xx_y_y, g_x_x_0_0_y_xx_y_z, g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xxx_y_x, g_xy_xxx_y_y, g_xy_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_y_x[i] = -4.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_y_y[i] = -4.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_y_z[i] = -4.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_z_x, g_x_x_0_0_y_xx_z_y, g_x_x_0_0_y_xx_z_z, g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xxx_z_x, g_xy_xxx_z_y, g_xy_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_z_x[i] = -4.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_z_y[i] = -4.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_z_z[i] = -4.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_x_0_0_y_xy_x_x, g_x_x_0_0_y_xy_x_y, g_x_x_0_0_y_xy_x_z, g_xy_xxy_x_x, g_xy_xxy_x_y, g_xy_xxy_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xy_x_x[i] = -2.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_x_y[i] = -2.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_x_z[i] = -2.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_x_0_0_y_xy_y_x, g_x_x_0_0_y_xy_y_y, g_x_x_0_0_y_xy_y_z, g_xy_xxy_y_x, g_xy_xxy_y_y, g_xy_xxy_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xy_y_x[i] = -2.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_y_y[i] = -2.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_y_z[i] = -2.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_x_0_0_y_xy_z_x, g_x_x_0_0_y_xy_z_y, g_x_x_0_0_y_xy_z_z, g_xy_xxy_z_x, g_xy_xxy_z_y, g_xy_xxy_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xy_z_x[i] = -2.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_z_y[i] = -2.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_z_z[i] = -2.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_x_0_0_y_xz_x_x, g_x_x_0_0_y_xz_x_y, g_x_x_0_0_y_xz_x_z, g_xy_xxz_x_x, g_xy_xxz_x_y, g_xy_xxz_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xz_x_x[i] = -2.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_x_y[i] = -2.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_x_z[i] = -2.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_x_0_0_y_xz_y_x, g_x_x_0_0_y_xz_y_y, g_x_x_0_0_y_xz_y_z, g_xy_xxz_y_x, g_xy_xxz_y_y, g_xy_xxz_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xz_y_x[i] = -2.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_y_y[i] = -2.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_y_z[i] = -2.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_x_0_0_y_xz_z_x, g_x_x_0_0_y_xz_z_y, g_x_x_0_0_y_xz_z_z, g_xy_xxz_z_x, g_xy_xxz_z_y, g_xy_xxz_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xz_z_x[i] = -2.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_z_y[i] = -2.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_z_z[i] = -2.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_x_0_0_y_yy_x_x, g_x_x_0_0_y_yy_x_y, g_x_x_0_0_y_yy_x_z, g_xy_xyy_x_x, g_xy_xyy_x_y, g_xy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yy_x_x[i] = 4.0 * g_xy_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_x_y[i] = 4.0 * g_xy_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_x_z[i] = 4.0 * g_xy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_x_0_0_y_yy_y_x, g_x_x_0_0_y_yy_y_y, g_x_x_0_0_y_yy_y_z, g_xy_xyy_y_x, g_xy_xyy_y_y, g_xy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yy_y_x[i] = 4.0 * g_xy_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_y_y[i] = 4.0 * g_xy_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_y_z[i] = 4.0 * g_xy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_x_0_0_y_yy_z_x, g_x_x_0_0_y_yy_z_y, g_x_x_0_0_y_yy_z_z, g_xy_xyy_z_x, g_xy_xyy_z_y, g_xy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yy_z_x[i] = 4.0 * g_xy_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_z_y[i] = 4.0 * g_xy_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_z_z[i] = 4.0 * g_xy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_x_0_0_y_yz_x_x, g_x_x_0_0_y_yz_x_y, g_x_x_0_0_y_yz_x_z, g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yz_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_x_0_0_y_yz_y_x, g_x_x_0_0_y_yz_y_y, g_x_x_0_0_y_yz_y_z, g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yz_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_x_0_0_y_yz_z_x, g_x_x_0_0_y_yz_z_y, g_x_x_0_0_y_yz_z_z, g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yz_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_x_0_0_y_zz_x_x, g_x_x_0_0_y_zz_x_y, g_x_x_0_0_y_zz_x_z, g_xy_xzz_x_x, g_xy_xzz_x_y, g_xy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_zz_x_x[i] = 4.0 * g_xy_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_x_y[i] = 4.0 * g_xy_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_x_z[i] = 4.0 * g_xy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_x_0_0_y_zz_y_x, g_x_x_0_0_y_zz_y_y, g_x_x_0_0_y_zz_y_z, g_xy_xzz_y_x, g_xy_xzz_y_y, g_xy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_zz_y_x[i] = 4.0 * g_xy_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_y_y[i] = 4.0 * g_xy_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_y_z[i] = 4.0 * g_xy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_x_0_0_y_zz_z_x, g_x_x_0_0_y_zz_z_y, g_x_x_0_0_y_zz_z_z, g_xy_xzz_z_x, g_xy_xzz_z_y, g_xy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_zz_z_x[i] = 4.0 * g_xy_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_z_y[i] = 4.0 * g_xy_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_z_z[i] = 4.0 * g_xy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_x_x, g_x_x_0_0_z_xx_x_y, g_x_x_0_0_z_xx_x_z, g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xxx_x_x, g_xz_xxx_x_y, g_xz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_x_x[i] = -4.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xxx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_x_y[i] = -4.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xxx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_x_z[i] = -4.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_y_x, g_x_x_0_0_z_xx_y_y, g_x_x_0_0_z_xx_y_z, g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xxx_y_x, g_xz_xxx_y_y, g_xz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_y_x[i] = -4.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xxx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_y_y[i] = -4.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xxx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_y_z[i] = -4.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_z_x, g_x_x_0_0_z_xx_z_y, g_x_x_0_0_z_xx_z_z, g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xxx_z_x, g_xz_xxx_z_y, g_xz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_z_x[i] = -4.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xxx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_z_y[i] = -4.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xxx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_z_z[i] = -4.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_x_0_0_z_xy_x_x, g_x_x_0_0_z_xy_x_y, g_x_x_0_0_z_xy_x_z, g_xz_xxy_x_x, g_xz_xxy_x_y, g_xz_xxy_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xy_x_x[i] = -2.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_xxy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_x_y[i] = -2.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_xxy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_x_z[i] = -2.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_x_0_0_z_xy_y_x, g_x_x_0_0_z_xy_y_y, g_x_x_0_0_z_xy_y_z, g_xz_xxy_y_x, g_xz_xxy_y_y, g_xz_xxy_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xy_y_x[i] = -2.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_xxy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_y_y[i] = -2.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_xxy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_y_z[i] = -2.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_x_0_0_z_xy_z_x, g_x_x_0_0_z_xy_z_y, g_x_x_0_0_z_xy_z_z, g_xz_xxy_z_x, g_xz_xxy_z_y, g_xz_xxy_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xy_z_x[i] = -2.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_xxy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_z_y[i] = -2.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_xxy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_z_z[i] = -2.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_x_0_0_z_xz_x_x, g_x_x_0_0_z_xz_x_y, g_x_x_0_0_z_xz_x_z, g_xz_xxz_x_x, g_xz_xxz_x_y, g_xz_xxz_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xz_x_x[i] = -2.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_xxz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_x_y[i] = -2.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_xxz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_x_z[i] = -2.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_x_0_0_z_xz_y_x, g_x_x_0_0_z_xz_y_y, g_x_x_0_0_z_xz_y_z, g_xz_xxz_y_x, g_xz_xxz_y_y, g_xz_xxz_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xz_y_x[i] = -2.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_xxz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_y_y[i] = -2.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_xxz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_y_z[i] = -2.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_x_0_0_z_xz_z_x, g_x_x_0_0_z_xz_z_y, g_x_x_0_0_z_xz_z_z, g_xz_xxz_z_x, g_xz_xxz_z_y, g_xz_xxz_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xz_z_x[i] = -2.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_xxz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_z_y[i] = -2.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_xxz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_z_z[i] = -2.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_x_0_0_z_yy_x_x, g_x_x_0_0_z_yy_x_y, g_x_x_0_0_z_yy_x_z, g_xz_xyy_x_x, g_xz_xyy_x_y, g_xz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yy_x_x[i] = 4.0 * g_xz_xyy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_x_y[i] = 4.0 * g_xz_xyy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_x_z[i] = 4.0 * g_xz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_x_0_0_z_yy_y_x, g_x_x_0_0_z_yy_y_y, g_x_x_0_0_z_yy_y_z, g_xz_xyy_y_x, g_xz_xyy_y_y, g_xz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yy_y_x[i] = 4.0 * g_xz_xyy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_y_y[i] = 4.0 * g_xz_xyy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_y_z[i] = 4.0 * g_xz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_x_0_0_z_yy_z_x, g_x_x_0_0_z_yy_z_y, g_x_x_0_0_z_yy_z_z, g_xz_xyy_z_x, g_xz_xyy_z_y, g_xz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yy_z_x[i] = 4.0 * g_xz_xyy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_z_y[i] = 4.0 * g_xz_xyy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_z_z[i] = 4.0 * g_xz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_x_0_0_z_yz_x_x, g_x_x_0_0_z_yz_x_y, g_x_x_0_0_z_yz_x_z, g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yz_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_x_0_0_z_yz_y_x, g_x_x_0_0_z_yz_y_y, g_x_x_0_0_z_yz_y_z, g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yz_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_x_0_0_z_yz_z_x, g_x_x_0_0_z_yz_z_y, g_x_x_0_0_z_yz_z_z, g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yz_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_x_0_0_z_zz_x_x, g_x_x_0_0_z_zz_x_y, g_x_x_0_0_z_zz_x_z, g_xz_xzz_x_x, g_xz_xzz_x_y, g_xz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_zz_x_x[i] = 4.0 * g_xz_xzz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_x_y[i] = 4.0 * g_xz_xzz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_x_z[i] = 4.0 * g_xz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_x_0_0_z_zz_y_x, g_x_x_0_0_z_zz_y_y, g_x_x_0_0_z_zz_y_z, g_xz_xzz_y_x, g_xz_xzz_y_y, g_xz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_zz_y_x[i] = 4.0 * g_xz_xzz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_y_y[i] = 4.0 * g_xz_xzz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_y_z[i] = 4.0 * g_xz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_x_0_0_z_zz_z_x, g_x_x_0_0_z_zz_z_y, g_x_x_0_0_z_zz_z_z, g_xz_xzz_z_x, g_xz_xzz_z_y, g_xz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_zz_z_x[i] = 4.0 * g_xz_xzz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_z_y[i] = 4.0 * g_xz_xzz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_z_z[i] = 4.0 * g_xz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_x_y_0_0_x_xx_x_x, g_x_y_0_0_x_xx_x_y, g_x_y_0_0_x_xx_x_z, g_xx_xxy_x_x, g_xx_xxy_x_y, g_xx_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_x_x[i] = -2.0 * g_0_xxy_x_x[i] * b_exp + 4.0 * g_xx_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_x_y[i] = -2.0 * g_0_xxy_x_y[i] * b_exp + 4.0 * g_xx_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_x_z[i] = -2.0 * g_0_xxy_x_z[i] * b_exp + 4.0 * g_xx_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_x_y_0_0_x_xx_y_x, g_x_y_0_0_x_xx_y_y, g_x_y_0_0_x_xx_y_z, g_xx_xxy_y_x, g_xx_xxy_y_y, g_xx_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_y_x[i] = -2.0 * g_0_xxy_y_x[i] * b_exp + 4.0 * g_xx_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_y_y[i] = -2.0 * g_0_xxy_y_y[i] * b_exp + 4.0 * g_xx_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_y_z[i] = -2.0 * g_0_xxy_y_z[i] * b_exp + 4.0 * g_xx_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_x_y_0_0_x_xx_z_x, g_x_y_0_0_x_xx_z_y, g_x_y_0_0_x_xx_z_z, g_xx_xxy_z_x, g_xx_xxy_z_y, g_xx_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_z_x[i] = -2.0 * g_0_xxy_z_x[i] * b_exp + 4.0 * g_xx_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_z_y[i] = -2.0 * g_0_xxy_z_y[i] * b_exp + 4.0 * g_xx_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_z_z[i] = -2.0 * g_0_xxy_z_z[i] * b_exp + 4.0 * g_xx_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_x_y_0_0_x_xy_x_x, g_x_y_0_0_x_xy_x_y, g_x_y_0_0_x_xy_x_z, g_xx_x_x_x, g_xx_x_x_y, g_xx_x_x_z, g_xx_xyy_x_x, g_xx_xyy_x_y, g_xx_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xy_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xyy_x_x[i] * b_exp - 2.0 * g_xx_x_x_x[i] * a_exp + 4.0 * g_xx_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xyy_x_y[i] * b_exp - 2.0 * g_xx_x_x_y[i] * a_exp + 4.0 * g_xx_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xyy_x_z[i] * b_exp - 2.0 * g_xx_x_x_z[i] * a_exp + 4.0 * g_xx_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_x_y_0_0_x_xy_y_x, g_x_y_0_0_x_xy_y_y, g_x_y_0_0_x_xy_y_z, g_xx_x_y_x, g_xx_x_y_y, g_xx_x_y_z, g_xx_xyy_y_x, g_xx_xyy_y_y, g_xx_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xy_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xyy_y_x[i] * b_exp - 2.0 * g_xx_x_y_x[i] * a_exp + 4.0 * g_xx_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xyy_y_y[i] * b_exp - 2.0 * g_xx_x_y_y[i] * a_exp + 4.0 * g_xx_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xyy_y_z[i] * b_exp - 2.0 * g_xx_x_y_z[i] * a_exp + 4.0 * g_xx_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_x_y_0_0_x_xy_z_x, g_x_y_0_0_x_xy_z_y, g_x_y_0_0_x_xy_z_z, g_xx_x_z_x, g_xx_x_z_y, g_xx_x_z_z, g_xx_xyy_z_x, g_xx_xyy_z_y, g_xx_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xy_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xyy_z_x[i] * b_exp - 2.0 * g_xx_x_z_x[i] * a_exp + 4.0 * g_xx_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xyy_z_y[i] * b_exp - 2.0 * g_xx_x_z_y[i] * a_exp + 4.0 * g_xx_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xyy_z_z[i] * b_exp - 2.0 * g_xx_x_z_z[i] * a_exp + 4.0 * g_xx_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_x_y_0_0_x_xz_x_x, g_x_y_0_0_x_xz_x_y, g_x_y_0_0_x_xz_x_z, g_xx_xyz_x_x, g_xx_xyz_x_y, g_xx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_xx_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_xx_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_xx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_x_y_0_0_x_xz_y_x, g_x_y_0_0_x_xz_y_y, g_x_y_0_0_x_xz_y_z, g_xx_xyz_y_x, g_xx_xyz_y_y, g_xx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_xx_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_xx_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_xx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_x_y_0_0_x_xz_z_x, g_x_y_0_0_x_xz_z_y, g_x_y_0_0_x_xz_z_z, g_xx_xyz_z_x, g_xx_xyz_z_y, g_xx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_xx_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_xx_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_xx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yyy_x_x, g_0_yyy_x_y, g_0_yyy_x_z, g_x_y_0_0_x_yy_x_x, g_x_y_0_0_x_yy_x_y, g_x_y_0_0_x_yy_x_z, g_xx_y_x_x, g_xx_y_x_y, g_xx_y_x_z, g_xx_yyy_x_x, g_xx_yyy_x_y, g_xx_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yy_x_x[i] = 2.0 * g_0_y_x_x[i] - 2.0 * g_0_yyy_x_x[i] * b_exp - 4.0 * g_xx_y_x_x[i] * a_exp + 4.0 * g_xx_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_x_y[i] = 2.0 * g_0_y_x_y[i] - 2.0 * g_0_yyy_x_y[i] * b_exp - 4.0 * g_xx_y_x_y[i] * a_exp + 4.0 * g_xx_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_x_z[i] = 2.0 * g_0_y_x_z[i] - 2.0 * g_0_yyy_x_z[i] * b_exp - 4.0 * g_xx_y_x_z[i] * a_exp + 4.0 * g_xx_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yyy_y_x, g_0_yyy_y_y, g_0_yyy_y_z, g_x_y_0_0_x_yy_y_x, g_x_y_0_0_x_yy_y_y, g_x_y_0_0_x_yy_y_z, g_xx_y_y_x, g_xx_y_y_y, g_xx_y_y_z, g_xx_yyy_y_x, g_xx_yyy_y_y, g_xx_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yy_y_x[i] = 2.0 * g_0_y_y_x[i] - 2.0 * g_0_yyy_y_x[i] * b_exp - 4.0 * g_xx_y_y_x[i] * a_exp + 4.0 * g_xx_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_y_y[i] = 2.0 * g_0_y_y_y[i] - 2.0 * g_0_yyy_y_y[i] * b_exp - 4.0 * g_xx_y_y_y[i] * a_exp + 4.0 * g_xx_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_y_z[i] = 2.0 * g_0_y_y_z[i] - 2.0 * g_0_yyy_y_z[i] * b_exp - 4.0 * g_xx_y_y_z[i] * a_exp + 4.0 * g_xx_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yyy_z_x, g_0_yyy_z_y, g_0_yyy_z_z, g_x_y_0_0_x_yy_z_x, g_x_y_0_0_x_yy_z_y, g_x_y_0_0_x_yy_z_z, g_xx_y_z_x, g_xx_y_z_y, g_xx_y_z_z, g_xx_yyy_z_x, g_xx_yyy_z_y, g_xx_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yy_z_x[i] = 2.0 * g_0_y_z_x[i] - 2.0 * g_0_yyy_z_x[i] * b_exp - 4.0 * g_xx_y_z_x[i] * a_exp + 4.0 * g_xx_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_z_y[i] = 2.0 * g_0_y_z_y[i] - 2.0 * g_0_yyy_z_y[i] * b_exp - 4.0 * g_xx_y_z_y[i] * a_exp + 4.0 * g_xx_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_z_z[i] = 2.0 * g_0_y_z_z[i] - 2.0 * g_0_yyy_z_z[i] * b_exp - 4.0 * g_xx_y_z_z[i] * a_exp + 4.0 * g_xx_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_x_y_0_0_x_yz_x_x, g_x_y_0_0_x_yz_x_y, g_x_y_0_0_x_yz_x_z, g_xx_yyz_x_x, g_xx_yyz_x_y, g_xx_yyz_x_z, g_xx_z_x_x, g_xx_z_x_y, g_xx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_yyz_x_x[i] * b_exp - 2.0 * g_xx_z_x_x[i] * a_exp + 4.0 * g_xx_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_yyz_x_y[i] * b_exp - 2.0 * g_xx_z_x_y[i] * a_exp + 4.0 * g_xx_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_yyz_x_z[i] * b_exp - 2.0 * g_xx_z_x_z[i] * a_exp + 4.0 * g_xx_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_x_y_0_0_x_yz_y_x, g_x_y_0_0_x_yz_y_y, g_x_y_0_0_x_yz_y_z, g_xx_yyz_y_x, g_xx_yyz_y_y, g_xx_yyz_y_z, g_xx_z_y_x, g_xx_z_y_y, g_xx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_yyz_y_x[i] * b_exp - 2.0 * g_xx_z_y_x[i] * a_exp + 4.0 * g_xx_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_yyz_y_y[i] * b_exp - 2.0 * g_xx_z_y_y[i] * a_exp + 4.0 * g_xx_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_yyz_y_z[i] * b_exp - 2.0 * g_xx_z_y_z[i] * a_exp + 4.0 * g_xx_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_x_y_0_0_x_yz_z_x, g_x_y_0_0_x_yz_z_y, g_x_y_0_0_x_yz_z_z, g_xx_yyz_z_x, g_xx_yyz_z_y, g_xx_yyz_z_z, g_xx_z_z_x, g_xx_z_z_y, g_xx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_yyz_z_x[i] * b_exp - 2.0 * g_xx_z_z_x[i] * a_exp + 4.0 * g_xx_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_yyz_z_y[i] * b_exp - 2.0 * g_xx_z_z_y[i] * a_exp + 4.0 * g_xx_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_yyz_z_z[i] * b_exp - 2.0 * g_xx_z_z_z[i] * a_exp + 4.0 * g_xx_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_x_y_0_0_x_zz_x_x, g_x_y_0_0_x_zz_x_y, g_x_y_0_0_x_zz_x_z, g_xx_yzz_x_x, g_xx_yzz_x_y, g_xx_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_zz_x_x[i] = -2.0 * g_0_yzz_x_x[i] * b_exp + 4.0 * g_xx_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_x_y[i] = -2.0 * g_0_yzz_x_y[i] * b_exp + 4.0 * g_xx_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_x_z[i] = -2.0 * g_0_yzz_x_z[i] * b_exp + 4.0 * g_xx_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_x_y_0_0_x_zz_y_x, g_x_y_0_0_x_zz_y_y, g_x_y_0_0_x_zz_y_z, g_xx_yzz_y_x, g_xx_yzz_y_y, g_xx_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_zz_y_x[i] = -2.0 * g_0_yzz_y_x[i] * b_exp + 4.0 * g_xx_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_y_y[i] = -2.0 * g_0_yzz_y_y[i] * b_exp + 4.0 * g_xx_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_y_z[i] = -2.0 * g_0_yzz_y_z[i] * b_exp + 4.0 * g_xx_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_x_y_0_0_x_zz_z_x, g_x_y_0_0_x_zz_z_y, g_x_y_0_0_x_zz_z_z, g_xx_yzz_z_x, g_xx_yzz_z_y, g_xx_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_zz_z_x[i] = -2.0 * g_0_yzz_z_x[i] * b_exp + 4.0 * g_xx_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_z_y[i] = -2.0 * g_0_yzz_z_y[i] * b_exp + 4.0 * g_xx_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_z_z[i] = -2.0 * g_0_yzz_z_z[i] * b_exp + 4.0 * g_xx_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_x_x, g_x_y_0_0_y_xx_x_y, g_x_y_0_0_y_xx_x_z, g_xy_xxy_x_x, g_xy_xxy_x_y, g_xy_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_x_x[i] = 4.0 * g_xy_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_x_y[i] = 4.0 * g_xy_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_x_z[i] = 4.0 * g_xy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_y_x, g_x_y_0_0_y_xx_y_y, g_x_y_0_0_y_xx_y_z, g_xy_xxy_y_x, g_xy_xxy_y_y, g_xy_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_y_x[i] = 4.0 * g_xy_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_y_y[i] = 4.0 * g_xy_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_y_z[i] = 4.0 * g_xy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_z_x, g_x_y_0_0_y_xx_z_y, g_x_y_0_0_y_xx_z_z, g_xy_xxy_z_x, g_xy_xxy_z_y, g_xy_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_z_x[i] = 4.0 * g_xy_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_z_y[i] = 4.0 * g_xy_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_z_z[i] = 4.0 * g_xy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_y_0_0_y_xy_x_x, g_x_y_0_0_y_xy_x_y, g_x_y_0_0_y_xy_x_z, g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xyy_x_x, g_xy_xyy_x_y, g_xy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xy_x_x[i] = -2.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_x_y[i] = -2.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_x_z[i] = -2.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_y_0_0_y_xy_y_x, g_x_y_0_0_y_xy_y_y, g_x_y_0_0_y_xy_y_z, g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xyy_y_x, g_xy_xyy_y_y, g_xy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xy_y_x[i] = -2.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_y_y[i] = -2.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_y_z[i] = -2.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_y_0_0_y_xy_z_x, g_x_y_0_0_y_xy_z_y, g_x_y_0_0_y_xy_z_z, g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xyy_z_x, g_xy_xyy_z_y, g_xy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xy_z_x[i] = -2.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_z_y[i] = -2.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_z_z[i] = -2.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_y_0_0_y_xz_x_x, g_x_y_0_0_y_xz_x_y, g_x_y_0_0_y_xz_x_z, g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xz_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_y_0_0_y_xz_y_x, g_x_y_0_0_y_xz_y_y, g_x_y_0_0_y_xz_y_z, g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xz_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_y_0_0_y_xz_z_x, g_x_y_0_0_y_xz_z_y, g_x_y_0_0_y_xz_z_z, g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xz_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_y_0_0_y_yy_x_x, g_x_y_0_0_y_yy_x_y, g_x_y_0_0_y_yy_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_xy_yyy_x_x, g_xy_yyy_x_y, g_xy_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yy_x_x[i] = -4.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_x_y[i] = -4.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_x_z[i] = -4.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_y_0_0_y_yy_y_x, g_x_y_0_0_y_yy_y_y, g_x_y_0_0_y_yy_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_xy_yyy_y_x, g_xy_yyy_y_y, g_xy_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yy_y_x[i] = -4.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_y_y[i] = -4.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_y_z[i] = -4.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_y_0_0_y_yy_z_x, g_x_y_0_0_y_yy_z_y, g_x_y_0_0_y_yy_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_xy_yyy_z_x, g_xy_yyy_z_y, g_xy_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yy_z_x[i] = -4.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_z_y[i] = -4.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_z_z[i] = -4.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_y_0_0_y_yz_x_x, g_x_y_0_0_y_yz_x_y, g_x_y_0_0_y_yz_x_z, g_xy_yyz_x_x, g_xy_yyz_x_y, g_xy_yyz_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yz_x_x[i] = -2.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_x_y[i] = -2.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_x_z[i] = -2.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_y_0_0_y_yz_y_x, g_x_y_0_0_y_yz_y_y, g_x_y_0_0_y_yz_y_z, g_xy_yyz_y_x, g_xy_yyz_y_y, g_xy_yyz_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yz_y_x[i] = -2.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_y_y[i] = -2.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_y_z[i] = -2.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_y_0_0_y_yz_z_x, g_x_y_0_0_y_yz_z_y, g_x_y_0_0_y_yz_z_z, g_xy_yyz_z_x, g_xy_yyz_z_y, g_xy_yyz_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yz_z_x[i] = -2.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_z_y[i] = -2.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_z_z[i] = -2.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_y_0_0_y_zz_x_x, g_x_y_0_0_y_zz_x_y, g_x_y_0_0_y_zz_x_z, g_xy_yzz_x_x, g_xy_yzz_x_y, g_xy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_zz_x_x[i] = 4.0 * g_xy_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_x_y[i] = 4.0 * g_xy_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_x_z[i] = 4.0 * g_xy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_y_0_0_y_zz_y_x, g_x_y_0_0_y_zz_y_y, g_x_y_0_0_y_zz_y_z, g_xy_yzz_y_x, g_xy_yzz_y_y, g_xy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_zz_y_x[i] = 4.0 * g_xy_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_y_y[i] = 4.0 * g_xy_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_y_z[i] = 4.0 * g_xy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_y_0_0_y_zz_z_x, g_x_y_0_0_y_zz_z_y, g_x_y_0_0_y_zz_z_z, g_xy_yzz_z_x, g_xy_yzz_z_y, g_xy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_zz_z_x[i] = 4.0 * g_xy_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_z_y[i] = 4.0 * g_xy_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_z_z[i] = 4.0 * g_xy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_x_x, g_x_y_0_0_z_xx_x_y, g_x_y_0_0_z_xx_x_z, g_xz_xxy_x_x, g_xz_xxy_x_y, g_xz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_x_x[i] = 4.0 * g_xz_xxy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_x_y[i] = 4.0 * g_xz_xxy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_x_z[i] = 4.0 * g_xz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_y_x, g_x_y_0_0_z_xx_y_y, g_x_y_0_0_z_xx_y_z, g_xz_xxy_y_x, g_xz_xxy_y_y, g_xz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_y_x[i] = 4.0 * g_xz_xxy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_y_y[i] = 4.0 * g_xz_xxy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_y_z[i] = 4.0 * g_xz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_z_x, g_x_y_0_0_z_xx_z_y, g_x_y_0_0_z_xx_z_z, g_xz_xxy_z_x, g_xz_xxy_z_y, g_xz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_z_x[i] = 4.0 * g_xz_xxy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_z_y[i] = 4.0 * g_xz_xxy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_z_z[i] = 4.0 * g_xz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_y_0_0_z_xy_x_x, g_x_y_0_0_z_xy_x_y, g_x_y_0_0_z_xy_x_z, g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xyy_x_x, g_xz_xyy_x_y, g_xz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xy_x_x[i] = -2.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_x_y[i] = -2.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_x_z[i] = -2.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_y_0_0_z_xy_y_x, g_x_y_0_0_z_xy_y_y, g_x_y_0_0_z_xy_y_z, g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xyy_y_x, g_xz_xyy_y_y, g_xz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xy_y_x[i] = -2.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_y_y[i] = -2.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_y_z[i] = -2.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_y_0_0_z_xy_z_x, g_x_y_0_0_z_xy_z_y, g_x_y_0_0_z_xy_z_z, g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xyy_z_x, g_xz_xyy_z_y, g_xz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xy_z_x[i] = -2.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_z_y[i] = -2.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_z_z[i] = -2.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_x_y_0_0_z_xz_x_x, g_x_y_0_0_z_xz_x_y, g_x_y_0_0_z_xz_x_z, g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xz_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_x_y_0_0_z_xz_y_x, g_x_y_0_0_z_xz_y_y, g_x_y_0_0_z_xz_y_z, g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xz_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_x_y_0_0_z_xz_z_x, g_x_y_0_0_z_xz_z_y, g_x_y_0_0_z_xz_z_z, g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xz_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_x_y_0_0_z_yy_x_x, g_x_y_0_0_z_yy_x_y, g_x_y_0_0_z_yy_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_xz_yyy_x_x, g_xz_yyy_x_y, g_xz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yy_x_x[i] = -4.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_yyy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_x_y[i] = -4.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_yyy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_x_z[i] = -4.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_x_y_0_0_z_yy_y_x, g_x_y_0_0_z_yy_y_y, g_x_y_0_0_z_yy_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_xz_yyy_y_x, g_xz_yyy_y_y, g_xz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yy_y_x[i] = -4.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_yyy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_y_y[i] = -4.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_yyy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_y_z[i] = -4.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_x_y_0_0_z_yy_z_x, g_x_y_0_0_z_yy_z_y, g_x_y_0_0_z_yy_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_xz_yyy_z_x, g_xz_yyy_z_y, g_xz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yy_z_x[i] = -4.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_yyy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_z_y[i] = -4.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_yyy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_z_z[i] = -4.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_x_y_0_0_z_yz_x_x, g_x_y_0_0_z_yz_x_y, g_x_y_0_0_z_yz_x_z, g_xz_yyz_x_x, g_xz_yyz_x_y, g_xz_yyz_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yz_x_x[i] = -2.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_yyz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_x_y[i] = -2.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_yyz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_x_z[i] = -2.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_x_y_0_0_z_yz_y_x, g_x_y_0_0_z_yz_y_y, g_x_y_0_0_z_yz_y_z, g_xz_yyz_y_x, g_xz_yyz_y_y, g_xz_yyz_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yz_y_x[i] = -2.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_yyz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_y_y[i] = -2.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_yyz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_y_z[i] = -2.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_x_y_0_0_z_yz_z_x, g_x_y_0_0_z_yz_z_y, g_x_y_0_0_z_yz_z_z, g_xz_yyz_z_x, g_xz_yyz_z_y, g_xz_yyz_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yz_z_x[i] = -2.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_yyz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_z_y[i] = -2.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_yyz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_z_z[i] = -2.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_x_y_0_0_z_zz_x_x, g_x_y_0_0_z_zz_x_y, g_x_y_0_0_z_zz_x_z, g_xz_yzz_x_x, g_xz_yzz_x_y, g_xz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_zz_x_x[i] = 4.0 * g_xz_yzz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_x_y[i] = 4.0 * g_xz_yzz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_x_z[i] = 4.0 * g_xz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_x_y_0_0_z_zz_y_x, g_x_y_0_0_z_zz_y_y, g_x_y_0_0_z_zz_y_z, g_xz_yzz_y_x, g_xz_yzz_y_y, g_xz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_zz_y_x[i] = 4.0 * g_xz_yzz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_y_y[i] = 4.0 * g_xz_yzz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_y_z[i] = 4.0 * g_xz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_x_y_0_0_z_zz_z_x, g_x_y_0_0_z_zz_z_y, g_x_y_0_0_z_zz_z_z, g_xz_yzz_z_x, g_xz_yzz_z_y, g_xz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_zz_z_x[i] = 4.0 * g_xz_yzz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_z_y[i] = 4.0 * g_xz_yzz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_z_z[i] = 4.0 * g_xz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_x_z_0_0_x_xx_x_x, g_x_z_0_0_x_xx_x_y, g_x_z_0_0_x_xx_x_z, g_xx_xxz_x_x, g_xx_xxz_x_y, g_xx_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_x_x[i] = -2.0 * g_0_xxz_x_x[i] * b_exp + 4.0 * g_xx_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_x_y[i] = -2.0 * g_0_xxz_x_y[i] * b_exp + 4.0 * g_xx_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_x_z[i] = -2.0 * g_0_xxz_x_z[i] * b_exp + 4.0 * g_xx_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_x_z_0_0_x_xx_y_x, g_x_z_0_0_x_xx_y_y, g_x_z_0_0_x_xx_y_z, g_xx_xxz_y_x, g_xx_xxz_y_y, g_xx_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_y_x[i] = -2.0 * g_0_xxz_y_x[i] * b_exp + 4.0 * g_xx_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_y_y[i] = -2.0 * g_0_xxz_y_y[i] * b_exp + 4.0 * g_xx_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_y_z[i] = -2.0 * g_0_xxz_y_z[i] * b_exp + 4.0 * g_xx_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_x_z_0_0_x_xx_z_x, g_x_z_0_0_x_xx_z_y, g_x_z_0_0_x_xx_z_z, g_xx_xxz_z_x, g_xx_xxz_z_y, g_xx_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_z_x[i] = -2.0 * g_0_xxz_z_x[i] * b_exp + 4.0 * g_xx_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_z_y[i] = -2.0 * g_0_xxz_z_y[i] * b_exp + 4.0 * g_xx_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_z_z[i] = -2.0 * g_0_xxz_z_z[i] * b_exp + 4.0 * g_xx_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_x_z_0_0_x_xy_x_x, g_x_z_0_0_x_xy_x_y, g_x_z_0_0_x_xy_x_z, g_xx_xyz_x_x, g_xx_xyz_x_y, g_xx_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xy_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_xx_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_xx_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_xx_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_x_z_0_0_x_xy_y_x, g_x_z_0_0_x_xy_y_y, g_x_z_0_0_x_xy_y_z, g_xx_xyz_y_x, g_xx_xyz_y_y, g_xx_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xy_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_xx_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_xx_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_xx_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_x_z_0_0_x_xy_z_x, g_x_z_0_0_x_xy_z_y, g_x_z_0_0_x_xy_z_z, g_xx_xyz_z_x, g_xx_xyz_z_y, g_xx_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xy_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_xx_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_xx_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_xx_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_x_z_0_0_x_xz_x_x, g_x_z_0_0_x_xz_x_y, g_x_z_0_0_x_xz_x_z, g_xx_x_x_x, g_xx_x_x_y, g_xx_x_x_z, g_xx_xzz_x_x, g_xx_xzz_x_y, g_xx_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xz_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xzz_x_x[i] * b_exp - 2.0 * g_xx_x_x_x[i] * a_exp + 4.0 * g_xx_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xzz_x_y[i] * b_exp - 2.0 * g_xx_x_x_y[i] * a_exp + 4.0 * g_xx_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xzz_x_z[i] * b_exp - 2.0 * g_xx_x_x_z[i] * a_exp + 4.0 * g_xx_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_x_z_0_0_x_xz_y_x, g_x_z_0_0_x_xz_y_y, g_x_z_0_0_x_xz_y_z, g_xx_x_y_x, g_xx_x_y_y, g_xx_x_y_z, g_xx_xzz_y_x, g_xx_xzz_y_y, g_xx_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xz_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xzz_y_x[i] * b_exp - 2.0 * g_xx_x_y_x[i] * a_exp + 4.0 * g_xx_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xzz_y_y[i] * b_exp - 2.0 * g_xx_x_y_y[i] * a_exp + 4.0 * g_xx_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xzz_y_z[i] * b_exp - 2.0 * g_xx_x_y_z[i] * a_exp + 4.0 * g_xx_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_x_z_0_0_x_xz_z_x, g_x_z_0_0_x_xz_z_y, g_x_z_0_0_x_xz_z_z, g_xx_x_z_x, g_xx_x_z_y, g_xx_x_z_z, g_xx_xzz_z_x, g_xx_xzz_z_y, g_xx_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xz_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xzz_z_x[i] * b_exp - 2.0 * g_xx_x_z_x[i] * a_exp + 4.0 * g_xx_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xzz_z_y[i] * b_exp - 2.0 * g_xx_x_z_y[i] * a_exp + 4.0 * g_xx_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xzz_z_z[i] * b_exp - 2.0 * g_xx_x_z_z[i] * a_exp + 4.0 * g_xx_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_x_z_0_0_x_yy_x_x, g_x_z_0_0_x_yy_x_y, g_x_z_0_0_x_yy_x_z, g_xx_yyz_x_x, g_xx_yyz_x_y, g_xx_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yy_x_x[i] = -2.0 * g_0_yyz_x_x[i] * b_exp + 4.0 * g_xx_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_x_y[i] = -2.0 * g_0_yyz_x_y[i] * b_exp + 4.0 * g_xx_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_x_z[i] = -2.0 * g_0_yyz_x_z[i] * b_exp + 4.0 * g_xx_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_x_z_0_0_x_yy_y_x, g_x_z_0_0_x_yy_y_y, g_x_z_0_0_x_yy_y_z, g_xx_yyz_y_x, g_xx_yyz_y_y, g_xx_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yy_y_x[i] = -2.0 * g_0_yyz_y_x[i] * b_exp + 4.0 * g_xx_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_y_y[i] = -2.0 * g_0_yyz_y_y[i] * b_exp + 4.0 * g_xx_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_y_z[i] = -2.0 * g_0_yyz_y_z[i] * b_exp + 4.0 * g_xx_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_x_z_0_0_x_yy_z_x, g_x_z_0_0_x_yy_z_y, g_x_z_0_0_x_yy_z_z, g_xx_yyz_z_x, g_xx_yyz_z_y, g_xx_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yy_z_x[i] = -2.0 * g_0_yyz_z_x[i] * b_exp + 4.0 * g_xx_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_z_y[i] = -2.0 * g_0_yyz_z_y[i] * b_exp + 4.0 * g_xx_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_z_z[i] = -2.0 * g_0_yyz_z_z[i] * b_exp + 4.0 * g_xx_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_x_z_0_0_x_yz_x_x, g_x_z_0_0_x_yz_x_y, g_x_z_0_0_x_yz_x_z, g_xx_y_x_x, g_xx_y_x_y, g_xx_y_x_z, g_xx_yzz_x_x, g_xx_yzz_x_y, g_xx_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yz_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_yzz_x_x[i] * b_exp - 2.0 * g_xx_y_x_x[i] * a_exp + 4.0 * g_xx_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_yzz_x_y[i] * b_exp - 2.0 * g_xx_y_x_y[i] * a_exp + 4.0 * g_xx_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_yzz_x_z[i] * b_exp - 2.0 * g_xx_y_x_z[i] * a_exp + 4.0 * g_xx_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_x_z_0_0_x_yz_y_x, g_x_z_0_0_x_yz_y_y, g_x_z_0_0_x_yz_y_z, g_xx_y_y_x, g_xx_y_y_y, g_xx_y_y_z, g_xx_yzz_y_x, g_xx_yzz_y_y, g_xx_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yz_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_yzz_y_x[i] * b_exp - 2.0 * g_xx_y_y_x[i] * a_exp + 4.0 * g_xx_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_yzz_y_y[i] * b_exp - 2.0 * g_xx_y_y_y[i] * a_exp + 4.0 * g_xx_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_yzz_y_z[i] * b_exp - 2.0 * g_xx_y_y_z[i] * a_exp + 4.0 * g_xx_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_x_z_0_0_x_yz_z_x, g_x_z_0_0_x_yz_z_y, g_x_z_0_0_x_yz_z_z, g_xx_y_z_x, g_xx_y_z_y, g_xx_y_z_z, g_xx_yzz_z_x, g_xx_yzz_z_y, g_xx_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yz_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_yzz_z_x[i] * b_exp - 2.0 * g_xx_y_z_x[i] * a_exp + 4.0 * g_xx_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_yzz_z_y[i] * b_exp - 2.0 * g_xx_y_z_y[i] * a_exp + 4.0 * g_xx_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_yzz_z_z[i] * b_exp - 2.0 * g_xx_y_z_z[i] * a_exp + 4.0 * g_xx_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_0_zzz_x_x, g_0_zzz_x_y, g_0_zzz_x_z, g_x_z_0_0_x_zz_x_x, g_x_z_0_0_x_zz_x_y, g_x_z_0_0_x_zz_x_z, g_xx_z_x_x, g_xx_z_x_y, g_xx_z_x_z, g_xx_zzz_x_x, g_xx_zzz_x_y, g_xx_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_zz_x_x[i] = 2.0 * g_0_z_x_x[i] - 2.0 * g_0_zzz_x_x[i] * b_exp - 4.0 * g_xx_z_x_x[i] * a_exp + 4.0 * g_xx_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_x_y[i] = 2.0 * g_0_z_x_y[i] - 2.0 * g_0_zzz_x_y[i] * b_exp - 4.0 * g_xx_z_x_y[i] * a_exp + 4.0 * g_xx_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_x_z[i] = 2.0 * g_0_z_x_z[i] - 2.0 * g_0_zzz_x_z[i] * b_exp - 4.0 * g_xx_z_x_z[i] * a_exp + 4.0 * g_xx_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_0_zzz_y_x, g_0_zzz_y_y, g_0_zzz_y_z, g_x_z_0_0_x_zz_y_x, g_x_z_0_0_x_zz_y_y, g_x_z_0_0_x_zz_y_z, g_xx_z_y_x, g_xx_z_y_y, g_xx_z_y_z, g_xx_zzz_y_x, g_xx_zzz_y_y, g_xx_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_zz_y_x[i] = 2.0 * g_0_z_y_x[i] - 2.0 * g_0_zzz_y_x[i] * b_exp - 4.0 * g_xx_z_y_x[i] * a_exp + 4.0 * g_xx_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_y_y[i] = 2.0 * g_0_z_y_y[i] - 2.0 * g_0_zzz_y_y[i] * b_exp - 4.0 * g_xx_z_y_y[i] * a_exp + 4.0 * g_xx_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_y_z[i] = 2.0 * g_0_z_y_z[i] - 2.0 * g_0_zzz_y_z[i] * b_exp - 4.0 * g_xx_z_y_z[i] * a_exp + 4.0 * g_xx_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_0_zzz_z_x, g_0_zzz_z_y, g_0_zzz_z_z, g_x_z_0_0_x_zz_z_x, g_x_z_0_0_x_zz_z_y, g_x_z_0_0_x_zz_z_z, g_xx_z_z_x, g_xx_z_z_y, g_xx_z_z_z, g_xx_zzz_z_x, g_xx_zzz_z_y, g_xx_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_zz_z_x[i] = 2.0 * g_0_z_z_x[i] - 2.0 * g_0_zzz_z_x[i] * b_exp - 4.0 * g_xx_z_z_x[i] * a_exp + 4.0 * g_xx_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_z_y[i] = 2.0 * g_0_z_z_y[i] - 2.0 * g_0_zzz_z_y[i] * b_exp - 4.0 * g_xx_z_z_y[i] * a_exp + 4.0 * g_xx_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_z_z[i] = 2.0 * g_0_z_z_z[i] - 2.0 * g_0_zzz_z_z[i] * b_exp - 4.0 * g_xx_z_z_z[i] * a_exp + 4.0 * g_xx_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_x_x, g_x_z_0_0_y_xx_x_y, g_x_z_0_0_y_xx_x_z, g_xy_xxz_x_x, g_xy_xxz_x_y, g_xy_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_x_x[i] = 4.0 * g_xy_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_x_y[i] = 4.0 * g_xy_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_x_z[i] = 4.0 * g_xy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_y_x, g_x_z_0_0_y_xx_y_y, g_x_z_0_0_y_xx_y_z, g_xy_xxz_y_x, g_xy_xxz_y_y, g_xy_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_y_x[i] = 4.0 * g_xy_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_y_y[i] = 4.0 * g_xy_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_y_z[i] = 4.0 * g_xy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_z_x, g_x_z_0_0_y_xx_z_y, g_x_z_0_0_y_xx_z_z, g_xy_xxz_z_x, g_xy_xxz_z_y, g_xy_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_z_x[i] = 4.0 * g_xy_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_z_y[i] = 4.0 * g_xy_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_z_z[i] = 4.0 * g_xy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_x_z_0_0_y_xy_x_x, g_x_z_0_0_y_xy_x_y, g_x_z_0_0_y_xy_x_z, g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xy_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_x_z_0_0_y_xy_y_x, g_x_z_0_0_y_xy_y_y, g_x_z_0_0_y_xy_y_z, g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xy_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_x_z_0_0_y_xy_z_x, g_x_z_0_0_y_xy_z_y, g_x_z_0_0_y_xy_z_z, g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xy_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_x_z_0_0_y_xz_x_x, g_x_z_0_0_y_xz_x_y, g_x_z_0_0_y_xz_x_z, g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xzz_x_x, g_xy_xzz_x_y, g_xy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xz_x_x[i] = -2.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_x_y[i] = -2.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_x_z[i] = -2.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_x_z_0_0_y_xz_y_x, g_x_z_0_0_y_xz_y_y, g_x_z_0_0_y_xz_y_z, g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xzz_y_x, g_xy_xzz_y_y, g_xy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xz_y_x[i] = -2.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_y_y[i] = -2.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_y_z[i] = -2.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_x_z_0_0_y_xz_z_x, g_x_z_0_0_y_xz_z_y, g_x_z_0_0_y_xz_z_z, g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xzz_z_x, g_xy_xzz_z_y, g_xy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xz_z_x[i] = -2.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_z_y[i] = -2.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_z_z[i] = -2.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_x_z_0_0_y_yy_x_x, g_x_z_0_0_y_yy_x_y, g_x_z_0_0_y_yy_x_z, g_xy_yyz_x_x, g_xy_yyz_x_y, g_xy_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yy_x_x[i] = 4.0 * g_xy_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_x_y[i] = 4.0 * g_xy_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_x_z[i] = 4.0 * g_xy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_x_z_0_0_y_yy_y_x, g_x_z_0_0_y_yy_y_y, g_x_z_0_0_y_yy_y_z, g_xy_yyz_y_x, g_xy_yyz_y_y, g_xy_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yy_y_x[i] = 4.0 * g_xy_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_y_y[i] = 4.0 * g_xy_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_y_z[i] = 4.0 * g_xy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_x_z_0_0_y_yy_z_x, g_x_z_0_0_y_yy_z_y, g_x_z_0_0_y_yy_z_z, g_xy_yyz_z_x, g_xy_yyz_z_y, g_xy_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yy_z_x[i] = 4.0 * g_xy_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_z_y[i] = 4.0 * g_xy_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_z_z[i] = 4.0 * g_xy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_x_z_0_0_y_yz_x_x, g_x_z_0_0_y_yz_x_y, g_x_z_0_0_y_yz_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_xy_yzz_x_x, g_xy_yzz_x_y, g_xy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yz_x_x[i] = -2.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_x_y[i] = -2.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_x_z[i] = -2.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_x_z_0_0_y_yz_y_x, g_x_z_0_0_y_yz_y_y, g_x_z_0_0_y_yz_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_xy_yzz_y_x, g_xy_yzz_y_y, g_xy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yz_y_x[i] = -2.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_y_y[i] = -2.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_y_z[i] = -2.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_x_z_0_0_y_yz_z_x, g_x_z_0_0_y_yz_z_y, g_x_z_0_0_y_yz_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_xy_yzz_z_x, g_xy_yzz_z_y, g_xy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yz_z_x[i] = -2.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_z_y[i] = -2.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_z_z[i] = -2.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_x_z_0_0_y_zz_x_x, g_x_z_0_0_y_zz_x_y, g_x_z_0_0_y_zz_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z, g_xy_zzz_x_x, g_xy_zzz_x_y, g_xy_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_zz_x_x[i] = -4.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_x_y[i] = -4.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_x_z[i] = -4.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_x_z_0_0_y_zz_y_x, g_x_z_0_0_y_zz_y_y, g_x_z_0_0_y_zz_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z, g_xy_zzz_y_x, g_xy_zzz_y_y, g_xy_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_zz_y_x[i] = -4.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_y_y[i] = -4.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_y_z[i] = -4.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_x_z_0_0_y_zz_z_x, g_x_z_0_0_y_zz_z_y, g_x_z_0_0_y_zz_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z, g_xy_zzz_z_x, g_xy_zzz_z_y, g_xy_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_zz_z_x[i] = -4.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_z_y[i] = -4.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_z_z[i] = -4.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_x_x, g_x_z_0_0_z_xx_x_y, g_x_z_0_0_z_xx_x_z, g_xz_xxz_x_x, g_xz_xxz_x_y, g_xz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_x_x[i] = 4.0 * g_xz_xxz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_x_y[i] = 4.0 * g_xz_xxz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_x_z[i] = 4.0 * g_xz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_y_x, g_x_z_0_0_z_xx_y_y, g_x_z_0_0_z_xx_y_z, g_xz_xxz_y_x, g_xz_xxz_y_y, g_xz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_y_x[i] = 4.0 * g_xz_xxz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_y_y[i] = 4.0 * g_xz_xxz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_y_z[i] = 4.0 * g_xz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_z_x, g_x_z_0_0_z_xx_z_y, g_x_z_0_0_z_xx_z_z, g_xz_xxz_z_x, g_xz_xxz_z_y, g_xz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_z_x[i] = 4.0 * g_xz_xxz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_z_y[i] = 4.0 * g_xz_xxz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_z_z[i] = 4.0 * g_xz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_x_z_0_0_z_xy_x_x, g_x_z_0_0_z_xy_x_y, g_x_z_0_0_z_xy_x_z, g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xy_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_x_z_0_0_z_xy_y_x, g_x_z_0_0_z_xy_y_y, g_x_z_0_0_z_xy_y_z, g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xy_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_x_z_0_0_z_xy_z_x, g_x_z_0_0_z_xy_z_y, g_x_z_0_0_z_xy_z_z, g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xy_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_x_z_0_0_z_xz_x_x, g_x_z_0_0_z_xz_x_y, g_x_z_0_0_z_xz_x_z, g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xzz_x_x, g_xz_xzz_x_y, g_xz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xz_x_x[i] = -2.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_x_y[i] = -2.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_x_z[i] = -2.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_x_z_0_0_z_xz_y_x, g_x_z_0_0_z_xz_y_y, g_x_z_0_0_z_xz_y_z, g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xzz_y_x, g_xz_xzz_y_y, g_xz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xz_y_x[i] = -2.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_y_y[i] = -2.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_y_z[i] = -2.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_x_z_0_0_z_xz_z_x, g_x_z_0_0_z_xz_z_y, g_x_z_0_0_z_xz_z_z, g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xzz_z_x, g_xz_xzz_z_y, g_xz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xz_z_x[i] = -2.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_z_y[i] = -2.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_z_z[i] = -2.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_x_z_0_0_z_yy_x_x, g_x_z_0_0_z_yy_x_y, g_x_z_0_0_z_yy_x_z, g_xz_yyz_x_x, g_xz_yyz_x_y, g_xz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yy_x_x[i] = 4.0 * g_xz_yyz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_x_y[i] = 4.0 * g_xz_yyz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_x_z[i] = 4.0 * g_xz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_x_z_0_0_z_yy_y_x, g_x_z_0_0_z_yy_y_y, g_x_z_0_0_z_yy_y_z, g_xz_yyz_y_x, g_xz_yyz_y_y, g_xz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yy_y_x[i] = 4.0 * g_xz_yyz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_y_y[i] = 4.0 * g_xz_yyz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_y_z[i] = 4.0 * g_xz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_x_z_0_0_z_yy_z_x, g_x_z_0_0_z_yy_z_y, g_x_z_0_0_z_yy_z_z, g_xz_yyz_z_x, g_xz_yyz_z_y, g_xz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yy_z_x[i] = 4.0 * g_xz_yyz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_z_y[i] = 4.0 * g_xz_yyz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_z_z[i] = 4.0 * g_xz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_x_z_0_0_z_yz_x_x, g_x_z_0_0_z_yz_x_y, g_x_z_0_0_z_yz_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_xz_yzz_x_x, g_xz_yzz_x_y, g_xz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yz_x_x[i] = -2.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_yzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_x_y[i] = -2.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_yzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_x_z[i] = -2.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_x_z_0_0_z_yz_y_x, g_x_z_0_0_z_yz_y_y, g_x_z_0_0_z_yz_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_xz_yzz_y_x, g_xz_yzz_y_y, g_xz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yz_y_x[i] = -2.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_yzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_y_y[i] = -2.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_yzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_y_z[i] = -2.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_x_z_0_0_z_yz_z_x, g_x_z_0_0_z_yz_z_y, g_x_z_0_0_z_yz_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_xz_yzz_z_x, g_xz_yzz_z_y, g_xz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yz_z_x[i] = -2.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_yzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_z_y[i] = -2.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_yzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_z_z[i] = -2.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_x_z_0_0_z_zz_x_x, g_x_z_0_0_z_zz_x_y, g_x_z_0_0_z_zz_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z, g_xz_zzz_x_x, g_xz_zzz_x_y, g_xz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_zz_x_x[i] = -4.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_zzz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_x_y[i] = -4.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_zzz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_x_z[i] = -4.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_x_z_0_0_z_zz_y_x, g_x_z_0_0_z_zz_y_y, g_x_z_0_0_z_zz_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z, g_xz_zzz_y_x, g_xz_zzz_y_y, g_xz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_zz_y_x[i] = -4.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_zzz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_y_y[i] = -4.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_zzz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_y_z[i] = -4.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_x_z_0_0_z_zz_z_x, g_x_z_0_0_z_zz_z_y, g_x_z_0_0_z_zz_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z, g_xz_zzz_z_x, g_xz_zzz_z_y, g_xz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_zz_z_x[i] = -4.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_zzz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_z_y[i] = -4.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_zzz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_z_z[i] = -4.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xxx_x_x, g_xy_xxx_x_y, g_xy_xxx_x_z, g_y_x_0_0_x_xx_x_x, g_y_x_0_0_x_xx_x_y, g_y_x_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_x_x[i] = -4.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_x_y[i] = -4.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_x_z[i] = -4.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xxx_y_x, g_xy_xxx_y_y, g_xy_xxx_y_z, g_y_x_0_0_x_xx_y_x, g_y_x_0_0_x_xx_y_y, g_y_x_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_y_x[i] = -4.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_y_y[i] = -4.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_y_z[i] = -4.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xxx_z_x, g_xy_xxx_z_y, g_xy_xxx_z_z, g_y_x_0_0_x_xx_z_x, g_y_x_0_0_x_xx_z_y, g_y_x_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_z_x[i] = -4.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_z_y[i] = -4.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_z_z[i] = -4.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_xy_xxy_x_x, g_xy_xxy_x_y, g_xy_xxy_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_y_x_0_0_x_xy_x_x, g_y_x_0_0_x_xy_x_y, g_y_x_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xy_x_x[i] = -2.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_x_y[i] = -2.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_x_z[i] = -2.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_xy_xxy_y_x, g_xy_xxy_y_y, g_xy_xxy_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_y_x_0_0_x_xy_y_x, g_y_x_0_0_x_xy_y_y, g_y_x_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xy_y_x[i] = -2.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_y_y[i] = -2.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_y_z[i] = -2.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_xy_xxy_z_x, g_xy_xxy_z_y, g_xy_xxy_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_y_x_0_0_x_xy_z_x, g_y_x_0_0_x_xy_z_y, g_y_x_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xy_z_x[i] = -2.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_z_y[i] = -2.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_z_z[i] = -2.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_xy_xxz_x_x, g_xy_xxz_x_y, g_xy_xxz_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z, g_y_x_0_0_x_xz_x_x, g_y_x_0_0_x_xz_x_y, g_y_x_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xz_x_x[i] = -2.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_x_y[i] = -2.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_x_z[i] = -2.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_xy_xxz_y_x, g_xy_xxz_y_y, g_xy_xxz_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z, g_y_x_0_0_x_xz_y_x, g_y_x_0_0_x_xz_y_y, g_y_x_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xz_y_x[i] = -2.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_y_y[i] = -2.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_y_z[i] = -2.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_xy_xxz_z_x, g_xy_xxz_z_y, g_xy_xxz_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z, g_y_x_0_0_x_xz_z_x, g_y_x_0_0_x_xz_z_y, g_y_x_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xz_z_x[i] = -2.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_z_y[i] = -2.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_z_z[i] = -2.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_xy_xyy_x_x, g_xy_xyy_x_y, g_xy_xyy_x_z, g_y_x_0_0_x_yy_x_x, g_y_x_0_0_x_yy_x_y, g_y_x_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yy_x_x[i] = 4.0 * g_xy_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_x_y[i] = 4.0 * g_xy_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_x_z[i] = 4.0 * g_xy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_xy_xyy_y_x, g_xy_xyy_y_y, g_xy_xyy_y_z, g_y_x_0_0_x_yy_y_x, g_y_x_0_0_x_yy_y_y, g_y_x_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yy_y_x[i] = 4.0 * g_xy_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_y_y[i] = 4.0 * g_xy_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_y_z[i] = 4.0 * g_xy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_xy_xyy_z_x, g_xy_xyy_z_y, g_xy_xyy_z_z, g_y_x_0_0_x_yy_z_x, g_y_x_0_0_x_yy_z_y, g_y_x_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yy_z_x[i] = 4.0 * g_xy_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_z_y[i] = 4.0 * g_xy_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_z_z[i] = 4.0 * g_xy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z, g_y_x_0_0_x_yz_x_x, g_y_x_0_0_x_yz_x_y, g_y_x_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yz_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z, g_y_x_0_0_x_yz_y_x, g_y_x_0_0_x_yz_y_y, g_y_x_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yz_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z, g_y_x_0_0_x_yz_z_x, g_y_x_0_0_x_yz_z_y, g_y_x_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yz_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_xy_xzz_x_x, g_xy_xzz_x_y, g_xy_xzz_x_z, g_y_x_0_0_x_zz_x_x, g_y_x_0_0_x_zz_x_y, g_y_x_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_zz_x_x[i] = 4.0 * g_xy_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_x_y[i] = 4.0 * g_xy_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_x_z[i] = 4.0 * g_xy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_xy_xzz_y_x, g_xy_xzz_y_y, g_xy_xzz_y_z, g_y_x_0_0_x_zz_y_x, g_y_x_0_0_x_zz_y_y, g_y_x_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_zz_y_x[i] = 4.0 * g_xy_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_y_y[i] = 4.0 * g_xy_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_y_z[i] = 4.0 * g_xy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_xy_xzz_z_x, g_xy_xzz_z_y, g_xy_xzz_z_z, g_y_x_0_0_x_zz_z_x, g_y_x_0_0_x_zz_z_y, g_y_x_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_zz_z_x[i] = 4.0 * g_xy_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_z_y[i] = 4.0 * g_xy_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_z_z[i] = 4.0 * g_xy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xxx_x_x, g_0_xxx_x_y, g_0_xxx_x_z, g_y_x_0_0_y_xx_x_x, g_y_x_0_0_y_xx_x_y, g_y_x_0_0_y_xx_x_z, g_yy_x_x_x, g_yy_x_x_y, g_yy_x_x_z, g_yy_xxx_x_x, g_yy_xxx_x_y, g_yy_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_x_x[i] = 2.0 * g_0_x_x_x[i] - 2.0 * g_0_xxx_x_x[i] * b_exp - 4.0 * g_yy_x_x_x[i] * a_exp + 4.0 * g_yy_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_x_y[i] = 2.0 * g_0_x_x_y[i] - 2.0 * g_0_xxx_x_y[i] * b_exp - 4.0 * g_yy_x_x_y[i] * a_exp + 4.0 * g_yy_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_x_z[i] = 2.0 * g_0_x_x_z[i] - 2.0 * g_0_xxx_x_z[i] * b_exp - 4.0 * g_yy_x_x_z[i] * a_exp + 4.0 * g_yy_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xxx_y_x, g_0_xxx_y_y, g_0_xxx_y_z, g_y_x_0_0_y_xx_y_x, g_y_x_0_0_y_xx_y_y, g_y_x_0_0_y_xx_y_z, g_yy_x_y_x, g_yy_x_y_y, g_yy_x_y_z, g_yy_xxx_y_x, g_yy_xxx_y_y, g_yy_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_y_x[i] = 2.0 * g_0_x_y_x[i] - 2.0 * g_0_xxx_y_x[i] * b_exp - 4.0 * g_yy_x_y_x[i] * a_exp + 4.0 * g_yy_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_y_y[i] = 2.0 * g_0_x_y_y[i] - 2.0 * g_0_xxx_y_y[i] * b_exp - 4.0 * g_yy_x_y_y[i] * a_exp + 4.0 * g_yy_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_y_z[i] = 2.0 * g_0_x_y_z[i] - 2.0 * g_0_xxx_y_z[i] * b_exp - 4.0 * g_yy_x_y_z[i] * a_exp + 4.0 * g_yy_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xxx_z_x, g_0_xxx_z_y, g_0_xxx_z_z, g_y_x_0_0_y_xx_z_x, g_y_x_0_0_y_xx_z_y, g_y_x_0_0_y_xx_z_z, g_yy_x_z_x, g_yy_x_z_y, g_yy_x_z_z, g_yy_xxx_z_x, g_yy_xxx_z_y, g_yy_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_z_x[i] = 2.0 * g_0_x_z_x[i] - 2.0 * g_0_xxx_z_x[i] * b_exp - 4.0 * g_yy_x_z_x[i] * a_exp + 4.0 * g_yy_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_z_y[i] = 2.0 * g_0_x_z_y[i] - 2.0 * g_0_xxx_z_y[i] * b_exp - 4.0 * g_yy_x_z_y[i] * a_exp + 4.0 * g_yy_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_z_z[i] = 2.0 * g_0_x_z_z[i] - 2.0 * g_0_xxx_z_z[i] * b_exp - 4.0 * g_yy_x_z_z[i] * a_exp + 4.0 * g_yy_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_y_x_0_0_y_xy_x_x, g_y_x_0_0_y_xy_x_y, g_y_x_0_0_y_xy_x_z, g_yy_xxy_x_x, g_yy_xxy_x_y, g_yy_xxy_x_z, g_yy_y_x_x, g_yy_y_x_y, g_yy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xy_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_xxy_x_x[i] * b_exp - 2.0 * g_yy_y_x_x[i] * a_exp + 4.0 * g_yy_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_xxy_x_y[i] * b_exp - 2.0 * g_yy_y_x_y[i] * a_exp + 4.0 * g_yy_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_xxy_x_z[i] * b_exp - 2.0 * g_yy_y_x_z[i] * a_exp + 4.0 * g_yy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_y_x_0_0_y_xy_y_x, g_y_x_0_0_y_xy_y_y, g_y_x_0_0_y_xy_y_z, g_yy_xxy_y_x, g_yy_xxy_y_y, g_yy_xxy_y_z, g_yy_y_y_x, g_yy_y_y_y, g_yy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xy_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_xxy_y_x[i] * b_exp - 2.0 * g_yy_y_y_x[i] * a_exp + 4.0 * g_yy_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_xxy_y_y[i] * b_exp - 2.0 * g_yy_y_y_y[i] * a_exp + 4.0 * g_yy_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_xxy_y_z[i] * b_exp - 2.0 * g_yy_y_y_z[i] * a_exp + 4.0 * g_yy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_y_x_0_0_y_xy_z_x, g_y_x_0_0_y_xy_z_y, g_y_x_0_0_y_xy_z_z, g_yy_xxy_z_x, g_yy_xxy_z_y, g_yy_xxy_z_z, g_yy_y_z_x, g_yy_y_z_y, g_yy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xy_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_xxy_z_x[i] * b_exp - 2.0 * g_yy_y_z_x[i] * a_exp + 4.0 * g_yy_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_xxy_z_y[i] * b_exp - 2.0 * g_yy_y_z_y[i] * a_exp + 4.0 * g_yy_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_xxy_z_z[i] * b_exp - 2.0 * g_yy_y_z_z[i] * a_exp + 4.0 * g_yy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_y_x_0_0_y_xz_x_x, g_y_x_0_0_y_xz_x_y, g_y_x_0_0_y_xz_x_z, g_yy_xxz_x_x, g_yy_xxz_x_y, g_yy_xxz_x_z, g_yy_z_x_x, g_yy_z_x_y, g_yy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_xxz_x_x[i] * b_exp - 2.0 * g_yy_z_x_x[i] * a_exp + 4.0 * g_yy_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_xxz_x_y[i] * b_exp - 2.0 * g_yy_z_x_y[i] * a_exp + 4.0 * g_yy_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_xxz_x_z[i] * b_exp - 2.0 * g_yy_z_x_z[i] * a_exp + 4.0 * g_yy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_y_x_0_0_y_xz_y_x, g_y_x_0_0_y_xz_y_y, g_y_x_0_0_y_xz_y_z, g_yy_xxz_y_x, g_yy_xxz_y_y, g_yy_xxz_y_z, g_yy_z_y_x, g_yy_z_y_y, g_yy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_xxz_y_x[i] * b_exp - 2.0 * g_yy_z_y_x[i] * a_exp + 4.0 * g_yy_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_xxz_y_y[i] * b_exp - 2.0 * g_yy_z_y_y[i] * a_exp + 4.0 * g_yy_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_xxz_y_z[i] * b_exp - 2.0 * g_yy_z_y_z[i] * a_exp + 4.0 * g_yy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_y_x_0_0_y_xz_z_x, g_y_x_0_0_y_xz_z_y, g_y_x_0_0_y_xz_z_z, g_yy_xxz_z_x, g_yy_xxz_z_y, g_yy_xxz_z_z, g_yy_z_z_x, g_yy_z_z_y, g_yy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_xxz_z_x[i] * b_exp - 2.0 * g_yy_z_z_x[i] * a_exp + 4.0 * g_yy_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_xxz_z_y[i] * b_exp - 2.0 * g_yy_z_z_y[i] * a_exp + 4.0 * g_yy_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_xxz_z_z[i] * b_exp - 2.0 * g_yy_z_z_z[i] * a_exp + 4.0 * g_yy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_y_x_0_0_y_yy_x_x, g_y_x_0_0_y_yy_x_y, g_y_x_0_0_y_yy_x_z, g_yy_xyy_x_x, g_yy_xyy_x_y, g_yy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yy_x_x[i] = -2.0 * g_0_xyy_x_x[i] * b_exp + 4.0 * g_yy_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_x_y[i] = -2.0 * g_0_xyy_x_y[i] * b_exp + 4.0 * g_yy_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_x_z[i] = -2.0 * g_0_xyy_x_z[i] * b_exp + 4.0 * g_yy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_y_x_0_0_y_yy_y_x, g_y_x_0_0_y_yy_y_y, g_y_x_0_0_y_yy_y_z, g_yy_xyy_y_x, g_yy_xyy_y_y, g_yy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yy_y_x[i] = -2.0 * g_0_xyy_y_x[i] * b_exp + 4.0 * g_yy_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_y_y[i] = -2.0 * g_0_xyy_y_y[i] * b_exp + 4.0 * g_yy_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_y_z[i] = -2.0 * g_0_xyy_y_z[i] * b_exp + 4.0 * g_yy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_y_x_0_0_y_yy_z_x, g_y_x_0_0_y_yy_z_y, g_y_x_0_0_y_yy_z_z, g_yy_xyy_z_x, g_yy_xyy_z_y, g_yy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yy_z_x[i] = -2.0 * g_0_xyy_z_x[i] * b_exp + 4.0 * g_yy_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_z_y[i] = -2.0 * g_0_xyy_z_y[i] * b_exp + 4.0 * g_yy_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_z_z[i] = -2.0 * g_0_xyy_z_z[i] * b_exp + 4.0 * g_yy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_y_x_0_0_y_yz_x_x, g_y_x_0_0_y_yz_x_y, g_y_x_0_0_y_yz_x_z, g_yy_xyz_x_x, g_yy_xyz_x_y, g_yy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_yy_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_yy_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_yy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_y_x_0_0_y_yz_y_x, g_y_x_0_0_y_yz_y_y, g_y_x_0_0_y_yz_y_z, g_yy_xyz_y_x, g_yy_xyz_y_y, g_yy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_yy_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_yy_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_yy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_y_x_0_0_y_yz_z_x, g_y_x_0_0_y_yz_z_y, g_y_x_0_0_y_yz_z_z, g_yy_xyz_z_x, g_yy_xyz_z_y, g_yy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_yy_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_yy_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_yy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_y_x_0_0_y_zz_x_x, g_y_x_0_0_y_zz_x_y, g_y_x_0_0_y_zz_x_z, g_yy_xzz_x_x, g_yy_xzz_x_y, g_yy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_zz_x_x[i] = -2.0 * g_0_xzz_x_x[i] * b_exp + 4.0 * g_yy_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_x_y[i] = -2.0 * g_0_xzz_x_y[i] * b_exp + 4.0 * g_yy_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_x_z[i] = -2.0 * g_0_xzz_x_z[i] * b_exp + 4.0 * g_yy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_y_x_0_0_y_zz_y_x, g_y_x_0_0_y_zz_y_y, g_y_x_0_0_y_zz_y_z, g_yy_xzz_y_x, g_yy_xzz_y_y, g_yy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_zz_y_x[i] = -2.0 * g_0_xzz_y_x[i] * b_exp + 4.0 * g_yy_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_y_y[i] = -2.0 * g_0_xzz_y_y[i] * b_exp + 4.0 * g_yy_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_y_z[i] = -2.0 * g_0_xzz_y_z[i] * b_exp + 4.0 * g_yy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_y_x_0_0_y_zz_z_x, g_y_x_0_0_y_zz_z_y, g_y_x_0_0_y_zz_z_z, g_yy_xzz_z_x, g_yy_xzz_z_y, g_yy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_zz_z_x[i] = -2.0 * g_0_xzz_z_x[i] * b_exp + 4.0 * g_yy_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_z_y[i] = -2.0 * g_0_xzz_z_y[i] * b_exp + 4.0 * g_yy_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_z_z[i] = -2.0 * g_0_xzz_z_z[i] * b_exp + 4.0 * g_yy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_x_x, g_y_x_0_0_z_xx_x_y, g_y_x_0_0_z_xx_x_z, g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xxx_x_x, g_yz_xxx_x_y, g_yz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_x_x[i] = -4.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xxx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_x_y[i] = -4.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xxx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_x_z[i] = -4.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_y_x, g_y_x_0_0_z_xx_y_y, g_y_x_0_0_z_xx_y_z, g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xxx_y_x, g_yz_xxx_y_y, g_yz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_y_x[i] = -4.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xxx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_y_y[i] = -4.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xxx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_y_z[i] = -4.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_z_x, g_y_x_0_0_z_xx_z_y, g_y_x_0_0_z_xx_z_z, g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xxx_z_x, g_yz_xxx_z_y, g_yz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_z_x[i] = -4.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xxx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_z_y[i] = -4.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xxx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_z_z[i] = -4.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_y_x_0_0_z_xy_x_x, g_y_x_0_0_z_xy_x_y, g_y_x_0_0_z_xy_x_z, g_yz_xxy_x_x, g_yz_xxy_x_y, g_yz_xxy_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xy_x_x[i] = -2.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_xxy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_x_y[i] = -2.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_xxy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_x_z[i] = -2.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_y_x_0_0_z_xy_y_x, g_y_x_0_0_z_xy_y_y, g_y_x_0_0_z_xy_y_z, g_yz_xxy_y_x, g_yz_xxy_y_y, g_yz_xxy_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xy_y_x[i] = -2.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_xxy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_y_y[i] = -2.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_xxy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_y_z[i] = -2.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_y_x_0_0_z_xy_z_x, g_y_x_0_0_z_xy_z_y, g_y_x_0_0_z_xy_z_z, g_yz_xxy_z_x, g_yz_xxy_z_y, g_yz_xxy_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xy_z_x[i] = -2.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_xxy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_z_y[i] = -2.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_xxy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_z_z[i] = -2.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_y_x_0_0_z_xz_x_x, g_y_x_0_0_z_xz_x_y, g_y_x_0_0_z_xz_x_z, g_yz_xxz_x_x, g_yz_xxz_x_y, g_yz_xxz_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xz_x_x[i] = -2.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_xxz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_x_y[i] = -2.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_xxz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_x_z[i] = -2.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_y_x_0_0_z_xz_y_x, g_y_x_0_0_z_xz_y_y, g_y_x_0_0_z_xz_y_z, g_yz_xxz_y_x, g_yz_xxz_y_y, g_yz_xxz_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xz_y_x[i] = -2.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_xxz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_y_y[i] = -2.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_xxz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_y_z[i] = -2.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_y_x_0_0_z_xz_z_x, g_y_x_0_0_z_xz_z_y, g_y_x_0_0_z_xz_z_z, g_yz_xxz_z_x, g_yz_xxz_z_y, g_yz_xxz_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xz_z_x[i] = -2.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_xxz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_z_y[i] = -2.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_xxz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_z_z[i] = -2.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_y_x_0_0_z_yy_x_x, g_y_x_0_0_z_yy_x_y, g_y_x_0_0_z_yy_x_z, g_yz_xyy_x_x, g_yz_xyy_x_y, g_yz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yy_x_x[i] = 4.0 * g_yz_xyy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_x_y[i] = 4.0 * g_yz_xyy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_x_z[i] = 4.0 * g_yz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_y_x_0_0_z_yy_y_x, g_y_x_0_0_z_yy_y_y, g_y_x_0_0_z_yy_y_z, g_yz_xyy_y_x, g_yz_xyy_y_y, g_yz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yy_y_x[i] = 4.0 * g_yz_xyy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_y_y[i] = 4.0 * g_yz_xyy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_y_z[i] = 4.0 * g_yz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_y_x_0_0_z_yy_z_x, g_y_x_0_0_z_yy_z_y, g_y_x_0_0_z_yy_z_z, g_yz_xyy_z_x, g_yz_xyy_z_y, g_yz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yy_z_x[i] = 4.0 * g_yz_xyy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_z_y[i] = 4.0 * g_yz_xyy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_z_z[i] = 4.0 * g_yz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_y_x_0_0_z_yz_x_x, g_y_x_0_0_z_yz_x_y, g_y_x_0_0_z_yz_x_z, g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yz_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_y_x_0_0_z_yz_y_x, g_y_x_0_0_z_yz_y_y, g_y_x_0_0_z_yz_y_z, g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yz_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_y_x_0_0_z_yz_z_x, g_y_x_0_0_z_yz_z_y, g_y_x_0_0_z_yz_z_z, g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yz_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_y_x_0_0_z_zz_x_x, g_y_x_0_0_z_zz_x_y, g_y_x_0_0_z_zz_x_z, g_yz_xzz_x_x, g_yz_xzz_x_y, g_yz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_zz_x_x[i] = 4.0 * g_yz_xzz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_x_y[i] = 4.0 * g_yz_xzz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_x_z[i] = 4.0 * g_yz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_y_x_0_0_z_zz_y_x, g_y_x_0_0_z_zz_y_y, g_y_x_0_0_z_zz_y_z, g_yz_xzz_y_x, g_yz_xzz_y_y, g_yz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_zz_y_x[i] = 4.0 * g_yz_xzz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_y_y[i] = 4.0 * g_yz_xzz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_y_z[i] = 4.0 * g_yz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_y_x_0_0_z_zz_z_x, g_y_x_0_0_z_zz_z_y, g_y_x_0_0_z_zz_z_z, g_yz_xzz_z_x, g_yz_xzz_z_y, g_yz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_zz_z_x[i] = 4.0 * g_yz_xzz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_z_y[i] = 4.0 * g_yz_xzz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_z_z[i] = 4.0 * g_yz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xy_xxy_x_x, g_xy_xxy_x_y, g_xy_xxy_x_z, g_y_y_0_0_x_xx_x_x, g_y_y_0_0_x_xx_x_y, g_y_y_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_x_x[i] = 4.0 * g_xy_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_x_y[i] = 4.0 * g_xy_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_x_z[i] = 4.0 * g_xy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xy_xxy_y_x, g_xy_xxy_y_y, g_xy_xxy_y_z, g_y_y_0_0_x_xx_y_x, g_y_y_0_0_x_xx_y_y, g_y_y_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_y_x[i] = 4.0 * g_xy_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_y_y[i] = 4.0 * g_xy_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_y_z[i] = 4.0 * g_xy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xy_xxy_z_x, g_xy_xxy_z_y, g_xy_xxy_z_z, g_y_y_0_0_x_xx_z_x, g_y_y_0_0_x_xx_z_y, g_y_y_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_z_x[i] = 4.0 * g_xy_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_z_y[i] = 4.0 * g_xy_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_z_z[i] = 4.0 * g_xy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xyy_x_x, g_xy_xyy_x_y, g_xy_xyy_x_z, g_y_y_0_0_x_xy_x_x, g_y_y_0_0_x_xy_x_y, g_y_y_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xy_x_x[i] = -2.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_x_y[i] = -2.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_x_z[i] = -2.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xyy_y_x, g_xy_xyy_y_y, g_xy_xyy_y_z, g_y_y_0_0_x_xy_y_x, g_y_y_0_0_x_xy_y_y, g_y_y_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xy_y_x[i] = -2.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_y_y[i] = -2.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_y_z[i] = -2.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xyy_z_x, g_xy_xyy_z_y, g_xy_xyy_z_z, g_y_y_0_0_x_xy_z_x, g_y_y_0_0_x_xy_z_y, g_y_y_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xy_z_x[i] = -2.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_z_y[i] = -2.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_z_z[i] = -2.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z, g_y_y_0_0_x_xz_x_x, g_y_y_0_0_x_xz_x_y, g_y_y_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xz_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z, g_y_y_0_0_x_xz_y_x, g_y_y_0_0_x_xz_y_y, g_y_y_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xz_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z, g_y_y_0_0_x_xz_z_x, g_y_y_0_0_x_xz_z_y, g_y_y_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xz_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_xy_yyy_x_x, g_xy_yyy_x_y, g_xy_yyy_x_z, g_y_y_0_0_x_yy_x_x, g_y_y_0_0_x_yy_x_y, g_y_y_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yy_x_x[i] = -4.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_x_y[i] = -4.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_x_z[i] = -4.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_xy_yyy_y_x, g_xy_yyy_y_y, g_xy_yyy_y_z, g_y_y_0_0_x_yy_y_x, g_y_y_0_0_x_yy_y_y, g_y_y_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yy_y_x[i] = -4.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_y_y[i] = -4.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_y_z[i] = -4.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_xy_yyy_z_x, g_xy_yyy_z_y, g_xy_yyy_z_z, g_y_y_0_0_x_yy_z_x, g_y_y_0_0_x_yy_z_y, g_y_y_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yy_z_x[i] = -4.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_z_y[i] = -4.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_z_z[i] = -4.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_xy_yyz_x_x, g_xy_yyz_x_y, g_xy_yyz_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z, g_y_y_0_0_x_yz_x_x, g_y_y_0_0_x_yz_x_y, g_y_y_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yz_x_x[i] = -2.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_x_y[i] = -2.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_x_z[i] = -2.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_xy_yyz_y_x, g_xy_yyz_y_y, g_xy_yyz_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z, g_y_y_0_0_x_yz_y_x, g_y_y_0_0_x_yz_y_y, g_y_y_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yz_y_x[i] = -2.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_y_y[i] = -2.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_y_z[i] = -2.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_xy_yyz_z_x, g_xy_yyz_z_y, g_xy_yyz_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z, g_y_y_0_0_x_yz_z_x, g_y_y_0_0_x_yz_z_y, g_y_y_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yz_z_x[i] = -2.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_z_y[i] = -2.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_z_z[i] = -2.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_xy_yzz_x_x, g_xy_yzz_x_y, g_xy_yzz_x_z, g_y_y_0_0_x_zz_x_x, g_y_y_0_0_x_zz_x_y, g_y_y_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_zz_x_x[i] = 4.0 * g_xy_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_x_y[i] = 4.0 * g_xy_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_x_z[i] = 4.0 * g_xy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_xy_yzz_y_x, g_xy_yzz_y_y, g_xy_yzz_y_z, g_y_y_0_0_x_zz_y_x, g_y_y_0_0_x_zz_y_y, g_y_y_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_zz_y_x[i] = 4.0 * g_xy_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_y_y[i] = 4.0 * g_xy_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_y_z[i] = 4.0 * g_xy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_xy_yzz_z_x, g_xy_yzz_z_y, g_xy_yzz_z_z, g_y_y_0_0_x_zz_z_x, g_y_y_0_0_x_zz_z_y, g_y_y_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_zz_z_x[i] = 4.0 * g_xy_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_z_y[i] = 4.0 * g_xy_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_z_z[i] = 4.0 * g_xy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_y_y_0_0_y_xx_x_x, g_y_y_0_0_y_xx_x_y, g_y_y_0_0_y_xx_x_z, g_yy_xxy_x_x, g_yy_xxy_x_y, g_yy_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_x_x[i] = -2.0 * g_0_xxy_x_x[i] * b_exp + 4.0 * g_yy_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_x_y[i] = -2.0 * g_0_xxy_x_y[i] * b_exp + 4.0 * g_yy_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_x_z[i] = -2.0 * g_0_xxy_x_z[i] * b_exp + 4.0 * g_yy_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_y_y_0_0_y_xx_y_x, g_y_y_0_0_y_xx_y_y, g_y_y_0_0_y_xx_y_z, g_yy_xxy_y_x, g_yy_xxy_y_y, g_yy_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_y_x[i] = -2.0 * g_0_xxy_y_x[i] * b_exp + 4.0 * g_yy_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_y_y[i] = -2.0 * g_0_xxy_y_y[i] * b_exp + 4.0 * g_yy_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_y_z[i] = -2.0 * g_0_xxy_y_z[i] * b_exp + 4.0 * g_yy_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_y_y_0_0_y_xx_z_x, g_y_y_0_0_y_xx_z_y, g_y_y_0_0_y_xx_z_z, g_yy_xxy_z_x, g_yy_xxy_z_y, g_yy_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_z_x[i] = -2.0 * g_0_xxy_z_x[i] * b_exp + 4.0 * g_yy_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_z_y[i] = -2.0 * g_0_xxy_z_y[i] * b_exp + 4.0 * g_yy_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_z_z[i] = -2.0 * g_0_xxy_z_z[i] * b_exp + 4.0 * g_yy_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_y_y_0_0_y_xy_x_x, g_y_y_0_0_y_xy_x_y, g_y_y_0_0_y_xy_x_z, g_yy_x_x_x, g_yy_x_x_y, g_yy_x_x_z, g_yy_xyy_x_x, g_yy_xyy_x_y, g_yy_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xy_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xyy_x_x[i] * b_exp - 2.0 * g_yy_x_x_x[i] * a_exp + 4.0 * g_yy_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xyy_x_y[i] * b_exp - 2.0 * g_yy_x_x_y[i] * a_exp + 4.0 * g_yy_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xyy_x_z[i] * b_exp - 2.0 * g_yy_x_x_z[i] * a_exp + 4.0 * g_yy_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_y_y_0_0_y_xy_y_x, g_y_y_0_0_y_xy_y_y, g_y_y_0_0_y_xy_y_z, g_yy_x_y_x, g_yy_x_y_y, g_yy_x_y_z, g_yy_xyy_y_x, g_yy_xyy_y_y, g_yy_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xy_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xyy_y_x[i] * b_exp - 2.0 * g_yy_x_y_x[i] * a_exp + 4.0 * g_yy_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xyy_y_y[i] * b_exp - 2.0 * g_yy_x_y_y[i] * a_exp + 4.0 * g_yy_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xyy_y_z[i] * b_exp - 2.0 * g_yy_x_y_z[i] * a_exp + 4.0 * g_yy_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_y_y_0_0_y_xy_z_x, g_y_y_0_0_y_xy_z_y, g_y_y_0_0_y_xy_z_z, g_yy_x_z_x, g_yy_x_z_y, g_yy_x_z_z, g_yy_xyy_z_x, g_yy_xyy_z_y, g_yy_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xy_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xyy_z_x[i] * b_exp - 2.0 * g_yy_x_z_x[i] * a_exp + 4.0 * g_yy_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xyy_z_y[i] * b_exp - 2.0 * g_yy_x_z_y[i] * a_exp + 4.0 * g_yy_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xyy_z_z[i] * b_exp - 2.0 * g_yy_x_z_z[i] * a_exp + 4.0 * g_yy_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_y_y_0_0_y_xz_x_x, g_y_y_0_0_y_xz_x_y, g_y_y_0_0_y_xz_x_z, g_yy_xyz_x_x, g_yy_xyz_x_y, g_yy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_yy_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_yy_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_yy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_y_y_0_0_y_xz_y_x, g_y_y_0_0_y_xz_y_y, g_y_y_0_0_y_xz_y_z, g_yy_xyz_y_x, g_yy_xyz_y_y, g_yy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_yy_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_yy_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_yy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_y_y_0_0_y_xz_z_x, g_y_y_0_0_y_xz_z_y, g_y_y_0_0_y_xz_z_z, g_yy_xyz_z_x, g_yy_xyz_z_y, g_yy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_yy_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_yy_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_yy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yyy_x_x, g_0_yyy_x_y, g_0_yyy_x_z, g_y_y_0_0_y_yy_x_x, g_y_y_0_0_y_yy_x_y, g_y_y_0_0_y_yy_x_z, g_yy_y_x_x, g_yy_y_x_y, g_yy_y_x_z, g_yy_yyy_x_x, g_yy_yyy_x_y, g_yy_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yy_x_x[i] = 2.0 * g_0_y_x_x[i] - 2.0 * g_0_yyy_x_x[i] * b_exp - 4.0 * g_yy_y_x_x[i] * a_exp + 4.0 * g_yy_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_x_y[i] = 2.0 * g_0_y_x_y[i] - 2.0 * g_0_yyy_x_y[i] * b_exp - 4.0 * g_yy_y_x_y[i] * a_exp + 4.0 * g_yy_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_x_z[i] = 2.0 * g_0_y_x_z[i] - 2.0 * g_0_yyy_x_z[i] * b_exp - 4.0 * g_yy_y_x_z[i] * a_exp + 4.0 * g_yy_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yyy_y_x, g_0_yyy_y_y, g_0_yyy_y_z, g_y_y_0_0_y_yy_y_x, g_y_y_0_0_y_yy_y_y, g_y_y_0_0_y_yy_y_z, g_yy_y_y_x, g_yy_y_y_y, g_yy_y_y_z, g_yy_yyy_y_x, g_yy_yyy_y_y, g_yy_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yy_y_x[i] = 2.0 * g_0_y_y_x[i] - 2.0 * g_0_yyy_y_x[i] * b_exp - 4.0 * g_yy_y_y_x[i] * a_exp + 4.0 * g_yy_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_y_y[i] = 2.0 * g_0_y_y_y[i] - 2.0 * g_0_yyy_y_y[i] * b_exp - 4.0 * g_yy_y_y_y[i] * a_exp + 4.0 * g_yy_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_y_z[i] = 2.0 * g_0_y_y_z[i] - 2.0 * g_0_yyy_y_z[i] * b_exp - 4.0 * g_yy_y_y_z[i] * a_exp + 4.0 * g_yy_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yyy_z_x, g_0_yyy_z_y, g_0_yyy_z_z, g_y_y_0_0_y_yy_z_x, g_y_y_0_0_y_yy_z_y, g_y_y_0_0_y_yy_z_z, g_yy_y_z_x, g_yy_y_z_y, g_yy_y_z_z, g_yy_yyy_z_x, g_yy_yyy_z_y, g_yy_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yy_z_x[i] = 2.0 * g_0_y_z_x[i] - 2.0 * g_0_yyy_z_x[i] * b_exp - 4.0 * g_yy_y_z_x[i] * a_exp + 4.0 * g_yy_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_z_y[i] = 2.0 * g_0_y_z_y[i] - 2.0 * g_0_yyy_z_y[i] * b_exp - 4.0 * g_yy_y_z_y[i] * a_exp + 4.0 * g_yy_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_z_z[i] = 2.0 * g_0_y_z_z[i] - 2.0 * g_0_yyy_z_z[i] * b_exp - 4.0 * g_yy_y_z_z[i] * a_exp + 4.0 * g_yy_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_y_y_0_0_y_yz_x_x, g_y_y_0_0_y_yz_x_y, g_y_y_0_0_y_yz_x_z, g_yy_yyz_x_x, g_yy_yyz_x_y, g_yy_yyz_x_z, g_yy_z_x_x, g_yy_z_x_y, g_yy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_yyz_x_x[i] * b_exp - 2.0 * g_yy_z_x_x[i] * a_exp + 4.0 * g_yy_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_yyz_x_y[i] * b_exp - 2.0 * g_yy_z_x_y[i] * a_exp + 4.0 * g_yy_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_yyz_x_z[i] * b_exp - 2.0 * g_yy_z_x_z[i] * a_exp + 4.0 * g_yy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_y_y_0_0_y_yz_y_x, g_y_y_0_0_y_yz_y_y, g_y_y_0_0_y_yz_y_z, g_yy_yyz_y_x, g_yy_yyz_y_y, g_yy_yyz_y_z, g_yy_z_y_x, g_yy_z_y_y, g_yy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_yyz_y_x[i] * b_exp - 2.0 * g_yy_z_y_x[i] * a_exp + 4.0 * g_yy_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_yyz_y_y[i] * b_exp - 2.0 * g_yy_z_y_y[i] * a_exp + 4.0 * g_yy_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_yyz_y_z[i] * b_exp - 2.0 * g_yy_z_y_z[i] * a_exp + 4.0 * g_yy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_y_y_0_0_y_yz_z_x, g_y_y_0_0_y_yz_z_y, g_y_y_0_0_y_yz_z_z, g_yy_yyz_z_x, g_yy_yyz_z_y, g_yy_yyz_z_z, g_yy_z_z_x, g_yy_z_z_y, g_yy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_yyz_z_x[i] * b_exp - 2.0 * g_yy_z_z_x[i] * a_exp + 4.0 * g_yy_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_yyz_z_y[i] * b_exp - 2.0 * g_yy_z_z_y[i] * a_exp + 4.0 * g_yy_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_yyz_z_z[i] * b_exp - 2.0 * g_yy_z_z_z[i] * a_exp + 4.0 * g_yy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_y_y_0_0_y_zz_x_x, g_y_y_0_0_y_zz_x_y, g_y_y_0_0_y_zz_x_z, g_yy_yzz_x_x, g_yy_yzz_x_y, g_yy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_zz_x_x[i] = -2.0 * g_0_yzz_x_x[i] * b_exp + 4.0 * g_yy_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_x_y[i] = -2.0 * g_0_yzz_x_y[i] * b_exp + 4.0 * g_yy_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_x_z[i] = -2.0 * g_0_yzz_x_z[i] * b_exp + 4.0 * g_yy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_y_y_0_0_y_zz_y_x, g_y_y_0_0_y_zz_y_y, g_y_y_0_0_y_zz_y_z, g_yy_yzz_y_x, g_yy_yzz_y_y, g_yy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_zz_y_x[i] = -2.0 * g_0_yzz_y_x[i] * b_exp + 4.0 * g_yy_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_y_y[i] = -2.0 * g_0_yzz_y_y[i] * b_exp + 4.0 * g_yy_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_y_z[i] = -2.0 * g_0_yzz_y_z[i] * b_exp + 4.0 * g_yy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_y_y_0_0_y_zz_z_x, g_y_y_0_0_y_zz_z_y, g_y_y_0_0_y_zz_z_z, g_yy_yzz_z_x, g_yy_yzz_z_y, g_yy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_zz_z_x[i] = -2.0 * g_0_yzz_z_x[i] * b_exp + 4.0 * g_yy_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_z_y[i] = -2.0 * g_0_yzz_z_y[i] * b_exp + 4.0 * g_yy_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_z_z[i] = -2.0 * g_0_yzz_z_z[i] * b_exp + 4.0 * g_yy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_x_x, g_y_y_0_0_z_xx_x_y, g_y_y_0_0_z_xx_x_z, g_yz_xxy_x_x, g_yz_xxy_x_y, g_yz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_x_x[i] = 4.0 * g_yz_xxy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_x_y[i] = 4.0 * g_yz_xxy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_x_z[i] = 4.0 * g_yz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_y_x, g_y_y_0_0_z_xx_y_y, g_y_y_0_0_z_xx_y_z, g_yz_xxy_y_x, g_yz_xxy_y_y, g_yz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_y_x[i] = 4.0 * g_yz_xxy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_y_y[i] = 4.0 * g_yz_xxy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_y_z[i] = 4.0 * g_yz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_z_x, g_y_y_0_0_z_xx_z_y, g_y_y_0_0_z_xx_z_z, g_yz_xxy_z_x, g_yz_xxy_z_y, g_yz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_z_x[i] = 4.0 * g_yz_xxy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_z_y[i] = 4.0 * g_yz_xxy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_z_z[i] = 4.0 * g_yz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_y_y_0_0_z_xy_x_x, g_y_y_0_0_z_xy_x_y, g_y_y_0_0_z_xy_x_z, g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xyy_x_x, g_yz_xyy_x_y, g_yz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xy_x_x[i] = -2.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_x_y[i] = -2.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_x_z[i] = -2.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_y_y_0_0_z_xy_y_x, g_y_y_0_0_z_xy_y_y, g_y_y_0_0_z_xy_y_z, g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xyy_y_x, g_yz_xyy_y_y, g_yz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xy_y_x[i] = -2.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_y_y[i] = -2.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_y_z[i] = -2.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_y_y_0_0_z_xy_z_x, g_y_y_0_0_z_xy_z_y, g_y_y_0_0_z_xy_z_z, g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xyy_z_x, g_yz_xyy_z_y, g_yz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xy_z_x[i] = -2.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_z_y[i] = -2.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_z_z[i] = -2.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_y_y_0_0_z_xz_x_x, g_y_y_0_0_z_xz_x_y, g_y_y_0_0_z_xz_x_z, g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xz_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_y_y_0_0_z_xz_y_x, g_y_y_0_0_z_xz_y_y, g_y_y_0_0_z_xz_y_z, g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xz_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_y_y_0_0_z_xz_z_x, g_y_y_0_0_z_xz_z_y, g_y_y_0_0_z_xz_z_z, g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xz_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_y_y_0_0_z_yy_x_x, g_y_y_0_0_z_yy_x_y, g_y_y_0_0_z_yy_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_yz_yyy_x_x, g_yz_yyy_x_y, g_yz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yy_x_x[i] = -4.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_yyy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_x_y[i] = -4.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_yyy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_x_z[i] = -4.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_y_y_0_0_z_yy_y_x, g_y_y_0_0_z_yy_y_y, g_y_y_0_0_z_yy_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_yz_yyy_y_x, g_yz_yyy_y_y, g_yz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yy_y_x[i] = -4.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_yyy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_y_y[i] = -4.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_yyy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_y_z[i] = -4.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_y_y_0_0_z_yy_z_x, g_y_y_0_0_z_yy_z_y, g_y_y_0_0_z_yy_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_yz_yyy_z_x, g_yz_yyy_z_y, g_yz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yy_z_x[i] = -4.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_yyy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_z_y[i] = -4.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_yyy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_z_z[i] = -4.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_y_y_0_0_z_yz_x_x, g_y_y_0_0_z_yz_x_y, g_y_y_0_0_z_yz_x_z, g_yz_yyz_x_x, g_yz_yyz_x_y, g_yz_yyz_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yz_x_x[i] = -2.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_yyz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_x_y[i] = -2.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_yyz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_x_z[i] = -2.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_y_y_0_0_z_yz_y_x, g_y_y_0_0_z_yz_y_y, g_y_y_0_0_z_yz_y_z, g_yz_yyz_y_x, g_yz_yyz_y_y, g_yz_yyz_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yz_y_x[i] = -2.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_yyz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_y_y[i] = -2.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_yyz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_y_z[i] = -2.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_y_y_0_0_z_yz_z_x, g_y_y_0_0_z_yz_z_y, g_y_y_0_0_z_yz_z_z, g_yz_yyz_z_x, g_yz_yyz_z_y, g_yz_yyz_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yz_z_x[i] = -2.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_yyz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_z_y[i] = -2.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_yyz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_z_z[i] = -2.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_y_y_0_0_z_zz_x_x, g_y_y_0_0_z_zz_x_y, g_y_y_0_0_z_zz_x_z, g_yz_yzz_x_x, g_yz_yzz_x_y, g_yz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_zz_x_x[i] = 4.0 * g_yz_yzz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_x_y[i] = 4.0 * g_yz_yzz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_x_z[i] = 4.0 * g_yz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_y_y_0_0_z_zz_y_x, g_y_y_0_0_z_zz_y_y, g_y_y_0_0_z_zz_y_z, g_yz_yzz_y_x, g_yz_yzz_y_y, g_yz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_zz_y_x[i] = 4.0 * g_yz_yzz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_y_y[i] = 4.0 * g_yz_yzz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_y_z[i] = 4.0 * g_yz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_y_y_0_0_z_zz_z_x, g_y_y_0_0_z_zz_z_y, g_y_y_0_0_z_zz_z_z, g_yz_yzz_z_x, g_yz_yzz_z_y, g_yz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_zz_z_x[i] = 4.0 * g_yz_yzz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_z_y[i] = 4.0 * g_yz_yzz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_z_z[i] = 4.0 * g_yz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_xy_xxz_x_x, g_xy_xxz_x_y, g_xy_xxz_x_z, g_y_z_0_0_x_xx_x_x, g_y_z_0_0_x_xx_x_y, g_y_z_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_x_x[i] = 4.0 * g_xy_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_x_y[i] = 4.0 * g_xy_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_x_z[i] = 4.0 * g_xy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_xy_xxz_y_x, g_xy_xxz_y_y, g_xy_xxz_y_z, g_y_z_0_0_x_xx_y_x, g_y_z_0_0_x_xx_y_y, g_y_z_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_y_x[i] = 4.0 * g_xy_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_y_y[i] = 4.0 * g_xy_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_y_z[i] = 4.0 * g_xy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_xy_xxz_z_x, g_xy_xxz_z_y, g_xy_xxz_z_z, g_y_z_0_0_x_xx_z_x, g_y_z_0_0_x_xx_z_y, g_y_z_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_z_x[i] = 4.0 * g_xy_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_z_y[i] = 4.0 * g_xy_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_z_z[i] = 4.0 * g_xy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_xy_xyz_x_x, g_xy_xyz_x_y, g_xy_xyz_x_z, g_y_z_0_0_x_xy_x_x, g_y_z_0_0_x_xy_x_y, g_y_z_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xy_x_x[i] = 4.0 * g_xy_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_x_y[i] = 4.0 * g_xy_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_x_z[i] = 4.0 * g_xy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_xy_xyz_y_x, g_xy_xyz_y_y, g_xy_xyz_y_z, g_y_z_0_0_x_xy_y_x, g_y_z_0_0_x_xy_y_y, g_y_z_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xy_y_x[i] = 4.0 * g_xy_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_y_y[i] = 4.0 * g_xy_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_y_z[i] = 4.0 * g_xy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_xy_xyz_z_x, g_xy_xyz_z_y, g_xy_xyz_z_z, g_y_z_0_0_x_xy_z_x, g_y_z_0_0_x_xy_z_y, g_y_z_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xy_z_x[i] = 4.0 * g_xy_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_z_y[i] = 4.0 * g_xy_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_z_z[i] = 4.0 * g_xy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_xy_xzz_x_x, g_xy_xzz_x_y, g_xy_xzz_x_z, g_y_z_0_0_x_xz_x_x, g_y_z_0_0_x_xz_x_y, g_y_z_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xz_x_x[i] = -2.0 * g_xy_x_x_x[i] * a_exp + 4.0 * g_xy_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_x_y[i] = -2.0 * g_xy_x_x_y[i] * a_exp + 4.0 * g_xy_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_x_z[i] = -2.0 * g_xy_x_x_z[i] * a_exp + 4.0 * g_xy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_xy_xzz_y_x, g_xy_xzz_y_y, g_xy_xzz_y_z, g_y_z_0_0_x_xz_y_x, g_y_z_0_0_x_xz_y_y, g_y_z_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xz_y_x[i] = -2.0 * g_xy_x_y_x[i] * a_exp + 4.0 * g_xy_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_y_y[i] = -2.0 * g_xy_x_y_y[i] * a_exp + 4.0 * g_xy_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_y_z[i] = -2.0 * g_xy_x_y_z[i] * a_exp + 4.0 * g_xy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_xy_xzz_z_x, g_xy_xzz_z_y, g_xy_xzz_z_z, g_y_z_0_0_x_xz_z_x, g_y_z_0_0_x_xz_z_y, g_y_z_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xz_z_x[i] = -2.0 * g_xy_x_z_x[i] * a_exp + 4.0 * g_xy_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_z_y[i] = -2.0 * g_xy_x_z_y[i] * a_exp + 4.0 * g_xy_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_z_z[i] = -2.0 * g_xy_x_z_z[i] * a_exp + 4.0 * g_xy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_xy_yyz_x_x, g_xy_yyz_x_y, g_xy_yyz_x_z, g_y_z_0_0_x_yy_x_x, g_y_z_0_0_x_yy_x_y, g_y_z_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yy_x_x[i] = 4.0 * g_xy_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_x_y[i] = 4.0 * g_xy_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_x_z[i] = 4.0 * g_xy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_xy_yyz_y_x, g_xy_yyz_y_y, g_xy_yyz_y_z, g_y_z_0_0_x_yy_y_x, g_y_z_0_0_x_yy_y_y, g_y_z_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yy_y_x[i] = 4.0 * g_xy_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_y_y[i] = 4.0 * g_xy_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_y_z[i] = 4.0 * g_xy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_xy_yyz_z_x, g_xy_yyz_z_y, g_xy_yyz_z_z, g_y_z_0_0_x_yy_z_x, g_y_z_0_0_x_yy_z_y, g_y_z_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yy_z_x[i] = 4.0 * g_xy_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_z_y[i] = 4.0 * g_xy_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_z_z[i] = 4.0 * g_xy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_xy_yzz_x_x, g_xy_yzz_x_y, g_xy_yzz_x_z, g_y_z_0_0_x_yz_x_x, g_y_z_0_0_x_yz_x_y, g_y_z_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yz_x_x[i] = -2.0 * g_xy_y_x_x[i] * a_exp + 4.0 * g_xy_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_x_y[i] = -2.0 * g_xy_y_x_y[i] * a_exp + 4.0 * g_xy_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_x_z[i] = -2.0 * g_xy_y_x_z[i] * a_exp + 4.0 * g_xy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_xy_yzz_y_x, g_xy_yzz_y_y, g_xy_yzz_y_z, g_y_z_0_0_x_yz_y_x, g_y_z_0_0_x_yz_y_y, g_y_z_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yz_y_x[i] = -2.0 * g_xy_y_y_x[i] * a_exp + 4.0 * g_xy_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_y_y[i] = -2.0 * g_xy_y_y_y[i] * a_exp + 4.0 * g_xy_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_y_z[i] = -2.0 * g_xy_y_y_z[i] * a_exp + 4.0 * g_xy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_xy_yzz_z_x, g_xy_yzz_z_y, g_xy_yzz_z_z, g_y_z_0_0_x_yz_z_x, g_y_z_0_0_x_yz_z_y, g_y_z_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yz_z_x[i] = -2.0 * g_xy_y_z_x[i] * a_exp + 4.0 * g_xy_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_z_y[i] = -2.0 * g_xy_y_z_y[i] * a_exp + 4.0 * g_xy_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_z_z[i] = -2.0 * g_xy_y_z_z[i] * a_exp + 4.0 * g_xy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z, g_xy_zzz_x_x, g_xy_zzz_x_y, g_xy_zzz_x_z, g_y_z_0_0_x_zz_x_x, g_y_z_0_0_x_zz_x_y, g_y_z_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_zz_x_x[i] = -4.0 * g_xy_z_x_x[i] * a_exp + 4.0 * g_xy_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_x_y[i] = -4.0 * g_xy_z_x_y[i] * a_exp + 4.0 * g_xy_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_x_z[i] = -4.0 * g_xy_z_x_z[i] * a_exp + 4.0 * g_xy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z, g_xy_zzz_y_x, g_xy_zzz_y_y, g_xy_zzz_y_z, g_y_z_0_0_x_zz_y_x, g_y_z_0_0_x_zz_y_y, g_y_z_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_zz_y_x[i] = -4.0 * g_xy_z_y_x[i] * a_exp + 4.0 * g_xy_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_y_y[i] = -4.0 * g_xy_z_y_y[i] * a_exp + 4.0 * g_xy_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_y_z[i] = -4.0 * g_xy_z_y_z[i] * a_exp + 4.0 * g_xy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z, g_xy_zzz_z_x, g_xy_zzz_z_y, g_xy_zzz_z_z, g_y_z_0_0_x_zz_z_x, g_y_z_0_0_x_zz_z_y, g_y_z_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_zz_z_x[i] = -4.0 * g_xy_z_z_x[i] * a_exp + 4.0 * g_xy_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_z_y[i] = -4.0 * g_xy_z_z_y[i] * a_exp + 4.0 * g_xy_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_z_z[i] = -4.0 * g_xy_z_z_z[i] * a_exp + 4.0 * g_xy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_y_z_0_0_y_xx_x_x, g_y_z_0_0_y_xx_x_y, g_y_z_0_0_y_xx_x_z, g_yy_xxz_x_x, g_yy_xxz_x_y, g_yy_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_x_x[i] = -2.0 * g_0_xxz_x_x[i] * b_exp + 4.0 * g_yy_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_x_y[i] = -2.0 * g_0_xxz_x_y[i] * b_exp + 4.0 * g_yy_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_x_z[i] = -2.0 * g_0_xxz_x_z[i] * b_exp + 4.0 * g_yy_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_y_z_0_0_y_xx_y_x, g_y_z_0_0_y_xx_y_y, g_y_z_0_0_y_xx_y_z, g_yy_xxz_y_x, g_yy_xxz_y_y, g_yy_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_y_x[i] = -2.0 * g_0_xxz_y_x[i] * b_exp + 4.0 * g_yy_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_y_y[i] = -2.0 * g_0_xxz_y_y[i] * b_exp + 4.0 * g_yy_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_y_z[i] = -2.0 * g_0_xxz_y_z[i] * b_exp + 4.0 * g_yy_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_y_z_0_0_y_xx_z_x, g_y_z_0_0_y_xx_z_y, g_y_z_0_0_y_xx_z_z, g_yy_xxz_z_x, g_yy_xxz_z_y, g_yy_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_z_x[i] = -2.0 * g_0_xxz_z_x[i] * b_exp + 4.0 * g_yy_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_z_y[i] = -2.0 * g_0_xxz_z_y[i] * b_exp + 4.0 * g_yy_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_z_z[i] = -2.0 * g_0_xxz_z_z[i] * b_exp + 4.0 * g_yy_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_y_z_0_0_y_xy_x_x, g_y_z_0_0_y_xy_x_y, g_y_z_0_0_y_xy_x_z, g_yy_xyz_x_x, g_yy_xyz_x_y, g_yy_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xy_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_yy_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_yy_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_yy_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_y_z_0_0_y_xy_y_x, g_y_z_0_0_y_xy_y_y, g_y_z_0_0_y_xy_y_z, g_yy_xyz_y_x, g_yy_xyz_y_y, g_yy_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xy_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_yy_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_yy_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_yy_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_y_z_0_0_y_xy_z_x, g_y_z_0_0_y_xy_z_y, g_y_z_0_0_y_xy_z_z, g_yy_xyz_z_x, g_yy_xyz_z_y, g_yy_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xy_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_yy_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_yy_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_yy_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_y_z_0_0_y_xz_x_x, g_y_z_0_0_y_xz_x_y, g_y_z_0_0_y_xz_x_z, g_yy_x_x_x, g_yy_x_x_y, g_yy_x_x_z, g_yy_xzz_x_x, g_yy_xzz_x_y, g_yy_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xz_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xzz_x_x[i] * b_exp - 2.0 * g_yy_x_x_x[i] * a_exp + 4.0 * g_yy_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xzz_x_y[i] * b_exp - 2.0 * g_yy_x_x_y[i] * a_exp + 4.0 * g_yy_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xzz_x_z[i] * b_exp - 2.0 * g_yy_x_x_z[i] * a_exp + 4.0 * g_yy_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_y_z_0_0_y_xz_y_x, g_y_z_0_0_y_xz_y_y, g_y_z_0_0_y_xz_y_z, g_yy_x_y_x, g_yy_x_y_y, g_yy_x_y_z, g_yy_xzz_y_x, g_yy_xzz_y_y, g_yy_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xz_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xzz_y_x[i] * b_exp - 2.0 * g_yy_x_y_x[i] * a_exp + 4.0 * g_yy_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xzz_y_y[i] * b_exp - 2.0 * g_yy_x_y_y[i] * a_exp + 4.0 * g_yy_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xzz_y_z[i] * b_exp - 2.0 * g_yy_x_y_z[i] * a_exp + 4.0 * g_yy_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_y_z_0_0_y_xz_z_x, g_y_z_0_0_y_xz_z_y, g_y_z_0_0_y_xz_z_z, g_yy_x_z_x, g_yy_x_z_y, g_yy_x_z_z, g_yy_xzz_z_x, g_yy_xzz_z_y, g_yy_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xz_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xzz_z_x[i] * b_exp - 2.0 * g_yy_x_z_x[i] * a_exp + 4.0 * g_yy_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xzz_z_y[i] * b_exp - 2.0 * g_yy_x_z_y[i] * a_exp + 4.0 * g_yy_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xzz_z_z[i] * b_exp - 2.0 * g_yy_x_z_z[i] * a_exp + 4.0 * g_yy_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_y_z_0_0_y_yy_x_x, g_y_z_0_0_y_yy_x_y, g_y_z_0_0_y_yy_x_z, g_yy_yyz_x_x, g_yy_yyz_x_y, g_yy_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yy_x_x[i] = -2.0 * g_0_yyz_x_x[i] * b_exp + 4.0 * g_yy_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_x_y[i] = -2.0 * g_0_yyz_x_y[i] * b_exp + 4.0 * g_yy_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_x_z[i] = -2.0 * g_0_yyz_x_z[i] * b_exp + 4.0 * g_yy_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_y_z_0_0_y_yy_y_x, g_y_z_0_0_y_yy_y_y, g_y_z_0_0_y_yy_y_z, g_yy_yyz_y_x, g_yy_yyz_y_y, g_yy_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yy_y_x[i] = -2.0 * g_0_yyz_y_x[i] * b_exp + 4.0 * g_yy_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_y_y[i] = -2.0 * g_0_yyz_y_y[i] * b_exp + 4.0 * g_yy_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_y_z[i] = -2.0 * g_0_yyz_y_z[i] * b_exp + 4.0 * g_yy_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_y_z_0_0_y_yy_z_x, g_y_z_0_0_y_yy_z_y, g_y_z_0_0_y_yy_z_z, g_yy_yyz_z_x, g_yy_yyz_z_y, g_yy_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yy_z_x[i] = -2.0 * g_0_yyz_z_x[i] * b_exp + 4.0 * g_yy_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_z_y[i] = -2.0 * g_0_yyz_z_y[i] * b_exp + 4.0 * g_yy_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_z_z[i] = -2.0 * g_0_yyz_z_z[i] * b_exp + 4.0 * g_yy_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_y_z_0_0_y_yz_x_x, g_y_z_0_0_y_yz_x_y, g_y_z_0_0_y_yz_x_z, g_yy_y_x_x, g_yy_y_x_y, g_yy_y_x_z, g_yy_yzz_x_x, g_yy_yzz_x_y, g_yy_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yz_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_yzz_x_x[i] * b_exp - 2.0 * g_yy_y_x_x[i] * a_exp + 4.0 * g_yy_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_yzz_x_y[i] * b_exp - 2.0 * g_yy_y_x_y[i] * a_exp + 4.0 * g_yy_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_yzz_x_z[i] * b_exp - 2.0 * g_yy_y_x_z[i] * a_exp + 4.0 * g_yy_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_y_z_0_0_y_yz_y_x, g_y_z_0_0_y_yz_y_y, g_y_z_0_0_y_yz_y_z, g_yy_y_y_x, g_yy_y_y_y, g_yy_y_y_z, g_yy_yzz_y_x, g_yy_yzz_y_y, g_yy_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yz_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_yzz_y_x[i] * b_exp - 2.0 * g_yy_y_y_x[i] * a_exp + 4.0 * g_yy_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_yzz_y_y[i] * b_exp - 2.0 * g_yy_y_y_y[i] * a_exp + 4.0 * g_yy_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_yzz_y_z[i] * b_exp - 2.0 * g_yy_y_y_z[i] * a_exp + 4.0 * g_yy_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_y_z_0_0_y_yz_z_x, g_y_z_0_0_y_yz_z_y, g_y_z_0_0_y_yz_z_z, g_yy_y_z_x, g_yy_y_z_y, g_yy_y_z_z, g_yy_yzz_z_x, g_yy_yzz_z_y, g_yy_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yz_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_yzz_z_x[i] * b_exp - 2.0 * g_yy_y_z_x[i] * a_exp + 4.0 * g_yy_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_yzz_z_y[i] * b_exp - 2.0 * g_yy_y_z_y[i] * a_exp + 4.0 * g_yy_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_yzz_z_z[i] * b_exp - 2.0 * g_yy_y_z_z[i] * a_exp + 4.0 * g_yy_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_0_zzz_x_x, g_0_zzz_x_y, g_0_zzz_x_z, g_y_z_0_0_y_zz_x_x, g_y_z_0_0_y_zz_x_y, g_y_z_0_0_y_zz_x_z, g_yy_z_x_x, g_yy_z_x_y, g_yy_z_x_z, g_yy_zzz_x_x, g_yy_zzz_x_y, g_yy_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_zz_x_x[i] = 2.0 * g_0_z_x_x[i] - 2.0 * g_0_zzz_x_x[i] * b_exp - 4.0 * g_yy_z_x_x[i] * a_exp + 4.0 * g_yy_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_x_y[i] = 2.0 * g_0_z_x_y[i] - 2.0 * g_0_zzz_x_y[i] * b_exp - 4.0 * g_yy_z_x_y[i] * a_exp + 4.0 * g_yy_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_x_z[i] = 2.0 * g_0_z_x_z[i] - 2.0 * g_0_zzz_x_z[i] * b_exp - 4.0 * g_yy_z_x_z[i] * a_exp + 4.0 * g_yy_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_0_zzz_y_x, g_0_zzz_y_y, g_0_zzz_y_z, g_y_z_0_0_y_zz_y_x, g_y_z_0_0_y_zz_y_y, g_y_z_0_0_y_zz_y_z, g_yy_z_y_x, g_yy_z_y_y, g_yy_z_y_z, g_yy_zzz_y_x, g_yy_zzz_y_y, g_yy_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_zz_y_x[i] = 2.0 * g_0_z_y_x[i] - 2.0 * g_0_zzz_y_x[i] * b_exp - 4.0 * g_yy_z_y_x[i] * a_exp + 4.0 * g_yy_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_y_y[i] = 2.0 * g_0_z_y_y[i] - 2.0 * g_0_zzz_y_y[i] * b_exp - 4.0 * g_yy_z_y_y[i] * a_exp + 4.0 * g_yy_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_y_z[i] = 2.0 * g_0_z_y_z[i] - 2.0 * g_0_zzz_y_z[i] * b_exp - 4.0 * g_yy_z_y_z[i] * a_exp + 4.0 * g_yy_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_0_zzz_z_x, g_0_zzz_z_y, g_0_zzz_z_z, g_y_z_0_0_y_zz_z_x, g_y_z_0_0_y_zz_z_y, g_y_z_0_0_y_zz_z_z, g_yy_z_z_x, g_yy_z_z_y, g_yy_z_z_z, g_yy_zzz_z_x, g_yy_zzz_z_y, g_yy_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_zz_z_x[i] = 2.0 * g_0_z_z_x[i] - 2.0 * g_0_zzz_z_x[i] * b_exp - 4.0 * g_yy_z_z_x[i] * a_exp + 4.0 * g_yy_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_z_y[i] = 2.0 * g_0_z_z_y[i] - 2.0 * g_0_zzz_z_y[i] * b_exp - 4.0 * g_yy_z_z_y[i] * a_exp + 4.0 * g_yy_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_z_z[i] = 2.0 * g_0_z_z_z[i] - 2.0 * g_0_zzz_z_z[i] * b_exp - 4.0 * g_yy_z_z_z[i] * a_exp + 4.0 * g_yy_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_x_x, g_y_z_0_0_z_xx_x_y, g_y_z_0_0_z_xx_x_z, g_yz_xxz_x_x, g_yz_xxz_x_y, g_yz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_x_x[i] = 4.0 * g_yz_xxz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_x_y[i] = 4.0 * g_yz_xxz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_x_z[i] = 4.0 * g_yz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_y_x, g_y_z_0_0_z_xx_y_y, g_y_z_0_0_z_xx_y_z, g_yz_xxz_y_x, g_yz_xxz_y_y, g_yz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_y_x[i] = 4.0 * g_yz_xxz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_y_y[i] = 4.0 * g_yz_xxz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_y_z[i] = 4.0 * g_yz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_z_x, g_y_z_0_0_z_xx_z_y, g_y_z_0_0_z_xx_z_z, g_yz_xxz_z_x, g_yz_xxz_z_y, g_yz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_z_x[i] = 4.0 * g_yz_xxz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_z_y[i] = 4.0 * g_yz_xxz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_z_z[i] = 4.0 * g_yz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_y_z_0_0_z_xy_x_x, g_y_z_0_0_z_xy_x_y, g_y_z_0_0_z_xy_x_z, g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xy_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_y_z_0_0_z_xy_y_x, g_y_z_0_0_z_xy_y_y, g_y_z_0_0_z_xy_y_z, g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xy_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_y_z_0_0_z_xy_z_x, g_y_z_0_0_z_xy_z_y, g_y_z_0_0_z_xy_z_z, g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xy_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_y_z_0_0_z_xz_x_x, g_y_z_0_0_z_xz_x_y, g_y_z_0_0_z_xz_x_z, g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xzz_x_x, g_yz_xzz_x_y, g_yz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xz_x_x[i] = -2.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_x_y[i] = -2.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_x_z[i] = -2.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_y_z_0_0_z_xz_y_x, g_y_z_0_0_z_xz_y_y, g_y_z_0_0_z_xz_y_z, g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xzz_y_x, g_yz_xzz_y_y, g_yz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xz_y_x[i] = -2.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_y_y[i] = -2.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_y_z[i] = -2.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_y_z_0_0_z_xz_z_x, g_y_z_0_0_z_xz_z_y, g_y_z_0_0_z_xz_z_z, g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xzz_z_x, g_yz_xzz_z_y, g_yz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xz_z_x[i] = -2.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_z_y[i] = -2.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_z_z[i] = -2.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_y_z_0_0_z_yy_x_x, g_y_z_0_0_z_yy_x_y, g_y_z_0_0_z_yy_x_z, g_yz_yyz_x_x, g_yz_yyz_x_y, g_yz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yy_x_x[i] = 4.0 * g_yz_yyz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_x_y[i] = 4.0 * g_yz_yyz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_x_z[i] = 4.0 * g_yz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_y_z_0_0_z_yy_y_x, g_y_z_0_0_z_yy_y_y, g_y_z_0_0_z_yy_y_z, g_yz_yyz_y_x, g_yz_yyz_y_y, g_yz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yy_y_x[i] = 4.0 * g_yz_yyz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_y_y[i] = 4.0 * g_yz_yyz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_y_z[i] = 4.0 * g_yz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_y_z_0_0_z_yy_z_x, g_y_z_0_0_z_yy_z_y, g_y_z_0_0_z_yy_z_z, g_yz_yyz_z_x, g_yz_yyz_z_y, g_yz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yy_z_x[i] = 4.0 * g_yz_yyz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_z_y[i] = 4.0 * g_yz_yyz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_z_z[i] = 4.0 * g_yz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_y_z_0_0_z_yz_x_x, g_y_z_0_0_z_yz_x_y, g_y_z_0_0_z_yz_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_yz_yzz_x_x, g_yz_yzz_x_y, g_yz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yz_x_x[i] = -2.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_yzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_x_y[i] = -2.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_yzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_x_z[i] = -2.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_y_z_0_0_z_yz_y_x, g_y_z_0_0_z_yz_y_y, g_y_z_0_0_z_yz_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_yz_yzz_y_x, g_yz_yzz_y_y, g_yz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yz_y_x[i] = -2.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_yzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_y_y[i] = -2.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_yzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_y_z[i] = -2.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_y_z_0_0_z_yz_z_x, g_y_z_0_0_z_yz_z_y, g_y_z_0_0_z_yz_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_yz_yzz_z_x, g_yz_yzz_z_y, g_yz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yz_z_x[i] = -2.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_yzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_z_y[i] = -2.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_yzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_z_z[i] = -2.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_y_z_0_0_z_zz_x_x, g_y_z_0_0_z_zz_x_y, g_y_z_0_0_z_zz_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z, g_yz_zzz_x_x, g_yz_zzz_x_y, g_yz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_zz_x_x[i] = -4.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_zzz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_x_y[i] = -4.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_zzz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_x_z[i] = -4.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_y_z_0_0_z_zz_y_x, g_y_z_0_0_z_zz_y_y, g_y_z_0_0_z_zz_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z, g_yz_zzz_y_x, g_yz_zzz_y_y, g_yz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_zz_y_x[i] = -4.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_zzz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_y_y[i] = -4.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_zzz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_y_z[i] = -4.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_y_z_0_0_z_zz_z_x, g_y_z_0_0_z_zz_z_y, g_y_z_0_0_z_zz_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z, g_yz_zzz_z_x, g_yz_zzz_z_y, g_yz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_zz_z_x[i] = -4.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_zzz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_z_y[i] = -4.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_zzz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_z_z[i] = -4.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (972-975)

    #pragma omp simd aligned(g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xxx_x_x, g_xz_xxx_x_y, g_xz_xxx_x_z, g_z_x_0_0_x_xx_x_x, g_z_x_0_0_x_xx_x_y, g_z_x_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_x_x[i] = -4.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_x_y[i] = -4.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_x_z[i] = -4.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (975-978)

    #pragma omp simd aligned(g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xxx_y_x, g_xz_xxx_y_y, g_xz_xxx_y_z, g_z_x_0_0_x_xx_y_x, g_z_x_0_0_x_xx_y_y, g_z_x_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_y_x[i] = -4.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_y_y[i] = -4.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_y_z[i] = -4.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (978-981)

    #pragma omp simd aligned(g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xxx_z_x, g_xz_xxx_z_y, g_xz_xxx_z_z, g_z_x_0_0_x_xx_z_x, g_z_x_0_0_x_xx_z_y, g_z_x_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_z_x[i] = -4.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_z_y[i] = -4.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_z_z[i] = -4.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (981-984)

    #pragma omp simd aligned(g_xz_xxy_x_x, g_xz_xxy_x_y, g_xz_xxy_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_z_x_0_0_x_xy_x_x, g_z_x_0_0_x_xy_x_y, g_z_x_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xy_x_x[i] = -2.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_x_y[i] = -2.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_x_z[i] = -2.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (984-987)

    #pragma omp simd aligned(g_xz_xxy_y_x, g_xz_xxy_y_y, g_xz_xxy_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_z_x_0_0_x_xy_y_x, g_z_x_0_0_x_xy_y_y, g_z_x_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xy_y_x[i] = -2.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_y_y[i] = -2.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_y_z[i] = -2.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (987-990)

    #pragma omp simd aligned(g_xz_xxy_z_x, g_xz_xxy_z_y, g_xz_xxy_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_z_x_0_0_x_xy_z_x, g_z_x_0_0_x_xy_z_y, g_z_x_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xy_z_x[i] = -2.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_z_y[i] = -2.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_z_z[i] = -2.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (990-993)

    #pragma omp simd aligned(g_xz_xxz_x_x, g_xz_xxz_x_y, g_xz_xxz_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z, g_z_x_0_0_x_xz_x_x, g_z_x_0_0_x_xz_x_y, g_z_x_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xz_x_x[i] = -2.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_x_y[i] = -2.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_x_z[i] = -2.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (993-996)

    #pragma omp simd aligned(g_xz_xxz_y_x, g_xz_xxz_y_y, g_xz_xxz_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z, g_z_x_0_0_x_xz_y_x, g_z_x_0_0_x_xz_y_y, g_z_x_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xz_y_x[i] = -2.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_y_y[i] = -2.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_y_z[i] = -2.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (996-999)

    #pragma omp simd aligned(g_xz_xxz_z_x, g_xz_xxz_z_y, g_xz_xxz_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z, g_z_x_0_0_x_xz_z_x, g_z_x_0_0_x_xz_z_y, g_z_x_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xz_z_x[i] = -2.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_z_y[i] = -2.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_z_z[i] = -2.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (999-1002)

    #pragma omp simd aligned(g_xz_xyy_x_x, g_xz_xyy_x_y, g_xz_xyy_x_z, g_z_x_0_0_x_yy_x_x, g_z_x_0_0_x_yy_x_y, g_z_x_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yy_x_x[i] = 4.0 * g_xz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_x_y[i] = 4.0 * g_xz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_x_z[i] = 4.0 * g_xz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1002-1005)

    #pragma omp simd aligned(g_xz_xyy_y_x, g_xz_xyy_y_y, g_xz_xyy_y_z, g_z_x_0_0_x_yy_y_x, g_z_x_0_0_x_yy_y_y, g_z_x_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yy_y_x[i] = 4.0 * g_xz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_y_y[i] = 4.0 * g_xz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_y_z[i] = 4.0 * g_xz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1005-1008)

    #pragma omp simd aligned(g_xz_xyy_z_x, g_xz_xyy_z_y, g_xz_xyy_z_z, g_z_x_0_0_x_yy_z_x, g_z_x_0_0_x_yy_z_y, g_z_x_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yy_z_x[i] = 4.0 * g_xz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_z_y[i] = 4.0 * g_xz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_z_z[i] = 4.0 * g_xz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1008-1011)

    #pragma omp simd aligned(g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z, g_z_x_0_0_x_yz_x_x, g_z_x_0_0_x_yz_x_y, g_z_x_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yz_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1011-1014)

    #pragma omp simd aligned(g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z, g_z_x_0_0_x_yz_y_x, g_z_x_0_0_x_yz_y_y, g_z_x_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yz_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1014-1017)

    #pragma omp simd aligned(g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z, g_z_x_0_0_x_yz_z_x, g_z_x_0_0_x_yz_z_y, g_z_x_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yz_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1017-1020)

    #pragma omp simd aligned(g_xz_xzz_x_x, g_xz_xzz_x_y, g_xz_xzz_x_z, g_z_x_0_0_x_zz_x_x, g_z_x_0_0_x_zz_x_y, g_z_x_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_zz_x_x[i] = 4.0 * g_xz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_x_y[i] = 4.0 * g_xz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_x_z[i] = 4.0 * g_xz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1020-1023)

    #pragma omp simd aligned(g_xz_xzz_y_x, g_xz_xzz_y_y, g_xz_xzz_y_z, g_z_x_0_0_x_zz_y_x, g_z_x_0_0_x_zz_y_y, g_z_x_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_zz_y_x[i] = 4.0 * g_xz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_y_y[i] = 4.0 * g_xz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_y_z[i] = 4.0 * g_xz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1023-1026)

    #pragma omp simd aligned(g_xz_xzz_z_x, g_xz_xzz_z_y, g_xz_xzz_z_z, g_z_x_0_0_x_zz_z_x, g_z_x_0_0_x_zz_z_y, g_z_x_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_zz_z_x[i] = 4.0 * g_xz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_z_y[i] = 4.0 * g_xz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_z_z[i] = 4.0 * g_xz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1026-1029)

    #pragma omp simd aligned(g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xxx_x_x, g_yz_xxx_x_y, g_yz_xxx_x_z, g_z_x_0_0_y_xx_x_x, g_z_x_0_0_y_xx_x_y, g_z_x_0_0_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_x_x[i] = -4.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_x_y[i] = -4.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_x_z[i] = -4.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1029-1032)

    #pragma omp simd aligned(g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xxx_y_x, g_yz_xxx_y_y, g_yz_xxx_y_z, g_z_x_0_0_y_xx_y_x, g_z_x_0_0_y_xx_y_y, g_z_x_0_0_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_y_x[i] = -4.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_y_y[i] = -4.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_y_z[i] = -4.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1032-1035)

    #pragma omp simd aligned(g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xxx_z_x, g_yz_xxx_z_y, g_yz_xxx_z_z, g_z_x_0_0_y_xx_z_x, g_z_x_0_0_y_xx_z_y, g_z_x_0_0_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_z_x[i] = -4.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_z_y[i] = -4.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_z_z[i] = -4.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1035-1038)

    #pragma omp simd aligned(g_yz_xxy_x_x, g_yz_xxy_x_y, g_yz_xxy_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_z_x_0_0_y_xy_x_x, g_z_x_0_0_y_xy_x_y, g_z_x_0_0_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xy_x_x[i] = -2.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_x_y[i] = -2.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_x_z[i] = -2.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1038-1041)

    #pragma omp simd aligned(g_yz_xxy_y_x, g_yz_xxy_y_y, g_yz_xxy_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_z_x_0_0_y_xy_y_x, g_z_x_0_0_y_xy_y_y, g_z_x_0_0_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xy_y_x[i] = -2.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_y_y[i] = -2.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_y_z[i] = -2.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1041-1044)

    #pragma omp simd aligned(g_yz_xxy_z_x, g_yz_xxy_z_y, g_yz_xxy_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_z_x_0_0_y_xy_z_x, g_z_x_0_0_y_xy_z_y, g_z_x_0_0_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xy_z_x[i] = -2.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_z_y[i] = -2.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_z_z[i] = -2.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1044-1047)

    #pragma omp simd aligned(g_yz_xxz_x_x, g_yz_xxz_x_y, g_yz_xxz_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z, g_z_x_0_0_y_xz_x_x, g_z_x_0_0_y_xz_x_y, g_z_x_0_0_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xz_x_x[i] = -2.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_x_y[i] = -2.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_x_z[i] = -2.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1047-1050)

    #pragma omp simd aligned(g_yz_xxz_y_x, g_yz_xxz_y_y, g_yz_xxz_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z, g_z_x_0_0_y_xz_y_x, g_z_x_0_0_y_xz_y_y, g_z_x_0_0_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xz_y_x[i] = -2.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_y_y[i] = -2.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_y_z[i] = -2.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1050-1053)

    #pragma omp simd aligned(g_yz_xxz_z_x, g_yz_xxz_z_y, g_yz_xxz_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z, g_z_x_0_0_y_xz_z_x, g_z_x_0_0_y_xz_z_y, g_z_x_0_0_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xz_z_x[i] = -2.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_z_y[i] = -2.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_z_z[i] = -2.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1053-1056)

    #pragma omp simd aligned(g_yz_xyy_x_x, g_yz_xyy_x_y, g_yz_xyy_x_z, g_z_x_0_0_y_yy_x_x, g_z_x_0_0_y_yy_x_y, g_z_x_0_0_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yy_x_x[i] = 4.0 * g_yz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_x_y[i] = 4.0 * g_yz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_x_z[i] = 4.0 * g_yz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1056-1059)

    #pragma omp simd aligned(g_yz_xyy_y_x, g_yz_xyy_y_y, g_yz_xyy_y_z, g_z_x_0_0_y_yy_y_x, g_z_x_0_0_y_yy_y_y, g_z_x_0_0_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yy_y_x[i] = 4.0 * g_yz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_y_y[i] = 4.0 * g_yz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_y_z[i] = 4.0 * g_yz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1059-1062)

    #pragma omp simd aligned(g_yz_xyy_z_x, g_yz_xyy_z_y, g_yz_xyy_z_z, g_z_x_0_0_y_yy_z_x, g_z_x_0_0_y_yy_z_y, g_z_x_0_0_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yy_z_x[i] = 4.0 * g_yz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_z_y[i] = 4.0 * g_yz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_z_z[i] = 4.0 * g_yz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1062-1065)

    #pragma omp simd aligned(g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z, g_z_x_0_0_y_yz_x_x, g_z_x_0_0_y_yz_x_y, g_z_x_0_0_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yz_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1065-1068)

    #pragma omp simd aligned(g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z, g_z_x_0_0_y_yz_y_x, g_z_x_0_0_y_yz_y_y, g_z_x_0_0_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yz_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1068-1071)

    #pragma omp simd aligned(g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z, g_z_x_0_0_y_yz_z_x, g_z_x_0_0_y_yz_z_y, g_z_x_0_0_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yz_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1071-1074)

    #pragma omp simd aligned(g_yz_xzz_x_x, g_yz_xzz_x_y, g_yz_xzz_x_z, g_z_x_0_0_y_zz_x_x, g_z_x_0_0_y_zz_x_y, g_z_x_0_0_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_zz_x_x[i] = 4.0 * g_yz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_x_y[i] = 4.0 * g_yz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_x_z[i] = 4.0 * g_yz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1074-1077)

    #pragma omp simd aligned(g_yz_xzz_y_x, g_yz_xzz_y_y, g_yz_xzz_y_z, g_z_x_0_0_y_zz_y_x, g_z_x_0_0_y_zz_y_y, g_z_x_0_0_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_zz_y_x[i] = 4.0 * g_yz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_y_y[i] = 4.0 * g_yz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_y_z[i] = 4.0 * g_yz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1077-1080)

    #pragma omp simd aligned(g_yz_xzz_z_x, g_yz_xzz_z_y, g_yz_xzz_z_z, g_z_x_0_0_y_zz_z_x, g_z_x_0_0_y_zz_z_y, g_z_x_0_0_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_zz_z_x[i] = 4.0 * g_yz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_z_y[i] = 4.0 * g_yz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_z_z[i] = 4.0 * g_yz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1080-1083)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xxx_x_x, g_0_xxx_x_y, g_0_xxx_x_z, g_z_x_0_0_z_xx_x_x, g_z_x_0_0_z_xx_x_y, g_z_x_0_0_z_xx_x_z, g_zz_x_x_x, g_zz_x_x_y, g_zz_x_x_z, g_zz_xxx_x_x, g_zz_xxx_x_y, g_zz_xxx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_x_x[i] = 2.0 * g_0_x_x_x[i] - 2.0 * g_0_xxx_x_x[i] * b_exp - 4.0 * g_zz_x_x_x[i] * a_exp + 4.0 * g_zz_xxx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_x_y[i] = 2.0 * g_0_x_x_y[i] - 2.0 * g_0_xxx_x_y[i] * b_exp - 4.0 * g_zz_x_x_y[i] * a_exp + 4.0 * g_zz_xxx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_x_z[i] = 2.0 * g_0_x_x_z[i] - 2.0 * g_0_xxx_x_z[i] * b_exp - 4.0 * g_zz_x_x_z[i] * a_exp + 4.0 * g_zz_xxx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1083-1086)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xxx_y_x, g_0_xxx_y_y, g_0_xxx_y_z, g_z_x_0_0_z_xx_y_x, g_z_x_0_0_z_xx_y_y, g_z_x_0_0_z_xx_y_z, g_zz_x_y_x, g_zz_x_y_y, g_zz_x_y_z, g_zz_xxx_y_x, g_zz_xxx_y_y, g_zz_xxx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_y_x[i] = 2.0 * g_0_x_y_x[i] - 2.0 * g_0_xxx_y_x[i] * b_exp - 4.0 * g_zz_x_y_x[i] * a_exp + 4.0 * g_zz_xxx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_y_y[i] = 2.0 * g_0_x_y_y[i] - 2.0 * g_0_xxx_y_y[i] * b_exp - 4.0 * g_zz_x_y_y[i] * a_exp + 4.0 * g_zz_xxx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_y_z[i] = 2.0 * g_0_x_y_z[i] - 2.0 * g_0_xxx_y_z[i] * b_exp - 4.0 * g_zz_x_y_z[i] * a_exp + 4.0 * g_zz_xxx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1086-1089)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xxx_z_x, g_0_xxx_z_y, g_0_xxx_z_z, g_z_x_0_0_z_xx_z_x, g_z_x_0_0_z_xx_z_y, g_z_x_0_0_z_xx_z_z, g_zz_x_z_x, g_zz_x_z_y, g_zz_x_z_z, g_zz_xxx_z_x, g_zz_xxx_z_y, g_zz_xxx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_z_x[i] = 2.0 * g_0_x_z_x[i] - 2.0 * g_0_xxx_z_x[i] * b_exp - 4.0 * g_zz_x_z_x[i] * a_exp + 4.0 * g_zz_xxx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_z_y[i] = 2.0 * g_0_x_z_y[i] - 2.0 * g_0_xxx_z_y[i] * b_exp - 4.0 * g_zz_x_z_y[i] * a_exp + 4.0 * g_zz_xxx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_z_z[i] = 2.0 * g_0_x_z_z[i] - 2.0 * g_0_xxx_z_z[i] * b_exp - 4.0 * g_zz_x_z_z[i] * a_exp + 4.0 * g_zz_xxx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1089-1092)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_z_x_0_0_z_xy_x_x, g_z_x_0_0_z_xy_x_y, g_z_x_0_0_z_xy_x_z, g_zz_xxy_x_x, g_zz_xxy_x_y, g_zz_xxy_x_z, g_zz_y_x_x, g_zz_y_x_y, g_zz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xy_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_xxy_x_x[i] * b_exp - 2.0 * g_zz_y_x_x[i] * a_exp + 4.0 * g_zz_xxy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_xxy_x_y[i] * b_exp - 2.0 * g_zz_y_x_y[i] * a_exp + 4.0 * g_zz_xxy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_xxy_x_z[i] * b_exp - 2.0 * g_zz_y_x_z[i] * a_exp + 4.0 * g_zz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1092-1095)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_z_x_0_0_z_xy_y_x, g_z_x_0_0_z_xy_y_y, g_z_x_0_0_z_xy_y_z, g_zz_xxy_y_x, g_zz_xxy_y_y, g_zz_xxy_y_z, g_zz_y_y_x, g_zz_y_y_y, g_zz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xy_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_xxy_y_x[i] * b_exp - 2.0 * g_zz_y_y_x[i] * a_exp + 4.0 * g_zz_xxy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_xxy_y_y[i] * b_exp - 2.0 * g_zz_y_y_y[i] * a_exp + 4.0 * g_zz_xxy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_xxy_y_z[i] * b_exp - 2.0 * g_zz_y_y_z[i] * a_exp + 4.0 * g_zz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1095-1098)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_z_x_0_0_z_xy_z_x, g_z_x_0_0_z_xy_z_y, g_z_x_0_0_z_xy_z_z, g_zz_xxy_z_x, g_zz_xxy_z_y, g_zz_xxy_z_z, g_zz_y_z_x, g_zz_y_z_y, g_zz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xy_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_xxy_z_x[i] * b_exp - 2.0 * g_zz_y_z_x[i] * a_exp + 4.0 * g_zz_xxy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_xxy_z_y[i] * b_exp - 2.0 * g_zz_y_z_y[i] * a_exp + 4.0 * g_zz_xxy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_xxy_z_z[i] * b_exp - 2.0 * g_zz_y_z_z[i] * a_exp + 4.0 * g_zz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1098-1101)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_z_x_0_0_z_xz_x_x, g_z_x_0_0_z_xz_x_y, g_z_x_0_0_z_xz_x_z, g_zz_xxz_x_x, g_zz_xxz_x_y, g_zz_xxz_x_z, g_zz_z_x_x, g_zz_z_x_y, g_zz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_xxz_x_x[i] * b_exp - 2.0 * g_zz_z_x_x[i] * a_exp + 4.0 * g_zz_xxz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_xxz_x_y[i] * b_exp - 2.0 * g_zz_z_x_y[i] * a_exp + 4.0 * g_zz_xxz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_xxz_x_z[i] * b_exp - 2.0 * g_zz_z_x_z[i] * a_exp + 4.0 * g_zz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1101-1104)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_z_x_0_0_z_xz_y_x, g_z_x_0_0_z_xz_y_y, g_z_x_0_0_z_xz_y_z, g_zz_xxz_y_x, g_zz_xxz_y_y, g_zz_xxz_y_z, g_zz_z_y_x, g_zz_z_y_y, g_zz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_xxz_y_x[i] * b_exp - 2.0 * g_zz_z_y_x[i] * a_exp + 4.0 * g_zz_xxz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_xxz_y_y[i] * b_exp - 2.0 * g_zz_z_y_y[i] * a_exp + 4.0 * g_zz_xxz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_xxz_y_z[i] * b_exp - 2.0 * g_zz_z_y_z[i] * a_exp + 4.0 * g_zz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1104-1107)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_z_x_0_0_z_xz_z_x, g_z_x_0_0_z_xz_z_y, g_z_x_0_0_z_xz_z_z, g_zz_xxz_z_x, g_zz_xxz_z_y, g_zz_xxz_z_z, g_zz_z_z_x, g_zz_z_z_y, g_zz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_xxz_z_x[i] * b_exp - 2.0 * g_zz_z_z_x[i] * a_exp + 4.0 * g_zz_xxz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_xxz_z_y[i] * b_exp - 2.0 * g_zz_z_z_y[i] * a_exp + 4.0 * g_zz_xxz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_xxz_z_z[i] * b_exp - 2.0 * g_zz_z_z_z[i] * a_exp + 4.0 * g_zz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1107-1110)

    #pragma omp simd aligned(g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_z_x_0_0_z_yy_x_x, g_z_x_0_0_z_yy_x_y, g_z_x_0_0_z_yy_x_z, g_zz_xyy_x_x, g_zz_xyy_x_y, g_zz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yy_x_x[i] = -2.0 * g_0_xyy_x_x[i] * b_exp + 4.0 * g_zz_xyy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_x_y[i] = -2.0 * g_0_xyy_x_y[i] * b_exp + 4.0 * g_zz_xyy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_x_z[i] = -2.0 * g_0_xyy_x_z[i] * b_exp + 4.0 * g_zz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1110-1113)

    #pragma omp simd aligned(g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_z_x_0_0_z_yy_y_x, g_z_x_0_0_z_yy_y_y, g_z_x_0_0_z_yy_y_z, g_zz_xyy_y_x, g_zz_xyy_y_y, g_zz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yy_y_x[i] = -2.0 * g_0_xyy_y_x[i] * b_exp + 4.0 * g_zz_xyy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_y_y[i] = -2.0 * g_0_xyy_y_y[i] * b_exp + 4.0 * g_zz_xyy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_y_z[i] = -2.0 * g_0_xyy_y_z[i] * b_exp + 4.0 * g_zz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1113-1116)

    #pragma omp simd aligned(g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_z_x_0_0_z_yy_z_x, g_z_x_0_0_z_yy_z_y, g_z_x_0_0_z_yy_z_z, g_zz_xyy_z_x, g_zz_xyy_z_y, g_zz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yy_z_x[i] = -2.0 * g_0_xyy_z_x[i] * b_exp + 4.0 * g_zz_xyy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_z_y[i] = -2.0 * g_0_xyy_z_y[i] * b_exp + 4.0 * g_zz_xyy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_z_z[i] = -2.0 * g_0_xyy_z_z[i] * b_exp + 4.0 * g_zz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1116-1119)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_z_x_0_0_z_yz_x_x, g_z_x_0_0_z_yz_x_y, g_z_x_0_0_z_yz_x_z, g_zz_xyz_x_x, g_zz_xyz_x_y, g_zz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_zz_xyz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_zz_xyz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_zz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1119-1122)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_z_x_0_0_z_yz_y_x, g_z_x_0_0_z_yz_y_y, g_z_x_0_0_z_yz_y_z, g_zz_xyz_y_x, g_zz_xyz_y_y, g_zz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_zz_xyz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_zz_xyz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_zz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1122-1125)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_z_x_0_0_z_yz_z_x, g_z_x_0_0_z_yz_z_y, g_z_x_0_0_z_yz_z_z, g_zz_xyz_z_x, g_zz_xyz_z_y, g_zz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_zz_xyz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_zz_xyz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_zz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1125-1128)

    #pragma omp simd aligned(g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_z_x_0_0_z_zz_x_x, g_z_x_0_0_z_zz_x_y, g_z_x_0_0_z_zz_x_z, g_zz_xzz_x_x, g_zz_xzz_x_y, g_zz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_zz_x_x[i] = -2.0 * g_0_xzz_x_x[i] * b_exp + 4.0 * g_zz_xzz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_x_y[i] = -2.0 * g_0_xzz_x_y[i] * b_exp + 4.0 * g_zz_xzz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_x_z[i] = -2.0 * g_0_xzz_x_z[i] * b_exp + 4.0 * g_zz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1128-1131)

    #pragma omp simd aligned(g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_z_x_0_0_z_zz_y_x, g_z_x_0_0_z_zz_y_y, g_z_x_0_0_z_zz_y_z, g_zz_xzz_y_x, g_zz_xzz_y_y, g_zz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_zz_y_x[i] = -2.0 * g_0_xzz_y_x[i] * b_exp + 4.0 * g_zz_xzz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_y_y[i] = -2.0 * g_0_xzz_y_y[i] * b_exp + 4.0 * g_zz_xzz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_y_z[i] = -2.0 * g_0_xzz_y_z[i] * b_exp + 4.0 * g_zz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1131-1134)

    #pragma omp simd aligned(g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_z_x_0_0_z_zz_z_x, g_z_x_0_0_z_zz_z_y, g_z_x_0_0_z_zz_z_z, g_zz_xzz_z_x, g_zz_xzz_z_y, g_zz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_zz_z_x[i] = -2.0 * g_0_xzz_z_x[i] * b_exp + 4.0 * g_zz_xzz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_z_y[i] = -2.0 * g_0_xzz_z_y[i] * b_exp + 4.0 * g_zz_xzz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_z_z[i] = -2.0 * g_0_xzz_z_z[i] * b_exp + 4.0 * g_zz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1134-1137)

    #pragma omp simd aligned(g_xz_xxy_x_x, g_xz_xxy_x_y, g_xz_xxy_x_z, g_z_y_0_0_x_xx_x_x, g_z_y_0_0_x_xx_x_y, g_z_y_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_x_x[i] = 4.0 * g_xz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_x_y[i] = 4.0 * g_xz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_x_z[i] = 4.0 * g_xz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1137-1140)

    #pragma omp simd aligned(g_xz_xxy_y_x, g_xz_xxy_y_y, g_xz_xxy_y_z, g_z_y_0_0_x_xx_y_x, g_z_y_0_0_x_xx_y_y, g_z_y_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_y_x[i] = 4.0 * g_xz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_y_y[i] = 4.0 * g_xz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_y_z[i] = 4.0 * g_xz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1140-1143)

    #pragma omp simd aligned(g_xz_xxy_z_x, g_xz_xxy_z_y, g_xz_xxy_z_z, g_z_y_0_0_x_xx_z_x, g_z_y_0_0_x_xx_z_y, g_z_y_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_z_x[i] = 4.0 * g_xz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_z_y[i] = 4.0 * g_xz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_z_z[i] = 4.0 * g_xz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1143-1146)

    #pragma omp simd aligned(g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xyy_x_x, g_xz_xyy_x_y, g_xz_xyy_x_z, g_z_y_0_0_x_xy_x_x, g_z_y_0_0_x_xy_x_y, g_z_y_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xy_x_x[i] = -2.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_x_y[i] = -2.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_x_z[i] = -2.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1146-1149)

    #pragma omp simd aligned(g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xyy_y_x, g_xz_xyy_y_y, g_xz_xyy_y_z, g_z_y_0_0_x_xy_y_x, g_z_y_0_0_x_xy_y_y, g_z_y_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xy_y_x[i] = -2.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_y_y[i] = -2.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_y_z[i] = -2.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1149-1152)

    #pragma omp simd aligned(g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xyy_z_x, g_xz_xyy_z_y, g_xz_xyy_z_z, g_z_y_0_0_x_xy_z_x, g_z_y_0_0_x_xy_z_y, g_z_y_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xy_z_x[i] = -2.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_z_y[i] = -2.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_z_z[i] = -2.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1152-1155)

    #pragma omp simd aligned(g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z, g_z_y_0_0_x_xz_x_x, g_z_y_0_0_x_xz_x_y, g_z_y_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xz_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1155-1158)

    #pragma omp simd aligned(g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z, g_z_y_0_0_x_xz_y_x, g_z_y_0_0_x_xz_y_y, g_z_y_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xz_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1158-1161)

    #pragma omp simd aligned(g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z, g_z_y_0_0_x_xz_z_x, g_z_y_0_0_x_xz_z_y, g_z_y_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xz_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1161-1164)

    #pragma omp simd aligned(g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_xz_yyy_x_x, g_xz_yyy_x_y, g_xz_yyy_x_z, g_z_y_0_0_x_yy_x_x, g_z_y_0_0_x_yy_x_y, g_z_y_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yy_x_x[i] = -4.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_x_y[i] = -4.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_x_z[i] = -4.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1164-1167)

    #pragma omp simd aligned(g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_xz_yyy_y_x, g_xz_yyy_y_y, g_xz_yyy_y_z, g_z_y_0_0_x_yy_y_x, g_z_y_0_0_x_yy_y_y, g_z_y_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yy_y_x[i] = -4.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_y_y[i] = -4.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_y_z[i] = -4.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1167-1170)

    #pragma omp simd aligned(g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_xz_yyy_z_x, g_xz_yyy_z_y, g_xz_yyy_z_z, g_z_y_0_0_x_yy_z_x, g_z_y_0_0_x_yy_z_y, g_z_y_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yy_z_x[i] = -4.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_z_y[i] = -4.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_z_z[i] = -4.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1170-1173)

    #pragma omp simd aligned(g_xz_yyz_x_x, g_xz_yyz_x_y, g_xz_yyz_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z, g_z_y_0_0_x_yz_x_x, g_z_y_0_0_x_yz_x_y, g_z_y_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yz_x_x[i] = -2.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_x_y[i] = -2.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_x_z[i] = -2.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1173-1176)

    #pragma omp simd aligned(g_xz_yyz_y_x, g_xz_yyz_y_y, g_xz_yyz_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z, g_z_y_0_0_x_yz_y_x, g_z_y_0_0_x_yz_y_y, g_z_y_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yz_y_x[i] = -2.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_y_y[i] = -2.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_y_z[i] = -2.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1176-1179)

    #pragma omp simd aligned(g_xz_yyz_z_x, g_xz_yyz_z_y, g_xz_yyz_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z, g_z_y_0_0_x_yz_z_x, g_z_y_0_0_x_yz_z_y, g_z_y_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yz_z_x[i] = -2.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_z_y[i] = -2.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_z_z[i] = -2.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1179-1182)

    #pragma omp simd aligned(g_xz_yzz_x_x, g_xz_yzz_x_y, g_xz_yzz_x_z, g_z_y_0_0_x_zz_x_x, g_z_y_0_0_x_zz_x_y, g_z_y_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_zz_x_x[i] = 4.0 * g_xz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_x_y[i] = 4.0 * g_xz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_x_z[i] = 4.0 * g_xz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1182-1185)

    #pragma omp simd aligned(g_xz_yzz_y_x, g_xz_yzz_y_y, g_xz_yzz_y_z, g_z_y_0_0_x_zz_y_x, g_z_y_0_0_x_zz_y_y, g_z_y_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_zz_y_x[i] = 4.0 * g_xz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_y_y[i] = 4.0 * g_xz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_y_z[i] = 4.0 * g_xz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1185-1188)

    #pragma omp simd aligned(g_xz_yzz_z_x, g_xz_yzz_z_y, g_xz_yzz_z_z, g_z_y_0_0_x_zz_z_x, g_z_y_0_0_x_zz_z_y, g_z_y_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_zz_z_x[i] = 4.0 * g_xz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_z_y[i] = 4.0 * g_xz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_z_z[i] = 4.0 * g_xz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1188-1191)

    #pragma omp simd aligned(g_yz_xxy_x_x, g_yz_xxy_x_y, g_yz_xxy_x_z, g_z_y_0_0_y_xx_x_x, g_z_y_0_0_y_xx_x_y, g_z_y_0_0_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_x_x[i] = 4.0 * g_yz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_x_y[i] = 4.0 * g_yz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_x_z[i] = 4.0 * g_yz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1191-1194)

    #pragma omp simd aligned(g_yz_xxy_y_x, g_yz_xxy_y_y, g_yz_xxy_y_z, g_z_y_0_0_y_xx_y_x, g_z_y_0_0_y_xx_y_y, g_z_y_0_0_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_y_x[i] = 4.0 * g_yz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_y_y[i] = 4.0 * g_yz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_y_z[i] = 4.0 * g_yz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1194-1197)

    #pragma omp simd aligned(g_yz_xxy_z_x, g_yz_xxy_z_y, g_yz_xxy_z_z, g_z_y_0_0_y_xx_z_x, g_z_y_0_0_y_xx_z_y, g_z_y_0_0_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_z_x[i] = 4.0 * g_yz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_z_y[i] = 4.0 * g_yz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_z_z[i] = 4.0 * g_yz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1197-1200)

    #pragma omp simd aligned(g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xyy_x_x, g_yz_xyy_x_y, g_yz_xyy_x_z, g_z_y_0_0_y_xy_x_x, g_z_y_0_0_y_xy_x_y, g_z_y_0_0_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xy_x_x[i] = -2.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_x_y[i] = -2.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_x_z[i] = -2.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1200-1203)

    #pragma omp simd aligned(g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xyy_y_x, g_yz_xyy_y_y, g_yz_xyy_y_z, g_z_y_0_0_y_xy_y_x, g_z_y_0_0_y_xy_y_y, g_z_y_0_0_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xy_y_x[i] = -2.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_y_y[i] = -2.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_y_z[i] = -2.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1203-1206)

    #pragma omp simd aligned(g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xyy_z_x, g_yz_xyy_z_y, g_yz_xyy_z_z, g_z_y_0_0_y_xy_z_x, g_z_y_0_0_y_xy_z_y, g_z_y_0_0_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xy_z_x[i] = -2.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_z_y[i] = -2.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_z_z[i] = -2.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1206-1209)

    #pragma omp simd aligned(g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z, g_z_y_0_0_y_xz_x_x, g_z_y_0_0_y_xz_x_y, g_z_y_0_0_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xz_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1209-1212)

    #pragma omp simd aligned(g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z, g_z_y_0_0_y_xz_y_x, g_z_y_0_0_y_xz_y_y, g_z_y_0_0_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xz_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1212-1215)

    #pragma omp simd aligned(g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z, g_z_y_0_0_y_xz_z_x, g_z_y_0_0_y_xz_z_y, g_z_y_0_0_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xz_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1215-1218)

    #pragma omp simd aligned(g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_yz_yyy_x_x, g_yz_yyy_x_y, g_yz_yyy_x_z, g_z_y_0_0_y_yy_x_x, g_z_y_0_0_y_yy_x_y, g_z_y_0_0_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yy_x_x[i] = -4.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_x_y[i] = -4.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_x_z[i] = -4.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1218-1221)

    #pragma omp simd aligned(g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_yz_yyy_y_x, g_yz_yyy_y_y, g_yz_yyy_y_z, g_z_y_0_0_y_yy_y_x, g_z_y_0_0_y_yy_y_y, g_z_y_0_0_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yy_y_x[i] = -4.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_y_y[i] = -4.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_y_z[i] = -4.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1221-1224)

    #pragma omp simd aligned(g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_yz_yyy_z_x, g_yz_yyy_z_y, g_yz_yyy_z_z, g_z_y_0_0_y_yy_z_x, g_z_y_0_0_y_yy_z_y, g_z_y_0_0_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yy_z_x[i] = -4.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_z_y[i] = -4.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_z_z[i] = -4.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1224-1227)

    #pragma omp simd aligned(g_yz_yyz_x_x, g_yz_yyz_x_y, g_yz_yyz_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z, g_z_y_0_0_y_yz_x_x, g_z_y_0_0_y_yz_x_y, g_z_y_0_0_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yz_x_x[i] = -2.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_x_y[i] = -2.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_x_z[i] = -2.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1227-1230)

    #pragma omp simd aligned(g_yz_yyz_y_x, g_yz_yyz_y_y, g_yz_yyz_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z, g_z_y_0_0_y_yz_y_x, g_z_y_0_0_y_yz_y_y, g_z_y_0_0_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yz_y_x[i] = -2.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_y_y[i] = -2.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_y_z[i] = -2.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1230-1233)

    #pragma omp simd aligned(g_yz_yyz_z_x, g_yz_yyz_z_y, g_yz_yyz_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z, g_z_y_0_0_y_yz_z_x, g_z_y_0_0_y_yz_z_y, g_z_y_0_0_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yz_z_x[i] = -2.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_z_y[i] = -2.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_z_z[i] = -2.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1233-1236)

    #pragma omp simd aligned(g_yz_yzz_x_x, g_yz_yzz_x_y, g_yz_yzz_x_z, g_z_y_0_0_y_zz_x_x, g_z_y_0_0_y_zz_x_y, g_z_y_0_0_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_zz_x_x[i] = 4.0 * g_yz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_x_y[i] = 4.0 * g_yz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_x_z[i] = 4.0 * g_yz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1236-1239)

    #pragma omp simd aligned(g_yz_yzz_y_x, g_yz_yzz_y_y, g_yz_yzz_y_z, g_z_y_0_0_y_zz_y_x, g_z_y_0_0_y_zz_y_y, g_z_y_0_0_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_zz_y_x[i] = 4.0 * g_yz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_y_y[i] = 4.0 * g_yz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_y_z[i] = 4.0 * g_yz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1239-1242)

    #pragma omp simd aligned(g_yz_yzz_z_x, g_yz_yzz_z_y, g_yz_yzz_z_z, g_z_y_0_0_y_zz_z_x, g_z_y_0_0_y_zz_z_y, g_z_y_0_0_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_zz_z_x[i] = 4.0 * g_yz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_z_y[i] = 4.0 * g_yz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_z_z[i] = 4.0 * g_yz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1242-1245)

    #pragma omp simd aligned(g_0_xxy_x_x, g_0_xxy_x_y, g_0_xxy_x_z, g_z_y_0_0_z_xx_x_x, g_z_y_0_0_z_xx_x_y, g_z_y_0_0_z_xx_x_z, g_zz_xxy_x_x, g_zz_xxy_x_y, g_zz_xxy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_x_x[i] = -2.0 * g_0_xxy_x_x[i] * b_exp + 4.0 * g_zz_xxy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_x_y[i] = -2.0 * g_0_xxy_x_y[i] * b_exp + 4.0 * g_zz_xxy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_x_z[i] = -2.0 * g_0_xxy_x_z[i] * b_exp + 4.0 * g_zz_xxy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1245-1248)

    #pragma omp simd aligned(g_0_xxy_y_x, g_0_xxy_y_y, g_0_xxy_y_z, g_z_y_0_0_z_xx_y_x, g_z_y_0_0_z_xx_y_y, g_z_y_0_0_z_xx_y_z, g_zz_xxy_y_x, g_zz_xxy_y_y, g_zz_xxy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_y_x[i] = -2.0 * g_0_xxy_y_x[i] * b_exp + 4.0 * g_zz_xxy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_y_y[i] = -2.0 * g_0_xxy_y_y[i] * b_exp + 4.0 * g_zz_xxy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_y_z[i] = -2.0 * g_0_xxy_y_z[i] * b_exp + 4.0 * g_zz_xxy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1248-1251)

    #pragma omp simd aligned(g_0_xxy_z_x, g_0_xxy_z_y, g_0_xxy_z_z, g_z_y_0_0_z_xx_z_x, g_z_y_0_0_z_xx_z_y, g_z_y_0_0_z_xx_z_z, g_zz_xxy_z_x, g_zz_xxy_z_y, g_zz_xxy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_z_x[i] = -2.0 * g_0_xxy_z_x[i] * b_exp + 4.0 * g_zz_xxy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_z_y[i] = -2.0 * g_0_xxy_z_y[i] * b_exp + 4.0 * g_zz_xxy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_z_z[i] = -2.0 * g_0_xxy_z_z[i] * b_exp + 4.0 * g_zz_xxy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1251-1254)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xyy_x_x, g_0_xyy_x_y, g_0_xyy_x_z, g_z_y_0_0_z_xy_x_x, g_z_y_0_0_z_xy_x_y, g_z_y_0_0_z_xy_x_z, g_zz_x_x_x, g_zz_x_x_y, g_zz_x_x_z, g_zz_xyy_x_x, g_zz_xyy_x_y, g_zz_xyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xy_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xyy_x_x[i] * b_exp - 2.0 * g_zz_x_x_x[i] * a_exp + 4.0 * g_zz_xyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xyy_x_y[i] * b_exp - 2.0 * g_zz_x_x_y[i] * a_exp + 4.0 * g_zz_xyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xyy_x_z[i] * b_exp - 2.0 * g_zz_x_x_z[i] * a_exp + 4.0 * g_zz_xyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1254-1257)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xyy_y_x, g_0_xyy_y_y, g_0_xyy_y_z, g_z_y_0_0_z_xy_y_x, g_z_y_0_0_z_xy_y_y, g_z_y_0_0_z_xy_y_z, g_zz_x_y_x, g_zz_x_y_y, g_zz_x_y_z, g_zz_xyy_y_x, g_zz_xyy_y_y, g_zz_xyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xy_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xyy_y_x[i] * b_exp - 2.0 * g_zz_x_y_x[i] * a_exp + 4.0 * g_zz_xyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xyy_y_y[i] * b_exp - 2.0 * g_zz_x_y_y[i] * a_exp + 4.0 * g_zz_xyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xyy_y_z[i] * b_exp - 2.0 * g_zz_x_y_z[i] * a_exp + 4.0 * g_zz_xyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1257-1260)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xyy_z_x, g_0_xyy_z_y, g_0_xyy_z_z, g_z_y_0_0_z_xy_z_x, g_z_y_0_0_z_xy_z_y, g_z_y_0_0_z_xy_z_z, g_zz_x_z_x, g_zz_x_z_y, g_zz_x_z_z, g_zz_xyy_z_x, g_zz_xyy_z_y, g_zz_xyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xy_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xyy_z_x[i] * b_exp - 2.0 * g_zz_x_z_x[i] * a_exp + 4.0 * g_zz_xyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xyy_z_y[i] * b_exp - 2.0 * g_zz_x_z_y[i] * a_exp + 4.0 * g_zz_xyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xyy_z_z[i] * b_exp - 2.0 * g_zz_x_z_z[i] * a_exp + 4.0 * g_zz_xyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1260-1263)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_z_y_0_0_z_xz_x_x, g_z_y_0_0_z_xz_x_y, g_z_y_0_0_z_xz_x_z, g_zz_xyz_x_x, g_zz_xyz_x_y, g_zz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xz_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_zz_xyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_zz_xyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_zz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1263-1266)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_z_y_0_0_z_xz_y_x, g_z_y_0_0_z_xz_y_y, g_z_y_0_0_z_xz_y_z, g_zz_xyz_y_x, g_zz_xyz_y_y, g_zz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xz_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_zz_xyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_zz_xyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_zz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1266-1269)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_z_y_0_0_z_xz_z_x, g_z_y_0_0_z_xz_z_y, g_z_y_0_0_z_xz_z_z, g_zz_xyz_z_x, g_zz_xyz_z_y, g_zz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xz_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_zz_xyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_zz_xyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_zz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1269-1272)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yyy_x_x, g_0_yyy_x_y, g_0_yyy_x_z, g_z_y_0_0_z_yy_x_x, g_z_y_0_0_z_yy_x_y, g_z_y_0_0_z_yy_x_z, g_zz_y_x_x, g_zz_y_x_y, g_zz_y_x_z, g_zz_yyy_x_x, g_zz_yyy_x_y, g_zz_yyy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yy_x_x[i] = 2.0 * g_0_y_x_x[i] - 2.0 * g_0_yyy_x_x[i] * b_exp - 4.0 * g_zz_y_x_x[i] * a_exp + 4.0 * g_zz_yyy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_x_y[i] = 2.0 * g_0_y_x_y[i] - 2.0 * g_0_yyy_x_y[i] * b_exp - 4.0 * g_zz_y_x_y[i] * a_exp + 4.0 * g_zz_yyy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_x_z[i] = 2.0 * g_0_y_x_z[i] - 2.0 * g_0_yyy_x_z[i] * b_exp - 4.0 * g_zz_y_x_z[i] * a_exp + 4.0 * g_zz_yyy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1272-1275)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yyy_y_x, g_0_yyy_y_y, g_0_yyy_y_z, g_z_y_0_0_z_yy_y_x, g_z_y_0_0_z_yy_y_y, g_z_y_0_0_z_yy_y_z, g_zz_y_y_x, g_zz_y_y_y, g_zz_y_y_z, g_zz_yyy_y_x, g_zz_yyy_y_y, g_zz_yyy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yy_y_x[i] = 2.0 * g_0_y_y_x[i] - 2.0 * g_0_yyy_y_x[i] * b_exp - 4.0 * g_zz_y_y_x[i] * a_exp + 4.0 * g_zz_yyy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_y_y[i] = 2.0 * g_0_y_y_y[i] - 2.0 * g_0_yyy_y_y[i] * b_exp - 4.0 * g_zz_y_y_y[i] * a_exp + 4.0 * g_zz_yyy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_y_z[i] = 2.0 * g_0_y_y_z[i] - 2.0 * g_0_yyy_y_z[i] * b_exp - 4.0 * g_zz_y_y_z[i] * a_exp + 4.0 * g_zz_yyy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1275-1278)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yyy_z_x, g_0_yyy_z_y, g_0_yyy_z_z, g_z_y_0_0_z_yy_z_x, g_z_y_0_0_z_yy_z_y, g_z_y_0_0_z_yy_z_z, g_zz_y_z_x, g_zz_y_z_y, g_zz_y_z_z, g_zz_yyy_z_x, g_zz_yyy_z_y, g_zz_yyy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yy_z_x[i] = 2.0 * g_0_y_z_x[i] - 2.0 * g_0_yyy_z_x[i] * b_exp - 4.0 * g_zz_y_z_x[i] * a_exp + 4.0 * g_zz_yyy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_z_y[i] = 2.0 * g_0_y_z_y[i] - 2.0 * g_0_yyy_z_y[i] * b_exp - 4.0 * g_zz_y_z_y[i] * a_exp + 4.0 * g_zz_yyy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_z_z[i] = 2.0 * g_0_y_z_z[i] - 2.0 * g_0_yyy_z_z[i] * b_exp - 4.0 * g_zz_y_z_z[i] * a_exp + 4.0 * g_zz_yyy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1278-1281)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_z_y_0_0_z_yz_x_x, g_z_y_0_0_z_yz_x_y, g_z_y_0_0_z_yz_x_z, g_zz_yyz_x_x, g_zz_yyz_x_y, g_zz_yyz_x_z, g_zz_z_x_x, g_zz_z_x_y, g_zz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yz_x_x[i] = g_0_z_x_x[i] - 2.0 * g_0_yyz_x_x[i] * b_exp - 2.0 * g_zz_z_x_x[i] * a_exp + 4.0 * g_zz_yyz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_x_y[i] = g_0_z_x_y[i] - 2.0 * g_0_yyz_x_y[i] * b_exp - 2.0 * g_zz_z_x_y[i] * a_exp + 4.0 * g_zz_yyz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_x_z[i] = g_0_z_x_z[i] - 2.0 * g_0_yyz_x_z[i] * b_exp - 2.0 * g_zz_z_x_z[i] * a_exp + 4.0 * g_zz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1281-1284)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_z_y_0_0_z_yz_y_x, g_z_y_0_0_z_yz_y_y, g_z_y_0_0_z_yz_y_z, g_zz_yyz_y_x, g_zz_yyz_y_y, g_zz_yyz_y_z, g_zz_z_y_x, g_zz_z_y_y, g_zz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yz_y_x[i] = g_0_z_y_x[i] - 2.0 * g_0_yyz_y_x[i] * b_exp - 2.0 * g_zz_z_y_x[i] * a_exp + 4.0 * g_zz_yyz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_y_y[i] = g_0_z_y_y[i] - 2.0 * g_0_yyz_y_y[i] * b_exp - 2.0 * g_zz_z_y_y[i] * a_exp + 4.0 * g_zz_yyz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_y_z[i] = g_0_z_y_z[i] - 2.0 * g_0_yyz_y_z[i] * b_exp - 2.0 * g_zz_z_y_z[i] * a_exp + 4.0 * g_zz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1284-1287)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_z_y_0_0_z_yz_z_x, g_z_y_0_0_z_yz_z_y, g_z_y_0_0_z_yz_z_z, g_zz_yyz_z_x, g_zz_yyz_z_y, g_zz_yyz_z_z, g_zz_z_z_x, g_zz_z_z_y, g_zz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yz_z_x[i] = g_0_z_z_x[i] - 2.0 * g_0_yyz_z_x[i] * b_exp - 2.0 * g_zz_z_z_x[i] * a_exp + 4.0 * g_zz_yyz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_z_y[i] = g_0_z_z_y[i] - 2.0 * g_0_yyz_z_y[i] * b_exp - 2.0 * g_zz_z_z_y[i] * a_exp + 4.0 * g_zz_yyz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_z_z[i] = g_0_z_z_z[i] - 2.0 * g_0_yyz_z_z[i] * b_exp - 2.0 * g_zz_z_z_z[i] * a_exp + 4.0 * g_zz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1287-1290)

    #pragma omp simd aligned(g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_z_y_0_0_z_zz_x_x, g_z_y_0_0_z_zz_x_y, g_z_y_0_0_z_zz_x_z, g_zz_yzz_x_x, g_zz_yzz_x_y, g_zz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_zz_x_x[i] = -2.0 * g_0_yzz_x_x[i] * b_exp + 4.0 * g_zz_yzz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_x_y[i] = -2.0 * g_0_yzz_x_y[i] * b_exp + 4.0 * g_zz_yzz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_x_z[i] = -2.0 * g_0_yzz_x_z[i] * b_exp + 4.0 * g_zz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1290-1293)

    #pragma omp simd aligned(g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_z_y_0_0_z_zz_y_x, g_z_y_0_0_z_zz_y_y, g_z_y_0_0_z_zz_y_z, g_zz_yzz_y_x, g_zz_yzz_y_y, g_zz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_zz_y_x[i] = -2.0 * g_0_yzz_y_x[i] * b_exp + 4.0 * g_zz_yzz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_y_y[i] = -2.0 * g_0_yzz_y_y[i] * b_exp + 4.0 * g_zz_yzz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_y_z[i] = -2.0 * g_0_yzz_y_z[i] * b_exp + 4.0 * g_zz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1293-1296)

    #pragma omp simd aligned(g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_z_y_0_0_z_zz_z_x, g_z_y_0_0_z_zz_z_y, g_z_y_0_0_z_zz_z_z, g_zz_yzz_z_x, g_zz_yzz_z_y, g_zz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_zz_z_x[i] = -2.0 * g_0_yzz_z_x[i] * b_exp + 4.0 * g_zz_yzz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_z_y[i] = -2.0 * g_0_yzz_z_y[i] * b_exp + 4.0 * g_zz_yzz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_z_z[i] = -2.0 * g_0_yzz_z_z[i] * b_exp + 4.0 * g_zz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1296-1299)

    #pragma omp simd aligned(g_xz_xxz_x_x, g_xz_xxz_x_y, g_xz_xxz_x_z, g_z_z_0_0_x_xx_x_x, g_z_z_0_0_x_xx_x_y, g_z_z_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_x_x[i] = 4.0 * g_xz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_x_y[i] = 4.0 * g_xz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_x_z[i] = 4.0 * g_xz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1299-1302)

    #pragma omp simd aligned(g_xz_xxz_y_x, g_xz_xxz_y_y, g_xz_xxz_y_z, g_z_z_0_0_x_xx_y_x, g_z_z_0_0_x_xx_y_y, g_z_z_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_y_x[i] = 4.0 * g_xz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_y_y[i] = 4.0 * g_xz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_y_z[i] = 4.0 * g_xz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1302-1305)

    #pragma omp simd aligned(g_xz_xxz_z_x, g_xz_xxz_z_y, g_xz_xxz_z_z, g_z_z_0_0_x_xx_z_x, g_z_z_0_0_x_xx_z_y, g_z_z_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_z_x[i] = 4.0 * g_xz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_z_y[i] = 4.0 * g_xz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_z_z[i] = 4.0 * g_xz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1305-1308)

    #pragma omp simd aligned(g_xz_xyz_x_x, g_xz_xyz_x_y, g_xz_xyz_x_z, g_z_z_0_0_x_xy_x_x, g_z_z_0_0_x_xy_x_y, g_z_z_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xy_x_x[i] = 4.0 * g_xz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_x_y[i] = 4.0 * g_xz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_x_z[i] = 4.0 * g_xz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1308-1311)

    #pragma omp simd aligned(g_xz_xyz_y_x, g_xz_xyz_y_y, g_xz_xyz_y_z, g_z_z_0_0_x_xy_y_x, g_z_z_0_0_x_xy_y_y, g_z_z_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xy_y_x[i] = 4.0 * g_xz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_y_y[i] = 4.0 * g_xz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_y_z[i] = 4.0 * g_xz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1311-1314)

    #pragma omp simd aligned(g_xz_xyz_z_x, g_xz_xyz_z_y, g_xz_xyz_z_z, g_z_z_0_0_x_xy_z_x, g_z_z_0_0_x_xy_z_y, g_z_z_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xy_z_x[i] = 4.0 * g_xz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_z_y[i] = 4.0 * g_xz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_z_z[i] = 4.0 * g_xz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1314-1317)

    #pragma omp simd aligned(g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_xz_xzz_x_x, g_xz_xzz_x_y, g_xz_xzz_x_z, g_z_z_0_0_x_xz_x_x, g_z_z_0_0_x_xz_x_y, g_z_z_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xz_x_x[i] = -2.0 * g_xz_x_x_x[i] * a_exp + 4.0 * g_xz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_x_y[i] = -2.0 * g_xz_x_x_y[i] * a_exp + 4.0 * g_xz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_x_z[i] = -2.0 * g_xz_x_x_z[i] * a_exp + 4.0 * g_xz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1317-1320)

    #pragma omp simd aligned(g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_xz_xzz_y_x, g_xz_xzz_y_y, g_xz_xzz_y_z, g_z_z_0_0_x_xz_y_x, g_z_z_0_0_x_xz_y_y, g_z_z_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xz_y_x[i] = -2.0 * g_xz_x_y_x[i] * a_exp + 4.0 * g_xz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_y_y[i] = -2.0 * g_xz_x_y_y[i] * a_exp + 4.0 * g_xz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_y_z[i] = -2.0 * g_xz_x_y_z[i] * a_exp + 4.0 * g_xz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1320-1323)

    #pragma omp simd aligned(g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_xz_xzz_z_x, g_xz_xzz_z_y, g_xz_xzz_z_z, g_z_z_0_0_x_xz_z_x, g_z_z_0_0_x_xz_z_y, g_z_z_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xz_z_x[i] = -2.0 * g_xz_x_z_x[i] * a_exp + 4.0 * g_xz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_z_y[i] = -2.0 * g_xz_x_z_y[i] * a_exp + 4.0 * g_xz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_z_z[i] = -2.0 * g_xz_x_z_z[i] * a_exp + 4.0 * g_xz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1323-1326)

    #pragma omp simd aligned(g_xz_yyz_x_x, g_xz_yyz_x_y, g_xz_yyz_x_z, g_z_z_0_0_x_yy_x_x, g_z_z_0_0_x_yy_x_y, g_z_z_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yy_x_x[i] = 4.0 * g_xz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_x_y[i] = 4.0 * g_xz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_x_z[i] = 4.0 * g_xz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1326-1329)

    #pragma omp simd aligned(g_xz_yyz_y_x, g_xz_yyz_y_y, g_xz_yyz_y_z, g_z_z_0_0_x_yy_y_x, g_z_z_0_0_x_yy_y_y, g_z_z_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yy_y_x[i] = 4.0 * g_xz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_y_y[i] = 4.0 * g_xz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_y_z[i] = 4.0 * g_xz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1329-1332)

    #pragma omp simd aligned(g_xz_yyz_z_x, g_xz_yyz_z_y, g_xz_yyz_z_z, g_z_z_0_0_x_yy_z_x, g_z_z_0_0_x_yy_z_y, g_z_z_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yy_z_x[i] = 4.0 * g_xz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_z_y[i] = 4.0 * g_xz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_z_z[i] = 4.0 * g_xz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1332-1335)

    #pragma omp simd aligned(g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_xz_yzz_x_x, g_xz_yzz_x_y, g_xz_yzz_x_z, g_z_z_0_0_x_yz_x_x, g_z_z_0_0_x_yz_x_y, g_z_z_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yz_x_x[i] = -2.0 * g_xz_y_x_x[i] * a_exp + 4.0 * g_xz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_x_y[i] = -2.0 * g_xz_y_x_y[i] * a_exp + 4.0 * g_xz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_x_z[i] = -2.0 * g_xz_y_x_z[i] * a_exp + 4.0 * g_xz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1335-1338)

    #pragma omp simd aligned(g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_xz_yzz_y_x, g_xz_yzz_y_y, g_xz_yzz_y_z, g_z_z_0_0_x_yz_y_x, g_z_z_0_0_x_yz_y_y, g_z_z_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yz_y_x[i] = -2.0 * g_xz_y_y_x[i] * a_exp + 4.0 * g_xz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_y_y[i] = -2.0 * g_xz_y_y_y[i] * a_exp + 4.0 * g_xz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_y_z[i] = -2.0 * g_xz_y_y_z[i] * a_exp + 4.0 * g_xz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1338-1341)

    #pragma omp simd aligned(g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_xz_yzz_z_x, g_xz_yzz_z_y, g_xz_yzz_z_z, g_z_z_0_0_x_yz_z_x, g_z_z_0_0_x_yz_z_y, g_z_z_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yz_z_x[i] = -2.0 * g_xz_y_z_x[i] * a_exp + 4.0 * g_xz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_z_y[i] = -2.0 * g_xz_y_z_y[i] * a_exp + 4.0 * g_xz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_z_z[i] = -2.0 * g_xz_y_z_z[i] * a_exp + 4.0 * g_xz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1341-1344)

    #pragma omp simd aligned(g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z, g_xz_zzz_x_x, g_xz_zzz_x_y, g_xz_zzz_x_z, g_z_z_0_0_x_zz_x_x, g_z_z_0_0_x_zz_x_y, g_z_z_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_zz_x_x[i] = -4.0 * g_xz_z_x_x[i] * a_exp + 4.0 * g_xz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_x_y[i] = -4.0 * g_xz_z_x_y[i] * a_exp + 4.0 * g_xz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_x_z[i] = -4.0 * g_xz_z_x_z[i] * a_exp + 4.0 * g_xz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1344-1347)

    #pragma omp simd aligned(g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z, g_xz_zzz_y_x, g_xz_zzz_y_y, g_xz_zzz_y_z, g_z_z_0_0_x_zz_y_x, g_z_z_0_0_x_zz_y_y, g_z_z_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_zz_y_x[i] = -4.0 * g_xz_z_y_x[i] * a_exp + 4.0 * g_xz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_y_y[i] = -4.0 * g_xz_z_y_y[i] * a_exp + 4.0 * g_xz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_y_z[i] = -4.0 * g_xz_z_y_z[i] * a_exp + 4.0 * g_xz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1347-1350)

    #pragma omp simd aligned(g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z, g_xz_zzz_z_x, g_xz_zzz_z_y, g_xz_zzz_z_z, g_z_z_0_0_x_zz_z_x, g_z_z_0_0_x_zz_z_y, g_z_z_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_zz_z_x[i] = -4.0 * g_xz_z_z_x[i] * a_exp + 4.0 * g_xz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_z_y[i] = -4.0 * g_xz_z_z_y[i] * a_exp + 4.0 * g_xz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_z_z[i] = -4.0 * g_xz_z_z_z[i] * a_exp + 4.0 * g_xz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1350-1353)

    #pragma omp simd aligned(g_yz_xxz_x_x, g_yz_xxz_x_y, g_yz_xxz_x_z, g_z_z_0_0_y_xx_x_x, g_z_z_0_0_y_xx_x_y, g_z_z_0_0_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_x_x[i] = 4.0 * g_yz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_x_y[i] = 4.0 * g_yz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_x_z[i] = 4.0 * g_yz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1353-1356)

    #pragma omp simd aligned(g_yz_xxz_y_x, g_yz_xxz_y_y, g_yz_xxz_y_z, g_z_z_0_0_y_xx_y_x, g_z_z_0_0_y_xx_y_y, g_z_z_0_0_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_y_x[i] = 4.0 * g_yz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_y_y[i] = 4.0 * g_yz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_y_z[i] = 4.0 * g_yz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1356-1359)

    #pragma omp simd aligned(g_yz_xxz_z_x, g_yz_xxz_z_y, g_yz_xxz_z_z, g_z_z_0_0_y_xx_z_x, g_z_z_0_0_y_xx_z_y, g_z_z_0_0_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_z_x[i] = 4.0 * g_yz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_z_y[i] = 4.0 * g_yz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_z_z[i] = 4.0 * g_yz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1359-1362)

    #pragma omp simd aligned(g_yz_xyz_x_x, g_yz_xyz_x_y, g_yz_xyz_x_z, g_z_z_0_0_y_xy_x_x, g_z_z_0_0_y_xy_x_y, g_z_z_0_0_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xy_x_x[i] = 4.0 * g_yz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_x_y[i] = 4.0 * g_yz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_x_z[i] = 4.0 * g_yz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1362-1365)

    #pragma omp simd aligned(g_yz_xyz_y_x, g_yz_xyz_y_y, g_yz_xyz_y_z, g_z_z_0_0_y_xy_y_x, g_z_z_0_0_y_xy_y_y, g_z_z_0_0_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xy_y_x[i] = 4.0 * g_yz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_y_y[i] = 4.0 * g_yz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_y_z[i] = 4.0 * g_yz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1365-1368)

    #pragma omp simd aligned(g_yz_xyz_z_x, g_yz_xyz_z_y, g_yz_xyz_z_z, g_z_z_0_0_y_xy_z_x, g_z_z_0_0_y_xy_z_y, g_z_z_0_0_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xy_z_x[i] = 4.0 * g_yz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_z_y[i] = 4.0 * g_yz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_z_z[i] = 4.0 * g_yz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1368-1371)

    #pragma omp simd aligned(g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_yz_xzz_x_x, g_yz_xzz_x_y, g_yz_xzz_x_z, g_z_z_0_0_y_xz_x_x, g_z_z_0_0_y_xz_x_y, g_z_z_0_0_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xz_x_x[i] = -2.0 * g_yz_x_x_x[i] * a_exp + 4.0 * g_yz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_x_y[i] = -2.0 * g_yz_x_x_y[i] * a_exp + 4.0 * g_yz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_x_z[i] = -2.0 * g_yz_x_x_z[i] * a_exp + 4.0 * g_yz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1371-1374)

    #pragma omp simd aligned(g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_yz_xzz_y_x, g_yz_xzz_y_y, g_yz_xzz_y_z, g_z_z_0_0_y_xz_y_x, g_z_z_0_0_y_xz_y_y, g_z_z_0_0_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xz_y_x[i] = -2.0 * g_yz_x_y_x[i] * a_exp + 4.0 * g_yz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_y_y[i] = -2.0 * g_yz_x_y_y[i] * a_exp + 4.0 * g_yz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_y_z[i] = -2.0 * g_yz_x_y_z[i] * a_exp + 4.0 * g_yz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1374-1377)

    #pragma omp simd aligned(g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_yz_xzz_z_x, g_yz_xzz_z_y, g_yz_xzz_z_z, g_z_z_0_0_y_xz_z_x, g_z_z_0_0_y_xz_z_y, g_z_z_0_0_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xz_z_x[i] = -2.0 * g_yz_x_z_x[i] * a_exp + 4.0 * g_yz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_z_y[i] = -2.0 * g_yz_x_z_y[i] * a_exp + 4.0 * g_yz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_z_z[i] = -2.0 * g_yz_x_z_z[i] * a_exp + 4.0 * g_yz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1377-1380)

    #pragma omp simd aligned(g_yz_yyz_x_x, g_yz_yyz_x_y, g_yz_yyz_x_z, g_z_z_0_0_y_yy_x_x, g_z_z_0_0_y_yy_x_y, g_z_z_0_0_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yy_x_x[i] = 4.0 * g_yz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_x_y[i] = 4.0 * g_yz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_x_z[i] = 4.0 * g_yz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1380-1383)

    #pragma omp simd aligned(g_yz_yyz_y_x, g_yz_yyz_y_y, g_yz_yyz_y_z, g_z_z_0_0_y_yy_y_x, g_z_z_0_0_y_yy_y_y, g_z_z_0_0_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yy_y_x[i] = 4.0 * g_yz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_y_y[i] = 4.0 * g_yz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_y_z[i] = 4.0 * g_yz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1383-1386)

    #pragma omp simd aligned(g_yz_yyz_z_x, g_yz_yyz_z_y, g_yz_yyz_z_z, g_z_z_0_0_y_yy_z_x, g_z_z_0_0_y_yy_z_y, g_z_z_0_0_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yy_z_x[i] = 4.0 * g_yz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_z_y[i] = 4.0 * g_yz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_z_z[i] = 4.0 * g_yz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1386-1389)

    #pragma omp simd aligned(g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_yz_yzz_x_x, g_yz_yzz_x_y, g_yz_yzz_x_z, g_z_z_0_0_y_yz_x_x, g_z_z_0_0_y_yz_x_y, g_z_z_0_0_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yz_x_x[i] = -2.0 * g_yz_y_x_x[i] * a_exp + 4.0 * g_yz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_x_y[i] = -2.0 * g_yz_y_x_y[i] * a_exp + 4.0 * g_yz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_x_z[i] = -2.0 * g_yz_y_x_z[i] * a_exp + 4.0 * g_yz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1389-1392)

    #pragma omp simd aligned(g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_yz_yzz_y_x, g_yz_yzz_y_y, g_yz_yzz_y_z, g_z_z_0_0_y_yz_y_x, g_z_z_0_0_y_yz_y_y, g_z_z_0_0_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yz_y_x[i] = -2.0 * g_yz_y_y_x[i] * a_exp + 4.0 * g_yz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_y_y[i] = -2.0 * g_yz_y_y_y[i] * a_exp + 4.0 * g_yz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_y_z[i] = -2.0 * g_yz_y_y_z[i] * a_exp + 4.0 * g_yz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1392-1395)

    #pragma omp simd aligned(g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_yz_yzz_z_x, g_yz_yzz_z_y, g_yz_yzz_z_z, g_z_z_0_0_y_yz_z_x, g_z_z_0_0_y_yz_z_y, g_z_z_0_0_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yz_z_x[i] = -2.0 * g_yz_y_z_x[i] * a_exp + 4.0 * g_yz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_z_y[i] = -2.0 * g_yz_y_z_y[i] * a_exp + 4.0 * g_yz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_z_z[i] = -2.0 * g_yz_y_z_z[i] * a_exp + 4.0 * g_yz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1395-1398)

    #pragma omp simd aligned(g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z, g_yz_zzz_x_x, g_yz_zzz_x_y, g_yz_zzz_x_z, g_z_z_0_0_y_zz_x_x, g_z_z_0_0_y_zz_x_y, g_z_z_0_0_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_zz_x_x[i] = -4.0 * g_yz_z_x_x[i] * a_exp + 4.0 * g_yz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_x_y[i] = -4.0 * g_yz_z_x_y[i] * a_exp + 4.0 * g_yz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_x_z[i] = -4.0 * g_yz_z_x_z[i] * a_exp + 4.0 * g_yz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1398-1401)

    #pragma omp simd aligned(g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z, g_yz_zzz_y_x, g_yz_zzz_y_y, g_yz_zzz_y_z, g_z_z_0_0_y_zz_y_x, g_z_z_0_0_y_zz_y_y, g_z_z_0_0_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_zz_y_x[i] = -4.0 * g_yz_z_y_x[i] * a_exp + 4.0 * g_yz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_y_y[i] = -4.0 * g_yz_z_y_y[i] * a_exp + 4.0 * g_yz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_y_z[i] = -4.0 * g_yz_z_y_z[i] * a_exp + 4.0 * g_yz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1401-1404)

    #pragma omp simd aligned(g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z, g_yz_zzz_z_x, g_yz_zzz_z_y, g_yz_zzz_z_z, g_z_z_0_0_y_zz_z_x, g_z_z_0_0_y_zz_z_y, g_z_z_0_0_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_zz_z_x[i] = -4.0 * g_yz_z_z_x[i] * a_exp + 4.0 * g_yz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_z_y[i] = -4.0 * g_yz_z_z_y[i] * a_exp + 4.0 * g_yz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_z_z[i] = -4.0 * g_yz_z_z_z[i] * a_exp + 4.0 * g_yz_zzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1404-1407)

    #pragma omp simd aligned(g_0_xxz_x_x, g_0_xxz_x_y, g_0_xxz_x_z, g_z_z_0_0_z_xx_x_x, g_z_z_0_0_z_xx_x_y, g_z_z_0_0_z_xx_x_z, g_zz_xxz_x_x, g_zz_xxz_x_y, g_zz_xxz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_x_x[i] = -2.0 * g_0_xxz_x_x[i] * b_exp + 4.0 * g_zz_xxz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_x_y[i] = -2.0 * g_0_xxz_x_y[i] * b_exp + 4.0 * g_zz_xxz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_x_z[i] = -2.0 * g_0_xxz_x_z[i] * b_exp + 4.0 * g_zz_xxz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1407-1410)

    #pragma omp simd aligned(g_0_xxz_y_x, g_0_xxz_y_y, g_0_xxz_y_z, g_z_z_0_0_z_xx_y_x, g_z_z_0_0_z_xx_y_y, g_z_z_0_0_z_xx_y_z, g_zz_xxz_y_x, g_zz_xxz_y_y, g_zz_xxz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_y_x[i] = -2.0 * g_0_xxz_y_x[i] * b_exp + 4.0 * g_zz_xxz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_y_y[i] = -2.0 * g_0_xxz_y_y[i] * b_exp + 4.0 * g_zz_xxz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_y_z[i] = -2.0 * g_0_xxz_y_z[i] * b_exp + 4.0 * g_zz_xxz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1410-1413)

    #pragma omp simd aligned(g_0_xxz_z_x, g_0_xxz_z_y, g_0_xxz_z_z, g_z_z_0_0_z_xx_z_x, g_z_z_0_0_z_xx_z_y, g_z_z_0_0_z_xx_z_z, g_zz_xxz_z_x, g_zz_xxz_z_y, g_zz_xxz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_z_x[i] = -2.0 * g_0_xxz_z_x[i] * b_exp + 4.0 * g_zz_xxz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_z_y[i] = -2.0 * g_0_xxz_z_y[i] * b_exp + 4.0 * g_zz_xxz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_z_z[i] = -2.0 * g_0_xxz_z_z[i] * b_exp + 4.0 * g_zz_xxz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1413-1416)

    #pragma omp simd aligned(g_0_xyz_x_x, g_0_xyz_x_y, g_0_xyz_x_z, g_z_z_0_0_z_xy_x_x, g_z_z_0_0_z_xy_x_y, g_z_z_0_0_z_xy_x_z, g_zz_xyz_x_x, g_zz_xyz_x_y, g_zz_xyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xy_x_x[i] = -2.0 * g_0_xyz_x_x[i] * b_exp + 4.0 * g_zz_xyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_x_y[i] = -2.0 * g_0_xyz_x_y[i] * b_exp + 4.0 * g_zz_xyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_x_z[i] = -2.0 * g_0_xyz_x_z[i] * b_exp + 4.0 * g_zz_xyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1416-1419)

    #pragma omp simd aligned(g_0_xyz_y_x, g_0_xyz_y_y, g_0_xyz_y_z, g_z_z_0_0_z_xy_y_x, g_z_z_0_0_z_xy_y_y, g_z_z_0_0_z_xy_y_z, g_zz_xyz_y_x, g_zz_xyz_y_y, g_zz_xyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xy_y_x[i] = -2.0 * g_0_xyz_y_x[i] * b_exp + 4.0 * g_zz_xyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_y_y[i] = -2.0 * g_0_xyz_y_y[i] * b_exp + 4.0 * g_zz_xyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_y_z[i] = -2.0 * g_0_xyz_y_z[i] * b_exp + 4.0 * g_zz_xyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1419-1422)

    #pragma omp simd aligned(g_0_xyz_z_x, g_0_xyz_z_y, g_0_xyz_z_z, g_z_z_0_0_z_xy_z_x, g_z_z_0_0_z_xy_z_y, g_z_z_0_0_z_xy_z_z, g_zz_xyz_z_x, g_zz_xyz_z_y, g_zz_xyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xy_z_x[i] = -2.0 * g_0_xyz_z_x[i] * b_exp + 4.0 * g_zz_xyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_z_y[i] = -2.0 * g_0_xyz_z_y[i] * b_exp + 4.0 * g_zz_xyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_z_z[i] = -2.0 * g_0_xyz_z_z[i] * b_exp + 4.0 * g_zz_xyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1422-1425)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_0_xzz_x_x, g_0_xzz_x_y, g_0_xzz_x_z, g_z_z_0_0_z_xz_x_x, g_z_z_0_0_z_xz_x_y, g_z_z_0_0_z_xz_x_z, g_zz_x_x_x, g_zz_x_x_y, g_zz_x_x_z, g_zz_xzz_x_x, g_zz_xzz_x_y, g_zz_xzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xz_x_x[i] = g_0_x_x_x[i] - 2.0 * g_0_xzz_x_x[i] * b_exp - 2.0 * g_zz_x_x_x[i] * a_exp + 4.0 * g_zz_xzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_x_y[i] = g_0_x_x_y[i] - 2.0 * g_0_xzz_x_y[i] * b_exp - 2.0 * g_zz_x_x_y[i] * a_exp + 4.0 * g_zz_xzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_x_z[i] = g_0_x_x_z[i] - 2.0 * g_0_xzz_x_z[i] * b_exp - 2.0 * g_zz_x_x_z[i] * a_exp + 4.0 * g_zz_xzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1425-1428)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_0_xzz_y_x, g_0_xzz_y_y, g_0_xzz_y_z, g_z_z_0_0_z_xz_y_x, g_z_z_0_0_z_xz_y_y, g_z_z_0_0_z_xz_y_z, g_zz_x_y_x, g_zz_x_y_y, g_zz_x_y_z, g_zz_xzz_y_x, g_zz_xzz_y_y, g_zz_xzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xz_y_x[i] = g_0_x_y_x[i] - 2.0 * g_0_xzz_y_x[i] * b_exp - 2.0 * g_zz_x_y_x[i] * a_exp + 4.0 * g_zz_xzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_y_y[i] = g_0_x_y_y[i] - 2.0 * g_0_xzz_y_y[i] * b_exp - 2.0 * g_zz_x_y_y[i] * a_exp + 4.0 * g_zz_xzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_y_z[i] = g_0_x_y_z[i] - 2.0 * g_0_xzz_y_z[i] * b_exp - 2.0 * g_zz_x_y_z[i] * a_exp + 4.0 * g_zz_xzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1428-1431)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_0_xzz_z_x, g_0_xzz_z_y, g_0_xzz_z_z, g_z_z_0_0_z_xz_z_x, g_z_z_0_0_z_xz_z_y, g_z_z_0_0_z_xz_z_z, g_zz_x_z_x, g_zz_x_z_y, g_zz_x_z_z, g_zz_xzz_z_x, g_zz_xzz_z_y, g_zz_xzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xz_z_x[i] = g_0_x_z_x[i] - 2.0 * g_0_xzz_z_x[i] * b_exp - 2.0 * g_zz_x_z_x[i] * a_exp + 4.0 * g_zz_xzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_z_y[i] = g_0_x_z_y[i] - 2.0 * g_0_xzz_z_y[i] * b_exp - 2.0 * g_zz_x_z_y[i] * a_exp + 4.0 * g_zz_xzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_z_z[i] = g_0_x_z_z[i] - 2.0 * g_0_xzz_z_z[i] * b_exp - 2.0 * g_zz_x_z_z[i] * a_exp + 4.0 * g_zz_xzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1431-1434)

    #pragma omp simd aligned(g_0_yyz_x_x, g_0_yyz_x_y, g_0_yyz_x_z, g_z_z_0_0_z_yy_x_x, g_z_z_0_0_z_yy_x_y, g_z_z_0_0_z_yy_x_z, g_zz_yyz_x_x, g_zz_yyz_x_y, g_zz_yyz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yy_x_x[i] = -2.0 * g_0_yyz_x_x[i] * b_exp + 4.0 * g_zz_yyz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_x_y[i] = -2.0 * g_0_yyz_x_y[i] * b_exp + 4.0 * g_zz_yyz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_x_z[i] = -2.0 * g_0_yyz_x_z[i] * b_exp + 4.0 * g_zz_yyz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1434-1437)

    #pragma omp simd aligned(g_0_yyz_y_x, g_0_yyz_y_y, g_0_yyz_y_z, g_z_z_0_0_z_yy_y_x, g_z_z_0_0_z_yy_y_y, g_z_z_0_0_z_yy_y_z, g_zz_yyz_y_x, g_zz_yyz_y_y, g_zz_yyz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yy_y_x[i] = -2.0 * g_0_yyz_y_x[i] * b_exp + 4.0 * g_zz_yyz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_y_y[i] = -2.0 * g_0_yyz_y_y[i] * b_exp + 4.0 * g_zz_yyz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_y_z[i] = -2.0 * g_0_yyz_y_z[i] * b_exp + 4.0 * g_zz_yyz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1437-1440)

    #pragma omp simd aligned(g_0_yyz_z_x, g_0_yyz_z_y, g_0_yyz_z_z, g_z_z_0_0_z_yy_z_x, g_z_z_0_0_z_yy_z_y, g_z_z_0_0_z_yy_z_z, g_zz_yyz_z_x, g_zz_yyz_z_y, g_zz_yyz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yy_z_x[i] = -2.0 * g_0_yyz_z_x[i] * b_exp + 4.0 * g_zz_yyz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_z_y[i] = -2.0 * g_0_yyz_z_y[i] * b_exp + 4.0 * g_zz_yyz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_z_z[i] = -2.0 * g_0_yyz_z_z[i] * b_exp + 4.0 * g_zz_yyz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1440-1443)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_0_yzz_x_x, g_0_yzz_x_y, g_0_yzz_x_z, g_z_z_0_0_z_yz_x_x, g_z_z_0_0_z_yz_x_y, g_z_z_0_0_z_yz_x_z, g_zz_y_x_x, g_zz_y_x_y, g_zz_y_x_z, g_zz_yzz_x_x, g_zz_yzz_x_y, g_zz_yzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yz_x_x[i] = g_0_y_x_x[i] - 2.0 * g_0_yzz_x_x[i] * b_exp - 2.0 * g_zz_y_x_x[i] * a_exp + 4.0 * g_zz_yzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_x_y[i] = g_0_y_x_y[i] - 2.0 * g_0_yzz_x_y[i] * b_exp - 2.0 * g_zz_y_x_y[i] * a_exp + 4.0 * g_zz_yzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_x_z[i] = g_0_y_x_z[i] - 2.0 * g_0_yzz_x_z[i] * b_exp - 2.0 * g_zz_y_x_z[i] * a_exp + 4.0 * g_zz_yzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1443-1446)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_0_yzz_y_x, g_0_yzz_y_y, g_0_yzz_y_z, g_z_z_0_0_z_yz_y_x, g_z_z_0_0_z_yz_y_y, g_z_z_0_0_z_yz_y_z, g_zz_y_y_x, g_zz_y_y_y, g_zz_y_y_z, g_zz_yzz_y_x, g_zz_yzz_y_y, g_zz_yzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yz_y_x[i] = g_0_y_y_x[i] - 2.0 * g_0_yzz_y_x[i] * b_exp - 2.0 * g_zz_y_y_x[i] * a_exp + 4.0 * g_zz_yzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_y_y[i] = g_0_y_y_y[i] - 2.0 * g_0_yzz_y_y[i] * b_exp - 2.0 * g_zz_y_y_y[i] * a_exp + 4.0 * g_zz_yzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_y_z[i] = g_0_y_y_z[i] - 2.0 * g_0_yzz_y_z[i] * b_exp - 2.0 * g_zz_y_y_z[i] * a_exp + 4.0 * g_zz_yzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1446-1449)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_0_yzz_z_x, g_0_yzz_z_y, g_0_yzz_z_z, g_z_z_0_0_z_yz_z_x, g_z_z_0_0_z_yz_z_y, g_z_z_0_0_z_yz_z_z, g_zz_y_z_x, g_zz_y_z_y, g_zz_y_z_z, g_zz_yzz_z_x, g_zz_yzz_z_y, g_zz_yzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yz_z_x[i] = g_0_y_z_x[i] - 2.0 * g_0_yzz_z_x[i] * b_exp - 2.0 * g_zz_y_z_x[i] * a_exp + 4.0 * g_zz_yzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_z_y[i] = g_0_y_z_y[i] - 2.0 * g_0_yzz_z_y[i] * b_exp - 2.0 * g_zz_y_z_y[i] * a_exp + 4.0 * g_zz_yzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_z_z[i] = g_0_y_z_z[i] - 2.0 * g_0_yzz_z_z[i] * b_exp - 2.0 * g_zz_y_z_z[i] * a_exp + 4.0 * g_zz_yzz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (1449-1452)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_0_zzz_x_x, g_0_zzz_x_y, g_0_zzz_x_z, g_z_z_0_0_z_zz_x_x, g_z_z_0_0_z_zz_x_y, g_z_z_0_0_z_zz_x_z, g_zz_z_x_x, g_zz_z_x_y, g_zz_z_x_z, g_zz_zzz_x_x, g_zz_zzz_x_y, g_zz_zzz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_zz_x_x[i] = 2.0 * g_0_z_x_x[i] - 2.0 * g_0_zzz_x_x[i] * b_exp - 4.0 * g_zz_z_x_x[i] * a_exp + 4.0 * g_zz_zzz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_x_y[i] = 2.0 * g_0_z_x_y[i] - 2.0 * g_0_zzz_x_y[i] * b_exp - 4.0 * g_zz_z_x_y[i] * a_exp + 4.0 * g_zz_zzz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_x_z[i] = 2.0 * g_0_z_x_z[i] - 2.0 * g_0_zzz_x_z[i] * b_exp - 4.0 * g_zz_z_x_z[i] * a_exp + 4.0 * g_zz_zzz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (1452-1455)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_0_zzz_y_x, g_0_zzz_y_y, g_0_zzz_y_z, g_z_z_0_0_z_zz_y_x, g_z_z_0_0_z_zz_y_y, g_z_z_0_0_z_zz_y_z, g_zz_z_y_x, g_zz_z_y_y, g_zz_z_y_z, g_zz_zzz_y_x, g_zz_zzz_y_y, g_zz_zzz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_zz_y_x[i] = 2.0 * g_0_z_y_x[i] - 2.0 * g_0_zzz_y_x[i] * b_exp - 4.0 * g_zz_z_y_x[i] * a_exp + 4.0 * g_zz_zzz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_y_y[i] = 2.0 * g_0_z_y_y[i] - 2.0 * g_0_zzz_y_y[i] * b_exp - 4.0 * g_zz_z_y_y[i] * a_exp + 4.0 * g_zz_zzz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_y_z[i] = 2.0 * g_0_z_y_z[i] - 2.0 * g_0_zzz_y_z[i] * b_exp - 4.0 * g_zz_z_y_z[i] * a_exp + 4.0 * g_zz_zzz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (1455-1458)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_0_zzz_z_x, g_0_zzz_z_y, g_0_zzz_z_z, g_z_z_0_0_z_zz_z_x, g_z_z_0_0_z_zz_z_y, g_z_z_0_0_z_zz_z_z, g_zz_z_z_x, g_zz_z_z_y, g_zz_z_z_z, g_zz_zzz_z_x, g_zz_zzz_z_y, g_zz_zzz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_zz_z_x[i] = 2.0 * g_0_z_z_x[i] - 2.0 * g_0_zzz_z_x[i] * b_exp - 4.0 * g_zz_z_z_x[i] * a_exp + 4.0 * g_zz_zzz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_z_y[i] = 2.0 * g_0_z_z_y[i] - 2.0 * g_0_zzz_z_y[i] * b_exp - 4.0 * g_zz_z_z_y[i] * a_exp + 4.0 * g_zz_zzz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_z_z[i] = 2.0 * g_0_z_z_z[i] - 2.0 * g_0_zzz_z_z[i] * b_exp - 4.0 * g_zz_z_z_z[i] * a_exp + 4.0 * g_zz_zzz_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


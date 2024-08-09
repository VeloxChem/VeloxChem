#include "GeomDeriv1100OfScalarForSDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sdsd_0(CSimdArray<double>& buffer_1100_sdsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_pfsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sdsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsd

    auto g_x_x_0_xx = buffer_ppsd[0];

    auto g_x_x_0_xy = buffer_ppsd[1];

    auto g_x_x_0_xz = buffer_ppsd[2];

    auto g_x_x_0_yy = buffer_ppsd[3];

    auto g_x_x_0_yz = buffer_ppsd[4];

    auto g_x_x_0_zz = buffer_ppsd[5];

    auto g_x_y_0_xx = buffer_ppsd[6];

    auto g_x_y_0_xy = buffer_ppsd[7];

    auto g_x_y_0_xz = buffer_ppsd[8];

    auto g_x_y_0_yy = buffer_ppsd[9];

    auto g_x_y_0_yz = buffer_ppsd[10];

    auto g_x_y_0_zz = buffer_ppsd[11];

    auto g_x_z_0_xx = buffer_ppsd[12];

    auto g_x_z_0_xy = buffer_ppsd[13];

    auto g_x_z_0_xz = buffer_ppsd[14];

    auto g_x_z_0_yy = buffer_ppsd[15];

    auto g_x_z_0_yz = buffer_ppsd[16];

    auto g_x_z_0_zz = buffer_ppsd[17];

    auto g_y_x_0_xx = buffer_ppsd[18];

    auto g_y_x_0_xy = buffer_ppsd[19];

    auto g_y_x_0_xz = buffer_ppsd[20];

    auto g_y_x_0_yy = buffer_ppsd[21];

    auto g_y_x_0_yz = buffer_ppsd[22];

    auto g_y_x_0_zz = buffer_ppsd[23];

    auto g_y_y_0_xx = buffer_ppsd[24];

    auto g_y_y_0_xy = buffer_ppsd[25];

    auto g_y_y_0_xz = buffer_ppsd[26];

    auto g_y_y_0_yy = buffer_ppsd[27];

    auto g_y_y_0_yz = buffer_ppsd[28];

    auto g_y_y_0_zz = buffer_ppsd[29];

    auto g_y_z_0_xx = buffer_ppsd[30];

    auto g_y_z_0_xy = buffer_ppsd[31];

    auto g_y_z_0_xz = buffer_ppsd[32];

    auto g_y_z_0_yy = buffer_ppsd[33];

    auto g_y_z_0_yz = buffer_ppsd[34];

    auto g_y_z_0_zz = buffer_ppsd[35];

    auto g_z_x_0_xx = buffer_ppsd[36];

    auto g_z_x_0_xy = buffer_ppsd[37];

    auto g_z_x_0_xz = buffer_ppsd[38];

    auto g_z_x_0_yy = buffer_ppsd[39];

    auto g_z_x_0_yz = buffer_ppsd[40];

    auto g_z_x_0_zz = buffer_ppsd[41];

    auto g_z_y_0_xx = buffer_ppsd[42];

    auto g_z_y_0_xy = buffer_ppsd[43];

    auto g_z_y_0_xz = buffer_ppsd[44];

    auto g_z_y_0_yy = buffer_ppsd[45];

    auto g_z_y_0_yz = buffer_ppsd[46];

    auto g_z_y_0_zz = buffer_ppsd[47];

    auto g_z_z_0_xx = buffer_ppsd[48];

    auto g_z_z_0_xy = buffer_ppsd[49];

    auto g_z_z_0_xz = buffer_ppsd[50];

    auto g_z_z_0_yy = buffer_ppsd[51];

    auto g_z_z_0_yz = buffer_ppsd[52];

    auto g_z_z_0_zz = buffer_ppsd[53];

    /// Set up components of auxilary buffer : buffer_pfsd

    auto g_x_xxx_0_xx = buffer_pfsd[0];

    auto g_x_xxx_0_xy = buffer_pfsd[1];

    auto g_x_xxx_0_xz = buffer_pfsd[2];

    auto g_x_xxx_0_yy = buffer_pfsd[3];

    auto g_x_xxx_0_yz = buffer_pfsd[4];

    auto g_x_xxx_0_zz = buffer_pfsd[5];

    auto g_x_xxy_0_xx = buffer_pfsd[6];

    auto g_x_xxy_0_xy = buffer_pfsd[7];

    auto g_x_xxy_0_xz = buffer_pfsd[8];

    auto g_x_xxy_0_yy = buffer_pfsd[9];

    auto g_x_xxy_0_yz = buffer_pfsd[10];

    auto g_x_xxy_0_zz = buffer_pfsd[11];

    auto g_x_xxz_0_xx = buffer_pfsd[12];

    auto g_x_xxz_0_xy = buffer_pfsd[13];

    auto g_x_xxz_0_xz = buffer_pfsd[14];

    auto g_x_xxz_0_yy = buffer_pfsd[15];

    auto g_x_xxz_0_yz = buffer_pfsd[16];

    auto g_x_xxz_0_zz = buffer_pfsd[17];

    auto g_x_xyy_0_xx = buffer_pfsd[18];

    auto g_x_xyy_0_xy = buffer_pfsd[19];

    auto g_x_xyy_0_xz = buffer_pfsd[20];

    auto g_x_xyy_0_yy = buffer_pfsd[21];

    auto g_x_xyy_0_yz = buffer_pfsd[22];

    auto g_x_xyy_0_zz = buffer_pfsd[23];

    auto g_x_xyz_0_xx = buffer_pfsd[24];

    auto g_x_xyz_0_xy = buffer_pfsd[25];

    auto g_x_xyz_0_xz = buffer_pfsd[26];

    auto g_x_xyz_0_yy = buffer_pfsd[27];

    auto g_x_xyz_0_yz = buffer_pfsd[28];

    auto g_x_xyz_0_zz = buffer_pfsd[29];

    auto g_x_xzz_0_xx = buffer_pfsd[30];

    auto g_x_xzz_0_xy = buffer_pfsd[31];

    auto g_x_xzz_0_xz = buffer_pfsd[32];

    auto g_x_xzz_0_yy = buffer_pfsd[33];

    auto g_x_xzz_0_yz = buffer_pfsd[34];

    auto g_x_xzz_0_zz = buffer_pfsd[35];

    auto g_x_yyy_0_xx = buffer_pfsd[36];

    auto g_x_yyy_0_xy = buffer_pfsd[37];

    auto g_x_yyy_0_xz = buffer_pfsd[38];

    auto g_x_yyy_0_yy = buffer_pfsd[39];

    auto g_x_yyy_0_yz = buffer_pfsd[40];

    auto g_x_yyy_0_zz = buffer_pfsd[41];

    auto g_x_yyz_0_xx = buffer_pfsd[42];

    auto g_x_yyz_0_xy = buffer_pfsd[43];

    auto g_x_yyz_0_xz = buffer_pfsd[44];

    auto g_x_yyz_0_yy = buffer_pfsd[45];

    auto g_x_yyz_0_yz = buffer_pfsd[46];

    auto g_x_yyz_0_zz = buffer_pfsd[47];

    auto g_x_yzz_0_xx = buffer_pfsd[48];

    auto g_x_yzz_0_xy = buffer_pfsd[49];

    auto g_x_yzz_0_xz = buffer_pfsd[50];

    auto g_x_yzz_0_yy = buffer_pfsd[51];

    auto g_x_yzz_0_yz = buffer_pfsd[52];

    auto g_x_yzz_0_zz = buffer_pfsd[53];

    auto g_x_zzz_0_xx = buffer_pfsd[54];

    auto g_x_zzz_0_xy = buffer_pfsd[55];

    auto g_x_zzz_0_xz = buffer_pfsd[56];

    auto g_x_zzz_0_yy = buffer_pfsd[57];

    auto g_x_zzz_0_yz = buffer_pfsd[58];

    auto g_x_zzz_0_zz = buffer_pfsd[59];

    auto g_y_xxx_0_xx = buffer_pfsd[60];

    auto g_y_xxx_0_xy = buffer_pfsd[61];

    auto g_y_xxx_0_xz = buffer_pfsd[62];

    auto g_y_xxx_0_yy = buffer_pfsd[63];

    auto g_y_xxx_0_yz = buffer_pfsd[64];

    auto g_y_xxx_0_zz = buffer_pfsd[65];

    auto g_y_xxy_0_xx = buffer_pfsd[66];

    auto g_y_xxy_0_xy = buffer_pfsd[67];

    auto g_y_xxy_0_xz = buffer_pfsd[68];

    auto g_y_xxy_0_yy = buffer_pfsd[69];

    auto g_y_xxy_0_yz = buffer_pfsd[70];

    auto g_y_xxy_0_zz = buffer_pfsd[71];

    auto g_y_xxz_0_xx = buffer_pfsd[72];

    auto g_y_xxz_0_xy = buffer_pfsd[73];

    auto g_y_xxz_0_xz = buffer_pfsd[74];

    auto g_y_xxz_0_yy = buffer_pfsd[75];

    auto g_y_xxz_0_yz = buffer_pfsd[76];

    auto g_y_xxz_0_zz = buffer_pfsd[77];

    auto g_y_xyy_0_xx = buffer_pfsd[78];

    auto g_y_xyy_0_xy = buffer_pfsd[79];

    auto g_y_xyy_0_xz = buffer_pfsd[80];

    auto g_y_xyy_0_yy = buffer_pfsd[81];

    auto g_y_xyy_0_yz = buffer_pfsd[82];

    auto g_y_xyy_0_zz = buffer_pfsd[83];

    auto g_y_xyz_0_xx = buffer_pfsd[84];

    auto g_y_xyz_0_xy = buffer_pfsd[85];

    auto g_y_xyz_0_xz = buffer_pfsd[86];

    auto g_y_xyz_0_yy = buffer_pfsd[87];

    auto g_y_xyz_0_yz = buffer_pfsd[88];

    auto g_y_xyz_0_zz = buffer_pfsd[89];

    auto g_y_xzz_0_xx = buffer_pfsd[90];

    auto g_y_xzz_0_xy = buffer_pfsd[91];

    auto g_y_xzz_0_xz = buffer_pfsd[92];

    auto g_y_xzz_0_yy = buffer_pfsd[93];

    auto g_y_xzz_0_yz = buffer_pfsd[94];

    auto g_y_xzz_0_zz = buffer_pfsd[95];

    auto g_y_yyy_0_xx = buffer_pfsd[96];

    auto g_y_yyy_0_xy = buffer_pfsd[97];

    auto g_y_yyy_0_xz = buffer_pfsd[98];

    auto g_y_yyy_0_yy = buffer_pfsd[99];

    auto g_y_yyy_0_yz = buffer_pfsd[100];

    auto g_y_yyy_0_zz = buffer_pfsd[101];

    auto g_y_yyz_0_xx = buffer_pfsd[102];

    auto g_y_yyz_0_xy = buffer_pfsd[103];

    auto g_y_yyz_0_xz = buffer_pfsd[104];

    auto g_y_yyz_0_yy = buffer_pfsd[105];

    auto g_y_yyz_0_yz = buffer_pfsd[106];

    auto g_y_yyz_0_zz = buffer_pfsd[107];

    auto g_y_yzz_0_xx = buffer_pfsd[108];

    auto g_y_yzz_0_xy = buffer_pfsd[109];

    auto g_y_yzz_0_xz = buffer_pfsd[110];

    auto g_y_yzz_0_yy = buffer_pfsd[111];

    auto g_y_yzz_0_yz = buffer_pfsd[112];

    auto g_y_yzz_0_zz = buffer_pfsd[113];

    auto g_y_zzz_0_xx = buffer_pfsd[114];

    auto g_y_zzz_0_xy = buffer_pfsd[115];

    auto g_y_zzz_0_xz = buffer_pfsd[116];

    auto g_y_zzz_0_yy = buffer_pfsd[117];

    auto g_y_zzz_0_yz = buffer_pfsd[118];

    auto g_y_zzz_0_zz = buffer_pfsd[119];

    auto g_z_xxx_0_xx = buffer_pfsd[120];

    auto g_z_xxx_0_xy = buffer_pfsd[121];

    auto g_z_xxx_0_xz = buffer_pfsd[122];

    auto g_z_xxx_0_yy = buffer_pfsd[123];

    auto g_z_xxx_0_yz = buffer_pfsd[124];

    auto g_z_xxx_0_zz = buffer_pfsd[125];

    auto g_z_xxy_0_xx = buffer_pfsd[126];

    auto g_z_xxy_0_xy = buffer_pfsd[127];

    auto g_z_xxy_0_xz = buffer_pfsd[128];

    auto g_z_xxy_0_yy = buffer_pfsd[129];

    auto g_z_xxy_0_yz = buffer_pfsd[130];

    auto g_z_xxy_0_zz = buffer_pfsd[131];

    auto g_z_xxz_0_xx = buffer_pfsd[132];

    auto g_z_xxz_0_xy = buffer_pfsd[133];

    auto g_z_xxz_0_xz = buffer_pfsd[134];

    auto g_z_xxz_0_yy = buffer_pfsd[135];

    auto g_z_xxz_0_yz = buffer_pfsd[136];

    auto g_z_xxz_0_zz = buffer_pfsd[137];

    auto g_z_xyy_0_xx = buffer_pfsd[138];

    auto g_z_xyy_0_xy = buffer_pfsd[139];

    auto g_z_xyy_0_xz = buffer_pfsd[140];

    auto g_z_xyy_0_yy = buffer_pfsd[141];

    auto g_z_xyy_0_yz = buffer_pfsd[142];

    auto g_z_xyy_0_zz = buffer_pfsd[143];

    auto g_z_xyz_0_xx = buffer_pfsd[144];

    auto g_z_xyz_0_xy = buffer_pfsd[145];

    auto g_z_xyz_0_xz = buffer_pfsd[146];

    auto g_z_xyz_0_yy = buffer_pfsd[147];

    auto g_z_xyz_0_yz = buffer_pfsd[148];

    auto g_z_xyz_0_zz = buffer_pfsd[149];

    auto g_z_xzz_0_xx = buffer_pfsd[150];

    auto g_z_xzz_0_xy = buffer_pfsd[151];

    auto g_z_xzz_0_xz = buffer_pfsd[152];

    auto g_z_xzz_0_yy = buffer_pfsd[153];

    auto g_z_xzz_0_yz = buffer_pfsd[154];

    auto g_z_xzz_0_zz = buffer_pfsd[155];

    auto g_z_yyy_0_xx = buffer_pfsd[156];

    auto g_z_yyy_0_xy = buffer_pfsd[157];

    auto g_z_yyy_0_xz = buffer_pfsd[158];

    auto g_z_yyy_0_yy = buffer_pfsd[159];

    auto g_z_yyy_0_yz = buffer_pfsd[160];

    auto g_z_yyy_0_zz = buffer_pfsd[161];

    auto g_z_yyz_0_xx = buffer_pfsd[162];

    auto g_z_yyz_0_xy = buffer_pfsd[163];

    auto g_z_yyz_0_xz = buffer_pfsd[164];

    auto g_z_yyz_0_yy = buffer_pfsd[165];

    auto g_z_yyz_0_yz = buffer_pfsd[166];

    auto g_z_yyz_0_zz = buffer_pfsd[167];

    auto g_z_yzz_0_xx = buffer_pfsd[168];

    auto g_z_yzz_0_xy = buffer_pfsd[169];

    auto g_z_yzz_0_xz = buffer_pfsd[170];

    auto g_z_yzz_0_yy = buffer_pfsd[171];

    auto g_z_yzz_0_yz = buffer_pfsd[172];

    auto g_z_yzz_0_zz = buffer_pfsd[173];

    auto g_z_zzz_0_xx = buffer_pfsd[174];

    auto g_z_zzz_0_xy = buffer_pfsd[175];

    auto g_z_zzz_0_xz = buffer_pfsd[176];

    auto g_z_zzz_0_yy = buffer_pfsd[177];

    auto g_z_zzz_0_yz = buffer_pfsd[178];

    auto g_z_zzz_0_zz = buffer_pfsd[179];

    /// Set up components of integrals buffer : buffer_1100_sdsd

    auto g_x_x_0_0_0_xx_0_xx = buffer_1100_sdsd[0];

    auto g_x_x_0_0_0_xx_0_xy = buffer_1100_sdsd[1];

    auto g_x_x_0_0_0_xx_0_xz = buffer_1100_sdsd[2];

    auto g_x_x_0_0_0_xx_0_yy = buffer_1100_sdsd[3];

    auto g_x_x_0_0_0_xx_0_yz = buffer_1100_sdsd[4];

    auto g_x_x_0_0_0_xx_0_zz = buffer_1100_sdsd[5];

    auto g_x_x_0_0_0_xy_0_xx = buffer_1100_sdsd[6];

    auto g_x_x_0_0_0_xy_0_xy = buffer_1100_sdsd[7];

    auto g_x_x_0_0_0_xy_0_xz = buffer_1100_sdsd[8];

    auto g_x_x_0_0_0_xy_0_yy = buffer_1100_sdsd[9];

    auto g_x_x_0_0_0_xy_0_yz = buffer_1100_sdsd[10];

    auto g_x_x_0_0_0_xy_0_zz = buffer_1100_sdsd[11];

    auto g_x_x_0_0_0_xz_0_xx = buffer_1100_sdsd[12];

    auto g_x_x_0_0_0_xz_0_xy = buffer_1100_sdsd[13];

    auto g_x_x_0_0_0_xz_0_xz = buffer_1100_sdsd[14];

    auto g_x_x_0_0_0_xz_0_yy = buffer_1100_sdsd[15];

    auto g_x_x_0_0_0_xz_0_yz = buffer_1100_sdsd[16];

    auto g_x_x_0_0_0_xz_0_zz = buffer_1100_sdsd[17];

    auto g_x_x_0_0_0_yy_0_xx = buffer_1100_sdsd[18];

    auto g_x_x_0_0_0_yy_0_xy = buffer_1100_sdsd[19];

    auto g_x_x_0_0_0_yy_0_xz = buffer_1100_sdsd[20];

    auto g_x_x_0_0_0_yy_0_yy = buffer_1100_sdsd[21];

    auto g_x_x_0_0_0_yy_0_yz = buffer_1100_sdsd[22];

    auto g_x_x_0_0_0_yy_0_zz = buffer_1100_sdsd[23];

    auto g_x_x_0_0_0_yz_0_xx = buffer_1100_sdsd[24];

    auto g_x_x_0_0_0_yz_0_xy = buffer_1100_sdsd[25];

    auto g_x_x_0_0_0_yz_0_xz = buffer_1100_sdsd[26];

    auto g_x_x_0_0_0_yz_0_yy = buffer_1100_sdsd[27];

    auto g_x_x_0_0_0_yz_0_yz = buffer_1100_sdsd[28];

    auto g_x_x_0_0_0_yz_0_zz = buffer_1100_sdsd[29];

    auto g_x_x_0_0_0_zz_0_xx = buffer_1100_sdsd[30];

    auto g_x_x_0_0_0_zz_0_xy = buffer_1100_sdsd[31];

    auto g_x_x_0_0_0_zz_0_xz = buffer_1100_sdsd[32];

    auto g_x_x_0_0_0_zz_0_yy = buffer_1100_sdsd[33];

    auto g_x_x_0_0_0_zz_0_yz = buffer_1100_sdsd[34];

    auto g_x_x_0_0_0_zz_0_zz = buffer_1100_sdsd[35];

    auto g_x_y_0_0_0_xx_0_xx = buffer_1100_sdsd[36];

    auto g_x_y_0_0_0_xx_0_xy = buffer_1100_sdsd[37];

    auto g_x_y_0_0_0_xx_0_xz = buffer_1100_sdsd[38];

    auto g_x_y_0_0_0_xx_0_yy = buffer_1100_sdsd[39];

    auto g_x_y_0_0_0_xx_0_yz = buffer_1100_sdsd[40];

    auto g_x_y_0_0_0_xx_0_zz = buffer_1100_sdsd[41];

    auto g_x_y_0_0_0_xy_0_xx = buffer_1100_sdsd[42];

    auto g_x_y_0_0_0_xy_0_xy = buffer_1100_sdsd[43];

    auto g_x_y_0_0_0_xy_0_xz = buffer_1100_sdsd[44];

    auto g_x_y_0_0_0_xy_0_yy = buffer_1100_sdsd[45];

    auto g_x_y_0_0_0_xy_0_yz = buffer_1100_sdsd[46];

    auto g_x_y_0_0_0_xy_0_zz = buffer_1100_sdsd[47];

    auto g_x_y_0_0_0_xz_0_xx = buffer_1100_sdsd[48];

    auto g_x_y_0_0_0_xz_0_xy = buffer_1100_sdsd[49];

    auto g_x_y_0_0_0_xz_0_xz = buffer_1100_sdsd[50];

    auto g_x_y_0_0_0_xz_0_yy = buffer_1100_sdsd[51];

    auto g_x_y_0_0_0_xz_0_yz = buffer_1100_sdsd[52];

    auto g_x_y_0_0_0_xz_0_zz = buffer_1100_sdsd[53];

    auto g_x_y_0_0_0_yy_0_xx = buffer_1100_sdsd[54];

    auto g_x_y_0_0_0_yy_0_xy = buffer_1100_sdsd[55];

    auto g_x_y_0_0_0_yy_0_xz = buffer_1100_sdsd[56];

    auto g_x_y_0_0_0_yy_0_yy = buffer_1100_sdsd[57];

    auto g_x_y_0_0_0_yy_0_yz = buffer_1100_sdsd[58];

    auto g_x_y_0_0_0_yy_0_zz = buffer_1100_sdsd[59];

    auto g_x_y_0_0_0_yz_0_xx = buffer_1100_sdsd[60];

    auto g_x_y_0_0_0_yz_0_xy = buffer_1100_sdsd[61];

    auto g_x_y_0_0_0_yz_0_xz = buffer_1100_sdsd[62];

    auto g_x_y_0_0_0_yz_0_yy = buffer_1100_sdsd[63];

    auto g_x_y_0_0_0_yz_0_yz = buffer_1100_sdsd[64];

    auto g_x_y_0_0_0_yz_0_zz = buffer_1100_sdsd[65];

    auto g_x_y_0_0_0_zz_0_xx = buffer_1100_sdsd[66];

    auto g_x_y_0_0_0_zz_0_xy = buffer_1100_sdsd[67];

    auto g_x_y_0_0_0_zz_0_xz = buffer_1100_sdsd[68];

    auto g_x_y_0_0_0_zz_0_yy = buffer_1100_sdsd[69];

    auto g_x_y_0_0_0_zz_0_yz = buffer_1100_sdsd[70];

    auto g_x_y_0_0_0_zz_0_zz = buffer_1100_sdsd[71];

    auto g_x_z_0_0_0_xx_0_xx = buffer_1100_sdsd[72];

    auto g_x_z_0_0_0_xx_0_xy = buffer_1100_sdsd[73];

    auto g_x_z_0_0_0_xx_0_xz = buffer_1100_sdsd[74];

    auto g_x_z_0_0_0_xx_0_yy = buffer_1100_sdsd[75];

    auto g_x_z_0_0_0_xx_0_yz = buffer_1100_sdsd[76];

    auto g_x_z_0_0_0_xx_0_zz = buffer_1100_sdsd[77];

    auto g_x_z_0_0_0_xy_0_xx = buffer_1100_sdsd[78];

    auto g_x_z_0_0_0_xy_0_xy = buffer_1100_sdsd[79];

    auto g_x_z_0_0_0_xy_0_xz = buffer_1100_sdsd[80];

    auto g_x_z_0_0_0_xy_0_yy = buffer_1100_sdsd[81];

    auto g_x_z_0_0_0_xy_0_yz = buffer_1100_sdsd[82];

    auto g_x_z_0_0_0_xy_0_zz = buffer_1100_sdsd[83];

    auto g_x_z_0_0_0_xz_0_xx = buffer_1100_sdsd[84];

    auto g_x_z_0_0_0_xz_0_xy = buffer_1100_sdsd[85];

    auto g_x_z_0_0_0_xz_0_xz = buffer_1100_sdsd[86];

    auto g_x_z_0_0_0_xz_0_yy = buffer_1100_sdsd[87];

    auto g_x_z_0_0_0_xz_0_yz = buffer_1100_sdsd[88];

    auto g_x_z_0_0_0_xz_0_zz = buffer_1100_sdsd[89];

    auto g_x_z_0_0_0_yy_0_xx = buffer_1100_sdsd[90];

    auto g_x_z_0_0_0_yy_0_xy = buffer_1100_sdsd[91];

    auto g_x_z_0_0_0_yy_0_xz = buffer_1100_sdsd[92];

    auto g_x_z_0_0_0_yy_0_yy = buffer_1100_sdsd[93];

    auto g_x_z_0_0_0_yy_0_yz = buffer_1100_sdsd[94];

    auto g_x_z_0_0_0_yy_0_zz = buffer_1100_sdsd[95];

    auto g_x_z_0_0_0_yz_0_xx = buffer_1100_sdsd[96];

    auto g_x_z_0_0_0_yz_0_xy = buffer_1100_sdsd[97];

    auto g_x_z_0_0_0_yz_0_xz = buffer_1100_sdsd[98];

    auto g_x_z_0_0_0_yz_0_yy = buffer_1100_sdsd[99];

    auto g_x_z_0_0_0_yz_0_yz = buffer_1100_sdsd[100];

    auto g_x_z_0_0_0_yz_0_zz = buffer_1100_sdsd[101];

    auto g_x_z_0_0_0_zz_0_xx = buffer_1100_sdsd[102];

    auto g_x_z_0_0_0_zz_0_xy = buffer_1100_sdsd[103];

    auto g_x_z_0_0_0_zz_0_xz = buffer_1100_sdsd[104];

    auto g_x_z_0_0_0_zz_0_yy = buffer_1100_sdsd[105];

    auto g_x_z_0_0_0_zz_0_yz = buffer_1100_sdsd[106];

    auto g_x_z_0_0_0_zz_0_zz = buffer_1100_sdsd[107];

    auto g_y_x_0_0_0_xx_0_xx = buffer_1100_sdsd[108];

    auto g_y_x_0_0_0_xx_0_xy = buffer_1100_sdsd[109];

    auto g_y_x_0_0_0_xx_0_xz = buffer_1100_sdsd[110];

    auto g_y_x_0_0_0_xx_0_yy = buffer_1100_sdsd[111];

    auto g_y_x_0_0_0_xx_0_yz = buffer_1100_sdsd[112];

    auto g_y_x_0_0_0_xx_0_zz = buffer_1100_sdsd[113];

    auto g_y_x_0_0_0_xy_0_xx = buffer_1100_sdsd[114];

    auto g_y_x_0_0_0_xy_0_xy = buffer_1100_sdsd[115];

    auto g_y_x_0_0_0_xy_0_xz = buffer_1100_sdsd[116];

    auto g_y_x_0_0_0_xy_0_yy = buffer_1100_sdsd[117];

    auto g_y_x_0_0_0_xy_0_yz = buffer_1100_sdsd[118];

    auto g_y_x_0_0_0_xy_0_zz = buffer_1100_sdsd[119];

    auto g_y_x_0_0_0_xz_0_xx = buffer_1100_sdsd[120];

    auto g_y_x_0_0_0_xz_0_xy = buffer_1100_sdsd[121];

    auto g_y_x_0_0_0_xz_0_xz = buffer_1100_sdsd[122];

    auto g_y_x_0_0_0_xz_0_yy = buffer_1100_sdsd[123];

    auto g_y_x_0_0_0_xz_0_yz = buffer_1100_sdsd[124];

    auto g_y_x_0_0_0_xz_0_zz = buffer_1100_sdsd[125];

    auto g_y_x_0_0_0_yy_0_xx = buffer_1100_sdsd[126];

    auto g_y_x_0_0_0_yy_0_xy = buffer_1100_sdsd[127];

    auto g_y_x_0_0_0_yy_0_xz = buffer_1100_sdsd[128];

    auto g_y_x_0_0_0_yy_0_yy = buffer_1100_sdsd[129];

    auto g_y_x_0_0_0_yy_0_yz = buffer_1100_sdsd[130];

    auto g_y_x_0_0_0_yy_0_zz = buffer_1100_sdsd[131];

    auto g_y_x_0_0_0_yz_0_xx = buffer_1100_sdsd[132];

    auto g_y_x_0_0_0_yz_0_xy = buffer_1100_sdsd[133];

    auto g_y_x_0_0_0_yz_0_xz = buffer_1100_sdsd[134];

    auto g_y_x_0_0_0_yz_0_yy = buffer_1100_sdsd[135];

    auto g_y_x_0_0_0_yz_0_yz = buffer_1100_sdsd[136];

    auto g_y_x_0_0_0_yz_0_zz = buffer_1100_sdsd[137];

    auto g_y_x_0_0_0_zz_0_xx = buffer_1100_sdsd[138];

    auto g_y_x_0_0_0_zz_0_xy = buffer_1100_sdsd[139];

    auto g_y_x_0_0_0_zz_0_xz = buffer_1100_sdsd[140];

    auto g_y_x_0_0_0_zz_0_yy = buffer_1100_sdsd[141];

    auto g_y_x_0_0_0_zz_0_yz = buffer_1100_sdsd[142];

    auto g_y_x_0_0_0_zz_0_zz = buffer_1100_sdsd[143];

    auto g_y_y_0_0_0_xx_0_xx = buffer_1100_sdsd[144];

    auto g_y_y_0_0_0_xx_0_xy = buffer_1100_sdsd[145];

    auto g_y_y_0_0_0_xx_0_xz = buffer_1100_sdsd[146];

    auto g_y_y_0_0_0_xx_0_yy = buffer_1100_sdsd[147];

    auto g_y_y_0_0_0_xx_0_yz = buffer_1100_sdsd[148];

    auto g_y_y_0_0_0_xx_0_zz = buffer_1100_sdsd[149];

    auto g_y_y_0_0_0_xy_0_xx = buffer_1100_sdsd[150];

    auto g_y_y_0_0_0_xy_0_xy = buffer_1100_sdsd[151];

    auto g_y_y_0_0_0_xy_0_xz = buffer_1100_sdsd[152];

    auto g_y_y_0_0_0_xy_0_yy = buffer_1100_sdsd[153];

    auto g_y_y_0_0_0_xy_0_yz = buffer_1100_sdsd[154];

    auto g_y_y_0_0_0_xy_0_zz = buffer_1100_sdsd[155];

    auto g_y_y_0_0_0_xz_0_xx = buffer_1100_sdsd[156];

    auto g_y_y_0_0_0_xz_0_xy = buffer_1100_sdsd[157];

    auto g_y_y_0_0_0_xz_0_xz = buffer_1100_sdsd[158];

    auto g_y_y_0_0_0_xz_0_yy = buffer_1100_sdsd[159];

    auto g_y_y_0_0_0_xz_0_yz = buffer_1100_sdsd[160];

    auto g_y_y_0_0_0_xz_0_zz = buffer_1100_sdsd[161];

    auto g_y_y_0_0_0_yy_0_xx = buffer_1100_sdsd[162];

    auto g_y_y_0_0_0_yy_0_xy = buffer_1100_sdsd[163];

    auto g_y_y_0_0_0_yy_0_xz = buffer_1100_sdsd[164];

    auto g_y_y_0_0_0_yy_0_yy = buffer_1100_sdsd[165];

    auto g_y_y_0_0_0_yy_0_yz = buffer_1100_sdsd[166];

    auto g_y_y_0_0_0_yy_0_zz = buffer_1100_sdsd[167];

    auto g_y_y_0_0_0_yz_0_xx = buffer_1100_sdsd[168];

    auto g_y_y_0_0_0_yz_0_xy = buffer_1100_sdsd[169];

    auto g_y_y_0_0_0_yz_0_xz = buffer_1100_sdsd[170];

    auto g_y_y_0_0_0_yz_0_yy = buffer_1100_sdsd[171];

    auto g_y_y_0_0_0_yz_0_yz = buffer_1100_sdsd[172];

    auto g_y_y_0_0_0_yz_0_zz = buffer_1100_sdsd[173];

    auto g_y_y_0_0_0_zz_0_xx = buffer_1100_sdsd[174];

    auto g_y_y_0_0_0_zz_0_xy = buffer_1100_sdsd[175];

    auto g_y_y_0_0_0_zz_0_xz = buffer_1100_sdsd[176];

    auto g_y_y_0_0_0_zz_0_yy = buffer_1100_sdsd[177];

    auto g_y_y_0_0_0_zz_0_yz = buffer_1100_sdsd[178];

    auto g_y_y_0_0_0_zz_0_zz = buffer_1100_sdsd[179];

    auto g_y_z_0_0_0_xx_0_xx = buffer_1100_sdsd[180];

    auto g_y_z_0_0_0_xx_0_xy = buffer_1100_sdsd[181];

    auto g_y_z_0_0_0_xx_0_xz = buffer_1100_sdsd[182];

    auto g_y_z_0_0_0_xx_0_yy = buffer_1100_sdsd[183];

    auto g_y_z_0_0_0_xx_0_yz = buffer_1100_sdsd[184];

    auto g_y_z_0_0_0_xx_0_zz = buffer_1100_sdsd[185];

    auto g_y_z_0_0_0_xy_0_xx = buffer_1100_sdsd[186];

    auto g_y_z_0_0_0_xy_0_xy = buffer_1100_sdsd[187];

    auto g_y_z_0_0_0_xy_0_xz = buffer_1100_sdsd[188];

    auto g_y_z_0_0_0_xy_0_yy = buffer_1100_sdsd[189];

    auto g_y_z_0_0_0_xy_0_yz = buffer_1100_sdsd[190];

    auto g_y_z_0_0_0_xy_0_zz = buffer_1100_sdsd[191];

    auto g_y_z_0_0_0_xz_0_xx = buffer_1100_sdsd[192];

    auto g_y_z_0_0_0_xz_0_xy = buffer_1100_sdsd[193];

    auto g_y_z_0_0_0_xz_0_xz = buffer_1100_sdsd[194];

    auto g_y_z_0_0_0_xz_0_yy = buffer_1100_sdsd[195];

    auto g_y_z_0_0_0_xz_0_yz = buffer_1100_sdsd[196];

    auto g_y_z_0_0_0_xz_0_zz = buffer_1100_sdsd[197];

    auto g_y_z_0_0_0_yy_0_xx = buffer_1100_sdsd[198];

    auto g_y_z_0_0_0_yy_0_xy = buffer_1100_sdsd[199];

    auto g_y_z_0_0_0_yy_0_xz = buffer_1100_sdsd[200];

    auto g_y_z_0_0_0_yy_0_yy = buffer_1100_sdsd[201];

    auto g_y_z_0_0_0_yy_0_yz = buffer_1100_sdsd[202];

    auto g_y_z_0_0_0_yy_0_zz = buffer_1100_sdsd[203];

    auto g_y_z_0_0_0_yz_0_xx = buffer_1100_sdsd[204];

    auto g_y_z_0_0_0_yz_0_xy = buffer_1100_sdsd[205];

    auto g_y_z_0_0_0_yz_0_xz = buffer_1100_sdsd[206];

    auto g_y_z_0_0_0_yz_0_yy = buffer_1100_sdsd[207];

    auto g_y_z_0_0_0_yz_0_yz = buffer_1100_sdsd[208];

    auto g_y_z_0_0_0_yz_0_zz = buffer_1100_sdsd[209];

    auto g_y_z_0_0_0_zz_0_xx = buffer_1100_sdsd[210];

    auto g_y_z_0_0_0_zz_0_xy = buffer_1100_sdsd[211];

    auto g_y_z_0_0_0_zz_0_xz = buffer_1100_sdsd[212];

    auto g_y_z_0_0_0_zz_0_yy = buffer_1100_sdsd[213];

    auto g_y_z_0_0_0_zz_0_yz = buffer_1100_sdsd[214];

    auto g_y_z_0_0_0_zz_0_zz = buffer_1100_sdsd[215];

    auto g_z_x_0_0_0_xx_0_xx = buffer_1100_sdsd[216];

    auto g_z_x_0_0_0_xx_0_xy = buffer_1100_sdsd[217];

    auto g_z_x_0_0_0_xx_0_xz = buffer_1100_sdsd[218];

    auto g_z_x_0_0_0_xx_0_yy = buffer_1100_sdsd[219];

    auto g_z_x_0_0_0_xx_0_yz = buffer_1100_sdsd[220];

    auto g_z_x_0_0_0_xx_0_zz = buffer_1100_sdsd[221];

    auto g_z_x_0_0_0_xy_0_xx = buffer_1100_sdsd[222];

    auto g_z_x_0_0_0_xy_0_xy = buffer_1100_sdsd[223];

    auto g_z_x_0_0_0_xy_0_xz = buffer_1100_sdsd[224];

    auto g_z_x_0_0_0_xy_0_yy = buffer_1100_sdsd[225];

    auto g_z_x_0_0_0_xy_0_yz = buffer_1100_sdsd[226];

    auto g_z_x_0_0_0_xy_0_zz = buffer_1100_sdsd[227];

    auto g_z_x_0_0_0_xz_0_xx = buffer_1100_sdsd[228];

    auto g_z_x_0_0_0_xz_0_xy = buffer_1100_sdsd[229];

    auto g_z_x_0_0_0_xz_0_xz = buffer_1100_sdsd[230];

    auto g_z_x_0_0_0_xz_0_yy = buffer_1100_sdsd[231];

    auto g_z_x_0_0_0_xz_0_yz = buffer_1100_sdsd[232];

    auto g_z_x_0_0_0_xz_0_zz = buffer_1100_sdsd[233];

    auto g_z_x_0_0_0_yy_0_xx = buffer_1100_sdsd[234];

    auto g_z_x_0_0_0_yy_0_xy = buffer_1100_sdsd[235];

    auto g_z_x_0_0_0_yy_0_xz = buffer_1100_sdsd[236];

    auto g_z_x_0_0_0_yy_0_yy = buffer_1100_sdsd[237];

    auto g_z_x_0_0_0_yy_0_yz = buffer_1100_sdsd[238];

    auto g_z_x_0_0_0_yy_0_zz = buffer_1100_sdsd[239];

    auto g_z_x_0_0_0_yz_0_xx = buffer_1100_sdsd[240];

    auto g_z_x_0_0_0_yz_0_xy = buffer_1100_sdsd[241];

    auto g_z_x_0_0_0_yz_0_xz = buffer_1100_sdsd[242];

    auto g_z_x_0_0_0_yz_0_yy = buffer_1100_sdsd[243];

    auto g_z_x_0_0_0_yz_0_yz = buffer_1100_sdsd[244];

    auto g_z_x_0_0_0_yz_0_zz = buffer_1100_sdsd[245];

    auto g_z_x_0_0_0_zz_0_xx = buffer_1100_sdsd[246];

    auto g_z_x_0_0_0_zz_0_xy = buffer_1100_sdsd[247];

    auto g_z_x_0_0_0_zz_0_xz = buffer_1100_sdsd[248];

    auto g_z_x_0_0_0_zz_0_yy = buffer_1100_sdsd[249];

    auto g_z_x_0_0_0_zz_0_yz = buffer_1100_sdsd[250];

    auto g_z_x_0_0_0_zz_0_zz = buffer_1100_sdsd[251];

    auto g_z_y_0_0_0_xx_0_xx = buffer_1100_sdsd[252];

    auto g_z_y_0_0_0_xx_0_xy = buffer_1100_sdsd[253];

    auto g_z_y_0_0_0_xx_0_xz = buffer_1100_sdsd[254];

    auto g_z_y_0_0_0_xx_0_yy = buffer_1100_sdsd[255];

    auto g_z_y_0_0_0_xx_0_yz = buffer_1100_sdsd[256];

    auto g_z_y_0_0_0_xx_0_zz = buffer_1100_sdsd[257];

    auto g_z_y_0_0_0_xy_0_xx = buffer_1100_sdsd[258];

    auto g_z_y_0_0_0_xy_0_xy = buffer_1100_sdsd[259];

    auto g_z_y_0_0_0_xy_0_xz = buffer_1100_sdsd[260];

    auto g_z_y_0_0_0_xy_0_yy = buffer_1100_sdsd[261];

    auto g_z_y_0_0_0_xy_0_yz = buffer_1100_sdsd[262];

    auto g_z_y_0_0_0_xy_0_zz = buffer_1100_sdsd[263];

    auto g_z_y_0_0_0_xz_0_xx = buffer_1100_sdsd[264];

    auto g_z_y_0_0_0_xz_0_xy = buffer_1100_sdsd[265];

    auto g_z_y_0_0_0_xz_0_xz = buffer_1100_sdsd[266];

    auto g_z_y_0_0_0_xz_0_yy = buffer_1100_sdsd[267];

    auto g_z_y_0_0_0_xz_0_yz = buffer_1100_sdsd[268];

    auto g_z_y_0_0_0_xz_0_zz = buffer_1100_sdsd[269];

    auto g_z_y_0_0_0_yy_0_xx = buffer_1100_sdsd[270];

    auto g_z_y_0_0_0_yy_0_xy = buffer_1100_sdsd[271];

    auto g_z_y_0_0_0_yy_0_xz = buffer_1100_sdsd[272];

    auto g_z_y_0_0_0_yy_0_yy = buffer_1100_sdsd[273];

    auto g_z_y_0_0_0_yy_0_yz = buffer_1100_sdsd[274];

    auto g_z_y_0_0_0_yy_0_zz = buffer_1100_sdsd[275];

    auto g_z_y_0_0_0_yz_0_xx = buffer_1100_sdsd[276];

    auto g_z_y_0_0_0_yz_0_xy = buffer_1100_sdsd[277];

    auto g_z_y_0_0_0_yz_0_xz = buffer_1100_sdsd[278];

    auto g_z_y_0_0_0_yz_0_yy = buffer_1100_sdsd[279];

    auto g_z_y_0_0_0_yz_0_yz = buffer_1100_sdsd[280];

    auto g_z_y_0_0_0_yz_0_zz = buffer_1100_sdsd[281];

    auto g_z_y_0_0_0_zz_0_xx = buffer_1100_sdsd[282];

    auto g_z_y_0_0_0_zz_0_xy = buffer_1100_sdsd[283];

    auto g_z_y_0_0_0_zz_0_xz = buffer_1100_sdsd[284];

    auto g_z_y_0_0_0_zz_0_yy = buffer_1100_sdsd[285];

    auto g_z_y_0_0_0_zz_0_yz = buffer_1100_sdsd[286];

    auto g_z_y_0_0_0_zz_0_zz = buffer_1100_sdsd[287];

    auto g_z_z_0_0_0_xx_0_xx = buffer_1100_sdsd[288];

    auto g_z_z_0_0_0_xx_0_xy = buffer_1100_sdsd[289];

    auto g_z_z_0_0_0_xx_0_xz = buffer_1100_sdsd[290];

    auto g_z_z_0_0_0_xx_0_yy = buffer_1100_sdsd[291];

    auto g_z_z_0_0_0_xx_0_yz = buffer_1100_sdsd[292];

    auto g_z_z_0_0_0_xx_0_zz = buffer_1100_sdsd[293];

    auto g_z_z_0_0_0_xy_0_xx = buffer_1100_sdsd[294];

    auto g_z_z_0_0_0_xy_0_xy = buffer_1100_sdsd[295];

    auto g_z_z_0_0_0_xy_0_xz = buffer_1100_sdsd[296];

    auto g_z_z_0_0_0_xy_0_yy = buffer_1100_sdsd[297];

    auto g_z_z_0_0_0_xy_0_yz = buffer_1100_sdsd[298];

    auto g_z_z_0_0_0_xy_0_zz = buffer_1100_sdsd[299];

    auto g_z_z_0_0_0_xz_0_xx = buffer_1100_sdsd[300];

    auto g_z_z_0_0_0_xz_0_xy = buffer_1100_sdsd[301];

    auto g_z_z_0_0_0_xz_0_xz = buffer_1100_sdsd[302];

    auto g_z_z_0_0_0_xz_0_yy = buffer_1100_sdsd[303];

    auto g_z_z_0_0_0_xz_0_yz = buffer_1100_sdsd[304];

    auto g_z_z_0_0_0_xz_0_zz = buffer_1100_sdsd[305];

    auto g_z_z_0_0_0_yy_0_xx = buffer_1100_sdsd[306];

    auto g_z_z_0_0_0_yy_0_xy = buffer_1100_sdsd[307];

    auto g_z_z_0_0_0_yy_0_xz = buffer_1100_sdsd[308];

    auto g_z_z_0_0_0_yy_0_yy = buffer_1100_sdsd[309];

    auto g_z_z_0_0_0_yy_0_yz = buffer_1100_sdsd[310];

    auto g_z_z_0_0_0_yy_0_zz = buffer_1100_sdsd[311];

    auto g_z_z_0_0_0_yz_0_xx = buffer_1100_sdsd[312];

    auto g_z_z_0_0_0_yz_0_xy = buffer_1100_sdsd[313];

    auto g_z_z_0_0_0_yz_0_xz = buffer_1100_sdsd[314];

    auto g_z_z_0_0_0_yz_0_yy = buffer_1100_sdsd[315];

    auto g_z_z_0_0_0_yz_0_yz = buffer_1100_sdsd[316];

    auto g_z_z_0_0_0_yz_0_zz = buffer_1100_sdsd[317];

    auto g_z_z_0_0_0_zz_0_xx = buffer_1100_sdsd[318];

    auto g_z_z_0_0_0_zz_0_xy = buffer_1100_sdsd[319];

    auto g_z_z_0_0_0_zz_0_xz = buffer_1100_sdsd[320];

    auto g_z_z_0_0_0_zz_0_yy = buffer_1100_sdsd[321];

    auto g_z_z_0_0_0_zz_0_yz = buffer_1100_sdsd[322];

    auto g_z_z_0_0_0_zz_0_zz = buffer_1100_sdsd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_0_xx, g_x_x_0_0_0_xx_0_xy, g_x_x_0_0_0_xx_0_xz, g_x_x_0_0_0_xx_0_yy, g_x_x_0_0_0_xx_0_yz, g_x_x_0_0_0_xx_0_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xxx_0_xx, g_x_xxx_0_xy, g_x_xxx_0_xz, g_x_xxx_0_yy, g_x_xxx_0_yz, g_x_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_0_xx[i] = -4.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_xy[i] = -4.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_xz[i] = -4.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_yy[i] = -4.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_yz[i] = -4.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_zz[i] = -4.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_0_xx, g_x_x_0_0_0_xy_0_xy, g_x_x_0_0_0_xy_0_xz, g_x_x_0_0_0_xy_0_yy, g_x_x_0_0_0_xy_0_yz, g_x_x_0_0_0_xy_0_zz, g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_0_xx, g_x_x_0_0_0_xz_0_xy, g_x_x_0_0_0_xz_0_xz, g_x_x_0_0_0_xz_0_yy, g_x_x_0_0_0_xz_0_yz, g_x_x_0_0_0_xz_0_zz, g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_0_xx, g_x_x_0_0_0_yy_0_xy, g_x_x_0_0_0_yy_0_xz, g_x_x_0_0_0_yy_0_yy, g_x_x_0_0_0_yy_0_yz, g_x_x_0_0_0_yy_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_0_xx[i] = 4.0 * g_x_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_xy[i] = 4.0 * g_x_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_xz[i] = 4.0 * g_x_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_yy[i] = 4.0 * g_x_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_yz[i] = 4.0 * g_x_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_zz[i] = 4.0 * g_x_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_0_xx, g_x_x_0_0_0_yz_0_xy, g_x_x_0_0_0_yz_0_xz, g_x_x_0_0_0_yz_0_yy, g_x_x_0_0_0_yz_0_yz, g_x_x_0_0_0_yz_0_zz, g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_0_xx[i] = 4.0 * g_x_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_xy[i] = 4.0 * g_x_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_xz[i] = 4.0 * g_x_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_yy[i] = 4.0 * g_x_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_yz[i] = 4.0 * g_x_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_zz[i] = 4.0 * g_x_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_0_xx, g_x_x_0_0_0_zz_0_xy, g_x_x_0_0_0_zz_0_xz, g_x_x_0_0_0_zz_0_yy, g_x_x_0_0_0_zz_0_yz, g_x_x_0_0_0_zz_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_0_xx[i] = 4.0 * g_x_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_xy[i] = 4.0 * g_x_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_xz[i] = 4.0 * g_x_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_yy[i] = 4.0 * g_x_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_yz[i] = 4.0 * g_x_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_zz[i] = 4.0 * g_x_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_0_0_xx_0_xx, g_x_y_0_0_0_xx_0_xy, g_x_y_0_0_0_xx_0_xz, g_x_y_0_0_0_xx_0_yy, g_x_y_0_0_0_xx_0_yz, g_x_y_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_0_xx[i] = 4.0 * g_x_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_xy[i] = 4.0 * g_x_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_xz[i] = 4.0 * g_x_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_yy[i] = 4.0 * g_x_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_yz[i] = 4.0 * g_x_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_zz[i] = 4.0 * g_x_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_x_y_0_0_0_xy_0_xx, g_x_y_0_0_0_xy_0_xy, g_x_y_0_0_0_xy_0_xz, g_x_y_0_0_0_xy_0_yy, g_x_y_0_0_0_xy_0_yz, g_x_y_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_x_y_0_0_0_xz_0_xx, g_x_y_0_0_0_xz_0_xy, g_x_y_0_0_0_xz_0_xz, g_x_y_0_0_0_xz_0_yy, g_x_y_0_0_0_xz_0_yz, g_x_y_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_0_xx[i] = 4.0 * g_x_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_xy[i] = 4.0 * g_x_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_xz[i] = 4.0 * g_x_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_yy[i] = 4.0 * g_x_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_yz[i] = 4.0 * g_x_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_zz[i] = 4.0 * g_x_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_0_xx, g_x_y_0_0_0_yy_0_xy, g_x_y_0_0_0_yy_0_xz, g_x_y_0_0_0_yy_0_yy, g_x_y_0_0_0_yy_0_yz, g_x_y_0_0_0_yy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yyy_0_xx, g_x_yyy_0_xy, g_x_yyy_0_xz, g_x_yyy_0_yy, g_x_yyy_0_yz, g_x_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_0_xx[i] = -4.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_xy[i] = -4.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_xz[i] = -4.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_yy[i] = -4.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_yz[i] = -4.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_zz[i] = -4.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_0_xx, g_x_y_0_0_0_yz_0_xy, g_x_y_0_0_0_yz_0_xz, g_x_y_0_0_0_yz_0_yy, g_x_y_0_0_0_yz_0_yz, g_x_y_0_0_0_yz_0_zz, g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_0_xx, g_x_y_0_0_0_zz_0_xy, g_x_y_0_0_0_zz_0_xz, g_x_y_0_0_0_zz_0_yy, g_x_y_0_0_0_zz_0_yz, g_x_y_0_0_0_zz_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_0_xx[i] = 4.0 * g_x_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_xy[i] = 4.0 * g_x_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_xz[i] = 4.0 * g_x_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_yy[i] = 4.0 * g_x_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_yz[i] = 4.0 * g_x_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_zz[i] = 4.0 * g_x_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_0_0_xx_0_xx, g_x_z_0_0_0_xx_0_xy, g_x_z_0_0_0_xx_0_xz, g_x_z_0_0_0_xx_0_yy, g_x_z_0_0_0_xx_0_yz, g_x_z_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_0_xx[i] = 4.0 * g_x_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_xy[i] = 4.0 * g_x_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_xz[i] = 4.0 * g_x_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_yy[i] = 4.0 * g_x_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_yz[i] = 4.0 * g_x_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_zz[i] = 4.0 * g_x_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_x_z_0_0_0_xy_0_xx, g_x_z_0_0_0_xy_0_xy, g_x_z_0_0_0_xy_0_xz, g_x_z_0_0_0_xy_0_yy, g_x_z_0_0_0_xy_0_yz, g_x_z_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_0_xx[i] = 4.0 * g_x_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_xy[i] = 4.0 * g_x_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_xz[i] = 4.0 * g_x_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_yy[i] = 4.0 * g_x_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_yz[i] = 4.0 * g_x_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_zz[i] = 4.0 * g_x_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_x_z_0_0_0_xz_0_xx, g_x_z_0_0_0_xz_0_xy, g_x_z_0_0_0_xz_0_xz, g_x_z_0_0_0_xz_0_yy, g_x_z_0_0_0_xz_0_yz, g_x_z_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_0_0_yy_0_xx, g_x_z_0_0_0_yy_0_xy, g_x_z_0_0_0_yy_0_xz, g_x_z_0_0_0_yy_0_yy, g_x_z_0_0_0_yy_0_yz, g_x_z_0_0_0_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_0_xx[i] = 4.0 * g_x_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_xy[i] = 4.0 * g_x_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_xz[i] = 4.0 * g_x_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_yy[i] = 4.0 * g_x_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_yz[i] = 4.0 * g_x_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_zz[i] = 4.0 * g_x_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_x_z_0_0_0_yz_0_xx, g_x_z_0_0_0_yz_0_xy, g_x_z_0_0_0_yz_0_xz, g_x_z_0_0_0_yz_0_yy, g_x_z_0_0_0_yz_0_yz, g_x_z_0_0_0_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_0_xx, g_x_z_0_0_0_zz_0_xy, g_x_z_0_0_0_zz_0_xz, g_x_z_0_0_0_zz_0_yy, g_x_z_0_0_0_zz_0_yz, g_x_z_0_0_0_zz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_zzz_0_xx, g_x_zzz_0_xy, g_x_zzz_0_xz, g_x_zzz_0_yy, g_x_zzz_0_yz, g_x_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_0_xx[i] = -4.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_xy[i] = -4.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_xz[i] = -4.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_yy[i] = -4.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_yz[i] = -4.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_zz[i] = -4.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_0_xx, g_y_x_0_0_0_xx_0_xy, g_y_x_0_0_0_xx_0_xz, g_y_x_0_0_0_xx_0_yy, g_y_x_0_0_0_xx_0_yz, g_y_x_0_0_0_xx_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xxx_0_xx, g_y_xxx_0_xy, g_y_xxx_0_xz, g_y_xxx_0_yy, g_y_xxx_0_yz, g_y_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_0_xx[i] = -4.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_xy[i] = -4.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_xz[i] = -4.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_yy[i] = -4.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_yz[i] = -4.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_zz[i] = -4.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_0_xx, g_y_x_0_0_0_xy_0_xy, g_y_x_0_0_0_xy_0_xz, g_y_x_0_0_0_xy_0_yy, g_y_x_0_0_0_xy_0_yz, g_y_x_0_0_0_xy_0_zz, g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_0_xx, g_y_x_0_0_0_xz_0_xy, g_y_x_0_0_0_xz_0_xz, g_y_x_0_0_0_xz_0_yy, g_y_x_0_0_0_xz_0_yz, g_y_x_0_0_0_xz_0_zz, g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_0_xx, g_y_x_0_0_0_yy_0_xy, g_y_x_0_0_0_yy_0_xz, g_y_x_0_0_0_yy_0_yy, g_y_x_0_0_0_yy_0_yz, g_y_x_0_0_0_yy_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_0_xx[i] = 4.0 * g_y_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_xy[i] = 4.0 * g_y_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_xz[i] = 4.0 * g_y_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_yy[i] = 4.0 * g_y_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_yz[i] = 4.0 * g_y_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_zz[i] = 4.0 * g_y_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_0_xx, g_y_x_0_0_0_yz_0_xy, g_y_x_0_0_0_yz_0_xz, g_y_x_0_0_0_yz_0_yy, g_y_x_0_0_0_yz_0_yz, g_y_x_0_0_0_yz_0_zz, g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_0_xx[i] = 4.0 * g_y_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_xy[i] = 4.0 * g_y_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_xz[i] = 4.0 * g_y_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_yy[i] = 4.0 * g_y_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_yz[i] = 4.0 * g_y_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_zz[i] = 4.0 * g_y_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_0_xx, g_y_x_0_0_0_zz_0_xy, g_y_x_0_0_0_zz_0_xz, g_y_x_0_0_0_zz_0_yy, g_y_x_0_0_0_zz_0_yz, g_y_x_0_0_0_zz_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_0_xx[i] = 4.0 * g_y_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_xy[i] = 4.0 * g_y_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_xz[i] = 4.0 * g_y_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_yy[i] = 4.0 * g_y_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_yz[i] = 4.0 * g_y_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_zz[i] = 4.0 * g_y_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_0_0_xx_0_xx, g_y_y_0_0_0_xx_0_xy, g_y_y_0_0_0_xx_0_xz, g_y_y_0_0_0_xx_0_yy, g_y_y_0_0_0_xx_0_yz, g_y_y_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_0_xx[i] = 4.0 * g_y_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_xy[i] = 4.0 * g_y_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_xz[i] = 4.0 * g_y_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_yy[i] = 4.0 * g_y_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_yz[i] = 4.0 * g_y_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_zz[i] = 4.0 * g_y_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz, g_y_y_0_0_0_xy_0_xx, g_y_y_0_0_0_xy_0_xy, g_y_y_0_0_0_xy_0_xz, g_y_y_0_0_0_xy_0_yy, g_y_y_0_0_0_xy_0_yz, g_y_y_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_y_y_0_0_0_xz_0_xx, g_y_y_0_0_0_xz_0_xy, g_y_y_0_0_0_xz_0_xz, g_y_y_0_0_0_xz_0_yy, g_y_y_0_0_0_xz_0_yz, g_y_y_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_0_xx[i] = 4.0 * g_y_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_xy[i] = 4.0 * g_y_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_xz[i] = 4.0 * g_y_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_yy[i] = 4.0 * g_y_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_yz[i] = 4.0 * g_y_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_zz[i] = 4.0 * g_y_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_0_xx, g_y_y_0_0_0_yy_0_xy, g_y_y_0_0_0_yy_0_xz, g_y_y_0_0_0_yy_0_yy, g_y_y_0_0_0_yy_0_yz, g_y_y_0_0_0_yy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yyy_0_xx, g_y_yyy_0_xy, g_y_yyy_0_xz, g_y_yyy_0_yy, g_y_yyy_0_yz, g_y_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_0_xx[i] = -4.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_xy[i] = -4.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_xz[i] = -4.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_yy[i] = -4.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_yz[i] = -4.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_zz[i] = -4.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_0_xx, g_y_y_0_0_0_yz_0_xy, g_y_y_0_0_0_yz_0_xz, g_y_y_0_0_0_yz_0_yy, g_y_y_0_0_0_yz_0_yz, g_y_y_0_0_0_yz_0_zz, g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_0_xx, g_y_y_0_0_0_zz_0_xy, g_y_y_0_0_0_zz_0_xz, g_y_y_0_0_0_zz_0_yy, g_y_y_0_0_0_zz_0_yz, g_y_y_0_0_0_zz_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_0_xx[i] = 4.0 * g_y_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_xy[i] = 4.0 * g_y_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_xz[i] = 4.0 * g_y_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_yy[i] = 4.0 * g_y_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_yz[i] = 4.0 * g_y_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_zz[i] = 4.0 * g_y_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_0_0_xx_0_xx, g_y_z_0_0_0_xx_0_xy, g_y_z_0_0_0_xx_0_xz, g_y_z_0_0_0_xx_0_yy, g_y_z_0_0_0_xx_0_yz, g_y_z_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_0_xx[i] = 4.0 * g_y_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_xy[i] = 4.0 * g_y_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_xz[i] = 4.0 * g_y_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_yy[i] = 4.0 * g_y_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_yz[i] = 4.0 * g_y_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_zz[i] = 4.0 * g_y_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_y_z_0_0_0_xy_0_xx, g_y_z_0_0_0_xy_0_xy, g_y_z_0_0_0_xy_0_xz, g_y_z_0_0_0_xy_0_yy, g_y_z_0_0_0_xy_0_yz, g_y_z_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_0_xx[i] = 4.0 * g_y_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_xy[i] = 4.0 * g_y_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_xz[i] = 4.0 * g_y_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_yy[i] = 4.0 * g_y_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_yz[i] = 4.0 * g_y_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_zz[i] = 4.0 * g_y_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz, g_y_z_0_0_0_xz_0_xx, g_y_z_0_0_0_xz_0_xy, g_y_z_0_0_0_xz_0_xz, g_y_z_0_0_0_xz_0_yy, g_y_z_0_0_0_xz_0_yz, g_y_z_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_0_0_yy_0_xx, g_y_z_0_0_0_yy_0_xy, g_y_z_0_0_0_yy_0_xz, g_y_z_0_0_0_yy_0_yy, g_y_z_0_0_0_yy_0_yz, g_y_z_0_0_0_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_0_xx[i] = 4.0 * g_y_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_xy[i] = 4.0 * g_y_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_xz[i] = 4.0 * g_y_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_yy[i] = 4.0 * g_y_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_yz[i] = 4.0 * g_y_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_zz[i] = 4.0 * g_y_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz, g_y_z_0_0_0_yz_0_xx, g_y_z_0_0_0_yz_0_xy, g_y_z_0_0_0_yz_0_xz, g_y_z_0_0_0_yz_0_yy, g_y_z_0_0_0_yz_0_yz, g_y_z_0_0_0_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_0_xx, g_y_z_0_0_0_zz_0_xy, g_y_z_0_0_0_zz_0_xz, g_y_z_0_0_0_zz_0_yy, g_y_z_0_0_0_zz_0_yz, g_y_z_0_0_0_zz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_zzz_0_xx, g_y_zzz_0_xy, g_y_zzz_0_xz, g_y_zzz_0_yy, g_y_zzz_0_yz, g_y_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_0_xx[i] = -4.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_xy[i] = -4.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_xz[i] = -4.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_yy[i] = -4.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_yz[i] = -4.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_zz[i] = -4.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_0_xx, g_z_x_0_0_0_xx_0_xy, g_z_x_0_0_0_xx_0_xz, g_z_x_0_0_0_xx_0_yy, g_z_x_0_0_0_xx_0_yz, g_z_x_0_0_0_xx_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xxx_0_xx, g_z_xxx_0_xy, g_z_xxx_0_xz, g_z_xxx_0_yy, g_z_xxx_0_yz, g_z_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_0_xx[i] = -4.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_xy[i] = -4.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_xz[i] = -4.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_yy[i] = -4.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_yz[i] = -4.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_zz[i] = -4.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_0_xx, g_z_x_0_0_0_xy_0_xy, g_z_x_0_0_0_xy_0_xz, g_z_x_0_0_0_xy_0_yy, g_z_x_0_0_0_xy_0_yz, g_z_x_0_0_0_xy_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_0_xx, g_z_x_0_0_0_xz_0_xy, g_z_x_0_0_0_xz_0_xz, g_z_x_0_0_0_xz_0_yy, g_z_x_0_0_0_xz_0_yz, g_z_x_0_0_0_xz_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_0_xx, g_z_x_0_0_0_yy_0_xy, g_z_x_0_0_0_yy_0_xz, g_z_x_0_0_0_yy_0_yy, g_z_x_0_0_0_yy_0_yz, g_z_x_0_0_0_yy_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_0_xx[i] = 4.0 * g_z_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_xy[i] = 4.0 * g_z_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_xz[i] = 4.0 * g_z_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_yy[i] = 4.0 * g_z_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_yz[i] = 4.0 * g_z_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_zz[i] = 4.0 * g_z_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_0_xx, g_z_x_0_0_0_yz_0_xy, g_z_x_0_0_0_yz_0_xz, g_z_x_0_0_0_yz_0_yy, g_z_x_0_0_0_yz_0_yz, g_z_x_0_0_0_yz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_0_xx[i] = 4.0 * g_z_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_xy[i] = 4.0 * g_z_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_xz[i] = 4.0 * g_z_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_yy[i] = 4.0 * g_z_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_yz[i] = 4.0 * g_z_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_zz[i] = 4.0 * g_z_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_0_xx, g_z_x_0_0_0_zz_0_xy, g_z_x_0_0_0_zz_0_xz, g_z_x_0_0_0_zz_0_yy, g_z_x_0_0_0_zz_0_yz, g_z_x_0_0_0_zz_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_0_xx[i] = 4.0 * g_z_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_xy[i] = 4.0 * g_z_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_xz[i] = 4.0 * g_z_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_yy[i] = 4.0 * g_z_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_yz[i] = 4.0 * g_z_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_zz[i] = 4.0 * g_z_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_0_0_xx_0_xx, g_z_y_0_0_0_xx_0_xy, g_z_y_0_0_0_xx_0_xz, g_z_y_0_0_0_xx_0_yy, g_z_y_0_0_0_xx_0_yz, g_z_y_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_0_xx[i] = 4.0 * g_z_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_xy[i] = 4.0 * g_z_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_xz[i] = 4.0 * g_z_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_yy[i] = 4.0 * g_z_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_yz[i] = 4.0 * g_z_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_zz[i] = 4.0 * g_z_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz, g_z_y_0_0_0_xy_0_xx, g_z_y_0_0_0_xy_0_xy, g_z_y_0_0_0_xy_0_xz, g_z_y_0_0_0_xy_0_yy, g_z_y_0_0_0_xy_0_yz, g_z_y_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz, g_z_y_0_0_0_xz_0_xx, g_z_y_0_0_0_xz_0_xy, g_z_y_0_0_0_xz_0_xz, g_z_y_0_0_0_xz_0_yy, g_z_y_0_0_0_xz_0_yz, g_z_y_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_0_xx[i] = 4.0 * g_z_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_xy[i] = 4.0 * g_z_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_xz[i] = 4.0 * g_z_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_yy[i] = 4.0 * g_z_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_yz[i] = 4.0 * g_z_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_zz[i] = 4.0 * g_z_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_0_xx, g_z_y_0_0_0_yy_0_xy, g_z_y_0_0_0_yy_0_xz, g_z_y_0_0_0_yy_0_yy, g_z_y_0_0_0_yy_0_yz, g_z_y_0_0_0_yy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yyy_0_xx, g_z_yyy_0_xy, g_z_yyy_0_xz, g_z_yyy_0_yy, g_z_yyy_0_yz, g_z_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_0_xx[i] = -4.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_xy[i] = -4.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_xz[i] = -4.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_yy[i] = -4.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_yz[i] = -4.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_zz[i] = -4.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_0_xx, g_z_y_0_0_0_yz_0_xy, g_z_y_0_0_0_yz_0_xz, g_z_y_0_0_0_yz_0_yy, g_z_y_0_0_0_yz_0_yz, g_z_y_0_0_0_yz_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_0_xx, g_z_y_0_0_0_zz_0_xy, g_z_y_0_0_0_zz_0_xz, g_z_y_0_0_0_zz_0_yy, g_z_y_0_0_0_zz_0_yz, g_z_y_0_0_0_zz_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_0_xx[i] = 4.0 * g_z_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_xy[i] = 4.0 * g_z_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_xz[i] = 4.0 * g_z_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_yy[i] = 4.0 * g_z_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_yz[i] = 4.0 * g_z_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_zz[i] = 4.0 * g_z_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_0_0_xx_0_xx, g_z_z_0_0_0_xx_0_xy, g_z_z_0_0_0_xx_0_xz, g_z_z_0_0_0_xx_0_yy, g_z_z_0_0_0_xx_0_yz, g_z_z_0_0_0_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_0_xx[i] = 4.0 * g_z_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_xy[i] = 4.0 * g_z_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_xz[i] = 4.0 * g_z_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_yy[i] = 4.0 * g_z_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_yz[i] = 4.0 * g_z_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_zz[i] = 4.0 * g_z_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz, g_z_z_0_0_0_xy_0_xx, g_z_z_0_0_0_xy_0_xy, g_z_z_0_0_0_xy_0_xz, g_z_z_0_0_0_xy_0_yy, g_z_z_0_0_0_xy_0_yz, g_z_z_0_0_0_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_0_xx[i] = 4.0 * g_z_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_xy[i] = 4.0 * g_z_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_xz[i] = 4.0 * g_z_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_yy[i] = 4.0 * g_z_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_yz[i] = 4.0 * g_z_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_zz[i] = 4.0 * g_z_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz, g_z_z_0_0_0_xz_0_xx, g_z_z_0_0_0_xz_0_xy, g_z_z_0_0_0_xz_0_xz, g_z_z_0_0_0_xz_0_yy, g_z_z_0_0_0_xz_0_yz, g_z_z_0_0_0_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_0_0_yy_0_xx, g_z_z_0_0_0_yy_0_xy, g_z_z_0_0_0_yy_0_xz, g_z_z_0_0_0_yy_0_yy, g_z_z_0_0_0_yy_0_yz, g_z_z_0_0_0_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_0_xx[i] = 4.0 * g_z_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_xy[i] = 4.0 * g_z_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_xz[i] = 4.0 * g_z_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_yy[i] = 4.0 * g_z_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_yz[i] = 4.0 * g_z_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_zz[i] = 4.0 * g_z_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz, g_z_z_0_0_0_yz_0_xx, g_z_z_0_0_0_yz_0_xy, g_z_z_0_0_0_yz_0_xz, g_z_z_0_0_0_yz_0_yy, g_z_z_0_0_0_yz_0_yz, g_z_z_0_0_0_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_0_xx, g_z_z_0_0_0_zz_0_xy, g_z_z_0_0_0_zz_0_xz, g_z_z_0_0_0_zz_0_yy, g_z_z_0_0_0_zz_0_yz, g_z_z_0_0_0_zz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_zzz_0_xx, g_z_zzz_0_xy, g_z_zzz_0_xz, g_z_zzz_0_yy, g_z_zzz_0_yz, g_z_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_0_xx[i] = -4.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_xy[i] = -4.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_xz[i] = -4.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_yy[i] = -4.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_yz[i] = -4.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_zz[i] = -4.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_zzz_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


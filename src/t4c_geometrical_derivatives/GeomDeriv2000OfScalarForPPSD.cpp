#include "GeomDeriv2000OfScalarForPPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ppsd_0(CSimdArray<double>& buffer_2000_ppsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_fpsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ppsd.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_fpsd

    auto g_xxx_x_0_xx = buffer_fpsd[0];

    auto g_xxx_x_0_xy = buffer_fpsd[1];

    auto g_xxx_x_0_xz = buffer_fpsd[2];

    auto g_xxx_x_0_yy = buffer_fpsd[3];

    auto g_xxx_x_0_yz = buffer_fpsd[4];

    auto g_xxx_x_0_zz = buffer_fpsd[5];

    auto g_xxx_y_0_xx = buffer_fpsd[6];

    auto g_xxx_y_0_xy = buffer_fpsd[7];

    auto g_xxx_y_0_xz = buffer_fpsd[8];

    auto g_xxx_y_0_yy = buffer_fpsd[9];

    auto g_xxx_y_0_yz = buffer_fpsd[10];

    auto g_xxx_y_0_zz = buffer_fpsd[11];

    auto g_xxx_z_0_xx = buffer_fpsd[12];

    auto g_xxx_z_0_xy = buffer_fpsd[13];

    auto g_xxx_z_0_xz = buffer_fpsd[14];

    auto g_xxx_z_0_yy = buffer_fpsd[15];

    auto g_xxx_z_0_yz = buffer_fpsd[16];

    auto g_xxx_z_0_zz = buffer_fpsd[17];

    auto g_xxy_x_0_xx = buffer_fpsd[18];

    auto g_xxy_x_0_xy = buffer_fpsd[19];

    auto g_xxy_x_0_xz = buffer_fpsd[20];

    auto g_xxy_x_0_yy = buffer_fpsd[21];

    auto g_xxy_x_0_yz = buffer_fpsd[22];

    auto g_xxy_x_0_zz = buffer_fpsd[23];

    auto g_xxy_y_0_xx = buffer_fpsd[24];

    auto g_xxy_y_0_xy = buffer_fpsd[25];

    auto g_xxy_y_0_xz = buffer_fpsd[26];

    auto g_xxy_y_0_yy = buffer_fpsd[27];

    auto g_xxy_y_0_yz = buffer_fpsd[28];

    auto g_xxy_y_0_zz = buffer_fpsd[29];

    auto g_xxy_z_0_xx = buffer_fpsd[30];

    auto g_xxy_z_0_xy = buffer_fpsd[31];

    auto g_xxy_z_0_xz = buffer_fpsd[32];

    auto g_xxy_z_0_yy = buffer_fpsd[33];

    auto g_xxy_z_0_yz = buffer_fpsd[34];

    auto g_xxy_z_0_zz = buffer_fpsd[35];

    auto g_xxz_x_0_xx = buffer_fpsd[36];

    auto g_xxz_x_0_xy = buffer_fpsd[37];

    auto g_xxz_x_0_xz = buffer_fpsd[38];

    auto g_xxz_x_0_yy = buffer_fpsd[39];

    auto g_xxz_x_0_yz = buffer_fpsd[40];

    auto g_xxz_x_0_zz = buffer_fpsd[41];

    auto g_xxz_y_0_xx = buffer_fpsd[42];

    auto g_xxz_y_0_xy = buffer_fpsd[43];

    auto g_xxz_y_0_xz = buffer_fpsd[44];

    auto g_xxz_y_0_yy = buffer_fpsd[45];

    auto g_xxz_y_0_yz = buffer_fpsd[46];

    auto g_xxz_y_0_zz = buffer_fpsd[47];

    auto g_xxz_z_0_xx = buffer_fpsd[48];

    auto g_xxz_z_0_xy = buffer_fpsd[49];

    auto g_xxz_z_0_xz = buffer_fpsd[50];

    auto g_xxz_z_0_yy = buffer_fpsd[51];

    auto g_xxz_z_0_yz = buffer_fpsd[52];

    auto g_xxz_z_0_zz = buffer_fpsd[53];

    auto g_xyy_x_0_xx = buffer_fpsd[54];

    auto g_xyy_x_0_xy = buffer_fpsd[55];

    auto g_xyy_x_0_xz = buffer_fpsd[56];

    auto g_xyy_x_0_yy = buffer_fpsd[57];

    auto g_xyy_x_0_yz = buffer_fpsd[58];

    auto g_xyy_x_0_zz = buffer_fpsd[59];

    auto g_xyy_y_0_xx = buffer_fpsd[60];

    auto g_xyy_y_0_xy = buffer_fpsd[61];

    auto g_xyy_y_0_xz = buffer_fpsd[62];

    auto g_xyy_y_0_yy = buffer_fpsd[63];

    auto g_xyy_y_0_yz = buffer_fpsd[64];

    auto g_xyy_y_0_zz = buffer_fpsd[65];

    auto g_xyy_z_0_xx = buffer_fpsd[66];

    auto g_xyy_z_0_xy = buffer_fpsd[67];

    auto g_xyy_z_0_xz = buffer_fpsd[68];

    auto g_xyy_z_0_yy = buffer_fpsd[69];

    auto g_xyy_z_0_yz = buffer_fpsd[70];

    auto g_xyy_z_0_zz = buffer_fpsd[71];

    auto g_xyz_x_0_xx = buffer_fpsd[72];

    auto g_xyz_x_0_xy = buffer_fpsd[73];

    auto g_xyz_x_0_xz = buffer_fpsd[74];

    auto g_xyz_x_0_yy = buffer_fpsd[75];

    auto g_xyz_x_0_yz = buffer_fpsd[76];

    auto g_xyz_x_0_zz = buffer_fpsd[77];

    auto g_xyz_y_0_xx = buffer_fpsd[78];

    auto g_xyz_y_0_xy = buffer_fpsd[79];

    auto g_xyz_y_0_xz = buffer_fpsd[80];

    auto g_xyz_y_0_yy = buffer_fpsd[81];

    auto g_xyz_y_0_yz = buffer_fpsd[82];

    auto g_xyz_y_0_zz = buffer_fpsd[83];

    auto g_xyz_z_0_xx = buffer_fpsd[84];

    auto g_xyz_z_0_xy = buffer_fpsd[85];

    auto g_xyz_z_0_xz = buffer_fpsd[86];

    auto g_xyz_z_0_yy = buffer_fpsd[87];

    auto g_xyz_z_0_yz = buffer_fpsd[88];

    auto g_xyz_z_0_zz = buffer_fpsd[89];

    auto g_xzz_x_0_xx = buffer_fpsd[90];

    auto g_xzz_x_0_xy = buffer_fpsd[91];

    auto g_xzz_x_0_xz = buffer_fpsd[92];

    auto g_xzz_x_0_yy = buffer_fpsd[93];

    auto g_xzz_x_0_yz = buffer_fpsd[94];

    auto g_xzz_x_0_zz = buffer_fpsd[95];

    auto g_xzz_y_0_xx = buffer_fpsd[96];

    auto g_xzz_y_0_xy = buffer_fpsd[97];

    auto g_xzz_y_0_xz = buffer_fpsd[98];

    auto g_xzz_y_0_yy = buffer_fpsd[99];

    auto g_xzz_y_0_yz = buffer_fpsd[100];

    auto g_xzz_y_0_zz = buffer_fpsd[101];

    auto g_xzz_z_0_xx = buffer_fpsd[102];

    auto g_xzz_z_0_xy = buffer_fpsd[103];

    auto g_xzz_z_0_xz = buffer_fpsd[104];

    auto g_xzz_z_0_yy = buffer_fpsd[105];

    auto g_xzz_z_0_yz = buffer_fpsd[106];

    auto g_xzz_z_0_zz = buffer_fpsd[107];

    auto g_yyy_x_0_xx = buffer_fpsd[108];

    auto g_yyy_x_0_xy = buffer_fpsd[109];

    auto g_yyy_x_0_xz = buffer_fpsd[110];

    auto g_yyy_x_0_yy = buffer_fpsd[111];

    auto g_yyy_x_0_yz = buffer_fpsd[112];

    auto g_yyy_x_0_zz = buffer_fpsd[113];

    auto g_yyy_y_0_xx = buffer_fpsd[114];

    auto g_yyy_y_0_xy = buffer_fpsd[115];

    auto g_yyy_y_0_xz = buffer_fpsd[116];

    auto g_yyy_y_0_yy = buffer_fpsd[117];

    auto g_yyy_y_0_yz = buffer_fpsd[118];

    auto g_yyy_y_0_zz = buffer_fpsd[119];

    auto g_yyy_z_0_xx = buffer_fpsd[120];

    auto g_yyy_z_0_xy = buffer_fpsd[121];

    auto g_yyy_z_0_xz = buffer_fpsd[122];

    auto g_yyy_z_0_yy = buffer_fpsd[123];

    auto g_yyy_z_0_yz = buffer_fpsd[124];

    auto g_yyy_z_0_zz = buffer_fpsd[125];

    auto g_yyz_x_0_xx = buffer_fpsd[126];

    auto g_yyz_x_0_xy = buffer_fpsd[127];

    auto g_yyz_x_0_xz = buffer_fpsd[128];

    auto g_yyz_x_0_yy = buffer_fpsd[129];

    auto g_yyz_x_0_yz = buffer_fpsd[130];

    auto g_yyz_x_0_zz = buffer_fpsd[131];

    auto g_yyz_y_0_xx = buffer_fpsd[132];

    auto g_yyz_y_0_xy = buffer_fpsd[133];

    auto g_yyz_y_0_xz = buffer_fpsd[134];

    auto g_yyz_y_0_yy = buffer_fpsd[135];

    auto g_yyz_y_0_yz = buffer_fpsd[136];

    auto g_yyz_y_0_zz = buffer_fpsd[137];

    auto g_yyz_z_0_xx = buffer_fpsd[138];

    auto g_yyz_z_0_xy = buffer_fpsd[139];

    auto g_yyz_z_0_xz = buffer_fpsd[140];

    auto g_yyz_z_0_yy = buffer_fpsd[141];

    auto g_yyz_z_0_yz = buffer_fpsd[142];

    auto g_yyz_z_0_zz = buffer_fpsd[143];

    auto g_yzz_x_0_xx = buffer_fpsd[144];

    auto g_yzz_x_0_xy = buffer_fpsd[145];

    auto g_yzz_x_0_xz = buffer_fpsd[146];

    auto g_yzz_x_0_yy = buffer_fpsd[147];

    auto g_yzz_x_0_yz = buffer_fpsd[148];

    auto g_yzz_x_0_zz = buffer_fpsd[149];

    auto g_yzz_y_0_xx = buffer_fpsd[150];

    auto g_yzz_y_0_xy = buffer_fpsd[151];

    auto g_yzz_y_0_xz = buffer_fpsd[152];

    auto g_yzz_y_0_yy = buffer_fpsd[153];

    auto g_yzz_y_0_yz = buffer_fpsd[154];

    auto g_yzz_y_0_zz = buffer_fpsd[155];

    auto g_yzz_z_0_xx = buffer_fpsd[156];

    auto g_yzz_z_0_xy = buffer_fpsd[157];

    auto g_yzz_z_0_xz = buffer_fpsd[158];

    auto g_yzz_z_0_yy = buffer_fpsd[159];

    auto g_yzz_z_0_yz = buffer_fpsd[160];

    auto g_yzz_z_0_zz = buffer_fpsd[161];

    auto g_zzz_x_0_xx = buffer_fpsd[162];

    auto g_zzz_x_0_xy = buffer_fpsd[163];

    auto g_zzz_x_0_xz = buffer_fpsd[164];

    auto g_zzz_x_0_yy = buffer_fpsd[165];

    auto g_zzz_x_0_yz = buffer_fpsd[166];

    auto g_zzz_x_0_zz = buffer_fpsd[167];

    auto g_zzz_y_0_xx = buffer_fpsd[168];

    auto g_zzz_y_0_xy = buffer_fpsd[169];

    auto g_zzz_y_0_xz = buffer_fpsd[170];

    auto g_zzz_y_0_yy = buffer_fpsd[171];

    auto g_zzz_y_0_yz = buffer_fpsd[172];

    auto g_zzz_y_0_zz = buffer_fpsd[173];

    auto g_zzz_z_0_xx = buffer_fpsd[174];

    auto g_zzz_z_0_xy = buffer_fpsd[175];

    auto g_zzz_z_0_xz = buffer_fpsd[176];

    auto g_zzz_z_0_yy = buffer_fpsd[177];

    auto g_zzz_z_0_yz = buffer_fpsd[178];

    auto g_zzz_z_0_zz = buffer_fpsd[179];

    /// Set up components of integrals buffer : buffer_2000_ppsd

    auto g_xx_0_0_0_x_x_0_xx = buffer_2000_ppsd[0];

    auto g_xx_0_0_0_x_x_0_xy = buffer_2000_ppsd[1];

    auto g_xx_0_0_0_x_x_0_xz = buffer_2000_ppsd[2];

    auto g_xx_0_0_0_x_x_0_yy = buffer_2000_ppsd[3];

    auto g_xx_0_0_0_x_x_0_yz = buffer_2000_ppsd[4];

    auto g_xx_0_0_0_x_x_0_zz = buffer_2000_ppsd[5];

    auto g_xx_0_0_0_x_y_0_xx = buffer_2000_ppsd[6];

    auto g_xx_0_0_0_x_y_0_xy = buffer_2000_ppsd[7];

    auto g_xx_0_0_0_x_y_0_xz = buffer_2000_ppsd[8];

    auto g_xx_0_0_0_x_y_0_yy = buffer_2000_ppsd[9];

    auto g_xx_0_0_0_x_y_0_yz = buffer_2000_ppsd[10];

    auto g_xx_0_0_0_x_y_0_zz = buffer_2000_ppsd[11];

    auto g_xx_0_0_0_x_z_0_xx = buffer_2000_ppsd[12];

    auto g_xx_0_0_0_x_z_0_xy = buffer_2000_ppsd[13];

    auto g_xx_0_0_0_x_z_0_xz = buffer_2000_ppsd[14];

    auto g_xx_0_0_0_x_z_0_yy = buffer_2000_ppsd[15];

    auto g_xx_0_0_0_x_z_0_yz = buffer_2000_ppsd[16];

    auto g_xx_0_0_0_x_z_0_zz = buffer_2000_ppsd[17];

    auto g_xx_0_0_0_y_x_0_xx = buffer_2000_ppsd[18];

    auto g_xx_0_0_0_y_x_0_xy = buffer_2000_ppsd[19];

    auto g_xx_0_0_0_y_x_0_xz = buffer_2000_ppsd[20];

    auto g_xx_0_0_0_y_x_0_yy = buffer_2000_ppsd[21];

    auto g_xx_0_0_0_y_x_0_yz = buffer_2000_ppsd[22];

    auto g_xx_0_0_0_y_x_0_zz = buffer_2000_ppsd[23];

    auto g_xx_0_0_0_y_y_0_xx = buffer_2000_ppsd[24];

    auto g_xx_0_0_0_y_y_0_xy = buffer_2000_ppsd[25];

    auto g_xx_0_0_0_y_y_0_xz = buffer_2000_ppsd[26];

    auto g_xx_0_0_0_y_y_0_yy = buffer_2000_ppsd[27];

    auto g_xx_0_0_0_y_y_0_yz = buffer_2000_ppsd[28];

    auto g_xx_0_0_0_y_y_0_zz = buffer_2000_ppsd[29];

    auto g_xx_0_0_0_y_z_0_xx = buffer_2000_ppsd[30];

    auto g_xx_0_0_0_y_z_0_xy = buffer_2000_ppsd[31];

    auto g_xx_0_0_0_y_z_0_xz = buffer_2000_ppsd[32];

    auto g_xx_0_0_0_y_z_0_yy = buffer_2000_ppsd[33];

    auto g_xx_0_0_0_y_z_0_yz = buffer_2000_ppsd[34];

    auto g_xx_0_0_0_y_z_0_zz = buffer_2000_ppsd[35];

    auto g_xx_0_0_0_z_x_0_xx = buffer_2000_ppsd[36];

    auto g_xx_0_0_0_z_x_0_xy = buffer_2000_ppsd[37];

    auto g_xx_0_0_0_z_x_0_xz = buffer_2000_ppsd[38];

    auto g_xx_0_0_0_z_x_0_yy = buffer_2000_ppsd[39];

    auto g_xx_0_0_0_z_x_0_yz = buffer_2000_ppsd[40];

    auto g_xx_0_0_0_z_x_0_zz = buffer_2000_ppsd[41];

    auto g_xx_0_0_0_z_y_0_xx = buffer_2000_ppsd[42];

    auto g_xx_0_0_0_z_y_0_xy = buffer_2000_ppsd[43];

    auto g_xx_0_0_0_z_y_0_xz = buffer_2000_ppsd[44];

    auto g_xx_0_0_0_z_y_0_yy = buffer_2000_ppsd[45];

    auto g_xx_0_0_0_z_y_0_yz = buffer_2000_ppsd[46];

    auto g_xx_0_0_0_z_y_0_zz = buffer_2000_ppsd[47];

    auto g_xx_0_0_0_z_z_0_xx = buffer_2000_ppsd[48];

    auto g_xx_0_0_0_z_z_0_xy = buffer_2000_ppsd[49];

    auto g_xx_0_0_0_z_z_0_xz = buffer_2000_ppsd[50];

    auto g_xx_0_0_0_z_z_0_yy = buffer_2000_ppsd[51];

    auto g_xx_0_0_0_z_z_0_yz = buffer_2000_ppsd[52];

    auto g_xx_0_0_0_z_z_0_zz = buffer_2000_ppsd[53];

    auto g_xy_0_0_0_x_x_0_xx = buffer_2000_ppsd[54];

    auto g_xy_0_0_0_x_x_0_xy = buffer_2000_ppsd[55];

    auto g_xy_0_0_0_x_x_0_xz = buffer_2000_ppsd[56];

    auto g_xy_0_0_0_x_x_0_yy = buffer_2000_ppsd[57];

    auto g_xy_0_0_0_x_x_0_yz = buffer_2000_ppsd[58];

    auto g_xy_0_0_0_x_x_0_zz = buffer_2000_ppsd[59];

    auto g_xy_0_0_0_x_y_0_xx = buffer_2000_ppsd[60];

    auto g_xy_0_0_0_x_y_0_xy = buffer_2000_ppsd[61];

    auto g_xy_0_0_0_x_y_0_xz = buffer_2000_ppsd[62];

    auto g_xy_0_0_0_x_y_0_yy = buffer_2000_ppsd[63];

    auto g_xy_0_0_0_x_y_0_yz = buffer_2000_ppsd[64];

    auto g_xy_0_0_0_x_y_0_zz = buffer_2000_ppsd[65];

    auto g_xy_0_0_0_x_z_0_xx = buffer_2000_ppsd[66];

    auto g_xy_0_0_0_x_z_0_xy = buffer_2000_ppsd[67];

    auto g_xy_0_0_0_x_z_0_xz = buffer_2000_ppsd[68];

    auto g_xy_0_0_0_x_z_0_yy = buffer_2000_ppsd[69];

    auto g_xy_0_0_0_x_z_0_yz = buffer_2000_ppsd[70];

    auto g_xy_0_0_0_x_z_0_zz = buffer_2000_ppsd[71];

    auto g_xy_0_0_0_y_x_0_xx = buffer_2000_ppsd[72];

    auto g_xy_0_0_0_y_x_0_xy = buffer_2000_ppsd[73];

    auto g_xy_0_0_0_y_x_0_xz = buffer_2000_ppsd[74];

    auto g_xy_0_0_0_y_x_0_yy = buffer_2000_ppsd[75];

    auto g_xy_0_0_0_y_x_0_yz = buffer_2000_ppsd[76];

    auto g_xy_0_0_0_y_x_0_zz = buffer_2000_ppsd[77];

    auto g_xy_0_0_0_y_y_0_xx = buffer_2000_ppsd[78];

    auto g_xy_0_0_0_y_y_0_xy = buffer_2000_ppsd[79];

    auto g_xy_0_0_0_y_y_0_xz = buffer_2000_ppsd[80];

    auto g_xy_0_0_0_y_y_0_yy = buffer_2000_ppsd[81];

    auto g_xy_0_0_0_y_y_0_yz = buffer_2000_ppsd[82];

    auto g_xy_0_0_0_y_y_0_zz = buffer_2000_ppsd[83];

    auto g_xy_0_0_0_y_z_0_xx = buffer_2000_ppsd[84];

    auto g_xy_0_0_0_y_z_0_xy = buffer_2000_ppsd[85];

    auto g_xy_0_0_0_y_z_0_xz = buffer_2000_ppsd[86];

    auto g_xy_0_0_0_y_z_0_yy = buffer_2000_ppsd[87];

    auto g_xy_0_0_0_y_z_0_yz = buffer_2000_ppsd[88];

    auto g_xy_0_0_0_y_z_0_zz = buffer_2000_ppsd[89];

    auto g_xy_0_0_0_z_x_0_xx = buffer_2000_ppsd[90];

    auto g_xy_0_0_0_z_x_0_xy = buffer_2000_ppsd[91];

    auto g_xy_0_0_0_z_x_0_xz = buffer_2000_ppsd[92];

    auto g_xy_0_0_0_z_x_0_yy = buffer_2000_ppsd[93];

    auto g_xy_0_0_0_z_x_0_yz = buffer_2000_ppsd[94];

    auto g_xy_0_0_0_z_x_0_zz = buffer_2000_ppsd[95];

    auto g_xy_0_0_0_z_y_0_xx = buffer_2000_ppsd[96];

    auto g_xy_0_0_0_z_y_0_xy = buffer_2000_ppsd[97];

    auto g_xy_0_0_0_z_y_0_xz = buffer_2000_ppsd[98];

    auto g_xy_0_0_0_z_y_0_yy = buffer_2000_ppsd[99];

    auto g_xy_0_0_0_z_y_0_yz = buffer_2000_ppsd[100];

    auto g_xy_0_0_0_z_y_0_zz = buffer_2000_ppsd[101];

    auto g_xy_0_0_0_z_z_0_xx = buffer_2000_ppsd[102];

    auto g_xy_0_0_0_z_z_0_xy = buffer_2000_ppsd[103];

    auto g_xy_0_0_0_z_z_0_xz = buffer_2000_ppsd[104];

    auto g_xy_0_0_0_z_z_0_yy = buffer_2000_ppsd[105];

    auto g_xy_0_0_0_z_z_0_yz = buffer_2000_ppsd[106];

    auto g_xy_0_0_0_z_z_0_zz = buffer_2000_ppsd[107];

    auto g_xz_0_0_0_x_x_0_xx = buffer_2000_ppsd[108];

    auto g_xz_0_0_0_x_x_0_xy = buffer_2000_ppsd[109];

    auto g_xz_0_0_0_x_x_0_xz = buffer_2000_ppsd[110];

    auto g_xz_0_0_0_x_x_0_yy = buffer_2000_ppsd[111];

    auto g_xz_0_0_0_x_x_0_yz = buffer_2000_ppsd[112];

    auto g_xz_0_0_0_x_x_0_zz = buffer_2000_ppsd[113];

    auto g_xz_0_0_0_x_y_0_xx = buffer_2000_ppsd[114];

    auto g_xz_0_0_0_x_y_0_xy = buffer_2000_ppsd[115];

    auto g_xz_0_0_0_x_y_0_xz = buffer_2000_ppsd[116];

    auto g_xz_0_0_0_x_y_0_yy = buffer_2000_ppsd[117];

    auto g_xz_0_0_0_x_y_0_yz = buffer_2000_ppsd[118];

    auto g_xz_0_0_0_x_y_0_zz = buffer_2000_ppsd[119];

    auto g_xz_0_0_0_x_z_0_xx = buffer_2000_ppsd[120];

    auto g_xz_0_0_0_x_z_0_xy = buffer_2000_ppsd[121];

    auto g_xz_0_0_0_x_z_0_xz = buffer_2000_ppsd[122];

    auto g_xz_0_0_0_x_z_0_yy = buffer_2000_ppsd[123];

    auto g_xz_0_0_0_x_z_0_yz = buffer_2000_ppsd[124];

    auto g_xz_0_0_0_x_z_0_zz = buffer_2000_ppsd[125];

    auto g_xz_0_0_0_y_x_0_xx = buffer_2000_ppsd[126];

    auto g_xz_0_0_0_y_x_0_xy = buffer_2000_ppsd[127];

    auto g_xz_0_0_0_y_x_0_xz = buffer_2000_ppsd[128];

    auto g_xz_0_0_0_y_x_0_yy = buffer_2000_ppsd[129];

    auto g_xz_0_0_0_y_x_0_yz = buffer_2000_ppsd[130];

    auto g_xz_0_0_0_y_x_0_zz = buffer_2000_ppsd[131];

    auto g_xz_0_0_0_y_y_0_xx = buffer_2000_ppsd[132];

    auto g_xz_0_0_0_y_y_0_xy = buffer_2000_ppsd[133];

    auto g_xz_0_0_0_y_y_0_xz = buffer_2000_ppsd[134];

    auto g_xz_0_0_0_y_y_0_yy = buffer_2000_ppsd[135];

    auto g_xz_0_0_0_y_y_0_yz = buffer_2000_ppsd[136];

    auto g_xz_0_0_0_y_y_0_zz = buffer_2000_ppsd[137];

    auto g_xz_0_0_0_y_z_0_xx = buffer_2000_ppsd[138];

    auto g_xz_0_0_0_y_z_0_xy = buffer_2000_ppsd[139];

    auto g_xz_0_0_0_y_z_0_xz = buffer_2000_ppsd[140];

    auto g_xz_0_0_0_y_z_0_yy = buffer_2000_ppsd[141];

    auto g_xz_0_0_0_y_z_0_yz = buffer_2000_ppsd[142];

    auto g_xz_0_0_0_y_z_0_zz = buffer_2000_ppsd[143];

    auto g_xz_0_0_0_z_x_0_xx = buffer_2000_ppsd[144];

    auto g_xz_0_0_0_z_x_0_xy = buffer_2000_ppsd[145];

    auto g_xz_0_0_0_z_x_0_xz = buffer_2000_ppsd[146];

    auto g_xz_0_0_0_z_x_0_yy = buffer_2000_ppsd[147];

    auto g_xz_0_0_0_z_x_0_yz = buffer_2000_ppsd[148];

    auto g_xz_0_0_0_z_x_0_zz = buffer_2000_ppsd[149];

    auto g_xz_0_0_0_z_y_0_xx = buffer_2000_ppsd[150];

    auto g_xz_0_0_0_z_y_0_xy = buffer_2000_ppsd[151];

    auto g_xz_0_0_0_z_y_0_xz = buffer_2000_ppsd[152];

    auto g_xz_0_0_0_z_y_0_yy = buffer_2000_ppsd[153];

    auto g_xz_0_0_0_z_y_0_yz = buffer_2000_ppsd[154];

    auto g_xz_0_0_0_z_y_0_zz = buffer_2000_ppsd[155];

    auto g_xz_0_0_0_z_z_0_xx = buffer_2000_ppsd[156];

    auto g_xz_0_0_0_z_z_0_xy = buffer_2000_ppsd[157];

    auto g_xz_0_0_0_z_z_0_xz = buffer_2000_ppsd[158];

    auto g_xz_0_0_0_z_z_0_yy = buffer_2000_ppsd[159];

    auto g_xz_0_0_0_z_z_0_yz = buffer_2000_ppsd[160];

    auto g_xz_0_0_0_z_z_0_zz = buffer_2000_ppsd[161];

    auto g_yy_0_0_0_x_x_0_xx = buffer_2000_ppsd[162];

    auto g_yy_0_0_0_x_x_0_xy = buffer_2000_ppsd[163];

    auto g_yy_0_0_0_x_x_0_xz = buffer_2000_ppsd[164];

    auto g_yy_0_0_0_x_x_0_yy = buffer_2000_ppsd[165];

    auto g_yy_0_0_0_x_x_0_yz = buffer_2000_ppsd[166];

    auto g_yy_0_0_0_x_x_0_zz = buffer_2000_ppsd[167];

    auto g_yy_0_0_0_x_y_0_xx = buffer_2000_ppsd[168];

    auto g_yy_0_0_0_x_y_0_xy = buffer_2000_ppsd[169];

    auto g_yy_0_0_0_x_y_0_xz = buffer_2000_ppsd[170];

    auto g_yy_0_0_0_x_y_0_yy = buffer_2000_ppsd[171];

    auto g_yy_0_0_0_x_y_0_yz = buffer_2000_ppsd[172];

    auto g_yy_0_0_0_x_y_0_zz = buffer_2000_ppsd[173];

    auto g_yy_0_0_0_x_z_0_xx = buffer_2000_ppsd[174];

    auto g_yy_0_0_0_x_z_0_xy = buffer_2000_ppsd[175];

    auto g_yy_0_0_0_x_z_0_xz = buffer_2000_ppsd[176];

    auto g_yy_0_0_0_x_z_0_yy = buffer_2000_ppsd[177];

    auto g_yy_0_0_0_x_z_0_yz = buffer_2000_ppsd[178];

    auto g_yy_0_0_0_x_z_0_zz = buffer_2000_ppsd[179];

    auto g_yy_0_0_0_y_x_0_xx = buffer_2000_ppsd[180];

    auto g_yy_0_0_0_y_x_0_xy = buffer_2000_ppsd[181];

    auto g_yy_0_0_0_y_x_0_xz = buffer_2000_ppsd[182];

    auto g_yy_0_0_0_y_x_0_yy = buffer_2000_ppsd[183];

    auto g_yy_0_0_0_y_x_0_yz = buffer_2000_ppsd[184];

    auto g_yy_0_0_0_y_x_0_zz = buffer_2000_ppsd[185];

    auto g_yy_0_0_0_y_y_0_xx = buffer_2000_ppsd[186];

    auto g_yy_0_0_0_y_y_0_xy = buffer_2000_ppsd[187];

    auto g_yy_0_0_0_y_y_0_xz = buffer_2000_ppsd[188];

    auto g_yy_0_0_0_y_y_0_yy = buffer_2000_ppsd[189];

    auto g_yy_0_0_0_y_y_0_yz = buffer_2000_ppsd[190];

    auto g_yy_0_0_0_y_y_0_zz = buffer_2000_ppsd[191];

    auto g_yy_0_0_0_y_z_0_xx = buffer_2000_ppsd[192];

    auto g_yy_0_0_0_y_z_0_xy = buffer_2000_ppsd[193];

    auto g_yy_0_0_0_y_z_0_xz = buffer_2000_ppsd[194];

    auto g_yy_0_0_0_y_z_0_yy = buffer_2000_ppsd[195];

    auto g_yy_0_0_0_y_z_0_yz = buffer_2000_ppsd[196];

    auto g_yy_0_0_0_y_z_0_zz = buffer_2000_ppsd[197];

    auto g_yy_0_0_0_z_x_0_xx = buffer_2000_ppsd[198];

    auto g_yy_0_0_0_z_x_0_xy = buffer_2000_ppsd[199];

    auto g_yy_0_0_0_z_x_0_xz = buffer_2000_ppsd[200];

    auto g_yy_0_0_0_z_x_0_yy = buffer_2000_ppsd[201];

    auto g_yy_0_0_0_z_x_0_yz = buffer_2000_ppsd[202];

    auto g_yy_0_0_0_z_x_0_zz = buffer_2000_ppsd[203];

    auto g_yy_0_0_0_z_y_0_xx = buffer_2000_ppsd[204];

    auto g_yy_0_0_0_z_y_0_xy = buffer_2000_ppsd[205];

    auto g_yy_0_0_0_z_y_0_xz = buffer_2000_ppsd[206];

    auto g_yy_0_0_0_z_y_0_yy = buffer_2000_ppsd[207];

    auto g_yy_0_0_0_z_y_0_yz = buffer_2000_ppsd[208];

    auto g_yy_0_0_0_z_y_0_zz = buffer_2000_ppsd[209];

    auto g_yy_0_0_0_z_z_0_xx = buffer_2000_ppsd[210];

    auto g_yy_0_0_0_z_z_0_xy = buffer_2000_ppsd[211];

    auto g_yy_0_0_0_z_z_0_xz = buffer_2000_ppsd[212];

    auto g_yy_0_0_0_z_z_0_yy = buffer_2000_ppsd[213];

    auto g_yy_0_0_0_z_z_0_yz = buffer_2000_ppsd[214];

    auto g_yy_0_0_0_z_z_0_zz = buffer_2000_ppsd[215];

    auto g_yz_0_0_0_x_x_0_xx = buffer_2000_ppsd[216];

    auto g_yz_0_0_0_x_x_0_xy = buffer_2000_ppsd[217];

    auto g_yz_0_0_0_x_x_0_xz = buffer_2000_ppsd[218];

    auto g_yz_0_0_0_x_x_0_yy = buffer_2000_ppsd[219];

    auto g_yz_0_0_0_x_x_0_yz = buffer_2000_ppsd[220];

    auto g_yz_0_0_0_x_x_0_zz = buffer_2000_ppsd[221];

    auto g_yz_0_0_0_x_y_0_xx = buffer_2000_ppsd[222];

    auto g_yz_0_0_0_x_y_0_xy = buffer_2000_ppsd[223];

    auto g_yz_0_0_0_x_y_0_xz = buffer_2000_ppsd[224];

    auto g_yz_0_0_0_x_y_0_yy = buffer_2000_ppsd[225];

    auto g_yz_0_0_0_x_y_0_yz = buffer_2000_ppsd[226];

    auto g_yz_0_0_0_x_y_0_zz = buffer_2000_ppsd[227];

    auto g_yz_0_0_0_x_z_0_xx = buffer_2000_ppsd[228];

    auto g_yz_0_0_0_x_z_0_xy = buffer_2000_ppsd[229];

    auto g_yz_0_0_0_x_z_0_xz = buffer_2000_ppsd[230];

    auto g_yz_0_0_0_x_z_0_yy = buffer_2000_ppsd[231];

    auto g_yz_0_0_0_x_z_0_yz = buffer_2000_ppsd[232];

    auto g_yz_0_0_0_x_z_0_zz = buffer_2000_ppsd[233];

    auto g_yz_0_0_0_y_x_0_xx = buffer_2000_ppsd[234];

    auto g_yz_0_0_0_y_x_0_xy = buffer_2000_ppsd[235];

    auto g_yz_0_0_0_y_x_0_xz = buffer_2000_ppsd[236];

    auto g_yz_0_0_0_y_x_0_yy = buffer_2000_ppsd[237];

    auto g_yz_0_0_0_y_x_0_yz = buffer_2000_ppsd[238];

    auto g_yz_0_0_0_y_x_0_zz = buffer_2000_ppsd[239];

    auto g_yz_0_0_0_y_y_0_xx = buffer_2000_ppsd[240];

    auto g_yz_0_0_0_y_y_0_xy = buffer_2000_ppsd[241];

    auto g_yz_0_0_0_y_y_0_xz = buffer_2000_ppsd[242];

    auto g_yz_0_0_0_y_y_0_yy = buffer_2000_ppsd[243];

    auto g_yz_0_0_0_y_y_0_yz = buffer_2000_ppsd[244];

    auto g_yz_0_0_0_y_y_0_zz = buffer_2000_ppsd[245];

    auto g_yz_0_0_0_y_z_0_xx = buffer_2000_ppsd[246];

    auto g_yz_0_0_0_y_z_0_xy = buffer_2000_ppsd[247];

    auto g_yz_0_0_0_y_z_0_xz = buffer_2000_ppsd[248];

    auto g_yz_0_0_0_y_z_0_yy = buffer_2000_ppsd[249];

    auto g_yz_0_0_0_y_z_0_yz = buffer_2000_ppsd[250];

    auto g_yz_0_0_0_y_z_0_zz = buffer_2000_ppsd[251];

    auto g_yz_0_0_0_z_x_0_xx = buffer_2000_ppsd[252];

    auto g_yz_0_0_0_z_x_0_xy = buffer_2000_ppsd[253];

    auto g_yz_0_0_0_z_x_0_xz = buffer_2000_ppsd[254];

    auto g_yz_0_0_0_z_x_0_yy = buffer_2000_ppsd[255];

    auto g_yz_0_0_0_z_x_0_yz = buffer_2000_ppsd[256];

    auto g_yz_0_0_0_z_x_0_zz = buffer_2000_ppsd[257];

    auto g_yz_0_0_0_z_y_0_xx = buffer_2000_ppsd[258];

    auto g_yz_0_0_0_z_y_0_xy = buffer_2000_ppsd[259];

    auto g_yz_0_0_0_z_y_0_xz = buffer_2000_ppsd[260];

    auto g_yz_0_0_0_z_y_0_yy = buffer_2000_ppsd[261];

    auto g_yz_0_0_0_z_y_0_yz = buffer_2000_ppsd[262];

    auto g_yz_0_0_0_z_y_0_zz = buffer_2000_ppsd[263];

    auto g_yz_0_0_0_z_z_0_xx = buffer_2000_ppsd[264];

    auto g_yz_0_0_0_z_z_0_xy = buffer_2000_ppsd[265];

    auto g_yz_0_0_0_z_z_0_xz = buffer_2000_ppsd[266];

    auto g_yz_0_0_0_z_z_0_yy = buffer_2000_ppsd[267];

    auto g_yz_0_0_0_z_z_0_yz = buffer_2000_ppsd[268];

    auto g_yz_0_0_0_z_z_0_zz = buffer_2000_ppsd[269];

    auto g_zz_0_0_0_x_x_0_xx = buffer_2000_ppsd[270];

    auto g_zz_0_0_0_x_x_0_xy = buffer_2000_ppsd[271];

    auto g_zz_0_0_0_x_x_0_xz = buffer_2000_ppsd[272];

    auto g_zz_0_0_0_x_x_0_yy = buffer_2000_ppsd[273];

    auto g_zz_0_0_0_x_x_0_yz = buffer_2000_ppsd[274];

    auto g_zz_0_0_0_x_x_0_zz = buffer_2000_ppsd[275];

    auto g_zz_0_0_0_x_y_0_xx = buffer_2000_ppsd[276];

    auto g_zz_0_0_0_x_y_0_xy = buffer_2000_ppsd[277];

    auto g_zz_0_0_0_x_y_0_xz = buffer_2000_ppsd[278];

    auto g_zz_0_0_0_x_y_0_yy = buffer_2000_ppsd[279];

    auto g_zz_0_0_0_x_y_0_yz = buffer_2000_ppsd[280];

    auto g_zz_0_0_0_x_y_0_zz = buffer_2000_ppsd[281];

    auto g_zz_0_0_0_x_z_0_xx = buffer_2000_ppsd[282];

    auto g_zz_0_0_0_x_z_0_xy = buffer_2000_ppsd[283];

    auto g_zz_0_0_0_x_z_0_xz = buffer_2000_ppsd[284];

    auto g_zz_0_0_0_x_z_0_yy = buffer_2000_ppsd[285];

    auto g_zz_0_0_0_x_z_0_yz = buffer_2000_ppsd[286];

    auto g_zz_0_0_0_x_z_0_zz = buffer_2000_ppsd[287];

    auto g_zz_0_0_0_y_x_0_xx = buffer_2000_ppsd[288];

    auto g_zz_0_0_0_y_x_0_xy = buffer_2000_ppsd[289];

    auto g_zz_0_0_0_y_x_0_xz = buffer_2000_ppsd[290];

    auto g_zz_0_0_0_y_x_0_yy = buffer_2000_ppsd[291];

    auto g_zz_0_0_0_y_x_0_yz = buffer_2000_ppsd[292];

    auto g_zz_0_0_0_y_x_0_zz = buffer_2000_ppsd[293];

    auto g_zz_0_0_0_y_y_0_xx = buffer_2000_ppsd[294];

    auto g_zz_0_0_0_y_y_0_xy = buffer_2000_ppsd[295];

    auto g_zz_0_0_0_y_y_0_xz = buffer_2000_ppsd[296];

    auto g_zz_0_0_0_y_y_0_yy = buffer_2000_ppsd[297];

    auto g_zz_0_0_0_y_y_0_yz = buffer_2000_ppsd[298];

    auto g_zz_0_0_0_y_y_0_zz = buffer_2000_ppsd[299];

    auto g_zz_0_0_0_y_z_0_xx = buffer_2000_ppsd[300];

    auto g_zz_0_0_0_y_z_0_xy = buffer_2000_ppsd[301];

    auto g_zz_0_0_0_y_z_0_xz = buffer_2000_ppsd[302];

    auto g_zz_0_0_0_y_z_0_yy = buffer_2000_ppsd[303];

    auto g_zz_0_0_0_y_z_0_yz = buffer_2000_ppsd[304];

    auto g_zz_0_0_0_y_z_0_zz = buffer_2000_ppsd[305];

    auto g_zz_0_0_0_z_x_0_xx = buffer_2000_ppsd[306];

    auto g_zz_0_0_0_z_x_0_xy = buffer_2000_ppsd[307];

    auto g_zz_0_0_0_z_x_0_xz = buffer_2000_ppsd[308];

    auto g_zz_0_0_0_z_x_0_yy = buffer_2000_ppsd[309];

    auto g_zz_0_0_0_z_x_0_yz = buffer_2000_ppsd[310];

    auto g_zz_0_0_0_z_x_0_zz = buffer_2000_ppsd[311];

    auto g_zz_0_0_0_z_y_0_xx = buffer_2000_ppsd[312];

    auto g_zz_0_0_0_z_y_0_xy = buffer_2000_ppsd[313];

    auto g_zz_0_0_0_z_y_0_xz = buffer_2000_ppsd[314];

    auto g_zz_0_0_0_z_y_0_yy = buffer_2000_ppsd[315];

    auto g_zz_0_0_0_z_y_0_yz = buffer_2000_ppsd[316];

    auto g_zz_0_0_0_z_y_0_zz = buffer_2000_ppsd[317];

    auto g_zz_0_0_0_z_z_0_xx = buffer_2000_ppsd[318];

    auto g_zz_0_0_0_z_z_0_xy = buffer_2000_ppsd[319];

    auto g_zz_0_0_0_z_z_0_xz = buffer_2000_ppsd[320];

    auto g_zz_0_0_0_z_z_0_yy = buffer_2000_ppsd[321];

    auto g_zz_0_0_0_z_z_0_yz = buffer_2000_ppsd[322];

    auto g_zz_0_0_0_z_z_0_zz = buffer_2000_ppsd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_xx_0_0_0_x_x_0_xx, g_xx_0_0_0_x_x_0_xy, g_xx_0_0_0_x_x_0_xz, g_xx_0_0_0_x_x_0_yy, g_xx_0_0_0_x_x_0_yz, g_xx_0_0_0_x_x_0_zz, g_xxx_x_0_xx, g_xxx_x_0_xy, g_xxx_x_0_xz, g_xxx_x_0_yy, g_xxx_x_0_yz, g_xxx_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_0_xx[i] = -6.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_xxx_x_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_xy[i] = -6.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_xxx_x_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_xz[i] = -6.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_xxx_x_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_yy[i] = -6.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_xxx_x_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_yz[i] = -6.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_xxx_x_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_zz[i] = -6.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_xxx_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xx_0_0_0_x_y_0_xx, g_xx_0_0_0_x_y_0_xy, g_xx_0_0_0_x_y_0_xz, g_xx_0_0_0_x_y_0_yy, g_xx_0_0_0_x_y_0_yz, g_xx_0_0_0_x_y_0_zz, g_xxx_y_0_xx, g_xxx_y_0_xy, g_xxx_y_0_xz, g_xxx_y_0_yy, g_xxx_y_0_yz, g_xxx_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_0_xx[i] = -6.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_xxx_y_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_xy[i] = -6.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_xxx_y_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_xz[i] = -6.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_xxx_y_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_yy[i] = -6.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_xxx_y_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_yz[i] = -6.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_xxx_y_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_zz[i] = -6.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_xxx_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xx_0_0_0_x_z_0_xx, g_xx_0_0_0_x_z_0_xy, g_xx_0_0_0_x_z_0_xz, g_xx_0_0_0_x_z_0_yy, g_xx_0_0_0_x_z_0_yz, g_xx_0_0_0_x_z_0_zz, g_xxx_z_0_xx, g_xxx_z_0_xy, g_xxx_z_0_xz, g_xxx_z_0_yy, g_xxx_z_0_yz, g_xxx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_0_xx[i] = -6.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_xxx_z_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_xy[i] = -6.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_xxx_z_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_xz[i] = -6.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_xxx_z_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_yy[i] = -6.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_xxx_z_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_yz[i] = -6.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_xxx_z_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_zz[i] = -6.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_xxx_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_0_xx, g_xx_0_0_0_y_x_0_xy, g_xx_0_0_0_y_x_0_xz, g_xx_0_0_0_y_x_0_yy, g_xx_0_0_0_y_x_0_yz, g_xx_0_0_0_y_x_0_zz, g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_xxy_x_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_xxy_x_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_xxy_x_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_xxy_x_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_xxy_x_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_xxy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_0_xx, g_xx_0_0_0_y_y_0_xy, g_xx_0_0_0_y_y_0_xz, g_xx_0_0_0_y_y_0_yy, g_xx_0_0_0_y_y_0_yz, g_xx_0_0_0_y_y_0_zz, g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_xxy_y_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_xxy_y_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_xxy_y_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_xxy_y_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_xxy_y_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_xxy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_0_xx, g_xx_0_0_0_y_z_0_xy, g_xx_0_0_0_y_z_0_xz, g_xx_0_0_0_y_z_0_yy, g_xx_0_0_0_y_z_0_yz, g_xx_0_0_0_y_z_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_xxy_z_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_xxy_z_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_xxy_z_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_xxy_z_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_xxy_z_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_xxy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_0_xx, g_xx_0_0_0_z_x_0_xy, g_xx_0_0_0_z_x_0_xz, g_xx_0_0_0_z_x_0_yy, g_xx_0_0_0_z_x_0_yz, g_xx_0_0_0_z_x_0_zz, g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_xxz_x_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_xxz_x_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_xxz_x_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_xxz_x_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_xxz_x_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_xxz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_0_xx, g_xx_0_0_0_z_y_0_xy, g_xx_0_0_0_z_y_0_xz, g_xx_0_0_0_z_y_0_yy, g_xx_0_0_0_z_y_0_yz, g_xx_0_0_0_z_y_0_zz, g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_xxz_y_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_xxz_y_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_xxz_y_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_xxz_y_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_xxz_y_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_xxz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_0_xx, g_xx_0_0_0_z_z_0_xy, g_xx_0_0_0_z_z_0_xz, g_xx_0_0_0_z_z_0_yy, g_xx_0_0_0_z_z_0_yz, g_xx_0_0_0_z_z_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_xxz_z_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_xxz_z_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_xxz_z_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_xxz_z_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_xxz_z_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_xxz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xy_0_0_0_x_x_0_xx, g_xy_0_0_0_x_x_0_xy, g_xy_0_0_0_x_x_0_xz, g_xy_0_0_0_x_x_0_yy, g_xy_0_0_0_x_x_0_yz, g_xy_0_0_0_x_x_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_xxy_x_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_xxy_x_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_xxy_x_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_xxy_x_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_xxy_x_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_xxy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_xy_0_0_0_x_y_0_xx, g_xy_0_0_0_x_y_0_xy, g_xy_0_0_0_x_y_0_xz, g_xy_0_0_0_x_y_0_yy, g_xy_0_0_0_x_y_0_yz, g_xy_0_0_0_x_y_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_xxy_y_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_xxy_y_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_xxy_y_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_xxy_y_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_xxy_y_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_xxy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_xy_0_0_0_x_z_0_xx, g_xy_0_0_0_x_z_0_xy, g_xy_0_0_0_x_z_0_xz, g_xy_0_0_0_x_z_0_yy, g_xy_0_0_0_x_z_0_yz, g_xy_0_0_0_x_z_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_xxy_z_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_xxy_z_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_xxy_z_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_xxy_z_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_xxy_z_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_xxy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_xy_0_0_0_y_x_0_xx, g_xy_0_0_0_y_x_0_xy, g_xy_0_0_0_y_x_0_xz, g_xy_0_0_0_y_x_0_yy, g_xy_0_0_0_y_x_0_yz, g_xy_0_0_0_y_x_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_xyy_x_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_xyy_x_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_xyy_x_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_xyy_x_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_xyy_x_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_xyy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xy_0_0_0_y_y_0_xx, g_xy_0_0_0_y_y_0_xy, g_xy_0_0_0_y_y_0_xz, g_xy_0_0_0_y_y_0_yy, g_xy_0_0_0_y_y_0_yz, g_xy_0_0_0_y_y_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_xyy_y_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_xyy_y_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_xyy_y_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_xyy_y_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_xyy_y_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_xyy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xy_0_0_0_y_z_0_xx, g_xy_0_0_0_y_z_0_xy, g_xy_0_0_0_y_z_0_xz, g_xy_0_0_0_y_z_0_yy, g_xy_0_0_0_y_z_0_yz, g_xy_0_0_0_y_z_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_xyy_z_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_xyy_z_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_xyy_z_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_xyy_z_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_xyy_z_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_xyy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_0_xx, g_xy_0_0_0_z_x_0_xy, g_xy_0_0_0_z_x_0_xz, g_xy_0_0_0_z_x_0_yy, g_xy_0_0_0_z_x_0_yz, g_xy_0_0_0_z_x_0_zz, g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_0_xx[i] = 4.0 * g_xyz_x_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_xy[i] = 4.0 * g_xyz_x_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_xz[i] = 4.0 * g_xyz_x_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_yy[i] = 4.0 * g_xyz_x_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_yz[i] = 4.0 * g_xyz_x_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_zz[i] = 4.0 * g_xyz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_0_xx, g_xy_0_0_0_z_y_0_xy, g_xy_0_0_0_z_y_0_xz, g_xy_0_0_0_z_y_0_yy, g_xy_0_0_0_z_y_0_yz, g_xy_0_0_0_z_y_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_0_xx[i] = 4.0 * g_xyz_y_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_xy[i] = 4.0 * g_xyz_y_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_xz[i] = 4.0 * g_xyz_y_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_yy[i] = 4.0 * g_xyz_y_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_yz[i] = 4.0 * g_xyz_y_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_zz[i] = 4.0 * g_xyz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_0_xx, g_xy_0_0_0_z_z_0_xy, g_xy_0_0_0_z_z_0_xz, g_xy_0_0_0_z_z_0_yy, g_xy_0_0_0_z_z_0_yz, g_xy_0_0_0_z_z_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_0_xx[i] = 4.0 * g_xyz_z_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_xy[i] = 4.0 * g_xyz_z_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_xz[i] = 4.0 * g_xyz_z_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_yy[i] = 4.0 * g_xyz_z_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_yz[i] = 4.0 * g_xyz_z_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_zz[i] = 4.0 * g_xyz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xz_0_0_0_x_x_0_xx, g_xz_0_0_0_x_x_0_xy, g_xz_0_0_0_x_x_0_xz, g_xz_0_0_0_x_x_0_yy, g_xz_0_0_0_x_x_0_yz, g_xz_0_0_0_x_x_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_xxz_x_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_xxz_x_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_xxz_x_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_xxz_x_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_xxz_x_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_xxz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_xz_0_0_0_x_y_0_xx, g_xz_0_0_0_x_y_0_xy, g_xz_0_0_0_x_y_0_xz, g_xz_0_0_0_x_y_0_yy, g_xz_0_0_0_x_y_0_yz, g_xz_0_0_0_x_y_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_xxz_y_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_xxz_y_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_xxz_y_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_xxz_y_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_xxz_y_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_xxz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_xz_0_0_0_x_z_0_xx, g_xz_0_0_0_x_z_0_xy, g_xz_0_0_0_x_z_0_xz, g_xz_0_0_0_x_z_0_yy, g_xz_0_0_0_x_z_0_yz, g_xz_0_0_0_x_z_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_xxz_z_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_xxz_z_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_xxz_z_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_xxz_z_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_xxz_z_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_xxz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xz_0_0_0_y_x_0_xx, g_xz_0_0_0_y_x_0_xy, g_xz_0_0_0_y_x_0_xz, g_xz_0_0_0_y_x_0_yy, g_xz_0_0_0_y_x_0_yz, g_xz_0_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_0_xx[i] = 4.0 * g_xyz_x_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_xy[i] = 4.0 * g_xyz_x_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_xz[i] = 4.0 * g_xyz_x_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_yy[i] = 4.0 * g_xyz_x_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_yz[i] = 4.0 * g_xyz_x_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_zz[i] = 4.0 * g_xyz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xz_0_0_0_y_y_0_xx, g_xz_0_0_0_y_y_0_xy, g_xz_0_0_0_y_y_0_xz, g_xz_0_0_0_y_y_0_yy, g_xz_0_0_0_y_y_0_yz, g_xz_0_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_0_xx[i] = 4.0 * g_xyz_y_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_xy[i] = 4.0 * g_xyz_y_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_xz[i] = 4.0 * g_xyz_y_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_yy[i] = 4.0 * g_xyz_y_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_yz[i] = 4.0 * g_xyz_y_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_zz[i] = 4.0 * g_xyz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_xz_0_0_0_y_z_0_xx, g_xz_0_0_0_y_z_0_xy, g_xz_0_0_0_y_z_0_xz, g_xz_0_0_0_y_z_0_yy, g_xz_0_0_0_y_z_0_yz, g_xz_0_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_0_xx[i] = 4.0 * g_xyz_z_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_xy[i] = 4.0 * g_xyz_z_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_xz[i] = 4.0 * g_xyz_z_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_yy[i] = 4.0 * g_xyz_z_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_yz[i] = 4.0 * g_xyz_z_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_zz[i] = 4.0 * g_xyz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_xz_0_0_0_z_x_0_xx, g_xz_0_0_0_z_x_0_xy, g_xz_0_0_0_z_x_0_xz, g_xz_0_0_0_z_x_0_yy, g_xz_0_0_0_z_x_0_yz, g_xz_0_0_0_z_x_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_xzz_x_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_xzz_x_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_xzz_x_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_xzz_x_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_xzz_x_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_xzz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xz_0_0_0_z_y_0_xx, g_xz_0_0_0_z_y_0_xy, g_xz_0_0_0_z_y_0_xz, g_xz_0_0_0_z_y_0_yy, g_xz_0_0_0_z_y_0_yz, g_xz_0_0_0_z_y_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_xzz_y_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_xzz_y_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_xzz_y_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_xzz_y_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_xzz_y_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_xzz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xz_0_0_0_z_z_0_xx, g_xz_0_0_0_z_z_0_xy, g_xz_0_0_0_z_z_0_xz, g_xz_0_0_0_z_z_0_yy, g_xz_0_0_0_z_z_0_yz, g_xz_0_0_0_z_z_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_xzz_z_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_xzz_z_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_xzz_z_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_xzz_z_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_xzz_z_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_xzz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_yy_0_0_0_x_x_0_xx, g_yy_0_0_0_x_x_0_xy, g_yy_0_0_0_x_x_0_xz, g_yy_0_0_0_x_x_0_yy, g_yy_0_0_0_x_x_0_yz, g_yy_0_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_xyy_x_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_xyy_x_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_xyy_x_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_xyy_x_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_xyy_x_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_xyy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_yy_0_0_0_x_y_0_xx, g_yy_0_0_0_x_y_0_xy, g_yy_0_0_0_x_y_0_xz, g_yy_0_0_0_x_y_0_yy, g_yy_0_0_0_x_y_0_yz, g_yy_0_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_xyy_y_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_xyy_y_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_xyy_y_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_xyy_y_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_xyy_y_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_xyy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz, g_yy_0_0_0_x_z_0_xx, g_yy_0_0_0_x_z_0_xy, g_yy_0_0_0_x_z_0_xz, g_yy_0_0_0_x_z_0_yy, g_yy_0_0_0_x_z_0_yz, g_yy_0_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_xyy_z_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_xyy_z_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_xyy_z_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_xyy_z_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_xyy_z_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_xyy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_yy_0_0_0_y_x_0_xx, g_yy_0_0_0_y_x_0_xy, g_yy_0_0_0_y_x_0_xz, g_yy_0_0_0_y_x_0_yy, g_yy_0_0_0_y_x_0_yz, g_yy_0_0_0_y_x_0_zz, g_yyy_x_0_xx, g_yyy_x_0_xy, g_yyy_x_0_xz, g_yyy_x_0_yy, g_yyy_x_0_yz, g_yyy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_0_xx[i] = -6.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_yyy_x_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_xy[i] = -6.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_yyy_x_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_xz[i] = -6.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_yyy_x_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_yy[i] = -6.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_yyy_x_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_yz[i] = -6.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_yyy_x_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_zz[i] = -6.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_yyy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_yy_0_0_0_y_y_0_xx, g_yy_0_0_0_y_y_0_xy, g_yy_0_0_0_y_y_0_xz, g_yy_0_0_0_y_y_0_yy, g_yy_0_0_0_y_y_0_yz, g_yy_0_0_0_y_y_0_zz, g_yyy_y_0_xx, g_yyy_y_0_xy, g_yyy_y_0_xz, g_yyy_y_0_yy, g_yyy_y_0_yz, g_yyy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_0_xx[i] = -6.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_yyy_y_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_xy[i] = -6.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_yyy_y_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_xz[i] = -6.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_yyy_y_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_yy[i] = -6.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_yyy_y_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_yz[i] = -6.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_yyy_y_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_zz[i] = -6.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_yyy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yy_0_0_0_y_z_0_xx, g_yy_0_0_0_y_z_0_xy, g_yy_0_0_0_y_z_0_xz, g_yy_0_0_0_y_z_0_yy, g_yy_0_0_0_y_z_0_yz, g_yy_0_0_0_y_z_0_zz, g_yyy_z_0_xx, g_yyy_z_0_xy, g_yyy_z_0_xz, g_yyy_z_0_yy, g_yyy_z_0_yz, g_yyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_0_xx[i] = -6.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_yyy_z_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_xy[i] = -6.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_yyy_z_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_xz[i] = -6.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_yyy_z_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_yy[i] = -6.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_yyy_z_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_yz[i] = -6.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_yyy_z_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_zz[i] = -6.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_yyy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_0_xx, g_yy_0_0_0_z_x_0_xy, g_yy_0_0_0_z_x_0_xz, g_yy_0_0_0_z_x_0_yy, g_yy_0_0_0_z_x_0_yz, g_yy_0_0_0_z_x_0_zz, g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_yyz_x_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_yyz_x_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_yyz_x_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_yyz_x_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_yyz_x_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_yyz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_0_xx, g_yy_0_0_0_z_y_0_xy, g_yy_0_0_0_z_y_0_xz, g_yy_0_0_0_z_y_0_yy, g_yy_0_0_0_z_y_0_yz, g_yy_0_0_0_z_y_0_zz, g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_yyz_y_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_yyz_y_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_yyz_y_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_yyz_y_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_yyz_y_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_yyz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_0_xx, g_yy_0_0_0_z_z_0_xy, g_yy_0_0_0_z_z_0_xz, g_yy_0_0_0_z_z_0_yy, g_yy_0_0_0_z_z_0_yz, g_yy_0_0_0_z_z_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_yyz_z_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_yyz_z_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_yyz_z_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_yyz_z_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_yyz_z_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_yyz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_yz_0_0_0_x_x_0_xx, g_yz_0_0_0_x_x_0_xy, g_yz_0_0_0_x_x_0_xz, g_yz_0_0_0_x_x_0_yy, g_yz_0_0_0_x_x_0_yz, g_yz_0_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_0_xx[i] = 4.0 * g_xyz_x_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_xy[i] = 4.0 * g_xyz_x_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_xz[i] = 4.0 * g_xyz_x_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_yy[i] = 4.0 * g_xyz_x_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_yz[i] = 4.0 * g_xyz_x_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_zz[i] = 4.0 * g_xyz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_yz_0_0_0_x_y_0_xx, g_yz_0_0_0_x_y_0_xy, g_yz_0_0_0_x_y_0_xz, g_yz_0_0_0_x_y_0_yy, g_yz_0_0_0_x_y_0_yz, g_yz_0_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_0_xx[i] = 4.0 * g_xyz_y_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_xy[i] = 4.0 * g_xyz_y_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_xz[i] = 4.0 * g_xyz_y_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_yy[i] = 4.0 * g_xyz_y_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_yz[i] = 4.0 * g_xyz_y_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_zz[i] = 4.0 * g_xyz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_yz_0_0_0_x_z_0_xx, g_yz_0_0_0_x_z_0_xy, g_yz_0_0_0_x_z_0_xz, g_yz_0_0_0_x_z_0_yy, g_yz_0_0_0_x_z_0_yz, g_yz_0_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_0_xx[i] = 4.0 * g_xyz_z_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_xy[i] = 4.0 * g_xyz_z_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_xz[i] = 4.0 * g_xyz_z_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_yy[i] = 4.0 * g_xyz_z_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_yz[i] = 4.0 * g_xyz_z_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_zz[i] = 4.0 * g_xyz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yz_0_0_0_y_x_0_xx, g_yz_0_0_0_y_x_0_xy, g_yz_0_0_0_y_x_0_xz, g_yz_0_0_0_y_x_0_yy, g_yz_0_0_0_y_x_0_yz, g_yz_0_0_0_y_x_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_0_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_yyz_x_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_yyz_x_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_yyz_x_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_yyz_x_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_yyz_x_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_yyz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_yz_0_0_0_y_y_0_xx, g_yz_0_0_0_y_y_0_xy, g_yz_0_0_0_y_y_0_xz, g_yz_0_0_0_y_y_0_yy, g_yz_0_0_0_y_y_0_yz, g_yz_0_0_0_y_y_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_0_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_yyz_y_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_yyz_y_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_yyz_y_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_yyz_y_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_yyz_y_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_yyz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_yz_0_0_0_y_z_0_xx, g_yz_0_0_0_y_z_0_xy, g_yz_0_0_0_y_z_0_xz, g_yz_0_0_0_y_z_0_yy, g_yz_0_0_0_y_z_0_yz, g_yz_0_0_0_y_z_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_0_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_yyz_z_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_yyz_z_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_yyz_z_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_yyz_z_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_yyz_z_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_yyz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_yz_0_0_0_z_x_0_xx, g_yz_0_0_0_z_x_0_xy, g_yz_0_0_0_z_x_0_xz, g_yz_0_0_0_z_x_0_yy, g_yz_0_0_0_z_x_0_yz, g_yz_0_0_0_z_x_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_yzz_x_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_yzz_x_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_yzz_x_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_yzz_x_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_yzz_x_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_yzz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_yz_0_0_0_z_y_0_xx, g_yz_0_0_0_z_y_0_xy, g_yz_0_0_0_z_y_0_xz, g_yz_0_0_0_z_y_0_yy, g_yz_0_0_0_z_y_0_yz, g_yz_0_0_0_z_y_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_yzz_y_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_yzz_y_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_yzz_y_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_yzz_y_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_yzz_y_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_yzz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yz_0_0_0_z_z_0_xx, g_yz_0_0_0_z_z_0_xy, g_yz_0_0_0_z_z_0_xz, g_yz_0_0_0_z_z_0_yy, g_yz_0_0_0_z_z_0_yz, g_yz_0_0_0_z_z_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_yzz_z_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_yzz_z_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_yzz_z_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_yzz_z_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_yzz_z_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_yzz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_zz_0_0_0_x_x_0_xx, g_zz_0_0_0_x_x_0_xy, g_zz_0_0_0_x_x_0_xz, g_zz_0_0_0_x_x_0_yy, g_zz_0_0_0_x_x_0_yz, g_zz_0_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_0_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_xzz_x_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_xzz_x_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_xzz_x_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_xzz_x_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_xzz_x_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_xzz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_zz_0_0_0_x_y_0_xx, g_zz_0_0_0_x_y_0_xy, g_zz_0_0_0_x_y_0_xz, g_zz_0_0_0_x_y_0_yy, g_zz_0_0_0_x_y_0_yz, g_zz_0_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_0_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_xzz_y_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_xzz_y_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_xzz_y_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_xzz_y_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_xzz_y_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_xzz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz, g_zz_0_0_0_x_z_0_xx, g_zz_0_0_0_x_z_0_xy, g_zz_0_0_0_x_z_0_xz, g_zz_0_0_0_x_z_0_yy, g_zz_0_0_0_x_z_0_yz, g_zz_0_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_0_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_xzz_z_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_xzz_z_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_xzz_z_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_xzz_z_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_xzz_z_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_xzz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_zz_0_0_0_y_x_0_xx, g_zz_0_0_0_y_x_0_xy, g_zz_0_0_0_y_x_0_xz, g_zz_0_0_0_y_x_0_yy, g_zz_0_0_0_y_x_0_yz, g_zz_0_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_0_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_yzz_x_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_yzz_x_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_yzz_x_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_yzz_x_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_yzz_x_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_yzz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_zz_0_0_0_y_y_0_xx, g_zz_0_0_0_y_y_0_xy, g_zz_0_0_0_y_y_0_xz, g_zz_0_0_0_y_y_0_yy, g_zz_0_0_0_y_y_0_yz, g_zz_0_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_0_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_yzz_y_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_yzz_y_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_yzz_y_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_yzz_y_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_yzz_y_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_yzz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz, g_zz_0_0_0_y_z_0_xx, g_zz_0_0_0_y_z_0_xy, g_zz_0_0_0_y_z_0_xz, g_zz_0_0_0_y_z_0_yy, g_zz_0_0_0_y_z_0_yz, g_zz_0_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_0_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_yzz_z_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_yzz_z_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_yzz_z_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_yzz_z_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_yzz_z_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_yzz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_zz_0_0_0_z_x_0_xx, g_zz_0_0_0_z_x_0_xy, g_zz_0_0_0_z_x_0_xz, g_zz_0_0_0_z_x_0_yy, g_zz_0_0_0_z_x_0_yz, g_zz_0_0_0_z_x_0_zz, g_zzz_x_0_xx, g_zzz_x_0_xy, g_zzz_x_0_xz, g_zzz_x_0_yy, g_zzz_x_0_yz, g_zzz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_0_xx[i] = -6.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_zzz_x_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_xy[i] = -6.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_zzz_x_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_xz[i] = -6.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_zzz_x_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_yy[i] = -6.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_zzz_x_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_yz[i] = -6.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_zzz_x_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_zz[i] = -6.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_zzz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_zz_0_0_0_z_y_0_xx, g_zz_0_0_0_z_y_0_xy, g_zz_0_0_0_z_y_0_xz, g_zz_0_0_0_z_y_0_yy, g_zz_0_0_0_z_y_0_yz, g_zz_0_0_0_z_y_0_zz, g_zzz_y_0_xx, g_zzz_y_0_xy, g_zzz_y_0_xz, g_zzz_y_0_yy, g_zzz_y_0_yz, g_zzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_0_xx[i] = -6.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_zzz_y_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_xy[i] = -6.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_zzz_y_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_xz[i] = -6.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_zzz_y_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_yy[i] = -6.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_zzz_y_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_yz[i] = -6.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_zzz_y_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_zz[i] = -6.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_zzz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_zz_0_0_0_z_z_0_xx, g_zz_0_0_0_z_z_0_xy, g_zz_0_0_0_z_z_0_xz, g_zz_0_0_0_z_z_0_yy, g_zz_0_0_0_z_z_0_yz, g_zz_0_0_0_z_z_0_zz, g_zzz_z_0_xx, g_zzz_z_0_xy, g_zzz_z_0_xz, g_zzz_z_0_yy, g_zzz_z_0_yz, g_zzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_0_xx[i] = -6.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_zzz_z_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_xy[i] = -6.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_zzz_z_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_xz[i] = -6.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_zzz_z_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_yy[i] = -6.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_zzz_z_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_yz[i] = -6.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_zzz_z_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_zz[i] = -6.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_zzz_z_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


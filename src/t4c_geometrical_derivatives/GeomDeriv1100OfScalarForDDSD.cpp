#include "GeomDeriv1100OfScalarForDDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ddsd_0(CSimdArray<double>& buffer_1100_ddsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_pfsd,
                     const CSimdArray<double>& buffer_fpsd,
                     const CSimdArray<double>& buffer_ffsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ddsd.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_ffsd

    auto g_xxx_xxx_0_xx = buffer_ffsd[0];

    auto g_xxx_xxx_0_xy = buffer_ffsd[1];

    auto g_xxx_xxx_0_xz = buffer_ffsd[2];

    auto g_xxx_xxx_0_yy = buffer_ffsd[3];

    auto g_xxx_xxx_0_yz = buffer_ffsd[4];

    auto g_xxx_xxx_0_zz = buffer_ffsd[5];

    auto g_xxx_xxy_0_xx = buffer_ffsd[6];

    auto g_xxx_xxy_0_xy = buffer_ffsd[7];

    auto g_xxx_xxy_0_xz = buffer_ffsd[8];

    auto g_xxx_xxy_0_yy = buffer_ffsd[9];

    auto g_xxx_xxy_0_yz = buffer_ffsd[10];

    auto g_xxx_xxy_0_zz = buffer_ffsd[11];

    auto g_xxx_xxz_0_xx = buffer_ffsd[12];

    auto g_xxx_xxz_0_xy = buffer_ffsd[13];

    auto g_xxx_xxz_0_xz = buffer_ffsd[14];

    auto g_xxx_xxz_0_yy = buffer_ffsd[15];

    auto g_xxx_xxz_0_yz = buffer_ffsd[16];

    auto g_xxx_xxz_0_zz = buffer_ffsd[17];

    auto g_xxx_xyy_0_xx = buffer_ffsd[18];

    auto g_xxx_xyy_0_xy = buffer_ffsd[19];

    auto g_xxx_xyy_0_xz = buffer_ffsd[20];

    auto g_xxx_xyy_0_yy = buffer_ffsd[21];

    auto g_xxx_xyy_0_yz = buffer_ffsd[22];

    auto g_xxx_xyy_0_zz = buffer_ffsd[23];

    auto g_xxx_xyz_0_xx = buffer_ffsd[24];

    auto g_xxx_xyz_0_xy = buffer_ffsd[25];

    auto g_xxx_xyz_0_xz = buffer_ffsd[26];

    auto g_xxx_xyz_0_yy = buffer_ffsd[27];

    auto g_xxx_xyz_0_yz = buffer_ffsd[28];

    auto g_xxx_xyz_0_zz = buffer_ffsd[29];

    auto g_xxx_xzz_0_xx = buffer_ffsd[30];

    auto g_xxx_xzz_0_xy = buffer_ffsd[31];

    auto g_xxx_xzz_0_xz = buffer_ffsd[32];

    auto g_xxx_xzz_0_yy = buffer_ffsd[33];

    auto g_xxx_xzz_0_yz = buffer_ffsd[34];

    auto g_xxx_xzz_0_zz = buffer_ffsd[35];

    auto g_xxx_yyy_0_xx = buffer_ffsd[36];

    auto g_xxx_yyy_0_xy = buffer_ffsd[37];

    auto g_xxx_yyy_0_xz = buffer_ffsd[38];

    auto g_xxx_yyy_0_yy = buffer_ffsd[39];

    auto g_xxx_yyy_0_yz = buffer_ffsd[40];

    auto g_xxx_yyy_0_zz = buffer_ffsd[41];

    auto g_xxx_yyz_0_xx = buffer_ffsd[42];

    auto g_xxx_yyz_0_xy = buffer_ffsd[43];

    auto g_xxx_yyz_0_xz = buffer_ffsd[44];

    auto g_xxx_yyz_0_yy = buffer_ffsd[45];

    auto g_xxx_yyz_0_yz = buffer_ffsd[46];

    auto g_xxx_yyz_0_zz = buffer_ffsd[47];

    auto g_xxx_yzz_0_xx = buffer_ffsd[48];

    auto g_xxx_yzz_0_xy = buffer_ffsd[49];

    auto g_xxx_yzz_0_xz = buffer_ffsd[50];

    auto g_xxx_yzz_0_yy = buffer_ffsd[51];

    auto g_xxx_yzz_0_yz = buffer_ffsd[52];

    auto g_xxx_yzz_0_zz = buffer_ffsd[53];

    auto g_xxx_zzz_0_xx = buffer_ffsd[54];

    auto g_xxx_zzz_0_xy = buffer_ffsd[55];

    auto g_xxx_zzz_0_xz = buffer_ffsd[56];

    auto g_xxx_zzz_0_yy = buffer_ffsd[57];

    auto g_xxx_zzz_0_yz = buffer_ffsd[58];

    auto g_xxx_zzz_0_zz = buffer_ffsd[59];

    auto g_xxy_xxx_0_xx = buffer_ffsd[60];

    auto g_xxy_xxx_0_xy = buffer_ffsd[61];

    auto g_xxy_xxx_0_xz = buffer_ffsd[62];

    auto g_xxy_xxx_0_yy = buffer_ffsd[63];

    auto g_xxy_xxx_0_yz = buffer_ffsd[64];

    auto g_xxy_xxx_0_zz = buffer_ffsd[65];

    auto g_xxy_xxy_0_xx = buffer_ffsd[66];

    auto g_xxy_xxy_0_xy = buffer_ffsd[67];

    auto g_xxy_xxy_0_xz = buffer_ffsd[68];

    auto g_xxy_xxy_0_yy = buffer_ffsd[69];

    auto g_xxy_xxy_0_yz = buffer_ffsd[70];

    auto g_xxy_xxy_0_zz = buffer_ffsd[71];

    auto g_xxy_xxz_0_xx = buffer_ffsd[72];

    auto g_xxy_xxz_0_xy = buffer_ffsd[73];

    auto g_xxy_xxz_0_xz = buffer_ffsd[74];

    auto g_xxy_xxz_0_yy = buffer_ffsd[75];

    auto g_xxy_xxz_0_yz = buffer_ffsd[76];

    auto g_xxy_xxz_0_zz = buffer_ffsd[77];

    auto g_xxy_xyy_0_xx = buffer_ffsd[78];

    auto g_xxy_xyy_0_xy = buffer_ffsd[79];

    auto g_xxy_xyy_0_xz = buffer_ffsd[80];

    auto g_xxy_xyy_0_yy = buffer_ffsd[81];

    auto g_xxy_xyy_0_yz = buffer_ffsd[82];

    auto g_xxy_xyy_0_zz = buffer_ffsd[83];

    auto g_xxy_xyz_0_xx = buffer_ffsd[84];

    auto g_xxy_xyz_0_xy = buffer_ffsd[85];

    auto g_xxy_xyz_0_xz = buffer_ffsd[86];

    auto g_xxy_xyz_0_yy = buffer_ffsd[87];

    auto g_xxy_xyz_0_yz = buffer_ffsd[88];

    auto g_xxy_xyz_0_zz = buffer_ffsd[89];

    auto g_xxy_xzz_0_xx = buffer_ffsd[90];

    auto g_xxy_xzz_0_xy = buffer_ffsd[91];

    auto g_xxy_xzz_0_xz = buffer_ffsd[92];

    auto g_xxy_xzz_0_yy = buffer_ffsd[93];

    auto g_xxy_xzz_0_yz = buffer_ffsd[94];

    auto g_xxy_xzz_0_zz = buffer_ffsd[95];

    auto g_xxy_yyy_0_xx = buffer_ffsd[96];

    auto g_xxy_yyy_0_xy = buffer_ffsd[97];

    auto g_xxy_yyy_0_xz = buffer_ffsd[98];

    auto g_xxy_yyy_0_yy = buffer_ffsd[99];

    auto g_xxy_yyy_0_yz = buffer_ffsd[100];

    auto g_xxy_yyy_0_zz = buffer_ffsd[101];

    auto g_xxy_yyz_0_xx = buffer_ffsd[102];

    auto g_xxy_yyz_0_xy = buffer_ffsd[103];

    auto g_xxy_yyz_0_xz = buffer_ffsd[104];

    auto g_xxy_yyz_0_yy = buffer_ffsd[105];

    auto g_xxy_yyz_0_yz = buffer_ffsd[106];

    auto g_xxy_yyz_0_zz = buffer_ffsd[107];

    auto g_xxy_yzz_0_xx = buffer_ffsd[108];

    auto g_xxy_yzz_0_xy = buffer_ffsd[109];

    auto g_xxy_yzz_0_xz = buffer_ffsd[110];

    auto g_xxy_yzz_0_yy = buffer_ffsd[111];

    auto g_xxy_yzz_0_yz = buffer_ffsd[112];

    auto g_xxy_yzz_0_zz = buffer_ffsd[113];

    auto g_xxy_zzz_0_xx = buffer_ffsd[114];

    auto g_xxy_zzz_0_xy = buffer_ffsd[115];

    auto g_xxy_zzz_0_xz = buffer_ffsd[116];

    auto g_xxy_zzz_0_yy = buffer_ffsd[117];

    auto g_xxy_zzz_0_yz = buffer_ffsd[118];

    auto g_xxy_zzz_0_zz = buffer_ffsd[119];

    auto g_xxz_xxx_0_xx = buffer_ffsd[120];

    auto g_xxz_xxx_0_xy = buffer_ffsd[121];

    auto g_xxz_xxx_0_xz = buffer_ffsd[122];

    auto g_xxz_xxx_0_yy = buffer_ffsd[123];

    auto g_xxz_xxx_0_yz = buffer_ffsd[124];

    auto g_xxz_xxx_0_zz = buffer_ffsd[125];

    auto g_xxz_xxy_0_xx = buffer_ffsd[126];

    auto g_xxz_xxy_0_xy = buffer_ffsd[127];

    auto g_xxz_xxy_0_xz = buffer_ffsd[128];

    auto g_xxz_xxy_0_yy = buffer_ffsd[129];

    auto g_xxz_xxy_0_yz = buffer_ffsd[130];

    auto g_xxz_xxy_0_zz = buffer_ffsd[131];

    auto g_xxz_xxz_0_xx = buffer_ffsd[132];

    auto g_xxz_xxz_0_xy = buffer_ffsd[133];

    auto g_xxz_xxz_0_xz = buffer_ffsd[134];

    auto g_xxz_xxz_0_yy = buffer_ffsd[135];

    auto g_xxz_xxz_0_yz = buffer_ffsd[136];

    auto g_xxz_xxz_0_zz = buffer_ffsd[137];

    auto g_xxz_xyy_0_xx = buffer_ffsd[138];

    auto g_xxz_xyy_0_xy = buffer_ffsd[139];

    auto g_xxz_xyy_0_xz = buffer_ffsd[140];

    auto g_xxz_xyy_0_yy = buffer_ffsd[141];

    auto g_xxz_xyy_0_yz = buffer_ffsd[142];

    auto g_xxz_xyy_0_zz = buffer_ffsd[143];

    auto g_xxz_xyz_0_xx = buffer_ffsd[144];

    auto g_xxz_xyz_0_xy = buffer_ffsd[145];

    auto g_xxz_xyz_0_xz = buffer_ffsd[146];

    auto g_xxz_xyz_0_yy = buffer_ffsd[147];

    auto g_xxz_xyz_0_yz = buffer_ffsd[148];

    auto g_xxz_xyz_0_zz = buffer_ffsd[149];

    auto g_xxz_xzz_0_xx = buffer_ffsd[150];

    auto g_xxz_xzz_0_xy = buffer_ffsd[151];

    auto g_xxz_xzz_0_xz = buffer_ffsd[152];

    auto g_xxz_xzz_0_yy = buffer_ffsd[153];

    auto g_xxz_xzz_0_yz = buffer_ffsd[154];

    auto g_xxz_xzz_0_zz = buffer_ffsd[155];

    auto g_xxz_yyy_0_xx = buffer_ffsd[156];

    auto g_xxz_yyy_0_xy = buffer_ffsd[157];

    auto g_xxz_yyy_0_xz = buffer_ffsd[158];

    auto g_xxz_yyy_0_yy = buffer_ffsd[159];

    auto g_xxz_yyy_0_yz = buffer_ffsd[160];

    auto g_xxz_yyy_0_zz = buffer_ffsd[161];

    auto g_xxz_yyz_0_xx = buffer_ffsd[162];

    auto g_xxz_yyz_0_xy = buffer_ffsd[163];

    auto g_xxz_yyz_0_xz = buffer_ffsd[164];

    auto g_xxz_yyz_0_yy = buffer_ffsd[165];

    auto g_xxz_yyz_0_yz = buffer_ffsd[166];

    auto g_xxz_yyz_0_zz = buffer_ffsd[167];

    auto g_xxz_yzz_0_xx = buffer_ffsd[168];

    auto g_xxz_yzz_0_xy = buffer_ffsd[169];

    auto g_xxz_yzz_0_xz = buffer_ffsd[170];

    auto g_xxz_yzz_0_yy = buffer_ffsd[171];

    auto g_xxz_yzz_0_yz = buffer_ffsd[172];

    auto g_xxz_yzz_0_zz = buffer_ffsd[173];

    auto g_xxz_zzz_0_xx = buffer_ffsd[174];

    auto g_xxz_zzz_0_xy = buffer_ffsd[175];

    auto g_xxz_zzz_0_xz = buffer_ffsd[176];

    auto g_xxz_zzz_0_yy = buffer_ffsd[177];

    auto g_xxz_zzz_0_yz = buffer_ffsd[178];

    auto g_xxz_zzz_0_zz = buffer_ffsd[179];

    auto g_xyy_xxx_0_xx = buffer_ffsd[180];

    auto g_xyy_xxx_0_xy = buffer_ffsd[181];

    auto g_xyy_xxx_0_xz = buffer_ffsd[182];

    auto g_xyy_xxx_0_yy = buffer_ffsd[183];

    auto g_xyy_xxx_0_yz = buffer_ffsd[184];

    auto g_xyy_xxx_0_zz = buffer_ffsd[185];

    auto g_xyy_xxy_0_xx = buffer_ffsd[186];

    auto g_xyy_xxy_0_xy = buffer_ffsd[187];

    auto g_xyy_xxy_0_xz = buffer_ffsd[188];

    auto g_xyy_xxy_0_yy = buffer_ffsd[189];

    auto g_xyy_xxy_0_yz = buffer_ffsd[190];

    auto g_xyy_xxy_0_zz = buffer_ffsd[191];

    auto g_xyy_xxz_0_xx = buffer_ffsd[192];

    auto g_xyy_xxz_0_xy = buffer_ffsd[193];

    auto g_xyy_xxz_0_xz = buffer_ffsd[194];

    auto g_xyy_xxz_0_yy = buffer_ffsd[195];

    auto g_xyy_xxz_0_yz = buffer_ffsd[196];

    auto g_xyy_xxz_0_zz = buffer_ffsd[197];

    auto g_xyy_xyy_0_xx = buffer_ffsd[198];

    auto g_xyy_xyy_0_xy = buffer_ffsd[199];

    auto g_xyy_xyy_0_xz = buffer_ffsd[200];

    auto g_xyy_xyy_0_yy = buffer_ffsd[201];

    auto g_xyy_xyy_0_yz = buffer_ffsd[202];

    auto g_xyy_xyy_0_zz = buffer_ffsd[203];

    auto g_xyy_xyz_0_xx = buffer_ffsd[204];

    auto g_xyy_xyz_0_xy = buffer_ffsd[205];

    auto g_xyy_xyz_0_xz = buffer_ffsd[206];

    auto g_xyy_xyz_0_yy = buffer_ffsd[207];

    auto g_xyy_xyz_0_yz = buffer_ffsd[208];

    auto g_xyy_xyz_0_zz = buffer_ffsd[209];

    auto g_xyy_xzz_0_xx = buffer_ffsd[210];

    auto g_xyy_xzz_0_xy = buffer_ffsd[211];

    auto g_xyy_xzz_0_xz = buffer_ffsd[212];

    auto g_xyy_xzz_0_yy = buffer_ffsd[213];

    auto g_xyy_xzz_0_yz = buffer_ffsd[214];

    auto g_xyy_xzz_0_zz = buffer_ffsd[215];

    auto g_xyy_yyy_0_xx = buffer_ffsd[216];

    auto g_xyy_yyy_0_xy = buffer_ffsd[217];

    auto g_xyy_yyy_0_xz = buffer_ffsd[218];

    auto g_xyy_yyy_0_yy = buffer_ffsd[219];

    auto g_xyy_yyy_0_yz = buffer_ffsd[220];

    auto g_xyy_yyy_0_zz = buffer_ffsd[221];

    auto g_xyy_yyz_0_xx = buffer_ffsd[222];

    auto g_xyy_yyz_0_xy = buffer_ffsd[223];

    auto g_xyy_yyz_0_xz = buffer_ffsd[224];

    auto g_xyy_yyz_0_yy = buffer_ffsd[225];

    auto g_xyy_yyz_0_yz = buffer_ffsd[226];

    auto g_xyy_yyz_0_zz = buffer_ffsd[227];

    auto g_xyy_yzz_0_xx = buffer_ffsd[228];

    auto g_xyy_yzz_0_xy = buffer_ffsd[229];

    auto g_xyy_yzz_0_xz = buffer_ffsd[230];

    auto g_xyy_yzz_0_yy = buffer_ffsd[231];

    auto g_xyy_yzz_0_yz = buffer_ffsd[232];

    auto g_xyy_yzz_0_zz = buffer_ffsd[233];

    auto g_xyy_zzz_0_xx = buffer_ffsd[234];

    auto g_xyy_zzz_0_xy = buffer_ffsd[235];

    auto g_xyy_zzz_0_xz = buffer_ffsd[236];

    auto g_xyy_zzz_0_yy = buffer_ffsd[237];

    auto g_xyy_zzz_0_yz = buffer_ffsd[238];

    auto g_xyy_zzz_0_zz = buffer_ffsd[239];

    auto g_xyz_xxx_0_xx = buffer_ffsd[240];

    auto g_xyz_xxx_0_xy = buffer_ffsd[241];

    auto g_xyz_xxx_0_xz = buffer_ffsd[242];

    auto g_xyz_xxx_0_yy = buffer_ffsd[243];

    auto g_xyz_xxx_0_yz = buffer_ffsd[244];

    auto g_xyz_xxx_0_zz = buffer_ffsd[245];

    auto g_xyz_xxy_0_xx = buffer_ffsd[246];

    auto g_xyz_xxy_0_xy = buffer_ffsd[247];

    auto g_xyz_xxy_0_xz = buffer_ffsd[248];

    auto g_xyz_xxy_0_yy = buffer_ffsd[249];

    auto g_xyz_xxy_0_yz = buffer_ffsd[250];

    auto g_xyz_xxy_0_zz = buffer_ffsd[251];

    auto g_xyz_xxz_0_xx = buffer_ffsd[252];

    auto g_xyz_xxz_0_xy = buffer_ffsd[253];

    auto g_xyz_xxz_0_xz = buffer_ffsd[254];

    auto g_xyz_xxz_0_yy = buffer_ffsd[255];

    auto g_xyz_xxz_0_yz = buffer_ffsd[256];

    auto g_xyz_xxz_0_zz = buffer_ffsd[257];

    auto g_xyz_xyy_0_xx = buffer_ffsd[258];

    auto g_xyz_xyy_0_xy = buffer_ffsd[259];

    auto g_xyz_xyy_0_xz = buffer_ffsd[260];

    auto g_xyz_xyy_0_yy = buffer_ffsd[261];

    auto g_xyz_xyy_0_yz = buffer_ffsd[262];

    auto g_xyz_xyy_0_zz = buffer_ffsd[263];

    auto g_xyz_xyz_0_xx = buffer_ffsd[264];

    auto g_xyz_xyz_0_xy = buffer_ffsd[265];

    auto g_xyz_xyz_0_xz = buffer_ffsd[266];

    auto g_xyz_xyz_0_yy = buffer_ffsd[267];

    auto g_xyz_xyz_0_yz = buffer_ffsd[268];

    auto g_xyz_xyz_0_zz = buffer_ffsd[269];

    auto g_xyz_xzz_0_xx = buffer_ffsd[270];

    auto g_xyz_xzz_0_xy = buffer_ffsd[271];

    auto g_xyz_xzz_0_xz = buffer_ffsd[272];

    auto g_xyz_xzz_0_yy = buffer_ffsd[273];

    auto g_xyz_xzz_0_yz = buffer_ffsd[274];

    auto g_xyz_xzz_0_zz = buffer_ffsd[275];

    auto g_xyz_yyy_0_xx = buffer_ffsd[276];

    auto g_xyz_yyy_0_xy = buffer_ffsd[277];

    auto g_xyz_yyy_0_xz = buffer_ffsd[278];

    auto g_xyz_yyy_0_yy = buffer_ffsd[279];

    auto g_xyz_yyy_0_yz = buffer_ffsd[280];

    auto g_xyz_yyy_0_zz = buffer_ffsd[281];

    auto g_xyz_yyz_0_xx = buffer_ffsd[282];

    auto g_xyz_yyz_0_xy = buffer_ffsd[283];

    auto g_xyz_yyz_0_xz = buffer_ffsd[284];

    auto g_xyz_yyz_0_yy = buffer_ffsd[285];

    auto g_xyz_yyz_0_yz = buffer_ffsd[286];

    auto g_xyz_yyz_0_zz = buffer_ffsd[287];

    auto g_xyz_yzz_0_xx = buffer_ffsd[288];

    auto g_xyz_yzz_0_xy = buffer_ffsd[289];

    auto g_xyz_yzz_0_xz = buffer_ffsd[290];

    auto g_xyz_yzz_0_yy = buffer_ffsd[291];

    auto g_xyz_yzz_0_yz = buffer_ffsd[292];

    auto g_xyz_yzz_0_zz = buffer_ffsd[293];

    auto g_xyz_zzz_0_xx = buffer_ffsd[294];

    auto g_xyz_zzz_0_xy = buffer_ffsd[295];

    auto g_xyz_zzz_0_xz = buffer_ffsd[296];

    auto g_xyz_zzz_0_yy = buffer_ffsd[297];

    auto g_xyz_zzz_0_yz = buffer_ffsd[298];

    auto g_xyz_zzz_0_zz = buffer_ffsd[299];

    auto g_xzz_xxx_0_xx = buffer_ffsd[300];

    auto g_xzz_xxx_0_xy = buffer_ffsd[301];

    auto g_xzz_xxx_0_xz = buffer_ffsd[302];

    auto g_xzz_xxx_0_yy = buffer_ffsd[303];

    auto g_xzz_xxx_0_yz = buffer_ffsd[304];

    auto g_xzz_xxx_0_zz = buffer_ffsd[305];

    auto g_xzz_xxy_0_xx = buffer_ffsd[306];

    auto g_xzz_xxy_0_xy = buffer_ffsd[307];

    auto g_xzz_xxy_0_xz = buffer_ffsd[308];

    auto g_xzz_xxy_0_yy = buffer_ffsd[309];

    auto g_xzz_xxy_0_yz = buffer_ffsd[310];

    auto g_xzz_xxy_0_zz = buffer_ffsd[311];

    auto g_xzz_xxz_0_xx = buffer_ffsd[312];

    auto g_xzz_xxz_0_xy = buffer_ffsd[313];

    auto g_xzz_xxz_0_xz = buffer_ffsd[314];

    auto g_xzz_xxz_0_yy = buffer_ffsd[315];

    auto g_xzz_xxz_0_yz = buffer_ffsd[316];

    auto g_xzz_xxz_0_zz = buffer_ffsd[317];

    auto g_xzz_xyy_0_xx = buffer_ffsd[318];

    auto g_xzz_xyy_0_xy = buffer_ffsd[319];

    auto g_xzz_xyy_0_xz = buffer_ffsd[320];

    auto g_xzz_xyy_0_yy = buffer_ffsd[321];

    auto g_xzz_xyy_0_yz = buffer_ffsd[322];

    auto g_xzz_xyy_0_zz = buffer_ffsd[323];

    auto g_xzz_xyz_0_xx = buffer_ffsd[324];

    auto g_xzz_xyz_0_xy = buffer_ffsd[325];

    auto g_xzz_xyz_0_xz = buffer_ffsd[326];

    auto g_xzz_xyz_0_yy = buffer_ffsd[327];

    auto g_xzz_xyz_0_yz = buffer_ffsd[328];

    auto g_xzz_xyz_0_zz = buffer_ffsd[329];

    auto g_xzz_xzz_0_xx = buffer_ffsd[330];

    auto g_xzz_xzz_0_xy = buffer_ffsd[331];

    auto g_xzz_xzz_0_xz = buffer_ffsd[332];

    auto g_xzz_xzz_0_yy = buffer_ffsd[333];

    auto g_xzz_xzz_0_yz = buffer_ffsd[334];

    auto g_xzz_xzz_0_zz = buffer_ffsd[335];

    auto g_xzz_yyy_0_xx = buffer_ffsd[336];

    auto g_xzz_yyy_0_xy = buffer_ffsd[337];

    auto g_xzz_yyy_0_xz = buffer_ffsd[338];

    auto g_xzz_yyy_0_yy = buffer_ffsd[339];

    auto g_xzz_yyy_0_yz = buffer_ffsd[340];

    auto g_xzz_yyy_0_zz = buffer_ffsd[341];

    auto g_xzz_yyz_0_xx = buffer_ffsd[342];

    auto g_xzz_yyz_0_xy = buffer_ffsd[343];

    auto g_xzz_yyz_0_xz = buffer_ffsd[344];

    auto g_xzz_yyz_0_yy = buffer_ffsd[345];

    auto g_xzz_yyz_0_yz = buffer_ffsd[346];

    auto g_xzz_yyz_0_zz = buffer_ffsd[347];

    auto g_xzz_yzz_0_xx = buffer_ffsd[348];

    auto g_xzz_yzz_0_xy = buffer_ffsd[349];

    auto g_xzz_yzz_0_xz = buffer_ffsd[350];

    auto g_xzz_yzz_0_yy = buffer_ffsd[351];

    auto g_xzz_yzz_0_yz = buffer_ffsd[352];

    auto g_xzz_yzz_0_zz = buffer_ffsd[353];

    auto g_xzz_zzz_0_xx = buffer_ffsd[354];

    auto g_xzz_zzz_0_xy = buffer_ffsd[355];

    auto g_xzz_zzz_0_xz = buffer_ffsd[356];

    auto g_xzz_zzz_0_yy = buffer_ffsd[357];

    auto g_xzz_zzz_0_yz = buffer_ffsd[358];

    auto g_xzz_zzz_0_zz = buffer_ffsd[359];

    auto g_yyy_xxx_0_xx = buffer_ffsd[360];

    auto g_yyy_xxx_0_xy = buffer_ffsd[361];

    auto g_yyy_xxx_0_xz = buffer_ffsd[362];

    auto g_yyy_xxx_0_yy = buffer_ffsd[363];

    auto g_yyy_xxx_0_yz = buffer_ffsd[364];

    auto g_yyy_xxx_0_zz = buffer_ffsd[365];

    auto g_yyy_xxy_0_xx = buffer_ffsd[366];

    auto g_yyy_xxy_0_xy = buffer_ffsd[367];

    auto g_yyy_xxy_0_xz = buffer_ffsd[368];

    auto g_yyy_xxy_0_yy = buffer_ffsd[369];

    auto g_yyy_xxy_0_yz = buffer_ffsd[370];

    auto g_yyy_xxy_0_zz = buffer_ffsd[371];

    auto g_yyy_xxz_0_xx = buffer_ffsd[372];

    auto g_yyy_xxz_0_xy = buffer_ffsd[373];

    auto g_yyy_xxz_0_xz = buffer_ffsd[374];

    auto g_yyy_xxz_0_yy = buffer_ffsd[375];

    auto g_yyy_xxz_0_yz = buffer_ffsd[376];

    auto g_yyy_xxz_0_zz = buffer_ffsd[377];

    auto g_yyy_xyy_0_xx = buffer_ffsd[378];

    auto g_yyy_xyy_0_xy = buffer_ffsd[379];

    auto g_yyy_xyy_0_xz = buffer_ffsd[380];

    auto g_yyy_xyy_0_yy = buffer_ffsd[381];

    auto g_yyy_xyy_0_yz = buffer_ffsd[382];

    auto g_yyy_xyy_0_zz = buffer_ffsd[383];

    auto g_yyy_xyz_0_xx = buffer_ffsd[384];

    auto g_yyy_xyz_0_xy = buffer_ffsd[385];

    auto g_yyy_xyz_0_xz = buffer_ffsd[386];

    auto g_yyy_xyz_0_yy = buffer_ffsd[387];

    auto g_yyy_xyz_0_yz = buffer_ffsd[388];

    auto g_yyy_xyz_0_zz = buffer_ffsd[389];

    auto g_yyy_xzz_0_xx = buffer_ffsd[390];

    auto g_yyy_xzz_0_xy = buffer_ffsd[391];

    auto g_yyy_xzz_0_xz = buffer_ffsd[392];

    auto g_yyy_xzz_0_yy = buffer_ffsd[393];

    auto g_yyy_xzz_0_yz = buffer_ffsd[394];

    auto g_yyy_xzz_0_zz = buffer_ffsd[395];

    auto g_yyy_yyy_0_xx = buffer_ffsd[396];

    auto g_yyy_yyy_0_xy = buffer_ffsd[397];

    auto g_yyy_yyy_0_xz = buffer_ffsd[398];

    auto g_yyy_yyy_0_yy = buffer_ffsd[399];

    auto g_yyy_yyy_0_yz = buffer_ffsd[400];

    auto g_yyy_yyy_0_zz = buffer_ffsd[401];

    auto g_yyy_yyz_0_xx = buffer_ffsd[402];

    auto g_yyy_yyz_0_xy = buffer_ffsd[403];

    auto g_yyy_yyz_0_xz = buffer_ffsd[404];

    auto g_yyy_yyz_0_yy = buffer_ffsd[405];

    auto g_yyy_yyz_0_yz = buffer_ffsd[406];

    auto g_yyy_yyz_0_zz = buffer_ffsd[407];

    auto g_yyy_yzz_0_xx = buffer_ffsd[408];

    auto g_yyy_yzz_0_xy = buffer_ffsd[409];

    auto g_yyy_yzz_0_xz = buffer_ffsd[410];

    auto g_yyy_yzz_0_yy = buffer_ffsd[411];

    auto g_yyy_yzz_0_yz = buffer_ffsd[412];

    auto g_yyy_yzz_0_zz = buffer_ffsd[413];

    auto g_yyy_zzz_0_xx = buffer_ffsd[414];

    auto g_yyy_zzz_0_xy = buffer_ffsd[415];

    auto g_yyy_zzz_0_xz = buffer_ffsd[416];

    auto g_yyy_zzz_0_yy = buffer_ffsd[417];

    auto g_yyy_zzz_0_yz = buffer_ffsd[418];

    auto g_yyy_zzz_0_zz = buffer_ffsd[419];

    auto g_yyz_xxx_0_xx = buffer_ffsd[420];

    auto g_yyz_xxx_0_xy = buffer_ffsd[421];

    auto g_yyz_xxx_0_xz = buffer_ffsd[422];

    auto g_yyz_xxx_0_yy = buffer_ffsd[423];

    auto g_yyz_xxx_0_yz = buffer_ffsd[424];

    auto g_yyz_xxx_0_zz = buffer_ffsd[425];

    auto g_yyz_xxy_0_xx = buffer_ffsd[426];

    auto g_yyz_xxy_0_xy = buffer_ffsd[427];

    auto g_yyz_xxy_0_xz = buffer_ffsd[428];

    auto g_yyz_xxy_0_yy = buffer_ffsd[429];

    auto g_yyz_xxy_0_yz = buffer_ffsd[430];

    auto g_yyz_xxy_0_zz = buffer_ffsd[431];

    auto g_yyz_xxz_0_xx = buffer_ffsd[432];

    auto g_yyz_xxz_0_xy = buffer_ffsd[433];

    auto g_yyz_xxz_0_xz = buffer_ffsd[434];

    auto g_yyz_xxz_0_yy = buffer_ffsd[435];

    auto g_yyz_xxz_0_yz = buffer_ffsd[436];

    auto g_yyz_xxz_0_zz = buffer_ffsd[437];

    auto g_yyz_xyy_0_xx = buffer_ffsd[438];

    auto g_yyz_xyy_0_xy = buffer_ffsd[439];

    auto g_yyz_xyy_0_xz = buffer_ffsd[440];

    auto g_yyz_xyy_0_yy = buffer_ffsd[441];

    auto g_yyz_xyy_0_yz = buffer_ffsd[442];

    auto g_yyz_xyy_0_zz = buffer_ffsd[443];

    auto g_yyz_xyz_0_xx = buffer_ffsd[444];

    auto g_yyz_xyz_0_xy = buffer_ffsd[445];

    auto g_yyz_xyz_0_xz = buffer_ffsd[446];

    auto g_yyz_xyz_0_yy = buffer_ffsd[447];

    auto g_yyz_xyz_0_yz = buffer_ffsd[448];

    auto g_yyz_xyz_0_zz = buffer_ffsd[449];

    auto g_yyz_xzz_0_xx = buffer_ffsd[450];

    auto g_yyz_xzz_0_xy = buffer_ffsd[451];

    auto g_yyz_xzz_0_xz = buffer_ffsd[452];

    auto g_yyz_xzz_0_yy = buffer_ffsd[453];

    auto g_yyz_xzz_0_yz = buffer_ffsd[454];

    auto g_yyz_xzz_0_zz = buffer_ffsd[455];

    auto g_yyz_yyy_0_xx = buffer_ffsd[456];

    auto g_yyz_yyy_0_xy = buffer_ffsd[457];

    auto g_yyz_yyy_0_xz = buffer_ffsd[458];

    auto g_yyz_yyy_0_yy = buffer_ffsd[459];

    auto g_yyz_yyy_0_yz = buffer_ffsd[460];

    auto g_yyz_yyy_0_zz = buffer_ffsd[461];

    auto g_yyz_yyz_0_xx = buffer_ffsd[462];

    auto g_yyz_yyz_0_xy = buffer_ffsd[463];

    auto g_yyz_yyz_0_xz = buffer_ffsd[464];

    auto g_yyz_yyz_0_yy = buffer_ffsd[465];

    auto g_yyz_yyz_0_yz = buffer_ffsd[466];

    auto g_yyz_yyz_0_zz = buffer_ffsd[467];

    auto g_yyz_yzz_0_xx = buffer_ffsd[468];

    auto g_yyz_yzz_0_xy = buffer_ffsd[469];

    auto g_yyz_yzz_0_xz = buffer_ffsd[470];

    auto g_yyz_yzz_0_yy = buffer_ffsd[471];

    auto g_yyz_yzz_0_yz = buffer_ffsd[472];

    auto g_yyz_yzz_0_zz = buffer_ffsd[473];

    auto g_yyz_zzz_0_xx = buffer_ffsd[474];

    auto g_yyz_zzz_0_xy = buffer_ffsd[475];

    auto g_yyz_zzz_0_xz = buffer_ffsd[476];

    auto g_yyz_zzz_0_yy = buffer_ffsd[477];

    auto g_yyz_zzz_0_yz = buffer_ffsd[478];

    auto g_yyz_zzz_0_zz = buffer_ffsd[479];

    auto g_yzz_xxx_0_xx = buffer_ffsd[480];

    auto g_yzz_xxx_0_xy = buffer_ffsd[481];

    auto g_yzz_xxx_0_xz = buffer_ffsd[482];

    auto g_yzz_xxx_0_yy = buffer_ffsd[483];

    auto g_yzz_xxx_0_yz = buffer_ffsd[484];

    auto g_yzz_xxx_0_zz = buffer_ffsd[485];

    auto g_yzz_xxy_0_xx = buffer_ffsd[486];

    auto g_yzz_xxy_0_xy = buffer_ffsd[487];

    auto g_yzz_xxy_0_xz = buffer_ffsd[488];

    auto g_yzz_xxy_0_yy = buffer_ffsd[489];

    auto g_yzz_xxy_0_yz = buffer_ffsd[490];

    auto g_yzz_xxy_0_zz = buffer_ffsd[491];

    auto g_yzz_xxz_0_xx = buffer_ffsd[492];

    auto g_yzz_xxz_0_xy = buffer_ffsd[493];

    auto g_yzz_xxz_0_xz = buffer_ffsd[494];

    auto g_yzz_xxz_0_yy = buffer_ffsd[495];

    auto g_yzz_xxz_0_yz = buffer_ffsd[496];

    auto g_yzz_xxz_0_zz = buffer_ffsd[497];

    auto g_yzz_xyy_0_xx = buffer_ffsd[498];

    auto g_yzz_xyy_0_xy = buffer_ffsd[499];

    auto g_yzz_xyy_0_xz = buffer_ffsd[500];

    auto g_yzz_xyy_0_yy = buffer_ffsd[501];

    auto g_yzz_xyy_0_yz = buffer_ffsd[502];

    auto g_yzz_xyy_0_zz = buffer_ffsd[503];

    auto g_yzz_xyz_0_xx = buffer_ffsd[504];

    auto g_yzz_xyz_0_xy = buffer_ffsd[505];

    auto g_yzz_xyz_0_xz = buffer_ffsd[506];

    auto g_yzz_xyz_0_yy = buffer_ffsd[507];

    auto g_yzz_xyz_0_yz = buffer_ffsd[508];

    auto g_yzz_xyz_0_zz = buffer_ffsd[509];

    auto g_yzz_xzz_0_xx = buffer_ffsd[510];

    auto g_yzz_xzz_0_xy = buffer_ffsd[511];

    auto g_yzz_xzz_0_xz = buffer_ffsd[512];

    auto g_yzz_xzz_0_yy = buffer_ffsd[513];

    auto g_yzz_xzz_0_yz = buffer_ffsd[514];

    auto g_yzz_xzz_0_zz = buffer_ffsd[515];

    auto g_yzz_yyy_0_xx = buffer_ffsd[516];

    auto g_yzz_yyy_0_xy = buffer_ffsd[517];

    auto g_yzz_yyy_0_xz = buffer_ffsd[518];

    auto g_yzz_yyy_0_yy = buffer_ffsd[519];

    auto g_yzz_yyy_0_yz = buffer_ffsd[520];

    auto g_yzz_yyy_0_zz = buffer_ffsd[521];

    auto g_yzz_yyz_0_xx = buffer_ffsd[522];

    auto g_yzz_yyz_0_xy = buffer_ffsd[523];

    auto g_yzz_yyz_0_xz = buffer_ffsd[524];

    auto g_yzz_yyz_0_yy = buffer_ffsd[525];

    auto g_yzz_yyz_0_yz = buffer_ffsd[526];

    auto g_yzz_yyz_0_zz = buffer_ffsd[527];

    auto g_yzz_yzz_0_xx = buffer_ffsd[528];

    auto g_yzz_yzz_0_xy = buffer_ffsd[529];

    auto g_yzz_yzz_0_xz = buffer_ffsd[530];

    auto g_yzz_yzz_0_yy = buffer_ffsd[531];

    auto g_yzz_yzz_0_yz = buffer_ffsd[532];

    auto g_yzz_yzz_0_zz = buffer_ffsd[533];

    auto g_yzz_zzz_0_xx = buffer_ffsd[534];

    auto g_yzz_zzz_0_xy = buffer_ffsd[535];

    auto g_yzz_zzz_0_xz = buffer_ffsd[536];

    auto g_yzz_zzz_0_yy = buffer_ffsd[537];

    auto g_yzz_zzz_0_yz = buffer_ffsd[538];

    auto g_yzz_zzz_0_zz = buffer_ffsd[539];

    auto g_zzz_xxx_0_xx = buffer_ffsd[540];

    auto g_zzz_xxx_0_xy = buffer_ffsd[541];

    auto g_zzz_xxx_0_xz = buffer_ffsd[542];

    auto g_zzz_xxx_0_yy = buffer_ffsd[543];

    auto g_zzz_xxx_0_yz = buffer_ffsd[544];

    auto g_zzz_xxx_0_zz = buffer_ffsd[545];

    auto g_zzz_xxy_0_xx = buffer_ffsd[546];

    auto g_zzz_xxy_0_xy = buffer_ffsd[547];

    auto g_zzz_xxy_0_xz = buffer_ffsd[548];

    auto g_zzz_xxy_0_yy = buffer_ffsd[549];

    auto g_zzz_xxy_0_yz = buffer_ffsd[550];

    auto g_zzz_xxy_0_zz = buffer_ffsd[551];

    auto g_zzz_xxz_0_xx = buffer_ffsd[552];

    auto g_zzz_xxz_0_xy = buffer_ffsd[553];

    auto g_zzz_xxz_0_xz = buffer_ffsd[554];

    auto g_zzz_xxz_0_yy = buffer_ffsd[555];

    auto g_zzz_xxz_0_yz = buffer_ffsd[556];

    auto g_zzz_xxz_0_zz = buffer_ffsd[557];

    auto g_zzz_xyy_0_xx = buffer_ffsd[558];

    auto g_zzz_xyy_0_xy = buffer_ffsd[559];

    auto g_zzz_xyy_0_xz = buffer_ffsd[560];

    auto g_zzz_xyy_0_yy = buffer_ffsd[561];

    auto g_zzz_xyy_0_yz = buffer_ffsd[562];

    auto g_zzz_xyy_0_zz = buffer_ffsd[563];

    auto g_zzz_xyz_0_xx = buffer_ffsd[564];

    auto g_zzz_xyz_0_xy = buffer_ffsd[565];

    auto g_zzz_xyz_0_xz = buffer_ffsd[566];

    auto g_zzz_xyz_0_yy = buffer_ffsd[567];

    auto g_zzz_xyz_0_yz = buffer_ffsd[568];

    auto g_zzz_xyz_0_zz = buffer_ffsd[569];

    auto g_zzz_xzz_0_xx = buffer_ffsd[570];

    auto g_zzz_xzz_0_xy = buffer_ffsd[571];

    auto g_zzz_xzz_0_xz = buffer_ffsd[572];

    auto g_zzz_xzz_0_yy = buffer_ffsd[573];

    auto g_zzz_xzz_0_yz = buffer_ffsd[574];

    auto g_zzz_xzz_0_zz = buffer_ffsd[575];

    auto g_zzz_yyy_0_xx = buffer_ffsd[576];

    auto g_zzz_yyy_0_xy = buffer_ffsd[577];

    auto g_zzz_yyy_0_xz = buffer_ffsd[578];

    auto g_zzz_yyy_0_yy = buffer_ffsd[579];

    auto g_zzz_yyy_0_yz = buffer_ffsd[580];

    auto g_zzz_yyy_0_zz = buffer_ffsd[581];

    auto g_zzz_yyz_0_xx = buffer_ffsd[582];

    auto g_zzz_yyz_0_xy = buffer_ffsd[583];

    auto g_zzz_yyz_0_xz = buffer_ffsd[584];

    auto g_zzz_yyz_0_yy = buffer_ffsd[585];

    auto g_zzz_yyz_0_yz = buffer_ffsd[586];

    auto g_zzz_yyz_0_zz = buffer_ffsd[587];

    auto g_zzz_yzz_0_xx = buffer_ffsd[588];

    auto g_zzz_yzz_0_xy = buffer_ffsd[589];

    auto g_zzz_yzz_0_xz = buffer_ffsd[590];

    auto g_zzz_yzz_0_yy = buffer_ffsd[591];

    auto g_zzz_yzz_0_yz = buffer_ffsd[592];

    auto g_zzz_yzz_0_zz = buffer_ffsd[593];

    auto g_zzz_zzz_0_xx = buffer_ffsd[594];

    auto g_zzz_zzz_0_xy = buffer_ffsd[595];

    auto g_zzz_zzz_0_xz = buffer_ffsd[596];

    auto g_zzz_zzz_0_yy = buffer_ffsd[597];

    auto g_zzz_zzz_0_yz = buffer_ffsd[598];

    auto g_zzz_zzz_0_zz = buffer_ffsd[599];

    /// Set up components of integrals buffer : buffer_1100_ddsd

    auto g_x_x_0_0_xx_xx_0_xx = buffer_1100_ddsd[0];

    auto g_x_x_0_0_xx_xx_0_xy = buffer_1100_ddsd[1];

    auto g_x_x_0_0_xx_xx_0_xz = buffer_1100_ddsd[2];

    auto g_x_x_0_0_xx_xx_0_yy = buffer_1100_ddsd[3];

    auto g_x_x_0_0_xx_xx_0_yz = buffer_1100_ddsd[4];

    auto g_x_x_0_0_xx_xx_0_zz = buffer_1100_ddsd[5];

    auto g_x_x_0_0_xx_xy_0_xx = buffer_1100_ddsd[6];

    auto g_x_x_0_0_xx_xy_0_xy = buffer_1100_ddsd[7];

    auto g_x_x_0_0_xx_xy_0_xz = buffer_1100_ddsd[8];

    auto g_x_x_0_0_xx_xy_0_yy = buffer_1100_ddsd[9];

    auto g_x_x_0_0_xx_xy_0_yz = buffer_1100_ddsd[10];

    auto g_x_x_0_0_xx_xy_0_zz = buffer_1100_ddsd[11];

    auto g_x_x_0_0_xx_xz_0_xx = buffer_1100_ddsd[12];

    auto g_x_x_0_0_xx_xz_0_xy = buffer_1100_ddsd[13];

    auto g_x_x_0_0_xx_xz_0_xz = buffer_1100_ddsd[14];

    auto g_x_x_0_0_xx_xz_0_yy = buffer_1100_ddsd[15];

    auto g_x_x_0_0_xx_xz_0_yz = buffer_1100_ddsd[16];

    auto g_x_x_0_0_xx_xz_0_zz = buffer_1100_ddsd[17];

    auto g_x_x_0_0_xx_yy_0_xx = buffer_1100_ddsd[18];

    auto g_x_x_0_0_xx_yy_0_xy = buffer_1100_ddsd[19];

    auto g_x_x_0_0_xx_yy_0_xz = buffer_1100_ddsd[20];

    auto g_x_x_0_0_xx_yy_0_yy = buffer_1100_ddsd[21];

    auto g_x_x_0_0_xx_yy_0_yz = buffer_1100_ddsd[22];

    auto g_x_x_0_0_xx_yy_0_zz = buffer_1100_ddsd[23];

    auto g_x_x_0_0_xx_yz_0_xx = buffer_1100_ddsd[24];

    auto g_x_x_0_0_xx_yz_0_xy = buffer_1100_ddsd[25];

    auto g_x_x_0_0_xx_yz_0_xz = buffer_1100_ddsd[26];

    auto g_x_x_0_0_xx_yz_0_yy = buffer_1100_ddsd[27];

    auto g_x_x_0_0_xx_yz_0_yz = buffer_1100_ddsd[28];

    auto g_x_x_0_0_xx_yz_0_zz = buffer_1100_ddsd[29];

    auto g_x_x_0_0_xx_zz_0_xx = buffer_1100_ddsd[30];

    auto g_x_x_0_0_xx_zz_0_xy = buffer_1100_ddsd[31];

    auto g_x_x_0_0_xx_zz_0_xz = buffer_1100_ddsd[32];

    auto g_x_x_0_0_xx_zz_0_yy = buffer_1100_ddsd[33];

    auto g_x_x_0_0_xx_zz_0_yz = buffer_1100_ddsd[34];

    auto g_x_x_0_0_xx_zz_0_zz = buffer_1100_ddsd[35];

    auto g_x_x_0_0_xy_xx_0_xx = buffer_1100_ddsd[36];

    auto g_x_x_0_0_xy_xx_0_xy = buffer_1100_ddsd[37];

    auto g_x_x_0_0_xy_xx_0_xz = buffer_1100_ddsd[38];

    auto g_x_x_0_0_xy_xx_0_yy = buffer_1100_ddsd[39];

    auto g_x_x_0_0_xy_xx_0_yz = buffer_1100_ddsd[40];

    auto g_x_x_0_0_xy_xx_0_zz = buffer_1100_ddsd[41];

    auto g_x_x_0_0_xy_xy_0_xx = buffer_1100_ddsd[42];

    auto g_x_x_0_0_xy_xy_0_xy = buffer_1100_ddsd[43];

    auto g_x_x_0_0_xy_xy_0_xz = buffer_1100_ddsd[44];

    auto g_x_x_0_0_xy_xy_0_yy = buffer_1100_ddsd[45];

    auto g_x_x_0_0_xy_xy_0_yz = buffer_1100_ddsd[46];

    auto g_x_x_0_0_xy_xy_0_zz = buffer_1100_ddsd[47];

    auto g_x_x_0_0_xy_xz_0_xx = buffer_1100_ddsd[48];

    auto g_x_x_0_0_xy_xz_0_xy = buffer_1100_ddsd[49];

    auto g_x_x_0_0_xy_xz_0_xz = buffer_1100_ddsd[50];

    auto g_x_x_0_0_xy_xz_0_yy = buffer_1100_ddsd[51];

    auto g_x_x_0_0_xy_xz_0_yz = buffer_1100_ddsd[52];

    auto g_x_x_0_0_xy_xz_0_zz = buffer_1100_ddsd[53];

    auto g_x_x_0_0_xy_yy_0_xx = buffer_1100_ddsd[54];

    auto g_x_x_0_0_xy_yy_0_xy = buffer_1100_ddsd[55];

    auto g_x_x_0_0_xy_yy_0_xz = buffer_1100_ddsd[56];

    auto g_x_x_0_0_xy_yy_0_yy = buffer_1100_ddsd[57];

    auto g_x_x_0_0_xy_yy_0_yz = buffer_1100_ddsd[58];

    auto g_x_x_0_0_xy_yy_0_zz = buffer_1100_ddsd[59];

    auto g_x_x_0_0_xy_yz_0_xx = buffer_1100_ddsd[60];

    auto g_x_x_0_0_xy_yz_0_xy = buffer_1100_ddsd[61];

    auto g_x_x_0_0_xy_yz_0_xz = buffer_1100_ddsd[62];

    auto g_x_x_0_0_xy_yz_0_yy = buffer_1100_ddsd[63];

    auto g_x_x_0_0_xy_yz_0_yz = buffer_1100_ddsd[64];

    auto g_x_x_0_0_xy_yz_0_zz = buffer_1100_ddsd[65];

    auto g_x_x_0_0_xy_zz_0_xx = buffer_1100_ddsd[66];

    auto g_x_x_0_0_xy_zz_0_xy = buffer_1100_ddsd[67];

    auto g_x_x_0_0_xy_zz_0_xz = buffer_1100_ddsd[68];

    auto g_x_x_0_0_xy_zz_0_yy = buffer_1100_ddsd[69];

    auto g_x_x_0_0_xy_zz_0_yz = buffer_1100_ddsd[70];

    auto g_x_x_0_0_xy_zz_0_zz = buffer_1100_ddsd[71];

    auto g_x_x_0_0_xz_xx_0_xx = buffer_1100_ddsd[72];

    auto g_x_x_0_0_xz_xx_0_xy = buffer_1100_ddsd[73];

    auto g_x_x_0_0_xz_xx_0_xz = buffer_1100_ddsd[74];

    auto g_x_x_0_0_xz_xx_0_yy = buffer_1100_ddsd[75];

    auto g_x_x_0_0_xz_xx_0_yz = buffer_1100_ddsd[76];

    auto g_x_x_0_0_xz_xx_0_zz = buffer_1100_ddsd[77];

    auto g_x_x_0_0_xz_xy_0_xx = buffer_1100_ddsd[78];

    auto g_x_x_0_0_xz_xy_0_xy = buffer_1100_ddsd[79];

    auto g_x_x_0_0_xz_xy_0_xz = buffer_1100_ddsd[80];

    auto g_x_x_0_0_xz_xy_0_yy = buffer_1100_ddsd[81];

    auto g_x_x_0_0_xz_xy_0_yz = buffer_1100_ddsd[82];

    auto g_x_x_0_0_xz_xy_0_zz = buffer_1100_ddsd[83];

    auto g_x_x_0_0_xz_xz_0_xx = buffer_1100_ddsd[84];

    auto g_x_x_0_0_xz_xz_0_xy = buffer_1100_ddsd[85];

    auto g_x_x_0_0_xz_xz_0_xz = buffer_1100_ddsd[86];

    auto g_x_x_0_0_xz_xz_0_yy = buffer_1100_ddsd[87];

    auto g_x_x_0_0_xz_xz_0_yz = buffer_1100_ddsd[88];

    auto g_x_x_0_0_xz_xz_0_zz = buffer_1100_ddsd[89];

    auto g_x_x_0_0_xz_yy_0_xx = buffer_1100_ddsd[90];

    auto g_x_x_0_0_xz_yy_0_xy = buffer_1100_ddsd[91];

    auto g_x_x_0_0_xz_yy_0_xz = buffer_1100_ddsd[92];

    auto g_x_x_0_0_xz_yy_0_yy = buffer_1100_ddsd[93];

    auto g_x_x_0_0_xz_yy_0_yz = buffer_1100_ddsd[94];

    auto g_x_x_0_0_xz_yy_0_zz = buffer_1100_ddsd[95];

    auto g_x_x_0_0_xz_yz_0_xx = buffer_1100_ddsd[96];

    auto g_x_x_0_0_xz_yz_0_xy = buffer_1100_ddsd[97];

    auto g_x_x_0_0_xz_yz_0_xz = buffer_1100_ddsd[98];

    auto g_x_x_0_0_xz_yz_0_yy = buffer_1100_ddsd[99];

    auto g_x_x_0_0_xz_yz_0_yz = buffer_1100_ddsd[100];

    auto g_x_x_0_0_xz_yz_0_zz = buffer_1100_ddsd[101];

    auto g_x_x_0_0_xz_zz_0_xx = buffer_1100_ddsd[102];

    auto g_x_x_0_0_xz_zz_0_xy = buffer_1100_ddsd[103];

    auto g_x_x_0_0_xz_zz_0_xz = buffer_1100_ddsd[104];

    auto g_x_x_0_0_xz_zz_0_yy = buffer_1100_ddsd[105];

    auto g_x_x_0_0_xz_zz_0_yz = buffer_1100_ddsd[106];

    auto g_x_x_0_0_xz_zz_0_zz = buffer_1100_ddsd[107];

    auto g_x_x_0_0_yy_xx_0_xx = buffer_1100_ddsd[108];

    auto g_x_x_0_0_yy_xx_0_xy = buffer_1100_ddsd[109];

    auto g_x_x_0_0_yy_xx_0_xz = buffer_1100_ddsd[110];

    auto g_x_x_0_0_yy_xx_0_yy = buffer_1100_ddsd[111];

    auto g_x_x_0_0_yy_xx_0_yz = buffer_1100_ddsd[112];

    auto g_x_x_0_0_yy_xx_0_zz = buffer_1100_ddsd[113];

    auto g_x_x_0_0_yy_xy_0_xx = buffer_1100_ddsd[114];

    auto g_x_x_0_0_yy_xy_0_xy = buffer_1100_ddsd[115];

    auto g_x_x_0_0_yy_xy_0_xz = buffer_1100_ddsd[116];

    auto g_x_x_0_0_yy_xy_0_yy = buffer_1100_ddsd[117];

    auto g_x_x_0_0_yy_xy_0_yz = buffer_1100_ddsd[118];

    auto g_x_x_0_0_yy_xy_0_zz = buffer_1100_ddsd[119];

    auto g_x_x_0_0_yy_xz_0_xx = buffer_1100_ddsd[120];

    auto g_x_x_0_0_yy_xz_0_xy = buffer_1100_ddsd[121];

    auto g_x_x_0_0_yy_xz_0_xz = buffer_1100_ddsd[122];

    auto g_x_x_0_0_yy_xz_0_yy = buffer_1100_ddsd[123];

    auto g_x_x_0_0_yy_xz_0_yz = buffer_1100_ddsd[124];

    auto g_x_x_0_0_yy_xz_0_zz = buffer_1100_ddsd[125];

    auto g_x_x_0_0_yy_yy_0_xx = buffer_1100_ddsd[126];

    auto g_x_x_0_0_yy_yy_0_xy = buffer_1100_ddsd[127];

    auto g_x_x_0_0_yy_yy_0_xz = buffer_1100_ddsd[128];

    auto g_x_x_0_0_yy_yy_0_yy = buffer_1100_ddsd[129];

    auto g_x_x_0_0_yy_yy_0_yz = buffer_1100_ddsd[130];

    auto g_x_x_0_0_yy_yy_0_zz = buffer_1100_ddsd[131];

    auto g_x_x_0_0_yy_yz_0_xx = buffer_1100_ddsd[132];

    auto g_x_x_0_0_yy_yz_0_xy = buffer_1100_ddsd[133];

    auto g_x_x_0_0_yy_yz_0_xz = buffer_1100_ddsd[134];

    auto g_x_x_0_0_yy_yz_0_yy = buffer_1100_ddsd[135];

    auto g_x_x_0_0_yy_yz_0_yz = buffer_1100_ddsd[136];

    auto g_x_x_0_0_yy_yz_0_zz = buffer_1100_ddsd[137];

    auto g_x_x_0_0_yy_zz_0_xx = buffer_1100_ddsd[138];

    auto g_x_x_0_0_yy_zz_0_xy = buffer_1100_ddsd[139];

    auto g_x_x_0_0_yy_zz_0_xz = buffer_1100_ddsd[140];

    auto g_x_x_0_0_yy_zz_0_yy = buffer_1100_ddsd[141];

    auto g_x_x_0_0_yy_zz_0_yz = buffer_1100_ddsd[142];

    auto g_x_x_0_0_yy_zz_0_zz = buffer_1100_ddsd[143];

    auto g_x_x_0_0_yz_xx_0_xx = buffer_1100_ddsd[144];

    auto g_x_x_0_0_yz_xx_0_xy = buffer_1100_ddsd[145];

    auto g_x_x_0_0_yz_xx_0_xz = buffer_1100_ddsd[146];

    auto g_x_x_0_0_yz_xx_0_yy = buffer_1100_ddsd[147];

    auto g_x_x_0_0_yz_xx_0_yz = buffer_1100_ddsd[148];

    auto g_x_x_0_0_yz_xx_0_zz = buffer_1100_ddsd[149];

    auto g_x_x_0_0_yz_xy_0_xx = buffer_1100_ddsd[150];

    auto g_x_x_0_0_yz_xy_0_xy = buffer_1100_ddsd[151];

    auto g_x_x_0_0_yz_xy_0_xz = buffer_1100_ddsd[152];

    auto g_x_x_0_0_yz_xy_0_yy = buffer_1100_ddsd[153];

    auto g_x_x_0_0_yz_xy_0_yz = buffer_1100_ddsd[154];

    auto g_x_x_0_0_yz_xy_0_zz = buffer_1100_ddsd[155];

    auto g_x_x_0_0_yz_xz_0_xx = buffer_1100_ddsd[156];

    auto g_x_x_0_0_yz_xz_0_xy = buffer_1100_ddsd[157];

    auto g_x_x_0_0_yz_xz_0_xz = buffer_1100_ddsd[158];

    auto g_x_x_0_0_yz_xz_0_yy = buffer_1100_ddsd[159];

    auto g_x_x_0_0_yz_xz_0_yz = buffer_1100_ddsd[160];

    auto g_x_x_0_0_yz_xz_0_zz = buffer_1100_ddsd[161];

    auto g_x_x_0_0_yz_yy_0_xx = buffer_1100_ddsd[162];

    auto g_x_x_0_0_yz_yy_0_xy = buffer_1100_ddsd[163];

    auto g_x_x_0_0_yz_yy_0_xz = buffer_1100_ddsd[164];

    auto g_x_x_0_0_yz_yy_0_yy = buffer_1100_ddsd[165];

    auto g_x_x_0_0_yz_yy_0_yz = buffer_1100_ddsd[166];

    auto g_x_x_0_0_yz_yy_0_zz = buffer_1100_ddsd[167];

    auto g_x_x_0_0_yz_yz_0_xx = buffer_1100_ddsd[168];

    auto g_x_x_0_0_yz_yz_0_xy = buffer_1100_ddsd[169];

    auto g_x_x_0_0_yz_yz_0_xz = buffer_1100_ddsd[170];

    auto g_x_x_0_0_yz_yz_0_yy = buffer_1100_ddsd[171];

    auto g_x_x_0_0_yz_yz_0_yz = buffer_1100_ddsd[172];

    auto g_x_x_0_0_yz_yz_0_zz = buffer_1100_ddsd[173];

    auto g_x_x_0_0_yz_zz_0_xx = buffer_1100_ddsd[174];

    auto g_x_x_0_0_yz_zz_0_xy = buffer_1100_ddsd[175];

    auto g_x_x_0_0_yz_zz_0_xz = buffer_1100_ddsd[176];

    auto g_x_x_0_0_yz_zz_0_yy = buffer_1100_ddsd[177];

    auto g_x_x_0_0_yz_zz_0_yz = buffer_1100_ddsd[178];

    auto g_x_x_0_0_yz_zz_0_zz = buffer_1100_ddsd[179];

    auto g_x_x_0_0_zz_xx_0_xx = buffer_1100_ddsd[180];

    auto g_x_x_0_0_zz_xx_0_xy = buffer_1100_ddsd[181];

    auto g_x_x_0_0_zz_xx_0_xz = buffer_1100_ddsd[182];

    auto g_x_x_0_0_zz_xx_0_yy = buffer_1100_ddsd[183];

    auto g_x_x_0_0_zz_xx_0_yz = buffer_1100_ddsd[184];

    auto g_x_x_0_0_zz_xx_0_zz = buffer_1100_ddsd[185];

    auto g_x_x_0_0_zz_xy_0_xx = buffer_1100_ddsd[186];

    auto g_x_x_0_0_zz_xy_0_xy = buffer_1100_ddsd[187];

    auto g_x_x_0_0_zz_xy_0_xz = buffer_1100_ddsd[188];

    auto g_x_x_0_0_zz_xy_0_yy = buffer_1100_ddsd[189];

    auto g_x_x_0_0_zz_xy_0_yz = buffer_1100_ddsd[190];

    auto g_x_x_0_0_zz_xy_0_zz = buffer_1100_ddsd[191];

    auto g_x_x_0_0_zz_xz_0_xx = buffer_1100_ddsd[192];

    auto g_x_x_0_0_zz_xz_0_xy = buffer_1100_ddsd[193];

    auto g_x_x_0_0_zz_xz_0_xz = buffer_1100_ddsd[194];

    auto g_x_x_0_0_zz_xz_0_yy = buffer_1100_ddsd[195];

    auto g_x_x_0_0_zz_xz_0_yz = buffer_1100_ddsd[196];

    auto g_x_x_0_0_zz_xz_0_zz = buffer_1100_ddsd[197];

    auto g_x_x_0_0_zz_yy_0_xx = buffer_1100_ddsd[198];

    auto g_x_x_0_0_zz_yy_0_xy = buffer_1100_ddsd[199];

    auto g_x_x_0_0_zz_yy_0_xz = buffer_1100_ddsd[200];

    auto g_x_x_0_0_zz_yy_0_yy = buffer_1100_ddsd[201];

    auto g_x_x_0_0_zz_yy_0_yz = buffer_1100_ddsd[202];

    auto g_x_x_0_0_zz_yy_0_zz = buffer_1100_ddsd[203];

    auto g_x_x_0_0_zz_yz_0_xx = buffer_1100_ddsd[204];

    auto g_x_x_0_0_zz_yz_0_xy = buffer_1100_ddsd[205];

    auto g_x_x_0_0_zz_yz_0_xz = buffer_1100_ddsd[206];

    auto g_x_x_0_0_zz_yz_0_yy = buffer_1100_ddsd[207];

    auto g_x_x_0_0_zz_yz_0_yz = buffer_1100_ddsd[208];

    auto g_x_x_0_0_zz_yz_0_zz = buffer_1100_ddsd[209];

    auto g_x_x_0_0_zz_zz_0_xx = buffer_1100_ddsd[210];

    auto g_x_x_0_0_zz_zz_0_xy = buffer_1100_ddsd[211];

    auto g_x_x_0_0_zz_zz_0_xz = buffer_1100_ddsd[212];

    auto g_x_x_0_0_zz_zz_0_yy = buffer_1100_ddsd[213];

    auto g_x_x_0_0_zz_zz_0_yz = buffer_1100_ddsd[214];

    auto g_x_x_0_0_zz_zz_0_zz = buffer_1100_ddsd[215];

    auto g_x_y_0_0_xx_xx_0_xx = buffer_1100_ddsd[216];

    auto g_x_y_0_0_xx_xx_0_xy = buffer_1100_ddsd[217];

    auto g_x_y_0_0_xx_xx_0_xz = buffer_1100_ddsd[218];

    auto g_x_y_0_0_xx_xx_0_yy = buffer_1100_ddsd[219];

    auto g_x_y_0_0_xx_xx_0_yz = buffer_1100_ddsd[220];

    auto g_x_y_0_0_xx_xx_0_zz = buffer_1100_ddsd[221];

    auto g_x_y_0_0_xx_xy_0_xx = buffer_1100_ddsd[222];

    auto g_x_y_0_0_xx_xy_0_xy = buffer_1100_ddsd[223];

    auto g_x_y_0_0_xx_xy_0_xz = buffer_1100_ddsd[224];

    auto g_x_y_0_0_xx_xy_0_yy = buffer_1100_ddsd[225];

    auto g_x_y_0_0_xx_xy_0_yz = buffer_1100_ddsd[226];

    auto g_x_y_0_0_xx_xy_0_zz = buffer_1100_ddsd[227];

    auto g_x_y_0_0_xx_xz_0_xx = buffer_1100_ddsd[228];

    auto g_x_y_0_0_xx_xz_0_xy = buffer_1100_ddsd[229];

    auto g_x_y_0_0_xx_xz_0_xz = buffer_1100_ddsd[230];

    auto g_x_y_0_0_xx_xz_0_yy = buffer_1100_ddsd[231];

    auto g_x_y_0_0_xx_xz_0_yz = buffer_1100_ddsd[232];

    auto g_x_y_0_0_xx_xz_0_zz = buffer_1100_ddsd[233];

    auto g_x_y_0_0_xx_yy_0_xx = buffer_1100_ddsd[234];

    auto g_x_y_0_0_xx_yy_0_xy = buffer_1100_ddsd[235];

    auto g_x_y_0_0_xx_yy_0_xz = buffer_1100_ddsd[236];

    auto g_x_y_0_0_xx_yy_0_yy = buffer_1100_ddsd[237];

    auto g_x_y_0_0_xx_yy_0_yz = buffer_1100_ddsd[238];

    auto g_x_y_0_0_xx_yy_0_zz = buffer_1100_ddsd[239];

    auto g_x_y_0_0_xx_yz_0_xx = buffer_1100_ddsd[240];

    auto g_x_y_0_0_xx_yz_0_xy = buffer_1100_ddsd[241];

    auto g_x_y_0_0_xx_yz_0_xz = buffer_1100_ddsd[242];

    auto g_x_y_0_0_xx_yz_0_yy = buffer_1100_ddsd[243];

    auto g_x_y_0_0_xx_yz_0_yz = buffer_1100_ddsd[244];

    auto g_x_y_0_0_xx_yz_0_zz = buffer_1100_ddsd[245];

    auto g_x_y_0_0_xx_zz_0_xx = buffer_1100_ddsd[246];

    auto g_x_y_0_0_xx_zz_0_xy = buffer_1100_ddsd[247];

    auto g_x_y_0_0_xx_zz_0_xz = buffer_1100_ddsd[248];

    auto g_x_y_0_0_xx_zz_0_yy = buffer_1100_ddsd[249];

    auto g_x_y_0_0_xx_zz_0_yz = buffer_1100_ddsd[250];

    auto g_x_y_0_0_xx_zz_0_zz = buffer_1100_ddsd[251];

    auto g_x_y_0_0_xy_xx_0_xx = buffer_1100_ddsd[252];

    auto g_x_y_0_0_xy_xx_0_xy = buffer_1100_ddsd[253];

    auto g_x_y_0_0_xy_xx_0_xz = buffer_1100_ddsd[254];

    auto g_x_y_0_0_xy_xx_0_yy = buffer_1100_ddsd[255];

    auto g_x_y_0_0_xy_xx_0_yz = buffer_1100_ddsd[256];

    auto g_x_y_0_0_xy_xx_0_zz = buffer_1100_ddsd[257];

    auto g_x_y_0_0_xy_xy_0_xx = buffer_1100_ddsd[258];

    auto g_x_y_0_0_xy_xy_0_xy = buffer_1100_ddsd[259];

    auto g_x_y_0_0_xy_xy_0_xz = buffer_1100_ddsd[260];

    auto g_x_y_0_0_xy_xy_0_yy = buffer_1100_ddsd[261];

    auto g_x_y_0_0_xy_xy_0_yz = buffer_1100_ddsd[262];

    auto g_x_y_0_0_xy_xy_0_zz = buffer_1100_ddsd[263];

    auto g_x_y_0_0_xy_xz_0_xx = buffer_1100_ddsd[264];

    auto g_x_y_0_0_xy_xz_0_xy = buffer_1100_ddsd[265];

    auto g_x_y_0_0_xy_xz_0_xz = buffer_1100_ddsd[266];

    auto g_x_y_0_0_xy_xz_0_yy = buffer_1100_ddsd[267];

    auto g_x_y_0_0_xy_xz_0_yz = buffer_1100_ddsd[268];

    auto g_x_y_0_0_xy_xz_0_zz = buffer_1100_ddsd[269];

    auto g_x_y_0_0_xy_yy_0_xx = buffer_1100_ddsd[270];

    auto g_x_y_0_0_xy_yy_0_xy = buffer_1100_ddsd[271];

    auto g_x_y_0_0_xy_yy_0_xz = buffer_1100_ddsd[272];

    auto g_x_y_0_0_xy_yy_0_yy = buffer_1100_ddsd[273];

    auto g_x_y_0_0_xy_yy_0_yz = buffer_1100_ddsd[274];

    auto g_x_y_0_0_xy_yy_0_zz = buffer_1100_ddsd[275];

    auto g_x_y_0_0_xy_yz_0_xx = buffer_1100_ddsd[276];

    auto g_x_y_0_0_xy_yz_0_xy = buffer_1100_ddsd[277];

    auto g_x_y_0_0_xy_yz_0_xz = buffer_1100_ddsd[278];

    auto g_x_y_0_0_xy_yz_0_yy = buffer_1100_ddsd[279];

    auto g_x_y_0_0_xy_yz_0_yz = buffer_1100_ddsd[280];

    auto g_x_y_0_0_xy_yz_0_zz = buffer_1100_ddsd[281];

    auto g_x_y_0_0_xy_zz_0_xx = buffer_1100_ddsd[282];

    auto g_x_y_0_0_xy_zz_0_xy = buffer_1100_ddsd[283];

    auto g_x_y_0_0_xy_zz_0_xz = buffer_1100_ddsd[284];

    auto g_x_y_0_0_xy_zz_0_yy = buffer_1100_ddsd[285];

    auto g_x_y_0_0_xy_zz_0_yz = buffer_1100_ddsd[286];

    auto g_x_y_0_0_xy_zz_0_zz = buffer_1100_ddsd[287];

    auto g_x_y_0_0_xz_xx_0_xx = buffer_1100_ddsd[288];

    auto g_x_y_0_0_xz_xx_0_xy = buffer_1100_ddsd[289];

    auto g_x_y_0_0_xz_xx_0_xz = buffer_1100_ddsd[290];

    auto g_x_y_0_0_xz_xx_0_yy = buffer_1100_ddsd[291];

    auto g_x_y_0_0_xz_xx_0_yz = buffer_1100_ddsd[292];

    auto g_x_y_0_0_xz_xx_0_zz = buffer_1100_ddsd[293];

    auto g_x_y_0_0_xz_xy_0_xx = buffer_1100_ddsd[294];

    auto g_x_y_0_0_xz_xy_0_xy = buffer_1100_ddsd[295];

    auto g_x_y_0_0_xz_xy_0_xz = buffer_1100_ddsd[296];

    auto g_x_y_0_0_xz_xy_0_yy = buffer_1100_ddsd[297];

    auto g_x_y_0_0_xz_xy_0_yz = buffer_1100_ddsd[298];

    auto g_x_y_0_0_xz_xy_0_zz = buffer_1100_ddsd[299];

    auto g_x_y_0_0_xz_xz_0_xx = buffer_1100_ddsd[300];

    auto g_x_y_0_0_xz_xz_0_xy = buffer_1100_ddsd[301];

    auto g_x_y_0_0_xz_xz_0_xz = buffer_1100_ddsd[302];

    auto g_x_y_0_0_xz_xz_0_yy = buffer_1100_ddsd[303];

    auto g_x_y_0_0_xz_xz_0_yz = buffer_1100_ddsd[304];

    auto g_x_y_0_0_xz_xz_0_zz = buffer_1100_ddsd[305];

    auto g_x_y_0_0_xz_yy_0_xx = buffer_1100_ddsd[306];

    auto g_x_y_0_0_xz_yy_0_xy = buffer_1100_ddsd[307];

    auto g_x_y_0_0_xz_yy_0_xz = buffer_1100_ddsd[308];

    auto g_x_y_0_0_xz_yy_0_yy = buffer_1100_ddsd[309];

    auto g_x_y_0_0_xz_yy_0_yz = buffer_1100_ddsd[310];

    auto g_x_y_0_0_xz_yy_0_zz = buffer_1100_ddsd[311];

    auto g_x_y_0_0_xz_yz_0_xx = buffer_1100_ddsd[312];

    auto g_x_y_0_0_xz_yz_0_xy = buffer_1100_ddsd[313];

    auto g_x_y_0_0_xz_yz_0_xz = buffer_1100_ddsd[314];

    auto g_x_y_0_0_xz_yz_0_yy = buffer_1100_ddsd[315];

    auto g_x_y_0_0_xz_yz_0_yz = buffer_1100_ddsd[316];

    auto g_x_y_0_0_xz_yz_0_zz = buffer_1100_ddsd[317];

    auto g_x_y_0_0_xz_zz_0_xx = buffer_1100_ddsd[318];

    auto g_x_y_0_0_xz_zz_0_xy = buffer_1100_ddsd[319];

    auto g_x_y_0_0_xz_zz_0_xz = buffer_1100_ddsd[320];

    auto g_x_y_0_0_xz_zz_0_yy = buffer_1100_ddsd[321];

    auto g_x_y_0_0_xz_zz_0_yz = buffer_1100_ddsd[322];

    auto g_x_y_0_0_xz_zz_0_zz = buffer_1100_ddsd[323];

    auto g_x_y_0_0_yy_xx_0_xx = buffer_1100_ddsd[324];

    auto g_x_y_0_0_yy_xx_0_xy = buffer_1100_ddsd[325];

    auto g_x_y_0_0_yy_xx_0_xz = buffer_1100_ddsd[326];

    auto g_x_y_0_0_yy_xx_0_yy = buffer_1100_ddsd[327];

    auto g_x_y_0_0_yy_xx_0_yz = buffer_1100_ddsd[328];

    auto g_x_y_0_0_yy_xx_0_zz = buffer_1100_ddsd[329];

    auto g_x_y_0_0_yy_xy_0_xx = buffer_1100_ddsd[330];

    auto g_x_y_0_0_yy_xy_0_xy = buffer_1100_ddsd[331];

    auto g_x_y_0_0_yy_xy_0_xz = buffer_1100_ddsd[332];

    auto g_x_y_0_0_yy_xy_0_yy = buffer_1100_ddsd[333];

    auto g_x_y_0_0_yy_xy_0_yz = buffer_1100_ddsd[334];

    auto g_x_y_0_0_yy_xy_0_zz = buffer_1100_ddsd[335];

    auto g_x_y_0_0_yy_xz_0_xx = buffer_1100_ddsd[336];

    auto g_x_y_0_0_yy_xz_0_xy = buffer_1100_ddsd[337];

    auto g_x_y_0_0_yy_xz_0_xz = buffer_1100_ddsd[338];

    auto g_x_y_0_0_yy_xz_0_yy = buffer_1100_ddsd[339];

    auto g_x_y_0_0_yy_xz_0_yz = buffer_1100_ddsd[340];

    auto g_x_y_0_0_yy_xz_0_zz = buffer_1100_ddsd[341];

    auto g_x_y_0_0_yy_yy_0_xx = buffer_1100_ddsd[342];

    auto g_x_y_0_0_yy_yy_0_xy = buffer_1100_ddsd[343];

    auto g_x_y_0_0_yy_yy_0_xz = buffer_1100_ddsd[344];

    auto g_x_y_0_0_yy_yy_0_yy = buffer_1100_ddsd[345];

    auto g_x_y_0_0_yy_yy_0_yz = buffer_1100_ddsd[346];

    auto g_x_y_0_0_yy_yy_0_zz = buffer_1100_ddsd[347];

    auto g_x_y_0_0_yy_yz_0_xx = buffer_1100_ddsd[348];

    auto g_x_y_0_0_yy_yz_0_xy = buffer_1100_ddsd[349];

    auto g_x_y_0_0_yy_yz_0_xz = buffer_1100_ddsd[350];

    auto g_x_y_0_0_yy_yz_0_yy = buffer_1100_ddsd[351];

    auto g_x_y_0_0_yy_yz_0_yz = buffer_1100_ddsd[352];

    auto g_x_y_0_0_yy_yz_0_zz = buffer_1100_ddsd[353];

    auto g_x_y_0_0_yy_zz_0_xx = buffer_1100_ddsd[354];

    auto g_x_y_0_0_yy_zz_0_xy = buffer_1100_ddsd[355];

    auto g_x_y_0_0_yy_zz_0_xz = buffer_1100_ddsd[356];

    auto g_x_y_0_0_yy_zz_0_yy = buffer_1100_ddsd[357];

    auto g_x_y_0_0_yy_zz_0_yz = buffer_1100_ddsd[358];

    auto g_x_y_0_0_yy_zz_0_zz = buffer_1100_ddsd[359];

    auto g_x_y_0_0_yz_xx_0_xx = buffer_1100_ddsd[360];

    auto g_x_y_0_0_yz_xx_0_xy = buffer_1100_ddsd[361];

    auto g_x_y_0_0_yz_xx_0_xz = buffer_1100_ddsd[362];

    auto g_x_y_0_0_yz_xx_0_yy = buffer_1100_ddsd[363];

    auto g_x_y_0_0_yz_xx_0_yz = buffer_1100_ddsd[364];

    auto g_x_y_0_0_yz_xx_0_zz = buffer_1100_ddsd[365];

    auto g_x_y_0_0_yz_xy_0_xx = buffer_1100_ddsd[366];

    auto g_x_y_0_0_yz_xy_0_xy = buffer_1100_ddsd[367];

    auto g_x_y_0_0_yz_xy_0_xz = buffer_1100_ddsd[368];

    auto g_x_y_0_0_yz_xy_0_yy = buffer_1100_ddsd[369];

    auto g_x_y_0_0_yz_xy_0_yz = buffer_1100_ddsd[370];

    auto g_x_y_0_0_yz_xy_0_zz = buffer_1100_ddsd[371];

    auto g_x_y_0_0_yz_xz_0_xx = buffer_1100_ddsd[372];

    auto g_x_y_0_0_yz_xz_0_xy = buffer_1100_ddsd[373];

    auto g_x_y_0_0_yz_xz_0_xz = buffer_1100_ddsd[374];

    auto g_x_y_0_0_yz_xz_0_yy = buffer_1100_ddsd[375];

    auto g_x_y_0_0_yz_xz_0_yz = buffer_1100_ddsd[376];

    auto g_x_y_0_0_yz_xz_0_zz = buffer_1100_ddsd[377];

    auto g_x_y_0_0_yz_yy_0_xx = buffer_1100_ddsd[378];

    auto g_x_y_0_0_yz_yy_0_xy = buffer_1100_ddsd[379];

    auto g_x_y_0_0_yz_yy_0_xz = buffer_1100_ddsd[380];

    auto g_x_y_0_0_yz_yy_0_yy = buffer_1100_ddsd[381];

    auto g_x_y_0_0_yz_yy_0_yz = buffer_1100_ddsd[382];

    auto g_x_y_0_0_yz_yy_0_zz = buffer_1100_ddsd[383];

    auto g_x_y_0_0_yz_yz_0_xx = buffer_1100_ddsd[384];

    auto g_x_y_0_0_yz_yz_0_xy = buffer_1100_ddsd[385];

    auto g_x_y_0_0_yz_yz_0_xz = buffer_1100_ddsd[386];

    auto g_x_y_0_0_yz_yz_0_yy = buffer_1100_ddsd[387];

    auto g_x_y_0_0_yz_yz_0_yz = buffer_1100_ddsd[388];

    auto g_x_y_0_0_yz_yz_0_zz = buffer_1100_ddsd[389];

    auto g_x_y_0_0_yz_zz_0_xx = buffer_1100_ddsd[390];

    auto g_x_y_0_0_yz_zz_0_xy = buffer_1100_ddsd[391];

    auto g_x_y_0_0_yz_zz_0_xz = buffer_1100_ddsd[392];

    auto g_x_y_0_0_yz_zz_0_yy = buffer_1100_ddsd[393];

    auto g_x_y_0_0_yz_zz_0_yz = buffer_1100_ddsd[394];

    auto g_x_y_0_0_yz_zz_0_zz = buffer_1100_ddsd[395];

    auto g_x_y_0_0_zz_xx_0_xx = buffer_1100_ddsd[396];

    auto g_x_y_0_0_zz_xx_0_xy = buffer_1100_ddsd[397];

    auto g_x_y_0_0_zz_xx_0_xz = buffer_1100_ddsd[398];

    auto g_x_y_0_0_zz_xx_0_yy = buffer_1100_ddsd[399];

    auto g_x_y_0_0_zz_xx_0_yz = buffer_1100_ddsd[400];

    auto g_x_y_0_0_zz_xx_0_zz = buffer_1100_ddsd[401];

    auto g_x_y_0_0_zz_xy_0_xx = buffer_1100_ddsd[402];

    auto g_x_y_0_0_zz_xy_0_xy = buffer_1100_ddsd[403];

    auto g_x_y_0_0_zz_xy_0_xz = buffer_1100_ddsd[404];

    auto g_x_y_0_0_zz_xy_0_yy = buffer_1100_ddsd[405];

    auto g_x_y_0_0_zz_xy_0_yz = buffer_1100_ddsd[406];

    auto g_x_y_0_0_zz_xy_0_zz = buffer_1100_ddsd[407];

    auto g_x_y_0_0_zz_xz_0_xx = buffer_1100_ddsd[408];

    auto g_x_y_0_0_zz_xz_0_xy = buffer_1100_ddsd[409];

    auto g_x_y_0_0_zz_xz_0_xz = buffer_1100_ddsd[410];

    auto g_x_y_0_0_zz_xz_0_yy = buffer_1100_ddsd[411];

    auto g_x_y_0_0_zz_xz_0_yz = buffer_1100_ddsd[412];

    auto g_x_y_0_0_zz_xz_0_zz = buffer_1100_ddsd[413];

    auto g_x_y_0_0_zz_yy_0_xx = buffer_1100_ddsd[414];

    auto g_x_y_0_0_zz_yy_0_xy = buffer_1100_ddsd[415];

    auto g_x_y_0_0_zz_yy_0_xz = buffer_1100_ddsd[416];

    auto g_x_y_0_0_zz_yy_0_yy = buffer_1100_ddsd[417];

    auto g_x_y_0_0_zz_yy_0_yz = buffer_1100_ddsd[418];

    auto g_x_y_0_0_zz_yy_0_zz = buffer_1100_ddsd[419];

    auto g_x_y_0_0_zz_yz_0_xx = buffer_1100_ddsd[420];

    auto g_x_y_0_0_zz_yz_0_xy = buffer_1100_ddsd[421];

    auto g_x_y_0_0_zz_yz_0_xz = buffer_1100_ddsd[422];

    auto g_x_y_0_0_zz_yz_0_yy = buffer_1100_ddsd[423];

    auto g_x_y_0_0_zz_yz_0_yz = buffer_1100_ddsd[424];

    auto g_x_y_0_0_zz_yz_0_zz = buffer_1100_ddsd[425];

    auto g_x_y_0_0_zz_zz_0_xx = buffer_1100_ddsd[426];

    auto g_x_y_0_0_zz_zz_0_xy = buffer_1100_ddsd[427];

    auto g_x_y_0_0_zz_zz_0_xz = buffer_1100_ddsd[428];

    auto g_x_y_0_0_zz_zz_0_yy = buffer_1100_ddsd[429];

    auto g_x_y_0_0_zz_zz_0_yz = buffer_1100_ddsd[430];

    auto g_x_y_0_0_zz_zz_0_zz = buffer_1100_ddsd[431];

    auto g_x_z_0_0_xx_xx_0_xx = buffer_1100_ddsd[432];

    auto g_x_z_0_0_xx_xx_0_xy = buffer_1100_ddsd[433];

    auto g_x_z_0_0_xx_xx_0_xz = buffer_1100_ddsd[434];

    auto g_x_z_0_0_xx_xx_0_yy = buffer_1100_ddsd[435];

    auto g_x_z_0_0_xx_xx_0_yz = buffer_1100_ddsd[436];

    auto g_x_z_0_0_xx_xx_0_zz = buffer_1100_ddsd[437];

    auto g_x_z_0_0_xx_xy_0_xx = buffer_1100_ddsd[438];

    auto g_x_z_0_0_xx_xy_0_xy = buffer_1100_ddsd[439];

    auto g_x_z_0_0_xx_xy_0_xz = buffer_1100_ddsd[440];

    auto g_x_z_0_0_xx_xy_0_yy = buffer_1100_ddsd[441];

    auto g_x_z_0_0_xx_xy_0_yz = buffer_1100_ddsd[442];

    auto g_x_z_0_0_xx_xy_0_zz = buffer_1100_ddsd[443];

    auto g_x_z_0_0_xx_xz_0_xx = buffer_1100_ddsd[444];

    auto g_x_z_0_0_xx_xz_0_xy = buffer_1100_ddsd[445];

    auto g_x_z_0_0_xx_xz_0_xz = buffer_1100_ddsd[446];

    auto g_x_z_0_0_xx_xz_0_yy = buffer_1100_ddsd[447];

    auto g_x_z_0_0_xx_xz_0_yz = buffer_1100_ddsd[448];

    auto g_x_z_0_0_xx_xz_0_zz = buffer_1100_ddsd[449];

    auto g_x_z_0_0_xx_yy_0_xx = buffer_1100_ddsd[450];

    auto g_x_z_0_0_xx_yy_0_xy = buffer_1100_ddsd[451];

    auto g_x_z_0_0_xx_yy_0_xz = buffer_1100_ddsd[452];

    auto g_x_z_0_0_xx_yy_0_yy = buffer_1100_ddsd[453];

    auto g_x_z_0_0_xx_yy_0_yz = buffer_1100_ddsd[454];

    auto g_x_z_0_0_xx_yy_0_zz = buffer_1100_ddsd[455];

    auto g_x_z_0_0_xx_yz_0_xx = buffer_1100_ddsd[456];

    auto g_x_z_0_0_xx_yz_0_xy = buffer_1100_ddsd[457];

    auto g_x_z_0_0_xx_yz_0_xz = buffer_1100_ddsd[458];

    auto g_x_z_0_0_xx_yz_0_yy = buffer_1100_ddsd[459];

    auto g_x_z_0_0_xx_yz_0_yz = buffer_1100_ddsd[460];

    auto g_x_z_0_0_xx_yz_0_zz = buffer_1100_ddsd[461];

    auto g_x_z_0_0_xx_zz_0_xx = buffer_1100_ddsd[462];

    auto g_x_z_0_0_xx_zz_0_xy = buffer_1100_ddsd[463];

    auto g_x_z_0_0_xx_zz_0_xz = buffer_1100_ddsd[464];

    auto g_x_z_0_0_xx_zz_0_yy = buffer_1100_ddsd[465];

    auto g_x_z_0_0_xx_zz_0_yz = buffer_1100_ddsd[466];

    auto g_x_z_0_0_xx_zz_0_zz = buffer_1100_ddsd[467];

    auto g_x_z_0_0_xy_xx_0_xx = buffer_1100_ddsd[468];

    auto g_x_z_0_0_xy_xx_0_xy = buffer_1100_ddsd[469];

    auto g_x_z_0_0_xy_xx_0_xz = buffer_1100_ddsd[470];

    auto g_x_z_0_0_xy_xx_0_yy = buffer_1100_ddsd[471];

    auto g_x_z_0_0_xy_xx_0_yz = buffer_1100_ddsd[472];

    auto g_x_z_0_0_xy_xx_0_zz = buffer_1100_ddsd[473];

    auto g_x_z_0_0_xy_xy_0_xx = buffer_1100_ddsd[474];

    auto g_x_z_0_0_xy_xy_0_xy = buffer_1100_ddsd[475];

    auto g_x_z_0_0_xy_xy_0_xz = buffer_1100_ddsd[476];

    auto g_x_z_0_0_xy_xy_0_yy = buffer_1100_ddsd[477];

    auto g_x_z_0_0_xy_xy_0_yz = buffer_1100_ddsd[478];

    auto g_x_z_0_0_xy_xy_0_zz = buffer_1100_ddsd[479];

    auto g_x_z_0_0_xy_xz_0_xx = buffer_1100_ddsd[480];

    auto g_x_z_0_0_xy_xz_0_xy = buffer_1100_ddsd[481];

    auto g_x_z_0_0_xy_xz_0_xz = buffer_1100_ddsd[482];

    auto g_x_z_0_0_xy_xz_0_yy = buffer_1100_ddsd[483];

    auto g_x_z_0_0_xy_xz_0_yz = buffer_1100_ddsd[484];

    auto g_x_z_0_0_xy_xz_0_zz = buffer_1100_ddsd[485];

    auto g_x_z_0_0_xy_yy_0_xx = buffer_1100_ddsd[486];

    auto g_x_z_0_0_xy_yy_0_xy = buffer_1100_ddsd[487];

    auto g_x_z_0_0_xy_yy_0_xz = buffer_1100_ddsd[488];

    auto g_x_z_0_0_xy_yy_0_yy = buffer_1100_ddsd[489];

    auto g_x_z_0_0_xy_yy_0_yz = buffer_1100_ddsd[490];

    auto g_x_z_0_0_xy_yy_0_zz = buffer_1100_ddsd[491];

    auto g_x_z_0_0_xy_yz_0_xx = buffer_1100_ddsd[492];

    auto g_x_z_0_0_xy_yz_0_xy = buffer_1100_ddsd[493];

    auto g_x_z_0_0_xy_yz_0_xz = buffer_1100_ddsd[494];

    auto g_x_z_0_0_xy_yz_0_yy = buffer_1100_ddsd[495];

    auto g_x_z_0_0_xy_yz_0_yz = buffer_1100_ddsd[496];

    auto g_x_z_0_0_xy_yz_0_zz = buffer_1100_ddsd[497];

    auto g_x_z_0_0_xy_zz_0_xx = buffer_1100_ddsd[498];

    auto g_x_z_0_0_xy_zz_0_xy = buffer_1100_ddsd[499];

    auto g_x_z_0_0_xy_zz_0_xz = buffer_1100_ddsd[500];

    auto g_x_z_0_0_xy_zz_0_yy = buffer_1100_ddsd[501];

    auto g_x_z_0_0_xy_zz_0_yz = buffer_1100_ddsd[502];

    auto g_x_z_0_0_xy_zz_0_zz = buffer_1100_ddsd[503];

    auto g_x_z_0_0_xz_xx_0_xx = buffer_1100_ddsd[504];

    auto g_x_z_0_0_xz_xx_0_xy = buffer_1100_ddsd[505];

    auto g_x_z_0_0_xz_xx_0_xz = buffer_1100_ddsd[506];

    auto g_x_z_0_0_xz_xx_0_yy = buffer_1100_ddsd[507];

    auto g_x_z_0_0_xz_xx_0_yz = buffer_1100_ddsd[508];

    auto g_x_z_0_0_xz_xx_0_zz = buffer_1100_ddsd[509];

    auto g_x_z_0_0_xz_xy_0_xx = buffer_1100_ddsd[510];

    auto g_x_z_0_0_xz_xy_0_xy = buffer_1100_ddsd[511];

    auto g_x_z_0_0_xz_xy_0_xz = buffer_1100_ddsd[512];

    auto g_x_z_0_0_xz_xy_0_yy = buffer_1100_ddsd[513];

    auto g_x_z_0_0_xz_xy_0_yz = buffer_1100_ddsd[514];

    auto g_x_z_0_0_xz_xy_0_zz = buffer_1100_ddsd[515];

    auto g_x_z_0_0_xz_xz_0_xx = buffer_1100_ddsd[516];

    auto g_x_z_0_0_xz_xz_0_xy = buffer_1100_ddsd[517];

    auto g_x_z_0_0_xz_xz_0_xz = buffer_1100_ddsd[518];

    auto g_x_z_0_0_xz_xz_0_yy = buffer_1100_ddsd[519];

    auto g_x_z_0_0_xz_xz_0_yz = buffer_1100_ddsd[520];

    auto g_x_z_0_0_xz_xz_0_zz = buffer_1100_ddsd[521];

    auto g_x_z_0_0_xz_yy_0_xx = buffer_1100_ddsd[522];

    auto g_x_z_0_0_xz_yy_0_xy = buffer_1100_ddsd[523];

    auto g_x_z_0_0_xz_yy_0_xz = buffer_1100_ddsd[524];

    auto g_x_z_0_0_xz_yy_0_yy = buffer_1100_ddsd[525];

    auto g_x_z_0_0_xz_yy_0_yz = buffer_1100_ddsd[526];

    auto g_x_z_0_0_xz_yy_0_zz = buffer_1100_ddsd[527];

    auto g_x_z_0_0_xz_yz_0_xx = buffer_1100_ddsd[528];

    auto g_x_z_0_0_xz_yz_0_xy = buffer_1100_ddsd[529];

    auto g_x_z_0_0_xz_yz_0_xz = buffer_1100_ddsd[530];

    auto g_x_z_0_0_xz_yz_0_yy = buffer_1100_ddsd[531];

    auto g_x_z_0_0_xz_yz_0_yz = buffer_1100_ddsd[532];

    auto g_x_z_0_0_xz_yz_0_zz = buffer_1100_ddsd[533];

    auto g_x_z_0_0_xz_zz_0_xx = buffer_1100_ddsd[534];

    auto g_x_z_0_0_xz_zz_0_xy = buffer_1100_ddsd[535];

    auto g_x_z_0_0_xz_zz_0_xz = buffer_1100_ddsd[536];

    auto g_x_z_0_0_xz_zz_0_yy = buffer_1100_ddsd[537];

    auto g_x_z_0_0_xz_zz_0_yz = buffer_1100_ddsd[538];

    auto g_x_z_0_0_xz_zz_0_zz = buffer_1100_ddsd[539];

    auto g_x_z_0_0_yy_xx_0_xx = buffer_1100_ddsd[540];

    auto g_x_z_0_0_yy_xx_0_xy = buffer_1100_ddsd[541];

    auto g_x_z_0_0_yy_xx_0_xz = buffer_1100_ddsd[542];

    auto g_x_z_0_0_yy_xx_0_yy = buffer_1100_ddsd[543];

    auto g_x_z_0_0_yy_xx_0_yz = buffer_1100_ddsd[544];

    auto g_x_z_0_0_yy_xx_0_zz = buffer_1100_ddsd[545];

    auto g_x_z_0_0_yy_xy_0_xx = buffer_1100_ddsd[546];

    auto g_x_z_0_0_yy_xy_0_xy = buffer_1100_ddsd[547];

    auto g_x_z_0_0_yy_xy_0_xz = buffer_1100_ddsd[548];

    auto g_x_z_0_0_yy_xy_0_yy = buffer_1100_ddsd[549];

    auto g_x_z_0_0_yy_xy_0_yz = buffer_1100_ddsd[550];

    auto g_x_z_0_0_yy_xy_0_zz = buffer_1100_ddsd[551];

    auto g_x_z_0_0_yy_xz_0_xx = buffer_1100_ddsd[552];

    auto g_x_z_0_0_yy_xz_0_xy = buffer_1100_ddsd[553];

    auto g_x_z_0_0_yy_xz_0_xz = buffer_1100_ddsd[554];

    auto g_x_z_0_0_yy_xz_0_yy = buffer_1100_ddsd[555];

    auto g_x_z_0_0_yy_xz_0_yz = buffer_1100_ddsd[556];

    auto g_x_z_0_0_yy_xz_0_zz = buffer_1100_ddsd[557];

    auto g_x_z_0_0_yy_yy_0_xx = buffer_1100_ddsd[558];

    auto g_x_z_0_0_yy_yy_0_xy = buffer_1100_ddsd[559];

    auto g_x_z_0_0_yy_yy_0_xz = buffer_1100_ddsd[560];

    auto g_x_z_0_0_yy_yy_0_yy = buffer_1100_ddsd[561];

    auto g_x_z_0_0_yy_yy_0_yz = buffer_1100_ddsd[562];

    auto g_x_z_0_0_yy_yy_0_zz = buffer_1100_ddsd[563];

    auto g_x_z_0_0_yy_yz_0_xx = buffer_1100_ddsd[564];

    auto g_x_z_0_0_yy_yz_0_xy = buffer_1100_ddsd[565];

    auto g_x_z_0_0_yy_yz_0_xz = buffer_1100_ddsd[566];

    auto g_x_z_0_0_yy_yz_0_yy = buffer_1100_ddsd[567];

    auto g_x_z_0_0_yy_yz_0_yz = buffer_1100_ddsd[568];

    auto g_x_z_0_0_yy_yz_0_zz = buffer_1100_ddsd[569];

    auto g_x_z_0_0_yy_zz_0_xx = buffer_1100_ddsd[570];

    auto g_x_z_0_0_yy_zz_0_xy = buffer_1100_ddsd[571];

    auto g_x_z_0_0_yy_zz_0_xz = buffer_1100_ddsd[572];

    auto g_x_z_0_0_yy_zz_0_yy = buffer_1100_ddsd[573];

    auto g_x_z_0_0_yy_zz_0_yz = buffer_1100_ddsd[574];

    auto g_x_z_0_0_yy_zz_0_zz = buffer_1100_ddsd[575];

    auto g_x_z_0_0_yz_xx_0_xx = buffer_1100_ddsd[576];

    auto g_x_z_0_0_yz_xx_0_xy = buffer_1100_ddsd[577];

    auto g_x_z_0_0_yz_xx_0_xz = buffer_1100_ddsd[578];

    auto g_x_z_0_0_yz_xx_0_yy = buffer_1100_ddsd[579];

    auto g_x_z_0_0_yz_xx_0_yz = buffer_1100_ddsd[580];

    auto g_x_z_0_0_yz_xx_0_zz = buffer_1100_ddsd[581];

    auto g_x_z_0_0_yz_xy_0_xx = buffer_1100_ddsd[582];

    auto g_x_z_0_0_yz_xy_0_xy = buffer_1100_ddsd[583];

    auto g_x_z_0_0_yz_xy_0_xz = buffer_1100_ddsd[584];

    auto g_x_z_0_0_yz_xy_0_yy = buffer_1100_ddsd[585];

    auto g_x_z_0_0_yz_xy_0_yz = buffer_1100_ddsd[586];

    auto g_x_z_0_0_yz_xy_0_zz = buffer_1100_ddsd[587];

    auto g_x_z_0_0_yz_xz_0_xx = buffer_1100_ddsd[588];

    auto g_x_z_0_0_yz_xz_0_xy = buffer_1100_ddsd[589];

    auto g_x_z_0_0_yz_xz_0_xz = buffer_1100_ddsd[590];

    auto g_x_z_0_0_yz_xz_0_yy = buffer_1100_ddsd[591];

    auto g_x_z_0_0_yz_xz_0_yz = buffer_1100_ddsd[592];

    auto g_x_z_0_0_yz_xz_0_zz = buffer_1100_ddsd[593];

    auto g_x_z_0_0_yz_yy_0_xx = buffer_1100_ddsd[594];

    auto g_x_z_0_0_yz_yy_0_xy = buffer_1100_ddsd[595];

    auto g_x_z_0_0_yz_yy_0_xz = buffer_1100_ddsd[596];

    auto g_x_z_0_0_yz_yy_0_yy = buffer_1100_ddsd[597];

    auto g_x_z_0_0_yz_yy_0_yz = buffer_1100_ddsd[598];

    auto g_x_z_0_0_yz_yy_0_zz = buffer_1100_ddsd[599];

    auto g_x_z_0_0_yz_yz_0_xx = buffer_1100_ddsd[600];

    auto g_x_z_0_0_yz_yz_0_xy = buffer_1100_ddsd[601];

    auto g_x_z_0_0_yz_yz_0_xz = buffer_1100_ddsd[602];

    auto g_x_z_0_0_yz_yz_0_yy = buffer_1100_ddsd[603];

    auto g_x_z_0_0_yz_yz_0_yz = buffer_1100_ddsd[604];

    auto g_x_z_0_0_yz_yz_0_zz = buffer_1100_ddsd[605];

    auto g_x_z_0_0_yz_zz_0_xx = buffer_1100_ddsd[606];

    auto g_x_z_0_0_yz_zz_0_xy = buffer_1100_ddsd[607];

    auto g_x_z_0_0_yz_zz_0_xz = buffer_1100_ddsd[608];

    auto g_x_z_0_0_yz_zz_0_yy = buffer_1100_ddsd[609];

    auto g_x_z_0_0_yz_zz_0_yz = buffer_1100_ddsd[610];

    auto g_x_z_0_0_yz_zz_0_zz = buffer_1100_ddsd[611];

    auto g_x_z_0_0_zz_xx_0_xx = buffer_1100_ddsd[612];

    auto g_x_z_0_0_zz_xx_0_xy = buffer_1100_ddsd[613];

    auto g_x_z_0_0_zz_xx_0_xz = buffer_1100_ddsd[614];

    auto g_x_z_0_0_zz_xx_0_yy = buffer_1100_ddsd[615];

    auto g_x_z_0_0_zz_xx_0_yz = buffer_1100_ddsd[616];

    auto g_x_z_0_0_zz_xx_0_zz = buffer_1100_ddsd[617];

    auto g_x_z_0_0_zz_xy_0_xx = buffer_1100_ddsd[618];

    auto g_x_z_0_0_zz_xy_0_xy = buffer_1100_ddsd[619];

    auto g_x_z_0_0_zz_xy_0_xz = buffer_1100_ddsd[620];

    auto g_x_z_0_0_zz_xy_0_yy = buffer_1100_ddsd[621];

    auto g_x_z_0_0_zz_xy_0_yz = buffer_1100_ddsd[622];

    auto g_x_z_0_0_zz_xy_0_zz = buffer_1100_ddsd[623];

    auto g_x_z_0_0_zz_xz_0_xx = buffer_1100_ddsd[624];

    auto g_x_z_0_0_zz_xz_0_xy = buffer_1100_ddsd[625];

    auto g_x_z_0_0_zz_xz_0_xz = buffer_1100_ddsd[626];

    auto g_x_z_0_0_zz_xz_0_yy = buffer_1100_ddsd[627];

    auto g_x_z_0_0_zz_xz_0_yz = buffer_1100_ddsd[628];

    auto g_x_z_0_0_zz_xz_0_zz = buffer_1100_ddsd[629];

    auto g_x_z_0_0_zz_yy_0_xx = buffer_1100_ddsd[630];

    auto g_x_z_0_0_zz_yy_0_xy = buffer_1100_ddsd[631];

    auto g_x_z_0_0_zz_yy_0_xz = buffer_1100_ddsd[632];

    auto g_x_z_0_0_zz_yy_0_yy = buffer_1100_ddsd[633];

    auto g_x_z_0_0_zz_yy_0_yz = buffer_1100_ddsd[634];

    auto g_x_z_0_0_zz_yy_0_zz = buffer_1100_ddsd[635];

    auto g_x_z_0_0_zz_yz_0_xx = buffer_1100_ddsd[636];

    auto g_x_z_0_0_zz_yz_0_xy = buffer_1100_ddsd[637];

    auto g_x_z_0_0_zz_yz_0_xz = buffer_1100_ddsd[638];

    auto g_x_z_0_0_zz_yz_0_yy = buffer_1100_ddsd[639];

    auto g_x_z_0_0_zz_yz_0_yz = buffer_1100_ddsd[640];

    auto g_x_z_0_0_zz_yz_0_zz = buffer_1100_ddsd[641];

    auto g_x_z_0_0_zz_zz_0_xx = buffer_1100_ddsd[642];

    auto g_x_z_0_0_zz_zz_0_xy = buffer_1100_ddsd[643];

    auto g_x_z_0_0_zz_zz_0_xz = buffer_1100_ddsd[644];

    auto g_x_z_0_0_zz_zz_0_yy = buffer_1100_ddsd[645];

    auto g_x_z_0_0_zz_zz_0_yz = buffer_1100_ddsd[646];

    auto g_x_z_0_0_zz_zz_0_zz = buffer_1100_ddsd[647];

    auto g_y_x_0_0_xx_xx_0_xx = buffer_1100_ddsd[648];

    auto g_y_x_0_0_xx_xx_0_xy = buffer_1100_ddsd[649];

    auto g_y_x_0_0_xx_xx_0_xz = buffer_1100_ddsd[650];

    auto g_y_x_0_0_xx_xx_0_yy = buffer_1100_ddsd[651];

    auto g_y_x_0_0_xx_xx_0_yz = buffer_1100_ddsd[652];

    auto g_y_x_0_0_xx_xx_0_zz = buffer_1100_ddsd[653];

    auto g_y_x_0_0_xx_xy_0_xx = buffer_1100_ddsd[654];

    auto g_y_x_0_0_xx_xy_0_xy = buffer_1100_ddsd[655];

    auto g_y_x_0_0_xx_xy_0_xz = buffer_1100_ddsd[656];

    auto g_y_x_0_0_xx_xy_0_yy = buffer_1100_ddsd[657];

    auto g_y_x_0_0_xx_xy_0_yz = buffer_1100_ddsd[658];

    auto g_y_x_0_0_xx_xy_0_zz = buffer_1100_ddsd[659];

    auto g_y_x_0_0_xx_xz_0_xx = buffer_1100_ddsd[660];

    auto g_y_x_0_0_xx_xz_0_xy = buffer_1100_ddsd[661];

    auto g_y_x_0_0_xx_xz_0_xz = buffer_1100_ddsd[662];

    auto g_y_x_0_0_xx_xz_0_yy = buffer_1100_ddsd[663];

    auto g_y_x_0_0_xx_xz_0_yz = buffer_1100_ddsd[664];

    auto g_y_x_0_0_xx_xz_0_zz = buffer_1100_ddsd[665];

    auto g_y_x_0_0_xx_yy_0_xx = buffer_1100_ddsd[666];

    auto g_y_x_0_0_xx_yy_0_xy = buffer_1100_ddsd[667];

    auto g_y_x_0_0_xx_yy_0_xz = buffer_1100_ddsd[668];

    auto g_y_x_0_0_xx_yy_0_yy = buffer_1100_ddsd[669];

    auto g_y_x_0_0_xx_yy_0_yz = buffer_1100_ddsd[670];

    auto g_y_x_0_0_xx_yy_0_zz = buffer_1100_ddsd[671];

    auto g_y_x_0_0_xx_yz_0_xx = buffer_1100_ddsd[672];

    auto g_y_x_0_0_xx_yz_0_xy = buffer_1100_ddsd[673];

    auto g_y_x_0_0_xx_yz_0_xz = buffer_1100_ddsd[674];

    auto g_y_x_0_0_xx_yz_0_yy = buffer_1100_ddsd[675];

    auto g_y_x_0_0_xx_yz_0_yz = buffer_1100_ddsd[676];

    auto g_y_x_0_0_xx_yz_0_zz = buffer_1100_ddsd[677];

    auto g_y_x_0_0_xx_zz_0_xx = buffer_1100_ddsd[678];

    auto g_y_x_0_0_xx_zz_0_xy = buffer_1100_ddsd[679];

    auto g_y_x_0_0_xx_zz_0_xz = buffer_1100_ddsd[680];

    auto g_y_x_0_0_xx_zz_0_yy = buffer_1100_ddsd[681];

    auto g_y_x_0_0_xx_zz_0_yz = buffer_1100_ddsd[682];

    auto g_y_x_0_0_xx_zz_0_zz = buffer_1100_ddsd[683];

    auto g_y_x_0_0_xy_xx_0_xx = buffer_1100_ddsd[684];

    auto g_y_x_0_0_xy_xx_0_xy = buffer_1100_ddsd[685];

    auto g_y_x_0_0_xy_xx_0_xz = buffer_1100_ddsd[686];

    auto g_y_x_0_0_xy_xx_0_yy = buffer_1100_ddsd[687];

    auto g_y_x_0_0_xy_xx_0_yz = buffer_1100_ddsd[688];

    auto g_y_x_0_0_xy_xx_0_zz = buffer_1100_ddsd[689];

    auto g_y_x_0_0_xy_xy_0_xx = buffer_1100_ddsd[690];

    auto g_y_x_0_0_xy_xy_0_xy = buffer_1100_ddsd[691];

    auto g_y_x_0_0_xy_xy_0_xz = buffer_1100_ddsd[692];

    auto g_y_x_0_0_xy_xy_0_yy = buffer_1100_ddsd[693];

    auto g_y_x_0_0_xy_xy_0_yz = buffer_1100_ddsd[694];

    auto g_y_x_0_0_xy_xy_0_zz = buffer_1100_ddsd[695];

    auto g_y_x_0_0_xy_xz_0_xx = buffer_1100_ddsd[696];

    auto g_y_x_0_0_xy_xz_0_xy = buffer_1100_ddsd[697];

    auto g_y_x_0_0_xy_xz_0_xz = buffer_1100_ddsd[698];

    auto g_y_x_0_0_xy_xz_0_yy = buffer_1100_ddsd[699];

    auto g_y_x_0_0_xy_xz_0_yz = buffer_1100_ddsd[700];

    auto g_y_x_0_0_xy_xz_0_zz = buffer_1100_ddsd[701];

    auto g_y_x_0_0_xy_yy_0_xx = buffer_1100_ddsd[702];

    auto g_y_x_0_0_xy_yy_0_xy = buffer_1100_ddsd[703];

    auto g_y_x_0_0_xy_yy_0_xz = buffer_1100_ddsd[704];

    auto g_y_x_0_0_xy_yy_0_yy = buffer_1100_ddsd[705];

    auto g_y_x_0_0_xy_yy_0_yz = buffer_1100_ddsd[706];

    auto g_y_x_0_0_xy_yy_0_zz = buffer_1100_ddsd[707];

    auto g_y_x_0_0_xy_yz_0_xx = buffer_1100_ddsd[708];

    auto g_y_x_0_0_xy_yz_0_xy = buffer_1100_ddsd[709];

    auto g_y_x_0_0_xy_yz_0_xz = buffer_1100_ddsd[710];

    auto g_y_x_0_0_xy_yz_0_yy = buffer_1100_ddsd[711];

    auto g_y_x_0_0_xy_yz_0_yz = buffer_1100_ddsd[712];

    auto g_y_x_0_0_xy_yz_0_zz = buffer_1100_ddsd[713];

    auto g_y_x_0_0_xy_zz_0_xx = buffer_1100_ddsd[714];

    auto g_y_x_0_0_xy_zz_0_xy = buffer_1100_ddsd[715];

    auto g_y_x_0_0_xy_zz_0_xz = buffer_1100_ddsd[716];

    auto g_y_x_0_0_xy_zz_0_yy = buffer_1100_ddsd[717];

    auto g_y_x_0_0_xy_zz_0_yz = buffer_1100_ddsd[718];

    auto g_y_x_0_0_xy_zz_0_zz = buffer_1100_ddsd[719];

    auto g_y_x_0_0_xz_xx_0_xx = buffer_1100_ddsd[720];

    auto g_y_x_0_0_xz_xx_0_xy = buffer_1100_ddsd[721];

    auto g_y_x_0_0_xz_xx_0_xz = buffer_1100_ddsd[722];

    auto g_y_x_0_0_xz_xx_0_yy = buffer_1100_ddsd[723];

    auto g_y_x_0_0_xz_xx_0_yz = buffer_1100_ddsd[724];

    auto g_y_x_0_0_xz_xx_0_zz = buffer_1100_ddsd[725];

    auto g_y_x_0_0_xz_xy_0_xx = buffer_1100_ddsd[726];

    auto g_y_x_0_0_xz_xy_0_xy = buffer_1100_ddsd[727];

    auto g_y_x_0_0_xz_xy_0_xz = buffer_1100_ddsd[728];

    auto g_y_x_0_0_xz_xy_0_yy = buffer_1100_ddsd[729];

    auto g_y_x_0_0_xz_xy_0_yz = buffer_1100_ddsd[730];

    auto g_y_x_0_0_xz_xy_0_zz = buffer_1100_ddsd[731];

    auto g_y_x_0_0_xz_xz_0_xx = buffer_1100_ddsd[732];

    auto g_y_x_0_0_xz_xz_0_xy = buffer_1100_ddsd[733];

    auto g_y_x_0_0_xz_xz_0_xz = buffer_1100_ddsd[734];

    auto g_y_x_0_0_xz_xz_0_yy = buffer_1100_ddsd[735];

    auto g_y_x_0_0_xz_xz_0_yz = buffer_1100_ddsd[736];

    auto g_y_x_0_0_xz_xz_0_zz = buffer_1100_ddsd[737];

    auto g_y_x_0_0_xz_yy_0_xx = buffer_1100_ddsd[738];

    auto g_y_x_0_0_xz_yy_0_xy = buffer_1100_ddsd[739];

    auto g_y_x_0_0_xz_yy_0_xz = buffer_1100_ddsd[740];

    auto g_y_x_0_0_xz_yy_0_yy = buffer_1100_ddsd[741];

    auto g_y_x_0_0_xz_yy_0_yz = buffer_1100_ddsd[742];

    auto g_y_x_0_0_xz_yy_0_zz = buffer_1100_ddsd[743];

    auto g_y_x_0_0_xz_yz_0_xx = buffer_1100_ddsd[744];

    auto g_y_x_0_0_xz_yz_0_xy = buffer_1100_ddsd[745];

    auto g_y_x_0_0_xz_yz_0_xz = buffer_1100_ddsd[746];

    auto g_y_x_0_0_xz_yz_0_yy = buffer_1100_ddsd[747];

    auto g_y_x_0_0_xz_yz_0_yz = buffer_1100_ddsd[748];

    auto g_y_x_0_0_xz_yz_0_zz = buffer_1100_ddsd[749];

    auto g_y_x_0_0_xz_zz_0_xx = buffer_1100_ddsd[750];

    auto g_y_x_0_0_xz_zz_0_xy = buffer_1100_ddsd[751];

    auto g_y_x_0_0_xz_zz_0_xz = buffer_1100_ddsd[752];

    auto g_y_x_0_0_xz_zz_0_yy = buffer_1100_ddsd[753];

    auto g_y_x_0_0_xz_zz_0_yz = buffer_1100_ddsd[754];

    auto g_y_x_0_0_xz_zz_0_zz = buffer_1100_ddsd[755];

    auto g_y_x_0_0_yy_xx_0_xx = buffer_1100_ddsd[756];

    auto g_y_x_0_0_yy_xx_0_xy = buffer_1100_ddsd[757];

    auto g_y_x_0_0_yy_xx_0_xz = buffer_1100_ddsd[758];

    auto g_y_x_0_0_yy_xx_0_yy = buffer_1100_ddsd[759];

    auto g_y_x_0_0_yy_xx_0_yz = buffer_1100_ddsd[760];

    auto g_y_x_0_0_yy_xx_0_zz = buffer_1100_ddsd[761];

    auto g_y_x_0_0_yy_xy_0_xx = buffer_1100_ddsd[762];

    auto g_y_x_0_0_yy_xy_0_xy = buffer_1100_ddsd[763];

    auto g_y_x_0_0_yy_xy_0_xz = buffer_1100_ddsd[764];

    auto g_y_x_0_0_yy_xy_0_yy = buffer_1100_ddsd[765];

    auto g_y_x_0_0_yy_xy_0_yz = buffer_1100_ddsd[766];

    auto g_y_x_0_0_yy_xy_0_zz = buffer_1100_ddsd[767];

    auto g_y_x_0_0_yy_xz_0_xx = buffer_1100_ddsd[768];

    auto g_y_x_0_0_yy_xz_0_xy = buffer_1100_ddsd[769];

    auto g_y_x_0_0_yy_xz_0_xz = buffer_1100_ddsd[770];

    auto g_y_x_0_0_yy_xz_0_yy = buffer_1100_ddsd[771];

    auto g_y_x_0_0_yy_xz_0_yz = buffer_1100_ddsd[772];

    auto g_y_x_0_0_yy_xz_0_zz = buffer_1100_ddsd[773];

    auto g_y_x_0_0_yy_yy_0_xx = buffer_1100_ddsd[774];

    auto g_y_x_0_0_yy_yy_0_xy = buffer_1100_ddsd[775];

    auto g_y_x_0_0_yy_yy_0_xz = buffer_1100_ddsd[776];

    auto g_y_x_0_0_yy_yy_0_yy = buffer_1100_ddsd[777];

    auto g_y_x_0_0_yy_yy_0_yz = buffer_1100_ddsd[778];

    auto g_y_x_0_0_yy_yy_0_zz = buffer_1100_ddsd[779];

    auto g_y_x_0_0_yy_yz_0_xx = buffer_1100_ddsd[780];

    auto g_y_x_0_0_yy_yz_0_xy = buffer_1100_ddsd[781];

    auto g_y_x_0_0_yy_yz_0_xz = buffer_1100_ddsd[782];

    auto g_y_x_0_0_yy_yz_0_yy = buffer_1100_ddsd[783];

    auto g_y_x_0_0_yy_yz_0_yz = buffer_1100_ddsd[784];

    auto g_y_x_0_0_yy_yz_0_zz = buffer_1100_ddsd[785];

    auto g_y_x_0_0_yy_zz_0_xx = buffer_1100_ddsd[786];

    auto g_y_x_0_0_yy_zz_0_xy = buffer_1100_ddsd[787];

    auto g_y_x_0_0_yy_zz_0_xz = buffer_1100_ddsd[788];

    auto g_y_x_0_0_yy_zz_0_yy = buffer_1100_ddsd[789];

    auto g_y_x_0_0_yy_zz_0_yz = buffer_1100_ddsd[790];

    auto g_y_x_0_0_yy_zz_0_zz = buffer_1100_ddsd[791];

    auto g_y_x_0_0_yz_xx_0_xx = buffer_1100_ddsd[792];

    auto g_y_x_0_0_yz_xx_0_xy = buffer_1100_ddsd[793];

    auto g_y_x_0_0_yz_xx_0_xz = buffer_1100_ddsd[794];

    auto g_y_x_0_0_yz_xx_0_yy = buffer_1100_ddsd[795];

    auto g_y_x_0_0_yz_xx_0_yz = buffer_1100_ddsd[796];

    auto g_y_x_0_0_yz_xx_0_zz = buffer_1100_ddsd[797];

    auto g_y_x_0_0_yz_xy_0_xx = buffer_1100_ddsd[798];

    auto g_y_x_0_0_yz_xy_0_xy = buffer_1100_ddsd[799];

    auto g_y_x_0_0_yz_xy_0_xz = buffer_1100_ddsd[800];

    auto g_y_x_0_0_yz_xy_0_yy = buffer_1100_ddsd[801];

    auto g_y_x_0_0_yz_xy_0_yz = buffer_1100_ddsd[802];

    auto g_y_x_0_0_yz_xy_0_zz = buffer_1100_ddsd[803];

    auto g_y_x_0_0_yz_xz_0_xx = buffer_1100_ddsd[804];

    auto g_y_x_0_0_yz_xz_0_xy = buffer_1100_ddsd[805];

    auto g_y_x_0_0_yz_xz_0_xz = buffer_1100_ddsd[806];

    auto g_y_x_0_0_yz_xz_0_yy = buffer_1100_ddsd[807];

    auto g_y_x_0_0_yz_xz_0_yz = buffer_1100_ddsd[808];

    auto g_y_x_0_0_yz_xz_0_zz = buffer_1100_ddsd[809];

    auto g_y_x_0_0_yz_yy_0_xx = buffer_1100_ddsd[810];

    auto g_y_x_0_0_yz_yy_0_xy = buffer_1100_ddsd[811];

    auto g_y_x_0_0_yz_yy_0_xz = buffer_1100_ddsd[812];

    auto g_y_x_0_0_yz_yy_0_yy = buffer_1100_ddsd[813];

    auto g_y_x_0_0_yz_yy_0_yz = buffer_1100_ddsd[814];

    auto g_y_x_0_0_yz_yy_0_zz = buffer_1100_ddsd[815];

    auto g_y_x_0_0_yz_yz_0_xx = buffer_1100_ddsd[816];

    auto g_y_x_0_0_yz_yz_0_xy = buffer_1100_ddsd[817];

    auto g_y_x_0_0_yz_yz_0_xz = buffer_1100_ddsd[818];

    auto g_y_x_0_0_yz_yz_0_yy = buffer_1100_ddsd[819];

    auto g_y_x_0_0_yz_yz_0_yz = buffer_1100_ddsd[820];

    auto g_y_x_0_0_yz_yz_0_zz = buffer_1100_ddsd[821];

    auto g_y_x_0_0_yz_zz_0_xx = buffer_1100_ddsd[822];

    auto g_y_x_0_0_yz_zz_0_xy = buffer_1100_ddsd[823];

    auto g_y_x_0_0_yz_zz_0_xz = buffer_1100_ddsd[824];

    auto g_y_x_0_0_yz_zz_0_yy = buffer_1100_ddsd[825];

    auto g_y_x_0_0_yz_zz_0_yz = buffer_1100_ddsd[826];

    auto g_y_x_0_0_yz_zz_0_zz = buffer_1100_ddsd[827];

    auto g_y_x_0_0_zz_xx_0_xx = buffer_1100_ddsd[828];

    auto g_y_x_0_0_zz_xx_0_xy = buffer_1100_ddsd[829];

    auto g_y_x_0_0_zz_xx_0_xz = buffer_1100_ddsd[830];

    auto g_y_x_0_0_zz_xx_0_yy = buffer_1100_ddsd[831];

    auto g_y_x_0_0_zz_xx_0_yz = buffer_1100_ddsd[832];

    auto g_y_x_0_0_zz_xx_0_zz = buffer_1100_ddsd[833];

    auto g_y_x_0_0_zz_xy_0_xx = buffer_1100_ddsd[834];

    auto g_y_x_0_0_zz_xy_0_xy = buffer_1100_ddsd[835];

    auto g_y_x_0_0_zz_xy_0_xz = buffer_1100_ddsd[836];

    auto g_y_x_0_0_zz_xy_0_yy = buffer_1100_ddsd[837];

    auto g_y_x_0_0_zz_xy_0_yz = buffer_1100_ddsd[838];

    auto g_y_x_0_0_zz_xy_0_zz = buffer_1100_ddsd[839];

    auto g_y_x_0_0_zz_xz_0_xx = buffer_1100_ddsd[840];

    auto g_y_x_0_0_zz_xz_0_xy = buffer_1100_ddsd[841];

    auto g_y_x_0_0_zz_xz_0_xz = buffer_1100_ddsd[842];

    auto g_y_x_0_0_zz_xz_0_yy = buffer_1100_ddsd[843];

    auto g_y_x_0_0_zz_xz_0_yz = buffer_1100_ddsd[844];

    auto g_y_x_0_0_zz_xz_0_zz = buffer_1100_ddsd[845];

    auto g_y_x_0_0_zz_yy_0_xx = buffer_1100_ddsd[846];

    auto g_y_x_0_0_zz_yy_0_xy = buffer_1100_ddsd[847];

    auto g_y_x_0_0_zz_yy_0_xz = buffer_1100_ddsd[848];

    auto g_y_x_0_0_zz_yy_0_yy = buffer_1100_ddsd[849];

    auto g_y_x_0_0_zz_yy_0_yz = buffer_1100_ddsd[850];

    auto g_y_x_0_0_zz_yy_0_zz = buffer_1100_ddsd[851];

    auto g_y_x_0_0_zz_yz_0_xx = buffer_1100_ddsd[852];

    auto g_y_x_0_0_zz_yz_0_xy = buffer_1100_ddsd[853];

    auto g_y_x_0_0_zz_yz_0_xz = buffer_1100_ddsd[854];

    auto g_y_x_0_0_zz_yz_0_yy = buffer_1100_ddsd[855];

    auto g_y_x_0_0_zz_yz_0_yz = buffer_1100_ddsd[856];

    auto g_y_x_0_0_zz_yz_0_zz = buffer_1100_ddsd[857];

    auto g_y_x_0_0_zz_zz_0_xx = buffer_1100_ddsd[858];

    auto g_y_x_0_0_zz_zz_0_xy = buffer_1100_ddsd[859];

    auto g_y_x_0_0_zz_zz_0_xz = buffer_1100_ddsd[860];

    auto g_y_x_0_0_zz_zz_0_yy = buffer_1100_ddsd[861];

    auto g_y_x_0_0_zz_zz_0_yz = buffer_1100_ddsd[862];

    auto g_y_x_0_0_zz_zz_0_zz = buffer_1100_ddsd[863];

    auto g_y_y_0_0_xx_xx_0_xx = buffer_1100_ddsd[864];

    auto g_y_y_0_0_xx_xx_0_xy = buffer_1100_ddsd[865];

    auto g_y_y_0_0_xx_xx_0_xz = buffer_1100_ddsd[866];

    auto g_y_y_0_0_xx_xx_0_yy = buffer_1100_ddsd[867];

    auto g_y_y_0_0_xx_xx_0_yz = buffer_1100_ddsd[868];

    auto g_y_y_0_0_xx_xx_0_zz = buffer_1100_ddsd[869];

    auto g_y_y_0_0_xx_xy_0_xx = buffer_1100_ddsd[870];

    auto g_y_y_0_0_xx_xy_0_xy = buffer_1100_ddsd[871];

    auto g_y_y_0_0_xx_xy_0_xz = buffer_1100_ddsd[872];

    auto g_y_y_0_0_xx_xy_0_yy = buffer_1100_ddsd[873];

    auto g_y_y_0_0_xx_xy_0_yz = buffer_1100_ddsd[874];

    auto g_y_y_0_0_xx_xy_0_zz = buffer_1100_ddsd[875];

    auto g_y_y_0_0_xx_xz_0_xx = buffer_1100_ddsd[876];

    auto g_y_y_0_0_xx_xz_0_xy = buffer_1100_ddsd[877];

    auto g_y_y_0_0_xx_xz_0_xz = buffer_1100_ddsd[878];

    auto g_y_y_0_0_xx_xz_0_yy = buffer_1100_ddsd[879];

    auto g_y_y_0_0_xx_xz_0_yz = buffer_1100_ddsd[880];

    auto g_y_y_0_0_xx_xz_0_zz = buffer_1100_ddsd[881];

    auto g_y_y_0_0_xx_yy_0_xx = buffer_1100_ddsd[882];

    auto g_y_y_0_0_xx_yy_0_xy = buffer_1100_ddsd[883];

    auto g_y_y_0_0_xx_yy_0_xz = buffer_1100_ddsd[884];

    auto g_y_y_0_0_xx_yy_0_yy = buffer_1100_ddsd[885];

    auto g_y_y_0_0_xx_yy_0_yz = buffer_1100_ddsd[886];

    auto g_y_y_0_0_xx_yy_0_zz = buffer_1100_ddsd[887];

    auto g_y_y_0_0_xx_yz_0_xx = buffer_1100_ddsd[888];

    auto g_y_y_0_0_xx_yz_0_xy = buffer_1100_ddsd[889];

    auto g_y_y_0_0_xx_yz_0_xz = buffer_1100_ddsd[890];

    auto g_y_y_0_0_xx_yz_0_yy = buffer_1100_ddsd[891];

    auto g_y_y_0_0_xx_yz_0_yz = buffer_1100_ddsd[892];

    auto g_y_y_0_0_xx_yz_0_zz = buffer_1100_ddsd[893];

    auto g_y_y_0_0_xx_zz_0_xx = buffer_1100_ddsd[894];

    auto g_y_y_0_0_xx_zz_0_xy = buffer_1100_ddsd[895];

    auto g_y_y_0_0_xx_zz_0_xz = buffer_1100_ddsd[896];

    auto g_y_y_0_0_xx_zz_0_yy = buffer_1100_ddsd[897];

    auto g_y_y_0_0_xx_zz_0_yz = buffer_1100_ddsd[898];

    auto g_y_y_0_0_xx_zz_0_zz = buffer_1100_ddsd[899];

    auto g_y_y_0_0_xy_xx_0_xx = buffer_1100_ddsd[900];

    auto g_y_y_0_0_xy_xx_0_xy = buffer_1100_ddsd[901];

    auto g_y_y_0_0_xy_xx_0_xz = buffer_1100_ddsd[902];

    auto g_y_y_0_0_xy_xx_0_yy = buffer_1100_ddsd[903];

    auto g_y_y_0_0_xy_xx_0_yz = buffer_1100_ddsd[904];

    auto g_y_y_0_0_xy_xx_0_zz = buffer_1100_ddsd[905];

    auto g_y_y_0_0_xy_xy_0_xx = buffer_1100_ddsd[906];

    auto g_y_y_0_0_xy_xy_0_xy = buffer_1100_ddsd[907];

    auto g_y_y_0_0_xy_xy_0_xz = buffer_1100_ddsd[908];

    auto g_y_y_0_0_xy_xy_0_yy = buffer_1100_ddsd[909];

    auto g_y_y_0_0_xy_xy_0_yz = buffer_1100_ddsd[910];

    auto g_y_y_0_0_xy_xy_0_zz = buffer_1100_ddsd[911];

    auto g_y_y_0_0_xy_xz_0_xx = buffer_1100_ddsd[912];

    auto g_y_y_0_0_xy_xz_0_xy = buffer_1100_ddsd[913];

    auto g_y_y_0_0_xy_xz_0_xz = buffer_1100_ddsd[914];

    auto g_y_y_0_0_xy_xz_0_yy = buffer_1100_ddsd[915];

    auto g_y_y_0_0_xy_xz_0_yz = buffer_1100_ddsd[916];

    auto g_y_y_0_0_xy_xz_0_zz = buffer_1100_ddsd[917];

    auto g_y_y_0_0_xy_yy_0_xx = buffer_1100_ddsd[918];

    auto g_y_y_0_0_xy_yy_0_xy = buffer_1100_ddsd[919];

    auto g_y_y_0_0_xy_yy_0_xz = buffer_1100_ddsd[920];

    auto g_y_y_0_0_xy_yy_0_yy = buffer_1100_ddsd[921];

    auto g_y_y_0_0_xy_yy_0_yz = buffer_1100_ddsd[922];

    auto g_y_y_0_0_xy_yy_0_zz = buffer_1100_ddsd[923];

    auto g_y_y_0_0_xy_yz_0_xx = buffer_1100_ddsd[924];

    auto g_y_y_0_0_xy_yz_0_xy = buffer_1100_ddsd[925];

    auto g_y_y_0_0_xy_yz_0_xz = buffer_1100_ddsd[926];

    auto g_y_y_0_0_xy_yz_0_yy = buffer_1100_ddsd[927];

    auto g_y_y_0_0_xy_yz_0_yz = buffer_1100_ddsd[928];

    auto g_y_y_0_0_xy_yz_0_zz = buffer_1100_ddsd[929];

    auto g_y_y_0_0_xy_zz_0_xx = buffer_1100_ddsd[930];

    auto g_y_y_0_0_xy_zz_0_xy = buffer_1100_ddsd[931];

    auto g_y_y_0_0_xy_zz_0_xz = buffer_1100_ddsd[932];

    auto g_y_y_0_0_xy_zz_0_yy = buffer_1100_ddsd[933];

    auto g_y_y_0_0_xy_zz_0_yz = buffer_1100_ddsd[934];

    auto g_y_y_0_0_xy_zz_0_zz = buffer_1100_ddsd[935];

    auto g_y_y_0_0_xz_xx_0_xx = buffer_1100_ddsd[936];

    auto g_y_y_0_0_xz_xx_0_xy = buffer_1100_ddsd[937];

    auto g_y_y_0_0_xz_xx_0_xz = buffer_1100_ddsd[938];

    auto g_y_y_0_0_xz_xx_0_yy = buffer_1100_ddsd[939];

    auto g_y_y_0_0_xz_xx_0_yz = buffer_1100_ddsd[940];

    auto g_y_y_0_0_xz_xx_0_zz = buffer_1100_ddsd[941];

    auto g_y_y_0_0_xz_xy_0_xx = buffer_1100_ddsd[942];

    auto g_y_y_0_0_xz_xy_0_xy = buffer_1100_ddsd[943];

    auto g_y_y_0_0_xz_xy_0_xz = buffer_1100_ddsd[944];

    auto g_y_y_0_0_xz_xy_0_yy = buffer_1100_ddsd[945];

    auto g_y_y_0_0_xz_xy_0_yz = buffer_1100_ddsd[946];

    auto g_y_y_0_0_xz_xy_0_zz = buffer_1100_ddsd[947];

    auto g_y_y_0_0_xz_xz_0_xx = buffer_1100_ddsd[948];

    auto g_y_y_0_0_xz_xz_0_xy = buffer_1100_ddsd[949];

    auto g_y_y_0_0_xz_xz_0_xz = buffer_1100_ddsd[950];

    auto g_y_y_0_0_xz_xz_0_yy = buffer_1100_ddsd[951];

    auto g_y_y_0_0_xz_xz_0_yz = buffer_1100_ddsd[952];

    auto g_y_y_0_0_xz_xz_0_zz = buffer_1100_ddsd[953];

    auto g_y_y_0_0_xz_yy_0_xx = buffer_1100_ddsd[954];

    auto g_y_y_0_0_xz_yy_0_xy = buffer_1100_ddsd[955];

    auto g_y_y_0_0_xz_yy_0_xz = buffer_1100_ddsd[956];

    auto g_y_y_0_0_xz_yy_0_yy = buffer_1100_ddsd[957];

    auto g_y_y_0_0_xz_yy_0_yz = buffer_1100_ddsd[958];

    auto g_y_y_0_0_xz_yy_0_zz = buffer_1100_ddsd[959];

    auto g_y_y_0_0_xz_yz_0_xx = buffer_1100_ddsd[960];

    auto g_y_y_0_0_xz_yz_0_xy = buffer_1100_ddsd[961];

    auto g_y_y_0_0_xz_yz_0_xz = buffer_1100_ddsd[962];

    auto g_y_y_0_0_xz_yz_0_yy = buffer_1100_ddsd[963];

    auto g_y_y_0_0_xz_yz_0_yz = buffer_1100_ddsd[964];

    auto g_y_y_0_0_xz_yz_0_zz = buffer_1100_ddsd[965];

    auto g_y_y_0_0_xz_zz_0_xx = buffer_1100_ddsd[966];

    auto g_y_y_0_0_xz_zz_0_xy = buffer_1100_ddsd[967];

    auto g_y_y_0_0_xz_zz_0_xz = buffer_1100_ddsd[968];

    auto g_y_y_0_0_xz_zz_0_yy = buffer_1100_ddsd[969];

    auto g_y_y_0_0_xz_zz_0_yz = buffer_1100_ddsd[970];

    auto g_y_y_0_0_xz_zz_0_zz = buffer_1100_ddsd[971];

    auto g_y_y_0_0_yy_xx_0_xx = buffer_1100_ddsd[972];

    auto g_y_y_0_0_yy_xx_0_xy = buffer_1100_ddsd[973];

    auto g_y_y_0_0_yy_xx_0_xz = buffer_1100_ddsd[974];

    auto g_y_y_0_0_yy_xx_0_yy = buffer_1100_ddsd[975];

    auto g_y_y_0_0_yy_xx_0_yz = buffer_1100_ddsd[976];

    auto g_y_y_0_0_yy_xx_0_zz = buffer_1100_ddsd[977];

    auto g_y_y_0_0_yy_xy_0_xx = buffer_1100_ddsd[978];

    auto g_y_y_0_0_yy_xy_0_xy = buffer_1100_ddsd[979];

    auto g_y_y_0_0_yy_xy_0_xz = buffer_1100_ddsd[980];

    auto g_y_y_0_0_yy_xy_0_yy = buffer_1100_ddsd[981];

    auto g_y_y_0_0_yy_xy_0_yz = buffer_1100_ddsd[982];

    auto g_y_y_0_0_yy_xy_0_zz = buffer_1100_ddsd[983];

    auto g_y_y_0_0_yy_xz_0_xx = buffer_1100_ddsd[984];

    auto g_y_y_0_0_yy_xz_0_xy = buffer_1100_ddsd[985];

    auto g_y_y_0_0_yy_xz_0_xz = buffer_1100_ddsd[986];

    auto g_y_y_0_0_yy_xz_0_yy = buffer_1100_ddsd[987];

    auto g_y_y_0_0_yy_xz_0_yz = buffer_1100_ddsd[988];

    auto g_y_y_0_0_yy_xz_0_zz = buffer_1100_ddsd[989];

    auto g_y_y_0_0_yy_yy_0_xx = buffer_1100_ddsd[990];

    auto g_y_y_0_0_yy_yy_0_xy = buffer_1100_ddsd[991];

    auto g_y_y_0_0_yy_yy_0_xz = buffer_1100_ddsd[992];

    auto g_y_y_0_0_yy_yy_0_yy = buffer_1100_ddsd[993];

    auto g_y_y_0_0_yy_yy_0_yz = buffer_1100_ddsd[994];

    auto g_y_y_0_0_yy_yy_0_zz = buffer_1100_ddsd[995];

    auto g_y_y_0_0_yy_yz_0_xx = buffer_1100_ddsd[996];

    auto g_y_y_0_0_yy_yz_0_xy = buffer_1100_ddsd[997];

    auto g_y_y_0_0_yy_yz_0_xz = buffer_1100_ddsd[998];

    auto g_y_y_0_0_yy_yz_0_yy = buffer_1100_ddsd[999];

    auto g_y_y_0_0_yy_yz_0_yz = buffer_1100_ddsd[1000];

    auto g_y_y_0_0_yy_yz_0_zz = buffer_1100_ddsd[1001];

    auto g_y_y_0_0_yy_zz_0_xx = buffer_1100_ddsd[1002];

    auto g_y_y_0_0_yy_zz_0_xy = buffer_1100_ddsd[1003];

    auto g_y_y_0_0_yy_zz_0_xz = buffer_1100_ddsd[1004];

    auto g_y_y_0_0_yy_zz_0_yy = buffer_1100_ddsd[1005];

    auto g_y_y_0_0_yy_zz_0_yz = buffer_1100_ddsd[1006];

    auto g_y_y_0_0_yy_zz_0_zz = buffer_1100_ddsd[1007];

    auto g_y_y_0_0_yz_xx_0_xx = buffer_1100_ddsd[1008];

    auto g_y_y_0_0_yz_xx_0_xy = buffer_1100_ddsd[1009];

    auto g_y_y_0_0_yz_xx_0_xz = buffer_1100_ddsd[1010];

    auto g_y_y_0_0_yz_xx_0_yy = buffer_1100_ddsd[1011];

    auto g_y_y_0_0_yz_xx_0_yz = buffer_1100_ddsd[1012];

    auto g_y_y_0_0_yz_xx_0_zz = buffer_1100_ddsd[1013];

    auto g_y_y_0_0_yz_xy_0_xx = buffer_1100_ddsd[1014];

    auto g_y_y_0_0_yz_xy_0_xy = buffer_1100_ddsd[1015];

    auto g_y_y_0_0_yz_xy_0_xz = buffer_1100_ddsd[1016];

    auto g_y_y_0_0_yz_xy_0_yy = buffer_1100_ddsd[1017];

    auto g_y_y_0_0_yz_xy_0_yz = buffer_1100_ddsd[1018];

    auto g_y_y_0_0_yz_xy_0_zz = buffer_1100_ddsd[1019];

    auto g_y_y_0_0_yz_xz_0_xx = buffer_1100_ddsd[1020];

    auto g_y_y_0_0_yz_xz_0_xy = buffer_1100_ddsd[1021];

    auto g_y_y_0_0_yz_xz_0_xz = buffer_1100_ddsd[1022];

    auto g_y_y_0_0_yz_xz_0_yy = buffer_1100_ddsd[1023];

    auto g_y_y_0_0_yz_xz_0_yz = buffer_1100_ddsd[1024];

    auto g_y_y_0_0_yz_xz_0_zz = buffer_1100_ddsd[1025];

    auto g_y_y_0_0_yz_yy_0_xx = buffer_1100_ddsd[1026];

    auto g_y_y_0_0_yz_yy_0_xy = buffer_1100_ddsd[1027];

    auto g_y_y_0_0_yz_yy_0_xz = buffer_1100_ddsd[1028];

    auto g_y_y_0_0_yz_yy_0_yy = buffer_1100_ddsd[1029];

    auto g_y_y_0_0_yz_yy_0_yz = buffer_1100_ddsd[1030];

    auto g_y_y_0_0_yz_yy_0_zz = buffer_1100_ddsd[1031];

    auto g_y_y_0_0_yz_yz_0_xx = buffer_1100_ddsd[1032];

    auto g_y_y_0_0_yz_yz_0_xy = buffer_1100_ddsd[1033];

    auto g_y_y_0_0_yz_yz_0_xz = buffer_1100_ddsd[1034];

    auto g_y_y_0_0_yz_yz_0_yy = buffer_1100_ddsd[1035];

    auto g_y_y_0_0_yz_yz_0_yz = buffer_1100_ddsd[1036];

    auto g_y_y_0_0_yz_yz_0_zz = buffer_1100_ddsd[1037];

    auto g_y_y_0_0_yz_zz_0_xx = buffer_1100_ddsd[1038];

    auto g_y_y_0_0_yz_zz_0_xy = buffer_1100_ddsd[1039];

    auto g_y_y_0_0_yz_zz_0_xz = buffer_1100_ddsd[1040];

    auto g_y_y_0_0_yz_zz_0_yy = buffer_1100_ddsd[1041];

    auto g_y_y_0_0_yz_zz_0_yz = buffer_1100_ddsd[1042];

    auto g_y_y_0_0_yz_zz_0_zz = buffer_1100_ddsd[1043];

    auto g_y_y_0_0_zz_xx_0_xx = buffer_1100_ddsd[1044];

    auto g_y_y_0_0_zz_xx_0_xy = buffer_1100_ddsd[1045];

    auto g_y_y_0_0_zz_xx_0_xz = buffer_1100_ddsd[1046];

    auto g_y_y_0_0_zz_xx_0_yy = buffer_1100_ddsd[1047];

    auto g_y_y_0_0_zz_xx_0_yz = buffer_1100_ddsd[1048];

    auto g_y_y_0_0_zz_xx_0_zz = buffer_1100_ddsd[1049];

    auto g_y_y_0_0_zz_xy_0_xx = buffer_1100_ddsd[1050];

    auto g_y_y_0_0_zz_xy_0_xy = buffer_1100_ddsd[1051];

    auto g_y_y_0_0_zz_xy_0_xz = buffer_1100_ddsd[1052];

    auto g_y_y_0_0_zz_xy_0_yy = buffer_1100_ddsd[1053];

    auto g_y_y_0_0_zz_xy_0_yz = buffer_1100_ddsd[1054];

    auto g_y_y_0_0_zz_xy_0_zz = buffer_1100_ddsd[1055];

    auto g_y_y_0_0_zz_xz_0_xx = buffer_1100_ddsd[1056];

    auto g_y_y_0_0_zz_xz_0_xy = buffer_1100_ddsd[1057];

    auto g_y_y_0_0_zz_xz_0_xz = buffer_1100_ddsd[1058];

    auto g_y_y_0_0_zz_xz_0_yy = buffer_1100_ddsd[1059];

    auto g_y_y_0_0_zz_xz_0_yz = buffer_1100_ddsd[1060];

    auto g_y_y_0_0_zz_xz_0_zz = buffer_1100_ddsd[1061];

    auto g_y_y_0_0_zz_yy_0_xx = buffer_1100_ddsd[1062];

    auto g_y_y_0_0_zz_yy_0_xy = buffer_1100_ddsd[1063];

    auto g_y_y_0_0_zz_yy_0_xz = buffer_1100_ddsd[1064];

    auto g_y_y_0_0_zz_yy_0_yy = buffer_1100_ddsd[1065];

    auto g_y_y_0_0_zz_yy_0_yz = buffer_1100_ddsd[1066];

    auto g_y_y_0_0_zz_yy_0_zz = buffer_1100_ddsd[1067];

    auto g_y_y_0_0_zz_yz_0_xx = buffer_1100_ddsd[1068];

    auto g_y_y_0_0_zz_yz_0_xy = buffer_1100_ddsd[1069];

    auto g_y_y_0_0_zz_yz_0_xz = buffer_1100_ddsd[1070];

    auto g_y_y_0_0_zz_yz_0_yy = buffer_1100_ddsd[1071];

    auto g_y_y_0_0_zz_yz_0_yz = buffer_1100_ddsd[1072];

    auto g_y_y_0_0_zz_yz_0_zz = buffer_1100_ddsd[1073];

    auto g_y_y_0_0_zz_zz_0_xx = buffer_1100_ddsd[1074];

    auto g_y_y_0_0_zz_zz_0_xy = buffer_1100_ddsd[1075];

    auto g_y_y_0_0_zz_zz_0_xz = buffer_1100_ddsd[1076];

    auto g_y_y_0_0_zz_zz_0_yy = buffer_1100_ddsd[1077];

    auto g_y_y_0_0_zz_zz_0_yz = buffer_1100_ddsd[1078];

    auto g_y_y_0_0_zz_zz_0_zz = buffer_1100_ddsd[1079];

    auto g_y_z_0_0_xx_xx_0_xx = buffer_1100_ddsd[1080];

    auto g_y_z_0_0_xx_xx_0_xy = buffer_1100_ddsd[1081];

    auto g_y_z_0_0_xx_xx_0_xz = buffer_1100_ddsd[1082];

    auto g_y_z_0_0_xx_xx_0_yy = buffer_1100_ddsd[1083];

    auto g_y_z_0_0_xx_xx_0_yz = buffer_1100_ddsd[1084];

    auto g_y_z_0_0_xx_xx_0_zz = buffer_1100_ddsd[1085];

    auto g_y_z_0_0_xx_xy_0_xx = buffer_1100_ddsd[1086];

    auto g_y_z_0_0_xx_xy_0_xy = buffer_1100_ddsd[1087];

    auto g_y_z_0_0_xx_xy_0_xz = buffer_1100_ddsd[1088];

    auto g_y_z_0_0_xx_xy_0_yy = buffer_1100_ddsd[1089];

    auto g_y_z_0_0_xx_xy_0_yz = buffer_1100_ddsd[1090];

    auto g_y_z_0_0_xx_xy_0_zz = buffer_1100_ddsd[1091];

    auto g_y_z_0_0_xx_xz_0_xx = buffer_1100_ddsd[1092];

    auto g_y_z_0_0_xx_xz_0_xy = buffer_1100_ddsd[1093];

    auto g_y_z_0_0_xx_xz_0_xz = buffer_1100_ddsd[1094];

    auto g_y_z_0_0_xx_xz_0_yy = buffer_1100_ddsd[1095];

    auto g_y_z_0_0_xx_xz_0_yz = buffer_1100_ddsd[1096];

    auto g_y_z_0_0_xx_xz_0_zz = buffer_1100_ddsd[1097];

    auto g_y_z_0_0_xx_yy_0_xx = buffer_1100_ddsd[1098];

    auto g_y_z_0_0_xx_yy_0_xy = buffer_1100_ddsd[1099];

    auto g_y_z_0_0_xx_yy_0_xz = buffer_1100_ddsd[1100];

    auto g_y_z_0_0_xx_yy_0_yy = buffer_1100_ddsd[1101];

    auto g_y_z_0_0_xx_yy_0_yz = buffer_1100_ddsd[1102];

    auto g_y_z_0_0_xx_yy_0_zz = buffer_1100_ddsd[1103];

    auto g_y_z_0_0_xx_yz_0_xx = buffer_1100_ddsd[1104];

    auto g_y_z_0_0_xx_yz_0_xy = buffer_1100_ddsd[1105];

    auto g_y_z_0_0_xx_yz_0_xz = buffer_1100_ddsd[1106];

    auto g_y_z_0_0_xx_yz_0_yy = buffer_1100_ddsd[1107];

    auto g_y_z_0_0_xx_yz_0_yz = buffer_1100_ddsd[1108];

    auto g_y_z_0_0_xx_yz_0_zz = buffer_1100_ddsd[1109];

    auto g_y_z_0_0_xx_zz_0_xx = buffer_1100_ddsd[1110];

    auto g_y_z_0_0_xx_zz_0_xy = buffer_1100_ddsd[1111];

    auto g_y_z_0_0_xx_zz_0_xz = buffer_1100_ddsd[1112];

    auto g_y_z_0_0_xx_zz_0_yy = buffer_1100_ddsd[1113];

    auto g_y_z_0_0_xx_zz_0_yz = buffer_1100_ddsd[1114];

    auto g_y_z_0_0_xx_zz_0_zz = buffer_1100_ddsd[1115];

    auto g_y_z_0_0_xy_xx_0_xx = buffer_1100_ddsd[1116];

    auto g_y_z_0_0_xy_xx_0_xy = buffer_1100_ddsd[1117];

    auto g_y_z_0_0_xy_xx_0_xz = buffer_1100_ddsd[1118];

    auto g_y_z_0_0_xy_xx_0_yy = buffer_1100_ddsd[1119];

    auto g_y_z_0_0_xy_xx_0_yz = buffer_1100_ddsd[1120];

    auto g_y_z_0_0_xy_xx_0_zz = buffer_1100_ddsd[1121];

    auto g_y_z_0_0_xy_xy_0_xx = buffer_1100_ddsd[1122];

    auto g_y_z_0_0_xy_xy_0_xy = buffer_1100_ddsd[1123];

    auto g_y_z_0_0_xy_xy_0_xz = buffer_1100_ddsd[1124];

    auto g_y_z_0_0_xy_xy_0_yy = buffer_1100_ddsd[1125];

    auto g_y_z_0_0_xy_xy_0_yz = buffer_1100_ddsd[1126];

    auto g_y_z_0_0_xy_xy_0_zz = buffer_1100_ddsd[1127];

    auto g_y_z_0_0_xy_xz_0_xx = buffer_1100_ddsd[1128];

    auto g_y_z_0_0_xy_xz_0_xy = buffer_1100_ddsd[1129];

    auto g_y_z_0_0_xy_xz_0_xz = buffer_1100_ddsd[1130];

    auto g_y_z_0_0_xy_xz_0_yy = buffer_1100_ddsd[1131];

    auto g_y_z_0_0_xy_xz_0_yz = buffer_1100_ddsd[1132];

    auto g_y_z_0_0_xy_xz_0_zz = buffer_1100_ddsd[1133];

    auto g_y_z_0_0_xy_yy_0_xx = buffer_1100_ddsd[1134];

    auto g_y_z_0_0_xy_yy_0_xy = buffer_1100_ddsd[1135];

    auto g_y_z_0_0_xy_yy_0_xz = buffer_1100_ddsd[1136];

    auto g_y_z_0_0_xy_yy_0_yy = buffer_1100_ddsd[1137];

    auto g_y_z_0_0_xy_yy_0_yz = buffer_1100_ddsd[1138];

    auto g_y_z_0_0_xy_yy_0_zz = buffer_1100_ddsd[1139];

    auto g_y_z_0_0_xy_yz_0_xx = buffer_1100_ddsd[1140];

    auto g_y_z_0_0_xy_yz_0_xy = buffer_1100_ddsd[1141];

    auto g_y_z_0_0_xy_yz_0_xz = buffer_1100_ddsd[1142];

    auto g_y_z_0_0_xy_yz_0_yy = buffer_1100_ddsd[1143];

    auto g_y_z_0_0_xy_yz_0_yz = buffer_1100_ddsd[1144];

    auto g_y_z_0_0_xy_yz_0_zz = buffer_1100_ddsd[1145];

    auto g_y_z_0_0_xy_zz_0_xx = buffer_1100_ddsd[1146];

    auto g_y_z_0_0_xy_zz_0_xy = buffer_1100_ddsd[1147];

    auto g_y_z_0_0_xy_zz_0_xz = buffer_1100_ddsd[1148];

    auto g_y_z_0_0_xy_zz_0_yy = buffer_1100_ddsd[1149];

    auto g_y_z_0_0_xy_zz_0_yz = buffer_1100_ddsd[1150];

    auto g_y_z_0_0_xy_zz_0_zz = buffer_1100_ddsd[1151];

    auto g_y_z_0_0_xz_xx_0_xx = buffer_1100_ddsd[1152];

    auto g_y_z_0_0_xz_xx_0_xy = buffer_1100_ddsd[1153];

    auto g_y_z_0_0_xz_xx_0_xz = buffer_1100_ddsd[1154];

    auto g_y_z_0_0_xz_xx_0_yy = buffer_1100_ddsd[1155];

    auto g_y_z_0_0_xz_xx_0_yz = buffer_1100_ddsd[1156];

    auto g_y_z_0_0_xz_xx_0_zz = buffer_1100_ddsd[1157];

    auto g_y_z_0_0_xz_xy_0_xx = buffer_1100_ddsd[1158];

    auto g_y_z_0_0_xz_xy_0_xy = buffer_1100_ddsd[1159];

    auto g_y_z_0_0_xz_xy_0_xz = buffer_1100_ddsd[1160];

    auto g_y_z_0_0_xz_xy_0_yy = buffer_1100_ddsd[1161];

    auto g_y_z_0_0_xz_xy_0_yz = buffer_1100_ddsd[1162];

    auto g_y_z_0_0_xz_xy_0_zz = buffer_1100_ddsd[1163];

    auto g_y_z_0_0_xz_xz_0_xx = buffer_1100_ddsd[1164];

    auto g_y_z_0_0_xz_xz_0_xy = buffer_1100_ddsd[1165];

    auto g_y_z_0_0_xz_xz_0_xz = buffer_1100_ddsd[1166];

    auto g_y_z_0_0_xz_xz_0_yy = buffer_1100_ddsd[1167];

    auto g_y_z_0_0_xz_xz_0_yz = buffer_1100_ddsd[1168];

    auto g_y_z_0_0_xz_xz_0_zz = buffer_1100_ddsd[1169];

    auto g_y_z_0_0_xz_yy_0_xx = buffer_1100_ddsd[1170];

    auto g_y_z_0_0_xz_yy_0_xy = buffer_1100_ddsd[1171];

    auto g_y_z_0_0_xz_yy_0_xz = buffer_1100_ddsd[1172];

    auto g_y_z_0_0_xz_yy_0_yy = buffer_1100_ddsd[1173];

    auto g_y_z_0_0_xz_yy_0_yz = buffer_1100_ddsd[1174];

    auto g_y_z_0_0_xz_yy_0_zz = buffer_1100_ddsd[1175];

    auto g_y_z_0_0_xz_yz_0_xx = buffer_1100_ddsd[1176];

    auto g_y_z_0_0_xz_yz_0_xy = buffer_1100_ddsd[1177];

    auto g_y_z_0_0_xz_yz_0_xz = buffer_1100_ddsd[1178];

    auto g_y_z_0_0_xz_yz_0_yy = buffer_1100_ddsd[1179];

    auto g_y_z_0_0_xz_yz_0_yz = buffer_1100_ddsd[1180];

    auto g_y_z_0_0_xz_yz_0_zz = buffer_1100_ddsd[1181];

    auto g_y_z_0_0_xz_zz_0_xx = buffer_1100_ddsd[1182];

    auto g_y_z_0_0_xz_zz_0_xy = buffer_1100_ddsd[1183];

    auto g_y_z_0_0_xz_zz_0_xz = buffer_1100_ddsd[1184];

    auto g_y_z_0_0_xz_zz_0_yy = buffer_1100_ddsd[1185];

    auto g_y_z_0_0_xz_zz_0_yz = buffer_1100_ddsd[1186];

    auto g_y_z_0_0_xz_zz_0_zz = buffer_1100_ddsd[1187];

    auto g_y_z_0_0_yy_xx_0_xx = buffer_1100_ddsd[1188];

    auto g_y_z_0_0_yy_xx_0_xy = buffer_1100_ddsd[1189];

    auto g_y_z_0_0_yy_xx_0_xz = buffer_1100_ddsd[1190];

    auto g_y_z_0_0_yy_xx_0_yy = buffer_1100_ddsd[1191];

    auto g_y_z_0_0_yy_xx_0_yz = buffer_1100_ddsd[1192];

    auto g_y_z_0_0_yy_xx_0_zz = buffer_1100_ddsd[1193];

    auto g_y_z_0_0_yy_xy_0_xx = buffer_1100_ddsd[1194];

    auto g_y_z_0_0_yy_xy_0_xy = buffer_1100_ddsd[1195];

    auto g_y_z_0_0_yy_xy_0_xz = buffer_1100_ddsd[1196];

    auto g_y_z_0_0_yy_xy_0_yy = buffer_1100_ddsd[1197];

    auto g_y_z_0_0_yy_xy_0_yz = buffer_1100_ddsd[1198];

    auto g_y_z_0_0_yy_xy_0_zz = buffer_1100_ddsd[1199];

    auto g_y_z_0_0_yy_xz_0_xx = buffer_1100_ddsd[1200];

    auto g_y_z_0_0_yy_xz_0_xy = buffer_1100_ddsd[1201];

    auto g_y_z_0_0_yy_xz_0_xz = buffer_1100_ddsd[1202];

    auto g_y_z_0_0_yy_xz_0_yy = buffer_1100_ddsd[1203];

    auto g_y_z_0_0_yy_xz_0_yz = buffer_1100_ddsd[1204];

    auto g_y_z_0_0_yy_xz_0_zz = buffer_1100_ddsd[1205];

    auto g_y_z_0_0_yy_yy_0_xx = buffer_1100_ddsd[1206];

    auto g_y_z_0_0_yy_yy_0_xy = buffer_1100_ddsd[1207];

    auto g_y_z_0_0_yy_yy_0_xz = buffer_1100_ddsd[1208];

    auto g_y_z_0_0_yy_yy_0_yy = buffer_1100_ddsd[1209];

    auto g_y_z_0_0_yy_yy_0_yz = buffer_1100_ddsd[1210];

    auto g_y_z_0_0_yy_yy_0_zz = buffer_1100_ddsd[1211];

    auto g_y_z_0_0_yy_yz_0_xx = buffer_1100_ddsd[1212];

    auto g_y_z_0_0_yy_yz_0_xy = buffer_1100_ddsd[1213];

    auto g_y_z_0_0_yy_yz_0_xz = buffer_1100_ddsd[1214];

    auto g_y_z_0_0_yy_yz_0_yy = buffer_1100_ddsd[1215];

    auto g_y_z_0_0_yy_yz_0_yz = buffer_1100_ddsd[1216];

    auto g_y_z_0_0_yy_yz_0_zz = buffer_1100_ddsd[1217];

    auto g_y_z_0_0_yy_zz_0_xx = buffer_1100_ddsd[1218];

    auto g_y_z_0_0_yy_zz_0_xy = buffer_1100_ddsd[1219];

    auto g_y_z_0_0_yy_zz_0_xz = buffer_1100_ddsd[1220];

    auto g_y_z_0_0_yy_zz_0_yy = buffer_1100_ddsd[1221];

    auto g_y_z_0_0_yy_zz_0_yz = buffer_1100_ddsd[1222];

    auto g_y_z_0_0_yy_zz_0_zz = buffer_1100_ddsd[1223];

    auto g_y_z_0_0_yz_xx_0_xx = buffer_1100_ddsd[1224];

    auto g_y_z_0_0_yz_xx_0_xy = buffer_1100_ddsd[1225];

    auto g_y_z_0_0_yz_xx_0_xz = buffer_1100_ddsd[1226];

    auto g_y_z_0_0_yz_xx_0_yy = buffer_1100_ddsd[1227];

    auto g_y_z_0_0_yz_xx_0_yz = buffer_1100_ddsd[1228];

    auto g_y_z_0_0_yz_xx_0_zz = buffer_1100_ddsd[1229];

    auto g_y_z_0_0_yz_xy_0_xx = buffer_1100_ddsd[1230];

    auto g_y_z_0_0_yz_xy_0_xy = buffer_1100_ddsd[1231];

    auto g_y_z_0_0_yz_xy_0_xz = buffer_1100_ddsd[1232];

    auto g_y_z_0_0_yz_xy_0_yy = buffer_1100_ddsd[1233];

    auto g_y_z_0_0_yz_xy_0_yz = buffer_1100_ddsd[1234];

    auto g_y_z_0_0_yz_xy_0_zz = buffer_1100_ddsd[1235];

    auto g_y_z_0_0_yz_xz_0_xx = buffer_1100_ddsd[1236];

    auto g_y_z_0_0_yz_xz_0_xy = buffer_1100_ddsd[1237];

    auto g_y_z_0_0_yz_xz_0_xz = buffer_1100_ddsd[1238];

    auto g_y_z_0_0_yz_xz_0_yy = buffer_1100_ddsd[1239];

    auto g_y_z_0_0_yz_xz_0_yz = buffer_1100_ddsd[1240];

    auto g_y_z_0_0_yz_xz_0_zz = buffer_1100_ddsd[1241];

    auto g_y_z_0_0_yz_yy_0_xx = buffer_1100_ddsd[1242];

    auto g_y_z_0_0_yz_yy_0_xy = buffer_1100_ddsd[1243];

    auto g_y_z_0_0_yz_yy_0_xz = buffer_1100_ddsd[1244];

    auto g_y_z_0_0_yz_yy_0_yy = buffer_1100_ddsd[1245];

    auto g_y_z_0_0_yz_yy_0_yz = buffer_1100_ddsd[1246];

    auto g_y_z_0_0_yz_yy_0_zz = buffer_1100_ddsd[1247];

    auto g_y_z_0_0_yz_yz_0_xx = buffer_1100_ddsd[1248];

    auto g_y_z_0_0_yz_yz_0_xy = buffer_1100_ddsd[1249];

    auto g_y_z_0_0_yz_yz_0_xz = buffer_1100_ddsd[1250];

    auto g_y_z_0_0_yz_yz_0_yy = buffer_1100_ddsd[1251];

    auto g_y_z_0_0_yz_yz_0_yz = buffer_1100_ddsd[1252];

    auto g_y_z_0_0_yz_yz_0_zz = buffer_1100_ddsd[1253];

    auto g_y_z_0_0_yz_zz_0_xx = buffer_1100_ddsd[1254];

    auto g_y_z_0_0_yz_zz_0_xy = buffer_1100_ddsd[1255];

    auto g_y_z_0_0_yz_zz_0_xz = buffer_1100_ddsd[1256];

    auto g_y_z_0_0_yz_zz_0_yy = buffer_1100_ddsd[1257];

    auto g_y_z_0_0_yz_zz_0_yz = buffer_1100_ddsd[1258];

    auto g_y_z_0_0_yz_zz_0_zz = buffer_1100_ddsd[1259];

    auto g_y_z_0_0_zz_xx_0_xx = buffer_1100_ddsd[1260];

    auto g_y_z_0_0_zz_xx_0_xy = buffer_1100_ddsd[1261];

    auto g_y_z_0_0_zz_xx_0_xz = buffer_1100_ddsd[1262];

    auto g_y_z_0_0_zz_xx_0_yy = buffer_1100_ddsd[1263];

    auto g_y_z_0_0_zz_xx_0_yz = buffer_1100_ddsd[1264];

    auto g_y_z_0_0_zz_xx_0_zz = buffer_1100_ddsd[1265];

    auto g_y_z_0_0_zz_xy_0_xx = buffer_1100_ddsd[1266];

    auto g_y_z_0_0_zz_xy_0_xy = buffer_1100_ddsd[1267];

    auto g_y_z_0_0_zz_xy_0_xz = buffer_1100_ddsd[1268];

    auto g_y_z_0_0_zz_xy_0_yy = buffer_1100_ddsd[1269];

    auto g_y_z_0_0_zz_xy_0_yz = buffer_1100_ddsd[1270];

    auto g_y_z_0_0_zz_xy_0_zz = buffer_1100_ddsd[1271];

    auto g_y_z_0_0_zz_xz_0_xx = buffer_1100_ddsd[1272];

    auto g_y_z_0_0_zz_xz_0_xy = buffer_1100_ddsd[1273];

    auto g_y_z_0_0_zz_xz_0_xz = buffer_1100_ddsd[1274];

    auto g_y_z_0_0_zz_xz_0_yy = buffer_1100_ddsd[1275];

    auto g_y_z_0_0_zz_xz_0_yz = buffer_1100_ddsd[1276];

    auto g_y_z_0_0_zz_xz_0_zz = buffer_1100_ddsd[1277];

    auto g_y_z_0_0_zz_yy_0_xx = buffer_1100_ddsd[1278];

    auto g_y_z_0_0_zz_yy_0_xy = buffer_1100_ddsd[1279];

    auto g_y_z_0_0_zz_yy_0_xz = buffer_1100_ddsd[1280];

    auto g_y_z_0_0_zz_yy_0_yy = buffer_1100_ddsd[1281];

    auto g_y_z_0_0_zz_yy_0_yz = buffer_1100_ddsd[1282];

    auto g_y_z_0_0_zz_yy_0_zz = buffer_1100_ddsd[1283];

    auto g_y_z_0_0_zz_yz_0_xx = buffer_1100_ddsd[1284];

    auto g_y_z_0_0_zz_yz_0_xy = buffer_1100_ddsd[1285];

    auto g_y_z_0_0_zz_yz_0_xz = buffer_1100_ddsd[1286];

    auto g_y_z_0_0_zz_yz_0_yy = buffer_1100_ddsd[1287];

    auto g_y_z_0_0_zz_yz_0_yz = buffer_1100_ddsd[1288];

    auto g_y_z_0_0_zz_yz_0_zz = buffer_1100_ddsd[1289];

    auto g_y_z_0_0_zz_zz_0_xx = buffer_1100_ddsd[1290];

    auto g_y_z_0_0_zz_zz_0_xy = buffer_1100_ddsd[1291];

    auto g_y_z_0_0_zz_zz_0_xz = buffer_1100_ddsd[1292];

    auto g_y_z_0_0_zz_zz_0_yy = buffer_1100_ddsd[1293];

    auto g_y_z_0_0_zz_zz_0_yz = buffer_1100_ddsd[1294];

    auto g_y_z_0_0_zz_zz_0_zz = buffer_1100_ddsd[1295];

    auto g_z_x_0_0_xx_xx_0_xx = buffer_1100_ddsd[1296];

    auto g_z_x_0_0_xx_xx_0_xy = buffer_1100_ddsd[1297];

    auto g_z_x_0_0_xx_xx_0_xz = buffer_1100_ddsd[1298];

    auto g_z_x_0_0_xx_xx_0_yy = buffer_1100_ddsd[1299];

    auto g_z_x_0_0_xx_xx_0_yz = buffer_1100_ddsd[1300];

    auto g_z_x_0_0_xx_xx_0_zz = buffer_1100_ddsd[1301];

    auto g_z_x_0_0_xx_xy_0_xx = buffer_1100_ddsd[1302];

    auto g_z_x_0_0_xx_xy_0_xy = buffer_1100_ddsd[1303];

    auto g_z_x_0_0_xx_xy_0_xz = buffer_1100_ddsd[1304];

    auto g_z_x_0_0_xx_xy_0_yy = buffer_1100_ddsd[1305];

    auto g_z_x_0_0_xx_xy_0_yz = buffer_1100_ddsd[1306];

    auto g_z_x_0_0_xx_xy_0_zz = buffer_1100_ddsd[1307];

    auto g_z_x_0_0_xx_xz_0_xx = buffer_1100_ddsd[1308];

    auto g_z_x_0_0_xx_xz_0_xy = buffer_1100_ddsd[1309];

    auto g_z_x_0_0_xx_xz_0_xz = buffer_1100_ddsd[1310];

    auto g_z_x_0_0_xx_xz_0_yy = buffer_1100_ddsd[1311];

    auto g_z_x_0_0_xx_xz_0_yz = buffer_1100_ddsd[1312];

    auto g_z_x_0_0_xx_xz_0_zz = buffer_1100_ddsd[1313];

    auto g_z_x_0_0_xx_yy_0_xx = buffer_1100_ddsd[1314];

    auto g_z_x_0_0_xx_yy_0_xy = buffer_1100_ddsd[1315];

    auto g_z_x_0_0_xx_yy_0_xz = buffer_1100_ddsd[1316];

    auto g_z_x_0_0_xx_yy_0_yy = buffer_1100_ddsd[1317];

    auto g_z_x_0_0_xx_yy_0_yz = buffer_1100_ddsd[1318];

    auto g_z_x_0_0_xx_yy_0_zz = buffer_1100_ddsd[1319];

    auto g_z_x_0_0_xx_yz_0_xx = buffer_1100_ddsd[1320];

    auto g_z_x_0_0_xx_yz_0_xy = buffer_1100_ddsd[1321];

    auto g_z_x_0_0_xx_yz_0_xz = buffer_1100_ddsd[1322];

    auto g_z_x_0_0_xx_yz_0_yy = buffer_1100_ddsd[1323];

    auto g_z_x_0_0_xx_yz_0_yz = buffer_1100_ddsd[1324];

    auto g_z_x_0_0_xx_yz_0_zz = buffer_1100_ddsd[1325];

    auto g_z_x_0_0_xx_zz_0_xx = buffer_1100_ddsd[1326];

    auto g_z_x_0_0_xx_zz_0_xy = buffer_1100_ddsd[1327];

    auto g_z_x_0_0_xx_zz_0_xz = buffer_1100_ddsd[1328];

    auto g_z_x_0_0_xx_zz_0_yy = buffer_1100_ddsd[1329];

    auto g_z_x_0_0_xx_zz_0_yz = buffer_1100_ddsd[1330];

    auto g_z_x_0_0_xx_zz_0_zz = buffer_1100_ddsd[1331];

    auto g_z_x_0_0_xy_xx_0_xx = buffer_1100_ddsd[1332];

    auto g_z_x_0_0_xy_xx_0_xy = buffer_1100_ddsd[1333];

    auto g_z_x_0_0_xy_xx_0_xz = buffer_1100_ddsd[1334];

    auto g_z_x_0_0_xy_xx_0_yy = buffer_1100_ddsd[1335];

    auto g_z_x_0_0_xy_xx_0_yz = buffer_1100_ddsd[1336];

    auto g_z_x_0_0_xy_xx_0_zz = buffer_1100_ddsd[1337];

    auto g_z_x_0_0_xy_xy_0_xx = buffer_1100_ddsd[1338];

    auto g_z_x_0_0_xy_xy_0_xy = buffer_1100_ddsd[1339];

    auto g_z_x_0_0_xy_xy_0_xz = buffer_1100_ddsd[1340];

    auto g_z_x_0_0_xy_xy_0_yy = buffer_1100_ddsd[1341];

    auto g_z_x_0_0_xy_xy_0_yz = buffer_1100_ddsd[1342];

    auto g_z_x_0_0_xy_xy_0_zz = buffer_1100_ddsd[1343];

    auto g_z_x_0_0_xy_xz_0_xx = buffer_1100_ddsd[1344];

    auto g_z_x_0_0_xy_xz_0_xy = buffer_1100_ddsd[1345];

    auto g_z_x_0_0_xy_xz_0_xz = buffer_1100_ddsd[1346];

    auto g_z_x_0_0_xy_xz_0_yy = buffer_1100_ddsd[1347];

    auto g_z_x_0_0_xy_xz_0_yz = buffer_1100_ddsd[1348];

    auto g_z_x_0_0_xy_xz_0_zz = buffer_1100_ddsd[1349];

    auto g_z_x_0_0_xy_yy_0_xx = buffer_1100_ddsd[1350];

    auto g_z_x_0_0_xy_yy_0_xy = buffer_1100_ddsd[1351];

    auto g_z_x_0_0_xy_yy_0_xz = buffer_1100_ddsd[1352];

    auto g_z_x_0_0_xy_yy_0_yy = buffer_1100_ddsd[1353];

    auto g_z_x_0_0_xy_yy_0_yz = buffer_1100_ddsd[1354];

    auto g_z_x_0_0_xy_yy_0_zz = buffer_1100_ddsd[1355];

    auto g_z_x_0_0_xy_yz_0_xx = buffer_1100_ddsd[1356];

    auto g_z_x_0_0_xy_yz_0_xy = buffer_1100_ddsd[1357];

    auto g_z_x_0_0_xy_yz_0_xz = buffer_1100_ddsd[1358];

    auto g_z_x_0_0_xy_yz_0_yy = buffer_1100_ddsd[1359];

    auto g_z_x_0_0_xy_yz_0_yz = buffer_1100_ddsd[1360];

    auto g_z_x_0_0_xy_yz_0_zz = buffer_1100_ddsd[1361];

    auto g_z_x_0_0_xy_zz_0_xx = buffer_1100_ddsd[1362];

    auto g_z_x_0_0_xy_zz_0_xy = buffer_1100_ddsd[1363];

    auto g_z_x_0_0_xy_zz_0_xz = buffer_1100_ddsd[1364];

    auto g_z_x_0_0_xy_zz_0_yy = buffer_1100_ddsd[1365];

    auto g_z_x_0_0_xy_zz_0_yz = buffer_1100_ddsd[1366];

    auto g_z_x_0_0_xy_zz_0_zz = buffer_1100_ddsd[1367];

    auto g_z_x_0_0_xz_xx_0_xx = buffer_1100_ddsd[1368];

    auto g_z_x_0_0_xz_xx_0_xy = buffer_1100_ddsd[1369];

    auto g_z_x_0_0_xz_xx_0_xz = buffer_1100_ddsd[1370];

    auto g_z_x_0_0_xz_xx_0_yy = buffer_1100_ddsd[1371];

    auto g_z_x_0_0_xz_xx_0_yz = buffer_1100_ddsd[1372];

    auto g_z_x_0_0_xz_xx_0_zz = buffer_1100_ddsd[1373];

    auto g_z_x_0_0_xz_xy_0_xx = buffer_1100_ddsd[1374];

    auto g_z_x_0_0_xz_xy_0_xy = buffer_1100_ddsd[1375];

    auto g_z_x_0_0_xz_xy_0_xz = buffer_1100_ddsd[1376];

    auto g_z_x_0_0_xz_xy_0_yy = buffer_1100_ddsd[1377];

    auto g_z_x_0_0_xz_xy_0_yz = buffer_1100_ddsd[1378];

    auto g_z_x_0_0_xz_xy_0_zz = buffer_1100_ddsd[1379];

    auto g_z_x_0_0_xz_xz_0_xx = buffer_1100_ddsd[1380];

    auto g_z_x_0_0_xz_xz_0_xy = buffer_1100_ddsd[1381];

    auto g_z_x_0_0_xz_xz_0_xz = buffer_1100_ddsd[1382];

    auto g_z_x_0_0_xz_xz_0_yy = buffer_1100_ddsd[1383];

    auto g_z_x_0_0_xz_xz_0_yz = buffer_1100_ddsd[1384];

    auto g_z_x_0_0_xz_xz_0_zz = buffer_1100_ddsd[1385];

    auto g_z_x_0_0_xz_yy_0_xx = buffer_1100_ddsd[1386];

    auto g_z_x_0_0_xz_yy_0_xy = buffer_1100_ddsd[1387];

    auto g_z_x_0_0_xz_yy_0_xz = buffer_1100_ddsd[1388];

    auto g_z_x_0_0_xz_yy_0_yy = buffer_1100_ddsd[1389];

    auto g_z_x_0_0_xz_yy_0_yz = buffer_1100_ddsd[1390];

    auto g_z_x_0_0_xz_yy_0_zz = buffer_1100_ddsd[1391];

    auto g_z_x_0_0_xz_yz_0_xx = buffer_1100_ddsd[1392];

    auto g_z_x_0_0_xz_yz_0_xy = buffer_1100_ddsd[1393];

    auto g_z_x_0_0_xz_yz_0_xz = buffer_1100_ddsd[1394];

    auto g_z_x_0_0_xz_yz_0_yy = buffer_1100_ddsd[1395];

    auto g_z_x_0_0_xz_yz_0_yz = buffer_1100_ddsd[1396];

    auto g_z_x_0_0_xz_yz_0_zz = buffer_1100_ddsd[1397];

    auto g_z_x_0_0_xz_zz_0_xx = buffer_1100_ddsd[1398];

    auto g_z_x_0_0_xz_zz_0_xy = buffer_1100_ddsd[1399];

    auto g_z_x_0_0_xz_zz_0_xz = buffer_1100_ddsd[1400];

    auto g_z_x_0_0_xz_zz_0_yy = buffer_1100_ddsd[1401];

    auto g_z_x_0_0_xz_zz_0_yz = buffer_1100_ddsd[1402];

    auto g_z_x_0_0_xz_zz_0_zz = buffer_1100_ddsd[1403];

    auto g_z_x_0_0_yy_xx_0_xx = buffer_1100_ddsd[1404];

    auto g_z_x_0_0_yy_xx_0_xy = buffer_1100_ddsd[1405];

    auto g_z_x_0_0_yy_xx_0_xz = buffer_1100_ddsd[1406];

    auto g_z_x_0_0_yy_xx_0_yy = buffer_1100_ddsd[1407];

    auto g_z_x_0_0_yy_xx_0_yz = buffer_1100_ddsd[1408];

    auto g_z_x_0_0_yy_xx_0_zz = buffer_1100_ddsd[1409];

    auto g_z_x_0_0_yy_xy_0_xx = buffer_1100_ddsd[1410];

    auto g_z_x_0_0_yy_xy_0_xy = buffer_1100_ddsd[1411];

    auto g_z_x_0_0_yy_xy_0_xz = buffer_1100_ddsd[1412];

    auto g_z_x_0_0_yy_xy_0_yy = buffer_1100_ddsd[1413];

    auto g_z_x_0_0_yy_xy_0_yz = buffer_1100_ddsd[1414];

    auto g_z_x_0_0_yy_xy_0_zz = buffer_1100_ddsd[1415];

    auto g_z_x_0_0_yy_xz_0_xx = buffer_1100_ddsd[1416];

    auto g_z_x_0_0_yy_xz_0_xy = buffer_1100_ddsd[1417];

    auto g_z_x_0_0_yy_xz_0_xz = buffer_1100_ddsd[1418];

    auto g_z_x_0_0_yy_xz_0_yy = buffer_1100_ddsd[1419];

    auto g_z_x_0_0_yy_xz_0_yz = buffer_1100_ddsd[1420];

    auto g_z_x_0_0_yy_xz_0_zz = buffer_1100_ddsd[1421];

    auto g_z_x_0_0_yy_yy_0_xx = buffer_1100_ddsd[1422];

    auto g_z_x_0_0_yy_yy_0_xy = buffer_1100_ddsd[1423];

    auto g_z_x_0_0_yy_yy_0_xz = buffer_1100_ddsd[1424];

    auto g_z_x_0_0_yy_yy_0_yy = buffer_1100_ddsd[1425];

    auto g_z_x_0_0_yy_yy_0_yz = buffer_1100_ddsd[1426];

    auto g_z_x_0_0_yy_yy_0_zz = buffer_1100_ddsd[1427];

    auto g_z_x_0_0_yy_yz_0_xx = buffer_1100_ddsd[1428];

    auto g_z_x_0_0_yy_yz_0_xy = buffer_1100_ddsd[1429];

    auto g_z_x_0_0_yy_yz_0_xz = buffer_1100_ddsd[1430];

    auto g_z_x_0_0_yy_yz_0_yy = buffer_1100_ddsd[1431];

    auto g_z_x_0_0_yy_yz_0_yz = buffer_1100_ddsd[1432];

    auto g_z_x_0_0_yy_yz_0_zz = buffer_1100_ddsd[1433];

    auto g_z_x_0_0_yy_zz_0_xx = buffer_1100_ddsd[1434];

    auto g_z_x_0_0_yy_zz_0_xy = buffer_1100_ddsd[1435];

    auto g_z_x_0_0_yy_zz_0_xz = buffer_1100_ddsd[1436];

    auto g_z_x_0_0_yy_zz_0_yy = buffer_1100_ddsd[1437];

    auto g_z_x_0_0_yy_zz_0_yz = buffer_1100_ddsd[1438];

    auto g_z_x_0_0_yy_zz_0_zz = buffer_1100_ddsd[1439];

    auto g_z_x_0_0_yz_xx_0_xx = buffer_1100_ddsd[1440];

    auto g_z_x_0_0_yz_xx_0_xy = buffer_1100_ddsd[1441];

    auto g_z_x_0_0_yz_xx_0_xz = buffer_1100_ddsd[1442];

    auto g_z_x_0_0_yz_xx_0_yy = buffer_1100_ddsd[1443];

    auto g_z_x_0_0_yz_xx_0_yz = buffer_1100_ddsd[1444];

    auto g_z_x_0_0_yz_xx_0_zz = buffer_1100_ddsd[1445];

    auto g_z_x_0_0_yz_xy_0_xx = buffer_1100_ddsd[1446];

    auto g_z_x_0_0_yz_xy_0_xy = buffer_1100_ddsd[1447];

    auto g_z_x_0_0_yz_xy_0_xz = buffer_1100_ddsd[1448];

    auto g_z_x_0_0_yz_xy_0_yy = buffer_1100_ddsd[1449];

    auto g_z_x_0_0_yz_xy_0_yz = buffer_1100_ddsd[1450];

    auto g_z_x_0_0_yz_xy_0_zz = buffer_1100_ddsd[1451];

    auto g_z_x_0_0_yz_xz_0_xx = buffer_1100_ddsd[1452];

    auto g_z_x_0_0_yz_xz_0_xy = buffer_1100_ddsd[1453];

    auto g_z_x_0_0_yz_xz_0_xz = buffer_1100_ddsd[1454];

    auto g_z_x_0_0_yz_xz_0_yy = buffer_1100_ddsd[1455];

    auto g_z_x_0_0_yz_xz_0_yz = buffer_1100_ddsd[1456];

    auto g_z_x_0_0_yz_xz_0_zz = buffer_1100_ddsd[1457];

    auto g_z_x_0_0_yz_yy_0_xx = buffer_1100_ddsd[1458];

    auto g_z_x_0_0_yz_yy_0_xy = buffer_1100_ddsd[1459];

    auto g_z_x_0_0_yz_yy_0_xz = buffer_1100_ddsd[1460];

    auto g_z_x_0_0_yz_yy_0_yy = buffer_1100_ddsd[1461];

    auto g_z_x_0_0_yz_yy_0_yz = buffer_1100_ddsd[1462];

    auto g_z_x_0_0_yz_yy_0_zz = buffer_1100_ddsd[1463];

    auto g_z_x_0_0_yz_yz_0_xx = buffer_1100_ddsd[1464];

    auto g_z_x_0_0_yz_yz_0_xy = buffer_1100_ddsd[1465];

    auto g_z_x_0_0_yz_yz_0_xz = buffer_1100_ddsd[1466];

    auto g_z_x_0_0_yz_yz_0_yy = buffer_1100_ddsd[1467];

    auto g_z_x_0_0_yz_yz_0_yz = buffer_1100_ddsd[1468];

    auto g_z_x_0_0_yz_yz_0_zz = buffer_1100_ddsd[1469];

    auto g_z_x_0_0_yz_zz_0_xx = buffer_1100_ddsd[1470];

    auto g_z_x_0_0_yz_zz_0_xy = buffer_1100_ddsd[1471];

    auto g_z_x_0_0_yz_zz_0_xz = buffer_1100_ddsd[1472];

    auto g_z_x_0_0_yz_zz_0_yy = buffer_1100_ddsd[1473];

    auto g_z_x_0_0_yz_zz_0_yz = buffer_1100_ddsd[1474];

    auto g_z_x_0_0_yz_zz_0_zz = buffer_1100_ddsd[1475];

    auto g_z_x_0_0_zz_xx_0_xx = buffer_1100_ddsd[1476];

    auto g_z_x_0_0_zz_xx_0_xy = buffer_1100_ddsd[1477];

    auto g_z_x_0_0_zz_xx_0_xz = buffer_1100_ddsd[1478];

    auto g_z_x_0_0_zz_xx_0_yy = buffer_1100_ddsd[1479];

    auto g_z_x_0_0_zz_xx_0_yz = buffer_1100_ddsd[1480];

    auto g_z_x_0_0_zz_xx_0_zz = buffer_1100_ddsd[1481];

    auto g_z_x_0_0_zz_xy_0_xx = buffer_1100_ddsd[1482];

    auto g_z_x_0_0_zz_xy_0_xy = buffer_1100_ddsd[1483];

    auto g_z_x_0_0_zz_xy_0_xz = buffer_1100_ddsd[1484];

    auto g_z_x_0_0_zz_xy_0_yy = buffer_1100_ddsd[1485];

    auto g_z_x_0_0_zz_xy_0_yz = buffer_1100_ddsd[1486];

    auto g_z_x_0_0_zz_xy_0_zz = buffer_1100_ddsd[1487];

    auto g_z_x_0_0_zz_xz_0_xx = buffer_1100_ddsd[1488];

    auto g_z_x_0_0_zz_xz_0_xy = buffer_1100_ddsd[1489];

    auto g_z_x_0_0_zz_xz_0_xz = buffer_1100_ddsd[1490];

    auto g_z_x_0_0_zz_xz_0_yy = buffer_1100_ddsd[1491];

    auto g_z_x_0_0_zz_xz_0_yz = buffer_1100_ddsd[1492];

    auto g_z_x_0_0_zz_xz_0_zz = buffer_1100_ddsd[1493];

    auto g_z_x_0_0_zz_yy_0_xx = buffer_1100_ddsd[1494];

    auto g_z_x_0_0_zz_yy_0_xy = buffer_1100_ddsd[1495];

    auto g_z_x_0_0_zz_yy_0_xz = buffer_1100_ddsd[1496];

    auto g_z_x_0_0_zz_yy_0_yy = buffer_1100_ddsd[1497];

    auto g_z_x_0_0_zz_yy_0_yz = buffer_1100_ddsd[1498];

    auto g_z_x_0_0_zz_yy_0_zz = buffer_1100_ddsd[1499];

    auto g_z_x_0_0_zz_yz_0_xx = buffer_1100_ddsd[1500];

    auto g_z_x_0_0_zz_yz_0_xy = buffer_1100_ddsd[1501];

    auto g_z_x_0_0_zz_yz_0_xz = buffer_1100_ddsd[1502];

    auto g_z_x_0_0_zz_yz_0_yy = buffer_1100_ddsd[1503];

    auto g_z_x_0_0_zz_yz_0_yz = buffer_1100_ddsd[1504];

    auto g_z_x_0_0_zz_yz_0_zz = buffer_1100_ddsd[1505];

    auto g_z_x_0_0_zz_zz_0_xx = buffer_1100_ddsd[1506];

    auto g_z_x_0_0_zz_zz_0_xy = buffer_1100_ddsd[1507];

    auto g_z_x_0_0_zz_zz_0_xz = buffer_1100_ddsd[1508];

    auto g_z_x_0_0_zz_zz_0_yy = buffer_1100_ddsd[1509];

    auto g_z_x_0_0_zz_zz_0_yz = buffer_1100_ddsd[1510];

    auto g_z_x_0_0_zz_zz_0_zz = buffer_1100_ddsd[1511];

    auto g_z_y_0_0_xx_xx_0_xx = buffer_1100_ddsd[1512];

    auto g_z_y_0_0_xx_xx_0_xy = buffer_1100_ddsd[1513];

    auto g_z_y_0_0_xx_xx_0_xz = buffer_1100_ddsd[1514];

    auto g_z_y_0_0_xx_xx_0_yy = buffer_1100_ddsd[1515];

    auto g_z_y_0_0_xx_xx_0_yz = buffer_1100_ddsd[1516];

    auto g_z_y_0_0_xx_xx_0_zz = buffer_1100_ddsd[1517];

    auto g_z_y_0_0_xx_xy_0_xx = buffer_1100_ddsd[1518];

    auto g_z_y_0_0_xx_xy_0_xy = buffer_1100_ddsd[1519];

    auto g_z_y_0_0_xx_xy_0_xz = buffer_1100_ddsd[1520];

    auto g_z_y_0_0_xx_xy_0_yy = buffer_1100_ddsd[1521];

    auto g_z_y_0_0_xx_xy_0_yz = buffer_1100_ddsd[1522];

    auto g_z_y_0_0_xx_xy_0_zz = buffer_1100_ddsd[1523];

    auto g_z_y_0_0_xx_xz_0_xx = buffer_1100_ddsd[1524];

    auto g_z_y_0_0_xx_xz_0_xy = buffer_1100_ddsd[1525];

    auto g_z_y_0_0_xx_xz_0_xz = buffer_1100_ddsd[1526];

    auto g_z_y_0_0_xx_xz_0_yy = buffer_1100_ddsd[1527];

    auto g_z_y_0_0_xx_xz_0_yz = buffer_1100_ddsd[1528];

    auto g_z_y_0_0_xx_xz_0_zz = buffer_1100_ddsd[1529];

    auto g_z_y_0_0_xx_yy_0_xx = buffer_1100_ddsd[1530];

    auto g_z_y_0_0_xx_yy_0_xy = buffer_1100_ddsd[1531];

    auto g_z_y_0_0_xx_yy_0_xz = buffer_1100_ddsd[1532];

    auto g_z_y_0_0_xx_yy_0_yy = buffer_1100_ddsd[1533];

    auto g_z_y_0_0_xx_yy_0_yz = buffer_1100_ddsd[1534];

    auto g_z_y_0_0_xx_yy_0_zz = buffer_1100_ddsd[1535];

    auto g_z_y_0_0_xx_yz_0_xx = buffer_1100_ddsd[1536];

    auto g_z_y_0_0_xx_yz_0_xy = buffer_1100_ddsd[1537];

    auto g_z_y_0_0_xx_yz_0_xz = buffer_1100_ddsd[1538];

    auto g_z_y_0_0_xx_yz_0_yy = buffer_1100_ddsd[1539];

    auto g_z_y_0_0_xx_yz_0_yz = buffer_1100_ddsd[1540];

    auto g_z_y_0_0_xx_yz_0_zz = buffer_1100_ddsd[1541];

    auto g_z_y_0_0_xx_zz_0_xx = buffer_1100_ddsd[1542];

    auto g_z_y_0_0_xx_zz_0_xy = buffer_1100_ddsd[1543];

    auto g_z_y_0_0_xx_zz_0_xz = buffer_1100_ddsd[1544];

    auto g_z_y_0_0_xx_zz_0_yy = buffer_1100_ddsd[1545];

    auto g_z_y_0_0_xx_zz_0_yz = buffer_1100_ddsd[1546];

    auto g_z_y_0_0_xx_zz_0_zz = buffer_1100_ddsd[1547];

    auto g_z_y_0_0_xy_xx_0_xx = buffer_1100_ddsd[1548];

    auto g_z_y_0_0_xy_xx_0_xy = buffer_1100_ddsd[1549];

    auto g_z_y_0_0_xy_xx_0_xz = buffer_1100_ddsd[1550];

    auto g_z_y_0_0_xy_xx_0_yy = buffer_1100_ddsd[1551];

    auto g_z_y_0_0_xy_xx_0_yz = buffer_1100_ddsd[1552];

    auto g_z_y_0_0_xy_xx_0_zz = buffer_1100_ddsd[1553];

    auto g_z_y_0_0_xy_xy_0_xx = buffer_1100_ddsd[1554];

    auto g_z_y_0_0_xy_xy_0_xy = buffer_1100_ddsd[1555];

    auto g_z_y_0_0_xy_xy_0_xz = buffer_1100_ddsd[1556];

    auto g_z_y_0_0_xy_xy_0_yy = buffer_1100_ddsd[1557];

    auto g_z_y_0_0_xy_xy_0_yz = buffer_1100_ddsd[1558];

    auto g_z_y_0_0_xy_xy_0_zz = buffer_1100_ddsd[1559];

    auto g_z_y_0_0_xy_xz_0_xx = buffer_1100_ddsd[1560];

    auto g_z_y_0_0_xy_xz_0_xy = buffer_1100_ddsd[1561];

    auto g_z_y_0_0_xy_xz_0_xz = buffer_1100_ddsd[1562];

    auto g_z_y_0_0_xy_xz_0_yy = buffer_1100_ddsd[1563];

    auto g_z_y_0_0_xy_xz_0_yz = buffer_1100_ddsd[1564];

    auto g_z_y_0_0_xy_xz_0_zz = buffer_1100_ddsd[1565];

    auto g_z_y_0_0_xy_yy_0_xx = buffer_1100_ddsd[1566];

    auto g_z_y_0_0_xy_yy_0_xy = buffer_1100_ddsd[1567];

    auto g_z_y_0_0_xy_yy_0_xz = buffer_1100_ddsd[1568];

    auto g_z_y_0_0_xy_yy_0_yy = buffer_1100_ddsd[1569];

    auto g_z_y_0_0_xy_yy_0_yz = buffer_1100_ddsd[1570];

    auto g_z_y_0_0_xy_yy_0_zz = buffer_1100_ddsd[1571];

    auto g_z_y_0_0_xy_yz_0_xx = buffer_1100_ddsd[1572];

    auto g_z_y_0_0_xy_yz_0_xy = buffer_1100_ddsd[1573];

    auto g_z_y_0_0_xy_yz_0_xz = buffer_1100_ddsd[1574];

    auto g_z_y_0_0_xy_yz_0_yy = buffer_1100_ddsd[1575];

    auto g_z_y_0_0_xy_yz_0_yz = buffer_1100_ddsd[1576];

    auto g_z_y_0_0_xy_yz_0_zz = buffer_1100_ddsd[1577];

    auto g_z_y_0_0_xy_zz_0_xx = buffer_1100_ddsd[1578];

    auto g_z_y_0_0_xy_zz_0_xy = buffer_1100_ddsd[1579];

    auto g_z_y_0_0_xy_zz_0_xz = buffer_1100_ddsd[1580];

    auto g_z_y_0_0_xy_zz_0_yy = buffer_1100_ddsd[1581];

    auto g_z_y_0_0_xy_zz_0_yz = buffer_1100_ddsd[1582];

    auto g_z_y_0_0_xy_zz_0_zz = buffer_1100_ddsd[1583];

    auto g_z_y_0_0_xz_xx_0_xx = buffer_1100_ddsd[1584];

    auto g_z_y_0_0_xz_xx_0_xy = buffer_1100_ddsd[1585];

    auto g_z_y_0_0_xz_xx_0_xz = buffer_1100_ddsd[1586];

    auto g_z_y_0_0_xz_xx_0_yy = buffer_1100_ddsd[1587];

    auto g_z_y_0_0_xz_xx_0_yz = buffer_1100_ddsd[1588];

    auto g_z_y_0_0_xz_xx_0_zz = buffer_1100_ddsd[1589];

    auto g_z_y_0_0_xz_xy_0_xx = buffer_1100_ddsd[1590];

    auto g_z_y_0_0_xz_xy_0_xy = buffer_1100_ddsd[1591];

    auto g_z_y_0_0_xz_xy_0_xz = buffer_1100_ddsd[1592];

    auto g_z_y_0_0_xz_xy_0_yy = buffer_1100_ddsd[1593];

    auto g_z_y_0_0_xz_xy_0_yz = buffer_1100_ddsd[1594];

    auto g_z_y_0_0_xz_xy_0_zz = buffer_1100_ddsd[1595];

    auto g_z_y_0_0_xz_xz_0_xx = buffer_1100_ddsd[1596];

    auto g_z_y_0_0_xz_xz_0_xy = buffer_1100_ddsd[1597];

    auto g_z_y_0_0_xz_xz_0_xz = buffer_1100_ddsd[1598];

    auto g_z_y_0_0_xz_xz_0_yy = buffer_1100_ddsd[1599];

    auto g_z_y_0_0_xz_xz_0_yz = buffer_1100_ddsd[1600];

    auto g_z_y_0_0_xz_xz_0_zz = buffer_1100_ddsd[1601];

    auto g_z_y_0_0_xz_yy_0_xx = buffer_1100_ddsd[1602];

    auto g_z_y_0_0_xz_yy_0_xy = buffer_1100_ddsd[1603];

    auto g_z_y_0_0_xz_yy_0_xz = buffer_1100_ddsd[1604];

    auto g_z_y_0_0_xz_yy_0_yy = buffer_1100_ddsd[1605];

    auto g_z_y_0_0_xz_yy_0_yz = buffer_1100_ddsd[1606];

    auto g_z_y_0_0_xz_yy_0_zz = buffer_1100_ddsd[1607];

    auto g_z_y_0_0_xz_yz_0_xx = buffer_1100_ddsd[1608];

    auto g_z_y_0_0_xz_yz_0_xy = buffer_1100_ddsd[1609];

    auto g_z_y_0_0_xz_yz_0_xz = buffer_1100_ddsd[1610];

    auto g_z_y_0_0_xz_yz_0_yy = buffer_1100_ddsd[1611];

    auto g_z_y_0_0_xz_yz_0_yz = buffer_1100_ddsd[1612];

    auto g_z_y_0_0_xz_yz_0_zz = buffer_1100_ddsd[1613];

    auto g_z_y_0_0_xz_zz_0_xx = buffer_1100_ddsd[1614];

    auto g_z_y_0_0_xz_zz_0_xy = buffer_1100_ddsd[1615];

    auto g_z_y_0_0_xz_zz_0_xz = buffer_1100_ddsd[1616];

    auto g_z_y_0_0_xz_zz_0_yy = buffer_1100_ddsd[1617];

    auto g_z_y_0_0_xz_zz_0_yz = buffer_1100_ddsd[1618];

    auto g_z_y_0_0_xz_zz_0_zz = buffer_1100_ddsd[1619];

    auto g_z_y_0_0_yy_xx_0_xx = buffer_1100_ddsd[1620];

    auto g_z_y_0_0_yy_xx_0_xy = buffer_1100_ddsd[1621];

    auto g_z_y_0_0_yy_xx_0_xz = buffer_1100_ddsd[1622];

    auto g_z_y_0_0_yy_xx_0_yy = buffer_1100_ddsd[1623];

    auto g_z_y_0_0_yy_xx_0_yz = buffer_1100_ddsd[1624];

    auto g_z_y_0_0_yy_xx_0_zz = buffer_1100_ddsd[1625];

    auto g_z_y_0_0_yy_xy_0_xx = buffer_1100_ddsd[1626];

    auto g_z_y_0_0_yy_xy_0_xy = buffer_1100_ddsd[1627];

    auto g_z_y_0_0_yy_xy_0_xz = buffer_1100_ddsd[1628];

    auto g_z_y_0_0_yy_xy_0_yy = buffer_1100_ddsd[1629];

    auto g_z_y_0_0_yy_xy_0_yz = buffer_1100_ddsd[1630];

    auto g_z_y_0_0_yy_xy_0_zz = buffer_1100_ddsd[1631];

    auto g_z_y_0_0_yy_xz_0_xx = buffer_1100_ddsd[1632];

    auto g_z_y_0_0_yy_xz_0_xy = buffer_1100_ddsd[1633];

    auto g_z_y_0_0_yy_xz_0_xz = buffer_1100_ddsd[1634];

    auto g_z_y_0_0_yy_xz_0_yy = buffer_1100_ddsd[1635];

    auto g_z_y_0_0_yy_xz_0_yz = buffer_1100_ddsd[1636];

    auto g_z_y_0_0_yy_xz_0_zz = buffer_1100_ddsd[1637];

    auto g_z_y_0_0_yy_yy_0_xx = buffer_1100_ddsd[1638];

    auto g_z_y_0_0_yy_yy_0_xy = buffer_1100_ddsd[1639];

    auto g_z_y_0_0_yy_yy_0_xz = buffer_1100_ddsd[1640];

    auto g_z_y_0_0_yy_yy_0_yy = buffer_1100_ddsd[1641];

    auto g_z_y_0_0_yy_yy_0_yz = buffer_1100_ddsd[1642];

    auto g_z_y_0_0_yy_yy_0_zz = buffer_1100_ddsd[1643];

    auto g_z_y_0_0_yy_yz_0_xx = buffer_1100_ddsd[1644];

    auto g_z_y_0_0_yy_yz_0_xy = buffer_1100_ddsd[1645];

    auto g_z_y_0_0_yy_yz_0_xz = buffer_1100_ddsd[1646];

    auto g_z_y_0_0_yy_yz_0_yy = buffer_1100_ddsd[1647];

    auto g_z_y_0_0_yy_yz_0_yz = buffer_1100_ddsd[1648];

    auto g_z_y_0_0_yy_yz_0_zz = buffer_1100_ddsd[1649];

    auto g_z_y_0_0_yy_zz_0_xx = buffer_1100_ddsd[1650];

    auto g_z_y_0_0_yy_zz_0_xy = buffer_1100_ddsd[1651];

    auto g_z_y_0_0_yy_zz_0_xz = buffer_1100_ddsd[1652];

    auto g_z_y_0_0_yy_zz_0_yy = buffer_1100_ddsd[1653];

    auto g_z_y_0_0_yy_zz_0_yz = buffer_1100_ddsd[1654];

    auto g_z_y_0_0_yy_zz_0_zz = buffer_1100_ddsd[1655];

    auto g_z_y_0_0_yz_xx_0_xx = buffer_1100_ddsd[1656];

    auto g_z_y_0_0_yz_xx_0_xy = buffer_1100_ddsd[1657];

    auto g_z_y_0_0_yz_xx_0_xz = buffer_1100_ddsd[1658];

    auto g_z_y_0_0_yz_xx_0_yy = buffer_1100_ddsd[1659];

    auto g_z_y_0_0_yz_xx_0_yz = buffer_1100_ddsd[1660];

    auto g_z_y_0_0_yz_xx_0_zz = buffer_1100_ddsd[1661];

    auto g_z_y_0_0_yz_xy_0_xx = buffer_1100_ddsd[1662];

    auto g_z_y_0_0_yz_xy_0_xy = buffer_1100_ddsd[1663];

    auto g_z_y_0_0_yz_xy_0_xz = buffer_1100_ddsd[1664];

    auto g_z_y_0_0_yz_xy_0_yy = buffer_1100_ddsd[1665];

    auto g_z_y_0_0_yz_xy_0_yz = buffer_1100_ddsd[1666];

    auto g_z_y_0_0_yz_xy_0_zz = buffer_1100_ddsd[1667];

    auto g_z_y_0_0_yz_xz_0_xx = buffer_1100_ddsd[1668];

    auto g_z_y_0_0_yz_xz_0_xy = buffer_1100_ddsd[1669];

    auto g_z_y_0_0_yz_xz_0_xz = buffer_1100_ddsd[1670];

    auto g_z_y_0_0_yz_xz_0_yy = buffer_1100_ddsd[1671];

    auto g_z_y_0_0_yz_xz_0_yz = buffer_1100_ddsd[1672];

    auto g_z_y_0_0_yz_xz_0_zz = buffer_1100_ddsd[1673];

    auto g_z_y_0_0_yz_yy_0_xx = buffer_1100_ddsd[1674];

    auto g_z_y_0_0_yz_yy_0_xy = buffer_1100_ddsd[1675];

    auto g_z_y_0_0_yz_yy_0_xz = buffer_1100_ddsd[1676];

    auto g_z_y_0_0_yz_yy_0_yy = buffer_1100_ddsd[1677];

    auto g_z_y_0_0_yz_yy_0_yz = buffer_1100_ddsd[1678];

    auto g_z_y_0_0_yz_yy_0_zz = buffer_1100_ddsd[1679];

    auto g_z_y_0_0_yz_yz_0_xx = buffer_1100_ddsd[1680];

    auto g_z_y_0_0_yz_yz_0_xy = buffer_1100_ddsd[1681];

    auto g_z_y_0_0_yz_yz_0_xz = buffer_1100_ddsd[1682];

    auto g_z_y_0_0_yz_yz_0_yy = buffer_1100_ddsd[1683];

    auto g_z_y_0_0_yz_yz_0_yz = buffer_1100_ddsd[1684];

    auto g_z_y_0_0_yz_yz_0_zz = buffer_1100_ddsd[1685];

    auto g_z_y_0_0_yz_zz_0_xx = buffer_1100_ddsd[1686];

    auto g_z_y_0_0_yz_zz_0_xy = buffer_1100_ddsd[1687];

    auto g_z_y_0_0_yz_zz_0_xz = buffer_1100_ddsd[1688];

    auto g_z_y_0_0_yz_zz_0_yy = buffer_1100_ddsd[1689];

    auto g_z_y_0_0_yz_zz_0_yz = buffer_1100_ddsd[1690];

    auto g_z_y_0_0_yz_zz_0_zz = buffer_1100_ddsd[1691];

    auto g_z_y_0_0_zz_xx_0_xx = buffer_1100_ddsd[1692];

    auto g_z_y_0_0_zz_xx_0_xy = buffer_1100_ddsd[1693];

    auto g_z_y_0_0_zz_xx_0_xz = buffer_1100_ddsd[1694];

    auto g_z_y_0_0_zz_xx_0_yy = buffer_1100_ddsd[1695];

    auto g_z_y_0_0_zz_xx_0_yz = buffer_1100_ddsd[1696];

    auto g_z_y_0_0_zz_xx_0_zz = buffer_1100_ddsd[1697];

    auto g_z_y_0_0_zz_xy_0_xx = buffer_1100_ddsd[1698];

    auto g_z_y_0_0_zz_xy_0_xy = buffer_1100_ddsd[1699];

    auto g_z_y_0_0_zz_xy_0_xz = buffer_1100_ddsd[1700];

    auto g_z_y_0_0_zz_xy_0_yy = buffer_1100_ddsd[1701];

    auto g_z_y_0_0_zz_xy_0_yz = buffer_1100_ddsd[1702];

    auto g_z_y_0_0_zz_xy_0_zz = buffer_1100_ddsd[1703];

    auto g_z_y_0_0_zz_xz_0_xx = buffer_1100_ddsd[1704];

    auto g_z_y_0_0_zz_xz_0_xy = buffer_1100_ddsd[1705];

    auto g_z_y_0_0_zz_xz_0_xz = buffer_1100_ddsd[1706];

    auto g_z_y_0_0_zz_xz_0_yy = buffer_1100_ddsd[1707];

    auto g_z_y_0_0_zz_xz_0_yz = buffer_1100_ddsd[1708];

    auto g_z_y_0_0_zz_xz_0_zz = buffer_1100_ddsd[1709];

    auto g_z_y_0_0_zz_yy_0_xx = buffer_1100_ddsd[1710];

    auto g_z_y_0_0_zz_yy_0_xy = buffer_1100_ddsd[1711];

    auto g_z_y_0_0_zz_yy_0_xz = buffer_1100_ddsd[1712];

    auto g_z_y_0_0_zz_yy_0_yy = buffer_1100_ddsd[1713];

    auto g_z_y_0_0_zz_yy_0_yz = buffer_1100_ddsd[1714];

    auto g_z_y_0_0_zz_yy_0_zz = buffer_1100_ddsd[1715];

    auto g_z_y_0_0_zz_yz_0_xx = buffer_1100_ddsd[1716];

    auto g_z_y_0_0_zz_yz_0_xy = buffer_1100_ddsd[1717];

    auto g_z_y_0_0_zz_yz_0_xz = buffer_1100_ddsd[1718];

    auto g_z_y_0_0_zz_yz_0_yy = buffer_1100_ddsd[1719];

    auto g_z_y_0_0_zz_yz_0_yz = buffer_1100_ddsd[1720];

    auto g_z_y_0_0_zz_yz_0_zz = buffer_1100_ddsd[1721];

    auto g_z_y_0_0_zz_zz_0_xx = buffer_1100_ddsd[1722];

    auto g_z_y_0_0_zz_zz_0_xy = buffer_1100_ddsd[1723];

    auto g_z_y_0_0_zz_zz_0_xz = buffer_1100_ddsd[1724];

    auto g_z_y_0_0_zz_zz_0_yy = buffer_1100_ddsd[1725];

    auto g_z_y_0_0_zz_zz_0_yz = buffer_1100_ddsd[1726];

    auto g_z_y_0_0_zz_zz_0_zz = buffer_1100_ddsd[1727];

    auto g_z_z_0_0_xx_xx_0_xx = buffer_1100_ddsd[1728];

    auto g_z_z_0_0_xx_xx_0_xy = buffer_1100_ddsd[1729];

    auto g_z_z_0_0_xx_xx_0_xz = buffer_1100_ddsd[1730];

    auto g_z_z_0_0_xx_xx_0_yy = buffer_1100_ddsd[1731];

    auto g_z_z_0_0_xx_xx_0_yz = buffer_1100_ddsd[1732];

    auto g_z_z_0_0_xx_xx_0_zz = buffer_1100_ddsd[1733];

    auto g_z_z_0_0_xx_xy_0_xx = buffer_1100_ddsd[1734];

    auto g_z_z_0_0_xx_xy_0_xy = buffer_1100_ddsd[1735];

    auto g_z_z_0_0_xx_xy_0_xz = buffer_1100_ddsd[1736];

    auto g_z_z_0_0_xx_xy_0_yy = buffer_1100_ddsd[1737];

    auto g_z_z_0_0_xx_xy_0_yz = buffer_1100_ddsd[1738];

    auto g_z_z_0_0_xx_xy_0_zz = buffer_1100_ddsd[1739];

    auto g_z_z_0_0_xx_xz_0_xx = buffer_1100_ddsd[1740];

    auto g_z_z_0_0_xx_xz_0_xy = buffer_1100_ddsd[1741];

    auto g_z_z_0_0_xx_xz_0_xz = buffer_1100_ddsd[1742];

    auto g_z_z_0_0_xx_xz_0_yy = buffer_1100_ddsd[1743];

    auto g_z_z_0_0_xx_xz_0_yz = buffer_1100_ddsd[1744];

    auto g_z_z_0_0_xx_xz_0_zz = buffer_1100_ddsd[1745];

    auto g_z_z_0_0_xx_yy_0_xx = buffer_1100_ddsd[1746];

    auto g_z_z_0_0_xx_yy_0_xy = buffer_1100_ddsd[1747];

    auto g_z_z_0_0_xx_yy_0_xz = buffer_1100_ddsd[1748];

    auto g_z_z_0_0_xx_yy_0_yy = buffer_1100_ddsd[1749];

    auto g_z_z_0_0_xx_yy_0_yz = buffer_1100_ddsd[1750];

    auto g_z_z_0_0_xx_yy_0_zz = buffer_1100_ddsd[1751];

    auto g_z_z_0_0_xx_yz_0_xx = buffer_1100_ddsd[1752];

    auto g_z_z_0_0_xx_yz_0_xy = buffer_1100_ddsd[1753];

    auto g_z_z_0_0_xx_yz_0_xz = buffer_1100_ddsd[1754];

    auto g_z_z_0_0_xx_yz_0_yy = buffer_1100_ddsd[1755];

    auto g_z_z_0_0_xx_yz_0_yz = buffer_1100_ddsd[1756];

    auto g_z_z_0_0_xx_yz_0_zz = buffer_1100_ddsd[1757];

    auto g_z_z_0_0_xx_zz_0_xx = buffer_1100_ddsd[1758];

    auto g_z_z_0_0_xx_zz_0_xy = buffer_1100_ddsd[1759];

    auto g_z_z_0_0_xx_zz_0_xz = buffer_1100_ddsd[1760];

    auto g_z_z_0_0_xx_zz_0_yy = buffer_1100_ddsd[1761];

    auto g_z_z_0_0_xx_zz_0_yz = buffer_1100_ddsd[1762];

    auto g_z_z_0_0_xx_zz_0_zz = buffer_1100_ddsd[1763];

    auto g_z_z_0_0_xy_xx_0_xx = buffer_1100_ddsd[1764];

    auto g_z_z_0_0_xy_xx_0_xy = buffer_1100_ddsd[1765];

    auto g_z_z_0_0_xy_xx_0_xz = buffer_1100_ddsd[1766];

    auto g_z_z_0_0_xy_xx_0_yy = buffer_1100_ddsd[1767];

    auto g_z_z_0_0_xy_xx_0_yz = buffer_1100_ddsd[1768];

    auto g_z_z_0_0_xy_xx_0_zz = buffer_1100_ddsd[1769];

    auto g_z_z_0_0_xy_xy_0_xx = buffer_1100_ddsd[1770];

    auto g_z_z_0_0_xy_xy_0_xy = buffer_1100_ddsd[1771];

    auto g_z_z_0_0_xy_xy_0_xz = buffer_1100_ddsd[1772];

    auto g_z_z_0_0_xy_xy_0_yy = buffer_1100_ddsd[1773];

    auto g_z_z_0_0_xy_xy_0_yz = buffer_1100_ddsd[1774];

    auto g_z_z_0_0_xy_xy_0_zz = buffer_1100_ddsd[1775];

    auto g_z_z_0_0_xy_xz_0_xx = buffer_1100_ddsd[1776];

    auto g_z_z_0_0_xy_xz_0_xy = buffer_1100_ddsd[1777];

    auto g_z_z_0_0_xy_xz_0_xz = buffer_1100_ddsd[1778];

    auto g_z_z_0_0_xy_xz_0_yy = buffer_1100_ddsd[1779];

    auto g_z_z_0_0_xy_xz_0_yz = buffer_1100_ddsd[1780];

    auto g_z_z_0_0_xy_xz_0_zz = buffer_1100_ddsd[1781];

    auto g_z_z_0_0_xy_yy_0_xx = buffer_1100_ddsd[1782];

    auto g_z_z_0_0_xy_yy_0_xy = buffer_1100_ddsd[1783];

    auto g_z_z_0_0_xy_yy_0_xz = buffer_1100_ddsd[1784];

    auto g_z_z_0_0_xy_yy_0_yy = buffer_1100_ddsd[1785];

    auto g_z_z_0_0_xy_yy_0_yz = buffer_1100_ddsd[1786];

    auto g_z_z_0_0_xy_yy_0_zz = buffer_1100_ddsd[1787];

    auto g_z_z_0_0_xy_yz_0_xx = buffer_1100_ddsd[1788];

    auto g_z_z_0_0_xy_yz_0_xy = buffer_1100_ddsd[1789];

    auto g_z_z_0_0_xy_yz_0_xz = buffer_1100_ddsd[1790];

    auto g_z_z_0_0_xy_yz_0_yy = buffer_1100_ddsd[1791];

    auto g_z_z_0_0_xy_yz_0_yz = buffer_1100_ddsd[1792];

    auto g_z_z_0_0_xy_yz_0_zz = buffer_1100_ddsd[1793];

    auto g_z_z_0_0_xy_zz_0_xx = buffer_1100_ddsd[1794];

    auto g_z_z_0_0_xy_zz_0_xy = buffer_1100_ddsd[1795];

    auto g_z_z_0_0_xy_zz_0_xz = buffer_1100_ddsd[1796];

    auto g_z_z_0_0_xy_zz_0_yy = buffer_1100_ddsd[1797];

    auto g_z_z_0_0_xy_zz_0_yz = buffer_1100_ddsd[1798];

    auto g_z_z_0_0_xy_zz_0_zz = buffer_1100_ddsd[1799];

    auto g_z_z_0_0_xz_xx_0_xx = buffer_1100_ddsd[1800];

    auto g_z_z_0_0_xz_xx_0_xy = buffer_1100_ddsd[1801];

    auto g_z_z_0_0_xz_xx_0_xz = buffer_1100_ddsd[1802];

    auto g_z_z_0_0_xz_xx_0_yy = buffer_1100_ddsd[1803];

    auto g_z_z_0_0_xz_xx_0_yz = buffer_1100_ddsd[1804];

    auto g_z_z_0_0_xz_xx_0_zz = buffer_1100_ddsd[1805];

    auto g_z_z_0_0_xz_xy_0_xx = buffer_1100_ddsd[1806];

    auto g_z_z_0_0_xz_xy_0_xy = buffer_1100_ddsd[1807];

    auto g_z_z_0_0_xz_xy_0_xz = buffer_1100_ddsd[1808];

    auto g_z_z_0_0_xz_xy_0_yy = buffer_1100_ddsd[1809];

    auto g_z_z_0_0_xz_xy_0_yz = buffer_1100_ddsd[1810];

    auto g_z_z_0_0_xz_xy_0_zz = buffer_1100_ddsd[1811];

    auto g_z_z_0_0_xz_xz_0_xx = buffer_1100_ddsd[1812];

    auto g_z_z_0_0_xz_xz_0_xy = buffer_1100_ddsd[1813];

    auto g_z_z_0_0_xz_xz_0_xz = buffer_1100_ddsd[1814];

    auto g_z_z_0_0_xz_xz_0_yy = buffer_1100_ddsd[1815];

    auto g_z_z_0_0_xz_xz_0_yz = buffer_1100_ddsd[1816];

    auto g_z_z_0_0_xz_xz_0_zz = buffer_1100_ddsd[1817];

    auto g_z_z_0_0_xz_yy_0_xx = buffer_1100_ddsd[1818];

    auto g_z_z_0_0_xz_yy_0_xy = buffer_1100_ddsd[1819];

    auto g_z_z_0_0_xz_yy_0_xz = buffer_1100_ddsd[1820];

    auto g_z_z_0_0_xz_yy_0_yy = buffer_1100_ddsd[1821];

    auto g_z_z_0_0_xz_yy_0_yz = buffer_1100_ddsd[1822];

    auto g_z_z_0_0_xz_yy_0_zz = buffer_1100_ddsd[1823];

    auto g_z_z_0_0_xz_yz_0_xx = buffer_1100_ddsd[1824];

    auto g_z_z_0_0_xz_yz_0_xy = buffer_1100_ddsd[1825];

    auto g_z_z_0_0_xz_yz_0_xz = buffer_1100_ddsd[1826];

    auto g_z_z_0_0_xz_yz_0_yy = buffer_1100_ddsd[1827];

    auto g_z_z_0_0_xz_yz_0_yz = buffer_1100_ddsd[1828];

    auto g_z_z_0_0_xz_yz_0_zz = buffer_1100_ddsd[1829];

    auto g_z_z_0_0_xz_zz_0_xx = buffer_1100_ddsd[1830];

    auto g_z_z_0_0_xz_zz_0_xy = buffer_1100_ddsd[1831];

    auto g_z_z_0_0_xz_zz_0_xz = buffer_1100_ddsd[1832];

    auto g_z_z_0_0_xz_zz_0_yy = buffer_1100_ddsd[1833];

    auto g_z_z_0_0_xz_zz_0_yz = buffer_1100_ddsd[1834];

    auto g_z_z_0_0_xz_zz_0_zz = buffer_1100_ddsd[1835];

    auto g_z_z_0_0_yy_xx_0_xx = buffer_1100_ddsd[1836];

    auto g_z_z_0_0_yy_xx_0_xy = buffer_1100_ddsd[1837];

    auto g_z_z_0_0_yy_xx_0_xz = buffer_1100_ddsd[1838];

    auto g_z_z_0_0_yy_xx_0_yy = buffer_1100_ddsd[1839];

    auto g_z_z_0_0_yy_xx_0_yz = buffer_1100_ddsd[1840];

    auto g_z_z_0_0_yy_xx_0_zz = buffer_1100_ddsd[1841];

    auto g_z_z_0_0_yy_xy_0_xx = buffer_1100_ddsd[1842];

    auto g_z_z_0_0_yy_xy_0_xy = buffer_1100_ddsd[1843];

    auto g_z_z_0_0_yy_xy_0_xz = buffer_1100_ddsd[1844];

    auto g_z_z_0_0_yy_xy_0_yy = buffer_1100_ddsd[1845];

    auto g_z_z_0_0_yy_xy_0_yz = buffer_1100_ddsd[1846];

    auto g_z_z_0_0_yy_xy_0_zz = buffer_1100_ddsd[1847];

    auto g_z_z_0_0_yy_xz_0_xx = buffer_1100_ddsd[1848];

    auto g_z_z_0_0_yy_xz_0_xy = buffer_1100_ddsd[1849];

    auto g_z_z_0_0_yy_xz_0_xz = buffer_1100_ddsd[1850];

    auto g_z_z_0_0_yy_xz_0_yy = buffer_1100_ddsd[1851];

    auto g_z_z_0_0_yy_xz_0_yz = buffer_1100_ddsd[1852];

    auto g_z_z_0_0_yy_xz_0_zz = buffer_1100_ddsd[1853];

    auto g_z_z_0_0_yy_yy_0_xx = buffer_1100_ddsd[1854];

    auto g_z_z_0_0_yy_yy_0_xy = buffer_1100_ddsd[1855];

    auto g_z_z_0_0_yy_yy_0_xz = buffer_1100_ddsd[1856];

    auto g_z_z_0_0_yy_yy_0_yy = buffer_1100_ddsd[1857];

    auto g_z_z_0_0_yy_yy_0_yz = buffer_1100_ddsd[1858];

    auto g_z_z_0_0_yy_yy_0_zz = buffer_1100_ddsd[1859];

    auto g_z_z_0_0_yy_yz_0_xx = buffer_1100_ddsd[1860];

    auto g_z_z_0_0_yy_yz_0_xy = buffer_1100_ddsd[1861];

    auto g_z_z_0_0_yy_yz_0_xz = buffer_1100_ddsd[1862];

    auto g_z_z_0_0_yy_yz_0_yy = buffer_1100_ddsd[1863];

    auto g_z_z_0_0_yy_yz_0_yz = buffer_1100_ddsd[1864];

    auto g_z_z_0_0_yy_yz_0_zz = buffer_1100_ddsd[1865];

    auto g_z_z_0_0_yy_zz_0_xx = buffer_1100_ddsd[1866];

    auto g_z_z_0_0_yy_zz_0_xy = buffer_1100_ddsd[1867];

    auto g_z_z_0_0_yy_zz_0_xz = buffer_1100_ddsd[1868];

    auto g_z_z_0_0_yy_zz_0_yy = buffer_1100_ddsd[1869];

    auto g_z_z_0_0_yy_zz_0_yz = buffer_1100_ddsd[1870];

    auto g_z_z_0_0_yy_zz_0_zz = buffer_1100_ddsd[1871];

    auto g_z_z_0_0_yz_xx_0_xx = buffer_1100_ddsd[1872];

    auto g_z_z_0_0_yz_xx_0_xy = buffer_1100_ddsd[1873];

    auto g_z_z_0_0_yz_xx_0_xz = buffer_1100_ddsd[1874];

    auto g_z_z_0_0_yz_xx_0_yy = buffer_1100_ddsd[1875];

    auto g_z_z_0_0_yz_xx_0_yz = buffer_1100_ddsd[1876];

    auto g_z_z_0_0_yz_xx_0_zz = buffer_1100_ddsd[1877];

    auto g_z_z_0_0_yz_xy_0_xx = buffer_1100_ddsd[1878];

    auto g_z_z_0_0_yz_xy_0_xy = buffer_1100_ddsd[1879];

    auto g_z_z_0_0_yz_xy_0_xz = buffer_1100_ddsd[1880];

    auto g_z_z_0_0_yz_xy_0_yy = buffer_1100_ddsd[1881];

    auto g_z_z_0_0_yz_xy_0_yz = buffer_1100_ddsd[1882];

    auto g_z_z_0_0_yz_xy_0_zz = buffer_1100_ddsd[1883];

    auto g_z_z_0_0_yz_xz_0_xx = buffer_1100_ddsd[1884];

    auto g_z_z_0_0_yz_xz_0_xy = buffer_1100_ddsd[1885];

    auto g_z_z_0_0_yz_xz_0_xz = buffer_1100_ddsd[1886];

    auto g_z_z_0_0_yz_xz_0_yy = buffer_1100_ddsd[1887];

    auto g_z_z_0_0_yz_xz_0_yz = buffer_1100_ddsd[1888];

    auto g_z_z_0_0_yz_xz_0_zz = buffer_1100_ddsd[1889];

    auto g_z_z_0_0_yz_yy_0_xx = buffer_1100_ddsd[1890];

    auto g_z_z_0_0_yz_yy_0_xy = buffer_1100_ddsd[1891];

    auto g_z_z_0_0_yz_yy_0_xz = buffer_1100_ddsd[1892];

    auto g_z_z_0_0_yz_yy_0_yy = buffer_1100_ddsd[1893];

    auto g_z_z_0_0_yz_yy_0_yz = buffer_1100_ddsd[1894];

    auto g_z_z_0_0_yz_yy_0_zz = buffer_1100_ddsd[1895];

    auto g_z_z_0_0_yz_yz_0_xx = buffer_1100_ddsd[1896];

    auto g_z_z_0_0_yz_yz_0_xy = buffer_1100_ddsd[1897];

    auto g_z_z_0_0_yz_yz_0_xz = buffer_1100_ddsd[1898];

    auto g_z_z_0_0_yz_yz_0_yy = buffer_1100_ddsd[1899];

    auto g_z_z_0_0_yz_yz_0_yz = buffer_1100_ddsd[1900];

    auto g_z_z_0_0_yz_yz_0_zz = buffer_1100_ddsd[1901];

    auto g_z_z_0_0_yz_zz_0_xx = buffer_1100_ddsd[1902];

    auto g_z_z_0_0_yz_zz_0_xy = buffer_1100_ddsd[1903];

    auto g_z_z_0_0_yz_zz_0_xz = buffer_1100_ddsd[1904];

    auto g_z_z_0_0_yz_zz_0_yy = buffer_1100_ddsd[1905];

    auto g_z_z_0_0_yz_zz_0_yz = buffer_1100_ddsd[1906];

    auto g_z_z_0_0_yz_zz_0_zz = buffer_1100_ddsd[1907];

    auto g_z_z_0_0_zz_xx_0_xx = buffer_1100_ddsd[1908];

    auto g_z_z_0_0_zz_xx_0_xy = buffer_1100_ddsd[1909];

    auto g_z_z_0_0_zz_xx_0_xz = buffer_1100_ddsd[1910];

    auto g_z_z_0_0_zz_xx_0_yy = buffer_1100_ddsd[1911];

    auto g_z_z_0_0_zz_xx_0_yz = buffer_1100_ddsd[1912];

    auto g_z_z_0_0_zz_xx_0_zz = buffer_1100_ddsd[1913];

    auto g_z_z_0_0_zz_xy_0_xx = buffer_1100_ddsd[1914];

    auto g_z_z_0_0_zz_xy_0_xy = buffer_1100_ddsd[1915];

    auto g_z_z_0_0_zz_xy_0_xz = buffer_1100_ddsd[1916];

    auto g_z_z_0_0_zz_xy_0_yy = buffer_1100_ddsd[1917];

    auto g_z_z_0_0_zz_xy_0_yz = buffer_1100_ddsd[1918];

    auto g_z_z_0_0_zz_xy_0_zz = buffer_1100_ddsd[1919];

    auto g_z_z_0_0_zz_xz_0_xx = buffer_1100_ddsd[1920];

    auto g_z_z_0_0_zz_xz_0_xy = buffer_1100_ddsd[1921];

    auto g_z_z_0_0_zz_xz_0_xz = buffer_1100_ddsd[1922];

    auto g_z_z_0_0_zz_xz_0_yy = buffer_1100_ddsd[1923];

    auto g_z_z_0_0_zz_xz_0_yz = buffer_1100_ddsd[1924];

    auto g_z_z_0_0_zz_xz_0_zz = buffer_1100_ddsd[1925];

    auto g_z_z_0_0_zz_yy_0_xx = buffer_1100_ddsd[1926];

    auto g_z_z_0_0_zz_yy_0_xy = buffer_1100_ddsd[1927];

    auto g_z_z_0_0_zz_yy_0_xz = buffer_1100_ddsd[1928];

    auto g_z_z_0_0_zz_yy_0_yy = buffer_1100_ddsd[1929];

    auto g_z_z_0_0_zz_yy_0_yz = buffer_1100_ddsd[1930];

    auto g_z_z_0_0_zz_yy_0_zz = buffer_1100_ddsd[1931];

    auto g_z_z_0_0_zz_yz_0_xx = buffer_1100_ddsd[1932];

    auto g_z_z_0_0_zz_yz_0_xy = buffer_1100_ddsd[1933];

    auto g_z_z_0_0_zz_yz_0_xz = buffer_1100_ddsd[1934];

    auto g_z_z_0_0_zz_yz_0_yy = buffer_1100_ddsd[1935];

    auto g_z_z_0_0_zz_yz_0_yz = buffer_1100_ddsd[1936];

    auto g_z_z_0_0_zz_yz_0_zz = buffer_1100_ddsd[1937];

    auto g_z_z_0_0_zz_zz_0_xx = buffer_1100_ddsd[1938];

    auto g_z_z_0_0_zz_zz_0_xy = buffer_1100_ddsd[1939];

    auto g_z_z_0_0_zz_zz_0_xz = buffer_1100_ddsd[1940];

    auto g_z_z_0_0_zz_zz_0_yy = buffer_1100_ddsd[1941];

    auto g_z_z_0_0_zz_zz_0_yz = buffer_1100_ddsd[1942];

    auto g_z_z_0_0_zz_zz_0_zz = buffer_1100_ddsd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_xx_xx_0_xx, g_x_x_0_0_xx_xx_0_xy, g_x_x_0_0_xx_xx_0_xz, g_x_x_0_0_xx_xx_0_yy, g_x_x_0_0_xx_xx_0_yz, g_x_x_0_0_xx_xx_0_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xxx_0_xx, g_x_xxx_0_xy, g_x_xxx_0_xz, g_x_xxx_0_yy, g_x_xxx_0_yz, g_x_xxx_0_zz, g_xxx_x_0_xx, g_xxx_x_0_xy, g_xxx_x_0_xz, g_xxx_x_0_yy, g_xxx_x_0_yz, g_xxx_x_0_zz, g_xxx_xxx_0_xx, g_xxx_xxx_0_xy, g_xxx_xxx_0_xz, g_xxx_xxx_0_yy, g_xxx_xxx_0_yz, g_xxx_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_0_xx[i] = 4.0 * g_x_x_0_xx[i] - 4.0 * g_x_xxx_0_xx[i] * b_exp - 4.0 * g_xxx_x_0_xx[i] * a_exp + 4.0 * g_xxx_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_xy[i] = 4.0 * g_x_x_0_xy[i] - 4.0 * g_x_xxx_0_xy[i] * b_exp - 4.0 * g_xxx_x_0_xy[i] * a_exp + 4.0 * g_xxx_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_xz[i] = 4.0 * g_x_x_0_xz[i] - 4.0 * g_x_xxx_0_xz[i] * b_exp - 4.0 * g_xxx_x_0_xz[i] * a_exp + 4.0 * g_xxx_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_yy[i] = 4.0 * g_x_x_0_yy[i] - 4.0 * g_x_xxx_0_yy[i] * b_exp - 4.0 * g_xxx_x_0_yy[i] * a_exp + 4.0 * g_xxx_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_yz[i] = 4.0 * g_x_x_0_yz[i] - 4.0 * g_x_xxx_0_yz[i] * b_exp - 4.0 * g_xxx_x_0_yz[i] * a_exp + 4.0 * g_xxx_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_zz[i] = 4.0 * g_x_x_0_zz[i] - 4.0 * g_x_xxx_0_zz[i] * b_exp - 4.0 * g_xxx_x_0_zz[i] * a_exp + 4.0 * g_xxx_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_xx_xy_0_xx, g_x_x_0_0_xx_xy_0_xy, g_x_x_0_0_xx_xy_0_xz, g_x_x_0_0_xx_xy_0_yy, g_x_x_0_0_xx_xy_0_yz, g_x_x_0_0_xx_xy_0_zz, g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xxx_xxy_0_xx, g_xxx_xxy_0_xy, g_xxx_xxy_0_xz, g_xxx_xxy_0_yy, g_xxx_xxy_0_yz, g_xxx_xxy_0_zz, g_xxx_y_0_xx, g_xxx_y_0_xy, g_xxx_y_0_xz, g_xxx_y_0_yy, g_xxx_y_0_yz, g_xxx_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xy_0_xx[i] = 2.0 * g_x_y_0_xx[i] - 4.0 * g_x_xxy_0_xx[i] * b_exp - 2.0 * g_xxx_y_0_xx[i] * a_exp + 4.0 * g_xxx_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_xy[i] = 2.0 * g_x_y_0_xy[i] - 4.0 * g_x_xxy_0_xy[i] * b_exp - 2.0 * g_xxx_y_0_xy[i] * a_exp + 4.0 * g_xxx_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_xz[i] = 2.0 * g_x_y_0_xz[i] - 4.0 * g_x_xxy_0_xz[i] * b_exp - 2.0 * g_xxx_y_0_xz[i] * a_exp + 4.0 * g_xxx_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_yy[i] = 2.0 * g_x_y_0_yy[i] - 4.0 * g_x_xxy_0_yy[i] * b_exp - 2.0 * g_xxx_y_0_yy[i] * a_exp + 4.0 * g_xxx_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_yz[i] = 2.0 * g_x_y_0_yz[i] - 4.0 * g_x_xxy_0_yz[i] * b_exp - 2.0 * g_xxx_y_0_yz[i] * a_exp + 4.0 * g_xxx_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_zz[i] = 2.0 * g_x_y_0_zz[i] - 4.0 * g_x_xxy_0_zz[i] * b_exp - 2.0 * g_xxx_y_0_zz[i] * a_exp + 4.0 * g_xxx_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_xx_xz_0_xx, g_x_x_0_0_xx_xz_0_xy, g_x_x_0_0_xx_xz_0_xz, g_x_x_0_0_xx_xz_0_yy, g_x_x_0_0_xx_xz_0_yz, g_x_x_0_0_xx_xz_0_zz, g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xxx_xxz_0_xx, g_xxx_xxz_0_xy, g_xxx_xxz_0_xz, g_xxx_xxz_0_yy, g_xxx_xxz_0_yz, g_xxx_xxz_0_zz, g_xxx_z_0_xx, g_xxx_z_0_xy, g_xxx_z_0_xz, g_xxx_z_0_yy, g_xxx_z_0_yz, g_xxx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xz_0_xx[i] = 2.0 * g_x_z_0_xx[i] - 4.0 * g_x_xxz_0_xx[i] * b_exp - 2.0 * g_xxx_z_0_xx[i] * a_exp + 4.0 * g_xxx_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_xy[i] = 2.0 * g_x_z_0_xy[i] - 4.0 * g_x_xxz_0_xy[i] * b_exp - 2.0 * g_xxx_z_0_xy[i] * a_exp + 4.0 * g_xxx_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_xz[i] = 2.0 * g_x_z_0_xz[i] - 4.0 * g_x_xxz_0_xz[i] * b_exp - 2.0 * g_xxx_z_0_xz[i] * a_exp + 4.0 * g_xxx_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_yy[i] = 2.0 * g_x_z_0_yy[i] - 4.0 * g_x_xxz_0_yy[i] * b_exp - 2.0 * g_xxx_z_0_yy[i] * a_exp + 4.0 * g_xxx_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_yz[i] = 2.0 * g_x_z_0_yz[i] - 4.0 * g_x_xxz_0_yz[i] * b_exp - 2.0 * g_xxx_z_0_yz[i] * a_exp + 4.0 * g_xxx_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_zz[i] = 2.0 * g_x_z_0_zz[i] - 4.0 * g_x_xxz_0_zz[i] * b_exp - 2.0 * g_xxx_z_0_zz[i] * a_exp + 4.0 * g_xxx_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_xx_yy_0_xx, g_x_x_0_0_xx_yy_0_xy, g_x_x_0_0_xx_yy_0_xz, g_x_x_0_0_xx_yy_0_yy, g_x_x_0_0_xx_yy_0_yz, g_x_x_0_0_xx_yy_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_xxx_xyy_0_xx, g_xxx_xyy_0_xy, g_xxx_xyy_0_xz, g_xxx_xyy_0_yy, g_xxx_xyy_0_yz, g_xxx_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yy_0_xx[i] = -4.0 * g_x_xyy_0_xx[i] * b_exp + 4.0 * g_xxx_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_xy[i] = -4.0 * g_x_xyy_0_xy[i] * b_exp + 4.0 * g_xxx_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_xz[i] = -4.0 * g_x_xyy_0_xz[i] * b_exp + 4.0 * g_xxx_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_yy[i] = -4.0 * g_x_xyy_0_yy[i] * b_exp + 4.0 * g_xxx_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_yz[i] = -4.0 * g_x_xyy_0_yz[i] * b_exp + 4.0 * g_xxx_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_zz[i] = -4.0 * g_x_xyy_0_zz[i] * b_exp + 4.0 * g_xxx_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_xx_yz_0_xx, g_x_x_0_0_xx_yz_0_xy, g_x_x_0_0_xx_yz_0_xz, g_x_x_0_0_xx_yz_0_yy, g_x_x_0_0_xx_yz_0_yz, g_x_x_0_0_xx_yz_0_zz, g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xxx_xyz_0_xx, g_xxx_xyz_0_xy, g_xxx_xyz_0_xz, g_xxx_xyz_0_yy, g_xxx_xyz_0_yz, g_xxx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yz_0_xx[i] = -4.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xxx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_xy[i] = -4.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xxx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_xz[i] = -4.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xxx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_yy[i] = -4.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xxx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_yz[i] = -4.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xxx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_zz[i] = -4.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xxx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_xx_zz_0_xx, g_x_x_0_0_xx_zz_0_xy, g_x_x_0_0_xx_zz_0_xz, g_x_x_0_0_xx_zz_0_yy, g_x_x_0_0_xx_zz_0_yz, g_x_x_0_0_xx_zz_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_xxx_xzz_0_xx, g_xxx_xzz_0_xy, g_xxx_xzz_0_xz, g_xxx_xzz_0_yy, g_xxx_xzz_0_yz, g_xxx_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_zz_0_xx[i] = -4.0 * g_x_xzz_0_xx[i] * b_exp + 4.0 * g_xxx_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_xy[i] = -4.0 * g_x_xzz_0_xy[i] * b_exp + 4.0 * g_xxx_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_xz[i] = -4.0 * g_x_xzz_0_xz[i] * b_exp + 4.0 * g_xxx_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_yy[i] = -4.0 * g_x_xzz_0_yy[i] * b_exp + 4.0 * g_xxx_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_yz[i] = -4.0 * g_x_xzz_0_yz[i] * b_exp + 4.0 * g_xxx_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_zz[i] = -4.0 * g_x_xzz_0_zz[i] * b_exp + 4.0 * g_xxx_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_0_xx, g_x_x_0_0_xy_xx_0_xy, g_x_x_0_0_xy_xx_0_xz, g_x_x_0_0_xy_xx_0_yy, g_x_x_0_0_xy_xx_0_yz, g_x_x_0_0_xy_xx_0_zz, g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xxx_0_xx, g_xxy_xxx_0_xy, g_xxy_xxx_0_xz, g_xxy_xxx_0_yy, g_xxy_xxx_0_yz, g_xxy_xxx_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xxx_0_xx, g_y_xxx_0_xy, g_y_xxx_0_xz, g_y_xxx_0_yy, g_y_xxx_0_yz, g_y_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_0_xx[i] = 2.0 * g_y_x_0_xx[i] - 2.0 * g_y_xxx_0_xx[i] * b_exp - 4.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_xy[i] = 2.0 * g_y_x_0_xy[i] - 2.0 * g_y_xxx_0_xy[i] * b_exp - 4.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_xz[i] = 2.0 * g_y_x_0_xz[i] - 2.0 * g_y_xxx_0_xz[i] * b_exp - 4.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_yy[i] = 2.0 * g_y_x_0_yy[i] - 2.0 * g_y_xxx_0_yy[i] * b_exp - 4.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_yz[i] = 2.0 * g_y_x_0_yz[i] - 2.0 * g_y_xxx_0_yz[i] * b_exp - 4.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_zz[i] = 2.0 * g_y_x_0_zz[i] - 2.0 * g_y_xxx_0_zz[i] * b_exp - 4.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_xy_xy_0_xx, g_x_x_0_0_xy_xy_0_xy, g_x_x_0_0_xy_xy_0_xz, g_x_x_0_0_xy_xy_0_yy, g_x_x_0_0_xy_xy_0_yz, g_x_x_0_0_xy_xy_0_zz, g_xxy_xxy_0_xx, g_xxy_xxy_0_xy, g_xxy_xxy_0_xz, g_xxy_xxy_0_yy, g_xxy_xxy_0_yz, g_xxy_xxy_0_zz, g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xy_0_xx[i] = g_y_y_0_xx[i] - 2.0 * g_y_xxy_0_xx[i] * b_exp - 2.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_xy[i] = g_y_y_0_xy[i] - 2.0 * g_y_xxy_0_xy[i] * b_exp - 2.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_xz[i] = g_y_y_0_xz[i] - 2.0 * g_y_xxy_0_xz[i] * b_exp - 2.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_yy[i] = g_y_y_0_yy[i] - 2.0 * g_y_xxy_0_yy[i] * b_exp - 2.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_yz[i] = g_y_y_0_yz[i] - 2.0 * g_y_xxy_0_yz[i] * b_exp - 2.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_zz[i] = g_y_y_0_zz[i] - 2.0 * g_y_xxy_0_zz[i] * b_exp - 2.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_xy_xz_0_xx, g_x_x_0_0_xy_xz_0_xy, g_x_x_0_0_xy_xz_0_xz, g_x_x_0_0_xy_xz_0_yy, g_x_x_0_0_xy_xz_0_yz, g_x_x_0_0_xy_xz_0_zz, g_xxy_xxz_0_xx, g_xxy_xxz_0_xy, g_xxy_xxz_0_xz, g_xxy_xxz_0_yy, g_xxy_xxz_0_yz, g_xxy_xxz_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xz_0_xx[i] = g_y_z_0_xx[i] - 2.0 * g_y_xxz_0_xx[i] * b_exp - 2.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_xy[i] = g_y_z_0_xy[i] - 2.0 * g_y_xxz_0_xy[i] * b_exp - 2.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_xz[i] = g_y_z_0_xz[i] - 2.0 * g_y_xxz_0_xz[i] * b_exp - 2.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_yy[i] = g_y_z_0_yy[i] - 2.0 * g_y_xxz_0_yy[i] * b_exp - 2.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_yz[i] = g_y_z_0_yz[i] - 2.0 * g_y_xxz_0_yz[i] * b_exp - 2.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_zz[i] = g_y_z_0_zz[i] - 2.0 * g_y_xxz_0_zz[i] * b_exp - 2.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_xy_yy_0_xx, g_x_x_0_0_xy_yy_0_xy, g_x_x_0_0_xy_yy_0_xz, g_x_x_0_0_xy_yy_0_yy, g_x_x_0_0_xy_yy_0_yz, g_x_x_0_0_xy_yy_0_zz, g_xxy_xyy_0_xx, g_xxy_xyy_0_xy, g_xxy_xyy_0_xz, g_xxy_xyy_0_yy, g_xxy_xyy_0_yz, g_xxy_xyy_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yy_0_xx[i] = -2.0 * g_y_xyy_0_xx[i] * b_exp + 4.0 * g_xxy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_xy[i] = -2.0 * g_y_xyy_0_xy[i] * b_exp + 4.0 * g_xxy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_xz[i] = -2.0 * g_y_xyy_0_xz[i] * b_exp + 4.0 * g_xxy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_yy[i] = -2.0 * g_y_xyy_0_yy[i] * b_exp + 4.0 * g_xxy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_yz[i] = -2.0 * g_y_xyy_0_yz[i] * b_exp + 4.0 * g_xxy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_zz[i] = -2.0 * g_y_xyy_0_zz[i] * b_exp + 4.0 * g_xxy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_xy_yz_0_xx, g_x_x_0_0_xy_yz_0_xy, g_x_x_0_0_xy_yz_0_xz, g_x_x_0_0_xy_yz_0_yy, g_x_x_0_0_xy_yz_0_yz, g_x_x_0_0_xy_yz_0_zz, g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yz_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_xy_zz_0_xx, g_x_x_0_0_xy_zz_0_xy, g_x_x_0_0_xy_zz_0_xz, g_x_x_0_0_xy_zz_0_yy, g_x_x_0_0_xy_zz_0_yz, g_x_x_0_0_xy_zz_0_zz, g_xxy_xzz_0_xx, g_xxy_xzz_0_xy, g_xxy_xzz_0_xz, g_xxy_xzz_0_yy, g_xxy_xzz_0_yz, g_xxy_xzz_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_zz_0_xx[i] = -2.0 * g_y_xzz_0_xx[i] * b_exp + 4.0 * g_xxy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_xy[i] = -2.0 * g_y_xzz_0_xy[i] * b_exp + 4.0 * g_xxy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_xz[i] = -2.0 * g_y_xzz_0_xz[i] * b_exp + 4.0 * g_xxy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_yy[i] = -2.0 * g_y_xzz_0_yy[i] * b_exp + 4.0 * g_xxy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_yz[i] = -2.0 * g_y_xzz_0_yz[i] * b_exp + 4.0 * g_xxy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_zz[i] = -2.0 * g_y_xzz_0_zz[i] * b_exp + 4.0 * g_xxy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_0_xx, g_x_x_0_0_xz_xx_0_xy, g_x_x_0_0_xz_xx_0_xz, g_x_x_0_0_xz_xx_0_yy, g_x_x_0_0_xz_xx_0_yz, g_x_x_0_0_xz_xx_0_zz, g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xxx_0_xx, g_xxz_xxx_0_xy, g_xxz_xxx_0_xz, g_xxz_xxx_0_yy, g_xxz_xxx_0_yz, g_xxz_xxx_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xxx_0_xx, g_z_xxx_0_xy, g_z_xxx_0_xz, g_z_xxx_0_yy, g_z_xxx_0_yz, g_z_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_0_xx[i] = 2.0 * g_z_x_0_xx[i] - 2.0 * g_z_xxx_0_xx[i] * b_exp - 4.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_xy[i] = 2.0 * g_z_x_0_xy[i] - 2.0 * g_z_xxx_0_xy[i] * b_exp - 4.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_xz[i] = 2.0 * g_z_x_0_xz[i] - 2.0 * g_z_xxx_0_xz[i] * b_exp - 4.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_yy[i] = 2.0 * g_z_x_0_yy[i] - 2.0 * g_z_xxx_0_yy[i] * b_exp - 4.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_yz[i] = 2.0 * g_z_x_0_yz[i] - 2.0 * g_z_xxx_0_yz[i] * b_exp - 4.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_zz[i] = 2.0 * g_z_x_0_zz[i] - 2.0 * g_z_xxx_0_zz[i] * b_exp - 4.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_xz_xy_0_xx, g_x_x_0_0_xz_xy_0_xy, g_x_x_0_0_xz_xy_0_xz, g_x_x_0_0_xz_xy_0_yy, g_x_x_0_0_xz_xy_0_yz, g_x_x_0_0_xz_xy_0_zz, g_xxz_xxy_0_xx, g_xxz_xxy_0_xy, g_xxz_xxy_0_xz, g_xxz_xxy_0_yy, g_xxz_xxy_0_yz, g_xxz_xxy_0_zz, g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xy_0_xx[i] = g_z_y_0_xx[i] - 2.0 * g_z_xxy_0_xx[i] * b_exp - 2.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_xy[i] = g_z_y_0_xy[i] - 2.0 * g_z_xxy_0_xy[i] * b_exp - 2.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_xz[i] = g_z_y_0_xz[i] - 2.0 * g_z_xxy_0_xz[i] * b_exp - 2.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_yy[i] = g_z_y_0_yy[i] - 2.0 * g_z_xxy_0_yy[i] * b_exp - 2.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_yz[i] = g_z_y_0_yz[i] - 2.0 * g_z_xxy_0_yz[i] * b_exp - 2.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_zz[i] = g_z_y_0_zz[i] - 2.0 * g_z_xxy_0_zz[i] * b_exp - 2.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_xz_xz_0_xx, g_x_x_0_0_xz_xz_0_xy, g_x_x_0_0_xz_xz_0_xz, g_x_x_0_0_xz_xz_0_yy, g_x_x_0_0_xz_xz_0_yz, g_x_x_0_0_xz_xz_0_zz, g_xxz_xxz_0_xx, g_xxz_xxz_0_xy, g_xxz_xxz_0_xz, g_xxz_xxz_0_yy, g_xxz_xxz_0_yz, g_xxz_xxz_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xz_0_xx[i] = g_z_z_0_xx[i] - 2.0 * g_z_xxz_0_xx[i] * b_exp - 2.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_xy[i] = g_z_z_0_xy[i] - 2.0 * g_z_xxz_0_xy[i] * b_exp - 2.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_xz[i] = g_z_z_0_xz[i] - 2.0 * g_z_xxz_0_xz[i] * b_exp - 2.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_yy[i] = g_z_z_0_yy[i] - 2.0 * g_z_xxz_0_yy[i] * b_exp - 2.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_yz[i] = g_z_z_0_yz[i] - 2.0 * g_z_xxz_0_yz[i] * b_exp - 2.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_zz[i] = g_z_z_0_zz[i] - 2.0 * g_z_xxz_0_zz[i] * b_exp - 2.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_xz_yy_0_xx, g_x_x_0_0_xz_yy_0_xy, g_x_x_0_0_xz_yy_0_xz, g_x_x_0_0_xz_yy_0_yy, g_x_x_0_0_xz_yy_0_yz, g_x_x_0_0_xz_yy_0_zz, g_xxz_xyy_0_xx, g_xxz_xyy_0_xy, g_xxz_xyy_0_xz, g_xxz_xyy_0_yy, g_xxz_xyy_0_yz, g_xxz_xyy_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yy_0_xx[i] = -2.0 * g_z_xyy_0_xx[i] * b_exp + 4.0 * g_xxz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_xy[i] = -2.0 * g_z_xyy_0_xy[i] * b_exp + 4.0 * g_xxz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_xz[i] = -2.0 * g_z_xyy_0_xz[i] * b_exp + 4.0 * g_xxz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_yy[i] = -2.0 * g_z_xyy_0_yy[i] * b_exp + 4.0 * g_xxz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_yz[i] = -2.0 * g_z_xyy_0_yz[i] * b_exp + 4.0 * g_xxz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_zz[i] = -2.0 * g_z_xyy_0_zz[i] * b_exp + 4.0 * g_xxz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_xz_yz_0_xx, g_x_x_0_0_xz_yz_0_xy, g_x_x_0_0_xz_yz_0_xz, g_x_x_0_0_xz_yz_0_yy, g_x_x_0_0_xz_yz_0_yz, g_x_x_0_0_xz_yz_0_zz, g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yz_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_xz_zz_0_xx, g_x_x_0_0_xz_zz_0_xy, g_x_x_0_0_xz_zz_0_xz, g_x_x_0_0_xz_zz_0_yy, g_x_x_0_0_xz_zz_0_yz, g_x_x_0_0_xz_zz_0_zz, g_xxz_xzz_0_xx, g_xxz_xzz_0_xy, g_xxz_xzz_0_xz, g_xxz_xzz_0_yy, g_xxz_xzz_0_yz, g_xxz_xzz_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_zz_0_xx[i] = -2.0 * g_z_xzz_0_xx[i] * b_exp + 4.0 * g_xxz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_xy[i] = -2.0 * g_z_xzz_0_xy[i] * b_exp + 4.0 * g_xxz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_xz[i] = -2.0 * g_z_xzz_0_xz[i] * b_exp + 4.0 * g_xxz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_yy[i] = -2.0 * g_z_xzz_0_yy[i] * b_exp + 4.0 * g_xxz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_yz[i] = -2.0 * g_z_xzz_0_yz[i] * b_exp + 4.0 * g_xxz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_zz[i] = -2.0 * g_z_xzz_0_zz[i] * b_exp + 4.0 * g_xxz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_0_xx, g_x_x_0_0_yy_xx_0_xy, g_x_x_0_0_yy_xx_0_xz, g_x_x_0_0_yy_xx_0_yy, g_x_x_0_0_yy_xx_0_yz, g_x_x_0_0_yy_xx_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xxx_0_xx, g_xyy_xxx_0_xy, g_xyy_xxx_0_xz, g_xyy_xxx_0_yy, g_xyy_xxx_0_yz, g_xyy_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_0_xx[i] = -4.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_xy[i] = -4.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_xz[i] = -4.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_yy[i] = -4.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_yz[i] = -4.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_zz[i] = -4.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_x_0_0_yy_xy_0_xx, g_x_x_0_0_yy_xy_0_xy, g_x_x_0_0_yy_xy_0_xz, g_x_x_0_0_yy_xy_0_yy, g_x_x_0_0_yy_xy_0_yz, g_x_x_0_0_yy_xy_0_zz, g_xyy_xxy_0_xx, g_xyy_xxy_0_xy, g_xyy_xxy_0_xz, g_xyy_xxy_0_yy, g_xyy_xxy_0_yz, g_xyy_xxy_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xy_0_xx[i] = -2.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_xy[i] = -2.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_xz[i] = -2.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_yy[i] = -2.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_yz[i] = -2.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_zz[i] = -2.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_x_0_0_yy_xz_0_xx, g_x_x_0_0_yy_xz_0_xy, g_x_x_0_0_yy_xz_0_xz, g_x_x_0_0_yy_xz_0_yy, g_x_x_0_0_yy_xz_0_yz, g_x_x_0_0_yy_xz_0_zz, g_xyy_xxz_0_xx, g_xyy_xxz_0_xy, g_xyy_xxz_0_xz, g_xyy_xxz_0_yy, g_xyy_xxz_0_yz, g_xyy_xxz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xz_0_xx[i] = -2.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_xy[i] = -2.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_xz[i] = -2.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_yy[i] = -2.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_yz[i] = -2.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_zz[i] = -2.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_x_0_0_yy_yy_0_xx, g_x_x_0_0_yy_yy_0_xy, g_x_x_0_0_yy_yy_0_xz, g_x_x_0_0_yy_yy_0_yy, g_x_x_0_0_yy_yy_0_yz, g_x_x_0_0_yy_yy_0_zz, g_xyy_xyy_0_xx, g_xyy_xyy_0_xy, g_xyy_xyy_0_xz, g_xyy_xyy_0_yy, g_xyy_xyy_0_yz, g_xyy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yy_0_xx[i] = 4.0 * g_xyy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_xy[i] = 4.0 * g_xyy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_xz[i] = 4.0 * g_xyy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_yy[i] = 4.0 * g_xyy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_yz[i] = 4.0 * g_xyy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_zz[i] = 4.0 * g_xyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_x_0_0_yy_yz_0_xx, g_x_x_0_0_yy_yz_0_xy, g_x_x_0_0_yy_yz_0_xz, g_x_x_0_0_yy_yz_0_yy, g_x_x_0_0_yy_yz_0_yz, g_x_x_0_0_yy_yz_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yz_0_xx[i] = 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_xy[i] = 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_xz[i] = 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_yy[i] = 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_yz[i] = 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_zz[i] = 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_x_0_0_yy_zz_0_xx, g_x_x_0_0_yy_zz_0_xy, g_x_x_0_0_yy_zz_0_xz, g_x_x_0_0_yy_zz_0_yy, g_x_x_0_0_yy_zz_0_yz, g_x_x_0_0_yy_zz_0_zz, g_xyy_xzz_0_xx, g_xyy_xzz_0_xy, g_xyy_xzz_0_xz, g_xyy_xzz_0_yy, g_xyy_xzz_0_yz, g_xyy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_zz_0_xx[i] = 4.0 * g_xyy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_xy[i] = 4.0 * g_xyy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_xz[i] = 4.0 * g_xyy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_yy[i] = 4.0 * g_xyy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_yz[i] = 4.0 * g_xyy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_zz[i] = 4.0 * g_xyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_0_xx, g_x_x_0_0_yz_xx_0_xy, g_x_x_0_0_yz_xx_0_xz, g_x_x_0_0_yz_xx_0_yy, g_x_x_0_0_yz_xx_0_yz, g_x_x_0_0_yz_xx_0_zz, g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xxx_0_xx, g_xyz_xxx_0_xy, g_xyz_xxx_0_xz, g_xyz_xxx_0_yy, g_xyz_xxx_0_yz, g_xyz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_0_xx[i] = -4.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_xy[i] = -4.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_xz[i] = -4.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_yy[i] = -4.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_yz[i] = -4.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_zz[i] = -4.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_x_0_0_yz_xy_0_xx, g_x_x_0_0_yz_xy_0_xy, g_x_x_0_0_yz_xy_0_xz, g_x_x_0_0_yz_xy_0_yy, g_x_x_0_0_yz_xy_0_yz, g_x_x_0_0_yz_xy_0_zz, g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xy_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_x_0_0_yz_xz_0_xx, g_x_x_0_0_yz_xz_0_xy, g_x_x_0_0_yz_xz_0_xz, g_x_x_0_0_yz_xz_0_yy, g_x_x_0_0_yz_xz_0_yz, g_x_x_0_0_yz_xz_0_zz, g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_x_0_0_yz_yy_0_xx, g_x_x_0_0_yz_yy_0_xy, g_x_x_0_0_yz_yy_0_xz, g_x_x_0_0_yz_yy_0_yy, g_x_x_0_0_yz_yy_0_yz, g_x_x_0_0_yz_yy_0_zz, g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yy_0_xx[i] = 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_xy[i] = 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_xz[i] = 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_yy[i] = 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_yz[i] = 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_zz[i] = 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_x_0_0_yz_yz_0_xx, g_x_x_0_0_yz_yz_0_xy, g_x_x_0_0_yz_yz_0_xz, g_x_x_0_0_yz_yz_0_yy, g_x_x_0_0_yz_yz_0_yz, g_x_x_0_0_yz_yz_0_zz, g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_x_0_0_yz_zz_0_xx, g_x_x_0_0_yz_zz_0_xy, g_x_x_0_0_yz_zz_0_xz, g_x_x_0_0_yz_zz_0_yy, g_x_x_0_0_yz_zz_0_yz, g_x_x_0_0_yz_zz_0_zz, g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_zz_0_xx[i] = 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_xy[i] = 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_xz[i] = 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_yy[i] = 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_yz[i] = 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_zz[i] = 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_0_xx, g_x_x_0_0_zz_xx_0_xy, g_x_x_0_0_zz_xx_0_xz, g_x_x_0_0_zz_xx_0_yy, g_x_x_0_0_zz_xx_0_yz, g_x_x_0_0_zz_xx_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xxx_0_xx, g_xzz_xxx_0_xy, g_xzz_xxx_0_xz, g_xzz_xxx_0_yy, g_xzz_xxx_0_yz, g_xzz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_0_xx[i] = -4.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_xy[i] = -4.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_xz[i] = -4.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_yy[i] = -4.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_yz[i] = -4.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_zz[i] = -4.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_x_0_0_zz_xy_0_xx, g_x_x_0_0_zz_xy_0_xy, g_x_x_0_0_zz_xy_0_xz, g_x_x_0_0_zz_xy_0_yy, g_x_x_0_0_zz_xy_0_yz, g_x_x_0_0_zz_xy_0_zz, g_xzz_xxy_0_xx, g_xzz_xxy_0_xy, g_xzz_xxy_0_xz, g_xzz_xxy_0_yy, g_xzz_xxy_0_yz, g_xzz_xxy_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xy_0_xx[i] = -2.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_xy[i] = -2.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_xz[i] = -2.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_yy[i] = -2.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_yz[i] = -2.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_zz[i] = -2.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_x_0_0_zz_xz_0_xx, g_x_x_0_0_zz_xz_0_xy, g_x_x_0_0_zz_xz_0_xz, g_x_x_0_0_zz_xz_0_yy, g_x_x_0_0_zz_xz_0_yz, g_x_x_0_0_zz_xz_0_zz, g_xzz_xxz_0_xx, g_xzz_xxz_0_xy, g_xzz_xxz_0_xz, g_xzz_xxz_0_yy, g_xzz_xxz_0_yz, g_xzz_xxz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xz_0_xx[i] = -2.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_xy[i] = -2.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_xz[i] = -2.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_yy[i] = -2.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_yz[i] = -2.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_zz[i] = -2.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_x_0_0_zz_yy_0_xx, g_x_x_0_0_zz_yy_0_xy, g_x_x_0_0_zz_yy_0_xz, g_x_x_0_0_zz_yy_0_yy, g_x_x_0_0_zz_yy_0_yz, g_x_x_0_0_zz_yy_0_zz, g_xzz_xyy_0_xx, g_xzz_xyy_0_xy, g_xzz_xyy_0_xz, g_xzz_xyy_0_yy, g_xzz_xyy_0_yz, g_xzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yy_0_xx[i] = 4.0 * g_xzz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_xy[i] = 4.0 * g_xzz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_xz[i] = 4.0 * g_xzz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_yy[i] = 4.0 * g_xzz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_yz[i] = 4.0 * g_xzz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_zz[i] = 4.0 * g_xzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_x_0_0_zz_yz_0_xx, g_x_x_0_0_zz_yz_0_xy, g_x_x_0_0_zz_yz_0_xz, g_x_x_0_0_zz_yz_0_yy, g_x_x_0_0_zz_yz_0_yz, g_x_x_0_0_zz_yz_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yz_0_xx[i] = 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_xy[i] = 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_xz[i] = 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_yy[i] = 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_yz[i] = 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_zz[i] = 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_x_0_0_zz_zz_0_xx, g_x_x_0_0_zz_zz_0_xy, g_x_x_0_0_zz_zz_0_xz, g_x_x_0_0_zz_zz_0_yy, g_x_x_0_0_zz_zz_0_yz, g_x_x_0_0_zz_zz_0_zz, g_xzz_xzz_0_xx, g_xzz_xzz_0_xy, g_xzz_xzz_0_xz, g_xzz_xzz_0_yy, g_xzz_xzz_0_yz, g_xzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_zz_0_xx[i] = 4.0 * g_xzz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_xy[i] = 4.0 * g_xzz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_xz[i] = 4.0 * g_xzz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_yy[i] = 4.0 * g_xzz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_yz[i] = 4.0 * g_xzz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_zz[i] = 4.0 * g_xzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_0_xx_xx_0_xx, g_x_y_0_0_xx_xx_0_xy, g_x_y_0_0_xx_xx_0_xz, g_x_y_0_0_xx_xx_0_yy, g_x_y_0_0_xx_xx_0_yz, g_x_y_0_0_xx_xx_0_zz, g_xxx_xxy_0_xx, g_xxx_xxy_0_xy, g_xxx_xxy_0_xz, g_xxx_xxy_0_yy, g_xxx_xxy_0_yz, g_xxx_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_0_xx[i] = -4.0 * g_x_xxy_0_xx[i] * b_exp + 4.0 * g_xxx_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_xy[i] = -4.0 * g_x_xxy_0_xy[i] * b_exp + 4.0 * g_xxx_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_xz[i] = -4.0 * g_x_xxy_0_xz[i] * b_exp + 4.0 * g_xxx_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_yy[i] = -4.0 * g_x_xxy_0_yy[i] * b_exp + 4.0 * g_xxx_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_yz[i] = -4.0 * g_x_xxy_0_yz[i] * b_exp + 4.0 * g_xxx_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_zz[i] = -4.0 * g_x_xxy_0_zz[i] * b_exp + 4.0 * g_xxx_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_x_y_0_0_xx_xy_0_xx, g_x_y_0_0_xx_xy_0_xy, g_x_y_0_0_xx_xy_0_xz, g_x_y_0_0_xx_xy_0_yy, g_x_y_0_0_xx_xy_0_yz, g_x_y_0_0_xx_xy_0_zz, g_xxx_x_0_xx, g_xxx_x_0_xy, g_xxx_x_0_xz, g_xxx_x_0_yy, g_xxx_x_0_yz, g_xxx_x_0_zz, g_xxx_xyy_0_xx, g_xxx_xyy_0_xy, g_xxx_xyy_0_xz, g_xxx_xyy_0_yy, g_xxx_xyy_0_yz, g_xxx_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xy_0_xx[i] = 2.0 * g_x_x_0_xx[i] - 4.0 * g_x_xyy_0_xx[i] * b_exp - 2.0 * g_xxx_x_0_xx[i] * a_exp + 4.0 * g_xxx_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_xy[i] = 2.0 * g_x_x_0_xy[i] - 4.0 * g_x_xyy_0_xy[i] * b_exp - 2.0 * g_xxx_x_0_xy[i] * a_exp + 4.0 * g_xxx_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_xz[i] = 2.0 * g_x_x_0_xz[i] - 4.0 * g_x_xyy_0_xz[i] * b_exp - 2.0 * g_xxx_x_0_xz[i] * a_exp + 4.0 * g_xxx_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_yy[i] = 2.0 * g_x_x_0_yy[i] - 4.0 * g_x_xyy_0_yy[i] * b_exp - 2.0 * g_xxx_x_0_yy[i] * a_exp + 4.0 * g_xxx_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_yz[i] = 2.0 * g_x_x_0_yz[i] - 4.0 * g_x_xyy_0_yz[i] * b_exp - 2.0 * g_xxx_x_0_yz[i] * a_exp + 4.0 * g_xxx_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_zz[i] = 2.0 * g_x_x_0_zz[i] - 4.0 * g_x_xyy_0_zz[i] * b_exp - 2.0 * g_xxx_x_0_zz[i] * a_exp + 4.0 * g_xxx_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_x_y_0_0_xx_xz_0_xx, g_x_y_0_0_xx_xz_0_xy, g_x_y_0_0_xx_xz_0_xz, g_x_y_0_0_xx_xz_0_yy, g_x_y_0_0_xx_xz_0_yz, g_x_y_0_0_xx_xz_0_zz, g_xxx_xyz_0_xx, g_xxx_xyz_0_xy, g_xxx_xyz_0_xz, g_xxx_xyz_0_yy, g_xxx_xyz_0_yz, g_xxx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xz_0_xx[i] = -4.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xxx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_xy[i] = -4.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xxx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_xz[i] = -4.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xxx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_yy[i] = -4.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xxx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_yz[i] = -4.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xxx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_zz[i] = -4.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xxx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_y_0_0_xx_yy_0_xx, g_x_y_0_0_xx_yy_0_xy, g_x_y_0_0_xx_yy_0_xz, g_x_y_0_0_xx_yy_0_yy, g_x_y_0_0_xx_yy_0_yz, g_x_y_0_0_xx_yy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yyy_0_xx, g_x_yyy_0_xy, g_x_yyy_0_xz, g_x_yyy_0_yy, g_x_yyy_0_yz, g_x_yyy_0_zz, g_xxx_y_0_xx, g_xxx_y_0_xy, g_xxx_y_0_xz, g_xxx_y_0_yy, g_xxx_y_0_yz, g_xxx_y_0_zz, g_xxx_yyy_0_xx, g_xxx_yyy_0_xy, g_xxx_yyy_0_xz, g_xxx_yyy_0_yy, g_xxx_yyy_0_yz, g_xxx_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yy_0_xx[i] = 4.0 * g_x_y_0_xx[i] - 4.0 * g_x_yyy_0_xx[i] * b_exp - 4.0 * g_xxx_y_0_xx[i] * a_exp + 4.0 * g_xxx_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_xy[i] = 4.0 * g_x_y_0_xy[i] - 4.0 * g_x_yyy_0_xy[i] * b_exp - 4.0 * g_xxx_y_0_xy[i] * a_exp + 4.0 * g_xxx_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_xz[i] = 4.0 * g_x_y_0_xz[i] - 4.0 * g_x_yyy_0_xz[i] * b_exp - 4.0 * g_xxx_y_0_xz[i] * a_exp + 4.0 * g_xxx_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_yy[i] = 4.0 * g_x_y_0_yy[i] - 4.0 * g_x_yyy_0_yy[i] * b_exp - 4.0 * g_xxx_y_0_yy[i] * a_exp + 4.0 * g_xxx_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_yz[i] = 4.0 * g_x_y_0_yz[i] - 4.0 * g_x_yyy_0_yz[i] * b_exp - 4.0 * g_xxx_y_0_yz[i] * a_exp + 4.0 * g_xxx_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_zz[i] = 4.0 * g_x_y_0_zz[i] - 4.0 * g_x_yyy_0_zz[i] * b_exp - 4.0 * g_xxx_y_0_zz[i] * a_exp + 4.0 * g_xxx_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_y_0_0_xx_yz_0_xx, g_x_y_0_0_xx_yz_0_xy, g_x_y_0_0_xx_yz_0_xz, g_x_y_0_0_xx_yz_0_yy, g_x_y_0_0_xx_yz_0_yz, g_x_y_0_0_xx_yz_0_zz, g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xxx_yyz_0_xx, g_xxx_yyz_0_xy, g_xxx_yyz_0_xz, g_xxx_yyz_0_yy, g_xxx_yyz_0_yz, g_xxx_yyz_0_zz, g_xxx_z_0_xx, g_xxx_z_0_xy, g_xxx_z_0_xz, g_xxx_z_0_yy, g_xxx_z_0_yz, g_xxx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yz_0_xx[i] = 2.0 * g_x_z_0_xx[i] - 4.0 * g_x_yyz_0_xx[i] * b_exp - 2.0 * g_xxx_z_0_xx[i] * a_exp + 4.0 * g_xxx_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_xy[i] = 2.0 * g_x_z_0_xy[i] - 4.0 * g_x_yyz_0_xy[i] * b_exp - 2.0 * g_xxx_z_0_xy[i] * a_exp + 4.0 * g_xxx_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_xz[i] = 2.0 * g_x_z_0_xz[i] - 4.0 * g_x_yyz_0_xz[i] * b_exp - 2.0 * g_xxx_z_0_xz[i] * a_exp + 4.0 * g_xxx_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_yy[i] = 2.0 * g_x_z_0_yy[i] - 4.0 * g_x_yyz_0_yy[i] * b_exp - 2.0 * g_xxx_z_0_yy[i] * a_exp + 4.0 * g_xxx_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_yz[i] = 2.0 * g_x_z_0_yz[i] - 4.0 * g_x_yyz_0_yz[i] * b_exp - 2.0 * g_xxx_z_0_yz[i] * a_exp + 4.0 * g_xxx_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_zz[i] = 2.0 * g_x_z_0_zz[i] - 4.0 * g_x_yyz_0_zz[i] * b_exp - 2.0 * g_xxx_z_0_zz[i] * a_exp + 4.0 * g_xxx_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_y_0_0_xx_zz_0_xx, g_x_y_0_0_xx_zz_0_xy, g_x_y_0_0_xx_zz_0_xz, g_x_y_0_0_xx_zz_0_yy, g_x_y_0_0_xx_zz_0_yz, g_x_y_0_0_xx_zz_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_xxx_yzz_0_xx, g_xxx_yzz_0_xy, g_xxx_yzz_0_xz, g_xxx_yzz_0_yy, g_xxx_yzz_0_yz, g_xxx_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_zz_0_xx[i] = -4.0 * g_x_yzz_0_xx[i] * b_exp + 4.0 * g_xxx_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_xy[i] = -4.0 * g_x_yzz_0_xy[i] * b_exp + 4.0 * g_xxx_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_xz[i] = -4.0 * g_x_yzz_0_xz[i] * b_exp + 4.0 * g_xxx_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_yy[i] = -4.0 * g_x_yzz_0_yy[i] * b_exp + 4.0 * g_xxx_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_yz[i] = -4.0 * g_x_yzz_0_yz[i] * b_exp + 4.0 * g_xxx_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_zz[i] = -4.0 * g_x_yzz_0_zz[i] * b_exp + 4.0 * g_xxx_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_0_xx, g_x_y_0_0_xy_xx_0_xy, g_x_y_0_0_xy_xx_0_xz, g_x_y_0_0_xy_xx_0_yy, g_x_y_0_0_xy_xx_0_yz, g_x_y_0_0_xy_xx_0_zz, g_xxy_xxy_0_xx, g_xxy_xxy_0_xy, g_xxy_xxy_0_xz, g_xxy_xxy_0_yy, g_xxy_xxy_0_yz, g_xxy_xxy_0_zz, g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_0_xx[i] = -2.0 * g_y_xxy_0_xx[i] * b_exp + 4.0 * g_xxy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_xy[i] = -2.0 * g_y_xxy_0_xy[i] * b_exp + 4.0 * g_xxy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_xz[i] = -2.0 * g_y_xxy_0_xz[i] * b_exp + 4.0 * g_xxy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_yy[i] = -2.0 * g_y_xxy_0_yy[i] * b_exp + 4.0 * g_xxy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_yz[i] = -2.0 * g_y_xxy_0_yz[i] * b_exp + 4.0 * g_xxy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_zz[i] = -2.0 * g_y_xxy_0_zz[i] * b_exp + 4.0 * g_xxy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_y_0_0_xy_xy_0_xx, g_x_y_0_0_xy_xy_0_xy, g_x_y_0_0_xy_xy_0_xz, g_x_y_0_0_xy_xy_0_yy, g_x_y_0_0_xy_xy_0_yz, g_x_y_0_0_xy_xy_0_zz, g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xyy_0_xx, g_xxy_xyy_0_xy, g_xxy_xyy_0_xz, g_xxy_xyy_0_yy, g_xxy_xyy_0_yz, g_xxy_xyy_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xy_0_xx[i] = g_y_x_0_xx[i] - 2.0 * g_y_xyy_0_xx[i] * b_exp - 2.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_xy[i] = g_y_x_0_xy[i] - 2.0 * g_y_xyy_0_xy[i] * b_exp - 2.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_xz[i] = g_y_x_0_xz[i] - 2.0 * g_y_xyy_0_xz[i] * b_exp - 2.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_yy[i] = g_y_x_0_yy[i] - 2.0 * g_y_xyy_0_yy[i] * b_exp - 2.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_yz[i] = g_y_x_0_yz[i] - 2.0 * g_y_xyy_0_yz[i] * b_exp - 2.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_zz[i] = g_y_x_0_zz[i] - 2.0 * g_y_xyy_0_zz[i] * b_exp - 2.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_y_0_0_xy_xz_0_xx, g_x_y_0_0_xy_xz_0_xy, g_x_y_0_0_xy_xz_0_xz, g_x_y_0_0_xy_xz_0_yy, g_x_y_0_0_xy_xz_0_yz, g_x_y_0_0_xy_xz_0_zz, g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xz_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_y_0_0_xy_yy_0_xx, g_x_y_0_0_xy_yy_0_xy, g_x_y_0_0_xy_yy_0_xz, g_x_y_0_0_xy_yy_0_yy, g_x_y_0_0_xy_yy_0_yz, g_x_y_0_0_xy_yy_0_zz, g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_xxy_yyy_0_xx, g_xxy_yyy_0_xy, g_xxy_yyy_0_xz, g_xxy_yyy_0_yy, g_xxy_yyy_0_yz, g_xxy_yyy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yyy_0_xx, g_y_yyy_0_xy, g_y_yyy_0_xz, g_y_yyy_0_yy, g_y_yyy_0_yz, g_y_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yy_0_xx[i] = 2.0 * g_y_y_0_xx[i] - 2.0 * g_y_yyy_0_xx[i] * b_exp - 4.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_xy[i] = 2.0 * g_y_y_0_xy[i] - 2.0 * g_y_yyy_0_xy[i] * b_exp - 4.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_xz[i] = 2.0 * g_y_y_0_xz[i] - 2.0 * g_y_yyy_0_xz[i] * b_exp - 4.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_yy[i] = 2.0 * g_y_y_0_yy[i] - 2.0 * g_y_yyy_0_yy[i] * b_exp - 4.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_yz[i] = 2.0 * g_y_y_0_yz[i] - 2.0 * g_y_yyy_0_yz[i] * b_exp - 4.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_zz[i] = 2.0 * g_y_y_0_zz[i] - 2.0 * g_y_yyy_0_zz[i] * b_exp - 4.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_y_0_0_xy_yz_0_xx, g_x_y_0_0_xy_yz_0_xy, g_x_y_0_0_xy_yz_0_xz, g_x_y_0_0_xy_yz_0_yy, g_x_y_0_0_xy_yz_0_yz, g_x_y_0_0_xy_yz_0_zz, g_xxy_yyz_0_xx, g_xxy_yyz_0_xy, g_xxy_yyz_0_xz, g_xxy_yyz_0_yy, g_xxy_yyz_0_yz, g_xxy_yyz_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yz_0_xx[i] = g_y_z_0_xx[i] - 2.0 * g_y_yyz_0_xx[i] * b_exp - 2.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_xy[i] = g_y_z_0_xy[i] - 2.0 * g_y_yyz_0_xy[i] * b_exp - 2.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_xz[i] = g_y_z_0_xz[i] - 2.0 * g_y_yyz_0_xz[i] * b_exp - 2.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_yy[i] = g_y_z_0_yy[i] - 2.0 * g_y_yyz_0_yy[i] * b_exp - 2.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_yz[i] = g_y_z_0_yz[i] - 2.0 * g_y_yyz_0_yz[i] * b_exp - 2.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_zz[i] = g_y_z_0_zz[i] - 2.0 * g_y_yyz_0_zz[i] * b_exp - 2.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_y_0_0_xy_zz_0_xx, g_x_y_0_0_xy_zz_0_xy, g_x_y_0_0_xy_zz_0_xz, g_x_y_0_0_xy_zz_0_yy, g_x_y_0_0_xy_zz_0_yz, g_x_y_0_0_xy_zz_0_zz, g_xxy_yzz_0_xx, g_xxy_yzz_0_xy, g_xxy_yzz_0_xz, g_xxy_yzz_0_yy, g_xxy_yzz_0_yz, g_xxy_yzz_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_zz_0_xx[i] = -2.0 * g_y_yzz_0_xx[i] * b_exp + 4.0 * g_xxy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_xy[i] = -2.0 * g_y_yzz_0_xy[i] * b_exp + 4.0 * g_xxy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_xz[i] = -2.0 * g_y_yzz_0_xz[i] * b_exp + 4.0 * g_xxy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_yy[i] = -2.0 * g_y_yzz_0_yy[i] * b_exp + 4.0 * g_xxy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_yz[i] = -2.0 * g_y_yzz_0_yz[i] * b_exp + 4.0 * g_xxy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_zz[i] = -2.0 * g_y_yzz_0_zz[i] * b_exp + 4.0 * g_xxy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_0_xx, g_x_y_0_0_xz_xx_0_xy, g_x_y_0_0_xz_xx_0_xz, g_x_y_0_0_xz_xx_0_yy, g_x_y_0_0_xz_xx_0_yz, g_x_y_0_0_xz_xx_0_zz, g_xxz_xxy_0_xx, g_xxz_xxy_0_xy, g_xxz_xxy_0_xz, g_xxz_xxy_0_yy, g_xxz_xxy_0_yz, g_xxz_xxy_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_0_xx[i] = -2.0 * g_z_xxy_0_xx[i] * b_exp + 4.0 * g_xxz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_xy[i] = -2.0 * g_z_xxy_0_xy[i] * b_exp + 4.0 * g_xxz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_xz[i] = -2.0 * g_z_xxy_0_xz[i] * b_exp + 4.0 * g_xxz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_yy[i] = -2.0 * g_z_xxy_0_yy[i] * b_exp + 4.0 * g_xxz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_yz[i] = -2.0 * g_z_xxy_0_yz[i] * b_exp + 4.0 * g_xxz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_zz[i] = -2.0 * g_z_xxy_0_zz[i] * b_exp + 4.0 * g_xxz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_y_0_0_xz_xy_0_xx, g_x_y_0_0_xz_xy_0_xy, g_x_y_0_0_xz_xy_0_xz, g_x_y_0_0_xz_xy_0_yy, g_x_y_0_0_xz_xy_0_yz, g_x_y_0_0_xz_xy_0_zz, g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xyy_0_xx, g_xxz_xyy_0_xy, g_xxz_xyy_0_xz, g_xxz_xyy_0_yy, g_xxz_xyy_0_yz, g_xxz_xyy_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xy_0_xx[i] = g_z_x_0_xx[i] - 2.0 * g_z_xyy_0_xx[i] * b_exp - 2.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_xy[i] = g_z_x_0_xy[i] - 2.0 * g_z_xyy_0_xy[i] * b_exp - 2.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_xz[i] = g_z_x_0_xz[i] - 2.0 * g_z_xyy_0_xz[i] * b_exp - 2.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_yy[i] = g_z_x_0_yy[i] - 2.0 * g_z_xyy_0_yy[i] * b_exp - 2.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_yz[i] = g_z_x_0_yz[i] - 2.0 * g_z_xyy_0_yz[i] * b_exp - 2.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_zz[i] = g_z_x_0_zz[i] - 2.0 * g_z_xyy_0_zz[i] * b_exp - 2.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_y_0_0_xz_xz_0_xx, g_x_y_0_0_xz_xz_0_xy, g_x_y_0_0_xz_xz_0_xz, g_x_y_0_0_xz_xz_0_yy, g_x_y_0_0_xz_xz_0_yz, g_x_y_0_0_xz_xz_0_zz, g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xz_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_y_0_0_xz_yy_0_xx, g_x_y_0_0_xz_yy_0_xy, g_x_y_0_0_xz_yy_0_xz, g_x_y_0_0_xz_yy_0_yy, g_x_y_0_0_xz_yy_0_yz, g_x_y_0_0_xz_yy_0_zz, g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_xxz_yyy_0_xx, g_xxz_yyy_0_xy, g_xxz_yyy_0_xz, g_xxz_yyy_0_yy, g_xxz_yyy_0_yz, g_xxz_yyy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yyy_0_xx, g_z_yyy_0_xy, g_z_yyy_0_xz, g_z_yyy_0_yy, g_z_yyy_0_yz, g_z_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yy_0_xx[i] = 2.0 * g_z_y_0_xx[i] - 2.0 * g_z_yyy_0_xx[i] * b_exp - 4.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_xy[i] = 2.0 * g_z_y_0_xy[i] - 2.0 * g_z_yyy_0_xy[i] * b_exp - 4.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_xz[i] = 2.0 * g_z_y_0_xz[i] - 2.0 * g_z_yyy_0_xz[i] * b_exp - 4.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_yy[i] = 2.0 * g_z_y_0_yy[i] - 2.0 * g_z_yyy_0_yy[i] * b_exp - 4.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_yz[i] = 2.0 * g_z_y_0_yz[i] - 2.0 * g_z_yyy_0_yz[i] * b_exp - 4.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_zz[i] = 2.0 * g_z_y_0_zz[i] - 2.0 * g_z_yyy_0_zz[i] * b_exp - 4.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_y_0_0_xz_yz_0_xx, g_x_y_0_0_xz_yz_0_xy, g_x_y_0_0_xz_yz_0_xz, g_x_y_0_0_xz_yz_0_yy, g_x_y_0_0_xz_yz_0_yz, g_x_y_0_0_xz_yz_0_zz, g_xxz_yyz_0_xx, g_xxz_yyz_0_xy, g_xxz_yyz_0_xz, g_xxz_yyz_0_yy, g_xxz_yyz_0_yz, g_xxz_yyz_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yz_0_xx[i] = g_z_z_0_xx[i] - 2.0 * g_z_yyz_0_xx[i] * b_exp - 2.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_xy[i] = g_z_z_0_xy[i] - 2.0 * g_z_yyz_0_xy[i] * b_exp - 2.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_xz[i] = g_z_z_0_xz[i] - 2.0 * g_z_yyz_0_xz[i] * b_exp - 2.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_yy[i] = g_z_z_0_yy[i] - 2.0 * g_z_yyz_0_yy[i] * b_exp - 2.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_yz[i] = g_z_z_0_yz[i] - 2.0 * g_z_yyz_0_yz[i] * b_exp - 2.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_zz[i] = g_z_z_0_zz[i] - 2.0 * g_z_yyz_0_zz[i] * b_exp - 2.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_y_0_0_xz_zz_0_xx, g_x_y_0_0_xz_zz_0_xy, g_x_y_0_0_xz_zz_0_xz, g_x_y_0_0_xz_zz_0_yy, g_x_y_0_0_xz_zz_0_yz, g_x_y_0_0_xz_zz_0_zz, g_xxz_yzz_0_xx, g_xxz_yzz_0_xy, g_xxz_yzz_0_xz, g_xxz_yzz_0_yy, g_xxz_yzz_0_yz, g_xxz_yzz_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_zz_0_xx[i] = -2.0 * g_z_yzz_0_xx[i] * b_exp + 4.0 * g_xxz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_xy[i] = -2.0 * g_z_yzz_0_xy[i] * b_exp + 4.0 * g_xxz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_xz[i] = -2.0 * g_z_yzz_0_xz[i] * b_exp + 4.0 * g_xxz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_yy[i] = -2.0 * g_z_yzz_0_yy[i] * b_exp + 4.0 * g_xxz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_yz[i] = -2.0 * g_z_yzz_0_yz[i] * b_exp + 4.0 * g_xxz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_zz[i] = -2.0 * g_z_yzz_0_zz[i] * b_exp + 4.0 * g_xxz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_0_xx, g_x_y_0_0_yy_xx_0_xy, g_x_y_0_0_yy_xx_0_xz, g_x_y_0_0_yy_xx_0_yy, g_x_y_0_0_yy_xx_0_yz, g_x_y_0_0_yy_xx_0_zz, g_xyy_xxy_0_xx, g_xyy_xxy_0_xy, g_xyy_xxy_0_xz, g_xyy_xxy_0_yy, g_xyy_xxy_0_yz, g_xyy_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_0_xx[i] = 4.0 * g_xyy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_xy[i] = 4.0 * g_xyy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_xz[i] = 4.0 * g_xyy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_yy[i] = 4.0 * g_xyy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_yz[i] = 4.0 * g_xyy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_zz[i] = 4.0 * g_xyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_y_0_0_yy_xy_0_xx, g_x_y_0_0_yy_xy_0_xy, g_x_y_0_0_yy_xy_0_xz, g_x_y_0_0_yy_xy_0_yy, g_x_y_0_0_yy_xy_0_yz, g_x_y_0_0_yy_xy_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xyy_0_xx, g_xyy_xyy_0_xy, g_xyy_xyy_0_xz, g_xyy_xyy_0_yy, g_xyy_xyy_0_yz, g_xyy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xy_0_xx[i] = -2.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_xy[i] = -2.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_xz[i] = -2.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_yy[i] = -2.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_yz[i] = -2.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_zz[i] = -2.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_y_0_0_yy_xz_0_xx, g_x_y_0_0_yy_xz_0_xy, g_x_y_0_0_yy_xz_0_xz, g_x_y_0_0_yy_xz_0_yy, g_x_y_0_0_yy_xz_0_yz, g_x_y_0_0_yy_xz_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xz_0_xx[i] = 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_xy[i] = 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_xz[i] = 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_yy[i] = 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_yz[i] = 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_zz[i] = 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_y_0_0_yy_yy_0_xx, g_x_y_0_0_yy_yy_0_xy, g_x_y_0_0_yy_yy_0_xz, g_x_y_0_0_yy_yy_0_yy, g_x_y_0_0_yy_yy_0_yz, g_x_y_0_0_yy_yy_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_xyy_yyy_0_xx, g_xyy_yyy_0_xy, g_xyy_yyy_0_xz, g_xyy_yyy_0_yy, g_xyy_yyy_0_yz, g_xyy_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yy_0_xx[i] = -4.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_xy[i] = -4.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_xz[i] = -4.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_yy[i] = -4.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_yz[i] = -4.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_zz[i] = -4.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_y_0_0_yy_yz_0_xx, g_x_y_0_0_yy_yz_0_xy, g_x_y_0_0_yy_yz_0_xz, g_x_y_0_0_yy_yz_0_yy, g_x_y_0_0_yy_yz_0_yz, g_x_y_0_0_yy_yz_0_zz, g_xyy_yyz_0_xx, g_xyy_yyz_0_xy, g_xyy_yyz_0_xz, g_xyy_yyz_0_yy, g_xyy_yyz_0_yz, g_xyy_yyz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yz_0_xx[i] = -2.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_xy[i] = -2.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_xz[i] = -2.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_yy[i] = -2.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_yz[i] = -2.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_zz[i] = -2.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_y_0_0_yy_zz_0_xx, g_x_y_0_0_yy_zz_0_xy, g_x_y_0_0_yy_zz_0_xz, g_x_y_0_0_yy_zz_0_yy, g_x_y_0_0_yy_zz_0_yz, g_x_y_0_0_yy_zz_0_zz, g_xyy_yzz_0_xx, g_xyy_yzz_0_xy, g_xyy_yzz_0_xz, g_xyy_yzz_0_yy, g_xyy_yzz_0_yz, g_xyy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_zz_0_xx[i] = 4.0 * g_xyy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_xy[i] = 4.0 * g_xyy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_xz[i] = 4.0 * g_xyy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_yy[i] = 4.0 * g_xyy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_yz[i] = 4.0 * g_xyy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_zz[i] = 4.0 * g_xyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_0_xx, g_x_y_0_0_yz_xx_0_xy, g_x_y_0_0_yz_xx_0_xz, g_x_y_0_0_yz_xx_0_yy, g_x_y_0_0_yz_xx_0_yz, g_x_y_0_0_yz_xx_0_zz, g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_0_xx[i] = 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_xy[i] = 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_xz[i] = 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_yy[i] = 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_yz[i] = 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_zz[i] = 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_y_0_0_yz_xy_0_xx, g_x_y_0_0_yz_xy_0_xy, g_x_y_0_0_yz_xy_0_xz, g_x_y_0_0_yz_xy_0_yy, g_x_y_0_0_yz_xy_0_yz, g_x_y_0_0_yz_xy_0_zz, g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xy_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_y_0_0_yz_xz_0_xx, g_x_y_0_0_yz_xz_0_xy, g_x_y_0_0_yz_xz_0_xz, g_x_y_0_0_yz_xz_0_yy, g_x_y_0_0_yz_xz_0_yz, g_x_y_0_0_yz_xz_0_zz, g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_y_0_0_yz_yy_0_xx, g_x_y_0_0_yz_yy_0_xy, g_x_y_0_0_yz_yy_0_xz, g_x_y_0_0_yz_yy_0_yy, g_x_y_0_0_yz_yy_0_yz, g_x_y_0_0_yz_yy_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yyy_0_xx, g_xyz_yyy_0_xy, g_xyz_yyy_0_xz, g_xyz_yyy_0_yy, g_xyz_yyy_0_yz, g_xyz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yy_0_xx[i] = -4.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_xy[i] = -4.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_xz[i] = -4.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_yy[i] = -4.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_yz[i] = -4.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_zz[i] = -4.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_y_0_0_yz_yz_0_xx, g_x_y_0_0_yz_yz_0_xy, g_x_y_0_0_yz_yz_0_xz, g_x_y_0_0_yz_yz_0_yy, g_x_y_0_0_yz_yz_0_yz, g_x_y_0_0_yz_yz_0_zz, g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_y_0_0_yz_zz_0_xx, g_x_y_0_0_yz_zz_0_xy, g_x_y_0_0_yz_zz_0_xz, g_x_y_0_0_yz_zz_0_yy, g_x_y_0_0_yz_zz_0_yz, g_x_y_0_0_yz_zz_0_zz, g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_zz_0_xx[i] = 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_xy[i] = 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_xz[i] = 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_yy[i] = 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_yz[i] = 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_zz[i] = 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_0_xx, g_x_y_0_0_zz_xx_0_xy, g_x_y_0_0_zz_xx_0_xz, g_x_y_0_0_zz_xx_0_yy, g_x_y_0_0_zz_xx_0_yz, g_x_y_0_0_zz_xx_0_zz, g_xzz_xxy_0_xx, g_xzz_xxy_0_xy, g_xzz_xxy_0_xz, g_xzz_xxy_0_yy, g_xzz_xxy_0_yz, g_xzz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_0_xx[i] = 4.0 * g_xzz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_xy[i] = 4.0 * g_xzz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_xz[i] = 4.0 * g_xzz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_yy[i] = 4.0 * g_xzz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_yz[i] = 4.0 * g_xzz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_zz[i] = 4.0 * g_xzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_y_0_0_zz_xy_0_xx, g_x_y_0_0_zz_xy_0_xy, g_x_y_0_0_zz_xy_0_xz, g_x_y_0_0_zz_xy_0_yy, g_x_y_0_0_zz_xy_0_yz, g_x_y_0_0_zz_xy_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xyy_0_xx, g_xzz_xyy_0_xy, g_xzz_xyy_0_xz, g_xzz_xyy_0_yy, g_xzz_xyy_0_yz, g_xzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xy_0_xx[i] = -2.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_xy[i] = -2.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_xz[i] = -2.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_yy[i] = -2.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_yz[i] = -2.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_zz[i] = -2.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_y_0_0_zz_xz_0_xx, g_x_y_0_0_zz_xz_0_xy, g_x_y_0_0_zz_xz_0_xz, g_x_y_0_0_zz_xz_0_yy, g_x_y_0_0_zz_xz_0_yz, g_x_y_0_0_zz_xz_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xz_0_xx[i] = 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_xy[i] = 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_xz[i] = 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_yy[i] = 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_yz[i] = 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_zz[i] = 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_y_0_0_zz_yy_0_xx, g_x_y_0_0_zz_yy_0_xy, g_x_y_0_0_zz_yy_0_xz, g_x_y_0_0_zz_yy_0_yy, g_x_y_0_0_zz_yy_0_yz, g_x_y_0_0_zz_yy_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_xzz_yyy_0_xx, g_xzz_yyy_0_xy, g_xzz_yyy_0_xz, g_xzz_yyy_0_yy, g_xzz_yyy_0_yz, g_xzz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yy_0_xx[i] = -4.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_xy[i] = -4.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_xz[i] = -4.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_yy[i] = -4.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_yz[i] = -4.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_zz[i] = -4.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_y_0_0_zz_yz_0_xx, g_x_y_0_0_zz_yz_0_xy, g_x_y_0_0_zz_yz_0_xz, g_x_y_0_0_zz_yz_0_yy, g_x_y_0_0_zz_yz_0_yz, g_x_y_0_0_zz_yz_0_zz, g_xzz_yyz_0_xx, g_xzz_yyz_0_xy, g_xzz_yyz_0_xz, g_xzz_yyz_0_yy, g_xzz_yyz_0_yz, g_xzz_yyz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yz_0_xx[i] = -2.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_xy[i] = -2.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_xz[i] = -2.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_yy[i] = -2.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_yz[i] = -2.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_zz[i] = -2.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_y_0_0_zz_zz_0_xx, g_x_y_0_0_zz_zz_0_xy, g_x_y_0_0_zz_zz_0_xz, g_x_y_0_0_zz_zz_0_yy, g_x_y_0_0_zz_zz_0_yz, g_x_y_0_0_zz_zz_0_zz, g_xzz_yzz_0_xx, g_xzz_yzz_0_xy, g_xzz_yzz_0_xz, g_xzz_yzz_0_yy, g_xzz_yzz_0_yz, g_xzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_zz_0_xx[i] = 4.0 * g_xzz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_xy[i] = 4.0 * g_xzz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_xz[i] = 4.0 * g_xzz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_yy[i] = 4.0 * g_xzz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_yz[i] = 4.0 * g_xzz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_zz[i] = 4.0 * g_xzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_0_xx_xx_0_xx, g_x_z_0_0_xx_xx_0_xy, g_x_z_0_0_xx_xx_0_xz, g_x_z_0_0_xx_xx_0_yy, g_x_z_0_0_xx_xx_0_yz, g_x_z_0_0_xx_xx_0_zz, g_xxx_xxz_0_xx, g_xxx_xxz_0_xy, g_xxx_xxz_0_xz, g_xxx_xxz_0_yy, g_xxx_xxz_0_yz, g_xxx_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_0_xx[i] = -4.0 * g_x_xxz_0_xx[i] * b_exp + 4.0 * g_xxx_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_xy[i] = -4.0 * g_x_xxz_0_xy[i] * b_exp + 4.0 * g_xxx_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_xz[i] = -4.0 * g_x_xxz_0_xz[i] * b_exp + 4.0 * g_xxx_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_yy[i] = -4.0 * g_x_xxz_0_yy[i] * b_exp + 4.0 * g_xxx_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_yz[i] = -4.0 * g_x_xxz_0_yz[i] * b_exp + 4.0 * g_xxx_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_zz[i] = -4.0 * g_x_xxz_0_zz[i] * b_exp + 4.0 * g_xxx_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_x_z_0_0_xx_xy_0_xx, g_x_z_0_0_xx_xy_0_xy, g_x_z_0_0_xx_xy_0_xz, g_x_z_0_0_xx_xy_0_yy, g_x_z_0_0_xx_xy_0_yz, g_x_z_0_0_xx_xy_0_zz, g_xxx_xyz_0_xx, g_xxx_xyz_0_xy, g_xxx_xyz_0_xz, g_xxx_xyz_0_yy, g_xxx_xyz_0_yz, g_xxx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xy_0_xx[i] = -4.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xxx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_xy[i] = -4.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xxx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_xz[i] = -4.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xxx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_yy[i] = -4.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xxx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_yz[i] = -4.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xxx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_zz[i] = -4.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xxx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_x_z_0_0_xx_xz_0_xx, g_x_z_0_0_xx_xz_0_xy, g_x_z_0_0_xx_xz_0_xz, g_x_z_0_0_xx_xz_0_yy, g_x_z_0_0_xx_xz_0_yz, g_x_z_0_0_xx_xz_0_zz, g_xxx_x_0_xx, g_xxx_x_0_xy, g_xxx_x_0_xz, g_xxx_x_0_yy, g_xxx_x_0_yz, g_xxx_x_0_zz, g_xxx_xzz_0_xx, g_xxx_xzz_0_xy, g_xxx_xzz_0_xz, g_xxx_xzz_0_yy, g_xxx_xzz_0_yz, g_xxx_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xz_0_xx[i] = 2.0 * g_x_x_0_xx[i] - 4.0 * g_x_xzz_0_xx[i] * b_exp - 2.0 * g_xxx_x_0_xx[i] * a_exp + 4.0 * g_xxx_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_xy[i] = 2.0 * g_x_x_0_xy[i] - 4.0 * g_x_xzz_0_xy[i] * b_exp - 2.0 * g_xxx_x_0_xy[i] * a_exp + 4.0 * g_xxx_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_xz[i] = 2.0 * g_x_x_0_xz[i] - 4.0 * g_x_xzz_0_xz[i] * b_exp - 2.0 * g_xxx_x_0_xz[i] * a_exp + 4.0 * g_xxx_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_yy[i] = 2.0 * g_x_x_0_yy[i] - 4.0 * g_x_xzz_0_yy[i] * b_exp - 2.0 * g_xxx_x_0_yy[i] * a_exp + 4.0 * g_xxx_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_yz[i] = 2.0 * g_x_x_0_yz[i] - 4.0 * g_x_xzz_0_yz[i] * b_exp - 2.0 * g_xxx_x_0_yz[i] * a_exp + 4.0 * g_xxx_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_zz[i] = 2.0 * g_x_x_0_zz[i] - 4.0 * g_x_xzz_0_zz[i] * b_exp - 2.0 * g_xxx_x_0_zz[i] * a_exp + 4.0 * g_xxx_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_0_xx_yy_0_xx, g_x_z_0_0_xx_yy_0_xy, g_x_z_0_0_xx_yy_0_xz, g_x_z_0_0_xx_yy_0_yy, g_x_z_0_0_xx_yy_0_yz, g_x_z_0_0_xx_yy_0_zz, g_xxx_yyz_0_xx, g_xxx_yyz_0_xy, g_xxx_yyz_0_xz, g_xxx_yyz_0_yy, g_xxx_yyz_0_yz, g_xxx_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yy_0_xx[i] = -4.0 * g_x_yyz_0_xx[i] * b_exp + 4.0 * g_xxx_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_xy[i] = -4.0 * g_x_yyz_0_xy[i] * b_exp + 4.0 * g_xxx_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_xz[i] = -4.0 * g_x_yyz_0_xz[i] * b_exp + 4.0 * g_xxx_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_yy[i] = -4.0 * g_x_yyz_0_yy[i] * b_exp + 4.0 * g_xxx_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_yz[i] = -4.0 * g_x_yyz_0_yz[i] * b_exp + 4.0 * g_xxx_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_zz[i] = -4.0 * g_x_yyz_0_zz[i] * b_exp + 4.0 * g_xxx_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_x_z_0_0_xx_yz_0_xx, g_x_z_0_0_xx_yz_0_xy, g_x_z_0_0_xx_yz_0_xz, g_x_z_0_0_xx_yz_0_yy, g_x_z_0_0_xx_yz_0_yz, g_x_z_0_0_xx_yz_0_zz, g_xxx_y_0_xx, g_xxx_y_0_xy, g_xxx_y_0_xz, g_xxx_y_0_yy, g_xxx_y_0_yz, g_xxx_y_0_zz, g_xxx_yzz_0_xx, g_xxx_yzz_0_xy, g_xxx_yzz_0_xz, g_xxx_yzz_0_yy, g_xxx_yzz_0_yz, g_xxx_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yz_0_xx[i] = 2.0 * g_x_y_0_xx[i] - 4.0 * g_x_yzz_0_xx[i] * b_exp - 2.0 * g_xxx_y_0_xx[i] * a_exp + 4.0 * g_xxx_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_xy[i] = 2.0 * g_x_y_0_xy[i] - 4.0 * g_x_yzz_0_xy[i] * b_exp - 2.0 * g_xxx_y_0_xy[i] * a_exp + 4.0 * g_xxx_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_xz[i] = 2.0 * g_x_y_0_xz[i] - 4.0 * g_x_yzz_0_xz[i] * b_exp - 2.0 * g_xxx_y_0_xz[i] * a_exp + 4.0 * g_xxx_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_yy[i] = 2.0 * g_x_y_0_yy[i] - 4.0 * g_x_yzz_0_yy[i] * b_exp - 2.0 * g_xxx_y_0_yy[i] * a_exp + 4.0 * g_xxx_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_yz[i] = 2.0 * g_x_y_0_yz[i] - 4.0 * g_x_yzz_0_yz[i] * b_exp - 2.0 * g_xxx_y_0_yz[i] * a_exp + 4.0 * g_xxx_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_zz[i] = 2.0 * g_x_y_0_zz[i] - 4.0 * g_x_yzz_0_zz[i] * b_exp - 2.0 * g_xxx_y_0_zz[i] * a_exp + 4.0 * g_xxx_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_z_0_0_xx_zz_0_xx, g_x_z_0_0_xx_zz_0_xy, g_x_z_0_0_xx_zz_0_xz, g_x_z_0_0_xx_zz_0_yy, g_x_z_0_0_xx_zz_0_yz, g_x_z_0_0_xx_zz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_zzz_0_xx, g_x_zzz_0_xy, g_x_zzz_0_xz, g_x_zzz_0_yy, g_x_zzz_0_yz, g_x_zzz_0_zz, g_xxx_z_0_xx, g_xxx_z_0_xy, g_xxx_z_0_xz, g_xxx_z_0_yy, g_xxx_z_0_yz, g_xxx_z_0_zz, g_xxx_zzz_0_xx, g_xxx_zzz_0_xy, g_xxx_zzz_0_xz, g_xxx_zzz_0_yy, g_xxx_zzz_0_yz, g_xxx_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_zz_0_xx[i] = 4.0 * g_x_z_0_xx[i] - 4.0 * g_x_zzz_0_xx[i] * b_exp - 4.0 * g_xxx_z_0_xx[i] * a_exp + 4.0 * g_xxx_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_xy[i] = 4.0 * g_x_z_0_xy[i] - 4.0 * g_x_zzz_0_xy[i] * b_exp - 4.0 * g_xxx_z_0_xy[i] * a_exp + 4.0 * g_xxx_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_xz[i] = 4.0 * g_x_z_0_xz[i] - 4.0 * g_x_zzz_0_xz[i] * b_exp - 4.0 * g_xxx_z_0_xz[i] * a_exp + 4.0 * g_xxx_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_yy[i] = 4.0 * g_x_z_0_yy[i] - 4.0 * g_x_zzz_0_yy[i] * b_exp - 4.0 * g_xxx_z_0_yy[i] * a_exp + 4.0 * g_xxx_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_yz[i] = 4.0 * g_x_z_0_yz[i] - 4.0 * g_x_zzz_0_yz[i] * b_exp - 4.0 * g_xxx_z_0_yz[i] * a_exp + 4.0 * g_xxx_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_zz[i] = 4.0 * g_x_z_0_zz[i] - 4.0 * g_x_zzz_0_zz[i] * b_exp - 4.0 * g_xxx_z_0_zz[i] * a_exp + 4.0 * g_xxx_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_0_xx, g_x_z_0_0_xy_xx_0_xy, g_x_z_0_0_xy_xx_0_xz, g_x_z_0_0_xy_xx_0_yy, g_x_z_0_0_xy_xx_0_yz, g_x_z_0_0_xy_xx_0_zz, g_xxy_xxz_0_xx, g_xxy_xxz_0_xy, g_xxy_xxz_0_xz, g_xxy_xxz_0_yy, g_xxy_xxz_0_yz, g_xxy_xxz_0_zz, g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_0_xx[i] = -2.0 * g_y_xxz_0_xx[i] * b_exp + 4.0 * g_xxy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_xy[i] = -2.0 * g_y_xxz_0_xy[i] * b_exp + 4.0 * g_xxy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_xz[i] = -2.0 * g_y_xxz_0_xz[i] * b_exp + 4.0 * g_xxy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_yy[i] = -2.0 * g_y_xxz_0_yy[i] * b_exp + 4.0 * g_xxy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_yz[i] = -2.0 * g_y_xxz_0_yz[i] * b_exp + 4.0 * g_xxy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_zz[i] = -2.0 * g_y_xxz_0_zz[i] * b_exp + 4.0 * g_xxy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_z_0_0_xy_xy_0_xx, g_x_z_0_0_xy_xy_0_xy, g_x_z_0_0_xy_xy_0_xz, g_x_z_0_0_xy_xy_0_yy, g_x_z_0_0_xy_xy_0_yz, g_x_z_0_0_xy_xy_0_zz, g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xy_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_z_0_0_xy_xz_0_xx, g_x_z_0_0_xy_xz_0_xy, g_x_z_0_0_xy_xz_0_xz, g_x_z_0_0_xy_xz_0_yy, g_x_z_0_0_xy_xz_0_yz, g_x_z_0_0_xy_xz_0_zz, g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xzz_0_xx, g_xxy_xzz_0_xy, g_xxy_xzz_0_xz, g_xxy_xzz_0_yy, g_xxy_xzz_0_yz, g_xxy_xzz_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xz_0_xx[i] = g_y_x_0_xx[i] - 2.0 * g_y_xzz_0_xx[i] * b_exp - 2.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_xy[i] = g_y_x_0_xy[i] - 2.0 * g_y_xzz_0_xy[i] * b_exp - 2.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_xz[i] = g_y_x_0_xz[i] - 2.0 * g_y_xzz_0_xz[i] * b_exp - 2.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_yy[i] = g_y_x_0_yy[i] - 2.0 * g_y_xzz_0_yy[i] * b_exp - 2.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_yz[i] = g_y_x_0_yz[i] - 2.0 * g_y_xzz_0_yz[i] * b_exp - 2.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_zz[i] = g_y_x_0_zz[i] - 2.0 * g_y_xzz_0_zz[i] * b_exp - 2.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_z_0_0_xy_yy_0_xx, g_x_z_0_0_xy_yy_0_xy, g_x_z_0_0_xy_yy_0_xz, g_x_z_0_0_xy_yy_0_yy, g_x_z_0_0_xy_yy_0_yz, g_x_z_0_0_xy_yy_0_zz, g_xxy_yyz_0_xx, g_xxy_yyz_0_xy, g_xxy_yyz_0_xz, g_xxy_yyz_0_yy, g_xxy_yyz_0_yz, g_xxy_yyz_0_zz, g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yy_0_xx[i] = -2.0 * g_y_yyz_0_xx[i] * b_exp + 4.0 * g_xxy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_xy[i] = -2.0 * g_y_yyz_0_xy[i] * b_exp + 4.0 * g_xxy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_xz[i] = -2.0 * g_y_yyz_0_xz[i] * b_exp + 4.0 * g_xxy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_yy[i] = -2.0 * g_y_yyz_0_yy[i] * b_exp + 4.0 * g_xxy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_yz[i] = -2.0 * g_y_yyz_0_yz[i] * b_exp + 4.0 * g_xxy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_zz[i] = -2.0 * g_y_yyz_0_zz[i] * b_exp + 4.0 * g_xxy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_z_0_0_xy_yz_0_xx, g_x_z_0_0_xy_yz_0_xy, g_x_z_0_0_xy_yz_0_xz, g_x_z_0_0_xy_yz_0_yy, g_x_z_0_0_xy_yz_0_yz, g_x_z_0_0_xy_yz_0_zz, g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_xxy_yzz_0_xx, g_xxy_yzz_0_xy, g_xxy_yzz_0_xz, g_xxy_yzz_0_yy, g_xxy_yzz_0_yz, g_xxy_yzz_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yz_0_xx[i] = g_y_y_0_xx[i] - 2.0 * g_y_yzz_0_xx[i] * b_exp - 2.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_xy[i] = g_y_y_0_xy[i] - 2.0 * g_y_yzz_0_xy[i] * b_exp - 2.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_xz[i] = g_y_y_0_xz[i] - 2.0 * g_y_yzz_0_xz[i] * b_exp - 2.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_yy[i] = g_y_y_0_yy[i] - 2.0 * g_y_yzz_0_yy[i] * b_exp - 2.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_yz[i] = g_y_y_0_yz[i] - 2.0 * g_y_yzz_0_yz[i] * b_exp - 2.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_zz[i] = g_y_y_0_zz[i] - 2.0 * g_y_yzz_0_zz[i] * b_exp - 2.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_z_0_0_xy_zz_0_xx, g_x_z_0_0_xy_zz_0_xy, g_x_z_0_0_xy_zz_0_xz, g_x_z_0_0_xy_zz_0_yy, g_x_z_0_0_xy_zz_0_yz, g_x_z_0_0_xy_zz_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_xxy_zzz_0_xx, g_xxy_zzz_0_xy, g_xxy_zzz_0_xz, g_xxy_zzz_0_yy, g_xxy_zzz_0_yz, g_xxy_zzz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_zzz_0_xx, g_y_zzz_0_xy, g_y_zzz_0_xz, g_y_zzz_0_yy, g_y_zzz_0_yz, g_y_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_zz_0_xx[i] = 2.0 * g_y_z_0_xx[i] - 2.0 * g_y_zzz_0_xx[i] * b_exp - 4.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_xy[i] = 2.0 * g_y_z_0_xy[i] - 2.0 * g_y_zzz_0_xy[i] * b_exp - 4.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_xz[i] = 2.0 * g_y_z_0_xz[i] - 2.0 * g_y_zzz_0_xz[i] * b_exp - 4.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_yy[i] = 2.0 * g_y_z_0_yy[i] - 2.0 * g_y_zzz_0_yy[i] * b_exp - 4.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_yz[i] = 2.0 * g_y_z_0_yz[i] - 2.0 * g_y_zzz_0_yz[i] * b_exp - 4.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_zz[i] = 2.0 * g_y_z_0_zz[i] - 2.0 * g_y_zzz_0_zz[i] * b_exp - 4.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_0_xx, g_x_z_0_0_xz_xx_0_xy, g_x_z_0_0_xz_xx_0_xz, g_x_z_0_0_xz_xx_0_yy, g_x_z_0_0_xz_xx_0_yz, g_x_z_0_0_xz_xx_0_zz, g_xxz_xxz_0_xx, g_xxz_xxz_0_xy, g_xxz_xxz_0_xz, g_xxz_xxz_0_yy, g_xxz_xxz_0_yz, g_xxz_xxz_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_0_xx[i] = -2.0 * g_z_xxz_0_xx[i] * b_exp + 4.0 * g_xxz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_xy[i] = -2.0 * g_z_xxz_0_xy[i] * b_exp + 4.0 * g_xxz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_xz[i] = -2.0 * g_z_xxz_0_xz[i] * b_exp + 4.0 * g_xxz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_yy[i] = -2.0 * g_z_xxz_0_yy[i] * b_exp + 4.0 * g_xxz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_yz[i] = -2.0 * g_z_xxz_0_yz[i] * b_exp + 4.0 * g_xxz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_zz[i] = -2.0 * g_z_xxz_0_zz[i] * b_exp + 4.0 * g_xxz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_z_0_0_xz_xy_0_xx, g_x_z_0_0_xz_xy_0_xy, g_x_z_0_0_xz_xy_0_xz, g_x_z_0_0_xz_xy_0_yy, g_x_z_0_0_xz_xy_0_yz, g_x_z_0_0_xz_xy_0_zz, g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xy_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_z_0_0_xz_xz_0_xx, g_x_z_0_0_xz_xz_0_xy, g_x_z_0_0_xz_xz_0_xz, g_x_z_0_0_xz_xz_0_yy, g_x_z_0_0_xz_xz_0_yz, g_x_z_0_0_xz_xz_0_zz, g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xzz_0_xx, g_xxz_xzz_0_xy, g_xxz_xzz_0_xz, g_xxz_xzz_0_yy, g_xxz_xzz_0_yz, g_xxz_xzz_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xz_0_xx[i] = g_z_x_0_xx[i] - 2.0 * g_z_xzz_0_xx[i] * b_exp - 2.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_xy[i] = g_z_x_0_xy[i] - 2.0 * g_z_xzz_0_xy[i] * b_exp - 2.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_xz[i] = g_z_x_0_xz[i] - 2.0 * g_z_xzz_0_xz[i] * b_exp - 2.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_yy[i] = g_z_x_0_yy[i] - 2.0 * g_z_xzz_0_yy[i] * b_exp - 2.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_yz[i] = g_z_x_0_yz[i] - 2.0 * g_z_xzz_0_yz[i] * b_exp - 2.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_zz[i] = g_z_x_0_zz[i] - 2.0 * g_z_xzz_0_zz[i] * b_exp - 2.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_z_0_0_xz_yy_0_xx, g_x_z_0_0_xz_yy_0_xy, g_x_z_0_0_xz_yy_0_xz, g_x_z_0_0_xz_yy_0_yy, g_x_z_0_0_xz_yy_0_yz, g_x_z_0_0_xz_yy_0_zz, g_xxz_yyz_0_xx, g_xxz_yyz_0_xy, g_xxz_yyz_0_xz, g_xxz_yyz_0_yy, g_xxz_yyz_0_yz, g_xxz_yyz_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yy_0_xx[i] = -2.0 * g_z_yyz_0_xx[i] * b_exp + 4.0 * g_xxz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_xy[i] = -2.0 * g_z_yyz_0_xy[i] * b_exp + 4.0 * g_xxz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_xz[i] = -2.0 * g_z_yyz_0_xz[i] * b_exp + 4.0 * g_xxz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_yy[i] = -2.0 * g_z_yyz_0_yy[i] * b_exp + 4.0 * g_xxz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_yz[i] = -2.0 * g_z_yyz_0_yz[i] * b_exp + 4.0 * g_xxz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_zz[i] = -2.0 * g_z_yyz_0_zz[i] * b_exp + 4.0 * g_xxz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_z_0_0_xz_yz_0_xx, g_x_z_0_0_xz_yz_0_xy, g_x_z_0_0_xz_yz_0_xz, g_x_z_0_0_xz_yz_0_yy, g_x_z_0_0_xz_yz_0_yz, g_x_z_0_0_xz_yz_0_zz, g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_xxz_yzz_0_xx, g_xxz_yzz_0_xy, g_xxz_yzz_0_xz, g_xxz_yzz_0_yy, g_xxz_yzz_0_yz, g_xxz_yzz_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yz_0_xx[i] = g_z_y_0_xx[i] - 2.0 * g_z_yzz_0_xx[i] * b_exp - 2.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_xy[i] = g_z_y_0_xy[i] - 2.0 * g_z_yzz_0_xy[i] * b_exp - 2.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_xz[i] = g_z_y_0_xz[i] - 2.0 * g_z_yzz_0_xz[i] * b_exp - 2.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_yy[i] = g_z_y_0_yy[i] - 2.0 * g_z_yzz_0_yy[i] * b_exp - 2.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_yz[i] = g_z_y_0_yz[i] - 2.0 * g_z_yzz_0_yz[i] * b_exp - 2.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_zz[i] = g_z_y_0_zz[i] - 2.0 * g_z_yzz_0_zz[i] * b_exp - 2.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_z_0_0_xz_zz_0_xx, g_x_z_0_0_xz_zz_0_xy, g_x_z_0_0_xz_zz_0_xz, g_x_z_0_0_xz_zz_0_yy, g_x_z_0_0_xz_zz_0_yz, g_x_z_0_0_xz_zz_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_xxz_zzz_0_xx, g_xxz_zzz_0_xy, g_xxz_zzz_0_xz, g_xxz_zzz_0_yy, g_xxz_zzz_0_yz, g_xxz_zzz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_zzz_0_xx, g_z_zzz_0_xy, g_z_zzz_0_xz, g_z_zzz_0_yy, g_z_zzz_0_yz, g_z_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_zz_0_xx[i] = 2.0 * g_z_z_0_xx[i] - 2.0 * g_z_zzz_0_xx[i] * b_exp - 4.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_xy[i] = 2.0 * g_z_z_0_xy[i] - 2.0 * g_z_zzz_0_xy[i] * b_exp - 4.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_xz[i] = 2.0 * g_z_z_0_xz[i] - 2.0 * g_z_zzz_0_xz[i] * b_exp - 4.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_yy[i] = 2.0 * g_z_z_0_yy[i] - 2.0 * g_z_zzz_0_yy[i] * b_exp - 4.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_yz[i] = 2.0 * g_z_z_0_yz[i] - 2.0 * g_z_zzz_0_yz[i] * b_exp - 4.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_zz[i] = 2.0 * g_z_z_0_zz[i] - 2.0 * g_z_zzz_0_zz[i] * b_exp - 4.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_0_xx, g_x_z_0_0_yy_xx_0_xy, g_x_z_0_0_yy_xx_0_xz, g_x_z_0_0_yy_xx_0_yy, g_x_z_0_0_yy_xx_0_yz, g_x_z_0_0_yy_xx_0_zz, g_xyy_xxz_0_xx, g_xyy_xxz_0_xy, g_xyy_xxz_0_xz, g_xyy_xxz_0_yy, g_xyy_xxz_0_yz, g_xyy_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_0_xx[i] = 4.0 * g_xyy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_xy[i] = 4.0 * g_xyy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_xz[i] = 4.0 * g_xyy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_yy[i] = 4.0 * g_xyy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_yz[i] = 4.0 * g_xyy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_zz[i] = 4.0 * g_xyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_z_0_0_yy_xy_0_xx, g_x_z_0_0_yy_xy_0_xy, g_x_z_0_0_yy_xy_0_xz, g_x_z_0_0_yy_xy_0_yy, g_x_z_0_0_yy_xy_0_yz, g_x_z_0_0_yy_xy_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xy_0_xx[i] = 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_xy[i] = 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_xz[i] = 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_yy[i] = 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_yz[i] = 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_zz[i] = 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_z_0_0_yy_xz_0_xx, g_x_z_0_0_yy_xz_0_xy, g_x_z_0_0_yy_xz_0_xz, g_x_z_0_0_yy_xz_0_yy, g_x_z_0_0_yy_xz_0_yz, g_x_z_0_0_yy_xz_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xzz_0_xx, g_xyy_xzz_0_xy, g_xyy_xzz_0_xz, g_xyy_xzz_0_yy, g_xyy_xzz_0_yz, g_xyy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xz_0_xx[i] = -2.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_xy[i] = -2.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_xz[i] = -2.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_yy[i] = -2.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_yz[i] = -2.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_zz[i] = -2.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_z_0_0_yy_yy_0_xx, g_x_z_0_0_yy_yy_0_xy, g_x_z_0_0_yy_yy_0_xz, g_x_z_0_0_yy_yy_0_yy, g_x_z_0_0_yy_yy_0_yz, g_x_z_0_0_yy_yy_0_zz, g_xyy_yyz_0_xx, g_xyy_yyz_0_xy, g_xyy_yyz_0_xz, g_xyy_yyz_0_yy, g_xyy_yyz_0_yz, g_xyy_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yy_0_xx[i] = 4.0 * g_xyy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_xy[i] = 4.0 * g_xyy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_xz[i] = 4.0 * g_xyy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_yy[i] = 4.0 * g_xyy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_yz[i] = 4.0 * g_xyy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_zz[i] = 4.0 * g_xyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_z_0_0_yy_yz_0_xx, g_x_z_0_0_yy_yz_0_xy, g_x_z_0_0_yy_yz_0_xz, g_x_z_0_0_yy_yz_0_yy, g_x_z_0_0_yy_yz_0_yz, g_x_z_0_0_yy_yz_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_xyy_yzz_0_xx, g_xyy_yzz_0_xy, g_xyy_yzz_0_xz, g_xyy_yzz_0_yy, g_xyy_yzz_0_yz, g_xyy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yz_0_xx[i] = -2.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_xy[i] = -2.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_xz[i] = -2.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_yy[i] = -2.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_yz[i] = -2.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_zz[i] = -2.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_z_0_0_yy_zz_0_xx, g_x_z_0_0_yy_zz_0_xy, g_x_z_0_0_yy_zz_0_xz, g_x_z_0_0_yy_zz_0_yy, g_x_z_0_0_yy_zz_0_yz, g_x_z_0_0_yy_zz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz, g_xyy_zzz_0_xx, g_xyy_zzz_0_xy, g_xyy_zzz_0_xz, g_xyy_zzz_0_yy, g_xyy_zzz_0_yz, g_xyy_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_zz_0_xx[i] = -4.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_xy[i] = -4.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_xz[i] = -4.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_yy[i] = -4.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_yz[i] = -4.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_zz[i] = -4.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_0_xx, g_x_z_0_0_yz_xx_0_xy, g_x_z_0_0_yz_xx_0_xz, g_x_z_0_0_yz_xx_0_yy, g_x_z_0_0_yz_xx_0_yz, g_x_z_0_0_yz_xx_0_zz, g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_0_xx[i] = 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_xy[i] = 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_xz[i] = 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_yy[i] = 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_yz[i] = 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_zz[i] = 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_z_0_0_yz_xy_0_xx, g_x_z_0_0_yz_xy_0_xy, g_x_z_0_0_yz_xy_0_xz, g_x_z_0_0_yz_xy_0_yy, g_x_z_0_0_yz_xy_0_yz, g_x_z_0_0_yz_xy_0_zz, g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xy_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_z_0_0_yz_xz_0_xx, g_x_z_0_0_yz_xz_0_xy, g_x_z_0_0_yz_xz_0_xz, g_x_z_0_0_yz_xz_0_yy, g_x_z_0_0_yz_xz_0_yz, g_x_z_0_0_yz_xz_0_zz, g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xz_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_z_0_0_yz_yy_0_xx, g_x_z_0_0_yz_yy_0_xy, g_x_z_0_0_yz_yy_0_xz, g_x_z_0_0_yz_yy_0_yy, g_x_z_0_0_yz_yy_0_yz, g_x_z_0_0_yz_yy_0_zz, g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yy_0_xx[i] = 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_xy[i] = 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_xz[i] = 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_yy[i] = 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_yz[i] = 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_zz[i] = 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_z_0_0_yz_yz_0_xx, g_x_z_0_0_yz_yz_0_xy, g_x_z_0_0_yz_yz_0_xz, g_x_z_0_0_yz_yz_0_yy, g_x_z_0_0_yz_yz_0_yz, g_x_z_0_0_yz_yz_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yz_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_z_0_0_yz_zz_0_xx, g_x_z_0_0_yz_zz_0_xy, g_x_z_0_0_yz_zz_0_xz, g_x_z_0_0_yz_zz_0_yy, g_x_z_0_0_yz_zz_0_yz, g_x_z_0_0_yz_zz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_xyz_zzz_0_xx, g_xyz_zzz_0_xy, g_xyz_zzz_0_xz, g_xyz_zzz_0_yy, g_xyz_zzz_0_yz, g_xyz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_zz_0_xx[i] = -4.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_xy[i] = -4.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_xz[i] = -4.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_yy[i] = -4.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_yz[i] = -4.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_zz[i] = -4.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_0_xx, g_x_z_0_0_zz_xx_0_xy, g_x_z_0_0_zz_xx_0_xz, g_x_z_0_0_zz_xx_0_yy, g_x_z_0_0_zz_xx_0_yz, g_x_z_0_0_zz_xx_0_zz, g_xzz_xxz_0_xx, g_xzz_xxz_0_xy, g_xzz_xxz_0_xz, g_xzz_xxz_0_yy, g_xzz_xxz_0_yz, g_xzz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_0_xx[i] = 4.0 * g_xzz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_xy[i] = 4.0 * g_xzz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_xz[i] = 4.0 * g_xzz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_yy[i] = 4.0 * g_xzz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_yz[i] = 4.0 * g_xzz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_zz[i] = 4.0 * g_xzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_z_0_0_zz_xy_0_xx, g_x_z_0_0_zz_xy_0_xy, g_x_z_0_0_zz_xy_0_xz, g_x_z_0_0_zz_xy_0_yy, g_x_z_0_0_zz_xy_0_yz, g_x_z_0_0_zz_xy_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xy_0_xx[i] = 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_xy[i] = 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_xz[i] = 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_yy[i] = 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_yz[i] = 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_zz[i] = 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_z_0_0_zz_xz_0_xx, g_x_z_0_0_zz_xz_0_xy, g_x_z_0_0_zz_xz_0_xz, g_x_z_0_0_zz_xz_0_yy, g_x_z_0_0_zz_xz_0_yz, g_x_z_0_0_zz_xz_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xzz_0_xx, g_xzz_xzz_0_xy, g_xzz_xzz_0_xz, g_xzz_xzz_0_yy, g_xzz_xzz_0_yz, g_xzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xz_0_xx[i] = -2.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_xy[i] = -2.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_xz[i] = -2.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_yy[i] = -2.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_yz[i] = -2.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_zz[i] = -2.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_z_0_0_zz_yy_0_xx, g_x_z_0_0_zz_yy_0_xy, g_x_z_0_0_zz_yy_0_xz, g_x_z_0_0_zz_yy_0_yy, g_x_z_0_0_zz_yy_0_yz, g_x_z_0_0_zz_yy_0_zz, g_xzz_yyz_0_xx, g_xzz_yyz_0_xy, g_xzz_yyz_0_xz, g_xzz_yyz_0_yy, g_xzz_yyz_0_yz, g_xzz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yy_0_xx[i] = 4.0 * g_xzz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_xy[i] = 4.0 * g_xzz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_xz[i] = 4.0 * g_xzz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_yy[i] = 4.0 * g_xzz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_yz[i] = 4.0 * g_xzz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_zz[i] = 4.0 * g_xzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_z_0_0_zz_yz_0_xx, g_x_z_0_0_zz_yz_0_xy, g_x_z_0_0_zz_yz_0_xz, g_x_z_0_0_zz_yz_0_yy, g_x_z_0_0_zz_yz_0_yz, g_x_z_0_0_zz_yz_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_xzz_yzz_0_xx, g_xzz_yzz_0_xy, g_xzz_yzz_0_xz, g_xzz_yzz_0_yy, g_xzz_yzz_0_yz, g_xzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yz_0_xx[i] = -2.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_xy[i] = -2.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_xz[i] = -2.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_yy[i] = -2.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_yz[i] = -2.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_zz[i] = -2.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_z_0_0_zz_zz_0_xx, g_x_z_0_0_zz_zz_0_xy, g_x_z_0_0_zz_zz_0_xz, g_x_z_0_0_zz_zz_0_yy, g_x_z_0_0_zz_zz_0_yz, g_x_z_0_0_zz_zz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz, g_xzz_zzz_0_xx, g_xzz_zzz_0_xy, g_xzz_zzz_0_xz, g_xzz_zzz_0_yy, g_xzz_zzz_0_yz, g_xzz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_zz_0_xx[i] = -4.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_xy[i] = -4.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_xz[i] = -4.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_yy[i] = -4.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_yz[i] = -4.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_zz[i] = -4.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xxx_0_xx, g_xxy_xxx_0_xy, g_xxy_xxx_0_xz, g_xxy_xxx_0_yy, g_xxy_xxx_0_yz, g_xxy_xxx_0_zz, g_y_x_0_0_xx_xx_0_xx, g_y_x_0_0_xx_xx_0_xy, g_y_x_0_0_xx_xx_0_xz, g_y_x_0_0_xx_xx_0_yy, g_y_x_0_0_xx_xx_0_yz, g_y_x_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_0_xx[i] = -4.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_xy[i] = -4.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_xz[i] = -4.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_yy[i] = -4.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_yz[i] = -4.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_zz[i] = -4.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xxy_xxy_0_xx, g_xxy_xxy_0_xy, g_xxy_xxy_0_xz, g_xxy_xxy_0_yy, g_xxy_xxy_0_yz, g_xxy_xxy_0_zz, g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_y_x_0_0_xx_xy_0_xx, g_y_x_0_0_xx_xy_0_xy, g_y_x_0_0_xx_xy_0_xz, g_y_x_0_0_xx_xy_0_yy, g_y_x_0_0_xx_xy_0_yz, g_y_x_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xy_0_xx[i] = -2.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_xy[i] = -2.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_xz[i] = -2.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_yy[i] = -2.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_yz[i] = -2.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_zz[i] = -2.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xxy_xxz_0_xx, g_xxy_xxz_0_xy, g_xxy_xxz_0_xz, g_xxy_xxz_0_yy, g_xxy_xxz_0_yz, g_xxy_xxz_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_y_x_0_0_xx_xz_0_xx, g_y_x_0_0_xx_xz_0_xy, g_y_x_0_0_xx_xz_0_xz, g_y_x_0_0_xx_xz_0_yy, g_y_x_0_0_xx_xz_0_yz, g_y_x_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xz_0_xx[i] = -2.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_xy[i] = -2.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_xz[i] = -2.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_yy[i] = -2.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_yz[i] = -2.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_zz[i] = -2.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xxy_xyy_0_xx, g_xxy_xyy_0_xy, g_xxy_xyy_0_xz, g_xxy_xyy_0_yy, g_xxy_xyy_0_yz, g_xxy_xyy_0_zz, g_y_x_0_0_xx_yy_0_xx, g_y_x_0_0_xx_yy_0_xy, g_y_x_0_0_xx_yy_0_xz, g_y_x_0_0_xx_yy_0_yy, g_y_x_0_0_xx_yy_0_yz, g_y_x_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yy_0_xx[i] = 4.0 * g_xxy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_xy[i] = 4.0 * g_xxy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_xz[i] = 4.0 * g_xxy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_yy[i] = 4.0 * g_xxy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_yz[i] = 4.0 * g_xxy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_zz[i] = 4.0 * g_xxy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_x_0_0_xx_yz_0_xx, g_y_x_0_0_xx_yz_0_xy, g_y_x_0_0_xx_yz_0_xz, g_y_x_0_0_xx_yz_0_yy, g_y_x_0_0_xx_yz_0_yz, g_y_x_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yz_0_xx[i] = 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_xy[i] = 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_xz[i] = 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_yy[i] = 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_yz[i] = 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_zz[i] = 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xxy_xzz_0_xx, g_xxy_xzz_0_xy, g_xxy_xzz_0_xz, g_xxy_xzz_0_yy, g_xxy_xzz_0_yz, g_xxy_xzz_0_zz, g_y_x_0_0_xx_zz_0_xx, g_y_x_0_0_xx_zz_0_xy, g_y_x_0_0_xx_zz_0_xz, g_y_x_0_0_xx_zz_0_yy, g_y_x_0_0_xx_zz_0_yz, g_y_x_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_zz_0_xx[i] = 4.0 * g_xxy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_xy[i] = 4.0 * g_xxy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_xz[i] = 4.0 * g_xxy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_yy[i] = 4.0 * g_xxy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_yz[i] = 4.0 * g_xxy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_zz[i] = 4.0 * g_xxy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xxx_0_xx, g_x_xxx_0_xy, g_x_xxx_0_xz, g_x_xxx_0_yy, g_x_xxx_0_yz, g_x_xxx_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xxx_0_xx, g_xyy_xxx_0_xy, g_xyy_xxx_0_xz, g_xyy_xxx_0_yy, g_xyy_xxx_0_yz, g_xyy_xxx_0_zz, g_y_x_0_0_xy_xx_0_xx, g_y_x_0_0_xy_xx_0_xy, g_y_x_0_0_xy_xx_0_xz, g_y_x_0_0_xy_xx_0_yy, g_y_x_0_0_xy_xx_0_yz, g_y_x_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_0_xx[i] = 2.0 * g_x_x_0_xx[i] - 2.0 * g_x_xxx_0_xx[i] * b_exp - 4.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_xy[i] = 2.0 * g_x_x_0_xy[i] - 2.0 * g_x_xxx_0_xy[i] * b_exp - 4.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_xz[i] = 2.0 * g_x_x_0_xz[i] - 2.0 * g_x_xxx_0_xz[i] * b_exp - 4.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_yy[i] = 2.0 * g_x_x_0_yy[i] - 2.0 * g_x_xxx_0_yy[i] * b_exp - 4.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_yz[i] = 2.0 * g_x_x_0_yz[i] - 2.0 * g_x_xxx_0_yz[i] * b_exp - 4.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_zz[i] = 2.0 * g_x_x_0_zz[i] - 2.0 * g_x_xxx_0_zz[i] * b_exp - 4.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xyy_xxy_0_xx, g_xyy_xxy_0_xy, g_xyy_xxy_0_xz, g_xyy_xxy_0_yy, g_xyy_xxy_0_yz, g_xyy_xxy_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_y_x_0_0_xy_xy_0_xx, g_y_x_0_0_xy_xy_0_xy, g_y_x_0_0_xy_xy_0_xz, g_y_x_0_0_xy_xy_0_yy, g_y_x_0_0_xy_xy_0_yz, g_y_x_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xy_0_xx[i] = g_x_y_0_xx[i] - 2.0 * g_x_xxy_0_xx[i] * b_exp - 2.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_xy[i] = g_x_y_0_xy[i] - 2.0 * g_x_xxy_0_xy[i] * b_exp - 2.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_xz[i] = g_x_y_0_xz[i] - 2.0 * g_x_xxy_0_xz[i] * b_exp - 2.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_yy[i] = g_x_y_0_yy[i] - 2.0 * g_x_xxy_0_yy[i] * b_exp - 2.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_yz[i] = g_x_y_0_yz[i] - 2.0 * g_x_xxy_0_yz[i] * b_exp - 2.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_zz[i] = g_x_y_0_zz[i] - 2.0 * g_x_xxy_0_zz[i] * b_exp - 2.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xyy_xxz_0_xx, g_xyy_xxz_0_xy, g_xyy_xxz_0_xz, g_xyy_xxz_0_yy, g_xyy_xxz_0_yz, g_xyy_xxz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz, g_y_x_0_0_xy_xz_0_xx, g_y_x_0_0_xy_xz_0_xy, g_y_x_0_0_xy_xz_0_xz, g_y_x_0_0_xy_xz_0_yy, g_y_x_0_0_xy_xz_0_yz, g_y_x_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xz_0_xx[i] = g_x_z_0_xx[i] - 2.0 * g_x_xxz_0_xx[i] * b_exp - 2.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_xy[i] = g_x_z_0_xy[i] - 2.0 * g_x_xxz_0_xy[i] * b_exp - 2.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_xz[i] = g_x_z_0_xz[i] - 2.0 * g_x_xxz_0_xz[i] * b_exp - 2.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_yy[i] = g_x_z_0_yy[i] - 2.0 * g_x_xxz_0_yy[i] * b_exp - 2.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_yz[i] = g_x_z_0_yz[i] - 2.0 * g_x_xxz_0_yz[i] * b_exp - 2.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_zz[i] = g_x_z_0_zz[i] - 2.0 * g_x_xxz_0_zz[i] * b_exp - 2.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_xyy_xyy_0_xx, g_xyy_xyy_0_xy, g_xyy_xyy_0_xz, g_xyy_xyy_0_yy, g_xyy_xyy_0_yz, g_xyy_xyy_0_zz, g_y_x_0_0_xy_yy_0_xx, g_y_x_0_0_xy_yy_0_xy, g_y_x_0_0_xy_yy_0_xz, g_y_x_0_0_xy_yy_0_yy, g_y_x_0_0_xy_yy_0_yz, g_y_x_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yy_0_xx[i] = -2.0 * g_x_xyy_0_xx[i] * b_exp + 4.0 * g_xyy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_xy[i] = -2.0 * g_x_xyy_0_xy[i] * b_exp + 4.0 * g_xyy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_xz[i] = -2.0 * g_x_xyy_0_xz[i] * b_exp + 4.0 * g_xyy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_yy[i] = -2.0 * g_x_xyy_0_yy[i] * b_exp + 4.0 * g_xyy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_yz[i] = -2.0 * g_x_xyy_0_yz[i] * b_exp + 4.0 * g_xyy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_zz[i] = -2.0 * g_x_xyy_0_zz[i] * b_exp + 4.0 * g_xyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz, g_y_x_0_0_xy_yz_0_xx, g_y_x_0_0_xy_yz_0_xy, g_y_x_0_0_xy_yz_0_xz, g_y_x_0_0_xy_yz_0_yy, g_y_x_0_0_xy_yz_0_yz, g_y_x_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yz_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_xyy_xzz_0_xx, g_xyy_xzz_0_xy, g_xyy_xzz_0_xz, g_xyy_xzz_0_yy, g_xyy_xzz_0_yz, g_xyy_xzz_0_zz, g_y_x_0_0_xy_zz_0_xx, g_y_x_0_0_xy_zz_0_xy, g_y_x_0_0_xy_zz_0_xz, g_y_x_0_0_xy_zz_0_yy, g_y_x_0_0_xy_zz_0_yz, g_y_x_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_zz_0_xx[i] = -2.0 * g_x_xzz_0_xx[i] * b_exp + 4.0 * g_xyy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_xy[i] = -2.0 * g_x_xzz_0_xy[i] * b_exp + 4.0 * g_xyy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_xz[i] = -2.0 * g_x_xzz_0_xz[i] * b_exp + 4.0 * g_xyy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_yy[i] = -2.0 * g_x_xzz_0_yy[i] * b_exp + 4.0 * g_xyy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_yz[i] = -2.0 * g_x_xzz_0_yz[i] * b_exp + 4.0 * g_xyy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_zz[i] = -2.0 * g_x_xzz_0_zz[i] * b_exp + 4.0 * g_xyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xxx_0_xx, g_xyz_xxx_0_xy, g_xyz_xxx_0_xz, g_xyz_xxx_0_yy, g_xyz_xxx_0_yz, g_xyz_xxx_0_zz, g_y_x_0_0_xz_xx_0_xx, g_y_x_0_0_xz_xx_0_xy, g_y_x_0_0_xz_xx_0_xz, g_y_x_0_0_xz_xx_0_yy, g_y_x_0_0_xz_xx_0_yz, g_y_x_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_0_xx[i] = -4.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_xy[i] = -4.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_xz[i] = -4.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_yy[i] = -4.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_yz[i] = -4.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_zz[i] = -4.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_y_x_0_0_xz_xy_0_xx, g_y_x_0_0_xz_xy_0_xy, g_y_x_0_0_xz_xy_0_xz, g_y_x_0_0_xz_xy_0_yy, g_y_x_0_0_xz_xy_0_yz, g_y_x_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xy_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_y_x_0_0_xz_xz_0_xx, g_y_x_0_0_xz_xz_0_xy, g_y_x_0_0_xz_xz_0_xz, g_y_x_0_0_xz_xz_0_yy, g_y_x_0_0_xz_xz_0_yz, g_y_x_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz, g_y_x_0_0_xz_yy_0_xx, g_y_x_0_0_xz_yy_0_xy, g_y_x_0_0_xz_yy_0_xz, g_y_x_0_0_xz_yy_0_yy, g_y_x_0_0_xz_yy_0_yz, g_y_x_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yy_0_xx[i] = 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_xy[i] = 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_xz[i] = 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_yy[i] = 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_yz[i] = 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_zz[i] = 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_y_x_0_0_xz_yz_0_xx, g_y_x_0_0_xz_yz_0_xy, g_y_x_0_0_xz_yz_0_xz, g_y_x_0_0_xz_yz_0_yy, g_y_x_0_0_xz_yz_0_yz, g_y_x_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz, g_y_x_0_0_xz_zz_0_xx, g_y_x_0_0_xz_zz_0_xy, g_y_x_0_0_xz_zz_0_xz, g_y_x_0_0_xz_zz_0_yy, g_y_x_0_0_xz_zz_0_yz, g_y_x_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_zz_0_xx[i] = 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_xy[i] = 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_xz[i] = 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_yy[i] = 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_yz[i] = 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_zz[i] = 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_x_0_0_yy_xx_0_xx, g_y_x_0_0_yy_xx_0_xy, g_y_x_0_0_yy_xx_0_xz, g_y_x_0_0_yy_xx_0_yy, g_y_x_0_0_yy_xx_0_yz, g_y_x_0_0_yy_xx_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xxx_0_xx, g_y_xxx_0_xy, g_y_xxx_0_xz, g_y_xxx_0_yy, g_y_xxx_0_yz, g_y_xxx_0_zz, g_yyy_x_0_xx, g_yyy_x_0_xy, g_yyy_x_0_xz, g_yyy_x_0_yy, g_yyy_x_0_yz, g_yyy_x_0_zz, g_yyy_xxx_0_xx, g_yyy_xxx_0_xy, g_yyy_xxx_0_xz, g_yyy_xxx_0_yy, g_yyy_xxx_0_yz, g_yyy_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_0_xx[i] = 4.0 * g_y_x_0_xx[i] - 4.0 * g_y_xxx_0_xx[i] * b_exp - 4.0 * g_yyy_x_0_xx[i] * a_exp + 4.0 * g_yyy_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_xy[i] = 4.0 * g_y_x_0_xy[i] - 4.0 * g_y_xxx_0_xy[i] * b_exp - 4.0 * g_yyy_x_0_xy[i] * a_exp + 4.0 * g_yyy_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_xz[i] = 4.0 * g_y_x_0_xz[i] - 4.0 * g_y_xxx_0_xz[i] * b_exp - 4.0 * g_yyy_x_0_xz[i] * a_exp + 4.0 * g_yyy_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_yy[i] = 4.0 * g_y_x_0_yy[i] - 4.0 * g_y_xxx_0_yy[i] * b_exp - 4.0 * g_yyy_x_0_yy[i] * a_exp + 4.0 * g_yyy_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_yz[i] = 4.0 * g_y_x_0_yz[i] - 4.0 * g_y_xxx_0_yz[i] * b_exp - 4.0 * g_yyy_x_0_yz[i] * a_exp + 4.0 * g_yyy_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_zz[i] = 4.0 * g_y_x_0_zz[i] - 4.0 * g_y_xxx_0_zz[i] * b_exp - 4.0 * g_yyy_x_0_zz[i] * a_exp + 4.0 * g_yyy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_x_0_0_yy_xy_0_xx, g_y_x_0_0_yy_xy_0_xy, g_y_x_0_0_yy_xy_0_xz, g_y_x_0_0_yy_xy_0_yy, g_y_x_0_0_yy_xy_0_yz, g_y_x_0_0_yy_xy_0_zz, g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_yyy_xxy_0_xx, g_yyy_xxy_0_xy, g_yyy_xxy_0_xz, g_yyy_xxy_0_yy, g_yyy_xxy_0_yz, g_yyy_xxy_0_zz, g_yyy_y_0_xx, g_yyy_y_0_xy, g_yyy_y_0_xz, g_yyy_y_0_yy, g_yyy_y_0_yz, g_yyy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xy_0_xx[i] = 2.0 * g_y_y_0_xx[i] - 4.0 * g_y_xxy_0_xx[i] * b_exp - 2.0 * g_yyy_y_0_xx[i] * a_exp + 4.0 * g_yyy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_xy[i] = 2.0 * g_y_y_0_xy[i] - 4.0 * g_y_xxy_0_xy[i] * b_exp - 2.0 * g_yyy_y_0_xy[i] * a_exp + 4.0 * g_yyy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_xz[i] = 2.0 * g_y_y_0_xz[i] - 4.0 * g_y_xxy_0_xz[i] * b_exp - 2.0 * g_yyy_y_0_xz[i] * a_exp + 4.0 * g_yyy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_yy[i] = 2.0 * g_y_y_0_yy[i] - 4.0 * g_y_xxy_0_yy[i] * b_exp - 2.0 * g_yyy_y_0_yy[i] * a_exp + 4.0 * g_yyy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_yz[i] = 2.0 * g_y_y_0_yz[i] - 4.0 * g_y_xxy_0_yz[i] * b_exp - 2.0 * g_yyy_y_0_yz[i] * a_exp + 4.0 * g_yyy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_zz[i] = 2.0 * g_y_y_0_zz[i] - 4.0 * g_y_xxy_0_zz[i] * b_exp - 2.0 * g_yyy_y_0_zz[i] * a_exp + 4.0 * g_yyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_x_0_0_yy_xz_0_xx, g_y_x_0_0_yy_xz_0_xy, g_y_x_0_0_yy_xz_0_xz, g_y_x_0_0_yy_xz_0_yy, g_y_x_0_0_yy_xz_0_yz, g_y_x_0_0_yy_xz_0_zz, g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yyy_xxz_0_xx, g_yyy_xxz_0_xy, g_yyy_xxz_0_xz, g_yyy_xxz_0_yy, g_yyy_xxz_0_yz, g_yyy_xxz_0_zz, g_yyy_z_0_xx, g_yyy_z_0_xy, g_yyy_z_0_xz, g_yyy_z_0_yy, g_yyy_z_0_yz, g_yyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xz_0_xx[i] = 2.0 * g_y_z_0_xx[i] - 4.0 * g_y_xxz_0_xx[i] * b_exp - 2.0 * g_yyy_z_0_xx[i] * a_exp + 4.0 * g_yyy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_xy[i] = 2.0 * g_y_z_0_xy[i] - 4.0 * g_y_xxz_0_xy[i] * b_exp - 2.0 * g_yyy_z_0_xy[i] * a_exp + 4.0 * g_yyy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_xz[i] = 2.0 * g_y_z_0_xz[i] - 4.0 * g_y_xxz_0_xz[i] * b_exp - 2.0 * g_yyy_z_0_xz[i] * a_exp + 4.0 * g_yyy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_yy[i] = 2.0 * g_y_z_0_yy[i] - 4.0 * g_y_xxz_0_yy[i] * b_exp - 2.0 * g_yyy_z_0_yy[i] * a_exp + 4.0 * g_yyy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_yz[i] = 2.0 * g_y_z_0_yz[i] - 4.0 * g_y_xxz_0_yz[i] * b_exp - 2.0 * g_yyy_z_0_yz[i] * a_exp + 4.0 * g_yyy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_zz[i] = 2.0 * g_y_z_0_zz[i] - 4.0 * g_y_xxz_0_zz[i] * b_exp - 2.0 * g_yyy_z_0_zz[i] * a_exp + 4.0 * g_yyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_x_0_0_yy_yy_0_xx, g_y_x_0_0_yy_yy_0_xy, g_y_x_0_0_yy_yy_0_xz, g_y_x_0_0_yy_yy_0_yy, g_y_x_0_0_yy_yy_0_yz, g_y_x_0_0_yy_yy_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz, g_yyy_xyy_0_xx, g_yyy_xyy_0_xy, g_yyy_xyy_0_xz, g_yyy_xyy_0_yy, g_yyy_xyy_0_yz, g_yyy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yy_0_xx[i] = -4.0 * g_y_xyy_0_xx[i] * b_exp + 4.0 * g_yyy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_xy[i] = -4.0 * g_y_xyy_0_xy[i] * b_exp + 4.0 * g_yyy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_xz[i] = -4.0 * g_y_xyy_0_xz[i] * b_exp + 4.0 * g_yyy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_yy[i] = -4.0 * g_y_xyy_0_yy[i] * b_exp + 4.0 * g_yyy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_yz[i] = -4.0 * g_y_xyy_0_yz[i] * b_exp + 4.0 * g_yyy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_zz[i] = -4.0 * g_y_xyy_0_zz[i] * b_exp + 4.0 * g_yyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_x_0_0_yy_yz_0_xx, g_y_x_0_0_yy_yz_0_xy, g_y_x_0_0_yy_yz_0_xz, g_y_x_0_0_yy_yz_0_yy, g_y_x_0_0_yy_yz_0_yz, g_y_x_0_0_yy_yz_0_zz, g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_yyy_xyz_0_xx, g_yyy_xyz_0_xy, g_yyy_xyz_0_xz, g_yyy_xyz_0_yy, g_yyy_xyz_0_yz, g_yyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yz_0_xx[i] = -4.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_xy[i] = -4.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_xz[i] = -4.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_yy[i] = -4.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_yz[i] = -4.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_zz[i] = -4.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_x_0_0_yy_zz_0_xx, g_y_x_0_0_yy_zz_0_xy, g_y_x_0_0_yy_zz_0_xz, g_y_x_0_0_yy_zz_0_yy, g_y_x_0_0_yy_zz_0_yz, g_y_x_0_0_yy_zz_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz, g_yyy_xzz_0_xx, g_yyy_xzz_0_xy, g_yyy_xzz_0_xz, g_yyy_xzz_0_yy, g_yyy_xzz_0_yz, g_yyy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_zz_0_xx[i] = -4.0 * g_y_xzz_0_xx[i] * b_exp + 4.0 * g_yyy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_xy[i] = -4.0 * g_y_xzz_0_xy[i] * b_exp + 4.0 * g_yyy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_xz[i] = -4.0 * g_y_xzz_0_xz[i] * b_exp + 4.0 * g_yyy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_yy[i] = -4.0 * g_y_xzz_0_yy[i] * b_exp + 4.0 * g_yyy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_yz[i] = -4.0 * g_y_xzz_0_yz[i] * b_exp + 4.0 * g_yyy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_zz[i] = -4.0 * g_y_xzz_0_zz[i] * b_exp + 4.0 * g_yyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_0_xx, g_y_x_0_0_yz_xx_0_xy, g_y_x_0_0_yz_xx_0_xz, g_y_x_0_0_yz_xx_0_yy, g_y_x_0_0_yz_xx_0_yz, g_y_x_0_0_yz_xx_0_zz, g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xxx_0_xx, g_yyz_xxx_0_xy, g_yyz_xxx_0_xz, g_yyz_xxx_0_yy, g_yyz_xxx_0_yz, g_yyz_xxx_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xxx_0_xx, g_z_xxx_0_xy, g_z_xxx_0_xz, g_z_xxx_0_yy, g_z_xxx_0_yz, g_z_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_0_xx[i] = 2.0 * g_z_x_0_xx[i] - 2.0 * g_z_xxx_0_xx[i] * b_exp - 4.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_xy[i] = 2.0 * g_z_x_0_xy[i] - 2.0 * g_z_xxx_0_xy[i] * b_exp - 4.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_xz[i] = 2.0 * g_z_x_0_xz[i] - 2.0 * g_z_xxx_0_xz[i] * b_exp - 4.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_yy[i] = 2.0 * g_z_x_0_yy[i] - 2.0 * g_z_xxx_0_yy[i] * b_exp - 4.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_yz[i] = 2.0 * g_z_x_0_yz[i] - 2.0 * g_z_xxx_0_yz[i] * b_exp - 4.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_zz[i] = 2.0 * g_z_x_0_zz[i] - 2.0 * g_z_xxx_0_zz[i] * b_exp - 4.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_x_0_0_yz_xy_0_xx, g_y_x_0_0_yz_xy_0_xy, g_y_x_0_0_yz_xy_0_xz, g_y_x_0_0_yz_xy_0_yy, g_y_x_0_0_yz_xy_0_yz, g_y_x_0_0_yz_xy_0_zz, g_yyz_xxy_0_xx, g_yyz_xxy_0_xy, g_yyz_xxy_0_xz, g_yyz_xxy_0_yy, g_yyz_xxy_0_yz, g_yyz_xxy_0_zz, g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xy_0_xx[i] = g_z_y_0_xx[i] - 2.0 * g_z_xxy_0_xx[i] * b_exp - 2.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_xy[i] = g_z_y_0_xy[i] - 2.0 * g_z_xxy_0_xy[i] * b_exp - 2.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_xz[i] = g_z_y_0_xz[i] - 2.0 * g_z_xxy_0_xz[i] * b_exp - 2.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_yy[i] = g_z_y_0_yy[i] - 2.0 * g_z_xxy_0_yy[i] * b_exp - 2.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_yz[i] = g_z_y_0_yz[i] - 2.0 * g_z_xxy_0_yz[i] * b_exp - 2.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_zz[i] = g_z_y_0_zz[i] - 2.0 * g_z_xxy_0_zz[i] * b_exp - 2.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_x_0_0_yz_xz_0_xx, g_y_x_0_0_yz_xz_0_xy, g_y_x_0_0_yz_xz_0_xz, g_y_x_0_0_yz_xz_0_yy, g_y_x_0_0_yz_xz_0_yz, g_y_x_0_0_yz_xz_0_zz, g_yyz_xxz_0_xx, g_yyz_xxz_0_xy, g_yyz_xxz_0_xz, g_yyz_xxz_0_yy, g_yyz_xxz_0_yz, g_yyz_xxz_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xz_0_xx[i] = g_z_z_0_xx[i] - 2.0 * g_z_xxz_0_xx[i] * b_exp - 2.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_xy[i] = g_z_z_0_xy[i] - 2.0 * g_z_xxz_0_xy[i] * b_exp - 2.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_xz[i] = g_z_z_0_xz[i] - 2.0 * g_z_xxz_0_xz[i] * b_exp - 2.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_yy[i] = g_z_z_0_yy[i] - 2.0 * g_z_xxz_0_yy[i] * b_exp - 2.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_yz[i] = g_z_z_0_yz[i] - 2.0 * g_z_xxz_0_yz[i] * b_exp - 2.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_zz[i] = g_z_z_0_zz[i] - 2.0 * g_z_xxz_0_zz[i] * b_exp - 2.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_y_x_0_0_yz_yy_0_xx, g_y_x_0_0_yz_yy_0_xy, g_y_x_0_0_yz_yy_0_xz, g_y_x_0_0_yz_yy_0_yy, g_y_x_0_0_yz_yy_0_yz, g_y_x_0_0_yz_yy_0_zz, g_yyz_xyy_0_xx, g_yyz_xyy_0_xy, g_yyz_xyy_0_xz, g_yyz_xyy_0_yy, g_yyz_xyy_0_yz, g_yyz_xyy_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yy_0_xx[i] = -2.0 * g_z_xyy_0_xx[i] * b_exp + 4.0 * g_yyz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_xy[i] = -2.0 * g_z_xyy_0_xy[i] * b_exp + 4.0 * g_yyz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_xz[i] = -2.0 * g_z_xyy_0_xz[i] * b_exp + 4.0 * g_yyz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_yy[i] = -2.0 * g_z_xyy_0_yy[i] * b_exp + 4.0 * g_yyz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_yz[i] = -2.0 * g_z_xyy_0_yz[i] * b_exp + 4.0 * g_yyz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_zz[i] = -2.0 * g_z_xyy_0_zz[i] * b_exp + 4.0 * g_yyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_y_x_0_0_yz_yz_0_xx, g_y_x_0_0_yz_yz_0_xy, g_y_x_0_0_yz_yz_0_xz, g_y_x_0_0_yz_yz_0_yy, g_y_x_0_0_yz_yz_0_yz, g_y_x_0_0_yz_yz_0_zz, g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yz_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_y_x_0_0_yz_zz_0_xx, g_y_x_0_0_yz_zz_0_xy, g_y_x_0_0_yz_zz_0_xz, g_y_x_0_0_yz_zz_0_yy, g_y_x_0_0_yz_zz_0_yz, g_y_x_0_0_yz_zz_0_zz, g_yyz_xzz_0_xx, g_yyz_xzz_0_xy, g_yyz_xzz_0_xz, g_yyz_xzz_0_yy, g_yyz_xzz_0_yz, g_yyz_xzz_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_zz_0_xx[i] = -2.0 * g_z_xzz_0_xx[i] * b_exp + 4.0 * g_yyz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_xy[i] = -2.0 * g_z_xzz_0_xy[i] * b_exp + 4.0 * g_yyz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_xz[i] = -2.0 * g_z_xzz_0_xz[i] * b_exp + 4.0 * g_yyz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_yy[i] = -2.0 * g_z_xzz_0_yy[i] * b_exp + 4.0 * g_yyz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_yz[i] = -2.0 * g_z_xzz_0_yz[i] * b_exp + 4.0 * g_yyz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_zz[i] = -2.0 * g_z_xzz_0_zz[i] * b_exp + 4.0 * g_yyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_0_xx, g_y_x_0_0_zz_xx_0_xy, g_y_x_0_0_zz_xx_0_xz, g_y_x_0_0_zz_xx_0_yy, g_y_x_0_0_zz_xx_0_yz, g_y_x_0_0_zz_xx_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xxx_0_xx, g_yzz_xxx_0_xy, g_yzz_xxx_0_xz, g_yzz_xxx_0_yy, g_yzz_xxx_0_yz, g_yzz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_0_xx[i] = -4.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_xy[i] = -4.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_xz[i] = -4.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_yy[i] = -4.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_yz[i] = -4.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_zz[i] = -4.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_y_x_0_0_zz_xy_0_xx, g_y_x_0_0_zz_xy_0_xy, g_y_x_0_0_zz_xy_0_xz, g_y_x_0_0_zz_xy_0_yy, g_y_x_0_0_zz_xy_0_yz, g_y_x_0_0_zz_xy_0_zz, g_yzz_xxy_0_xx, g_yzz_xxy_0_xy, g_yzz_xxy_0_xz, g_yzz_xxy_0_yy, g_yzz_xxy_0_yz, g_yzz_xxy_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xy_0_xx[i] = -2.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_xy[i] = -2.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_xz[i] = -2.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_yy[i] = -2.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_yz[i] = -2.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_zz[i] = -2.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_y_x_0_0_zz_xz_0_xx, g_y_x_0_0_zz_xz_0_xy, g_y_x_0_0_zz_xz_0_xz, g_y_x_0_0_zz_xz_0_yy, g_y_x_0_0_zz_xz_0_yz, g_y_x_0_0_zz_xz_0_zz, g_yzz_xxz_0_xx, g_yzz_xxz_0_xy, g_yzz_xxz_0_xz, g_yzz_xxz_0_yy, g_yzz_xxz_0_yz, g_yzz_xxz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xz_0_xx[i] = -2.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_xy[i] = -2.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_xz[i] = -2.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_yy[i] = -2.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_yz[i] = -2.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_zz[i] = -2.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_y_x_0_0_zz_yy_0_xx, g_y_x_0_0_zz_yy_0_xy, g_y_x_0_0_zz_yy_0_xz, g_y_x_0_0_zz_yy_0_yy, g_y_x_0_0_zz_yy_0_yz, g_y_x_0_0_zz_yy_0_zz, g_yzz_xyy_0_xx, g_yzz_xyy_0_xy, g_yzz_xyy_0_xz, g_yzz_xyy_0_yy, g_yzz_xyy_0_yz, g_yzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yy_0_xx[i] = 4.0 * g_yzz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_xy[i] = 4.0 * g_yzz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_xz[i] = 4.0 * g_yzz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_yy[i] = 4.0 * g_yzz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_yz[i] = 4.0 * g_yzz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_zz[i] = 4.0 * g_yzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_y_x_0_0_zz_yz_0_xx, g_y_x_0_0_zz_yz_0_xy, g_y_x_0_0_zz_yz_0_xz, g_y_x_0_0_zz_yz_0_yy, g_y_x_0_0_zz_yz_0_yz, g_y_x_0_0_zz_yz_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yz_0_xx[i] = 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_xy[i] = 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_xz[i] = 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_yy[i] = 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_yz[i] = 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_zz[i] = 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_y_x_0_0_zz_zz_0_xx, g_y_x_0_0_zz_zz_0_xy, g_y_x_0_0_zz_zz_0_xz, g_y_x_0_0_zz_zz_0_yy, g_y_x_0_0_zz_zz_0_yz, g_y_x_0_0_zz_zz_0_zz, g_yzz_xzz_0_xx, g_yzz_xzz_0_xy, g_yzz_xzz_0_xz, g_yzz_xzz_0_yy, g_yzz_xzz_0_yz, g_yzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_zz_0_xx[i] = 4.0 * g_yzz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_xy[i] = 4.0 * g_yzz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_xz[i] = 4.0 * g_yzz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_yy[i] = 4.0 * g_yzz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_yz[i] = 4.0 * g_yzz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_zz[i] = 4.0 * g_yzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xxy_xxy_0_xx, g_xxy_xxy_0_xy, g_xxy_xxy_0_xz, g_xxy_xxy_0_yy, g_xxy_xxy_0_yz, g_xxy_xxy_0_zz, g_y_y_0_0_xx_xx_0_xx, g_y_y_0_0_xx_xx_0_xy, g_y_y_0_0_xx_xx_0_xz, g_y_y_0_0_xx_xx_0_yy, g_y_y_0_0_xx_xx_0_yz, g_y_y_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_0_xx[i] = 4.0 * g_xxy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_xy[i] = 4.0 * g_xxy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_xz[i] = 4.0 * g_xxy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_yy[i] = 4.0 * g_xxy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_yz[i] = 4.0 * g_xxy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_zz[i] = 4.0 * g_xxy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xyy_0_xx, g_xxy_xyy_0_xy, g_xxy_xyy_0_xz, g_xxy_xyy_0_yy, g_xxy_xyy_0_yz, g_xxy_xyy_0_zz, g_y_y_0_0_xx_xy_0_xx, g_y_y_0_0_xx_xy_0_xy, g_y_y_0_0_xx_xy_0_xz, g_y_y_0_0_xx_xy_0_yy, g_y_y_0_0_xx_xy_0_yz, g_y_y_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xy_0_xx[i] = -2.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_xy[i] = -2.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_xz[i] = -2.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_yy[i] = -2.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_yz[i] = -2.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_zz[i] = -2.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_y_0_0_xx_xz_0_xx, g_y_y_0_0_xx_xz_0_xy, g_y_y_0_0_xx_xz_0_xz, g_y_y_0_0_xx_xz_0_yy, g_y_y_0_0_xx_xz_0_yz, g_y_y_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xz_0_xx[i] = 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_xy[i] = 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_xz[i] = 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_yy[i] = 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_yz[i] = 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_zz[i] = 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_xxy_yyy_0_xx, g_xxy_yyy_0_xy, g_xxy_yyy_0_xz, g_xxy_yyy_0_yy, g_xxy_yyy_0_yz, g_xxy_yyy_0_zz, g_y_y_0_0_xx_yy_0_xx, g_y_y_0_0_xx_yy_0_xy, g_y_y_0_0_xx_yy_0_xz, g_y_y_0_0_xx_yy_0_yy, g_y_y_0_0_xx_yy_0_yz, g_y_y_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yy_0_xx[i] = -4.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_xy[i] = -4.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_xz[i] = -4.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_yy[i] = -4.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_yz[i] = -4.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_zz[i] = -4.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xxy_yyz_0_xx, g_xxy_yyz_0_xy, g_xxy_yyz_0_xz, g_xxy_yyz_0_yy, g_xxy_yyz_0_yz, g_xxy_yyz_0_zz, g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_y_y_0_0_xx_yz_0_xx, g_y_y_0_0_xx_yz_0_xy, g_y_y_0_0_xx_yz_0_xz, g_y_y_0_0_xx_yz_0_yy, g_y_y_0_0_xx_yz_0_yz, g_y_y_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yz_0_xx[i] = -2.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_xy[i] = -2.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_xz[i] = -2.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_yy[i] = -2.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_yz[i] = -2.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_zz[i] = -2.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xxy_yzz_0_xx, g_xxy_yzz_0_xy, g_xxy_yzz_0_xz, g_xxy_yzz_0_yy, g_xxy_yzz_0_yz, g_xxy_yzz_0_zz, g_y_y_0_0_xx_zz_0_xx, g_y_y_0_0_xx_zz_0_xy, g_y_y_0_0_xx_zz_0_xz, g_y_y_0_0_xx_zz_0_yy, g_y_y_0_0_xx_zz_0_yz, g_y_y_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_zz_0_xx[i] = 4.0 * g_xxy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_xy[i] = 4.0 * g_xxy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_xz[i] = 4.0 * g_xxy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_yy[i] = 4.0 * g_xxy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_yz[i] = 4.0 * g_xxy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_zz[i] = 4.0 * g_xxy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_xyy_xxy_0_xx, g_xyy_xxy_0_xy, g_xyy_xxy_0_xz, g_xyy_xxy_0_yy, g_xyy_xxy_0_yz, g_xyy_xxy_0_zz, g_y_y_0_0_xy_xx_0_xx, g_y_y_0_0_xy_xx_0_xy, g_y_y_0_0_xy_xx_0_xz, g_y_y_0_0_xy_xx_0_yy, g_y_y_0_0_xy_xx_0_yz, g_y_y_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_0_xx[i] = -2.0 * g_x_xxy_0_xx[i] * b_exp + 4.0 * g_xyy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_xy[i] = -2.0 * g_x_xxy_0_xy[i] * b_exp + 4.0 * g_xyy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_xz[i] = -2.0 * g_x_xxy_0_xz[i] * b_exp + 4.0 * g_xyy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_yy[i] = -2.0 * g_x_xxy_0_yy[i] * b_exp + 4.0 * g_xyy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_yz[i] = -2.0 * g_x_xxy_0_yz[i] * b_exp + 4.0 * g_xyy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_zz[i] = -2.0 * g_x_xxy_0_zz[i] * b_exp + 4.0 * g_xyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xyy_0_xx, g_xyy_xyy_0_xy, g_xyy_xyy_0_xz, g_xyy_xyy_0_yy, g_xyy_xyy_0_yz, g_xyy_xyy_0_zz, g_y_y_0_0_xy_xy_0_xx, g_y_y_0_0_xy_xy_0_xy, g_y_y_0_0_xy_xy_0_xz, g_y_y_0_0_xy_xy_0_yy, g_y_y_0_0_xy_xy_0_yz, g_y_y_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xy_0_xx[i] = g_x_x_0_xx[i] - 2.0 * g_x_xyy_0_xx[i] * b_exp - 2.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_xy[i] = g_x_x_0_xy[i] - 2.0 * g_x_xyy_0_xy[i] * b_exp - 2.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_xz[i] = g_x_x_0_xz[i] - 2.0 * g_x_xyy_0_xz[i] * b_exp - 2.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_yy[i] = g_x_x_0_yy[i] - 2.0 * g_x_xyy_0_yy[i] * b_exp - 2.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_yz[i] = g_x_x_0_yz[i] - 2.0 * g_x_xyy_0_yz[i] * b_exp - 2.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_zz[i] = g_x_x_0_zz[i] - 2.0 * g_x_xyy_0_zz[i] * b_exp - 2.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz, g_y_y_0_0_xy_xz_0_xx, g_y_y_0_0_xy_xz_0_xy, g_y_y_0_0_xy_xz_0_xz, g_y_y_0_0_xy_xz_0_yy, g_y_y_0_0_xy_xz_0_yz, g_y_y_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xz_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yyy_0_xx, g_x_yyy_0_xy, g_x_yyy_0_xz, g_x_yyy_0_yy, g_x_yyy_0_yz, g_x_yyy_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_xyy_yyy_0_xx, g_xyy_yyy_0_xy, g_xyy_yyy_0_xz, g_xyy_yyy_0_yy, g_xyy_yyy_0_yz, g_xyy_yyy_0_zz, g_y_y_0_0_xy_yy_0_xx, g_y_y_0_0_xy_yy_0_xy, g_y_y_0_0_xy_yy_0_xz, g_y_y_0_0_xy_yy_0_yy, g_y_y_0_0_xy_yy_0_yz, g_y_y_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yy_0_xx[i] = 2.0 * g_x_y_0_xx[i] - 2.0 * g_x_yyy_0_xx[i] * b_exp - 4.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_xy[i] = 2.0 * g_x_y_0_xy[i] - 2.0 * g_x_yyy_0_xy[i] * b_exp - 4.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_xz[i] = 2.0 * g_x_y_0_xz[i] - 2.0 * g_x_yyy_0_xz[i] * b_exp - 4.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_yy[i] = 2.0 * g_x_y_0_yy[i] - 2.0 * g_x_yyy_0_yy[i] * b_exp - 4.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_yz[i] = 2.0 * g_x_y_0_yz[i] - 2.0 * g_x_yyy_0_yz[i] * b_exp - 4.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_zz[i] = 2.0 * g_x_y_0_zz[i] - 2.0 * g_x_yyy_0_zz[i] * b_exp - 4.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xyy_yyz_0_xx, g_xyy_yyz_0_xy, g_xyy_yyz_0_xz, g_xyy_yyz_0_yy, g_xyy_yyz_0_yz, g_xyy_yyz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz, g_y_y_0_0_xy_yz_0_xx, g_y_y_0_0_xy_yz_0_xy, g_y_y_0_0_xy_yz_0_xz, g_y_y_0_0_xy_yz_0_yy, g_y_y_0_0_xy_yz_0_yz, g_y_y_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yz_0_xx[i] = g_x_z_0_xx[i] - 2.0 * g_x_yyz_0_xx[i] * b_exp - 2.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_xy[i] = g_x_z_0_xy[i] - 2.0 * g_x_yyz_0_xy[i] * b_exp - 2.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_xz[i] = g_x_z_0_xz[i] - 2.0 * g_x_yyz_0_xz[i] * b_exp - 2.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_yy[i] = g_x_z_0_yy[i] - 2.0 * g_x_yyz_0_yy[i] * b_exp - 2.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_yz[i] = g_x_z_0_yz[i] - 2.0 * g_x_yyz_0_yz[i] * b_exp - 2.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_zz[i] = g_x_z_0_zz[i] - 2.0 * g_x_yyz_0_zz[i] * b_exp - 2.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_xyy_yzz_0_xx, g_xyy_yzz_0_xy, g_xyy_yzz_0_xz, g_xyy_yzz_0_yy, g_xyy_yzz_0_yz, g_xyy_yzz_0_zz, g_y_y_0_0_xy_zz_0_xx, g_y_y_0_0_xy_zz_0_xy, g_y_y_0_0_xy_zz_0_xz, g_y_y_0_0_xy_zz_0_yy, g_y_y_0_0_xy_zz_0_yz, g_y_y_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_zz_0_xx[i] = -2.0 * g_x_yzz_0_xx[i] * b_exp + 4.0 * g_xyy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_xy[i] = -2.0 * g_x_yzz_0_xy[i] * b_exp + 4.0 * g_xyy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_xz[i] = -2.0 * g_x_yzz_0_xz[i] * b_exp + 4.0 * g_xyy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_yy[i] = -2.0 * g_x_yzz_0_yy[i] * b_exp + 4.0 * g_xyy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_yz[i] = -2.0 * g_x_yzz_0_yz[i] * b_exp + 4.0 * g_xyy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_zz[i] = -2.0 * g_x_yzz_0_zz[i] * b_exp + 4.0 * g_xyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz, g_y_y_0_0_xz_xx_0_xx, g_y_y_0_0_xz_xx_0_xy, g_y_y_0_0_xz_xx_0_xz, g_y_y_0_0_xz_xx_0_yy, g_y_y_0_0_xz_xx_0_yz, g_y_y_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_0_xx[i] = 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_xy[i] = 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_xz[i] = 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_yy[i] = 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_yz[i] = 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_zz[i] = 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz, g_y_y_0_0_xz_xy_0_xx, g_y_y_0_0_xz_xy_0_xy, g_y_y_0_0_xz_xy_0_xz, g_y_y_0_0_xz_xy_0_yy, g_y_y_0_0_xz_xy_0_yz, g_y_y_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xy_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_y_y_0_0_xz_xz_0_xx, g_y_y_0_0_xz_xz_0_xy, g_y_y_0_0_xz_xz_0_xz, g_y_y_0_0_xz_xz_0_yy, g_y_y_0_0_xz_xz_0_yz, g_y_y_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yyy_0_xx, g_xyz_yyy_0_xy, g_xyz_yyy_0_xz, g_xyz_yyy_0_yy, g_xyz_yyy_0_yz, g_xyz_yyy_0_zz, g_y_y_0_0_xz_yy_0_xx, g_y_y_0_0_xz_yy_0_xy, g_y_y_0_0_xz_yy_0_xz, g_y_y_0_0_xz_yy_0_yy, g_y_y_0_0_xz_yy_0_yz, g_y_y_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yy_0_xx[i] = -4.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_xy[i] = -4.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_xz[i] = -4.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_yy[i] = -4.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_yz[i] = -4.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_zz[i] = -4.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_y_y_0_0_xz_yz_0_xx, g_y_y_0_0_xz_yz_0_xy, g_y_y_0_0_xz_yz_0_xz, g_y_y_0_0_xz_yz_0_yy, g_y_y_0_0_xz_yz_0_yz, g_y_y_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz, g_y_y_0_0_xz_zz_0_xx, g_y_y_0_0_xz_zz_0_xy, g_y_y_0_0_xz_zz_0_xz, g_y_y_0_0_xz_zz_0_yy, g_y_y_0_0_xz_zz_0_yz, g_y_y_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_zz_0_xx[i] = 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_xy[i] = 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_xz[i] = 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_yy[i] = 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_yz[i] = 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_zz[i] = 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_0_yy_xx_0_xx, g_y_y_0_0_yy_xx_0_xy, g_y_y_0_0_yy_xx_0_xz, g_y_y_0_0_yy_xx_0_yy, g_y_y_0_0_yy_xx_0_yz, g_y_y_0_0_yy_xx_0_zz, g_yyy_xxy_0_xx, g_yyy_xxy_0_xy, g_yyy_xxy_0_xz, g_yyy_xxy_0_yy, g_yyy_xxy_0_yz, g_yyy_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_0_xx[i] = -4.0 * g_y_xxy_0_xx[i] * b_exp + 4.0 * g_yyy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_xy[i] = -4.0 * g_y_xxy_0_xy[i] * b_exp + 4.0 * g_yyy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_xz[i] = -4.0 * g_y_xxy_0_xz[i] * b_exp + 4.0 * g_yyy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_yy[i] = -4.0 * g_y_xxy_0_yy[i] * b_exp + 4.0 * g_yyy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_yz[i] = -4.0 * g_y_xxy_0_yz[i] * b_exp + 4.0 * g_yyy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_zz[i] = -4.0 * g_y_xxy_0_zz[i] * b_exp + 4.0 * g_yyy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz, g_y_y_0_0_yy_xy_0_xx, g_y_y_0_0_yy_xy_0_xy, g_y_y_0_0_yy_xy_0_xz, g_y_y_0_0_yy_xy_0_yy, g_y_y_0_0_yy_xy_0_yz, g_y_y_0_0_yy_xy_0_zz, g_yyy_x_0_xx, g_yyy_x_0_xy, g_yyy_x_0_xz, g_yyy_x_0_yy, g_yyy_x_0_yz, g_yyy_x_0_zz, g_yyy_xyy_0_xx, g_yyy_xyy_0_xy, g_yyy_xyy_0_xz, g_yyy_xyy_0_yy, g_yyy_xyy_0_yz, g_yyy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xy_0_xx[i] = 2.0 * g_y_x_0_xx[i] - 4.0 * g_y_xyy_0_xx[i] * b_exp - 2.0 * g_yyy_x_0_xx[i] * a_exp + 4.0 * g_yyy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_xy[i] = 2.0 * g_y_x_0_xy[i] - 4.0 * g_y_xyy_0_xy[i] * b_exp - 2.0 * g_yyy_x_0_xy[i] * a_exp + 4.0 * g_yyy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_xz[i] = 2.0 * g_y_x_0_xz[i] - 4.0 * g_y_xyy_0_xz[i] * b_exp - 2.0 * g_yyy_x_0_xz[i] * a_exp + 4.0 * g_yyy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_yy[i] = 2.0 * g_y_x_0_yy[i] - 4.0 * g_y_xyy_0_yy[i] * b_exp - 2.0 * g_yyy_x_0_yy[i] * a_exp + 4.0 * g_yyy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_yz[i] = 2.0 * g_y_x_0_yz[i] - 4.0 * g_y_xyy_0_yz[i] * b_exp - 2.0 * g_yyy_x_0_yz[i] * a_exp + 4.0 * g_yyy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_zz[i] = 2.0 * g_y_x_0_zz[i] - 4.0 * g_y_xyy_0_zz[i] * b_exp - 2.0 * g_yyy_x_0_zz[i] * a_exp + 4.0 * g_yyy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_y_y_0_0_yy_xz_0_xx, g_y_y_0_0_yy_xz_0_xy, g_y_y_0_0_yy_xz_0_xz, g_y_y_0_0_yy_xz_0_yy, g_y_y_0_0_yy_xz_0_yz, g_y_y_0_0_yy_xz_0_zz, g_yyy_xyz_0_xx, g_yyy_xyz_0_xy, g_yyy_xyz_0_xz, g_yyy_xyz_0_yy, g_yyy_xyz_0_yz, g_yyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xz_0_xx[i] = -4.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_xy[i] = -4.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_xz[i] = -4.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_yy[i] = -4.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_yz[i] = -4.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_zz[i] = -4.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_y_y_0_0_yy_yy_0_xx, g_y_y_0_0_yy_yy_0_xy, g_y_y_0_0_yy_yy_0_xz, g_y_y_0_0_yy_yy_0_yy, g_y_y_0_0_yy_yy_0_yz, g_y_y_0_0_yy_yy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yyy_0_xx, g_y_yyy_0_xy, g_y_yyy_0_xz, g_y_yyy_0_yy, g_y_yyy_0_yz, g_y_yyy_0_zz, g_yyy_y_0_xx, g_yyy_y_0_xy, g_yyy_y_0_xz, g_yyy_y_0_yy, g_yyy_y_0_yz, g_yyy_y_0_zz, g_yyy_yyy_0_xx, g_yyy_yyy_0_xy, g_yyy_yyy_0_xz, g_yyy_yyy_0_yy, g_yyy_yyy_0_yz, g_yyy_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yy_0_xx[i] = 4.0 * g_y_y_0_xx[i] - 4.0 * g_y_yyy_0_xx[i] * b_exp - 4.0 * g_yyy_y_0_xx[i] * a_exp + 4.0 * g_yyy_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_xy[i] = 4.0 * g_y_y_0_xy[i] - 4.0 * g_y_yyy_0_xy[i] * b_exp - 4.0 * g_yyy_y_0_xy[i] * a_exp + 4.0 * g_yyy_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_xz[i] = 4.0 * g_y_y_0_xz[i] - 4.0 * g_y_yyy_0_xz[i] * b_exp - 4.0 * g_yyy_y_0_xz[i] * a_exp + 4.0 * g_yyy_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_yy[i] = 4.0 * g_y_y_0_yy[i] - 4.0 * g_y_yyy_0_yy[i] * b_exp - 4.0 * g_yyy_y_0_yy[i] * a_exp + 4.0 * g_yyy_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_yz[i] = 4.0 * g_y_y_0_yz[i] - 4.0 * g_y_yyy_0_yz[i] * b_exp - 4.0 * g_yyy_y_0_yz[i] * a_exp + 4.0 * g_yyy_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_zz[i] = 4.0 * g_y_y_0_zz[i] - 4.0 * g_y_yyy_0_zz[i] * b_exp - 4.0 * g_yyy_y_0_zz[i] * a_exp + 4.0 * g_yyy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_y_y_0_0_yy_yz_0_xx, g_y_y_0_0_yy_yz_0_xy, g_y_y_0_0_yy_yz_0_xz, g_y_y_0_0_yy_yz_0_yy, g_y_y_0_0_yy_yz_0_yz, g_y_y_0_0_yy_yz_0_zz, g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yyy_yyz_0_xx, g_yyy_yyz_0_xy, g_yyy_yyz_0_xz, g_yyy_yyz_0_yy, g_yyy_yyz_0_yz, g_yyy_yyz_0_zz, g_yyy_z_0_xx, g_yyy_z_0_xy, g_yyy_z_0_xz, g_yyy_z_0_yy, g_yyy_z_0_yz, g_yyy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yz_0_xx[i] = 2.0 * g_y_z_0_xx[i] - 4.0 * g_y_yyz_0_xx[i] * b_exp - 2.0 * g_yyy_z_0_xx[i] * a_exp + 4.0 * g_yyy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_xy[i] = 2.0 * g_y_z_0_xy[i] - 4.0 * g_y_yyz_0_xy[i] * b_exp - 2.0 * g_yyy_z_0_xy[i] * a_exp + 4.0 * g_yyy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_xz[i] = 2.0 * g_y_z_0_xz[i] - 4.0 * g_y_yyz_0_xz[i] * b_exp - 2.0 * g_yyy_z_0_xz[i] * a_exp + 4.0 * g_yyy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_yy[i] = 2.0 * g_y_z_0_yy[i] - 4.0 * g_y_yyz_0_yy[i] * b_exp - 2.0 * g_yyy_z_0_yy[i] * a_exp + 4.0 * g_yyy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_yz[i] = 2.0 * g_y_z_0_yz[i] - 4.0 * g_y_yyz_0_yz[i] * b_exp - 2.0 * g_yyy_z_0_yz[i] * a_exp + 4.0 * g_yyy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_zz[i] = 2.0 * g_y_z_0_zz[i] - 4.0 * g_y_yyz_0_zz[i] * b_exp - 2.0 * g_yyy_z_0_zz[i] * a_exp + 4.0 * g_yyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_y_y_0_0_yy_zz_0_xx, g_y_y_0_0_yy_zz_0_xy, g_y_y_0_0_yy_zz_0_xz, g_y_y_0_0_yy_zz_0_yy, g_y_y_0_0_yy_zz_0_yz, g_y_y_0_0_yy_zz_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz, g_yyy_yzz_0_xx, g_yyy_yzz_0_xy, g_yyy_yzz_0_xz, g_yyy_yzz_0_yy, g_yyy_yzz_0_yz, g_yyy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_zz_0_xx[i] = -4.0 * g_y_yzz_0_xx[i] * b_exp + 4.0 * g_yyy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_xy[i] = -4.0 * g_y_yzz_0_xy[i] * b_exp + 4.0 * g_yyy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_xz[i] = -4.0 * g_y_yzz_0_xz[i] * b_exp + 4.0 * g_yyy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_yy[i] = -4.0 * g_y_yzz_0_yy[i] * b_exp + 4.0 * g_yyy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_yz[i] = -4.0 * g_y_yzz_0_yz[i] * b_exp + 4.0 * g_yyy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_zz[i] = -4.0 * g_y_yzz_0_zz[i] * b_exp + 4.0 * g_yyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_0_xx, g_y_y_0_0_yz_xx_0_xy, g_y_y_0_0_yz_xx_0_xz, g_y_y_0_0_yz_xx_0_yy, g_y_y_0_0_yz_xx_0_yz, g_y_y_0_0_yz_xx_0_zz, g_yyz_xxy_0_xx, g_yyz_xxy_0_xy, g_yyz_xxy_0_xz, g_yyz_xxy_0_yy, g_yyz_xxy_0_yz, g_yyz_xxy_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_0_xx[i] = -2.0 * g_z_xxy_0_xx[i] * b_exp + 4.0 * g_yyz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_xy[i] = -2.0 * g_z_xxy_0_xy[i] * b_exp + 4.0 * g_yyz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_xz[i] = -2.0 * g_z_xxy_0_xz[i] * b_exp + 4.0 * g_yyz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_yy[i] = -2.0 * g_z_xxy_0_yy[i] * b_exp + 4.0 * g_yyz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_yz[i] = -2.0 * g_z_xxy_0_yz[i] * b_exp + 4.0 * g_yyz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_zz[i] = -2.0 * g_z_xxy_0_zz[i] * b_exp + 4.0 * g_yyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_y_y_0_0_yz_xy_0_xx, g_y_y_0_0_yz_xy_0_xy, g_y_y_0_0_yz_xy_0_xz, g_y_y_0_0_yz_xy_0_yy, g_y_y_0_0_yz_xy_0_yz, g_y_y_0_0_yz_xy_0_zz, g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xyy_0_xx, g_yyz_xyy_0_xy, g_yyz_xyy_0_xz, g_yyz_xyy_0_yy, g_yyz_xyy_0_yz, g_yyz_xyy_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xy_0_xx[i] = g_z_x_0_xx[i] - 2.0 * g_z_xyy_0_xx[i] * b_exp - 2.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_xy[i] = g_z_x_0_xy[i] - 2.0 * g_z_xyy_0_xy[i] * b_exp - 2.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_xz[i] = g_z_x_0_xz[i] - 2.0 * g_z_xyy_0_xz[i] * b_exp - 2.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_yy[i] = g_z_x_0_yy[i] - 2.0 * g_z_xyy_0_yy[i] * b_exp - 2.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_yz[i] = g_z_x_0_yz[i] - 2.0 * g_z_xyy_0_yz[i] * b_exp - 2.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_zz[i] = g_z_x_0_zz[i] - 2.0 * g_z_xyy_0_zz[i] * b_exp - 2.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_y_y_0_0_yz_xz_0_xx, g_y_y_0_0_yz_xz_0_xy, g_y_y_0_0_yz_xz_0_xz, g_y_y_0_0_yz_xz_0_yy, g_y_y_0_0_yz_xz_0_yz, g_y_y_0_0_yz_xz_0_zz, g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xz_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_y_y_0_0_yz_yy_0_xx, g_y_y_0_0_yz_yy_0_xy, g_y_y_0_0_yz_yy_0_xz, g_y_y_0_0_yz_yy_0_yy, g_y_y_0_0_yz_yy_0_yz, g_y_y_0_0_yz_yy_0_zz, g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_yyz_yyy_0_xx, g_yyz_yyy_0_xy, g_yyz_yyy_0_xz, g_yyz_yyy_0_yy, g_yyz_yyy_0_yz, g_yyz_yyy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yyy_0_xx, g_z_yyy_0_xy, g_z_yyy_0_xz, g_z_yyy_0_yy, g_z_yyy_0_yz, g_z_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yy_0_xx[i] = 2.0 * g_z_y_0_xx[i] - 2.0 * g_z_yyy_0_xx[i] * b_exp - 4.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_xy[i] = 2.0 * g_z_y_0_xy[i] - 2.0 * g_z_yyy_0_xy[i] * b_exp - 4.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_xz[i] = 2.0 * g_z_y_0_xz[i] - 2.0 * g_z_yyy_0_xz[i] * b_exp - 4.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_yy[i] = 2.0 * g_z_y_0_yy[i] - 2.0 * g_z_yyy_0_yy[i] * b_exp - 4.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_yz[i] = 2.0 * g_z_y_0_yz[i] - 2.0 * g_z_yyy_0_yz[i] * b_exp - 4.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_zz[i] = 2.0 * g_z_y_0_zz[i] - 2.0 * g_z_yyy_0_zz[i] * b_exp - 4.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_y_y_0_0_yz_yz_0_xx, g_y_y_0_0_yz_yz_0_xy, g_y_y_0_0_yz_yz_0_xz, g_y_y_0_0_yz_yz_0_yy, g_y_y_0_0_yz_yz_0_yz, g_y_y_0_0_yz_yz_0_zz, g_yyz_yyz_0_xx, g_yyz_yyz_0_xy, g_yyz_yyz_0_xz, g_yyz_yyz_0_yy, g_yyz_yyz_0_yz, g_yyz_yyz_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yz_0_xx[i] = g_z_z_0_xx[i] - 2.0 * g_z_yyz_0_xx[i] * b_exp - 2.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_xy[i] = g_z_z_0_xy[i] - 2.0 * g_z_yyz_0_xy[i] * b_exp - 2.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_xz[i] = g_z_z_0_xz[i] - 2.0 * g_z_yyz_0_xz[i] * b_exp - 2.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_yy[i] = g_z_z_0_yy[i] - 2.0 * g_z_yyz_0_yy[i] * b_exp - 2.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_yz[i] = g_z_z_0_yz[i] - 2.0 * g_z_yyz_0_yz[i] * b_exp - 2.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_zz[i] = g_z_z_0_zz[i] - 2.0 * g_z_yyz_0_zz[i] * b_exp - 2.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_y_y_0_0_yz_zz_0_xx, g_y_y_0_0_yz_zz_0_xy, g_y_y_0_0_yz_zz_0_xz, g_y_y_0_0_yz_zz_0_yy, g_y_y_0_0_yz_zz_0_yz, g_y_y_0_0_yz_zz_0_zz, g_yyz_yzz_0_xx, g_yyz_yzz_0_xy, g_yyz_yzz_0_xz, g_yyz_yzz_0_yy, g_yyz_yzz_0_yz, g_yyz_yzz_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_zz_0_xx[i] = -2.0 * g_z_yzz_0_xx[i] * b_exp + 4.0 * g_yyz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_xy[i] = -2.0 * g_z_yzz_0_xy[i] * b_exp + 4.0 * g_yyz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_xz[i] = -2.0 * g_z_yzz_0_xz[i] * b_exp + 4.0 * g_yyz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_yy[i] = -2.0 * g_z_yzz_0_yy[i] * b_exp + 4.0 * g_yyz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_yz[i] = -2.0 * g_z_yzz_0_yz[i] * b_exp + 4.0 * g_yyz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_zz[i] = -2.0 * g_z_yzz_0_zz[i] * b_exp + 4.0 * g_yyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_0_xx, g_y_y_0_0_zz_xx_0_xy, g_y_y_0_0_zz_xx_0_xz, g_y_y_0_0_zz_xx_0_yy, g_y_y_0_0_zz_xx_0_yz, g_y_y_0_0_zz_xx_0_zz, g_yzz_xxy_0_xx, g_yzz_xxy_0_xy, g_yzz_xxy_0_xz, g_yzz_xxy_0_yy, g_yzz_xxy_0_yz, g_yzz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_0_xx[i] = 4.0 * g_yzz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_xy[i] = 4.0 * g_yzz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_xz[i] = 4.0 * g_yzz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_yy[i] = 4.0 * g_yzz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_yz[i] = 4.0 * g_yzz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_zz[i] = 4.0 * g_yzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_y_y_0_0_zz_xy_0_xx, g_y_y_0_0_zz_xy_0_xy, g_y_y_0_0_zz_xy_0_xz, g_y_y_0_0_zz_xy_0_yy, g_y_y_0_0_zz_xy_0_yz, g_y_y_0_0_zz_xy_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xyy_0_xx, g_yzz_xyy_0_xy, g_yzz_xyy_0_xz, g_yzz_xyy_0_yy, g_yzz_xyy_0_yz, g_yzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xy_0_xx[i] = -2.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_xy[i] = -2.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_xz[i] = -2.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_yy[i] = -2.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_yz[i] = -2.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_zz[i] = -2.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_y_y_0_0_zz_xz_0_xx, g_y_y_0_0_zz_xz_0_xy, g_y_y_0_0_zz_xz_0_xz, g_y_y_0_0_zz_xz_0_yy, g_y_y_0_0_zz_xz_0_yz, g_y_y_0_0_zz_xz_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xz_0_xx[i] = 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_xy[i] = 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_xz[i] = 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_yy[i] = 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_yz[i] = 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_zz[i] = 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_y_y_0_0_zz_yy_0_xx, g_y_y_0_0_zz_yy_0_xy, g_y_y_0_0_zz_yy_0_xz, g_y_y_0_0_zz_yy_0_yy, g_y_y_0_0_zz_yy_0_yz, g_y_y_0_0_zz_yy_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_yzz_yyy_0_xx, g_yzz_yyy_0_xy, g_yzz_yyy_0_xz, g_yzz_yyy_0_yy, g_yzz_yyy_0_yz, g_yzz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yy_0_xx[i] = -4.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_xy[i] = -4.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_xz[i] = -4.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_yy[i] = -4.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_yz[i] = -4.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_zz[i] = -4.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_y_y_0_0_zz_yz_0_xx, g_y_y_0_0_zz_yz_0_xy, g_y_y_0_0_zz_yz_0_xz, g_y_y_0_0_zz_yz_0_yy, g_y_y_0_0_zz_yz_0_yz, g_y_y_0_0_zz_yz_0_zz, g_yzz_yyz_0_xx, g_yzz_yyz_0_xy, g_yzz_yyz_0_xz, g_yzz_yyz_0_yy, g_yzz_yyz_0_yz, g_yzz_yyz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yz_0_xx[i] = -2.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_xy[i] = -2.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_xz[i] = -2.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_yy[i] = -2.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_yz[i] = -2.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_zz[i] = -2.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_y_y_0_0_zz_zz_0_xx, g_y_y_0_0_zz_zz_0_xy, g_y_y_0_0_zz_zz_0_xz, g_y_y_0_0_zz_zz_0_yy, g_y_y_0_0_zz_zz_0_yz, g_y_y_0_0_zz_zz_0_zz, g_yzz_yzz_0_xx, g_yzz_yzz_0_xy, g_yzz_yzz_0_xz, g_yzz_yzz_0_yy, g_yzz_yzz_0_yz, g_yzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_zz_0_xx[i] = 4.0 * g_yzz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_xy[i] = 4.0 * g_yzz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_xz[i] = 4.0 * g_yzz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_yy[i] = 4.0 * g_yzz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_yz[i] = 4.0 * g_yzz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_zz[i] = 4.0 * g_yzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_xxy_xxz_0_xx, g_xxy_xxz_0_xy, g_xxy_xxz_0_xz, g_xxy_xxz_0_yy, g_xxy_xxz_0_yz, g_xxy_xxz_0_zz, g_y_z_0_0_xx_xx_0_xx, g_y_z_0_0_xx_xx_0_xy, g_y_z_0_0_xx_xx_0_xz, g_y_z_0_0_xx_xx_0_yy, g_y_z_0_0_xx_xx_0_yz, g_y_z_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_0_xx[i] = 4.0 * g_xxy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_xy[i] = 4.0 * g_xxy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_xz[i] = 4.0 * g_xxy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_yy[i] = 4.0 * g_xxy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_yz[i] = 4.0 * g_xxy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_zz[i] = 4.0 * g_xxy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_xxy_xyz_0_xx, g_xxy_xyz_0_xy, g_xxy_xyz_0_xz, g_xxy_xyz_0_yy, g_xxy_xyz_0_yz, g_xxy_xyz_0_zz, g_y_z_0_0_xx_xy_0_xx, g_y_z_0_0_xx_xy_0_xy, g_y_z_0_0_xx_xy_0_xz, g_y_z_0_0_xx_xy_0_yy, g_y_z_0_0_xx_xy_0_yz, g_y_z_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xy_0_xx[i] = 4.0 * g_xxy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_xy[i] = 4.0 * g_xxy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_xz[i] = 4.0 * g_xxy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_yy[i] = 4.0 * g_xxy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_yz[i] = 4.0 * g_xxy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_zz[i] = 4.0 * g_xxy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_xxy_x_0_xx, g_xxy_x_0_xy, g_xxy_x_0_xz, g_xxy_x_0_yy, g_xxy_x_0_yz, g_xxy_x_0_zz, g_xxy_xzz_0_xx, g_xxy_xzz_0_xy, g_xxy_xzz_0_xz, g_xxy_xzz_0_yy, g_xxy_xzz_0_yz, g_xxy_xzz_0_zz, g_y_z_0_0_xx_xz_0_xx, g_y_z_0_0_xx_xz_0_xy, g_y_z_0_0_xx_xz_0_xz, g_y_z_0_0_xx_xz_0_yy, g_y_z_0_0_xx_xz_0_yz, g_y_z_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xz_0_xx[i] = -2.0 * g_xxy_x_0_xx[i] * a_exp + 4.0 * g_xxy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_xy[i] = -2.0 * g_xxy_x_0_xy[i] * a_exp + 4.0 * g_xxy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_xz[i] = -2.0 * g_xxy_x_0_xz[i] * a_exp + 4.0 * g_xxy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_yy[i] = -2.0 * g_xxy_x_0_yy[i] * a_exp + 4.0 * g_xxy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_yz[i] = -2.0 * g_xxy_x_0_yz[i] * a_exp + 4.0 * g_xxy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_zz[i] = -2.0 * g_xxy_x_0_zz[i] * a_exp + 4.0 * g_xxy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_xxy_yyz_0_xx, g_xxy_yyz_0_xy, g_xxy_yyz_0_xz, g_xxy_yyz_0_yy, g_xxy_yyz_0_yz, g_xxy_yyz_0_zz, g_y_z_0_0_xx_yy_0_xx, g_y_z_0_0_xx_yy_0_xy, g_y_z_0_0_xx_yy_0_xz, g_y_z_0_0_xx_yy_0_yy, g_y_z_0_0_xx_yy_0_yz, g_y_z_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yy_0_xx[i] = 4.0 * g_xxy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_xy[i] = 4.0 * g_xxy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_xz[i] = 4.0 * g_xxy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_yy[i] = 4.0 * g_xxy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_yz[i] = 4.0 * g_xxy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_zz[i] = 4.0 * g_xxy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_xxy_y_0_xx, g_xxy_y_0_xy, g_xxy_y_0_xz, g_xxy_y_0_yy, g_xxy_y_0_yz, g_xxy_y_0_zz, g_xxy_yzz_0_xx, g_xxy_yzz_0_xy, g_xxy_yzz_0_xz, g_xxy_yzz_0_yy, g_xxy_yzz_0_yz, g_xxy_yzz_0_zz, g_y_z_0_0_xx_yz_0_xx, g_y_z_0_0_xx_yz_0_xy, g_y_z_0_0_xx_yz_0_xz, g_y_z_0_0_xx_yz_0_yy, g_y_z_0_0_xx_yz_0_yz, g_y_z_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yz_0_xx[i] = -2.0 * g_xxy_y_0_xx[i] * a_exp + 4.0 * g_xxy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_xy[i] = -2.0 * g_xxy_y_0_xy[i] * a_exp + 4.0 * g_xxy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_xz[i] = -2.0 * g_xxy_y_0_xz[i] * a_exp + 4.0 * g_xxy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_yy[i] = -2.0 * g_xxy_y_0_yy[i] * a_exp + 4.0 * g_xxy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_yz[i] = -2.0 * g_xxy_y_0_yz[i] * a_exp + 4.0 * g_xxy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_zz[i] = -2.0 * g_xxy_y_0_zz[i] * a_exp + 4.0 * g_xxy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_xxy_z_0_xx, g_xxy_z_0_xy, g_xxy_z_0_xz, g_xxy_z_0_yy, g_xxy_z_0_yz, g_xxy_z_0_zz, g_xxy_zzz_0_xx, g_xxy_zzz_0_xy, g_xxy_zzz_0_xz, g_xxy_zzz_0_yy, g_xxy_zzz_0_yz, g_xxy_zzz_0_zz, g_y_z_0_0_xx_zz_0_xx, g_y_z_0_0_xx_zz_0_xy, g_y_z_0_0_xx_zz_0_xz, g_y_z_0_0_xx_zz_0_yy, g_y_z_0_0_xx_zz_0_yz, g_y_z_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_zz_0_xx[i] = -4.0 * g_xxy_z_0_xx[i] * a_exp + 4.0 * g_xxy_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_xy[i] = -4.0 * g_xxy_z_0_xy[i] * a_exp + 4.0 * g_xxy_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_xz[i] = -4.0 * g_xxy_z_0_xz[i] * a_exp + 4.0 * g_xxy_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_yy[i] = -4.0 * g_xxy_z_0_yy[i] * a_exp + 4.0 * g_xxy_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_yz[i] = -4.0 * g_xxy_z_0_yz[i] * a_exp + 4.0 * g_xxy_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_zz[i] = -4.0 * g_xxy_z_0_zz[i] * a_exp + 4.0 * g_xxy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_xyy_xxz_0_xx, g_xyy_xxz_0_xy, g_xyy_xxz_0_xz, g_xyy_xxz_0_yy, g_xyy_xxz_0_yz, g_xyy_xxz_0_zz, g_y_z_0_0_xy_xx_0_xx, g_y_z_0_0_xy_xx_0_xy, g_y_z_0_0_xy_xx_0_xz, g_y_z_0_0_xy_xx_0_yy, g_y_z_0_0_xy_xx_0_yz, g_y_z_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_0_xx[i] = -2.0 * g_x_xxz_0_xx[i] * b_exp + 4.0 * g_xyy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_xy[i] = -2.0 * g_x_xxz_0_xy[i] * b_exp + 4.0 * g_xyy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_xz[i] = -2.0 * g_x_xxz_0_xz[i] * b_exp + 4.0 * g_xyy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_yy[i] = -2.0 * g_x_xxz_0_yy[i] * b_exp + 4.0 * g_xyy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_yz[i] = -2.0 * g_x_xxz_0_yz[i] * b_exp + 4.0 * g_xyy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_zz[i] = -2.0 * g_x_xxz_0_zz[i] * b_exp + 4.0 * g_xyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xyy_xyz_0_xx, g_xyy_xyz_0_xy, g_xyy_xyz_0_xz, g_xyy_xyz_0_yy, g_xyy_xyz_0_yz, g_xyy_xyz_0_zz, g_y_z_0_0_xy_xy_0_xx, g_y_z_0_0_xy_xy_0_xy, g_y_z_0_0_xy_xy_0_xz, g_y_z_0_0_xy_xy_0_yy, g_y_z_0_0_xy_xy_0_yz, g_y_z_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xy_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_xyy_x_0_xx, g_xyy_x_0_xy, g_xyy_x_0_xz, g_xyy_x_0_yy, g_xyy_x_0_yz, g_xyy_x_0_zz, g_xyy_xzz_0_xx, g_xyy_xzz_0_xy, g_xyy_xzz_0_xz, g_xyy_xzz_0_yy, g_xyy_xzz_0_yz, g_xyy_xzz_0_zz, g_y_z_0_0_xy_xz_0_xx, g_y_z_0_0_xy_xz_0_xy, g_y_z_0_0_xy_xz_0_xz, g_y_z_0_0_xy_xz_0_yy, g_y_z_0_0_xy_xz_0_yz, g_y_z_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xz_0_xx[i] = g_x_x_0_xx[i] - 2.0 * g_x_xzz_0_xx[i] * b_exp - 2.0 * g_xyy_x_0_xx[i] * a_exp + 4.0 * g_xyy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_xy[i] = g_x_x_0_xy[i] - 2.0 * g_x_xzz_0_xy[i] * b_exp - 2.0 * g_xyy_x_0_xy[i] * a_exp + 4.0 * g_xyy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_xz[i] = g_x_x_0_xz[i] - 2.0 * g_x_xzz_0_xz[i] * b_exp - 2.0 * g_xyy_x_0_xz[i] * a_exp + 4.0 * g_xyy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_yy[i] = g_x_x_0_yy[i] - 2.0 * g_x_xzz_0_yy[i] * b_exp - 2.0 * g_xyy_x_0_yy[i] * a_exp + 4.0 * g_xyy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_yz[i] = g_x_x_0_yz[i] - 2.0 * g_x_xzz_0_yz[i] * b_exp - 2.0 * g_xyy_x_0_yz[i] * a_exp + 4.0 * g_xyy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_zz[i] = g_x_x_0_zz[i] - 2.0 * g_x_xzz_0_zz[i] * b_exp - 2.0 * g_xyy_x_0_zz[i] * a_exp + 4.0 * g_xyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_xyy_yyz_0_xx, g_xyy_yyz_0_xy, g_xyy_yyz_0_xz, g_xyy_yyz_0_yy, g_xyy_yyz_0_yz, g_xyy_yyz_0_zz, g_y_z_0_0_xy_yy_0_xx, g_y_z_0_0_xy_yy_0_xy, g_y_z_0_0_xy_yy_0_xz, g_y_z_0_0_xy_yy_0_yy, g_y_z_0_0_xy_yy_0_yz, g_y_z_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yy_0_xx[i] = -2.0 * g_x_yyz_0_xx[i] * b_exp + 4.0 * g_xyy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_xy[i] = -2.0 * g_x_yyz_0_xy[i] * b_exp + 4.0 * g_xyy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_xz[i] = -2.0 * g_x_yyz_0_xz[i] * b_exp + 4.0 * g_xyy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_yy[i] = -2.0 * g_x_yyz_0_yy[i] * b_exp + 4.0 * g_xyy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_yz[i] = -2.0 * g_x_yyz_0_yz[i] * b_exp + 4.0 * g_xyy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_zz[i] = -2.0 * g_x_yyz_0_zz[i] * b_exp + 4.0 * g_xyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_xyy_y_0_xx, g_xyy_y_0_xy, g_xyy_y_0_xz, g_xyy_y_0_yy, g_xyy_y_0_yz, g_xyy_y_0_zz, g_xyy_yzz_0_xx, g_xyy_yzz_0_xy, g_xyy_yzz_0_xz, g_xyy_yzz_0_yy, g_xyy_yzz_0_yz, g_xyy_yzz_0_zz, g_y_z_0_0_xy_yz_0_xx, g_y_z_0_0_xy_yz_0_xy, g_y_z_0_0_xy_yz_0_xz, g_y_z_0_0_xy_yz_0_yy, g_y_z_0_0_xy_yz_0_yz, g_y_z_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yz_0_xx[i] = g_x_y_0_xx[i] - 2.0 * g_x_yzz_0_xx[i] * b_exp - 2.0 * g_xyy_y_0_xx[i] * a_exp + 4.0 * g_xyy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_xy[i] = g_x_y_0_xy[i] - 2.0 * g_x_yzz_0_xy[i] * b_exp - 2.0 * g_xyy_y_0_xy[i] * a_exp + 4.0 * g_xyy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_xz[i] = g_x_y_0_xz[i] - 2.0 * g_x_yzz_0_xz[i] * b_exp - 2.0 * g_xyy_y_0_xz[i] * a_exp + 4.0 * g_xyy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_yy[i] = g_x_y_0_yy[i] - 2.0 * g_x_yzz_0_yy[i] * b_exp - 2.0 * g_xyy_y_0_yy[i] * a_exp + 4.0 * g_xyy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_yz[i] = g_x_y_0_yz[i] - 2.0 * g_x_yzz_0_yz[i] * b_exp - 2.0 * g_xyy_y_0_yz[i] * a_exp + 4.0 * g_xyy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_zz[i] = g_x_y_0_zz[i] - 2.0 * g_x_yzz_0_zz[i] * b_exp - 2.0 * g_xyy_y_0_zz[i] * a_exp + 4.0 * g_xyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_zzz_0_xx, g_x_zzz_0_xy, g_x_zzz_0_xz, g_x_zzz_0_yy, g_x_zzz_0_yz, g_x_zzz_0_zz, g_xyy_z_0_xx, g_xyy_z_0_xy, g_xyy_z_0_xz, g_xyy_z_0_yy, g_xyy_z_0_yz, g_xyy_z_0_zz, g_xyy_zzz_0_xx, g_xyy_zzz_0_xy, g_xyy_zzz_0_xz, g_xyy_zzz_0_yy, g_xyy_zzz_0_yz, g_xyy_zzz_0_zz, g_y_z_0_0_xy_zz_0_xx, g_y_z_0_0_xy_zz_0_xy, g_y_z_0_0_xy_zz_0_xz, g_y_z_0_0_xy_zz_0_yy, g_y_z_0_0_xy_zz_0_yz, g_y_z_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_zz_0_xx[i] = 2.0 * g_x_z_0_xx[i] - 2.0 * g_x_zzz_0_xx[i] * b_exp - 4.0 * g_xyy_z_0_xx[i] * a_exp + 4.0 * g_xyy_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_xy[i] = 2.0 * g_x_z_0_xy[i] - 2.0 * g_x_zzz_0_xy[i] * b_exp - 4.0 * g_xyy_z_0_xy[i] * a_exp + 4.0 * g_xyy_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_xz[i] = 2.0 * g_x_z_0_xz[i] - 2.0 * g_x_zzz_0_xz[i] * b_exp - 4.0 * g_xyy_z_0_xz[i] * a_exp + 4.0 * g_xyy_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_yy[i] = 2.0 * g_x_z_0_yy[i] - 2.0 * g_x_zzz_0_yy[i] * b_exp - 4.0 * g_xyy_z_0_yy[i] * a_exp + 4.0 * g_xyy_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_yz[i] = 2.0 * g_x_z_0_yz[i] - 2.0 * g_x_zzz_0_yz[i] * b_exp - 4.0 * g_xyy_z_0_yz[i] * a_exp + 4.0 * g_xyy_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_zz[i] = 2.0 * g_x_z_0_zz[i] - 2.0 * g_x_zzz_0_zz[i] * b_exp - 4.0 * g_xyy_z_0_zz[i] * a_exp + 4.0 * g_xyy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz, g_y_z_0_0_xz_xx_0_xx, g_y_z_0_0_xz_xx_0_xy, g_y_z_0_0_xz_xx_0_xz, g_y_z_0_0_xz_xx_0_yy, g_y_z_0_0_xz_xx_0_yz, g_y_z_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_0_xx[i] = 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_xy[i] = 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_xz[i] = 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_yy[i] = 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_yz[i] = 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_zz[i] = 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_y_z_0_0_xz_xy_0_xx, g_y_z_0_0_xz_xy_0_xy, g_y_z_0_0_xz_xy_0_xz, g_y_z_0_0_xz_xy_0_yy, g_y_z_0_0_xz_xy_0_yz, g_y_z_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xy_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz, g_y_z_0_0_xz_xz_0_xx, g_y_z_0_0_xz_xz_0_xy, g_y_z_0_0_xz_xz_0_xz, g_y_z_0_0_xz_xz_0_yy, g_y_z_0_0_xz_xz_0_yz, g_y_z_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xz_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz, g_y_z_0_0_xz_yy_0_xx, g_y_z_0_0_xz_yy_0_xy, g_y_z_0_0_xz_yy_0_xz, g_y_z_0_0_xz_yy_0_yy, g_y_z_0_0_xz_yy_0_yz, g_y_z_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yy_0_xx[i] = 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_xy[i] = 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_xz[i] = 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_yy[i] = 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_yz[i] = 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_zz[i] = 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz, g_y_z_0_0_xz_yz_0_xx, g_y_z_0_0_xz_yz_0_xy, g_y_z_0_0_xz_yz_0_xz, g_y_z_0_0_xz_yz_0_yy, g_y_z_0_0_xz_yz_0_yz, g_y_z_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yz_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_xyz_zzz_0_xx, g_xyz_zzz_0_xy, g_xyz_zzz_0_xz, g_xyz_zzz_0_yy, g_xyz_zzz_0_yz, g_xyz_zzz_0_zz, g_y_z_0_0_xz_zz_0_xx, g_y_z_0_0_xz_zz_0_xy, g_y_z_0_0_xz_zz_0_xz, g_y_z_0_0_xz_zz_0_yy, g_y_z_0_0_xz_zz_0_yz, g_y_z_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_zz_0_xx[i] = -4.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_xy[i] = -4.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_xz[i] = -4.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_yy[i] = -4.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_yz[i] = -4.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_zz[i] = -4.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_0_yy_xx_0_xx, g_y_z_0_0_yy_xx_0_xy, g_y_z_0_0_yy_xx_0_xz, g_y_z_0_0_yy_xx_0_yy, g_y_z_0_0_yy_xx_0_yz, g_y_z_0_0_yy_xx_0_zz, g_yyy_xxz_0_xx, g_yyy_xxz_0_xy, g_yyy_xxz_0_xz, g_yyy_xxz_0_yy, g_yyy_xxz_0_yz, g_yyy_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_0_xx[i] = -4.0 * g_y_xxz_0_xx[i] * b_exp + 4.0 * g_yyy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_xy[i] = -4.0 * g_y_xxz_0_xy[i] * b_exp + 4.0 * g_yyy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_xz[i] = -4.0 * g_y_xxz_0_xz[i] * b_exp + 4.0 * g_yyy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_yy[i] = -4.0 * g_y_xxz_0_yy[i] * b_exp + 4.0 * g_yyy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_yz[i] = -4.0 * g_y_xxz_0_yz[i] * b_exp + 4.0 * g_yyy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_zz[i] = -4.0 * g_y_xxz_0_zz[i] * b_exp + 4.0 * g_yyy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_y_z_0_0_yy_xy_0_xx, g_y_z_0_0_yy_xy_0_xy, g_y_z_0_0_yy_xy_0_xz, g_y_z_0_0_yy_xy_0_yy, g_y_z_0_0_yy_xy_0_yz, g_y_z_0_0_yy_xy_0_zz, g_yyy_xyz_0_xx, g_yyy_xyz_0_xy, g_yyy_xyz_0_xz, g_yyy_xyz_0_yy, g_yyy_xyz_0_yz, g_yyy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xy_0_xx[i] = -4.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yyy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_xy[i] = -4.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yyy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_xz[i] = -4.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yyy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_yy[i] = -4.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yyy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_yz[i] = -4.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yyy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_zz[i] = -4.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yyy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz, g_y_z_0_0_yy_xz_0_xx, g_y_z_0_0_yy_xz_0_xy, g_y_z_0_0_yy_xz_0_xz, g_y_z_0_0_yy_xz_0_yy, g_y_z_0_0_yy_xz_0_yz, g_y_z_0_0_yy_xz_0_zz, g_yyy_x_0_xx, g_yyy_x_0_xy, g_yyy_x_0_xz, g_yyy_x_0_yy, g_yyy_x_0_yz, g_yyy_x_0_zz, g_yyy_xzz_0_xx, g_yyy_xzz_0_xy, g_yyy_xzz_0_xz, g_yyy_xzz_0_yy, g_yyy_xzz_0_yz, g_yyy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xz_0_xx[i] = 2.0 * g_y_x_0_xx[i] - 4.0 * g_y_xzz_0_xx[i] * b_exp - 2.0 * g_yyy_x_0_xx[i] * a_exp + 4.0 * g_yyy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_xy[i] = 2.0 * g_y_x_0_xy[i] - 4.0 * g_y_xzz_0_xy[i] * b_exp - 2.0 * g_yyy_x_0_xy[i] * a_exp + 4.0 * g_yyy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_xz[i] = 2.0 * g_y_x_0_xz[i] - 4.0 * g_y_xzz_0_xz[i] * b_exp - 2.0 * g_yyy_x_0_xz[i] * a_exp + 4.0 * g_yyy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_yy[i] = 2.0 * g_y_x_0_yy[i] - 4.0 * g_y_xzz_0_yy[i] * b_exp - 2.0 * g_yyy_x_0_yy[i] * a_exp + 4.0 * g_yyy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_yz[i] = 2.0 * g_y_x_0_yz[i] - 4.0 * g_y_xzz_0_yz[i] * b_exp - 2.0 * g_yyy_x_0_yz[i] * a_exp + 4.0 * g_yyy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_zz[i] = 2.0 * g_y_x_0_zz[i] - 4.0 * g_y_xzz_0_zz[i] * b_exp - 2.0 * g_yyy_x_0_zz[i] * a_exp + 4.0 * g_yyy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_0_yy_yy_0_xx, g_y_z_0_0_yy_yy_0_xy, g_y_z_0_0_yy_yy_0_xz, g_y_z_0_0_yy_yy_0_yy, g_y_z_0_0_yy_yy_0_yz, g_y_z_0_0_yy_yy_0_zz, g_yyy_yyz_0_xx, g_yyy_yyz_0_xy, g_yyy_yyz_0_xz, g_yyy_yyz_0_yy, g_yyy_yyz_0_yz, g_yyy_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yy_0_xx[i] = -4.0 * g_y_yyz_0_xx[i] * b_exp + 4.0 * g_yyy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_xy[i] = -4.0 * g_y_yyz_0_xy[i] * b_exp + 4.0 * g_yyy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_xz[i] = -4.0 * g_y_yyz_0_xz[i] * b_exp + 4.0 * g_yyy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_yy[i] = -4.0 * g_y_yyz_0_yy[i] * b_exp + 4.0 * g_yyy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_yz[i] = -4.0 * g_y_yyz_0_yz[i] * b_exp + 4.0 * g_yyy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_zz[i] = -4.0 * g_y_yyz_0_zz[i] * b_exp + 4.0 * g_yyy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz, g_y_z_0_0_yy_yz_0_xx, g_y_z_0_0_yy_yz_0_xy, g_y_z_0_0_yy_yz_0_xz, g_y_z_0_0_yy_yz_0_yy, g_y_z_0_0_yy_yz_0_yz, g_y_z_0_0_yy_yz_0_zz, g_yyy_y_0_xx, g_yyy_y_0_xy, g_yyy_y_0_xz, g_yyy_y_0_yy, g_yyy_y_0_yz, g_yyy_y_0_zz, g_yyy_yzz_0_xx, g_yyy_yzz_0_xy, g_yyy_yzz_0_xz, g_yyy_yzz_0_yy, g_yyy_yzz_0_yz, g_yyy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yz_0_xx[i] = 2.0 * g_y_y_0_xx[i] - 4.0 * g_y_yzz_0_xx[i] * b_exp - 2.0 * g_yyy_y_0_xx[i] * a_exp + 4.0 * g_yyy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_xy[i] = 2.0 * g_y_y_0_xy[i] - 4.0 * g_y_yzz_0_xy[i] * b_exp - 2.0 * g_yyy_y_0_xy[i] * a_exp + 4.0 * g_yyy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_xz[i] = 2.0 * g_y_y_0_xz[i] - 4.0 * g_y_yzz_0_xz[i] * b_exp - 2.0 * g_yyy_y_0_xz[i] * a_exp + 4.0 * g_yyy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_yy[i] = 2.0 * g_y_y_0_yy[i] - 4.0 * g_y_yzz_0_yy[i] * b_exp - 2.0 * g_yyy_y_0_yy[i] * a_exp + 4.0 * g_yyy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_yz[i] = 2.0 * g_y_y_0_yz[i] - 4.0 * g_y_yzz_0_yz[i] * b_exp - 2.0 * g_yyy_y_0_yz[i] * a_exp + 4.0 * g_yyy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_zz[i] = 2.0 * g_y_y_0_zz[i] - 4.0 * g_y_yzz_0_zz[i] * b_exp - 2.0 * g_yyy_y_0_zz[i] * a_exp + 4.0 * g_yyy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_z_0_0_yy_zz_0_xx, g_y_z_0_0_yy_zz_0_xy, g_y_z_0_0_yy_zz_0_xz, g_y_z_0_0_yy_zz_0_yy, g_y_z_0_0_yy_zz_0_yz, g_y_z_0_0_yy_zz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_zzz_0_xx, g_y_zzz_0_xy, g_y_zzz_0_xz, g_y_zzz_0_yy, g_y_zzz_0_yz, g_y_zzz_0_zz, g_yyy_z_0_xx, g_yyy_z_0_xy, g_yyy_z_0_xz, g_yyy_z_0_yy, g_yyy_z_0_yz, g_yyy_z_0_zz, g_yyy_zzz_0_xx, g_yyy_zzz_0_xy, g_yyy_zzz_0_xz, g_yyy_zzz_0_yy, g_yyy_zzz_0_yz, g_yyy_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_zz_0_xx[i] = 4.0 * g_y_z_0_xx[i] - 4.0 * g_y_zzz_0_xx[i] * b_exp - 4.0 * g_yyy_z_0_xx[i] * a_exp + 4.0 * g_yyy_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_xy[i] = 4.0 * g_y_z_0_xy[i] - 4.0 * g_y_zzz_0_xy[i] * b_exp - 4.0 * g_yyy_z_0_xy[i] * a_exp + 4.0 * g_yyy_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_xz[i] = 4.0 * g_y_z_0_xz[i] - 4.0 * g_y_zzz_0_xz[i] * b_exp - 4.0 * g_yyy_z_0_xz[i] * a_exp + 4.0 * g_yyy_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_yy[i] = 4.0 * g_y_z_0_yy[i] - 4.0 * g_y_zzz_0_yy[i] * b_exp - 4.0 * g_yyy_z_0_yy[i] * a_exp + 4.0 * g_yyy_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_yz[i] = 4.0 * g_y_z_0_yz[i] - 4.0 * g_y_zzz_0_yz[i] * b_exp - 4.0 * g_yyy_z_0_yz[i] * a_exp + 4.0 * g_yyy_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_zz[i] = 4.0 * g_y_z_0_zz[i] - 4.0 * g_y_zzz_0_zz[i] * b_exp - 4.0 * g_yyy_z_0_zz[i] * a_exp + 4.0 * g_yyy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_0_xx, g_y_z_0_0_yz_xx_0_xy, g_y_z_0_0_yz_xx_0_xz, g_y_z_0_0_yz_xx_0_yy, g_y_z_0_0_yz_xx_0_yz, g_y_z_0_0_yz_xx_0_zz, g_yyz_xxz_0_xx, g_yyz_xxz_0_xy, g_yyz_xxz_0_xz, g_yyz_xxz_0_yy, g_yyz_xxz_0_yz, g_yyz_xxz_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_0_xx[i] = -2.0 * g_z_xxz_0_xx[i] * b_exp + 4.0 * g_yyz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_xy[i] = -2.0 * g_z_xxz_0_xy[i] * b_exp + 4.0 * g_yyz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_xz[i] = -2.0 * g_z_xxz_0_xz[i] * b_exp + 4.0 * g_yyz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_yy[i] = -2.0 * g_z_xxz_0_yy[i] * b_exp + 4.0 * g_yyz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_yz[i] = -2.0 * g_z_xxz_0_yz[i] * b_exp + 4.0 * g_yyz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_zz[i] = -2.0 * g_z_xxz_0_zz[i] * b_exp + 4.0 * g_yyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_z_0_0_yz_xy_0_xx, g_y_z_0_0_yz_xy_0_xy, g_y_z_0_0_yz_xy_0_xz, g_y_z_0_0_yz_xy_0_yy, g_y_z_0_0_yz_xy_0_yz, g_y_z_0_0_yz_xy_0_zz, g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xy_0_xx[i] = -2.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_xy[i] = -2.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_xz[i] = -2.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_yy[i] = -2.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_yz[i] = -2.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_zz[i] = -2.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_z_0_0_yz_xz_0_xx, g_y_z_0_0_yz_xz_0_xy, g_y_z_0_0_yz_xz_0_xz, g_y_z_0_0_yz_xz_0_yy, g_y_z_0_0_yz_xz_0_yz, g_y_z_0_0_yz_xz_0_zz, g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xzz_0_xx, g_yyz_xzz_0_xy, g_yyz_xzz_0_xz, g_yyz_xzz_0_yy, g_yyz_xzz_0_yz, g_yyz_xzz_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xz_0_xx[i] = g_z_x_0_xx[i] - 2.0 * g_z_xzz_0_xx[i] * b_exp - 2.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_xy[i] = g_z_x_0_xy[i] - 2.0 * g_z_xzz_0_xy[i] * b_exp - 2.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_xz[i] = g_z_x_0_xz[i] - 2.0 * g_z_xzz_0_xz[i] * b_exp - 2.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_yy[i] = g_z_x_0_yy[i] - 2.0 * g_z_xzz_0_yy[i] * b_exp - 2.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_yz[i] = g_z_x_0_yz[i] - 2.0 * g_z_xzz_0_yz[i] * b_exp - 2.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_zz[i] = g_z_x_0_zz[i] - 2.0 * g_z_xzz_0_zz[i] * b_exp - 2.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_z_0_0_yz_yy_0_xx, g_y_z_0_0_yz_yy_0_xy, g_y_z_0_0_yz_yy_0_xz, g_y_z_0_0_yz_yy_0_yy, g_y_z_0_0_yz_yy_0_yz, g_y_z_0_0_yz_yy_0_zz, g_yyz_yyz_0_xx, g_yyz_yyz_0_xy, g_yyz_yyz_0_xz, g_yyz_yyz_0_yy, g_yyz_yyz_0_yz, g_yyz_yyz_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yy_0_xx[i] = -2.0 * g_z_yyz_0_xx[i] * b_exp + 4.0 * g_yyz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_xy[i] = -2.0 * g_z_yyz_0_xy[i] * b_exp + 4.0 * g_yyz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_xz[i] = -2.0 * g_z_yyz_0_xz[i] * b_exp + 4.0 * g_yyz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_yy[i] = -2.0 * g_z_yyz_0_yy[i] * b_exp + 4.0 * g_yyz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_yz[i] = -2.0 * g_z_yyz_0_yz[i] * b_exp + 4.0 * g_yyz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_zz[i] = -2.0 * g_z_yyz_0_zz[i] * b_exp + 4.0 * g_yyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_z_0_0_yz_yz_0_xx, g_y_z_0_0_yz_yz_0_xy, g_y_z_0_0_yz_yz_0_xz, g_y_z_0_0_yz_yz_0_yy, g_y_z_0_0_yz_yz_0_yz, g_y_z_0_0_yz_yz_0_zz, g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_yyz_yzz_0_xx, g_yyz_yzz_0_xy, g_yyz_yzz_0_xz, g_yyz_yzz_0_yy, g_yyz_yzz_0_yz, g_yyz_yzz_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yz_0_xx[i] = g_z_y_0_xx[i] - 2.0 * g_z_yzz_0_xx[i] * b_exp - 2.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_xy[i] = g_z_y_0_xy[i] - 2.0 * g_z_yzz_0_xy[i] * b_exp - 2.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_xz[i] = g_z_y_0_xz[i] - 2.0 * g_z_yzz_0_xz[i] * b_exp - 2.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_yy[i] = g_z_y_0_yy[i] - 2.0 * g_z_yzz_0_yy[i] * b_exp - 2.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_yz[i] = g_z_y_0_yz[i] - 2.0 * g_z_yzz_0_yz[i] * b_exp - 2.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_zz[i] = g_z_y_0_zz[i] - 2.0 * g_z_yzz_0_zz[i] * b_exp - 2.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_z_0_0_yz_zz_0_xx, g_y_z_0_0_yz_zz_0_xy, g_y_z_0_0_yz_zz_0_xz, g_y_z_0_0_yz_zz_0_yy, g_y_z_0_0_yz_zz_0_yz, g_y_z_0_0_yz_zz_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_yyz_zzz_0_xx, g_yyz_zzz_0_xy, g_yyz_zzz_0_xz, g_yyz_zzz_0_yy, g_yyz_zzz_0_yz, g_yyz_zzz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_zzz_0_xx, g_z_zzz_0_xy, g_z_zzz_0_xz, g_z_zzz_0_yy, g_z_zzz_0_yz, g_z_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_zz_0_xx[i] = 2.0 * g_z_z_0_xx[i] - 2.0 * g_z_zzz_0_xx[i] * b_exp - 4.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_xy[i] = 2.0 * g_z_z_0_xy[i] - 2.0 * g_z_zzz_0_xy[i] * b_exp - 4.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_xz[i] = 2.0 * g_z_z_0_xz[i] - 2.0 * g_z_zzz_0_xz[i] * b_exp - 4.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_yy[i] = 2.0 * g_z_z_0_yy[i] - 2.0 * g_z_zzz_0_yy[i] * b_exp - 4.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_yz[i] = 2.0 * g_z_z_0_yz[i] - 2.0 * g_z_zzz_0_yz[i] * b_exp - 4.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_zz[i] = 2.0 * g_z_z_0_zz[i] - 2.0 * g_z_zzz_0_zz[i] * b_exp - 4.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_0_xx, g_y_z_0_0_zz_xx_0_xy, g_y_z_0_0_zz_xx_0_xz, g_y_z_0_0_zz_xx_0_yy, g_y_z_0_0_zz_xx_0_yz, g_y_z_0_0_zz_xx_0_zz, g_yzz_xxz_0_xx, g_yzz_xxz_0_xy, g_yzz_xxz_0_xz, g_yzz_xxz_0_yy, g_yzz_xxz_0_yz, g_yzz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_0_xx[i] = 4.0 * g_yzz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_xy[i] = 4.0 * g_yzz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_xz[i] = 4.0 * g_yzz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_yy[i] = 4.0 * g_yzz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_yz[i] = 4.0 * g_yzz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_zz[i] = 4.0 * g_yzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_z_0_0_zz_xy_0_xx, g_y_z_0_0_zz_xy_0_xy, g_y_z_0_0_zz_xy_0_xz, g_y_z_0_0_zz_xy_0_yy, g_y_z_0_0_zz_xy_0_yz, g_y_z_0_0_zz_xy_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xy_0_xx[i] = 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_xy[i] = 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_xz[i] = 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_yy[i] = 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_yz[i] = 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_zz[i] = 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_z_0_0_zz_xz_0_xx, g_y_z_0_0_zz_xz_0_xy, g_y_z_0_0_zz_xz_0_xz, g_y_z_0_0_zz_xz_0_yy, g_y_z_0_0_zz_xz_0_yz, g_y_z_0_0_zz_xz_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xzz_0_xx, g_yzz_xzz_0_xy, g_yzz_xzz_0_xz, g_yzz_xzz_0_yy, g_yzz_xzz_0_yz, g_yzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xz_0_xx[i] = -2.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_xy[i] = -2.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_xz[i] = -2.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_yy[i] = -2.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_yz[i] = -2.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_zz[i] = -2.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_z_0_0_zz_yy_0_xx, g_y_z_0_0_zz_yy_0_xy, g_y_z_0_0_zz_yy_0_xz, g_y_z_0_0_zz_yy_0_yy, g_y_z_0_0_zz_yy_0_yz, g_y_z_0_0_zz_yy_0_zz, g_yzz_yyz_0_xx, g_yzz_yyz_0_xy, g_yzz_yyz_0_xz, g_yzz_yyz_0_yy, g_yzz_yyz_0_yz, g_yzz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yy_0_xx[i] = 4.0 * g_yzz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_xy[i] = 4.0 * g_yzz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_xz[i] = 4.0 * g_yzz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_yy[i] = 4.0 * g_yzz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_yz[i] = 4.0 * g_yzz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_zz[i] = 4.0 * g_yzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_z_0_0_zz_yz_0_xx, g_y_z_0_0_zz_yz_0_xy, g_y_z_0_0_zz_yz_0_xz, g_y_z_0_0_zz_yz_0_yy, g_y_z_0_0_zz_yz_0_yz, g_y_z_0_0_zz_yz_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_yzz_yzz_0_xx, g_yzz_yzz_0_xy, g_yzz_yzz_0_xz, g_yzz_yzz_0_yy, g_yzz_yzz_0_yz, g_yzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yz_0_xx[i] = -2.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_xy[i] = -2.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_xz[i] = -2.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_yy[i] = -2.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_yz[i] = -2.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_zz[i] = -2.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_z_0_0_zz_zz_0_xx, g_y_z_0_0_zz_zz_0_xy, g_y_z_0_0_zz_zz_0_xz, g_y_z_0_0_zz_zz_0_yy, g_y_z_0_0_zz_zz_0_yz, g_y_z_0_0_zz_zz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz, g_yzz_zzz_0_xx, g_yzz_zzz_0_xy, g_yzz_zzz_0_xz, g_yzz_zzz_0_yy, g_yzz_zzz_0_yz, g_yzz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_zz_0_xx[i] = -4.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_xy[i] = -4.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_xz[i] = -4.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_yy[i] = -4.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_yz[i] = -4.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_zz[i] = -4.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xxx_0_xx, g_xxz_xxx_0_xy, g_xxz_xxx_0_xz, g_xxz_xxx_0_yy, g_xxz_xxx_0_yz, g_xxz_xxx_0_zz, g_z_x_0_0_xx_xx_0_xx, g_z_x_0_0_xx_xx_0_xy, g_z_x_0_0_xx_xx_0_xz, g_z_x_0_0_xx_xx_0_yy, g_z_x_0_0_xx_xx_0_yz, g_z_x_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_0_xx[i] = -4.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_xy[i] = -4.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_xz[i] = -4.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_yy[i] = -4.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_yz[i] = -4.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_zz[i] = -4.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xxz_xxy_0_xx, g_xxz_xxy_0_xy, g_xxz_xxy_0_xz, g_xxz_xxy_0_yy, g_xxz_xxy_0_yz, g_xxz_xxy_0_zz, g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_z_x_0_0_xx_xy_0_xx, g_z_x_0_0_xx_xy_0_xy, g_z_x_0_0_xx_xy_0_xz, g_z_x_0_0_xx_xy_0_yy, g_z_x_0_0_xx_xy_0_yz, g_z_x_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xy_0_xx[i] = -2.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_xy[i] = -2.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_xz[i] = -2.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_yy[i] = -2.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_yz[i] = -2.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_zz[i] = -2.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xxz_xxz_0_xx, g_xxz_xxz_0_xy, g_xxz_xxz_0_xz, g_xxz_xxz_0_yy, g_xxz_xxz_0_yz, g_xxz_xxz_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_z_x_0_0_xx_xz_0_xx, g_z_x_0_0_xx_xz_0_xy, g_z_x_0_0_xx_xz_0_xz, g_z_x_0_0_xx_xz_0_yy, g_z_x_0_0_xx_xz_0_yz, g_z_x_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xz_0_xx[i] = -2.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_xy[i] = -2.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_xz[i] = -2.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_yy[i] = -2.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_yz[i] = -2.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_zz[i] = -2.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xxz_xyy_0_xx, g_xxz_xyy_0_xy, g_xxz_xyy_0_xz, g_xxz_xyy_0_yy, g_xxz_xyy_0_yz, g_xxz_xyy_0_zz, g_z_x_0_0_xx_yy_0_xx, g_z_x_0_0_xx_yy_0_xy, g_z_x_0_0_xx_yy_0_xz, g_z_x_0_0_xx_yy_0_yy, g_z_x_0_0_xx_yy_0_yz, g_z_x_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yy_0_xx[i] = 4.0 * g_xxz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_xy[i] = 4.0 * g_xxz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_xz[i] = 4.0 * g_xxz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_yy[i] = 4.0 * g_xxz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_yz[i] = 4.0 * g_xxz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_zz[i] = 4.0 * g_xxz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_x_0_0_xx_yz_0_xx, g_z_x_0_0_xx_yz_0_xy, g_z_x_0_0_xx_yz_0_xz, g_z_x_0_0_xx_yz_0_yy, g_z_x_0_0_xx_yz_0_yz, g_z_x_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yz_0_xx[i] = 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_xy[i] = 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_xz[i] = 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_yy[i] = 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_yz[i] = 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_zz[i] = 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xxz_xzz_0_xx, g_xxz_xzz_0_xy, g_xxz_xzz_0_xz, g_xxz_xzz_0_yy, g_xxz_xzz_0_yz, g_xxz_xzz_0_zz, g_z_x_0_0_xx_zz_0_xx, g_z_x_0_0_xx_zz_0_xy, g_z_x_0_0_xx_zz_0_xz, g_z_x_0_0_xx_zz_0_yy, g_z_x_0_0_xx_zz_0_yz, g_z_x_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_zz_0_xx[i] = 4.0 * g_xxz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_xy[i] = 4.0 * g_xxz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_xz[i] = 4.0 * g_xxz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_yy[i] = 4.0 * g_xxz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_yz[i] = 4.0 * g_xxz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_zz[i] = 4.0 * g_xxz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xxx_0_xx, g_xyz_xxx_0_xy, g_xyz_xxx_0_xz, g_xyz_xxx_0_yy, g_xyz_xxx_0_yz, g_xyz_xxx_0_zz, g_z_x_0_0_xy_xx_0_xx, g_z_x_0_0_xy_xx_0_xy, g_z_x_0_0_xy_xx_0_xz, g_z_x_0_0_xy_xx_0_yy, g_z_x_0_0_xy_xx_0_yz, g_z_x_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_0_xx[i] = -4.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_xy[i] = -4.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_xz[i] = -4.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_yy[i] = -4.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_yz[i] = -4.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_zz[i] = -4.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz, g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_z_x_0_0_xy_xy_0_xx, g_z_x_0_0_xy_xy_0_xy, g_z_x_0_0_xy_xy_0_xz, g_z_x_0_0_xy_xy_0_yy, g_z_x_0_0_xy_xy_0_yz, g_z_x_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xy_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_z_x_0_0_xy_xz_0_xx, g_z_x_0_0_xy_xz_0_xy, g_z_x_0_0_xy_xz_0_xz, g_z_x_0_0_xy_xz_0_yy, g_z_x_0_0_xy_xz_0_yz, g_z_x_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz, g_z_x_0_0_xy_yy_0_xx, g_z_x_0_0_xy_yy_0_xy, g_z_x_0_0_xy_yy_0_xz, g_z_x_0_0_xy_yy_0_yy, g_z_x_0_0_xy_yy_0_yz, g_z_x_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yy_0_xx[i] = 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_xy[i] = 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_xz[i] = 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_yy[i] = 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_yz[i] = 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_zz[i] = 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_z_x_0_0_xy_yz_0_xx, g_z_x_0_0_xy_yz_0_xy, g_z_x_0_0_xy_yz_0_xz, g_z_x_0_0_xy_yz_0_yy, g_z_x_0_0_xy_yz_0_yz, g_z_x_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz, g_z_x_0_0_xy_zz_0_xx, g_z_x_0_0_xy_zz_0_xy, g_z_x_0_0_xy_zz_0_xz, g_z_x_0_0_xy_zz_0_yy, g_z_x_0_0_xy_zz_0_yz, g_z_x_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_zz_0_xx[i] = 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_xy[i] = 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_xz[i] = 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_yy[i] = 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_yz[i] = 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_zz[i] = 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xxx_0_xx, g_x_xxx_0_xy, g_x_xxx_0_xz, g_x_xxx_0_yy, g_x_xxx_0_yz, g_x_xxx_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xxx_0_xx, g_xzz_xxx_0_xy, g_xzz_xxx_0_xz, g_xzz_xxx_0_yy, g_xzz_xxx_0_yz, g_xzz_xxx_0_zz, g_z_x_0_0_xz_xx_0_xx, g_z_x_0_0_xz_xx_0_xy, g_z_x_0_0_xz_xx_0_xz, g_z_x_0_0_xz_xx_0_yy, g_z_x_0_0_xz_xx_0_yz, g_z_x_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_0_xx[i] = 2.0 * g_x_x_0_xx[i] - 2.0 * g_x_xxx_0_xx[i] * b_exp - 4.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_xy[i] = 2.0 * g_x_x_0_xy[i] - 2.0 * g_x_xxx_0_xy[i] * b_exp - 4.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_xz[i] = 2.0 * g_x_x_0_xz[i] - 2.0 * g_x_xxx_0_xz[i] * b_exp - 4.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_yy[i] = 2.0 * g_x_x_0_yy[i] - 2.0 * g_x_xxx_0_yy[i] * b_exp - 4.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_yz[i] = 2.0 * g_x_x_0_yz[i] - 2.0 * g_x_xxx_0_yz[i] * b_exp - 4.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_zz[i] = 2.0 * g_x_x_0_zz[i] - 2.0 * g_x_xxx_0_zz[i] * b_exp - 4.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_xzz_xxy_0_xx, g_xzz_xxy_0_xy, g_xzz_xxy_0_xz, g_xzz_xxy_0_yy, g_xzz_xxy_0_yz, g_xzz_xxy_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_z_x_0_0_xz_xy_0_xx, g_z_x_0_0_xz_xy_0_xy, g_z_x_0_0_xz_xy_0_xz, g_z_x_0_0_xz_xy_0_yy, g_z_x_0_0_xz_xy_0_yz, g_z_x_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xy_0_xx[i] = g_x_y_0_xx[i] - 2.0 * g_x_xxy_0_xx[i] * b_exp - 2.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_xy[i] = g_x_y_0_xy[i] - 2.0 * g_x_xxy_0_xy[i] * b_exp - 2.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_xz[i] = g_x_y_0_xz[i] - 2.0 * g_x_xxy_0_xz[i] * b_exp - 2.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_yy[i] = g_x_y_0_yy[i] - 2.0 * g_x_xxy_0_yy[i] * b_exp - 2.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_yz[i] = g_x_y_0_yz[i] - 2.0 * g_x_xxy_0_yz[i] * b_exp - 2.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_zz[i] = g_x_y_0_zz[i] - 2.0 * g_x_xxy_0_zz[i] * b_exp - 2.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xzz_xxz_0_xx, g_xzz_xxz_0_xy, g_xzz_xxz_0_xz, g_xzz_xxz_0_yy, g_xzz_xxz_0_yz, g_xzz_xxz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz, g_z_x_0_0_xz_xz_0_xx, g_z_x_0_0_xz_xz_0_xy, g_z_x_0_0_xz_xz_0_xz, g_z_x_0_0_xz_xz_0_yy, g_z_x_0_0_xz_xz_0_yz, g_z_x_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xz_0_xx[i] = g_x_z_0_xx[i] - 2.0 * g_x_xxz_0_xx[i] * b_exp - 2.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_xy[i] = g_x_z_0_xy[i] - 2.0 * g_x_xxz_0_xy[i] * b_exp - 2.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_xz[i] = g_x_z_0_xz[i] - 2.0 * g_x_xxz_0_xz[i] * b_exp - 2.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_yy[i] = g_x_z_0_yy[i] - 2.0 * g_x_xxz_0_yy[i] * b_exp - 2.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_yz[i] = g_x_z_0_yz[i] - 2.0 * g_x_xxz_0_yz[i] * b_exp - 2.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_zz[i] = g_x_z_0_zz[i] - 2.0 * g_x_xxz_0_zz[i] * b_exp - 2.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_xzz_xyy_0_xx, g_xzz_xyy_0_xy, g_xzz_xyy_0_xz, g_xzz_xyy_0_yy, g_xzz_xyy_0_yz, g_xzz_xyy_0_zz, g_z_x_0_0_xz_yy_0_xx, g_z_x_0_0_xz_yy_0_xy, g_z_x_0_0_xz_yy_0_xz, g_z_x_0_0_xz_yy_0_yy, g_z_x_0_0_xz_yy_0_yz, g_z_x_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yy_0_xx[i] = -2.0 * g_x_xyy_0_xx[i] * b_exp + 4.0 * g_xzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_xy[i] = -2.0 * g_x_xyy_0_xy[i] * b_exp + 4.0 * g_xzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_xz[i] = -2.0 * g_x_xyy_0_xz[i] * b_exp + 4.0 * g_xzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_yy[i] = -2.0 * g_x_xyy_0_yy[i] * b_exp + 4.0 * g_xzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_yz[i] = -2.0 * g_x_xyy_0_yz[i] * b_exp + 4.0 * g_xzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_zz[i] = -2.0 * g_x_xyy_0_zz[i] * b_exp + 4.0 * g_xzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz, g_z_x_0_0_xz_yz_0_xx, g_z_x_0_0_xz_yz_0_xy, g_z_x_0_0_xz_yz_0_xz, g_z_x_0_0_xz_yz_0_yy, g_z_x_0_0_xz_yz_0_yz, g_z_x_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yz_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_xzz_xzz_0_xx, g_xzz_xzz_0_xy, g_xzz_xzz_0_xz, g_xzz_xzz_0_yy, g_xzz_xzz_0_yz, g_xzz_xzz_0_zz, g_z_x_0_0_xz_zz_0_xx, g_z_x_0_0_xz_zz_0_xy, g_z_x_0_0_xz_zz_0_xz, g_z_x_0_0_xz_zz_0_yy, g_z_x_0_0_xz_zz_0_yz, g_z_x_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_zz_0_xx[i] = -2.0 * g_x_xzz_0_xx[i] * b_exp + 4.0 * g_xzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_xy[i] = -2.0 * g_x_xzz_0_xy[i] * b_exp + 4.0 * g_xzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_xz[i] = -2.0 * g_x_xzz_0_xz[i] * b_exp + 4.0 * g_xzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_yy[i] = -2.0 * g_x_xzz_0_yy[i] * b_exp + 4.0 * g_xzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_yz[i] = -2.0 * g_x_xzz_0_yz[i] * b_exp + 4.0 * g_xzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_zz[i] = -2.0 * g_x_xzz_0_zz[i] * b_exp + 4.0 * g_xzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xxx_0_xx, g_yyz_xxx_0_xy, g_yyz_xxx_0_xz, g_yyz_xxx_0_yy, g_yyz_xxx_0_yz, g_yyz_xxx_0_zz, g_z_x_0_0_yy_xx_0_xx, g_z_x_0_0_yy_xx_0_xy, g_z_x_0_0_yy_xx_0_xz, g_z_x_0_0_yy_xx_0_yy, g_z_x_0_0_yy_xx_0_yz, g_z_x_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_0_xx[i] = -4.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_xy[i] = -4.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_xz[i] = -4.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_yy[i] = -4.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_yz[i] = -4.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_zz[i] = -4.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_yyz_xxy_0_xx, g_yyz_xxy_0_xy, g_yyz_xxy_0_xz, g_yyz_xxy_0_yy, g_yyz_xxy_0_yz, g_yyz_xxy_0_zz, g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_z_x_0_0_yy_xy_0_xx, g_z_x_0_0_yy_xy_0_xy, g_z_x_0_0_yy_xy_0_xz, g_z_x_0_0_yy_xy_0_yy, g_z_x_0_0_yy_xy_0_yz, g_z_x_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xy_0_xx[i] = -2.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_xy[i] = -2.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_xz[i] = -2.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_yy[i] = -2.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_yz[i] = -2.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_zz[i] = -2.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_yyz_xxz_0_xx, g_yyz_xxz_0_xy, g_yyz_xxz_0_xz, g_yyz_xxz_0_yy, g_yyz_xxz_0_yz, g_yyz_xxz_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_z_x_0_0_yy_xz_0_xx, g_z_x_0_0_yy_xz_0_xy, g_z_x_0_0_yy_xz_0_xz, g_z_x_0_0_yy_xz_0_yy, g_z_x_0_0_yy_xz_0_yz, g_z_x_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xz_0_xx[i] = -2.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_xy[i] = -2.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_xz[i] = -2.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_yy[i] = -2.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_yz[i] = -2.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_zz[i] = -2.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_yyz_xyy_0_xx, g_yyz_xyy_0_xy, g_yyz_xyy_0_xz, g_yyz_xyy_0_yy, g_yyz_xyy_0_yz, g_yyz_xyy_0_zz, g_z_x_0_0_yy_yy_0_xx, g_z_x_0_0_yy_yy_0_xy, g_z_x_0_0_yy_yy_0_xz, g_z_x_0_0_yy_yy_0_yy, g_z_x_0_0_yy_yy_0_yz, g_z_x_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yy_0_xx[i] = 4.0 * g_yyz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_xy[i] = 4.0 * g_yyz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_xz[i] = 4.0 * g_yyz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_yy[i] = 4.0 * g_yyz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_yz[i] = 4.0 * g_yyz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_zz[i] = 4.0 * g_yyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_x_0_0_yy_yz_0_xx, g_z_x_0_0_yy_yz_0_xy, g_z_x_0_0_yy_yz_0_xz, g_z_x_0_0_yy_yz_0_yy, g_z_x_0_0_yy_yz_0_yz, g_z_x_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yz_0_xx[i] = 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_xy[i] = 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_xz[i] = 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_yy[i] = 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_yz[i] = 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_zz[i] = 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_yyz_xzz_0_xx, g_yyz_xzz_0_xy, g_yyz_xzz_0_xz, g_yyz_xzz_0_yy, g_yyz_xzz_0_yz, g_yyz_xzz_0_zz, g_z_x_0_0_yy_zz_0_xx, g_z_x_0_0_yy_zz_0_xy, g_z_x_0_0_yy_zz_0_xz, g_z_x_0_0_yy_zz_0_yy, g_z_x_0_0_yy_zz_0_yz, g_z_x_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_zz_0_xx[i] = 4.0 * g_yyz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_xy[i] = 4.0 * g_yyz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_xz[i] = 4.0 * g_yyz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_yy[i] = 4.0 * g_yyz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_yz[i] = 4.0 * g_yyz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_zz[i] = 4.0 * g_yyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xxx_0_xx, g_y_xxx_0_xy, g_y_xxx_0_xz, g_y_xxx_0_yy, g_y_xxx_0_yz, g_y_xxx_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xxx_0_xx, g_yzz_xxx_0_xy, g_yzz_xxx_0_xz, g_yzz_xxx_0_yy, g_yzz_xxx_0_yz, g_yzz_xxx_0_zz, g_z_x_0_0_yz_xx_0_xx, g_z_x_0_0_yz_xx_0_xy, g_z_x_0_0_yz_xx_0_xz, g_z_x_0_0_yz_xx_0_yy, g_z_x_0_0_yz_xx_0_yz, g_z_x_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_0_xx[i] = 2.0 * g_y_x_0_xx[i] - 2.0 * g_y_xxx_0_xx[i] * b_exp - 4.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_xy[i] = 2.0 * g_y_x_0_xy[i] - 2.0 * g_y_xxx_0_xy[i] * b_exp - 4.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_xz[i] = 2.0 * g_y_x_0_xz[i] - 2.0 * g_y_xxx_0_xz[i] * b_exp - 4.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_yy[i] = 2.0 * g_y_x_0_yy[i] - 2.0 * g_y_xxx_0_yy[i] * b_exp - 4.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_yz[i] = 2.0 * g_y_x_0_yz[i] - 2.0 * g_y_xxx_0_yz[i] * b_exp - 4.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_zz[i] = 2.0 * g_y_x_0_zz[i] - 2.0 * g_y_xxx_0_zz[i] * b_exp - 4.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_yzz_xxy_0_xx, g_yzz_xxy_0_xy, g_yzz_xxy_0_xz, g_yzz_xxy_0_yy, g_yzz_xxy_0_yz, g_yzz_xxy_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_z_x_0_0_yz_xy_0_xx, g_z_x_0_0_yz_xy_0_xy, g_z_x_0_0_yz_xy_0_xz, g_z_x_0_0_yz_xy_0_yy, g_z_x_0_0_yz_xy_0_yz, g_z_x_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xy_0_xx[i] = g_y_y_0_xx[i] - 2.0 * g_y_xxy_0_xx[i] * b_exp - 2.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_xy[i] = g_y_y_0_xy[i] - 2.0 * g_y_xxy_0_xy[i] * b_exp - 2.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_xz[i] = g_y_y_0_xz[i] - 2.0 * g_y_xxy_0_xz[i] * b_exp - 2.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_yy[i] = g_y_y_0_yy[i] - 2.0 * g_y_xxy_0_yy[i] * b_exp - 2.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_yz[i] = g_y_y_0_yz[i] - 2.0 * g_y_xxy_0_yz[i] * b_exp - 2.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_zz[i] = g_y_y_0_zz[i] - 2.0 * g_y_xxy_0_zz[i] * b_exp - 2.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yzz_xxz_0_xx, g_yzz_xxz_0_xy, g_yzz_xxz_0_xz, g_yzz_xxz_0_yy, g_yzz_xxz_0_yz, g_yzz_xxz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz, g_z_x_0_0_yz_xz_0_xx, g_z_x_0_0_yz_xz_0_xy, g_z_x_0_0_yz_xz_0_xz, g_z_x_0_0_yz_xz_0_yy, g_z_x_0_0_yz_xz_0_yz, g_z_x_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xz_0_xx[i] = g_y_z_0_xx[i] - 2.0 * g_y_xxz_0_xx[i] * b_exp - 2.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_xy[i] = g_y_z_0_xy[i] - 2.0 * g_y_xxz_0_xy[i] * b_exp - 2.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_xz[i] = g_y_z_0_xz[i] - 2.0 * g_y_xxz_0_xz[i] * b_exp - 2.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_yy[i] = g_y_z_0_yy[i] - 2.0 * g_y_xxz_0_yy[i] * b_exp - 2.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_yz[i] = g_y_z_0_yz[i] - 2.0 * g_y_xxz_0_yz[i] * b_exp - 2.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_zz[i] = g_y_z_0_zz[i] - 2.0 * g_y_xxz_0_zz[i] * b_exp - 2.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz, g_yzz_xyy_0_xx, g_yzz_xyy_0_xy, g_yzz_xyy_0_xz, g_yzz_xyy_0_yy, g_yzz_xyy_0_yz, g_yzz_xyy_0_zz, g_z_x_0_0_yz_yy_0_xx, g_z_x_0_0_yz_yy_0_xy, g_z_x_0_0_yz_yy_0_xz, g_z_x_0_0_yz_yy_0_yy, g_z_x_0_0_yz_yy_0_yz, g_z_x_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yy_0_xx[i] = -2.0 * g_y_xyy_0_xx[i] * b_exp + 4.0 * g_yzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_xy[i] = -2.0 * g_y_xyy_0_xy[i] * b_exp + 4.0 * g_yzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_xz[i] = -2.0 * g_y_xyy_0_xz[i] * b_exp + 4.0 * g_yzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_yy[i] = -2.0 * g_y_xyy_0_yy[i] * b_exp + 4.0 * g_yzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_yz[i] = -2.0 * g_y_xyy_0_yz[i] * b_exp + 4.0 * g_yzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_zz[i] = -2.0 * g_y_xyy_0_zz[i] * b_exp + 4.0 * g_yzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz, g_z_x_0_0_yz_yz_0_xx, g_z_x_0_0_yz_yz_0_xy, g_z_x_0_0_yz_yz_0_xz, g_z_x_0_0_yz_yz_0_yy, g_z_x_0_0_yz_yz_0_yz, g_z_x_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yz_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz, g_yzz_xzz_0_xx, g_yzz_xzz_0_xy, g_yzz_xzz_0_xz, g_yzz_xzz_0_yy, g_yzz_xzz_0_yz, g_yzz_xzz_0_zz, g_z_x_0_0_yz_zz_0_xx, g_z_x_0_0_yz_zz_0_xy, g_z_x_0_0_yz_zz_0_xz, g_z_x_0_0_yz_zz_0_yy, g_z_x_0_0_yz_zz_0_yz, g_z_x_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_zz_0_xx[i] = -2.0 * g_y_xzz_0_xx[i] * b_exp + 4.0 * g_yzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_xy[i] = -2.0 * g_y_xzz_0_xy[i] * b_exp + 4.0 * g_yzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_xz[i] = -2.0 * g_y_xzz_0_xz[i] * b_exp + 4.0 * g_yzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_yy[i] = -2.0 * g_y_xzz_0_yy[i] * b_exp + 4.0 * g_yzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_yz[i] = -2.0 * g_y_xzz_0_yz[i] * b_exp + 4.0 * g_yzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_zz[i] = -2.0 * g_y_xzz_0_zz[i] * b_exp + 4.0 * g_yzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_z_x_0_0_zz_xx_0_xx, g_z_x_0_0_zz_xx_0_xy, g_z_x_0_0_zz_xx_0_xz, g_z_x_0_0_zz_xx_0_yy, g_z_x_0_0_zz_xx_0_yz, g_z_x_0_0_zz_xx_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xxx_0_xx, g_z_xxx_0_xy, g_z_xxx_0_xz, g_z_xxx_0_yy, g_z_xxx_0_yz, g_z_xxx_0_zz, g_zzz_x_0_xx, g_zzz_x_0_xy, g_zzz_x_0_xz, g_zzz_x_0_yy, g_zzz_x_0_yz, g_zzz_x_0_zz, g_zzz_xxx_0_xx, g_zzz_xxx_0_xy, g_zzz_xxx_0_xz, g_zzz_xxx_0_yy, g_zzz_xxx_0_yz, g_zzz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_0_xx[i] = 4.0 * g_z_x_0_xx[i] - 4.0 * g_z_xxx_0_xx[i] * b_exp - 4.0 * g_zzz_x_0_xx[i] * a_exp + 4.0 * g_zzz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_xy[i] = 4.0 * g_z_x_0_xy[i] - 4.0 * g_z_xxx_0_xy[i] * b_exp - 4.0 * g_zzz_x_0_xy[i] * a_exp + 4.0 * g_zzz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_xz[i] = 4.0 * g_z_x_0_xz[i] - 4.0 * g_z_xxx_0_xz[i] * b_exp - 4.0 * g_zzz_x_0_xz[i] * a_exp + 4.0 * g_zzz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_yy[i] = 4.0 * g_z_x_0_yy[i] - 4.0 * g_z_xxx_0_yy[i] * b_exp - 4.0 * g_zzz_x_0_yy[i] * a_exp + 4.0 * g_zzz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_yz[i] = 4.0 * g_z_x_0_yz[i] - 4.0 * g_z_xxx_0_yz[i] * b_exp - 4.0 * g_zzz_x_0_yz[i] * a_exp + 4.0 * g_zzz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_zz[i] = 4.0 * g_z_x_0_zz[i] - 4.0 * g_z_xxx_0_zz[i] * b_exp - 4.0 * g_zzz_x_0_zz[i] * a_exp + 4.0 * g_zzz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_z_x_0_0_zz_xy_0_xx, g_z_x_0_0_zz_xy_0_xy, g_z_x_0_0_zz_xy_0_xz, g_z_x_0_0_zz_xy_0_yy, g_z_x_0_0_zz_xy_0_yz, g_z_x_0_0_zz_xy_0_zz, g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_zzz_xxy_0_xx, g_zzz_xxy_0_xy, g_zzz_xxy_0_xz, g_zzz_xxy_0_yy, g_zzz_xxy_0_yz, g_zzz_xxy_0_zz, g_zzz_y_0_xx, g_zzz_y_0_xy, g_zzz_y_0_xz, g_zzz_y_0_yy, g_zzz_y_0_yz, g_zzz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xy_0_xx[i] = 2.0 * g_z_y_0_xx[i] - 4.0 * g_z_xxy_0_xx[i] * b_exp - 2.0 * g_zzz_y_0_xx[i] * a_exp + 4.0 * g_zzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_xy[i] = 2.0 * g_z_y_0_xy[i] - 4.0 * g_z_xxy_0_xy[i] * b_exp - 2.0 * g_zzz_y_0_xy[i] * a_exp + 4.0 * g_zzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_xz[i] = 2.0 * g_z_y_0_xz[i] - 4.0 * g_z_xxy_0_xz[i] * b_exp - 2.0 * g_zzz_y_0_xz[i] * a_exp + 4.0 * g_zzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_yy[i] = 2.0 * g_z_y_0_yy[i] - 4.0 * g_z_xxy_0_yy[i] * b_exp - 2.0 * g_zzz_y_0_yy[i] * a_exp + 4.0 * g_zzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_yz[i] = 2.0 * g_z_y_0_yz[i] - 4.0 * g_z_xxy_0_yz[i] * b_exp - 2.0 * g_zzz_y_0_yz[i] * a_exp + 4.0 * g_zzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_zz[i] = 2.0 * g_z_y_0_zz[i] - 4.0 * g_z_xxy_0_zz[i] * b_exp - 2.0 * g_zzz_y_0_zz[i] * a_exp + 4.0 * g_zzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_z_x_0_0_zz_xz_0_xx, g_z_x_0_0_zz_xz_0_xy, g_z_x_0_0_zz_xz_0_xz, g_z_x_0_0_zz_xz_0_yy, g_z_x_0_0_zz_xz_0_yz, g_z_x_0_0_zz_xz_0_zz, g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_zzz_xxz_0_xx, g_zzz_xxz_0_xy, g_zzz_xxz_0_xz, g_zzz_xxz_0_yy, g_zzz_xxz_0_yz, g_zzz_xxz_0_zz, g_zzz_z_0_xx, g_zzz_z_0_xy, g_zzz_z_0_xz, g_zzz_z_0_yy, g_zzz_z_0_yz, g_zzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xz_0_xx[i] = 2.0 * g_z_z_0_xx[i] - 4.0 * g_z_xxz_0_xx[i] * b_exp - 2.0 * g_zzz_z_0_xx[i] * a_exp + 4.0 * g_zzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_xy[i] = 2.0 * g_z_z_0_xy[i] - 4.0 * g_z_xxz_0_xy[i] * b_exp - 2.0 * g_zzz_z_0_xy[i] * a_exp + 4.0 * g_zzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_xz[i] = 2.0 * g_z_z_0_xz[i] - 4.0 * g_z_xxz_0_xz[i] * b_exp - 2.0 * g_zzz_z_0_xz[i] * a_exp + 4.0 * g_zzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_yy[i] = 2.0 * g_z_z_0_yy[i] - 4.0 * g_z_xxz_0_yy[i] * b_exp - 2.0 * g_zzz_z_0_yy[i] * a_exp + 4.0 * g_zzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_yz[i] = 2.0 * g_z_z_0_yz[i] - 4.0 * g_z_xxz_0_yz[i] * b_exp - 2.0 * g_zzz_z_0_yz[i] * a_exp + 4.0 * g_zzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_zz[i] = 2.0 * g_z_z_0_zz[i] - 4.0 * g_z_xxz_0_zz[i] * b_exp - 2.0 * g_zzz_z_0_zz[i] * a_exp + 4.0 * g_zzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_z_x_0_0_zz_yy_0_xx, g_z_x_0_0_zz_yy_0_xy, g_z_x_0_0_zz_yy_0_xz, g_z_x_0_0_zz_yy_0_yy, g_z_x_0_0_zz_yy_0_yz, g_z_x_0_0_zz_yy_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz, g_zzz_xyy_0_xx, g_zzz_xyy_0_xy, g_zzz_xyy_0_xz, g_zzz_xyy_0_yy, g_zzz_xyy_0_yz, g_zzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yy_0_xx[i] = -4.0 * g_z_xyy_0_xx[i] * b_exp + 4.0 * g_zzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_xy[i] = -4.0 * g_z_xyy_0_xy[i] * b_exp + 4.0 * g_zzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_xz[i] = -4.0 * g_z_xyy_0_xz[i] * b_exp + 4.0 * g_zzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_yy[i] = -4.0 * g_z_xyy_0_yy[i] * b_exp + 4.0 * g_zzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_yz[i] = -4.0 * g_z_xyy_0_yz[i] * b_exp + 4.0 * g_zzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_zz[i] = -4.0 * g_z_xyy_0_zz[i] * b_exp + 4.0 * g_zzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_z_x_0_0_zz_yz_0_xx, g_z_x_0_0_zz_yz_0_xy, g_z_x_0_0_zz_yz_0_xz, g_z_x_0_0_zz_yz_0_yy, g_z_x_0_0_zz_yz_0_yz, g_z_x_0_0_zz_yz_0_zz, g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz, g_zzz_xyz_0_xx, g_zzz_xyz_0_xy, g_zzz_xyz_0_xz, g_zzz_xyz_0_yy, g_zzz_xyz_0_yz, g_zzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yz_0_xx[i] = -4.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_zzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_xy[i] = -4.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_zzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_xz[i] = -4.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_zzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_yy[i] = -4.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_zzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_yz[i] = -4.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_zzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_zz[i] = -4.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_zzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_z_x_0_0_zz_zz_0_xx, g_z_x_0_0_zz_zz_0_xy, g_z_x_0_0_zz_zz_0_xz, g_z_x_0_0_zz_zz_0_yy, g_z_x_0_0_zz_zz_0_yz, g_z_x_0_0_zz_zz_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz, g_zzz_xzz_0_xx, g_zzz_xzz_0_xy, g_zzz_xzz_0_xz, g_zzz_xzz_0_yy, g_zzz_xzz_0_yz, g_zzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_zz_0_xx[i] = -4.0 * g_z_xzz_0_xx[i] * b_exp + 4.0 * g_zzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_xy[i] = -4.0 * g_z_xzz_0_xy[i] * b_exp + 4.0 * g_zzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_xz[i] = -4.0 * g_z_xzz_0_xz[i] * b_exp + 4.0 * g_zzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_yy[i] = -4.0 * g_z_xzz_0_yy[i] * b_exp + 4.0 * g_zzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_yz[i] = -4.0 * g_z_xzz_0_yz[i] * b_exp + 4.0 * g_zzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_zz[i] = -4.0 * g_z_xzz_0_zz[i] * b_exp + 4.0 * g_zzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_xxz_xxy_0_xx, g_xxz_xxy_0_xy, g_xxz_xxy_0_xz, g_xxz_xxy_0_yy, g_xxz_xxy_0_yz, g_xxz_xxy_0_zz, g_z_y_0_0_xx_xx_0_xx, g_z_y_0_0_xx_xx_0_xy, g_z_y_0_0_xx_xx_0_xz, g_z_y_0_0_xx_xx_0_yy, g_z_y_0_0_xx_xx_0_yz, g_z_y_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_0_xx[i] = 4.0 * g_xxz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_xy[i] = 4.0 * g_xxz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_xz[i] = 4.0 * g_xxz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_yy[i] = 4.0 * g_xxz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_yz[i] = 4.0 * g_xxz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_zz[i] = 4.0 * g_xxz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xyy_0_xx, g_xxz_xyy_0_xy, g_xxz_xyy_0_xz, g_xxz_xyy_0_yy, g_xxz_xyy_0_yz, g_xxz_xyy_0_zz, g_z_y_0_0_xx_xy_0_xx, g_z_y_0_0_xx_xy_0_xy, g_z_y_0_0_xx_xy_0_xz, g_z_y_0_0_xx_xy_0_yy, g_z_y_0_0_xx_xy_0_yz, g_z_y_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xy_0_xx[i] = -2.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_xy[i] = -2.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_xz[i] = -2.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_yy[i] = -2.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_yz[i] = -2.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_zz[i] = -2.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_y_0_0_xx_xz_0_xx, g_z_y_0_0_xx_xz_0_xy, g_z_y_0_0_xx_xz_0_xz, g_z_y_0_0_xx_xz_0_yy, g_z_y_0_0_xx_xz_0_yz, g_z_y_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xz_0_xx[i] = 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_xy[i] = 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_xz[i] = 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_yy[i] = 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_yz[i] = 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_zz[i] = 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_xxz_yyy_0_xx, g_xxz_yyy_0_xy, g_xxz_yyy_0_xz, g_xxz_yyy_0_yy, g_xxz_yyy_0_yz, g_xxz_yyy_0_zz, g_z_y_0_0_xx_yy_0_xx, g_z_y_0_0_xx_yy_0_xy, g_z_y_0_0_xx_yy_0_xz, g_z_y_0_0_xx_yy_0_yy, g_z_y_0_0_xx_yy_0_yz, g_z_y_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yy_0_xx[i] = -4.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_xy[i] = -4.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_xz[i] = -4.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_yy[i] = -4.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_yz[i] = -4.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_zz[i] = -4.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_xxz_yyz_0_xx, g_xxz_yyz_0_xy, g_xxz_yyz_0_xz, g_xxz_yyz_0_yy, g_xxz_yyz_0_yz, g_xxz_yyz_0_zz, g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_z_y_0_0_xx_yz_0_xx, g_z_y_0_0_xx_yz_0_xy, g_z_y_0_0_xx_yz_0_xz, g_z_y_0_0_xx_yz_0_yy, g_z_y_0_0_xx_yz_0_yz, g_z_y_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yz_0_xx[i] = -2.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_xy[i] = -2.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_xz[i] = -2.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_yy[i] = -2.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_yz[i] = -2.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_zz[i] = -2.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_xxz_yzz_0_xx, g_xxz_yzz_0_xy, g_xxz_yzz_0_xz, g_xxz_yzz_0_yy, g_xxz_yzz_0_yz, g_xxz_yzz_0_zz, g_z_y_0_0_xx_zz_0_xx, g_z_y_0_0_xx_zz_0_xy, g_z_y_0_0_xx_zz_0_xz, g_z_y_0_0_xx_zz_0_yy, g_z_y_0_0_xx_zz_0_yz, g_z_y_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_zz_0_xx[i] = 4.0 * g_xxz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_xy[i] = 4.0 * g_xxz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_xz[i] = 4.0 * g_xxz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_yy[i] = 4.0 * g_xxz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_yz[i] = 4.0 * g_xxz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_zz[i] = 4.0 * g_xxz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_xyz_xxy_0_xx, g_xyz_xxy_0_xy, g_xyz_xxy_0_xz, g_xyz_xxy_0_yy, g_xyz_xxy_0_yz, g_xyz_xxy_0_zz, g_z_y_0_0_xy_xx_0_xx, g_z_y_0_0_xy_xx_0_xy, g_z_y_0_0_xy_xx_0_xz, g_z_y_0_0_xy_xx_0_yy, g_z_y_0_0_xy_xx_0_yz, g_z_y_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_0_xx[i] = 4.0 * g_xyz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_xy[i] = 4.0 * g_xyz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_xz[i] = 4.0 * g_xyz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_yy[i] = 4.0 * g_xyz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_yz[i] = 4.0 * g_xyz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_zz[i] = 4.0 * g_xyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xyy_0_xx, g_xyz_xyy_0_xy, g_xyz_xyy_0_xz, g_xyz_xyy_0_yy, g_xyz_xyy_0_yz, g_xyz_xyy_0_zz, g_z_y_0_0_xy_xy_0_xx, g_z_y_0_0_xy_xy_0_xy, g_z_y_0_0_xy_xy_0_xz, g_z_y_0_0_xy_xy_0_yy, g_z_y_0_0_xy_xy_0_yz, g_z_y_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xy_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_z_y_0_0_xy_xz_0_xx, g_z_y_0_0_xy_xz_0_xy, g_z_y_0_0_xy_xz_0_xz, g_z_y_0_0_xy_xz_0_yy, g_z_y_0_0_xy_xz_0_yz, g_z_y_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xz_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yyy_0_xx, g_xyz_yyy_0_xy, g_xyz_yyy_0_xz, g_xyz_yyy_0_yy, g_xyz_yyy_0_yz, g_xyz_yyy_0_zz, g_z_y_0_0_xy_yy_0_xx, g_z_y_0_0_xy_yy_0_xy, g_z_y_0_0_xy_yy_0_xz, g_z_y_0_0_xy_yy_0_yy, g_z_y_0_0_xy_yy_0_yz, g_z_y_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yy_0_xx[i] = -4.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_xy[i] = -4.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_xz[i] = -4.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_yy[i] = -4.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_yz[i] = -4.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_zz[i] = -4.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz, g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_z_y_0_0_xy_yz_0_xx, g_z_y_0_0_xy_yz_0_xy, g_z_y_0_0_xy_yz_0_xz, g_z_y_0_0_xy_yz_0_yy, g_z_y_0_0_xy_yz_0_yz, g_z_y_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yz_0_xx[i] = -2.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_xy[i] = -2.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_xz[i] = -2.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_yy[i] = -2.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_yz[i] = -2.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_zz[i] = -2.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz, g_z_y_0_0_xy_zz_0_xx, g_z_y_0_0_xy_zz_0_xy, g_z_y_0_0_xy_zz_0_xz, g_z_y_0_0_xy_zz_0_yy, g_z_y_0_0_xy_zz_0_yz, g_z_y_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_zz_0_xx[i] = 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_xy[i] = 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_xz[i] = 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_yy[i] = 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_yz[i] = 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_zz[i] = 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_x_xxy_0_xx, g_x_xxy_0_xy, g_x_xxy_0_xz, g_x_xxy_0_yy, g_x_xxy_0_yz, g_x_xxy_0_zz, g_xzz_xxy_0_xx, g_xzz_xxy_0_xy, g_xzz_xxy_0_xz, g_xzz_xxy_0_yy, g_xzz_xxy_0_yz, g_xzz_xxy_0_zz, g_z_y_0_0_xz_xx_0_xx, g_z_y_0_0_xz_xx_0_xy, g_z_y_0_0_xz_xx_0_xz, g_z_y_0_0_xz_xx_0_yy, g_z_y_0_0_xz_xx_0_yz, g_z_y_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_0_xx[i] = -2.0 * g_x_xxy_0_xx[i] * b_exp + 4.0 * g_xzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_xy[i] = -2.0 * g_x_xxy_0_xy[i] * b_exp + 4.0 * g_xzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_xz[i] = -2.0 * g_x_xxy_0_xz[i] * b_exp + 4.0 * g_xzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_yy[i] = -2.0 * g_x_xxy_0_yy[i] * b_exp + 4.0 * g_xzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_yz[i] = -2.0 * g_x_xxy_0_yz[i] * b_exp + 4.0 * g_xzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_zz[i] = -2.0 * g_x_xxy_0_zz[i] * b_exp + 4.0 * g_xzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xyy_0_xx, g_x_xyy_0_xy, g_x_xyy_0_xz, g_x_xyy_0_yy, g_x_xyy_0_yz, g_x_xyy_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xyy_0_xx, g_xzz_xyy_0_xy, g_xzz_xyy_0_xz, g_xzz_xyy_0_yy, g_xzz_xyy_0_yz, g_xzz_xyy_0_zz, g_z_y_0_0_xz_xy_0_xx, g_z_y_0_0_xz_xy_0_xy, g_z_y_0_0_xz_xy_0_xz, g_z_y_0_0_xz_xy_0_yy, g_z_y_0_0_xz_xy_0_yz, g_z_y_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xy_0_xx[i] = g_x_x_0_xx[i] - 2.0 * g_x_xyy_0_xx[i] * b_exp - 2.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_xy[i] = g_x_x_0_xy[i] - 2.0 * g_x_xyy_0_xy[i] * b_exp - 2.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_xz[i] = g_x_x_0_xz[i] - 2.0 * g_x_xyy_0_xz[i] * b_exp - 2.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_yy[i] = g_x_x_0_yy[i] - 2.0 * g_x_xyy_0_yy[i] * b_exp - 2.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_yz[i] = g_x_x_0_yz[i] - 2.0 * g_x_xyy_0_yz[i] * b_exp - 2.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_zz[i] = g_x_x_0_zz[i] - 2.0 * g_x_xyy_0_zz[i] * b_exp - 2.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz, g_z_y_0_0_xz_xz_0_xx, g_z_y_0_0_xz_xz_0_xy, g_z_y_0_0_xz_xz_0_xz, g_z_y_0_0_xz_xz_0_yy, g_z_y_0_0_xz_xz_0_yz, g_z_y_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xz_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yyy_0_xx, g_x_yyy_0_xy, g_x_yyy_0_xz, g_x_yyy_0_yy, g_x_yyy_0_yz, g_x_yyy_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_xzz_yyy_0_xx, g_xzz_yyy_0_xy, g_xzz_yyy_0_xz, g_xzz_yyy_0_yy, g_xzz_yyy_0_yz, g_xzz_yyy_0_zz, g_z_y_0_0_xz_yy_0_xx, g_z_y_0_0_xz_yy_0_xy, g_z_y_0_0_xz_yy_0_xz, g_z_y_0_0_xz_yy_0_yy, g_z_y_0_0_xz_yy_0_yz, g_z_y_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yy_0_xx[i] = 2.0 * g_x_y_0_xx[i] - 2.0 * g_x_yyy_0_xx[i] * b_exp - 4.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_xy[i] = 2.0 * g_x_y_0_xy[i] - 2.0 * g_x_yyy_0_xy[i] * b_exp - 4.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_xz[i] = 2.0 * g_x_y_0_xz[i] - 2.0 * g_x_yyy_0_xz[i] * b_exp - 4.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_yy[i] = 2.0 * g_x_y_0_yy[i] - 2.0 * g_x_yyy_0_yy[i] * b_exp - 4.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_yz[i] = 2.0 * g_x_y_0_yz[i] - 2.0 * g_x_yyy_0_yz[i] * b_exp - 4.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_zz[i] = 2.0 * g_x_y_0_zz[i] - 2.0 * g_x_yyy_0_zz[i] * b_exp - 4.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_xzz_yyz_0_xx, g_xzz_yyz_0_xy, g_xzz_yyz_0_xz, g_xzz_yyz_0_yy, g_xzz_yyz_0_yz, g_xzz_yyz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz, g_z_y_0_0_xz_yz_0_xx, g_z_y_0_0_xz_yz_0_xy, g_z_y_0_0_xz_yz_0_xz, g_z_y_0_0_xz_yz_0_yy, g_z_y_0_0_xz_yz_0_yz, g_z_y_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yz_0_xx[i] = g_x_z_0_xx[i] - 2.0 * g_x_yyz_0_xx[i] * b_exp - 2.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_xy[i] = g_x_z_0_xy[i] - 2.0 * g_x_yyz_0_xy[i] * b_exp - 2.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_xz[i] = g_x_z_0_xz[i] - 2.0 * g_x_yyz_0_xz[i] * b_exp - 2.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_yy[i] = g_x_z_0_yy[i] - 2.0 * g_x_yyz_0_yy[i] * b_exp - 2.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_yz[i] = g_x_z_0_yz[i] - 2.0 * g_x_yyz_0_yz[i] * b_exp - 2.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_zz[i] = g_x_z_0_zz[i] - 2.0 * g_x_yyz_0_zz[i] * b_exp - 2.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_xzz_yzz_0_xx, g_xzz_yzz_0_xy, g_xzz_yzz_0_xz, g_xzz_yzz_0_yy, g_xzz_yzz_0_yz, g_xzz_yzz_0_zz, g_z_y_0_0_xz_zz_0_xx, g_z_y_0_0_xz_zz_0_xy, g_z_y_0_0_xz_zz_0_xz, g_z_y_0_0_xz_zz_0_yy, g_z_y_0_0_xz_zz_0_yz, g_z_y_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_zz_0_xx[i] = -2.0 * g_x_yzz_0_xx[i] * b_exp + 4.0 * g_xzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_xy[i] = -2.0 * g_x_yzz_0_xy[i] * b_exp + 4.0 * g_xzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_xz[i] = -2.0 * g_x_yzz_0_xz[i] * b_exp + 4.0 * g_xzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_yy[i] = -2.0 * g_x_yzz_0_yy[i] * b_exp + 4.0 * g_xzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_yz[i] = -2.0 * g_x_yzz_0_yz[i] * b_exp + 4.0 * g_xzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_zz[i] = -2.0 * g_x_yzz_0_zz[i] * b_exp + 4.0 * g_xzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_yyz_xxy_0_xx, g_yyz_xxy_0_xy, g_yyz_xxy_0_xz, g_yyz_xxy_0_yy, g_yyz_xxy_0_yz, g_yyz_xxy_0_zz, g_z_y_0_0_yy_xx_0_xx, g_z_y_0_0_yy_xx_0_xy, g_z_y_0_0_yy_xx_0_xz, g_z_y_0_0_yy_xx_0_yy, g_z_y_0_0_yy_xx_0_yz, g_z_y_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_0_xx[i] = 4.0 * g_yyz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_xy[i] = 4.0 * g_yyz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_xz[i] = 4.0 * g_yyz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_yy[i] = 4.0 * g_yyz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_yz[i] = 4.0 * g_yyz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_zz[i] = 4.0 * g_yyz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xyy_0_xx, g_yyz_xyy_0_xy, g_yyz_xyy_0_xz, g_yyz_xyy_0_yy, g_yyz_xyy_0_yz, g_yyz_xyy_0_zz, g_z_y_0_0_yy_xy_0_xx, g_z_y_0_0_yy_xy_0_xy, g_z_y_0_0_yy_xy_0_xz, g_z_y_0_0_yy_xy_0_yy, g_z_y_0_0_yy_xy_0_yz, g_z_y_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xy_0_xx[i] = -2.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_xy[i] = -2.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_xz[i] = -2.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_yy[i] = -2.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_yz[i] = -2.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_zz[i] = -2.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_y_0_0_yy_xz_0_xx, g_z_y_0_0_yy_xz_0_xy, g_z_y_0_0_yy_xz_0_xz, g_z_y_0_0_yy_xz_0_yy, g_z_y_0_0_yy_xz_0_yz, g_z_y_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xz_0_xx[i] = 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_xy[i] = 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_xz[i] = 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_yy[i] = 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_yz[i] = 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_zz[i] = 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_yyz_yyy_0_xx, g_yyz_yyy_0_xy, g_yyz_yyy_0_xz, g_yyz_yyy_0_yy, g_yyz_yyy_0_yz, g_yyz_yyy_0_zz, g_z_y_0_0_yy_yy_0_xx, g_z_y_0_0_yy_yy_0_xy, g_z_y_0_0_yy_yy_0_xz, g_z_y_0_0_yy_yy_0_yy, g_z_y_0_0_yy_yy_0_yz, g_z_y_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yy_0_xx[i] = -4.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_xy[i] = -4.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_xz[i] = -4.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_yy[i] = -4.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_yz[i] = -4.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_zz[i] = -4.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_yyz_yyz_0_xx, g_yyz_yyz_0_xy, g_yyz_yyz_0_xz, g_yyz_yyz_0_yy, g_yyz_yyz_0_yz, g_yyz_yyz_0_zz, g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_z_y_0_0_yy_yz_0_xx, g_z_y_0_0_yy_yz_0_xy, g_z_y_0_0_yy_yz_0_xz, g_z_y_0_0_yy_yz_0_yy, g_z_y_0_0_yy_yz_0_yz, g_z_y_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yz_0_xx[i] = -2.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_xy[i] = -2.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_xz[i] = -2.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_yy[i] = -2.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_yz[i] = -2.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_zz[i] = -2.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_yyz_yzz_0_xx, g_yyz_yzz_0_xy, g_yyz_yzz_0_xz, g_yyz_yzz_0_yy, g_yyz_yzz_0_yz, g_yyz_yzz_0_zz, g_z_y_0_0_yy_zz_0_xx, g_z_y_0_0_yy_zz_0_xy, g_z_y_0_0_yy_zz_0_xz, g_z_y_0_0_yy_zz_0_yy, g_z_y_0_0_yy_zz_0_yz, g_z_y_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_zz_0_xx[i] = 4.0 * g_yyz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_xy[i] = 4.0 * g_yyz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_xz[i] = 4.0 * g_yyz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_yy[i] = 4.0 * g_yyz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_yz[i] = 4.0 * g_yyz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_zz[i] = 4.0 * g_yyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_y_xxy_0_xx, g_y_xxy_0_xy, g_y_xxy_0_xz, g_y_xxy_0_yy, g_y_xxy_0_yz, g_y_xxy_0_zz, g_yzz_xxy_0_xx, g_yzz_xxy_0_xy, g_yzz_xxy_0_xz, g_yzz_xxy_0_yy, g_yzz_xxy_0_yz, g_yzz_xxy_0_zz, g_z_y_0_0_yz_xx_0_xx, g_z_y_0_0_yz_xx_0_xy, g_z_y_0_0_yz_xx_0_xz, g_z_y_0_0_yz_xx_0_yy, g_z_y_0_0_yz_xx_0_yz, g_z_y_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_0_xx[i] = -2.0 * g_y_xxy_0_xx[i] * b_exp + 4.0 * g_yzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_xy[i] = -2.0 * g_y_xxy_0_xy[i] * b_exp + 4.0 * g_yzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_xz[i] = -2.0 * g_y_xxy_0_xz[i] * b_exp + 4.0 * g_yzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_yy[i] = -2.0 * g_y_xxy_0_yy[i] * b_exp + 4.0 * g_yzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_yz[i] = -2.0 * g_y_xxy_0_yz[i] * b_exp + 4.0 * g_yzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_zz[i] = -2.0 * g_y_xxy_0_zz[i] * b_exp + 4.0 * g_yzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xyy_0_xx, g_y_xyy_0_xy, g_y_xyy_0_xz, g_y_xyy_0_yy, g_y_xyy_0_yz, g_y_xyy_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xyy_0_xx, g_yzz_xyy_0_xy, g_yzz_xyy_0_xz, g_yzz_xyy_0_yy, g_yzz_xyy_0_yz, g_yzz_xyy_0_zz, g_z_y_0_0_yz_xy_0_xx, g_z_y_0_0_yz_xy_0_xy, g_z_y_0_0_yz_xy_0_xz, g_z_y_0_0_yz_xy_0_yy, g_z_y_0_0_yz_xy_0_yz, g_z_y_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xy_0_xx[i] = g_y_x_0_xx[i] - 2.0 * g_y_xyy_0_xx[i] * b_exp - 2.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_xy[i] = g_y_x_0_xy[i] - 2.0 * g_y_xyy_0_xy[i] * b_exp - 2.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_xz[i] = g_y_x_0_xz[i] - 2.0 * g_y_xyy_0_xz[i] * b_exp - 2.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_yy[i] = g_y_x_0_yy[i] - 2.0 * g_y_xyy_0_yy[i] * b_exp - 2.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_yz[i] = g_y_x_0_yz[i] - 2.0 * g_y_xyy_0_yz[i] * b_exp - 2.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_zz[i] = g_y_x_0_zz[i] - 2.0 * g_y_xyy_0_zz[i] * b_exp - 2.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz, g_z_y_0_0_yz_xz_0_xx, g_z_y_0_0_yz_xz_0_xy, g_z_y_0_0_yz_xz_0_xz, g_z_y_0_0_yz_xz_0_yy, g_z_y_0_0_yz_xz_0_yz, g_z_y_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xz_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yyy_0_xx, g_y_yyy_0_xy, g_y_yyy_0_xz, g_y_yyy_0_yy, g_y_yyy_0_yz, g_y_yyy_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_yzz_yyy_0_xx, g_yzz_yyy_0_xy, g_yzz_yyy_0_xz, g_yzz_yyy_0_yy, g_yzz_yyy_0_yz, g_yzz_yyy_0_zz, g_z_y_0_0_yz_yy_0_xx, g_z_y_0_0_yz_yy_0_xy, g_z_y_0_0_yz_yy_0_xz, g_z_y_0_0_yz_yy_0_yy, g_z_y_0_0_yz_yy_0_yz, g_z_y_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yy_0_xx[i] = 2.0 * g_y_y_0_xx[i] - 2.0 * g_y_yyy_0_xx[i] * b_exp - 4.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_xy[i] = 2.0 * g_y_y_0_xy[i] - 2.0 * g_y_yyy_0_xy[i] * b_exp - 4.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_xz[i] = 2.0 * g_y_y_0_xz[i] - 2.0 * g_y_yyy_0_xz[i] * b_exp - 4.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_yy[i] = 2.0 * g_y_y_0_yy[i] - 2.0 * g_y_yyy_0_yy[i] * b_exp - 4.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_yz[i] = 2.0 * g_y_y_0_yz[i] - 2.0 * g_y_yyy_0_yz[i] * b_exp - 4.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_zz[i] = 2.0 * g_y_y_0_zz[i] - 2.0 * g_y_yyy_0_zz[i] * b_exp - 4.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_yzz_yyz_0_xx, g_yzz_yyz_0_xy, g_yzz_yyz_0_xz, g_yzz_yyz_0_yy, g_yzz_yyz_0_yz, g_yzz_yyz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz, g_z_y_0_0_yz_yz_0_xx, g_z_y_0_0_yz_yz_0_xy, g_z_y_0_0_yz_yz_0_xz, g_z_y_0_0_yz_yz_0_yy, g_z_y_0_0_yz_yz_0_yz, g_z_y_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yz_0_xx[i] = g_y_z_0_xx[i] - 2.0 * g_y_yyz_0_xx[i] * b_exp - 2.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_xy[i] = g_y_z_0_xy[i] - 2.0 * g_y_yyz_0_xy[i] * b_exp - 2.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_xz[i] = g_y_z_0_xz[i] - 2.0 * g_y_yyz_0_xz[i] * b_exp - 2.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_yy[i] = g_y_z_0_yy[i] - 2.0 * g_y_yyz_0_yy[i] * b_exp - 2.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_yz[i] = g_y_z_0_yz[i] - 2.0 * g_y_yyz_0_yz[i] * b_exp - 2.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_zz[i] = g_y_z_0_zz[i] - 2.0 * g_y_yyz_0_zz[i] * b_exp - 2.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz, g_yzz_yzz_0_xx, g_yzz_yzz_0_xy, g_yzz_yzz_0_xz, g_yzz_yzz_0_yy, g_yzz_yzz_0_yz, g_yzz_yzz_0_zz, g_z_y_0_0_yz_zz_0_xx, g_z_y_0_0_yz_zz_0_xy, g_z_y_0_0_yz_zz_0_xz, g_z_y_0_0_yz_zz_0_yy, g_z_y_0_0_yz_zz_0_yz, g_z_y_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_zz_0_xx[i] = -2.0 * g_y_yzz_0_xx[i] * b_exp + 4.0 * g_yzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_xy[i] = -2.0 * g_y_yzz_0_xy[i] * b_exp + 4.0 * g_yzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_xz[i] = -2.0 * g_y_yzz_0_xz[i] * b_exp + 4.0 * g_yzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_yy[i] = -2.0 * g_y_yzz_0_yy[i] * b_exp + 4.0 * g_yzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_yz[i] = -2.0 * g_y_yzz_0_yz[i] * b_exp + 4.0 * g_yzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_zz[i] = -2.0 * g_y_yzz_0_zz[i] * b_exp + 4.0 * g_yzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_z_xxy_0_xx, g_z_xxy_0_xy, g_z_xxy_0_xz, g_z_xxy_0_yy, g_z_xxy_0_yz, g_z_xxy_0_zz, g_z_y_0_0_zz_xx_0_xx, g_z_y_0_0_zz_xx_0_xy, g_z_y_0_0_zz_xx_0_xz, g_z_y_0_0_zz_xx_0_yy, g_z_y_0_0_zz_xx_0_yz, g_z_y_0_0_zz_xx_0_zz, g_zzz_xxy_0_xx, g_zzz_xxy_0_xy, g_zzz_xxy_0_xz, g_zzz_xxy_0_yy, g_zzz_xxy_0_yz, g_zzz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_0_xx[i] = -4.0 * g_z_xxy_0_xx[i] * b_exp + 4.0 * g_zzz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_xy[i] = -4.0 * g_z_xxy_0_xy[i] * b_exp + 4.0 * g_zzz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_xz[i] = -4.0 * g_z_xxy_0_xz[i] * b_exp + 4.0 * g_zzz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_yy[i] = -4.0 * g_z_xxy_0_yy[i] * b_exp + 4.0 * g_zzz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_yz[i] = -4.0 * g_z_xxy_0_yz[i] * b_exp + 4.0 * g_zzz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_zz[i] = -4.0 * g_z_xxy_0_zz[i] * b_exp + 4.0 * g_zzz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xyy_0_xx, g_z_xyy_0_xy, g_z_xyy_0_xz, g_z_xyy_0_yy, g_z_xyy_0_yz, g_z_xyy_0_zz, g_z_y_0_0_zz_xy_0_xx, g_z_y_0_0_zz_xy_0_xy, g_z_y_0_0_zz_xy_0_xz, g_z_y_0_0_zz_xy_0_yy, g_z_y_0_0_zz_xy_0_yz, g_z_y_0_0_zz_xy_0_zz, g_zzz_x_0_xx, g_zzz_x_0_xy, g_zzz_x_0_xz, g_zzz_x_0_yy, g_zzz_x_0_yz, g_zzz_x_0_zz, g_zzz_xyy_0_xx, g_zzz_xyy_0_xy, g_zzz_xyy_0_xz, g_zzz_xyy_0_yy, g_zzz_xyy_0_yz, g_zzz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xy_0_xx[i] = 2.0 * g_z_x_0_xx[i] - 4.0 * g_z_xyy_0_xx[i] * b_exp - 2.0 * g_zzz_x_0_xx[i] * a_exp + 4.0 * g_zzz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_xy[i] = 2.0 * g_z_x_0_xy[i] - 4.0 * g_z_xyy_0_xy[i] * b_exp - 2.0 * g_zzz_x_0_xy[i] * a_exp + 4.0 * g_zzz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_xz[i] = 2.0 * g_z_x_0_xz[i] - 4.0 * g_z_xyy_0_xz[i] * b_exp - 2.0 * g_zzz_x_0_xz[i] * a_exp + 4.0 * g_zzz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_yy[i] = 2.0 * g_z_x_0_yy[i] - 4.0 * g_z_xyy_0_yy[i] * b_exp - 2.0 * g_zzz_x_0_yy[i] * a_exp + 4.0 * g_zzz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_yz[i] = 2.0 * g_z_x_0_yz[i] - 4.0 * g_z_xyy_0_yz[i] * b_exp - 2.0 * g_zzz_x_0_yz[i] * a_exp + 4.0 * g_zzz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_zz[i] = 2.0 * g_z_x_0_zz[i] - 4.0 * g_z_xyy_0_zz[i] * b_exp - 2.0 * g_zzz_x_0_zz[i] * a_exp + 4.0 * g_zzz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz, g_z_y_0_0_zz_xz_0_xx, g_z_y_0_0_zz_xz_0_xy, g_z_y_0_0_zz_xz_0_xz, g_z_y_0_0_zz_xz_0_yy, g_z_y_0_0_zz_xz_0_yz, g_z_y_0_0_zz_xz_0_zz, g_zzz_xyz_0_xx, g_zzz_xyz_0_xy, g_zzz_xyz_0_xz, g_zzz_xyz_0_yy, g_zzz_xyz_0_yz, g_zzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xz_0_xx[i] = -4.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_zzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_xy[i] = -4.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_zzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_xz[i] = -4.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_zzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_yy[i] = -4.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_zzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_yz[i] = -4.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_zzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_zz[i] = -4.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_zzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_z_y_0_0_zz_yy_0_xx, g_z_y_0_0_zz_yy_0_xy, g_z_y_0_0_zz_yy_0_xz, g_z_y_0_0_zz_yy_0_yy, g_z_y_0_0_zz_yy_0_yz, g_z_y_0_0_zz_yy_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yyy_0_xx, g_z_yyy_0_xy, g_z_yyy_0_xz, g_z_yyy_0_yy, g_z_yyy_0_yz, g_z_yyy_0_zz, g_zzz_y_0_xx, g_zzz_y_0_xy, g_zzz_y_0_xz, g_zzz_y_0_yy, g_zzz_y_0_yz, g_zzz_y_0_zz, g_zzz_yyy_0_xx, g_zzz_yyy_0_xy, g_zzz_yyy_0_xz, g_zzz_yyy_0_yy, g_zzz_yyy_0_yz, g_zzz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yy_0_xx[i] = 4.0 * g_z_y_0_xx[i] - 4.0 * g_z_yyy_0_xx[i] * b_exp - 4.0 * g_zzz_y_0_xx[i] * a_exp + 4.0 * g_zzz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_xy[i] = 4.0 * g_z_y_0_xy[i] - 4.0 * g_z_yyy_0_xy[i] * b_exp - 4.0 * g_zzz_y_0_xy[i] * a_exp + 4.0 * g_zzz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_xz[i] = 4.0 * g_z_y_0_xz[i] - 4.0 * g_z_yyy_0_xz[i] * b_exp - 4.0 * g_zzz_y_0_xz[i] * a_exp + 4.0 * g_zzz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_yy[i] = 4.0 * g_z_y_0_yy[i] - 4.0 * g_z_yyy_0_yy[i] * b_exp - 4.0 * g_zzz_y_0_yy[i] * a_exp + 4.0 * g_zzz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_yz[i] = 4.0 * g_z_y_0_yz[i] - 4.0 * g_z_yyy_0_yz[i] * b_exp - 4.0 * g_zzz_y_0_yz[i] * a_exp + 4.0 * g_zzz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_zz[i] = 4.0 * g_z_y_0_zz[i] - 4.0 * g_z_yyy_0_zz[i] * b_exp - 4.0 * g_zzz_y_0_zz[i] * a_exp + 4.0 * g_zzz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_z_y_0_0_zz_yz_0_xx, g_z_y_0_0_zz_yz_0_xy, g_z_y_0_0_zz_yz_0_xz, g_z_y_0_0_zz_yz_0_yy, g_z_y_0_0_zz_yz_0_yz, g_z_y_0_0_zz_yz_0_zz, g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_zzz_yyz_0_xx, g_zzz_yyz_0_xy, g_zzz_yyz_0_xz, g_zzz_yyz_0_yy, g_zzz_yyz_0_yz, g_zzz_yyz_0_zz, g_zzz_z_0_xx, g_zzz_z_0_xy, g_zzz_z_0_xz, g_zzz_z_0_yy, g_zzz_z_0_yz, g_zzz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yz_0_xx[i] = 2.0 * g_z_z_0_xx[i] - 4.0 * g_z_yyz_0_xx[i] * b_exp - 2.0 * g_zzz_z_0_xx[i] * a_exp + 4.0 * g_zzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_xy[i] = 2.0 * g_z_z_0_xy[i] - 4.0 * g_z_yyz_0_xy[i] * b_exp - 2.0 * g_zzz_z_0_xy[i] * a_exp + 4.0 * g_zzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_xz[i] = 2.0 * g_z_z_0_xz[i] - 4.0 * g_z_yyz_0_xz[i] * b_exp - 2.0 * g_zzz_z_0_xz[i] * a_exp + 4.0 * g_zzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_yy[i] = 2.0 * g_z_z_0_yy[i] - 4.0 * g_z_yyz_0_yy[i] * b_exp - 2.0 * g_zzz_z_0_yy[i] * a_exp + 4.0 * g_zzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_yz[i] = 2.0 * g_z_z_0_yz[i] - 4.0 * g_z_yyz_0_yz[i] * b_exp - 2.0 * g_zzz_z_0_yz[i] * a_exp + 4.0 * g_zzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_zz[i] = 2.0 * g_z_z_0_zz[i] - 4.0 * g_z_yyz_0_zz[i] * b_exp - 2.0 * g_zzz_z_0_zz[i] * a_exp + 4.0 * g_zzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_z_y_0_0_zz_zz_0_xx, g_z_y_0_0_zz_zz_0_xy, g_z_y_0_0_zz_zz_0_xz, g_z_y_0_0_zz_zz_0_yy, g_z_y_0_0_zz_zz_0_yz, g_z_y_0_0_zz_zz_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz, g_zzz_yzz_0_xx, g_zzz_yzz_0_xy, g_zzz_yzz_0_xz, g_zzz_yzz_0_yy, g_zzz_yzz_0_yz, g_zzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_zz_0_xx[i] = -4.0 * g_z_yzz_0_xx[i] * b_exp + 4.0 * g_zzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_xy[i] = -4.0 * g_z_yzz_0_xy[i] * b_exp + 4.0 * g_zzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_xz[i] = -4.0 * g_z_yzz_0_xz[i] * b_exp + 4.0 * g_zzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_yy[i] = -4.0 * g_z_yzz_0_yy[i] * b_exp + 4.0 * g_zzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_yz[i] = -4.0 * g_z_yzz_0_yz[i] * b_exp + 4.0 * g_zzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_zz[i] = -4.0 * g_z_yzz_0_zz[i] * b_exp + 4.0 * g_zzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_xxz_xxz_0_xx, g_xxz_xxz_0_xy, g_xxz_xxz_0_xz, g_xxz_xxz_0_yy, g_xxz_xxz_0_yz, g_xxz_xxz_0_zz, g_z_z_0_0_xx_xx_0_xx, g_z_z_0_0_xx_xx_0_xy, g_z_z_0_0_xx_xx_0_xz, g_z_z_0_0_xx_xx_0_yy, g_z_z_0_0_xx_xx_0_yz, g_z_z_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_0_xx[i] = 4.0 * g_xxz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_xy[i] = 4.0 * g_xxz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_xz[i] = 4.0 * g_xxz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_yy[i] = 4.0 * g_xxz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_yz[i] = 4.0 * g_xxz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_zz[i] = 4.0 * g_xxz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_xxz_xyz_0_xx, g_xxz_xyz_0_xy, g_xxz_xyz_0_xz, g_xxz_xyz_0_yy, g_xxz_xyz_0_yz, g_xxz_xyz_0_zz, g_z_z_0_0_xx_xy_0_xx, g_z_z_0_0_xx_xy_0_xy, g_z_z_0_0_xx_xy_0_xz, g_z_z_0_0_xx_xy_0_yy, g_z_z_0_0_xx_xy_0_yz, g_z_z_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xy_0_xx[i] = 4.0 * g_xxz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_xy[i] = 4.0 * g_xxz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_xz[i] = 4.0 * g_xxz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_yy[i] = 4.0 * g_xxz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_yz[i] = 4.0 * g_xxz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_zz[i] = 4.0 * g_xxz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_xxz_x_0_xx, g_xxz_x_0_xy, g_xxz_x_0_xz, g_xxz_x_0_yy, g_xxz_x_0_yz, g_xxz_x_0_zz, g_xxz_xzz_0_xx, g_xxz_xzz_0_xy, g_xxz_xzz_0_xz, g_xxz_xzz_0_yy, g_xxz_xzz_0_yz, g_xxz_xzz_0_zz, g_z_z_0_0_xx_xz_0_xx, g_z_z_0_0_xx_xz_0_xy, g_z_z_0_0_xx_xz_0_xz, g_z_z_0_0_xx_xz_0_yy, g_z_z_0_0_xx_xz_0_yz, g_z_z_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xz_0_xx[i] = -2.0 * g_xxz_x_0_xx[i] * a_exp + 4.0 * g_xxz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_xy[i] = -2.0 * g_xxz_x_0_xy[i] * a_exp + 4.0 * g_xxz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_xz[i] = -2.0 * g_xxz_x_0_xz[i] * a_exp + 4.0 * g_xxz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_yy[i] = -2.0 * g_xxz_x_0_yy[i] * a_exp + 4.0 * g_xxz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_yz[i] = -2.0 * g_xxz_x_0_yz[i] * a_exp + 4.0 * g_xxz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_zz[i] = -2.0 * g_xxz_x_0_zz[i] * a_exp + 4.0 * g_xxz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_xxz_yyz_0_xx, g_xxz_yyz_0_xy, g_xxz_yyz_0_xz, g_xxz_yyz_0_yy, g_xxz_yyz_0_yz, g_xxz_yyz_0_zz, g_z_z_0_0_xx_yy_0_xx, g_z_z_0_0_xx_yy_0_xy, g_z_z_0_0_xx_yy_0_xz, g_z_z_0_0_xx_yy_0_yy, g_z_z_0_0_xx_yy_0_yz, g_z_z_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yy_0_xx[i] = 4.0 * g_xxz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_xy[i] = 4.0 * g_xxz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_xz[i] = 4.0 * g_xxz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_yy[i] = 4.0 * g_xxz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_yz[i] = 4.0 * g_xxz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_zz[i] = 4.0 * g_xxz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_xxz_y_0_xx, g_xxz_y_0_xy, g_xxz_y_0_xz, g_xxz_y_0_yy, g_xxz_y_0_yz, g_xxz_y_0_zz, g_xxz_yzz_0_xx, g_xxz_yzz_0_xy, g_xxz_yzz_0_xz, g_xxz_yzz_0_yy, g_xxz_yzz_0_yz, g_xxz_yzz_0_zz, g_z_z_0_0_xx_yz_0_xx, g_z_z_0_0_xx_yz_0_xy, g_z_z_0_0_xx_yz_0_xz, g_z_z_0_0_xx_yz_0_yy, g_z_z_0_0_xx_yz_0_yz, g_z_z_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yz_0_xx[i] = -2.0 * g_xxz_y_0_xx[i] * a_exp + 4.0 * g_xxz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_xy[i] = -2.0 * g_xxz_y_0_xy[i] * a_exp + 4.0 * g_xxz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_xz[i] = -2.0 * g_xxz_y_0_xz[i] * a_exp + 4.0 * g_xxz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_yy[i] = -2.0 * g_xxz_y_0_yy[i] * a_exp + 4.0 * g_xxz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_yz[i] = -2.0 * g_xxz_y_0_yz[i] * a_exp + 4.0 * g_xxz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_zz[i] = -2.0 * g_xxz_y_0_zz[i] * a_exp + 4.0 * g_xxz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_xxz_z_0_xx, g_xxz_z_0_xy, g_xxz_z_0_xz, g_xxz_z_0_yy, g_xxz_z_0_yz, g_xxz_z_0_zz, g_xxz_zzz_0_xx, g_xxz_zzz_0_xy, g_xxz_zzz_0_xz, g_xxz_zzz_0_yy, g_xxz_zzz_0_yz, g_xxz_zzz_0_zz, g_z_z_0_0_xx_zz_0_xx, g_z_z_0_0_xx_zz_0_xy, g_z_z_0_0_xx_zz_0_xz, g_z_z_0_0_xx_zz_0_yy, g_z_z_0_0_xx_zz_0_yz, g_z_z_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_zz_0_xx[i] = -4.0 * g_xxz_z_0_xx[i] * a_exp + 4.0 * g_xxz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_xy[i] = -4.0 * g_xxz_z_0_xy[i] * a_exp + 4.0 * g_xxz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_xz[i] = -4.0 * g_xxz_z_0_xz[i] * a_exp + 4.0 * g_xxz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_yy[i] = -4.0 * g_xxz_z_0_yy[i] * a_exp + 4.0 * g_xxz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_yz[i] = -4.0 * g_xxz_z_0_yz[i] * a_exp + 4.0 * g_xxz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_zz[i] = -4.0 * g_xxz_z_0_zz[i] * a_exp + 4.0 * g_xxz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_xyz_xxz_0_xx, g_xyz_xxz_0_xy, g_xyz_xxz_0_xz, g_xyz_xxz_0_yy, g_xyz_xxz_0_yz, g_xyz_xxz_0_zz, g_z_z_0_0_xy_xx_0_xx, g_z_z_0_0_xy_xx_0_xy, g_z_z_0_0_xy_xx_0_xz, g_z_z_0_0_xy_xx_0_yy, g_z_z_0_0_xy_xx_0_yz, g_z_z_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_0_xx[i] = 4.0 * g_xyz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_xy[i] = 4.0 * g_xyz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_xz[i] = 4.0 * g_xyz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_yy[i] = 4.0 * g_xyz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_yz[i] = 4.0 * g_xyz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_zz[i] = 4.0 * g_xyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_xyz_xyz_0_xx, g_xyz_xyz_0_xy, g_xyz_xyz_0_xz, g_xyz_xyz_0_yy, g_xyz_xyz_0_yz, g_xyz_xyz_0_zz, g_z_z_0_0_xy_xy_0_xx, g_z_z_0_0_xy_xy_0_xy, g_z_z_0_0_xy_xy_0_xz, g_z_z_0_0_xy_xy_0_yy, g_z_z_0_0_xy_xy_0_yz, g_z_z_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xy_0_xx[i] = 4.0 * g_xyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_xy[i] = 4.0 * g_xyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_xz[i] = 4.0 * g_xyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_yy[i] = 4.0 * g_xyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_yz[i] = 4.0 * g_xyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_zz[i] = 4.0 * g_xyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_xyz_x_0_xx, g_xyz_x_0_xy, g_xyz_x_0_xz, g_xyz_x_0_yy, g_xyz_x_0_yz, g_xyz_x_0_zz, g_xyz_xzz_0_xx, g_xyz_xzz_0_xy, g_xyz_xzz_0_xz, g_xyz_xzz_0_yy, g_xyz_xzz_0_yz, g_xyz_xzz_0_zz, g_z_z_0_0_xy_xz_0_xx, g_z_z_0_0_xy_xz_0_xy, g_z_z_0_0_xy_xz_0_xz, g_z_z_0_0_xy_xz_0_yy, g_z_z_0_0_xy_xz_0_yz, g_z_z_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xz_0_xx[i] = -2.0 * g_xyz_x_0_xx[i] * a_exp + 4.0 * g_xyz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_xy[i] = -2.0 * g_xyz_x_0_xy[i] * a_exp + 4.0 * g_xyz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_xz[i] = -2.0 * g_xyz_x_0_xz[i] * a_exp + 4.0 * g_xyz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_yy[i] = -2.0 * g_xyz_x_0_yy[i] * a_exp + 4.0 * g_xyz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_yz[i] = -2.0 * g_xyz_x_0_yz[i] * a_exp + 4.0 * g_xyz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_zz[i] = -2.0 * g_xyz_x_0_zz[i] * a_exp + 4.0 * g_xyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_xyz_yyz_0_xx, g_xyz_yyz_0_xy, g_xyz_yyz_0_xz, g_xyz_yyz_0_yy, g_xyz_yyz_0_yz, g_xyz_yyz_0_zz, g_z_z_0_0_xy_yy_0_xx, g_z_z_0_0_xy_yy_0_xy, g_z_z_0_0_xy_yy_0_xz, g_z_z_0_0_xy_yy_0_yy, g_z_z_0_0_xy_yy_0_yz, g_z_z_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yy_0_xx[i] = 4.0 * g_xyz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_xy[i] = 4.0 * g_xyz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_xz[i] = 4.0 * g_xyz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_yy[i] = 4.0 * g_xyz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_yz[i] = 4.0 * g_xyz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_zz[i] = 4.0 * g_xyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_xyz_y_0_xx, g_xyz_y_0_xy, g_xyz_y_0_xz, g_xyz_y_0_yy, g_xyz_y_0_yz, g_xyz_y_0_zz, g_xyz_yzz_0_xx, g_xyz_yzz_0_xy, g_xyz_yzz_0_xz, g_xyz_yzz_0_yy, g_xyz_yzz_0_yz, g_xyz_yzz_0_zz, g_z_z_0_0_xy_yz_0_xx, g_z_z_0_0_xy_yz_0_xy, g_z_z_0_0_xy_yz_0_xz, g_z_z_0_0_xy_yz_0_yy, g_z_z_0_0_xy_yz_0_yz, g_z_z_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yz_0_xx[i] = -2.0 * g_xyz_y_0_xx[i] * a_exp + 4.0 * g_xyz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_xy[i] = -2.0 * g_xyz_y_0_xy[i] * a_exp + 4.0 * g_xyz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_xz[i] = -2.0 * g_xyz_y_0_xz[i] * a_exp + 4.0 * g_xyz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_yy[i] = -2.0 * g_xyz_y_0_yy[i] * a_exp + 4.0 * g_xyz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_yz[i] = -2.0 * g_xyz_y_0_yz[i] * a_exp + 4.0 * g_xyz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_zz[i] = -2.0 * g_xyz_y_0_zz[i] * a_exp + 4.0 * g_xyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_xyz_z_0_xx, g_xyz_z_0_xy, g_xyz_z_0_xz, g_xyz_z_0_yy, g_xyz_z_0_yz, g_xyz_z_0_zz, g_xyz_zzz_0_xx, g_xyz_zzz_0_xy, g_xyz_zzz_0_xz, g_xyz_zzz_0_yy, g_xyz_zzz_0_yz, g_xyz_zzz_0_zz, g_z_z_0_0_xy_zz_0_xx, g_z_z_0_0_xy_zz_0_xy, g_z_z_0_0_xy_zz_0_xz, g_z_z_0_0_xy_zz_0_yy, g_z_z_0_0_xy_zz_0_yz, g_z_z_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_zz_0_xx[i] = -4.0 * g_xyz_z_0_xx[i] * a_exp + 4.0 * g_xyz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_xy[i] = -4.0 * g_xyz_z_0_xy[i] * a_exp + 4.0 * g_xyz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_xz[i] = -4.0 * g_xyz_z_0_xz[i] * a_exp + 4.0 * g_xyz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_yy[i] = -4.0 * g_xyz_z_0_yy[i] * a_exp + 4.0 * g_xyz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_yz[i] = -4.0 * g_xyz_z_0_yz[i] * a_exp + 4.0 * g_xyz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_zz[i] = -4.0 * g_xyz_z_0_zz[i] * a_exp + 4.0 * g_xyz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_x_xxz_0_xx, g_x_xxz_0_xy, g_x_xxz_0_xz, g_x_xxz_0_yy, g_x_xxz_0_yz, g_x_xxz_0_zz, g_xzz_xxz_0_xx, g_xzz_xxz_0_xy, g_xzz_xxz_0_xz, g_xzz_xxz_0_yy, g_xzz_xxz_0_yz, g_xzz_xxz_0_zz, g_z_z_0_0_xz_xx_0_xx, g_z_z_0_0_xz_xx_0_xy, g_z_z_0_0_xz_xx_0_xz, g_z_z_0_0_xz_xx_0_yy, g_z_z_0_0_xz_xx_0_yz, g_z_z_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_0_xx[i] = -2.0 * g_x_xxz_0_xx[i] * b_exp + 4.0 * g_xzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_xy[i] = -2.0 * g_x_xxz_0_xy[i] * b_exp + 4.0 * g_xzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_xz[i] = -2.0 * g_x_xxz_0_xz[i] * b_exp + 4.0 * g_xzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_yy[i] = -2.0 * g_x_xxz_0_yy[i] * b_exp + 4.0 * g_xzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_yz[i] = -2.0 * g_x_xxz_0_yz[i] * b_exp + 4.0 * g_xzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_zz[i] = -2.0 * g_x_xxz_0_zz[i] * b_exp + 4.0 * g_xzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_x_xyz_0_xx, g_x_xyz_0_xy, g_x_xyz_0_xz, g_x_xyz_0_yy, g_x_xyz_0_yz, g_x_xyz_0_zz, g_xzz_xyz_0_xx, g_xzz_xyz_0_xy, g_xzz_xyz_0_xz, g_xzz_xyz_0_yy, g_xzz_xyz_0_yz, g_xzz_xyz_0_zz, g_z_z_0_0_xz_xy_0_xx, g_z_z_0_0_xz_xy_0_xy, g_z_z_0_0_xz_xy_0_xz, g_z_z_0_0_xz_xy_0_yy, g_z_z_0_0_xz_xy_0_yz, g_z_z_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xy_0_xx[i] = -2.0 * g_x_xyz_0_xx[i] * b_exp + 4.0 * g_xzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_xy[i] = -2.0 * g_x_xyz_0_xy[i] * b_exp + 4.0 * g_xzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_xz[i] = -2.0 * g_x_xyz_0_xz[i] * b_exp + 4.0 * g_xzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_yy[i] = -2.0 * g_x_xyz_0_yy[i] * b_exp + 4.0 * g_xzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_yz[i] = -2.0 * g_x_xyz_0_yz[i] * b_exp + 4.0 * g_xzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_zz[i] = -2.0 * g_x_xyz_0_zz[i] * b_exp + 4.0 * g_xzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_xzz_0_xx, g_x_xzz_0_xy, g_x_xzz_0_xz, g_x_xzz_0_yy, g_x_xzz_0_yz, g_x_xzz_0_zz, g_xzz_x_0_xx, g_xzz_x_0_xy, g_xzz_x_0_xz, g_xzz_x_0_yy, g_xzz_x_0_yz, g_xzz_x_0_zz, g_xzz_xzz_0_xx, g_xzz_xzz_0_xy, g_xzz_xzz_0_xz, g_xzz_xzz_0_yy, g_xzz_xzz_0_yz, g_xzz_xzz_0_zz, g_z_z_0_0_xz_xz_0_xx, g_z_z_0_0_xz_xz_0_xy, g_z_z_0_0_xz_xz_0_xz, g_z_z_0_0_xz_xz_0_yy, g_z_z_0_0_xz_xz_0_yz, g_z_z_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xz_0_xx[i] = g_x_x_0_xx[i] - 2.0 * g_x_xzz_0_xx[i] * b_exp - 2.0 * g_xzz_x_0_xx[i] * a_exp + 4.0 * g_xzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_xy[i] = g_x_x_0_xy[i] - 2.0 * g_x_xzz_0_xy[i] * b_exp - 2.0 * g_xzz_x_0_xy[i] * a_exp + 4.0 * g_xzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_xz[i] = g_x_x_0_xz[i] - 2.0 * g_x_xzz_0_xz[i] * b_exp - 2.0 * g_xzz_x_0_xz[i] * a_exp + 4.0 * g_xzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_yy[i] = g_x_x_0_yy[i] - 2.0 * g_x_xzz_0_yy[i] * b_exp - 2.0 * g_xzz_x_0_yy[i] * a_exp + 4.0 * g_xzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_yz[i] = g_x_x_0_yz[i] - 2.0 * g_x_xzz_0_yz[i] * b_exp - 2.0 * g_xzz_x_0_yz[i] * a_exp + 4.0 * g_xzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_zz[i] = g_x_x_0_zz[i] - 2.0 * g_x_xzz_0_zz[i] * b_exp - 2.0 * g_xzz_x_0_zz[i] * a_exp + 4.0 * g_xzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_x_yyz_0_xx, g_x_yyz_0_xy, g_x_yyz_0_xz, g_x_yyz_0_yy, g_x_yyz_0_yz, g_x_yyz_0_zz, g_xzz_yyz_0_xx, g_xzz_yyz_0_xy, g_xzz_yyz_0_xz, g_xzz_yyz_0_yy, g_xzz_yyz_0_yz, g_xzz_yyz_0_zz, g_z_z_0_0_xz_yy_0_xx, g_z_z_0_0_xz_yy_0_xy, g_z_z_0_0_xz_yy_0_xz, g_z_z_0_0_xz_yy_0_yy, g_z_z_0_0_xz_yy_0_yz, g_z_z_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yy_0_xx[i] = -2.0 * g_x_yyz_0_xx[i] * b_exp + 4.0 * g_xzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_xy[i] = -2.0 * g_x_yyz_0_xy[i] * b_exp + 4.0 * g_xzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_xz[i] = -2.0 * g_x_yyz_0_xz[i] * b_exp + 4.0 * g_xzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_yy[i] = -2.0 * g_x_yyz_0_yy[i] * b_exp + 4.0 * g_xzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_yz[i] = -2.0 * g_x_yyz_0_yz[i] * b_exp + 4.0 * g_xzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_zz[i] = -2.0 * g_x_yyz_0_zz[i] * b_exp + 4.0 * g_xzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_yzz_0_xx, g_x_yzz_0_xy, g_x_yzz_0_xz, g_x_yzz_0_yy, g_x_yzz_0_yz, g_x_yzz_0_zz, g_xzz_y_0_xx, g_xzz_y_0_xy, g_xzz_y_0_xz, g_xzz_y_0_yy, g_xzz_y_0_yz, g_xzz_y_0_zz, g_xzz_yzz_0_xx, g_xzz_yzz_0_xy, g_xzz_yzz_0_xz, g_xzz_yzz_0_yy, g_xzz_yzz_0_yz, g_xzz_yzz_0_zz, g_z_z_0_0_xz_yz_0_xx, g_z_z_0_0_xz_yz_0_xy, g_z_z_0_0_xz_yz_0_xz, g_z_z_0_0_xz_yz_0_yy, g_z_z_0_0_xz_yz_0_yz, g_z_z_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yz_0_xx[i] = g_x_y_0_xx[i] - 2.0 * g_x_yzz_0_xx[i] * b_exp - 2.0 * g_xzz_y_0_xx[i] * a_exp + 4.0 * g_xzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_xy[i] = g_x_y_0_xy[i] - 2.0 * g_x_yzz_0_xy[i] * b_exp - 2.0 * g_xzz_y_0_xy[i] * a_exp + 4.0 * g_xzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_xz[i] = g_x_y_0_xz[i] - 2.0 * g_x_yzz_0_xz[i] * b_exp - 2.0 * g_xzz_y_0_xz[i] * a_exp + 4.0 * g_xzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_yy[i] = g_x_y_0_yy[i] - 2.0 * g_x_yzz_0_yy[i] * b_exp - 2.0 * g_xzz_y_0_yy[i] * a_exp + 4.0 * g_xzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_yz[i] = g_x_y_0_yz[i] - 2.0 * g_x_yzz_0_yz[i] * b_exp - 2.0 * g_xzz_y_0_yz[i] * a_exp + 4.0 * g_xzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_zz[i] = g_x_y_0_zz[i] - 2.0 * g_x_yzz_0_zz[i] * b_exp - 2.0 * g_xzz_y_0_zz[i] * a_exp + 4.0 * g_xzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_zzz_0_xx, g_x_zzz_0_xy, g_x_zzz_0_xz, g_x_zzz_0_yy, g_x_zzz_0_yz, g_x_zzz_0_zz, g_xzz_z_0_xx, g_xzz_z_0_xy, g_xzz_z_0_xz, g_xzz_z_0_yy, g_xzz_z_0_yz, g_xzz_z_0_zz, g_xzz_zzz_0_xx, g_xzz_zzz_0_xy, g_xzz_zzz_0_xz, g_xzz_zzz_0_yy, g_xzz_zzz_0_yz, g_xzz_zzz_0_zz, g_z_z_0_0_xz_zz_0_xx, g_z_z_0_0_xz_zz_0_xy, g_z_z_0_0_xz_zz_0_xz, g_z_z_0_0_xz_zz_0_yy, g_z_z_0_0_xz_zz_0_yz, g_z_z_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_zz_0_xx[i] = 2.0 * g_x_z_0_xx[i] - 2.0 * g_x_zzz_0_xx[i] * b_exp - 4.0 * g_xzz_z_0_xx[i] * a_exp + 4.0 * g_xzz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_xy[i] = 2.0 * g_x_z_0_xy[i] - 2.0 * g_x_zzz_0_xy[i] * b_exp - 4.0 * g_xzz_z_0_xy[i] * a_exp + 4.0 * g_xzz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_xz[i] = 2.0 * g_x_z_0_xz[i] - 2.0 * g_x_zzz_0_xz[i] * b_exp - 4.0 * g_xzz_z_0_xz[i] * a_exp + 4.0 * g_xzz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_yy[i] = 2.0 * g_x_z_0_yy[i] - 2.0 * g_x_zzz_0_yy[i] * b_exp - 4.0 * g_xzz_z_0_yy[i] * a_exp + 4.0 * g_xzz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_yz[i] = 2.0 * g_x_z_0_yz[i] - 2.0 * g_x_zzz_0_yz[i] * b_exp - 4.0 * g_xzz_z_0_yz[i] * a_exp + 4.0 * g_xzz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_zz[i] = 2.0 * g_x_z_0_zz[i] - 2.0 * g_x_zzz_0_zz[i] * b_exp - 4.0 * g_xzz_z_0_zz[i] * a_exp + 4.0 * g_xzz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_yyz_xxz_0_xx, g_yyz_xxz_0_xy, g_yyz_xxz_0_xz, g_yyz_xxz_0_yy, g_yyz_xxz_0_yz, g_yyz_xxz_0_zz, g_z_z_0_0_yy_xx_0_xx, g_z_z_0_0_yy_xx_0_xy, g_z_z_0_0_yy_xx_0_xz, g_z_z_0_0_yy_xx_0_yy, g_z_z_0_0_yy_xx_0_yz, g_z_z_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_0_xx[i] = 4.0 * g_yyz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_xy[i] = 4.0 * g_yyz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_xz[i] = 4.0 * g_yyz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_yy[i] = 4.0 * g_yyz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_yz[i] = 4.0 * g_yyz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_zz[i] = 4.0 * g_yyz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_yyz_xyz_0_xx, g_yyz_xyz_0_xy, g_yyz_xyz_0_xz, g_yyz_xyz_0_yy, g_yyz_xyz_0_yz, g_yyz_xyz_0_zz, g_z_z_0_0_yy_xy_0_xx, g_z_z_0_0_yy_xy_0_xy, g_z_z_0_0_yy_xy_0_xz, g_z_z_0_0_yy_xy_0_yy, g_z_z_0_0_yy_xy_0_yz, g_z_z_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xy_0_xx[i] = 4.0 * g_yyz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_xy[i] = 4.0 * g_yyz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_xz[i] = 4.0 * g_yyz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_yy[i] = 4.0 * g_yyz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_yz[i] = 4.0 * g_yyz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_zz[i] = 4.0 * g_yyz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_yyz_x_0_xx, g_yyz_x_0_xy, g_yyz_x_0_xz, g_yyz_x_0_yy, g_yyz_x_0_yz, g_yyz_x_0_zz, g_yyz_xzz_0_xx, g_yyz_xzz_0_xy, g_yyz_xzz_0_xz, g_yyz_xzz_0_yy, g_yyz_xzz_0_yz, g_yyz_xzz_0_zz, g_z_z_0_0_yy_xz_0_xx, g_z_z_0_0_yy_xz_0_xy, g_z_z_0_0_yy_xz_0_xz, g_z_z_0_0_yy_xz_0_yy, g_z_z_0_0_yy_xz_0_yz, g_z_z_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xz_0_xx[i] = -2.0 * g_yyz_x_0_xx[i] * a_exp + 4.0 * g_yyz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_xy[i] = -2.0 * g_yyz_x_0_xy[i] * a_exp + 4.0 * g_yyz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_xz[i] = -2.0 * g_yyz_x_0_xz[i] * a_exp + 4.0 * g_yyz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_yy[i] = -2.0 * g_yyz_x_0_yy[i] * a_exp + 4.0 * g_yyz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_yz[i] = -2.0 * g_yyz_x_0_yz[i] * a_exp + 4.0 * g_yyz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_zz[i] = -2.0 * g_yyz_x_0_zz[i] * a_exp + 4.0 * g_yyz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_yyz_yyz_0_xx, g_yyz_yyz_0_xy, g_yyz_yyz_0_xz, g_yyz_yyz_0_yy, g_yyz_yyz_0_yz, g_yyz_yyz_0_zz, g_z_z_0_0_yy_yy_0_xx, g_z_z_0_0_yy_yy_0_xy, g_z_z_0_0_yy_yy_0_xz, g_z_z_0_0_yy_yy_0_yy, g_z_z_0_0_yy_yy_0_yz, g_z_z_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yy_0_xx[i] = 4.0 * g_yyz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_xy[i] = 4.0 * g_yyz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_xz[i] = 4.0 * g_yyz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_yy[i] = 4.0 * g_yyz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_yz[i] = 4.0 * g_yyz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_zz[i] = 4.0 * g_yyz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_yyz_y_0_xx, g_yyz_y_0_xy, g_yyz_y_0_xz, g_yyz_y_0_yy, g_yyz_y_0_yz, g_yyz_y_0_zz, g_yyz_yzz_0_xx, g_yyz_yzz_0_xy, g_yyz_yzz_0_xz, g_yyz_yzz_0_yy, g_yyz_yzz_0_yz, g_yyz_yzz_0_zz, g_z_z_0_0_yy_yz_0_xx, g_z_z_0_0_yy_yz_0_xy, g_z_z_0_0_yy_yz_0_xz, g_z_z_0_0_yy_yz_0_yy, g_z_z_0_0_yy_yz_0_yz, g_z_z_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yz_0_xx[i] = -2.0 * g_yyz_y_0_xx[i] * a_exp + 4.0 * g_yyz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_xy[i] = -2.0 * g_yyz_y_0_xy[i] * a_exp + 4.0 * g_yyz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_xz[i] = -2.0 * g_yyz_y_0_xz[i] * a_exp + 4.0 * g_yyz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_yy[i] = -2.0 * g_yyz_y_0_yy[i] * a_exp + 4.0 * g_yyz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_yz[i] = -2.0 * g_yyz_y_0_yz[i] * a_exp + 4.0 * g_yyz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_zz[i] = -2.0 * g_yyz_y_0_zz[i] * a_exp + 4.0 * g_yyz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_yyz_z_0_xx, g_yyz_z_0_xy, g_yyz_z_0_xz, g_yyz_z_0_yy, g_yyz_z_0_yz, g_yyz_z_0_zz, g_yyz_zzz_0_xx, g_yyz_zzz_0_xy, g_yyz_zzz_0_xz, g_yyz_zzz_0_yy, g_yyz_zzz_0_yz, g_yyz_zzz_0_zz, g_z_z_0_0_yy_zz_0_xx, g_z_z_0_0_yy_zz_0_xy, g_z_z_0_0_yy_zz_0_xz, g_z_z_0_0_yy_zz_0_yy, g_z_z_0_0_yy_zz_0_yz, g_z_z_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_zz_0_xx[i] = -4.0 * g_yyz_z_0_xx[i] * a_exp + 4.0 * g_yyz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_xy[i] = -4.0 * g_yyz_z_0_xy[i] * a_exp + 4.0 * g_yyz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_xz[i] = -4.0 * g_yyz_z_0_xz[i] * a_exp + 4.0 * g_yyz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_yy[i] = -4.0 * g_yyz_z_0_yy[i] * a_exp + 4.0 * g_yyz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_yz[i] = -4.0 * g_yyz_z_0_yz[i] * a_exp + 4.0 * g_yyz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_zz[i] = -4.0 * g_yyz_z_0_zz[i] * a_exp + 4.0 * g_yyz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_y_xxz_0_xx, g_y_xxz_0_xy, g_y_xxz_0_xz, g_y_xxz_0_yy, g_y_xxz_0_yz, g_y_xxz_0_zz, g_yzz_xxz_0_xx, g_yzz_xxz_0_xy, g_yzz_xxz_0_xz, g_yzz_xxz_0_yy, g_yzz_xxz_0_yz, g_yzz_xxz_0_zz, g_z_z_0_0_yz_xx_0_xx, g_z_z_0_0_yz_xx_0_xy, g_z_z_0_0_yz_xx_0_xz, g_z_z_0_0_yz_xx_0_yy, g_z_z_0_0_yz_xx_0_yz, g_z_z_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_0_xx[i] = -2.0 * g_y_xxz_0_xx[i] * b_exp + 4.0 * g_yzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_xy[i] = -2.0 * g_y_xxz_0_xy[i] * b_exp + 4.0 * g_yzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_xz[i] = -2.0 * g_y_xxz_0_xz[i] * b_exp + 4.0 * g_yzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_yy[i] = -2.0 * g_y_xxz_0_yy[i] * b_exp + 4.0 * g_yzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_yz[i] = -2.0 * g_y_xxz_0_yz[i] * b_exp + 4.0 * g_yzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_zz[i] = -2.0 * g_y_xxz_0_zz[i] * b_exp + 4.0 * g_yzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_y_xyz_0_xx, g_y_xyz_0_xy, g_y_xyz_0_xz, g_y_xyz_0_yy, g_y_xyz_0_yz, g_y_xyz_0_zz, g_yzz_xyz_0_xx, g_yzz_xyz_0_xy, g_yzz_xyz_0_xz, g_yzz_xyz_0_yy, g_yzz_xyz_0_yz, g_yzz_xyz_0_zz, g_z_z_0_0_yz_xy_0_xx, g_z_z_0_0_yz_xy_0_xy, g_z_z_0_0_yz_xy_0_xz, g_z_z_0_0_yz_xy_0_yy, g_z_z_0_0_yz_xy_0_yz, g_z_z_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xy_0_xx[i] = -2.0 * g_y_xyz_0_xx[i] * b_exp + 4.0 * g_yzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_xy[i] = -2.0 * g_y_xyz_0_xy[i] * b_exp + 4.0 * g_yzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_xz[i] = -2.0 * g_y_xyz_0_xz[i] * b_exp + 4.0 * g_yzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_yy[i] = -2.0 * g_y_xyz_0_yy[i] * b_exp + 4.0 * g_yzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_yz[i] = -2.0 * g_y_xyz_0_yz[i] * b_exp + 4.0 * g_yzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_zz[i] = -2.0 * g_y_xyz_0_zz[i] * b_exp + 4.0 * g_yzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_xzz_0_xx, g_y_xzz_0_xy, g_y_xzz_0_xz, g_y_xzz_0_yy, g_y_xzz_0_yz, g_y_xzz_0_zz, g_yzz_x_0_xx, g_yzz_x_0_xy, g_yzz_x_0_xz, g_yzz_x_0_yy, g_yzz_x_0_yz, g_yzz_x_0_zz, g_yzz_xzz_0_xx, g_yzz_xzz_0_xy, g_yzz_xzz_0_xz, g_yzz_xzz_0_yy, g_yzz_xzz_0_yz, g_yzz_xzz_0_zz, g_z_z_0_0_yz_xz_0_xx, g_z_z_0_0_yz_xz_0_xy, g_z_z_0_0_yz_xz_0_xz, g_z_z_0_0_yz_xz_0_yy, g_z_z_0_0_yz_xz_0_yz, g_z_z_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xz_0_xx[i] = g_y_x_0_xx[i] - 2.0 * g_y_xzz_0_xx[i] * b_exp - 2.0 * g_yzz_x_0_xx[i] * a_exp + 4.0 * g_yzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_xy[i] = g_y_x_0_xy[i] - 2.0 * g_y_xzz_0_xy[i] * b_exp - 2.0 * g_yzz_x_0_xy[i] * a_exp + 4.0 * g_yzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_xz[i] = g_y_x_0_xz[i] - 2.0 * g_y_xzz_0_xz[i] * b_exp - 2.0 * g_yzz_x_0_xz[i] * a_exp + 4.0 * g_yzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_yy[i] = g_y_x_0_yy[i] - 2.0 * g_y_xzz_0_yy[i] * b_exp - 2.0 * g_yzz_x_0_yy[i] * a_exp + 4.0 * g_yzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_yz[i] = g_y_x_0_yz[i] - 2.0 * g_y_xzz_0_yz[i] * b_exp - 2.0 * g_yzz_x_0_yz[i] * a_exp + 4.0 * g_yzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_zz[i] = g_y_x_0_zz[i] - 2.0 * g_y_xzz_0_zz[i] * b_exp - 2.0 * g_yzz_x_0_zz[i] * a_exp + 4.0 * g_yzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_y_yyz_0_xx, g_y_yyz_0_xy, g_y_yyz_0_xz, g_y_yyz_0_yy, g_y_yyz_0_yz, g_y_yyz_0_zz, g_yzz_yyz_0_xx, g_yzz_yyz_0_xy, g_yzz_yyz_0_xz, g_yzz_yyz_0_yy, g_yzz_yyz_0_yz, g_yzz_yyz_0_zz, g_z_z_0_0_yz_yy_0_xx, g_z_z_0_0_yz_yy_0_xy, g_z_z_0_0_yz_yy_0_xz, g_z_z_0_0_yz_yy_0_yy, g_z_z_0_0_yz_yy_0_yz, g_z_z_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yy_0_xx[i] = -2.0 * g_y_yyz_0_xx[i] * b_exp + 4.0 * g_yzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_xy[i] = -2.0 * g_y_yyz_0_xy[i] * b_exp + 4.0 * g_yzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_xz[i] = -2.0 * g_y_yyz_0_xz[i] * b_exp + 4.0 * g_yzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_yy[i] = -2.0 * g_y_yyz_0_yy[i] * b_exp + 4.0 * g_yzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_yz[i] = -2.0 * g_y_yyz_0_yz[i] * b_exp + 4.0 * g_yzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_zz[i] = -2.0 * g_y_yyz_0_zz[i] * b_exp + 4.0 * g_yzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_yzz_0_xx, g_y_yzz_0_xy, g_y_yzz_0_xz, g_y_yzz_0_yy, g_y_yzz_0_yz, g_y_yzz_0_zz, g_yzz_y_0_xx, g_yzz_y_0_xy, g_yzz_y_0_xz, g_yzz_y_0_yy, g_yzz_y_0_yz, g_yzz_y_0_zz, g_yzz_yzz_0_xx, g_yzz_yzz_0_xy, g_yzz_yzz_0_xz, g_yzz_yzz_0_yy, g_yzz_yzz_0_yz, g_yzz_yzz_0_zz, g_z_z_0_0_yz_yz_0_xx, g_z_z_0_0_yz_yz_0_xy, g_z_z_0_0_yz_yz_0_xz, g_z_z_0_0_yz_yz_0_yy, g_z_z_0_0_yz_yz_0_yz, g_z_z_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yz_0_xx[i] = g_y_y_0_xx[i] - 2.0 * g_y_yzz_0_xx[i] * b_exp - 2.0 * g_yzz_y_0_xx[i] * a_exp + 4.0 * g_yzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_xy[i] = g_y_y_0_xy[i] - 2.0 * g_y_yzz_0_xy[i] * b_exp - 2.0 * g_yzz_y_0_xy[i] * a_exp + 4.0 * g_yzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_xz[i] = g_y_y_0_xz[i] - 2.0 * g_y_yzz_0_xz[i] * b_exp - 2.0 * g_yzz_y_0_xz[i] * a_exp + 4.0 * g_yzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_yy[i] = g_y_y_0_yy[i] - 2.0 * g_y_yzz_0_yy[i] * b_exp - 2.0 * g_yzz_y_0_yy[i] * a_exp + 4.0 * g_yzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_yz[i] = g_y_y_0_yz[i] - 2.0 * g_y_yzz_0_yz[i] * b_exp - 2.0 * g_yzz_y_0_yz[i] * a_exp + 4.0 * g_yzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_zz[i] = g_y_y_0_zz[i] - 2.0 * g_y_yzz_0_zz[i] * b_exp - 2.0 * g_yzz_y_0_zz[i] * a_exp + 4.0 * g_yzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_zzz_0_xx, g_y_zzz_0_xy, g_y_zzz_0_xz, g_y_zzz_0_yy, g_y_zzz_0_yz, g_y_zzz_0_zz, g_yzz_z_0_xx, g_yzz_z_0_xy, g_yzz_z_0_xz, g_yzz_z_0_yy, g_yzz_z_0_yz, g_yzz_z_0_zz, g_yzz_zzz_0_xx, g_yzz_zzz_0_xy, g_yzz_zzz_0_xz, g_yzz_zzz_0_yy, g_yzz_zzz_0_yz, g_yzz_zzz_0_zz, g_z_z_0_0_yz_zz_0_xx, g_z_z_0_0_yz_zz_0_xy, g_z_z_0_0_yz_zz_0_xz, g_z_z_0_0_yz_zz_0_yy, g_z_z_0_0_yz_zz_0_yz, g_z_z_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_zz_0_xx[i] = 2.0 * g_y_z_0_xx[i] - 2.0 * g_y_zzz_0_xx[i] * b_exp - 4.0 * g_yzz_z_0_xx[i] * a_exp + 4.0 * g_yzz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_xy[i] = 2.0 * g_y_z_0_xy[i] - 2.0 * g_y_zzz_0_xy[i] * b_exp - 4.0 * g_yzz_z_0_xy[i] * a_exp + 4.0 * g_yzz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_xz[i] = 2.0 * g_y_z_0_xz[i] - 2.0 * g_y_zzz_0_xz[i] * b_exp - 4.0 * g_yzz_z_0_xz[i] * a_exp + 4.0 * g_yzz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_yy[i] = 2.0 * g_y_z_0_yy[i] - 2.0 * g_y_zzz_0_yy[i] * b_exp - 4.0 * g_yzz_z_0_yy[i] * a_exp + 4.0 * g_yzz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_yz[i] = 2.0 * g_y_z_0_yz[i] - 2.0 * g_y_zzz_0_yz[i] * b_exp - 4.0 * g_yzz_z_0_yz[i] * a_exp + 4.0 * g_yzz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_zz[i] = 2.0 * g_y_z_0_zz[i] - 2.0 * g_y_zzz_0_zz[i] * b_exp - 4.0 * g_yzz_z_0_zz[i] * a_exp + 4.0 * g_yzz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_xxz_0_xx, g_z_xxz_0_xy, g_z_xxz_0_xz, g_z_xxz_0_yy, g_z_xxz_0_yz, g_z_xxz_0_zz, g_z_z_0_0_zz_xx_0_xx, g_z_z_0_0_zz_xx_0_xy, g_z_z_0_0_zz_xx_0_xz, g_z_z_0_0_zz_xx_0_yy, g_z_z_0_0_zz_xx_0_yz, g_z_z_0_0_zz_xx_0_zz, g_zzz_xxz_0_xx, g_zzz_xxz_0_xy, g_zzz_xxz_0_xz, g_zzz_xxz_0_yy, g_zzz_xxz_0_yz, g_zzz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_0_xx[i] = -4.0 * g_z_xxz_0_xx[i] * b_exp + 4.0 * g_zzz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_xy[i] = -4.0 * g_z_xxz_0_xy[i] * b_exp + 4.0 * g_zzz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_xz[i] = -4.0 * g_z_xxz_0_xz[i] * b_exp + 4.0 * g_zzz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_yy[i] = -4.0 * g_z_xxz_0_yy[i] * b_exp + 4.0 * g_zzz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_yz[i] = -4.0 * g_z_xxz_0_yz[i] * b_exp + 4.0 * g_zzz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_zz[i] = -4.0 * g_z_xxz_0_zz[i] * b_exp + 4.0 * g_zzz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_xyz_0_xx, g_z_xyz_0_xy, g_z_xyz_0_xz, g_z_xyz_0_yy, g_z_xyz_0_yz, g_z_xyz_0_zz, g_z_z_0_0_zz_xy_0_xx, g_z_z_0_0_zz_xy_0_xy, g_z_z_0_0_zz_xy_0_xz, g_z_z_0_0_zz_xy_0_yy, g_z_z_0_0_zz_xy_0_yz, g_z_z_0_0_zz_xy_0_zz, g_zzz_xyz_0_xx, g_zzz_xyz_0_xy, g_zzz_xyz_0_xz, g_zzz_xyz_0_yy, g_zzz_xyz_0_yz, g_zzz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xy_0_xx[i] = -4.0 * g_z_xyz_0_xx[i] * b_exp + 4.0 * g_zzz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_xy[i] = -4.0 * g_z_xyz_0_xy[i] * b_exp + 4.0 * g_zzz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_xz[i] = -4.0 * g_z_xyz_0_xz[i] * b_exp + 4.0 * g_zzz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_yy[i] = -4.0 * g_z_xyz_0_yy[i] * b_exp + 4.0 * g_zzz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_yz[i] = -4.0 * g_z_xyz_0_yz[i] * b_exp + 4.0 * g_zzz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_zz[i] = -4.0 * g_z_xyz_0_zz[i] * b_exp + 4.0 * g_zzz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_xzz_0_xx, g_z_xzz_0_xy, g_z_xzz_0_xz, g_z_xzz_0_yy, g_z_xzz_0_yz, g_z_xzz_0_zz, g_z_z_0_0_zz_xz_0_xx, g_z_z_0_0_zz_xz_0_xy, g_z_z_0_0_zz_xz_0_xz, g_z_z_0_0_zz_xz_0_yy, g_z_z_0_0_zz_xz_0_yz, g_z_z_0_0_zz_xz_0_zz, g_zzz_x_0_xx, g_zzz_x_0_xy, g_zzz_x_0_xz, g_zzz_x_0_yy, g_zzz_x_0_yz, g_zzz_x_0_zz, g_zzz_xzz_0_xx, g_zzz_xzz_0_xy, g_zzz_xzz_0_xz, g_zzz_xzz_0_yy, g_zzz_xzz_0_yz, g_zzz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xz_0_xx[i] = 2.0 * g_z_x_0_xx[i] - 4.0 * g_z_xzz_0_xx[i] * b_exp - 2.0 * g_zzz_x_0_xx[i] * a_exp + 4.0 * g_zzz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_xy[i] = 2.0 * g_z_x_0_xy[i] - 4.0 * g_z_xzz_0_xy[i] * b_exp - 2.0 * g_zzz_x_0_xy[i] * a_exp + 4.0 * g_zzz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_xz[i] = 2.0 * g_z_x_0_xz[i] - 4.0 * g_z_xzz_0_xz[i] * b_exp - 2.0 * g_zzz_x_0_xz[i] * a_exp + 4.0 * g_zzz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_yy[i] = 2.0 * g_z_x_0_yy[i] - 4.0 * g_z_xzz_0_yy[i] * b_exp - 2.0 * g_zzz_x_0_yy[i] * a_exp + 4.0 * g_zzz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_yz[i] = 2.0 * g_z_x_0_yz[i] - 4.0 * g_z_xzz_0_yz[i] * b_exp - 2.0 * g_zzz_x_0_yz[i] * a_exp + 4.0 * g_zzz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_zz[i] = 2.0 * g_z_x_0_zz[i] - 4.0 * g_z_xzz_0_zz[i] * b_exp - 2.0 * g_zzz_x_0_zz[i] * a_exp + 4.0 * g_zzz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_yyz_0_xx, g_z_yyz_0_xy, g_z_yyz_0_xz, g_z_yyz_0_yy, g_z_yyz_0_yz, g_z_yyz_0_zz, g_z_z_0_0_zz_yy_0_xx, g_z_z_0_0_zz_yy_0_xy, g_z_z_0_0_zz_yy_0_xz, g_z_z_0_0_zz_yy_0_yy, g_z_z_0_0_zz_yy_0_yz, g_z_z_0_0_zz_yy_0_zz, g_zzz_yyz_0_xx, g_zzz_yyz_0_xy, g_zzz_yyz_0_xz, g_zzz_yyz_0_yy, g_zzz_yyz_0_yz, g_zzz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yy_0_xx[i] = -4.0 * g_z_yyz_0_xx[i] * b_exp + 4.0 * g_zzz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_xy[i] = -4.0 * g_z_yyz_0_xy[i] * b_exp + 4.0 * g_zzz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_xz[i] = -4.0 * g_z_yyz_0_xz[i] * b_exp + 4.0 * g_zzz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_yy[i] = -4.0 * g_z_yyz_0_yy[i] * b_exp + 4.0 * g_zzz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_yz[i] = -4.0 * g_z_yyz_0_yz[i] * b_exp + 4.0 * g_zzz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_zz[i] = -4.0 * g_z_yyz_0_zz[i] * b_exp + 4.0 * g_zzz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_yzz_0_xx, g_z_yzz_0_xy, g_z_yzz_0_xz, g_z_yzz_0_yy, g_z_yzz_0_yz, g_z_yzz_0_zz, g_z_z_0_0_zz_yz_0_xx, g_z_z_0_0_zz_yz_0_xy, g_z_z_0_0_zz_yz_0_xz, g_z_z_0_0_zz_yz_0_yy, g_z_z_0_0_zz_yz_0_yz, g_z_z_0_0_zz_yz_0_zz, g_zzz_y_0_xx, g_zzz_y_0_xy, g_zzz_y_0_xz, g_zzz_y_0_yy, g_zzz_y_0_yz, g_zzz_y_0_zz, g_zzz_yzz_0_xx, g_zzz_yzz_0_xy, g_zzz_yzz_0_xz, g_zzz_yzz_0_yy, g_zzz_yzz_0_yz, g_zzz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yz_0_xx[i] = 2.0 * g_z_y_0_xx[i] - 4.0 * g_z_yzz_0_xx[i] * b_exp - 2.0 * g_zzz_y_0_xx[i] * a_exp + 4.0 * g_zzz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_xy[i] = 2.0 * g_z_y_0_xy[i] - 4.0 * g_z_yzz_0_xy[i] * b_exp - 2.0 * g_zzz_y_0_xy[i] * a_exp + 4.0 * g_zzz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_xz[i] = 2.0 * g_z_y_0_xz[i] - 4.0 * g_z_yzz_0_xz[i] * b_exp - 2.0 * g_zzz_y_0_xz[i] * a_exp + 4.0 * g_zzz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_yy[i] = 2.0 * g_z_y_0_yy[i] - 4.0 * g_z_yzz_0_yy[i] * b_exp - 2.0 * g_zzz_y_0_yy[i] * a_exp + 4.0 * g_zzz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_yz[i] = 2.0 * g_z_y_0_yz[i] - 4.0 * g_z_yzz_0_yz[i] * b_exp - 2.0 * g_zzz_y_0_yz[i] * a_exp + 4.0 * g_zzz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_zz[i] = 2.0 * g_z_y_0_zz[i] - 4.0 * g_z_yzz_0_zz[i] * b_exp - 2.0 * g_zzz_y_0_zz[i] * a_exp + 4.0 * g_zzz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_z_0_0_zz_zz_0_xx, g_z_z_0_0_zz_zz_0_xy, g_z_z_0_0_zz_zz_0_xz, g_z_z_0_0_zz_zz_0_yy, g_z_z_0_0_zz_zz_0_yz, g_z_z_0_0_zz_zz_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_zzz_0_xx, g_z_zzz_0_xy, g_z_zzz_0_xz, g_z_zzz_0_yy, g_z_zzz_0_yz, g_z_zzz_0_zz, g_zzz_z_0_xx, g_zzz_z_0_xy, g_zzz_z_0_xz, g_zzz_z_0_yy, g_zzz_z_0_yz, g_zzz_z_0_zz, g_zzz_zzz_0_xx, g_zzz_zzz_0_xy, g_zzz_zzz_0_xz, g_zzz_zzz_0_yy, g_zzz_zzz_0_yz, g_zzz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_zz_0_xx[i] = 4.0 * g_z_z_0_xx[i] - 4.0 * g_z_zzz_0_xx[i] * b_exp - 4.0 * g_zzz_z_0_xx[i] * a_exp + 4.0 * g_zzz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_xy[i] = 4.0 * g_z_z_0_xy[i] - 4.0 * g_z_zzz_0_xy[i] * b_exp - 4.0 * g_zzz_z_0_xy[i] * a_exp + 4.0 * g_zzz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_xz[i] = 4.0 * g_z_z_0_xz[i] - 4.0 * g_z_zzz_0_xz[i] * b_exp - 4.0 * g_zzz_z_0_xz[i] * a_exp + 4.0 * g_zzz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_yy[i] = 4.0 * g_z_z_0_yy[i] - 4.0 * g_z_zzz_0_yy[i] * b_exp - 4.0 * g_zzz_z_0_yy[i] * a_exp + 4.0 * g_zzz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_yz[i] = 4.0 * g_z_z_0_yz[i] - 4.0 * g_z_zzz_0_yz[i] * b_exp - 4.0 * g_zzz_z_0_yz[i] * a_exp + 4.0 * g_zzz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_zz[i] = 4.0 * g_z_z_0_zz[i] - 4.0 * g_z_zzz_0_zz[i] * b_exp - 4.0 * g_zzz_z_0_zz[i] * a_exp + 4.0 * g_zzz_zzz_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


#include "GeomDeriv1100OfScalarForPDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pdsd_0(CSimdArray<double>& buffer_1100_pdsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_sfsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const CSimdArray<double>& buffer_dfsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pdsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spsd

    auto g_0_x_0_xx = buffer_spsd[0];

    auto g_0_x_0_xy = buffer_spsd[1];

    auto g_0_x_0_xz = buffer_spsd[2];

    auto g_0_x_0_yy = buffer_spsd[3];

    auto g_0_x_0_yz = buffer_spsd[4];

    auto g_0_x_0_zz = buffer_spsd[5];

    auto g_0_y_0_xx = buffer_spsd[6];

    auto g_0_y_0_xy = buffer_spsd[7];

    auto g_0_y_0_xz = buffer_spsd[8];

    auto g_0_y_0_yy = buffer_spsd[9];

    auto g_0_y_0_yz = buffer_spsd[10];

    auto g_0_y_0_zz = buffer_spsd[11];

    auto g_0_z_0_xx = buffer_spsd[12];

    auto g_0_z_0_xy = buffer_spsd[13];

    auto g_0_z_0_xz = buffer_spsd[14];

    auto g_0_z_0_yy = buffer_spsd[15];

    auto g_0_z_0_yz = buffer_spsd[16];

    auto g_0_z_0_zz = buffer_spsd[17];

    /// Set up components of auxilary buffer : buffer_sfsd

    auto g_0_xxx_0_xx = buffer_sfsd[0];

    auto g_0_xxx_0_xy = buffer_sfsd[1];

    auto g_0_xxx_0_xz = buffer_sfsd[2];

    auto g_0_xxx_0_yy = buffer_sfsd[3];

    auto g_0_xxx_0_yz = buffer_sfsd[4];

    auto g_0_xxx_0_zz = buffer_sfsd[5];

    auto g_0_xxy_0_xx = buffer_sfsd[6];

    auto g_0_xxy_0_xy = buffer_sfsd[7];

    auto g_0_xxy_0_xz = buffer_sfsd[8];

    auto g_0_xxy_0_yy = buffer_sfsd[9];

    auto g_0_xxy_0_yz = buffer_sfsd[10];

    auto g_0_xxy_0_zz = buffer_sfsd[11];

    auto g_0_xxz_0_xx = buffer_sfsd[12];

    auto g_0_xxz_0_xy = buffer_sfsd[13];

    auto g_0_xxz_0_xz = buffer_sfsd[14];

    auto g_0_xxz_0_yy = buffer_sfsd[15];

    auto g_0_xxz_0_yz = buffer_sfsd[16];

    auto g_0_xxz_0_zz = buffer_sfsd[17];

    auto g_0_xyy_0_xx = buffer_sfsd[18];

    auto g_0_xyy_0_xy = buffer_sfsd[19];

    auto g_0_xyy_0_xz = buffer_sfsd[20];

    auto g_0_xyy_0_yy = buffer_sfsd[21];

    auto g_0_xyy_0_yz = buffer_sfsd[22];

    auto g_0_xyy_0_zz = buffer_sfsd[23];

    auto g_0_xyz_0_xx = buffer_sfsd[24];

    auto g_0_xyz_0_xy = buffer_sfsd[25];

    auto g_0_xyz_0_xz = buffer_sfsd[26];

    auto g_0_xyz_0_yy = buffer_sfsd[27];

    auto g_0_xyz_0_yz = buffer_sfsd[28];

    auto g_0_xyz_0_zz = buffer_sfsd[29];

    auto g_0_xzz_0_xx = buffer_sfsd[30];

    auto g_0_xzz_0_xy = buffer_sfsd[31];

    auto g_0_xzz_0_xz = buffer_sfsd[32];

    auto g_0_xzz_0_yy = buffer_sfsd[33];

    auto g_0_xzz_0_yz = buffer_sfsd[34];

    auto g_0_xzz_0_zz = buffer_sfsd[35];

    auto g_0_yyy_0_xx = buffer_sfsd[36];

    auto g_0_yyy_0_xy = buffer_sfsd[37];

    auto g_0_yyy_0_xz = buffer_sfsd[38];

    auto g_0_yyy_0_yy = buffer_sfsd[39];

    auto g_0_yyy_0_yz = buffer_sfsd[40];

    auto g_0_yyy_0_zz = buffer_sfsd[41];

    auto g_0_yyz_0_xx = buffer_sfsd[42];

    auto g_0_yyz_0_xy = buffer_sfsd[43];

    auto g_0_yyz_0_xz = buffer_sfsd[44];

    auto g_0_yyz_0_yy = buffer_sfsd[45];

    auto g_0_yyz_0_yz = buffer_sfsd[46];

    auto g_0_yyz_0_zz = buffer_sfsd[47];

    auto g_0_yzz_0_xx = buffer_sfsd[48];

    auto g_0_yzz_0_xy = buffer_sfsd[49];

    auto g_0_yzz_0_xz = buffer_sfsd[50];

    auto g_0_yzz_0_yy = buffer_sfsd[51];

    auto g_0_yzz_0_yz = buffer_sfsd[52];

    auto g_0_yzz_0_zz = buffer_sfsd[53];

    auto g_0_zzz_0_xx = buffer_sfsd[54];

    auto g_0_zzz_0_xy = buffer_sfsd[55];

    auto g_0_zzz_0_xz = buffer_sfsd[56];

    auto g_0_zzz_0_yy = buffer_sfsd[57];

    auto g_0_zzz_0_yz = buffer_sfsd[58];

    auto g_0_zzz_0_zz = buffer_sfsd[59];

    /// Set up components of auxilary buffer : buffer_dpsd

    auto g_xx_x_0_xx = buffer_dpsd[0];

    auto g_xx_x_0_xy = buffer_dpsd[1];

    auto g_xx_x_0_xz = buffer_dpsd[2];

    auto g_xx_x_0_yy = buffer_dpsd[3];

    auto g_xx_x_0_yz = buffer_dpsd[4];

    auto g_xx_x_0_zz = buffer_dpsd[5];

    auto g_xx_y_0_xx = buffer_dpsd[6];

    auto g_xx_y_0_xy = buffer_dpsd[7];

    auto g_xx_y_0_xz = buffer_dpsd[8];

    auto g_xx_y_0_yy = buffer_dpsd[9];

    auto g_xx_y_0_yz = buffer_dpsd[10];

    auto g_xx_y_0_zz = buffer_dpsd[11];

    auto g_xx_z_0_xx = buffer_dpsd[12];

    auto g_xx_z_0_xy = buffer_dpsd[13];

    auto g_xx_z_0_xz = buffer_dpsd[14];

    auto g_xx_z_0_yy = buffer_dpsd[15];

    auto g_xx_z_0_yz = buffer_dpsd[16];

    auto g_xx_z_0_zz = buffer_dpsd[17];

    auto g_xy_x_0_xx = buffer_dpsd[18];

    auto g_xy_x_0_xy = buffer_dpsd[19];

    auto g_xy_x_0_xz = buffer_dpsd[20];

    auto g_xy_x_0_yy = buffer_dpsd[21];

    auto g_xy_x_0_yz = buffer_dpsd[22];

    auto g_xy_x_0_zz = buffer_dpsd[23];

    auto g_xy_y_0_xx = buffer_dpsd[24];

    auto g_xy_y_0_xy = buffer_dpsd[25];

    auto g_xy_y_0_xz = buffer_dpsd[26];

    auto g_xy_y_0_yy = buffer_dpsd[27];

    auto g_xy_y_0_yz = buffer_dpsd[28];

    auto g_xy_y_0_zz = buffer_dpsd[29];

    auto g_xy_z_0_xx = buffer_dpsd[30];

    auto g_xy_z_0_xy = buffer_dpsd[31];

    auto g_xy_z_0_xz = buffer_dpsd[32];

    auto g_xy_z_0_yy = buffer_dpsd[33];

    auto g_xy_z_0_yz = buffer_dpsd[34];

    auto g_xy_z_0_zz = buffer_dpsd[35];

    auto g_xz_x_0_xx = buffer_dpsd[36];

    auto g_xz_x_0_xy = buffer_dpsd[37];

    auto g_xz_x_0_xz = buffer_dpsd[38];

    auto g_xz_x_0_yy = buffer_dpsd[39];

    auto g_xz_x_0_yz = buffer_dpsd[40];

    auto g_xz_x_0_zz = buffer_dpsd[41];

    auto g_xz_y_0_xx = buffer_dpsd[42];

    auto g_xz_y_0_xy = buffer_dpsd[43];

    auto g_xz_y_0_xz = buffer_dpsd[44];

    auto g_xz_y_0_yy = buffer_dpsd[45];

    auto g_xz_y_0_yz = buffer_dpsd[46];

    auto g_xz_y_0_zz = buffer_dpsd[47];

    auto g_xz_z_0_xx = buffer_dpsd[48];

    auto g_xz_z_0_xy = buffer_dpsd[49];

    auto g_xz_z_0_xz = buffer_dpsd[50];

    auto g_xz_z_0_yy = buffer_dpsd[51];

    auto g_xz_z_0_yz = buffer_dpsd[52];

    auto g_xz_z_0_zz = buffer_dpsd[53];

    auto g_yy_x_0_xx = buffer_dpsd[54];

    auto g_yy_x_0_xy = buffer_dpsd[55];

    auto g_yy_x_0_xz = buffer_dpsd[56];

    auto g_yy_x_0_yy = buffer_dpsd[57];

    auto g_yy_x_0_yz = buffer_dpsd[58];

    auto g_yy_x_0_zz = buffer_dpsd[59];

    auto g_yy_y_0_xx = buffer_dpsd[60];

    auto g_yy_y_0_xy = buffer_dpsd[61];

    auto g_yy_y_0_xz = buffer_dpsd[62];

    auto g_yy_y_0_yy = buffer_dpsd[63];

    auto g_yy_y_0_yz = buffer_dpsd[64];

    auto g_yy_y_0_zz = buffer_dpsd[65];

    auto g_yy_z_0_xx = buffer_dpsd[66];

    auto g_yy_z_0_xy = buffer_dpsd[67];

    auto g_yy_z_0_xz = buffer_dpsd[68];

    auto g_yy_z_0_yy = buffer_dpsd[69];

    auto g_yy_z_0_yz = buffer_dpsd[70];

    auto g_yy_z_0_zz = buffer_dpsd[71];

    auto g_yz_x_0_xx = buffer_dpsd[72];

    auto g_yz_x_0_xy = buffer_dpsd[73];

    auto g_yz_x_0_xz = buffer_dpsd[74];

    auto g_yz_x_0_yy = buffer_dpsd[75];

    auto g_yz_x_0_yz = buffer_dpsd[76];

    auto g_yz_x_0_zz = buffer_dpsd[77];

    auto g_yz_y_0_xx = buffer_dpsd[78];

    auto g_yz_y_0_xy = buffer_dpsd[79];

    auto g_yz_y_0_xz = buffer_dpsd[80];

    auto g_yz_y_0_yy = buffer_dpsd[81];

    auto g_yz_y_0_yz = buffer_dpsd[82];

    auto g_yz_y_0_zz = buffer_dpsd[83];

    auto g_yz_z_0_xx = buffer_dpsd[84];

    auto g_yz_z_0_xy = buffer_dpsd[85];

    auto g_yz_z_0_xz = buffer_dpsd[86];

    auto g_yz_z_0_yy = buffer_dpsd[87];

    auto g_yz_z_0_yz = buffer_dpsd[88];

    auto g_yz_z_0_zz = buffer_dpsd[89];

    auto g_zz_x_0_xx = buffer_dpsd[90];

    auto g_zz_x_0_xy = buffer_dpsd[91];

    auto g_zz_x_0_xz = buffer_dpsd[92];

    auto g_zz_x_0_yy = buffer_dpsd[93];

    auto g_zz_x_0_yz = buffer_dpsd[94];

    auto g_zz_x_0_zz = buffer_dpsd[95];

    auto g_zz_y_0_xx = buffer_dpsd[96];

    auto g_zz_y_0_xy = buffer_dpsd[97];

    auto g_zz_y_0_xz = buffer_dpsd[98];

    auto g_zz_y_0_yy = buffer_dpsd[99];

    auto g_zz_y_0_yz = buffer_dpsd[100];

    auto g_zz_y_0_zz = buffer_dpsd[101];

    auto g_zz_z_0_xx = buffer_dpsd[102];

    auto g_zz_z_0_xy = buffer_dpsd[103];

    auto g_zz_z_0_xz = buffer_dpsd[104];

    auto g_zz_z_0_yy = buffer_dpsd[105];

    auto g_zz_z_0_yz = buffer_dpsd[106];

    auto g_zz_z_0_zz = buffer_dpsd[107];

    /// Set up components of auxilary buffer : buffer_dfsd

    auto g_xx_xxx_0_xx = buffer_dfsd[0];

    auto g_xx_xxx_0_xy = buffer_dfsd[1];

    auto g_xx_xxx_0_xz = buffer_dfsd[2];

    auto g_xx_xxx_0_yy = buffer_dfsd[3];

    auto g_xx_xxx_0_yz = buffer_dfsd[4];

    auto g_xx_xxx_0_zz = buffer_dfsd[5];

    auto g_xx_xxy_0_xx = buffer_dfsd[6];

    auto g_xx_xxy_0_xy = buffer_dfsd[7];

    auto g_xx_xxy_0_xz = buffer_dfsd[8];

    auto g_xx_xxy_0_yy = buffer_dfsd[9];

    auto g_xx_xxy_0_yz = buffer_dfsd[10];

    auto g_xx_xxy_0_zz = buffer_dfsd[11];

    auto g_xx_xxz_0_xx = buffer_dfsd[12];

    auto g_xx_xxz_0_xy = buffer_dfsd[13];

    auto g_xx_xxz_0_xz = buffer_dfsd[14];

    auto g_xx_xxz_0_yy = buffer_dfsd[15];

    auto g_xx_xxz_0_yz = buffer_dfsd[16];

    auto g_xx_xxz_0_zz = buffer_dfsd[17];

    auto g_xx_xyy_0_xx = buffer_dfsd[18];

    auto g_xx_xyy_0_xy = buffer_dfsd[19];

    auto g_xx_xyy_0_xz = buffer_dfsd[20];

    auto g_xx_xyy_0_yy = buffer_dfsd[21];

    auto g_xx_xyy_0_yz = buffer_dfsd[22];

    auto g_xx_xyy_0_zz = buffer_dfsd[23];

    auto g_xx_xyz_0_xx = buffer_dfsd[24];

    auto g_xx_xyz_0_xy = buffer_dfsd[25];

    auto g_xx_xyz_0_xz = buffer_dfsd[26];

    auto g_xx_xyz_0_yy = buffer_dfsd[27];

    auto g_xx_xyz_0_yz = buffer_dfsd[28];

    auto g_xx_xyz_0_zz = buffer_dfsd[29];

    auto g_xx_xzz_0_xx = buffer_dfsd[30];

    auto g_xx_xzz_0_xy = buffer_dfsd[31];

    auto g_xx_xzz_0_xz = buffer_dfsd[32];

    auto g_xx_xzz_0_yy = buffer_dfsd[33];

    auto g_xx_xzz_0_yz = buffer_dfsd[34];

    auto g_xx_xzz_0_zz = buffer_dfsd[35];

    auto g_xx_yyy_0_xx = buffer_dfsd[36];

    auto g_xx_yyy_0_xy = buffer_dfsd[37];

    auto g_xx_yyy_0_xz = buffer_dfsd[38];

    auto g_xx_yyy_0_yy = buffer_dfsd[39];

    auto g_xx_yyy_0_yz = buffer_dfsd[40];

    auto g_xx_yyy_0_zz = buffer_dfsd[41];

    auto g_xx_yyz_0_xx = buffer_dfsd[42];

    auto g_xx_yyz_0_xy = buffer_dfsd[43];

    auto g_xx_yyz_0_xz = buffer_dfsd[44];

    auto g_xx_yyz_0_yy = buffer_dfsd[45];

    auto g_xx_yyz_0_yz = buffer_dfsd[46];

    auto g_xx_yyz_0_zz = buffer_dfsd[47];

    auto g_xx_yzz_0_xx = buffer_dfsd[48];

    auto g_xx_yzz_0_xy = buffer_dfsd[49];

    auto g_xx_yzz_0_xz = buffer_dfsd[50];

    auto g_xx_yzz_0_yy = buffer_dfsd[51];

    auto g_xx_yzz_0_yz = buffer_dfsd[52];

    auto g_xx_yzz_0_zz = buffer_dfsd[53];

    auto g_xx_zzz_0_xx = buffer_dfsd[54];

    auto g_xx_zzz_0_xy = buffer_dfsd[55];

    auto g_xx_zzz_0_xz = buffer_dfsd[56];

    auto g_xx_zzz_0_yy = buffer_dfsd[57];

    auto g_xx_zzz_0_yz = buffer_dfsd[58];

    auto g_xx_zzz_0_zz = buffer_dfsd[59];

    auto g_xy_xxx_0_xx = buffer_dfsd[60];

    auto g_xy_xxx_0_xy = buffer_dfsd[61];

    auto g_xy_xxx_0_xz = buffer_dfsd[62];

    auto g_xy_xxx_0_yy = buffer_dfsd[63];

    auto g_xy_xxx_0_yz = buffer_dfsd[64];

    auto g_xy_xxx_0_zz = buffer_dfsd[65];

    auto g_xy_xxy_0_xx = buffer_dfsd[66];

    auto g_xy_xxy_0_xy = buffer_dfsd[67];

    auto g_xy_xxy_0_xz = buffer_dfsd[68];

    auto g_xy_xxy_0_yy = buffer_dfsd[69];

    auto g_xy_xxy_0_yz = buffer_dfsd[70];

    auto g_xy_xxy_0_zz = buffer_dfsd[71];

    auto g_xy_xxz_0_xx = buffer_dfsd[72];

    auto g_xy_xxz_0_xy = buffer_dfsd[73];

    auto g_xy_xxz_0_xz = buffer_dfsd[74];

    auto g_xy_xxz_0_yy = buffer_dfsd[75];

    auto g_xy_xxz_0_yz = buffer_dfsd[76];

    auto g_xy_xxz_0_zz = buffer_dfsd[77];

    auto g_xy_xyy_0_xx = buffer_dfsd[78];

    auto g_xy_xyy_0_xy = buffer_dfsd[79];

    auto g_xy_xyy_0_xz = buffer_dfsd[80];

    auto g_xy_xyy_0_yy = buffer_dfsd[81];

    auto g_xy_xyy_0_yz = buffer_dfsd[82];

    auto g_xy_xyy_0_zz = buffer_dfsd[83];

    auto g_xy_xyz_0_xx = buffer_dfsd[84];

    auto g_xy_xyz_0_xy = buffer_dfsd[85];

    auto g_xy_xyz_0_xz = buffer_dfsd[86];

    auto g_xy_xyz_0_yy = buffer_dfsd[87];

    auto g_xy_xyz_0_yz = buffer_dfsd[88];

    auto g_xy_xyz_0_zz = buffer_dfsd[89];

    auto g_xy_xzz_0_xx = buffer_dfsd[90];

    auto g_xy_xzz_0_xy = buffer_dfsd[91];

    auto g_xy_xzz_0_xz = buffer_dfsd[92];

    auto g_xy_xzz_0_yy = buffer_dfsd[93];

    auto g_xy_xzz_0_yz = buffer_dfsd[94];

    auto g_xy_xzz_0_zz = buffer_dfsd[95];

    auto g_xy_yyy_0_xx = buffer_dfsd[96];

    auto g_xy_yyy_0_xy = buffer_dfsd[97];

    auto g_xy_yyy_0_xz = buffer_dfsd[98];

    auto g_xy_yyy_0_yy = buffer_dfsd[99];

    auto g_xy_yyy_0_yz = buffer_dfsd[100];

    auto g_xy_yyy_0_zz = buffer_dfsd[101];

    auto g_xy_yyz_0_xx = buffer_dfsd[102];

    auto g_xy_yyz_0_xy = buffer_dfsd[103];

    auto g_xy_yyz_0_xz = buffer_dfsd[104];

    auto g_xy_yyz_0_yy = buffer_dfsd[105];

    auto g_xy_yyz_0_yz = buffer_dfsd[106];

    auto g_xy_yyz_0_zz = buffer_dfsd[107];

    auto g_xy_yzz_0_xx = buffer_dfsd[108];

    auto g_xy_yzz_0_xy = buffer_dfsd[109];

    auto g_xy_yzz_0_xz = buffer_dfsd[110];

    auto g_xy_yzz_0_yy = buffer_dfsd[111];

    auto g_xy_yzz_0_yz = buffer_dfsd[112];

    auto g_xy_yzz_0_zz = buffer_dfsd[113];

    auto g_xy_zzz_0_xx = buffer_dfsd[114];

    auto g_xy_zzz_0_xy = buffer_dfsd[115];

    auto g_xy_zzz_0_xz = buffer_dfsd[116];

    auto g_xy_zzz_0_yy = buffer_dfsd[117];

    auto g_xy_zzz_0_yz = buffer_dfsd[118];

    auto g_xy_zzz_0_zz = buffer_dfsd[119];

    auto g_xz_xxx_0_xx = buffer_dfsd[120];

    auto g_xz_xxx_0_xy = buffer_dfsd[121];

    auto g_xz_xxx_0_xz = buffer_dfsd[122];

    auto g_xz_xxx_0_yy = buffer_dfsd[123];

    auto g_xz_xxx_0_yz = buffer_dfsd[124];

    auto g_xz_xxx_0_zz = buffer_dfsd[125];

    auto g_xz_xxy_0_xx = buffer_dfsd[126];

    auto g_xz_xxy_0_xy = buffer_dfsd[127];

    auto g_xz_xxy_0_xz = buffer_dfsd[128];

    auto g_xz_xxy_0_yy = buffer_dfsd[129];

    auto g_xz_xxy_0_yz = buffer_dfsd[130];

    auto g_xz_xxy_0_zz = buffer_dfsd[131];

    auto g_xz_xxz_0_xx = buffer_dfsd[132];

    auto g_xz_xxz_0_xy = buffer_dfsd[133];

    auto g_xz_xxz_0_xz = buffer_dfsd[134];

    auto g_xz_xxz_0_yy = buffer_dfsd[135];

    auto g_xz_xxz_0_yz = buffer_dfsd[136];

    auto g_xz_xxz_0_zz = buffer_dfsd[137];

    auto g_xz_xyy_0_xx = buffer_dfsd[138];

    auto g_xz_xyy_0_xy = buffer_dfsd[139];

    auto g_xz_xyy_0_xz = buffer_dfsd[140];

    auto g_xz_xyy_0_yy = buffer_dfsd[141];

    auto g_xz_xyy_0_yz = buffer_dfsd[142];

    auto g_xz_xyy_0_zz = buffer_dfsd[143];

    auto g_xz_xyz_0_xx = buffer_dfsd[144];

    auto g_xz_xyz_0_xy = buffer_dfsd[145];

    auto g_xz_xyz_0_xz = buffer_dfsd[146];

    auto g_xz_xyz_0_yy = buffer_dfsd[147];

    auto g_xz_xyz_0_yz = buffer_dfsd[148];

    auto g_xz_xyz_0_zz = buffer_dfsd[149];

    auto g_xz_xzz_0_xx = buffer_dfsd[150];

    auto g_xz_xzz_0_xy = buffer_dfsd[151];

    auto g_xz_xzz_0_xz = buffer_dfsd[152];

    auto g_xz_xzz_0_yy = buffer_dfsd[153];

    auto g_xz_xzz_0_yz = buffer_dfsd[154];

    auto g_xz_xzz_0_zz = buffer_dfsd[155];

    auto g_xz_yyy_0_xx = buffer_dfsd[156];

    auto g_xz_yyy_0_xy = buffer_dfsd[157];

    auto g_xz_yyy_0_xz = buffer_dfsd[158];

    auto g_xz_yyy_0_yy = buffer_dfsd[159];

    auto g_xz_yyy_0_yz = buffer_dfsd[160];

    auto g_xz_yyy_0_zz = buffer_dfsd[161];

    auto g_xz_yyz_0_xx = buffer_dfsd[162];

    auto g_xz_yyz_0_xy = buffer_dfsd[163];

    auto g_xz_yyz_0_xz = buffer_dfsd[164];

    auto g_xz_yyz_0_yy = buffer_dfsd[165];

    auto g_xz_yyz_0_yz = buffer_dfsd[166];

    auto g_xz_yyz_0_zz = buffer_dfsd[167];

    auto g_xz_yzz_0_xx = buffer_dfsd[168];

    auto g_xz_yzz_0_xy = buffer_dfsd[169];

    auto g_xz_yzz_0_xz = buffer_dfsd[170];

    auto g_xz_yzz_0_yy = buffer_dfsd[171];

    auto g_xz_yzz_0_yz = buffer_dfsd[172];

    auto g_xz_yzz_0_zz = buffer_dfsd[173];

    auto g_xz_zzz_0_xx = buffer_dfsd[174];

    auto g_xz_zzz_0_xy = buffer_dfsd[175];

    auto g_xz_zzz_0_xz = buffer_dfsd[176];

    auto g_xz_zzz_0_yy = buffer_dfsd[177];

    auto g_xz_zzz_0_yz = buffer_dfsd[178];

    auto g_xz_zzz_0_zz = buffer_dfsd[179];

    auto g_yy_xxx_0_xx = buffer_dfsd[180];

    auto g_yy_xxx_0_xy = buffer_dfsd[181];

    auto g_yy_xxx_0_xz = buffer_dfsd[182];

    auto g_yy_xxx_0_yy = buffer_dfsd[183];

    auto g_yy_xxx_0_yz = buffer_dfsd[184];

    auto g_yy_xxx_0_zz = buffer_dfsd[185];

    auto g_yy_xxy_0_xx = buffer_dfsd[186];

    auto g_yy_xxy_0_xy = buffer_dfsd[187];

    auto g_yy_xxy_0_xz = buffer_dfsd[188];

    auto g_yy_xxy_0_yy = buffer_dfsd[189];

    auto g_yy_xxy_0_yz = buffer_dfsd[190];

    auto g_yy_xxy_0_zz = buffer_dfsd[191];

    auto g_yy_xxz_0_xx = buffer_dfsd[192];

    auto g_yy_xxz_0_xy = buffer_dfsd[193];

    auto g_yy_xxz_0_xz = buffer_dfsd[194];

    auto g_yy_xxz_0_yy = buffer_dfsd[195];

    auto g_yy_xxz_0_yz = buffer_dfsd[196];

    auto g_yy_xxz_0_zz = buffer_dfsd[197];

    auto g_yy_xyy_0_xx = buffer_dfsd[198];

    auto g_yy_xyy_0_xy = buffer_dfsd[199];

    auto g_yy_xyy_0_xz = buffer_dfsd[200];

    auto g_yy_xyy_0_yy = buffer_dfsd[201];

    auto g_yy_xyy_0_yz = buffer_dfsd[202];

    auto g_yy_xyy_0_zz = buffer_dfsd[203];

    auto g_yy_xyz_0_xx = buffer_dfsd[204];

    auto g_yy_xyz_0_xy = buffer_dfsd[205];

    auto g_yy_xyz_0_xz = buffer_dfsd[206];

    auto g_yy_xyz_0_yy = buffer_dfsd[207];

    auto g_yy_xyz_0_yz = buffer_dfsd[208];

    auto g_yy_xyz_0_zz = buffer_dfsd[209];

    auto g_yy_xzz_0_xx = buffer_dfsd[210];

    auto g_yy_xzz_0_xy = buffer_dfsd[211];

    auto g_yy_xzz_0_xz = buffer_dfsd[212];

    auto g_yy_xzz_0_yy = buffer_dfsd[213];

    auto g_yy_xzz_0_yz = buffer_dfsd[214];

    auto g_yy_xzz_0_zz = buffer_dfsd[215];

    auto g_yy_yyy_0_xx = buffer_dfsd[216];

    auto g_yy_yyy_0_xy = buffer_dfsd[217];

    auto g_yy_yyy_0_xz = buffer_dfsd[218];

    auto g_yy_yyy_0_yy = buffer_dfsd[219];

    auto g_yy_yyy_0_yz = buffer_dfsd[220];

    auto g_yy_yyy_0_zz = buffer_dfsd[221];

    auto g_yy_yyz_0_xx = buffer_dfsd[222];

    auto g_yy_yyz_0_xy = buffer_dfsd[223];

    auto g_yy_yyz_0_xz = buffer_dfsd[224];

    auto g_yy_yyz_0_yy = buffer_dfsd[225];

    auto g_yy_yyz_0_yz = buffer_dfsd[226];

    auto g_yy_yyz_0_zz = buffer_dfsd[227];

    auto g_yy_yzz_0_xx = buffer_dfsd[228];

    auto g_yy_yzz_0_xy = buffer_dfsd[229];

    auto g_yy_yzz_0_xz = buffer_dfsd[230];

    auto g_yy_yzz_0_yy = buffer_dfsd[231];

    auto g_yy_yzz_0_yz = buffer_dfsd[232];

    auto g_yy_yzz_0_zz = buffer_dfsd[233];

    auto g_yy_zzz_0_xx = buffer_dfsd[234];

    auto g_yy_zzz_0_xy = buffer_dfsd[235];

    auto g_yy_zzz_0_xz = buffer_dfsd[236];

    auto g_yy_zzz_0_yy = buffer_dfsd[237];

    auto g_yy_zzz_0_yz = buffer_dfsd[238];

    auto g_yy_zzz_0_zz = buffer_dfsd[239];

    auto g_yz_xxx_0_xx = buffer_dfsd[240];

    auto g_yz_xxx_0_xy = buffer_dfsd[241];

    auto g_yz_xxx_0_xz = buffer_dfsd[242];

    auto g_yz_xxx_0_yy = buffer_dfsd[243];

    auto g_yz_xxx_0_yz = buffer_dfsd[244];

    auto g_yz_xxx_0_zz = buffer_dfsd[245];

    auto g_yz_xxy_0_xx = buffer_dfsd[246];

    auto g_yz_xxy_0_xy = buffer_dfsd[247];

    auto g_yz_xxy_0_xz = buffer_dfsd[248];

    auto g_yz_xxy_0_yy = buffer_dfsd[249];

    auto g_yz_xxy_0_yz = buffer_dfsd[250];

    auto g_yz_xxy_0_zz = buffer_dfsd[251];

    auto g_yz_xxz_0_xx = buffer_dfsd[252];

    auto g_yz_xxz_0_xy = buffer_dfsd[253];

    auto g_yz_xxz_0_xz = buffer_dfsd[254];

    auto g_yz_xxz_0_yy = buffer_dfsd[255];

    auto g_yz_xxz_0_yz = buffer_dfsd[256];

    auto g_yz_xxz_0_zz = buffer_dfsd[257];

    auto g_yz_xyy_0_xx = buffer_dfsd[258];

    auto g_yz_xyy_0_xy = buffer_dfsd[259];

    auto g_yz_xyy_0_xz = buffer_dfsd[260];

    auto g_yz_xyy_0_yy = buffer_dfsd[261];

    auto g_yz_xyy_0_yz = buffer_dfsd[262];

    auto g_yz_xyy_0_zz = buffer_dfsd[263];

    auto g_yz_xyz_0_xx = buffer_dfsd[264];

    auto g_yz_xyz_0_xy = buffer_dfsd[265];

    auto g_yz_xyz_0_xz = buffer_dfsd[266];

    auto g_yz_xyz_0_yy = buffer_dfsd[267];

    auto g_yz_xyz_0_yz = buffer_dfsd[268];

    auto g_yz_xyz_0_zz = buffer_dfsd[269];

    auto g_yz_xzz_0_xx = buffer_dfsd[270];

    auto g_yz_xzz_0_xy = buffer_dfsd[271];

    auto g_yz_xzz_0_xz = buffer_dfsd[272];

    auto g_yz_xzz_0_yy = buffer_dfsd[273];

    auto g_yz_xzz_0_yz = buffer_dfsd[274];

    auto g_yz_xzz_0_zz = buffer_dfsd[275];

    auto g_yz_yyy_0_xx = buffer_dfsd[276];

    auto g_yz_yyy_0_xy = buffer_dfsd[277];

    auto g_yz_yyy_0_xz = buffer_dfsd[278];

    auto g_yz_yyy_0_yy = buffer_dfsd[279];

    auto g_yz_yyy_0_yz = buffer_dfsd[280];

    auto g_yz_yyy_0_zz = buffer_dfsd[281];

    auto g_yz_yyz_0_xx = buffer_dfsd[282];

    auto g_yz_yyz_0_xy = buffer_dfsd[283];

    auto g_yz_yyz_0_xz = buffer_dfsd[284];

    auto g_yz_yyz_0_yy = buffer_dfsd[285];

    auto g_yz_yyz_0_yz = buffer_dfsd[286];

    auto g_yz_yyz_0_zz = buffer_dfsd[287];

    auto g_yz_yzz_0_xx = buffer_dfsd[288];

    auto g_yz_yzz_0_xy = buffer_dfsd[289];

    auto g_yz_yzz_0_xz = buffer_dfsd[290];

    auto g_yz_yzz_0_yy = buffer_dfsd[291];

    auto g_yz_yzz_0_yz = buffer_dfsd[292];

    auto g_yz_yzz_0_zz = buffer_dfsd[293];

    auto g_yz_zzz_0_xx = buffer_dfsd[294];

    auto g_yz_zzz_0_xy = buffer_dfsd[295];

    auto g_yz_zzz_0_xz = buffer_dfsd[296];

    auto g_yz_zzz_0_yy = buffer_dfsd[297];

    auto g_yz_zzz_0_yz = buffer_dfsd[298];

    auto g_yz_zzz_0_zz = buffer_dfsd[299];

    auto g_zz_xxx_0_xx = buffer_dfsd[300];

    auto g_zz_xxx_0_xy = buffer_dfsd[301];

    auto g_zz_xxx_0_xz = buffer_dfsd[302];

    auto g_zz_xxx_0_yy = buffer_dfsd[303];

    auto g_zz_xxx_0_yz = buffer_dfsd[304];

    auto g_zz_xxx_0_zz = buffer_dfsd[305];

    auto g_zz_xxy_0_xx = buffer_dfsd[306];

    auto g_zz_xxy_0_xy = buffer_dfsd[307];

    auto g_zz_xxy_0_xz = buffer_dfsd[308];

    auto g_zz_xxy_0_yy = buffer_dfsd[309];

    auto g_zz_xxy_0_yz = buffer_dfsd[310];

    auto g_zz_xxy_0_zz = buffer_dfsd[311];

    auto g_zz_xxz_0_xx = buffer_dfsd[312];

    auto g_zz_xxz_0_xy = buffer_dfsd[313];

    auto g_zz_xxz_0_xz = buffer_dfsd[314];

    auto g_zz_xxz_0_yy = buffer_dfsd[315];

    auto g_zz_xxz_0_yz = buffer_dfsd[316];

    auto g_zz_xxz_0_zz = buffer_dfsd[317];

    auto g_zz_xyy_0_xx = buffer_dfsd[318];

    auto g_zz_xyy_0_xy = buffer_dfsd[319];

    auto g_zz_xyy_0_xz = buffer_dfsd[320];

    auto g_zz_xyy_0_yy = buffer_dfsd[321];

    auto g_zz_xyy_0_yz = buffer_dfsd[322];

    auto g_zz_xyy_0_zz = buffer_dfsd[323];

    auto g_zz_xyz_0_xx = buffer_dfsd[324];

    auto g_zz_xyz_0_xy = buffer_dfsd[325];

    auto g_zz_xyz_0_xz = buffer_dfsd[326];

    auto g_zz_xyz_0_yy = buffer_dfsd[327];

    auto g_zz_xyz_0_yz = buffer_dfsd[328];

    auto g_zz_xyz_0_zz = buffer_dfsd[329];

    auto g_zz_xzz_0_xx = buffer_dfsd[330];

    auto g_zz_xzz_0_xy = buffer_dfsd[331];

    auto g_zz_xzz_0_xz = buffer_dfsd[332];

    auto g_zz_xzz_0_yy = buffer_dfsd[333];

    auto g_zz_xzz_0_yz = buffer_dfsd[334];

    auto g_zz_xzz_0_zz = buffer_dfsd[335];

    auto g_zz_yyy_0_xx = buffer_dfsd[336];

    auto g_zz_yyy_0_xy = buffer_dfsd[337];

    auto g_zz_yyy_0_xz = buffer_dfsd[338];

    auto g_zz_yyy_0_yy = buffer_dfsd[339];

    auto g_zz_yyy_0_yz = buffer_dfsd[340];

    auto g_zz_yyy_0_zz = buffer_dfsd[341];

    auto g_zz_yyz_0_xx = buffer_dfsd[342];

    auto g_zz_yyz_0_xy = buffer_dfsd[343];

    auto g_zz_yyz_0_xz = buffer_dfsd[344];

    auto g_zz_yyz_0_yy = buffer_dfsd[345];

    auto g_zz_yyz_0_yz = buffer_dfsd[346];

    auto g_zz_yyz_0_zz = buffer_dfsd[347];

    auto g_zz_yzz_0_xx = buffer_dfsd[348];

    auto g_zz_yzz_0_xy = buffer_dfsd[349];

    auto g_zz_yzz_0_xz = buffer_dfsd[350];

    auto g_zz_yzz_0_yy = buffer_dfsd[351];

    auto g_zz_yzz_0_yz = buffer_dfsd[352];

    auto g_zz_yzz_0_zz = buffer_dfsd[353];

    auto g_zz_zzz_0_xx = buffer_dfsd[354];

    auto g_zz_zzz_0_xy = buffer_dfsd[355];

    auto g_zz_zzz_0_xz = buffer_dfsd[356];

    auto g_zz_zzz_0_yy = buffer_dfsd[357];

    auto g_zz_zzz_0_yz = buffer_dfsd[358];

    auto g_zz_zzz_0_zz = buffer_dfsd[359];

    /// Set up components of integrals buffer : buffer_1100_pdsd

    auto g_x_x_0_0_x_xx_0_xx = buffer_1100_pdsd[0];

    auto g_x_x_0_0_x_xx_0_xy = buffer_1100_pdsd[1];

    auto g_x_x_0_0_x_xx_0_xz = buffer_1100_pdsd[2];

    auto g_x_x_0_0_x_xx_0_yy = buffer_1100_pdsd[3];

    auto g_x_x_0_0_x_xx_0_yz = buffer_1100_pdsd[4];

    auto g_x_x_0_0_x_xx_0_zz = buffer_1100_pdsd[5];

    auto g_x_x_0_0_x_xy_0_xx = buffer_1100_pdsd[6];

    auto g_x_x_0_0_x_xy_0_xy = buffer_1100_pdsd[7];

    auto g_x_x_0_0_x_xy_0_xz = buffer_1100_pdsd[8];

    auto g_x_x_0_0_x_xy_0_yy = buffer_1100_pdsd[9];

    auto g_x_x_0_0_x_xy_0_yz = buffer_1100_pdsd[10];

    auto g_x_x_0_0_x_xy_0_zz = buffer_1100_pdsd[11];

    auto g_x_x_0_0_x_xz_0_xx = buffer_1100_pdsd[12];

    auto g_x_x_0_0_x_xz_0_xy = buffer_1100_pdsd[13];

    auto g_x_x_0_0_x_xz_0_xz = buffer_1100_pdsd[14];

    auto g_x_x_0_0_x_xz_0_yy = buffer_1100_pdsd[15];

    auto g_x_x_0_0_x_xz_0_yz = buffer_1100_pdsd[16];

    auto g_x_x_0_0_x_xz_0_zz = buffer_1100_pdsd[17];

    auto g_x_x_0_0_x_yy_0_xx = buffer_1100_pdsd[18];

    auto g_x_x_0_0_x_yy_0_xy = buffer_1100_pdsd[19];

    auto g_x_x_0_0_x_yy_0_xz = buffer_1100_pdsd[20];

    auto g_x_x_0_0_x_yy_0_yy = buffer_1100_pdsd[21];

    auto g_x_x_0_0_x_yy_0_yz = buffer_1100_pdsd[22];

    auto g_x_x_0_0_x_yy_0_zz = buffer_1100_pdsd[23];

    auto g_x_x_0_0_x_yz_0_xx = buffer_1100_pdsd[24];

    auto g_x_x_0_0_x_yz_0_xy = buffer_1100_pdsd[25];

    auto g_x_x_0_0_x_yz_0_xz = buffer_1100_pdsd[26];

    auto g_x_x_0_0_x_yz_0_yy = buffer_1100_pdsd[27];

    auto g_x_x_0_0_x_yz_0_yz = buffer_1100_pdsd[28];

    auto g_x_x_0_0_x_yz_0_zz = buffer_1100_pdsd[29];

    auto g_x_x_0_0_x_zz_0_xx = buffer_1100_pdsd[30];

    auto g_x_x_0_0_x_zz_0_xy = buffer_1100_pdsd[31];

    auto g_x_x_0_0_x_zz_0_xz = buffer_1100_pdsd[32];

    auto g_x_x_0_0_x_zz_0_yy = buffer_1100_pdsd[33];

    auto g_x_x_0_0_x_zz_0_yz = buffer_1100_pdsd[34];

    auto g_x_x_0_0_x_zz_0_zz = buffer_1100_pdsd[35];

    auto g_x_x_0_0_y_xx_0_xx = buffer_1100_pdsd[36];

    auto g_x_x_0_0_y_xx_0_xy = buffer_1100_pdsd[37];

    auto g_x_x_0_0_y_xx_0_xz = buffer_1100_pdsd[38];

    auto g_x_x_0_0_y_xx_0_yy = buffer_1100_pdsd[39];

    auto g_x_x_0_0_y_xx_0_yz = buffer_1100_pdsd[40];

    auto g_x_x_0_0_y_xx_0_zz = buffer_1100_pdsd[41];

    auto g_x_x_0_0_y_xy_0_xx = buffer_1100_pdsd[42];

    auto g_x_x_0_0_y_xy_0_xy = buffer_1100_pdsd[43];

    auto g_x_x_0_0_y_xy_0_xz = buffer_1100_pdsd[44];

    auto g_x_x_0_0_y_xy_0_yy = buffer_1100_pdsd[45];

    auto g_x_x_0_0_y_xy_0_yz = buffer_1100_pdsd[46];

    auto g_x_x_0_0_y_xy_0_zz = buffer_1100_pdsd[47];

    auto g_x_x_0_0_y_xz_0_xx = buffer_1100_pdsd[48];

    auto g_x_x_0_0_y_xz_0_xy = buffer_1100_pdsd[49];

    auto g_x_x_0_0_y_xz_0_xz = buffer_1100_pdsd[50];

    auto g_x_x_0_0_y_xz_0_yy = buffer_1100_pdsd[51];

    auto g_x_x_0_0_y_xz_0_yz = buffer_1100_pdsd[52];

    auto g_x_x_0_0_y_xz_0_zz = buffer_1100_pdsd[53];

    auto g_x_x_0_0_y_yy_0_xx = buffer_1100_pdsd[54];

    auto g_x_x_0_0_y_yy_0_xy = buffer_1100_pdsd[55];

    auto g_x_x_0_0_y_yy_0_xz = buffer_1100_pdsd[56];

    auto g_x_x_0_0_y_yy_0_yy = buffer_1100_pdsd[57];

    auto g_x_x_0_0_y_yy_0_yz = buffer_1100_pdsd[58];

    auto g_x_x_0_0_y_yy_0_zz = buffer_1100_pdsd[59];

    auto g_x_x_0_0_y_yz_0_xx = buffer_1100_pdsd[60];

    auto g_x_x_0_0_y_yz_0_xy = buffer_1100_pdsd[61];

    auto g_x_x_0_0_y_yz_0_xz = buffer_1100_pdsd[62];

    auto g_x_x_0_0_y_yz_0_yy = buffer_1100_pdsd[63];

    auto g_x_x_0_0_y_yz_0_yz = buffer_1100_pdsd[64];

    auto g_x_x_0_0_y_yz_0_zz = buffer_1100_pdsd[65];

    auto g_x_x_0_0_y_zz_0_xx = buffer_1100_pdsd[66];

    auto g_x_x_0_0_y_zz_0_xy = buffer_1100_pdsd[67];

    auto g_x_x_0_0_y_zz_0_xz = buffer_1100_pdsd[68];

    auto g_x_x_0_0_y_zz_0_yy = buffer_1100_pdsd[69];

    auto g_x_x_0_0_y_zz_0_yz = buffer_1100_pdsd[70];

    auto g_x_x_0_0_y_zz_0_zz = buffer_1100_pdsd[71];

    auto g_x_x_0_0_z_xx_0_xx = buffer_1100_pdsd[72];

    auto g_x_x_0_0_z_xx_0_xy = buffer_1100_pdsd[73];

    auto g_x_x_0_0_z_xx_0_xz = buffer_1100_pdsd[74];

    auto g_x_x_0_0_z_xx_0_yy = buffer_1100_pdsd[75];

    auto g_x_x_0_0_z_xx_0_yz = buffer_1100_pdsd[76];

    auto g_x_x_0_0_z_xx_0_zz = buffer_1100_pdsd[77];

    auto g_x_x_0_0_z_xy_0_xx = buffer_1100_pdsd[78];

    auto g_x_x_0_0_z_xy_0_xy = buffer_1100_pdsd[79];

    auto g_x_x_0_0_z_xy_0_xz = buffer_1100_pdsd[80];

    auto g_x_x_0_0_z_xy_0_yy = buffer_1100_pdsd[81];

    auto g_x_x_0_0_z_xy_0_yz = buffer_1100_pdsd[82];

    auto g_x_x_0_0_z_xy_0_zz = buffer_1100_pdsd[83];

    auto g_x_x_0_0_z_xz_0_xx = buffer_1100_pdsd[84];

    auto g_x_x_0_0_z_xz_0_xy = buffer_1100_pdsd[85];

    auto g_x_x_0_0_z_xz_0_xz = buffer_1100_pdsd[86];

    auto g_x_x_0_0_z_xz_0_yy = buffer_1100_pdsd[87];

    auto g_x_x_0_0_z_xz_0_yz = buffer_1100_pdsd[88];

    auto g_x_x_0_0_z_xz_0_zz = buffer_1100_pdsd[89];

    auto g_x_x_0_0_z_yy_0_xx = buffer_1100_pdsd[90];

    auto g_x_x_0_0_z_yy_0_xy = buffer_1100_pdsd[91];

    auto g_x_x_0_0_z_yy_0_xz = buffer_1100_pdsd[92];

    auto g_x_x_0_0_z_yy_0_yy = buffer_1100_pdsd[93];

    auto g_x_x_0_0_z_yy_0_yz = buffer_1100_pdsd[94];

    auto g_x_x_0_0_z_yy_0_zz = buffer_1100_pdsd[95];

    auto g_x_x_0_0_z_yz_0_xx = buffer_1100_pdsd[96];

    auto g_x_x_0_0_z_yz_0_xy = buffer_1100_pdsd[97];

    auto g_x_x_0_0_z_yz_0_xz = buffer_1100_pdsd[98];

    auto g_x_x_0_0_z_yz_0_yy = buffer_1100_pdsd[99];

    auto g_x_x_0_0_z_yz_0_yz = buffer_1100_pdsd[100];

    auto g_x_x_0_0_z_yz_0_zz = buffer_1100_pdsd[101];

    auto g_x_x_0_0_z_zz_0_xx = buffer_1100_pdsd[102];

    auto g_x_x_0_0_z_zz_0_xy = buffer_1100_pdsd[103];

    auto g_x_x_0_0_z_zz_0_xz = buffer_1100_pdsd[104];

    auto g_x_x_0_0_z_zz_0_yy = buffer_1100_pdsd[105];

    auto g_x_x_0_0_z_zz_0_yz = buffer_1100_pdsd[106];

    auto g_x_x_0_0_z_zz_0_zz = buffer_1100_pdsd[107];

    auto g_x_y_0_0_x_xx_0_xx = buffer_1100_pdsd[108];

    auto g_x_y_0_0_x_xx_0_xy = buffer_1100_pdsd[109];

    auto g_x_y_0_0_x_xx_0_xz = buffer_1100_pdsd[110];

    auto g_x_y_0_0_x_xx_0_yy = buffer_1100_pdsd[111];

    auto g_x_y_0_0_x_xx_0_yz = buffer_1100_pdsd[112];

    auto g_x_y_0_0_x_xx_0_zz = buffer_1100_pdsd[113];

    auto g_x_y_0_0_x_xy_0_xx = buffer_1100_pdsd[114];

    auto g_x_y_0_0_x_xy_0_xy = buffer_1100_pdsd[115];

    auto g_x_y_0_0_x_xy_0_xz = buffer_1100_pdsd[116];

    auto g_x_y_0_0_x_xy_0_yy = buffer_1100_pdsd[117];

    auto g_x_y_0_0_x_xy_0_yz = buffer_1100_pdsd[118];

    auto g_x_y_0_0_x_xy_0_zz = buffer_1100_pdsd[119];

    auto g_x_y_0_0_x_xz_0_xx = buffer_1100_pdsd[120];

    auto g_x_y_0_0_x_xz_0_xy = buffer_1100_pdsd[121];

    auto g_x_y_0_0_x_xz_0_xz = buffer_1100_pdsd[122];

    auto g_x_y_0_0_x_xz_0_yy = buffer_1100_pdsd[123];

    auto g_x_y_0_0_x_xz_0_yz = buffer_1100_pdsd[124];

    auto g_x_y_0_0_x_xz_0_zz = buffer_1100_pdsd[125];

    auto g_x_y_0_0_x_yy_0_xx = buffer_1100_pdsd[126];

    auto g_x_y_0_0_x_yy_0_xy = buffer_1100_pdsd[127];

    auto g_x_y_0_0_x_yy_0_xz = buffer_1100_pdsd[128];

    auto g_x_y_0_0_x_yy_0_yy = buffer_1100_pdsd[129];

    auto g_x_y_0_0_x_yy_0_yz = buffer_1100_pdsd[130];

    auto g_x_y_0_0_x_yy_0_zz = buffer_1100_pdsd[131];

    auto g_x_y_0_0_x_yz_0_xx = buffer_1100_pdsd[132];

    auto g_x_y_0_0_x_yz_0_xy = buffer_1100_pdsd[133];

    auto g_x_y_0_0_x_yz_0_xz = buffer_1100_pdsd[134];

    auto g_x_y_0_0_x_yz_0_yy = buffer_1100_pdsd[135];

    auto g_x_y_0_0_x_yz_0_yz = buffer_1100_pdsd[136];

    auto g_x_y_0_0_x_yz_0_zz = buffer_1100_pdsd[137];

    auto g_x_y_0_0_x_zz_0_xx = buffer_1100_pdsd[138];

    auto g_x_y_0_0_x_zz_0_xy = buffer_1100_pdsd[139];

    auto g_x_y_0_0_x_zz_0_xz = buffer_1100_pdsd[140];

    auto g_x_y_0_0_x_zz_0_yy = buffer_1100_pdsd[141];

    auto g_x_y_0_0_x_zz_0_yz = buffer_1100_pdsd[142];

    auto g_x_y_0_0_x_zz_0_zz = buffer_1100_pdsd[143];

    auto g_x_y_0_0_y_xx_0_xx = buffer_1100_pdsd[144];

    auto g_x_y_0_0_y_xx_0_xy = buffer_1100_pdsd[145];

    auto g_x_y_0_0_y_xx_0_xz = buffer_1100_pdsd[146];

    auto g_x_y_0_0_y_xx_0_yy = buffer_1100_pdsd[147];

    auto g_x_y_0_0_y_xx_0_yz = buffer_1100_pdsd[148];

    auto g_x_y_0_0_y_xx_0_zz = buffer_1100_pdsd[149];

    auto g_x_y_0_0_y_xy_0_xx = buffer_1100_pdsd[150];

    auto g_x_y_0_0_y_xy_0_xy = buffer_1100_pdsd[151];

    auto g_x_y_0_0_y_xy_0_xz = buffer_1100_pdsd[152];

    auto g_x_y_0_0_y_xy_0_yy = buffer_1100_pdsd[153];

    auto g_x_y_0_0_y_xy_0_yz = buffer_1100_pdsd[154];

    auto g_x_y_0_0_y_xy_0_zz = buffer_1100_pdsd[155];

    auto g_x_y_0_0_y_xz_0_xx = buffer_1100_pdsd[156];

    auto g_x_y_0_0_y_xz_0_xy = buffer_1100_pdsd[157];

    auto g_x_y_0_0_y_xz_0_xz = buffer_1100_pdsd[158];

    auto g_x_y_0_0_y_xz_0_yy = buffer_1100_pdsd[159];

    auto g_x_y_0_0_y_xz_0_yz = buffer_1100_pdsd[160];

    auto g_x_y_0_0_y_xz_0_zz = buffer_1100_pdsd[161];

    auto g_x_y_0_0_y_yy_0_xx = buffer_1100_pdsd[162];

    auto g_x_y_0_0_y_yy_0_xy = buffer_1100_pdsd[163];

    auto g_x_y_0_0_y_yy_0_xz = buffer_1100_pdsd[164];

    auto g_x_y_0_0_y_yy_0_yy = buffer_1100_pdsd[165];

    auto g_x_y_0_0_y_yy_0_yz = buffer_1100_pdsd[166];

    auto g_x_y_0_0_y_yy_0_zz = buffer_1100_pdsd[167];

    auto g_x_y_0_0_y_yz_0_xx = buffer_1100_pdsd[168];

    auto g_x_y_0_0_y_yz_0_xy = buffer_1100_pdsd[169];

    auto g_x_y_0_0_y_yz_0_xz = buffer_1100_pdsd[170];

    auto g_x_y_0_0_y_yz_0_yy = buffer_1100_pdsd[171];

    auto g_x_y_0_0_y_yz_0_yz = buffer_1100_pdsd[172];

    auto g_x_y_0_0_y_yz_0_zz = buffer_1100_pdsd[173];

    auto g_x_y_0_0_y_zz_0_xx = buffer_1100_pdsd[174];

    auto g_x_y_0_0_y_zz_0_xy = buffer_1100_pdsd[175];

    auto g_x_y_0_0_y_zz_0_xz = buffer_1100_pdsd[176];

    auto g_x_y_0_0_y_zz_0_yy = buffer_1100_pdsd[177];

    auto g_x_y_0_0_y_zz_0_yz = buffer_1100_pdsd[178];

    auto g_x_y_0_0_y_zz_0_zz = buffer_1100_pdsd[179];

    auto g_x_y_0_0_z_xx_0_xx = buffer_1100_pdsd[180];

    auto g_x_y_0_0_z_xx_0_xy = buffer_1100_pdsd[181];

    auto g_x_y_0_0_z_xx_0_xz = buffer_1100_pdsd[182];

    auto g_x_y_0_0_z_xx_0_yy = buffer_1100_pdsd[183];

    auto g_x_y_0_0_z_xx_0_yz = buffer_1100_pdsd[184];

    auto g_x_y_0_0_z_xx_0_zz = buffer_1100_pdsd[185];

    auto g_x_y_0_0_z_xy_0_xx = buffer_1100_pdsd[186];

    auto g_x_y_0_0_z_xy_0_xy = buffer_1100_pdsd[187];

    auto g_x_y_0_0_z_xy_0_xz = buffer_1100_pdsd[188];

    auto g_x_y_0_0_z_xy_0_yy = buffer_1100_pdsd[189];

    auto g_x_y_0_0_z_xy_0_yz = buffer_1100_pdsd[190];

    auto g_x_y_0_0_z_xy_0_zz = buffer_1100_pdsd[191];

    auto g_x_y_0_0_z_xz_0_xx = buffer_1100_pdsd[192];

    auto g_x_y_0_0_z_xz_0_xy = buffer_1100_pdsd[193];

    auto g_x_y_0_0_z_xz_0_xz = buffer_1100_pdsd[194];

    auto g_x_y_0_0_z_xz_0_yy = buffer_1100_pdsd[195];

    auto g_x_y_0_0_z_xz_0_yz = buffer_1100_pdsd[196];

    auto g_x_y_0_0_z_xz_0_zz = buffer_1100_pdsd[197];

    auto g_x_y_0_0_z_yy_0_xx = buffer_1100_pdsd[198];

    auto g_x_y_0_0_z_yy_0_xy = buffer_1100_pdsd[199];

    auto g_x_y_0_0_z_yy_0_xz = buffer_1100_pdsd[200];

    auto g_x_y_0_0_z_yy_0_yy = buffer_1100_pdsd[201];

    auto g_x_y_0_0_z_yy_0_yz = buffer_1100_pdsd[202];

    auto g_x_y_0_0_z_yy_0_zz = buffer_1100_pdsd[203];

    auto g_x_y_0_0_z_yz_0_xx = buffer_1100_pdsd[204];

    auto g_x_y_0_0_z_yz_0_xy = buffer_1100_pdsd[205];

    auto g_x_y_0_0_z_yz_0_xz = buffer_1100_pdsd[206];

    auto g_x_y_0_0_z_yz_0_yy = buffer_1100_pdsd[207];

    auto g_x_y_0_0_z_yz_0_yz = buffer_1100_pdsd[208];

    auto g_x_y_0_0_z_yz_0_zz = buffer_1100_pdsd[209];

    auto g_x_y_0_0_z_zz_0_xx = buffer_1100_pdsd[210];

    auto g_x_y_0_0_z_zz_0_xy = buffer_1100_pdsd[211];

    auto g_x_y_0_0_z_zz_0_xz = buffer_1100_pdsd[212];

    auto g_x_y_0_0_z_zz_0_yy = buffer_1100_pdsd[213];

    auto g_x_y_0_0_z_zz_0_yz = buffer_1100_pdsd[214];

    auto g_x_y_0_0_z_zz_0_zz = buffer_1100_pdsd[215];

    auto g_x_z_0_0_x_xx_0_xx = buffer_1100_pdsd[216];

    auto g_x_z_0_0_x_xx_0_xy = buffer_1100_pdsd[217];

    auto g_x_z_0_0_x_xx_0_xz = buffer_1100_pdsd[218];

    auto g_x_z_0_0_x_xx_0_yy = buffer_1100_pdsd[219];

    auto g_x_z_0_0_x_xx_0_yz = buffer_1100_pdsd[220];

    auto g_x_z_0_0_x_xx_0_zz = buffer_1100_pdsd[221];

    auto g_x_z_0_0_x_xy_0_xx = buffer_1100_pdsd[222];

    auto g_x_z_0_0_x_xy_0_xy = buffer_1100_pdsd[223];

    auto g_x_z_0_0_x_xy_0_xz = buffer_1100_pdsd[224];

    auto g_x_z_0_0_x_xy_0_yy = buffer_1100_pdsd[225];

    auto g_x_z_0_0_x_xy_0_yz = buffer_1100_pdsd[226];

    auto g_x_z_0_0_x_xy_0_zz = buffer_1100_pdsd[227];

    auto g_x_z_0_0_x_xz_0_xx = buffer_1100_pdsd[228];

    auto g_x_z_0_0_x_xz_0_xy = buffer_1100_pdsd[229];

    auto g_x_z_0_0_x_xz_0_xz = buffer_1100_pdsd[230];

    auto g_x_z_0_0_x_xz_0_yy = buffer_1100_pdsd[231];

    auto g_x_z_0_0_x_xz_0_yz = buffer_1100_pdsd[232];

    auto g_x_z_0_0_x_xz_0_zz = buffer_1100_pdsd[233];

    auto g_x_z_0_0_x_yy_0_xx = buffer_1100_pdsd[234];

    auto g_x_z_0_0_x_yy_0_xy = buffer_1100_pdsd[235];

    auto g_x_z_0_0_x_yy_0_xz = buffer_1100_pdsd[236];

    auto g_x_z_0_0_x_yy_0_yy = buffer_1100_pdsd[237];

    auto g_x_z_0_0_x_yy_0_yz = buffer_1100_pdsd[238];

    auto g_x_z_0_0_x_yy_0_zz = buffer_1100_pdsd[239];

    auto g_x_z_0_0_x_yz_0_xx = buffer_1100_pdsd[240];

    auto g_x_z_0_0_x_yz_0_xy = buffer_1100_pdsd[241];

    auto g_x_z_0_0_x_yz_0_xz = buffer_1100_pdsd[242];

    auto g_x_z_0_0_x_yz_0_yy = buffer_1100_pdsd[243];

    auto g_x_z_0_0_x_yz_0_yz = buffer_1100_pdsd[244];

    auto g_x_z_0_0_x_yz_0_zz = buffer_1100_pdsd[245];

    auto g_x_z_0_0_x_zz_0_xx = buffer_1100_pdsd[246];

    auto g_x_z_0_0_x_zz_0_xy = buffer_1100_pdsd[247];

    auto g_x_z_0_0_x_zz_0_xz = buffer_1100_pdsd[248];

    auto g_x_z_0_0_x_zz_0_yy = buffer_1100_pdsd[249];

    auto g_x_z_0_0_x_zz_0_yz = buffer_1100_pdsd[250];

    auto g_x_z_0_0_x_zz_0_zz = buffer_1100_pdsd[251];

    auto g_x_z_0_0_y_xx_0_xx = buffer_1100_pdsd[252];

    auto g_x_z_0_0_y_xx_0_xy = buffer_1100_pdsd[253];

    auto g_x_z_0_0_y_xx_0_xz = buffer_1100_pdsd[254];

    auto g_x_z_0_0_y_xx_0_yy = buffer_1100_pdsd[255];

    auto g_x_z_0_0_y_xx_0_yz = buffer_1100_pdsd[256];

    auto g_x_z_0_0_y_xx_0_zz = buffer_1100_pdsd[257];

    auto g_x_z_0_0_y_xy_0_xx = buffer_1100_pdsd[258];

    auto g_x_z_0_0_y_xy_0_xy = buffer_1100_pdsd[259];

    auto g_x_z_0_0_y_xy_0_xz = buffer_1100_pdsd[260];

    auto g_x_z_0_0_y_xy_0_yy = buffer_1100_pdsd[261];

    auto g_x_z_0_0_y_xy_0_yz = buffer_1100_pdsd[262];

    auto g_x_z_0_0_y_xy_0_zz = buffer_1100_pdsd[263];

    auto g_x_z_0_0_y_xz_0_xx = buffer_1100_pdsd[264];

    auto g_x_z_0_0_y_xz_0_xy = buffer_1100_pdsd[265];

    auto g_x_z_0_0_y_xz_0_xz = buffer_1100_pdsd[266];

    auto g_x_z_0_0_y_xz_0_yy = buffer_1100_pdsd[267];

    auto g_x_z_0_0_y_xz_0_yz = buffer_1100_pdsd[268];

    auto g_x_z_0_0_y_xz_0_zz = buffer_1100_pdsd[269];

    auto g_x_z_0_0_y_yy_0_xx = buffer_1100_pdsd[270];

    auto g_x_z_0_0_y_yy_0_xy = buffer_1100_pdsd[271];

    auto g_x_z_0_0_y_yy_0_xz = buffer_1100_pdsd[272];

    auto g_x_z_0_0_y_yy_0_yy = buffer_1100_pdsd[273];

    auto g_x_z_0_0_y_yy_0_yz = buffer_1100_pdsd[274];

    auto g_x_z_0_0_y_yy_0_zz = buffer_1100_pdsd[275];

    auto g_x_z_0_0_y_yz_0_xx = buffer_1100_pdsd[276];

    auto g_x_z_0_0_y_yz_0_xy = buffer_1100_pdsd[277];

    auto g_x_z_0_0_y_yz_0_xz = buffer_1100_pdsd[278];

    auto g_x_z_0_0_y_yz_0_yy = buffer_1100_pdsd[279];

    auto g_x_z_0_0_y_yz_0_yz = buffer_1100_pdsd[280];

    auto g_x_z_0_0_y_yz_0_zz = buffer_1100_pdsd[281];

    auto g_x_z_0_0_y_zz_0_xx = buffer_1100_pdsd[282];

    auto g_x_z_0_0_y_zz_0_xy = buffer_1100_pdsd[283];

    auto g_x_z_0_0_y_zz_0_xz = buffer_1100_pdsd[284];

    auto g_x_z_0_0_y_zz_0_yy = buffer_1100_pdsd[285];

    auto g_x_z_0_0_y_zz_0_yz = buffer_1100_pdsd[286];

    auto g_x_z_0_0_y_zz_0_zz = buffer_1100_pdsd[287];

    auto g_x_z_0_0_z_xx_0_xx = buffer_1100_pdsd[288];

    auto g_x_z_0_0_z_xx_0_xy = buffer_1100_pdsd[289];

    auto g_x_z_0_0_z_xx_0_xz = buffer_1100_pdsd[290];

    auto g_x_z_0_0_z_xx_0_yy = buffer_1100_pdsd[291];

    auto g_x_z_0_0_z_xx_0_yz = buffer_1100_pdsd[292];

    auto g_x_z_0_0_z_xx_0_zz = buffer_1100_pdsd[293];

    auto g_x_z_0_0_z_xy_0_xx = buffer_1100_pdsd[294];

    auto g_x_z_0_0_z_xy_0_xy = buffer_1100_pdsd[295];

    auto g_x_z_0_0_z_xy_0_xz = buffer_1100_pdsd[296];

    auto g_x_z_0_0_z_xy_0_yy = buffer_1100_pdsd[297];

    auto g_x_z_0_0_z_xy_0_yz = buffer_1100_pdsd[298];

    auto g_x_z_0_0_z_xy_0_zz = buffer_1100_pdsd[299];

    auto g_x_z_0_0_z_xz_0_xx = buffer_1100_pdsd[300];

    auto g_x_z_0_0_z_xz_0_xy = buffer_1100_pdsd[301];

    auto g_x_z_0_0_z_xz_0_xz = buffer_1100_pdsd[302];

    auto g_x_z_0_0_z_xz_0_yy = buffer_1100_pdsd[303];

    auto g_x_z_0_0_z_xz_0_yz = buffer_1100_pdsd[304];

    auto g_x_z_0_0_z_xz_0_zz = buffer_1100_pdsd[305];

    auto g_x_z_0_0_z_yy_0_xx = buffer_1100_pdsd[306];

    auto g_x_z_0_0_z_yy_0_xy = buffer_1100_pdsd[307];

    auto g_x_z_0_0_z_yy_0_xz = buffer_1100_pdsd[308];

    auto g_x_z_0_0_z_yy_0_yy = buffer_1100_pdsd[309];

    auto g_x_z_0_0_z_yy_0_yz = buffer_1100_pdsd[310];

    auto g_x_z_0_0_z_yy_0_zz = buffer_1100_pdsd[311];

    auto g_x_z_0_0_z_yz_0_xx = buffer_1100_pdsd[312];

    auto g_x_z_0_0_z_yz_0_xy = buffer_1100_pdsd[313];

    auto g_x_z_0_0_z_yz_0_xz = buffer_1100_pdsd[314];

    auto g_x_z_0_0_z_yz_0_yy = buffer_1100_pdsd[315];

    auto g_x_z_0_0_z_yz_0_yz = buffer_1100_pdsd[316];

    auto g_x_z_0_0_z_yz_0_zz = buffer_1100_pdsd[317];

    auto g_x_z_0_0_z_zz_0_xx = buffer_1100_pdsd[318];

    auto g_x_z_0_0_z_zz_0_xy = buffer_1100_pdsd[319];

    auto g_x_z_0_0_z_zz_0_xz = buffer_1100_pdsd[320];

    auto g_x_z_0_0_z_zz_0_yy = buffer_1100_pdsd[321];

    auto g_x_z_0_0_z_zz_0_yz = buffer_1100_pdsd[322];

    auto g_x_z_0_0_z_zz_0_zz = buffer_1100_pdsd[323];

    auto g_y_x_0_0_x_xx_0_xx = buffer_1100_pdsd[324];

    auto g_y_x_0_0_x_xx_0_xy = buffer_1100_pdsd[325];

    auto g_y_x_0_0_x_xx_0_xz = buffer_1100_pdsd[326];

    auto g_y_x_0_0_x_xx_0_yy = buffer_1100_pdsd[327];

    auto g_y_x_0_0_x_xx_0_yz = buffer_1100_pdsd[328];

    auto g_y_x_0_0_x_xx_0_zz = buffer_1100_pdsd[329];

    auto g_y_x_0_0_x_xy_0_xx = buffer_1100_pdsd[330];

    auto g_y_x_0_0_x_xy_0_xy = buffer_1100_pdsd[331];

    auto g_y_x_0_0_x_xy_0_xz = buffer_1100_pdsd[332];

    auto g_y_x_0_0_x_xy_0_yy = buffer_1100_pdsd[333];

    auto g_y_x_0_0_x_xy_0_yz = buffer_1100_pdsd[334];

    auto g_y_x_0_0_x_xy_0_zz = buffer_1100_pdsd[335];

    auto g_y_x_0_0_x_xz_0_xx = buffer_1100_pdsd[336];

    auto g_y_x_0_0_x_xz_0_xy = buffer_1100_pdsd[337];

    auto g_y_x_0_0_x_xz_0_xz = buffer_1100_pdsd[338];

    auto g_y_x_0_0_x_xz_0_yy = buffer_1100_pdsd[339];

    auto g_y_x_0_0_x_xz_0_yz = buffer_1100_pdsd[340];

    auto g_y_x_0_0_x_xz_0_zz = buffer_1100_pdsd[341];

    auto g_y_x_0_0_x_yy_0_xx = buffer_1100_pdsd[342];

    auto g_y_x_0_0_x_yy_0_xy = buffer_1100_pdsd[343];

    auto g_y_x_0_0_x_yy_0_xz = buffer_1100_pdsd[344];

    auto g_y_x_0_0_x_yy_0_yy = buffer_1100_pdsd[345];

    auto g_y_x_0_0_x_yy_0_yz = buffer_1100_pdsd[346];

    auto g_y_x_0_0_x_yy_0_zz = buffer_1100_pdsd[347];

    auto g_y_x_0_0_x_yz_0_xx = buffer_1100_pdsd[348];

    auto g_y_x_0_0_x_yz_0_xy = buffer_1100_pdsd[349];

    auto g_y_x_0_0_x_yz_0_xz = buffer_1100_pdsd[350];

    auto g_y_x_0_0_x_yz_0_yy = buffer_1100_pdsd[351];

    auto g_y_x_0_0_x_yz_0_yz = buffer_1100_pdsd[352];

    auto g_y_x_0_0_x_yz_0_zz = buffer_1100_pdsd[353];

    auto g_y_x_0_0_x_zz_0_xx = buffer_1100_pdsd[354];

    auto g_y_x_0_0_x_zz_0_xy = buffer_1100_pdsd[355];

    auto g_y_x_0_0_x_zz_0_xz = buffer_1100_pdsd[356];

    auto g_y_x_0_0_x_zz_0_yy = buffer_1100_pdsd[357];

    auto g_y_x_0_0_x_zz_0_yz = buffer_1100_pdsd[358];

    auto g_y_x_0_0_x_zz_0_zz = buffer_1100_pdsd[359];

    auto g_y_x_0_0_y_xx_0_xx = buffer_1100_pdsd[360];

    auto g_y_x_0_0_y_xx_0_xy = buffer_1100_pdsd[361];

    auto g_y_x_0_0_y_xx_0_xz = buffer_1100_pdsd[362];

    auto g_y_x_0_0_y_xx_0_yy = buffer_1100_pdsd[363];

    auto g_y_x_0_0_y_xx_0_yz = buffer_1100_pdsd[364];

    auto g_y_x_0_0_y_xx_0_zz = buffer_1100_pdsd[365];

    auto g_y_x_0_0_y_xy_0_xx = buffer_1100_pdsd[366];

    auto g_y_x_0_0_y_xy_0_xy = buffer_1100_pdsd[367];

    auto g_y_x_0_0_y_xy_0_xz = buffer_1100_pdsd[368];

    auto g_y_x_0_0_y_xy_0_yy = buffer_1100_pdsd[369];

    auto g_y_x_0_0_y_xy_0_yz = buffer_1100_pdsd[370];

    auto g_y_x_0_0_y_xy_0_zz = buffer_1100_pdsd[371];

    auto g_y_x_0_0_y_xz_0_xx = buffer_1100_pdsd[372];

    auto g_y_x_0_0_y_xz_0_xy = buffer_1100_pdsd[373];

    auto g_y_x_0_0_y_xz_0_xz = buffer_1100_pdsd[374];

    auto g_y_x_0_0_y_xz_0_yy = buffer_1100_pdsd[375];

    auto g_y_x_0_0_y_xz_0_yz = buffer_1100_pdsd[376];

    auto g_y_x_0_0_y_xz_0_zz = buffer_1100_pdsd[377];

    auto g_y_x_0_0_y_yy_0_xx = buffer_1100_pdsd[378];

    auto g_y_x_0_0_y_yy_0_xy = buffer_1100_pdsd[379];

    auto g_y_x_0_0_y_yy_0_xz = buffer_1100_pdsd[380];

    auto g_y_x_0_0_y_yy_0_yy = buffer_1100_pdsd[381];

    auto g_y_x_0_0_y_yy_0_yz = buffer_1100_pdsd[382];

    auto g_y_x_0_0_y_yy_0_zz = buffer_1100_pdsd[383];

    auto g_y_x_0_0_y_yz_0_xx = buffer_1100_pdsd[384];

    auto g_y_x_0_0_y_yz_0_xy = buffer_1100_pdsd[385];

    auto g_y_x_0_0_y_yz_0_xz = buffer_1100_pdsd[386];

    auto g_y_x_0_0_y_yz_0_yy = buffer_1100_pdsd[387];

    auto g_y_x_0_0_y_yz_0_yz = buffer_1100_pdsd[388];

    auto g_y_x_0_0_y_yz_0_zz = buffer_1100_pdsd[389];

    auto g_y_x_0_0_y_zz_0_xx = buffer_1100_pdsd[390];

    auto g_y_x_0_0_y_zz_0_xy = buffer_1100_pdsd[391];

    auto g_y_x_0_0_y_zz_0_xz = buffer_1100_pdsd[392];

    auto g_y_x_0_0_y_zz_0_yy = buffer_1100_pdsd[393];

    auto g_y_x_0_0_y_zz_0_yz = buffer_1100_pdsd[394];

    auto g_y_x_0_0_y_zz_0_zz = buffer_1100_pdsd[395];

    auto g_y_x_0_0_z_xx_0_xx = buffer_1100_pdsd[396];

    auto g_y_x_0_0_z_xx_0_xy = buffer_1100_pdsd[397];

    auto g_y_x_0_0_z_xx_0_xz = buffer_1100_pdsd[398];

    auto g_y_x_0_0_z_xx_0_yy = buffer_1100_pdsd[399];

    auto g_y_x_0_0_z_xx_0_yz = buffer_1100_pdsd[400];

    auto g_y_x_0_0_z_xx_0_zz = buffer_1100_pdsd[401];

    auto g_y_x_0_0_z_xy_0_xx = buffer_1100_pdsd[402];

    auto g_y_x_0_0_z_xy_0_xy = buffer_1100_pdsd[403];

    auto g_y_x_0_0_z_xy_0_xz = buffer_1100_pdsd[404];

    auto g_y_x_0_0_z_xy_0_yy = buffer_1100_pdsd[405];

    auto g_y_x_0_0_z_xy_0_yz = buffer_1100_pdsd[406];

    auto g_y_x_0_0_z_xy_0_zz = buffer_1100_pdsd[407];

    auto g_y_x_0_0_z_xz_0_xx = buffer_1100_pdsd[408];

    auto g_y_x_0_0_z_xz_0_xy = buffer_1100_pdsd[409];

    auto g_y_x_0_0_z_xz_0_xz = buffer_1100_pdsd[410];

    auto g_y_x_0_0_z_xz_0_yy = buffer_1100_pdsd[411];

    auto g_y_x_0_0_z_xz_0_yz = buffer_1100_pdsd[412];

    auto g_y_x_0_0_z_xz_0_zz = buffer_1100_pdsd[413];

    auto g_y_x_0_0_z_yy_0_xx = buffer_1100_pdsd[414];

    auto g_y_x_0_0_z_yy_0_xy = buffer_1100_pdsd[415];

    auto g_y_x_0_0_z_yy_0_xz = buffer_1100_pdsd[416];

    auto g_y_x_0_0_z_yy_0_yy = buffer_1100_pdsd[417];

    auto g_y_x_0_0_z_yy_0_yz = buffer_1100_pdsd[418];

    auto g_y_x_0_0_z_yy_0_zz = buffer_1100_pdsd[419];

    auto g_y_x_0_0_z_yz_0_xx = buffer_1100_pdsd[420];

    auto g_y_x_0_0_z_yz_0_xy = buffer_1100_pdsd[421];

    auto g_y_x_0_0_z_yz_0_xz = buffer_1100_pdsd[422];

    auto g_y_x_0_0_z_yz_0_yy = buffer_1100_pdsd[423];

    auto g_y_x_0_0_z_yz_0_yz = buffer_1100_pdsd[424];

    auto g_y_x_0_0_z_yz_0_zz = buffer_1100_pdsd[425];

    auto g_y_x_0_0_z_zz_0_xx = buffer_1100_pdsd[426];

    auto g_y_x_0_0_z_zz_0_xy = buffer_1100_pdsd[427];

    auto g_y_x_0_0_z_zz_0_xz = buffer_1100_pdsd[428];

    auto g_y_x_0_0_z_zz_0_yy = buffer_1100_pdsd[429];

    auto g_y_x_0_0_z_zz_0_yz = buffer_1100_pdsd[430];

    auto g_y_x_0_0_z_zz_0_zz = buffer_1100_pdsd[431];

    auto g_y_y_0_0_x_xx_0_xx = buffer_1100_pdsd[432];

    auto g_y_y_0_0_x_xx_0_xy = buffer_1100_pdsd[433];

    auto g_y_y_0_0_x_xx_0_xz = buffer_1100_pdsd[434];

    auto g_y_y_0_0_x_xx_0_yy = buffer_1100_pdsd[435];

    auto g_y_y_0_0_x_xx_0_yz = buffer_1100_pdsd[436];

    auto g_y_y_0_0_x_xx_0_zz = buffer_1100_pdsd[437];

    auto g_y_y_0_0_x_xy_0_xx = buffer_1100_pdsd[438];

    auto g_y_y_0_0_x_xy_0_xy = buffer_1100_pdsd[439];

    auto g_y_y_0_0_x_xy_0_xz = buffer_1100_pdsd[440];

    auto g_y_y_0_0_x_xy_0_yy = buffer_1100_pdsd[441];

    auto g_y_y_0_0_x_xy_0_yz = buffer_1100_pdsd[442];

    auto g_y_y_0_0_x_xy_0_zz = buffer_1100_pdsd[443];

    auto g_y_y_0_0_x_xz_0_xx = buffer_1100_pdsd[444];

    auto g_y_y_0_0_x_xz_0_xy = buffer_1100_pdsd[445];

    auto g_y_y_0_0_x_xz_0_xz = buffer_1100_pdsd[446];

    auto g_y_y_0_0_x_xz_0_yy = buffer_1100_pdsd[447];

    auto g_y_y_0_0_x_xz_0_yz = buffer_1100_pdsd[448];

    auto g_y_y_0_0_x_xz_0_zz = buffer_1100_pdsd[449];

    auto g_y_y_0_0_x_yy_0_xx = buffer_1100_pdsd[450];

    auto g_y_y_0_0_x_yy_0_xy = buffer_1100_pdsd[451];

    auto g_y_y_0_0_x_yy_0_xz = buffer_1100_pdsd[452];

    auto g_y_y_0_0_x_yy_0_yy = buffer_1100_pdsd[453];

    auto g_y_y_0_0_x_yy_0_yz = buffer_1100_pdsd[454];

    auto g_y_y_0_0_x_yy_0_zz = buffer_1100_pdsd[455];

    auto g_y_y_0_0_x_yz_0_xx = buffer_1100_pdsd[456];

    auto g_y_y_0_0_x_yz_0_xy = buffer_1100_pdsd[457];

    auto g_y_y_0_0_x_yz_0_xz = buffer_1100_pdsd[458];

    auto g_y_y_0_0_x_yz_0_yy = buffer_1100_pdsd[459];

    auto g_y_y_0_0_x_yz_0_yz = buffer_1100_pdsd[460];

    auto g_y_y_0_0_x_yz_0_zz = buffer_1100_pdsd[461];

    auto g_y_y_0_0_x_zz_0_xx = buffer_1100_pdsd[462];

    auto g_y_y_0_0_x_zz_0_xy = buffer_1100_pdsd[463];

    auto g_y_y_0_0_x_zz_0_xz = buffer_1100_pdsd[464];

    auto g_y_y_0_0_x_zz_0_yy = buffer_1100_pdsd[465];

    auto g_y_y_0_0_x_zz_0_yz = buffer_1100_pdsd[466];

    auto g_y_y_0_0_x_zz_0_zz = buffer_1100_pdsd[467];

    auto g_y_y_0_0_y_xx_0_xx = buffer_1100_pdsd[468];

    auto g_y_y_0_0_y_xx_0_xy = buffer_1100_pdsd[469];

    auto g_y_y_0_0_y_xx_0_xz = buffer_1100_pdsd[470];

    auto g_y_y_0_0_y_xx_0_yy = buffer_1100_pdsd[471];

    auto g_y_y_0_0_y_xx_0_yz = buffer_1100_pdsd[472];

    auto g_y_y_0_0_y_xx_0_zz = buffer_1100_pdsd[473];

    auto g_y_y_0_0_y_xy_0_xx = buffer_1100_pdsd[474];

    auto g_y_y_0_0_y_xy_0_xy = buffer_1100_pdsd[475];

    auto g_y_y_0_0_y_xy_0_xz = buffer_1100_pdsd[476];

    auto g_y_y_0_0_y_xy_0_yy = buffer_1100_pdsd[477];

    auto g_y_y_0_0_y_xy_0_yz = buffer_1100_pdsd[478];

    auto g_y_y_0_0_y_xy_0_zz = buffer_1100_pdsd[479];

    auto g_y_y_0_0_y_xz_0_xx = buffer_1100_pdsd[480];

    auto g_y_y_0_0_y_xz_0_xy = buffer_1100_pdsd[481];

    auto g_y_y_0_0_y_xz_0_xz = buffer_1100_pdsd[482];

    auto g_y_y_0_0_y_xz_0_yy = buffer_1100_pdsd[483];

    auto g_y_y_0_0_y_xz_0_yz = buffer_1100_pdsd[484];

    auto g_y_y_0_0_y_xz_0_zz = buffer_1100_pdsd[485];

    auto g_y_y_0_0_y_yy_0_xx = buffer_1100_pdsd[486];

    auto g_y_y_0_0_y_yy_0_xy = buffer_1100_pdsd[487];

    auto g_y_y_0_0_y_yy_0_xz = buffer_1100_pdsd[488];

    auto g_y_y_0_0_y_yy_0_yy = buffer_1100_pdsd[489];

    auto g_y_y_0_0_y_yy_0_yz = buffer_1100_pdsd[490];

    auto g_y_y_0_0_y_yy_0_zz = buffer_1100_pdsd[491];

    auto g_y_y_0_0_y_yz_0_xx = buffer_1100_pdsd[492];

    auto g_y_y_0_0_y_yz_0_xy = buffer_1100_pdsd[493];

    auto g_y_y_0_0_y_yz_0_xz = buffer_1100_pdsd[494];

    auto g_y_y_0_0_y_yz_0_yy = buffer_1100_pdsd[495];

    auto g_y_y_0_0_y_yz_0_yz = buffer_1100_pdsd[496];

    auto g_y_y_0_0_y_yz_0_zz = buffer_1100_pdsd[497];

    auto g_y_y_0_0_y_zz_0_xx = buffer_1100_pdsd[498];

    auto g_y_y_0_0_y_zz_0_xy = buffer_1100_pdsd[499];

    auto g_y_y_0_0_y_zz_0_xz = buffer_1100_pdsd[500];

    auto g_y_y_0_0_y_zz_0_yy = buffer_1100_pdsd[501];

    auto g_y_y_0_0_y_zz_0_yz = buffer_1100_pdsd[502];

    auto g_y_y_0_0_y_zz_0_zz = buffer_1100_pdsd[503];

    auto g_y_y_0_0_z_xx_0_xx = buffer_1100_pdsd[504];

    auto g_y_y_0_0_z_xx_0_xy = buffer_1100_pdsd[505];

    auto g_y_y_0_0_z_xx_0_xz = buffer_1100_pdsd[506];

    auto g_y_y_0_0_z_xx_0_yy = buffer_1100_pdsd[507];

    auto g_y_y_0_0_z_xx_0_yz = buffer_1100_pdsd[508];

    auto g_y_y_0_0_z_xx_0_zz = buffer_1100_pdsd[509];

    auto g_y_y_0_0_z_xy_0_xx = buffer_1100_pdsd[510];

    auto g_y_y_0_0_z_xy_0_xy = buffer_1100_pdsd[511];

    auto g_y_y_0_0_z_xy_0_xz = buffer_1100_pdsd[512];

    auto g_y_y_0_0_z_xy_0_yy = buffer_1100_pdsd[513];

    auto g_y_y_0_0_z_xy_0_yz = buffer_1100_pdsd[514];

    auto g_y_y_0_0_z_xy_0_zz = buffer_1100_pdsd[515];

    auto g_y_y_0_0_z_xz_0_xx = buffer_1100_pdsd[516];

    auto g_y_y_0_0_z_xz_0_xy = buffer_1100_pdsd[517];

    auto g_y_y_0_0_z_xz_0_xz = buffer_1100_pdsd[518];

    auto g_y_y_0_0_z_xz_0_yy = buffer_1100_pdsd[519];

    auto g_y_y_0_0_z_xz_0_yz = buffer_1100_pdsd[520];

    auto g_y_y_0_0_z_xz_0_zz = buffer_1100_pdsd[521];

    auto g_y_y_0_0_z_yy_0_xx = buffer_1100_pdsd[522];

    auto g_y_y_0_0_z_yy_0_xy = buffer_1100_pdsd[523];

    auto g_y_y_0_0_z_yy_0_xz = buffer_1100_pdsd[524];

    auto g_y_y_0_0_z_yy_0_yy = buffer_1100_pdsd[525];

    auto g_y_y_0_0_z_yy_0_yz = buffer_1100_pdsd[526];

    auto g_y_y_0_0_z_yy_0_zz = buffer_1100_pdsd[527];

    auto g_y_y_0_0_z_yz_0_xx = buffer_1100_pdsd[528];

    auto g_y_y_0_0_z_yz_0_xy = buffer_1100_pdsd[529];

    auto g_y_y_0_0_z_yz_0_xz = buffer_1100_pdsd[530];

    auto g_y_y_0_0_z_yz_0_yy = buffer_1100_pdsd[531];

    auto g_y_y_0_0_z_yz_0_yz = buffer_1100_pdsd[532];

    auto g_y_y_0_0_z_yz_0_zz = buffer_1100_pdsd[533];

    auto g_y_y_0_0_z_zz_0_xx = buffer_1100_pdsd[534];

    auto g_y_y_0_0_z_zz_0_xy = buffer_1100_pdsd[535];

    auto g_y_y_0_0_z_zz_0_xz = buffer_1100_pdsd[536];

    auto g_y_y_0_0_z_zz_0_yy = buffer_1100_pdsd[537];

    auto g_y_y_0_0_z_zz_0_yz = buffer_1100_pdsd[538];

    auto g_y_y_0_0_z_zz_0_zz = buffer_1100_pdsd[539];

    auto g_y_z_0_0_x_xx_0_xx = buffer_1100_pdsd[540];

    auto g_y_z_0_0_x_xx_0_xy = buffer_1100_pdsd[541];

    auto g_y_z_0_0_x_xx_0_xz = buffer_1100_pdsd[542];

    auto g_y_z_0_0_x_xx_0_yy = buffer_1100_pdsd[543];

    auto g_y_z_0_0_x_xx_0_yz = buffer_1100_pdsd[544];

    auto g_y_z_0_0_x_xx_0_zz = buffer_1100_pdsd[545];

    auto g_y_z_0_0_x_xy_0_xx = buffer_1100_pdsd[546];

    auto g_y_z_0_0_x_xy_0_xy = buffer_1100_pdsd[547];

    auto g_y_z_0_0_x_xy_0_xz = buffer_1100_pdsd[548];

    auto g_y_z_0_0_x_xy_0_yy = buffer_1100_pdsd[549];

    auto g_y_z_0_0_x_xy_0_yz = buffer_1100_pdsd[550];

    auto g_y_z_0_0_x_xy_0_zz = buffer_1100_pdsd[551];

    auto g_y_z_0_0_x_xz_0_xx = buffer_1100_pdsd[552];

    auto g_y_z_0_0_x_xz_0_xy = buffer_1100_pdsd[553];

    auto g_y_z_0_0_x_xz_0_xz = buffer_1100_pdsd[554];

    auto g_y_z_0_0_x_xz_0_yy = buffer_1100_pdsd[555];

    auto g_y_z_0_0_x_xz_0_yz = buffer_1100_pdsd[556];

    auto g_y_z_0_0_x_xz_0_zz = buffer_1100_pdsd[557];

    auto g_y_z_0_0_x_yy_0_xx = buffer_1100_pdsd[558];

    auto g_y_z_0_0_x_yy_0_xy = buffer_1100_pdsd[559];

    auto g_y_z_0_0_x_yy_0_xz = buffer_1100_pdsd[560];

    auto g_y_z_0_0_x_yy_0_yy = buffer_1100_pdsd[561];

    auto g_y_z_0_0_x_yy_0_yz = buffer_1100_pdsd[562];

    auto g_y_z_0_0_x_yy_0_zz = buffer_1100_pdsd[563];

    auto g_y_z_0_0_x_yz_0_xx = buffer_1100_pdsd[564];

    auto g_y_z_0_0_x_yz_0_xy = buffer_1100_pdsd[565];

    auto g_y_z_0_0_x_yz_0_xz = buffer_1100_pdsd[566];

    auto g_y_z_0_0_x_yz_0_yy = buffer_1100_pdsd[567];

    auto g_y_z_0_0_x_yz_0_yz = buffer_1100_pdsd[568];

    auto g_y_z_0_0_x_yz_0_zz = buffer_1100_pdsd[569];

    auto g_y_z_0_0_x_zz_0_xx = buffer_1100_pdsd[570];

    auto g_y_z_0_0_x_zz_0_xy = buffer_1100_pdsd[571];

    auto g_y_z_0_0_x_zz_0_xz = buffer_1100_pdsd[572];

    auto g_y_z_0_0_x_zz_0_yy = buffer_1100_pdsd[573];

    auto g_y_z_0_0_x_zz_0_yz = buffer_1100_pdsd[574];

    auto g_y_z_0_0_x_zz_0_zz = buffer_1100_pdsd[575];

    auto g_y_z_0_0_y_xx_0_xx = buffer_1100_pdsd[576];

    auto g_y_z_0_0_y_xx_0_xy = buffer_1100_pdsd[577];

    auto g_y_z_0_0_y_xx_0_xz = buffer_1100_pdsd[578];

    auto g_y_z_0_0_y_xx_0_yy = buffer_1100_pdsd[579];

    auto g_y_z_0_0_y_xx_0_yz = buffer_1100_pdsd[580];

    auto g_y_z_0_0_y_xx_0_zz = buffer_1100_pdsd[581];

    auto g_y_z_0_0_y_xy_0_xx = buffer_1100_pdsd[582];

    auto g_y_z_0_0_y_xy_0_xy = buffer_1100_pdsd[583];

    auto g_y_z_0_0_y_xy_0_xz = buffer_1100_pdsd[584];

    auto g_y_z_0_0_y_xy_0_yy = buffer_1100_pdsd[585];

    auto g_y_z_0_0_y_xy_0_yz = buffer_1100_pdsd[586];

    auto g_y_z_0_0_y_xy_0_zz = buffer_1100_pdsd[587];

    auto g_y_z_0_0_y_xz_0_xx = buffer_1100_pdsd[588];

    auto g_y_z_0_0_y_xz_0_xy = buffer_1100_pdsd[589];

    auto g_y_z_0_0_y_xz_0_xz = buffer_1100_pdsd[590];

    auto g_y_z_0_0_y_xz_0_yy = buffer_1100_pdsd[591];

    auto g_y_z_0_0_y_xz_0_yz = buffer_1100_pdsd[592];

    auto g_y_z_0_0_y_xz_0_zz = buffer_1100_pdsd[593];

    auto g_y_z_0_0_y_yy_0_xx = buffer_1100_pdsd[594];

    auto g_y_z_0_0_y_yy_0_xy = buffer_1100_pdsd[595];

    auto g_y_z_0_0_y_yy_0_xz = buffer_1100_pdsd[596];

    auto g_y_z_0_0_y_yy_0_yy = buffer_1100_pdsd[597];

    auto g_y_z_0_0_y_yy_0_yz = buffer_1100_pdsd[598];

    auto g_y_z_0_0_y_yy_0_zz = buffer_1100_pdsd[599];

    auto g_y_z_0_0_y_yz_0_xx = buffer_1100_pdsd[600];

    auto g_y_z_0_0_y_yz_0_xy = buffer_1100_pdsd[601];

    auto g_y_z_0_0_y_yz_0_xz = buffer_1100_pdsd[602];

    auto g_y_z_0_0_y_yz_0_yy = buffer_1100_pdsd[603];

    auto g_y_z_0_0_y_yz_0_yz = buffer_1100_pdsd[604];

    auto g_y_z_0_0_y_yz_0_zz = buffer_1100_pdsd[605];

    auto g_y_z_0_0_y_zz_0_xx = buffer_1100_pdsd[606];

    auto g_y_z_0_0_y_zz_0_xy = buffer_1100_pdsd[607];

    auto g_y_z_0_0_y_zz_0_xz = buffer_1100_pdsd[608];

    auto g_y_z_0_0_y_zz_0_yy = buffer_1100_pdsd[609];

    auto g_y_z_0_0_y_zz_0_yz = buffer_1100_pdsd[610];

    auto g_y_z_0_0_y_zz_0_zz = buffer_1100_pdsd[611];

    auto g_y_z_0_0_z_xx_0_xx = buffer_1100_pdsd[612];

    auto g_y_z_0_0_z_xx_0_xy = buffer_1100_pdsd[613];

    auto g_y_z_0_0_z_xx_0_xz = buffer_1100_pdsd[614];

    auto g_y_z_0_0_z_xx_0_yy = buffer_1100_pdsd[615];

    auto g_y_z_0_0_z_xx_0_yz = buffer_1100_pdsd[616];

    auto g_y_z_0_0_z_xx_0_zz = buffer_1100_pdsd[617];

    auto g_y_z_0_0_z_xy_0_xx = buffer_1100_pdsd[618];

    auto g_y_z_0_0_z_xy_0_xy = buffer_1100_pdsd[619];

    auto g_y_z_0_0_z_xy_0_xz = buffer_1100_pdsd[620];

    auto g_y_z_0_0_z_xy_0_yy = buffer_1100_pdsd[621];

    auto g_y_z_0_0_z_xy_0_yz = buffer_1100_pdsd[622];

    auto g_y_z_0_0_z_xy_0_zz = buffer_1100_pdsd[623];

    auto g_y_z_0_0_z_xz_0_xx = buffer_1100_pdsd[624];

    auto g_y_z_0_0_z_xz_0_xy = buffer_1100_pdsd[625];

    auto g_y_z_0_0_z_xz_0_xz = buffer_1100_pdsd[626];

    auto g_y_z_0_0_z_xz_0_yy = buffer_1100_pdsd[627];

    auto g_y_z_0_0_z_xz_0_yz = buffer_1100_pdsd[628];

    auto g_y_z_0_0_z_xz_0_zz = buffer_1100_pdsd[629];

    auto g_y_z_0_0_z_yy_0_xx = buffer_1100_pdsd[630];

    auto g_y_z_0_0_z_yy_0_xy = buffer_1100_pdsd[631];

    auto g_y_z_0_0_z_yy_0_xz = buffer_1100_pdsd[632];

    auto g_y_z_0_0_z_yy_0_yy = buffer_1100_pdsd[633];

    auto g_y_z_0_0_z_yy_0_yz = buffer_1100_pdsd[634];

    auto g_y_z_0_0_z_yy_0_zz = buffer_1100_pdsd[635];

    auto g_y_z_0_0_z_yz_0_xx = buffer_1100_pdsd[636];

    auto g_y_z_0_0_z_yz_0_xy = buffer_1100_pdsd[637];

    auto g_y_z_0_0_z_yz_0_xz = buffer_1100_pdsd[638];

    auto g_y_z_0_0_z_yz_0_yy = buffer_1100_pdsd[639];

    auto g_y_z_0_0_z_yz_0_yz = buffer_1100_pdsd[640];

    auto g_y_z_0_0_z_yz_0_zz = buffer_1100_pdsd[641];

    auto g_y_z_0_0_z_zz_0_xx = buffer_1100_pdsd[642];

    auto g_y_z_0_0_z_zz_0_xy = buffer_1100_pdsd[643];

    auto g_y_z_0_0_z_zz_0_xz = buffer_1100_pdsd[644];

    auto g_y_z_0_0_z_zz_0_yy = buffer_1100_pdsd[645];

    auto g_y_z_0_0_z_zz_0_yz = buffer_1100_pdsd[646];

    auto g_y_z_0_0_z_zz_0_zz = buffer_1100_pdsd[647];

    auto g_z_x_0_0_x_xx_0_xx = buffer_1100_pdsd[648];

    auto g_z_x_0_0_x_xx_0_xy = buffer_1100_pdsd[649];

    auto g_z_x_0_0_x_xx_0_xz = buffer_1100_pdsd[650];

    auto g_z_x_0_0_x_xx_0_yy = buffer_1100_pdsd[651];

    auto g_z_x_0_0_x_xx_0_yz = buffer_1100_pdsd[652];

    auto g_z_x_0_0_x_xx_0_zz = buffer_1100_pdsd[653];

    auto g_z_x_0_0_x_xy_0_xx = buffer_1100_pdsd[654];

    auto g_z_x_0_0_x_xy_0_xy = buffer_1100_pdsd[655];

    auto g_z_x_0_0_x_xy_0_xz = buffer_1100_pdsd[656];

    auto g_z_x_0_0_x_xy_0_yy = buffer_1100_pdsd[657];

    auto g_z_x_0_0_x_xy_0_yz = buffer_1100_pdsd[658];

    auto g_z_x_0_0_x_xy_0_zz = buffer_1100_pdsd[659];

    auto g_z_x_0_0_x_xz_0_xx = buffer_1100_pdsd[660];

    auto g_z_x_0_0_x_xz_0_xy = buffer_1100_pdsd[661];

    auto g_z_x_0_0_x_xz_0_xz = buffer_1100_pdsd[662];

    auto g_z_x_0_0_x_xz_0_yy = buffer_1100_pdsd[663];

    auto g_z_x_0_0_x_xz_0_yz = buffer_1100_pdsd[664];

    auto g_z_x_0_0_x_xz_0_zz = buffer_1100_pdsd[665];

    auto g_z_x_0_0_x_yy_0_xx = buffer_1100_pdsd[666];

    auto g_z_x_0_0_x_yy_0_xy = buffer_1100_pdsd[667];

    auto g_z_x_0_0_x_yy_0_xz = buffer_1100_pdsd[668];

    auto g_z_x_0_0_x_yy_0_yy = buffer_1100_pdsd[669];

    auto g_z_x_0_0_x_yy_0_yz = buffer_1100_pdsd[670];

    auto g_z_x_0_0_x_yy_0_zz = buffer_1100_pdsd[671];

    auto g_z_x_0_0_x_yz_0_xx = buffer_1100_pdsd[672];

    auto g_z_x_0_0_x_yz_0_xy = buffer_1100_pdsd[673];

    auto g_z_x_0_0_x_yz_0_xz = buffer_1100_pdsd[674];

    auto g_z_x_0_0_x_yz_0_yy = buffer_1100_pdsd[675];

    auto g_z_x_0_0_x_yz_0_yz = buffer_1100_pdsd[676];

    auto g_z_x_0_0_x_yz_0_zz = buffer_1100_pdsd[677];

    auto g_z_x_0_0_x_zz_0_xx = buffer_1100_pdsd[678];

    auto g_z_x_0_0_x_zz_0_xy = buffer_1100_pdsd[679];

    auto g_z_x_0_0_x_zz_0_xz = buffer_1100_pdsd[680];

    auto g_z_x_0_0_x_zz_0_yy = buffer_1100_pdsd[681];

    auto g_z_x_0_0_x_zz_0_yz = buffer_1100_pdsd[682];

    auto g_z_x_0_0_x_zz_0_zz = buffer_1100_pdsd[683];

    auto g_z_x_0_0_y_xx_0_xx = buffer_1100_pdsd[684];

    auto g_z_x_0_0_y_xx_0_xy = buffer_1100_pdsd[685];

    auto g_z_x_0_0_y_xx_0_xz = buffer_1100_pdsd[686];

    auto g_z_x_0_0_y_xx_0_yy = buffer_1100_pdsd[687];

    auto g_z_x_0_0_y_xx_0_yz = buffer_1100_pdsd[688];

    auto g_z_x_0_0_y_xx_0_zz = buffer_1100_pdsd[689];

    auto g_z_x_0_0_y_xy_0_xx = buffer_1100_pdsd[690];

    auto g_z_x_0_0_y_xy_0_xy = buffer_1100_pdsd[691];

    auto g_z_x_0_0_y_xy_0_xz = buffer_1100_pdsd[692];

    auto g_z_x_0_0_y_xy_0_yy = buffer_1100_pdsd[693];

    auto g_z_x_0_0_y_xy_0_yz = buffer_1100_pdsd[694];

    auto g_z_x_0_0_y_xy_0_zz = buffer_1100_pdsd[695];

    auto g_z_x_0_0_y_xz_0_xx = buffer_1100_pdsd[696];

    auto g_z_x_0_0_y_xz_0_xy = buffer_1100_pdsd[697];

    auto g_z_x_0_0_y_xz_0_xz = buffer_1100_pdsd[698];

    auto g_z_x_0_0_y_xz_0_yy = buffer_1100_pdsd[699];

    auto g_z_x_0_0_y_xz_0_yz = buffer_1100_pdsd[700];

    auto g_z_x_0_0_y_xz_0_zz = buffer_1100_pdsd[701];

    auto g_z_x_0_0_y_yy_0_xx = buffer_1100_pdsd[702];

    auto g_z_x_0_0_y_yy_0_xy = buffer_1100_pdsd[703];

    auto g_z_x_0_0_y_yy_0_xz = buffer_1100_pdsd[704];

    auto g_z_x_0_0_y_yy_0_yy = buffer_1100_pdsd[705];

    auto g_z_x_0_0_y_yy_0_yz = buffer_1100_pdsd[706];

    auto g_z_x_0_0_y_yy_0_zz = buffer_1100_pdsd[707];

    auto g_z_x_0_0_y_yz_0_xx = buffer_1100_pdsd[708];

    auto g_z_x_0_0_y_yz_0_xy = buffer_1100_pdsd[709];

    auto g_z_x_0_0_y_yz_0_xz = buffer_1100_pdsd[710];

    auto g_z_x_0_0_y_yz_0_yy = buffer_1100_pdsd[711];

    auto g_z_x_0_0_y_yz_0_yz = buffer_1100_pdsd[712];

    auto g_z_x_0_0_y_yz_0_zz = buffer_1100_pdsd[713];

    auto g_z_x_0_0_y_zz_0_xx = buffer_1100_pdsd[714];

    auto g_z_x_0_0_y_zz_0_xy = buffer_1100_pdsd[715];

    auto g_z_x_0_0_y_zz_0_xz = buffer_1100_pdsd[716];

    auto g_z_x_0_0_y_zz_0_yy = buffer_1100_pdsd[717];

    auto g_z_x_0_0_y_zz_0_yz = buffer_1100_pdsd[718];

    auto g_z_x_0_0_y_zz_0_zz = buffer_1100_pdsd[719];

    auto g_z_x_0_0_z_xx_0_xx = buffer_1100_pdsd[720];

    auto g_z_x_0_0_z_xx_0_xy = buffer_1100_pdsd[721];

    auto g_z_x_0_0_z_xx_0_xz = buffer_1100_pdsd[722];

    auto g_z_x_0_0_z_xx_0_yy = buffer_1100_pdsd[723];

    auto g_z_x_0_0_z_xx_0_yz = buffer_1100_pdsd[724];

    auto g_z_x_0_0_z_xx_0_zz = buffer_1100_pdsd[725];

    auto g_z_x_0_0_z_xy_0_xx = buffer_1100_pdsd[726];

    auto g_z_x_0_0_z_xy_0_xy = buffer_1100_pdsd[727];

    auto g_z_x_0_0_z_xy_0_xz = buffer_1100_pdsd[728];

    auto g_z_x_0_0_z_xy_0_yy = buffer_1100_pdsd[729];

    auto g_z_x_0_0_z_xy_0_yz = buffer_1100_pdsd[730];

    auto g_z_x_0_0_z_xy_0_zz = buffer_1100_pdsd[731];

    auto g_z_x_0_0_z_xz_0_xx = buffer_1100_pdsd[732];

    auto g_z_x_0_0_z_xz_0_xy = buffer_1100_pdsd[733];

    auto g_z_x_0_0_z_xz_0_xz = buffer_1100_pdsd[734];

    auto g_z_x_0_0_z_xz_0_yy = buffer_1100_pdsd[735];

    auto g_z_x_0_0_z_xz_0_yz = buffer_1100_pdsd[736];

    auto g_z_x_0_0_z_xz_0_zz = buffer_1100_pdsd[737];

    auto g_z_x_0_0_z_yy_0_xx = buffer_1100_pdsd[738];

    auto g_z_x_0_0_z_yy_0_xy = buffer_1100_pdsd[739];

    auto g_z_x_0_0_z_yy_0_xz = buffer_1100_pdsd[740];

    auto g_z_x_0_0_z_yy_0_yy = buffer_1100_pdsd[741];

    auto g_z_x_0_0_z_yy_0_yz = buffer_1100_pdsd[742];

    auto g_z_x_0_0_z_yy_0_zz = buffer_1100_pdsd[743];

    auto g_z_x_0_0_z_yz_0_xx = buffer_1100_pdsd[744];

    auto g_z_x_0_0_z_yz_0_xy = buffer_1100_pdsd[745];

    auto g_z_x_0_0_z_yz_0_xz = buffer_1100_pdsd[746];

    auto g_z_x_0_0_z_yz_0_yy = buffer_1100_pdsd[747];

    auto g_z_x_0_0_z_yz_0_yz = buffer_1100_pdsd[748];

    auto g_z_x_0_0_z_yz_0_zz = buffer_1100_pdsd[749];

    auto g_z_x_0_0_z_zz_0_xx = buffer_1100_pdsd[750];

    auto g_z_x_0_0_z_zz_0_xy = buffer_1100_pdsd[751];

    auto g_z_x_0_0_z_zz_0_xz = buffer_1100_pdsd[752];

    auto g_z_x_0_0_z_zz_0_yy = buffer_1100_pdsd[753];

    auto g_z_x_0_0_z_zz_0_yz = buffer_1100_pdsd[754];

    auto g_z_x_0_0_z_zz_0_zz = buffer_1100_pdsd[755];

    auto g_z_y_0_0_x_xx_0_xx = buffer_1100_pdsd[756];

    auto g_z_y_0_0_x_xx_0_xy = buffer_1100_pdsd[757];

    auto g_z_y_0_0_x_xx_0_xz = buffer_1100_pdsd[758];

    auto g_z_y_0_0_x_xx_0_yy = buffer_1100_pdsd[759];

    auto g_z_y_0_0_x_xx_0_yz = buffer_1100_pdsd[760];

    auto g_z_y_0_0_x_xx_0_zz = buffer_1100_pdsd[761];

    auto g_z_y_0_0_x_xy_0_xx = buffer_1100_pdsd[762];

    auto g_z_y_0_0_x_xy_0_xy = buffer_1100_pdsd[763];

    auto g_z_y_0_0_x_xy_0_xz = buffer_1100_pdsd[764];

    auto g_z_y_0_0_x_xy_0_yy = buffer_1100_pdsd[765];

    auto g_z_y_0_0_x_xy_0_yz = buffer_1100_pdsd[766];

    auto g_z_y_0_0_x_xy_0_zz = buffer_1100_pdsd[767];

    auto g_z_y_0_0_x_xz_0_xx = buffer_1100_pdsd[768];

    auto g_z_y_0_0_x_xz_0_xy = buffer_1100_pdsd[769];

    auto g_z_y_0_0_x_xz_0_xz = buffer_1100_pdsd[770];

    auto g_z_y_0_0_x_xz_0_yy = buffer_1100_pdsd[771];

    auto g_z_y_0_0_x_xz_0_yz = buffer_1100_pdsd[772];

    auto g_z_y_0_0_x_xz_0_zz = buffer_1100_pdsd[773];

    auto g_z_y_0_0_x_yy_0_xx = buffer_1100_pdsd[774];

    auto g_z_y_0_0_x_yy_0_xy = buffer_1100_pdsd[775];

    auto g_z_y_0_0_x_yy_0_xz = buffer_1100_pdsd[776];

    auto g_z_y_0_0_x_yy_0_yy = buffer_1100_pdsd[777];

    auto g_z_y_0_0_x_yy_0_yz = buffer_1100_pdsd[778];

    auto g_z_y_0_0_x_yy_0_zz = buffer_1100_pdsd[779];

    auto g_z_y_0_0_x_yz_0_xx = buffer_1100_pdsd[780];

    auto g_z_y_0_0_x_yz_0_xy = buffer_1100_pdsd[781];

    auto g_z_y_0_0_x_yz_0_xz = buffer_1100_pdsd[782];

    auto g_z_y_0_0_x_yz_0_yy = buffer_1100_pdsd[783];

    auto g_z_y_0_0_x_yz_0_yz = buffer_1100_pdsd[784];

    auto g_z_y_0_0_x_yz_0_zz = buffer_1100_pdsd[785];

    auto g_z_y_0_0_x_zz_0_xx = buffer_1100_pdsd[786];

    auto g_z_y_0_0_x_zz_0_xy = buffer_1100_pdsd[787];

    auto g_z_y_0_0_x_zz_0_xz = buffer_1100_pdsd[788];

    auto g_z_y_0_0_x_zz_0_yy = buffer_1100_pdsd[789];

    auto g_z_y_0_0_x_zz_0_yz = buffer_1100_pdsd[790];

    auto g_z_y_0_0_x_zz_0_zz = buffer_1100_pdsd[791];

    auto g_z_y_0_0_y_xx_0_xx = buffer_1100_pdsd[792];

    auto g_z_y_0_0_y_xx_0_xy = buffer_1100_pdsd[793];

    auto g_z_y_0_0_y_xx_0_xz = buffer_1100_pdsd[794];

    auto g_z_y_0_0_y_xx_0_yy = buffer_1100_pdsd[795];

    auto g_z_y_0_0_y_xx_0_yz = buffer_1100_pdsd[796];

    auto g_z_y_0_0_y_xx_0_zz = buffer_1100_pdsd[797];

    auto g_z_y_0_0_y_xy_0_xx = buffer_1100_pdsd[798];

    auto g_z_y_0_0_y_xy_0_xy = buffer_1100_pdsd[799];

    auto g_z_y_0_0_y_xy_0_xz = buffer_1100_pdsd[800];

    auto g_z_y_0_0_y_xy_0_yy = buffer_1100_pdsd[801];

    auto g_z_y_0_0_y_xy_0_yz = buffer_1100_pdsd[802];

    auto g_z_y_0_0_y_xy_0_zz = buffer_1100_pdsd[803];

    auto g_z_y_0_0_y_xz_0_xx = buffer_1100_pdsd[804];

    auto g_z_y_0_0_y_xz_0_xy = buffer_1100_pdsd[805];

    auto g_z_y_0_0_y_xz_0_xz = buffer_1100_pdsd[806];

    auto g_z_y_0_0_y_xz_0_yy = buffer_1100_pdsd[807];

    auto g_z_y_0_0_y_xz_0_yz = buffer_1100_pdsd[808];

    auto g_z_y_0_0_y_xz_0_zz = buffer_1100_pdsd[809];

    auto g_z_y_0_0_y_yy_0_xx = buffer_1100_pdsd[810];

    auto g_z_y_0_0_y_yy_0_xy = buffer_1100_pdsd[811];

    auto g_z_y_0_0_y_yy_0_xz = buffer_1100_pdsd[812];

    auto g_z_y_0_0_y_yy_0_yy = buffer_1100_pdsd[813];

    auto g_z_y_0_0_y_yy_0_yz = buffer_1100_pdsd[814];

    auto g_z_y_0_0_y_yy_0_zz = buffer_1100_pdsd[815];

    auto g_z_y_0_0_y_yz_0_xx = buffer_1100_pdsd[816];

    auto g_z_y_0_0_y_yz_0_xy = buffer_1100_pdsd[817];

    auto g_z_y_0_0_y_yz_0_xz = buffer_1100_pdsd[818];

    auto g_z_y_0_0_y_yz_0_yy = buffer_1100_pdsd[819];

    auto g_z_y_0_0_y_yz_0_yz = buffer_1100_pdsd[820];

    auto g_z_y_0_0_y_yz_0_zz = buffer_1100_pdsd[821];

    auto g_z_y_0_0_y_zz_0_xx = buffer_1100_pdsd[822];

    auto g_z_y_0_0_y_zz_0_xy = buffer_1100_pdsd[823];

    auto g_z_y_0_0_y_zz_0_xz = buffer_1100_pdsd[824];

    auto g_z_y_0_0_y_zz_0_yy = buffer_1100_pdsd[825];

    auto g_z_y_0_0_y_zz_0_yz = buffer_1100_pdsd[826];

    auto g_z_y_0_0_y_zz_0_zz = buffer_1100_pdsd[827];

    auto g_z_y_0_0_z_xx_0_xx = buffer_1100_pdsd[828];

    auto g_z_y_0_0_z_xx_0_xy = buffer_1100_pdsd[829];

    auto g_z_y_0_0_z_xx_0_xz = buffer_1100_pdsd[830];

    auto g_z_y_0_0_z_xx_0_yy = buffer_1100_pdsd[831];

    auto g_z_y_0_0_z_xx_0_yz = buffer_1100_pdsd[832];

    auto g_z_y_0_0_z_xx_0_zz = buffer_1100_pdsd[833];

    auto g_z_y_0_0_z_xy_0_xx = buffer_1100_pdsd[834];

    auto g_z_y_0_0_z_xy_0_xy = buffer_1100_pdsd[835];

    auto g_z_y_0_0_z_xy_0_xz = buffer_1100_pdsd[836];

    auto g_z_y_0_0_z_xy_0_yy = buffer_1100_pdsd[837];

    auto g_z_y_0_0_z_xy_0_yz = buffer_1100_pdsd[838];

    auto g_z_y_0_0_z_xy_0_zz = buffer_1100_pdsd[839];

    auto g_z_y_0_0_z_xz_0_xx = buffer_1100_pdsd[840];

    auto g_z_y_0_0_z_xz_0_xy = buffer_1100_pdsd[841];

    auto g_z_y_0_0_z_xz_0_xz = buffer_1100_pdsd[842];

    auto g_z_y_0_0_z_xz_0_yy = buffer_1100_pdsd[843];

    auto g_z_y_0_0_z_xz_0_yz = buffer_1100_pdsd[844];

    auto g_z_y_0_0_z_xz_0_zz = buffer_1100_pdsd[845];

    auto g_z_y_0_0_z_yy_0_xx = buffer_1100_pdsd[846];

    auto g_z_y_0_0_z_yy_0_xy = buffer_1100_pdsd[847];

    auto g_z_y_0_0_z_yy_0_xz = buffer_1100_pdsd[848];

    auto g_z_y_0_0_z_yy_0_yy = buffer_1100_pdsd[849];

    auto g_z_y_0_0_z_yy_0_yz = buffer_1100_pdsd[850];

    auto g_z_y_0_0_z_yy_0_zz = buffer_1100_pdsd[851];

    auto g_z_y_0_0_z_yz_0_xx = buffer_1100_pdsd[852];

    auto g_z_y_0_0_z_yz_0_xy = buffer_1100_pdsd[853];

    auto g_z_y_0_0_z_yz_0_xz = buffer_1100_pdsd[854];

    auto g_z_y_0_0_z_yz_0_yy = buffer_1100_pdsd[855];

    auto g_z_y_0_0_z_yz_0_yz = buffer_1100_pdsd[856];

    auto g_z_y_0_0_z_yz_0_zz = buffer_1100_pdsd[857];

    auto g_z_y_0_0_z_zz_0_xx = buffer_1100_pdsd[858];

    auto g_z_y_0_0_z_zz_0_xy = buffer_1100_pdsd[859];

    auto g_z_y_0_0_z_zz_0_xz = buffer_1100_pdsd[860];

    auto g_z_y_0_0_z_zz_0_yy = buffer_1100_pdsd[861];

    auto g_z_y_0_0_z_zz_0_yz = buffer_1100_pdsd[862];

    auto g_z_y_0_0_z_zz_0_zz = buffer_1100_pdsd[863];

    auto g_z_z_0_0_x_xx_0_xx = buffer_1100_pdsd[864];

    auto g_z_z_0_0_x_xx_0_xy = buffer_1100_pdsd[865];

    auto g_z_z_0_0_x_xx_0_xz = buffer_1100_pdsd[866];

    auto g_z_z_0_0_x_xx_0_yy = buffer_1100_pdsd[867];

    auto g_z_z_0_0_x_xx_0_yz = buffer_1100_pdsd[868];

    auto g_z_z_0_0_x_xx_0_zz = buffer_1100_pdsd[869];

    auto g_z_z_0_0_x_xy_0_xx = buffer_1100_pdsd[870];

    auto g_z_z_0_0_x_xy_0_xy = buffer_1100_pdsd[871];

    auto g_z_z_0_0_x_xy_0_xz = buffer_1100_pdsd[872];

    auto g_z_z_0_0_x_xy_0_yy = buffer_1100_pdsd[873];

    auto g_z_z_0_0_x_xy_0_yz = buffer_1100_pdsd[874];

    auto g_z_z_0_0_x_xy_0_zz = buffer_1100_pdsd[875];

    auto g_z_z_0_0_x_xz_0_xx = buffer_1100_pdsd[876];

    auto g_z_z_0_0_x_xz_0_xy = buffer_1100_pdsd[877];

    auto g_z_z_0_0_x_xz_0_xz = buffer_1100_pdsd[878];

    auto g_z_z_0_0_x_xz_0_yy = buffer_1100_pdsd[879];

    auto g_z_z_0_0_x_xz_0_yz = buffer_1100_pdsd[880];

    auto g_z_z_0_0_x_xz_0_zz = buffer_1100_pdsd[881];

    auto g_z_z_0_0_x_yy_0_xx = buffer_1100_pdsd[882];

    auto g_z_z_0_0_x_yy_0_xy = buffer_1100_pdsd[883];

    auto g_z_z_0_0_x_yy_0_xz = buffer_1100_pdsd[884];

    auto g_z_z_0_0_x_yy_0_yy = buffer_1100_pdsd[885];

    auto g_z_z_0_0_x_yy_0_yz = buffer_1100_pdsd[886];

    auto g_z_z_0_0_x_yy_0_zz = buffer_1100_pdsd[887];

    auto g_z_z_0_0_x_yz_0_xx = buffer_1100_pdsd[888];

    auto g_z_z_0_0_x_yz_0_xy = buffer_1100_pdsd[889];

    auto g_z_z_0_0_x_yz_0_xz = buffer_1100_pdsd[890];

    auto g_z_z_0_0_x_yz_0_yy = buffer_1100_pdsd[891];

    auto g_z_z_0_0_x_yz_0_yz = buffer_1100_pdsd[892];

    auto g_z_z_0_0_x_yz_0_zz = buffer_1100_pdsd[893];

    auto g_z_z_0_0_x_zz_0_xx = buffer_1100_pdsd[894];

    auto g_z_z_0_0_x_zz_0_xy = buffer_1100_pdsd[895];

    auto g_z_z_0_0_x_zz_0_xz = buffer_1100_pdsd[896];

    auto g_z_z_0_0_x_zz_0_yy = buffer_1100_pdsd[897];

    auto g_z_z_0_0_x_zz_0_yz = buffer_1100_pdsd[898];

    auto g_z_z_0_0_x_zz_0_zz = buffer_1100_pdsd[899];

    auto g_z_z_0_0_y_xx_0_xx = buffer_1100_pdsd[900];

    auto g_z_z_0_0_y_xx_0_xy = buffer_1100_pdsd[901];

    auto g_z_z_0_0_y_xx_0_xz = buffer_1100_pdsd[902];

    auto g_z_z_0_0_y_xx_0_yy = buffer_1100_pdsd[903];

    auto g_z_z_0_0_y_xx_0_yz = buffer_1100_pdsd[904];

    auto g_z_z_0_0_y_xx_0_zz = buffer_1100_pdsd[905];

    auto g_z_z_0_0_y_xy_0_xx = buffer_1100_pdsd[906];

    auto g_z_z_0_0_y_xy_0_xy = buffer_1100_pdsd[907];

    auto g_z_z_0_0_y_xy_0_xz = buffer_1100_pdsd[908];

    auto g_z_z_0_0_y_xy_0_yy = buffer_1100_pdsd[909];

    auto g_z_z_0_0_y_xy_0_yz = buffer_1100_pdsd[910];

    auto g_z_z_0_0_y_xy_0_zz = buffer_1100_pdsd[911];

    auto g_z_z_0_0_y_xz_0_xx = buffer_1100_pdsd[912];

    auto g_z_z_0_0_y_xz_0_xy = buffer_1100_pdsd[913];

    auto g_z_z_0_0_y_xz_0_xz = buffer_1100_pdsd[914];

    auto g_z_z_0_0_y_xz_0_yy = buffer_1100_pdsd[915];

    auto g_z_z_0_0_y_xz_0_yz = buffer_1100_pdsd[916];

    auto g_z_z_0_0_y_xz_0_zz = buffer_1100_pdsd[917];

    auto g_z_z_0_0_y_yy_0_xx = buffer_1100_pdsd[918];

    auto g_z_z_0_0_y_yy_0_xy = buffer_1100_pdsd[919];

    auto g_z_z_0_0_y_yy_0_xz = buffer_1100_pdsd[920];

    auto g_z_z_0_0_y_yy_0_yy = buffer_1100_pdsd[921];

    auto g_z_z_0_0_y_yy_0_yz = buffer_1100_pdsd[922];

    auto g_z_z_0_0_y_yy_0_zz = buffer_1100_pdsd[923];

    auto g_z_z_0_0_y_yz_0_xx = buffer_1100_pdsd[924];

    auto g_z_z_0_0_y_yz_0_xy = buffer_1100_pdsd[925];

    auto g_z_z_0_0_y_yz_0_xz = buffer_1100_pdsd[926];

    auto g_z_z_0_0_y_yz_0_yy = buffer_1100_pdsd[927];

    auto g_z_z_0_0_y_yz_0_yz = buffer_1100_pdsd[928];

    auto g_z_z_0_0_y_yz_0_zz = buffer_1100_pdsd[929];

    auto g_z_z_0_0_y_zz_0_xx = buffer_1100_pdsd[930];

    auto g_z_z_0_0_y_zz_0_xy = buffer_1100_pdsd[931];

    auto g_z_z_0_0_y_zz_0_xz = buffer_1100_pdsd[932];

    auto g_z_z_0_0_y_zz_0_yy = buffer_1100_pdsd[933];

    auto g_z_z_0_0_y_zz_0_yz = buffer_1100_pdsd[934];

    auto g_z_z_0_0_y_zz_0_zz = buffer_1100_pdsd[935];

    auto g_z_z_0_0_z_xx_0_xx = buffer_1100_pdsd[936];

    auto g_z_z_0_0_z_xx_0_xy = buffer_1100_pdsd[937];

    auto g_z_z_0_0_z_xx_0_xz = buffer_1100_pdsd[938];

    auto g_z_z_0_0_z_xx_0_yy = buffer_1100_pdsd[939];

    auto g_z_z_0_0_z_xx_0_yz = buffer_1100_pdsd[940];

    auto g_z_z_0_0_z_xx_0_zz = buffer_1100_pdsd[941];

    auto g_z_z_0_0_z_xy_0_xx = buffer_1100_pdsd[942];

    auto g_z_z_0_0_z_xy_0_xy = buffer_1100_pdsd[943];

    auto g_z_z_0_0_z_xy_0_xz = buffer_1100_pdsd[944];

    auto g_z_z_0_0_z_xy_0_yy = buffer_1100_pdsd[945];

    auto g_z_z_0_0_z_xy_0_yz = buffer_1100_pdsd[946];

    auto g_z_z_0_0_z_xy_0_zz = buffer_1100_pdsd[947];

    auto g_z_z_0_0_z_xz_0_xx = buffer_1100_pdsd[948];

    auto g_z_z_0_0_z_xz_0_xy = buffer_1100_pdsd[949];

    auto g_z_z_0_0_z_xz_0_xz = buffer_1100_pdsd[950];

    auto g_z_z_0_0_z_xz_0_yy = buffer_1100_pdsd[951];

    auto g_z_z_0_0_z_xz_0_yz = buffer_1100_pdsd[952];

    auto g_z_z_0_0_z_xz_0_zz = buffer_1100_pdsd[953];

    auto g_z_z_0_0_z_yy_0_xx = buffer_1100_pdsd[954];

    auto g_z_z_0_0_z_yy_0_xy = buffer_1100_pdsd[955];

    auto g_z_z_0_0_z_yy_0_xz = buffer_1100_pdsd[956];

    auto g_z_z_0_0_z_yy_0_yy = buffer_1100_pdsd[957];

    auto g_z_z_0_0_z_yy_0_yz = buffer_1100_pdsd[958];

    auto g_z_z_0_0_z_yy_0_zz = buffer_1100_pdsd[959];

    auto g_z_z_0_0_z_yz_0_xx = buffer_1100_pdsd[960];

    auto g_z_z_0_0_z_yz_0_xy = buffer_1100_pdsd[961];

    auto g_z_z_0_0_z_yz_0_xz = buffer_1100_pdsd[962];

    auto g_z_z_0_0_z_yz_0_yy = buffer_1100_pdsd[963];

    auto g_z_z_0_0_z_yz_0_yz = buffer_1100_pdsd[964];

    auto g_z_z_0_0_z_yz_0_zz = buffer_1100_pdsd[965];

    auto g_z_z_0_0_z_zz_0_xx = buffer_1100_pdsd[966];

    auto g_z_z_0_0_z_zz_0_xy = buffer_1100_pdsd[967];

    auto g_z_z_0_0_z_zz_0_xz = buffer_1100_pdsd[968];

    auto g_z_z_0_0_z_zz_0_yy = buffer_1100_pdsd[969];

    auto g_z_z_0_0_z_zz_0_yz = buffer_1100_pdsd[970];

    auto g_z_z_0_0_z_zz_0_zz = buffer_1100_pdsd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xxx_0_xx, g_0_xxx_0_xy, g_0_xxx_0_xz, g_0_xxx_0_yy, g_0_xxx_0_yz, g_0_xxx_0_zz, g_x_x_0_0_x_xx_0_xx, g_x_x_0_0_x_xx_0_xy, g_x_x_0_0_x_xx_0_xz, g_x_x_0_0_x_xx_0_yy, g_x_x_0_0_x_xx_0_yz, g_x_x_0_0_x_xx_0_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_xxx_0_xx, g_xx_xxx_0_xy, g_xx_xxx_0_xz, g_xx_xxx_0_yy, g_xx_xxx_0_yz, g_xx_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_0_xx[i] = 2.0 * g_0_x_0_xx[i] - 2.0 * g_0_xxx_0_xx[i] * b_exp - 4.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_xy[i] = 2.0 * g_0_x_0_xy[i] - 2.0 * g_0_xxx_0_xy[i] * b_exp - 4.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_xz[i] = 2.0 * g_0_x_0_xz[i] - 2.0 * g_0_xxx_0_xz[i] * b_exp - 4.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_yy[i] = 2.0 * g_0_x_0_yy[i] - 2.0 * g_0_xxx_0_yy[i] * b_exp - 4.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_yz[i] = 2.0 * g_0_x_0_yz[i] - 2.0 * g_0_xxx_0_yz[i] * b_exp - 4.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_zz[i] = 2.0 * g_0_x_0_zz[i] - 2.0 * g_0_xxx_0_zz[i] * b_exp - 4.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_x_x_0_0_x_xy_0_xx, g_x_x_0_0_x_xy_0_xy, g_x_x_0_0_x_xy_0_xz, g_x_x_0_0_x_xy_0_yy, g_x_x_0_0_x_xy_0_yz, g_x_x_0_0_x_xy_0_zz, g_xx_xxy_0_xx, g_xx_xxy_0_xy, g_xx_xxy_0_xz, g_xx_xxy_0_yy, g_xx_xxy_0_yz, g_xx_xxy_0_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xy_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_xxy_0_xx[i] * b_exp - 2.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_xxy_0_xy[i] * b_exp - 2.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_xxy_0_xz[i] * b_exp - 2.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_xxy_0_yy[i] * b_exp - 2.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_xxy_0_yz[i] * b_exp - 2.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_xxy_0_zz[i] * b_exp - 2.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_x_x_0_0_x_xz_0_xx, g_x_x_0_0_x_xz_0_xy, g_x_x_0_0_x_xz_0_xz, g_x_x_0_0_x_xz_0_yy, g_x_x_0_0_x_xz_0_yz, g_x_x_0_0_x_xz_0_zz, g_xx_xxz_0_xx, g_xx_xxz_0_xy, g_xx_xxz_0_xz, g_xx_xxz_0_yy, g_xx_xxz_0_yz, g_xx_xxz_0_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_xxz_0_xx[i] * b_exp - 2.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_xxz_0_xy[i] * b_exp - 2.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_xxz_0_xz[i] * b_exp - 2.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_xxz_0_yy[i] * b_exp - 2.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_xxz_0_yz[i] * b_exp - 2.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_xxz_0_zz[i] * b_exp - 2.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_x_x_0_0_x_yy_0_xx, g_x_x_0_0_x_yy_0_xy, g_x_x_0_0_x_yy_0_xz, g_x_x_0_0_x_yy_0_yy, g_x_x_0_0_x_yy_0_yz, g_x_x_0_0_x_yy_0_zz, g_xx_xyy_0_xx, g_xx_xyy_0_xy, g_xx_xyy_0_xz, g_xx_xyy_0_yy, g_xx_xyy_0_yz, g_xx_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yy_0_xx[i] = -2.0 * g_0_xyy_0_xx[i] * b_exp + 4.0 * g_xx_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_xy[i] = -2.0 * g_0_xyy_0_xy[i] * b_exp + 4.0 * g_xx_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_xz[i] = -2.0 * g_0_xyy_0_xz[i] * b_exp + 4.0 * g_xx_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_yy[i] = -2.0 * g_0_xyy_0_yy[i] * b_exp + 4.0 * g_xx_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_yz[i] = -2.0 * g_0_xyy_0_yz[i] * b_exp + 4.0 * g_xx_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_zz[i] = -2.0 * g_0_xyy_0_zz[i] * b_exp + 4.0 * g_xx_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_x_x_0_0_x_yz_0_xx, g_x_x_0_0_x_yz_0_xy, g_x_x_0_0_x_yz_0_xz, g_x_x_0_0_x_yz_0_yy, g_x_x_0_0_x_yz_0_yz, g_x_x_0_0_x_yz_0_zz, g_xx_xyz_0_xx, g_xx_xyz_0_xy, g_xx_xyz_0_xz, g_xx_xyz_0_yy, g_xx_xyz_0_yz, g_xx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_xx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_xx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_xx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_xx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_xx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_xx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_x_x_0_0_x_zz_0_xx, g_x_x_0_0_x_zz_0_xy, g_x_x_0_0_x_zz_0_xz, g_x_x_0_0_x_zz_0_yy, g_x_x_0_0_x_zz_0_yz, g_x_x_0_0_x_zz_0_zz, g_xx_xzz_0_xx, g_xx_xzz_0_xy, g_xx_xzz_0_xz, g_xx_xzz_0_yy, g_xx_xzz_0_yz, g_xx_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_zz_0_xx[i] = -2.0 * g_0_xzz_0_xx[i] * b_exp + 4.0 * g_xx_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_xy[i] = -2.0 * g_0_xzz_0_xy[i] * b_exp + 4.0 * g_xx_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_xz[i] = -2.0 * g_0_xzz_0_xz[i] * b_exp + 4.0 * g_xx_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_yy[i] = -2.0 * g_0_xzz_0_yy[i] * b_exp + 4.0 * g_xx_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_yz[i] = -2.0 * g_0_xzz_0_yz[i] * b_exp + 4.0 * g_xx_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_zz[i] = -2.0 * g_0_xzz_0_zz[i] * b_exp + 4.0 * g_xx_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_0_xx, g_x_x_0_0_y_xx_0_xy, g_x_x_0_0_y_xx_0_xz, g_x_x_0_0_y_xx_0_yy, g_x_x_0_0_y_xx_0_yz, g_x_x_0_0_y_xx_0_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xxx_0_xx, g_xy_xxx_0_xy, g_xy_xxx_0_xz, g_xy_xxx_0_yy, g_xy_xxx_0_yz, g_xy_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_0_xx[i] = -4.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_xy[i] = -4.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_xz[i] = -4.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_yy[i] = -4.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_yz[i] = -4.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_zz[i] = -4.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_y_xy_0_xx, g_x_x_0_0_y_xy_0_xy, g_x_x_0_0_y_xy_0_xz, g_x_x_0_0_y_xy_0_yy, g_x_x_0_0_y_xy_0_yz, g_x_x_0_0_y_xy_0_zz, g_xy_xxy_0_xx, g_xy_xxy_0_xy, g_xy_xxy_0_xz, g_xy_xxy_0_yy, g_xy_xxy_0_yz, g_xy_xxy_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xy_0_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_y_xz_0_xx, g_x_x_0_0_y_xz_0_xy, g_x_x_0_0_y_xz_0_xz, g_x_x_0_0_y_xz_0_yy, g_x_x_0_0_y_xz_0_yz, g_x_x_0_0_y_xz_0_zz, g_xy_xxz_0_xx, g_xy_xxz_0_xy, g_xy_xxz_0_xz, g_xy_xxz_0_yy, g_xy_xxz_0_yz, g_xy_xxz_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xz_0_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_y_yy_0_xx, g_x_x_0_0_y_yy_0_xy, g_x_x_0_0_y_yy_0_xz, g_x_x_0_0_y_yy_0_yy, g_x_x_0_0_y_yy_0_yz, g_x_x_0_0_y_yy_0_zz, g_xy_xyy_0_xx, g_xy_xyy_0_xy, g_xy_xyy_0_xz, g_xy_xyy_0_yy, g_xy_xyy_0_yz, g_xy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yy_0_xx[i] = 4.0 * g_xy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_xy[i] = 4.0 * g_xy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_xz[i] = 4.0 * g_xy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_yy[i] = 4.0 * g_xy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_yz[i] = 4.0 * g_xy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_zz[i] = 4.0 * g_xy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_y_yz_0_xx, g_x_x_0_0_y_yz_0_xy, g_x_x_0_0_y_yz_0_xz, g_x_x_0_0_y_yz_0_yy, g_x_x_0_0_y_yz_0_yz, g_x_x_0_0_y_yz_0_zz, g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yz_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_y_zz_0_xx, g_x_x_0_0_y_zz_0_xy, g_x_x_0_0_y_zz_0_xz, g_x_x_0_0_y_zz_0_yy, g_x_x_0_0_y_zz_0_yz, g_x_x_0_0_y_zz_0_zz, g_xy_xzz_0_xx, g_xy_xzz_0_xy, g_xy_xzz_0_xz, g_xy_xzz_0_yy, g_xy_xzz_0_yz, g_xy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_zz_0_xx[i] = 4.0 * g_xy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_xy[i] = 4.0 * g_xy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_xz[i] = 4.0 * g_xy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_yy[i] = 4.0 * g_xy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_yz[i] = 4.0 * g_xy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_zz[i] = 4.0 * g_xy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_0_xx, g_x_x_0_0_z_xx_0_xy, g_x_x_0_0_z_xx_0_xz, g_x_x_0_0_z_xx_0_yy, g_x_x_0_0_z_xx_0_yz, g_x_x_0_0_z_xx_0_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xxx_0_xx, g_xz_xxx_0_xy, g_xz_xxx_0_xz, g_xz_xxx_0_yy, g_xz_xxx_0_yz, g_xz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_0_xx[i] = -4.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xxx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_xy[i] = -4.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xxx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_xz[i] = -4.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xxx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_yy[i] = -4.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xxx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_yz[i] = -4.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xxx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_zz[i] = -4.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_z_xy_0_xx, g_x_x_0_0_z_xy_0_xy, g_x_x_0_0_z_xy_0_xz, g_x_x_0_0_z_xy_0_yy, g_x_x_0_0_z_xy_0_yz, g_x_x_0_0_z_xy_0_zz, g_xz_xxy_0_xx, g_xz_xxy_0_xy, g_xz_xxy_0_xz, g_xz_xxy_0_yy, g_xz_xxy_0_yz, g_xz_xxy_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xy_0_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_z_xz_0_xx, g_x_x_0_0_z_xz_0_xy, g_x_x_0_0_z_xz_0_xz, g_x_x_0_0_z_xz_0_yy, g_x_x_0_0_z_xz_0_yz, g_x_x_0_0_z_xz_0_zz, g_xz_xxz_0_xx, g_xz_xxz_0_xy, g_xz_xxz_0_xz, g_xz_xxz_0_yy, g_xz_xxz_0_yz, g_xz_xxz_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xz_0_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_z_yy_0_xx, g_x_x_0_0_z_yy_0_xy, g_x_x_0_0_z_yy_0_xz, g_x_x_0_0_z_yy_0_yy, g_x_x_0_0_z_yy_0_yz, g_x_x_0_0_z_yy_0_zz, g_xz_xyy_0_xx, g_xz_xyy_0_xy, g_xz_xyy_0_xz, g_xz_xyy_0_yy, g_xz_xyy_0_yz, g_xz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yy_0_xx[i] = 4.0 * g_xz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_xy[i] = 4.0 * g_xz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_xz[i] = 4.0 * g_xz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_yy[i] = 4.0 * g_xz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_yz[i] = 4.0 * g_xz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_zz[i] = 4.0 * g_xz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_z_yz_0_xx, g_x_x_0_0_z_yz_0_xy, g_x_x_0_0_z_yz_0_xz, g_x_x_0_0_z_yz_0_yy, g_x_x_0_0_z_yz_0_yz, g_x_x_0_0_z_yz_0_zz, g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yz_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_z_zz_0_xx, g_x_x_0_0_z_zz_0_xy, g_x_x_0_0_z_zz_0_xz, g_x_x_0_0_z_zz_0_yy, g_x_x_0_0_z_zz_0_yz, g_x_x_0_0_z_zz_0_zz, g_xz_xzz_0_xx, g_xz_xzz_0_xy, g_xz_xzz_0_xz, g_xz_xzz_0_yy, g_xz_xzz_0_yz, g_xz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_zz_0_xx[i] = 4.0 * g_xz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_xy[i] = 4.0 * g_xz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_xz[i] = 4.0 * g_xz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_yy[i] = 4.0 * g_xz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_yz[i] = 4.0 * g_xz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_zz[i] = 4.0 * g_xz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_x_y_0_0_x_xx_0_xx, g_x_y_0_0_x_xx_0_xy, g_x_y_0_0_x_xx_0_xz, g_x_y_0_0_x_xx_0_yy, g_x_y_0_0_x_xx_0_yz, g_x_y_0_0_x_xx_0_zz, g_xx_xxy_0_xx, g_xx_xxy_0_xy, g_xx_xxy_0_xz, g_xx_xxy_0_yy, g_xx_xxy_0_yz, g_xx_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_0_xx[i] = -2.0 * g_0_xxy_0_xx[i] * b_exp + 4.0 * g_xx_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_xy[i] = -2.0 * g_0_xxy_0_xy[i] * b_exp + 4.0 * g_xx_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_xz[i] = -2.0 * g_0_xxy_0_xz[i] * b_exp + 4.0 * g_xx_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_yy[i] = -2.0 * g_0_xxy_0_yy[i] * b_exp + 4.0 * g_xx_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_yz[i] = -2.0 * g_0_xxy_0_yz[i] * b_exp + 4.0 * g_xx_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_zz[i] = -2.0 * g_0_xxy_0_zz[i] * b_exp + 4.0 * g_xx_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_x_y_0_0_x_xy_0_xx, g_x_y_0_0_x_xy_0_xy, g_x_y_0_0_x_xy_0_xz, g_x_y_0_0_x_xy_0_yy, g_x_y_0_0_x_xy_0_yz, g_x_y_0_0_x_xy_0_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_xyy_0_xx, g_xx_xyy_0_xy, g_xx_xyy_0_xz, g_xx_xyy_0_yy, g_xx_xyy_0_yz, g_xx_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xy_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xyy_0_xx[i] * b_exp - 2.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xyy_0_xy[i] * b_exp - 2.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xyy_0_xz[i] * b_exp - 2.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xyy_0_yy[i] * b_exp - 2.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xyy_0_yz[i] * b_exp - 2.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xyy_0_zz[i] * b_exp - 2.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_x_y_0_0_x_xz_0_xx, g_x_y_0_0_x_xz_0_xy, g_x_y_0_0_x_xz_0_xz, g_x_y_0_0_x_xz_0_yy, g_x_y_0_0_x_xz_0_yz, g_x_y_0_0_x_xz_0_zz, g_xx_xyz_0_xx, g_xx_xyz_0_xy, g_xx_xyz_0_xz, g_xx_xyz_0_yy, g_xx_xyz_0_yz, g_xx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_xx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_xx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_xx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_xx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_xx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_xx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yyy_0_xx, g_0_yyy_0_xy, g_0_yyy_0_xz, g_0_yyy_0_yy, g_0_yyy_0_yz, g_0_yyy_0_zz, g_x_y_0_0_x_yy_0_xx, g_x_y_0_0_x_yy_0_xy, g_x_y_0_0_x_yy_0_xz, g_x_y_0_0_x_yy_0_yy, g_x_y_0_0_x_yy_0_yz, g_x_y_0_0_x_yy_0_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz, g_xx_yyy_0_xx, g_xx_yyy_0_xy, g_xx_yyy_0_xz, g_xx_yyy_0_yy, g_xx_yyy_0_yz, g_xx_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yy_0_xx[i] = 2.0 * g_0_y_0_xx[i] - 2.0 * g_0_yyy_0_xx[i] * b_exp - 4.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_xy[i] = 2.0 * g_0_y_0_xy[i] - 2.0 * g_0_yyy_0_xy[i] * b_exp - 4.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_xz[i] = 2.0 * g_0_y_0_xz[i] - 2.0 * g_0_yyy_0_xz[i] * b_exp - 4.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_yy[i] = 2.0 * g_0_y_0_yy[i] - 2.0 * g_0_yyy_0_yy[i] * b_exp - 4.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_yz[i] = 2.0 * g_0_y_0_yz[i] - 2.0 * g_0_yyy_0_yz[i] * b_exp - 4.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_zz[i] = 2.0 * g_0_y_0_zz[i] - 2.0 * g_0_yyy_0_zz[i] * b_exp - 4.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_x_y_0_0_x_yz_0_xx, g_x_y_0_0_x_yz_0_xy, g_x_y_0_0_x_yz_0_xz, g_x_y_0_0_x_yz_0_yy, g_x_y_0_0_x_yz_0_yz, g_x_y_0_0_x_yz_0_zz, g_xx_yyz_0_xx, g_xx_yyz_0_xy, g_xx_yyz_0_xz, g_xx_yyz_0_yy, g_xx_yyz_0_yz, g_xx_yyz_0_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_yyz_0_xx[i] * b_exp - 2.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_yyz_0_xy[i] * b_exp - 2.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_yyz_0_xz[i] * b_exp - 2.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_yyz_0_yy[i] * b_exp - 2.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_yyz_0_yz[i] * b_exp - 2.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_yyz_0_zz[i] * b_exp - 2.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_x_y_0_0_x_zz_0_xx, g_x_y_0_0_x_zz_0_xy, g_x_y_0_0_x_zz_0_xz, g_x_y_0_0_x_zz_0_yy, g_x_y_0_0_x_zz_0_yz, g_x_y_0_0_x_zz_0_zz, g_xx_yzz_0_xx, g_xx_yzz_0_xy, g_xx_yzz_0_xz, g_xx_yzz_0_yy, g_xx_yzz_0_yz, g_xx_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_zz_0_xx[i] = -2.0 * g_0_yzz_0_xx[i] * b_exp + 4.0 * g_xx_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_xy[i] = -2.0 * g_0_yzz_0_xy[i] * b_exp + 4.0 * g_xx_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_xz[i] = -2.0 * g_0_yzz_0_xz[i] * b_exp + 4.0 * g_xx_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_yy[i] = -2.0 * g_0_yzz_0_yy[i] * b_exp + 4.0 * g_xx_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_yz[i] = -2.0 * g_0_yzz_0_yz[i] * b_exp + 4.0 * g_xx_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_zz[i] = -2.0 * g_0_yzz_0_zz[i] * b_exp + 4.0 * g_xx_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_0_xx, g_x_y_0_0_y_xx_0_xy, g_x_y_0_0_y_xx_0_xz, g_x_y_0_0_y_xx_0_yy, g_x_y_0_0_y_xx_0_yz, g_x_y_0_0_y_xx_0_zz, g_xy_xxy_0_xx, g_xy_xxy_0_xy, g_xy_xxy_0_xz, g_xy_xxy_0_yy, g_xy_xxy_0_yz, g_xy_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_0_xx[i] = 4.0 * g_xy_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_xy[i] = 4.0 * g_xy_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_xz[i] = 4.0 * g_xy_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_yy[i] = 4.0 * g_xy_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_yz[i] = 4.0 * g_xy_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_zz[i] = 4.0 * g_xy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_y_0_0_y_xy_0_xx, g_x_y_0_0_y_xy_0_xy, g_x_y_0_0_y_xy_0_xz, g_x_y_0_0_y_xy_0_yy, g_x_y_0_0_y_xy_0_yz, g_x_y_0_0_y_xy_0_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xyy_0_xx, g_xy_xyy_0_xy, g_xy_xyy_0_xz, g_xy_xyy_0_yy, g_xy_xyy_0_yz, g_xy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xy_0_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_y_0_0_y_xz_0_xx, g_x_y_0_0_y_xz_0_xy, g_x_y_0_0_y_xz_0_xz, g_x_y_0_0_y_xz_0_yy, g_x_y_0_0_y_xz_0_yz, g_x_y_0_0_y_xz_0_zz, g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xz_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_y_0_0_y_yy_0_xx, g_x_y_0_0_y_yy_0_xy, g_x_y_0_0_y_yy_0_xz, g_x_y_0_0_y_yy_0_yy, g_x_y_0_0_y_yy_0_yz, g_x_y_0_0_y_yy_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_yyy_0_xx, g_xy_yyy_0_xy, g_xy_yyy_0_xz, g_xy_yyy_0_yy, g_xy_yyy_0_yz, g_xy_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yy_0_xx[i] = -4.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_xy[i] = -4.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_xz[i] = -4.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_yy[i] = -4.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_yz[i] = -4.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_zz[i] = -4.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_y_0_0_y_yz_0_xx, g_x_y_0_0_y_yz_0_xy, g_x_y_0_0_y_yz_0_xz, g_x_y_0_0_y_yz_0_yy, g_x_y_0_0_y_yz_0_yz, g_x_y_0_0_y_yz_0_zz, g_xy_yyz_0_xx, g_xy_yyz_0_xy, g_xy_yyz_0_xz, g_xy_yyz_0_yy, g_xy_yyz_0_yz, g_xy_yyz_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yz_0_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_y_0_0_y_zz_0_xx, g_x_y_0_0_y_zz_0_xy, g_x_y_0_0_y_zz_0_xz, g_x_y_0_0_y_zz_0_yy, g_x_y_0_0_y_zz_0_yz, g_x_y_0_0_y_zz_0_zz, g_xy_yzz_0_xx, g_xy_yzz_0_xy, g_xy_yzz_0_xz, g_xy_yzz_0_yy, g_xy_yzz_0_yz, g_xy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_zz_0_xx[i] = 4.0 * g_xy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_xy[i] = 4.0 * g_xy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_xz[i] = 4.0 * g_xy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_yy[i] = 4.0 * g_xy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_yz[i] = 4.0 * g_xy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_zz[i] = 4.0 * g_xy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_0_xx, g_x_y_0_0_z_xx_0_xy, g_x_y_0_0_z_xx_0_xz, g_x_y_0_0_z_xx_0_yy, g_x_y_0_0_z_xx_0_yz, g_x_y_0_0_z_xx_0_zz, g_xz_xxy_0_xx, g_xz_xxy_0_xy, g_xz_xxy_0_xz, g_xz_xxy_0_yy, g_xz_xxy_0_yz, g_xz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_0_xx[i] = 4.0 * g_xz_xxy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_xy[i] = 4.0 * g_xz_xxy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_xz[i] = 4.0 * g_xz_xxy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_yy[i] = 4.0 * g_xz_xxy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_yz[i] = 4.0 * g_xz_xxy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_zz[i] = 4.0 * g_xz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_y_0_0_z_xy_0_xx, g_x_y_0_0_z_xy_0_xy, g_x_y_0_0_z_xy_0_xz, g_x_y_0_0_z_xy_0_yy, g_x_y_0_0_z_xy_0_yz, g_x_y_0_0_z_xy_0_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xyy_0_xx, g_xz_xyy_0_xy, g_xz_xyy_0_xz, g_xz_xyy_0_yy, g_xz_xyy_0_yz, g_xz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xy_0_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_y_0_0_z_xz_0_xx, g_x_y_0_0_z_xz_0_xy, g_x_y_0_0_z_xz_0_xz, g_x_y_0_0_z_xz_0_yy, g_x_y_0_0_z_xz_0_yz, g_x_y_0_0_z_xz_0_zz, g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xz_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_y_0_0_z_yy_0_xx, g_x_y_0_0_z_yy_0_xy, g_x_y_0_0_z_yy_0_xz, g_x_y_0_0_z_yy_0_yy, g_x_y_0_0_z_yy_0_yz, g_x_y_0_0_z_yy_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_yyy_0_xx, g_xz_yyy_0_xy, g_xz_yyy_0_xz, g_xz_yyy_0_yy, g_xz_yyy_0_yz, g_xz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yy_0_xx[i] = -4.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_yyy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_xy[i] = -4.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_yyy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_xz[i] = -4.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_yyy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_yy[i] = -4.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_yyy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_yz[i] = -4.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_yyy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_zz[i] = -4.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_y_0_0_z_yz_0_xx, g_x_y_0_0_z_yz_0_xy, g_x_y_0_0_z_yz_0_xz, g_x_y_0_0_z_yz_0_yy, g_x_y_0_0_z_yz_0_yz, g_x_y_0_0_z_yz_0_zz, g_xz_yyz_0_xx, g_xz_yyz_0_xy, g_xz_yyz_0_xz, g_xz_yyz_0_yy, g_xz_yyz_0_yz, g_xz_yyz_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yz_0_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_y_0_0_z_zz_0_xx, g_x_y_0_0_z_zz_0_xy, g_x_y_0_0_z_zz_0_xz, g_x_y_0_0_z_zz_0_yy, g_x_y_0_0_z_zz_0_yz, g_x_y_0_0_z_zz_0_zz, g_xz_yzz_0_xx, g_xz_yzz_0_xy, g_xz_yzz_0_xz, g_xz_yzz_0_yy, g_xz_yzz_0_yz, g_xz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_zz_0_xx[i] = 4.0 * g_xz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_xy[i] = 4.0 * g_xz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_xz[i] = 4.0 * g_xz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_yy[i] = 4.0 * g_xz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_yz[i] = 4.0 * g_xz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_zz[i] = 4.0 * g_xz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_x_z_0_0_x_xx_0_xx, g_x_z_0_0_x_xx_0_xy, g_x_z_0_0_x_xx_0_xz, g_x_z_0_0_x_xx_0_yy, g_x_z_0_0_x_xx_0_yz, g_x_z_0_0_x_xx_0_zz, g_xx_xxz_0_xx, g_xx_xxz_0_xy, g_xx_xxz_0_xz, g_xx_xxz_0_yy, g_xx_xxz_0_yz, g_xx_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_0_xx[i] = -2.0 * g_0_xxz_0_xx[i] * b_exp + 4.0 * g_xx_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_xy[i] = -2.0 * g_0_xxz_0_xy[i] * b_exp + 4.0 * g_xx_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_xz[i] = -2.0 * g_0_xxz_0_xz[i] * b_exp + 4.0 * g_xx_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_yy[i] = -2.0 * g_0_xxz_0_yy[i] * b_exp + 4.0 * g_xx_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_yz[i] = -2.0 * g_0_xxz_0_yz[i] * b_exp + 4.0 * g_xx_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_zz[i] = -2.0 * g_0_xxz_0_zz[i] * b_exp + 4.0 * g_xx_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_x_z_0_0_x_xy_0_xx, g_x_z_0_0_x_xy_0_xy, g_x_z_0_0_x_xy_0_xz, g_x_z_0_0_x_xy_0_yy, g_x_z_0_0_x_xy_0_yz, g_x_z_0_0_x_xy_0_zz, g_xx_xyz_0_xx, g_xx_xyz_0_xy, g_xx_xyz_0_xz, g_xx_xyz_0_yy, g_xx_xyz_0_yz, g_xx_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xy_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_xx_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_xx_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_xx_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_xx_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_xx_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_xx_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_x_z_0_0_x_xz_0_xx, g_x_z_0_0_x_xz_0_xy, g_x_z_0_0_x_xz_0_xz, g_x_z_0_0_x_xz_0_yy, g_x_z_0_0_x_xz_0_yz, g_x_z_0_0_x_xz_0_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_xzz_0_xx, g_xx_xzz_0_xy, g_xx_xzz_0_xz, g_xx_xzz_0_yy, g_xx_xzz_0_yz, g_xx_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xz_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xzz_0_xx[i] * b_exp - 2.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xzz_0_xy[i] * b_exp - 2.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xzz_0_xz[i] * b_exp - 2.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xzz_0_yy[i] * b_exp - 2.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xzz_0_yz[i] * b_exp - 2.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xzz_0_zz[i] * b_exp - 2.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_x_z_0_0_x_yy_0_xx, g_x_z_0_0_x_yy_0_xy, g_x_z_0_0_x_yy_0_xz, g_x_z_0_0_x_yy_0_yy, g_x_z_0_0_x_yy_0_yz, g_x_z_0_0_x_yy_0_zz, g_xx_yyz_0_xx, g_xx_yyz_0_xy, g_xx_yyz_0_xz, g_xx_yyz_0_yy, g_xx_yyz_0_yz, g_xx_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yy_0_xx[i] = -2.0 * g_0_yyz_0_xx[i] * b_exp + 4.0 * g_xx_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_xy[i] = -2.0 * g_0_yyz_0_xy[i] * b_exp + 4.0 * g_xx_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_xz[i] = -2.0 * g_0_yyz_0_xz[i] * b_exp + 4.0 * g_xx_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_yy[i] = -2.0 * g_0_yyz_0_yy[i] * b_exp + 4.0 * g_xx_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_yz[i] = -2.0 * g_0_yyz_0_yz[i] * b_exp + 4.0 * g_xx_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_zz[i] = -2.0 * g_0_yyz_0_zz[i] * b_exp + 4.0 * g_xx_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_x_z_0_0_x_yz_0_xx, g_x_z_0_0_x_yz_0_xy, g_x_z_0_0_x_yz_0_xz, g_x_z_0_0_x_yz_0_yy, g_x_z_0_0_x_yz_0_yz, g_x_z_0_0_x_yz_0_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz, g_xx_yzz_0_xx, g_xx_yzz_0_xy, g_xx_yzz_0_xz, g_xx_yzz_0_yy, g_xx_yzz_0_yz, g_xx_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yz_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_yzz_0_xx[i] * b_exp - 2.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_yzz_0_xy[i] * b_exp - 2.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_yzz_0_xz[i] * b_exp - 2.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_yzz_0_yy[i] * b_exp - 2.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_yzz_0_yz[i] * b_exp - 2.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_yzz_0_zz[i] * b_exp - 2.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_zzz_0_xx, g_0_zzz_0_xy, g_0_zzz_0_xz, g_0_zzz_0_yy, g_0_zzz_0_yz, g_0_zzz_0_zz, g_x_z_0_0_x_zz_0_xx, g_x_z_0_0_x_zz_0_xy, g_x_z_0_0_x_zz_0_xz, g_x_z_0_0_x_zz_0_yy, g_x_z_0_0_x_zz_0_yz, g_x_z_0_0_x_zz_0_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz, g_xx_zzz_0_xx, g_xx_zzz_0_xy, g_xx_zzz_0_xz, g_xx_zzz_0_yy, g_xx_zzz_0_yz, g_xx_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_zz_0_xx[i] = 2.0 * g_0_z_0_xx[i] - 2.0 * g_0_zzz_0_xx[i] * b_exp - 4.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_xy[i] = 2.0 * g_0_z_0_xy[i] - 2.0 * g_0_zzz_0_xy[i] * b_exp - 4.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_xz[i] = 2.0 * g_0_z_0_xz[i] - 2.0 * g_0_zzz_0_xz[i] * b_exp - 4.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_yy[i] = 2.0 * g_0_z_0_yy[i] - 2.0 * g_0_zzz_0_yy[i] * b_exp - 4.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_yz[i] = 2.0 * g_0_z_0_yz[i] - 2.0 * g_0_zzz_0_yz[i] * b_exp - 4.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_zz[i] = 2.0 * g_0_z_0_zz[i] - 2.0 * g_0_zzz_0_zz[i] * b_exp - 4.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_0_xx, g_x_z_0_0_y_xx_0_xy, g_x_z_0_0_y_xx_0_xz, g_x_z_0_0_y_xx_0_yy, g_x_z_0_0_y_xx_0_yz, g_x_z_0_0_y_xx_0_zz, g_xy_xxz_0_xx, g_xy_xxz_0_xy, g_xy_xxz_0_xz, g_xy_xxz_0_yy, g_xy_xxz_0_yz, g_xy_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_0_xx[i] = 4.0 * g_xy_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_xy[i] = 4.0 * g_xy_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_xz[i] = 4.0 * g_xy_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_yy[i] = 4.0 * g_xy_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_yz[i] = 4.0 * g_xy_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_zz[i] = 4.0 * g_xy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_z_0_0_y_xy_0_xx, g_x_z_0_0_y_xy_0_xy, g_x_z_0_0_y_xy_0_xz, g_x_z_0_0_y_xy_0_yy, g_x_z_0_0_y_xy_0_yz, g_x_z_0_0_y_xy_0_zz, g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xy_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_z_0_0_y_xz_0_xx, g_x_z_0_0_y_xz_0_xy, g_x_z_0_0_y_xz_0_xz, g_x_z_0_0_y_xz_0_yy, g_x_z_0_0_y_xz_0_yz, g_x_z_0_0_y_xz_0_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xzz_0_xx, g_xy_xzz_0_xy, g_xy_xzz_0_xz, g_xy_xzz_0_yy, g_xy_xzz_0_yz, g_xy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xz_0_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_z_0_0_y_yy_0_xx, g_x_z_0_0_y_yy_0_xy, g_x_z_0_0_y_yy_0_xz, g_x_z_0_0_y_yy_0_yy, g_x_z_0_0_y_yy_0_yz, g_x_z_0_0_y_yy_0_zz, g_xy_yyz_0_xx, g_xy_yyz_0_xy, g_xy_yyz_0_xz, g_xy_yyz_0_yy, g_xy_yyz_0_yz, g_xy_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yy_0_xx[i] = 4.0 * g_xy_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_xy[i] = 4.0 * g_xy_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_xz[i] = 4.0 * g_xy_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_yy[i] = 4.0 * g_xy_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_yz[i] = 4.0 * g_xy_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_zz[i] = 4.0 * g_xy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_z_0_0_y_yz_0_xx, g_x_z_0_0_y_yz_0_xy, g_x_z_0_0_y_yz_0_xz, g_x_z_0_0_y_yz_0_yy, g_x_z_0_0_y_yz_0_yz, g_x_z_0_0_y_yz_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_yzz_0_xx, g_xy_yzz_0_xy, g_xy_yzz_0_xz, g_xy_yzz_0_yy, g_xy_yzz_0_yz, g_xy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yz_0_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_z_0_0_y_zz_0_xx, g_x_z_0_0_y_zz_0_xy, g_x_z_0_0_y_zz_0_xz, g_x_z_0_0_y_zz_0_yy, g_x_z_0_0_y_zz_0_yz, g_x_z_0_0_y_zz_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_zzz_0_xx, g_xy_zzz_0_xy, g_xy_zzz_0_xz, g_xy_zzz_0_yy, g_xy_zzz_0_yz, g_xy_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_zz_0_xx[i] = -4.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_xy[i] = -4.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_xz[i] = -4.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_yy[i] = -4.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_yz[i] = -4.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_zz[i] = -4.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_0_xx, g_x_z_0_0_z_xx_0_xy, g_x_z_0_0_z_xx_0_xz, g_x_z_0_0_z_xx_0_yy, g_x_z_0_0_z_xx_0_yz, g_x_z_0_0_z_xx_0_zz, g_xz_xxz_0_xx, g_xz_xxz_0_xy, g_xz_xxz_0_xz, g_xz_xxz_0_yy, g_xz_xxz_0_yz, g_xz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_0_xx[i] = 4.0 * g_xz_xxz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_xy[i] = 4.0 * g_xz_xxz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_xz[i] = 4.0 * g_xz_xxz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_yy[i] = 4.0 * g_xz_xxz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_yz[i] = 4.0 * g_xz_xxz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_zz[i] = 4.0 * g_xz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_z_0_0_z_xy_0_xx, g_x_z_0_0_z_xy_0_xy, g_x_z_0_0_z_xy_0_xz, g_x_z_0_0_z_xy_0_yy, g_x_z_0_0_z_xy_0_yz, g_x_z_0_0_z_xy_0_zz, g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xy_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_z_0_0_z_xz_0_xx, g_x_z_0_0_z_xz_0_xy, g_x_z_0_0_z_xz_0_xz, g_x_z_0_0_z_xz_0_yy, g_x_z_0_0_z_xz_0_yz, g_x_z_0_0_z_xz_0_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xzz_0_xx, g_xz_xzz_0_xy, g_xz_xzz_0_xz, g_xz_xzz_0_yy, g_xz_xzz_0_yz, g_xz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xz_0_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_z_0_0_z_yy_0_xx, g_x_z_0_0_z_yy_0_xy, g_x_z_0_0_z_yy_0_xz, g_x_z_0_0_z_yy_0_yy, g_x_z_0_0_z_yy_0_yz, g_x_z_0_0_z_yy_0_zz, g_xz_yyz_0_xx, g_xz_yyz_0_xy, g_xz_yyz_0_xz, g_xz_yyz_0_yy, g_xz_yyz_0_yz, g_xz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yy_0_xx[i] = 4.0 * g_xz_yyz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_xy[i] = 4.0 * g_xz_yyz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_xz[i] = 4.0 * g_xz_yyz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_yy[i] = 4.0 * g_xz_yyz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_yz[i] = 4.0 * g_xz_yyz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_zz[i] = 4.0 * g_xz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_z_0_0_z_yz_0_xx, g_x_z_0_0_z_yz_0_xy, g_x_z_0_0_z_yz_0_xz, g_x_z_0_0_z_yz_0_yy, g_x_z_0_0_z_yz_0_yz, g_x_z_0_0_z_yz_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_yzz_0_xx, g_xz_yzz_0_xy, g_xz_yzz_0_xz, g_xz_yzz_0_yy, g_xz_yzz_0_yz, g_xz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yz_0_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_yzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_yzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_yzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_yzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_yzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_z_0_0_z_zz_0_xx, g_x_z_0_0_z_zz_0_xy, g_x_z_0_0_z_zz_0_xz, g_x_z_0_0_z_zz_0_yy, g_x_z_0_0_z_zz_0_yz, g_x_z_0_0_z_zz_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_zzz_0_xx, g_xz_zzz_0_xy, g_xz_zzz_0_xz, g_xz_zzz_0_yy, g_xz_zzz_0_yz, g_xz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_zz_0_xx[i] = -4.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_zzz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_xy[i] = -4.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_zzz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_xz[i] = -4.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_zzz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_yy[i] = -4.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_zzz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_yz[i] = -4.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_zzz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_zz[i] = -4.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xxx_0_xx, g_xy_xxx_0_xy, g_xy_xxx_0_xz, g_xy_xxx_0_yy, g_xy_xxx_0_yz, g_xy_xxx_0_zz, g_y_x_0_0_x_xx_0_xx, g_y_x_0_0_x_xx_0_xy, g_y_x_0_0_x_xx_0_xz, g_y_x_0_0_x_xx_0_yy, g_y_x_0_0_x_xx_0_yz, g_y_x_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_0_xx[i] = -4.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_xy[i] = -4.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_xz[i] = -4.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_yy[i] = -4.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_yz[i] = -4.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_zz[i] = -4.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_xxy_0_xx, g_xy_xxy_0_xy, g_xy_xxy_0_xz, g_xy_xxy_0_yy, g_xy_xxy_0_yz, g_xy_xxy_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_y_x_0_0_x_xy_0_xx, g_y_x_0_0_x_xy_0_xy, g_y_x_0_0_x_xy_0_xz, g_y_x_0_0_x_xy_0_yy, g_y_x_0_0_x_xy_0_yz, g_y_x_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xy_0_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_xxz_0_xx, g_xy_xxz_0_xy, g_xy_xxz_0_xz, g_xy_xxz_0_yy, g_xy_xxz_0_yz, g_xy_xxz_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_y_x_0_0_x_xz_0_xx, g_y_x_0_0_x_xz_0_xy, g_y_x_0_0_x_xz_0_xz, g_y_x_0_0_x_xz_0_yy, g_y_x_0_0_x_xz_0_yz, g_y_x_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xz_0_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_xyy_0_xx, g_xy_xyy_0_xy, g_xy_xyy_0_xz, g_xy_xyy_0_yy, g_xy_xyy_0_yz, g_xy_xyy_0_zz, g_y_x_0_0_x_yy_0_xx, g_y_x_0_0_x_yy_0_xy, g_y_x_0_0_x_yy_0_xz, g_y_x_0_0_x_yy_0_yy, g_y_x_0_0_x_yy_0_yz, g_y_x_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yy_0_xx[i] = 4.0 * g_xy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_xy[i] = 4.0 * g_xy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_xz[i] = 4.0 * g_xy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_yy[i] = 4.0 * g_xy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_yz[i] = 4.0 * g_xy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_zz[i] = 4.0 * g_xy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz, g_y_x_0_0_x_yz_0_xx, g_y_x_0_0_x_yz_0_xy, g_y_x_0_0_x_yz_0_xz, g_y_x_0_0_x_yz_0_yy, g_y_x_0_0_x_yz_0_yz, g_y_x_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yz_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_xzz_0_xx, g_xy_xzz_0_xy, g_xy_xzz_0_xz, g_xy_xzz_0_yy, g_xy_xzz_0_yz, g_xy_xzz_0_zz, g_y_x_0_0_x_zz_0_xx, g_y_x_0_0_x_zz_0_xy, g_y_x_0_0_x_zz_0_xz, g_y_x_0_0_x_zz_0_yy, g_y_x_0_0_x_zz_0_yz, g_y_x_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_zz_0_xx[i] = 4.0 * g_xy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_xy[i] = 4.0 * g_xy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_xz[i] = 4.0 * g_xy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_yy[i] = 4.0 * g_xy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_yz[i] = 4.0 * g_xy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_zz[i] = 4.0 * g_xy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xxx_0_xx, g_0_xxx_0_xy, g_0_xxx_0_xz, g_0_xxx_0_yy, g_0_xxx_0_yz, g_0_xxx_0_zz, g_y_x_0_0_y_xx_0_xx, g_y_x_0_0_y_xx_0_xy, g_y_x_0_0_y_xx_0_xz, g_y_x_0_0_y_xx_0_yy, g_y_x_0_0_y_xx_0_yz, g_y_x_0_0_y_xx_0_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_xxx_0_xx, g_yy_xxx_0_xy, g_yy_xxx_0_xz, g_yy_xxx_0_yy, g_yy_xxx_0_yz, g_yy_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_0_xx[i] = 2.0 * g_0_x_0_xx[i] - 2.0 * g_0_xxx_0_xx[i] * b_exp - 4.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_xy[i] = 2.0 * g_0_x_0_xy[i] - 2.0 * g_0_xxx_0_xy[i] * b_exp - 4.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_xz[i] = 2.0 * g_0_x_0_xz[i] - 2.0 * g_0_xxx_0_xz[i] * b_exp - 4.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_yy[i] = 2.0 * g_0_x_0_yy[i] - 2.0 * g_0_xxx_0_yy[i] * b_exp - 4.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_yz[i] = 2.0 * g_0_x_0_yz[i] - 2.0 * g_0_xxx_0_yz[i] * b_exp - 4.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_zz[i] = 2.0 * g_0_x_0_zz[i] - 2.0 * g_0_xxx_0_zz[i] * b_exp - 4.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_y_x_0_0_y_xy_0_xx, g_y_x_0_0_y_xy_0_xy, g_y_x_0_0_y_xy_0_xz, g_y_x_0_0_y_xy_0_yy, g_y_x_0_0_y_xy_0_yz, g_y_x_0_0_y_xy_0_zz, g_yy_xxy_0_xx, g_yy_xxy_0_xy, g_yy_xxy_0_xz, g_yy_xxy_0_yy, g_yy_xxy_0_yz, g_yy_xxy_0_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xy_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_xxy_0_xx[i] * b_exp - 2.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_xxy_0_xy[i] * b_exp - 2.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_xxy_0_xz[i] * b_exp - 2.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_xxy_0_yy[i] * b_exp - 2.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_xxy_0_yz[i] * b_exp - 2.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_xxy_0_zz[i] * b_exp - 2.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_y_x_0_0_y_xz_0_xx, g_y_x_0_0_y_xz_0_xy, g_y_x_0_0_y_xz_0_xz, g_y_x_0_0_y_xz_0_yy, g_y_x_0_0_y_xz_0_yz, g_y_x_0_0_y_xz_0_zz, g_yy_xxz_0_xx, g_yy_xxz_0_xy, g_yy_xxz_0_xz, g_yy_xxz_0_yy, g_yy_xxz_0_yz, g_yy_xxz_0_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_xxz_0_xx[i] * b_exp - 2.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_xxz_0_xy[i] * b_exp - 2.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_xxz_0_xz[i] * b_exp - 2.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_xxz_0_yy[i] * b_exp - 2.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_xxz_0_yz[i] * b_exp - 2.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_xxz_0_zz[i] * b_exp - 2.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_y_x_0_0_y_yy_0_xx, g_y_x_0_0_y_yy_0_xy, g_y_x_0_0_y_yy_0_xz, g_y_x_0_0_y_yy_0_yy, g_y_x_0_0_y_yy_0_yz, g_y_x_0_0_y_yy_0_zz, g_yy_xyy_0_xx, g_yy_xyy_0_xy, g_yy_xyy_0_xz, g_yy_xyy_0_yy, g_yy_xyy_0_yz, g_yy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yy_0_xx[i] = -2.0 * g_0_xyy_0_xx[i] * b_exp + 4.0 * g_yy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_xy[i] = -2.0 * g_0_xyy_0_xy[i] * b_exp + 4.0 * g_yy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_xz[i] = -2.0 * g_0_xyy_0_xz[i] * b_exp + 4.0 * g_yy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_yy[i] = -2.0 * g_0_xyy_0_yy[i] * b_exp + 4.0 * g_yy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_yz[i] = -2.0 * g_0_xyy_0_yz[i] * b_exp + 4.0 * g_yy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_zz[i] = -2.0 * g_0_xyy_0_zz[i] * b_exp + 4.0 * g_yy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_y_x_0_0_y_yz_0_xx, g_y_x_0_0_y_yz_0_xy, g_y_x_0_0_y_yz_0_xz, g_y_x_0_0_y_yz_0_yy, g_y_x_0_0_y_yz_0_yz, g_y_x_0_0_y_yz_0_zz, g_yy_xyz_0_xx, g_yy_xyz_0_xy, g_yy_xyz_0_xz, g_yy_xyz_0_yy, g_yy_xyz_0_yz, g_yy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_yy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_yy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_yy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_yy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_yy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_yy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_y_x_0_0_y_zz_0_xx, g_y_x_0_0_y_zz_0_xy, g_y_x_0_0_y_zz_0_xz, g_y_x_0_0_y_zz_0_yy, g_y_x_0_0_y_zz_0_yz, g_y_x_0_0_y_zz_0_zz, g_yy_xzz_0_xx, g_yy_xzz_0_xy, g_yy_xzz_0_xz, g_yy_xzz_0_yy, g_yy_xzz_0_yz, g_yy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_zz_0_xx[i] = -2.0 * g_0_xzz_0_xx[i] * b_exp + 4.0 * g_yy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_xy[i] = -2.0 * g_0_xzz_0_xy[i] * b_exp + 4.0 * g_yy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_xz[i] = -2.0 * g_0_xzz_0_xz[i] * b_exp + 4.0 * g_yy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_yy[i] = -2.0 * g_0_xzz_0_yy[i] * b_exp + 4.0 * g_yy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_yz[i] = -2.0 * g_0_xzz_0_yz[i] * b_exp + 4.0 * g_yy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_zz[i] = -2.0 * g_0_xzz_0_zz[i] * b_exp + 4.0 * g_yy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_0_xx, g_y_x_0_0_z_xx_0_xy, g_y_x_0_0_z_xx_0_xz, g_y_x_0_0_z_xx_0_yy, g_y_x_0_0_z_xx_0_yz, g_y_x_0_0_z_xx_0_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xxx_0_xx, g_yz_xxx_0_xy, g_yz_xxx_0_xz, g_yz_xxx_0_yy, g_yz_xxx_0_yz, g_yz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_0_xx[i] = -4.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xxx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_xy[i] = -4.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xxx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_xz[i] = -4.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xxx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_yy[i] = -4.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xxx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_yz[i] = -4.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xxx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_zz[i] = -4.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_x_0_0_z_xy_0_xx, g_y_x_0_0_z_xy_0_xy, g_y_x_0_0_z_xy_0_xz, g_y_x_0_0_z_xy_0_yy, g_y_x_0_0_z_xy_0_yz, g_y_x_0_0_z_xy_0_zz, g_yz_xxy_0_xx, g_yz_xxy_0_xy, g_yz_xxy_0_xz, g_yz_xxy_0_yy, g_yz_xxy_0_yz, g_yz_xxy_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xy_0_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_x_0_0_z_xz_0_xx, g_y_x_0_0_z_xz_0_xy, g_y_x_0_0_z_xz_0_xz, g_y_x_0_0_z_xz_0_yy, g_y_x_0_0_z_xz_0_yz, g_y_x_0_0_z_xz_0_zz, g_yz_xxz_0_xx, g_yz_xxz_0_xy, g_yz_xxz_0_xz, g_yz_xxz_0_yy, g_yz_xxz_0_yz, g_yz_xxz_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xz_0_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_x_0_0_z_yy_0_xx, g_y_x_0_0_z_yy_0_xy, g_y_x_0_0_z_yy_0_xz, g_y_x_0_0_z_yy_0_yy, g_y_x_0_0_z_yy_0_yz, g_y_x_0_0_z_yy_0_zz, g_yz_xyy_0_xx, g_yz_xyy_0_xy, g_yz_xyy_0_xz, g_yz_xyy_0_yy, g_yz_xyy_0_yz, g_yz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yy_0_xx[i] = 4.0 * g_yz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_xy[i] = 4.0 * g_yz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_xz[i] = 4.0 * g_yz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_yy[i] = 4.0 * g_yz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_yz[i] = 4.0 * g_yz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_zz[i] = 4.0 * g_yz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_x_0_0_z_yz_0_xx, g_y_x_0_0_z_yz_0_xy, g_y_x_0_0_z_yz_0_xz, g_y_x_0_0_z_yz_0_yy, g_y_x_0_0_z_yz_0_yz, g_y_x_0_0_z_yz_0_zz, g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yz_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_x_0_0_z_zz_0_xx, g_y_x_0_0_z_zz_0_xy, g_y_x_0_0_z_zz_0_xz, g_y_x_0_0_z_zz_0_yy, g_y_x_0_0_z_zz_0_yz, g_y_x_0_0_z_zz_0_zz, g_yz_xzz_0_xx, g_yz_xzz_0_xy, g_yz_xzz_0_xz, g_yz_xzz_0_yy, g_yz_xzz_0_yz, g_yz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_zz_0_xx[i] = 4.0 * g_yz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_xy[i] = 4.0 * g_yz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_xz[i] = 4.0 * g_yz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_yy[i] = 4.0 * g_yz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_yz[i] = 4.0 * g_yz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_zz[i] = 4.0 * g_yz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xy_xxy_0_xx, g_xy_xxy_0_xy, g_xy_xxy_0_xz, g_xy_xxy_0_yy, g_xy_xxy_0_yz, g_xy_xxy_0_zz, g_y_y_0_0_x_xx_0_xx, g_y_y_0_0_x_xx_0_xy, g_y_y_0_0_x_xx_0_xz, g_y_y_0_0_x_xx_0_yy, g_y_y_0_0_x_xx_0_yz, g_y_y_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_0_xx[i] = 4.0 * g_xy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_xy[i] = 4.0 * g_xy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_xz[i] = 4.0 * g_xy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_yy[i] = 4.0 * g_xy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_yz[i] = 4.0 * g_xy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_zz[i] = 4.0 * g_xy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xyy_0_xx, g_xy_xyy_0_xy, g_xy_xyy_0_xz, g_xy_xyy_0_yy, g_xy_xyy_0_yz, g_xy_xyy_0_zz, g_y_y_0_0_x_xy_0_xx, g_y_y_0_0_x_xy_0_xy, g_y_y_0_0_x_xy_0_xz, g_y_y_0_0_x_xy_0_yy, g_y_y_0_0_x_xy_0_yz, g_y_y_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xy_0_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz, g_y_y_0_0_x_xz_0_xx, g_y_y_0_0_x_xz_0_xy, g_y_y_0_0_x_xz_0_xz, g_y_y_0_0_x_xz_0_yy, g_y_y_0_0_x_xz_0_yz, g_y_y_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xz_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_yyy_0_xx, g_xy_yyy_0_xy, g_xy_yyy_0_xz, g_xy_yyy_0_yy, g_xy_yyy_0_yz, g_xy_yyy_0_zz, g_y_y_0_0_x_yy_0_xx, g_y_y_0_0_x_yy_0_xy, g_y_y_0_0_x_yy_0_xz, g_y_y_0_0_x_yy_0_yy, g_y_y_0_0_x_yy_0_yz, g_y_y_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yy_0_xx[i] = -4.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_xy[i] = -4.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_xz[i] = -4.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_yy[i] = -4.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_yz[i] = -4.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_zz[i] = -4.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xy_yyz_0_xx, g_xy_yyz_0_xy, g_xy_yyz_0_xz, g_xy_yyz_0_yy, g_xy_yyz_0_yz, g_xy_yyz_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_y_y_0_0_x_yz_0_xx, g_y_y_0_0_x_yz_0_xy, g_y_y_0_0_x_yz_0_xz, g_y_y_0_0_x_yz_0_yy, g_y_y_0_0_x_yz_0_yz, g_y_y_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yz_0_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xy_yzz_0_xx, g_xy_yzz_0_xy, g_xy_yzz_0_xz, g_xy_yzz_0_yy, g_xy_yzz_0_yz, g_xy_yzz_0_zz, g_y_y_0_0_x_zz_0_xx, g_y_y_0_0_x_zz_0_xy, g_y_y_0_0_x_zz_0_xz, g_y_y_0_0_x_zz_0_yy, g_y_y_0_0_x_zz_0_yz, g_y_y_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_zz_0_xx[i] = 4.0 * g_xy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_xy[i] = 4.0 * g_xy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_xz[i] = 4.0 * g_xy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_yy[i] = 4.0 * g_xy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_yz[i] = 4.0 * g_xy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_zz[i] = 4.0 * g_xy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_y_y_0_0_y_xx_0_xx, g_y_y_0_0_y_xx_0_xy, g_y_y_0_0_y_xx_0_xz, g_y_y_0_0_y_xx_0_yy, g_y_y_0_0_y_xx_0_yz, g_y_y_0_0_y_xx_0_zz, g_yy_xxy_0_xx, g_yy_xxy_0_xy, g_yy_xxy_0_xz, g_yy_xxy_0_yy, g_yy_xxy_0_yz, g_yy_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_0_xx[i] = -2.0 * g_0_xxy_0_xx[i] * b_exp + 4.0 * g_yy_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_xy[i] = -2.0 * g_0_xxy_0_xy[i] * b_exp + 4.0 * g_yy_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_xz[i] = -2.0 * g_0_xxy_0_xz[i] * b_exp + 4.0 * g_yy_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_yy[i] = -2.0 * g_0_xxy_0_yy[i] * b_exp + 4.0 * g_yy_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_yz[i] = -2.0 * g_0_xxy_0_yz[i] * b_exp + 4.0 * g_yy_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_zz[i] = -2.0 * g_0_xxy_0_zz[i] * b_exp + 4.0 * g_yy_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_y_y_0_0_y_xy_0_xx, g_y_y_0_0_y_xy_0_xy, g_y_y_0_0_y_xy_0_xz, g_y_y_0_0_y_xy_0_yy, g_y_y_0_0_y_xy_0_yz, g_y_y_0_0_y_xy_0_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_xyy_0_xx, g_yy_xyy_0_xy, g_yy_xyy_0_xz, g_yy_xyy_0_yy, g_yy_xyy_0_yz, g_yy_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xy_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xyy_0_xx[i] * b_exp - 2.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xyy_0_xy[i] * b_exp - 2.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xyy_0_xz[i] * b_exp - 2.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xyy_0_yy[i] * b_exp - 2.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xyy_0_yz[i] * b_exp - 2.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xyy_0_zz[i] * b_exp - 2.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_y_y_0_0_y_xz_0_xx, g_y_y_0_0_y_xz_0_xy, g_y_y_0_0_y_xz_0_xz, g_y_y_0_0_y_xz_0_yy, g_y_y_0_0_y_xz_0_yz, g_y_y_0_0_y_xz_0_zz, g_yy_xyz_0_xx, g_yy_xyz_0_xy, g_yy_xyz_0_xz, g_yy_xyz_0_yy, g_yy_xyz_0_yz, g_yy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_yy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_yy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_yy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_yy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_yy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_yy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yyy_0_xx, g_0_yyy_0_xy, g_0_yyy_0_xz, g_0_yyy_0_yy, g_0_yyy_0_yz, g_0_yyy_0_zz, g_y_y_0_0_y_yy_0_xx, g_y_y_0_0_y_yy_0_xy, g_y_y_0_0_y_yy_0_xz, g_y_y_0_0_y_yy_0_yy, g_y_y_0_0_y_yy_0_yz, g_y_y_0_0_y_yy_0_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz, g_yy_yyy_0_xx, g_yy_yyy_0_xy, g_yy_yyy_0_xz, g_yy_yyy_0_yy, g_yy_yyy_0_yz, g_yy_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yy_0_xx[i] = 2.0 * g_0_y_0_xx[i] - 2.0 * g_0_yyy_0_xx[i] * b_exp - 4.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_xy[i] = 2.0 * g_0_y_0_xy[i] - 2.0 * g_0_yyy_0_xy[i] * b_exp - 4.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_xz[i] = 2.0 * g_0_y_0_xz[i] - 2.0 * g_0_yyy_0_xz[i] * b_exp - 4.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_yy[i] = 2.0 * g_0_y_0_yy[i] - 2.0 * g_0_yyy_0_yy[i] * b_exp - 4.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_yz[i] = 2.0 * g_0_y_0_yz[i] - 2.0 * g_0_yyy_0_yz[i] * b_exp - 4.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_zz[i] = 2.0 * g_0_y_0_zz[i] - 2.0 * g_0_yyy_0_zz[i] * b_exp - 4.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_y_y_0_0_y_yz_0_xx, g_y_y_0_0_y_yz_0_xy, g_y_y_0_0_y_yz_0_xz, g_y_y_0_0_y_yz_0_yy, g_y_y_0_0_y_yz_0_yz, g_y_y_0_0_y_yz_0_zz, g_yy_yyz_0_xx, g_yy_yyz_0_xy, g_yy_yyz_0_xz, g_yy_yyz_0_yy, g_yy_yyz_0_yz, g_yy_yyz_0_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_yyz_0_xx[i] * b_exp - 2.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_yyz_0_xy[i] * b_exp - 2.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_yyz_0_xz[i] * b_exp - 2.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_yyz_0_yy[i] * b_exp - 2.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_yyz_0_yz[i] * b_exp - 2.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_yyz_0_zz[i] * b_exp - 2.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_y_y_0_0_y_zz_0_xx, g_y_y_0_0_y_zz_0_xy, g_y_y_0_0_y_zz_0_xz, g_y_y_0_0_y_zz_0_yy, g_y_y_0_0_y_zz_0_yz, g_y_y_0_0_y_zz_0_zz, g_yy_yzz_0_xx, g_yy_yzz_0_xy, g_yy_yzz_0_xz, g_yy_yzz_0_yy, g_yy_yzz_0_yz, g_yy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_zz_0_xx[i] = -2.0 * g_0_yzz_0_xx[i] * b_exp + 4.0 * g_yy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_xy[i] = -2.0 * g_0_yzz_0_xy[i] * b_exp + 4.0 * g_yy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_xz[i] = -2.0 * g_0_yzz_0_xz[i] * b_exp + 4.0 * g_yy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_yy[i] = -2.0 * g_0_yzz_0_yy[i] * b_exp + 4.0 * g_yy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_yz[i] = -2.0 * g_0_yzz_0_yz[i] * b_exp + 4.0 * g_yy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_zz[i] = -2.0 * g_0_yzz_0_zz[i] * b_exp + 4.0 * g_yy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_0_xx, g_y_y_0_0_z_xx_0_xy, g_y_y_0_0_z_xx_0_xz, g_y_y_0_0_z_xx_0_yy, g_y_y_0_0_z_xx_0_yz, g_y_y_0_0_z_xx_0_zz, g_yz_xxy_0_xx, g_yz_xxy_0_xy, g_yz_xxy_0_xz, g_yz_xxy_0_yy, g_yz_xxy_0_yz, g_yz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_0_xx[i] = 4.0 * g_yz_xxy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_xy[i] = 4.0 * g_yz_xxy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_xz[i] = 4.0 * g_yz_xxy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_yy[i] = 4.0 * g_yz_xxy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_yz[i] = 4.0 * g_yz_xxy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_zz[i] = 4.0 * g_yz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_y_0_0_z_xy_0_xx, g_y_y_0_0_z_xy_0_xy, g_y_y_0_0_z_xy_0_xz, g_y_y_0_0_z_xy_0_yy, g_y_y_0_0_z_xy_0_yz, g_y_y_0_0_z_xy_0_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xyy_0_xx, g_yz_xyy_0_xy, g_yz_xyy_0_xz, g_yz_xyy_0_yy, g_yz_xyy_0_yz, g_yz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xy_0_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_y_0_0_z_xz_0_xx, g_y_y_0_0_z_xz_0_xy, g_y_y_0_0_z_xz_0_xz, g_y_y_0_0_z_xz_0_yy, g_y_y_0_0_z_xz_0_yz, g_y_y_0_0_z_xz_0_zz, g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xz_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_y_0_0_z_yy_0_xx, g_y_y_0_0_z_yy_0_xy, g_y_y_0_0_z_yy_0_xz, g_y_y_0_0_z_yy_0_yy, g_y_y_0_0_z_yy_0_yz, g_y_y_0_0_z_yy_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_yyy_0_xx, g_yz_yyy_0_xy, g_yz_yyy_0_xz, g_yz_yyy_0_yy, g_yz_yyy_0_yz, g_yz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yy_0_xx[i] = -4.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_yyy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_xy[i] = -4.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_yyy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_xz[i] = -4.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_yyy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_yy[i] = -4.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_yyy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_yz[i] = -4.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_yyy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_zz[i] = -4.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_y_0_0_z_yz_0_xx, g_y_y_0_0_z_yz_0_xy, g_y_y_0_0_z_yz_0_xz, g_y_y_0_0_z_yz_0_yy, g_y_y_0_0_z_yz_0_yz, g_y_y_0_0_z_yz_0_zz, g_yz_yyz_0_xx, g_yz_yyz_0_xy, g_yz_yyz_0_xz, g_yz_yyz_0_yy, g_yz_yyz_0_yz, g_yz_yyz_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yz_0_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_y_0_0_z_zz_0_xx, g_y_y_0_0_z_zz_0_xy, g_y_y_0_0_z_zz_0_xz, g_y_y_0_0_z_zz_0_yy, g_y_y_0_0_z_zz_0_yz, g_y_y_0_0_z_zz_0_zz, g_yz_yzz_0_xx, g_yz_yzz_0_xy, g_yz_yzz_0_xz, g_yz_yzz_0_yy, g_yz_yzz_0_yz, g_yz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_zz_0_xx[i] = 4.0 * g_yz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_xy[i] = 4.0 * g_yz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_xz[i] = 4.0 * g_yz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_yy[i] = 4.0 * g_yz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_yz[i] = 4.0 * g_yz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_zz[i] = 4.0 * g_yz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xy_xxz_0_xx, g_xy_xxz_0_xy, g_xy_xxz_0_xz, g_xy_xxz_0_yy, g_xy_xxz_0_yz, g_xy_xxz_0_zz, g_y_z_0_0_x_xx_0_xx, g_y_z_0_0_x_xx_0_xy, g_y_z_0_0_x_xx_0_xz, g_y_z_0_0_x_xx_0_yy, g_y_z_0_0_x_xx_0_yz, g_y_z_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_0_xx[i] = 4.0 * g_xy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_xy[i] = 4.0 * g_xy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_xz[i] = 4.0 * g_xy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_yy[i] = 4.0 * g_xy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_yz[i] = 4.0 * g_xy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_zz[i] = 4.0 * g_xy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xy_xyz_0_xx, g_xy_xyz_0_xy, g_xy_xyz_0_xz, g_xy_xyz_0_yy, g_xy_xyz_0_yz, g_xy_xyz_0_zz, g_y_z_0_0_x_xy_0_xx, g_y_z_0_0_x_xy_0_xy, g_y_z_0_0_x_xy_0_xz, g_y_z_0_0_x_xy_0_yy, g_y_z_0_0_x_xy_0_yz, g_y_z_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xy_0_xx[i] = 4.0 * g_xy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_xy[i] = 4.0 * g_xy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_xz[i] = 4.0 * g_xy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_yy[i] = 4.0 * g_xy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_yz[i] = 4.0 * g_xy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_zz[i] = 4.0 * g_xy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_xzz_0_xx, g_xy_xzz_0_xy, g_xy_xzz_0_xz, g_xy_xzz_0_yy, g_xy_xzz_0_yz, g_xy_xzz_0_zz, g_y_z_0_0_x_xz_0_xx, g_y_z_0_0_x_xz_0_xy, g_y_z_0_0_x_xz_0_xz, g_y_z_0_0_x_xz_0_yy, g_y_z_0_0_x_xz_0_yz, g_y_z_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xz_0_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xy_yyz_0_xx, g_xy_yyz_0_xy, g_xy_yyz_0_xz, g_xy_yyz_0_yy, g_xy_yyz_0_yz, g_xy_yyz_0_zz, g_y_z_0_0_x_yy_0_xx, g_y_z_0_0_x_yy_0_xy, g_y_z_0_0_x_yy_0_xz, g_y_z_0_0_x_yy_0_yy, g_y_z_0_0_x_yy_0_yz, g_y_z_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yy_0_xx[i] = 4.0 * g_xy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_xy[i] = 4.0 * g_xy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_xz[i] = 4.0 * g_xy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_yy[i] = 4.0 * g_xy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_yz[i] = 4.0 * g_xy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_zz[i] = 4.0 * g_xy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_yzz_0_xx, g_xy_yzz_0_xy, g_xy_yzz_0_xz, g_xy_yzz_0_yy, g_xy_yzz_0_yz, g_xy_yzz_0_zz, g_y_z_0_0_x_yz_0_xx, g_y_z_0_0_x_yz_0_xy, g_y_z_0_0_x_yz_0_xz, g_y_z_0_0_x_yz_0_yy, g_y_z_0_0_x_yz_0_yz, g_y_z_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yz_0_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_zzz_0_xx, g_xy_zzz_0_xy, g_xy_zzz_0_xz, g_xy_zzz_0_yy, g_xy_zzz_0_yz, g_xy_zzz_0_zz, g_y_z_0_0_x_zz_0_xx, g_y_z_0_0_x_zz_0_xy, g_y_z_0_0_x_zz_0_xz, g_y_z_0_0_x_zz_0_yy, g_y_z_0_0_x_zz_0_yz, g_y_z_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_zz_0_xx[i] = -4.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_xy[i] = -4.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_xz[i] = -4.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_yy[i] = -4.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_yz[i] = -4.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_zz[i] = -4.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_y_z_0_0_y_xx_0_xx, g_y_z_0_0_y_xx_0_xy, g_y_z_0_0_y_xx_0_xz, g_y_z_0_0_y_xx_0_yy, g_y_z_0_0_y_xx_0_yz, g_y_z_0_0_y_xx_0_zz, g_yy_xxz_0_xx, g_yy_xxz_0_xy, g_yy_xxz_0_xz, g_yy_xxz_0_yy, g_yy_xxz_0_yz, g_yy_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_0_xx[i] = -2.0 * g_0_xxz_0_xx[i] * b_exp + 4.0 * g_yy_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_xy[i] = -2.0 * g_0_xxz_0_xy[i] * b_exp + 4.0 * g_yy_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_xz[i] = -2.0 * g_0_xxz_0_xz[i] * b_exp + 4.0 * g_yy_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_yy[i] = -2.0 * g_0_xxz_0_yy[i] * b_exp + 4.0 * g_yy_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_yz[i] = -2.0 * g_0_xxz_0_yz[i] * b_exp + 4.0 * g_yy_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_zz[i] = -2.0 * g_0_xxz_0_zz[i] * b_exp + 4.0 * g_yy_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_y_z_0_0_y_xy_0_xx, g_y_z_0_0_y_xy_0_xy, g_y_z_0_0_y_xy_0_xz, g_y_z_0_0_y_xy_0_yy, g_y_z_0_0_y_xy_0_yz, g_y_z_0_0_y_xy_0_zz, g_yy_xyz_0_xx, g_yy_xyz_0_xy, g_yy_xyz_0_xz, g_yy_xyz_0_yy, g_yy_xyz_0_yz, g_yy_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xy_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_yy_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_yy_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_yy_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_yy_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_yy_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_yy_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_y_z_0_0_y_xz_0_xx, g_y_z_0_0_y_xz_0_xy, g_y_z_0_0_y_xz_0_xz, g_y_z_0_0_y_xz_0_yy, g_y_z_0_0_y_xz_0_yz, g_y_z_0_0_y_xz_0_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_xzz_0_xx, g_yy_xzz_0_xy, g_yy_xzz_0_xz, g_yy_xzz_0_yy, g_yy_xzz_0_yz, g_yy_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xz_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xzz_0_xx[i] * b_exp - 2.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xzz_0_xy[i] * b_exp - 2.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xzz_0_xz[i] * b_exp - 2.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xzz_0_yy[i] * b_exp - 2.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xzz_0_yz[i] * b_exp - 2.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xzz_0_zz[i] * b_exp - 2.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_y_z_0_0_y_yy_0_xx, g_y_z_0_0_y_yy_0_xy, g_y_z_0_0_y_yy_0_xz, g_y_z_0_0_y_yy_0_yy, g_y_z_0_0_y_yy_0_yz, g_y_z_0_0_y_yy_0_zz, g_yy_yyz_0_xx, g_yy_yyz_0_xy, g_yy_yyz_0_xz, g_yy_yyz_0_yy, g_yy_yyz_0_yz, g_yy_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yy_0_xx[i] = -2.0 * g_0_yyz_0_xx[i] * b_exp + 4.0 * g_yy_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_xy[i] = -2.0 * g_0_yyz_0_xy[i] * b_exp + 4.0 * g_yy_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_xz[i] = -2.0 * g_0_yyz_0_xz[i] * b_exp + 4.0 * g_yy_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_yy[i] = -2.0 * g_0_yyz_0_yy[i] * b_exp + 4.0 * g_yy_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_yz[i] = -2.0 * g_0_yyz_0_yz[i] * b_exp + 4.0 * g_yy_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_zz[i] = -2.0 * g_0_yyz_0_zz[i] * b_exp + 4.0 * g_yy_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_y_z_0_0_y_yz_0_xx, g_y_z_0_0_y_yz_0_xy, g_y_z_0_0_y_yz_0_xz, g_y_z_0_0_y_yz_0_yy, g_y_z_0_0_y_yz_0_yz, g_y_z_0_0_y_yz_0_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz, g_yy_yzz_0_xx, g_yy_yzz_0_xy, g_yy_yzz_0_xz, g_yy_yzz_0_yy, g_yy_yzz_0_yz, g_yy_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yz_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_yzz_0_xx[i] * b_exp - 2.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_yzz_0_xy[i] * b_exp - 2.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_yzz_0_xz[i] * b_exp - 2.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_yzz_0_yy[i] * b_exp - 2.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_yzz_0_yz[i] * b_exp - 2.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_yzz_0_zz[i] * b_exp - 2.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_zzz_0_xx, g_0_zzz_0_xy, g_0_zzz_0_xz, g_0_zzz_0_yy, g_0_zzz_0_yz, g_0_zzz_0_zz, g_y_z_0_0_y_zz_0_xx, g_y_z_0_0_y_zz_0_xy, g_y_z_0_0_y_zz_0_xz, g_y_z_0_0_y_zz_0_yy, g_y_z_0_0_y_zz_0_yz, g_y_z_0_0_y_zz_0_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz, g_yy_zzz_0_xx, g_yy_zzz_0_xy, g_yy_zzz_0_xz, g_yy_zzz_0_yy, g_yy_zzz_0_yz, g_yy_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_zz_0_xx[i] = 2.0 * g_0_z_0_xx[i] - 2.0 * g_0_zzz_0_xx[i] * b_exp - 4.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_xy[i] = 2.0 * g_0_z_0_xy[i] - 2.0 * g_0_zzz_0_xy[i] * b_exp - 4.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_xz[i] = 2.0 * g_0_z_0_xz[i] - 2.0 * g_0_zzz_0_xz[i] * b_exp - 4.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_yy[i] = 2.0 * g_0_z_0_yy[i] - 2.0 * g_0_zzz_0_yy[i] * b_exp - 4.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_yz[i] = 2.0 * g_0_z_0_yz[i] - 2.0 * g_0_zzz_0_yz[i] * b_exp - 4.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_zz[i] = 2.0 * g_0_z_0_zz[i] - 2.0 * g_0_zzz_0_zz[i] * b_exp - 4.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_0_xx, g_y_z_0_0_z_xx_0_xy, g_y_z_0_0_z_xx_0_xz, g_y_z_0_0_z_xx_0_yy, g_y_z_0_0_z_xx_0_yz, g_y_z_0_0_z_xx_0_zz, g_yz_xxz_0_xx, g_yz_xxz_0_xy, g_yz_xxz_0_xz, g_yz_xxz_0_yy, g_yz_xxz_0_yz, g_yz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_0_xx[i] = 4.0 * g_yz_xxz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_xy[i] = 4.0 * g_yz_xxz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_xz[i] = 4.0 * g_yz_xxz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_yy[i] = 4.0 * g_yz_xxz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_yz[i] = 4.0 * g_yz_xxz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_zz[i] = 4.0 * g_yz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_z_0_0_z_xy_0_xx, g_y_z_0_0_z_xy_0_xy, g_y_z_0_0_z_xy_0_xz, g_y_z_0_0_z_xy_0_yy, g_y_z_0_0_z_xy_0_yz, g_y_z_0_0_z_xy_0_zz, g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xy_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_z_0_0_z_xz_0_xx, g_y_z_0_0_z_xz_0_xy, g_y_z_0_0_z_xz_0_xz, g_y_z_0_0_z_xz_0_yy, g_y_z_0_0_z_xz_0_yz, g_y_z_0_0_z_xz_0_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xzz_0_xx, g_yz_xzz_0_xy, g_yz_xzz_0_xz, g_yz_xzz_0_yy, g_yz_xzz_0_yz, g_yz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xz_0_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_z_0_0_z_yy_0_xx, g_y_z_0_0_z_yy_0_xy, g_y_z_0_0_z_yy_0_xz, g_y_z_0_0_z_yy_0_yy, g_y_z_0_0_z_yy_0_yz, g_y_z_0_0_z_yy_0_zz, g_yz_yyz_0_xx, g_yz_yyz_0_xy, g_yz_yyz_0_xz, g_yz_yyz_0_yy, g_yz_yyz_0_yz, g_yz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yy_0_xx[i] = 4.0 * g_yz_yyz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_xy[i] = 4.0 * g_yz_yyz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_xz[i] = 4.0 * g_yz_yyz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_yy[i] = 4.0 * g_yz_yyz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_yz[i] = 4.0 * g_yz_yyz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_zz[i] = 4.0 * g_yz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_z_0_0_z_yz_0_xx, g_y_z_0_0_z_yz_0_xy, g_y_z_0_0_z_yz_0_xz, g_y_z_0_0_z_yz_0_yy, g_y_z_0_0_z_yz_0_yz, g_y_z_0_0_z_yz_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_yzz_0_xx, g_yz_yzz_0_xy, g_yz_yzz_0_xz, g_yz_yzz_0_yy, g_yz_yzz_0_yz, g_yz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yz_0_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_yzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_yzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_yzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_yzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_yzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_z_0_0_z_zz_0_xx, g_y_z_0_0_z_zz_0_xy, g_y_z_0_0_z_zz_0_xz, g_y_z_0_0_z_zz_0_yy, g_y_z_0_0_z_zz_0_yz, g_y_z_0_0_z_zz_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_zzz_0_xx, g_yz_zzz_0_xy, g_yz_zzz_0_xz, g_yz_zzz_0_yy, g_yz_zzz_0_yz, g_yz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_zz_0_xx[i] = -4.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_zzz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_xy[i] = -4.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_zzz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_xz[i] = -4.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_zzz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_yy[i] = -4.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_zzz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_yz[i] = -4.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_zzz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_zz[i] = -4.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xxx_0_xx, g_xz_xxx_0_xy, g_xz_xxx_0_xz, g_xz_xxx_0_yy, g_xz_xxx_0_yz, g_xz_xxx_0_zz, g_z_x_0_0_x_xx_0_xx, g_z_x_0_0_x_xx_0_xy, g_z_x_0_0_x_xx_0_xz, g_z_x_0_0_x_xx_0_yy, g_z_x_0_0_x_xx_0_yz, g_z_x_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_0_xx[i] = -4.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_xy[i] = -4.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_xz[i] = -4.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_yy[i] = -4.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_yz[i] = -4.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_zz[i] = -4.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xz_xxy_0_xx, g_xz_xxy_0_xy, g_xz_xxy_0_xz, g_xz_xxy_0_yy, g_xz_xxy_0_yz, g_xz_xxy_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_z_x_0_0_x_xy_0_xx, g_z_x_0_0_x_xy_0_xy, g_z_x_0_0_x_xy_0_xz, g_z_x_0_0_x_xy_0_yy, g_z_x_0_0_x_xy_0_yz, g_z_x_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xy_0_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xz_xxz_0_xx, g_xz_xxz_0_xy, g_xz_xxz_0_xz, g_xz_xxz_0_yy, g_xz_xxz_0_yz, g_xz_xxz_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_z_x_0_0_x_xz_0_xx, g_z_x_0_0_x_xz_0_xy, g_z_x_0_0_x_xz_0_xz, g_z_x_0_0_x_xz_0_yy, g_z_x_0_0_x_xz_0_yz, g_z_x_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xz_0_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xz_xyy_0_xx, g_xz_xyy_0_xy, g_xz_xyy_0_xz, g_xz_xyy_0_yy, g_xz_xyy_0_yz, g_xz_xyy_0_zz, g_z_x_0_0_x_yy_0_xx, g_z_x_0_0_x_yy_0_xy, g_z_x_0_0_x_yy_0_xz, g_z_x_0_0_x_yy_0_yy, g_z_x_0_0_x_yy_0_yz, g_z_x_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yy_0_xx[i] = 4.0 * g_xz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_xy[i] = 4.0 * g_xz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_xz[i] = 4.0 * g_xz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_yy[i] = 4.0 * g_xz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_yz[i] = 4.0 * g_xz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_zz[i] = 4.0 * g_xz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz, g_z_x_0_0_x_yz_0_xx, g_z_x_0_0_x_yz_0_xy, g_z_x_0_0_x_yz_0_xz, g_z_x_0_0_x_yz_0_yy, g_z_x_0_0_x_yz_0_yz, g_z_x_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yz_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xz_xzz_0_xx, g_xz_xzz_0_xy, g_xz_xzz_0_xz, g_xz_xzz_0_yy, g_xz_xzz_0_yz, g_xz_xzz_0_zz, g_z_x_0_0_x_zz_0_xx, g_z_x_0_0_x_zz_0_xy, g_z_x_0_0_x_zz_0_xz, g_z_x_0_0_x_zz_0_yy, g_z_x_0_0_x_zz_0_yz, g_z_x_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_zz_0_xx[i] = 4.0 * g_xz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_xy[i] = 4.0 * g_xz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_xz[i] = 4.0 * g_xz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_yy[i] = 4.0 * g_xz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_yz[i] = 4.0 * g_xz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_zz[i] = 4.0 * g_xz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xxx_0_xx, g_yz_xxx_0_xy, g_yz_xxx_0_xz, g_yz_xxx_0_yy, g_yz_xxx_0_yz, g_yz_xxx_0_zz, g_z_x_0_0_y_xx_0_xx, g_z_x_0_0_y_xx_0_xy, g_z_x_0_0_y_xx_0_xz, g_z_x_0_0_y_xx_0_yy, g_z_x_0_0_y_xx_0_yz, g_z_x_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_0_xx[i] = -4.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_xy[i] = -4.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_xz[i] = -4.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_yy[i] = -4.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_yz[i] = -4.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_zz[i] = -4.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_yz_xxy_0_xx, g_yz_xxy_0_xy, g_yz_xxy_0_xz, g_yz_xxy_0_yy, g_yz_xxy_0_yz, g_yz_xxy_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_z_x_0_0_y_xy_0_xx, g_z_x_0_0_y_xy_0_xy, g_z_x_0_0_y_xy_0_xz, g_z_x_0_0_y_xy_0_yy, g_z_x_0_0_y_xy_0_yz, g_z_x_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xy_0_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_yz_xxz_0_xx, g_yz_xxz_0_xy, g_yz_xxz_0_xz, g_yz_xxz_0_yy, g_yz_xxz_0_yz, g_yz_xxz_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_z_x_0_0_y_xz_0_xx, g_z_x_0_0_y_xz_0_xy, g_z_x_0_0_y_xz_0_xz, g_z_x_0_0_y_xz_0_yy, g_z_x_0_0_y_xz_0_yz, g_z_x_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xz_0_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_yz_xyy_0_xx, g_yz_xyy_0_xy, g_yz_xyy_0_xz, g_yz_xyy_0_yy, g_yz_xyy_0_yz, g_yz_xyy_0_zz, g_z_x_0_0_y_yy_0_xx, g_z_x_0_0_y_yy_0_xy, g_z_x_0_0_y_yy_0_xz, g_z_x_0_0_y_yy_0_yy, g_z_x_0_0_y_yy_0_yz, g_z_x_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yy_0_xx[i] = 4.0 * g_yz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_xy[i] = 4.0 * g_yz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_xz[i] = 4.0 * g_yz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_yy[i] = 4.0 * g_yz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_yz[i] = 4.0 * g_yz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_zz[i] = 4.0 * g_yz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz, g_z_x_0_0_y_yz_0_xx, g_z_x_0_0_y_yz_0_xy, g_z_x_0_0_y_yz_0_xz, g_z_x_0_0_y_yz_0_yy, g_z_x_0_0_y_yz_0_yz, g_z_x_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yz_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_yz_xzz_0_xx, g_yz_xzz_0_xy, g_yz_xzz_0_xz, g_yz_xzz_0_yy, g_yz_xzz_0_yz, g_yz_xzz_0_zz, g_z_x_0_0_y_zz_0_xx, g_z_x_0_0_y_zz_0_xy, g_z_x_0_0_y_zz_0_xz, g_z_x_0_0_y_zz_0_yy, g_z_x_0_0_y_zz_0_yz, g_z_x_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_zz_0_xx[i] = 4.0 * g_yz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_xy[i] = 4.0 * g_yz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_xz[i] = 4.0 * g_yz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_yy[i] = 4.0 * g_yz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_yz[i] = 4.0 * g_yz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_zz[i] = 4.0 * g_yz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xxx_0_xx, g_0_xxx_0_xy, g_0_xxx_0_xz, g_0_xxx_0_yy, g_0_xxx_0_yz, g_0_xxx_0_zz, g_z_x_0_0_z_xx_0_xx, g_z_x_0_0_z_xx_0_xy, g_z_x_0_0_z_xx_0_xz, g_z_x_0_0_z_xx_0_yy, g_z_x_0_0_z_xx_0_yz, g_z_x_0_0_z_xx_0_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_xxx_0_xx, g_zz_xxx_0_xy, g_zz_xxx_0_xz, g_zz_xxx_0_yy, g_zz_xxx_0_yz, g_zz_xxx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_0_xx[i] = 2.0 * g_0_x_0_xx[i] - 2.0 * g_0_xxx_0_xx[i] * b_exp - 4.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_xxx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_xy[i] = 2.0 * g_0_x_0_xy[i] - 2.0 * g_0_xxx_0_xy[i] * b_exp - 4.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_xxx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_xz[i] = 2.0 * g_0_x_0_xz[i] - 2.0 * g_0_xxx_0_xz[i] * b_exp - 4.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_xxx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_yy[i] = 2.0 * g_0_x_0_yy[i] - 2.0 * g_0_xxx_0_yy[i] * b_exp - 4.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_xxx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_yz[i] = 2.0 * g_0_x_0_yz[i] - 2.0 * g_0_xxx_0_yz[i] * b_exp - 4.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_xxx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_zz[i] = 2.0 * g_0_x_0_zz[i] - 2.0 * g_0_xxx_0_zz[i] * b_exp - 4.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_xxx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_z_x_0_0_z_xy_0_xx, g_z_x_0_0_z_xy_0_xy, g_z_x_0_0_z_xy_0_xz, g_z_x_0_0_z_xy_0_yy, g_z_x_0_0_z_xy_0_yz, g_z_x_0_0_z_xy_0_zz, g_zz_xxy_0_xx, g_zz_xxy_0_xy, g_zz_xxy_0_xz, g_zz_xxy_0_yy, g_zz_xxy_0_yz, g_zz_xxy_0_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xy_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_xxy_0_xx[i] * b_exp - 2.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_xxy_0_xy[i] * b_exp - 2.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_xxy_0_xz[i] * b_exp - 2.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_xxy_0_yy[i] * b_exp - 2.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_xxy_0_yz[i] * b_exp - 2.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_xxy_0_zz[i] * b_exp - 2.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_z_x_0_0_z_xz_0_xx, g_z_x_0_0_z_xz_0_xy, g_z_x_0_0_z_xz_0_xz, g_z_x_0_0_z_xz_0_yy, g_z_x_0_0_z_xz_0_yz, g_z_x_0_0_z_xz_0_zz, g_zz_xxz_0_xx, g_zz_xxz_0_xy, g_zz_xxz_0_xz, g_zz_xxz_0_yy, g_zz_xxz_0_yz, g_zz_xxz_0_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_xxz_0_xx[i] * b_exp - 2.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_xxz_0_xy[i] * b_exp - 2.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_xxz_0_xz[i] * b_exp - 2.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_xxz_0_yy[i] * b_exp - 2.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_xxz_0_yz[i] * b_exp - 2.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_xxz_0_zz[i] * b_exp - 2.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_z_x_0_0_z_yy_0_xx, g_z_x_0_0_z_yy_0_xy, g_z_x_0_0_z_yy_0_xz, g_z_x_0_0_z_yy_0_yy, g_z_x_0_0_z_yy_0_yz, g_z_x_0_0_z_yy_0_zz, g_zz_xyy_0_xx, g_zz_xyy_0_xy, g_zz_xyy_0_xz, g_zz_xyy_0_yy, g_zz_xyy_0_yz, g_zz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yy_0_xx[i] = -2.0 * g_0_xyy_0_xx[i] * b_exp + 4.0 * g_zz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_xy[i] = -2.0 * g_0_xyy_0_xy[i] * b_exp + 4.0 * g_zz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_xz[i] = -2.0 * g_0_xyy_0_xz[i] * b_exp + 4.0 * g_zz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_yy[i] = -2.0 * g_0_xyy_0_yy[i] * b_exp + 4.0 * g_zz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_yz[i] = -2.0 * g_0_xyy_0_yz[i] * b_exp + 4.0 * g_zz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_zz[i] = -2.0 * g_0_xyy_0_zz[i] * b_exp + 4.0 * g_zz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_z_x_0_0_z_yz_0_xx, g_z_x_0_0_z_yz_0_xy, g_z_x_0_0_z_yz_0_xz, g_z_x_0_0_z_yz_0_yy, g_z_x_0_0_z_yz_0_yz, g_z_x_0_0_z_yz_0_zz, g_zz_xyz_0_xx, g_zz_xyz_0_xy, g_zz_xyz_0_xz, g_zz_xyz_0_yy, g_zz_xyz_0_yz, g_zz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_zz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_zz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_zz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_zz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_zz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_zz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_z_x_0_0_z_zz_0_xx, g_z_x_0_0_z_zz_0_xy, g_z_x_0_0_z_zz_0_xz, g_z_x_0_0_z_zz_0_yy, g_z_x_0_0_z_zz_0_yz, g_z_x_0_0_z_zz_0_zz, g_zz_xzz_0_xx, g_zz_xzz_0_xy, g_zz_xzz_0_xz, g_zz_xzz_0_yy, g_zz_xzz_0_yz, g_zz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_zz_0_xx[i] = -2.0 * g_0_xzz_0_xx[i] * b_exp + 4.0 * g_zz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_xy[i] = -2.0 * g_0_xzz_0_xy[i] * b_exp + 4.0 * g_zz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_xz[i] = -2.0 * g_0_xzz_0_xz[i] * b_exp + 4.0 * g_zz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_yy[i] = -2.0 * g_0_xzz_0_yy[i] * b_exp + 4.0 * g_zz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_yz[i] = -2.0 * g_0_xzz_0_yz[i] * b_exp + 4.0 * g_zz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_zz[i] = -2.0 * g_0_xzz_0_zz[i] * b_exp + 4.0 * g_zz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_xz_xxy_0_xx, g_xz_xxy_0_xy, g_xz_xxy_0_xz, g_xz_xxy_0_yy, g_xz_xxy_0_yz, g_xz_xxy_0_zz, g_z_y_0_0_x_xx_0_xx, g_z_y_0_0_x_xx_0_xy, g_z_y_0_0_x_xx_0_xz, g_z_y_0_0_x_xx_0_yy, g_z_y_0_0_x_xx_0_yz, g_z_y_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_0_xx[i] = 4.0 * g_xz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_xy[i] = 4.0 * g_xz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_xz[i] = 4.0 * g_xz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_yy[i] = 4.0 * g_xz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_yz[i] = 4.0 * g_xz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_zz[i] = 4.0 * g_xz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xyy_0_xx, g_xz_xyy_0_xy, g_xz_xyy_0_xz, g_xz_xyy_0_yy, g_xz_xyy_0_yz, g_xz_xyy_0_zz, g_z_y_0_0_x_xy_0_xx, g_z_y_0_0_x_xy_0_xy, g_z_y_0_0_x_xy_0_xz, g_z_y_0_0_x_xy_0_yy, g_z_y_0_0_x_xy_0_yz, g_z_y_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xy_0_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz, g_z_y_0_0_x_xz_0_xx, g_z_y_0_0_x_xz_0_xy, g_z_y_0_0_x_xz_0_xz, g_z_y_0_0_x_xz_0_yy, g_z_y_0_0_x_xz_0_yz, g_z_y_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xz_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_yyy_0_xx, g_xz_yyy_0_xy, g_xz_yyy_0_xz, g_xz_yyy_0_yy, g_xz_yyy_0_yz, g_xz_yyy_0_zz, g_z_y_0_0_x_yy_0_xx, g_z_y_0_0_x_yy_0_xy, g_z_y_0_0_x_yy_0_xz, g_z_y_0_0_x_yy_0_yy, g_z_y_0_0_x_yy_0_yz, g_z_y_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yy_0_xx[i] = -4.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_xy[i] = -4.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_xz[i] = -4.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_yy[i] = -4.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_yz[i] = -4.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_zz[i] = -4.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_xz_yyz_0_xx, g_xz_yyz_0_xy, g_xz_yyz_0_xz, g_xz_yyz_0_yy, g_xz_yyz_0_yz, g_xz_yyz_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_z_y_0_0_x_yz_0_xx, g_z_y_0_0_x_yz_0_xy, g_z_y_0_0_x_yz_0_xz, g_z_y_0_0_x_yz_0_yy, g_z_y_0_0_x_yz_0_yz, g_z_y_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yz_0_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_xz_yzz_0_xx, g_xz_yzz_0_xy, g_xz_yzz_0_xz, g_xz_yzz_0_yy, g_xz_yzz_0_yz, g_xz_yzz_0_zz, g_z_y_0_0_x_zz_0_xx, g_z_y_0_0_x_zz_0_xy, g_z_y_0_0_x_zz_0_xz, g_z_y_0_0_x_zz_0_yy, g_z_y_0_0_x_zz_0_yz, g_z_y_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_zz_0_xx[i] = 4.0 * g_xz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_xy[i] = 4.0 * g_xz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_xz[i] = 4.0 * g_xz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_yy[i] = 4.0 * g_xz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_yz[i] = 4.0 * g_xz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_zz[i] = 4.0 * g_xz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_yz_xxy_0_xx, g_yz_xxy_0_xy, g_yz_xxy_0_xz, g_yz_xxy_0_yy, g_yz_xxy_0_yz, g_yz_xxy_0_zz, g_z_y_0_0_y_xx_0_xx, g_z_y_0_0_y_xx_0_xy, g_z_y_0_0_y_xx_0_xz, g_z_y_0_0_y_xx_0_yy, g_z_y_0_0_y_xx_0_yz, g_z_y_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_0_xx[i] = 4.0 * g_yz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_xy[i] = 4.0 * g_yz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_xz[i] = 4.0 * g_yz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_yy[i] = 4.0 * g_yz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_yz[i] = 4.0 * g_yz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_zz[i] = 4.0 * g_yz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xyy_0_xx, g_yz_xyy_0_xy, g_yz_xyy_0_xz, g_yz_xyy_0_yy, g_yz_xyy_0_yz, g_yz_xyy_0_zz, g_z_y_0_0_y_xy_0_xx, g_z_y_0_0_y_xy_0_xy, g_z_y_0_0_y_xy_0_xz, g_z_y_0_0_y_xy_0_yy, g_z_y_0_0_y_xy_0_yz, g_z_y_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xy_0_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz, g_z_y_0_0_y_xz_0_xx, g_z_y_0_0_y_xz_0_xy, g_z_y_0_0_y_xz_0_xz, g_z_y_0_0_y_xz_0_yy, g_z_y_0_0_y_xz_0_yz, g_z_y_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xz_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_yyy_0_xx, g_yz_yyy_0_xy, g_yz_yyy_0_xz, g_yz_yyy_0_yy, g_yz_yyy_0_yz, g_yz_yyy_0_zz, g_z_y_0_0_y_yy_0_xx, g_z_y_0_0_y_yy_0_xy, g_z_y_0_0_y_yy_0_xz, g_z_y_0_0_y_yy_0_yy, g_z_y_0_0_y_yy_0_yz, g_z_y_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yy_0_xx[i] = -4.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_xy[i] = -4.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_xz[i] = -4.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_yy[i] = -4.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_yz[i] = -4.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_zz[i] = -4.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_yz_yyz_0_xx, g_yz_yyz_0_xy, g_yz_yyz_0_xz, g_yz_yyz_0_yy, g_yz_yyz_0_yz, g_yz_yyz_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_z_y_0_0_y_yz_0_xx, g_z_y_0_0_y_yz_0_xy, g_z_y_0_0_y_yz_0_xz, g_z_y_0_0_y_yz_0_yy, g_z_y_0_0_y_yz_0_yz, g_z_y_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yz_0_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_yz_yzz_0_xx, g_yz_yzz_0_xy, g_yz_yzz_0_xz, g_yz_yzz_0_yy, g_yz_yzz_0_yz, g_yz_yzz_0_zz, g_z_y_0_0_y_zz_0_xx, g_z_y_0_0_y_zz_0_xy, g_z_y_0_0_y_zz_0_xz, g_z_y_0_0_y_zz_0_yy, g_z_y_0_0_y_zz_0_yz, g_z_y_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_zz_0_xx[i] = 4.0 * g_yz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_xy[i] = 4.0 * g_yz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_xz[i] = 4.0 * g_yz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_yy[i] = 4.0 * g_yz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_yz[i] = 4.0 * g_yz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_zz[i] = 4.0 * g_yz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_0_xxy_0_xx, g_0_xxy_0_xy, g_0_xxy_0_xz, g_0_xxy_0_yy, g_0_xxy_0_yz, g_0_xxy_0_zz, g_z_y_0_0_z_xx_0_xx, g_z_y_0_0_z_xx_0_xy, g_z_y_0_0_z_xx_0_xz, g_z_y_0_0_z_xx_0_yy, g_z_y_0_0_z_xx_0_yz, g_z_y_0_0_z_xx_0_zz, g_zz_xxy_0_xx, g_zz_xxy_0_xy, g_zz_xxy_0_xz, g_zz_xxy_0_yy, g_zz_xxy_0_yz, g_zz_xxy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_0_xx[i] = -2.0 * g_0_xxy_0_xx[i] * b_exp + 4.0 * g_zz_xxy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_xy[i] = -2.0 * g_0_xxy_0_xy[i] * b_exp + 4.0 * g_zz_xxy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_xz[i] = -2.0 * g_0_xxy_0_xz[i] * b_exp + 4.0 * g_zz_xxy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_yy[i] = -2.0 * g_0_xxy_0_yy[i] * b_exp + 4.0 * g_zz_xxy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_yz[i] = -2.0 * g_0_xxy_0_yz[i] * b_exp + 4.0 * g_zz_xxy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_zz[i] = -2.0 * g_0_xxy_0_zz[i] * b_exp + 4.0 * g_zz_xxy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xyy_0_xx, g_0_xyy_0_xy, g_0_xyy_0_xz, g_0_xyy_0_yy, g_0_xyy_0_yz, g_0_xyy_0_zz, g_z_y_0_0_z_xy_0_xx, g_z_y_0_0_z_xy_0_xy, g_z_y_0_0_z_xy_0_xz, g_z_y_0_0_z_xy_0_yy, g_z_y_0_0_z_xy_0_yz, g_z_y_0_0_z_xy_0_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_xyy_0_xx, g_zz_xyy_0_xy, g_zz_xyy_0_xz, g_zz_xyy_0_yy, g_zz_xyy_0_yz, g_zz_xyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xy_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xyy_0_xx[i] * b_exp - 2.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_xyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xyy_0_xy[i] * b_exp - 2.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_xyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xyy_0_xz[i] * b_exp - 2.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_xyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xyy_0_yy[i] * b_exp - 2.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_xyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xyy_0_yz[i] * b_exp - 2.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_xyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xyy_0_zz[i] * b_exp - 2.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_xyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_z_y_0_0_z_xz_0_xx, g_z_y_0_0_z_xz_0_xy, g_z_y_0_0_z_xz_0_xz, g_z_y_0_0_z_xz_0_yy, g_z_y_0_0_z_xz_0_yz, g_z_y_0_0_z_xz_0_zz, g_zz_xyz_0_xx, g_zz_xyz_0_xy, g_zz_xyz_0_xz, g_zz_xyz_0_yy, g_zz_xyz_0_yz, g_zz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xz_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_zz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_zz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_zz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_zz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_zz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_zz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yyy_0_xx, g_0_yyy_0_xy, g_0_yyy_0_xz, g_0_yyy_0_yy, g_0_yyy_0_yz, g_0_yyy_0_zz, g_z_y_0_0_z_yy_0_xx, g_z_y_0_0_z_yy_0_xy, g_z_y_0_0_z_yy_0_xz, g_z_y_0_0_z_yy_0_yy, g_z_y_0_0_z_yy_0_yz, g_z_y_0_0_z_yy_0_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz, g_zz_yyy_0_xx, g_zz_yyy_0_xy, g_zz_yyy_0_xz, g_zz_yyy_0_yy, g_zz_yyy_0_yz, g_zz_yyy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yy_0_xx[i] = 2.0 * g_0_y_0_xx[i] - 2.0 * g_0_yyy_0_xx[i] * b_exp - 4.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_yyy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_xy[i] = 2.0 * g_0_y_0_xy[i] - 2.0 * g_0_yyy_0_xy[i] * b_exp - 4.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_yyy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_xz[i] = 2.0 * g_0_y_0_xz[i] - 2.0 * g_0_yyy_0_xz[i] * b_exp - 4.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_yyy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_yy[i] = 2.0 * g_0_y_0_yy[i] - 2.0 * g_0_yyy_0_yy[i] * b_exp - 4.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_yyy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_yz[i] = 2.0 * g_0_y_0_yz[i] - 2.0 * g_0_yyy_0_yz[i] * b_exp - 4.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_yyy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_zz[i] = 2.0 * g_0_y_0_zz[i] - 2.0 * g_0_yyy_0_zz[i] * b_exp - 4.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_yyy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_z_y_0_0_z_yz_0_xx, g_z_y_0_0_z_yz_0_xy, g_z_y_0_0_z_yz_0_xz, g_z_y_0_0_z_yz_0_yy, g_z_y_0_0_z_yz_0_yz, g_z_y_0_0_z_yz_0_zz, g_zz_yyz_0_xx, g_zz_yyz_0_xy, g_zz_yyz_0_xz, g_zz_yyz_0_yy, g_zz_yyz_0_yz, g_zz_yyz_0_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yz_0_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_yyz_0_xx[i] * b_exp - 2.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_yyz_0_xy[i] * b_exp - 2.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_yyz_0_xz[i] * b_exp - 2.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_yyz_0_yy[i] * b_exp - 2.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_yyz_0_yz[i] * b_exp - 2.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_yyz_0_zz[i] * b_exp - 2.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_z_y_0_0_z_zz_0_xx, g_z_y_0_0_z_zz_0_xy, g_z_y_0_0_z_zz_0_xz, g_z_y_0_0_z_zz_0_yy, g_z_y_0_0_z_zz_0_yz, g_z_y_0_0_z_zz_0_zz, g_zz_yzz_0_xx, g_zz_yzz_0_xy, g_zz_yzz_0_xz, g_zz_yzz_0_yy, g_zz_yzz_0_yz, g_zz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_zz_0_xx[i] = -2.0 * g_0_yzz_0_xx[i] * b_exp + 4.0 * g_zz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_xy[i] = -2.0 * g_0_yzz_0_xy[i] * b_exp + 4.0 * g_zz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_xz[i] = -2.0 * g_0_yzz_0_xz[i] * b_exp + 4.0 * g_zz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_yy[i] = -2.0 * g_0_yzz_0_yy[i] * b_exp + 4.0 * g_zz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_yz[i] = -2.0 * g_0_yzz_0_yz[i] * b_exp + 4.0 * g_zz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_zz[i] = -2.0 * g_0_yzz_0_zz[i] * b_exp + 4.0 * g_zz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xz_xxz_0_xx, g_xz_xxz_0_xy, g_xz_xxz_0_xz, g_xz_xxz_0_yy, g_xz_xxz_0_yz, g_xz_xxz_0_zz, g_z_z_0_0_x_xx_0_xx, g_z_z_0_0_x_xx_0_xy, g_z_z_0_0_x_xx_0_xz, g_z_z_0_0_x_xx_0_yy, g_z_z_0_0_x_xx_0_yz, g_z_z_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_0_xx[i] = 4.0 * g_xz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_xy[i] = 4.0 * g_xz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_xz[i] = 4.0 * g_xz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_yy[i] = 4.0 * g_xz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_yz[i] = 4.0 * g_xz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_zz[i] = 4.0 * g_xz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xz_xyz_0_xx, g_xz_xyz_0_xy, g_xz_xyz_0_xz, g_xz_xyz_0_yy, g_xz_xyz_0_yz, g_xz_xyz_0_zz, g_z_z_0_0_x_xy_0_xx, g_z_z_0_0_x_xy_0_xy, g_z_z_0_0_x_xy_0_xz, g_z_z_0_0_x_xy_0_yy, g_z_z_0_0_x_xy_0_yz, g_z_z_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xy_0_xx[i] = 4.0 * g_xz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_xy[i] = 4.0 * g_xz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_xz[i] = 4.0 * g_xz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_yy[i] = 4.0 * g_xz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_yz[i] = 4.0 * g_xz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_zz[i] = 4.0 * g_xz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_xzz_0_xx, g_xz_xzz_0_xy, g_xz_xzz_0_xz, g_xz_xzz_0_yy, g_xz_xzz_0_yz, g_xz_xzz_0_zz, g_z_z_0_0_x_xz_0_xx, g_z_z_0_0_x_xz_0_xy, g_z_z_0_0_x_xz_0_xz, g_z_z_0_0_x_xz_0_yy, g_z_z_0_0_x_xz_0_yz, g_z_z_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xz_0_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xz_yyz_0_xx, g_xz_yyz_0_xy, g_xz_yyz_0_xz, g_xz_yyz_0_yy, g_xz_yyz_0_yz, g_xz_yyz_0_zz, g_z_z_0_0_x_yy_0_xx, g_z_z_0_0_x_yy_0_xy, g_z_z_0_0_x_yy_0_xz, g_z_z_0_0_x_yy_0_yy, g_z_z_0_0_x_yy_0_yz, g_z_z_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yy_0_xx[i] = 4.0 * g_xz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_xy[i] = 4.0 * g_xz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_xz[i] = 4.0 * g_xz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_yy[i] = 4.0 * g_xz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_yz[i] = 4.0 * g_xz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_zz[i] = 4.0 * g_xz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_yzz_0_xx, g_xz_yzz_0_xy, g_xz_yzz_0_xz, g_xz_yzz_0_yy, g_xz_yzz_0_yz, g_xz_yzz_0_zz, g_z_z_0_0_x_yz_0_xx, g_z_z_0_0_x_yz_0_xy, g_z_z_0_0_x_yz_0_xz, g_z_z_0_0_x_yz_0_yy, g_z_z_0_0_x_yz_0_yz, g_z_z_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yz_0_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_zzz_0_xx, g_xz_zzz_0_xy, g_xz_zzz_0_xz, g_xz_zzz_0_yy, g_xz_zzz_0_yz, g_xz_zzz_0_zz, g_z_z_0_0_x_zz_0_xx, g_z_z_0_0_x_zz_0_xy, g_z_z_0_0_x_zz_0_xz, g_z_z_0_0_x_zz_0_yy, g_z_z_0_0_x_zz_0_yz, g_z_z_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_zz_0_xx[i] = -4.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_xy[i] = -4.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_xz[i] = -4.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_yy[i] = -4.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_yz[i] = -4.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_zz[i] = -4.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_yz_xxz_0_xx, g_yz_xxz_0_xy, g_yz_xxz_0_xz, g_yz_xxz_0_yy, g_yz_xxz_0_yz, g_yz_xxz_0_zz, g_z_z_0_0_y_xx_0_xx, g_z_z_0_0_y_xx_0_xy, g_z_z_0_0_y_xx_0_xz, g_z_z_0_0_y_xx_0_yy, g_z_z_0_0_y_xx_0_yz, g_z_z_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_0_xx[i] = 4.0 * g_yz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_xy[i] = 4.0 * g_yz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_xz[i] = 4.0 * g_yz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_yy[i] = 4.0 * g_yz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_yz[i] = 4.0 * g_yz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_zz[i] = 4.0 * g_yz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_yz_xyz_0_xx, g_yz_xyz_0_xy, g_yz_xyz_0_xz, g_yz_xyz_0_yy, g_yz_xyz_0_yz, g_yz_xyz_0_zz, g_z_z_0_0_y_xy_0_xx, g_z_z_0_0_y_xy_0_xy, g_z_z_0_0_y_xy_0_xz, g_z_z_0_0_y_xy_0_yy, g_z_z_0_0_y_xy_0_yz, g_z_z_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xy_0_xx[i] = 4.0 * g_yz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_xy[i] = 4.0 * g_yz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_xz[i] = 4.0 * g_yz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_yy[i] = 4.0 * g_yz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_yz[i] = 4.0 * g_yz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_zz[i] = 4.0 * g_yz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_xzz_0_xx, g_yz_xzz_0_xy, g_yz_xzz_0_xz, g_yz_xzz_0_yy, g_yz_xzz_0_yz, g_yz_xzz_0_zz, g_z_z_0_0_y_xz_0_xx, g_z_z_0_0_y_xz_0_xy, g_z_z_0_0_y_xz_0_xz, g_z_z_0_0_y_xz_0_yy, g_z_z_0_0_y_xz_0_yz, g_z_z_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xz_0_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_yz_yyz_0_xx, g_yz_yyz_0_xy, g_yz_yyz_0_xz, g_yz_yyz_0_yy, g_yz_yyz_0_yz, g_yz_yyz_0_zz, g_z_z_0_0_y_yy_0_xx, g_z_z_0_0_y_yy_0_xy, g_z_z_0_0_y_yy_0_xz, g_z_z_0_0_y_yy_0_yy, g_z_z_0_0_y_yy_0_yz, g_z_z_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yy_0_xx[i] = 4.0 * g_yz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_xy[i] = 4.0 * g_yz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_xz[i] = 4.0 * g_yz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_yy[i] = 4.0 * g_yz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_yz[i] = 4.0 * g_yz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_zz[i] = 4.0 * g_yz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_yzz_0_xx, g_yz_yzz_0_xy, g_yz_yzz_0_xz, g_yz_yzz_0_yy, g_yz_yzz_0_yz, g_yz_yzz_0_zz, g_z_z_0_0_y_yz_0_xx, g_z_z_0_0_y_yz_0_xy, g_z_z_0_0_y_yz_0_xz, g_z_z_0_0_y_yz_0_yy, g_z_z_0_0_y_yz_0_yz, g_z_z_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yz_0_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_zzz_0_xx, g_yz_zzz_0_xy, g_yz_zzz_0_xz, g_yz_zzz_0_yy, g_yz_zzz_0_yz, g_yz_zzz_0_zz, g_z_z_0_0_y_zz_0_xx, g_z_z_0_0_y_zz_0_xy, g_z_z_0_0_y_zz_0_xz, g_z_z_0_0_y_zz_0_yy, g_z_z_0_0_y_zz_0_yz, g_z_z_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_zz_0_xx[i] = -4.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_xy[i] = -4.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_xz[i] = -4.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_yy[i] = -4.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_yz[i] = -4.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_zz[i] = -4.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_zzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_0_xxz_0_xx, g_0_xxz_0_xy, g_0_xxz_0_xz, g_0_xxz_0_yy, g_0_xxz_0_yz, g_0_xxz_0_zz, g_z_z_0_0_z_xx_0_xx, g_z_z_0_0_z_xx_0_xy, g_z_z_0_0_z_xx_0_xz, g_z_z_0_0_z_xx_0_yy, g_z_z_0_0_z_xx_0_yz, g_z_z_0_0_z_xx_0_zz, g_zz_xxz_0_xx, g_zz_xxz_0_xy, g_zz_xxz_0_xz, g_zz_xxz_0_yy, g_zz_xxz_0_yz, g_zz_xxz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_0_xx[i] = -2.0 * g_0_xxz_0_xx[i] * b_exp + 4.0 * g_zz_xxz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_xy[i] = -2.0 * g_0_xxz_0_xy[i] * b_exp + 4.0 * g_zz_xxz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_xz[i] = -2.0 * g_0_xxz_0_xz[i] * b_exp + 4.0 * g_zz_xxz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_yy[i] = -2.0 * g_0_xxz_0_yy[i] * b_exp + 4.0 * g_zz_xxz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_yz[i] = -2.0 * g_0_xxz_0_yz[i] * b_exp + 4.0 * g_zz_xxz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_zz[i] = -2.0 * g_0_xxz_0_zz[i] * b_exp + 4.0 * g_zz_xxz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_0_xyz_0_xx, g_0_xyz_0_xy, g_0_xyz_0_xz, g_0_xyz_0_yy, g_0_xyz_0_yz, g_0_xyz_0_zz, g_z_z_0_0_z_xy_0_xx, g_z_z_0_0_z_xy_0_xy, g_z_z_0_0_z_xy_0_xz, g_z_z_0_0_z_xy_0_yy, g_z_z_0_0_z_xy_0_yz, g_z_z_0_0_z_xy_0_zz, g_zz_xyz_0_xx, g_zz_xyz_0_xy, g_zz_xyz_0_xz, g_zz_xyz_0_yy, g_zz_xyz_0_yz, g_zz_xyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xy_0_xx[i] = -2.0 * g_0_xyz_0_xx[i] * b_exp + 4.0 * g_zz_xyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_xy[i] = -2.0 * g_0_xyz_0_xy[i] * b_exp + 4.0 * g_zz_xyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_xz[i] = -2.0 * g_0_xyz_0_xz[i] * b_exp + 4.0 * g_zz_xyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_yy[i] = -2.0 * g_0_xyz_0_yy[i] * b_exp + 4.0 * g_zz_xyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_yz[i] = -2.0 * g_0_xyz_0_yz[i] * b_exp + 4.0 * g_zz_xyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_zz[i] = -2.0 * g_0_xyz_0_zz[i] * b_exp + 4.0 * g_zz_xyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xzz_0_xx, g_0_xzz_0_xy, g_0_xzz_0_xz, g_0_xzz_0_yy, g_0_xzz_0_yz, g_0_xzz_0_zz, g_z_z_0_0_z_xz_0_xx, g_z_z_0_0_z_xz_0_xy, g_z_z_0_0_z_xz_0_xz, g_z_z_0_0_z_xz_0_yy, g_z_z_0_0_z_xz_0_yz, g_z_z_0_0_z_xz_0_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_xzz_0_xx, g_zz_xzz_0_xy, g_zz_xzz_0_xz, g_zz_xzz_0_yy, g_zz_xzz_0_yz, g_zz_xzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xz_0_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_xzz_0_xx[i] * b_exp - 2.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_xzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_xzz_0_xy[i] * b_exp - 2.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_xzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_xzz_0_xz[i] * b_exp - 2.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_xzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_xzz_0_yy[i] * b_exp - 2.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_xzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_xzz_0_yz[i] * b_exp - 2.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_xzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_xzz_0_zz[i] * b_exp - 2.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_xzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_0_yyz_0_xx, g_0_yyz_0_xy, g_0_yyz_0_xz, g_0_yyz_0_yy, g_0_yyz_0_yz, g_0_yyz_0_zz, g_z_z_0_0_z_yy_0_xx, g_z_z_0_0_z_yy_0_xy, g_z_z_0_0_z_yy_0_xz, g_z_z_0_0_z_yy_0_yy, g_z_z_0_0_z_yy_0_yz, g_z_z_0_0_z_yy_0_zz, g_zz_yyz_0_xx, g_zz_yyz_0_xy, g_zz_yyz_0_xz, g_zz_yyz_0_yy, g_zz_yyz_0_yz, g_zz_yyz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yy_0_xx[i] = -2.0 * g_0_yyz_0_xx[i] * b_exp + 4.0 * g_zz_yyz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_xy[i] = -2.0 * g_0_yyz_0_xy[i] * b_exp + 4.0 * g_zz_yyz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_xz[i] = -2.0 * g_0_yyz_0_xz[i] * b_exp + 4.0 * g_zz_yyz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_yy[i] = -2.0 * g_0_yyz_0_yy[i] * b_exp + 4.0 * g_zz_yyz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_yz[i] = -2.0 * g_0_yyz_0_yz[i] * b_exp + 4.0 * g_zz_yyz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_zz[i] = -2.0 * g_0_yyz_0_zz[i] * b_exp + 4.0 * g_zz_yyz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yzz_0_xx, g_0_yzz_0_xy, g_0_yzz_0_xz, g_0_yzz_0_yy, g_0_yzz_0_yz, g_0_yzz_0_zz, g_z_z_0_0_z_yz_0_xx, g_z_z_0_0_z_yz_0_xy, g_z_z_0_0_z_yz_0_xz, g_z_z_0_0_z_yz_0_yy, g_z_z_0_0_z_yz_0_yz, g_z_z_0_0_z_yz_0_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz, g_zz_yzz_0_xx, g_zz_yzz_0_xy, g_zz_yzz_0_xz, g_zz_yzz_0_yy, g_zz_yzz_0_yz, g_zz_yzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yz_0_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_yzz_0_xx[i] * b_exp - 2.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_yzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_yzz_0_xy[i] * b_exp - 2.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_yzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_yzz_0_xz[i] * b_exp - 2.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_yzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_yzz_0_yy[i] * b_exp - 2.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_yzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_yzz_0_yz[i] * b_exp - 2.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_yzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_yzz_0_zz[i] * b_exp - 2.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_yzz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_zzz_0_xx, g_0_zzz_0_xy, g_0_zzz_0_xz, g_0_zzz_0_yy, g_0_zzz_0_yz, g_0_zzz_0_zz, g_z_z_0_0_z_zz_0_xx, g_z_z_0_0_z_zz_0_xy, g_z_z_0_0_z_zz_0_xz, g_z_z_0_0_z_zz_0_yy, g_z_z_0_0_z_zz_0_yz, g_z_z_0_0_z_zz_0_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz, g_zz_zzz_0_xx, g_zz_zzz_0_xy, g_zz_zzz_0_xz, g_zz_zzz_0_yy, g_zz_zzz_0_yz, g_zz_zzz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_zz_0_xx[i] = 2.0 * g_0_z_0_xx[i] - 2.0 * g_0_zzz_0_xx[i] * b_exp - 4.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_zzz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_xy[i] = 2.0 * g_0_z_0_xy[i] - 2.0 * g_0_zzz_0_xy[i] * b_exp - 4.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_zzz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_xz[i] = 2.0 * g_0_z_0_xz[i] - 2.0 * g_0_zzz_0_xz[i] * b_exp - 4.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_zzz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_yy[i] = 2.0 * g_0_z_0_yy[i] - 2.0 * g_0_zzz_0_yy[i] * b_exp - 4.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_zzz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_yz[i] = 2.0 * g_0_z_0_yz[i] - 2.0 * g_0_zzz_0_yz[i] * b_exp - 4.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_zzz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_zz[i] = 2.0 * g_0_z_0_zz[i] - 2.0 * g_0_zzz_0_zz[i] * b_exp - 4.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_zzz_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


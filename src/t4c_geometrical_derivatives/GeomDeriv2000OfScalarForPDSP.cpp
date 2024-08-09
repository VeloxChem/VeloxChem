#include "GeomDeriv2000OfScalarForPDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pdsp_0(CSimdArray<double>& buffer_2000_pdsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_fdsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pdsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsp

    auto g_x_xx_0_x = buffer_pdsp[0];

    auto g_x_xx_0_y = buffer_pdsp[1];

    auto g_x_xx_0_z = buffer_pdsp[2];

    auto g_x_xy_0_x = buffer_pdsp[3];

    auto g_x_xy_0_y = buffer_pdsp[4];

    auto g_x_xy_0_z = buffer_pdsp[5];

    auto g_x_xz_0_x = buffer_pdsp[6];

    auto g_x_xz_0_y = buffer_pdsp[7];

    auto g_x_xz_0_z = buffer_pdsp[8];

    auto g_x_yy_0_x = buffer_pdsp[9];

    auto g_x_yy_0_y = buffer_pdsp[10];

    auto g_x_yy_0_z = buffer_pdsp[11];

    auto g_x_yz_0_x = buffer_pdsp[12];

    auto g_x_yz_0_y = buffer_pdsp[13];

    auto g_x_yz_0_z = buffer_pdsp[14];

    auto g_x_zz_0_x = buffer_pdsp[15];

    auto g_x_zz_0_y = buffer_pdsp[16];

    auto g_x_zz_0_z = buffer_pdsp[17];

    auto g_y_xx_0_x = buffer_pdsp[18];

    auto g_y_xx_0_y = buffer_pdsp[19];

    auto g_y_xx_0_z = buffer_pdsp[20];

    auto g_y_xy_0_x = buffer_pdsp[21];

    auto g_y_xy_0_y = buffer_pdsp[22];

    auto g_y_xy_0_z = buffer_pdsp[23];

    auto g_y_xz_0_x = buffer_pdsp[24];

    auto g_y_xz_0_y = buffer_pdsp[25];

    auto g_y_xz_0_z = buffer_pdsp[26];

    auto g_y_yy_0_x = buffer_pdsp[27];

    auto g_y_yy_0_y = buffer_pdsp[28];

    auto g_y_yy_0_z = buffer_pdsp[29];

    auto g_y_yz_0_x = buffer_pdsp[30];

    auto g_y_yz_0_y = buffer_pdsp[31];

    auto g_y_yz_0_z = buffer_pdsp[32];

    auto g_y_zz_0_x = buffer_pdsp[33];

    auto g_y_zz_0_y = buffer_pdsp[34];

    auto g_y_zz_0_z = buffer_pdsp[35];

    auto g_z_xx_0_x = buffer_pdsp[36];

    auto g_z_xx_0_y = buffer_pdsp[37];

    auto g_z_xx_0_z = buffer_pdsp[38];

    auto g_z_xy_0_x = buffer_pdsp[39];

    auto g_z_xy_0_y = buffer_pdsp[40];

    auto g_z_xy_0_z = buffer_pdsp[41];

    auto g_z_xz_0_x = buffer_pdsp[42];

    auto g_z_xz_0_y = buffer_pdsp[43];

    auto g_z_xz_0_z = buffer_pdsp[44];

    auto g_z_yy_0_x = buffer_pdsp[45];

    auto g_z_yy_0_y = buffer_pdsp[46];

    auto g_z_yy_0_z = buffer_pdsp[47];

    auto g_z_yz_0_x = buffer_pdsp[48];

    auto g_z_yz_0_y = buffer_pdsp[49];

    auto g_z_yz_0_z = buffer_pdsp[50];

    auto g_z_zz_0_x = buffer_pdsp[51];

    auto g_z_zz_0_y = buffer_pdsp[52];

    auto g_z_zz_0_z = buffer_pdsp[53];

    /// Set up components of auxilary buffer : buffer_fdsp

    auto g_xxx_xx_0_x = buffer_fdsp[0];

    auto g_xxx_xx_0_y = buffer_fdsp[1];

    auto g_xxx_xx_0_z = buffer_fdsp[2];

    auto g_xxx_xy_0_x = buffer_fdsp[3];

    auto g_xxx_xy_0_y = buffer_fdsp[4];

    auto g_xxx_xy_0_z = buffer_fdsp[5];

    auto g_xxx_xz_0_x = buffer_fdsp[6];

    auto g_xxx_xz_0_y = buffer_fdsp[7];

    auto g_xxx_xz_0_z = buffer_fdsp[8];

    auto g_xxx_yy_0_x = buffer_fdsp[9];

    auto g_xxx_yy_0_y = buffer_fdsp[10];

    auto g_xxx_yy_0_z = buffer_fdsp[11];

    auto g_xxx_yz_0_x = buffer_fdsp[12];

    auto g_xxx_yz_0_y = buffer_fdsp[13];

    auto g_xxx_yz_0_z = buffer_fdsp[14];

    auto g_xxx_zz_0_x = buffer_fdsp[15];

    auto g_xxx_zz_0_y = buffer_fdsp[16];

    auto g_xxx_zz_0_z = buffer_fdsp[17];

    auto g_xxy_xx_0_x = buffer_fdsp[18];

    auto g_xxy_xx_0_y = buffer_fdsp[19];

    auto g_xxy_xx_0_z = buffer_fdsp[20];

    auto g_xxy_xy_0_x = buffer_fdsp[21];

    auto g_xxy_xy_0_y = buffer_fdsp[22];

    auto g_xxy_xy_0_z = buffer_fdsp[23];

    auto g_xxy_xz_0_x = buffer_fdsp[24];

    auto g_xxy_xz_0_y = buffer_fdsp[25];

    auto g_xxy_xz_0_z = buffer_fdsp[26];

    auto g_xxy_yy_0_x = buffer_fdsp[27];

    auto g_xxy_yy_0_y = buffer_fdsp[28];

    auto g_xxy_yy_0_z = buffer_fdsp[29];

    auto g_xxy_yz_0_x = buffer_fdsp[30];

    auto g_xxy_yz_0_y = buffer_fdsp[31];

    auto g_xxy_yz_0_z = buffer_fdsp[32];

    auto g_xxy_zz_0_x = buffer_fdsp[33];

    auto g_xxy_zz_0_y = buffer_fdsp[34];

    auto g_xxy_zz_0_z = buffer_fdsp[35];

    auto g_xxz_xx_0_x = buffer_fdsp[36];

    auto g_xxz_xx_0_y = buffer_fdsp[37];

    auto g_xxz_xx_0_z = buffer_fdsp[38];

    auto g_xxz_xy_0_x = buffer_fdsp[39];

    auto g_xxz_xy_0_y = buffer_fdsp[40];

    auto g_xxz_xy_0_z = buffer_fdsp[41];

    auto g_xxz_xz_0_x = buffer_fdsp[42];

    auto g_xxz_xz_0_y = buffer_fdsp[43];

    auto g_xxz_xz_0_z = buffer_fdsp[44];

    auto g_xxz_yy_0_x = buffer_fdsp[45];

    auto g_xxz_yy_0_y = buffer_fdsp[46];

    auto g_xxz_yy_0_z = buffer_fdsp[47];

    auto g_xxz_yz_0_x = buffer_fdsp[48];

    auto g_xxz_yz_0_y = buffer_fdsp[49];

    auto g_xxz_yz_0_z = buffer_fdsp[50];

    auto g_xxz_zz_0_x = buffer_fdsp[51];

    auto g_xxz_zz_0_y = buffer_fdsp[52];

    auto g_xxz_zz_0_z = buffer_fdsp[53];

    auto g_xyy_xx_0_x = buffer_fdsp[54];

    auto g_xyy_xx_0_y = buffer_fdsp[55];

    auto g_xyy_xx_0_z = buffer_fdsp[56];

    auto g_xyy_xy_0_x = buffer_fdsp[57];

    auto g_xyy_xy_0_y = buffer_fdsp[58];

    auto g_xyy_xy_0_z = buffer_fdsp[59];

    auto g_xyy_xz_0_x = buffer_fdsp[60];

    auto g_xyy_xz_0_y = buffer_fdsp[61];

    auto g_xyy_xz_0_z = buffer_fdsp[62];

    auto g_xyy_yy_0_x = buffer_fdsp[63];

    auto g_xyy_yy_0_y = buffer_fdsp[64];

    auto g_xyy_yy_0_z = buffer_fdsp[65];

    auto g_xyy_yz_0_x = buffer_fdsp[66];

    auto g_xyy_yz_0_y = buffer_fdsp[67];

    auto g_xyy_yz_0_z = buffer_fdsp[68];

    auto g_xyy_zz_0_x = buffer_fdsp[69];

    auto g_xyy_zz_0_y = buffer_fdsp[70];

    auto g_xyy_zz_0_z = buffer_fdsp[71];

    auto g_xyz_xx_0_x = buffer_fdsp[72];

    auto g_xyz_xx_0_y = buffer_fdsp[73];

    auto g_xyz_xx_0_z = buffer_fdsp[74];

    auto g_xyz_xy_0_x = buffer_fdsp[75];

    auto g_xyz_xy_0_y = buffer_fdsp[76];

    auto g_xyz_xy_0_z = buffer_fdsp[77];

    auto g_xyz_xz_0_x = buffer_fdsp[78];

    auto g_xyz_xz_0_y = buffer_fdsp[79];

    auto g_xyz_xz_0_z = buffer_fdsp[80];

    auto g_xyz_yy_0_x = buffer_fdsp[81];

    auto g_xyz_yy_0_y = buffer_fdsp[82];

    auto g_xyz_yy_0_z = buffer_fdsp[83];

    auto g_xyz_yz_0_x = buffer_fdsp[84];

    auto g_xyz_yz_0_y = buffer_fdsp[85];

    auto g_xyz_yz_0_z = buffer_fdsp[86];

    auto g_xyz_zz_0_x = buffer_fdsp[87];

    auto g_xyz_zz_0_y = buffer_fdsp[88];

    auto g_xyz_zz_0_z = buffer_fdsp[89];

    auto g_xzz_xx_0_x = buffer_fdsp[90];

    auto g_xzz_xx_0_y = buffer_fdsp[91];

    auto g_xzz_xx_0_z = buffer_fdsp[92];

    auto g_xzz_xy_0_x = buffer_fdsp[93];

    auto g_xzz_xy_0_y = buffer_fdsp[94];

    auto g_xzz_xy_0_z = buffer_fdsp[95];

    auto g_xzz_xz_0_x = buffer_fdsp[96];

    auto g_xzz_xz_0_y = buffer_fdsp[97];

    auto g_xzz_xz_0_z = buffer_fdsp[98];

    auto g_xzz_yy_0_x = buffer_fdsp[99];

    auto g_xzz_yy_0_y = buffer_fdsp[100];

    auto g_xzz_yy_0_z = buffer_fdsp[101];

    auto g_xzz_yz_0_x = buffer_fdsp[102];

    auto g_xzz_yz_0_y = buffer_fdsp[103];

    auto g_xzz_yz_0_z = buffer_fdsp[104];

    auto g_xzz_zz_0_x = buffer_fdsp[105];

    auto g_xzz_zz_0_y = buffer_fdsp[106];

    auto g_xzz_zz_0_z = buffer_fdsp[107];

    auto g_yyy_xx_0_x = buffer_fdsp[108];

    auto g_yyy_xx_0_y = buffer_fdsp[109];

    auto g_yyy_xx_0_z = buffer_fdsp[110];

    auto g_yyy_xy_0_x = buffer_fdsp[111];

    auto g_yyy_xy_0_y = buffer_fdsp[112];

    auto g_yyy_xy_0_z = buffer_fdsp[113];

    auto g_yyy_xz_0_x = buffer_fdsp[114];

    auto g_yyy_xz_0_y = buffer_fdsp[115];

    auto g_yyy_xz_0_z = buffer_fdsp[116];

    auto g_yyy_yy_0_x = buffer_fdsp[117];

    auto g_yyy_yy_0_y = buffer_fdsp[118];

    auto g_yyy_yy_0_z = buffer_fdsp[119];

    auto g_yyy_yz_0_x = buffer_fdsp[120];

    auto g_yyy_yz_0_y = buffer_fdsp[121];

    auto g_yyy_yz_0_z = buffer_fdsp[122];

    auto g_yyy_zz_0_x = buffer_fdsp[123];

    auto g_yyy_zz_0_y = buffer_fdsp[124];

    auto g_yyy_zz_0_z = buffer_fdsp[125];

    auto g_yyz_xx_0_x = buffer_fdsp[126];

    auto g_yyz_xx_0_y = buffer_fdsp[127];

    auto g_yyz_xx_0_z = buffer_fdsp[128];

    auto g_yyz_xy_0_x = buffer_fdsp[129];

    auto g_yyz_xy_0_y = buffer_fdsp[130];

    auto g_yyz_xy_0_z = buffer_fdsp[131];

    auto g_yyz_xz_0_x = buffer_fdsp[132];

    auto g_yyz_xz_0_y = buffer_fdsp[133];

    auto g_yyz_xz_0_z = buffer_fdsp[134];

    auto g_yyz_yy_0_x = buffer_fdsp[135];

    auto g_yyz_yy_0_y = buffer_fdsp[136];

    auto g_yyz_yy_0_z = buffer_fdsp[137];

    auto g_yyz_yz_0_x = buffer_fdsp[138];

    auto g_yyz_yz_0_y = buffer_fdsp[139];

    auto g_yyz_yz_0_z = buffer_fdsp[140];

    auto g_yyz_zz_0_x = buffer_fdsp[141];

    auto g_yyz_zz_0_y = buffer_fdsp[142];

    auto g_yyz_zz_0_z = buffer_fdsp[143];

    auto g_yzz_xx_0_x = buffer_fdsp[144];

    auto g_yzz_xx_0_y = buffer_fdsp[145];

    auto g_yzz_xx_0_z = buffer_fdsp[146];

    auto g_yzz_xy_0_x = buffer_fdsp[147];

    auto g_yzz_xy_0_y = buffer_fdsp[148];

    auto g_yzz_xy_0_z = buffer_fdsp[149];

    auto g_yzz_xz_0_x = buffer_fdsp[150];

    auto g_yzz_xz_0_y = buffer_fdsp[151];

    auto g_yzz_xz_0_z = buffer_fdsp[152];

    auto g_yzz_yy_0_x = buffer_fdsp[153];

    auto g_yzz_yy_0_y = buffer_fdsp[154];

    auto g_yzz_yy_0_z = buffer_fdsp[155];

    auto g_yzz_yz_0_x = buffer_fdsp[156];

    auto g_yzz_yz_0_y = buffer_fdsp[157];

    auto g_yzz_yz_0_z = buffer_fdsp[158];

    auto g_yzz_zz_0_x = buffer_fdsp[159];

    auto g_yzz_zz_0_y = buffer_fdsp[160];

    auto g_yzz_zz_0_z = buffer_fdsp[161];

    auto g_zzz_xx_0_x = buffer_fdsp[162];

    auto g_zzz_xx_0_y = buffer_fdsp[163];

    auto g_zzz_xx_0_z = buffer_fdsp[164];

    auto g_zzz_xy_0_x = buffer_fdsp[165];

    auto g_zzz_xy_0_y = buffer_fdsp[166];

    auto g_zzz_xy_0_z = buffer_fdsp[167];

    auto g_zzz_xz_0_x = buffer_fdsp[168];

    auto g_zzz_xz_0_y = buffer_fdsp[169];

    auto g_zzz_xz_0_z = buffer_fdsp[170];

    auto g_zzz_yy_0_x = buffer_fdsp[171];

    auto g_zzz_yy_0_y = buffer_fdsp[172];

    auto g_zzz_yy_0_z = buffer_fdsp[173];

    auto g_zzz_yz_0_x = buffer_fdsp[174];

    auto g_zzz_yz_0_y = buffer_fdsp[175];

    auto g_zzz_yz_0_z = buffer_fdsp[176];

    auto g_zzz_zz_0_x = buffer_fdsp[177];

    auto g_zzz_zz_0_y = buffer_fdsp[178];

    auto g_zzz_zz_0_z = buffer_fdsp[179];

    /// Set up components of integrals buffer : buffer_2000_pdsp

    auto g_xx_0_0_0_x_xx_0_x = buffer_2000_pdsp[0];

    auto g_xx_0_0_0_x_xx_0_y = buffer_2000_pdsp[1];

    auto g_xx_0_0_0_x_xx_0_z = buffer_2000_pdsp[2];

    auto g_xx_0_0_0_x_xy_0_x = buffer_2000_pdsp[3];

    auto g_xx_0_0_0_x_xy_0_y = buffer_2000_pdsp[4];

    auto g_xx_0_0_0_x_xy_0_z = buffer_2000_pdsp[5];

    auto g_xx_0_0_0_x_xz_0_x = buffer_2000_pdsp[6];

    auto g_xx_0_0_0_x_xz_0_y = buffer_2000_pdsp[7];

    auto g_xx_0_0_0_x_xz_0_z = buffer_2000_pdsp[8];

    auto g_xx_0_0_0_x_yy_0_x = buffer_2000_pdsp[9];

    auto g_xx_0_0_0_x_yy_0_y = buffer_2000_pdsp[10];

    auto g_xx_0_0_0_x_yy_0_z = buffer_2000_pdsp[11];

    auto g_xx_0_0_0_x_yz_0_x = buffer_2000_pdsp[12];

    auto g_xx_0_0_0_x_yz_0_y = buffer_2000_pdsp[13];

    auto g_xx_0_0_0_x_yz_0_z = buffer_2000_pdsp[14];

    auto g_xx_0_0_0_x_zz_0_x = buffer_2000_pdsp[15];

    auto g_xx_0_0_0_x_zz_0_y = buffer_2000_pdsp[16];

    auto g_xx_0_0_0_x_zz_0_z = buffer_2000_pdsp[17];

    auto g_xx_0_0_0_y_xx_0_x = buffer_2000_pdsp[18];

    auto g_xx_0_0_0_y_xx_0_y = buffer_2000_pdsp[19];

    auto g_xx_0_0_0_y_xx_0_z = buffer_2000_pdsp[20];

    auto g_xx_0_0_0_y_xy_0_x = buffer_2000_pdsp[21];

    auto g_xx_0_0_0_y_xy_0_y = buffer_2000_pdsp[22];

    auto g_xx_0_0_0_y_xy_0_z = buffer_2000_pdsp[23];

    auto g_xx_0_0_0_y_xz_0_x = buffer_2000_pdsp[24];

    auto g_xx_0_0_0_y_xz_0_y = buffer_2000_pdsp[25];

    auto g_xx_0_0_0_y_xz_0_z = buffer_2000_pdsp[26];

    auto g_xx_0_0_0_y_yy_0_x = buffer_2000_pdsp[27];

    auto g_xx_0_0_0_y_yy_0_y = buffer_2000_pdsp[28];

    auto g_xx_0_0_0_y_yy_0_z = buffer_2000_pdsp[29];

    auto g_xx_0_0_0_y_yz_0_x = buffer_2000_pdsp[30];

    auto g_xx_0_0_0_y_yz_0_y = buffer_2000_pdsp[31];

    auto g_xx_0_0_0_y_yz_0_z = buffer_2000_pdsp[32];

    auto g_xx_0_0_0_y_zz_0_x = buffer_2000_pdsp[33];

    auto g_xx_0_0_0_y_zz_0_y = buffer_2000_pdsp[34];

    auto g_xx_0_0_0_y_zz_0_z = buffer_2000_pdsp[35];

    auto g_xx_0_0_0_z_xx_0_x = buffer_2000_pdsp[36];

    auto g_xx_0_0_0_z_xx_0_y = buffer_2000_pdsp[37];

    auto g_xx_0_0_0_z_xx_0_z = buffer_2000_pdsp[38];

    auto g_xx_0_0_0_z_xy_0_x = buffer_2000_pdsp[39];

    auto g_xx_0_0_0_z_xy_0_y = buffer_2000_pdsp[40];

    auto g_xx_0_0_0_z_xy_0_z = buffer_2000_pdsp[41];

    auto g_xx_0_0_0_z_xz_0_x = buffer_2000_pdsp[42];

    auto g_xx_0_0_0_z_xz_0_y = buffer_2000_pdsp[43];

    auto g_xx_0_0_0_z_xz_0_z = buffer_2000_pdsp[44];

    auto g_xx_0_0_0_z_yy_0_x = buffer_2000_pdsp[45];

    auto g_xx_0_0_0_z_yy_0_y = buffer_2000_pdsp[46];

    auto g_xx_0_0_0_z_yy_0_z = buffer_2000_pdsp[47];

    auto g_xx_0_0_0_z_yz_0_x = buffer_2000_pdsp[48];

    auto g_xx_0_0_0_z_yz_0_y = buffer_2000_pdsp[49];

    auto g_xx_0_0_0_z_yz_0_z = buffer_2000_pdsp[50];

    auto g_xx_0_0_0_z_zz_0_x = buffer_2000_pdsp[51];

    auto g_xx_0_0_0_z_zz_0_y = buffer_2000_pdsp[52];

    auto g_xx_0_0_0_z_zz_0_z = buffer_2000_pdsp[53];

    auto g_xy_0_0_0_x_xx_0_x = buffer_2000_pdsp[54];

    auto g_xy_0_0_0_x_xx_0_y = buffer_2000_pdsp[55];

    auto g_xy_0_0_0_x_xx_0_z = buffer_2000_pdsp[56];

    auto g_xy_0_0_0_x_xy_0_x = buffer_2000_pdsp[57];

    auto g_xy_0_0_0_x_xy_0_y = buffer_2000_pdsp[58];

    auto g_xy_0_0_0_x_xy_0_z = buffer_2000_pdsp[59];

    auto g_xy_0_0_0_x_xz_0_x = buffer_2000_pdsp[60];

    auto g_xy_0_0_0_x_xz_0_y = buffer_2000_pdsp[61];

    auto g_xy_0_0_0_x_xz_0_z = buffer_2000_pdsp[62];

    auto g_xy_0_0_0_x_yy_0_x = buffer_2000_pdsp[63];

    auto g_xy_0_0_0_x_yy_0_y = buffer_2000_pdsp[64];

    auto g_xy_0_0_0_x_yy_0_z = buffer_2000_pdsp[65];

    auto g_xy_0_0_0_x_yz_0_x = buffer_2000_pdsp[66];

    auto g_xy_0_0_0_x_yz_0_y = buffer_2000_pdsp[67];

    auto g_xy_0_0_0_x_yz_0_z = buffer_2000_pdsp[68];

    auto g_xy_0_0_0_x_zz_0_x = buffer_2000_pdsp[69];

    auto g_xy_0_0_0_x_zz_0_y = buffer_2000_pdsp[70];

    auto g_xy_0_0_0_x_zz_0_z = buffer_2000_pdsp[71];

    auto g_xy_0_0_0_y_xx_0_x = buffer_2000_pdsp[72];

    auto g_xy_0_0_0_y_xx_0_y = buffer_2000_pdsp[73];

    auto g_xy_0_0_0_y_xx_0_z = buffer_2000_pdsp[74];

    auto g_xy_0_0_0_y_xy_0_x = buffer_2000_pdsp[75];

    auto g_xy_0_0_0_y_xy_0_y = buffer_2000_pdsp[76];

    auto g_xy_0_0_0_y_xy_0_z = buffer_2000_pdsp[77];

    auto g_xy_0_0_0_y_xz_0_x = buffer_2000_pdsp[78];

    auto g_xy_0_0_0_y_xz_0_y = buffer_2000_pdsp[79];

    auto g_xy_0_0_0_y_xz_0_z = buffer_2000_pdsp[80];

    auto g_xy_0_0_0_y_yy_0_x = buffer_2000_pdsp[81];

    auto g_xy_0_0_0_y_yy_0_y = buffer_2000_pdsp[82];

    auto g_xy_0_0_0_y_yy_0_z = buffer_2000_pdsp[83];

    auto g_xy_0_0_0_y_yz_0_x = buffer_2000_pdsp[84];

    auto g_xy_0_0_0_y_yz_0_y = buffer_2000_pdsp[85];

    auto g_xy_0_0_0_y_yz_0_z = buffer_2000_pdsp[86];

    auto g_xy_0_0_0_y_zz_0_x = buffer_2000_pdsp[87];

    auto g_xy_0_0_0_y_zz_0_y = buffer_2000_pdsp[88];

    auto g_xy_0_0_0_y_zz_0_z = buffer_2000_pdsp[89];

    auto g_xy_0_0_0_z_xx_0_x = buffer_2000_pdsp[90];

    auto g_xy_0_0_0_z_xx_0_y = buffer_2000_pdsp[91];

    auto g_xy_0_0_0_z_xx_0_z = buffer_2000_pdsp[92];

    auto g_xy_0_0_0_z_xy_0_x = buffer_2000_pdsp[93];

    auto g_xy_0_0_0_z_xy_0_y = buffer_2000_pdsp[94];

    auto g_xy_0_0_0_z_xy_0_z = buffer_2000_pdsp[95];

    auto g_xy_0_0_0_z_xz_0_x = buffer_2000_pdsp[96];

    auto g_xy_0_0_0_z_xz_0_y = buffer_2000_pdsp[97];

    auto g_xy_0_0_0_z_xz_0_z = buffer_2000_pdsp[98];

    auto g_xy_0_0_0_z_yy_0_x = buffer_2000_pdsp[99];

    auto g_xy_0_0_0_z_yy_0_y = buffer_2000_pdsp[100];

    auto g_xy_0_0_0_z_yy_0_z = buffer_2000_pdsp[101];

    auto g_xy_0_0_0_z_yz_0_x = buffer_2000_pdsp[102];

    auto g_xy_0_0_0_z_yz_0_y = buffer_2000_pdsp[103];

    auto g_xy_0_0_0_z_yz_0_z = buffer_2000_pdsp[104];

    auto g_xy_0_0_0_z_zz_0_x = buffer_2000_pdsp[105];

    auto g_xy_0_0_0_z_zz_0_y = buffer_2000_pdsp[106];

    auto g_xy_0_0_0_z_zz_0_z = buffer_2000_pdsp[107];

    auto g_xz_0_0_0_x_xx_0_x = buffer_2000_pdsp[108];

    auto g_xz_0_0_0_x_xx_0_y = buffer_2000_pdsp[109];

    auto g_xz_0_0_0_x_xx_0_z = buffer_2000_pdsp[110];

    auto g_xz_0_0_0_x_xy_0_x = buffer_2000_pdsp[111];

    auto g_xz_0_0_0_x_xy_0_y = buffer_2000_pdsp[112];

    auto g_xz_0_0_0_x_xy_0_z = buffer_2000_pdsp[113];

    auto g_xz_0_0_0_x_xz_0_x = buffer_2000_pdsp[114];

    auto g_xz_0_0_0_x_xz_0_y = buffer_2000_pdsp[115];

    auto g_xz_0_0_0_x_xz_0_z = buffer_2000_pdsp[116];

    auto g_xz_0_0_0_x_yy_0_x = buffer_2000_pdsp[117];

    auto g_xz_0_0_0_x_yy_0_y = buffer_2000_pdsp[118];

    auto g_xz_0_0_0_x_yy_0_z = buffer_2000_pdsp[119];

    auto g_xz_0_0_0_x_yz_0_x = buffer_2000_pdsp[120];

    auto g_xz_0_0_0_x_yz_0_y = buffer_2000_pdsp[121];

    auto g_xz_0_0_0_x_yz_0_z = buffer_2000_pdsp[122];

    auto g_xz_0_0_0_x_zz_0_x = buffer_2000_pdsp[123];

    auto g_xz_0_0_0_x_zz_0_y = buffer_2000_pdsp[124];

    auto g_xz_0_0_0_x_zz_0_z = buffer_2000_pdsp[125];

    auto g_xz_0_0_0_y_xx_0_x = buffer_2000_pdsp[126];

    auto g_xz_0_0_0_y_xx_0_y = buffer_2000_pdsp[127];

    auto g_xz_0_0_0_y_xx_0_z = buffer_2000_pdsp[128];

    auto g_xz_0_0_0_y_xy_0_x = buffer_2000_pdsp[129];

    auto g_xz_0_0_0_y_xy_0_y = buffer_2000_pdsp[130];

    auto g_xz_0_0_0_y_xy_0_z = buffer_2000_pdsp[131];

    auto g_xz_0_0_0_y_xz_0_x = buffer_2000_pdsp[132];

    auto g_xz_0_0_0_y_xz_0_y = buffer_2000_pdsp[133];

    auto g_xz_0_0_0_y_xz_0_z = buffer_2000_pdsp[134];

    auto g_xz_0_0_0_y_yy_0_x = buffer_2000_pdsp[135];

    auto g_xz_0_0_0_y_yy_0_y = buffer_2000_pdsp[136];

    auto g_xz_0_0_0_y_yy_0_z = buffer_2000_pdsp[137];

    auto g_xz_0_0_0_y_yz_0_x = buffer_2000_pdsp[138];

    auto g_xz_0_0_0_y_yz_0_y = buffer_2000_pdsp[139];

    auto g_xz_0_0_0_y_yz_0_z = buffer_2000_pdsp[140];

    auto g_xz_0_0_0_y_zz_0_x = buffer_2000_pdsp[141];

    auto g_xz_0_0_0_y_zz_0_y = buffer_2000_pdsp[142];

    auto g_xz_0_0_0_y_zz_0_z = buffer_2000_pdsp[143];

    auto g_xz_0_0_0_z_xx_0_x = buffer_2000_pdsp[144];

    auto g_xz_0_0_0_z_xx_0_y = buffer_2000_pdsp[145];

    auto g_xz_0_0_0_z_xx_0_z = buffer_2000_pdsp[146];

    auto g_xz_0_0_0_z_xy_0_x = buffer_2000_pdsp[147];

    auto g_xz_0_0_0_z_xy_0_y = buffer_2000_pdsp[148];

    auto g_xz_0_0_0_z_xy_0_z = buffer_2000_pdsp[149];

    auto g_xz_0_0_0_z_xz_0_x = buffer_2000_pdsp[150];

    auto g_xz_0_0_0_z_xz_0_y = buffer_2000_pdsp[151];

    auto g_xz_0_0_0_z_xz_0_z = buffer_2000_pdsp[152];

    auto g_xz_0_0_0_z_yy_0_x = buffer_2000_pdsp[153];

    auto g_xz_0_0_0_z_yy_0_y = buffer_2000_pdsp[154];

    auto g_xz_0_0_0_z_yy_0_z = buffer_2000_pdsp[155];

    auto g_xz_0_0_0_z_yz_0_x = buffer_2000_pdsp[156];

    auto g_xz_0_0_0_z_yz_0_y = buffer_2000_pdsp[157];

    auto g_xz_0_0_0_z_yz_0_z = buffer_2000_pdsp[158];

    auto g_xz_0_0_0_z_zz_0_x = buffer_2000_pdsp[159];

    auto g_xz_0_0_0_z_zz_0_y = buffer_2000_pdsp[160];

    auto g_xz_0_0_0_z_zz_0_z = buffer_2000_pdsp[161];

    auto g_yy_0_0_0_x_xx_0_x = buffer_2000_pdsp[162];

    auto g_yy_0_0_0_x_xx_0_y = buffer_2000_pdsp[163];

    auto g_yy_0_0_0_x_xx_0_z = buffer_2000_pdsp[164];

    auto g_yy_0_0_0_x_xy_0_x = buffer_2000_pdsp[165];

    auto g_yy_0_0_0_x_xy_0_y = buffer_2000_pdsp[166];

    auto g_yy_0_0_0_x_xy_0_z = buffer_2000_pdsp[167];

    auto g_yy_0_0_0_x_xz_0_x = buffer_2000_pdsp[168];

    auto g_yy_0_0_0_x_xz_0_y = buffer_2000_pdsp[169];

    auto g_yy_0_0_0_x_xz_0_z = buffer_2000_pdsp[170];

    auto g_yy_0_0_0_x_yy_0_x = buffer_2000_pdsp[171];

    auto g_yy_0_0_0_x_yy_0_y = buffer_2000_pdsp[172];

    auto g_yy_0_0_0_x_yy_0_z = buffer_2000_pdsp[173];

    auto g_yy_0_0_0_x_yz_0_x = buffer_2000_pdsp[174];

    auto g_yy_0_0_0_x_yz_0_y = buffer_2000_pdsp[175];

    auto g_yy_0_0_0_x_yz_0_z = buffer_2000_pdsp[176];

    auto g_yy_0_0_0_x_zz_0_x = buffer_2000_pdsp[177];

    auto g_yy_0_0_0_x_zz_0_y = buffer_2000_pdsp[178];

    auto g_yy_0_0_0_x_zz_0_z = buffer_2000_pdsp[179];

    auto g_yy_0_0_0_y_xx_0_x = buffer_2000_pdsp[180];

    auto g_yy_0_0_0_y_xx_0_y = buffer_2000_pdsp[181];

    auto g_yy_0_0_0_y_xx_0_z = buffer_2000_pdsp[182];

    auto g_yy_0_0_0_y_xy_0_x = buffer_2000_pdsp[183];

    auto g_yy_0_0_0_y_xy_0_y = buffer_2000_pdsp[184];

    auto g_yy_0_0_0_y_xy_0_z = buffer_2000_pdsp[185];

    auto g_yy_0_0_0_y_xz_0_x = buffer_2000_pdsp[186];

    auto g_yy_0_0_0_y_xz_0_y = buffer_2000_pdsp[187];

    auto g_yy_0_0_0_y_xz_0_z = buffer_2000_pdsp[188];

    auto g_yy_0_0_0_y_yy_0_x = buffer_2000_pdsp[189];

    auto g_yy_0_0_0_y_yy_0_y = buffer_2000_pdsp[190];

    auto g_yy_0_0_0_y_yy_0_z = buffer_2000_pdsp[191];

    auto g_yy_0_0_0_y_yz_0_x = buffer_2000_pdsp[192];

    auto g_yy_0_0_0_y_yz_0_y = buffer_2000_pdsp[193];

    auto g_yy_0_0_0_y_yz_0_z = buffer_2000_pdsp[194];

    auto g_yy_0_0_0_y_zz_0_x = buffer_2000_pdsp[195];

    auto g_yy_0_0_0_y_zz_0_y = buffer_2000_pdsp[196];

    auto g_yy_0_0_0_y_zz_0_z = buffer_2000_pdsp[197];

    auto g_yy_0_0_0_z_xx_0_x = buffer_2000_pdsp[198];

    auto g_yy_0_0_0_z_xx_0_y = buffer_2000_pdsp[199];

    auto g_yy_0_0_0_z_xx_0_z = buffer_2000_pdsp[200];

    auto g_yy_0_0_0_z_xy_0_x = buffer_2000_pdsp[201];

    auto g_yy_0_0_0_z_xy_0_y = buffer_2000_pdsp[202];

    auto g_yy_0_0_0_z_xy_0_z = buffer_2000_pdsp[203];

    auto g_yy_0_0_0_z_xz_0_x = buffer_2000_pdsp[204];

    auto g_yy_0_0_0_z_xz_0_y = buffer_2000_pdsp[205];

    auto g_yy_0_0_0_z_xz_0_z = buffer_2000_pdsp[206];

    auto g_yy_0_0_0_z_yy_0_x = buffer_2000_pdsp[207];

    auto g_yy_0_0_0_z_yy_0_y = buffer_2000_pdsp[208];

    auto g_yy_0_0_0_z_yy_0_z = buffer_2000_pdsp[209];

    auto g_yy_0_0_0_z_yz_0_x = buffer_2000_pdsp[210];

    auto g_yy_0_0_0_z_yz_0_y = buffer_2000_pdsp[211];

    auto g_yy_0_0_0_z_yz_0_z = buffer_2000_pdsp[212];

    auto g_yy_0_0_0_z_zz_0_x = buffer_2000_pdsp[213];

    auto g_yy_0_0_0_z_zz_0_y = buffer_2000_pdsp[214];

    auto g_yy_0_0_0_z_zz_0_z = buffer_2000_pdsp[215];

    auto g_yz_0_0_0_x_xx_0_x = buffer_2000_pdsp[216];

    auto g_yz_0_0_0_x_xx_0_y = buffer_2000_pdsp[217];

    auto g_yz_0_0_0_x_xx_0_z = buffer_2000_pdsp[218];

    auto g_yz_0_0_0_x_xy_0_x = buffer_2000_pdsp[219];

    auto g_yz_0_0_0_x_xy_0_y = buffer_2000_pdsp[220];

    auto g_yz_0_0_0_x_xy_0_z = buffer_2000_pdsp[221];

    auto g_yz_0_0_0_x_xz_0_x = buffer_2000_pdsp[222];

    auto g_yz_0_0_0_x_xz_0_y = buffer_2000_pdsp[223];

    auto g_yz_0_0_0_x_xz_0_z = buffer_2000_pdsp[224];

    auto g_yz_0_0_0_x_yy_0_x = buffer_2000_pdsp[225];

    auto g_yz_0_0_0_x_yy_0_y = buffer_2000_pdsp[226];

    auto g_yz_0_0_0_x_yy_0_z = buffer_2000_pdsp[227];

    auto g_yz_0_0_0_x_yz_0_x = buffer_2000_pdsp[228];

    auto g_yz_0_0_0_x_yz_0_y = buffer_2000_pdsp[229];

    auto g_yz_0_0_0_x_yz_0_z = buffer_2000_pdsp[230];

    auto g_yz_0_0_0_x_zz_0_x = buffer_2000_pdsp[231];

    auto g_yz_0_0_0_x_zz_0_y = buffer_2000_pdsp[232];

    auto g_yz_0_0_0_x_zz_0_z = buffer_2000_pdsp[233];

    auto g_yz_0_0_0_y_xx_0_x = buffer_2000_pdsp[234];

    auto g_yz_0_0_0_y_xx_0_y = buffer_2000_pdsp[235];

    auto g_yz_0_0_0_y_xx_0_z = buffer_2000_pdsp[236];

    auto g_yz_0_0_0_y_xy_0_x = buffer_2000_pdsp[237];

    auto g_yz_0_0_0_y_xy_0_y = buffer_2000_pdsp[238];

    auto g_yz_0_0_0_y_xy_0_z = buffer_2000_pdsp[239];

    auto g_yz_0_0_0_y_xz_0_x = buffer_2000_pdsp[240];

    auto g_yz_0_0_0_y_xz_0_y = buffer_2000_pdsp[241];

    auto g_yz_0_0_0_y_xz_0_z = buffer_2000_pdsp[242];

    auto g_yz_0_0_0_y_yy_0_x = buffer_2000_pdsp[243];

    auto g_yz_0_0_0_y_yy_0_y = buffer_2000_pdsp[244];

    auto g_yz_0_0_0_y_yy_0_z = buffer_2000_pdsp[245];

    auto g_yz_0_0_0_y_yz_0_x = buffer_2000_pdsp[246];

    auto g_yz_0_0_0_y_yz_0_y = buffer_2000_pdsp[247];

    auto g_yz_0_0_0_y_yz_0_z = buffer_2000_pdsp[248];

    auto g_yz_0_0_0_y_zz_0_x = buffer_2000_pdsp[249];

    auto g_yz_0_0_0_y_zz_0_y = buffer_2000_pdsp[250];

    auto g_yz_0_0_0_y_zz_0_z = buffer_2000_pdsp[251];

    auto g_yz_0_0_0_z_xx_0_x = buffer_2000_pdsp[252];

    auto g_yz_0_0_0_z_xx_0_y = buffer_2000_pdsp[253];

    auto g_yz_0_0_0_z_xx_0_z = buffer_2000_pdsp[254];

    auto g_yz_0_0_0_z_xy_0_x = buffer_2000_pdsp[255];

    auto g_yz_0_0_0_z_xy_0_y = buffer_2000_pdsp[256];

    auto g_yz_0_0_0_z_xy_0_z = buffer_2000_pdsp[257];

    auto g_yz_0_0_0_z_xz_0_x = buffer_2000_pdsp[258];

    auto g_yz_0_0_0_z_xz_0_y = buffer_2000_pdsp[259];

    auto g_yz_0_0_0_z_xz_0_z = buffer_2000_pdsp[260];

    auto g_yz_0_0_0_z_yy_0_x = buffer_2000_pdsp[261];

    auto g_yz_0_0_0_z_yy_0_y = buffer_2000_pdsp[262];

    auto g_yz_0_0_0_z_yy_0_z = buffer_2000_pdsp[263];

    auto g_yz_0_0_0_z_yz_0_x = buffer_2000_pdsp[264];

    auto g_yz_0_0_0_z_yz_0_y = buffer_2000_pdsp[265];

    auto g_yz_0_0_0_z_yz_0_z = buffer_2000_pdsp[266];

    auto g_yz_0_0_0_z_zz_0_x = buffer_2000_pdsp[267];

    auto g_yz_0_0_0_z_zz_0_y = buffer_2000_pdsp[268];

    auto g_yz_0_0_0_z_zz_0_z = buffer_2000_pdsp[269];

    auto g_zz_0_0_0_x_xx_0_x = buffer_2000_pdsp[270];

    auto g_zz_0_0_0_x_xx_0_y = buffer_2000_pdsp[271];

    auto g_zz_0_0_0_x_xx_0_z = buffer_2000_pdsp[272];

    auto g_zz_0_0_0_x_xy_0_x = buffer_2000_pdsp[273];

    auto g_zz_0_0_0_x_xy_0_y = buffer_2000_pdsp[274];

    auto g_zz_0_0_0_x_xy_0_z = buffer_2000_pdsp[275];

    auto g_zz_0_0_0_x_xz_0_x = buffer_2000_pdsp[276];

    auto g_zz_0_0_0_x_xz_0_y = buffer_2000_pdsp[277];

    auto g_zz_0_0_0_x_xz_0_z = buffer_2000_pdsp[278];

    auto g_zz_0_0_0_x_yy_0_x = buffer_2000_pdsp[279];

    auto g_zz_0_0_0_x_yy_0_y = buffer_2000_pdsp[280];

    auto g_zz_0_0_0_x_yy_0_z = buffer_2000_pdsp[281];

    auto g_zz_0_0_0_x_yz_0_x = buffer_2000_pdsp[282];

    auto g_zz_0_0_0_x_yz_0_y = buffer_2000_pdsp[283];

    auto g_zz_0_0_0_x_yz_0_z = buffer_2000_pdsp[284];

    auto g_zz_0_0_0_x_zz_0_x = buffer_2000_pdsp[285];

    auto g_zz_0_0_0_x_zz_0_y = buffer_2000_pdsp[286];

    auto g_zz_0_0_0_x_zz_0_z = buffer_2000_pdsp[287];

    auto g_zz_0_0_0_y_xx_0_x = buffer_2000_pdsp[288];

    auto g_zz_0_0_0_y_xx_0_y = buffer_2000_pdsp[289];

    auto g_zz_0_0_0_y_xx_0_z = buffer_2000_pdsp[290];

    auto g_zz_0_0_0_y_xy_0_x = buffer_2000_pdsp[291];

    auto g_zz_0_0_0_y_xy_0_y = buffer_2000_pdsp[292];

    auto g_zz_0_0_0_y_xy_0_z = buffer_2000_pdsp[293];

    auto g_zz_0_0_0_y_xz_0_x = buffer_2000_pdsp[294];

    auto g_zz_0_0_0_y_xz_0_y = buffer_2000_pdsp[295];

    auto g_zz_0_0_0_y_xz_0_z = buffer_2000_pdsp[296];

    auto g_zz_0_0_0_y_yy_0_x = buffer_2000_pdsp[297];

    auto g_zz_0_0_0_y_yy_0_y = buffer_2000_pdsp[298];

    auto g_zz_0_0_0_y_yy_0_z = buffer_2000_pdsp[299];

    auto g_zz_0_0_0_y_yz_0_x = buffer_2000_pdsp[300];

    auto g_zz_0_0_0_y_yz_0_y = buffer_2000_pdsp[301];

    auto g_zz_0_0_0_y_yz_0_z = buffer_2000_pdsp[302];

    auto g_zz_0_0_0_y_zz_0_x = buffer_2000_pdsp[303];

    auto g_zz_0_0_0_y_zz_0_y = buffer_2000_pdsp[304];

    auto g_zz_0_0_0_y_zz_0_z = buffer_2000_pdsp[305];

    auto g_zz_0_0_0_z_xx_0_x = buffer_2000_pdsp[306];

    auto g_zz_0_0_0_z_xx_0_y = buffer_2000_pdsp[307];

    auto g_zz_0_0_0_z_xx_0_z = buffer_2000_pdsp[308];

    auto g_zz_0_0_0_z_xy_0_x = buffer_2000_pdsp[309];

    auto g_zz_0_0_0_z_xy_0_y = buffer_2000_pdsp[310];

    auto g_zz_0_0_0_z_xy_0_z = buffer_2000_pdsp[311];

    auto g_zz_0_0_0_z_xz_0_x = buffer_2000_pdsp[312];

    auto g_zz_0_0_0_z_xz_0_y = buffer_2000_pdsp[313];

    auto g_zz_0_0_0_z_xz_0_z = buffer_2000_pdsp[314];

    auto g_zz_0_0_0_z_yy_0_x = buffer_2000_pdsp[315];

    auto g_zz_0_0_0_z_yy_0_y = buffer_2000_pdsp[316];

    auto g_zz_0_0_0_z_yy_0_z = buffer_2000_pdsp[317];

    auto g_zz_0_0_0_z_yz_0_x = buffer_2000_pdsp[318];

    auto g_zz_0_0_0_z_yz_0_y = buffer_2000_pdsp[319];

    auto g_zz_0_0_0_z_yz_0_z = buffer_2000_pdsp[320];

    auto g_zz_0_0_0_z_zz_0_x = buffer_2000_pdsp[321];

    auto g_zz_0_0_0_z_zz_0_y = buffer_2000_pdsp[322];

    auto g_zz_0_0_0_z_zz_0_z = buffer_2000_pdsp[323];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xx_0_0_0_x_xx_0_x, g_xx_0_0_0_x_xx_0_y, g_xx_0_0_0_x_xx_0_z, g_xxx_xx_0_x, g_xxx_xx_0_y, g_xxx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_0_x[i] = -6.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_xxx_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_y[i] = -6.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_xxx_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_0_z[i] = -6.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_xxx_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xx_0_0_0_x_xy_0_x, g_xx_0_0_0_x_xy_0_y, g_xx_0_0_0_x_xy_0_z, g_xxx_xy_0_x, g_xxx_xy_0_y, g_xxx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xy_0_x[i] = -6.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_xxx_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_y[i] = -6.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_xxx_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_z[i] = -6.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_xxx_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xx_0_0_0_x_xz_0_x, g_xx_0_0_0_x_xz_0_y, g_xx_0_0_0_x_xz_0_z, g_xxx_xz_0_x, g_xxx_xz_0_y, g_xxx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xz_0_x[i] = -6.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_xxx_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_y[i] = -6.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_xxx_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_z[i] = -6.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_xxx_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xx_0_0_0_x_yy_0_x, g_xx_0_0_0_x_yy_0_y, g_xx_0_0_0_x_yy_0_z, g_xxx_yy_0_x, g_xxx_yy_0_y, g_xxx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yy_0_x[i] = -6.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_xxx_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_y[i] = -6.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_xxx_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_z[i] = -6.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_xxx_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xx_0_0_0_x_yz_0_x, g_xx_0_0_0_x_yz_0_y, g_xx_0_0_0_x_yz_0_z, g_xxx_yz_0_x, g_xxx_yz_0_y, g_xxx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yz_0_x[i] = -6.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_xxx_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_y[i] = -6.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_xxx_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_z[i] = -6.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_xxx_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xx_0_0_0_x_zz_0_x, g_xx_0_0_0_x_zz_0_y, g_xx_0_0_0_x_zz_0_z, g_xxx_zz_0_x, g_xxx_zz_0_y, g_xxx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_zz_0_x[i] = -6.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_xxx_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_y[i] = -6.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_xxx_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_z[i] = -6.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_xxx_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_0_x, g_xx_0_0_0_y_xx_0_y, g_xx_0_0_0_y_xx_0_z, g_xxy_xx_0_x, g_xxy_xx_0_y, g_xxy_xx_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_0_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_xxy_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_xxy_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_0_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_xxy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xx_0_0_0_y_xy_0_x, g_xx_0_0_0_y_xy_0_y, g_xx_0_0_0_y_xy_0_z, g_xxy_xy_0_x, g_xxy_xy_0_y, g_xxy_xy_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xy_0_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_xxy_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_xxy_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_xxy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_xx_0_0_0_y_xz_0_x, g_xx_0_0_0_y_xz_0_y, g_xx_0_0_0_y_xz_0_z, g_xxy_xz_0_x, g_xxy_xz_0_y, g_xxy_xz_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xz_0_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_xxy_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_xxy_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_xxy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xx_0_0_0_y_yy_0_x, g_xx_0_0_0_y_yy_0_y, g_xx_0_0_0_y_yy_0_z, g_xxy_yy_0_x, g_xxy_yy_0_y, g_xxy_yy_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yy_0_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_xxy_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_xxy_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_xxy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xx_0_0_0_y_yz_0_x, g_xx_0_0_0_y_yz_0_y, g_xx_0_0_0_y_yz_0_z, g_xxy_yz_0_x, g_xxy_yz_0_y, g_xxy_yz_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yz_0_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_xxy_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_xxy_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_xxy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xx_0_0_0_y_zz_0_x, g_xx_0_0_0_y_zz_0_y, g_xx_0_0_0_y_zz_0_z, g_xxy_zz_0_x, g_xxy_zz_0_y, g_xxy_zz_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_zz_0_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_xxy_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_xxy_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_xxy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_0_x, g_xx_0_0_0_z_xx_0_y, g_xx_0_0_0_z_xx_0_z, g_xxz_xx_0_x, g_xxz_xx_0_y, g_xxz_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_0_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_xxz_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_xxz_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_0_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_xxz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_xx_0_0_0_z_xy_0_x, g_xx_0_0_0_z_xy_0_y, g_xx_0_0_0_z_xy_0_z, g_xxz_xy_0_x, g_xxz_xy_0_y, g_xxz_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xy_0_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_xxz_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_xxz_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_xxz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_xx_0_0_0_z_xz_0_x, g_xx_0_0_0_z_xz_0_y, g_xx_0_0_0_z_xz_0_z, g_xxz_xz_0_x, g_xxz_xz_0_y, g_xxz_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xz_0_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_xxz_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_xxz_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_xxz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xx_0_0_0_z_yy_0_x, g_xx_0_0_0_z_yy_0_y, g_xx_0_0_0_z_yy_0_z, g_xxz_yy_0_x, g_xxz_yy_0_y, g_xxz_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yy_0_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_xxz_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_xxz_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_xxz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_xx_0_0_0_z_yz_0_x, g_xx_0_0_0_z_yz_0_y, g_xx_0_0_0_z_yz_0_z, g_xxz_yz_0_x, g_xxz_yz_0_y, g_xxz_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yz_0_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_xxz_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_xxz_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_xxz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_xx_0_0_0_z_zz_0_x, g_xx_0_0_0_z_zz_0_y, g_xx_0_0_0_z_zz_0_z, g_xxz_zz_0_x, g_xxz_zz_0_y, g_xxz_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_zz_0_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_xxz_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_xxz_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_xxz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xxy_xx_0_x, g_xxy_xx_0_y, g_xxy_xx_0_z, g_xy_0_0_0_x_xx_0_x, g_xy_0_0_0_x_xx_0_y, g_xy_0_0_0_x_xx_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_0_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_xxy_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_xxy_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_0_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_xxy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xxy_xy_0_x, g_xxy_xy_0_y, g_xxy_xy_0_z, g_xy_0_0_0_x_xy_0_x, g_xy_0_0_0_x_xy_0_y, g_xy_0_0_0_x_xy_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xy_0_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_xxy_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_xxy_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_xxy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xxy_xz_0_x, g_xxy_xz_0_y, g_xxy_xz_0_z, g_xy_0_0_0_x_xz_0_x, g_xy_0_0_0_x_xz_0_y, g_xy_0_0_0_x_xz_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xz_0_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_xxy_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_xxy_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_xxy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xxy_yy_0_x, g_xxy_yy_0_y, g_xxy_yy_0_z, g_xy_0_0_0_x_yy_0_x, g_xy_0_0_0_x_yy_0_y, g_xy_0_0_0_x_yy_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yy_0_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_xxy_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_xxy_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_xxy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xxy_yz_0_x, g_xxy_yz_0_y, g_xxy_yz_0_z, g_xy_0_0_0_x_yz_0_x, g_xy_0_0_0_x_yz_0_y, g_xy_0_0_0_x_yz_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yz_0_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_xxy_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_xxy_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_xxy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xxy_zz_0_x, g_xxy_zz_0_y, g_xxy_zz_0_z, g_xy_0_0_0_x_zz_0_x, g_xy_0_0_0_x_zz_0_y, g_xy_0_0_0_x_zz_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_zz_0_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_xxy_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_xxy_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_xxy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xy_0_0_0_y_xx_0_x, g_xy_0_0_0_y_xx_0_y, g_xy_0_0_0_y_xx_0_z, g_xyy_xx_0_x, g_xyy_xx_0_y, g_xyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_0_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_xyy_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_xyy_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_0_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_xyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xy_0_0_0_y_xy_0_x, g_xy_0_0_0_y_xy_0_y, g_xy_0_0_0_y_xy_0_z, g_xyy_xy_0_x, g_xyy_xy_0_y, g_xyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xy_0_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_xyy_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_xyy_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_xyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xy_0_0_0_y_xz_0_x, g_xy_0_0_0_y_xz_0_y, g_xy_0_0_0_y_xz_0_z, g_xyy_xz_0_x, g_xyy_xz_0_y, g_xyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xz_0_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_xyy_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_xyy_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_xyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xy_0_0_0_y_yy_0_x, g_xy_0_0_0_y_yy_0_y, g_xy_0_0_0_y_yy_0_z, g_xyy_yy_0_x, g_xyy_yy_0_y, g_xyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yy_0_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_xyy_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_xyy_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_xyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xy_0_0_0_y_yz_0_x, g_xy_0_0_0_y_yz_0_y, g_xy_0_0_0_y_yz_0_z, g_xyy_yz_0_x, g_xyy_yz_0_y, g_xyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yz_0_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_xyy_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_xyy_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_xyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xy_0_0_0_y_zz_0_x, g_xy_0_0_0_y_zz_0_y, g_xy_0_0_0_y_zz_0_z, g_xyy_zz_0_x, g_xyy_zz_0_y, g_xyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_zz_0_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_xyy_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_xyy_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_xyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_0_x, g_xy_0_0_0_z_xx_0_y, g_xy_0_0_0_z_xx_0_z, g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_0_x[i] = 4.0 * g_xyz_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_y[i] = 4.0 * g_xyz_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_0_z[i] = 4.0 * g_xyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xy_0_0_0_z_xy_0_x, g_xy_0_0_0_z_xy_0_y, g_xy_0_0_0_z_xy_0_z, g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xy_0_x[i] = 4.0 * g_xyz_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_y[i] = 4.0 * g_xyz_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_z[i] = 4.0 * g_xyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xy_0_0_0_z_xz_0_x, g_xy_0_0_0_z_xz_0_y, g_xy_0_0_0_z_xz_0_z, g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xz_0_x[i] = 4.0 * g_xyz_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_y[i] = 4.0 * g_xyz_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_z[i] = 4.0 * g_xyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xy_0_0_0_z_yy_0_x, g_xy_0_0_0_z_yy_0_y, g_xy_0_0_0_z_yy_0_z, g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yy_0_x[i] = 4.0 * g_xyz_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_y[i] = 4.0 * g_xyz_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_z[i] = 4.0 * g_xyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xy_0_0_0_z_yz_0_x, g_xy_0_0_0_z_yz_0_y, g_xy_0_0_0_z_yz_0_z, g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yz_0_x[i] = 4.0 * g_xyz_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_y[i] = 4.0 * g_xyz_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_z[i] = 4.0 * g_xyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xy_0_0_0_z_zz_0_x, g_xy_0_0_0_z_zz_0_y, g_xy_0_0_0_z_zz_0_z, g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_zz_0_x[i] = 4.0 * g_xyz_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_y[i] = 4.0 * g_xyz_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_z[i] = 4.0 * g_xyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xxz_xx_0_x, g_xxz_xx_0_y, g_xxz_xx_0_z, g_xz_0_0_0_x_xx_0_x, g_xz_0_0_0_x_xx_0_y, g_xz_0_0_0_x_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_0_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_xxz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_xxz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_0_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_xxz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xxz_xy_0_x, g_xxz_xy_0_y, g_xxz_xy_0_z, g_xz_0_0_0_x_xy_0_x, g_xz_0_0_0_x_xy_0_y, g_xz_0_0_0_x_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xy_0_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_xxz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_xxz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_xxz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xxz_xz_0_x, g_xxz_xz_0_y, g_xxz_xz_0_z, g_xz_0_0_0_x_xz_0_x, g_xz_0_0_0_x_xz_0_y, g_xz_0_0_0_x_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xz_0_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_xxz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_xxz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_xxz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xxz_yy_0_x, g_xxz_yy_0_y, g_xxz_yy_0_z, g_xz_0_0_0_x_yy_0_x, g_xz_0_0_0_x_yy_0_y, g_xz_0_0_0_x_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yy_0_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_xxz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_xxz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_xxz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xxz_yz_0_x, g_xxz_yz_0_y, g_xxz_yz_0_z, g_xz_0_0_0_x_yz_0_x, g_xz_0_0_0_x_yz_0_y, g_xz_0_0_0_x_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yz_0_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_xxz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_xxz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_xxz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xxz_zz_0_x, g_xxz_zz_0_y, g_xxz_zz_0_z, g_xz_0_0_0_x_zz_0_x, g_xz_0_0_0_x_zz_0_y, g_xz_0_0_0_x_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_zz_0_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_xxz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_xxz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_xxz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z, g_xz_0_0_0_y_xx_0_x, g_xz_0_0_0_y_xx_0_y, g_xz_0_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_0_x[i] = 4.0 * g_xyz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_y[i] = 4.0 * g_xyz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_0_z[i] = 4.0 * g_xyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z, g_xz_0_0_0_y_xy_0_x, g_xz_0_0_0_y_xy_0_y, g_xz_0_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xy_0_x[i] = 4.0 * g_xyz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_y[i] = 4.0 * g_xyz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_z[i] = 4.0 * g_xyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z, g_xz_0_0_0_y_xz_0_x, g_xz_0_0_0_y_xz_0_y, g_xz_0_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xz_0_x[i] = 4.0 * g_xyz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_y[i] = 4.0 * g_xyz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_z[i] = 4.0 * g_xyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z, g_xz_0_0_0_y_yy_0_x, g_xz_0_0_0_y_yy_0_y, g_xz_0_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yy_0_x[i] = 4.0 * g_xyz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_y[i] = 4.0 * g_xyz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_z[i] = 4.0 * g_xyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z, g_xz_0_0_0_y_yz_0_x, g_xz_0_0_0_y_yz_0_y, g_xz_0_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yz_0_x[i] = 4.0 * g_xyz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_y[i] = 4.0 * g_xyz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_z[i] = 4.0 * g_xyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z, g_xz_0_0_0_y_zz_0_x, g_xz_0_0_0_y_zz_0_y, g_xz_0_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_zz_0_x[i] = 4.0 * g_xyz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_y[i] = 4.0 * g_xyz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_z[i] = 4.0 * g_xyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xz_0_0_0_z_xx_0_x, g_xz_0_0_0_z_xx_0_y, g_xz_0_0_0_z_xx_0_z, g_xzz_xx_0_x, g_xzz_xx_0_y, g_xzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_0_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_xzz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_xzz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_0_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_xzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xz_0_0_0_z_xy_0_x, g_xz_0_0_0_z_xy_0_y, g_xz_0_0_0_z_xy_0_z, g_xzz_xy_0_x, g_xzz_xy_0_y, g_xzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xy_0_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_xzz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_xzz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_xzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xz_0_0_0_z_xz_0_x, g_xz_0_0_0_z_xz_0_y, g_xz_0_0_0_z_xz_0_z, g_xzz_xz_0_x, g_xzz_xz_0_y, g_xzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xz_0_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_xzz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_xzz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_xzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xz_0_0_0_z_yy_0_x, g_xz_0_0_0_z_yy_0_y, g_xz_0_0_0_z_yy_0_z, g_xzz_yy_0_x, g_xzz_yy_0_y, g_xzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yy_0_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_xzz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_xzz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_xzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xz_0_0_0_z_yz_0_x, g_xz_0_0_0_z_yz_0_y, g_xz_0_0_0_z_yz_0_z, g_xzz_yz_0_x, g_xzz_yz_0_y, g_xzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yz_0_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_xzz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_xzz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_xzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xz_0_0_0_z_zz_0_x, g_xz_0_0_0_z_zz_0_y, g_xz_0_0_0_z_zz_0_z, g_xzz_zz_0_x, g_xzz_zz_0_y, g_xzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_zz_0_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_xzz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_xzz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_xzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xyy_xx_0_x, g_xyy_xx_0_y, g_xyy_xx_0_z, g_yy_0_0_0_x_xx_0_x, g_yy_0_0_0_x_xx_0_y, g_yy_0_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_0_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_xyy_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_xyy_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_0_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_xyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xyy_xy_0_x, g_xyy_xy_0_y, g_xyy_xy_0_z, g_yy_0_0_0_x_xy_0_x, g_yy_0_0_0_x_xy_0_y, g_yy_0_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xy_0_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_xyy_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_xyy_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_xyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xyy_xz_0_x, g_xyy_xz_0_y, g_xyy_xz_0_z, g_yy_0_0_0_x_xz_0_x, g_yy_0_0_0_x_xz_0_y, g_yy_0_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xz_0_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_xyy_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_xyy_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_xyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xyy_yy_0_x, g_xyy_yy_0_y, g_xyy_yy_0_z, g_yy_0_0_0_x_yy_0_x, g_yy_0_0_0_x_yy_0_y, g_yy_0_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yy_0_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_xyy_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_xyy_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_xyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xyy_yz_0_x, g_xyy_yz_0_y, g_xyy_yz_0_z, g_yy_0_0_0_x_yz_0_x, g_yy_0_0_0_x_yz_0_y, g_yy_0_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yz_0_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_xyy_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_xyy_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_xyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xyy_zz_0_x, g_xyy_zz_0_y, g_xyy_zz_0_z, g_yy_0_0_0_x_zz_0_x, g_yy_0_0_0_x_zz_0_y, g_yy_0_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_zz_0_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_xyy_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_xyy_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_xyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_yy_0_0_0_y_xx_0_x, g_yy_0_0_0_y_xx_0_y, g_yy_0_0_0_y_xx_0_z, g_yyy_xx_0_x, g_yyy_xx_0_y, g_yyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_0_x[i] = -6.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_yyy_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_y[i] = -6.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_yyy_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_0_z[i] = -6.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_yyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_yy_0_0_0_y_xy_0_x, g_yy_0_0_0_y_xy_0_y, g_yy_0_0_0_y_xy_0_z, g_yyy_xy_0_x, g_yyy_xy_0_y, g_yyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xy_0_x[i] = -6.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_yyy_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_y[i] = -6.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_yyy_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_z[i] = -6.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_yyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_yy_0_0_0_y_xz_0_x, g_yy_0_0_0_y_xz_0_y, g_yy_0_0_0_y_xz_0_z, g_yyy_xz_0_x, g_yyy_xz_0_y, g_yyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xz_0_x[i] = -6.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_yyy_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_y[i] = -6.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_yyy_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_z[i] = -6.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_yyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_yy_0_0_0_y_yy_0_x, g_yy_0_0_0_y_yy_0_y, g_yy_0_0_0_y_yy_0_z, g_yyy_yy_0_x, g_yyy_yy_0_y, g_yyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yy_0_x[i] = -6.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_yyy_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_y[i] = -6.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_yyy_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_z[i] = -6.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_yyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_yy_0_0_0_y_yz_0_x, g_yy_0_0_0_y_yz_0_y, g_yy_0_0_0_y_yz_0_z, g_yyy_yz_0_x, g_yyy_yz_0_y, g_yyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yz_0_x[i] = -6.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_yyy_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_y[i] = -6.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_yyy_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_z[i] = -6.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_yyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_yy_0_0_0_y_zz_0_x, g_yy_0_0_0_y_zz_0_y, g_yy_0_0_0_y_zz_0_z, g_yyy_zz_0_x, g_yyy_zz_0_y, g_yyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_zz_0_x[i] = -6.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_yyy_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_y[i] = -6.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_yyy_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_z[i] = -6.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_yyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_0_x, g_yy_0_0_0_z_xx_0_y, g_yy_0_0_0_z_xx_0_z, g_yyz_xx_0_x, g_yyz_xx_0_y, g_yyz_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_0_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_yyz_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_yyz_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_0_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_yyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_yy_0_0_0_z_xy_0_x, g_yy_0_0_0_z_xy_0_y, g_yy_0_0_0_z_xy_0_z, g_yyz_xy_0_x, g_yyz_xy_0_y, g_yyz_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xy_0_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_yyz_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_yyz_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_yyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_yy_0_0_0_z_xz_0_x, g_yy_0_0_0_z_xz_0_y, g_yy_0_0_0_z_xz_0_z, g_yyz_xz_0_x, g_yyz_xz_0_y, g_yyz_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xz_0_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_yyz_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_yyz_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_yyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_yy_0_0_0_z_yy_0_x, g_yy_0_0_0_z_yy_0_y, g_yy_0_0_0_z_yy_0_z, g_yyz_yy_0_x, g_yyz_yy_0_y, g_yyz_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yy_0_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_yyz_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_yyz_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_yyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_yy_0_0_0_z_yz_0_x, g_yy_0_0_0_z_yz_0_y, g_yy_0_0_0_z_yz_0_z, g_yyz_yz_0_x, g_yyz_yz_0_y, g_yyz_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yz_0_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_yyz_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_yyz_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_yyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_yy_0_0_0_z_zz_0_x, g_yy_0_0_0_z_zz_0_y, g_yy_0_0_0_z_zz_0_z, g_yyz_zz_0_x, g_yyz_zz_0_y, g_yyz_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_zz_0_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_yyz_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_yyz_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_yyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z, g_yz_0_0_0_x_xx_0_x, g_yz_0_0_0_x_xx_0_y, g_yz_0_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_0_x[i] = 4.0 * g_xyz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_y[i] = 4.0 * g_xyz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_0_z[i] = 4.0 * g_xyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z, g_yz_0_0_0_x_xy_0_x, g_yz_0_0_0_x_xy_0_y, g_yz_0_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xy_0_x[i] = 4.0 * g_xyz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_y[i] = 4.0 * g_xyz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_z[i] = 4.0 * g_xyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z, g_yz_0_0_0_x_xz_0_x, g_yz_0_0_0_x_xz_0_y, g_yz_0_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xz_0_x[i] = 4.0 * g_xyz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_y[i] = 4.0 * g_xyz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_z[i] = 4.0 * g_xyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z, g_yz_0_0_0_x_yy_0_x, g_yz_0_0_0_x_yy_0_y, g_yz_0_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yy_0_x[i] = 4.0 * g_xyz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_y[i] = 4.0 * g_xyz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_z[i] = 4.0 * g_xyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z, g_yz_0_0_0_x_yz_0_x, g_yz_0_0_0_x_yz_0_y, g_yz_0_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yz_0_x[i] = 4.0 * g_xyz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_y[i] = 4.0 * g_xyz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_z[i] = 4.0 * g_xyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z, g_yz_0_0_0_x_zz_0_x, g_yz_0_0_0_x_zz_0_y, g_yz_0_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_zz_0_x[i] = 4.0 * g_xyz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_y[i] = 4.0 * g_xyz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_z[i] = 4.0 * g_xyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_yyz_xx_0_x, g_yyz_xx_0_y, g_yyz_xx_0_z, g_yz_0_0_0_y_xx_0_x, g_yz_0_0_0_y_xx_0_y, g_yz_0_0_0_y_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_0_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_yyz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_yyz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_0_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_yyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_yyz_xy_0_x, g_yyz_xy_0_y, g_yyz_xy_0_z, g_yz_0_0_0_y_xy_0_x, g_yz_0_0_0_y_xy_0_y, g_yz_0_0_0_y_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xy_0_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_yyz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_yyz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_yyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_yyz_xz_0_x, g_yyz_xz_0_y, g_yyz_xz_0_z, g_yz_0_0_0_y_xz_0_x, g_yz_0_0_0_y_xz_0_y, g_yz_0_0_0_y_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xz_0_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_yyz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_yyz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_yyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_yyz_yy_0_x, g_yyz_yy_0_y, g_yyz_yy_0_z, g_yz_0_0_0_y_yy_0_x, g_yz_0_0_0_y_yy_0_y, g_yz_0_0_0_y_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yy_0_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_yyz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_yyz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_yyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_yyz_yz_0_x, g_yyz_yz_0_y, g_yyz_yz_0_z, g_yz_0_0_0_y_yz_0_x, g_yz_0_0_0_y_yz_0_y, g_yz_0_0_0_y_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yz_0_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_yyz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_yyz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_yyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_yyz_zz_0_x, g_yyz_zz_0_y, g_yyz_zz_0_z, g_yz_0_0_0_y_zz_0_x, g_yz_0_0_0_y_zz_0_y, g_yz_0_0_0_y_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_zz_0_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_yyz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_yyz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_yyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_yz_0_0_0_z_xx_0_x, g_yz_0_0_0_z_xx_0_y, g_yz_0_0_0_z_xx_0_z, g_yzz_xx_0_x, g_yzz_xx_0_y, g_yzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_0_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_yzz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_yzz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_0_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_yzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_yz_0_0_0_z_xy_0_x, g_yz_0_0_0_z_xy_0_y, g_yz_0_0_0_z_xy_0_z, g_yzz_xy_0_x, g_yzz_xy_0_y, g_yzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xy_0_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_yzz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_yzz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_yzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_yz_0_0_0_z_xz_0_x, g_yz_0_0_0_z_xz_0_y, g_yz_0_0_0_z_xz_0_z, g_yzz_xz_0_x, g_yzz_xz_0_y, g_yzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xz_0_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_yzz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_yzz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_yzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_yz_0_0_0_z_yy_0_x, g_yz_0_0_0_z_yy_0_y, g_yz_0_0_0_z_yy_0_z, g_yzz_yy_0_x, g_yzz_yy_0_y, g_yzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yy_0_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_yzz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_yzz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_yzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_yz_0_0_0_z_yz_0_x, g_yz_0_0_0_z_yz_0_y, g_yz_0_0_0_z_yz_0_z, g_yzz_yz_0_x, g_yzz_yz_0_y, g_yzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yz_0_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_yzz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_yzz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_yzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_yz_0_0_0_z_zz_0_x, g_yz_0_0_0_z_zz_0_y, g_yz_0_0_0_z_zz_0_z, g_yzz_zz_0_x, g_yzz_zz_0_y, g_yzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_zz_0_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_yzz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_yzz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_yzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xzz_xx_0_x, g_xzz_xx_0_y, g_xzz_xx_0_z, g_zz_0_0_0_x_xx_0_x, g_zz_0_0_0_x_xx_0_y, g_zz_0_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_0_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_xzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_xzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_0_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_xzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xzz_xy_0_x, g_xzz_xy_0_y, g_xzz_xy_0_z, g_zz_0_0_0_x_xy_0_x, g_zz_0_0_0_x_xy_0_y, g_zz_0_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xy_0_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_xzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_xzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_xzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xzz_xz_0_x, g_xzz_xz_0_y, g_xzz_xz_0_z, g_zz_0_0_0_x_xz_0_x, g_zz_0_0_0_x_xz_0_y, g_zz_0_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xz_0_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_xzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_xzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_xzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xzz_yy_0_x, g_xzz_yy_0_y, g_xzz_yy_0_z, g_zz_0_0_0_x_yy_0_x, g_zz_0_0_0_x_yy_0_y, g_zz_0_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yy_0_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_xzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_xzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_xzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xzz_yz_0_x, g_xzz_yz_0_y, g_xzz_yz_0_z, g_zz_0_0_0_x_yz_0_x, g_zz_0_0_0_x_yz_0_y, g_zz_0_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yz_0_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_xzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_xzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_xzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xzz_zz_0_x, g_xzz_zz_0_y, g_xzz_zz_0_z, g_zz_0_0_0_x_zz_0_x, g_zz_0_0_0_x_zz_0_y, g_zz_0_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_zz_0_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_xzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_xzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_xzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_yzz_xx_0_x, g_yzz_xx_0_y, g_yzz_xx_0_z, g_zz_0_0_0_y_xx_0_x, g_zz_0_0_0_y_xx_0_y, g_zz_0_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_0_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_yzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_yzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_0_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_yzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_yzz_xy_0_x, g_yzz_xy_0_y, g_yzz_xy_0_z, g_zz_0_0_0_y_xy_0_x, g_zz_0_0_0_y_xy_0_y, g_zz_0_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xy_0_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_yzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_yzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_yzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_yzz_xz_0_x, g_yzz_xz_0_y, g_yzz_xz_0_z, g_zz_0_0_0_y_xz_0_x, g_zz_0_0_0_y_xz_0_y, g_zz_0_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xz_0_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_yzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_yzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_yzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_yzz_yy_0_x, g_yzz_yy_0_y, g_yzz_yy_0_z, g_zz_0_0_0_y_yy_0_x, g_zz_0_0_0_y_yy_0_y, g_zz_0_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yy_0_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_yzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_yzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_yzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_yzz_yz_0_x, g_yzz_yz_0_y, g_yzz_yz_0_z, g_zz_0_0_0_y_yz_0_x, g_zz_0_0_0_y_yz_0_y, g_zz_0_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yz_0_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_yzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_yzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_yzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_yzz_zz_0_x, g_yzz_zz_0_y, g_yzz_zz_0_z, g_zz_0_0_0_y_zz_0_x, g_zz_0_0_0_y_zz_0_y, g_zz_0_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_zz_0_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_yzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_yzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_yzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z, g_zz_0_0_0_z_xx_0_x, g_zz_0_0_0_z_xx_0_y, g_zz_0_0_0_z_xx_0_z, g_zzz_xx_0_x, g_zzz_xx_0_y, g_zzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_0_x[i] = -6.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_zzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_y[i] = -6.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_zzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_0_z[i] = -6.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_zzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_zz_0_0_0_z_xy_0_x, g_zz_0_0_0_z_xy_0_y, g_zz_0_0_0_z_xy_0_z, g_zzz_xy_0_x, g_zzz_xy_0_y, g_zzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xy_0_x[i] = -6.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_zzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_y[i] = -6.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_zzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_z[i] = -6.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_zzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_zz_0_0_0_z_xz_0_x, g_zz_0_0_0_z_xz_0_y, g_zz_0_0_0_z_xz_0_z, g_zzz_xz_0_x, g_zzz_xz_0_y, g_zzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xz_0_x[i] = -6.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_zzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_y[i] = -6.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_zzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_z[i] = -6.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_zzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z, g_zz_0_0_0_z_yy_0_x, g_zz_0_0_0_z_yy_0_y, g_zz_0_0_0_z_yy_0_z, g_zzz_yy_0_x, g_zzz_yy_0_y, g_zzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yy_0_x[i] = -6.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_zzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_y[i] = -6.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_zzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_z[i] = -6.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_zzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_zz_0_0_0_z_yz_0_x, g_zz_0_0_0_z_yz_0_y, g_zz_0_0_0_z_yz_0_z, g_zzz_yz_0_x, g_zzz_yz_0_y, g_zzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yz_0_x[i] = -6.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_zzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_y[i] = -6.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_zzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_z[i] = -6.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_zzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z, g_zz_0_0_0_z_zz_0_x, g_zz_0_0_0_z_zz_0_y, g_zz_0_0_0_z_zz_0_z, g_zzz_zz_0_x, g_zzz_zz_0_y, g_zzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_zz_0_x[i] = -6.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_zzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_y[i] = -6.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_zzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_z[i] = -6.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_zzz_zz_0_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


#include "GeomDeriv1000OfScalarForDDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ddsp_0(CSimdArray<double>& buffer_1000_ddsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_fdsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ddsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ddsp

    auto g_x_0_0_0_xx_xx_0_x = buffer_1000_ddsp[0];

    auto g_x_0_0_0_xx_xx_0_y = buffer_1000_ddsp[1];

    auto g_x_0_0_0_xx_xx_0_z = buffer_1000_ddsp[2];

    auto g_x_0_0_0_xx_xy_0_x = buffer_1000_ddsp[3];

    auto g_x_0_0_0_xx_xy_0_y = buffer_1000_ddsp[4];

    auto g_x_0_0_0_xx_xy_0_z = buffer_1000_ddsp[5];

    auto g_x_0_0_0_xx_xz_0_x = buffer_1000_ddsp[6];

    auto g_x_0_0_0_xx_xz_0_y = buffer_1000_ddsp[7];

    auto g_x_0_0_0_xx_xz_0_z = buffer_1000_ddsp[8];

    auto g_x_0_0_0_xx_yy_0_x = buffer_1000_ddsp[9];

    auto g_x_0_0_0_xx_yy_0_y = buffer_1000_ddsp[10];

    auto g_x_0_0_0_xx_yy_0_z = buffer_1000_ddsp[11];

    auto g_x_0_0_0_xx_yz_0_x = buffer_1000_ddsp[12];

    auto g_x_0_0_0_xx_yz_0_y = buffer_1000_ddsp[13];

    auto g_x_0_0_0_xx_yz_0_z = buffer_1000_ddsp[14];

    auto g_x_0_0_0_xx_zz_0_x = buffer_1000_ddsp[15];

    auto g_x_0_0_0_xx_zz_0_y = buffer_1000_ddsp[16];

    auto g_x_0_0_0_xx_zz_0_z = buffer_1000_ddsp[17];

    auto g_x_0_0_0_xy_xx_0_x = buffer_1000_ddsp[18];

    auto g_x_0_0_0_xy_xx_0_y = buffer_1000_ddsp[19];

    auto g_x_0_0_0_xy_xx_0_z = buffer_1000_ddsp[20];

    auto g_x_0_0_0_xy_xy_0_x = buffer_1000_ddsp[21];

    auto g_x_0_0_0_xy_xy_0_y = buffer_1000_ddsp[22];

    auto g_x_0_0_0_xy_xy_0_z = buffer_1000_ddsp[23];

    auto g_x_0_0_0_xy_xz_0_x = buffer_1000_ddsp[24];

    auto g_x_0_0_0_xy_xz_0_y = buffer_1000_ddsp[25];

    auto g_x_0_0_0_xy_xz_0_z = buffer_1000_ddsp[26];

    auto g_x_0_0_0_xy_yy_0_x = buffer_1000_ddsp[27];

    auto g_x_0_0_0_xy_yy_0_y = buffer_1000_ddsp[28];

    auto g_x_0_0_0_xy_yy_0_z = buffer_1000_ddsp[29];

    auto g_x_0_0_0_xy_yz_0_x = buffer_1000_ddsp[30];

    auto g_x_0_0_0_xy_yz_0_y = buffer_1000_ddsp[31];

    auto g_x_0_0_0_xy_yz_0_z = buffer_1000_ddsp[32];

    auto g_x_0_0_0_xy_zz_0_x = buffer_1000_ddsp[33];

    auto g_x_0_0_0_xy_zz_0_y = buffer_1000_ddsp[34];

    auto g_x_0_0_0_xy_zz_0_z = buffer_1000_ddsp[35];

    auto g_x_0_0_0_xz_xx_0_x = buffer_1000_ddsp[36];

    auto g_x_0_0_0_xz_xx_0_y = buffer_1000_ddsp[37];

    auto g_x_0_0_0_xz_xx_0_z = buffer_1000_ddsp[38];

    auto g_x_0_0_0_xz_xy_0_x = buffer_1000_ddsp[39];

    auto g_x_0_0_0_xz_xy_0_y = buffer_1000_ddsp[40];

    auto g_x_0_0_0_xz_xy_0_z = buffer_1000_ddsp[41];

    auto g_x_0_0_0_xz_xz_0_x = buffer_1000_ddsp[42];

    auto g_x_0_0_0_xz_xz_0_y = buffer_1000_ddsp[43];

    auto g_x_0_0_0_xz_xz_0_z = buffer_1000_ddsp[44];

    auto g_x_0_0_0_xz_yy_0_x = buffer_1000_ddsp[45];

    auto g_x_0_0_0_xz_yy_0_y = buffer_1000_ddsp[46];

    auto g_x_0_0_0_xz_yy_0_z = buffer_1000_ddsp[47];

    auto g_x_0_0_0_xz_yz_0_x = buffer_1000_ddsp[48];

    auto g_x_0_0_0_xz_yz_0_y = buffer_1000_ddsp[49];

    auto g_x_0_0_0_xz_yz_0_z = buffer_1000_ddsp[50];

    auto g_x_0_0_0_xz_zz_0_x = buffer_1000_ddsp[51];

    auto g_x_0_0_0_xz_zz_0_y = buffer_1000_ddsp[52];

    auto g_x_0_0_0_xz_zz_0_z = buffer_1000_ddsp[53];

    auto g_x_0_0_0_yy_xx_0_x = buffer_1000_ddsp[54];

    auto g_x_0_0_0_yy_xx_0_y = buffer_1000_ddsp[55];

    auto g_x_0_0_0_yy_xx_0_z = buffer_1000_ddsp[56];

    auto g_x_0_0_0_yy_xy_0_x = buffer_1000_ddsp[57];

    auto g_x_0_0_0_yy_xy_0_y = buffer_1000_ddsp[58];

    auto g_x_0_0_0_yy_xy_0_z = buffer_1000_ddsp[59];

    auto g_x_0_0_0_yy_xz_0_x = buffer_1000_ddsp[60];

    auto g_x_0_0_0_yy_xz_0_y = buffer_1000_ddsp[61];

    auto g_x_0_0_0_yy_xz_0_z = buffer_1000_ddsp[62];

    auto g_x_0_0_0_yy_yy_0_x = buffer_1000_ddsp[63];

    auto g_x_0_0_0_yy_yy_0_y = buffer_1000_ddsp[64];

    auto g_x_0_0_0_yy_yy_0_z = buffer_1000_ddsp[65];

    auto g_x_0_0_0_yy_yz_0_x = buffer_1000_ddsp[66];

    auto g_x_0_0_0_yy_yz_0_y = buffer_1000_ddsp[67];

    auto g_x_0_0_0_yy_yz_0_z = buffer_1000_ddsp[68];

    auto g_x_0_0_0_yy_zz_0_x = buffer_1000_ddsp[69];

    auto g_x_0_0_0_yy_zz_0_y = buffer_1000_ddsp[70];

    auto g_x_0_0_0_yy_zz_0_z = buffer_1000_ddsp[71];

    auto g_x_0_0_0_yz_xx_0_x = buffer_1000_ddsp[72];

    auto g_x_0_0_0_yz_xx_0_y = buffer_1000_ddsp[73];

    auto g_x_0_0_0_yz_xx_0_z = buffer_1000_ddsp[74];

    auto g_x_0_0_0_yz_xy_0_x = buffer_1000_ddsp[75];

    auto g_x_0_0_0_yz_xy_0_y = buffer_1000_ddsp[76];

    auto g_x_0_0_0_yz_xy_0_z = buffer_1000_ddsp[77];

    auto g_x_0_0_0_yz_xz_0_x = buffer_1000_ddsp[78];

    auto g_x_0_0_0_yz_xz_0_y = buffer_1000_ddsp[79];

    auto g_x_0_0_0_yz_xz_0_z = buffer_1000_ddsp[80];

    auto g_x_0_0_0_yz_yy_0_x = buffer_1000_ddsp[81];

    auto g_x_0_0_0_yz_yy_0_y = buffer_1000_ddsp[82];

    auto g_x_0_0_0_yz_yy_0_z = buffer_1000_ddsp[83];

    auto g_x_0_0_0_yz_yz_0_x = buffer_1000_ddsp[84];

    auto g_x_0_0_0_yz_yz_0_y = buffer_1000_ddsp[85];

    auto g_x_0_0_0_yz_yz_0_z = buffer_1000_ddsp[86];

    auto g_x_0_0_0_yz_zz_0_x = buffer_1000_ddsp[87];

    auto g_x_0_0_0_yz_zz_0_y = buffer_1000_ddsp[88];

    auto g_x_0_0_0_yz_zz_0_z = buffer_1000_ddsp[89];

    auto g_x_0_0_0_zz_xx_0_x = buffer_1000_ddsp[90];

    auto g_x_0_0_0_zz_xx_0_y = buffer_1000_ddsp[91];

    auto g_x_0_0_0_zz_xx_0_z = buffer_1000_ddsp[92];

    auto g_x_0_0_0_zz_xy_0_x = buffer_1000_ddsp[93];

    auto g_x_0_0_0_zz_xy_0_y = buffer_1000_ddsp[94];

    auto g_x_0_0_0_zz_xy_0_z = buffer_1000_ddsp[95];

    auto g_x_0_0_0_zz_xz_0_x = buffer_1000_ddsp[96];

    auto g_x_0_0_0_zz_xz_0_y = buffer_1000_ddsp[97];

    auto g_x_0_0_0_zz_xz_0_z = buffer_1000_ddsp[98];

    auto g_x_0_0_0_zz_yy_0_x = buffer_1000_ddsp[99];

    auto g_x_0_0_0_zz_yy_0_y = buffer_1000_ddsp[100];

    auto g_x_0_0_0_zz_yy_0_z = buffer_1000_ddsp[101];

    auto g_x_0_0_0_zz_yz_0_x = buffer_1000_ddsp[102];

    auto g_x_0_0_0_zz_yz_0_y = buffer_1000_ddsp[103];

    auto g_x_0_0_0_zz_yz_0_z = buffer_1000_ddsp[104];

    auto g_x_0_0_0_zz_zz_0_x = buffer_1000_ddsp[105];

    auto g_x_0_0_0_zz_zz_0_y = buffer_1000_ddsp[106];

    auto g_x_0_0_0_zz_zz_0_z = buffer_1000_ddsp[107];

    auto g_y_0_0_0_xx_xx_0_x = buffer_1000_ddsp[108];

    auto g_y_0_0_0_xx_xx_0_y = buffer_1000_ddsp[109];

    auto g_y_0_0_0_xx_xx_0_z = buffer_1000_ddsp[110];

    auto g_y_0_0_0_xx_xy_0_x = buffer_1000_ddsp[111];

    auto g_y_0_0_0_xx_xy_0_y = buffer_1000_ddsp[112];

    auto g_y_0_0_0_xx_xy_0_z = buffer_1000_ddsp[113];

    auto g_y_0_0_0_xx_xz_0_x = buffer_1000_ddsp[114];

    auto g_y_0_0_0_xx_xz_0_y = buffer_1000_ddsp[115];

    auto g_y_0_0_0_xx_xz_0_z = buffer_1000_ddsp[116];

    auto g_y_0_0_0_xx_yy_0_x = buffer_1000_ddsp[117];

    auto g_y_0_0_0_xx_yy_0_y = buffer_1000_ddsp[118];

    auto g_y_0_0_0_xx_yy_0_z = buffer_1000_ddsp[119];

    auto g_y_0_0_0_xx_yz_0_x = buffer_1000_ddsp[120];

    auto g_y_0_0_0_xx_yz_0_y = buffer_1000_ddsp[121];

    auto g_y_0_0_0_xx_yz_0_z = buffer_1000_ddsp[122];

    auto g_y_0_0_0_xx_zz_0_x = buffer_1000_ddsp[123];

    auto g_y_0_0_0_xx_zz_0_y = buffer_1000_ddsp[124];

    auto g_y_0_0_0_xx_zz_0_z = buffer_1000_ddsp[125];

    auto g_y_0_0_0_xy_xx_0_x = buffer_1000_ddsp[126];

    auto g_y_0_0_0_xy_xx_0_y = buffer_1000_ddsp[127];

    auto g_y_0_0_0_xy_xx_0_z = buffer_1000_ddsp[128];

    auto g_y_0_0_0_xy_xy_0_x = buffer_1000_ddsp[129];

    auto g_y_0_0_0_xy_xy_0_y = buffer_1000_ddsp[130];

    auto g_y_0_0_0_xy_xy_0_z = buffer_1000_ddsp[131];

    auto g_y_0_0_0_xy_xz_0_x = buffer_1000_ddsp[132];

    auto g_y_0_0_0_xy_xz_0_y = buffer_1000_ddsp[133];

    auto g_y_0_0_0_xy_xz_0_z = buffer_1000_ddsp[134];

    auto g_y_0_0_0_xy_yy_0_x = buffer_1000_ddsp[135];

    auto g_y_0_0_0_xy_yy_0_y = buffer_1000_ddsp[136];

    auto g_y_0_0_0_xy_yy_0_z = buffer_1000_ddsp[137];

    auto g_y_0_0_0_xy_yz_0_x = buffer_1000_ddsp[138];

    auto g_y_0_0_0_xy_yz_0_y = buffer_1000_ddsp[139];

    auto g_y_0_0_0_xy_yz_0_z = buffer_1000_ddsp[140];

    auto g_y_0_0_0_xy_zz_0_x = buffer_1000_ddsp[141];

    auto g_y_0_0_0_xy_zz_0_y = buffer_1000_ddsp[142];

    auto g_y_0_0_0_xy_zz_0_z = buffer_1000_ddsp[143];

    auto g_y_0_0_0_xz_xx_0_x = buffer_1000_ddsp[144];

    auto g_y_0_0_0_xz_xx_0_y = buffer_1000_ddsp[145];

    auto g_y_0_0_0_xz_xx_0_z = buffer_1000_ddsp[146];

    auto g_y_0_0_0_xz_xy_0_x = buffer_1000_ddsp[147];

    auto g_y_0_0_0_xz_xy_0_y = buffer_1000_ddsp[148];

    auto g_y_0_0_0_xz_xy_0_z = buffer_1000_ddsp[149];

    auto g_y_0_0_0_xz_xz_0_x = buffer_1000_ddsp[150];

    auto g_y_0_0_0_xz_xz_0_y = buffer_1000_ddsp[151];

    auto g_y_0_0_0_xz_xz_0_z = buffer_1000_ddsp[152];

    auto g_y_0_0_0_xz_yy_0_x = buffer_1000_ddsp[153];

    auto g_y_0_0_0_xz_yy_0_y = buffer_1000_ddsp[154];

    auto g_y_0_0_0_xz_yy_0_z = buffer_1000_ddsp[155];

    auto g_y_0_0_0_xz_yz_0_x = buffer_1000_ddsp[156];

    auto g_y_0_0_0_xz_yz_0_y = buffer_1000_ddsp[157];

    auto g_y_0_0_0_xz_yz_0_z = buffer_1000_ddsp[158];

    auto g_y_0_0_0_xz_zz_0_x = buffer_1000_ddsp[159];

    auto g_y_0_0_0_xz_zz_0_y = buffer_1000_ddsp[160];

    auto g_y_0_0_0_xz_zz_0_z = buffer_1000_ddsp[161];

    auto g_y_0_0_0_yy_xx_0_x = buffer_1000_ddsp[162];

    auto g_y_0_0_0_yy_xx_0_y = buffer_1000_ddsp[163];

    auto g_y_0_0_0_yy_xx_0_z = buffer_1000_ddsp[164];

    auto g_y_0_0_0_yy_xy_0_x = buffer_1000_ddsp[165];

    auto g_y_0_0_0_yy_xy_0_y = buffer_1000_ddsp[166];

    auto g_y_0_0_0_yy_xy_0_z = buffer_1000_ddsp[167];

    auto g_y_0_0_0_yy_xz_0_x = buffer_1000_ddsp[168];

    auto g_y_0_0_0_yy_xz_0_y = buffer_1000_ddsp[169];

    auto g_y_0_0_0_yy_xz_0_z = buffer_1000_ddsp[170];

    auto g_y_0_0_0_yy_yy_0_x = buffer_1000_ddsp[171];

    auto g_y_0_0_0_yy_yy_0_y = buffer_1000_ddsp[172];

    auto g_y_0_0_0_yy_yy_0_z = buffer_1000_ddsp[173];

    auto g_y_0_0_0_yy_yz_0_x = buffer_1000_ddsp[174];

    auto g_y_0_0_0_yy_yz_0_y = buffer_1000_ddsp[175];

    auto g_y_0_0_0_yy_yz_0_z = buffer_1000_ddsp[176];

    auto g_y_0_0_0_yy_zz_0_x = buffer_1000_ddsp[177];

    auto g_y_0_0_0_yy_zz_0_y = buffer_1000_ddsp[178];

    auto g_y_0_0_0_yy_zz_0_z = buffer_1000_ddsp[179];

    auto g_y_0_0_0_yz_xx_0_x = buffer_1000_ddsp[180];

    auto g_y_0_0_0_yz_xx_0_y = buffer_1000_ddsp[181];

    auto g_y_0_0_0_yz_xx_0_z = buffer_1000_ddsp[182];

    auto g_y_0_0_0_yz_xy_0_x = buffer_1000_ddsp[183];

    auto g_y_0_0_0_yz_xy_0_y = buffer_1000_ddsp[184];

    auto g_y_0_0_0_yz_xy_0_z = buffer_1000_ddsp[185];

    auto g_y_0_0_0_yz_xz_0_x = buffer_1000_ddsp[186];

    auto g_y_0_0_0_yz_xz_0_y = buffer_1000_ddsp[187];

    auto g_y_0_0_0_yz_xz_0_z = buffer_1000_ddsp[188];

    auto g_y_0_0_0_yz_yy_0_x = buffer_1000_ddsp[189];

    auto g_y_0_0_0_yz_yy_0_y = buffer_1000_ddsp[190];

    auto g_y_0_0_0_yz_yy_0_z = buffer_1000_ddsp[191];

    auto g_y_0_0_0_yz_yz_0_x = buffer_1000_ddsp[192];

    auto g_y_0_0_0_yz_yz_0_y = buffer_1000_ddsp[193];

    auto g_y_0_0_0_yz_yz_0_z = buffer_1000_ddsp[194];

    auto g_y_0_0_0_yz_zz_0_x = buffer_1000_ddsp[195];

    auto g_y_0_0_0_yz_zz_0_y = buffer_1000_ddsp[196];

    auto g_y_0_0_0_yz_zz_0_z = buffer_1000_ddsp[197];

    auto g_y_0_0_0_zz_xx_0_x = buffer_1000_ddsp[198];

    auto g_y_0_0_0_zz_xx_0_y = buffer_1000_ddsp[199];

    auto g_y_0_0_0_zz_xx_0_z = buffer_1000_ddsp[200];

    auto g_y_0_0_0_zz_xy_0_x = buffer_1000_ddsp[201];

    auto g_y_0_0_0_zz_xy_0_y = buffer_1000_ddsp[202];

    auto g_y_0_0_0_zz_xy_0_z = buffer_1000_ddsp[203];

    auto g_y_0_0_0_zz_xz_0_x = buffer_1000_ddsp[204];

    auto g_y_0_0_0_zz_xz_0_y = buffer_1000_ddsp[205];

    auto g_y_0_0_0_zz_xz_0_z = buffer_1000_ddsp[206];

    auto g_y_0_0_0_zz_yy_0_x = buffer_1000_ddsp[207];

    auto g_y_0_0_0_zz_yy_0_y = buffer_1000_ddsp[208];

    auto g_y_0_0_0_zz_yy_0_z = buffer_1000_ddsp[209];

    auto g_y_0_0_0_zz_yz_0_x = buffer_1000_ddsp[210];

    auto g_y_0_0_0_zz_yz_0_y = buffer_1000_ddsp[211];

    auto g_y_0_0_0_zz_yz_0_z = buffer_1000_ddsp[212];

    auto g_y_0_0_0_zz_zz_0_x = buffer_1000_ddsp[213];

    auto g_y_0_0_0_zz_zz_0_y = buffer_1000_ddsp[214];

    auto g_y_0_0_0_zz_zz_0_z = buffer_1000_ddsp[215];

    auto g_z_0_0_0_xx_xx_0_x = buffer_1000_ddsp[216];

    auto g_z_0_0_0_xx_xx_0_y = buffer_1000_ddsp[217];

    auto g_z_0_0_0_xx_xx_0_z = buffer_1000_ddsp[218];

    auto g_z_0_0_0_xx_xy_0_x = buffer_1000_ddsp[219];

    auto g_z_0_0_0_xx_xy_0_y = buffer_1000_ddsp[220];

    auto g_z_0_0_0_xx_xy_0_z = buffer_1000_ddsp[221];

    auto g_z_0_0_0_xx_xz_0_x = buffer_1000_ddsp[222];

    auto g_z_0_0_0_xx_xz_0_y = buffer_1000_ddsp[223];

    auto g_z_0_0_0_xx_xz_0_z = buffer_1000_ddsp[224];

    auto g_z_0_0_0_xx_yy_0_x = buffer_1000_ddsp[225];

    auto g_z_0_0_0_xx_yy_0_y = buffer_1000_ddsp[226];

    auto g_z_0_0_0_xx_yy_0_z = buffer_1000_ddsp[227];

    auto g_z_0_0_0_xx_yz_0_x = buffer_1000_ddsp[228];

    auto g_z_0_0_0_xx_yz_0_y = buffer_1000_ddsp[229];

    auto g_z_0_0_0_xx_yz_0_z = buffer_1000_ddsp[230];

    auto g_z_0_0_0_xx_zz_0_x = buffer_1000_ddsp[231];

    auto g_z_0_0_0_xx_zz_0_y = buffer_1000_ddsp[232];

    auto g_z_0_0_0_xx_zz_0_z = buffer_1000_ddsp[233];

    auto g_z_0_0_0_xy_xx_0_x = buffer_1000_ddsp[234];

    auto g_z_0_0_0_xy_xx_0_y = buffer_1000_ddsp[235];

    auto g_z_0_0_0_xy_xx_0_z = buffer_1000_ddsp[236];

    auto g_z_0_0_0_xy_xy_0_x = buffer_1000_ddsp[237];

    auto g_z_0_0_0_xy_xy_0_y = buffer_1000_ddsp[238];

    auto g_z_0_0_0_xy_xy_0_z = buffer_1000_ddsp[239];

    auto g_z_0_0_0_xy_xz_0_x = buffer_1000_ddsp[240];

    auto g_z_0_0_0_xy_xz_0_y = buffer_1000_ddsp[241];

    auto g_z_0_0_0_xy_xz_0_z = buffer_1000_ddsp[242];

    auto g_z_0_0_0_xy_yy_0_x = buffer_1000_ddsp[243];

    auto g_z_0_0_0_xy_yy_0_y = buffer_1000_ddsp[244];

    auto g_z_0_0_0_xy_yy_0_z = buffer_1000_ddsp[245];

    auto g_z_0_0_0_xy_yz_0_x = buffer_1000_ddsp[246];

    auto g_z_0_0_0_xy_yz_0_y = buffer_1000_ddsp[247];

    auto g_z_0_0_0_xy_yz_0_z = buffer_1000_ddsp[248];

    auto g_z_0_0_0_xy_zz_0_x = buffer_1000_ddsp[249];

    auto g_z_0_0_0_xy_zz_0_y = buffer_1000_ddsp[250];

    auto g_z_0_0_0_xy_zz_0_z = buffer_1000_ddsp[251];

    auto g_z_0_0_0_xz_xx_0_x = buffer_1000_ddsp[252];

    auto g_z_0_0_0_xz_xx_0_y = buffer_1000_ddsp[253];

    auto g_z_0_0_0_xz_xx_0_z = buffer_1000_ddsp[254];

    auto g_z_0_0_0_xz_xy_0_x = buffer_1000_ddsp[255];

    auto g_z_0_0_0_xz_xy_0_y = buffer_1000_ddsp[256];

    auto g_z_0_0_0_xz_xy_0_z = buffer_1000_ddsp[257];

    auto g_z_0_0_0_xz_xz_0_x = buffer_1000_ddsp[258];

    auto g_z_0_0_0_xz_xz_0_y = buffer_1000_ddsp[259];

    auto g_z_0_0_0_xz_xz_0_z = buffer_1000_ddsp[260];

    auto g_z_0_0_0_xz_yy_0_x = buffer_1000_ddsp[261];

    auto g_z_0_0_0_xz_yy_0_y = buffer_1000_ddsp[262];

    auto g_z_0_0_0_xz_yy_0_z = buffer_1000_ddsp[263];

    auto g_z_0_0_0_xz_yz_0_x = buffer_1000_ddsp[264];

    auto g_z_0_0_0_xz_yz_0_y = buffer_1000_ddsp[265];

    auto g_z_0_0_0_xz_yz_0_z = buffer_1000_ddsp[266];

    auto g_z_0_0_0_xz_zz_0_x = buffer_1000_ddsp[267];

    auto g_z_0_0_0_xz_zz_0_y = buffer_1000_ddsp[268];

    auto g_z_0_0_0_xz_zz_0_z = buffer_1000_ddsp[269];

    auto g_z_0_0_0_yy_xx_0_x = buffer_1000_ddsp[270];

    auto g_z_0_0_0_yy_xx_0_y = buffer_1000_ddsp[271];

    auto g_z_0_0_0_yy_xx_0_z = buffer_1000_ddsp[272];

    auto g_z_0_0_0_yy_xy_0_x = buffer_1000_ddsp[273];

    auto g_z_0_0_0_yy_xy_0_y = buffer_1000_ddsp[274];

    auto g_z_0_0_0_yy_xy_0_z = buffer_1000_ddsp[275];

    auto g_z_0_0_0_yy_xz_0_x = buffer_1000_ddsp[276];

    auto g_z_0_0_0_yy_xz_0_y = buffer_1000_ddsp[277];

    auto g_z_0_0_0_yy_xz_0_z = buffer_1000_ddsp[278];

    auto g_z_0_0_0_yy_yy_0_x = buffer_1000_ddsp[279];

    auto g_z_0_0_0_yy_yy_0_y = buffer_1000_ddsp[280];

    auto g_z_0_0_0_yy_yy_0_z = buffer_1000_ddsp[281];

    auto g_z_0_0_0_yy_yz_0_x = buffer_1000_ddsp[282];

    auto g_z_0_0_0_yy_yz_0_y = buffer_1000_ddsp[283];

    auto g_z_0_0_0_yy_yz_0_z = buffer_1000_ddsp[284];

    auto g_z_0_0_0_yy_zz_0_x = buffer_1000_ddsp[285];

    auto g_z_0_0_0_yy_zz_0_y = buffer_1000_ddsp[286];

    auto g_z_0_0_0_yy_zz_0_z = buffer_1000_ddsp[287];

    auto g_z_0_0_0_yz_xx_0_x = buffer_1000_ddsp[288];

    auto g_z_0_0_0_yz_xx_0_y = buffer_1000_ddsp[289];

    auto g_z_0_0_0_yz_xx_0_z = buffer_1000_ddsp[290];

    auto g_z_0_0_0_yz_xy_0_x = buffer_1000_ddsp[291];

    auto g_z_0_0_0_yz_xy_0_y = buffer_1000_ddsp[292];

    auto g_z_0_0_0_yz_xy_0_z = buffer_1000_ddsp[293];

    auto g_z_0_0_0_yz_xz_0_x = buffer_1000_ddsp[294];

    auto g_z_0_0_0_yz_xz_0_y = buffer_1000_ddsp[295];

    auto g_z_0_0_0_yz_xz_0_z = buffer_1000_ddsp[296];

    auto g_z_0_0_0_yz_yy_0_x = buffer_1000_ddsp[297];

    auto g_z_0_0_0_yz_yy_0_y = buffer_1000_ddsp[298];

    auto g_z_0_0_0_yz_yy_0_z = buffer_1000_ddsp[299];

    auto g_z_0_0_0_yz_yz_0_x = buffer_1000_ddsp[300];

    auto g_z_0_0_0_yz_yz_0_y = buffer_1000_ddsp[301];

    auto g_z_0_0_0_yz_yz_0_z = buffer_1000_ddsp[302];

    auto g_z_0_0_0_yz_zz_0_x = buffer_1000_ddsp[303];

    auto g_z_0_0_0_yz_zz_0_y = buffer_1000_ddsp[304];

    auto g_z_0_0_0_yz_zz_0_z = buffer_1000_ddsp[305];

    auto g_z_0_0_0_zz_xx_0_x = buffer_1000_ddsp[306];

    auto g_z_0_0_0_zz_xx_0_y = buffer_1000_ddsp[307];

    auto g_z_0_0_0_zz_xx_0_z = buffer_1000_ddsp[308];

    auto g_z_0_0_0_zz_xy_0_x = buffer_1000_ddsp[309];

    auto g_z_0_0_0_zz_xy_0_y = buffer_1000_ddsp[310];

    auto g_z_0_0_0_zz_xy_0_z = buffer_1000_ddsp[311];

    auto g_z_0_0_0_zz_xz_0_x = buffer_1000_ddsp[312];

    auto g_z_0_0_0_zz_xz_0_y = buffer_1000_ddsp[313];

    auto g_z_0_0_0_zz_xz_0_z = buffer_1000_ddsp[314];

    auto g_z_0_0_0_zz_yy_0_x = buffer_1000_ddsp[315];

    auto g_z_0_0_0_zz_yy_0_y = buffer_1000_ddsp[316];

    auto g_z_0_0_0_zz_yy_0_z = buffer_1000_ddsp[317];

    auto g_z_0_0_0_zz_yz_0_x = buffer_1000_ddsp[318];

    auto g_z_0_0_0_zz_yz_0_y = buffer_1000_ddsp[319];

    auto g_z_0_0_0_zz_yz_0_z = buffer_1000_ddsp[320];

    auto g_z_0_0_0_zz_zz_0_x = buffer_1000_ddsp[321];

    auto g_z_0_0_0_zz_zz_0_y = buffer_1000_ddsp[322];

    auto g_z_0_0_0_zz_zz_0_z = buffer_1000_ddsp[323];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_0_x, g_x_0_0_0_xx_xx_0_y, g_x_0_0_0_xx_xx_0_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xxx_xx_0_x, g_xxx_xx_0_y, g_xxx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_0_x[i] = -2.0 * g_x_xx_0_x[i] + 2.0 * g_xxx_xx_0_x[i] * a_exp;

        g_x_0_0_0_xx_xx_0_y[i] = -2.0 * g_x_xx_0_y[i] + 2.0 * g_xxx_xx_0_y[i] * a_exp;

        g_x_0_0_0_xx_xx_0_z[i] = -2.0 * g_x_xx_0_z[i] + 2.0 * g_xxx_xx_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_0_x, g_x_0_0_0_xx_xy_0_y, g_x_0_0_0_xx_xy_0_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xxx_xy_0_x, g_xxx_xy_0_y, g_xxx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_0_x[i] = -2.0 * g_x_xy_0_x[i] + 2.0 * g_xxx_xy_0_x[i] * a_exp;

        g_x_0_0_0_xx_xy_0_y[i] = -2.0 * g_x_xy_0_y[i] + 2.0 * g_xxx_xy_0_y[i] * a_exp;

        g_x_0_0_0_xx_xy_0_z[i] = -2.0 * g_x_xy_0_z[i] + 2.0 * g_xxx_xy_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_0_x, g_x_0_0_0_xx_xz_0_y, g_x_0_0_0_xx_xz_0_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xxx_xz_0_x, g_xxx_xz_0_y, g_xxx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_0_x[i] = -2.0 * g_x_xz_0_x[i] + 2.0 * g_xxx_xz_0_x[i] * a_exp;

        g_x_0_0_0_xx_xz_0_y[i] = -2.0 * g_x_xz_0_y[i] + 2.0 * g_xxx_xz_0_y[i] * a_exp;

        g_x_0_0_0_xx_xz_0_z[i] = -2.0 * g_x_xz_0_z[i] + 2.0 * g_xxx_xz_0_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_0_x, g_x_0_0_0_xx_yy_0_y, g_x_0_0_0_xx_yy_0_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xxx_yy_0_x, g_xxx_yy_0_y, g_xxx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_0_x[i] = -2.0 * g_x_yy_0_x[i] + 2.0 * g_xxx_yy_0_x[i] * a_exp;

        g_x_0_0_0_xx_yy_0_y[i] = -2.0 * g_x_yy_0_y[i] + 2.0 * g_xxx_yy_0_y[i] * a_exp;

        g_x_0_0_0_xx_yy_0_z[i] = -2.0 * g_x_yy_0_z[i] + 2.0 * g_xxx_yy_0_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_0_x, g_x_0_0_0_xx_yz_0_y, g_x_0_0_0_xx_yz_0_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xxx_yz_0_x, g_xxx_yz_0_y, g_xxx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_0_x[i] = -2.0 * g_x_yz_0_x[i] + 2.0 * g_xxx_yz_0_x[i] * a_exp;

        g_x_0_0_0_xx_yz_0_y[i] = -2.0 * g_x_yz_0_y[i] + 2.0 * g_xxx_yz_0_y[i] * a_exp;

        g_x_0_0_0_xx_yz_0_z[i] = -2.0 * g_x_yz_0_z[i] + 2.0 * g_xxx_yz_0_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_0_x, g_x_0_0_0_xx_zz_0_y, g_x_0_0_0_xx_zz_0_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xxx_zz_0_x, g_xxx_zz_0_y, g_xxx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_0_x[i] = -2.0 * g_x_zz_0_x[i] + 2.0 * g_xxx_zz_0_x[i] * a_exp;

        g_x_0_0_0_xx_zz_0_y[i] = -2.0 * g_x_zz_0_y[i] + 2.0 * g_xxx_zz_0_y[i] * a_exp;

        g_x_0_0_0_xx_zz_0_z[i] = -2.0 * g_x_zz_0_z[i] + 2.0 * g_xxx_zz_0_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_0_x, g_x_0_0_0_xy_xx_0_y, g_x_0_0_0_xy_xx_0_z, g_xxy_xx_0_x, g_xxy_xx_0_y, g_xxy_xx_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_0_x[i] = -g_y_xx_0_x[i] + 2.0 * g_xxy_xx_0_x[i] * a_exp;

        g_x_0_0_0_xy_xx_0_y[i] = -g_y_xx_0_y[i] + 2.0 * g_xxy_xx_0_y[i] * a_exp;

        g_x_0_0_0_xy_xx_0_z[i] = -g_y_xx_0_z[i] + 2.0 * g_xxy_xx_0_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_0_x, g_x_0_0_0_xy_xy_0_y, g_x_0_0_0_xy_xy_0_z, g_xxy_xy_0_x, g_xxy_xy_0_y, g_xxy_xy_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_0_x[i] = -g_y_xy_0_x[i] + 2.0 * g_xxy_xy_0_x[i] * a_exp;

        g_x_0_0_0_xy_xy_0_y[i] = -g_y_xy_0_y[i] + 2.0 * g_xxy_xy_0_y[i] * a_exp;

        g_x_0_0_0_xy_xy_0_z[i] = -g_y_xy_0_z[i] + 2.0 * g_xxy_xy_0_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_0_x, g_x_0_0_0_xy_xz_0_y, g_x_0_0_0_xy_xz_0_z, g_xxy_xz_0_x, g_xxy_xz_0_y, g_xxy_xz_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_0_x[i] = -g_y_xz_0_x[i] + 2.0 * g_xxy_xz_0_x[i] * a_exp;

        g_x_0_0_0_xy_xz_0_y[i] = -g_y_xz_0_y[i] + 2.0 * g_xxy_xz_0_y[i] * a_exp;

        g_x_0_0_0_xy_xz_0_z[i] = -g_y_xz_0_z[i] + 2.0 * g_xxy_xz_0_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_0_x, g_x_0_0_0_xy_yy_0_y, g_x_0_0_0_xy_yy_0_z, g_xxy_yy_0_x, g_xxy_yy_0_y, g_xxy_yy_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_0_x[i] = -g_y_yy_0_x[i] + 2.0 * g_xxy_yy_0_x[i] * a_exp;

        g_x_0_0_0_xy_yy_0_y[i] = -g_y_yy_0_y[i] + 2.0 * g_xxy_yy_0_y[i] * a_exp;

        g_x_0_0_0_xy_yy_0_z[i] = -g_y_yy_0_z[i] + 2.0 * g_xxy_yy_0_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_0_x, g_x_0_0_0_xy_yz_0_y, g_x_0_0_0_xy_yz_0_z, g_xxy_yz_0_x, g_xxy_yz_0_y, g_xxy_yz_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_0_x[i] = -g_y_yz_0_x[i] + 2.0 * g_xxy_yz_0_x[i] * a_exp;

        g_x_0_0_0_xy_yz_0_y[i] = -g_y_yz_0_y[i] + 2.0 * g_xxy_yz_0_y[i] * a_exp;

        g_x_0_0_0_xy_yz_0_z[i] = -g_y_yz_0_z[i] + 2.0 * g_xxy_yz_0_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_0_x, g_x_0_0_0_xy_zz_0_y, g_x_0_0_0_xy_zz_0_z, g_xxy_zz_0_x, g_xxy_zz_0_y, g_xxy_zz_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_0_x[i] = -g_y_zz_0_x[i] + 2.0 * g_xxy_zz_0_x[i] * a_exp;

        g_x_0_0_0_xy_zz_0_y[i] = -g_y_zz_0_y[i] + 2.0 * g_xxy_zz_0_y[i] * a_exp;

        g_x_0_0_0_xy_zz_0_z[i] = -g_y_zz_0_z[i] + 2.0 * g_xxy_zz_0_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_0_x, g_x_0_0_0_xz_xx_0_y, g_x_0_0_0_xz_xx_0_z, g_xxz_xx_0_x, g_xxz_xx_0_y, g_xxz_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_0_x[i] = -g_z_xx_0_x[i] + 2.0 * g_xxz_xx_0_x[i] * a_exp;

        g_x_0_0_0_xz_xx_0_y[i] = -g_z_xx_0_y[i] + 2.0 * g_xxz_xx_0_y[i] * a_exp;

        g_x_0_0_0_xz_xx_0_z[i] = -g_z_xx_0_z[i] + 2.0 * g_xxz_xx_0_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_0_x, g_x_0_0_0_xz_xy_0_y, g_x_0_0_0_xz_xy_0_z, g_xxz_xy_0_x, g_xxz_xy_0_y, g_xxz_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_0_x[i] = -g_z_xy_0_x[i] + 2.0 * g_xxz_xy_0_x[i] * a_exp;

        g_x_0_0_0_xz_xy_0_y[i] = -g_z_xy_0_y[i] + 2.0 * g_xxz_xy_0_y[i] * a_exp;

        g_x_0_0_0_xz_xy_0_z[i] = -g_z_xy_0_z[i] + 2.0 * g_xxz_xy_0_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_0_x, g_x_0_0_0_xz_xz_0_y, g_x_0_0_0_xz_xz_0_z, g_xxz_xz_0_x, g_xxz_xz_0_y, g_xxz_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_0_x[i] = -g_z_xz_0_x[i] + 2.0 * g_xxz_xz_0_x[i] * a_exp;

        g_x_0_0_0_xz_xz_0_y[i] = -g_z_xz_0_y[i] + 2.0 * g_xxz_xz_0_y[i] * a_exp;

        g_x_0_0_0_xz_xz_0_z[i] = -g_z_xz_0_z[i] + 2.0 * g_xxz_xz_0_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_0_x, g_x_0_0_0_xz_yy_0_y, g_x_0_0_0_xz_yy_0_z, g_xxz_yy_0_x, g_xxz_yy_0_y, g_xxz_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_0_x[i] = -g_z_yy_0_x[i] + 2.0 * g_xxz_yy_0_x[i] * a_exp;

        g_x_0_0_0_xz_yy_0_y[i] = -g_z_yy_0_y[i] + 2.0 * g_xxz_yy_0_y[i] * a_exp;

        g_x_0_0_0_xz_yy_0_z[i] = -g_z_yy_0_z[i] + 2.0 * g_xxz_yy_0_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_0_x, g_x_0_0_0_xz_yz_0_y, g_x_0_0_0_xz_yz_0_z, g_xxz_yz_0_x, g_xxz_yz_0_y, g_xxz_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_0_x[i] = -g_z_yz_0_x[i] + 2.0 * g_xxz_yz_0_x[i] * a_exp;

        g_x_0_0_0_xz_yz_0_y[i] = -g_z_yz_0_y[i] + 2.0 * g_xxz_yz_0_y[i] * a_exp;

        g_x_0_0_0_xz_yz_0_z[i] = -g_z_yz_0_z[i] + 2.0 * g_xxz_yz_0_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_0_x, g_x_0_0_0_xz_zz_0_y, g_x_0_0_0_xz_zz_0_z, g_xxz_zz_0_x, g_xxz_zz_0_y, g_xxz_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_0_x[i] = -g_z_zz_0_x[i] + 2.0 * g_xxz_zz_0_x[i] * a_exp;

        g_x_0_0_0_xz_zz_0_y[i] = -g_z_zz_0_y[i] + 2.0 * g_xxz_zz_0_y[i] * a_exp;

        g_x_0_0_0_xz_zz_0_z[i] = -g_z_zz_0_z[i] + 2.0 * g_xxz_zz_0_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_0_x, g_x_0_0_0_yy_xx_0_y, g_x_0_0_0_yy_xx_0_z, g_xyy_xx_0_x, g_xyy_xx_0_y, g_xyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_0_x[i] = 2.0 * g_xyy_xx_0_x[i] * a_exp;

        g_x_0_0_0_yy_xx_0_y[i] = 2.0 * g_xyy_xx_0_y[i] * a_exp;

        g_x_0_0_0_yy_xx_0_z[i] = 2.0 * g_xyy_xx_0_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_0_x, g_x_0_0_0_yy_xy_0_y, g_x_0_0_0_yy_xy_0_z, g_xyy_xy_0_x, g_xyy_xy_0_y, g_xyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_0_x[i] = 2.0 * g_xyy_xy_0_x[i] * a_exp;

        g_x_0_0_0_yy_xy_0_y[i] = 2.0 * g_xyy_xy_0_y[i] * a_exp;

        g_x_0_0_0_yy_xy_0_z[i] = 2.0 * g_xyy_xy_0_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_0_x, g_x_0_0_0_yy_xz_0_y, g_x_0_0_0_yy_xz_0_z, g_xyy_xz_0_x, g_xyy_xz_0_y, g_xyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_0_x[i] = 2.0 * g_xyy_xz_0_x[i] * a_exp;

        g_x_0_0_0_yy_xz_0_y[i] = 2.0 * g_xyy_xz_0_y[i] * a_exp;

        g_x_0_0_0_yy_xz_0_z[i] = 2.0 * g_xyy_xz_0_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_0_x, g_x_0_0_0_yy_yy_0_y, g_x_0_0_0_yy_yy_0_z, g_xyy_yy_0_x, g_xyy_yy_0_y, g_xyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_0_x[i] = 2.0 * g_xyy_yy_0_x[i] * a_exp;

        g_x_0_0_0_yy_yy_0_y[i] = 2.0 * g_xyy_yy_0_y[i] * a_exp;

        g_x_0_0_0_yy_yy_0_z[i] = 2.0 * g_xyy_yy_0_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_0_x, g_x_0_0_0_yy_yz_0_y, g_x_0_0_0_yy_yz_0_z, g_xyy_yz_0_x, g_xyy_yz_0_y, g_xyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_0_x[i] = 2.0 * g_xyy_yz_0_x[i] * a_exp;

        g_x_0_0_0_yy_yz_0_y[i] = 2.0 * g_xyy_yz_0_y[i] * a_exp;

        g_x_0_0_0_yy_yz_0_z[i] = 2.0 * g_xyy_yz_0_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_0_x, g_x_0_0_0_yy_zz_0_y, g_x_0_0_0_yy_zz_0_z, g_xyy_zz_0_x, g_xyy_zz_0_y, g_xyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_0_x[i] = 2.0 * g_xyy_zz_0_x[i] * a_exp;

        g_x_0_0_0_yy_zz_0_y[i] = 2.0 * g_xyy_zz_0_y[i] * a_exp;

        g_x_0_0_0_yy_zz_0_z[i] = 2.0 * g_xyy_zz_0_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_0_x, g_x_0_0_0_yz_xx_0_y, g_x_0_0_0_yz_xx_0_z, g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_0_x[i] = 2.0 * g_xyz_xx_0_x[i] * a_exp;

        g_x_0_0_0_yz_xx_0_y[i] = 2.0 * g_xyz_xx_0_y[i] * a_exp;

        g_x_0_0_0_yz_xx_0_z[i] = 2.0 * g_xyz_xx_0_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_0_x, g_x_0_0_0_yz_xy_0_y, g_x_0_0_0_yz_xy_0_z, g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_0_x[i] = 2.0 * g_xyz_xy_0_x[i] * a_exp;

        g_x_0_0_0_yz_xy_0_y[i] = 2.0 * g_xyz_xy_0_y[i] * a_exp;

        g_x_0_0_0_yz_xy_0_z[i] = 2.0 * g_xyz_xy_0_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_0_x, g_x_0_0_0_yz_xz_0_y, g_x_0_0_0_yz_xz_0_z, g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_0_x[i] = 2.0 * g_xyz_xz_0_x[i] * a_exp;

        g_x_0_0_0_yz_xz_0_y[i] = 2.0 * g_xyz_xz_0_y[i] * a_exp;

        g_x_0_0_0_yz_xz_0_z[i] = 2.0 * g_xyz_xz_0_z[i] * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_0_x, g_x_0_0_0_yz_yy_0_y, g_x_0_0_0_yz_yy_0_z, g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_0_x[i] = 2.0 * g_xyz_yy_0_x[i] * a_exp;

        g_x_0_0_0_yz_yy_0_y[i] = 2.0 * g_xyz_yy_0_y[i] * a_exp;

        g_x_0_0_0_yz_yy_0_z[i] = 2.0 * g_xyz_yy_0_z[i] * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_0_x, g_x_0_0_0_yz_yz_0_y, g_x_0_0_0_yz_yz_0_z, g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_0_x[i] = 2.0 * g_xyz_yz_0_x[i] * a_exp;

        g_x_0_0_0_yz_yz_0_y[i] = 2.0 * g_xyz_yz_0_y[i] * a_exp;

        g_x_0_0_0_yz_yz_0_z[i] = 2.0 * g_xyz_yz_0_z[i] * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_0_x, g_x_0_0_0_yz_zz_0_y, g_x_0_0_0_yz_zz_0_z, g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_0_x[i] = 2.0 * g_xyz_zz_0_x[i] * a_exp;

        g_x_0_0_0_yz_zz_0_y[i] = 2.0 * g_xyz_zz_0_y[i] * a_exp;

        g_x_0_0_0_yz_zz_0_z[i] = 2.0 * g_xyz_zz_0_z[i] * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_0_x, g_x_0_0_0_zz_xx_0_y, g_x_0_0_0_zz_xx_0_z, g_xzz_xx_0_x, g_xzz_xx_0_y, g_xzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_0_x[i] = 2.0 * g_xzz_xx_0_x[i] * a_exp;

        g_x_0_0_0_zz_xx_0_y[i] = 2.0 * g_xzz_xx_0_y[i] * a_exp;

        g_x_0_0_0_zz_xx_0_z[i] = 2.0 * g_xzz_xx_0_z[i] * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_0_x, g_x_0_0_0_zz_xy_0_y, g_x_0_0_0_zz_xy_0_z, g_xzz_xy_0_x, g_xzz_xy_0_y, g_xzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_0_x[i] = 2.0 * g_xzz_xy_0_x[i] * a_exp;

        g_x_0_0_0_zz_xy_0_y[i] = 2.0 * g_xzz_xy_0_y[i] * a_exp;

        g_x_0_0_0_zz_xy_0_z[i] = 2.0 * g_xzz_xy_0_z[i] * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_0_x, g_x_0_0_0_zz_xz_0_y, g_x_0_0_0_zz_xz_0_z, g_xzz_xz_0_x, g_xzz_xz_0_y, g_xzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_0_x[i] = 2.0 * g_xzz_xz_0_x[i] * a_exp;

        g_x_0_0_0_zz_xz_0_y[i] = 2.0 * g_xzz_xz_0_y[i] * a_exp;

        g_x_0_0_0_zz_xz_0_z[i] = 2.0 * g_xzz_xz_0_z[i] * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_0_x, g_x_0_0_0_zz_yy_0_y, g_x_0_0_0_zz_yy_0_z, g_xzz_yy_0_x, g_xzz_yy_0_y, g_xzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_0_x[i] = 2.0 * g_xzz_yy_0_x[i] * a_exp;

        g_x_0_0_0_zz_yy_0_y[i] = 2.0 * g_xzz_yy_0_y[i] * a_exp;

        g_x_0_0_0_zz_yy_0_z[i] = 2.0 * g_xzz_yy_0_z[i] * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_0_x, g_x_0_0_0_zz_yz_0_y, g_x_0_0_0_zz_yz_0_z, g_xzz_yz_0_x, g_xzz_yz_0_y, g_xzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_0_x[i] = 2.0 * g_xzz_yz_0_x[i] * a_exp;

        g_x_0_0_0_zz_yz_0_y[i] = 2.0 * g_xzz_yz_0_y[i] * a_exp;

        g_x_0_0_0_zz_yz_0_z[i] = 2.0 * g_xzz_yz_0_z[i] * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_0_x, g_x_0_0_0_zz_zz_0_y, g_x_0_0_0_zz_zz_0_z, g_xzz_zz_0_x, g_xzz_zz_0_y, g_xzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_0_x[i] = 2.0 * g_xzz_zz_0_x[i] * a_exp;

        g_x_0_0_0_zz_zz_0_y[i] = 2.0 * g_xzz_zz_0_y[i] * a_exp;

        g_x_0_0_0_zz_zz_0_z[i] = 2.0 * g_xzz_zz_0_z[i] * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xxy_xx_0_x, g_xxy_xx_0_y, g_xxy_xx_0_z, g_y_0_0_0_xx_xx_0_x, g_y_0_0_0_xx_xx_0_y, g_y_0_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_0_x[i] = 2.0 * g_xxy_xx_0_x[i] * a_exp;

        g_y_0_0_0_xx_xx_0_y[i] = 2.0 * g_xxy_xx_0_y[i] * a_exp;

        g_y_0_0_0_xx_xx_0_z[i] = 2.0 * g_xxy_xx_0_z[i] * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xxy_xy_0_x, g_xxy_xy_0_y, g_xxy_xy_0_z, g_y_0_0_0_xx_xy_0_x, g_y_0_0_0_xx_xy_0_y, g_y_0_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_0_x[i] = 2.0 * g_xxy_xy_0_x[i] * a_exp;

        g_y_0_0_0_xx_xy_0_y[i] = 2.0 * g_xxy_xy_0_y[i] * a_exp;

        g_y_0_0_0_xx_xy_0_z[i] = 2.0 * g_xxy_xy_0_z[i] * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xxy_xz_0_x, g_xxy_xz_0_y, g_xxy_xz_0_z, g_y_0_0_0_xx_xz_0_x, g_y_0_0_0_xx_xz_0_y, g_y_0_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_0_x[i] = 2.0 * g_xxy_xz_0_x[i] * a_exp;

        g_y_0_0_0_xx_xz_0_y[i] = 2.0 * g_xxy_xz_0_y[i] * a_exp;

        g_y_0_0_0_xx_xz_0_z[i] = 2.0 * g_xxy_xz_0_z[i] * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xxy_yy_0_x, g_xxy_yy_0_y, g_xxy_yy_0_z, g_y_0_0_0_xx_yy_0_x, g_y_0_0_0_xx_yy_0_y, g_y_0_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_0_x[i] = 2.0 * g_xxy_yy_0_x[i] * a_exp;

        g_y_0_0_0_xx_yy_0_y[i] = 2.0 * g_xxy_yy_0_y[i] * a_exp;

        g_y_0_0_0_xx_yy_0_z[i] = 2.0 * g_xxy_yy_0_z[i] * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xxy_yz_0_x, g_xxy_yz_0_y, g_xxy_yz_0_z, g_y_0_0_0_xx_yz_0_x, g_y_0_0_0_xx_yz_0_y, g_y_0_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_0_x[i] = 2.0 * g_xxy_yz_0_x[i] * a_exp;

        g_y_0_0_0_xx_yz_0_y[i] = 2.0 * g_xxy_yz_0_y[i] * a_exp;

        g_y_0_0_0_xx_yz_0_z[i] = 2.0 * g_xxy_yz_0_z[i] * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xxy_zz_0_x, g_xxy_zz_0_y, g_xxy_zz_0_z, g_y_0_0_0_xx_zz_0_x, g_y_0_0_0_xx_zz_0_y, g_y_0_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_0_x[i] = 2.0 * g_xxy_zz_0_x[i] * a_exp;

        g_y_0_0_0_xx_zz_0_y[i] = 2.0 * g_xxy_zz_0_y[i] * a_exp;

        g_y_0_0_0_xx_zz_0_z[i] = 2.0 * g_xxy_zz_0_z[i] * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xyy_xx_0_x, g_xyy_xx_0_y, g_xyy_xx_0_z, g_y_0_0_0_xy_xx_0_x, g_y_0_0_0_xy_xx_0_y, g_y_0_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_0_x[i] = -g_x_xx_0_x[i] + 2.0 * g_xyy_xx_0_x[i] * a_exp;

        g_y_0_0_0_xy_xx_0_y[i] = -g_x_xx_0_y[i] + 2.0 * g_xyy_xx_0_y[i] * a_exp;

        g_y_0_0_0_xy_xx_0_z[i] = -g_x_xx_0_z[i] + 2.0 * g_xyy_xx_0_z[i] * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xyy_xy_0_x, g_xyy_xy_0_y, g_xyy_xy_0_z, g_y_0_0_0_xy_xy_0_x, g_y_0_0_0_xy_xy_0_y, g_y_0_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_0_x[i] = -g_x_xy_0_x[i] + 2.0 * g_xyy_xy_0_x[i] * a_exp;

        g_y_0_0_0_xy_xy_0_y[i] = -g_x_xy_0_y[i] + 2.0 * g_xyy_xy_0_y[i] * a_exp;

        g_y_0_0_0_xy_xy_0_z[i] = -g_x_xy_0_z[i] + 2.0 * g_xyy_xy_0_z[i] * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xyy_xz_0_x, g_xyy_xz_0_y, g_xyy_xz_0_z, g_y_0_0_0_xy_xz_0_x, g_y_0_0_0_xy_xz_0_y, g_y_0_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_0_x[i] = -g_x_xz_0_x[i] + 2.0 * g_xyy_xz_0_x[i] * a_exp;

        g_y_0_0_0_xy_xz_0_y[i] = -g_x_xz_0_y[i] + 2.0 * g_xyy_xz_0_y[i] * a_exp;

        g_y_0_0_0_xy_xz_0_z[i] = -g_x_xz_0_z[i] + 2.0 * g_xyy_xz_0_z[i] * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xyy_yy_0_x, g_xyy_yy_0_y, g_xyy_yy_0_z, g_y_0_0_0_xy_yy_0_x, g_y_0_0_0_xy_yy_0_y, g_y_0_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_0_x[i] = -g_x_yy_0_x[i] + 2.0 * g_xyy_yy_0_x[i] * a_exp;

        g_y_0_0_0_xy_yy_0_y[i] = -g_x_yy_0_y[i] + 2.0 * g_xyy_yy_0_y[i] * a_exp;

        g_y_0_0_0_xy_yy_0_z[i] = -g_x_yy_0_z[i] + 2.0 * g_xyy_yy_0_z[i] * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xyy_yz_0_x, g_xyy_yz_0_y, g_xyy_yz_0_z, g_y_0_0_0_xy_yz_0_x, g_y_0_0_0_xy_yz_0_y, g_y_0_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_0_x[i] = -g_x_yz_0_x[i] + 2.0 * g_xyy_yz_0_x[i] * a_exp;

        g_y_0_0_0_xy_yz_0_y[i] = -g_x_yz_0_y[i] + 2.0 * g_xyy_yz_0_y[i] * a_exp;

        g_y_0_0_0_xy_yz_0_z[i] = -g_x_yz_0_z[i] + 2.0 * g_xyy_yz_0_z[i] * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xyy_zz_0_x, g_xyy_zz_0_y, g_xyy_zz_0_z, g_y_0_0_0_xy_zz_0_x, g_y_0_0_0_xy_zz_0_y, g_y_0_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_0_x[i] = -g_x_zz_0_x[i] + 2.0 * g_xyy_zz_0_x[i] * a_exp;

        g_y_0_0_0_xy_zz_0_y[i] = -g_x_zz_0_y[i] + 2.0 * g_xyy_zz_0_y[i] * a_exp;

        g_y_0_0_0_xy_zz_0_z[i] = -g_x_zz_0_z[i] + 2.0 * g_xyy_zz_0_z[i] * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z, g_y_0_0_0_xz_xx_0_x, g_y_0_0_0_xz_xx_0_y, g_y_0_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_0_x[i] = 2.0 * g_xyz_xx_0_x[i] * a_exp;

        g_y_0_0_0_xz_xx_0_y[i] = 2.0 * g_xyz_xx_0_y[i] * a_exp;

        g_y_0_0_0_xz_xx_0_z[i] = 2.0 * g_xyz_xx_0_z[i] * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z, g_y_0_0_0_xz_xy_0_x, g_y_0_0_0_xz_xy_0_y, g_y_0_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_0_x[i] = 2.0 * g_xyz_xy_0_x[i] * a_exp;

        g_y_0_0_0_xz_xy_0_y[i] = 2.0 * g_xyz_xy_0_y[i] * a_exp;

        g_y_0_0_0_xz_xy_0_z[i] = 2.0 * g_xyz_xy_0_z[i] * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z, g_y_0_0_0_xz_xz_0_x, g_y_0_0_0_xz_xz_0_y, g_y_0_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_0_x[i] = 2.0 * g_xyz_xz_0_x[i] * a_exp;

        g_y_0_0_0_xz_xz_0_y[i] = 2.0 * g_xyz_xz_0_y[i] * a_exp;

        g_y_0_0_0_xz_xz_0_z[i] = 2.0 * g_xyz_xz_0_z[i] * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z, g_y_0_0_0_xz_yy_0_x, g_y_0_0_0_xz_yy_0_y, g_y_0_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_0_x[i] = 2.0 * g_xyz_yy_0_x[i] * a_exp;

        g_y_0_0_0_xz_yy_0_y[i] = 2.0 * g_xyz_yy_0_y[i] * a_exp;

        g_y_0_0_0_xz_yy_0_z[i] = 2.0 * g_xyz_yy_0_z[i] * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z, g_y_0_0_0_xz_yz_0_x, g_y_0_0_0_xz_yz_0_y, g_y_0_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_0_x[i] = 2.0 * g_xyz_yz_0_x[i] * a_exp;

        g_y_0_0_0_xz_yz_0_y[i] = 2.0 * g_xyz_yz_0_y[i] * a_exp;

        g_y_0_0_0_xz_yz_0_z[i] = 2.0 * g_xyz_yz_0_z[i] * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z, g_y_0_0_0_xz_zz_0_x, g_y_0_0_0_xz_zz_0_y, g_y_0_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_0_x[i] = 2.0 * g_xyz_zz_0_x[i] * a_exp;

        g_y_0_0_0_xz_zz_0_y[i] = 2.0 * g_xyz_zz_0_y[i] * a_exp;

        g_y_0_0_0_xz_zz_0_z[i] = 2.0 * g_xyz_zz_0_z[i] * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_0_x, g_y_0_0_0_yy_xx_0_y, g_y_0_0_0_yy_xx_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_yyy_xx_0_x, g_yyy_xx_0_y, g_yyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_0_x[i] = -2.0 * g_y_xx_0_x[i] + 2.0 * g_yyy_xx_0_x[i] * a_exp;

        g_y_0_0_0_yy_xx_0_y[i] = -2.0 * g_y_xx_0_y[i] + 2.0 * g_yyy_xx_0_y[i] * a_exp;

        g_y_0_0_0_yy_xx_0_z[i] = -2.0 * g_y_xx_0_z[i] + 2.0 * g_yyy_xx_0_z[i] * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_0_x, g_y_0_0_0_yy_xy_0_y, g_y_0_0_0_yy_xy_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_yyy_xy_0_x, g_yyy_xy_0_y, g_yyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_0_x[i] = -2.0 * g_y_xy_0_x[i] + 2.0 * g_yyy_xy_0_x[i] * a_exp;

        g_y_0_0_0_yy_xy_0_y[i] = -2.0 * g_y_xy_0_y[i] + 2.0 * g_yyy_xy_0_y[i] * a_exp;

        g_y_0_0_0_yy_xy_0_z[i] = -2.0 * g_y_xy_0_z[i] + 2.0 * g_yyy_xy_0_z[i] * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_0_x, g_y_0_0_0_yy_xz_0_y, g_y_0_0_0_yy_xz_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_yyy_xz_0_x, g_yyy_xz_0_y, g_yyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_0_x[i] = -2.0 * g_y_xz_0_x[i] + 2.0 * g_yyy_xz_0_x[i] * a_exp;

        g_y_0_0_0_yy_xz_0_y[i] = -2.0 * g_y_xz_0_y[i] + 2.0 * g_yyy_xz_0_y[i] * a_exp;

        g_y_0_0_0_yy_xz_0_z[i] = -2.0 * g_y_xz_0_z[i] + 2.0 * g_yyy_xz_0_z[i] * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_0_x, g_y_0_0_0_yy_yy_0_y, g_y_0_0_0_yy_yy_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_yyy_yy_0_x, g_yyy_yy_0_y, g_yyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_0_x[i] = -2.0 * g_y_yy_0_x[i] + 2.0 * g_yyy_yy_0_x[i] * a_exp;

        g_y_0_0_0_yy_yy_0_y[i] = -2.0 * g_y_yy_0_y[i] + 2.0 * g_yyy_yy_0_y[i] * a_exp;

        g_y_0_0_0_yy_yy_0_z[i] = -2.0 * g_y_yy_0_z[i] + 2.0 * g_yyy_yy_0_z[i] * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_0_x, g_y_0_0_0_yy_yz_0_y, g_y_0_0_0_yy_yz_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_yyy_yz_0_x, g_yyy_yz_0_y, g_yyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_0_x[i] = -2.0 * g_y_yz_0_x[i] + 2.0 * g_yyy_yz_0_x[i] * a_exp;

        g_y_0_0_0_yy_yz_0_y[i] = -2.0 * g_y_yz_0_y[i] + 2.0 * g_yyy_yz_0_y[i] * a_exp;

        g_y_0_0_0_yy_yz_0_z[i] = -2.0 * g_y_yz_0_z[i] + 2.0 * g_yyy_yz_0_z[i] * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_0_x, g_y_0_0_0_yy_zz_0_y, g_y_0_0_0_yy_zz_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_yyy_zz_0_x, g_yyy_zz_0_y, g_yyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_0_x[i] = -2.0 * g_y_zz_0_x[i] + 2.0 * g_yyy_zz_0_x[i] * a_exp;

        g_y_0_0_0_yy_zz_0_y[i] = -2.0 * g_y_zz_0_y[i] + 2.0 * g_yyy_zz_0_y[i] * a_exp;

        g_y_0_0_0_yy_zz_0_z[i] = -2.0 * g_y_zz_0_z[i] + 2.0 * g_yyy_zz_0_z[i] * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_0_x, g_y_0_0_0_yz_xx_0_y, g_y_0_0_0_yz_xx_0_z, g_yyz_xx_0_x, g_yyz_xx_0_y, g_yyz_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_0_x[i] = -g_z_xx_0_x[i] + 2.0 * g_yyz_xx_0_x[i] * a_exp;

        g_y_0_0_0_yz_xx_0_y[i] = -g_z_xx_0_y[i] + 2.0 * g_yyz_xx_0_y[i] * a_exp;

        g_y_0_0_0_yz_xx_0_z[i] = -g_z_xx_0_z[i] + 2.0 * g_yyz_xx_0_z[i] * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_0_x, g_y_0_0_0_yz_xy_0_y, g_y_0_0_0_yz_xy_0_z, g_yyz_xy_0_x, g_yyz_xy_0_y, g_yyz_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_0_x[i] = -g_z_xy_0_x[i] + 2.0 * g_yyz_xy_0_x[i] * a_exp;

        g_y_0_0_0_yz_xy_0_y[i] = -g_z_xy_0_y[i] + 2.0 * g_yyz_xy_0_y[i] * a_exp;

        g_y_0_0_0_yz_xy_0_z[i] = -g_z_xy_0_z[i] + 2.0 * g_yyz_xy_0_z[i] * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_0_x, g_y_0_0_0_yz_xz_0_y, g_y_0_0_0_yz_xz_0_z, g_yyz_xz_0_x, g_yyz_xz_0_y, g_yyz_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_0_x[i] = -g_z_xz_0_x[i] + 2.0 * g_yyz_xz_0_x[i] * a_exp;

        g_y_0_0_0_yz_xz_0_y[i] = -g_z_xz_0_y[i] + 2.0 * g_yyz_xz_0_y[i] * a_exp;

        g_y_0_0_0_yz_xz_0_z[i] = -g_z_xz_0_z[i] + 2.0 * g_yyz_xz_0_z[i] * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_0_x, g_y_0_0_0_yz_yy_0_y, g_y_0_0_0_yz_yy_0_z, g_yyz_yy_0_x, g_yyz_yy_0_y, g_yyz_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_0_x[i] = -g_z_yy_0_x[i] + 2.0 * g_yyz_yy_0_x[i] * a_exp;

        g_y_0_0_0_yz_yy_0_y[i] = -g_z_yy_0_y[i] + 2.0 * g_yyz_yy_0_y[i] * a_exp;

        g_y_0_0_0_yz_yy_0_z[i] = -g_z_yy_0_z[i] + 2.0 * g_yyz_yy_0_z[i] * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_0_x, g_y_0_0_0_yz_yz_0_y, g_y_0_0_0_yz_yz_0_z, g_yyz_yz_0_x, g_yyz_yz_0_y, g_yyz_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_0_x[i] = -g_z_yz_0_x[i] + 2.0 * g_yyz_yz_0_x[i] * a_exp;

        g_y_0_0_0_yz_yz_0_y[i] = -g_z_yz_0_y[i] + 2.0 * g_yyz_yz_0_y[i] * a_exp;

        g_y_0_0_0_yz_yz_0_z[i] = -g_z_yz_0_z[i] + 2.0 * g_yyz_yz_0_z[i] * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_0_x, g_y_0_0_0_yz_zz_0_y, g_y_0_0_0_yz_zz_0_z, g_yyz_zz_0_x, g_yyz_zz_0_y, g_yyz_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_0_x[i] = -g_z_zz_0_x[i] + 2.0 * g_yyz_zz_0_x[i] * a_exp;

        g_y_0_0_0_yz_zz_0_y[i] = -g_z_zz_0_y[i] + 2.0 * g_yyz_zz_0_y[i] * a_exp;

        g_y_0_0_0_yz_zz_0_z[i] = -g_z_zz_0_z[i] + 2.0 * g_yyz_zz_0_z[i] * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_0_x, g_y_0_0_0_zz_xx_0_y, g_y_0_0_0_zz_xx_0_z, g_yzz_xx_0_x, g_yzz_xx_0_y, g_yzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_0_x[i] = 2.0 * g_yzz_xx_0_x[i] * a_exp;

        g_y_0_0_0_zz_xx_0_y[i] = 2.0 * g_yzz_xx_0_y[i] * a_exp;

        g_y_0_0_0_zz_xx_0_z[i] = 2.0 * g_yzz_xx_0_z[i] * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_0_x, g_y_0_0_0_zz_xy_0_y, g_y_0_0_0_zz_xy_0_z, g_yzz_xy_0_x, g_yzz_xy_0_y, g_yzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_0_x[i] = 2.0 * g_yzz_xy_0_x[i] * a_exp;

        g_y_0_0_0_zz_xy_0_y[i] = 2.0 * g_yzz_xy_0_y[i] * a_exp;

        g_y_0_0_0_zz_xy_0_z[i] = 2.0 * g_yzz_xy_0_z[i] * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_0_x, g_y_0_0_0_zz_xz_0_y, g_y_0_0_0_zz_xz_0_z, g_yzz_xz_0_x, g_yzz_xz_0_y, g_yzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_0_x[i] = 2.0 * g_yzz_xz_0_x[i] * a_exp;

        g_y_0_0_0_zz_xz_0_y[i] = 2.0 * g_yzz_xz_0_y[i] * a_exp;

        g_y_0_0_0_zz_xz_0_z[i] = 2.0 * g_yzz_xz_0_z[i] * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_0_x, g_y_0_0_0_zz_yy_0_y, g_y_0_0_0_zz_yy_0_z, g_yzz_yy_0_x, g_yzz_yy_0_y, g_yzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_0_x[i] = 2.0 * g_yzz_yy_0_x[i] * a_exp;

        g_y_0_0_0_zz_yy_0_y[i] = 2.0 * g_yzz_yy_0_y[i] * a_exp;

        g_y_0_0_0_zz_yy_0_z[i] = 2.0 * g_yzz_yy_0_z[i] * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_0_x, g_y_0_0_0_zz_yz_0_y, g_y_0_0_0_zz_yz_0_z, g_yzz_yz_0_x, g_yzz_yz_0_y, g_yzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_0_x[i] = 2.0 * g_yzz_yz_0_x[i] * a_exp;

        g_y_0_0_0_zz_yz_0_y[i] = 2.0 * g_yzz_yz_0_y[i] * a_exp;

        g_y_0_0_0_zz_yz_0_z[i] = 2.0 * g_yzz_yz_0_z[i] * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_0_x, g_y_0_0_0_zz_zz_0_y, g_y_0_0_0_zz_zz_0_z, g_yzz_zz_0_x, g_yzz_zz_0_y, g_yzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_0_x[i] = 2.0 * g_yzz_zz_0_x[i] * a_exp;

        g_y_0_0_0_zz_zz_0_y[i] = 2.0 * g_yzz_zz_0_y[i] * a_exp;

        g_y_0_0_0_zz_zz_0_z[i] = 2.0 * g_yzz_zz_0_z[i] * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xxz_xx_0_x, g_xxz_xx_0_y, g_xxz_xx_0_z, g_z_0_0_0_xx_xx_0_x, g_z_0_0_0_xx_xx_0_y, g_z_0_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_0_x[i] = 2.0 * g_xxz_xx_0_x[i] * a_exp;

        g_z_0_0_0_xx_xx_0_y[i] = 2.0 * g_xxz_xx_0_y[i] * a_exp;

        g_z_0_0_0_xx_xx_0_z[i] = 2.0 * g_xxz_xx_0_z[i] * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xxz_xy_0_x, g_xxz_xy_0_y, g_xxz_xy_0_z, g_z_0_0_0_xx_xy_0_x, g_z_0_0_0_xx_xy_0_y, g_z_0_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_0_x[i] = 2.0 * g_xxz_xy_0_x[i] * a_exp;

        g_z_0_0_0_xx_xy_0_y[i] = 2.0 * g_xxz_xy_0_y[i] * a_exp;

        g_z_0_0_0_xx_xy_0_z[i] = 2.0 * g_xxz_xy_0_z[i] * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xxz_xz_0_x, g_xxz_xz_0_y, g_xxz_xz_0_z, g_z_0_0_0_xx_xz_0_x, g_z_0_0_0_xx_xz_0_y, g_z_0_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_0_x[i] = 2.0 * g_xxz_xz_0_x[i] * a_exp;

        g_z_0_0_0_xx_xz_0_y[i] = 2.0 * g_xxz_xz_0_y[i] * a_exp;

        g_z_0_0_0_xx_xz_0_z[i] = 2.0 * g_xxz_xz_0_z[i] * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xxz_yy_0_x, g_xxz_yy_0_y, g_xxz_yy_0_z, g_z_0_0_0_xx_yy_0_x, g_z_0_0_0_xx_yy_0_y, g_z_0_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_0_x[i] = 2.0 * g_xxz_yy_0_x[i] * a_exp;

        g_z_0_0_0_xx_yy_0_y[i] = 2.0 * g_xxz_yy_0_y[i] * a_exp;

        g_z_0_0_0_xx_yy_0_z[i] = 2.0 * g_xxz_yy_0_z[i] * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xxz_yz_0_x, g_xxz_yz_0_y, g_xxz_yz_0_z, g_z_0_0_0_xx_yz_0_x, g_z_0_0_0_xx_yz_0_y, g_z_0_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_0_x[i] = 2.0 * g_xxz_yz_0_x[i] * a_exp;

        g_z_0_0_0_xx_yz_0_y[i] = 2.0 * g_xxz_yz_0_y[i] * a_exp;

        g_z_0_0_0_xx_yz_0_z[i] = 2.0 * g_xxz_yz_0_z[i] * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xxz_zz_0_x, g_xxz_zz_0_y, g_xxz_zz_0_z, g_z_0_0_0_xx_zz_0_x, g_z_0_0_0_xx_zz_0_y, g_z_0_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_0_x[i] = 2.0 * g_xxz_zz_0_x[i] * a_exp;

        g_z_0_0_0_xx_zz_0_y[i] = 2.0 * g_xxz_zz_0_y[i] * a_exp;

        g_z_0_0_0_xx_zz_0_z[i] = 2.0 * g_xxz_zz_0_z[i] * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_xyz_xx_0_x, g_xyz_xx_0_y, g_xyz_xx_0_z, g_z_0_0_0_xy_xx_0_x, g_z_0_0_0_xy_xx_0_y, g_z_0_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_0_x[i] = 2.0 * g_xyz_xx_0_x[i] * a_exp;

        g_z_0_0_0_xy_xx_0_y[i] = 2.0 * g_xyz_xx_0_y[i] * a_exp;

        g_z_0_0_0_xy_xx_0_z[i] = 2.0 * g_xyz_xx_0_z[i] * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_xyz_xy_0_x, g_xyz_xy_0_y, g_xyz_xy_0_z, g_z_0_0_0_xy_xy_0_x, g_z_0_0_0_xy_xy_0_y, g_z_0_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_0_x[i] = 2.0 * g_xyz_xy_0_x[i] * a_exp;

        g_z_0_0_0_xy_xy_0_y[i] = 2.0 * g_xyz_xy_0_y[i] * a_exp;

        g_z_0_0_0_xy_xy_0_z[i] = 2.0 * g_xyz_xy_0_z[i] * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_xyz_xz_0_x, g_xyz_xz_0_y, g_xyz_xz_0_z, g_z_0_0_0_xy_xz_0_x, g_z_0_0_0_xy_xz_0_y, g_z_0_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_0_x[i] = 2.0 * g_xyz_xz_0_x[i] * a_exp;

        g_z_0_0_0_xy_xz_0_y[i] = 2.0 * g_xyz_xz_0_y[i] * a_exp;

        g_z_0_0_0_xy_xz_0_z[i] = 2.0 * g_xyz_xz_0_z[i] * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_xyz_yy_0_x, g_xyz_yy_0_y, g_xyz_yy_0_z, g_z_0_0_0_xy_yy_0_x, g_z_0_0_0_xy_yy_0_y, g_z_0_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_0_x[i] = 2.0 * g_xyz_yy_0_x[i] * a_exp;

        g_z_0_0_0_xy_yy_0_y[i] = 2.0 * g_xyz_yy_0_y[i] * a_exp;

        g_z_0_0_0_xy_yy_0_z[i] = 2.0 * g_xyz_yy_0_z[i] * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_xyz_yz_0_x, g_xyz_yz_0_y, g_xyz_yz_0_z, g_z_0_0_0_xy_yz_0_x, g_z_0_0_0_xy_yz_0_y, g_z_0_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_0_x[i] = 2.0 * g_xyz_yz_0_x[i] * a_exp;

        g_z_0_0_0_xy_yz_0_y[i] = 2.0 * g_xyz_yz_0_y[i] * a_exp;

        g_z_0_0_0_xy_yz_0_z[i] = 2.0 * g_xyz_yz_0_z[i] * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_xyz_zz_0_x, g_xyz_zz_0_y, g_xyz_zz_0_z, g_z_0_0_0_xy_zz_0_x, g_z_0_0_0_xy_zz_0_y, g_z_0_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_0_x[i] = 2.0 * g_xyz_zz_0_x[i] * a_exp;

        g_z_0_0_0_xy_zz_0_y[i] = 2.0 * g_xyz_zz_0_y[i] * a_exp;

        g_z_0_0_0_xy_zz_0_z[i] = 2.0 * g_xyz_zz_0_z[i] * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_xzz_xx_0_x, g_xzz_xx_0_y, g_xzz_xx_0_z, g_z_0_0_0_xz_xx_0_x, g_z_0_0_0_xz_xx_0_y, g_z_0_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_0_x[i] = -g_x_xx_0_x[i] + 2.0 * g_xzz_xx_0_x[i] * a_exp;

        g_z_0_0_0_xz_xx_0_y[i] = -g_x_xx_0_y[i] + 2.0 * g_xzz_xx_0_y[i] * a_exp;

        g_z_0_0_0_xz_xx_0_z[i] = -g_x_xx_0_z[i] + 2.0 * g_xzz_xx_0_z[i] * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_xzz_xy_0_x, g_xzz_xy_0_y, g_xzz_xy_0_z, g_z_0_0_0_xz_xy_0_x, g_z_0_0_0_xz_xy_0_y, g_z_0_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_0_x[i] = -g_x_xy_0_x[i] + 2.0 * g_xzz_xy_0_x[i] * a_exp;

        g_z_0_0_0_xz_xy_0_y[i] = -g_x_xy_0_y[i] + 2.0 * g_xzz_xy_0_y[i] * a_exp;

        g_z_0_0_0_xz_xy_0_z[i] = -g_x_xy_0_z[i] + 2.0 * g_xzz_xy_0_z[i] * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_xzz_xz_0_x, g_xzz_xz_0_y, g_xzz_xz_0_z, g_z_0_0_0_xz_xz_0_x, g_z_0_0_0_xz_xz_0_y, g_z_0_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_0_x[i] = -g_x_xz_0_x[i] + 2.0 * g_xzz_xz_0_x[i] * a_exp;

        g_z_0_0_0_xz_xz_0_y[i] = -g_x_xz_0_y[i] + 2.0 * g_xzz_xz_0_y[i] * a_exp;

        g_z_0_0_0_xz_xz_0_z[i] = -g_x_xz_0_z[i] + 2.0 * g_xzz_xz_0_z[i] * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_xzz_yy_0_x, g_xzz_yy_0_y, g_xzz_yy_0_z, g_z_0_0_0_xz_yy_0_x, g_z_0_0_0_xz_yy_0_y, g_z_0_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_0_x[i] = -g_x_yy_0_x[i] + 2.0 * g_xzz_yy_0_x[i] * a_exp;

        g_z_0_0_0_xz_yy_0_y[i] = -g_x_yy_0_y[i] + 2.0 * g_xzz_yy_0_y[i] * a_exp;

        g_z_0_0_0_xz_yy_0_z[i] = -g_x_yy_0_z[i] + 2.0 * g_xzz_yy_0_z[i] * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_xzz_yz_0_x, g_xzz_yz_0_y, g_xzz_yz_0_z, g_z_0_0_0_xz_yz_0_x, g_z_0_0_0_xz_yz_0_y, g_z_0_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_0_x[i] = -g_x_yz_0_x[i] + 2.0 * g_xzz_yz_0_x[i] * a_exp;

        g_z_0_0_0_xz_yz_0_y[i] = -g_x_yz_0_y[i] + 2.0 * g_xzz_yz_0_y[i] * a_exp;

        g_z_0_0_0_xz_yz_0_z[i] = -g_x_yz_0_z[i] + 2.0 * g_xzz_yz_0_z[i] * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_xzz_zz_0_x, g_xzz_zz_0_y, g_xzz_zz_0_z, g_z_0_0_0_xz_zz_0_x, g_z_0_0_0_xz_zz_0_y, g_z_0_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_0_x[i] = -g_x_zz_0_x[i] + 2.0 * g_xzz_zz_0_x[i] * a_exp;

        g_z_0_0_0_xz_zz_0_y[i] = -g_x_zz_0_y[i] + 2.0 * g_xzz_zz_0_y[i] * a_exp;

        g_z_0_0_0_xz_zz_0_z[i] = -g_x_zz_0_z[i] + 2.0 * g_xzz_zz_0_z[i] * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_yyz_xx_0_x, g_yyz_xx_0_y, g_yyz_xx_0_z, g_z_0_0_0_yy_xx_0_x, g_z_0_0_0_yy_xx_0_y, g_z_0_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_0_x[i] = 2.0 * g_yyz_xx_0_x[i] * a_exp;

        g_z_0_0_0_yy_xx_0_y[i] = 2.0 * g_yyz_xx_0_y[i] * a_exp;

        g_z_0_0_0_yy_xx_0_z[i] = 2.0 * g_yyz_xx_0_z[i] * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_yyz_xy_0_x, g_yyz_xy_0_y, g_yyz_xy_0_z, g_z_0_0_0_yy_xy_0_x, g_z_0_0_0_yy_xy_0_y, g_z_0_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_0_x[i] = 2.0 * g_yyz_xy_0_x[i] * a_exp;

        g_z_0_0_0_yy_xy_0_y[i] = 2.0 * g_yyz_xy_0_y[i] * a_exp;

        g_z_0_0_0_yy_xy_0_z[i] = 2.0 * g_yyz_xy_0_z[i] * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_yyz_xz_0_x, g_yyz_xz_0_y, g_yyz_xz_0_z, g_z_0_0_0_yy_xz_0_x, g_z_0_0_0_yy_xz_0_y, g_z_0_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_0_x[i] = 2.0 * g_yyz_xz_0_x[i] * a_exp;

        g_z_0_0_0_yy_xz_0_y[i] = 2.0 * g_yyz_xz_0_y[i] * a_exp;

        g_z_0_0_0_yy_xz_0_z[i] = 2.0 * g_yyz_xz_0_z[i] * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_yyz_yy_0_x, g_yyz_yy_0_y, g_yyz_yy_0_z, g_z_0_0_0_yy_yy_0_x, g_z_0_0_0_yy_yy_0_y, g_z_0_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_0_x[i] = 2.0 * g_yyz_yy_0_x[i] * a_exp;

        g_z_0_0_0_yy_yy_0_y[i] = 2.0 * g_yyz_yy_0_y[i] * a_exp;

        g_z_0_0_0_yy_yy_0_z[i] = 2.0 * g_yyz_yy_0_z[i] * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_yyz_yz_0_x, g_yyz_yz_0_y, g_yyz_yz_0_z, g_z_0_0_0_yy_yz_0_x, g_z_0_0_0_yy_yz_0_y, g_z_0_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_0_x[i] = 2.0 * g_yyz_yz_0_x[i] * a_exp;

        g_z_0_0_0_yy_yz_0_y[i] = 2.0 * g_yyz_yz_0_y[i] * a_exp;

        g_z_0_0_0_yy_yz_0_z[i] = 2.0 * g_yyz_yz_0_z[i] * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_yyz_zz_0_x, g_yyz_zz_0_y, g_yyz_zz_0_z, g_z_0_0_0_yy_zz_0_x, g_z_0_0_0_yy_zz_0_y, g_z_0_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_0_x[i] = 2.0 * g_yyz_zz_0_x[i] * a_exp;

        g_z_0_0_0_yy_zz_0_y[i] = 2.0 * g_yyz_zz_0_y[i] * a_exp;

        g_z_0_0_0_yy_zz_0_z[i] = 2.0 * g_yyz_zz_0_z[i] * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_yzz_xx_0_x, g_yzz_xx_0_y, g_yzz_xx_0_z, g_z_0_0_0_yz_xx_0_x, g_z_0_0_0_yz_xx_0_y, g_z_0_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_0_x[i] = -g_y_xx_0_x[i] + 2.0 * g_yzz_xx_0_x[i] * a_exp;

        g_z_0_0_0_yz_xx_0_y[i] = -g_y_xx_0_y[i] + 2.0 * g_yzz_xx_0_y[i] * a_exp;

        g_z_0_0_0_yz_xx_0_z[i] = -g_y_xx_0_z[i] + 2.0 * g_yzz_xx_0_z[i] * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_yzz_xy_0_x, g_yzz_xy_0_y, g_yzz_xy_0_z, g_z_0_0_0_yz_xy_0_x, g_z_0_0_0_yz_xy_0_y, g_z_0_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_0_x[i] = -g_y_xy_0_x[i] + 2.0 * g_yzz_xy_0_x[i] * a_exp;

        g_z_0_0_0_yz_xy_0_y[i] = -g_y_xy_0_y[i] + 2.0 * g_yzz_xy_0_y[i] * a_exp;

        g_z_0_0_0_yz_xy_0_z[i] = -g_y_xy_0_z[i] + 2.0 * g_yzz_xy_0_z[i] * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_yzz_xz_0_x, g_yzz_xz_0_y, g_yzz_xz_0_z, g_z_0_0_0_yz_xz_0_x, g_z_0_0_0_yz_xz_0_y, g_z_0_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_0_x[i] = -g_y_xz_0_x[i] + 2.0 * g_yzz_xz_0_x[i] * a_exp;

        g_z_0_0_0_yz_xz_0_y[i] = -g_y_xz_0_y[i] + 2.0 * g_yzz_xz_0_y[i] * a_exp;

        g_z_0_0_0_yz_xz_0_z[i] = -g_y_xz_0_z[i] + 2.0 * g_yzz_xz_0_z[i] * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_yzz_yy_0_x, g_yzz_yy_0_y, g_yzz_yy_0_z, g_z_0_0_0_yz_yy_0_x, g_z_0_0_0_yz_yy_0_y, g_z_0_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_0_x[i] = -g_y_yy_0_x[i] + 2.0 * g_yzz_yy_0_x[i] * a_exp;

        g_z_0_0_0_yz_yy_0_y[i] = -g_y_yy_0_y[i] + 2.0 * g_yzz_yy_0_y[i] * a_exp;

        g_z_0_0_0_yz_yy_0_z[i] = -g_y_yy_0_z[i] + 2.0 * g_yzz_yy_0_z[i] * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_yzz_yz_0_x, g_yzz_yz_0_y, g_yzz_yz_0_z, g_z_0_0_0_yz_yz_0_x, g_z_0_0_0_yz_yz_0_y, g_z_0_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_0_x[i] = -g_y_yz_0_x[i] + 2.0 * g_yzz_yz_0_x[i] * a_exp;

        g_z_0_0_0_yz_yz_0_y[i] = -g_y_yz_0_y[i] + 2.0 * g_yzz_yz_0_y[i] * a_exp;

        g_z_0_0_0_yz_yz_0_z[i] = -g_y_yz_0_z[i] + 2.0 * g_yzz_yz_0_z[i] * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_yzz_zz_0_x, g_yzz_zz_0_y, g_yzz_zz_0_z, g_z_0_0_0_yz_zz_0_x, g_z_0_0_0_yz_zz_0_y, g_z_0_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_0_x[i] = -g_y_zz_0_x[i] + 2.0 * g_yzz_zz_0_x[i] * a_exp;

        g_z_0_0_0_yz_zz_0_y[i] = -g_y_zz_0_y[i] + 2.0 * g_yzz_zz_0_y[i] * a_exp;

        g_z_0_0_0_yz_zz_0_z[i] = -g_y_zz_0_z[i] + 2.0 * g_yzz_zz_0_z[i] * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_0_x, g_z_0_0_0_zz_xx_0_y, g_z_0_0_0_zz_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z, g_zzz_xx_0_x, g_zzz_xx_0_y, g_zzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_0_x[i] = -2.0 * g_z_xx_0_x[i] + 2.0 * g_zzz_xx_0_x[i] * a_exp;

        g_z_0_0_0_zz_xx_0_y[i] = -2.0 * g_z_xx_0_y[i] + 2.0 * g_zzz_xx_0_y[i] * a_exp;

        g_z_0_0_0_zz_xx_0_z[i] = -2.0 * g_z_xx_0_z[i] + 2.0 * g_zzz_xx_0_z[i] * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_0_x, g_z_0_0_0_zz_xy_0_y, g_z_0_0_0_zz_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_zzz_xy_0_x, g_zzz_xy_0_y, g_zzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_0_x[i] = -2.0 * g_z_xy_0_x[i] + 2.0 * g_zzz_xy_0_x[i] * a_exp;

        g_z_0_0_0_zz_xy_0_y[i] = -2.0 * g_z_xy_0_y[i] + 2.0 * g_zzz_xy_0_y[i] * a_exp;

        g_z_0_0_0_zz_xy_0_z[i] = -2.0 * g_z_xy_0_z[i] + 2.0 * g_zzz_xy_0_z[i] * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_0_x, g_z_0_0_0_zz_xz_0_y, g_z_0_0_0_zz_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_zzz_xz_0_x, g_zzz_xz_0_y, g_zzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_0_x[i] = -2.0 * g_z_xz_0_x[i] + 2.0 * g_zzz_xz_0_x[i] * a_exp;

        g_z_0_0_0_zz_xz_0_y[i] = -2.0 * g_z_xz_0_y[i] + 2.0 * g_zzz_xz_0_y[i] * a_exp;

        g_z_0_0_0_zz_xz_0_z[i] = -2.0 * g_z_xz_0_z[i] + 2.0 * g_zzz_xz_0_z[i] * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_0_x, g_z_0_0_0_zz_yy_0_y, g_z_0_0_0_zz_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z, g_zzz_yy_0_x, g_zzz_yy_0_y, g_zzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_0_x[i] = -2.0 * g_z_yy_0_x[i] + 2.0 * g_zzz_yy_0_x[i] * a_exp;

        g_z_0_0_0_zz_yy_0_y[i] = -2.0 * g_z_yy_0_y[i] + 2.0 * g_zzz_yy_0_y[i] * a_exp;

        g_z_0_0_0_zz_yy_0_z[i] = -2.0 * g_z_yy_0_z[i] + 2.0 * g_zzz_yy_0_z[i] * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_0_x, g_z_0_0_0_zz_yz_0_y, g_z_0_0_0_zz_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_zzz_yz_0_x, g_zzz_yz_0_y, g_zzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_0_x[i] = -2.0 * g_z_yz_0_x[i] + 2.0 * g_zzz_yz_0_x[i] * a_exp;

        g_z_0_0_0_zz_yz_0_y[i] = -2.0 * g_z_yz_0_y[i] + 2.0 * g_zzz_yz_0_y[i] * a_exp;

        g_z_0_0_0_zz_yz_0_z[i] = -2.0 * g_z_yz_0_z[i] + 2.0 * g_zzz_yz_0_z[i] * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_0_x, g_z_0_0_0_zz_zz_0_y, g_z_0_0_0_zz_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z, g_zzz_zz_0_x, g_zzz_zz_0_y, g_zzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_0_x[i] = -2.0 * g_z_zz_0_x[i] + 2.0 * g_zzz_zz_0_x[i] * a_exp;

        g_z_0_0_0_zz_zz_0_y[i] = -2.0 * g_z_zz_0_y[i] + 2.0 * g_zzz_zz_0_y[i] * a_exp;

        g_z_0_0_0_zz_zz_0_z[i] = -2.0 * g_z_zz_0_z[i] + 2.0 * g_zzz_zz_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace


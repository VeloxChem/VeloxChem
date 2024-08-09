#include "GeomDeriv1100OfScalarForPDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pdsp_0(CSimdArray<double>& buffer_1100_pdsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_sfsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const CSimdArray<double>& buffer_dfsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pdsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spsp

    auto g_0_x_0_x = buffer_spsp[0];

    auto g_0_x_0_y = buffer_spsp[1];

    auto g_0_x_0_z = buffer_spsp[2];

    auto g_0_y_0_x = buffer_spsp[3];

    auto g_0_y_0_y = buffer_spsp[4];

    auto g_0_y_0_z = buffer_spsp[5];

    auto g_0_z_0_x = buffer_spsp[6];

    auto g_0_z_0_y = buffer_spsp[7];

    auto g_0_z_0_z = buffer_spsp[8];

    /// Set up components of auxilary buffer : buffer_sfsp

    auto g_0_xxx_0_x = buffer_sfsp[0];

    auto g_0_xxx_0_y = buffer_sfsp[1];

    auto g_0_xxx_0_z = buffer_sfsp[2];

    auto g_0_xxy_0_x = buffer_sfsp[3];

    auto g_0_xxy_0_y = buffer_sfsp[4];

    auto g_0_xxy_0_z = buffer_sfsp[5];

    auto g_0_xxz_0_x = buffer_sfsp[6];

    auto g_0_xxz_0_y = buffer_sfsp[7];

    auto g_0_xxz_0_z = buffer_sfsp[8];

    auto g_0_xyy_0_x = buffer_sfsp[9];

    auto g_0_xyy_0_y = buffer_sfsp[10];

    auto g_0_xyy_0_z = buffer_sfsp[11];

    auto g_0_xyz_0_x = buffer_sfsp[12];

    auto g_0_xyz_0_y = buffer_sfsp[13];

    auto g_0_xyz_0_z = buffer_sfsp[14];

    auto g_0_xzz_0_x = buffer_sfsp[15];

    auto g_0_xzz_0_y = buffer_sfsp[16];

    auto g_0_xzz_0_z = buffer_sfsp[17];

    auto g_0_yyy_0_x = buffer_sfsp[18];

    auto g_0_yyy_0_y = buffer_sfsp[19];

    auto g_0_yyy_0_z = buffer_sfsp[20];

    auto g_0_yyz_0_x = buffer_sfsp[21];

    auto g_0_yyz_0_y = buffer_sfsp[22];

    auto g_0_yyz_0_z = buffer_sfsp[23];

    auto g_0_yzz_0_x = buffer_sfsp[24];

    auto g_0_yzz_0_y = buffer_sfsp[25];

    auto g_0_yzz_0_z = buffer_sfsp[26];

    auto g_0_zzz_0_x = buffer_sfsp[27];

    auto g_0_zzz_0_y = buffer_sfsp[28];

    auto g_0_zzz_0_z = buffer_sfsp[29];

    /// Set up components of auxilary buffer : buffer_dpsp

    auto g_xx_x_0_x = buffer_dpsp[0];

    auto g_xx_x_0_y = buffer_dpsp[1];

    auto g_xx_x_0_z = buffer_dpsp[2];

    auto g_xx_y_0_x = buffer_dpsp[3];

    auto g_xx_y_0_y = buffer_dpsp[4];

    auto g_xx_y_0_z = buffer_dpsp[5];

    auto g_xx_z_0_x = buffer_dpsp[6];

    auto g_xx_z_0_y = buffer_dpsp[7];

    auto g_xx_z_0_z = buffer_dpsp[8];

    auto g_xy_x_0_x = buffer_dpsp[9];

    auto g_xy_x_0_y = buffer_dpsp[10];

    auto g_xy_x_0_z = buffer_dpsp[11];

    auto g_xy_y_0_x = buffer_dpsp[12];

    auto g_xy_y_0_y = buffer_dpsp[13];

    auto g_xy_y_0_z = buffer_dpsp[14];

    auto g_xy_z_0_x = buffer_dpsp[15];

    auto g_xy_z_0_y = buffer_dpsp[16];

    auto g_xy_z_0_z = buffer_dpsp[17];

    auto g_xz_x_0_x = buffer_dpsp[18];

    auto g_xz_x_0_y = buffer_dpsp[19];

    auto g_xz_x_0_z = buffer_dpsp[20];

    auto g_xz_y_0_x = buffer_dpsp[21];

    auto g_xz_y_0_y = buffer_dpsp[22];

    auto g_xz_y_0_z = buffer_dpsp[23];

    auto g_xz_z_0_x = buffer_dpsp[24];

    auto g_xz_z_0_y = buffer_dpsp[25];

    auto g_xz_z_0_z = buffer_dpsp[26];

    auto g_yy_x_0_x = buffer_dpsp[27];

    auto g_yy_x_0_y = buffer_dpsp[28];

    auto g_yy_x_0_z = buffer_dpsp[29];

    auto g_yy_y_0_x = buffer_dpsp[30];

    auto g_yy_y_0_y = buffer_dpsp[31];

    auto g_yy_y_0_z = buffer_dpsp[32];

    auto g_yy_z_0_x = buffer_dpsp[33];

    auto g_yy_z_0_y = buffer_dpsp[34];

    auto g_yy_z_0_z = buffer_dpsp[35];

    auto g_yz_x_0_x = buffer_dpsp[36];

    auto g_yz_x_0_y = buffer_dpsp[37];

    auto g_yz_x_0_z = buffer_dpsp[38];

    auto g_yz_y_0_x = buffer_dpsp[39];

    auto g_yz_y_0_y = buffer_dpsp[40];

    auto g_yz_y_0_z = buffer_dpsp[41];

    auto g_yz_z_0_x = buffer_dpsp[42];

    auto g_yz_z_0_y = buffer_dpsp[43];

    auto g_yz_z_0_z = buffer_dpsp[44];

    auto g_zz_x_0_x = buffer_dpsp[45];

    auto g_zz_x_0_y = buffer_dpsp[46];

    auto g_zz_x_0_z = buffer_dpsp[47];

    auto g_zz_y_0_x = buffer_dpsp[48];

    auto g_zz_y_0_y = buffer_dpsp[49];

    auto g_zz_y_0_z = buffer_dpsp[50];

    auto g_zz_z_0_x = buffer_dpsp[51];

    auto g_zz_z_0_y = buffer_dpsp[52];

    auto g_zz_z_0_z = buffer_dpsp[53];

    /// Set up components of auxilary buffer : buffer_dfsp

    auto g_xx_xxx_0_x = buffer_dfsp[0];

    auto g_xx_xxx_0_y = buffer_dfsp[1];

    auto g_xx_xxx_0_z = buffer_dfsp[2];

    auto g_xx_xxy_0_x = buffer_dfsp[3];

    auto g_xx_xxy_0_y = buffer_dfsp[4];

    auto g_xx_xxy_0_z = buffer_dfsp[5];

    auto g_xx_xxz_0_x = buffer_dfsp[6];

    auto g_xx_xxz_0_y = buffer_dfsp[7];

    auto g_xx_xxz_0_z = buffer_dfsp[8];

    auto g_xx_xyy_0_x = buffer_dfsp[9];

    auto g_xx_xyy_0_y = buffer_dfsp[10];

    auto g_xx_xyy_0_z = buffer_dfsp[11];

    auto g_xx_xyz_0_x = buffer_dfsp[12];

    auto g_xx_xyz_0_y = buffer_dfsp[13];

    auto g_xx_xyz_0_z = buffer_dfsp[14];

    auto g_xx_xzz_0_x = buffer_dfsp[15];

    auto g_xx_xzz_0_y = buffer_dfsp[16];

    auto g_xx_xzz_0_z = buffer_dfsp[17];

    auto g_xx_yyy_0_x = buffer_dfsp[18];

    auto g_xx_yyy_0_y = buffer_dfsp[19];

    auto g_xx_yyy_0_z = buffer_dfsp[20];

    auto g_xx_yyz_0_x = buffer_dfsp[21];

    auto g_xx_yyz_0_y = buffer_dfsp[22];

    auto g_xx_yyz_0_z = buffer_dfsp[23];

    auto g_xx_yzz_0_x = buffer_dfsp[24];

    auto g_xx_yzz_0_y = buffer_dfsp[25];

    auto g_xx_yzz_0_z = buffer_dfsp[26];

    auto g_xx_zzz_0_x = buffer_dfsp[27];

    auto g_xx_zzz_0_y = buffer_dfsp[28];

    auto g_xx_zzz_0_z = buffer_dfsp[29];

    auto g_xy_xxx_0_x = buffer_dfsp[30];

    auto g_xy_xxx_0_y = buffer_dfsp[31];

    auto g_xy_xxx_0_z = buffer_dfsp[32];

    auto g_xy_xxy_0_x = buffer_dfsp[33];

    auto g_xy_xxy_0_y = buffer_dfsp[34];

    auto g_xy_xxy_0_z = buffer_dfsp[35];

    auto g_xy_xxz_0_x = buffer_dfsp[36];

    auto g_xy_xxz_0_y = buffer_dfsp[37];

    auto g_xy_xxz_0_z = buffer_dfsp[38];

    auto g_xy_xyy_0_x = buffer_dfsp[39];

    auto g_xy_xyy_0_y = buffer_dfsp[40];

    auto g_xy_xyy_0_z = buffer_dfsp[41];

    auto g_xy_xyz_0_x = buffer_dfsp[42];

    auto g_xy_xyz_0_y = buffer_dfsp[43];

    auto g_xy_xyz_0_z = buffer_dfsp[44];

    auto g_xy_xzz_0_x = buffer_dfsp[45];

    auto g_xy_xzz_0_y = buffer_dfsp[46];

    auto g_xy_xzz_0_z = buffer_dfsp[47];

    auto g_xy_yyy_0_x = buffer_dfsp[48];

    auto g_xy_yyy_0_y = buffer_dfsp[49];

    auto g_xy_yyy_0_z = buffer_dfsp[50];

    auto g_xy_yyz_0_x = buffer_dfsp[51];

    auto g_xy_yyz_0_y = buffer_dfsp[52];

    auto g_xy_yyz_0_z = buffer_dfsp[53];

    auto g_xy_yzz_0_x = buffer_dfsp[54];

    auto g_xy_yzz_0_y = buffer_dfsp[55];

    auto g_xy_yzz_0_z = buffer_dfsp[56];

    auto g_xy_zzz_0_x = buffer_dfsp[57];

    auto g_xy_zzz_0_y = buffer_dfsp[58];

    auto g_xy_zzz_0_z = buffer_dfsp[59];

    auto g_xz_xxx_0_x = buffer_dfsp[60];

    auto g_xz_xxx_0_y = buffer_dfsp[61];

    auto g_xz_xxx_0_z = buffer_dfsp[62];

    auto g_xz_xxy_0_x = buffer_dfsp[63];

    auto g_xz_xxy_0_y = buffer_dfsp[64];

    auto g_xz_xxy_0_z = buffer_dfsp[65];

    auto g_xz_xxz_0_x = buffer_dfsp[66];

    auto g_xz_xxz_0_y = buffer_dfsp[67];

    auto g_xz_xxz_0_z = buffer_dfsp[68];

    auto g_xz_xyy_0_x = buffer_dfsp[69];

    auto g_xz_xyy_0_y = buffer_dfsp[70];

    auto g_xz_xyy_0_z = buffer_dfsp[71];

    auto g_xz_xyz_0_x = buffer_dfsp[72];

    auto g_xz_xyz_0_y = buffer_dfsp[73];

    auto g_xz_xyz_0_z = buffer_dfsp[74];

    auto g_xz_xzz_0_x = buffer_dfsp[75];

    auto g_xz_xzz_0_y = buffer_dfsp[76];

    auto g_xz_xzz_0_z = buffer_dfsp[77];

    auto g_xz_yyy_0_x = buffer_dfsp[78];

    auto g_xz_yyy_0_y = buffer_dfsp[79];

    auto g_xz_yyy_0_z = buffer_dfsp[80];

    auto g_xz_yyz_0_x = buffer_dfsp[81];

    auto g_xz_yyz_0_y = buffer_dfsp[82];

    auto g_xz_yyz_0_z = buffer_dfsp[83];

    auto g_xz_yzz_0_x = buffer_dfsp[84];

    auto g_xz_yzz_0_y = buffer_dfsp[85];

    auto g_xz_yzz_0_z = buffer_dfsp[86];

    auto g_xz_zzz_0_x = buffer_dfsp[87];

    auto g_xz_zzz_0_y = buffer_dfsp[88];

    auto g_xz_zzz_0_z = buffer_dfsp[89];

    auto g_yy_xxx_0_x = buffer_dfsp[90];

    auto g_yy_xxx_0_y = buffer_dfsp[91];

    auto g_yy_xxx_0_z = buffer_dfsp[92];

    auto g_yy_xxy_0_x = buffer_dfsp[93];

    auto g_yy_xxy_0_y = buffer_dfsp[94];

    auto g_yy_xxy_0_z = buffer_dfsp[95];

    auto g_yy_xxz_0_x = buffer_dfsp[96];

    auto g_yy_xxz_0_y = buffer_dfsp[97];

    auto g_yy_xxz_0_z = buffer_dfsp[98];

    auto g_yy_xyy_0_x = buffer_dfsp[99];

    auto g_yy_xyy_0_y = buffer_dfsp[100];

    auto g_yy_xyy_0_z = buffer_dfsp[101];

    auto g_yy_xyz_0_x = buffer_dfsp[102];

    auto g_yy_xyz_0_y = buffer_dfsp[103];

    auto g_yy_xyz_0_z = buffer_dfsp[104];

    auto g_yy_xzz_0_x = buffer_dfsp[105];

    auto g_yy_xzz_0_y = buffer_dfsp[106];

    auto g_yy_xzz_0_z = buffer_dfsp[107];

    auto g_yy_yyy_0_x = buffer_dfsp[108];

    auto g_yy_yyy_0_y = buffer_dfsp[109];

    auto g_yy_yyy_0_z = buffer_dfsp[110];

    auto g_yy_yyz_0_x = buffer_dfsp[111];

    auto g_yy_yyz_0_y = buffer_dfsp[112];

    auto g_yy_yyz_0_z = buffer_dfsp[113];

    auto g_yy_yzz_0_x = buffer_dfsp[114];

    auto g_yy_yzz_0_y = buffer_dfsp[115];

    auto g_yy_yzz_0_z = buffer_dfsp[116];

    auto g_yy_zzz_0_x = buffer_dfsp[117];

    auto g_yy_zzz_0_y = buffer_dfsp[118];

    auto g_yy_zzz_0_z = buffer_dfsp[119];

    auto g_yz_xxx_0_x = buffer_dfsp[120];

    auto g_yz_xxx_0_y = buffer_dfsp[121];

    auto g_yz_xxx_0_z = buffer_dfsp[122];

    auto g_yz_xxy_0_x = buffer_dfsp[123];

    auto g_yz_xxy_0_y = buffer_dfsp[124];

    auto g_yz_xxy_0_z = buffer_dfsp[125];

    auto g_yz_xxz_0_x = buffer_dfsp[126];

    auto g_yz_xxz_0_y = buffer_dfsp[127];

    auto g_yz_xxz_0_z = buffer_dfsp[128];

    auto g_yz_xyy_0_x = buffer_dfsp[129];

    auto g_yz_xyy_0_y = buffer_dfsp[130];

    auto g_yz_xyy_0_z = buffer_dfsp[131];

    auto g_yz_xyz_0_x = buffer_dfsp[132];

    auto g_yz_xyz_0_y = buffer_dfsp[133];

    auto g_yz_xyz_0_z = buffer_dfsp[134];

    auto g_yz_xzz_0_x = buffer_dfsp[135];

    auto g_yz_xzz_0_y = buffer_dfsp[136];

    auto g_yz_xzz_0_z = buffer_dfsp[137];

    auto g_yz_yyy_0_x = buffer_dfsp[138];

    auto g_yz_yyy_0_y = buffer_dfsp[139];

    auto g_yz_yyy_0_z = buffer_dfsp[140];

    auto g_yz_yyz_0_x = buffer_dfsp[141];

    auto g_yz_yyz_0_y = buffer_dfsp[142];

    auto g_yz_yyz_0_z = buffer_dfsp[143];

    auto g_yz_yzz_0_x = buffer_dfsp[144];

    auto g_yz_yzz_0_y = buffer_dfsp[145];

    auto g_yz_yzz_0_z = buffer_dfsp[146];

    auto g_yz_zzz_0_x = buffer_dfsp[147];

    auto g_yz_zzz_0_y = buffer_dfsp[148];

    auto g_yz_zzz_0_z = buffer_dfsp[149];

    auto g_zz_xxx_0_x = buffer_dfsp[150];

    auto g_zz_xxx_0_y = buffer_dfsp[151];

    auto g_zz_xxx_0_z = buffer_dfsp[152];

    auto g_zz_xxy_0_x = buffer_dfsp[153];

    auto g_zz_xxy_0_y = buffer_dfsp[154];

    auto g_zz_xxy_0_z = buffer_dfsp[155];

    auto g_zz_xxz_0_x = buffer_dfsp[156];

    auto g_zz_xxz_0_y = buffer_dfsp[157];

    auto g_zz_xxz_0_z = buffer_dfsp[158];

    auto g_zz_xyy_0_x = buffer_dfsp[159];

    auto g_zz_xyy_0_y = buffer_dfsp[160];

    auto g_zz_xyy_0_z = buffer_dfsp[161];

    auto g_zz_xyz_0_x = buffer_dfsp[162];

    auto g_zz_xyz_0_y = buffer_dfsp[163];

    auto g_zz_xyz_0_z = buffer_dfsp[164];

    auto g_zz_xzz_0_x = buffer_dfsp[165];

    auto g_zz_xzz_0_y = buffer_dfsp[166];

    auto g_zz_xzz_0_z = buffer_dfsp[167];

    auto g_zz_yyy_0_x = buffer_dfsp[168];

    auto g_zz_yyy_0_y = buffer_dfsp[169];

    auto g_zz_yyy_0_z = buffer_dfsp[170];

    auto g_zz_yyz_0_x = buffer_dfsp[171];

    auto g_zz_yyz_0_y = buffer_dfsp[172];

    auto g_zz_yyz_0_z = buffer_dfsp[173];

    auto g_zz_yzz_0_x = buffer_dfsp[174];

    auto g_zz_yzz_0_y = buffer_dfsp[175];

    auto g_zz_yzz_0_z = buffer_dfsp[176];

    auto g_zz_zzz_0_x = buffer_dfsp[177];

    auto g_zz_zzz_0_y = buffer_dfsp[178];

    auto g_zz_zzz_0_z = buffer_dfsp[179];

    /// Set up components of integrals buffer : buffer_1100_pdsp

    auto g_x_x_0_0_x_xx_0_x = buffer_1100_pdsp[0];

    auto g_x_x_0_0_x_xx_0_y = buffer_1100_pdsp[1];

    auto g_x_x_0_0_x_xx_0_z = buffer_1100_pdsp[2];

    auto g_x_x_0_0_x_xy_0_x = buffer_1100_pdsp[3];

    auto g_x_x_0_0_x_xy_0_y = buffer_1100_pdsp[4];

    auto g_x_x_0_0_x_xy_0_z = buffer_1100_pdsp[5];

    auto g_x_x_0_0_x_xz_0_x = buffer_1100_pdsp[6];

    auto g_x_x_0_0_x_xz_0_y = buffer_1100_pdsp[7];

    auto g_x_x_0_0_x_xz_0_z = buffer_1100_pdsp[8];

    auto g_x_x_0_0_x_yy_0_x = buffer_1100_pdsp[9];

    auto g_x_x_0_0_x_yy_0_y = buffer_1100_pdsp[10];

    auto g_x_x_0_0_x_yy_0_z = buffer_1100_pdsp[11];

    auto g_x_x_0_0_x_yz_0_x = buffer_1100_pdsp[12];

    auto g_x_x_0_0_x_yz_0_y = buffer_1100_pdsp[13];

    auto g_x_x_0_0_x_yz_0_z = buffer_1100_pdsp[14];

    auto g_x_x_0_0_x_zz_0_x = buffer_1100_pdsp[15];

    auto g_x_x_0_0_x_zz_0_y = buffer_1100_pdsp[16];

    auto g_x_x_0_0_x_zz_0_z = buffer_1100_pdsp[17];

    auto g_x_x_0_0_y_xx_0_x = buffer_1100_pdsp[18];

    auto g_x_x_0_0_y_xx_0_y = buffer_1100_pdsp[19];

    auto g_x_x_0_0_y_xx_0_z = buffer_1100_pdsp[20];

    auto g_x_x_0_0_y_xy_0_x = buffer_1100_pdsp[21];

    auto g_x_x_0_0_y_xy_0_y = buffer_1100_pdsp[22];

    auto g_x_x_0_0_y_xy_0_z = buffer_1100_pdsp[23];

    auto g_x_x_0_0_y_xz_0_x = buffer_1100_pdsp[24];

    auto g_x_x_0_0_y_xz_0_y = buffer_1100_pdsp[25];

    auto g_x_x_0_0_y_xz_0_z = buffer_1100_pdsp[26];

    auto g_x_x_0_0_y_yy_0_x = buffer_1100_pdsp[27];

    auto g_x_x_0_0_y_yy_0_y = buffer_1100_pdsp[28];

    auto g_x_x_0_0_y_yy_0_z = buffer_1100_pdsp[29];

    auto g_x_x_0_0_y_yz_0_x = buffer_1100_pdsp[30];

    auto g_x_x_0_0_y_yz_0_y = buffer_1100_pdsp[31];

    auto g_x_x_0_0_y_yz_0_z = buffer_1100_pdsp[32];

    auto g_x_x_0_0_y_zz_0_x = buffer_1100_pdsp[33];

    auto g_x_x_0_0_y_zz_0_y = buffer_1100_pdsp[34];

    auto g_x_x_0_0_y_zz_0_z = buffer_1100_pdsp[35];

    auto g_x_x_0_0_z_xx_0_x = buffer_1100_pdsp[36];

    auto g_x_x_0_0_z_xx_0_y = buffer_1100_pdsp[37];

    auto g_x_x_0_0_z_xx_0_z = buffer_1100_pdsp[38];

    auto g_x_x_0_0_z_xy_0_x = buffer_1100_pdsp[39];

    auto g_x_x_0_0_z_xy_0_y = buffer_1100_pdsp[40];

    auto g_x_x_0_0_z_xy_0_z = buffer_1100_pdsp[41];

    auto g_x_x_0_0_z_xz_0_x = buffer_1100_pdsp[42];

    auto g_x_x_0_0_z_xz_0_y = buffer_1100_pdsp[43];

    auto g_x_x_0_0_z_xz_0_z = buffer_1100_pdsp[44];

    auto g_x_x_0_0_z_yy_0_x = buffer_1100_pdsp[45];

    auto g_x_x_0_0_z_yy_0_y = buffer_1100_pdsp[46];

    auto g_x_x_0_0_z_yy_0_z = buffer_1100_pdsp[47];

    auto g_x_x_0_0_z_yz_0_x = buffer_1100_pdsp[48];

    auto g_x_x_0_0_z_yz_0_y = buffer_1100_pdsp[49];

    auto g_x_x_0_0_z_yz_0_z = buffer_1100_pdsp[50];

    auto g_x_x_0_0_z_zz_0_x = buffer_1100_pdsp[51];

    auto g_x_x_0_0_z_zz_0_y = buffer_1100_pdsp[52];

    auto g_x_x_0_0_z_zz_0_z = buffer_1100_pdsp[53];

    auto g_x_y_0_0_x_xx_0_x = buffer_1100_pdsp[54];

    auto g_x_y_0_0_x_xx_0_y = buffer_1100_pdsp[55];

    auto g_x_y_0_0_x_xx_0_z = buffer_1100_pdsp[56];

    auto g_x_y_0_0_x_xy_0_x = buffer_1100_pdsp[57];

    auto g_x_y_0_0_x_xy_0_y = buffer_1100_pdsp[58];

    auto g_x_y_0_0_x_xy_0_z = buffer_1100_pdsp[59];

    auto g_x_y_0_0_x_xz_0_x = buffer_1100_pdsp[60];

    auto g_x_y_0_0_x_xz_0_y = buffer_1100_pdsp[61];

    auto g_x_y_0_0_x_xz_0_z = buffer_1100_pdsp[62];

    auto g_x_y_0_0_x_yy_0_x = buffer_1100_pdsp[63];

    auto g_x_y_0_0_x_yy_0_y = buffer_1100_pdsp[64];

    auto g_x_y_0_0_x_yy_0_z = buffer_1100_pdsp[65];

    auto g_x_y_0_0_x_yz_0_x = buffer_1100_pdsp[66];

    auto g_x_y_0_0_x_yz_0_y = buffer_1100_pdsp[67];

    auto g_x_y_0_0_x_yz_0_z = buffer_1100_pdsp[68];

    auto g_x_y_0_0_x_zz_0_x = buffer_1100_pdsp[69];

    auto g_x_y_0_0_x_zz_0_y = buffer_1100_pdsp[70];

    auto g_x_y_0_0_x_zz_0_z = buffer_1100_pdsp[71];

    auto g_x_y_0_0_y_xx_0_x = buffer_1100_pdsp[72];

    auto g_x_y_0_0_y_xx_0_y = buffer_1100_pdsp[73];

    auto g_x_y_0_0_y_xx_0_z = buffer_1100_pdsp[74];

    auto g_x_y_0_0_y_xy_0_x = buffer_1100_pdsp[75];

    auto g_x_y_0_0_y_xy_0_y = buffer_1100_pdsp[76];

    auto g_x_y_0_0_y_xy_0_z = buffer_1100_pdsp[77];

    auto g_x_y_0_0_y_xz_0_x = buffer_1100_pdsp[78];

    auto g_x_y_0_0_y_xz_0_y = buffer_1100_pdsp[79];

    auto g_x_y_0_0_y_xz_0_z = buffer_1100_pdsp[80];

    auto g_x_y_0_0_y_yy_0_x = buffer_1100_pdsp[81];

    auto g_x_y_0_0_y_yy_0_y = buffer_1100_pdsp[82];

    auto g_x_y_0_0_y_yy_0_z = buffer_1100_pdsp[83];

    auto g_x_y_0_0_y_yz_0_x = buffer_1100_pdsp[84];

    auto g_x_y_0_0_y_yz_0_y = buffer_1100_pdsp[85];

    auto g_x_y_0_0_y_yz_0_z = buffer_1100_pdsp[86];

    auto g_x_y_0_0_y_zz_0_x = buffer_1100_pdsp[87];

    auto g_x_y_0_0_y_zz_0_y = buffer_1100_pdsp[88];

    auto g_x_y_0_0_y_zz_0_z = buffer_1100_pdsp[89];

    auto g_x_y_0_0_z_xx_0_x = buffer_1100_pdsp[90];

    auto g_x_y_0_0_z_xx_0_y = buffer_1100_pdsp[91];

    auto g_x_y_0_0_z_xx_0_z = buffer_1100_pdsp[92];

    auto g_x_y_0_0_z_xy_0_x = buffer_1100_pdsp[93];

    auto g_x_y_0_0_z_xy_0_y = buffer_1100_pdsp[94];

    auto g_x_y_0_0_z_xy_0_z = buffer_1100_pdsp[95];

    auto g_x_y_0_0_z_xz_0_x = buffer_1100_pdsp[96];

    auto g_x_y_0_0_z_xz_0_y = buffer_1100_pdsp[97];

    auto g_x_y_0_0_z_xz_0_z = buffer_1100_pdsp[98];

    auto g_x_y_0_0_z_yy_0_x = buffer_1100_pdsp[99];

    auto g_x_y_0_0_z_yy_0_y = buffer_1100_pdsp[100];

    auto g_x_y_0_0_z_yy_0_z = buffer_1100_pdsp[101];

    auto g_x_y_0_0_z_yz_0_x = buffer_1100_pdsp[102];

    auto g_x_y_0_0_z_yz_0_y = buffer_1100_pdsp[103];

    auto g_x_y_0_0_z_yz_0_z = buffer_1100_pdsp[104];

    auto g_x_y_0_0_z_zz_0_x = buffer_1100_pdsp[105];

    auto g_x_y_0_0_z_zz_0_y = buffer_1100_pdsp[106];

    auto g_x_y_0_0_z_zz_0_z = buffer_1100_pdsp[107];

    auto g_x_z_0_0_x_xx_0_x = buffer_1100_pdsp[108];

    auto g_x_z_0_0_x_xx_0_y = buffer_1100_pdsp[109];

    auto g_x_z_0_0_x_xx_0_z = buffer_1100_pdsp[110];

    auto g_x_z_0_0_x_xy_0_x = buffer_1100_pdsp[111];

    auto g_x_z_0_0_x_xy_0_y = buffer_1100_pdsp[112];

    auto g_x_z_0_0_x_xy_0_z = buffer_1100_pdsp[113];

    auto g_x_z_0_0_x_xz_0_x = buffer_1100_pdsp[114];

    auto g_x_z_0_0_x_xz_0_y = buffer_1100_pdsp[115];

    auto g_x_z_0_0_x_xz_0_z = buffer_1100_pdsp[116];

    auto g_x_z_0_0_x_yy_0_x = buffer_1100_pdsp[117];

    auto g_x_z_0_0_x_yy_0_y = buffer_1100_pdsp[118];

    auto g_x_z_0_0_x_yy_0_z = buffer_1100_pdsp[119];

    auto g_x_z_0_0_x_yz_0_x = buffer_1100_pdsp[120];

    auto g_x_z_0_0_x_yz_0_y = buffer_1100_pdsp[121];

    auto g_x_z_0_0_x_yz_0_z = buffer_1100_pdsp[122];

    auto g_x_z_0_0_x_zz_0_x = buffer_1100_pdsp[123];

    auto g_x_z_0_0_x_zz_0_y = buffer_1100_pdsp[124];

    auto g_x_z_0_0_x_zz_0_z = buffer_1100_pdsp[125];

    auto g_x_z_0_0_y_xx_0_x = buffer_1100_pdsp[126];

    auto g_x_z_0_0_y_xx_0_y = buffer_1100_pdsp[127];

    auto g_x_z_0_0_y_xx_0_z = buffer_1100_pdsp[128];

    auto g_x_z_0_0_y_xy_0_x = buffer_1100_pdsp[129];

    auto g_x_z_0_0_y_xy_0_y = buffer_1100_pdsp[130];

    auto g_x_z_0_0_y_xy_0_z = buffer_1100_pdsp[131];

    auto g_x_z_0_0_y_xz_0_x = buffer_1100_pdsp[132];

    auto g_x_z_0_0_y_xz_0_y = buffer_1100_pdsp[133];

    auto g_x_z_0_0_y_xz_0_z = buffer_1100_pdsp[134];

    auto g_x_z_0_0_y_yy_0_x = buffer_1100_pdsp[135];

    auto g_x_z_0_0_y_yy_0_y = buffer_1100_pdsp[136];

    auto g_x_z_0_0_y_yy_0_z = buffer_1100_pdsp[137];

    auto g_x_z_0_0_y_yz_0_x = buffer_1100_pdsp[138];

    auto g_x_z_0_0_y_yz_0_y = buffer_1100_pdsp[139];

    auto g_x_z_0_0_y_yz_0_z = buffer_1100_pdsp[140];

    auto g_x_z_0_0_y_zz_0_x = buffer_1100_pdsp[141];

    auto g_x_z_0_0_y_zz_0_y = buffer_1100_pdsp[142];

    auto g_x_z_0_0_y_zz_0_z = buffer_1100_pdsp[143];

    auto g_x_z_0_0_z_xx_0_x = buffer_1100_pdsp[144];

    auto g_x_z_0_0_z_xx_0_y = buffer_1100_pdsp[145];

    auto g_x_z_0_0_z_xx_0_z = buffer_1100_pdsp[146];

    auto g_x_z_0_0_z_xy_0_x = buffer_1100_pdsp[147];

    auto g_x_z_0_0_z_xy_0_y = buffer_1100_pdsp[148];

    auto g_x_z_0_0_z_xy_0_z = buffer_1100_pdsp[149];

    auto g_x_z_0_0_z_xz_0_x = buffer_1100_pdsp[150];

    auto g_x_z_0_0_z_xz_0_y = buffer_1100_pdsp[151];

    auto g_x_z_0_0_z_xz_0_z = buffer_1100_pdsp[152];

    auto g_x_z_0_0_z_yy_0_x = buffer_1100_pdsp[153];

    auto g_x_z_0_0_z_yy_0_y = buffer_1100_pdsp[154];

    auto g_x_z_0_0_z_yy_0_z = buffer_1100_pdsp[155];

    auto g_x_z_0_0_z_yz_0_x = buffer_1100_pdsp[156];

    auto g_x_z_0_0_z_yz_0_y = buffer_1100_pdsp[157];

    auto g_x_z_0_0_z_yz_0_z = buffer_1100_pdsp[158];

    auto g_x_z_0_0_z_zz_0_x = buffer_1100_pdsp[159];

    auto g_x_z_0_0_z_zz_0_y = buffer_1100_pdsp[160];

    auto g_x_z_0_0_z_zz_0_z = buffer_1100_pdsp[161];

    auto g_y_x_0_0_x_xx_0_x = buffer_1100_pdsp[162];

    auto g_y_x_0_0_x_xx_0_y = buffer_1100_pdsp[163];

    auto g_y_x_0_0_x_xx_0_z = buffer_1100_pdsp[164];

    auto g_y_x_0_0_x_xy_0_x = buffer_1100_pdsp[165];

    auto g_y_x_0_0_x_xy_0_y = buffer_1100_pdsp[166];

    auto g_y_x_0_0_x_xy_0_z = buffer_1100_pdsp[167];

    auto g_y_x_0_0_x_xz_0_x = buffer_1100_pdsp[168];

    auto g_y_x_0_0_x_xz_0_y = buffer_1100_pdsp[169];

    auto g_y_x_0_0_x_xz_0_z = buffer_1100_pdsp[170];

    auto g_y_x_0_0_x_yy_0_x = buffer_1100_pdsp[171];

    auto g_y_x_0_0_x_yy_0_y = buffer_1100_pdsp[172];

    auto g_y_x_0_0_x_yy_0_z = buffer_1100_pdsp[173];

    auto g_y_x_0_0_x_yz_0_x = buffer_1100_pdsp[174];

    auto g_y_x_0_0_x_yz_0_y = buffer_1100_pdsp[175];

    auto g_y_x_0_0_x_yz_0_z = buffer_1100_pdsp[176];

    auto g_y_x_0_0_x_zz_0_x = buffer_1100_pdsp[177];

    auto g_y_x_0_0_x_zz_0_y = buffer_1100_pdsp[178];

    auto g_y_x_0_0_x_zz_0_z = buffer_1100_pdsp[179];

    auto g_y_x_0_0_y_xx_0_x = buffer_1100_pdsp[180];

    auto g_y_x_0_0_y_xx_0_y = buffer_1100_pdsp[181];

    auto g_y_x_0_0_y_xx_0_z = buffer_1100_pdsp[182];

    auto g_y_x_0_0_y_xy_0_x = buffer_1100_pdsp[183];

    auto g_y_x_0_0_y_xy_0_y = buffer_1100_pdsp[184];

    auto g_y_x_0_0_y_xy_0_z = buffer_1100_pdsp[185];

    auto g_y_x_0_0_y_xz_0_x = buffer_1100_pdsp[186];

    auto g_y_x_0_0_y_xz_0_y = buffer_1100_pdsp[187];

    auto g_y_x_0_0_y_xz_0_z = buffer_1100_pdsp[188];

    auto g_y_x_0_0_y_yy_0_x = buffer_1100_pdsp[189];

    auto g_y_x_0_0_y_yy_0_y = buffer_1100_pdsp[190];

    auto g_y_x_0_0_y_yy_0_z = buffer_1100_pdsp[191];

    auto g_y_x_0_0_y_yz_0_x = buffer_1100_pdsp[192];

    auto g_y_x_0_0_y_yz_0_y = buffer_1100_pdsp[193];

    auto g_y_x_0_0_y_yz_0_z = buffer_1100_pdsp[194];

    auto g_y_x_0_0_y_zz_0_x = buffer_1100_pdsp[195];

    auto g_y_x_0_0_y_zz_0_y = buffer_1100_pdsp[196];

    auto g_y_x_0_0_y_zz_0_z = buffer_1100_pdsp[197];

    auto g_y_x_0_0_z_xx_0_x = buffer_1100_pdsp[198];

    auto g_y_x_0_0_z_xx_0_y = buffer_1100_pdsp[199];

    auto g_y_x_0_0_z_xx_0_z = buffer_1100_pdsp[200];

    auto g_y_x_0_0_z_xy_0_x = buffer_1100_pdsp[201];

    auto g_y_x_0_0_z_xy_0_y = buffer_1100_pdsp[202];

    auto g_y_x_0_0_z_xy_0_z = buffer_1100_pdsp[203];

    auto g_y_x_0_0_z_xz_0_x = buffer_1100_pdsp[204];

    auto g_y_x_0_0_z_xz_0_y = buffer_1100_pdsp[205];

    auto g_y_x_0_0_z_xz_0_z = buffer_1100_pdsp[206];

    auto g_y_x_0_0_z_yy_0_x = buffer_1100_pdsp[207];

    auto g_y_x_0_0_z_yy_0_y = buffer_1100_pdsp[208];

    auto g_y_x_0_0_z_yy_0_z = buffer_1100_pdsp[209];

    auto g_y_x_0_0_z_yz_0_x = buffer_1100_pdsp[210];

    auto g_y_x_0_0_z_yz_0_y = buffer_1100_pdsp[211];

    auto g_y_x_0_0_z_yz_0_z = buffer_1100_pdsp[212];

    auto g_y_x_0_0_z_zz_0_x = buffer_1100_pdsp[213];

    auto g_y_x_0_0_z_zz_0_y = buffer_1100_pdsp[214];

    auto g_y_x_0_0_z_zz_0_z = buffer_1100_pdsp[215];

    auto g_y_y_0_0_x_xx_0_x = buffer_1100_pdsp[216];

    auto g_y_y_0_0_x_xx_0_y = buffer_1100_pdsp[217];

    auto g_y_y_0_0_x_xx_0_z = buffer_1100_pdsp[218];

    auto g_y_y_0_0_x_xy_0_x = buffer_1100_pdsp[219];

    auto g_y_y_0_0_x_xy_0_y = buffer_1100_pdsp[220];

    auto g_y_y_0_0_x_xy_0_z = buffer_1100_pdsp[221];

    auto g_y_y_0_0_x_xz_0_x = buffer_1100_pdsp[222];

    auto g_y_y_0_0_x_xz_0_y = buffer_1100_pdsp[223];

    auto g_y_y_0_0_x_xz_0_z = buffer_1100_pdsp[224];

    auto g_y_y_0_0_x_yy_0_x = buffer_1100_pdsp[225];

    auto g_y_y_0_0_x_yy_0_y = buffer_1100_pdsp[226];

    auto g_y_y_0_0_x_yy_0_z = buffer_1100_pdsp[227];

    auto g_y_y_0_0_x_yz_0_x = buffer_1100_pdsp[228];

    auto g_y_y_0_0_x_yz_0_y = buffer_1100_pdsp[229];

    auto g_y_y_0_0_x_yz_0_z = buffer_1100_pdsp[230];

    auto g_y_y_0_0_x_zz_0_x = buffer_1100_pdsp[231];

    auto g_y_y_0_0_x_zz_0_y = buffer_1100_pdsp[232];

    auto g_y_y_0_0_x_zz_0_z = buffer_1100_pdsp[233];

    auto g_y_y_0_0_y_xx_0_x = buffer_1100_pdsp[234];

    auto g_y_y_0_0_y_xx_0_y = buffer_1100_pdsp[235];

    auto g_y_y_0_0_y_xx_0_z = buffer_1100_pdsp[236];

    auto g_y_y_0_0_y_xy_0_x = buffer_1100_pdsp[237];

    auto g_y_y_0_0_y_xy_0_y = buffer_1100_pdsp[238];

    auto g_y_y_0_0_y_xy_0_z = buffer_1100_pdsp[239];

    auto g_y_y_0_0_y_xz_0_x = buffer_1100_pdsp[240];

    auto g_y_y_0_0_y_xz_0_y = buffer_1100_pdsp[241];

    auto g_y_y_0_0_y_xz_0_z = buffer_1100_pdsp[242];

    auto g_y_y_0_0_y_yy_0_x = buffer_1100_pdsp[243];

    auto g_y_y_0_0_y_yy_0_y = buffer_1100_pdsp[244];

    auto g_y_y_0_0_y_yy_0_z = buffer_1100_pdsp[245];

    auto g_y_y_0_0_y_yz_0_x = buffer_1100_pdsp[246];

    auto g_y_y_0_0_y_yz_0_y = buffer_1100_pdsp[247];

    auto g_y_y_0_0_y_yz_0_z = buffer_1100_pdsp[248];

    auto g_y_y_0_0_y_zz_0_x = buffer_1100_pdsp[249];

    auto g_y_y_0_0_y_zz_0_y = buffer_1100_pdsp[250];

    auto g_y_y_0_0_y_zz_0_z = buffer_1100_pdsp[251];

    auto g_y_y_0_0_z_xx_0_x = buffer_1100_pdsp[252];

    auto g_y_y_0_0_z_xx_0_y = buffer_1100_pdsp[253];

    auto g_y_y_0_0_z_xx_0_z = buffer_1100_pdsp[254];

    auto g_y_y_0_0_z_xy_0_x = buffer_1100_pdsp[255];

    auto g_y_y_0_0_z_xy_0_y = buffer_1100_pdsp[256];

    auto g_y_y_0_0_z_xy_0_z = buffer_1100_pdsp[257];

    auto g_y_y_0_0_z_xz_0_x = buffer_1100_pdsp[258];

    auto g_y_y_0_0_z_xz_0_y = buffer_1100_pdsp[259];

    auto g_y_y_0_0_z_xz_0_z = buffer_1100_pdsp[260];

    auto g_y_y_0_0_z_yy_0_x = buffer_1100_pdsp[261];

    auto g_y_y_0_0_z_yy_0_y = buffer_1100_pdsp[262];

    auto g_y_y_0_0_z_yy_0_z = buffer_1100_pdsp[263];

    auto g_y_y_0_0_z_yz_0_x = buffer_1100_pdsp[264];

    auto g_y_y_0_0_z_yz_0_y = buffer_1100_pdsp[265];

    auto g_y_y_0_0_z_yz_0_z = buffer_1100_pdsp[266];

    auto g_y_y_0_0_z_zz_0_x = buffer_1100_pdsp[267];

    auto g_y_y_0_0_z_zz_0_y = buffer_1100_pdsp[268];

    auto g_y_y_0_0_z_zz_0_z = buffer_1100_pdsp[269];

    auto g_y_z_0_0_x_xx_0_x = buffer_1100_pdsp[270];

    auto g_y_z_0_0_x_xx_0_y = buffer_1100_pdsp[271];

    auto g_y_z_0_0_x_xx_0_z = buffer_1100_pdsp[272];

    auto g_y_z_0_0_x_xy_0_x = buffer_1100_pdsp[273];

    auto g_y_z_0_0_x_xy_0_y = buffer_1100_pdsp[274];

    auto g_y_z_0_0_x_xy_0_z = buffer_1100_pdsp[275];

    auto g_y_z_0_0_x_xz_0_x = buffer_1100_pdsp[276];

    auto g_y_z_0_0_x_xz_0_y = buffer_1100_pdsp[277];

    auto g_y_z_0_0_x_xz_0_z = buffer_1100_pdsp[278];

    auto g_y_z_0_0_x_yy_0_x = buffer_1100_pdsp[279];

    auto g_y_z_0_0_x_yy_0_y = buffer_1100_pdsp[280];

    auto g_y_z_0_0_x_yy_0_z = buffer_1100_pdsp[281];

    auto g_y_z_0_0_x_yz_0_x = buffer_1100_pdsp[282];

    auto g_y_z_0_0_x_yz_0_y = buffer_1100_pdsp[283];

    auto g_y_z_0_0_x_yz_0_z = buffer_1100_pdsp[284];

    auto g_y_z_0_0_x_zz_0_x = buffer_1100_pdsp[285];

    auto g_y_z_0_0_x_zz_0_y = buffer_1100_pdsp[286];

    auto g_y_z_0_0_x_zz_0_z = buffer_1100_pdsp[287];

    auto g_y_z_0_0_y_xx_0_x = buffer_1100_pdsp[288];

    auto g_y_z_0_0_y_xx_0_y = buffer_1100_pdsp[289];

    auto g_y_z_0_0_y_xx_0_z = buffer_1100_pdsp[290];

    auto g_y_z_0_0_y_xy_0_x = buffer_1100_pdsp[291];

    auto g_y_z_0_0_y_xy_0_y = buffer_1100_pdsp[292];

    auto g_y_z_0_0_y_xy_0_z = buffer_1100_pdsp[293];

    auto g_y_z_0_0_y_xz_0_x = buffer_1100_pdsp[294];

    auto g_y_z_0_0_y_xz_0_y = buffer_1100_pdsp[295];

    auto g_y_z_0_0_y_xz_0_z = buffer_1100_pdsp[296];

    auto g_y_z_0_0_y_yy_0_x = buffer_1100_pdsp[297];

    auto g_y_z_0_0_y_yy_0_y = buffer_1100_pdsp[298];

    auto g_y_z_0_0_y_yy_0_z = buffer_1100_pdsp[299];

    auto g_y_z_0_0_y_yz_0_x = buffer_1100_pdsp[300];

    auto g_y_z_0_0_y_yz_0_y = buffer_1100_pdsp[301];

    auto g_y_z_0_0_y_yz_0_z = buffer_1100_pdsp[302];

    auto g_y_z_0_0_y_zz_0_x = buffer_1100_pdsp[303];

    auto g_y_z_0_0_y_zz_0_y = buffer_1100_pdsp[304];

    auto g_y_z_0_0_y_zz_0_z = buffer_1100_pdsp[305];

    auto g_y_z_0_0_z_xx_0_x = buffer_1100_pdsp[306];

    auto g_y_z_0_0_z_xx_0_y = buffer_1100_pdsp[307];

    auto g_y_z_0_0_z_xx_0_z = buffer_1100_pdsp[308];

    auto g_y_z_0_0_z_xy_0_x = buffer_1100_pdsp[309];

    auto g_y_z_0_0_z_xy_0_y = buffer_1100_pdsp[310];

    auto g_y_z_0_0_z_xy_0_z = buffer_1100_pdsp[311];

    auto g_y_z_0_0_z_xz_0_x = buffer_1100_pdsp[312];

    auto g_y_z_0_0_z_xz_0_y = buffer_1100_pdsp[313];

    auto g_y_z_0_0_z_xz_0_z = buffer_1100_pdsp[314];

    auto g_y_z_0_0_z_yy_0_x = buffer_1100_pdsp[315];

    auto g_y_z_0_0_z_yy_0_y = buffer_1100_pdsp[316];

    auto g_y_z_0_0_z_yy_0_z = buffer_1100_pdsp[317];

    auto g_y_z_0_0_z_yz_0_x = buffer_1100_pdsp[318];

    auto g_y_z_0_0_z_yz_0_y = buffer_1100_pdsp[319];

    auto g_y_z_0_0_z_yz_0_z = buffer_1100_pdsp[320];

    auto g_y_z_0_0_z_zz_0_x = buffer_1100_pdsp[321];

    auto g_y_z_0_0_z_zz_0_y = buffer_1100_pdsp[322];

    auto g_y_z_0_0_z_zz_0_z = buffer_1100_pdsp[323];

    auto g_z_x_0_0_x_xx_0_x = buffer_1100_pdsp[324];

    auto g_z_x_0_0_x_xx_0_y = buffer_1100_pdsp[325];

    auto g_z_x_0_0_x_xx_0_z = buffer_1100_pdsp[326];

    auto g_z_x_0_0_x_xy_0_x = buffer_1100_pdsp[327];

    auto g_z_x_0_0_x_xy_0_y = buffer_1100_pdsp[328];

    auto g_z_x_0_0_x_xy_0_z = buffer_1100_pdsp[329];

    auto g_z_x_0_0_x_xz_0_x = buffer_1100_pdsp[330];

    auto g_z_x_0_0_x_xz_0_y = buffer_1100_pdsp[331];

    auto g_z_x_0_0_x_xz_0_z = buffer_1100_pdsp[332];

    auto g_z_x_0_0_x_yy_0_x = buffer_1100_pdsp[333];

    auto g_z_x_0_0_x_yy_0_y = buffer_1100_pdsp[334];

    auto g_z_x_0_0_x_yy_0_z = buffer_1100_pdsp[335];

    auto g_z_x_0_0_x_yz_0_x = buffer_1100_pdsp[336];

    auto g_z_x_0_0_x_yz_0_y = buffer_1100_pdsp[337];

    auto g_z_x_0_0_x_yz_0_z = buffer_1100_pdsp[338];

    auto g_z_x_0_0_x_zz_0_x = buffer_1100_pdsp[339];

    auto g_z_x_0_0_x_zz_0_y = buffer_1100_pdsp[340];

    auto g_z_x_0_0_x_zz_0_z = buffer_1100_pdsp[341];

    auto g_z_x_0_0_y_xx_0_x = buffer_1100_pdsp[342];

    auto g_z_x_0_0_y_xx_0_y = buffer_1100_pdsp[343];

    auto g_z_x_0_0_y_xx_0_z = buffer_1100_pdsp[344];

    auto g_z_x_0_0_y_xy_0_x = buffer_1100_pdsp[345];

    auto g_z_x_0_0_y_xy_0_y = buffer_1100_pdsp[346];

    auto g_z_x_0_0_y_xy_0_z = buffer_1100_pdsp[347];

    auto g_z_x_0_0_y_xz_0_x = buffer_1100_pdsp[348];

    auto g_z_x_0_0_y_xz_0_y = buffer_1100_pdsp[349];

    auto g_z_x_0_0_y_xz_0_z = buffer_1100_pdsp[350];

    auto g_z_x_0_0_y_yy_0_x = buffer_1100_pdsp[351];

    auto g_z_x_0_0_y_yy_0_y = buffer_1100_pdsp[352];

    auto g_z_x_0_0_y_yy_0_z = buffer_1100_pdsp[353];

    auto g_z_x_0_0_y_yz_0_x = buffer_1100_pdsp[354];

    auto g_z_x_0_0_y_yz_0_y = buffer_1100_pdsp[355];

    auto g_z_x_0_0_y_yz_0_z = buffer_1100_pdsp[356];

    auto g_z_x_0_0_y_zz_0_x = buffer_1100_pdsp[357];

    auto g_z_x_0_0_y_zz_0_y = buffer_1100_pdsp[358];

    auto g_z_x_0_0_y_zz_0_z = buffer_1100_pdsp[359];

    auto g_z_x_0_0_z_xx_0_x = buffer_1100_pdsp[360];

    auto g_z_x_0_0_z_xx_0_y = buffer_1100_pdsp[361];

    auto g_z_x_0_0_z_xx_0_z = buffer_1100_pdsp[362];

    auto g_z_x_0_0_z_xy_0_x = buffer_1100_pdsp[363];

    auto g_z_x_0_0_z_xy_0_y = buffer_1100_pdsp[364];

    auto g_z_x_0_0_z_xy_0_z = buffer_1100_pdsp[365];

    auto g_z_x_0_0_z_xz_0_x = buffer_1100_pdsp[366];

    auto g_z_x_0_0_z_xz_0_y = buffer_1100_pdsp[367];

    auto g_z_x_0_0_z_xz_0_z = buffer_1100_pdsp[368];

    auto g_z_x_0_0_z_yy_0_x = buffer_1100_pdsp[369];

    auto g_z_x_0_0_z_yy_0_y = buffer_1100_pdsp[370];

    auto g_z_x_0_0_z_yy_0_z = buffer_1100_pdsp[371];

    auto g_z_x_0_0_z_yz_0_x = buffer_1100_pdsp[372];

    auto g_z_x_0_0_z_yz_0_y = buffer_1100_pdsp[373];

    auto g_z_x_0_0_z_yz_0_z = buffer_1100_pdsp[374];

    auto g_z_x_0_0_z_zz_0_x = buffer_1100_pdsp[375];

    auto g_z_x_0_0_z_zz_0_y = buffer_1100_pdsp[376];

    auto g_z_x_0_0_z_zz_0_z = buffer_1100_pdsp[377];

    auto g_z_y_0_0_x_xx_0_x = buffer_1100_pdsp[378];

    auto g_z_y_0_0_x_xx_0_y = buffer_1100_pdsp[379];

    auto g_z_y_0_0_x_xx_0_z = buffer_1100_pdsp[380];

    auto g_z_y_0_0_x_xy_0_x = buffer_1100_pdsp[381];

    auto g_z_y_0_0_x_xy_0_y = buffer_1100_pdsp[382];

    auto g_z_y_0_0_x_xy_0_z = buffer_1100_pdsp[383];

    auto g_z_y_0_0_x_xz_0_x = buffer_1100_pdsp[384];

    auto g_z_y_0_0_x_xz_0_y = buffer_1100_pdsp[385];

    auto g_z_y_0_0_x_xz_0_z = buffer_1100_pdsp[386];

    auto g_z_y_0_0_x_yy_0_x = buffer_1100_pdsp[387];

    auto g_z_y_0_0_x_yy_0_y = buffer_1100_pdsp[388];

    auto g_z_y_0_0_x_yy_0_z = buffer_1100_pdsp[389];

    auto g_z_y_0_0_x_yz_0_x = buffer_1100_pdsp[390];

    auto g_z_y_0_0_x_yz_0_y = buffer_1100_pdsp[391];

    auto g_z_y_0_0_x_yz_0_z = buffer_1100_pdsp[392];

    auto g_z_y_0_0_x_zz_0_x = buffer_1100_pdsp[393];

    auto g_z_y_0_0_x_zz_0_y = buffer_1100_pdsp[394];

    auto g_z_y_0_0_x_zz_0_z = buffer_1100_pdsp[395];

    auto g_z_y_0_0_y_xx_0_x = buffer_1100_pdsp[396];

    auto g_z_y_0_0_y_xx_0_y = buffer_1100_pdsp[397];

    auto g_z_y_0_0_y_xx_0_z = buffer_1100_pdsp[398];

    auto g_z_y_0_0_y_xy_0_x = buffer_1100_pdsp[399];

    auto g_z_y_0_0_y_xy_0_y = buffer_1100_pdsp[400];

    auto g_z_y_0_0_y_xy_0_z = buffer_1100_pdsp[401];

    auto g_z_y_0_0_y_xz_0_x = buffer_1100_pdsp[402];

    auto g_z_y_0_0_y_xz_0_y = buffer_1100_pdsp[403];

    auto g_z_y_0_0_y_xz_0_z = buffer_1100_pdsp[404];

    auto g_z_y_0_0_y_yy_0_x = buffer_1100_pdsp[405];

    auto g_z_y_0_0_y_yy_0_y = buffer_1100_pdsp[406];

    auto g_z_y_0_0_y_yy_0_z = buffer_1100_pdsp[407];

    auto g_z_y_0_0_y_yz_0_x = buffer_1100_pdsp[408];

    auto g_z_y_0_0_y_yz_0_y = buffer_1100_pdsp[409];

    auto g_z_y_0_0_y_yz_0_z = buffer_1100_pdsp[410];

    auto g_z_y_0_0_y_zz_0_x = buffer_1100_pdsp[411];

    auto g_z_y_0_0_y_zz_0_y = buffer_1100_pdsp[412];

    auto g_z_y_0_0_y_zz_0_z = buffer_1100_pdsp[413];

    auto g_z_y_0_0_z_xx_0_x = buffer_1100_pdsp[414];

    auto g_z_y_0_0_z_xx_0_y = buffer_1100_pdsp[415];

    auto g_z_y_0_0_z_xx_0_z = buffer_1100_pdsp[416];

    auto g_z_y_0_0_z_xy_0_x = buffer_1100_pdsp[417];

    auto g_z_y_0_0_z_xy_0_y = buffer_1100_pdsp[418];

    auto g_z_y_0_0_z_xy_0_z = buffer_1100_pdsp[419];

    auto g_z_y_0_0_z_xz_0_x = buffer_1100_pdsp[420];

    auto g_z_y_0_0_z_xz_0_y = buffer_1100_pdsp[421];

    auto g_z_y_0_0_z_xz_0_z = buffer_1100_pdsp[422];

    auto g_z_y_0_0_z_yy_0_x = buffer_1100_pdsp[423];

    auto g_z_y_0_0_z_yy_0_y = buffer_1100_pdsp[424];

    auto g_z_y_0_0_z_yy_0_z = buffer_1100_pdsp[425];

    auto g_z_y_0_0_z_yz_0_x = buffer_1100_pdsp[426];

    auto g_z_y_0_0_z_yz_0_y = buffer_1100_pdsp[427];

    auto g_z_y_0_0_z_yz_0_z = buffer_1100_pdsp[428];

    auto g_z_y_0_0_z_zz_0_x = buffer_1100_pdsp[429];

    auto g_z_y_0_0_z_zz_0_y = buffer_1100_pdsp[430];

    auto g_z_y_0_0_z_zz_0_z = buffer_1100_pdsp[431];

    auto g_z_z_0_0_x_xx_0_x = buffer_1100_pdsp[432];

    auto g_z_z_0_0_x_xx_0_y = buffer_1100_pdsp[433];

    auto g_z_z_0_0_x_xx_0_z = buffer_1100_pdsp[434];

    auto g_z_z_0_0_x_xy_0_x = buffer_1100_pdsp[435];

    auto g_z_z_0_0_x_xy_0_y = buffer_1100_pdsp[436];

    auto g_z_z_0_0_x_xy_0_z = buffer_1100_pdsp[437];

    auto g_z_z_0_0_x_xz_0_x = buffer_1100_pdsp[438];

    auto g_z_z_0_0_x_xz_0_y = buffer_1100_pdsp[439];

    auto g_z_z_0_0_x_xz_0_z = buffer_1100_pdsp[440];

    auto g_z_z_0_0_x_yy_0_x = buffer_1100_pdsp[441];

    auto g_z_z_0_0_x_yy_0_y = buffer_1100_pdsp[442];

    auto g_z_z_0_0_x_yy_0_z = buffer_1100_pdsp[443];

    auto g_z_z_0_0_x_yz_0_x = buffer_1100_pdsp[444];

    auto g_z_z_0_0_x_yz_0_y = buffer_1100_pdsp[445];

    auto g_z_z_0_0_x_yz_0_z = buffer_1100_pdsp[446];

    auto g_z_z_0_0_x_zz_0_x = buffer_1100_pdsp[447];

    auto g_z_z_0_0_x_zz_0_y = buffer_1100_pdsp[448];

    auto g_z_z_0_0_x_zz_0_z = buffer_1100_pdsp[449];

    auto g_z_z_0_0_y_xx_0_x = buffer_1100_pdsp[450];

    auto g_z_z_0_0_y_xx_0_y = buffer_1100_pdsp[451];

    auto g_z_z_0_0_y_xx_0_z = buffer_1100_pdsp[452];

    auto g_z_z_0_0_y_xy_0_x = buffer_1100_pdsp[453];

    auto g_z_z_0_0_y_xy_0_y = buffer_1100_pdsp[454];

    auto g_z_z_0_0_y_xy_0_z = buffer_1100_pdsp[455];

    auto g_z_z_0_0_y_xz_0_x = buffer_1100_pdsp[456];

    auto g_z_z_0_0_y_xz_0_y = buffer_1100_pdsp[457];

    auto g_z_z_0_0_y_xz_0_z = buffer_1100_pdsp[458];

    auto g_z_z_0_0_y_yy_0_x = buffer_1100_pdsp[459];

    auto g_z_z_0_0_y_yy_0_y = buffer_1100_pdsp[460];

    auto g_z_z_0_0_y_yy_0_z = buffer_1100_pdsp[461];

    auto g_z_z_0_0_y_yz_0_x = buffer_1100_pdsp[462];

    auto g_z_z_0_0_y_yz_0_y = buffer_1100_pdsp[463];

    auto g_z_z_0_0_y_yz_0_z = buffer_1100_pdsp[464];

    auto g_z_z_0_0_y_zz_0_x = buffer_1100_pdsp[465];

    auto g_z_z_0_0_y_zz_0_y = buffer_1100_pdsp[466];

    auto g_z_z_0_0_y_zz_0_z = buffer_1100_pdsp[467];

    auto g_z_z_0_0_z_xx_0_x = buffer_1100_pdsp[468];

    auto g_z_z_0_0_z_xx_0_y = buffer_1100_pdsp[469];

    auto g_z_z_0_0_z_xx_0_z = buffer_1100_pdsp[470];

    auto g_z_z_0_0_z_xy_0_x = buffer_1100_pdsp[471];

    auto g_z_z_0_0_z_xy_0_y = buffer_1100_pdsp[472];

    auto g_z_z_0_0_z_xy_0_z = buffer_1100_pdsp[473];

    auto g_z_z_0_0_z_xz_0_x = buffer_1100_pdsp[474];

    auto g_z_z_0_0_z_xz_0_y = buffer_1100_pdsp[475];

    auto g_z_z_0_0_z_xz_0_z = buffer_1100_pdsp[476];

    auto g_z_z_0_0_z_yy_0_x = buffer_1100_pdsp[477];

    auto g_z_z_0_0_z_yy_0_y = buffer_1100_pdsp[478];

    auto g_z_z_0_0_z_yy_0_z = buffer_1100_pdsp[479];

    auto g_z_z_0_0_z_yz_0_x = buffer_1100_pdsp[480];

    auto g_z_z_0_0_z_yz_0_y = buffer_1100_pdsp[481];

    auto g_z_z_0_0_z_yz_0_z = buffer_1100_pdsp[482];

    auto g_z_z_0_0_z_zz_0_x = buffer_1100_pdsp[483];

    auto g_z_z_0_0_z_zz_0_y = buffer_1100_pdsp[484];

    auto g_z_z_0_0_z_zz_0_z = buffer_1100_pdsp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xxx_0_x, g_0_xxx_0_y, g_0_xxx_0_z, g_x_x_0_0_x_xx_0_x, g_x_x_0_0_x_xx_0_y, g_x_x_0_0_x_xx_0_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_xxx_0_x, g_xx_xxx_0_y, g_xx_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_0_x[i] = 2.0 * g_0_x_0_x[i] - 2.0 * g_0_xxx_0_x[i] * b_exp - 4.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_y[i] = 2.0 * g_0_x_0_y[i] - 2.0 * g_0_xxx_0_y[i] * b_exp - 4.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xx_0_z[i] = 2.0 * g_0_x_0_z[i] - 2.0 * g_0_xxx_0_z[i] * b_exp - 4.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_x_x_0_0_x_xy_0_x, g_x_x_0_0_x_xy_0_y, g_x_x_0_0_x_xy_0_z, g_xx_xxy_0_x, g_xx_xxy_0_y, g_xx_xxy_0_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xy_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_xxy_0_x[i] * b_exp - 2.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_xxy_0_y[i] * b_exp - 2.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_xxy_0_z[i] * b_exp - 2.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_x_x_0_0_x_xz_0_x, g_x_x_0_0_x_xz_0_y, g_x_x_0_0_x_xz_0_z, g_xx_xxz_0_x, g_xx_xxz_0_y, g_xx_xxz_0_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_xxz_0_x[i] * b_exp - 2.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_xxz_0_y[i] * b_exp - 2.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_xxz_0_z[i] * b_exp - 2.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_x_x_0_0_x_yy_0_x, g_x_x_0_0_x_yy_0_y, g_x_x_0_0_x_yy_0_z, g_xx_xyy_0_x, g_xx_xyy_0_y, g_xx_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yy_0_x[i] = -2.0 * g_0_xyy_0_x[i] * b_exp + 4.0 * g_xx_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_y[i] = -2.0 * g_0_xyy_0_y[i] * b_exp + 4.0 * g_xx_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_z[i] = -2.0 * g_0_xyy_0_z[i] * b_exp + 4.0 * g_xx_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_x_x_0_0_x_yz_0_x, g_x_x_0_0_x_yz_0_y, g_x_x_0_0_x_yz_0_z, g_xx_xyz_0_x, g_xx_xyz_0_y, g_xx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_yz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_xx_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_xx_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_xx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_x_x_0_0_x_zz_0_x, g_x_x_0_0_x_zz_0_y, g_x_x_0_0_x_zz_0_z, g_xx_xzz_0_x, g_xx_xzz_0_y, g_xx_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_zz_0_x[i] = -2.0 * g_0_xzz_0_x[i] * b_exp + 4.0 * g_xx_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_y[i] = -2.0 * g_0_xzz_0_y[i] * b_exp + 4.0 * g_xx_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_z[i] = -2.0 * g_0_xzz_0_z[i] * b_exp + 4.0 * g_xx_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_0_x, g_x_x_0_0_y_xx_0_y, g_x_x_0_0_y_xx_0_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xxx_0_x, g_xy_xxx_0_y, g_xy_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_0_x[i] = -4.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_y[i] = -4.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xx_0_z[i] = -4.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_y_xy_0_x, g_x_x_0_0_y_xy_0_y, g_x_x_0_0_y_xy_0_z, g_xy_xxy_0_x, g_xy_xxy_0_y, g_xy_xxy_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xy_0_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_y_xz_0_x, g_x_x_0_0_y_xz_0_y, g_x_x_0_0_y_xz_0_z, g_xy_xxz_0_x, g_xy_xxz_0_y, g_xy_xxz_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xz_0_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0_y_yy_0_x, g_x_x_0_0_y_yy_0_y, g_x_x_0_0_y_yy_0_z, g_xy_xyy_0_x, g_xy_xyy_0_y, g_xy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yy_0_x[i] = 4.0 * g_xy_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_y[i] = 4.0 * g_xy_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_z[i] = 4.0 * g_xy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_x_0_0_y_yz_0_x, g_x_x_0_0_y_yz_0_y, g_x_x_0_0_y_yz_0_z, g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_yz_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_x_0_0_y_zz_0_x, g_x_x_0_0_y_zz_0_y, g_x_x_0_0_y_zz_0_z, g_xy_xzz_0_x, g_xy_xzz_0_y, g_xy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_zz_0_x[i] = 4.0 * g_xy_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_y[i] = 4.0 * g_xy_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_z[i] = 4.0 * g_xy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_0_x, g_x_x_0_0_z_xx_0_y, g_x_x_0_0_z_xx_0_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xxx_0_x, g_xz_xxx_0_y, g_xz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_0_x[i] = -4.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_y[i] = -4.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xx_0_z[i] = -4.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_x_0_0_z_xy_0_x, g_x_x_0_0_z_xy_0_y, g_x_x_0_0_z_xy_0_z, g_xz_xxy_0_x, g_xz_xxy_0_y, g_xz_xxy_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xy_0_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_0_z_xz_0_x, g_x_x_0_0_z_xz_0_y, g_x_x_0_0_z_xz_0_z, g_xz_xxz_0_x, g_xz_xxz_0_y, g_xz_xxz_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xz_0_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0_z_yy_0_x, g_x_x_0_0_z_yy_0_y, g_x_x_0_0_z_yy_0_z, g_xz_xyy_0_x, g_xz_xyy_0_y, g_xz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yy_0_x[i] = 4.0 * g_xz_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_y[i] = 4.0 * g_xz_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_z[i] = 4.0 * g_xz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_x_0_0_z_yz_0_x, g_x_x_0_0_z_yz_0_y, g_x_x_0_0_z_yz_0_z, g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_yz_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_x_0_0_z_zz_0_x, g_x_x_0_0_z_zz_0_y, g_x_x_0_0_z_zz_0_z, g_xz_xzz_0_x, g_xz_xzz_0_y, g_xz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_zz_0_x[i] = 4.0 * g_xz_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_y[i] = 4.0 * g_xz_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_z[i] = 4.0 * g_xz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_x_y_0_0_x_xx_0_x, g_x_y_0_0_x_xx_0_y, g_x_y_0_0_x_xx_0_z, g_xx_xxy_0_x, g_xx_xxy_0_y, g_xx_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_0_x[i] = -2.0 * g_0_xxy_0_x[i] * b_exp + 4.0 * g_xx_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_y[i] = -2.0 * g_0_xxy_0_y[i] * b_exp + 4.0 * g_xx_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xx_0_z[i] = -2.0 * g_0_xxy_0_z[i] * b_exp + 4.0 * g_xx_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_x_y_0_0_x_xy_0_x, g_x_y_0_0_x_xy_0_y, g_x_y_0_0_x_xy_0_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_xyy_0_x, g_xx_xyy_0_y, g_xx_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xy_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xyy_0_x[i] * b_exp - 2.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xyy_0_y[i] * b_exp - 2.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xyy_0_z[i] * b_exp - 2.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_x_y_0_0_x_xz_0_x, g_x_y_0_0_x_xz_0_y, g_x_y_0_0_x_xz_0_z, g_xx_xyz_0_x, g_xx_xyz_0_y, g_xx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_xx_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_xx_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_xx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yyy_0_x, g_0_yyy_0_y, g_0_yyy_0_z, g_x_y_0_0_x_yy_0_x, g_x_y_0_0_x_yy_0_y, g_x_y_0_0_x_yy_0_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z, g_xx_yyy_0_x, g_xx_yyy_0_y, g_xx_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yy_0_x[i] = 2.0 * g_0_y_0_x[i] - 2.0 * g_0_yyy_0_x[i] * b_exp - 4.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_y[i] = 2.0 * g_0_y_0_y[i] - 2.0 * g_0_yyy_0_y[i] * b_exp - 4.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_z[i] = 2.0 * g_0_y_0_z[i] - 2.0 * g_0_yyy_0_z[i] * b_exp - 4.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_x_y_0_0_x_yz_0_x, g_x_y_0_0_x_yz_0_y, g_x_y_0_0_x_yz_0_z, g_xx_yyz_0_x, g_xx_yyz_0_y, g_xx_yyz_0_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_yz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_yyz_0_x[i] * b_exp - 2.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_yyz_0_y[i] * b_exp - 2.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_yyz_0_z[i] * b_exp - 2.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_x_y_0_0_x_zz_0_x, g_x_y_0_0_x_zz_0_y, g_x_y_0_0_x_zz_0_z, g_xx_yzz_0_x, g_xx_yzz_0_y, g_xx_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_zz_0_x[i] = -2.0 * g_0_yzz_0_x[i] * b_exp + 4.0 * g_xx_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_y[i] = -2.0 * g_0_yzz_0_y[i] * b_exp + 4.0 * g_xx_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_z[i] = -2.0 * g_0_yzz_0_z[i] * b_exp + 4.0 * g_xx_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_0_x, g_x_y_0_0_y_xx_0_y, g_x_y_0_0_y_xx_0_z, g_xy_xxy_0_x, g_xy_xxy_0_y, g_xy_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_0_x[i] = 4.0 * g_xy_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_y[i] = 4.0 * g_xy_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xx_0_z[i] = 4.0 * g_xy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_y_0_0_y_xy_0_x, g_x_y_0_0_y_xy_0_y, g_x_y_0_0_y_xy_0_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xyy_0_x, g_xy_xyy_0_y, g_xy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xy_0_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_y_0_0_y_xz_0_x, g_x_y_0_0_y_xz_0_y, g_x_y_0_0_y_xz_0_z, g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xz_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_y_0_0_y_yy_0_x, g_x_y_0_0_y_yy_0_y, g_x_y_0_0_y_yy_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_yyy_0_x, g_xy_yyy_0_y, g_xy_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yy_0_x[i] = -4.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_y[i] = -4.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_z[i] = -4.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_y_0_0_y_yz_0_x, g_x_y_0_0_y_yz_0_y, g_x_y_0_0_y_yz_0_z, g_xy_yyz_0_x, g_xy_yyz_0_y, g_xy_yyz_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_yz_0_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_y_0_0_y_zz_0_x, g_x_y_0_0_y_zz_0_y, g_x_y_0_0_y_zz_0_z, g_xy_yzz_0_x, g_xy_yzz_0_y, g_xy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_zz_0_x[i] = 4.0 * g_xy_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_y[i] = 4.0 * g_xy_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_z[i] = 4.0 * g_xy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_0_x, g_x_y_0_0_z_xx_0_y, g_x_y_0_0_z_xx_0_z, g_xz_xxy_0_x, g_xz_xxy_0_y, g_xz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_0_x[i] = 4.0 * g_xz_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_y[i] = 4.0 * g_xz_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xx_0_z[i] = 4.0 * g_xz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_y_0_0_z_xy_0_x, g_x_y_0_0_z_xy_0_y, g_x_y_0_0_z_xy_0_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xyy_0_x, g_xz_xyy_0_y, g_xz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xy_0_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_y_0_0_z_xz_0_x, g_x_y_0_0_z_xz_0_y, g_x_y_0_0_z_xz_0_z, g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xz_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_y_0_0_z_yy_0_x, g_x_y_0_0_z_yy_0_y, g_x_y_0_0_z_yy_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_yyy_0_x, g_xz_yyy_0_y, g_xz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yy_0_x[i] = -4.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_y[i] = -4.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_z[i] = -4.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_y_0_0_z_yz_0_x, g_x_y_0_0_z_yz_0_y, g_x_y_0_0_z_yz_0_z, g_xz_yyz_0_x, g_xz_yyz_0_y, g_xz_yyz_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_yz_0_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_y_0_0_z_zz_0_x, g_x_y_0_0_z_zz_0_y, g_x_y_0_0_z_zz_0_z, g_xz_yzz_0_x, g_xz_yzz_0_y, g_xz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_zz_0_x[i] = 4.0 * g_xz_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_y[i] = 4.0 * g_xz_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_z[i] = 4.0 * g_xz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_x_z_0_0_x_xx_0_x, g_x_z_0_0_x_xx_0_y, g_x_z_0_0_x_xx_0_z, g_xx_xxz_0_x, g_xx_xxz_0_y, g_xx_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_0_x[i] = -2.0 * g_0_xxz_0_x[i] * b_exp + 4.0 * g_xx_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_y[i] = -2.0 * g_0_xxz_0_y[i] * b_exp + 4.0 * g_xx_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xx_0_z[i] = -2.0 * g_0_xxz_0_z[i] * b_exp + 4.0 * g_xx_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_x_z_0_0_x_xy_0_x, g_x_z_0_0_x_xy_0_y, g_x_z_0_0_x_xy_0_z, g_xx_xyz_0_x, g_xx_xyz_0_y, g_xx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xy_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_xx_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_xx_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_xx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_x_z_0_0_x_xz_0_x, g_x_z_0_0_x_xz_0_y, g_x_z_0_0_x_xz_0_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_xzz_0_x, g_xx_xzz_0_y, g_xx_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xz_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xzz_0_x[i] * b_exp - 2.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xzz_0_y[i] * b_exp - 2.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xzz_0_z[i] * b_exp - 2.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_x_z_0_0_x_yy_0_x, g_x_z_0_0_x_yy_0_y, g_x_z_0_0_x_yy_0_z, g_xx_yyz_0_x, g_xx_yyz_0_y, g_xx_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yy_0_x[i] = -2.0 * g_0_yyz_0_x[i] * b_exp + 4.0 * g_xx_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_y[i] = -2.0 * g_0_yyz_0_y[i] * b_exp + 4.0 * g_xx_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_z[i] = -2.0 * g_0_yyz_0_z[i] * b_exp + 4.0 * g_xx_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_x_z_0_0_x_yz_0_x, g_x_z_0_0_x_yz_0_y, g_x_z_0_0_x_yz_0_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z, g_xx_yzz_0_x, g_xx_yzz_0_y, g_xx_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_yz_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_yzz_0_x[i] * b_exp - 2.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_yzz_0_y[i] * b_exp - 2.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_yzz_0_z[i] * b_exp - 2.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_zzz_0_x, g_0_zzz_0_y, g_0_zzz_0_z, g_x_z_0_0_x_zz_0_x, g_x_z_0_0_x_zz_0_y, g_x_z_0_0_x_zz_0_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z, g_xx_zzz_0_x, g_xx_zzz_0_y, g_xx_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_zz_0_x[i] = 2.0 * g_0_z_0_x[i] - 2.0 * g_0_zzz_0_x[i] * b_exp - 4.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_y[i] = 2.0 * g_0_z_0_y[i] - 2.0 * g_0_zzz_0_y[i] * b_exp - 4.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_z[i] = 2.0 * g_0_z_0_z[i] - 2.0 * g_0_zzz_0_z[i] * b_exp - 4.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_0_x, g_x_z_0_0_y_xx_0_y, g_x_z_0_0_y_xx_0_z, g_xy_xxz_0_x, g_xy_xxz_0_y, g_xy_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_0_x[i] = 4.0 * g_xy_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_y[i] = 4.0 * g_xy_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xx_0_z[i] = 4.0 * g_xy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_z_0_0_y_xy_0_x, g_x_z_0_0_y_xy_0_y, g_x_z_0_0_y_xy_0_z, g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xy_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_z_0_0_y_xz_0_x, g_x_z_0_0_y_xz_0_y, g_x_z_0_0_y_xz_0_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xzz_0_x, g_xy_xzz_0_y, g_xy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xz_0_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_z_0_0_y_yy_0_x, g_x_z_0_0_y_yy_0_y, g_x_z_0_0_y_yy_0_z, g_xy_yyz_0_x, g_xy_yyz_0_y, g_xy_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yy_0_x[i] = 4.0 * g_xy_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_y[i] = 4.0 * g_xy_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_z[i] = 4.0 * g_xy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_z_0_0_y_yz_0_x, g_x_z_0_0_y_yz_0_y, g_x_z_0_0_y_yz_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_yzz_0_x, g_xy_yzz_0_y, g_xy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_yz_0_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_z_0_0_y_zz_0_x, g_x_z_0_0_y_zz_0_y, g_x_z_0_0_y_zz_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_zzz_0_x, g_xy_zzz_0_y, g_xy_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_zz_0_x[i] = -4.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_y[i] = -4.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_z[i] = -4.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_0_x, g_x_z_0_0_z_xx_0_y, g_x_z_0_0_z_xx_0_z, g_xz_xxz_0_x, g_xz_xxz_0_y, g_xz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_0_x[i] = 4.0 * g_xz_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_y[i] = 4.0 * g_xz_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xx_0_z[i] = 4.0 * g_xz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_z_0_0_z_xy_0_x, g_x_z_0_0_z_xy_0_y, g_x_z_0_0_z_xy_0_z, g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xy_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_z_0_0_z_xz_0_x, g_x_z_0_0_z_xz_0_y, g_x_z_0_0_z_xz_0_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xzz_0_x, g_xz_xzz_0_y, g_xz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xz_0_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_z_0_0_z_yy_0_x, g_x_z_0_0_z_yy_0_y, g_x_z_0_0_z_yy_0_z, g_xz_yyz_0_x, g_xz_yyz_0_y, g_xz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yy_0_x[i] = 4.0 * g_xz_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_y[i] = 4.0 * g_xz_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_z[i] = 4.0 * g_xz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_z_0_0_z_yz_0_x, g_x_z_0_0_z_yz_0_y, g_x_z_0_0_z_yz_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_yzz_0_x, g_xz_yzz_0_y, g_xz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_yz_0_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_z_0_0_z_zz_0_x, g_x_z_0_0_z_zz_0_y, g_x_z_0_0_z_zz_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_zzz_0_x, g_xz_zzz_0_y, g_xz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_zz_0_x[i] = -4.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_y[i] = -4.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_z[i] = -4.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xxx_0_x, g_xy_xxx_0_y, g_xy_xxx_0_z, g_y_x_0_0_x_xx_0_x, g_y_x_0_0_x_xx_0_y, g_y_x_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_0_x[i] = -4.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_y[i] = -4.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xx_0_z[i] = -4.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xy_xxy_0_x, g_xy_xxy_0_y, g_xy_xxy_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_y_x_0_0_x_xy_0_x, g_y_x_0_0_x_xy_0_y, g_y_x_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xy_0_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xy_xxz_0_x, g_xy_xxz_0_y, g_xy_xxz_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_y_x_0_0_x_xz_0_x, g_y_x_0_0_x_xz_0_y, g_y_x_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xz_0_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xy_xyy_0_x, g_xy_xyy_0_y, g_xy_xyy_0_z, g_y_x_0_0_x_yy_0_x, g_y_x_0_0_x_yy_0_y, g_y_x_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yy_0_x[i] = 4.0 * g_xy_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_y[i] = 4.0 * g_xy_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_z[i] = 4.0 * g_xy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z, g_y_x_0_0_x_yz_0_x, g_y_x_0_0_x_yz_0_y, g_y_x_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_yz_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xy_xzz_0_x, g_xy_xzz_0_y, g_xy_xzz_0_z, g_y_x_0_0_x_zz_0_x, g_y_x_0_0_x_zz_0_y, g_y_x_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_zz_0_x[i] = 4.0 * g_xy_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_y[i] = 4.0 * g_xy_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_z[i] = 4.0 * g_xy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xxx_0_x, g_0_xxx_0_y, g_0_xxx_0_z, g_y_x_0_0_y_xx_0_x, g_y_x_0_0_y_xx_0_y, g_y_x_0_0_y_xx_0_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_xxx_0_x, g_yy_xxx_0_y, g_yy_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_0_x[i] = 2.0 * g_0_x_0_x[i] - 2.0 * g_0_xxx_0_x[i] * b_exp - 4.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_y[i] = 2.0 * g_0_x_0_y[i] - 2.0 * g_0_xxx_0_y[i] * b_exp - 4.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xx_0_z[i] = 2.0 * g_0_x_0_z[i] - 2.0 * g_0_xxx_0_z[i] * b_exp - 4.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_y_x_0_0_y_xy_0_x, g_y_x_0_0_y_xy_0_y, g_y_x_0_0_y_xy_0_z, g_yy_xxy_0_x, g_yy_xxy_0_y, g_yy_xxy_0_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xy_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_xxy_0_x[i] * b_exp - 2.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_xxy_0_y[i] * b_exp - 2.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_xxy_0_z[i] * b_exp - 2.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_y_x_0_0_y_xz_0_x, g_y_x_0_0_y_xz_0_y, g_y_x_0_0_y_xz_0_z, g_yy_xxz_0_x, g_yy_xxz_0_y, g_yy_xxz_0_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_xxz_0_x[i] * b_exp - 2.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_xxz_0_y[i] * b_exp - 2.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_xxz_0_z[i] * b_exp - 2.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_y_x_0_0_y_yy_0_x, g_y_x_0_0_y_yy_0_y, g_y_x_0_0_y_yy_0_z, g_yy_xyy_0_x, g_yy_xyy_0_y, g_yy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yy_0_x[i] = -2.0 * g_0_xyy_0_x[i] * b_exp + 4.0 * g_yy_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_y[i] = -2.0 * g_0_xyy_0_y[i] * b_exp + 4.0 * g_yy_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_z[i] = -2.0 * g_0_xyy_0_z[i] * b_exp + 4.0 * g_yy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_y_x_0_0_y_yz_0_x, g_y_x_0_0_y_yz_0_y, g_y_x_0_0_y_yz_0_z, g_yy_xyz_0_x, g_yy_xyz_0_y, g_yy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_yz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_yy_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_yy_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_yy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_y_x_0_0_y_zz_0_x, g_y_x_0_0_y_zz_0_y, g_y_x_0_0_y_zz_0_z, g_yy_xzz_0_x, g_yy_xzz_0_y, g_yy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_zz_0_x[i] = -2.0 * g_0_xzz_0_x[i] * b_exp + 4.0 * g_yy_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_y[i] = -2.0 * g_0_xzz_0_y[i] * b_exp + 4.0 * g_yy_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_z[i] = -2.0 * g_0_xzz_0_z[i] * b_exp + 4.0 * g_yy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_0_x, g_y_x_0_0_z_xx_0_y, g_y_x_0_0_z_xx_0_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xxx_0_x, g_yz_xxx_0_y, g_yz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_0_x[i] = -4.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_y[i] = -4.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xx_0_z[i] = -4.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_y_x_0_0_z_xy_0_x, g_y_x_0_0_z_xy_0_y, g_y_x_0_0_z_xy_0_z, g_yz_xxy_0_x, g_yz_xxy_0_y, g_yz_xxy_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xy_0_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_y_x_0_0_z_xz_0_x, g_y_x_0_0_z_xz_0_y, g_y_x_0_0_z_xz_0_z, g_yz_xxz_0_x, g_yz_xxz_0_y, g_yz_xxz_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xz_0_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_y_x_0_0_z_yy_0_x, g_y_x_0_0_z_yy_0_y, g_y_x_0_0_z_yy_0_z, g_yz_xyy_0_x, g_yz_xyy_0_y, g_yz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yy_0_x[i] = 4.0 * g_yz_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_y[i] = 4.0 * g_yz_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_z[i] = 4.0 * g_yz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_y_x_0_0_z_yz_0_x, g_y_x_0_0_z_yz_0_y, g_y_x_0_0_z_yz_0_z, g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_yz_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_y_x_0_0_z_zz_0_x, g_y_x_0_0_z_zz_0_y, g_y_x_0_0_z_zz_0_z, g_yz_xzz_0_x, g_yz_xzz_0_y, g_yz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_zz_0_x[i] = 4.0 * g_yz_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_y[i] = 4.0 * g_yz_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_z[i] = 4.0 * g_yz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xy_xxy_0_x, g_xy_xxy_0_y, g_xy_xxy_0_z, g_y_y_0_0_x_xx_0_x, g_y_y_0_0_x_xx_0_y, g_y_y_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_0_x[i] = 4.0 * g_xy_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_y[i] = 4.0 * g_xy_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xx_0_z[i] = 4.0 * g_xy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xyy_0_x, g_xy_xyy_0_y, g_xy_xyy_0_z, g_y_y_0_0_x_xy_0_x, g_y_y_0_0_x_xy_0_y, g_y_y_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xy_0_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z, g_y_y_0_0_x_xz_0_x, g_y_y_0_0_x_xz_0_y, g_y_y_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xz_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_yyy_0_x, g_xy_yyy_0_y, g_xy_yyy_0_z, g_y_y_0_0_x_yy_0_x, g_y_y_0_0_x_yy_0_y, g_y_y_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yy_0_x[i] = -4.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_y[i] = -4.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_z[i] = -4.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xy_yyz_0_x, g_xy_yyz_0_y, g_xy_yyz_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_y_y_0_0_x_yz_0_x, g_y_y_0_0_x_yz_0_y, g_y_y_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_yz_0_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xy_yzz_0_x, g_xy_yzz_0_y, g_xy_yzz_0_z, g_y_y_0_0_x_zz_0_x, g_y_y_0_0_x_zz_0_y, g_y_y_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_zz_0_x[i] = 4.0 * g_xy_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_y[i] = 4.0 * g_xy_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_z[i] = 4.0 * g_xy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_y_y_0_0_y_xx_0_x, g_y_y_0_0_y_xx_0_y, g_y_y_0_0_y_xx_0_z, g_yy_xxy_0_x, g_yy_xxy_0_y, g_yy_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_0_x[i] = -2.0 * g_0_xxy_0_x[i] * b_exp + 4.0 * g_yy_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_y[i] = -2.0 * g_0_xxy_0_y[i] * b_exp + 4.0 * g_yy_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xx_0_z[i] = -2.0 * g_0_xxy_0_z[i] * b_exp + 4.0 * g_yy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_y_y_0_0_y_xy_0_x, g_y_y_0_0_y_xy_0_y, g_y_y_0_0_y_xy_0_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_xyy_0_x, g_yy_xyy_0_y, g_yy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xy_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xyy_0_x[i] * b_exp - 2.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xyy_0_y[i] * b_exp - 2.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xyy_0_z[i] * b_exp - 2.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_y_y_0_0_y_xz_0_x, g_y_y_0_0_y_xz_0_y, g_y_y_0_0_y_xz_0_z, g_yy_xyz_0_x, g_yy_xyz_0_y, g_yy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_yy_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_yy_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_yy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yyy_0_x, g_0_yyy_0_y, g_0_yyy_0_z, g_y_y_0_0_y_yy_0_x, g_y_y_0_0_y_yy_0_y, g_y_y_0_0_y_yy_0_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z, g_yy_yyy_0_x, g_yy_yyy_0_y, g_yy_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yy_0_x[i] = 2.0 * g_0_y_0_x[i] - 2.0 * g_0_yyy_0_x[i] * b_exp - 4.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_y[i] = 2.0 * g_0_y_0_y[i] - 2.0 * g_0_yyy_0_y[i] * b_exp - 4.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_z[i] = 2.0 * g_0_y_0_z[i] - 2.0 * g_0_yyy_0_z[i] * b_exp - 4.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_y_y_0_0_y_yz_0_x, g_y_y_0_0_y_yz_0_y, g_y_y_0_0_y_yz_0_z, g_yy_yyz_0_x, g_yy_yyz_0_y, g_yy_yyz_0_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_yz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_yyz_0_x[i] * b_exp - 2.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_yyz_0_y[i] * b_exp - 2.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_yyz_0_z[i] * b_exp - 2.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_y_y_0_0_y_zz_0_x, g_y_y_0_0_y_zz_0_y, g_y_y_0_0_y_zz_0_z, g_yy_yzz_0_x, g_yy_yzz_0_y, g_yy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_zz_0_x[i] = -2.0 * g_0_yzz_0_x[i] * b_exp + 4.0 * g_yy_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_y[i] = -2.0 * g_0_yzz_0_y[i] * b_exp + 4.0 * g_yy_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_z[i] = -2.0 * g_0_yzz_0_z[i] * b_exp + 4.0 * g_yy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_0_x, g_y_y_0_0_z_xx_0_y, g_y_y_0_0_z_xx_0_z, g_yz_xxy_0_x, g_yz_xxy_0_y, g_yz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_0_x[i] = 4.0 * g_yz_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_y[i] = 4.0 * g_yz_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xx_0_z[i] = 4.0 * g_yz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_y_y_0_0_z_xy_0_x, g_y_y_0_0_z_xy_0_y, g_y_y_0_0_z_xy_0_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xyy_0_x, g_yz_xyy_0_y, g_yz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xy_0_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_y_y_0_0_z_xz_0_x, g_y_y_0_0_z_xz_0_y, g_y_y_0_0_z_xz_0_z, g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xz_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_y_y_0_0_z_yy_0_x, g_y_y_0_0_z_yy_0_y, g_y_y_0_0_z_yy_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_yyy_0_x, g_yz_yyy_0_y, g_yz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yy_0_x[i] = -4.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_y[i] = -4.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_z[i] = -4.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_y_y_0_0_z_yz_0_x, g_y_y_0_0_z_yz_0_y, g_y_y_0_0_z_yz_0_z, g_yz_yyz_0_x, g_yz_yyz_0_y, g_yz_yyz_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_yz_0_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_y_y_0_0_z_zz_0_x, g_y_y_0_0_z_zz_0_y, g_y_y_0_0_z_zz_0_z, g_yz_yzz_0_x, g_yz_yzz_0_y, g_yz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_zz_0_x[i] = 4.0 * g_yz_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_y[i] = 4.0 * g_yz_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_z[i] = 4.0 * g_yz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_xy_xxz_0_x, g_xy_xxz_0_y, g_xy_xxz_0_z, g_y_z_0_0_x_xx_0_x, g_y_z_0_0_x_xx_0_y, g_y_z_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_0_x[i] = 4.0 * g_xy_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_y[i] = 4.0 * g_xy_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xx_0_z[i] = 4.0 * g_xy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_xy_xyz_0_x, g_xy_xyz_0_y, g_xy_xyz_0_z, g_y_z_0_0_x_xy_0_x, g_y_z_0_0_x_xy_0_y, g_y_z_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xy_0_x[i] = 4.0 * g_xy_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_y[i] = 4.0 * g_xy_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_z[i] = 4.0 * g_xy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_xzz_0_x, g_xy_xzz_0_y, g_xy_xzz_0_z, g_y_z_0_0_x_xz_0_x, g_y_z_0_0_x_xz_0_y, g_y_z_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xz_0_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_xy_yyz_0_x, g_xy_yyz_0_y, g_xy_yyz_0_z, g_y_z_0_0_x_yy_0_x, g_y_z_0_0_x_yy_0_y, g_y_z_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yy_0_x[i] = 4.0 * g_xy_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_y[i] = 4.0 * g_xy_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_z[i] = 4.0 * g_xy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_yzz_0_x, g_xy_yzz_0_y, g_xy_yzz_0_z, g_y_z_0_0_x_yz_0_x, g_y_z_0_0_x_yz_0_y, g_y_z_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_yz_0_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_zzz_0_x, g_xy_zzz_0_y, g_xy_zzz_0_z, g_y_z_0_0_x_zz_0_x, g_y_z_0_0_x_zz_0_y, g_y_z_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_zz_0_x[i] = -4.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_y[i] = -4.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_z[i] = -4.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_y_z_0_0_y_xx_0_x, g_y_z_0_0_y_xx_0_y, g_y_z_0_0_y_xx_0_z, g_yy_xxz_0_x, g_yy_xxz_0_y, g_yy_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_0_x[i] = -2.0 * g_0_xxz_0_x[i] * b_exp + 4.0 * g_yy_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_y[i] = -2.0 * g_0_xxz_0_y[i] * b_exp + 4.0 * g_yy_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xx_0_z[i] = -2.0 * g_0_xxz_0_z[i] * b_exp + 4.0 * g_yy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_y_z_0_0_y_xy_0_x, g_y_z_0_0_y_xy_0_y, g_y_z_0_0_y_xy_0_z, g_yy_xyz_0_x, g_yy_xyz_0_y, g_yy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xy_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_yy_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_yy_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_yy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_y_z_0_0_y_xz_0_x, g_y_z_0_0_y_xz_0_y, g_y_z_0_0_y_xz_0_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_xzz_0_x, g_yy_xzz_0_y, g_yy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xz_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xzz_0_x[i] * b_exp - 2.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xzz_0_y[i] * b_exp - 2.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xzz_0_z[i] * b_exp - 2.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_y_z_0_0_y_yy_0_x, g_y_z_0_0_y_yy_0_y, g_y_z_0_0_y_yy_0_z, g_yy_yyz_0_x, g_yy_yyz_0_y, g_yy_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yy_0_x[i] = -2.0 * g_0_yyz_0_x[i] * b_exp + 4.0 * g_yy_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_y[i] = -2.0 * g_0_yyz_0_y[i] * b_exp + 4.0 * g_yy_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_z[i] = -2.0 * g_0_yyz_0_z[i] * b_exp + 4.0 * g_yy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_y_z_0_0_y_yz_0_x, g_y_z_0_0_y_yz_0_y, g_y_z_0_0_y_yz_0_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z, g_yy_yzz_0_x, g_yy_yzz_0_y, g_yy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_yz_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_yzz_0_x[i] * b_exp - 2.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_yzz_0_y[i] * b_exp - 2.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_yzz_0_z[i] * b_exp - 2.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_zzz_0_x, g_0_zzz_0_y, g_0_zzz_0_z, g_y_z_0_0_y_zz_0_x, g_y_z_0_0_y_zz_0_y, g_y_z_0_0_y_zz_0_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z, g_yy_zzz_0_x, g_yy_zzz_0_y, g_yy_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_zz_0_x[i] = 2.0 * g_0_z_0_x[i] - 2.0 * g_0_zzz_0_x[i] * b_exp - 4.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_y[i] = 2.0 * g_0_z_0_y[i] - 2.0 * g_0_zzz_0_y[i] * b_exp - 4.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_z[i] = 2.0 * g_0_z_0_z[i] - 2.0 * g_0_zzz_0_z[i] * b_exp - 4.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_0_x, g_y_z_0_0_z_xx_0_y, g_y_z_0_0_z_xx_0_z, g_yz_xxz_0_x, g_yz_xxz_0_y, g_yz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_0_x[i] = 4.0 * g_yz_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_y[i] = 4.0 * g_yz_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xx_0_z[i] = 4.0 * g_yz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_z_0_0_z_xy_0_x, g_y_z_0_0_z_xy_0_y, g_y_z_0_0_z_xy_0_z, g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xy_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_z_0_0_z_xz_0_x, g_y_z_0_0_z_xz_0_y, g_y_z_0_0_z_xz_0_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xzz_0_x, g_yz_xzz_0_y, g_yz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xz_0_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_z_0_0_z_yy_0_x, g_y_z_0_0_z_yy_0_y, g_y_z_0_0_z_yy_0_z, g_yz_yyz_0_x, g_yz_yyz_0_y, g_yz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yy_0_x[i] = 4.0 * g_yz_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_y[i] = 4.0 * g_yz_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_z[i] = 4.0 * g_yz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_z_0_0_z_yz_0_x, g_y_z_0_0_z_yz_0_y, g_y_z_0_0_z_yz_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_yzz_0_x, g_yz_yzz_0_y, g_yz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_yz_0_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_z_0_0_z_zz_0_x, g_y_z_0_0_z_zz_0_y, g_y_z_0_0_z_zz_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_zzz_0_x, g_yz_zzz_0_y, g_yz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_zz_0_x[i] = -4.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_y[i] = -4.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_z[i] = -4.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xxx_0_x, g_xz_xxx_0_y, g_xz_xxx_0_z, g_z_x_0_0_x_xx_0_x, g_z_x_0_0_x_xx_0_y, g_z_x_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_0_x[i] = -4.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_y[i] = -4.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xx_0_z[i] = -4.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xz_xxy_0_x, g_xz_xxy_0_y, g_xz_xxy_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_z_x_0_0_x_xy_0_x, g_z_x_0_0_x_xy_0_y, g_z_x_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xy_0_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xz_xxz_0_x, g_xz_xxz_0_y, g_xz_xxz_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_z_x_0_0_x_xz_0_x, g_z_x_0_0_x_xz_0_y, g_z_x_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xz_0_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xz_xyy_0_x, g_xz_xyy_0_y, g_xz_xyy_0_z, g_z_x_0_0_x_yy_0_x, g_z_x_0_0_x_yy_0_y, g_z_x_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yy_0_x[i] = 4.0 * g_xz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_y[i] = 4.0 * g_xz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_z[i] = 4.0 * g_xz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z, g_z_x_0_0_x_yz_0_x, g_z_x_0_0_x_yz_0_y, g_z_x_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_yz_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xz_xzz_0_x, g_xz_xzz_0_y, g_xz_xzz_0_z, g_z_x_0_0_x_zz_0_x, g_z_x_0_0_x_zz_0_y, g_z_x_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_zz_0_x[i] = 4.0 * g_xz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_y[i] = 4.0 * g_xz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_z[i] = 4.0 * g_xz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xxx_0_x, g_yz_xxx_0_y, g_yz_xxx_0_z, g_z_x_0_0_y_xx_0_x, g_z_x_0_0_y_xx_0_y, g_z_x_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_0_x[i] = -4.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_y[i] = -4.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xx_0_z[i] = -4.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_yz_xxy_0_x, g_yz_xxy_0_y, g_yz_xxy_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_z_x_0_0_y_xy_0_x, g_z_x_0_0_y_xy_0_y, g_z_x_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xy_0_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_yz_xxz_0_x, g_yz_xxz_0_y, g_yz_xxz_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_z_x_0_0_y_xz_0_x, g_z_x_0_0_y_xz_0_y, g_z_x_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xz_0_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_yz_xyy_0_x, g_yz_xyy_0_y, g_yz_xyy_0_z, g_z_x_0_0_y_yy_0_x, g_z_x_0_0_y_yy_0_y, g_z_x_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yy_0_x[i] = 4.0 * g_yz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_y[i] = 4.0 * g_yz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_z[i] = 4.0 * g_yz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z, g_z_x_0_0_y_yz_0_x, g_z_x_0_0_y_yz_0_y, g_z_x_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_yz_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_yz_xzz_0_x, g_yz_xzz_0_y, g_yz_xzz_0_z, g_z_x_0_0_y_zz_0_x, g_z_x_0_0_y_zz_0_y, g_z_x_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_zz_0_x[i] = 4.0 * g_yz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_y[i] = 4.0 * g_yz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_z[i] = 4.0 * g_yz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xxx_0_x, g_0_xxx_0_y, g_0_xxx_0_z, g_z_x_0_0_z_xx_0_x, g_z_x_0_0_z_xx_0_y, g_z_x_0_0_z_xx_0_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_xxx_0_x, g_zz_xxx_0_y, g_zz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_0_x[i] = 2.0 * g_0_x_0_x[i] - 2.0 * g_0_xxx_0_x[i] * b_exp - 4.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_y[i] = 2.0 * g_0_x_0_y[i] - 2.0 * g_0_xxx_0_y[i] * b_exp - 4.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xx_0_z[i] = 2.0 * g_0_x_0_z[i] - 2.0 * g_0_xxx_0_z[i] * b_exp - 4.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_z_x_0_0_z_xy_0_x, g_z_x_0_0_z_xy_0_y, g_z_x_0_0_z_xy_0_z, g_zz_xxy_0_x, g_zz_xxy_0_y, g_zz_xxy_0_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xy_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_xxy_0_x[i] * b_exp - 2.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_xxy_0_y[i] * b_exp - 2.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_xxy_0_z[i] * b_exp - 2.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_z_x_0_0_z_xz_0_x, g_z_x_0_0_z_xz_0_y, g_z_x_0_0_z_xz_0_z, g_zz_xxz_0_x, g_zz_xxz_0_y, g_zz_xxz_0_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_xxz_0_x[i] * b_exp - 2.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_xxz_0_y[i] * b_exp - 2.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_xxz_0_z[i] * b_exp - 2.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_z_x_0_0_z_yy_0_x, g_z_x_0_0_z_yy_0_y, g_z_x_0_0_z_yy_0_z, g_zz_xyy_0_x, g_zz_xyy_0_y, g_zz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yy_0_x[i] = -2.0 * g_0_xyy_0_x[i] * b_exp + 4.0 * g_zz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_y[i] = -2.0 * g_0_xyy_0_y[i] * b_exp + 4.0 * g_zz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_z[i] = -2.0 * g_0_xyy_0_z[i] * b_exp + 4.0 * g_zz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_z_x_0_0_z_yz_0_x, g_z_x_0_0_z_yz_0_y, g_z_x_0_0_z_yz_0_z, g_zz_xyz_0_x, g_zz_xyz_0_y, g_zz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_yz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_zz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_zz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_zz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_z_x_0_0_z_zz_0_x, g_z_x_0_0_z_zz_0_y, g_z_x_0_0_z_zz_0_z, g_zz_xzz_0_x, g_zz_xzz_0_y, g_zz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_zz_0_x[i] = -2.0 * g_0_xzz_0_x[i] * b_exp + 4.0 * g_zz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_y[i] = -2.0 * g_0_xzz_0_y[i] * b_exp + 4.0 * g_zz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_z[i] = -2.0 * g_0_xzz_0_z[i] * b_exp + 4.0 * g_zz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_xz_xxy_0_x, g_xz_xxy_0_y, g_xz_xxy_0_z, g_z_y_0_0_x_xx_0_x, g_z_y_0_0_x_xx_0_y, g_z_y_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_0_x[i] = 4.0 * g_xz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_y[i] = 4.0 * g_xz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xx_0_z[i] = 4.0 * g_xz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xyy_0_x, g_xz_xyy_0_y, g_xz_xyy_0_z, g_z_y_0_0_x_xy_0_x, g_z_y_0_0_x_xy_0_y, g_z_y_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xy_0_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z, g_z_y_0_0_x_xz_0_x, g_z_y_0_0_x_xz_0_y, g_z_y_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xz_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_yyy_0_x, g_xz_yyy_0_y, g_xz_yyy_0_z, g_z_y_0_0_x_yy_0_x, g_z_y_0_0_x_yy_0_y, g_z_y_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yy_0_x[i] = -4.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_y[i] = -4.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_z[i] = -4.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_xz_yyz_0_x, g_xz_yyz_0_y, g_xz_yyz_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_z_y_0_0_x_yz_0_x, g_z_y_0_0_x_yz_0_y, g_z_y_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_yz_0_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_xz_yzz_0_x, g_xz_yzz_0_y, g_xz_yzz_0_z, g_z_y_0_0_x_zz_0_x, g_z_y_0_0_x_zz_0_y, g_z_y_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_zz_0_x[i] = 4.0 * g_xz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_y[i] = 4.0 * g_xz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_z[i] = 4.0 * g_xz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_yz_xxy_0_x, g_yz_xxy_0_y, g_yz_xxy_0_z, g_z_y_0_0_y_xx_0_x, g_z_y_0_0_y_xx_0_y, g_z_y_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_0_x[i] = 4.0 * g_yz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_y[i] = 4.0 * g_yz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xx_0_z[i] = 4.0 * g_yz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xyy_0_x, g_yz_xyy_0_y, g_yz_xyy_0_z, g_z_y_0_0_y_xy_0_x, g_z_y_0_0_y_xy_0_y, g_z_y_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xy_0_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z, g_z_y_0_0_y_xz_0_x, g_z_y_0_0_y_xz_0_y, g_z_y_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xz_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_yyy_0_x, g_yz_yyy_0_y, g_yz_yyy_0_z, g_z_y_0_0_y_yy_0_x, g_z_y_0_0_y_yy_0_y, g_z_y_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yy_0_x[i] = -4.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_y[i] = -4.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_z[i] = -4.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_yz_yyz_0_x, g_yz_yyz_0_y, g_yz_yyz_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_z_y_0_0_y_yz_0_x, g_z_y_0_0_y_yz_0_y, g_z_y_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_yz_0_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_yz_yzz_0_x, g_yz_yzz_0_y, g_yz_yzz_0_z, g_z_y_0_0_y_zz_0_x, g_z_y_0_0_y_zz_0_y, g_z_y_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_zz_0_x[i] = 4.0 * g_yz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_y[i] = 4.0 * g_yz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_z[i] = 4.0 * g_yz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_0_xxy_0_x, g_0_xxy_0_y, g_0_xxy_0_z, g_z_y_0_0_z_xx_0_x, g_z_y_0_0_z_xx_0_y, g_z_y_0_0_z_xx_0_z, g_zz_xxy_0_x, g_zz_xxy_0_y, g_zz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_0_x[i] = -2.0 * g_0_xxy_0_x[i] * b_exp + 4.0 * g_zz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_y[i] = -2.0 * g_0_xxy_0_y[i] * b_exp + 4.0 * g_zz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xx_0_z[i] = -2.0 * g_0_xxy_0_z[i] * b_exp + 4.0 * g_zz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xyy_0_x, g_0_xyy_0_y, g_0_xyy_0_z, g_z_y_0_0_z_xy_0_x, g_z_y_0_0_z_xy_0_y, g_z_y_0_0_z_xy_0_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_xyy_0_x, g_zz_xyy_0_y, g_zz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xy_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xyy_0_x[i] * b_exp - 2.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xyy_0_y[i] * b_exp - 2.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xyy_0_z[i] * b_exp - 2.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_z_y_0_0_z_xz_0_x, g_z_y_0_0_z_xz_0_y, g_z_y_0_0_z_xz_0_z, g_zz_xyz_0_x, g_zz_xyz_0_y, g_zz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xz_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_zz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_zz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_zz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yyy_0_x, g_0_yyy_0_y, g_0_yyy_0_z, g_z_y_0_0_z_yy_0_x, g_z_y_0_0_z_yy_0_y, g_z_y_0_0_z_yy_0_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z, g_zz_yyy_0_x, g_zz_yyy_0_y, g_zz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yy_0_x[i] = 2.0 * g_0_y_0_x[i] - 2.0 * g_0_yyy_0_x[i] * b_exp - 4.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_y[i] = 2.0 * g_0_y_0_y[i] - 2.0 * g_0_yyy_0_y[i] * b_exp - 4.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_z[i] = 2.0 * g_0_y_0_z[i] - 2.0 * g_0_yyy_0_z[i] * b_exp - 4.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_z_y_0_0_z_yz_0_x, g_z_y_0_0_z_yz_0_y, g_z_y_0_0_z_yz_0_z, g_zz_yyz_0_x, g_zz_yyz_0_y, g_zz_yyz_0_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_yz_0_x[i] = g_0_z_0_x[i] - 2.0 * g_0_yyz_0_x[i] * b_exp - 2.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_y[i] = g_0_z_0_y[i] - 2.0 * g_0_yyz_0_y[i] * b_exp - 2.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_z[i] = g_0_z_0_z[i] - 2.0 * g_0_yyz_0_z[i] * b_exp - 2.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_z_y_0_0_z_zz_0_x, g_z_y_0_0_z_zz_0_y, g_z_y_0_0_z_zz_0_z, g_zz_yzz_0_x, g_zz_yzz_0_y, g_zz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_zz_0_x[i] = -2.0 * g_0_yzz_0_x[i] * b_exp + 4.0 * g_zz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_y[i] = -2.0 * g_0_yzz_0_y[i] * b_exp + 4.0 * g_zz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_z[i] = -2.0 * g_0_yzz_0_z[i] * b_exp + 4.0 * g_zz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xz_xxz_0_x, g_xz_xxz_0_y, g_xz_xxz_0_z, g_z_z_0_0_x_xx_0_x, g_z_z_0_0_x_xx_0_y, g_z_z_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_0_x[i] = 4.0 * g_xz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_y[i] = 4.0 * g_xz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xx_0_z[i] = 4.0 * g_xz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xz_xyz_0_x, g_xz_xyz_0_y, g_xz_xyz_0_z, g_z_z_0_0_x_xy_0_x, g_z_z_0_0_x_xy_0_y, g_z_z_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xy_0_x[i] = 4.0 * g_xz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_y[i] = 4.0 * g_xz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_z[i] = 4.0 * g_xz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_xzz_0_x, g_xz_xzz_0_y, g_xz_xzz_0_z, g_z_z_0_0_x_xz_0_x, g_z_z_0_0_x_xz_0_y, g_z_z_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xz_0_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xz_yyz_0_x, g_xz_yyz_0_y, g_xz_yyz_0_z, g_z_z_0_0_x_yy_0_x, g_z_z_0_0_x_yy_0_y, g_z_z_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yy_0_x[i] = 4.0 * g_xz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_y[i] = 4.0 * g_xz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_z[i] = 4.0 * g_xz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_yzz_0_x, g_xz_yzz_0_y, g_xz_yzz_0_z, g_z_z_0_0_x_yz_0_x, g_z_z_0_0_x_yz_0_y, g_z_z_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_yz_0_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_zzz_0_x, g_xz_zzz_0_y, g_xz_zzz_0_z, g_z_z_0_0_x_zz_0_x, g_z_z_0_0_x_zz_0_y, g_z_z_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_zz_0_x[i] = -4.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_y[i] = -4.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_z[i] = -4.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_yz_xxz_0_x, g_yz_xxz_0_y, g_yz_xxz_0_z, g_z_z_0_0_y_xx_0_x, g_z_z_0_0_y_xx_0_y, g_z_z_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_0_x[i] = 4.0 * g_yz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_y[i] = 4.0 * g_yz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xx_0_z[i] = 4.0 * g_yz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_yz_xyz_0_x, g_yz_xyz_0_y, g_yz_xyz_0_z, g_z_z_0_0_y_xy_0_x, g_z_z_0_0_y_xy_0_y, g_z_z_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xy_0_x[i] = 4.0 * g_yz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_y[i] = 4.0 * g_yz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_z[i] = 4.0 * g_yz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_xzz_0_x, g_yz_xzz_0_y, g_yz_xzz_0_z, g_z_z_0_0_y_xz_0_x, g_z_z_0_0_y_xz_0_y, g_z_z_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xz_0_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_yz_yyz_0_x, g_yz_yyz_0_y, g_yz_yyz_0_z, g_z_z_0_0_y_yy_0_x, g_z_z_0_0_y_yy_0_y, g_z_z_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yy_0_x[i] = 4.0 * g_yz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_y[i] = 4.0 * g_yz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_z[i] = 4.0 * g_yz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_yzz_0_x, g_yz_yzz_0_y, g_yz_yzz_0_z, g_z_z_0_0_y_yz_0_x, g_z_z_0_0_y_yz_0_y, g_z_z_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_yz_0_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_zzz_0_x, g_yz_zzz_0_y, g_yz_zzz_0_z, g_z_z_0_0_y_zz_0_x, g_z_z_0_0_y_zz_0_y, g_z_z_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_zz_0_x[i] = -4.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_y[i] = -4.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_z[i] = -4.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_0_xxz_0_x, g_0_xxz_0_y, g_0_xxz_0_z, g_z_z_0_0_z_xx_0_x, g_z_z_0_0_z_xx_0_y, g_z_z_0_0_z_xx_0_z, g_zz_xxz_0_x, g_zz_xxz_0_y, g_zz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_0_x[i] = -2.0 * g_0_xxz_0_x[i] * b_exp + 4.0 * g_zz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_y[i] = -2.0 * g_0_xxz_0_y[i] * b_exp + 4.0 * g_zz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xx_0_z[i] = -2.0 * g_0_xxz_0_z[i] * b_exp + 4.0 * g_zz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_0_xyz_0_x, g_0_xyz_0_y, g_0_xyz_0_z, g_z_z_0_0_z_xy_0_x, g_z_z_0_0_z_xy_0_y, g_z_z_0_0_z_xy_0_z, g_zz_xyz_0_x, g_zz_xyz_0_y, g_zz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xy_0_x[i] = -2.0 * g_0_xyz_0_x[i] * b_exp + 4.0 * g_zz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_y[i] = -2.0 * g_0_xyz_0_y[i] * b_exp + 4.0 * g_zz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_z[i] = -2.0 * g_0_xyz_0_z[i] * b_exp + 4.0 * g_zz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xzz_0_x, g_0_xzz_0_y, g_0_xzz_0_z, g_z_z_0_0_z_xz_0_x, g_z_z_0_0_z_xz_0_y, g_z_z_0_0_z_xz_0_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_xzz_0_x, g_zz_xzz_0_y, g_zz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xz_0_x[i] = g_0_x_0_x[i] - 2.0 * g_0_xzz_0_x[i] * b_exp - 2.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_y[i] = g_0_x_0_y[i] - 2.0 * g_0_xzz_0_y[i] * b_exp - 2.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_z[i] = g_0_x_0_z[i] - 2.0 * g_0_xzz_0_z[i] * b_exp - 2.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_0_yyz_0_x, g_0_yyz_0_y, g_0_yyz_0_z, g_z_z_0_0_z_yy_0_x, g_z_z_0_0_z_yy_0_y, g_z_z_0_0_z_yy_0_z, g_zz_yyz_0_x, g_zz_yyz_0_y, g_zz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yy_0_x[i] = -2.0 * g_0_yyz_0_x[i] * b_exp + 4.0 * g_zz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_y[i] = -2.0 * g_0_yyz_0_y[i] * b_exp + 4.0 * g_zz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_z[i] = -2.0 * g_0_yyz_0_z[i] * b_exp + 4.0 * g_zz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yzz_0_x, g_0_yzz_0_y, g_0_yzz_0_z, g_z_z_0_0_z_yz_0_x, g_z_z_0_0_z_yz_0_y, g_z_z_0_0_z_yz_0_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z, g_zz_yzz_0_x, g_zz_yzz_0_y, g_zz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_yz_0_x[i] = g_0_y_0_x[i] - 2.0 * g_0_yzz_0_x[i] * b_exp - 2.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_y[i] = g_0_y_0_y[i] - 2.0 * g_0_yzz_0_y[i] * b_exp - 2.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_z[i] = g_0_y_0_z[i] - 2.0 * g_0_yzz_0_z[i] * b_exp - 2.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_zzz_0_x, g_0_zzz_0_y, g_0_zzz_0_z, g_z_z_0_0_z_zz_0_x, g_z_z_0_0_z_zz_0_y, g_z_z_0_0_z_zz_0_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z, g_zz_zzz_0_x, g_zz_zzz_0_y, g_zz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_zz_0_x[i] = 2.0 * g_0_z_0_x[i] - 2.0 * g_0_zzz_0_x[i] * b_exp - 4.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_y[i] = 2.0 * g_0_z_0_y[i] - 2.0 * g_0_zzz_0_y[i] * b_exp - 4.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_z[i] = 2.0 * g_0_z_0_z[i] - 2.0 * g_0_zzz_0_z[i] * b_exp - 4.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_zzz_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


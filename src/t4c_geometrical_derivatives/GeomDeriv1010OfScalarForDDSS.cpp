#include "GeomDeriv1010OfScalarForDDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ddss_0(CSimdArray<double>& buffer_1010_ddss,
                     const CSimdArray<double>& buffer_pdps,
                     const CSimdArray<double>& buffer_fdps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ddss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdps

    auto g_x_xx_x_0 = buffer_pdps[0];

    auto g_x_xx_y_0 = buffer_pdps[1];

    auto g_x_xx_z_0 = buffer_pdps[2];

    auto g_x_xy_x_0 = buffer_pdps[3];

    auto g_x_xy_y_0 = buffer_pdps[4];

    auto g_x_xy_z_0 = buffer_pdps[5];

    auto g_x_xz_x_0 = buffer_pdps[6];

    auto g_x_xz_y_0 = buffer_pdps[7];

    auto g_x_xz_z_0 = buffer_pdps[8];

    auto g_x_yy_x_0 = buffer_pdps[9];

    auto g_x_yy_y_0 = buffer_pdps[10];

    auto g_x_yy_z_0 = buffer_pdps[11];

    auto g_x_yz_x_0 = buffer_pdps[12];

    auto g_x_yz_y_0 = buffer_pdps[13];

    auto g_x_yz_z_0 = buffer_pdps[14];

    auto g_x_zz_x_0 = buffer_pdps[15];

    auto g_x_zz_y_0 = buffer_pdps[16];

    auto g_x_zz_z_0 = buffer_pdps[17];

    auto g_y_xx_x_0 = buffer_pdps[18];

    auto g_y_xx_y_0 = buffer_pdps[19];

    auto g_y_xx_z_0 = buffer_pdps[20];

    auto g_y_xy_x_0 = buffer_pdps[21];

    auto g_y_xy_y_0 = buffer_pdps[22];

    auto g_y_xy_z_0 = buffer_pdps[23];

    auto g_y_xz_x_0 = buffer_pdps[24];

    auto g_y_xz_y_0 = buffer_pdps[25];

    auto g_y_xz_z_0 = buffer_pdps[26];

    auto g_y_yy_x_0 = buffer_pdps[27];

    auto g_y_yy_y_0 = buffer_pdps[28];

    auto g_y_yy_z_0 = buffer_pdps[29];

    auto g_y_yz_x_0 = buffer_pdps[30];

    auto g_y_yz_y_0 = buffer_pdps[31];

    auto g_y_yz_z_0 = buffer_pdps[32];

    auto g_y_zz_x_0 = buffer_pdps[33];

    auto g_y_zz_y_0 = buffer_pdps[34];

    auto g_y_zz_z_0 = buffer_pdps[35];

    auto g_z_xx_x_0 = buffer_pdps[36];

    auto g_z_xx_y_0 = buffer_pdps[37];

    auto g_z_xx_z_0 = buffer_pdps[38];

    auto g_z_xy_x_0 = buffer_pdps[39];

    auto g_z_xy_y_0 = buffer_pdps[40];

    auto g_z_xy_z_0 = buffer_pdps[41];

    auto g_z_xz_x_0 = buffer_pdps[42];

    auto g_z_xz_y_0 = buffer_pdps[43];

    auto g_z_xz_z_0 = buffer_pdps[44];

    auto g_z_yy_x_0 = buffer_pdps[45];

    auto g_z_yy_y_0 = buffer_pdps[46];

    auto g_z_yy_z_0 = buffer_pdps[47];

    auto g_z_yz_x_0 = buffer_pdps[48];

    auto g_z_yz_y_0 = buffer_pdps[49];

    auto g_z_yz_z_0 = buffer_pdps[50];

    auto g_z_zz_x_0 = buffer_pdps[51];

    auto g_z_zz_y_0 = buffer_pdps[52];

    auto g_z_zz_z_0 = buffer_pdps[53];

    /// Set up components of auxilary buffer : buffer_fdps

    auto g_xxx_xx_x_0 = buffer_fdps[0];

    auto g_xxx_xx_y_0 = buffer_fdps[1];

    auto g_xxx_xx_z_0 = buffer_fdps[2];

    auto g_xxx_xy_x_0 = buffer_fdps[3];

    auto g_xxx_xy_y_0 = buffer_fdps[4];

    auto g_xxx_xy_z_0 = buffer_fdps[5];

    auto g_xxx_xz_x_0 = buffer_fdps[6];

    auto g_xxx_xz_y_0 = buffer_fdps[7];

    auto g_xxx_xz_z_0 = buffer_fdps[8];

    auto g_xxx_yy_x_0 = buffer_fdps[9];

    auto g_xxx_yy_y_0 = buffer_fdps[10];

    auto g_xxx_yy_z_0 = buffer_fdps[11];

    auto g_xxx_yz_x_0 = buffer_fdps[12];

    auto g_xxx_yz_y_0 = buffer_fdps[13];

    auto g_xxx_yz_z_0 = buffer_fdps[14];

    auto g_xxx_zz_x_0 = buffer_fdps[15];

    auto g_xxx_zz_y_0 = buffer_fdps[16];

    auto g_xxx_zz_z_0 = buffer_fdps[17];

    auto g_xxy_xx_x_0 = buffer_fdps[18];

    auto g_xxy_xx_y_0 = buffer_fdps[19];

    auto g_xxy_xx_z_0 = buffer_fdps[20];

    auto g_xxy_xy_x_0 = buffer_fdps[21];

    auto g_xxy_xy_y_0 = buffer_fdps[22];

    auto g_xxy_xy_z_0 = buffer_fdps[23];

    auto g_xxy_xz_x_0 = buffer_fdps[24];

    auto g_xxy_xz_y_0 = buffer_fdps[25];

    auto g_xxy_xz_z_0 = buffer_fdps[26];

    auto g_xxy_yy_x_0 = buffer_fdps[27];

    auto g_xxy_yy_y_0 = buffer_fdps[28];

    auto g_xxy_yy_z_0 = buffer_fdps[29];

    auto g_xxy_yz_x_0 = buffer_fdps[30];

    auto g_xxy_yz_y_0 = buffer_fdps[31];

    auto g_xxy_yz_z_0 = buffer_fdps[32];

    auto g_xxy_zz_x_0 = buffer_fdps[33];

    auto g_xxy_zz_y_0 = buffer_fdps[34];

    auto g_xxy_zz_z_0 = buffer_fdps[35];

    auto g_xxz_xx_x_0 = buffer_fdps[36];

    auto g_xxz_xx_y_0 = buffer_fdps[37];

    auto g_xxz_xx_z_0 = buffer_fdps[38];

    auto g_xxz_xy_x_0 = buffer_fdps[39];

    auto g_xxz_xy_y_0 = buffer_fdps[40];

    auto g_xxz_xy_z_0 = buffer_fdps[41];

    auto g_xxz_xz_x_0 = buffer_fdps[42];

    auto g_xxz_xz_y_0 = buffer_fdps[43];

    auto g_xxz_xz_z_0 = buffer_fdps[44];

    auto g_xxz_yy_x_0 = buffer_fdps[45];

    auto g_xxz_yy_y_0 = buffer_fdps[46];

    auto g_xxz_yy_z_0 = buffer_fdps[47];

    auto g_xxz_yz_x_0 = buffer_fdps[48];

    auto g_xxz_yz_y_0 = buffer_fdps[49];

    auto g_xxz_yz_z_0 = buffer_fdps[50];

    auto g_xxz_zz_x_0 = buffer_fdps[51];

    auto g_xxz_zz_y_0 = buffer_fdps[52];

    auto g_xxz_zz_z_0 = buffer_fdps[53];

    auto g_xyy_xx_x_0 = buffer_fdps[54];

    auto g_xyy_xx_y_0 = buffer_fdps[55];

    auto g_xyy_xx_z_0 = buffer_fdps[56];

    auto g_xyy_xy_x_0 = buffer_fdps[57];

    auto g_xyy_xy_y_0 = buffer_fdps[58];

    auto g_xyy_xy_z_0 = buffer_fdps[59];

    auto g_xyy_xz_x_0 = buffer_fdps[60];

    auto g_xyy_xz_y_0 = buffer_fdps[61];

    auto g_xyy_xz_z_0 = buffer_fdps[62];

    auto g_xyy_yy_x_0 = buffer_fdps[63];

    auto g_xyy_yy_y_0 = buffer_fdps[64];

    auto g_xyy_yy_z_0 = buffer_fdps[65];

    auto g_xyy_yz_x_0 = buffer_fdps[66];

    auto g_xyy_yz_y_0 = buffer_fdps[67];

    auto g_xyy_yz_z_0 = buffer_fdps[68];

    auto g_xyy_zz_x_0 = buffer_fdps[69];

    auto g_xyy_zz_y_0 = buffer_fdps[70];

    auto g_xyy_zz_z_0 = buffer_fdps[71];

    auto g_xyz_xx_x_0 = buffer_fdps[72];

    auto g_xyz_xx_y_0 = buffer_fdps[73];

    auto g_xyz_xx_z_0 = buffer_fdps[74];

    auto g_xyz_xy_x_0 = buffer_fdps[75];

    auto g_xyz_xy_y_0 = buffer_fdps[76];

    auto g_xyz_xy_z_0 = buffer_fdps[77];

    auto g_xyz_xz_x_0 = buffer_fdps[78];

    auto g_xyz_xz_y_0 = buffer_fdps[79];

    auto g_xyz_xz_z_0 = buffer_fdps[80];

    auto g_xyz_yy_x_0 = buffer_fdps[81];

    auto g_xyz_yy_y_0 = buffer_fdps[82];

    auto g_xyz_yy_z_0 = buffer_fdps[83];

    auto g_xyz_yz_x_0 = buffer_fdps[84];

    auto g_xyz_yz_y_0 = buffer_fdps[85];

    auto g_xyz_yz_z_0 = buffer_fdps[86];

    auto g_xyz_zz_x_0 = buffer_fdps[87];

    auto g_xyz_zz_y_0 = buffer_fdps[88];

    auto g_xyz_zz_z_0 = buffer_fdps[89];

    auto g_xzz_xx_x_0 = buffer_fdps[90];

    auto g_xzz_xx_y_0 = buffer_fdps[91];

    auto g_xzz_xx_z_0 = buffer_fdps[92];

    auto g_xzz_xy_x_0 = buffer_fdps[93];

    auto g_xzz_xy_y_0 = buffer_fdps[94];

    auto g_xzz_xy_z_0 = buffer_fdps[95];

    auto g_xzz_xz_x_0 = buffer_fdps[96];

    auto g_xzz_xz_y_0 = buffer_fdps[97];

    auto g_xzz_xz_z_0 = buffer_fdps[98];

    auto g_xzz_yy_x_0 = buffer_fdps[99];

    auto g_xzz_yy_y_0 = buffer_fdps[100];

    auto g_xzz_yy_z_0 = buffer_fdps[101];

    auto g_xzz_yz_x_0 = buffer_fdps[102];

    auto g_xzz_yz_y_0 = buffer_fdps[103];

    auto g_xzz_yz_z_0 = buffer_fdps[104];

    auto g_xzz_zz_x_0 = buffer_fdps[105];

    auto g_xzz_zz_y_0 = buffer_fdps[106];

    auto g_xzz_zz_z_0 = buffer_fdps[107];

    auto g_yyy_xx_x_0 = buffer_fdps[108];

    auto g_yyy_xx_y_0 = buffer_fdps[109];

    auto g_yyy_xx_z_0 = buffer_fdps[110];

    auto g_yyy_xy_x_0 = buffer_fdps[111];

    auto g_yyy_xy_y_0 = buffer_fdps[112];

    auto g_yyy_xy_z_0 = buffer_fdps[113];

    auto g_yyy_xz_x_0 = buffer_fdps[114];

    auto g_yyy_xz_y_0 = buffer_fdps[115];

    auto g_yyy_xz_z_0 = buffer_fdps[116];

    auto g_yyy_yy_x_0 = buffer_fdps[117];

    auto g_yyy_yy_y_0 = buffer_fdps[118];

    auto g_yyy_yy_z_0 = buffer_fdps[119];

    auto g_yyy_yz_x_0 = buffer_fdps[120];

    auto g_yyy_yz_y_0 = buffer_fdps[121];

    auto g_yyy_yz_z_0 = buffer_fdps[122];

    auto g_yyy_zz_x_0 = buffer_fdps[123];

    auto g_yyy_zz_y_0 = buffer_fdps[124];

    auto g_yyy_zz_z_0 = buffer_fdps[125];

    auto g_yyz_xx_x_0 = buffer_fdps[126];

    auto g_yyz_xx_y_0 = buffer_fdps[127];

    auto g_yyz_xx_z_0 = buffer_fdps[128];

    auto g_yyz_xy_x_0 = buffer_fdps[129];

    auto g_yyz_xy_y_0 = buffer_fdps[130];

    auto g_yyz_xy_z_0 = buffer_fdps[131];

    auto g_yyz_xz_x_0 = buffer_fdps[132];

    auto g_yyz_xz_y_0 = buffer_fdps[133];

    auto g_yyz_xz_z_0 = buffer_fdps[134];

    auto g_yyz_yy_x_0 = buffer_fdps[135];

    auto g_yyz_yy_y_0 = buffer_fdps[136];

    auto g_yyz_yy_z_0 = buffer_fdps[137];

    auto g_yyz_yz_x_0 = buffer_fdps[138];

    auto g_yyz_yz_y_0 = buffer_fdps[139];

    auto g_yyz_yz_z_0 = buffer_fdps[140];

    auto g_yyz_zz_x_0 = buffer_fdps[141];

    auto g_yyz_zz_y_0 = buffer_fdps[142];

    auto g_yyz_zz_z_0 = buffer_fdps[143];

    auto g_yzz_xx_x_0 = buffer_fdps[144];

    auto g_yzz_xx_y_0 = buffer_fdps[145];

    auto g_yzz_xx_z_0 = buffer_fdps[146];

    auto g_yzz_xy_x_0 = buffer_fdps[147];

    auto g_yzz_xy_y_0 = buffer_fdps[148];

    auto g_yzz_xy_z_0 = buffer_fdps[149];

    auto g_yzz_xz_x_0 = buffer_fdps[150];

    auto g_yzz_xz_y_0 = buffer_fdps[151];

    auto g_yzz_xz_z_0 = buffer_fdps[152];

    auto g_yzz_yy_x_0 = buffer_fdps[153];

    auto g_yzz_yy_y_0 = buffer_fdps[154];

    auto g_yzz_yy_z_0 = buffer_fdps[155];

    auto g_yzz_yz_x_0 = buffer_fdps[156];

    auto g_yzz_yz_y_0 = buffer_fdps[157];

    auto g_yzz_yz_z_0 = buffer_fdps[158];

    auto g_yzz_zz_x_0 = buffer_fdps[159];

    auto g_yzz_zz_y_0 = buffer_fdps[160];

    auto g_yzz_zz_z_0 = buffer_fdps[161];

    auto g_zzz_xx_x_0 = buffer_fdps[162];

    auto g_zzz_xx_y_0 = buffer_fdps[163];

    auto g_zzz_xx_z_0 = buffer_fdps[164];

    auto g_zzz_xy_x_0 = buffer_fdps[165];

    auto g_zzz_xy_y_0 = buffer_fdps[166];

    auto g_zzz_xy_z_0 = buffer_fdps[167];

    auto g_zzz_xz_x_0 = buffer_fdps[168];

    auto g_zzz_xz_y_0 = buffer_fdps[169];

    auto g_zzz_xz_z_0 = buffer_fdps[170];

    auto g_zzz_yy_x_0 = buffer_fdps[171];

    auto g_zzz_yy_y_0 = buffer_fdps[172];

    auto g_zzz_yy_z_0 = buffer_fdps[173];

    auto g_zzz_yz_x_0 = buffer_fdps[174];

    auto g_zzz_yz_y_0 = buffer_fdps[175];

    auto g_zzz_yz_z_0 = buffer_fdps[176];

    auto g_zzz_zz_x_0 = buffer_fdps[177];

    auto g_zzz_zz_y_0 = buffer_fdps[178];

    auto g_zzz_zz_z_0 = buffer_fdps[179];

    /// Set up components of integrals buffer : buffer_1010_ddss

    auto g_x_0_x_0_xx_xx_0_0 = buffer_1010_ddss[0];

    auto g_x_0_x_0_xx_xy_0_0 = buffer_1010_ddss[1];

    auto g_x_0_x_0_xx_xz_0_0 = buffer_1010_ddss[2];

    auto g_x_0_x_0_xx_yy_0_0 = buffer_1010_ddss[3];

    auto g_x_0_x_0_xx_yz_0_0 = buffer_1010_ddss[4];

    auto g_x_0_x_0_xx_zz_0_0 = buffer_1010_ddss[5];

    auto g_x_0_x_0_xy_xx_0_0 = buffer_1010_ddss[6];

    auto g_x_0_x_0_xy_xy_0_0 = buffer_1010_ddss[7];

    auto g_x_0_x_0_xy_xz_0_0 = buffer_1010_ddss[8];

    auto g_x_0_x_0_xy_yy_0_0 = buffer_1010_ddss[9];

    auto g_x_0_x_0_xy_yz_0_0 = buffer_1010_ddss[10];

    auto g_x_0_x_0_xy_zz_0_0 = buffer_1010_ddss[11];

    auto g_x_0_x_0_xz_xx_0_0 = buffer_1010_ddss[12];

    auto g_x_0_x_0_xz_xy_0_0 = buffer_1010_ddss[13];

    auto g_x_0_x_0_xz_xz_0_0 = buffer_1010_ddss[14];

    auto g_x_0_x_0_xz_yy_0_0 = buffer_1010_ddss[15];

    auto g_x_0_x_0_xz_yz_0_0 = buffer_1010_ddss[16];

    auto g_x_0_x_0_xz_zz_0_0 = buffer_1010_ddss[17];

    auto g_x_0_x_0_yy_xx_0_0 = buffer_1010_ddss[18];

    auto g_x_0_x_0_yy_xy_0_0 = buffer_1010_ddss[19];

    auto g_x_0_x_0_yy_xz_0_0 = buffer_1010_ddss[20];

    auto g_x_0_x_0_yy_yy_0_0 = buffer_1010_ddss[21];

    auto g_x_0_x_0_yy_yz_0_0 = buffer_1010_ddss[22];

    auto g_x_0_x_0_yy_zz_0_0 = buffer_1010_ddss[23];

    auto g_x_0_x_0_yz_xx_0_0 = buffer_1010_ddss[24];

    auto g_x_0_x_0_yz_xy_0_0 = buffer_1010_ddss[25];

    auto g_x_0_x_0_yz_xz_0_0 = buffer_1010_ddss[26];

    auto g_x_0_x_0_yz_yy_0_0 = buffer_1010_ddss[27];

    auto g_x_0_x_0_yz_yz_0_0 = buffer_1010_ddss[28];

    auto g_x_0_x_0_yz_zz_0_0 = buffer_1010_ddss[29];

    auto g_x_0_x_0_zz_xx_0_0 = buffer_1010_ddss[30];

    auto g_x_0_x_0_zz_xy_0_0 = buffer_1010_ddss[31];

    auto g_x_0_x_0_zz_xz_0_0 = buffer_1010_ddss[32];

    auto g_x_0_x_0_zz_yy_0_0 = buffer_1010_ddss[33];

    auto g_x_0_x_0_zz_yz_0_0 = buffer_1010_ddss[34];

    auto g_x_0_x_0_zz_zz_0_0 = buffer_1010_ddss[35];

    auto g_x_0_y_0_xx_xx_0_0 = buffer_1010_ddss[36];

    auto g_x_0_y_0_xx_xy_0_0 = buffer_1010_ddss[37];

    auto g_x_0_y_0_xx_xz_0_0 = buffer_1010_ddss[38];

    auto g_x_0_y_0_xx_yy_0_0 = buffer_1010_ddss[39];

    auto g_x_0_y_0_xx_yz_0_0 = buffer_1010_ddss[40];

    auto g_x_0_y_0_xx_zz_0_0 = buffer_1010_ddss[41];

    auto g_x_0_y_0_xy_xx_0_0 = buffer_1010_ddss[42];

    auto g_x_0_y_0_xy_xy_0_0 = buffer_1010_ddss[43];

    auto g_x_0_y_0_xy_xz_0_0 = buffer_1010_ddss[44];

    auto g_x_0_y_0_xy_yy_0_0 = buffer_1010_ddss[45];

    auto g_x_0_y_0_xy_yz_0_0 = buffer_1010_ddss[46];

    auto g_x_0_y_0_xy_zz_0_0 = buffer_1010_ddss[47];

    auto g_x_0_y_0_xz_xx_0_0 = buffer_1010_ddss[48];

    auto g_x_0_y_0_xz_xy_0_0 = buffer_1010_ddss[49];

    auto g_x_0_y_0_xz_xz_0_0 = buffer_1010_ddss[50];

    auto g_x_0_y_0_xz_yy_0_0 = buffer_1010_ddss[51];

    auto g_x_0_y_0_xz_yz_0_0 = buffer_1010_ddss[52];

    auto g_x_0_y_0_xz_zz_0_0 = buffer_1010_ddss[53];

    auto g_x_0_y_0_yy_xx_0_0 = buffer_1010_ddss[54];

    auto g_x_0_y_0_yy_xy_0_0 = buffer_1010_ddss[55];

    auto g_x_0_y_0_yy_xz_0_0 = buffer_1010_ddss[56];

    auto g_x_0_y_0_yy_yy_0_0 = buffer_1010_ddss[57];

    auto g_x_0_y_0_yy_yz_0_0 = buffer_1010_ddss[58];

    auto g_x_0_y_0_yy_zz_0_0 = buffer_1010_ddss[59];

    auto g_x_0_y_0_yz_xx_0_0 = buffer_1010_ddss[60];

    auto g_x_0_y_0_yz_xy_0_0 = buffer_1010_ddss[61];

    auto g_x_0_y_0_yz_xz_0_0 = buffer_1010_ddss[62];

    auto g_x_0_y_0_yz_yy_0_0 = buffer_1010_ddss[63];

    auto g_x_0_y_0_yz_yz_0_0 = buffer_1010_ddss[64];

    auto g_x_0_y_0_yz_zz_0_0 = buffer_1010_ddss[65];

    auto g_x_0_y_0_zz_xx_0_0 = buffer_1010_ddss[66];

    auto g_x_0_y_0_zz_xy_0_0 = buffer_1010_ddss[67];

    auto g_x_0_y_0_zz_xz_0_0 = buffer_1010_ddss[68];

    auto g_x_0_y_0_zz_yy_0_0 = buffer_1010_ddss[69];

    auto g_x_0_y_0_zz_yz_0_0 = buffer_1010_ddss[70];

    auto g_x_0_y_0_zz_zz_0_0 = buffer_1010_ddss[71];

    auto g_x_0_z_0_xx_xx_0_0 = buffer_1010_ddss[72];

    auto g_x_0_z_0_xx_xy_0_0 = buffer_1010_ddss[73];

    auto g_x_0_z_0_xx_xz_0_0 = buffer_1010_ddss[74];

    auto g_x_0_z_0_xx_yy_0_0 = buffer_1010_ddss[75];

    auto g_x_0_z_0_xx_yz_0_0 = buffer_1010_ddss[76];

    auto g_x_0_z_0_xx_zz_0_0 = buffer_1010_ddss[77];

    auto g_x_0_z_0_xy_xx_0_0 = buffer_1010_ddss[78];

    auto g_x_0_z_0_xy_xy_0_0 = buffer_1010_ddss[79];

    auto g_x_0_z_0_xy_xz_0_0 = buffer_1010_ddss[80];

    auto g_x_0_z_0_xy_yy_0_0 = buffer_1010_ddss[81];

    auto g_x_0_z_0_xy_yz_0_0 = buffer_1010_ddss[82];

    auto g_x_0_z_0_xy_zz_0_0 = buffer_1010_ddss[83];

    auto g_x_0_z_0_xz_xx_0_0 = buffer_1010_ddss[84];

    auto g_x_0_z_0_xz_xy_0_0 = buffer_1010_ddss[85];

    auto g_x_0_z_0_xz_xz_0_0 = buffer_1010_ddss[86];

    auto g_x_0_z_0_xz_yy_0_0 = buffer_1010_ddss[87];

    auto g_x_0_z_0_xz_yz_0_0 = buffer_1010_ddss[88];

    auto g_x_0_z_0_xz_zz_0_0 = buffer_1010_ddss[89];

    auto g_x_0_z_0_yy_xx_0_0 = buffer_1010_ddss[90];

    auto g_x_0_z_0_yy_xy_0_0 = buffer_1010_ddss[91];

    auto g_x_0_z_0_yy_xz_0_0 = buffer_1010_ddss[92];

    auto g_x_0_z_0_yy_yy_0_0 = buffer_1010_ddss[93];

    auto g_x_0_z_0_yy_yz_0_0 = buffer_1010_ddss[94];

    auto g_x_0_z_0_yy_zz_0_0 = buffer_1010_ddss[95];

    auto g_x_0_z_0_yz_xx_0_0 = buffer_1010_ddss[96];

    auto g_x_0_z_0_yz_xy_0_0 = buffer_1010_ddss[97];

    auto g_x_0_z_0_yz_xz_0_0 = buffer_1010_ddss[98];

    auto g_x_0_z_0_yz_yy_0_0 = buffer_1010_ddss[99];

    auto g_x_0_z_0_yz_yz_0_0 = buffer_1010_ddss[100];

    auto g_x_0_z_0_yz_zz_0_0 = buffer_1010_ddss[101];

    auto g_x_0_z_0_zz_xx_0_0 = buffer_1010_ddss[102];

    auto g_x_0_z_0_zz_xy_0_0 = buffer_1010_ddss[103];

    auto g_x_0_z_0_zz_xz_0_0 = buffer_1010_ddss[104];

    auto g_x_0_z_0_zz_yy_0_0 = buffer_1010_ddss[105];

    auto g_x_0_z_0_zz_yz_0_0 = buffer_1010_ddss[106];

    auto g_x_0_z_0_zz_zz_0_0 = buffer_1010_ddss[107];

    auto g_y_0_x_0_xx_xx_0_0 = buffer_1010_ddss[108];

    auto g_y_0_x_0_xx_xy_0_0 = buffer_1010_ddss[109];

    auto g_y_0_x_0_xx_xz_0_0 = buffer_1010_ddss[110];

    auto g_y_0_x_0_xx_yy_0_0 = buffer_1010_ddss[111];

    auto g_y_0_x_0_xx_yz_0_0 = buffer_1010_ddss[112];

    auto g_y_0_x_0_xx_zz_0_0 = buffer_1010_ddss[113];

    auto g_y_0_x_0_xy_xx_0_0 = buffer_1010_ddss[114];

    auto g_y_0_x_0_xy_xy_0_0 = buffer_1010_ddss[115];

    auto g_y_0_x_0_xy_xz_0_0 = buffer_1010_ddss[116];

    auto g_y_0_x_0_xy_yy_0_0 = buffer_1010_ddss[117];

    auto g_y_0_x_0_xy_yz_0_0 = buffer_1010_ddss[118];

    auto g_y_0_x_0_xy_zz_0_0 = buffer_1010_ddss[119];

    auto g_y_0_x_0_xz_xx_0_0 = buffer_1010_ddss[120];

    auto g_y_0_x_0_xz_xy_0_0 = buffer_1010_ddss[121];

    auto g_y_0_x_0_xz_xz_0_0 = buffer_1010_ddss[122];

    auto g_y_0_x_0_xz_yy_0_0 = buffer_1010_ddss[123];

    auto g_y_0_x_0_xz_yz_0_0 = buffer_1010_ddss[124];

    auto g_y_0_x_0_xz_zz_0_0 = buffer_1010_ddss[125];

    auto g_y_0_x_0_yy_xx_0_0 = buffer_1010_ddss[126];

    auto g_y_0_x_0_yy_xy_0_0 = buffer_1010_ddss[127];

    auto g_y_0_x_0_yy_xz_0_0 = buffer_1010_ddss[128];

    auto g_y_0_x_0_yy_yy_0_0 = buffer_1010_ddss[129];

    auto g_y_0_x_0_yy_yz_0_0 = buffer_1010_ddss[130];

    auto g_y_0_x_0_yy_zz_0_0 = buffer_1010_ddss[131];

    auto g_y_0_x_0_yz_xx_0_0 = buffer_1010_ddss[132];

    auto g_y_0_x_0_yz_xy_0_0 = buffer_1010_ddss[133];

    auto g_y_0_x_0_yz_xz_0_0 = buffer_1010_ddss[134];

    auto g_y_0_x_0_yz_yy_0_0 = buffer_1010_ddss[135];

    auto g_y_0_x_0_yz_yz_0_0 = buffer_1010_ddss[136];

    auto g_y_0_x_0_yz_zz_0_0 = buffer_1010_ddss[137];

    auto g_y_0_x_0_zz_xx_0_0 = buffer_1010_ddss[138];

    auto g_y_0_x_0_zz_xy_0_0 = buffer_1010_ddss[139];

    auto g_y_0_x_0_zz_xz_0_0 = buffer_1010_ddss[140];

    auto g_y_0_x_0_zz_yy_0_0 = buffer_1010_ddss[141];

    auto g_y_0_x_0_zz_yz_0_0 = buffer_1010_ddss[142];

    auto g_y_0_x_0_zz_zz_0_0 = buffer_1010_ddss[143];

    auto g_y_0_y_0_xx_xx_0_0 = buffer_1010_ddss[144];

    auto g_y_0_y_0_xx_xy_0_0 = buffer_1010_ddss[145];

    auto g_y_0_y_0_xx_xz_0_0 = buffer_1010_ddss[146];

    auto g_y_0_y_0_xx_yy_0_0 = buffer_1010_ddss[147];

    auto g_y_0_y_0_xx_yz_0_0 = buffer_1010_ddss[148];

    auto g_y_0_y_0_xx_zz_0_0 = buffer_1010_ddss[149];

    auto g_y_0_y_0_xy_xx_0_0 = buffer_1010_ddss[150];

    auto g_y_0_y_0_xy_xy_0_0 = buffer_1010_ddss[151];

    auto g_y_0_y_0_xy_xz_0_0 = buffer_1010_ddss[152];

    auto g_y_0_y_0_xy_yy_0_0 = buffer_1010_ddss[153];

    auto g_y_0_y_0_xy_yz_0_0 = buffer_1010_ddss[154];

    auto g_y_0_y_0_xy_zz_0_0 = buffer_1010_ddss[155];

    auto g_y_0_y_0_xz_xx_0_0 = buffer_1010_ddss[156];

    auto g_y_0_y_0_xz_xy_0_0 = buffer_1010_ddss[157];

    auto g_y_0_y_0_xz_xz_0_0 = buffer_1010_ddss[158];

    auto g_y_0_y_0_xz_yy_0_0 = buffer_1010_ddss[159];

    auto g_y_0_y_0_xz_yz_0_0 = buffer_1010_ddss[160];

    auto g_y_0_y_0_xz_zz_0_0 = buffer_1010_ddss[161];

    auto g_y_0_y_0_yy_xx_0_0 = buffer_1010_ddss[162];

    auto g_y_0_y_0_yy_xy_0_0 = buffer_1010_ddss[163];

    auto g_y_0_y_0_yy_xz_0_0 = buffer_1010_ddss[164];

    auto g_y_0_y_0_yy_yy_0_0 = buffer_1010_ddss[165];

    auto g_y_0_y_0_yy_yz_0_0 = buffer_1010_ddss[166];

    auto g_y_0_y_0_yy_zz_0_0 = buffer_1010_ddss[167];

    auto g_y_0_y_0_yz_xx_0_0 = buffer_1010_ddss[168];

    auto g_y_0_y_0_yz_xy_0_0 = buffer_1010_ddss[169];

    auto g_y_0_y_0_yz_xz_0_0 = buffer_1010_ddss[170];

    auto g_y_0_y_0_yz_yy_0_0 = buffer_1010_ddss[171];

    auto g_y_0_y_0_yz_yz_0_0 = buffer_1010_ddss[172];

    auto g_y_0_y_0_yz_zz_0_0 = buffer_1010_ddss[173];

    auto g_y_0_y_0_zz_xx_0_0 = buffer_1010_ddss[174];

    auto g_y_0_y_0_zz_xy_0_0 = buffer_1010_ddss[175];

    auto g_y_0_y_0_zz_xz_0_0 = buffer_1010_ddss[176];

    auto g_y_0_y_0_zz_yy_0_0 = buffer_1010_ddss[177];

    auto g_y_0_y_0_zz_yz_0_0 = buffer_1010_ddss[178];

    auto g_y_0_y_0_zz_zz_0_0 = buffer_1010_ddss[179];

    auto g_y_0_z_0_xx_xx_0_0 = buffer_1010_ddss[180];

    auto g_y_0_z_0_xx_xy_0_0 = buffer_1010_ddss[181];

    auto g_y_0_z_0_xx_xz_0_0 = buffer_1010_ddss[182];

    auto g_y_0_z_0_xx_yy_0_0 = buffer_1010_ddss[183];

    auto g_y_0_z_0_xx_yz_0_0 = buffer_1010_ddss[184];

    auto g_y_0_z_0_xx_zz_0_0 = buffer_1010_ddss[185];

    auto g_y_0_z_0_xy_xx_0_0 = buffer_1010_ddss[186];

    auto g_y_0_z_0_xy_xy_0_0 = buffer_1010_ddss[187];

    auto g_y_0_z_0_xy_xz_0_0 = buffer_1010_ddss[188];

    auto g_y_0_z_0_xy_yy_0_0 = buffer_1010_ddss[189];

    auto g_y_0_z_0_xy_yz_0_0 = buffer_1010_ddss[190];

    auto g_y_0_z_0_xy_zz_0_0 = buffer_1010_ddss[191];

    auto g_y_0_z_0_xz_xx_0_0 = buffer_1010_ddss[192];

    auto g_y_0_z_0_xz_xy_0_0 = buffer_1010_ddss[193];

    auto g_y_0_z_0_xz_xz_0_0 = buffer_1010_ddss[194];

    auto g_y_0_z_0_xz_yy_0_0 = buffer_1010_ddss[195];

    auto g_y_0_z_0_xz_yz_0_0 = buffer_1010_ddss[196];

    auto g_y_0_z_0_xz_zz_0_0 = buffer_1010_ddss[197];

    auto g_y_0_z_0_yy_xx_0_0 = buffer_1010_ddss[198];

    auto g_y_0_z_0_yy_xy_0_0 = buffer_1010_ddss[199];

    auto g_y_0_z_0_yy_xz_0_0 = buffer_1010_ddss[200];

    auto g_y_0_z_0_yy_yy_0_0 = buffer_1010_ddss[201];

    auto g_y_0_z_0_yy_yz_0_0 = buffer_1010_ddss[202];

    auto g_y_0_z_0_yy_zz_0_0 = buffer_1010_ddss[203];

    auto g_y_0_z_0_yz_xx_0_0 = buffer_1010_ddss[204];

    auto g_y_0_z_0_yz_xy_0_0 = buffer_1010_ddss[205];

    auto g_y_0_z_0_yz_xz_0_0 = buffer_1010_ddss[206];

    auto g_y_0_z_0_yz_yy_0_0 = buffer_1010_ddss[207];

    auto g_y_0_z_0_yz_yz_0_0 = buffer_1010_ddss[208];

    auto g_y_0_z_0_yz_zz_0_0 = buffer_1010_ddss[209];

    auto g_y_0_z_0_zz_xx_0_0 = buffer_1010_ddss[210];

    auto g_y_0_z_0_zz_xy_0_0 = buffer_1010_ddss[211];

    auto g_y_0_z_0_zz_xz_0_0 = buffer_1010_ddss[212];

    auto g_y_0_z_0_zz_yy_0_0 = buffer_1010_ddss[213];

    auto g_y_0_z_0_zz_yz_0_0 = buffer_1010_ddss[214];

    auto g_y_0_z_0_zz_zz_0_0 = buffer_1010_ddss[215];

    auto g_z_0_x_0_xx_xx_0_0 = buffer_1010_ddss[216];

    auto g_z_0_x_0_xx_xy_0_0 = buffer_1010_ddss[217];

    auto g_z_0_x_0_xx_xz_0_0 = buffer_1010_ddss[218];

    auto g_z_0_x_0_xx_yy_0_0 = buffer_1010_ddss[219];

    auto g_z_0_x_0_xx_yz_0_0 = buffer_1010_ddss[220];

    auto g_z_0_x_0_xx_zz_0_0 = buffer_1010_ddss[221];

    auto g_z_0_x_0_xy_xx_0_0 = buffer_1010_ddss[222];

    auto g_z_0_x_0_xy_xy_0_0 = buffer_1010_ddss[223];

    auto g_z_0_x_0_xy_xz_0_0 = buffer_1010_ddss[224];

    auto g_z_0_x_0_xy_yy_0_0 = buffer_1010_ddss[225];

    auto g_z_0_x_0_xy_yz_0_0 = buffer_1010_ddss[226];

    auto g_z_0_x_0_xy_zz_0_0 = buffer_1010_ddss[227];

    auto g_z_0_x_0_xz_xx_0_0 = buffer_1010_ddss[228];

    auto g_z_0_x_0_xz_xy_0_0 = buffer_1010_ddss[229];

    auto g_z_0_x_0_xz_xz_0_0 = buffer_1010_ddss[230];

    auto g_z_0_x_0_xz_yy_0_0 = buffer_1010_ddss[231];

    auto g_z_0_x_0_xz_yz_0_0 = buffer_1010_ddss[232];

    auto g_z_0_x_0_xz_zz_0_0 = buffer_1010_ddss[233];

    auto g_z_0_x_0_yy_xx_0_0 = buffer_1010_ddss[234];

    auto g_z_0_x_0_yy_xy_0_0 = buffer_1010_ddss[235];

    auto g_z_0_x_0_yy_xz_0_0 = buffer_1010_ddss[236];

    auto g_z_0_x_0_yy_yy_0_0 = buffer_1010_ddss[237];

    auto g_z_0_x_0_yy_yz_0_0 = buffer_1010_ddss[238];

    auto g_z_0_x_0_yy_zz_0_0 = buffer_1010_ddss[239];

    auto g_z_0_x_0_yz_xx_0_0 = buffer_1010_ddss[240];

    auto g_z_0_x_0_yz_xy_0_0 = buffer_1010_ddss[241];

    auto g_z_0_x_0_yz_xz_0_0 = buffer_1010_ddss[242];

    auto g_z_0_x_0_yz_yy_0_0 = buffer_1010_ddss[243];

    auto g_z_0_x_0_yz_yz_0_0 = buffer_1010_ddss[244];

    auto g_z_0_x_0_yz_zz_0_0 = buffer_1010_ddss[245];

    auto g_z_0_x_0_zz_xx_0_0 = buffer_1010_ddss[246];

    auto g_z_0_x_0_zz_xy_0_0 = buffer_1010_ddss[247];

    auto g_z_0_x_0_zz_xz_0_0 = buffer_1010_ddss[248];

    auto g_z_0_x_0_zz_yy_0_0 = buffer_1010_ddss[249];

    auto g_z_0_x_0_zz_yz_0_0 = buffer_1010_ddss[250];

    auto g_z_0_x_0_zz_zz_0_0 = buffer_1010_ddss[251];

    auto g_z_0_y_0_xx_xx_0_0 = buffer_1010_ddss[252];

    auto g_z_0_y_0_xx_xy_0_0 = buffer_1010_ddss[253];

    auto g_z_0_y_0_xx_xz_0_0 = buffer_1010_ddss[254];

    auto g_z_0_y_0_xx_yy_0_0 = buffer_1010_ddss[255];

    auto g_z_0_y_0_xx_yz_0_0 = buffer_1010_ddss[256];

    auto g_z_0_y_0_xx_zz_0_0 = buffer_1010_ddss[257];

    auto g_z_0_y_0_xy_xx_0_0 = buffer_1010_ddss[258];

    auto g_z_0_y_0_xy_xy_0_0 = buffer_1010_ddss[259];

    auto g_z_0_y_0_xy_xz_0_0 = buffer_1010_ddss[260];

    auto g_z_0_y_0_xy_yy_0_0 = buffer_1010_ddss[261];

    auto g_z_0_y_0_xy_yz_0_0 = buffer_1010_ddss[262];

    auto g_z_0_y_0_xy_zz_0_0 = buffer_1010_ddss[263];

    auto g_z_0_y_0_xz_xx_0_0 = buffer_1010_ddss[264];

    auto g_z_0_y_0_xz_xy_0_0 = buffer_1010_ddss[265];

    auto g_z_0_y_0_xz_xz_0_0 = buffer_1010_ddss[266];

    auto g_z_0_y_0_xz_yy_0_0 = buffer_1010_ddss[267];

    auto g_z_0_y_0_xz_yz_0_0 = buffer_1010_ddss[268];

    auto g_z_0_y_0_xz_zz_0_0 = buffer_1010_ddss[269];

    auto g_z_0_y_0_yy_xx_0_0 = buffer_1010_ddss[270];

    auto g_z_0_y_0_yy_xy_0_0 = buffer_1010_ddss[271];

    auto g_z_0_y_0_yy_xz_0_0 = buffer_1010_ddss[272];

    auto g_z_0_y_0_yy_yy_0_0 = buffer_1010_ddss[273];

    auto g_z_0_y_0_yy_yz_0_0 = buffer_1010_ddss[274];

    auto g_z_0_y_0_yy_zz_0_0 = buffer_1010_ddss[275];

    auto g_z_0_y_0_yz_xx_0_0 = buffer_1010_ddss[276];

    auto g_z_0_y_0_yz_xy_0_0 = buffer_1010_ddss[277];

    auto g_z_0_y_0_yz_xz_0_0 = buffer_1010_ddss[278];

    auto g_z_0_y_0_yz_yy_0_0 = buffer_1010_ddss[279];

    auto g_z_0_y_0_yz_yz_0_0 = buffer_1010_ddss[280];

    auto g_z_0_y_0_yz_zz_0_0 = buffer_1010_ddss[281];

    auto g_z_0_y_0_zz_xx_0_0 = buffer_1010_ddss[282];

    auto g_z_0_y_0_zz_xy_0_0 = buffer_1010_ddss[283];

    auto g_z_0_y_0_zz_xz_0_0 = buffer_1010_ddss[284];

    auto g_z_0_y_0_zz_yy_0_0 = buffer_1010_ddss[285];

    auto g_z_0_y_0_zz_yz_0_0 = buffer_1010_ddss[286];

    auto g_z_0_y_0_zz_zz_0_0 = buffer_1010_ddss[287];

    auto g_z_0_z_0_xx_xx_0_0 = buffer_1010_ddss[288];

    auto g_z_0_z_0_xx_xy_0_0 = buffer_1010_ddss[289];

    auto g_z_0_z_0_xx_xz_0_0 = buffer_1010_ddss[290];

    auto g_z_0_z_0_xx_yy_0_0 = buffer_1010_ddss[291];

    auto g_z_0_z_0_xx_yz_0_0 = buffer_1010_ddss[292];

    auto g_z_0_z_0_xx_zz_0_0 = buffer_1010_ddss[293];

    auto g_z_0_z_0_xy_xx_0_0 = buffer_1010_ddss[294];

    auto g_z_0_z_0_xy_xy_0_0 = buffer_1010_ddss[295];

    auto g_z_0_z_0_xy_xz_0_0 = buffer_1010_ddss[296];

    auto g_z_0_z_0_xy_yy_0_0 = buffer_1010_ddss[297];

    auto g_z_0_z_0_xy_yz_0_0 = buffer_1010_ddss[298];

    auto g_z_0_z_0_xy_zz_0_0 = buffer_1010_ddss[299];

    auto g_z_0_z_0_xz_xx_0_0 = buffer_1010_ddss[300];

    auto g_z_0_z_0_xz_xy_0_0 = buffer_1010_ddss[301];

    auto g_z_0_z_0_xz_xz_0_0 = buffer_1010_ddss[302];

    auto g_z_0_z_0_xz_yy_0_0 = buffer_1010_ddss[303];

    auto g_z_0_z_0_xz_yz_0_0 = buffer_1010_ddss[304];

    auto g_z_0_z_0_xz_zz_0_0 = buffer_1010_ddss[305];

    auto g_z_0_z_0_yy_xx_0_0 = buffer_1010_ddss[306];

    auto g_z_0_z_0_yy_xy_0_0 = buffer_1010_ddss[307];

    auto g_z_0_z_0_yy_xz_0_0 = buffer_1010_ddss[308];

    auto g_z_0_z_0_yy_yy_0_0 = buffer_1010_ddss[309];

    auto g_z_0_z_0_yy_yz_0_0 = buffer_1010_ddss[310];

    auto g_z_0_z_0_yy_zz_0_0 = buffer_1010_ddss[311];

    auto g_z_0_z_0_yz_xx_0_0 = buffer_1010_ddss[312];

    auto g_z_0_z_0_yz_xy_0_0 = buffer_1010_ddss[313];

    auto g_z_0_z_0_yz_xz_0_0 = buffer_1010_ddss[314];

    auto g_z_0_z_0_yz_yy_0_0 = buffer_1010_ddss[315];

    auto g_z_0_z_0_yz_yz_0_0 = buffer_1010_ddss[316];

    auto g_z_0_z_0_yz_zz_0_0 = buffer_1010_ddss[317];

    auto g_z_0_z_0_zz_xx_0_0 = buffer_1010_ddss[318];

    auto g_z_0_z_0_zz_xy_0_0 = buffer_1010_ddss[319];

    auto g_z_0_z_0_zz_xz_0_0 = buffer_1010_ddss[320];

    auto g_z_0_z_0_zz_yy_0_0 = buffer_1010_ddss[321];

    auto g_z_0_z_0_zz_yz_0_0 = buffer_1010_ddss[322];

    auto g_z_0_z_0_zz_zz_0_0 = buffer_1010_ddss[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_xx_xx_0_0, g_x_0_x_0_xx_xy_0_0, g_x_0_x_0_xx_xz_0_0, g_x_0_x_0_xx_yy_0_0, g_x_0_x_0_xx_yz_0_0, g_x_0_x_0_xx_zz_0_0, g_x_xx_x_0, g_x_xy_x_0, g_x_xz_x_0, g_x_yy_x_0, g_x_yz_x_0, g_x_zz_x_0, g_xxx_xx_x_0, g_xxx_xy_x_0, g_xxx_xz_x_0, g_xxx_yy_x_0, g_xxx_yz_x_0, g_xxx_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xx_0_0[i] = -4.0 * g_x_xx_x_0[i] * c_exps[i] + 4.0 * g_xxx_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_0[i] = -4.0 * g_x_xy_x_0[i] * c_exps[i] + 4.0 * g_xxx_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_0[i] = -4.0 * g_x_xz_x_0[i] * c_exps[i] + 4.0 * g_xxx_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_0[i] = -4.0 * g_x_yy_x_0[i] * c_exps[i] + 4.0 * g_xxx_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_0[i] = -4.0 * g_x_yz_x_0[i] * c_exps[i] + 4.0 * g_xxx_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_0[i] = -4.0 * g_x_zz_x_0[i] * c_exps[i] + 4.0 * g_xxx_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_xy_xx_0_0, g_x_0_x_0_xy_xy_0_0, g_x_0_x_0_xy_xz_0_0, g_x_0_x_0_xy_yy_0_0, g_x_0_x_0_xy_yz_0_0, g_x_0_x_0_xy_zz_0_0, g_xxy_xx_x_0, g_xxy_xy_x_0, g_xxy_xz_x_0, g_xxy_yy_x_0, g_xxy_yz_x_0, g_xxy_zz_x_0, g_y_xx_x_0, g_y_xy_x_0, g_y_xz_x_0, g_y_yy_x_0, g_y_yz_x_0, g_y_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xx_0_0[i] = -2.0 * g_y_xx_x_0[i] * c_exps[i] + 4.0 * g_xxy_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_0[i] = -2.0 * g_y_xy_x_0[i] * c_exps[i] + 4.0 * g_xxy_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_0[i] = -2.0 * g_y_xz_x_0[i] * c_exps[i] + 4.0 * g_xxy_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_0[i] = -2.0 * g_y_yy_x_0[i] * c_exps[i] + 4.0 * g_xxy_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_0[i] = -2.0 * g_y_yz_x_0[i] * c_exps[i] + 4.0 * g_xxy_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_0[i] = -2.0 * g_y_zz_x_0[i] * c_exps[i] + 4.0 * g_xxy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_xz_xx_0_0, g_x_0_x_0_xz_xy_0_0, g_x_0_x_0_xz_xz_0_0, g_x_0_x_0_xz_yy_0_0, g_x_0_x_0_xz_yz_0_0, g_x_0_x_0_xz_zz_0_0, g_xxz_xx_x_0, g_xxz_xy_x_0, g_xxz_xz_x_0, g_xxz_yy_x_0, g_xxz_yz_x_0, g_xxz_zz_x_0, g_z_xx_x_0, g_z_xy_x_0, g_z_xz_x_0, g_z_yy_x_0, g_z_yz_x_0, g_z_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xx_0_0[i] = -2.0 * g_z_xx_x_0[i] * c_exps[i] + 4.0 * g_xxz_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_0[i] = -2.0 * g_z_xy_x_0[i] * c_exps[i] + 4.0 * g_xxz_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_0[i] = -2.0 * g_z_xz_x_0[i] * c_exps[i] + 4.0 * g_xxz_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_0[i] = -2.0 * g_z_yy_x_0[i] * c_exps[i] + 4.0 * g_xxz_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_0[i] = -2.0 * g_z_yz_x_0[i] * c_exps[i] + 4.0 * g_xxz_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_0[i] = -2.0 * g_z_zz_x_0[i] * c_exps[i] + 4.0 * g_xxz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_yy_xx_0_0, g_x_0_x_0_yy_xy_0_0, g_x_0_x_0_yy_xz_0_0, g_x_0_x_0_yy_yy_0_0, g_x_0_x_0_yy_yz_0_0, g_x_0_x_0_yy_zz_0_0, g_xyy_xx_x_0, g_xyy_xy_x_0, g_xyy_xz_x_0, g_xyy_yy_x_0, g_xyy_yz_x_0, g_xyy_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xx_0_0[i] = 4.0 * g_xyy_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_0[i] = 4.0 * g_xyy_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_0[i] = 4.0 * g_xyy_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_0[i] = 4.0 * g_xyy_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_0[i] = 4.0 * g_xyy_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_0[i] = 4.0 * g_xyy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_yz_xx_0_0, g_x_0_x_0_yz_xy_0_0, g_x_0_x_0_yz_xz_0_0, g_x_0_x_0_yz_yy_0_0, g_x_0_x_0_yz_yz_0_0, g_x_0_x_0_yz_zz_0_0, g_xyz_xx_x_0, g_xyz_xy_x_0, g_xyz_xz_x_0, g_xyz_yy_x_0, g_xyz_yz_x_0, g_xyz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xx_0_0[i] = 4.0 * g_xyz_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_0[i] = 4.0 * g_xyz_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_0[i] = 4.0 * g_xyz_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_0[i] = 4.0 * g_xyz_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_0[i] = 4.0 * g_xyz_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_0[i] = 4.0 * g_xyz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_zz_xx_0_0, g_x_0_x_0_zz_xy_0_0, g_x_0_x_0_zz_xz_0_0, g_x_0_x_0_zz_yy_0_0, g_x_0_x_0_zz_yz_0_0, g_x_0_x_0_zz_zz_0_0, g_xzz_xx_x_0, g_xzz_xy_x_0, g_xzz_xz_x_0, g_xzz_yy_x_0, g_xzz_yz_x_0, g_xzz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xx_0_0[i] = 4.0 * g_xzz_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_0[i] = 4.0 * g_xzz_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_0[i] = 4.0 * g_xzz_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_0[i] = 4.0 * g_xzz_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_0[i] = 4.0 * g_xzz_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_0[i] = 4.0 * g_xzz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_y_0_xx_xx_0_0, g_x_0_y_0_xx_xy_0_0, g_x_0_y_0_xx_xz_0_0, g_x_0_y_0_xx_yy_0_0, g_x_0_y_0_xx_yz_0_0, g_x_0_y_0_xx_zz_0_0, g_x_xx_y_0, g_x_xy_y_0, g_x_xz_y_0, g_x_yy_y_0, g_x_yz_y_0, g_x_zz_y_0, g_xxx_xx_y_0, g_xxx_xy_y_0, g_xxx_xz_y_0, g_xxx_yy_y_0, g_xxx_yz_y_0, g_xxx_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xx_0_0[i] = -4.0 * g_x_xx_y_0[i] * c_exps[i] + 4.0 * g_xxx_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_0[i] = -4.0 * g_x_xy_y_0[i] * c_exps[i] + 4.0 * g_xxx_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_0[i] = -4.0 * g_x_xz_y_0[i] * c_exps[i] + 4.0 * g_xxx_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_0[i] = -4.0 * g_x_yy_y_0[i] * c_exps[i] + 4.0 * g_xxx_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_0[i] = -4.0 * g_x_yz_y_0[i] * c_exps[i] + 4.0 * g_xxx_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_0[i] = -4.0 * g_x_zz_y_0[i] * c_exps[i] + 4.0 * g_xxx_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_y_0_xy_xx_0_0, g_x_0_y_0_xy_xy_0_0, g_x_0_y_0_xy_xz_0_0, g_x_0_y_0_xy_yy_0_0, g_x_0_y_0_xy_yz_0_0, g_x_0_y_0_xy_zz_0_0, g_xxy_xx_y_0, g_xxy_xy_y_0, g_xxy_xz_y_0, g_xxy_yy_y_0, g_xxy_yz_y_0, g_xxy_zz_y_0, g_y_xx_y_0, g_y_xy_y_0, g_y_xz_y_0, g_y_yy_y_0, g_y_yz_y_0, g_y_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xx_0_0[i] = -2.0 * g_y_xx_y_0[i] * c_exps[i] + 4.0 * g_xxy_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_0[i] = -2.0 * g_y_xy_y_0[i] * c_exps[i] + 4.0 * g_xxy_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_0[i] = -2.0 * g_y_xz_y_0[i] * c_exps[i] + 4.0 * g_xxy_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_0[i] = -2.0 * g_y_yy_y_0[i] * c_exps[i] + 4.0 * g_xxy_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_0[i] = -2.0 * g_y_yz_y_0[i] * c_exps[i] + 4.0 * g_xxy_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_0[i] = -2.0 * g_y_zz_y_0[i] * c_exps[i] + 4.0 * g_xxy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_y_0_xz_xx_0_0, g_x_0_y_0_xz_xy_0_0, g_x_0_y_0_xz_xz_0_0, g_x_0_y_0_xz_yy_0_0, g_x_0_y_0_xz_yz_0_0, g_x_0_y_0_xz_zz_0_0, g_xxz_xx_y_0, g_xxz_xy_y_0, g_xxz_xz_y_0, g_xxz_yy_y_0, g_xxz_yz_y_0, g_xxz_zz_y_0, g_z_xx_y_0, g_z_xy_y_0, g_z_xz_y_0, g_z_yy_y_0, g_z_yz_y_0, g_z_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xx_0_0[i] = -2.0 * g_z_xx_y_0[i] * c_exps[i] + 4.0 * g_xxz_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_0[i] = -2.0 * g_z_xy_y_0[i] * c_exps[i] + 4.0 * g_xxz_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_0[i] = -2.0 * g_z_xz_y_0[i] * c_exps[i] + 4.0 * g_xxz_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_0[i] = -2.0 * g_z_yy_y_0[i] * c_exps[i] + 4.0 * g_xxz_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_0[i] = -2.0 * g_z_yz_y_0[i] * c_exps[i] + 4.0 * g_xxz_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_0[i] = -2.0 * g_z_zz_y_0[i] * c_exps[i] + 4.0 * g_xxz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_y_0_yy_xx_0_0, g_x_0_y_0_yy_xy_0_0, g_x_0_y_0_yy_xz_0_0, g_x_0_y_0_yy_yy_0_0, g_x_0_y_0_yy_yz_0_0, g_x_0_y_0_yy_zz_0_0, g_xyy_xx_y_0, g_xyy_xy_y_0, g_xyy_xz_y_0, g_xyy_yy_y_0, g_xyy_yz_y_0, g_xyy_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xx_0_0[i] = 4.0 * g_xyy_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_0[i] = 4.0 * g_xyy_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_0[i] = 4.0 * g_xyy_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_0[i] = 4.0 * g_xyy_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_0[i] = 4.0 * g_xyy_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_0[i] = 4.0 * g_xyy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_y_0_yz_xx_0_0, g_x_0_y_0_yz_xy_0_0, g_x_0_y_0_yz_xz_0_0, g_x_0_y_0_yz_yy_0_0, g_x_0_y_0_yz_yz_0_0, g_x_0_y_0_yz_zz_0_0, g_xyz_xx_y_0, g_xyz_xy_y_0, g_xyz_xz_y_0, g_xyz_yy_y_0, g_xyz_yz_y_0, g_xyz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xx_0_0[i] = 4.0 * g_xyz_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_0[i] = 4.0 * g_xyz_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_0[i] = 4.0 * g_xyz_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_0[i] = 4.0 * g_xyz_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_0[i] = 4.0 * g_xyz_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_0[i] = 4.0 * g_xyz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_y_0_zz_xx_0_0, g_x_0_y_0_zz_xy_0_0, g_x_0_y_0_zz_xz_0_0, g_x_0_y_0_zz_yy_0_0, g_x_0_y_0_zz_yz_0_0, g_x_0_y_0_zz_zz_0_0, g_xzz_xx_y_0, g_xzz_xy_y_0, g_xzz_xz_y_0, g_xzz_yy_y_0, g_xzz_yz_y_0, g_xzz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xx_0_0[i] = 4.0 * g_xzz_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_0[i] = 4.0 * g_xzz_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_0[i] = 4.0 * g_xzz_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_0[i] = 4.0 * g_xzz_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_0[i] = 4.0 * g_xzz_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_0[i] = 4.0 * g_xzz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_z_0_xx_xx_0_0, g_x_0_z_0_xx_xy_0_0, g_x_0_z_0_xx_xz_0_0, g_x_0_z_0_xx_yy_0_0, g_x_0_z_0_xx_yz_0_0, g_x_0_z_0_xx_zz_0_0, g_x_xx_z_0, g_x_xy_z_0, g_x_xz_z_0, g_x_yy_z_0, g_x_yz_z_0, g_x_zz_z_0, g_xxx_xx_z_0, g_xxx_xy_z_0, g_xxx_xz_z_0, g_xxx_yy_z_0, g_xxx_yz_z_0, g_xxx_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xx_0_0[i] = -4.0 * g_x_xx_z_0[i] * c_exps[i] + 4.0 * g_xxx_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_0[i] = -4.0 * g_x_xy_z_0[i] * c_exps[i] + 4.0 * g_xxx_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_0[i] = -4.0 * g_x_xz_z_0[i] * c_exps[i] + 4.0 * g_xxx_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_0[i] = -4.0 * g_x_yy_z_0[i] * c_exps[i] + 4.0 * g_xxx_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_0[i] = -4.0 * g_x_yz_z_0[i] * c_exps[i] + 4.0 * g_xxx_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_0[i] = -4.0 * g_x_zz_z_0[i] * c_exps[i] + 4.0 * g_xxx_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_z_0_xy_xx_0_0, g_x_0_z_0_xy_xy_0_0, g_x_0_z_0_xy_xz_0_0, g_x_0_z_0_xy_yy_0_0, g_x_0_z_0_xy_yz_0_0, g_x_0_z_0_xy_zz_0_0, g_xxy_xx_z_0, g_xxy_xy_z_0, g_xxy_xz_z_0, g_xxy_yy_z_0, g_xxy_yz_z_0, g_xxy_zz_z_0, g_y_xx_z_0, g_y_xy_z_0, g_y_xz_z_0, g_y_yy_z_0, g_y_yz_z_0, g_y_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xx_0_0[i] = -2.0 * g_y_xx_z_0[i] * c_exps[i] + 4.0 * g_xxy_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_0[i] = -2.0 * g_y_xy_z_0[i] * c_exps[i] + 4.0 * g_xxy_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_0[i] = -2.0 * g_y_xz_z_0[i] * c_exps[i] + 4.0 * g_xxy_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_0[i] = -2.0 * g_y_yy_z_0[i] * c_exps[i] + 4.0 * g_xxy_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_0[i] = -2.0 * g_y_yz_z_0[i] * c_exps[i] + 4.0 * g_xxy_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_0[i] = -2.0 * g_y_zz_z_0[i] * c_exps[i] + 4.0 * g_xxy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_z_0_xz_xx_0_0, g_x_0_z_0_xz_xy_0_0, g_x_0_z_0_xz_xz_0_0, g_x_0_z_0_xz_yy_0_0, g_x_0_z_0_xz_yz_0_0, g_x_0_z_0_xz_zz_0_0, g_xxz_xx_z_0, g_xxz_xy_z_0, g_xxz_xz_z_0, g_xxz_yy_z_0, g_xxz_yz_z_0, g_xxz_zz_z_0, g_z_xx_z_0, g_z_xy_z_0, g_z_xz_z_0, g_z_yy_z_0, g_z_yz_z_0, g_z_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xx_0_0[i] = -2.0 * g_z_xx_z_0[i] * c_exps[i] + 4.0 * g_xxz_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_0[i] = -2.0 * g_z_xy_z_0[i] * c_exps[i] + 4.0 * g_xxz_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_0[i] = -2.0 * g_z_xz_z_0[i] * c_exps[i] + 4.0 * g_xxz_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_0[i] = -2.0 * g_z_yy_z_0[i] * c_exps[i] + 4.0 * g_xxz_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_0[i] = -2.0 * g_z_yz_z_0[i] * c_exps[i] + 4.0 * g_xxz_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_0[i] = -2.0 * g_z_zz_z_0[i] * c_exps[i] + 4.0 * g_xxz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_z_0_yy_xx_0_0, g_x_0_z_0_yy_xy_0_0, g_x_0_z_0_yy_xz_0_0, g_x_0_z_0_yy_yy_0_0, g_x_0_z_0_yy_yz_0_0, g_x_0_z_0_yy_zz_0_0, g_xyy_xx_z_0, g_xyy_xy_z_0, g_xyy_xz_z_0, g_xyy_yy_z_0, g_xyy_yz_z_0, g_xyy_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xx_0_0[i] = 4.0 * g_xyy_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_0[i] = 4.0 * g_xyy_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_0[i] = 4.0 * g_xyy_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_0[i] = 4.0 * g_xyy_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_0[i] = 4.0 * g_xyy_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_0[i] = 4.0 * g_xyy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_z_0_yz_xx_0_0, g_x_0_z_0_yz_xy_0_0, g_x_0_z_0_yz_xz_0_0, g_x_0_z_0_yz_yy_0_0, g_x_0_z_0_yz_yz_0_0, g_x_0_z_0_yz_zz_0_0, g_xyz_xx_z_0, g_xyz_xy_z_0, g_xyz_xz_z_0, g_xyz_yy_z_0, g_xyz_yz_z_0, g_xyz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xx_0_0[i] = 4.0 * g_xyz_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_0[i] = 4.0 * g_xyz_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_0[i] = 4.0 * g_xyz_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_0[i] = 4.0 * g_xyz_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_0[i] = 4.0 * g_xyz_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_0[i] = 4.0 * g_xyz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_z_0_zz_xx_0_0, g_x_0_z_0_zz_xy_0_0, g_x_0_z_0_zz_xz_0_0, g_x_0_z_0_zz_yy_0_0, g_x_0_z_0_zz_yz_0_0, g_x_0_z_0_zz_zz_0_0, g_xzz_xx_z_0, g_xzz_xy_z_0, g_xzz_xz_z_0, g_xzz_yy_z_0, g_xzz_yz_z_0, g_xzz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xx_0_0[i] = 4.0 * g_xzz_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_0[i] = 4.0 * g_xzz_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_0[i] = 4.0 * g_xzz_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_0[i] = 4.0 * g_xzz_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_0[i] = 4.0 * g_xzz_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_0[i] = 4.0 * g_xzz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xxy_xx_x_0, g_xxy_xy_x_0, g_xxy_xz_x_0, g_xxy_yy_x_0, g_xxy_yz_x_0, g_xxy_zz_x_0, g_y_0_x_0_xx_xx_0_0, g_y_0_x_0_xx_xy_0_0, g_y_0_x_0_xx_xz_0_0, g_y_0_x_0_xx_yy_0_0, g_y_0_x_0_xx_yz_0_0, g_y_0_x_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xx_0_0[i] = 4.0 * g_xxy_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_0[i] = 4.0 * g_xxy_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_0[i] = 4.0 * g_xxy_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_0[i] = 4.0 * g_xxy_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_0[i] = 4.0 * g_xxy_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_0[i] = 4.0 * g_xxy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_xx_x_0, g_x_xy_x_0, g_x_xz_x_0, g_x_yy_x_0, g_x_yz_x_0, g_x_zz_x_0, g_xyy_xx_x_0, g_xyy_xy_x_0, g_xyy_xz_x_0, g_xyy_yy_x_0, g_xyy_yz_x_0, g_xyy_zz_x_0, g_y_0_x_0_xy_xx_0_0, g_y_0_x_0_xy_xy_0_0, g_y_0_x_0_xy_xz_0_0, g_y_0_x_0_xy_yy_0_0, g_y_0_x_0_xy_yz_0_0, g_y_0_x_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xx_0_0[i] = -2.0 * g_x_xx_x_0[i] * c_exps[i] + 4.0 * g_xyy_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_0[i] = -2.0 * g_x_xy_x_0[i] * c_exps[i] + 4.0 * g_xyy_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_0[i] = -2.0 * g_x_xz_x_0[i] * c_exps[i] + 4.0 * g_xyy_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_0[i] = -2.0 * g_x_yy_x_0[i] * c_exps[i] + 4.0 * g_xyy_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_0[i] = -2.0 * g_x_yz_x_0[i] * c_exps[i] + 4.0 * g_xyy_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_0[i] = -2.0 * g_x_zz_x_0[i] * c_exps[i] + 4.0 * g_xyy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xyz_xx_x_0, g_xyz_xy_x_0, g_xyz_xz_x_0, g_xyz_yy_x_0, g_xyz_yz_x_0, g_xyz_zz_x_0, g_y_0_x_0_xz_xx_0_0, g_y_0_x_0_xz_xy_0_0, g_y_0_x_0_xz_xz_0_0, g_y_0_x_0_xz_yy_0_0, g_y_0_x_0_xz_yz_0_0, g_y_0_x_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xx_0_0[i] = 4.0 * g_xyz_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_0[i] = 4.0 * g_xyz_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_0[i] = 4.0 * g_xyz_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_0[i] = 4.0 * g_xyz_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_0[i] = 4.0 * g_xyz_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_0[i] = 4.0 * g_xyz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_0_x_0_yy_xx_0_0, g_y_0_x_0_yy_xy_0_0, g_y_0_x_0_yy_xz_0_0, g_y_0_x_0_yy_yy_0_0, g_y_0_x_0_yy_yz_0_0, g_y_0_x_0_yy_zz_0_0, g_y_xx_x_0, g_y_xy_x_0, g_y_xz_x_0, g_y_yy_x_0, g_y_yz_x_0, g_y_zz_x_0, g_yyy_xx_x_0, g_yyy_xy_x_0, g_yyy_xz_x_0, g_yyy_yy_x_0, g_yyy_yz_x_0, g_yyy_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xx_0_0[i] = -4.0 * g_y_xx_x_0[i] * c_exps[i] + 4.0 * g_yyy_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_0[i] = -4.0 * g_y_xy_x_0[i] * c_exps[i] + 4.0 * g_yyy_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_0[i] = -4.0 * g_y_xz_x_0[i] * c_exps[i] + 4.0 * g_yyy_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_0[i] = -4.0 * g_y_yy_x_0[i] * c_exps[i] + 4.0 * g_yyy_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_0[i] = -4.0 * g_y_yz_x_0[i] * c_exps[i] + 4.0 * g_yyy_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_0[i] = -4.0 * g_y_zz_x_0[i] * c_exps[i] + 4.0 * g_yyy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_0_x_0_yz_xx_0_0, g_y_0_x_0_yz_xy_0_0, g_y_0_x_0_yz_xz_0_0, g_y_0_x_0_yz_yy_0_0, g_y_0_x_0_yz_yz_0_0, g_y_0_x_0_yz_zz_0_0, g_yyz_xx_x_0, g_yyz_xy_x_0, g_yyz_xz_x_0, g_yyz_yy_x_0, g_yyz_yz_x_0, g_yyz_zz_x_0, g_z_xx_x_0, g_z_xy_x_0, g_z_xz_x_0, g_z_yy_x_0, g_z_yz_x_0, g_z_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xx_0_0[i] = -2.0 * g_z_xx_x_0[i] * c_exps[i] + 4.0 * g_yyz_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_0[i] = -2.0 * g_z_xy_x_0[i] * c_exps[i] + 4.0 * g_yyz_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_0[i] = -2.0 * g_z_xz_x_0[i] * c_exps[i] + 4.0 * g_yyz_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_0[i] = -2.0 * g_z_yy_x_0[i] * c_exps[i] + 4.0 * g_yyz_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_0[i] = -2.0 * g_z_yz_x_0[i] * c_exps[i] + 4.0 * g_yyz_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_0[i] = -2.0 * g_z_zz_x_0[i] * c_exps[i] + 4.0 * g_yyz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_0_x_0_zz_xx_0_0, g_y_0_x_0_zz_xy_0_0, g_y_0_x_0_zz_xz_0_0, g_y_0_x_0_zz_yy_0_0, g_y_0_x_0_zz_yz_0_0, g_y_0_x_0_zz_zz_0_0, g_yzz_xx_x_0, g_yzz_xy_x_0, g_yzz_xz_x_0, g_yzz_yy_x_0, g_yzz_yz_x_0, g_yzz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xx_0_0[i] = 4.0 * g_yzz_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_0[i] = 4.0 * g_yzz_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_0[i] = 4.0 * g_yzz_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_0[i] = 4.0 * g_yzz_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_0[i] = 4.0 * g_yzz_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_0[i] = 4.0 * g_yzz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xxy_xx_y_0, g_xxy_xy_y_0, g_xxy_xz_y_0, g_xxy_yy_y_0, g_xxy_yz_y_0, g_xxy_zz_y_0, g_y_0_y_0_xx_xx_0_0, g_y_0_y_0_xx_xy_0_0, g_y_0_y_0_xx_xz_0_0, g_y_0_y_0_xx_yy_0_0, g_y_0_y_0_xx_yz_0_0, g_y_0_y_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xx_0_0[i] = 4.0 * g_xxy_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_0[i] = 4.0 * g_xxy_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_0[i] = 4.0 * g_xxy_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_0[i] = 4.0 * g_xxy_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_0[i] = 4.0 * g_xxy_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_0[i] = 4.0 * g_xxy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_xx_y_0, g_x_xy_y_0, g_x_xz_y_0, g_x_yy_y_0, g_x_yz_y_0, g_x_zz_y_0, g_xyy_xx_y_0, g_xyy_xy_y_0, g_xyy_xz_y_0, g_xyy_yy_y_0, g_xyy_yz_y_0, g_xyy_zz_y_0, g_y_0_y_0_xy_xx_0_0, g_y_0_y_0_xy_xy_0_0, g_y_0_y_0_xy_xz_0_0, g_y_0_y_0_xy_yy_0_0, g_y_0_y_0_xy_yz_0_0, g_y_0_y_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xx_0_0[i] = -2.0 * g_x_xx_y_0[i] * c_exps[i] + 4.0 * g_xyy_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_0[i] = -2.0 * g_x_xy_y_0[i] * c_exps[i] + 4.0 * g_xyy_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_0[i] = -2.0 * g_x_xz_y_0[i] * c_exps[i] + 4.0 * g_xyy_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_0[i] = -2.0 * g_x_yy_y_0[i] * c_exps[i] + 4.0 * g_xyy_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_0[i] = -2.0 * g_x_yz_y_0[i] * c_exps[i] + 4.0 * g_xyy_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_0[i] = -2.0 * g_x_zz_y_0[i] * c_exps[i] + 4.0 * g_xyy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xyz_xx_y_0, g_xyz_xy_y_0, g_xyz_xz_y_0, g_xyz_yy_y_0, g_xyz_yz_y_0, g_xyz_zz_y_0, g_y_0_y_0_xz_xx_0_0, g_y_0_y_0_xz_xy_0_0, g_y_0_y_0_xz_xz_0_0, g_y_0_y_0_xz_yy_0_0, g_y_0_y_0_xz_yz_0_0, g_y_0_y_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xx_0_0[i] = 4.0 * g_xyz_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_0[i] = 4.0 * g_xyz_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_0[i] = 4.0 * g_xyz_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_0[i] = 4.0 * g_xyz_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_0[i] = 4.0 * g_xyz_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_0[i] = 4.0 * g_xyz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_y_0_yy_xx_0_0, g_y_0_y_0_yy_xy_0_0, g_y_0_y_0_yy_xz_0_0, g_y_0_y_0_yy_yy_0_0, g_y_0_y_0_yy_yz_0_0, g_y_0_y_0_yy_zz_0_0, g_y_xx_y_0, g_y_xy_y_0, g_y_xz_y_0, g_y_yy_y_0, g_y_yz_y_0, g_y_zz_y_0, g_yyy_xx_y_0, g_yyy_xy_y_0, g_yyy_xz_y_0, g_yyy_yy_y_0, g_yyy_yz_y_0, g_yyy_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xx_0_0[i] = -4.0 * g_y_xx_y_0[i] * c_exps[i] + 4.0 * g_yyy_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_0[i] = -4.0 * g_y_xy_y_0[i] * c_exps[i] + 4.0 * g_yyy_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_0[i] = -4.0 * g_y_xz_y_0[i] * c_exps[i] + 4.0 * g_yyy_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_0[i] = -4.0 * g_y_yy_y_0[i] * c_exps[i] + 4.0 * g_yyy_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_0[i] = -4.0 * g_y_yz_y_0[i] * c_exps[i] + 4.0 * g_yyy_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_0[i] = -4.0 * g_y_zz_y_0[i] * c_exps[i] + 4.0 * g_yyy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_y_0_yz_xx_0_0, g_y_0_y_0_yz_xy_0_0, g_y_0_y_0_yz_xz_0_0, g_y_0_y_0_yz_yy_0_0, g_y_0_y_0_yz_yz_0_0, g_y_0_y_0_yz_zz_0_0, g_yyz_xx_y_0, g_yyz_xy_y_0, g_yyz_xz_y_0, g_yyz_yy_y_0, g_yyz_yz_y_0, g_yyz_zz_y_0, g_z_xx_y_0, g_z_xy_y_0, g_z_xz_y_0, g_z_yy_y_0, g_z_yz_y_0, g_z_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xx_0_0[i] = -2.0 * g_z_xx_y_0[i] * c_exps[i] + 4.0 * g_yyz_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_0[i] = -2.0 * g_z_xy_y_0[i] * c_exps[i] + 4.0 * g_yyz_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_0[i] = -2.0 * g_z_xz_y_0[i] * c_exps[i] + 4.0 * g_yyz_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_0[i] = -2.0 * g_z_yy_y_0[i] * c_exps[i] + 4.0 * g_yyz_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_0[i] = -2.0 * g_z_yz_y_0[i] * c_exps[i] + 4.0 * g_yyz_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_0[i] = -2.0 * g_z_zz_y_0[i] * c_exps[i] + 4.0 * g_yyz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_y_0_zz_xx_0_0, g_y_0_y_0_zz_xy_0_0, g_y_0_y_0_zz_xz_0_0, g_y_0_y_0_zz_yy_0_0, g_y_0_y_0_zz_yz_0_0, g_y_0_y_0_zz_zz_0_0, g_yzz_xx_y_0, g_yzz_xy_y_0, g_yzz_xz_y_0, g_yzz_yy_y_0, g_yzz_yz_y_0, g_yzz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xx_0_0[i] = 4.0 * g_yzz_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_0[i] = 4.0 * g_yzz_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_0[i] = 4.0 * g_yzz_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_0[i] = 4.0 * g_yzz_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_0[i] = 4.0 * g_yzz_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_0[i] = 4.0 * g_yzz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xxy_xx_z_0, g_xxy_xy_z_0, g_xxy_xz_z_0, g_xxy_yy_z_0, g_xxy_yz_z_0, g_xxy_zz_z_0, g_y_0_z_0_xx_xx_0_0, g_y_0_z_0_xx_xy_0_0, g_y_0_z_0_xx_xz_0_0, g_y_0_z_0_xx_yy_0_0, g_y_0_z_0_xx_yz_0_0, g_y_0_z_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xx_0_0[i] = 4.0 * g_xxy_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_0[i] = 4.0 * g_xxy_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_0[i] = 4.0 * g_xxy_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_0[i] = 4.0 * g_xxy_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_0[i] = 4.0 * g_xxy_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_0[i] = 4.0 * g_xxy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_xx_z_0, g_x_xy_z_0, g_x_xz_z_0, g_x_yy_z_0, g_x_yz_z_0, g_x_zz_z_0, g_xyy_xx_z_0, g_xyy_xy_z_0, g_xyy_xz_z_0, g_xyy_yy_z_0, g_xyy_yz_z_0, g_xyy_zz_z_0, g_y_0_z_0_xy_xx_0_0, g_y_0_z_0_xy_xy_0_0, g_y_0_z_0_xy_xz_0_0, g_y_0_z_0_xy_yy_0_0, g_y_0_z_0_xy_yz_0_0, g_y_0_z_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xx_0_0[i] = -2.0 * g_x_xx_z_0[i] * c_exps[i] + 4.0 * g_xyy_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_0[i] = -2.0 * g_x_xy_z_0[i] * c_exps[i] + 4.0 * g_xyy_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_0[i] = -2.0 * g_x_xz_z_0[i] * c_exps[i] + 4.0 * g_xyy_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_0[i] = -2.0 * g_x_yy_z_0[i] * c_exps[i] + 4.0 * g_xyy_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_0[i] = -2.0 * g_x_yz_z_0[i] * c_exps[i] + 4.0 * g_xyy_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_0[i] = -2.0 * g_x_zz_z_0[i] * c_exps[i] + 4.0 * g_xyy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xyz_xx_z_0, g_xyz_xy_z_0, g_xyz_xz_z_0, g_xyz_yy_z_0, g_xyz_yz_z_0, g_xyz_zz_z_0, g_y_0_z_0_xz_xx_0_0, g_y_0_z_0_xz_xy_0_0, g_y_0_z_0_xz_xz_0_0, g_y_0_z_0_xz_yy_0_0, g_y_0_z_0_xz_yz_0_0, g_y_0_z_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xx_0_0[i] = 4.0 * g_xyz_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_0[i] = 4.0 * g_xyz_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_0[i] = 4.0 * g_xyz_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_0[i] = 4.0 * g_xyz_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_0[i] = 4.0 * g_xyz_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_0[i] = 4.0 * g_xyz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_z_0_yy_xx_0_0, g_y_0_z_0_yy_xy_0_0, g_y_0_z_0_yy_xz_0_0, g_y_0_z_0_yy_yy_0_0, g_y_0_z_0_yy_yz_0_0, g_y_0_z_0_yy_zz_0_0, g_y_xx_z_0, g_y_xy_z_0, g_y_xz_z_0, g_y_yy_z_0, g_y_yz_z_0, g_y_zz_z_0, g_yyy_xx_z_0, g_yyy_xy_z_0, g_yyy_xz_z_0, g_yyy_yy_z_0, g_yyy_yz_z_0, g_yyy_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xx_0_0[i] = -4.0 * g_y_xx_z_0[i] * c_exps[i] + 4.0 * g_yyy_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_0[i] = -4.0 * g_y_xy_z_0[i] * c_exps[i] + 4.0 * g_yyy_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_0[i] = -4.0 * g_y_xz_z_0[i] * c_exps[i] + 4.0 * g_yyy_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_0[i] = -4.0 * g_y_yy_z_0[i] * c_exps[i] + 4.0 * g_yyy_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_0[i] = -4.0 * g_y_yz_z_0[i] * c_exps[i] + 4.0 * g_yyy_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_0[i] = -4.0 * g_y_zz_z_0[i] * c_exps[i] + 4.0 * g_yyy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_z_0_yz_xx_0_0, g_y_0_z_0_yz_xy_0_0, g_y_0_z_0_yz_xz_0_0, g_y_0_z_0_yz_yy_0_0, g_y_0_z_0_yz_yz_0_0, g_y_0_z_0_yz_zz_0_0, g_yyz_xx_z_0, g_yyz_xy_z_0, g_yyz_xz_z_0, g_yyz_yy_z_0, g_yyz_yz_z_0, g_yyz_zz_z_0, g_z_xx_z_0, g_z_xy_z_0, g_z_xz_z_0, g_z_yy_z_0, g_z_yz_z_0, g_z_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xx_0_0[i] = -2.0 * g_z_xx_z_0[i] * c_exps[i] + 4.0 * g_yyz_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_0[i] = -2.0 * g_z_xy_z_0[i] * c_exps[i] + 4.0 * g_yyz_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_0[i] = -2.0 * g_z_xz_z_0[i] * c_exps[i] + 4.0 * g_yyz_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_0[i] = -2.0 * g_z_yy_z_0[i] * c_exps[i] + 4.0 * g_yyz_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_0[i] = -2.0 * g_z_yz_z_0[i] * c_exps[i] + 4.0 * g_yyz_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_0[i] = -2.0 * g_z_zz_z_0[i] * c_exps[i] + 4.0 * g_yyz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_z_0_zz_xx_0_0, g_y_0_z_0_zz_xy_0_0, g_y_0_z_0_zz_xz_0_0, g_y_0_z_0_zz_yy_0_0, g_y_0_z_0_zz_yz_0_0, g_y_0_z_0_zz_zz_0_0, g_yzz_xx_z_0, g_yzz_xy_z_0, g_yzz_xz_z_0, g_yzz_yy_z_0, g_yzz_yz_z_0, g_yzz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xx_0_0[i] = 4.0 * g_yzz_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_0[i] = 4.0 * g_yzz_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_0[i] = 4.0 * g_yzz_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_0[i] = 4.0 * g_yzz_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_0[i] = 4.0 * g_yzz_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_0[i] = 4.0 * g_yzz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xxz_xx_x_0, g_xxz_xy_x_0, g_xxz_xz_x_0, g_xxz_yy_x_0, g_xxz_yz_x_0, g_xxz_zz_x_0, g_z_0_x_0_xx_xx_0_0, g_z_0_x_0_xx_xy_0_0, g_z_0_x_0_xx_xz_0_0, g_z_0_x_0_xx_yy_0_0, g_z_0_x_0_xx_yz_0_0, g_z_0_x_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xx_0_0[i] = 4.0 * g_xxz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_0[i] = 4.0 * g_xxz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_0[i] = 4.0 * g_xxz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_0[i] = 4.0 * g_xxz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_0[i] = 4.0 * g_xxz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_0[i] = 4.0 * g_xxz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xyz_xx_x_0, g_xyz_xy_x_0, g_xyz_xz_x_0, g_xyz_yy_x_0, g_xyz_yz_x_0, g_xyz_zz_x_0, g_z_0_x_0_xy_xx_0_0, g_z_0_x_0_xy_xy_0_0, g_z_0_x_0_xy_xz_0_0, g_z_0_x_0_xy_yy_0_0, g_z_0_x_0_xy_yz_0_0, g_z_0_x_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xx_0_0[i] = 4.0 * g_xyz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_0[i] = 4.0 * g_xyz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_0[i] = 4.0 * g_xyz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_0[i] = 4.0 * g_xyz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_0[i] = 4.0 * g_xyz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_0[i] = 4.0 * g_xyz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_xx_x_0, g_x_xy_x_0, g_x_xz_x_0, g_x_yy_x_0, g_x_yz_x_0, g_x_zz_x_0, g_xzz_xx_x_0, g_xzz_xy_x_0, g_xzz_xz_x_0, g_xzz_yy_x_0, g_xzz_yz_x_0, g_xzz_zz_x_0, g_z_0_x_0_xz_xx_0_0, g_z_0_x_0_xz_xy_0_0, g_z_0_x_0_xz_xz_0_0, g_z_0_x_0_xz_yy_0_0, g_z_0_x_0_xz_yz_0_0, g_z_0_x_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xx_0_0[i] = -2.0 * g_x_xx_x_0[i] * c_exps[i] + 4.0 * g_xzz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_0[i] = -2.0 * g_x_xy_x_0[i] * c_exps[i] + 4.0 * g_xzz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_0[i] = -2.0 * g_x_xz_x_0[i] * c_exps[i] + 4.0 * g_xzz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_0[i] = -2.0 * g_x_yy_x_0[i] * c_exps[i] + 4.0 * g_xzz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_0[i] = -2.0 * g_x_yz_x_0[i] * c_exps[i] + 4.0 * g_xzz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_0[i] = -2.0 * g_x_zz_x_0[i] * c_exps[i] + 4.0 * g_xzz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_yyz_xx_x_0, g_yyz_xy_x_0, g_yyz_xz_x_0, g_yyz_yy_x_0, g_yyz_yz_x_0, g_yyz_zz_x_0, g_z_0_x_0_yy_xx_0_0, g_z_0_x_0_yy_xy_0_0, g_z_0_x_0_yy_xz_0_0, g_z_0_x_0_yy_yy_0_0, g_z_0_x_0_yy_yz_0_0, g_z_0_x_0_yy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xx_0_0[i] = 4.0 * g_yyz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_0[i] = 4.0 * g_yyz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_0[i] = 4.0 * g_yyz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_0[i] = 4.0 * g_yyz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_0[i] = 4.0 * g_yyz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_0[i] = 4.0 * g_yyz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_y_xx_x_0, g_y_xy_x_0, g_y_xz_x_0, g_y_yy_x_0, g_y_yz_x_0, g_y_zz_x_0, g_yzz_xx_x_0, g_yzz_xy_x_0, g_yzz_xz_x_0, g_yzz_yy_x_0, g_yzz_yz_x_0, g_yzz_zz_x_0, g_z_0_x_0_yz_xx_0_0, g_z_0_x_0_yz_xy_0_0, g_z_0_x_0_yz_xz_0_0, g_z_0_x_0_yz_yy_0_0, g_z_0_x_0_yz_yz_0_0, g_z_0_x_0_yz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xx_0_0[i] = -2.0 * g_y_xx_x_0[i] * c_exps[i] + 4.0 * g_yzz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_0[i] = -2.0 * g_y_xy_x_0[i] * c_exps[i] + 4.0 * g_yzz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_0[i] = -2.0 * g_y_xz_x_0[i] * c_exps[i] + 4.0 * g_yzz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_0[i] = -2.0 * g_y_yy_x_0[i] * c_exps[i] + 4.0 * g_yzz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_0[i] = -2.0 * g_y_yz_x_0[i] * c_exps[i] + 4.0 * g_yzz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_0[i] = -2.0 * g_y_zz_x_0[i] * c_exps[i] + 4.0 * g_yzz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_0_x_0_zz_xx_0_0, g_z_0_x_0_zz_xy_0_0, g_z_0_x_0_zz_xz_0_0, g_z_0_x_0_zz_yy_0_0, g_z_0_x_0_zz_yz_0_0, g_z_0_x_0_zz_zz_0_0, g_z_xx_x_0, g_z_xy_x_0, g_z_xz_x_0, g_z_yy_x_0, g_z_yz_x_0, g_z_zz_x_0, g_zzz_xx_x_0, g_zzz_xy_x_0, g_zzz_xz_x_0, g_zzz_yy_x_0, g_zzz_yz_x_0, g_zzz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xx_0_0[i] = -4.0 * g_z_xx_x_0[i] * c_exps[i] + 4.0 * g_zzz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_0[i] = -4.0 * g_z_xy_x_0[i] * c_exps[i] + 4.0 * g_zzz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_0[i] = -4.0 * g_z_xz_x_0[i] * c_exps[i] + 4.0 * g_zzz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_0[i] = -4.0 * g_z_yy_x_0[i] * c_exps[i] + 4.0 * g_zzz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_0[i] = -4.0 * g_z_yz_x_0[i] * c_exps[i] + 4.0 * g_zzz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_0[i] = -4.0 * g_z_zz_x_0[i] * c_exps[i] + 4.0 * g_zzz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xxz_xx_y_0, g_xxz_xy_y_0, g_xxz_xz_y_0, g_xxz_yy_y_0, g_xxz_yz_y_0, g_xxz_zz_y_0, g_z_0_y_0_xx_xx_0_0, g_z_0_y_0_xx_xy_0_0, g_z_0_y_0_xx_xz_0_0, g_z_0_y_0_xx_yy_0_0, g_z_0_y_0_xx_yz_0_0, g_z_0_y_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xx_0_0[i] = 4.0 * g_xxz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_0[i] = 4.0 * g_xxz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_0[i] = 4.0 * g_xxz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_0[i] = 4.0 * g_xxz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_0[i] = 4.0 * g_xxz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_0[i] = 4.0 * g_xxz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xyz_xx_y_0, g_xyz_xy_y_0, g_xyz_xz_y_0, g_xyz_yy_y_0, g_xyz_yz_y_0, g_xyz_zz_y_0, g_z_0_y_0_xy_xx_0_0, g_z_0_y_0_xy_xy_0_0, g_z_0_y_0_xy_xz_0_0, g_z_0_y_0_xy_yy_0_0, g_z_0_y_0_xy_yz_0_0, g_z_0_y_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xx_0_0[i] = 4.0 * g_xyz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_0[i] = 4.0 * g_xyz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_0[i] = 4.0 * g_xyz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_0[i] = 4.0 * g_xyz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_0[i] = 4.0 * g_xyz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_0[i] = 4.0 * g_xyz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_xx_y_0, g_x_xy_y_0, g_x_xz_y_0, g_x_yy_y_0, g_x_yz_y_0, g_x_zz_y_0, g_xzz_xx_y_0, g_xzz_xy_y_0, g_xzz_xz_y_0, g_xzz_yy_y_0, g_xzz_yz_y_0, g_xzz_zz_y_0, g_z_0_y_0_xz_xx_0_0, g_z_0_y_0_xz_xy_0_0, g_z_0_y_0_xz_xz_0_0, g_z_0_y_0_xz_yy_0_0, g_z_0_y_0_xz_yz_0_0, g_z_0_y_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xx_0_0[i] = -2.0 * g_x_xx_y_0[i] * c_exps[i] + 4.0 * g_xzz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_0[i] = -2.0 * g_x_xy_y_0[i] * c_exps[i] + 4.0 * g_xzz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_0[i] = -2.0 * g_x_xz_y_0[i] * c_exps[i] + 4.0 * g_xzz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_0[i] = -2.0 * g_x_yy_y_0[i] * c_exps[i] + 4.0 * g_xzz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_0[i] = -2.0 * g_x_yz_y_0[i] * c_exps[i] + 4.0 * g_xzz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_0[i] = -2.0 * g_x_zz_y_0[i] * c_exps[i] + 4.0 * g_xzz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_yyz_xx_y_0, g_yyz_xy_y_0, g_yyz_xz_y_0, g_yyz_yy_y_0, g_yyz_yz_y_0, g_yyz_zz_y_0, g_z_0_y_0_yy_xx_0_0, g_z_0_y_0_yy_xy_0_0, g_z_0_y_0_yy_xz_0_0, g_z_0_y_0_yy_yy_0_0, g_z_0_y_0_yy_yz_0_0, g_z_0_y_0_yy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xx_0_0[i] = 4.0 * g_yyz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_0[i] = 4.0 * g_yyz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_0[i] = 4.0 * g_yyz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_0[i] = 4.0 * g_yyz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_0[i] = 4.0 * g_yyz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_0[i] = 4.0 * g_yyz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_xx_y_0, g_y_xy_y_0, g_y_xz_y_0, g_y_yy_y_0, g_y_yz_y_0, g_y_zz_y_0, g_yzz_xx_y_0, g_yzz_xy_y_0, g_yzz_xz_y_0, g_yzz_yy_y_0, g_yzz_yz_y_0, g_yzz_zz_y_0, g_z_0_y_0_yz_xx_0_0, g_z_0_y_0_yz_xy_0_0, g_z_0_y_0_yz_xz_0_0, g_z_0_y_0_yz_yy_0_0, g_z_0_y_0_yz_yz_0_0, g_z_0_y_0_yz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xx_0_0[i] = -2.0 * g_y_xx_y_0[i] * c_exps[i] + 4.0 * g_yzz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_0[i] = -2.0 * g_y_xy_y_0[i] * c_exps[i] + 4.0 * g_yzz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_0[i] = -2.0 * g_y_xz_y_0[i] * c_exps[i] + 4.0 * g_yzz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_0[i] = -2.0 * g_y_yy_y_0[i] * c_exps[i] + 4.0 * g_yzz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_0[i] = -2.0 * g_y_yz_y_0[i] * c_exps[i] + 4.0 * g_yzz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_0[i] = -2.0 * g_y_zz_y_0[i] * c_exps[i] + 4.0 * g_yzz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_0_y_0_zz_xx_0_0, g_z_0_y_0_zz_xy_0_0, g_z_0_y_0_zz_xz_0_0, g_z_0_y_0_zz_yy_0_0, g_z_0_y_0_zz_yz_0_0, g_z_0_y_0_zz_zz_0_0, g_z_xx_y_0, g_z_xy_y_0, g_z_xz_y_0, g_z_yy_y_0, g_z_yz_y_0, g_z_zz_y_0, g_zzz_xx_y_0, g_zzz_xy_y_0, g_zzz_xz_y_0, g_zzz_yy_y_0, g_zzz_yz_y_0, g_zzz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xx_0_0[i] = -4.0 * g_z_xx_y_0[i] * c_exps[i] + 4.0 * g_zzz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_0[i] = -4.0 * g_z_xy_y_0[i] * c_exps[i] + 4.0 * g_zzz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_0[i] = -4.0 * g_z_xz_y_0[i] * c_exps[i] + 4.0 * g_zzz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_0[i] = -4.0 * g_z_yy_y_0[i] * c_exps[i] + 4.0 * g_zzz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_0[i] = -4.0 * g_z_yz_y_0[i] * c_exps[i] + 4.0 * g_zzz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_0[i] = -4.0 * g_z_zz_y_0[i] * c_exps[i] + 4.0 * g_zzz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xxz_xx_z_0, g_xxz_xy_z_0, g_xxz_xz_z_0, g_xxz_yy_z_0, g_xxz_yz_z_0, g_xxz_zz_z_0, g_z_0_z_0_xx_xx_0_0, g_z_0_z_0_xx_xy_0_0, g_z_0_z_0_xx_xz_0_0, g_z_0_z_0_xx_yy_0_0, g_z_0_z_0_xx_yz_0_0, g_z_0_z_0_xx_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xx_0_0[i] = 4.0 * g_xxz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_0[i] = 4.0 * g_xxz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_0[i] = 4.0 * g_xxz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_0[i] = 4.0 * g_xxz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_0[i] = 4.0 * g_xxz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_0[i] = 4.0 * g_xxz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xyz_xx_z_0, g_xyz_xy_z_0, g_xyz_xz_z_0, g_xyz_yy_z_0, g_xyz_yz_z_0, g_xyz_zz_z_0, g_z_0_z_0_xy_xx_0_0, g_z_0_z_0_xy_xy_0_0, g_z_0_z_0_xy_xz_0_0, g_z_0_z_0_xy_yy_0_0, g_z_0_z_0_xy_yz_0_0, g_z_0_z_0_xy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xx_0_0[i] = 4.0 * g_xyz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_0[i] = 4.0 * g_xyz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_0[i] = 4.0 * g_xyz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_0[i] = 4.0 * g_xyz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_0[i] = 4.0 * g_xyz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_0[i] = 4.0 * g_xyz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_xx_z_0, g_x_xy_z_0, g_x_xz_z_0, g_x_yy_z_0, g_x_yz_z_0, g_x_zz_z_0, g_xzz_xx_z_0, g_xzz_xy_z_0, g_xzz_xz_z_0, g_xzz_yy_z_0, g_xzz_yz_z_0, g_xzz_zz_z_0, g_z_0_z_0_xz_xx_0_0, g_z_0_z_0_xz_xy_0_0, g_z_0_z_0_xz_xz_0_0, g_z_0_z_0_xz_yy_0_0, g_z_0_z_0_xz_yz_0_0, g_z_0_z_0_xz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xx_0_0[i] = -2.0 * g_x_xx_z_0[i] * c_exps[i] + 4.0 * g_xzz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_0[i] = -2.0 * g_x_xy_z_0[i] * c_exps[i] + 4.0 * g_xzz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_0[i] = -2.0 * g_x_xz_z_0[i] * c_exps[i] + 4.0 * g_xzz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_0[i] = -2.0 * g_x_yy_z_0[i] * c_exps[i] + 4.0 * g_xzz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_0[i] = -2.0 * g_x_yz_z_0[i] * c_exps[i] + 4.0 * g_xzz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_0[i] = -2.0 * g_x_zz_z_0[i] * c_exps[i] + 4.0 * g_xzz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_yyz_xx_z_0, g_yyz_xy_z_0, g_yyz_xz_z_0, g_yyz_yy_z_0, g_yyz_yz_z_0, g_yyz_zz_z_0, g_z_0_z_0_yy_xx_0_0, g_z_0_z_0_yy_xy_0_0, g_z_0_z_0_yy_xz_0_0, g_z_0_z_0_yy_yy_0_0, g_z_0_z_0_yy_yz_0_0, g_z_0_z_0_yy_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xx_0_0[i] = 4.0 * g_yyz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_0[i] = 4.0 * g_yyz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_0[i] = 4.0 * g_yyz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_0[i] = 4.0 * g_yyz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_0[i] = 4.0 * g_yyz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_0[i] = 4.0 * g_yyz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_xx_z_0, g_y_xy_z_0, g_y_xz_z_0, g_y_yy_z_0, g_y_yz_z_0, g_y_zz_z_0, g_yzz_xx_z_0, g_yzz_xy_z_0, g_yzz_xz_z_0, g_yzz_yy_z_0, g_yzz_yz_z_0, g_yzz_zz_z_0, g_z_0_z_0_yz_xx_0_0, g_z_0_z_0_yz_xy_0_0, g_z_0_z_0_yz_xz_0_0, g_z_0_z_0_yz_yy_0_0, g_z_0_z_0_yz_yz_0_0, g_z_0_z_0_yz_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xx_0_0[i] = -2.0 * g_y_xx_z_0[i] * c_exps[i] + 4.0 * g_yzz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_0[i] = -2.0 * g_y_xy_z_0[i] * c_exps[i] + 4.0 * g_yzz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_0[i] = -2.0 * g_y_xz_z_0[i] * c_exps[i] + 4.0 * g_yzz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_0[i] = -2.0 * g_y_yy_z_0[i] * c_exps[i] + 4.0 * g_yzz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_0[i] = -2.0 * g_y_yz_z_0[i] * c_exps[i] + 4.0 * g_yzz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_0[i] = -2.0 * g_y_zz_z_0[i] * c_exps[i] + 4.0 * g_yzz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_0_z_0_zz_xx_0_0, g_z_0_z_0_zz_xy_0_0, g_z_0_z_0_zz_xz_0_0, g_z_0_z_0_zz_yy_0_0, g_z_0_z_0_zz_yz_0_0, g_z_0_z_0_zz_zz_0_0, g_z_xx_z_0, g_z_xy_z_0, g_z_xz_z_0, g_z_yy_z_0, g_z_yz_z_0, g_z_zz_z_0, g_zzz_xx_z_0, g_zzz_xy_z_0, g_zzz_xz_z_0, g_zzz_yy_z_0, g_zzz_yz_z_0, g_zzz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xx_0_0[i] = -4.0 * g_z_xx_z_0[i] * c_exps[i] + 4.0 * g_zzz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_0[i] = -4.0 * g_z_xy_z_0[i] * c_exps[i] + 4.0 * g_zzz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_0[i] = -4.0 * g_z_xz_z_0[i] * c_exps[i] + 4.0 * g_zzz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_0[i] = -4.0 * g_z_yy_z_0[i] * c_exps[i] + 4.0 * g_zzz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_0[i] = -4.0 * g_z_yz_z_0[i] * c_exps[i] + 4.0 * g_zzz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_0[i] = -4.0 * g_z_zz_z_0[i] * c_exps[i] + 4.0 * g_zzz_zz_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

